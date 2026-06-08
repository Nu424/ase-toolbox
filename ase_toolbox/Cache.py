"""
# Cache.py
ASEの計算結果や構造を再利用するための軽量キャッシュ機構。

## 関数・クラス一覧
- stable_json_dumps(): オブジェクトを安定したJSON文字列へ変換する
- stable_key(): オブジェクトからsha256キーを作成する
- atoms_key_data() / atoms_fingerprint(): ase.Atomsをキー化する
- KeyBuilder: 条件を段階的に追加してキーを作成する
- JsonValueStore / AtomsStore / PathStore: 保存形式の定義
- CacheStore: 保存・取得・get_or_computeを扱う

## 全体像
- CacheStore がメイン。キャッシュの読み書きの入口となる
- CacheStore に渡す store（StoreSpec）が、実際のファイル読み書きを担当する
  - キーに対応するディレクトリ内へペイロードを save / load する
- JsonValueStore / AtomsStore / PathStore は store の具体実装（通常は直接 save/load しない）

## キャッシュの流れ
- 前提: キー・バリュー形式で管理する
- キーの作り方
  - 条件を dict（Mapping）にまとめ、`stable_key()` で SHA256 ハッシュ（64文字）にする
  - `stable_key()` は辞書のキー順や float 表現のぶれを抑えるため、内部で JSON 正規化してからハッシュ化する
  - 構造は `atoms_fingerprint(atoms)` や `KeyBuilder` でキー構成に含める
- CacheStore の set/get/get_or_compute に渡す key
  - 条件 dict（Mapping）: 内部で `stable_key()` してディレクトリ名にする
  - 64文字の sha256 文字列: そのままディレクトリ名として使う（再ハッシュしない）
  - 注意: Atoms など任意の Python オブジェクトをそのまま渡すわけではない
- 書き込み（CacheStore.set）
  - namespace 配下の一時ディレクトリへ StoreSpec.save() と meta.json を書く
  - 完了後 `os.replace()` で原子的に差し替える（途中失敗で壊れたキャッシュが残りにくい）
- 読み込み（CacheStore.get）
  - キーに対応するディレクトリを特定し、StoreSpec.load() で値を返す
- よく使う入口は get_or_compute（キャッシュがあれば get、なければ計算して set）

## ディレクトリ構造の例
base_dir と namespace は CacheStore 作成時に指定する（例: base_dir=".cache/ase_toolbox", namespace="energy"）。

```
.cache/ase_toolbox/
  formation_energy/          # namespace（JsonValueStore）
    a1b2c3...（64文字）/
      meta.json              # キー情報・作成時刻・store種別など
      payload.json           # 実データ（エネルギーなど）
    d4e5f6.../
      meta.json
      payload.json
  optimized_structures/      # namespace（AtomsStore）
    a1b2c3.../
      meta.json
      structure.traj         # デフォルト。extxyz 等に変更可
  structure_paths/           # namespace（PathStore）
    a1b2c3.../
      meta.json
      path.json              # 外部ファイルへの参照パス
```
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field, is_dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import tempfile
from typing import Any, Generic, Protocol, TypeVar

import numpy as np
from ase import Atoms
from ase.io import read, write

try:
    from .util import sanitize_atoms_for_xyz_write
except ImportError:  # pragma: no cover - 直接実行時の保険
    sanitize_atoms_for_xyz_write = None


T = TypeVar("T")

# ----------
# ---定数
# ----------
_CACHE_SCHEMA_VERSION = 1
_SAFE_COMPONENT_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


# ----------
# ---キー生成（内部ヘルパー）
# ----------
def _round_float(value: float, decimals: int | None) -> float:
    """
    浮動小数点を指定桁数で丸める。

    Args:
        value (float): 丸め対象の値。
        decimals (int | None): 小数点以下の桁数。None の場合は丸めない。

    Returns:
        float: 丸め後の値。
    """
    if decimals is None:
        return float(value)
    return round(float(value), decimals)


def _to_jsonable(obj: Any, *, decimals: int | None = None) -> Any:
    """
    オブジェクトを JSON 化しやすい型へ再帰的に変換する。

    辞書のキーは文字列化してソートし、set は安定した順序のリストへ変換する。
    これにより、同じ内容なら常に同じ JSON 表現になる。

    Args:
        obj: 変換対象のオブジェクト。
        decimals (int | None, optional): float を丸める桁数。デフォルトは None。

    Returns:
        Any: JSON 互換の値。

    Raises:
        TypeError: 変換できない型が含まれる場合。
    """
    # --- スカラー型 ---
    if obj is None or isinstance(obj, (str, int, bool)):
        return obj
    if isinstance(obj, float):
        return _round_float(obj, decimals)

    # --- Path / NumPy ---
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, np.generic):
        return _to_jsonable(obj.item(), decimals=decimals)
    if isinstance(obj, np.ndarray):
        return _to_jsonable(obj.tolist(), decimals=decimals)

    # --- dataclass ---
    if is_dataclass(obj) and not isinstance(obj, type):
        return _to_jsonable(asdict(obj), decimals=decimals)

    # --- 辞書（キーをソートして順序を固定） ---
    if isinstance(obj, Mapping):
        return {
            str(key): _to_jsonable(value, decimals=decimals)
            for key, value in sorted(obj.items(), key=lambda item: str(item[0]))
        }

    # --- シーケンス ---
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(value, decimals=decimals) for value in obj]

    # --- 集合（内容でソートして順序を固定） ---
    if isinstance(obj, set):
        values = [_to_jsonable(value, decimals=decimals) for value in obj]
        return sorted(values, key=lambda value: stable_json_dumps(value))

    raise TypeError(f"JSONに変換できない型です: {type(obj).__name__}")


def stable_json_dumps(obj: Any, *, decimals: int | None = None) -> str:
    """
    オブジェクトを、順序と表現が安定した JSON 文字列へ変換する。

    キャッシュキー生成の前処理として使う。辞書のキー順や float の表現揺れを
    抑えるため、内部で正規化してから `json.dumps()` する。

    Args:
        obj: JSON 互換値、Path、NumPy 値、dataclass など。
        decimals (int | None, optional): float を丸める桁数。デフォルトは None。

    Returns:
        str: キー生成に使える正規化済み JSON 文字列。

    Examples:
        >>> stable_json_dumps({"b": 1, "a": 2})
        '{"a":2,"b":1}'
    """
    jsonable = _to_jsonable(obj, decimals=decimals)
    return json.dumps(
        jsonable,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def stable_key(obj: Any, *, decimals: int | None = None) -> str:
    """
    オブジェクトから sha256 ハッシュキーを作成する。オブジェクトをJSON化してから、sha256ハッシュ化する。

    Args:
        obj: キー化したい条件・設定など。
        decimals (int | None, optional): float を丸める桁数。デフォルトは None。

    Returns:
        str: 64 文字の sha256 ハッシュ（16進小文字）。

    Examples:
        >>> len(stable_key({"task": "test", "value": 1.0}))
        64
    """
    payload = stable_json_dumps(obj, decimals=decimals).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _constraint_key_data(constraint: Any) -> Any:
    """
    ASE の制約オブジェクト(constraint)をキー化用の辞書へ変換する。

    `todict()` が使える場合はそれを優先し、なければクラス名と index 属性から
    最低限の情報を抽出する。

    Args:
        constraint: ASE の Constraint オブジェクト。

    Returns:
        Any: `_to_jsonable()` に渡せる辞書相当の値。
    """
    # ----constraintに対して.todict()できるならそれを返す
    if hasattr(constraint, "todict"):
        try:
            return _to_jsonable(constraint.todict())
        except Exception:
            pass

    # ----それができない場合はクラス名と index 属性から最低限の情報を抽出する
    data: dict[str, Any] = {"name": type(constraint).__name__}
    if hasattr(constraint, "index"):
        try:
            data["index"] = _to_jsonable(getattr(constraint, "index"))
        except Exception:
            data["index"] = str(getattr(constraint, "index"))
    return data


def atoms_key_data(
    atoms: Atoms,
    *,
    decimals: int = 8,
    include_constraints: bool = True,
    include_info: bool = False,
    info_keys: list[str] | tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """
    ase.Atoms を、キー生成に使える辞書へ変換する。

    元素種・座標・セル・周期境界条件を基本として含める。calculator は含めない。

    Args:
        atoms (ase.Atoms): キー化したい構造。
        decimals (int, optional): 座標・セルを丸める桁数。デフォルトは 8。
        include_constraints (bool, optional): True の場合、制約条件もキーに含める。
            デフォルトは True。
        include_info (bool, optional): True の場合、`atoms.info` の単純な値も
            キーに含める。デフォルトは False。
        info_keys (list[str] | tuple[str, ...] | None, optional): info のうち
            含めるキー。None の場合は全キー候補を使う。

    Returns:
        dict[str, Any]: `stable_key()` に渡せる構造情報。

    Note:
        - calculator は fingerprint に含めない。計算条件は別途キーへ含めること。
        - info に JSON 化できない値が含まれる場合、そのキーはスキップされる。
    """
    # --- 基本構造情報 ---
    key_data: dict[str, Any] = {
        "symbols": atoms.get_chemical_symbols(),
        "positions": np.round(atoms.get_positions(), decimals).tolist(),
        "cell": np.round(atoms.get_cell().array, decimals).tolist(),
        "pbc": atoms.get_pbc().tolist(),
    }

    # --- 制約条件 ---
    if include_constraints:
        constraints = atoms.constraints
        if constraints:
            key_data["constraints"] = [
                _constraint_key_data(constraint) for constraint in constraints
            ]

    # --- atoms.info（任意） ---
    if include_info:
        selected_keys = info_keys if info_keys is not None else tuple(atoms.info.keys())
        info: dict[str, Any] = {}
        for key in selected_keys:
            if key not in atoms.info:
                continue
            try:
                info[str(key)] = _to_jsonable(atoms.info[key], decimals=decimals)
            except TypeError:
                continue
        if info:
            key_data["info"] = info

    return key_data


def atoms_fingerprint(
    atoms: Atoms,
    *,
    decimals: int = 8,
    include_constraints: bool = True,
    include_info: bool = False,
    info_keys: list[str] | tuple[str, ...] | None = None,
) -> str:
    """
    ase.Atoms から sha256 fingerprint を作成する。

    `atoms_key_data()` で構造情報を辞書化し、`stable_key()` でハッシュ化する。

    Args:
        atoms (ase.Atoms): キー化したい構造。
        decimals (int, optional): 座標・セルを丸める桁数。デフォルトは 8。
        include_constraints (bool, optional): 制約条件を含めるか。デフォルトは True。
        include_info (bool, optional): `atoms.info` を含めるか。デフォルトは False。
        info_keys (list[str] | tuple[str, ...] | None, optional): info のうち
            含めるキー。None の場合は全キー候補を使う。

    Returns:
        str: 64 文字の sha256 ハッシュ。

    Note:
        calculator は含めない。計算条件は別途 `stable_key()` の入力へ含めること。
    """
    return stable_key(
        atoms_key_data(
            atoms,
            decimals=decimals,
            include_constraints=include_constraints,
            include_info=include_info,
            info_keys=info_keys,
        )
    )


# ----------
# ---キー生成（ビルダー）
# ----------
@dataclass
class KeyBuilder:
    """
    条件を段階的に追加してキャッシュキーを作る補助クラス。

    複数の条件を名前付きで蓄積し、最後に `build()` で sha256 キーを生成する。

    Attributes:
        parts (dict[str, Any]): キー構成要素の辞書。
        decimals (int | None): float を丸める桁数。None の場合は丸めない。
    """

    parts: dict[str, Any] = field(default_factory=dict)
    decimals: int | None = None

    def add(self, name: str, value: Any) -> "KeyBuilder":
        """
        任意の値をキー構成要素として追加する。

        Args:
            name (str): 要素名（キー辞書のキーになる）。
            value: 追加する値。

        Returns:
            KeyBuilder: メソッドチェーン用に自身を返す。
        """
        self.parts[name] = value
        return self

    def add_atoms(
        self,
        name: str,
        atoms: Atoms,
        *,
        decimals: int = 8,
        include_constraints: bool = True,
        include_info: bool = False,
        info_keys: list[str] | tuple[str, ...] | None = None,
    ) -> "KeyBuilder":
        """
        ase.Atoms を fingerprint 化してキー構成要素として追加する。

        Args:
            name (str): 要素名。
            atoms (ase.Atoms): キー化したい構造。
            decimals (int, optional): 座標・セルを丸める桁数。デフォルトは 8。
            include_constraints (bool, optional): 制約条件を含めるか。デフォルトは True。
            include_info (bool, optional): `atoms.info` を含めるか。デフォルトは False。
            info_keys (list[str] | tuple[str, ...] | None, optional): info のうち
                含めるキー。

        Returns:
            KeyBuilder: メソッドチェーン用に自身を返す。
        """
        self.parts[name] = atoms_fingerprint(
            atoms,
            decimals=decimals,
            include_constraints=include_constraints,
            include_info=include_info,
            info_keys=info_keys,
        )
        return self

    def to_data(self) -> dict[str, Any]:
        """
        現状のKeyBuilderを辞書として返す。

        Returns:
            dict[str, Any]: `stable_key()` に渡せる辞書のコピー。
        """
        return dict(self.parts)

    def build(self) -> str:
        """
        蓄積した条件から sha256 キーを生成する。

        Returns:
            str: 64 文字の sha256 ハッシュ。
        """
        return stable_key(self.parts, decimals=self.decimals)


# ----------
# ---保存形式（StoreSpec）
# ----------
class StoreSpec(Protocol[T]):
    """
    CacheStore が利用する保存形式のインターフェース。

    値の型ごとに `save()` / `load()` の実装を切り替える。
    """

    name: str
    filename: str

    def save(self, directory: Path, value: T) -> None:
        """値をディレクトリへ保存する。"""
        ...

    def load(self, directory: Path) -> T:
        """ディレクトリから値を読み込む。"""
        ...


@dataclass
class JsonValueStore(Generic[T]):
    """
    JSON 互換値を保存する StoreSpec。

    float、int、str、dict、list などを `payload.json` として保存する。

    Attributes:
        filename (str): 保存ファイル名。デフォルトは ``payload.json``。
        name (str): 保存形式の識別名。デフォルトは ``json``。
        indent (int): JSON のインデント幅。デフォルトは 2。
    """

    filename: str = "payload.json"
    name: str = "json"
    indent: int = 2

    def save(self, directory: Path, value: T) -> None:
        """
        JSON 互換値をファイルへ保存する。

        Args:
            directory (Path): 保存先ディレクトリ。
            value: 保存する値。
        """
        path = directory / self.filename
        with path.open("w", encoding="utf-8") as handle:
            json.dump(
                _to_jsonable(value),
                handle,
                ensure_ascii=False,
                sort_keys=True,
                indent=self.indent,
            )
            handle.write("\n")

    def load(self, directory: Path) -> T:
        """
        JSON ファイルから値を読み込む。

        Args:
            directory (Path): 読み込み元ディレクトリ。

        Returns:
            T: 読み込んだ値。
        """
        path = directory / self.filename
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)


@dataclass
class AtomsStore:
    """
    ase.Atoms を ase.io で保存する StoreSpec。

    Attributes:
        filename (str): 保存ファイル名。デフォルトは ``structure.traj``。
        file_format (str | None): ASE のファイル形式。None の場合は拡張子から推定。
        name (str): 保存形式の識別名。デフォルトは ``atoms``。
        sanitize_xyz (bool): ``.xyz`` 保存時に `sanitize_atoms_for_xyz_write()` を
            使うか。デフォルトは False。

    Note:
        計算再利用を重視する場合は ``traj`` や ``extxyz`` を推奨する。
        ``xyz`` は情報落ちがあり得るため、デフォルトにはしない。
    """

    filename: str = "structure.traj"
    file_format: str | None = None
    name: str = "atoms"
    sanitize_xyz: bool = False

    def save(self, directory: Path, value: Atoms) -> None:
        """
        ase.Atoms をファイルへ保存する。

        Args:
            directory (Path): 保存先ディレクトリ。
            value (ase.Atoms): 保存する構造。

        Raises:
            RuntimeError: ``sanitize_xyz=True`` かつユーティリティを import できない場合。
        """
        atoms = value
        # --- xyz 書き出し時は安全な形式へ整形 ---
        if self.sanitize_xyz and self.filename.lower().endswith(".xyz"):
            if sanitize_atoms_for_xyz_write is None:
                raise RuntimeError("sanitize_atoms_for_xyz_write をimportできません。")
            atoms = sanitize_atoms_for_xyz_write(value)
        write(directory / self.filename, atoms, format=self.file_format)

    def load(self, directory: Path) -> Atoms:
        """
        ファイルから ase.Atoms を読み込む。

        Args:
            directory (Path): 読み込み元ディレクトリ。

        Returns:
            ase.Atoms: 読み込んだ構造。

        Raises:
            TypeError: 読み込み結果が ase.Atoms でない場合。
        """
        loaded = read(directory / self.filename, format=self.file_format)
        if not isinstance(loaded, Atoms):
            raise TypeError("保存された構造を ase.Atoms として読み込めませんでした。")
        return loaded


@dataclass
class PathStore:
    """
    ファイルやディレクトリのパスを JSON として保存する StoreSpec。

    計算済み構造ファイルなど、外部パスへの参照をキャッシュする用途向け。

    Attributes:
        filename (str): 保存ファイル名。デフォルトは ``path.json``。
        relative_to (str | Path | None): 相対パス化の基準ディレクトリ。
            None の場合は絶対/相対をそのまま記録する。
        name (str): 保存形式の識別名。デフォルトは ``path``。
    """

    filename: str = "path.json"
    relative_to: str | Path | None = None
    name: str = "path"

    def save(self, directory: Path, value: str | Path) -> None:
        """
        パスを JSON ファイルへ保存する。

        `relative_to` が指定されていれば、可能な限り相対パスで記録する。

        Args:
            directory (Path): 保存先ディレクトリ。
            value (str | Path): 保存するパス。
        """
        path_value = Path(value)
        # --- 相対パス化の試行 ---
        if self.relative_to is not None:
            try:
                path_text = str(path_value.relative_to(Path(self.relative_to)))
                is_relative = True
            except ValueError:
                path_text = str(path_value)
                is_relative = False
        else:
            path_text = str(path_value)
            is_relative = not path_value.is_absolute()

        data = {"path": path_text, "is_relative": is_relative}
        with (directory / self.filename).open("w", encoding="utf-8") as handle:
            json.dump(data, handle, ensure_ascii=False, sort_keys=True, indent=2)
            handle.write("\n")

    def load(self, directory: Path) -> Path:
        """
        JSON ファイルからパスを読み込む。

        Args:
            directory (Path): 読み込み元ディレクトリ。

        Returns:
            Path: 復元したパス。相対パス記録の場合は `relative_to` 基準で解決する。
        """
        with (directory / self.filename).open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        path = Path(data["path"])
        if data.get("is_relative") and self.relative_to is not None:
            return Path(self.relative_to) / path
        return path


# ----------
# ---キャッシュ本体
# ----------
@dataclass
class CacheStore(Generic[T]):
    """
    ファイルシステム上で計算結果を保存・取得するキャッシュ。

    1 キーにつき 1 ディレクトリを割り当て、`meta.json` とペイロードファイルを保存する。
    書き込みは一時ディレクトリ経由で行うことで、途中失敗で壊れたキャッシュが残りにくい。

    Attributes:
        base_dir (str | Path): キャッシュのルートディレクトリ。
        namespace (str): 用途ごとの名前空間（サブディレクトリ名）。
        store (StoreSpec[T]): 値の保存形式。
        schema_version (int): メタデータのスキーマ版。デフォルトは 1。

    Examples:
        >>> from tempfile import TemporaryDirectory
        >>> from ase_toolbox.Cache import CacheStore, JsonValueStore
        >>> td = TemporaryDirectory()
        >>> cache = CacheStore(td.name, "energy", JsonValueStore())
        >>> cache.set({"task": "demo"}, 1.23)
        1.23
        >>> cache.get({"task": "demo"})
        1.23
        >>> td.cleanup()
    """

    base_dir: str | Path
    namespace: str
    store: StoreSpec[T]
    schema_version: int = _CACHE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        """namespace の妥当性を検証し、base_dir を Path に正規化する。"""
        if not _SAFE_COMPONENT_RE.fullmatch(self.namespace):
            raise ValueError(
                "namespaceには英数字、アンダースコア、ハイフン、ドットのみ使用できます。"
            )
        self.base_dir = Path(self.base_dir)

    @property
    def namespace_dir(self) -> Path:
        """
        名前空間に対応するディレクトリパスを返す。

        Returns:
            Path: ``base_dir / namespace``。
        """
        return Path(self.base_dir) / self.namespace

    def _key_hash_and_data(
        self, key: str | Mapping[str, Any]
    ) -> tuple[str, Any | None]:
        """
        キー引数を、ハッシュ文字列・メタデータ用の元データへ変換する。

        既に 64 文字の sha256 文字列が渡された場合は、そのままハッシュとして扱う。

        Args:
            key (str | Mapping[str, Any]): キーまたはキー構成辞書。

        Returns:
            tuple[str, Any | None]: ``(key_hash, key_data)``。ハッシュ文字列入力時は
                ``key_data`` は None。
        """
        if isinstance(key, str) and _SHA256_RE.fullmatch(key):
            return key, None
        return stable_key(key), _to_jsonable(key)

    def path_for(self, key: str | Mapping[str, Any]) -> Path:
        """
        キーに対応するキャッシュエントリが保存されているディレクトリパスを返す。

        Args:
            key (str | Mapping[str, Any]): キーまたはキー構成辞書。

        Returns:
            Path: エントリディレクトリのパス。
        """
        key_hash, _ = self._key_hash_and_data(key)
        return self.namespace_dir / key_hash

    def exists(self, key: str | Mapping[str, Any]) -> bool:
        """
        キーに対応するキャッシュが存在するか判定する。

        ``meta.json`` と StoreSpec のペイロードファイルの両方が存在する場合に True。

        Args:
            key (str | Mapping[str, Any]): キーまたはキー構成辞書。

        Returns:
            bool: キャッシュが存在すれば True。
        """
        entry_dir = self.path_for(key)
        return (entry_dir / "meta.json").exists() and (
            entry_dir / self.store.filename
        ).exists()

    def metadata(self, key: str | Mapping[str, Any]) -> dict[str, Any]:
        """
        キーに対応する ``meta.json`` の内容を返す。

        Args:
            key (str | Mapping[str, Any]): キーまたはキー構成辞書。

        Returns:
            dict[str, Any]: メタデータ辞書。

        Raises:
            FileNotFoundError: キャッシュが存在しない場合。
        """
        with (self.path_for(key) / "meta.json").open("r", encoding="utf-8") as handle:
            return json.load(handle)

    def get(self, key: str | Mapping[str, Any]) -> T:
        """
        キャッシュから値を取得する。

        Args:
            key (str | Mapping[str, Any]): キーまたはキー構成辞書。

        Returns:
            T: 保存されていた値。

        Raises:
            KeyError: キャッシュが存在しない場合。
        """
        if not self.exists(key):
            raise KeyError(f"キャッシュが見つかりません: {self.path_for(key)}")
        # キーに対応するディレクトリを受け取り、そのディレクトリをStoreSpecのloadに渡して値を取得する。
        return self.store.load(self.path_for(key))

    def set(
        self,
        key: str | Mapping[str, Any],
        value: T,
        *,
        tags: Mapping[str, Any] | None = None,
    ) -> T:
        """
        値をキャッシュへ保存する。

        一時ディレクトリへ書き込んだ後、`os.replace()` で原子的に差し替える。

        Args:
            key (str | Mapping[str, Any]): キーまたはキー構成辞書。
            value (T): 保存する値。
            tags (Mapping[str, Any] | None, optional): メタデータに付与する任意タグ。

        Returns:
            T: 保存した値（引数 value と同じオブジェクト）。

        Raises:
            Exception: 保存処理中にエラーが発生した場合。一時ディレクトリは削除される。
        """
        # ---保存先名(=キー名)、保存するデータを用意する
        key_hash, key_data = self._key_hash_and_data(key)
        destination = self.namespace_dir / key_hash
        self.namespace_dir.mkdir(parents=True, exist_ok=True)

        # --- 一時ディレクトリへ書き込み ---
        tmp_dir = Path(
            tempfile.mkdtemp(prefix=f".tmp-{key_hash}-", dir=self.namespace_dir)
        )
        try:
            # ---StoreSpecを使って、値を保存する
            self.store.save(tmp_dir, value)

            # --- メタデータ作成 ---
            metadata = {
                "schema_version": self.schema_version,
                "created_at": datetime.now(timezone.utc).isoformat(),
                "namespace": self.namespace,
                "key_hash": key_hash,
                "key_data": key_data,
                "store": {
                    "name": self.store.name,
                    "filename": self.store.filename,
                },
                "tags": _to_jsonable(tags or {}),
            }
            with (tmp_dir / "meta.json").open("w", encoding="utf-8") as handle:
                json.dump(metadata, handle, ensure_ascii=False, sort_keys=True, indent=2)
                handle.write("\n")

            # --- 原子的に差し替え ---
            backup = destination.with_name(f".old-{destination.name}")
            if backup.exists():
                shutil.rmtree(backup)
            if destination.exists():
                os.replace(destination, backup)
            os.replace(tmp_dir, destination) # 一時ディレクトリを、保存先に差し替える(移動して置換みたいな)
            if backup.exists():
                shutil.rmtree(backup)
        except Exception:
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir)
            raise

        return value

    def get_or_compute(
        self,
        key: str | Mapping[str, Any],
        compute: Callable[[], T],
        *,
        tags: Mapping[str, Any] | None = None,
        overwrite: bool = False,
    ) -> T:
        """
        キャッシュがあれば取得し、なければ計算して保存する。

        Args:
            key (str | Mapping[str, Any]): キーまたはキー構成辞書。
            compute (Callable[[], T]): キャッシュが存在しないときに呼ぶ計算関数。
            tags (Mapping[str, Any] | None, optional): 保存時に付与する任意タグ。
            overwrite (bool, optional): True の場合、既存キャッシュを無視して再計算する。
                デフォルトは False。

        Returns:
            T: キャッシュまたは新規計算の結果。
        """
        if not overwrite and self.exists(key):
            return self.get(key)

        value = compute()
        return self.set(key, value, tags=tags)

    def invalidate(self, key: str | Mapping[str, Any]) -> bool:
        """
        指定キーのキャッシュエントリを削除する。

        Args:
            key (str | Mapping[str, Any]): キーまたはキー構成辞書。

        Returns:
            bool: 削除した場合は True。もともと存在しなければ False。
        """
        entry_dir = self.path_for(key)
        if not entry_dir.exists():
            return False
        shutil.rmtree(entry_dir)
        return True

    def clear(self) -> None:
        """
        名前空間配下のキャッシュをすべて削除する。

        ``namespace_dir`` 自体を削除する。存在しない場合は何もしない。
        """
        if self.namespace_dir.exists():
            shutil.rmtree(self.namespace_dir)
