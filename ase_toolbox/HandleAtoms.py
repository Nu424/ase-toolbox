"""
# HandleAtoms.py
ASEのAtomsオブジェクトを操作・処理するための関数をまとめたファイル

## 関数一覧
- smiles_to_atoms(): SMILES文字列からASE Atomsオブジェクトを生成する
- set_substrate_mask_all(): 全原子に基板マスクを設定する
- move_atoms(): 原子を指定方向に指定距離だけ移動する
- fix_layers(): (平面用)層を固定する
- substitute_elements(): 原子を置き換える
- compute_surface_normal(): (クラスター用)法線ベクトルを計算する
- place_adsorbate_along_normal(): (クラスター用)クラスターにくっつける
- place_adsorbate_on_surface(): (表面用)表面にくっつける
"""

from typing import Literal, Sequence, Mapping, Optional
from ase import Atoms, Atom
from ase.constraints import FixAtoms
import numpy as np
from numpy.typing import NDArray
import re
import warnings
from .FindAtoms import (
    find_central_atom,
    separate_layers,
    get_neighbors_with_coordination_condition,
    get_neighbors,
)
from .util import ConditionalLogger, ensure_logger, resolve_target_indices


# 例外クラス
class InvalidElementSymbolError(ValueError):
    """元素記号が不正な場合に投げる例外。"""


class LatticeConstantNotFoundError(KeyError):
    """格子定数が見つからない場合に投げる例外。"""


class CompositionSumError(ValueError):
    """組成の合計が1にならない等の不整合に対して投げる例外。"""


# 内部デフォルトの格子定数マップ（代表値, 単位: Å）
# 注: 室温代表値/相依存。実使用ではユーザが上書きしてください。
INTERNAL_LATTICE_MAP: dict[str, float] = {
    # fcc
    "Al": 4.049,
    "Cu": 3.615,
    "Au": 4.078,
    "Ag": 4.086,
    "Ni": 3.524,
    "Pd": 3.890,
    "Pt": 3.923,
    # bcc（立方）
    "Fe": 2.866,
    "V": 3.027,
    "Nb": 3.300,
    "Mo": 3.147,
    "W": 3.165,
    # hcp（a格）
    "Mg": 3.209,
    "Ti": 2.951,
}

# 体積混合法の注意喚起用：立方系ではない代表元素（hcp）と bcc/hcp の集合
_NON_CUBIC_ELEMENTS: set[str] = {"Mg", "Ti", "Zn", "Zr", "Co", "Cd", "Be", "Ru", "Os", "Re"}
_BCC_ELEMENTS: set[str] = {"Fe", "Cr", "W", "Mo", "V", "Nb", "Ta", "Ba"}


# SMILES文字列からASE Atomsオブジェクトを生成する
def smiles_to_atoms(
    smiles: str,
    *,
    optimize: Literal["UFF", "MMFF"] | None = "UFF",
    random_seed: int | None = None,
) -> Atoms:
    """
    SMILES文字列からASE Atomsオブジェクトを生成する。

    RDKitを使用してSMILES文字列を3D構造に変換します。
    ETKDG法で3D座標を埋め込み、オプションで力場最適化を実行します。

    Args:
        smiles (str): SMILES文字列（例: "CCO", "c1ccccc1"）。
        optimize (Literal["UFF", "MMFF"] | None, optional): 力場最適化の種類。
            "UFF": Universal Force Field で最適化。
            "MMFF": Merck Molecular Force Field で最適化。
            None: 最適化なし（ETKDG埋め込みのみ）。
            デフォルトは "UFF"。
        random_seed (int | None, optional): 3D座標埋め込み時の乱数シード。
            Noneの場合はランダム。デフォルトはNone。

    Returns:
        ase.Atoms: 生成された分子構造。単位はÅ。

    Raises:
        ImportError: RDKitがインストールされていない場合。
        ValueError: SMILES文字列が不正、3D座標埋め込みに失敗、または
                   最適化に失敗した場合。

    Examples:
        >>> # エタノールの生成
        >>> ethanol = smiles_to_atoms("CCO")
        >>> print(f"原子数: {len(ethanol)}")

        >>> # ベンゼンの生成（最適化なし）
        >>> benzene = smiles_to_atoms("c1ccccc1", optimize=None)

        >>> # 再現性のある生成
        >>> mol = smiles_to_atoms("CC(C)O", random_seed=42)

    Note:
        - この関数はRDKitに依存します。未導入の場合は明確なエラーメッセージを表示します。
        - 3D座標はETKDGv3アルゴリズムで生成されます。
        - 水素原子は自動的に付加されます。
    """
    # --- RDKitの遅延import ---
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise ImportError(
            "RDKitがインストールされていません。\n"
            "以下のコマンドでインストールしてください:\n"
            "  pip install rdkit\n"
            "または\n"
            "  conda install -c conda-forge rdkit"
        ) from e

    # --- SMILES文字列から分子オブジェクトを作成 ---
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(
            f"SMILES文字列の解析に失敗しました: '{smiles}'\n"
            "SMILES文字列が正しいか確認してください。"
        )

    # --- 水素原子を付加 ---
    mol = Chem.AddHs(mol)

    # --- 3D座標を埋め込む（ETKDG法） ---
    params = AllChem.ETKDGv3()
    if random_seed is not None:
        params.randomSeed = random_seed

    embed_result = AllChem.EmbedMolecule(mol, params)
    if embed_result == -1:
        raise ValueError(
            f"3D座標の埋め込みに失敗しました: '{smiles}'\n"
            "分子が大きすぎる、または構造的に不安定な可能性があります。"
        )

    # --- 力場最適化（オプション） ---
    if optimize == "UFF":
        converged = AllChem.UFFOptimizeMolecule(mol)
        if converged != 0:
            warnings.warn(
                f"UFF最適化が完全に収束しませんでした（収束コード: {converged}）。\n"
                "結果の精度が低い可能性があります。",
                RuntimeWarning,
            )
    elif optimize == "MMFF":
        # MMFFプロパティの取得
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is None:
            raise ValueError(
                f"MMFF力場のパラメータ取得に失敗しました: '{smiles}'\n"
                "この分子にはMMFFを適用できません。optimize='UFF'を試してください。"
            )
        converged = AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94")
        if converged != 0:
            warnings.warn(
                f"MMFF最適化が完全に収束しませんでした（収束コード: {converged}）。\n"
                "結果の精度が低い可能性があります。",
                RuntimeWarning,
            )
    elif optimize is not None:
        raise ValueError(
            f"optimizeは 'UFF'、'MMFF'、または None を指定してください。指定値: {optimize}"
        )

    # --- RDKit分子からASE Atomsへの変換 ---
    symbols = []
    positions = []
    conf = mol.GetConformer()

    for atom in mol.GetAtoms():
        # 元素記号を取得
        symbols.append(atom.GetSymbol())
        # 3D座標を取得（単位: Å）
        pos = conf.GetAtomPosition(atom.GetIdx())
        positions.append([pos.x, pos.y, pos.z])

    # --- ASE Atomsオブジェクトを作成 ---
    ase_atoms = Atoms(symbols=symbols, positions=positions)

    return ase_atoms


# 元素の格子定数を取得
def _get_element_lattice_constant(symbol: str, user_map: dict[str, float] | None = None) -> float:
    """
    元素の格子定数[a]（Å）を取得する。優先順は user_map > INTERNAL_LATTICE_MAP。

    Args:
        symbol: 元素記号（例: 'Cu'）。大文字小文字は自動正規化。
        user_map: ユーザ提供の格子定数マップ。キーは元素記号、値は格子定数[Å]。

    Returns:
        float: 格子定数[Å]。

    Raises:
        InvalidElementSymbolError: 元素記号が不正な形式の場合。
        LatticeConstantNotFoundError: 対象元素の格子定数が見つからない場合。
    """
    if not isinstance(symbol, str):
        raise InvalidElementSymbolError("元素記号 symbol は str を指定してください。")

    sym = symbol.capitalize()
    if re.fullmatch(r"^[A-Z][a-z]?$", sym) is None:
        raise InvalidElementSymbolError(f"元素記号の形式が不正です: {symbol}")

    # user_map を優先
    if user_map is not None:
        # キーの大小統一のために正規化
        normalized_user_map = {k.capitalize(): float(v) for k, v in user_map.items()}
        if sym in normalized_user_map:
            return float(normalized_user_map[sym])

    if sym in INTERNAL_LATTICE_MAP:
        return float(INTERNAL_LATTICE_MAP[sym])

    raise LatticeConstantNotFoundError(
        f"元素 '{sym}' の格子定数が見つかりません。user_map に追加してください。"
    )


# 組成から混合格子定数を計算
def mix_lattice_constant(
    composition: dict[str, float],
    lattice_map: dict[str, float] | None = None,
    method: Literal["vegard", "volume"] = "vegard",
    *,
    tol: float = 1e-6,
    return_detail: bool = False,
) -> float | tuple[float, dict]:
    """
    組成に基づいて混合格子定数 a を計算する。

    - vegard: a = Σ x_i a_i
    - volume: a = (Σ x_i a_i^3)^(1/3)  （立方晶を想定）

    Args:
        composition: {元素記号: 比率} の辞書。比率の合計はおよそ1。
        lattice_map: ユーザ提供の格子定数マップ。
        method: 'vegard'（線形）または 'volume'（体積基準）。
        tol: 合計1.0 の許容誤差。
        return_detail: True の場合、詳細情報の辞書も返す。

    Returns:
        float | tuple[float, dict]: 混合格子定数（Å）。return_detail=True の場合は (a, detail)。

    Raises:
        CompositionSumError: 比率合計が1から外れる、負値/ゼロが含まれる等。
        InvalidElementSymbolError: 不正な元素記号。
        LatticeConstantNotFoundError: 格子定数未定義の元素を含む場合。
        NotImplementedError: 未対応の method を指定した場合。
    """
    if not isinstance(composition, dict) or len(composition) == 0:
        raise CompositionSumError("composition は非空の dict[str, float] を指定してください。")

    # 検証と正規化
    symbols: list[str] = []
    fractions: list[float] = []
    for k, v in composition.items():
        sym = k.capitalize()
        if re.fullmatch(r"^[A-Z][a-z]?$", sym) is None:
            raise InvalidElementSymbolError(f"元素記号の形式が不正です: {k}")
        f = float(v)
        if f <= 0:
            raise CompositionSumError(f"組成比は正の値である必要があります: {k}={v}")
        symbols.append(sym)
        fractions.append(f)

    total = float(np.sum(fractions))
    if not np.isclose(total, 1.0, atol=tol):
        raise CompositionSumError(f"組成の合計が1ではありません（現在: {total:.6f}）。")

    # a_i の収集
    a_map: dict[str, float] = {}
    for sym in symbols:
        a_map[sym] = _get_element_lattice_constant(sym, user_map=lattice_map)

    # 計算
    x = np.array(fractions, dtype=float)
    a_list = np.array([a_map[sym] for sym in symbols], dtype=float)

    if method == "vegard":
        a_mixed = float(np.dot(x, a_list))
    elif method == "volume":
        # 注意喚起（立方晶想定）
        if any((s in _NON_CUBIC_ELEMENTS) for s in symbols) or any((s in _BCC_ELEMENTS) for s in symbols):
            warnings.warn(
                "method='volume' は立方晶を仮定しています。bcc/hcp を含む系では注意してください。",
                RuntimeWarning,
            )
        a_mixed = float(np.power(np.dot(x, np.power(a_list, 3.0)), 1.0 / 3.0))
    else:
        raise NotImplementedError(f"未対応の method です: {method}")

    if not return_detail:
        return a_mixed

    detail = {
        "method": method,
        "composition": {s: float(f) for s, f in zip(symbols, fractions)},
        "constants": {s: float(a_map[s]) for s in symbols},
        "a_mixed": float(a_mixed),
    }
    return a_mixed, detail


# 全原子に基板マスクを設定する
def set_substrate_mask_all(
    atoms: Atoms,
    *,
    is_substrate: bool = True,
    inplace: bool = True,
) -> Atoms:
    """
    全原子に is_substrate マスクを設定する。

    吸着分子を配置する前の基板構造に対して、全原子を「基板」としてマークするために使用します。
    これにより、後続の吸着分子配置時に層検出や高さ基準が基板のみに基づいて行われます。

    Args:
        atoms (ase.Atoms): マスクを設定する原子構造。
        is_substrate (bool, optional): 設定する値。True で基板、False で非基板。デフォルトは True。
        inplace (bool, optional): True の場合は atoms を直接変更。
            False の場合はコピーを作成して返します。デフォルトは True。

    Returns:
        ase.Atoms: マスクが設定された原子構造。inplace=True の場合は引数 atoms 自身が返されます。

    Examples:
        >>> from ase.build import bulk, surface
        >>> slab = surface(bulk('Cu'), (1,1,1), layers=3)
        >>> # 基板として初期化
        >>> slab = set_substrate_mask_all(slab, is_substrate=True)
        >>> # この後、place_adsorbate_on_surface を複数回呼んでも、層検出は常に基板のみで行われる
    """
    target = atoms if inplace else atoms.copy()
    mask = np.full(len(target), is_substrate, dtype=bool)
    target.set_array("is_substrate", mask)
    return target


# 手動で微妙に動かす
def move_atoms(
    base_structure: Atoms,
    target: int | Atom | Atoms | list[int] | list[Atom],
    direction: tuple[float, float, float] | Sequence[float] | NDArray[np.floating],
    distance: float,
    *,
    inplace: bool = False,
) -> Atoms:
    """
    指定した原子を指定方向に指定距離だけ移動する。

    複数の入力形式に対応し、原子を任意の方向に微小変位させることができます。
    方向ベクトルは自動的に正規化され、指定距離だけ移動されます。

    Args:
        base_structure (ase.Atoms): 操作対象の原子構造。
        target (int | Atom | Atoms | list[int] | list[Atom]): 移動させる原子の指定。以下の形式に対応：
            - int: 単一原子のインデックス
            - ase.Atom: 単一原子オブジェクト
            - ase.Atoms: 構造全体（全原子を移動）
            - list[int]: 複数原子のインデックスリスト
            - list[ase.Atom]: 複数原子オブジェクトのリスト
        direction (tuple[float, float, float] | Sequence[float] | NDArray[np.floating]): 移動方向を示すベクトル。
            tuple[float, float, float]、Sequence[float]、またはnumpy配列を指定可能。
            自動的に正規化されます。
        distance (float): 移動させる距離 [Å]。正の値で指定方向、負の値で逆方向。
        inplace (bool, optional): Trueの場合、base_structureを直接変更。
            Falseの場合、コピーを作成して返します。デフォルトはFalse。

    Returns:
        ase.Atoms: 原子が移動された構造。inplace=Trueの場合は引数と同じオブジェクト。

    Raises:
        ValueError: 方向ベクトルがゼロベクトルの場合、または指定されたAtomが
                   base_structure内に存在しない場合。
        TypeError: target、directionの型が不正な場合。
        IndexError: 指定されたインデックスが範囲外の場合。

    """
    # --- 入力検証：direction ---
    try:
        direction_vec = np.array(direction, dtype=float)
    except (ValueError, TypeError) as e:
        raise TypeError(
            f"directionは数値のsequenceまたはnumpy配列を指定してください: {e}"
        )

    if direction_vec.shape != (3,):
        raise ValueError(
            f"directionは3次元ベクトルを指定してください。現在の形状: {direction_vec.shape}"
        )

    # --- ゼロベクトルチェック ---
    norm = np.linalg.norm(direction_vec)
    if norm == 0:
        raise ValueError("移動方向ベクトルがゼロです。正しい方向を指定してください。")

    # --- 単位ベクトルと変位ベクトルを計算 ---
    unit_vector = direction_vec / norm
    displacement = unit_vector * distance

    # --- 作業用構造の準備 ---
    if inplace:
        structure = base_structure
    else:
        structure = base_structure.copy()

    # --- target から移動対象原子のインデックスリストを取得 ---
    target_indices = resolve_target_indices(base_structure, target)

    # --- 原子の移動を実行 ---
    for idx in target_indices:
        structure[idx].position += displacement

    return structure


# (平面用)層を固定する
def fix_layers(
    atoms: Atoms,
    fixed_layers: int,
    *,
    inplace: bool = False,
    decimals: int = 4,
    logger: Optional[ConditionalLogger] = None,
    enable_logging: bool = False,
    use_substrate_mask: Literal["auto", True, False] = "auto",
) -> Atoms:
    """
    指定された数の下層を固定する制約を追加する。

    z座標に基づいて層を特定し、下から指定された数の層に含まれる
    全ての原子に FixAtoms 制約を適用します。

    Args:
        atoms (ase.Atoms): 制約を適用する原子構造。
        fixed_layers (int): 固定する層数（下から数えて）。
            0 の場合は何も固定しません。
        inplace (bool, optional): True の場合は atoms を破壊的に更新。
            False の場合はコピーを作成して返します。デフォルトはFalse。
        decimals (int, optional): z座標の丸め精度（小数点以下の桁数）。
            デフォルトは4。層の判定精度に影響します。
        logger (Optional[ConditionalLogger], optional): ログ出力制御。
            Noneの場合は新規作成。
        enable_logging (bool, optional): ログ出力の有効/無効。デフォルトは True。
        use_substrate_mask (Literal["auto", True, False], optional): 基板マスクの使用設定。
            "auto": atoms.arrays に "is_substrate" が存在する場合、基板原子のみで層を検出。
            True: 基板マスクを使用（存在しない場合は全原子）。
            False: マスクを無視して全原子を対象。
            デフォルトは "auto"。

    Returns:
        ase.Atoms: 制約が適用された原子構造。
            inplace=True の場合は引数 atoms 自身が返されます。

    Raises:
        ValueError: fixed_layers が負の値の場合。
    """
    # --- ログ設定 ---
    logger = ensure_logger("fix_layers", enable_logging, logger)

    # --- 引数の検証 ---
    if fixed_layers < 0:
        raise ValueError("fixed_layers は0以上の値を指定してください。")

    # --- 作業用原子構造の準備 ---
    if inplace:
        result_atoms = atoms
    else:
        result_atoms = atoms.copy()

    # --- 固定層数が0の場合は何もしない ---
    if fixed_layers == 0:
        return result_atoms

    # --- 層別に分離して固定対象を特定 ---
    layers_indices = separate_layers(
        atoms,
        return_type="indices",
        decimals=decimals,
        sort_by_z=True,  # 下層から上層の順序
        use_substrate_mask=use_substrate_mask,
    )

    total_layers = len(layers_indices)

    # --- 固定層数が総層数以上の場合の警告 ---
    if fixed_layers >= total_layers:
        logger.warning(
            f"固定層数 ({fixed_layers}) が総層数 ({total_layers}) 以上です。"
        )
        logger.warning("全ての原子が固定されます。")
        fixed_layers = total_layers

    # --- 固定対象原子のインデックスを収集 ---
    fixed_indices: list[int] = []
    for layer_idx in range(fixed_layers):
        fixed_indices.extend(layers_indices[layer_idx])

    # --- FixAtoms制約を作成・適用 ---
    if fixed_indices:
        # マスク配列を作成（True = 固定する原子）
        mask = np.zeros(len(result_atoms), dtype=bool)
        mask[fixed_indices] = True

        constraint = FixAtoms(mask=mask)
        result_atoms.set_constraint(constraint)

        # 最上固定層のz座標を取得（報告用）
        fixed_z_coords = result_atoms.positions[fixed_indices, 2]
        max_fixed_z = np.max(np.round(fixed_z_coords, decimals=decimals))

        logger.info(
            f"Z≤{max_fixed_z:.{decimals}f} Å の原子 {len(fixed_indices)} 個を固定しました。"
        )
    else:
        logger.info("固定対象の原子が見つかりませんでした。")

    return result_atoms


# 置き換える
def substitute_elements(
    atoms: Atoms,
    target: int | Atom | Atoms | list[int] | list[Atom],
    new: str | Mapping[str, float],
    *,
    inplace: bool = False,
    seed: Optional[int] = None,
) -> Atoms:
    """
    指定した原子（または全原子）を指定した元素に置換します。

    置換後の元素は、単一の元素記号（例: "Pd"）または組成辞書
    （例: {"Cu": 0.8, "Pd": 0.2}）で指定できます。

    Args:
        atoms (ase.Atoms): 操作対象の構造。
        target (Union[int, ase.Atom, ase.Atoms, list[int], list[ase.Atom]]): 置換対象の指定。
            - int: 単一原子のインデックス。
            - list[int]: 複数原子のインデックス。
            - ase.Atom: 単一原子（atoms内に存在する必要あり）。
            - list[ase.Atom]: 複数原子（atoms内に存在する必要あり）。
            - ase.Atoms: 構造全体を対象（全原子を置換）。
        new (Union[str, Mapping[str, float]]): 置換後の指定。
            - str: 単一元素記号 (例: "Cu")。
            - Mapping[str, float]: 組成辞書 (例: {"Cu": 0.9, "Pd": 0.1})。合計は1。
        inplace (bool, optional): True の場合は atoms を破壊的に更新。
            False の場合はコピーを作成して返す。デフォルトは False。
        seed (Optional[int], optional): シャッフルの再現性確保のための乱数シード。

    Returns:
        ase.Atoms: 要求どおりに置換された Atoms オブジェクト。
                   inplace=True の場合は引数 atoms 自身が返ります。

    """
    # ----------
    # --- 置換対象インデックスの決定
    # ----------
    atom_indices = resolve_target_indices(atoms, target)

    n_targets = len(atom_indices)
    if n_targets == 0:
        # 置換対象なしの場合はそのまま返す
        return atoms if inplace else atoms.copy()

    # ----------
    # --- 置換後のシンボル列を生成
    # ----------
    replacement_symbols: list[str] = []
    # --- 単一シンボル指定の場合 ---
    if isinstance(new, str):
        replacement_symbols = [new] * n_targets

    # --- 組成辞書指定の場合 ---
    elif isinstance(new, Mapping):
        if n_targets < 0:
            raise ValueError("count は0以上である必要があります。")
        if n_targets == 0:
            replacement_symbols = []
        else:
            # 合計の検証（許容誤差内で1）
            total = float(sum(new.values()))
            if not np.isclose(total, 1.0):
                raise ValueError(f"組成の合計が1ではありません（現在: {total}）。")

            # 元素記号の検証
            symbols = list(new.keys())

            # --- 丸めベースで個数を算出し、最後の元素に差分を集約 ---
            counts = {s: int(round(n_targets * frac))
                      for s, frac in new.items()}

            # 丸め誤差の補正
            diff = n_targets - sum(counts.values())
            counts[symbols[-1]] += diff

            # --- シンボル列を生成し、シャッフル ---
            replacement_symbols = []
            for s, c in counts.items():
                replacement_symbols.extend([s] * c)

            # 万一長さが合わない場合は例外
            if len(replacement_symbols) != n_targets:
                raise RuntimeError(
                    f"置換シンボル列の長さが一致しません: {len(replacement_symbols)} != {n_targets}"
                )

            rng = np.random.default_rng(seed)
            rng.shuffle(replacement_symbols)
    else:
        raise TypeError("new は str または Mapping[str, float] を指定してください。")

    # ----------
    # --- 書き換え先の Atoms を準備
    # ----------
    if inplace:
        out = atoms
    else:
        out = atoms.copy()

    # --- インデックス順にシンボルを割り当て ---
    # 対象インデックスの順序に対応するようにする
    for idx, sym in zip(atom_indices, replacement_symbols):
        out[idx].symbol = sym

    return out


# 法線ベクトルを計算する
def compute_surface_normal(
    atoms: Atoms,
    target_atom: int | Atom,
    *,
    include_target: bool = True,
    reference_vector: NDArray[np.float64] | None = None,
    normalize: bool = True,
    return_plane: bool = False,
) -> NDArray[np.float64] | tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """
    指定した原子周辺の局所平面を主成分分析（PCA）で近似し、その法線ベクトルを計算する。

    指定した原子とその隣接原子の座標を用いてPCAを実行し、最小固有値に対応する
    固有ベクトルを表面の法線ベクトルとして算出する。

    Args:
        atoms (ase.Atoms): 原子構造を保持するASEのAtomsオブジェクト。
        target_atom (int | ase.Atom): 対象原子。インデックスまたはAtomオブジェクトを指定可能。
        include_target (bool, optional): PCA計算時に対象原子自身を含めるかどうか。
            デフォルトはTrue。
        reference_vector (NDArray[np.float64] | None, optional): 法線ベクトルの符号を
            決定するための参照ベクトル。Noneの場合は符号調整を行わない。
            デフォルトはNone。
        normalize (bool, optional): 法線ベクトルをユニットベクトルに正規化するか。
            デフォルトはTrue。
        return_plane (bool, optional): 平面情報（法線、重心、d係数）も返すか。
            Falseの場合は法線ベクトルのみ、Trueの場合は(normal, centroid, d)の
            タプルを返す。デフォルトはFalse。

    Returns:
        NDArray[np.float64] | tuple[NDArray[np.float64], NDArray[np.float64], float]:
            return_plane=Falseの場合: 法線ベクトル（形状: (3,)）
            return_plane=Trueの場合: (法線ベクトル, 重心, d係数)のタプル

            平面の方程式: normal · (r - centroid) = 0 または normal · r + d = 0

    Raises:
        TypeError: target_atomの型がintでもase.Atomでもない場合。
        ValueError: 指定されたAtomがatoms内に存在しない場合、または
                   平面を定義するのに十分な点数がない場合（3点未満）、または
                   点が共線で平面が定義できない場合。
        IndexError: 指定されたインデックスが範囲外の場合。

    Note:
        この関数は事前に定義されたget_neighbors関数を使用して隣接原子を取得する。
        PCAにより平面に最もよくフィットする法線ベクトルを求めるため、
        ノイズがある場合でも安定した結果が得られる。
    """
    # --- target_atom の型に応じてインデックスを取得 ---
    if isinstance(target_atom, int):
        index = target_atom
    elif isinstance(target_atom, Atom):
        try:
            # Atomオブジェクトがatoms内に存在する場合、そのインデックスを取得
            index = target_atom.index
        except ValueError:
            raise ValueError("指定されたAtomはatoms内に存在しません。")
    else:
        raise TypeError("target_atom は int または ase.Atom を指定してください。")

    # --- インデックスの範囲チェック ---
    if not (0 <= index < len(atoms)):
        raise IndexError(
            f"指定されたインデックス {index} は範囲外です（0-{len(atoms)-1}）。"
        )

    # --- 隣接原子のインデックスを取得 ---
    neighbor_indices = get_neighbors(atoms, index, return_type="indices")

    # --- PCA用の点群を構築 ---
    point_indices = list(neighbor_indices)
    if include_target:
        point_indices.append(index)

    # --- 点数の妥当性チェック ---
    if len(point_indices) < 3:
        raise ValueError(
            f"平面を一意に定義するための点数が不足しています。"
            f"必要: 3点以上、取得: {len(point_indices)}点"
        )

    # --- 座標データを取得 ---
    points = np.array([atoms[i].position for i in point_indices])

    # --- 1. 重心を計算 ---
    centroid = np.mean(points, axis=0)

    # --- 2. 中心化 ---
    centered_points = points - centroid

    # --- 3. 共分散行列を計算 ---
    n_points = len(points)
    H = np.dot(centered_points.T, centered_points) / n_points

    # --- 4. 固有値分解 ---
    eigvals, eigvecs = np.linalg.eigh(H)

    # --- 5. 最小固有値の固有ベクトルを法線ベクトルとする ---
    min_eigval_idx = np.argmin(eigvals)
    normal = eigvecs[:, min_eigval_idx].copy()

    # --- 共線性（退化）のチェック ---
    # 固有値を昇順にソート
    sorted_eigvals = np.sort(eigvals)
    min_eigval = sorted_eigvals[0]
    second_min_eigval = sorted_eigvals[1]

    # 最小固有値が非常に小さく、かつ2番目との比が小さい場合は共線と判定
    min_tolerance = 1e-12
    ratio_tolerance = 1e-6

    if min_eigval < min_tolerance:
        raise ValueError(
            "指定された点群が共線状態で、平面を一意に定義できません。"
            f"最小固有値: {min_eigval:.2e} < 許容値: {min_tolerance:.2e}"
        )

    # 固有値の比が小さすぎる場合も共線と判定
    if second_min_eigval > 0 and min_eigval / second_min_eigval < ratio_tolerance:
        raise ValueError(
            "指定された点群が近似的に共線状態で、平面を安定に定義できません。"
            f"固有値比: {min_eigval / second_min_eigval:.2e} < 許容値: {ratio_tolerance:.2e}"
        )

    # --- 6. 参照ベクトルによる符号調整 ---
    if reference_vector is not None:
        ref_vec = np.asarray(reference_vector)
        if np.allclose(ref_vec, 0):
            raise ValueError("reference_vectorはゼロベクトルにできません。")

        # 内積が負の場合は符号を反転
        if np.dot(normal, ref_vec) < 0:
            normal = -normal

    # --- 7. 正規化 ---
    if normalize:
        normal = normal / np.linalg.norm(normal)

    # --- 8. 返却形式に応じて出力 ---
    if return_plane:
        # 平面の方程式: normal · r + d = 0
        d = -np.dot(normal, centroid)
        return normal, centroid, d
    else:
        return normal


# (クラスター用)クラスターにくっつける
def place_adsorbate_along_normal(
    substrate: Atoms,
    adsorbate: Atoms,
    target_atom: int | Atom,
    distance: float,
    *,
    upper_tolerance: int = 1,
    lower_tolerance: int = 1,
) -> Atoms:
    """
    指定した基板の表面法線方向に、吸着分子を配置する。

    指定した基板の原子を基準として、その周辺の局所平面を主成分分析（PCA）で求め、
    法線方向に指定距離だけ離れた位置に吸着分子を配置する。吸着分子の+z軸が
    法線ベクトルに一致するように回転調整も行う。

    Args:
        substrate (ase.Atoms): 基板となる原子構造。
        adsorbate (ase.Atoms): 配置する吸着分子。
        target_atom (int | ase.Atom): 基準となる基板上の原子。
            インデックスまたはAtomオブジェクトを指定可能。
        distance (float): 法線方向に離す距離 [Å]。正の値を指定。
        upper_tolerance (int): 最小配位数から上方向への許容範囲。デフォルトは1。
        lower_tolerance (int): 最小配位数から下方向への許容範囲。デフォルトは1。

    Returns:
        ase.Atoms: 基板と配置済み吸着分子を結合した構造。

    Raises:
        ValueError: distance が負の値の場合、または法線計算に必要な点数が不足している場合。
        TypeError: target_atom の型が int でも ase.Atom でもない場合。
        IndexError: 指定されたインデックスが範囲外の場合。

    Note:
        - 法線ベクトルの向きは「重心→target_atom」方向を参照として決定されるため、
          常に基板表面から外側を向く。
        - 吸着分子の内部+z軸が法線に整列するように回転される。
        - この関数は compute_surface_normal と get_neighbors を利用する。
    """
    # --- 入力検証 ---
    if distance < 0:
        raise ValueError(f"distance は正の値を指定してください。指定値: {distance}")

    # --- target_atom の型に応じてインデックスを取得 ---
    if isinstance(target_atom, int):
        target_index = target_atom
    elif isinstance(target_atom, Atom):
        try:
            target_index = target_atom.index
        except ValueError:
            raise ValueError("指定されたAtomはsubstrate内に存在しません。")
    else:
        raise TypeError("target_atom は int または ase.Atom を指定してください。")

    # --- インデックスの範囲チェック ---
    if not (0 <= target_index < len(substrate)):
        raise IndexError(
            f"指定されたインデックス {target_index} は範囲外です（0-{len(substrate)-1}）。"
        )

    # --- 局所点群の重心を計算 ---
    neighbor_indices = get_neighbors_with_coordination_condition(
        substrate,
        target_index,
        return_type="indices",
        upper_tolerance=upper_tolerance,
        lower_tolerance=lower_tolerance,
    )

    point_indices = list(neighbor_indices) + \
        [target_index]  # include_target=True相当

    if len(point_indices) < 3:
        raise ValueError(
            f"法線計算に必要な点数が不足しています。"
            f"必要: 3点以上、取得: {len(point_indices)}点"
        )

    # points = np.array([substrate[i].position for i in point_indices])
    # centroid = np.mean(points, axis=0)  # 現在未使用だが将来の拡張用

    # --- 参照ベクトルの計算（重心→target_atom方向） ---
    target_pos = substrate[target_index].position

    substrate_centroid = substrate.get_center_of_mass()
    reference_vector = target_pos - substrate_centroid

    # --- 法線ベクトルを計算 ---
    normal = compute_surface_normal(
        substrate,
        target_index,
        include_target=True,
        reference_vector=reference_vector,
        normalize=True,
        return_plane=False,
    )

    # --- 吸着分子のコピーを作成 ---
    adsorbate_copy = adsorbate.copy()

    # --- 吸着分子の重心を計算 ---
    com = adsorbate_copy.get_center_of_mass()

    # --- z軸を法線に整列する回転行列を計算 ---
    z_axis = np.array([0.0, 0.0, 1.0])
    rotation_matrix = _compute_rotation_matrix(z_axis, normal)

    # --- 重心を中心とした回転を適用 ---
    for atom in adsorbate_copy:
        # 重心からの相対位置を計算
        rel_pos = atom.position - com
        # 回転を適用
        rotated_pos = rotation_matrix @ rel_pos
        # 重心に戻す
        atom.position = rotated_pos + com

    # --- 回転後の重心を再計算 ---
    com_rotated = adsorbate_copy.get_center_of_mass()

    # --- 目標位置への平行移動 ---
    target_position = target_pos + normal * distance
    translation = target_position - com_rotated

    for atom in adsorbate_copy:
        atom.position += translation

    # --- 基板と吸着分子を結合 ---
    combined = substrate + adsorbate_copy

    # --- 基板マスクの拡張（存在する場合） ---
    if "is_substrate" in substrate.arrays:
        old_mask = substrate.arrays["is_substrate"].astype(bool)
        new_mask = np.concatenate([old_mask, np.zeros(len(adsorbate_copy), dtype=bool)])
        combined.set_array("is_substrate", new_mask)

    return combined


def _compute_rotation_matrix(
    v_from: NDArray[np.float64], v_to: NDArray[np.float64]
) -> NDArray[np.float64]:
    """
    ベクトル v_from を v_to に回転させる回転行列を Rodrigues の回転公式で計算する。

    Args:
        v_from (NDArray[np.float64]): 回転前のベクトル（正規化済みを想定）。
        v_to (NDArray[np.float64]): 回転後のベクトル（正規化済みを想定）。

    Returns:
        NDArray[np.float64]: 3x3 回転行列。

    Note:
        v_from と v_to が同方向または反対方向の場合は特別に処理される。
    """
    # --- ベクトルの正規化 ---
    v_from = v_from / np.linalg.norm(v_from)
    v_to = v_to / np.linalg.norm(v_to)

    # --- 内積で角度関係をチェック ---
    dot_product = np.dot(v_from, v_to)

    # --- 同方向の場合（ほぼ平行） ---
    if np.isclose(dot_product, 1.0, atol=1e-10):
        return np.eye(3)

    # --- 反対方向の場合（ほぼ反平行） ---
    if np.isclose(dot_product, -1.0, atol=1e-10):
        # 任意の垂直軸を見つけて180度回転
        # v_from に垂直なベクトルを作成
        if abs(v_from[0]) < 0.9:
            perpendicular = np.array([1.0, 0.0, 0.0])
        else:
            perpendicular = np.array([0.0, 1.0, 0.0])

        # グラム・シュミット法で正規直交化
        perpendicular = perpendicular - np.dot(perpendicular, v_from) * v_from
        perpendicular = perpendicular / np.linalg.norm(perpendicular)

        # 180度回転行列: R = 2 * n * n^T - I （nは回転軸）
        return 2.0 * np.outer(perpendicular, perpendicular) - np.eye(3)

    # --- 一般的な場合：Rodrigues の回転公式 ---
    # 回転軸は外積で求める
    k = np.cross(v_from, v_to)
    k = k / np.linalg.norm(k)

    # 回転角度
    theta = np.arccos(np.clip(dot_product, -1.0, 1.0))

    # Rodrigues の回転公式: R = I + sin(θ) * [k]× + (1 - cos(θ)) * [k]×^2
    # [k]× は k の外積行列（歪対称行列）
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])

    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)

    return R


def place_adsorbate_on_surface(
    substrate: Atoms,
    adsorbate: Atoms,
    target_atom: int | Atom | None,
    height: float,
    position: Literal["top", "bridge", "hollow"],
    *,
    rotation_deg: tuple[float, float, float] | None = None,
    align_vector: Sequence[float] | None = None,
    rotate_about: Literal["com", "cog"] = "com",
    separate_layers_decimals: int = 4,
    allow_search_surface_atom: bool = True,
    inplace: bool = False,
    use_substrate_mask: Literal["auto", True, False] = "auto",
) -> Atoms:
    """
    指定した構造表面に、吸着分子を配置する。add_adsorbate()の高性能なラッパー関数。

    Args:
        substrate (ase.Atoms): ベースとなる原子構造。
        adsorbate (ase.Atoms): 配置する吸着分子。
        target_atom (int | ase.Atom | None): 基準となるベース上の原子。
            インデックスまたはAtomオブジェクトを指定可能。Noneの場合は重心に最も近い原子を探す。
        height (float): 吸着分子の高さ [Å]。
        position (Literal["top", "bridge", "hollow"]): 吸着分子の位置。
        rotation_deg (tuple[float, float, float] | None, optional): 吸着分子の回転角度 [度]。
            (rx, ry, rz) の形式で、X軸→Y軸→Z軸の順にオイラー角回転を適用します。
            Noneの場合は回転なし。デフォルトはNone。
        align_vector (Sequence[float] | None, optional): 吸着分子を整列させる方向ベクトル。
            指定した場合、このベクトルをグローバル +z 軸に整列させます。
            Noneの場合は整列なし。デフォルトはNone。
        rotate_about (Literal["com", "cog"], optional): 回転の中心。
            "com": 質量中心（Center of Mass）を中心に回転。
            "cog": 幾何中心（Center of Geometry）を中心に回転。
            デフォルトは "com"。
        separate_layers_decimals (int): 層を分割する際の小数点以下の桁数。
        allow_search_surface_atom (bool): target_atomが表面に存在しない場合、target_atomのxyに近い表面原子を探すかどうか。Falseの場合はエラーを返す。
        inplace (bool): もとの構造を置き換えるかどうか。
        use_substrate_mask (Literal["auto", True, False], optional): 基板マスクの使用設定。
            "auto": substrate.arrays に "is_substrate" が存在する場合、基板原子のみで層検出と高さ基準を決定。
            True: 基板マスクを使用（存在しない場合は全原子）。
            False: マスクを無視して全原子を対象。
            デフォルトは "auto"。複数の吸着分子を配置する場合、set_substrate_mask_all() で
            事前に基板マスクを設定し、このパラメータを "auto" または True にすることで、
            既存の吸着分子の影響を受けずに正しく配置できます。

    Returns:
        ase.Atoms: ベースと配置済み吸着分子を結合した構造。

    Raises:
        ValueError: position が 'top'、'bridge'、'hollow' 以外の場合。
        ValueError: target_atomが表面に存在しない場合、allow_search_surface_atomがFalseの場合。
        ValueError: 隣接原子が存在せず、bridgeまたはhollowの位置を決定できない場合。
        ValueError: 共通して隣接する原子が存在せず、hollowの位置を決定できない場合。
        ValueError: align_vectorがゼロベクトルの場合、またはrotate_aboutが不正な値の場合。
        TypeError: target_atom の型が不正な場合。
        IndexError: 指定されたインデックスが範囲外の場合。

    Note:
        - 回転は配置前に適用されます（整列回転 → オイラー角回転 → 位置決定）。
        - align_vectorとrotation_degは併用可能です。この場合、先にalign_vector整列、次にrotation_deg回転が適用されます。
    """
    # ---(準備)target_atomを、Atom|Noneにする
    if isinstance(target_atom, int):
        target_atom_index: int = target_atom
        target_atom: Atom | None = substrate[target_atom]
    elif isinstance(target_atom, Atom):
        target_atom_index: int = target_atom.index
        target_atom: Atom | None = target_atom
    elif target_atom is None:
        target_atom: Atom | None = None
        target_atom_index: int = None
    else:
        raise TypeError("target_atom は int または ase.Atom または None を指定してください。")

    # ---表面原子を取得する（基板マスク適用）
    layers: list[list[Atom]] = separate_layers(
        substrate,
        decimals=separate_layers_decimals,
        return_type="atoms",
        use_substrate_mask=use_substrate_mask,
    )
    top_layer: list[Atom] = layers[-1]
    top_layer_indices: list[int] = [atom.index for atom in top_layer]

    # ---target_atomを決定する
    if target_atom is None:
        # target_atomがNoneの場合、重心に最も近い原子を探す
        target_atom = find_central_atom(top_layer)
        target_atom_index = target_atom.index
    elif target_atom_index not in top_layer_indices:  # 含まれるか判定は、インデックスで行う
        # target_atomが表面に存在しない場合、target_atomのxyに近い表面原子を探す
        if allow_search_surface_atom:
            nearest_surface_atom = min(top_layer, key=lambda x: np.linalg.norm(
                x.position[:2] - target_atom.position[:2]))
            target_atom = nearest_surface_atom
        else:
            raise ValueError("指定された原子は表面に存在しません。")

    target_atom: Atom = target_atom  # ここまでで、target_atomはAtomオブジェクトになっている

    # ---positionから、xy座標を取得する
    if position == "top":
        position_xy = target_atom.position[:2]
    elif position in ["bridge", "hollow"]:
        # 表面上で、target_atomの隣接原子を探す
        target_atom_neighbors = get_neighbors(substrate,
                                              target_atom,
                                              return_type="atoms")
        # 表面上の隣接原子のみを残す
        target_atom_neighbors = list(
            filter(
                lambda x: np.isclose(
                    x.position[2],
                    target_atom.position[2]),
                target_atom_neighbors))
        if len(target_atom_neighbors) == 0:
            raise ValueError("隣接原子が存在せず、bridgeまたはhollowの位置を決定できません。")
        # ---bridgeの場合、1番目の隣接原子との中点を探す
        if position == "bridge":
            position_xy = (target_atom.position[:2] +
                           target_atom_neighbors[0].position[:2]) / 2
        elif position == "hollow":
            # ---hollowの場合、もう1つの原子を探す
            # target_atomも隣接原子も、共通して隣接する原子を探す
            # お互いの友だちだったら、そこは3人グループだよね、という考え方
            for neighbor in target_atom_neighbors:
                neighbor_neighbors = get_neighbors(substrate,
                                                   neighbor,
                                                   return_type="atoms")
                # (表面上の隣接原子のみを残す)
                neighbor_neighbors = list(
                    filter(
                        lambda x: np.isclose(x.position[2],
                                             neighbor.position[2]),
                        neighbor_neighbors))
                # 積集合を求める
                common_neighbor_indices = set([a.index for a in target_atom_neighbors]) & \
                    set([a.index for a in neighbor_neighbors])

                if len(common_neighbor_indices) > 0:
                    common_neighbor_index = list(common_neighbor_indices)[0]
                    common_neighbor = substrate[common_neighbor_index]
                    position_xy = (target_atom.position[:2] +
                                   neighbor.position[:2] +
                                   common_neighbor.position[:2]) / 3
                    break
            else:
                raise ValueError("共通して隣接する原子が存在せず、hollowの位置を決定できません。")
    else:
        raise ValueError("position は 'top'、'bridge'、'hollow' を指定してください。")

    # ---基板上面のz座標を計算（基板のみから）
    z_top = max(a.position[2] for a in top_layer)

    # ---出力構造の準備
    if inplace:
        out = substrate
    else:
        out = substrate.copy()

    # ---吸着分子を手動配置（高さ基準を基板のみに依存）
    ads = adsorbate.copy()
    
    # --- 回転処理の適用（配置前） ---
    # 回転中心の決定
    if rotate_about == "com":
        rotation_center = ads.get_center_of_mass()
    elif rotate_about == "cog":
        rotation_center = ads.get_positions().mean(axis=0)
    else:
        raise ValueError(f"rotate_aboutは 'com' または 'cog' を指定してください。指定値: {rotate_about}")
    
    # 1. align_vector による整列回転
    if align_vector is not None:
        align_vec = np.array(align_vector, dtype=float)
        if align_vec.shape != (3,):
            raise ValueError(f"align_vectorは3次元ベクトルを指定してください。現在の形状: {align_vec.shape}")
        if np.allclose(align_vec, 0):
            raise ValueError("align_vectorはゼロベクトルにできません。")
        
        # align_vector を +z 軸に整列させる回転行列を計算
        z_axis = np.array([0.0, 0.0, 1.0])
        R_align = _compute_rotation_matrix(align_vec, z_axis)
        
        # 回転中心を基準に回転を適用
        for atom in ads:
            rel_pos = atom.position - rotation_center
            rotated_pos = R_align @ rel_pos
            atom.position = rotated_pos + rotation_center
        
        # 回転後の中心を再計算
        if rotate_about == "com":
            rotation_center = ads.get_center_of_mass()
        else:
            rotation_center = ads.get_positions().mean(axis=0)
    
    # 2. rotation_deg によるオイラー角回転（XYZ順）
    if rotation_deg is not None:
        if len(rotation_deg) != 3:
            raise ValueError(f"rotation_degは3つの値 (rx, ry, rz) を指定してください。指定値: {rotation_deg}")
        
        rx_deg, ry_deg, rz_deg = rotation_deg
        rx_rad = np.deg2rad(rx_deg)
        ry_rad = np.deg2rad(ry_deg)
        rz_rad = np.deg2rad(rz_deg)
        
        # X軸回りの回転行列
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(rx_rad), -np.sin(rx_rad)],
            [0, np.sin(rx_rad), np.cos(rx_rad)]
        ])
        
        # Y軸回りの回転行列
        Ry = np.array([
            [np.cos(ry_rad), 0, np.sin(ry_rad)],
            [0, 1, 0],
            [-np.sin(ry_rad), 0, np.cos(ry_rad)]
        ])
        
        # Z軸回りの回転行列
        Rz = np.array([
            [np.cos(rz_rad), -np.sin(rz_rad), 0],
            [np.sin(rz_rad), np.cos(rz_rad), 0],
            [0, 0, 1]
        ])
        
        # 合成回転行列（Z * Y * X の順）
        R_euler = Rz @ Ry @ Rx
        
        # 回転中心を基準に回転を適用
        for atom in ads:
            rel_pos = atom.position - rotation_center
            rotated_pos = R_euler @ rel_pos
            atom.position = rotated_pos + rotation_center
    
    # XY位置合わせ: 吸着分子の重心XYを目標XY位置に移動
    com_xy = ads.get_center_of_mass()[:2]
    shift_xy = position_xy - com_xy
    ads.positions[:, :2] += shift_xy
    
    # Z位置合わせ: 吸着分子の最下点が z_top + height になるように移動
    z_min = ads.positions[:, 2].min()
    ads.positions[:, 2] += (z_top + height - z_min)

    # ---構造を結合
    out = out + ads

    # ---基板マスクの拡張（存在する場合）
    if "is_substrate" in substrate.arrays:
        old_mask = substrate.arrays["is_substrate"].astype(bool)
        new_mask = np.concatenate([old_mask, np.zeros(len(ads), dtype=bool)])
        out.set_array("is_substrate", new_mask)

    # ---返却する
    return out
