from datetime import datetime
import logging
from typing import Optional
from ase import Atoms, Atom
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer


class ConditionalLogger:
    """
    ログ出力を条件によって制御するラッパークラス。

    enabled=False の場合、すべてのログ出力を無効化する。
    try-delta-g.py から移植。
    """

    def __init__(self, base_logger, enabled: bool = True):
        self.base_logger = base_logger
        self.enabled = enabled

    def __getattr__(self, name):
        """
        logger.info(), logger.warning() などの呼び出しを透過的に処理する。
        enabled=False の場合は何もしない関数を返す。
        """
        if not self.enabled:
            # ログが無効な場合は何もしない関数を返す
            return lambda *args, **kwargs: None

        # ログが有効な場合は元のloggerのメソッドを返す
        return getattr(self.base_logger, name)


def setup_logger(prefix: str = "calc"):
    """デバッグ用ログの設定

    prefixごとに独立したLoggerを作成し、重複するハンドラ付与を防ぐ。
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{prefix}_{timestamp}.log"

    # prefixごとに独立したログ名を作成
    logger_name = f"ase_toolbox.{prefix}"
    base_logger = logging.getLogger(logger_name)

    # 既にハンドラが設定されている場合は重複設定を避ける
    if base_logger.handlers:
        return ConditionalLogger(base_logger, enabled=True)

    # ログレベル設定
    base_logger.setLevel(logging.DEBUG)

    # フォーマッタ作成
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # ファイルハンドラ作成
    file_handler = logging.FileHandler(log_filename, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    base_logger.addHandler(file_handler)

    # コンソールハンドラ作成
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    base_logger.addHandler(console_handler)

    # 親ロガーへの伝播を無効化（重複出力を防ぐ）
    base_logger.propagate = False

    base_logger.info(f"デバッグログファイル: {log_filename}")
    return ConditionalLogger(base_logger, enabled=True)


def ensure_logger(
    prefix: str, enable_logging: bool, logger: Optional[ConditionalLogger] = None
) -> ConditionalLogger:
    """
    ロガーの初期化を統一的に処理する。

    既存のロガーがあればそれを使用し、なければ新規作成する。

    Args:
        prefix (str): ログファイルのプレフィックス。
        enable_logging (bool): ログ出力を有効にするかどうか。
        logger (Optional[ConditionalLogger]): 既存のロガー。Noneの場合は新規作成。

    Returns:
        ConditionalLogger: 設定済みのロガー。
    """
    if logger is None:
        if enable_logging:
            logger = setup_logger(prefix=prefix)
        else:
            logger = ConditionalLogger(None, enabled=False)

    return logger


def optimize_and_get_energy(
    atoms: Atoms,
    calculator: Calculator,
    optimizer_cls: type[Optimizer],
    fmax: float,
    maxsteps: int,
    label: str,
    logger: ConditionalLogger,
    *,
    copy_atoms: bool = True,
) -> float:
    """
    構造最適化をおこない、最適化後のエネルギーを取得する。

    Args:
        atoms: 最適化する原子構造。
        calculator: 使用する計算機。
        optimizer_cls: 最適化アルゴリズムのクラス。
        fmax: 収束判定の力の閾値 [eV/Å]。
        maxsteps: 最大最適化ステップ数。
        label: ログ出力用のラベル。
        logger: ロガー。
        copy_atoms: Trueの場合、atomsをコピーして処理。Falseの場合は直接変更。

    Returns:
        float: 最適化後のポテンシャルエネルギー [eV]。
    """
    # 作業用コピーを作成
    if copy_atoms:
        work_atoms = atoms.copy()
    else:
        work_atoms = atoms

    # 計算機を設定
    work_atoms.calc = calculator

    logger.info(f"--- {label} 処理開始 ---")
    logger.info(f"原子数: {len(work_atoms)}")
    logger.info(f"組成: {work_atoms.symbols}")

    # 制約情報のログ出力
    constraints = work_atoms.constraints
    if constraints:
        logger.info(f"制約条件: {len(constraints)} 個")
        for i, constraint in enumerate(constraints):
            constraint_name = type(constraint).__name__
            logger.info(f"  制約{i}: {constraint_name}")
            # 固定原子数の情報（FixAtomsの場合）
            if hasattr(constraint, "index"):
                n_fixed = (
                    len(constraint.index) if hasattr(constraint.index, "__len__") else 1
                )
                logger.info(f"    固定原子数: {n_fixed}")
    else:
        logger.info("制約条件: なし")

    # 初期エネルギー取得
    e_initial = work_atoms.get_potential_energy()
    logger.info(f"初期エネルギー: {e_initial:.6f} eV")

    # 構造最適化実行
    logger.info("構造最適化開始")
    opt = optimizer_cls(work_atoms)
    opt.run(fmax=fmax, steps=maxsteps)
    logger.info("構造最適化完了")

    # 最適化後エネルギー取得
    e_final = work_atoms.get_potential_energy()
    e_change = e_final - e_initial
    logger.info(f"最適化後エネルギー: {e_final:.6f} eV")
    logger.info(f"エネルギー変化: {e_change:.6f} eV")
    logger.info(f"--- {label} 処理完了 ---")

    return e_final


def resolve_target_indices(
    base_atoms: Atoms, target: int | Atom | Atoms | list[int] | list[Atom]
) -> list[int]:
    """
    様々な形式のターゲット指定を、インデックスリストに正規化する。

    Args:
        base_atoms: 基準となる原子構造。
        target: ターゲット指定。以下の形式に対応：
            - int: 単一原子のインデックス
            - Atom: 単一原子オブジェクト
            - Atoms: 構造全体（全原子を選択）
            - list[int]: 複数原子のインデックスリスト
            - list[Atom]: 複数原子オブジェクトのリスト

    Returns:
        list[int]: ターゲットのインデックスリスト。

    Raises:
        TypeError: target の型が不正な場合。
        IndexError: 指定されたインデックスが範囲外の場合。
        ValueError: 指定されたAtomがbase_atoms内に存在しない場合。
    """
    target_indices: list[int] = []

    # --- int: 単一原子インデックス ---
    if isinstance(target, int):
        if not (0 <= target < len(base_atoms)):
            raise IndexError(
                f"インデックス{target}は範囲外です（0-{len(base_atoms)-1}）。"
            )
        target_indices = [target]

    # --- Atom: 単一原子オブジェクト ---
    elif isinstance(target, Atom):
        try:
            index = target.index
            if not (0 <= index < len(base_atoms)):
                raise IndexError(
                    f"インデックス{index}は範囲外です（0-{len(base_atoms)-1}）。"
                )
            target_indices = [index]
        except (AttributeError, ValueError):
            raise ValueError("指定されたAtomはbase_atoms内に存在しません。")

    # --- Atoms: 構造全体 ---
    elif isinstance(target, Atoms):
        target_indices = list(range(len(base_atoms)))

    # --- list: 複数指定 ---
    elif isinstance(target, list):
        if not target:  # 空リスト
            target_indices = []
        elif all(isinstance(item, int) for item in target):
            # list[int]の場合
            for idx in target:
                if not (0 <= idx < len(base_atoms)):
                    raise IndexError(
                        f"インデックス{idx}は範囲外です（0-{len(base_atoms)-1}）。"
                    )
            target_indices = list(target)
        elif all(isinstance(item, Atom) for item in target):
            # list[Atom]の場合
            for atom in target:
                try:
                    index = atom.index
                    if not (0 <= index < len(base_atoms)):
                        raise IndexError(
                            f"インデックス{index}は範囲外です（0-{len(base_atoms)-1}）。"
                        )
                    target_indices.append(index)
                except (AttributeError, ValueError):
                    raise ValueError("指定されたAtomはbase_atoms内に存在しません。")
        else:
            raise TypeError("リストの要素は全てintまたは全てAtomである必要があります。")

    else:
        raise TypeError(
            "targetはint、Atom、Atoms、list[int]、またはlist[Atom]を指定してください。"
        )

    return target_indices


def sanitize_atoms_for_xyz_write(atoms: Atoms) -> Atoms:
    """
    XYZ形式で安全に書き出すために、Atomsオブジェクトをクリーンアップする。

    長さが不一致なper-atom配列や、シリアライズできない情報を除去し、
    最小限の情報（元素記号、座標、セル、PBC）のみを保持した新しいAtomsオブジェクトを作成します。

    Args:
        atoms (ase.Atoms): クリーンアップ対象の原子構造。

    Returns:
        ase.Atoms: クリーンアップされた原子構造。以下の情報のみを保持：
            - 元素記号（symbols）
            - 座標（positions）
            - セル（cell）
            - 周期境界条件（pbc）
            - タグ（tags、長さが一致する場合のみ）
            - info辞書の単純型エントリ（str, int, float, bool のみ）

    Examples:
        >>> from ase import Atoms
        >>> import numpy as np
        >>> # 不一致な配列を持つ構造を作成（例：誤って設定）
        >>> atoms = Atoms("H2O", positions=[[0,0,0],[1,0,0],[0,1,0]])
        >>> # クリーンアップしてXYZ書き出し
        >>> clean_atoms = sanitize_atoms_for_xyz_write(atoms)
        >>> from ase.io import write
        >>> write("output.xyz", clean_atoms)

    Note:
        - この関数は元のAtomsオブジェクトを変更せず、新しいオブジェクトを返します。
        - per-atom配列（arrays）で長さが原子数と一致しないものは除外されます。
        - info辞書は単純な型（str, int, float, bool）のみが保存されます。
        - 計算機（calculator）やその他の複雑なオブジェクトは保存されません。
    """
    n = len(atoms)

    # --- 基本情報で新しいAtomsを作成 ---
    clean = Atoms(
        symbols=atoms.get_chemical_symbols(),
        positions=atoms.get_positions(),
        cell=atoms.get_cell(),
        pbc=atoms.get_pbc(),
    )

    # --- タグの保存（長さが一致する場合のみ） ---
    try:
        tags = atoms.get_tags()
        if tags is not None and len(tags) == n:
            clean.set_tags(tags)
    except Exception:
        # タグの取得や設定に失敗した場合は無視
        pass

    # --- info辞書の単純型エントリのみを保存 ---
    try:
        if hasattr(atoms, "info") and isinstance(atoms.info, dict):
            for k, v in atoms.info.items():
                if isinstance(v, (str, int, float, bool)):
                    clean.info[k] = v
    except Exception:
        # info処理に失敗した場合は無視
        pass

    return clean