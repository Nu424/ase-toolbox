"""
# HandleAtoms.py
ASEのAtomsオブジェクトを操作・処理するための関数をまとめたファイル

## 関数一覧
- move_atoms(): 原子を指定方向に指定距離だけ移動する
- fix_layers(): (平面用)層を固定する
- substitute_elements(): 原子を置き換える
- compute_surface_normal(): (クラスター用)法線ベクトルを計算する
- place_adsorbate_along_normal(): (クラスター用)クラスターにくっつける
- place_adsorbate_on_surface(): (表面用)表面にくっつける
"""

from typing import Literal, Sequence, Mapping, Optional
from ase import Atoms, Atom
from ase.build import add_adsorbate
from ase.constraints import FixAtoms
import numpy as np
from numpy.typing import NDArray
from .FindAtoms import (
    find_central_atom,
    separate_layers,
    get_neighbors_with_coordination_condition,
    get_neighbors,
)
from .util import ConditionalLogger, ensure_logger, resolve_target_indices


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
    separate_layers_decimals: int = 4,
    allow_search_surface_atom: bool = True,
    inplace: bool = False,
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
        separate_layers_decimals (int): 層を分割する際の小数点以下の桁数。
        allow_search_surface_atom (bool): target_atomが表面に存在しない場合、target_atomのxyに近い表面原子を探すかどうか。Falseの場合はエラーを返す。
        inplace (bool): もとの構造を置き換えるかどうか。

    Returns:
        ase.Atoms: ベースと配置済み吸着分子を結合した構造。

    Raises:
        ValueError: position が 'top'、'bridge'、'hollow' 以外の場合。
        ValueError: target_atomが表面に存在しない場合、allow_search_surface_atomがFalseの場合。
        ValueError: 隣接原子が存在せず、bridgeまたはhollowの位置を決定できない場合。
        ValueError: 共通して隣接する原子が存在せず、hollowの位置を決定できない場合。
        TypeError: target_atom の型が不正な場合。
        IndexError: 指定されたインデックスが範囲外の場合。
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

    # ---表面原子を取得する
    layers: list[list[Atom]] = separate_layers(substrate,
                                               decimals=separate_layers_decimals,
                                               return_type="atoms")
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

    # ---吸着させる
    if inplace:
        new_substrate = substrate
    else:
        new_substrate = substrate.copy()

    add_adsorbate(new_substrate, adsorbate,
                  height=height, position=position_xy)

    # ---返却する
    return new_substrate
