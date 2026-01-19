"""
# AlignSlab.py
スラブ構造の傾きを補正する関数群。

## 関数一覧
- fit_plane_normal_svd(): 点群から平面法線を推定
- align_slab(): 1段階/複数段階で傾きを補正（selectorは単体/リスト対応）
- select_points_seed_neighbors(): 中央付近のシード原子と近傍を選ぶ
- select_points_middle_layer(): 真ん中あたりの層を返す
- select_points_middle_layer_farthest_k(): 中間層から遠い点を選ぶ
"""

from typing import Callable, Iterable, Literal, Sequence

import numpy as np
from numpy.typing import NDArray
from ase import Atoms
from ase.neighborlist import neighbor_list

from .FindAtoms import separate_layers
from .util import ConditionalLogger, ensure_logger

PlanePointSelector = Callable[[Atoms], Sequence[int]]
_Z_AXIS = np.array([0.0, 0.0, 1.0])

# ----------
# ---内部関数
# ----------
def _normalize_indices(indices: Iterable[int], n_atoms: int) -> list[int]:
    """インデックス配列を正規化する。

    重複を除去し、範囲外インデックスを検出する。

    Args:
        indices (Iterable[int]): 入力インデックス列。
        n_atoms (int): 原子数。

    Returns:
        list[int]: 正規化されたインデックス。

    Raises:
        IndexError: 範囲外のインデックスが含まれる場合。
    """
    # --- 正規化と重複排除 ---
    unique_indices: list[int] = []
    seen: set[int] = set()
    for idx in indices:
        index = int(idx)
        # --- 範囲チェック ---
        if index < 0 or index >= n_atoms:
            raise IndexError(f"インデックス {index} は範囲外です（0-{n_atoms - 1}）。")
        if index in seen:
            continue
        seen.add(index)
        unique_indices.append(index)
    return unique_indices


def _normalize_selectors(
    selector: PlanePointSelector | Sequence[PlanePointSelector] | None,
    default_selectors: Sequence[PlanePointSelector],
) -> list[PlanePointSelector]:
    """selector入力を正規化し、セレクタのリストを返す。

    Args:
        selector (PlanePointSelector | Sequence[PlanePointSelector] | None):
            単体セレクタ、またはセレクタ列。
        default_selectors (Sequence[PlanePointSelector]): selectorがNoneのときに使う既定値。

    Returns:
        list[PlanePointSelector]: セレクタのリスト。

    Raises:
        TypeError: selectorが不正な型の場合。
        ValueError: セレクタが空リストの場合。
    """
    # --- selector の形を正規化 ---
    if selector is None:
        selectors = list(default_selectors)
    elif callable(selector):
        selectors = [selector]
    elif isinstance(selector, Sequence):
        selectors = list(selector)
    else:
        raise TypeError("selector は関数または関数のリストを指定してください。")

    if len(selectors) == 0:
        raise ValueError("selector は1つ以上の関数を指定してください。")
    if not all(callable(item) for item in selectors):
        raise TypeError("selector の各要素は呼び出し可能な関数である必要があります。")

    return selectors


def fit_plane_normal_svd(points: NDArray[np.floating]) -> NDArray[np.float64]:
    """点群に最小二乗でフィットする平面法線を推定する。

    Args:
        points (NDArray[np.floating]): 形状 (N, 3) の点群座標。

    Returns:
        NDArray[np.float64]: 正規化された法線ベクトル（形状: (3,)）。

    Raises:
        ValueError: 形状が不正、点数不足、または点群が共線で平面が定義できない場合。
    """
    # --- 入力検証 ---
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points は (N, 3) 形状の配列を指定してください。")
    if points.shape[0] < 3:
        raise ValueError("平面定義には3点以上が必要です。")

    # --- 重心で中心化 ---
    centroid = points.mean(axis=0)
    centered = points - centroid

    # --- 共線性チェック ---
    if np.linalg.matrix_rank(centered) < 2:
        raise ValueError("点群が共線のため平面を定義できません。")

    # --- SVD で法線を取得 ---
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]
    norm = float(np.linalg.norm(normal))
    if norm == 0.0:
        raise ValueError("法線ベクトルがゼロになりました。点群を確認してください。")
    return normal / norm


def _ensure_positive_z(normal: NDArray[np.float64]) -> NDArray[np.float64]:
    """法線ベクトルの向きを +z 方向に揃える。

    Args:
        normal (NDArray[np.float64]): 法線ベクトル。

    Returns:
        NDArray[np.float64]: +z 方向に揃えた法線ベクトル。
    """
    # --- z方向の符号を揃える ---
    if normal[2] < 0.0:
        return -normal
    return normal


def _angle_to_z(normal: NDArray[np.float64]) -> float:
    """法線ベクトルと +z 軸の角度（rad）を返す。

    Args:
        normal (NDArray[np.float64]): 法線ベクトル（正規化済みを想定）。

    Returns:
        float: +z 軸との角度（ラジアン）。
    """
    # --- 内積から角度を算出 ---
    # ベクトル→(dot)→内積→(arccos)→角度…という流れ
    dot = float(np.dot(normal, _Z_AXIS))
    return float(np.arccos(np.clip(dot, -1.0, 1.0)))


def _recenter_to_cell(atoms: Atoms) -> None:
    """原子の重心をセル中心に移動する。

    Args:
        atoms (Atoms): 対象構造。
    """
    # --- 重心とセル中心を計算 ---
    atom_center = atoms.positions.mean(axis=0)
    cell_center = atoms.get_cell().sum(axis=0) / 2.0
    # --- 平行移動 ---
    atoms.translate(cell_center - atom_center)


def _apply_alignment(
    atoms: Atoms,
    selector: PlanePointSelector,
    *,
    rotate_cell: bool,
    recenter_to_cell: bool,
    logger: ConditionalLogger,
    label: str,
) -> float | None:
    """点群選択→平面フィット→回転を適用する。

    失敗した場合は警告ログを出し、None を返す。

    Args:
        atoms (Atoms): 対象構造。
        selector (PlanePointSelector): フィット点の選択関数。
        rotate_cell (bool): セルも回転させるか。
        recenter_to_cell (bool): 回転後にセル中心へ平行移動するか。
        logger (ConditionalLogger): ロガー。
        label (str): ログ用ラベル。

    Returns:
        float | None: 回転角（rad）。失敗時は None。
    """
    # --- 点群の取得 ---
    try:
        indices = _normalize_indices(selector(atoms), len(atoms))
    except Exception as exc:
        logger.warning(f"{label}: 点群選択に失敗しました: {exc}")
        return None

    # --- 点数チェック ---
    if len(indices) < 3:
        logger.warning(
            f"{label}: 平面定義に必要な点数が不足しています（取得: {len(indices)}点）。"
        )
        return None

    # --- 平面フィット ---
    points = atoms.positions[indices]
    try:
        normal = fit_plane_normal_svd(points)
    except ValueError as exc:
        logger.warning(f"{label}: 平面フィットに失敗しました: {exc}")
        return None

    # --- 回転の適用 ---
    normal = _ensure_positive_z(normal)
    angle_rad = _angle_to_z(normal)

    if angle_rad > 0.0:
        atoms.rotate(normal, _Z_AXIS, center=(0.0, 0.0, 0.0), rotate_cell=rotate_cell)

    # --- 平行移動 ---
    if recenter_to_cell:
        _recenter_to_cell(atoms)

    # --- ログ出力 ---
    logger.info(
        f"{label}: 点数={len(indices)}, 回転角={np.degrees(angle_rad):.4f} deg"
    )
    return angle_rad

# ----------
# ---平面フィッティング用の原子選択に使用する関数
# ----------
def select_points_seed_neighbors(
    atoms: Atoms,
    *,
    cutoff: float = 3.0,
    z_diff_threshold: float = 1.0,
    center_xy_tolerance: float = 1.0,
    max_points: int | None = None,
) -> list[int]:
    """中央付近のシード原子と近傍から平面フィット点を選ぶ。

    Args:
        atoms (Atoms): 対象構造。
        cutoff (float): 隣接探索のカットオフ距離 [Å]。
        z_diff_threshold (float): シード原子とのz差の許容幅 [Å]。
        center_xy_tolerance (float): xy重心からの許容距離 [Å]。
        max_points (int | None): 選択点数の上限。Noneなら全候補を使用。

    Returns:
        list[int]: フィットに使う原子インデックス。

    Raises:
        ValueError: 各種パラメータが不正な場合。
    """
    # --- 入力検証 ---
    if cutoff <= 0:
        raise ValueError("cutoff は正の値を指定してください。")
    if z_diff_threshold < 0:
        raise ValueError("z_diff_threshold は0以上の値を指定してください。")
    if center_xy_tolerance < 0:
        raise ValueError("center_xy_tolerance は0以上の値を指定してください。")

    # --- 中央列（カラム）を探索 ---
    positions = atoms.positions
    pos_xy = positions[:, :2]
    center_xy = np.mean(pos_xy, axis=0)
    dists_sq = np.sum((pos_xy - center_xy) ** 2, axis=1)
    min_dist_sq = float(np.min(dists_sq)) # 最も中心に近い距離を取得

    # やや誤差を許容して、中心に近い原子を抽出する
    tolerance_sq = center_xy_tolerance**2
    central_column_indices = np.where(dists_sq <= (min_dist_sq + tolerance_sq))[0]
    if len(central_column_indices) == 0:
        return []

    # --- シード原子の決定（中央列の最上部） ---
    # 中央列の中で、最もZ座標が高い原子をシード原子とする
    best_local_idx = int(np.argmax(positions[central_column_indices, 2]))
    seed_idx = int(central_column_indices[best_local_idx]) # グローバルなインデックスに変換
    seed_z = float(positions[seed_idx, 2])

    # --- 近傍原子の取得 ---
    # シード原子(1つ)の近傍原子を取得する
    i_list, j_list = neighbor_list("ij", atoms, cutoff)
    neighbors_indices = j_list[i_list == seed_idx]

    # --- 近傍から候補を抽出 ---
    # 近くの原子との距離を改めてチェックし、良さそうであれば候補に入れる
    candidates: list[int] = [seed_idx]
    for n_idx in neighbors_indices:
        z_diff = abs(float(positions[n_idx, 2]) - seed_z)
        if z_diff <= z_diff_threshold:
            candidates.append(int(n_idx))

    # --- 重複除去と点数制限 ---
    # 重複除去と点数制限を行う
    candidates = _normalize_indices(candidates, len(atoms))
    if max_points is not None and len(candidates) > max_points:
        z_diffs = np.abs(positions[candidates, 2] - seed_z)
        sorted_order = np.argsort(z_diffs)
        candidates = [candidates[int(k)] for k in sorted_order[:max_points]]

    # 最終的に、シード+近傍原子処理で求められた、平面フィッティングに使用可能な原子インデックスを返す
    return candidates


def select_points_middle_layer(
    atoms: Atoms,
    *,
    z_tolerance: float = 0.3,
    use_substrate_mask: Literal["auto", True, False] = False,
) -> list[int]:
    """真ん中あたりの層の原子インデックスを返す。

    層数が偶数の場合は、上側寄りの中間層を返します。

    Args:
        atoms (Atoms): 対象構造。
        z_tolerance (float): 同一層とみなすz許容幅 [Å]。
        use_substrate_mask (Literal["auto", True, False]): 基板マスクの使用設定。

    Returns:
        list[int]: 中間層の原子インデックス。層が見つからない場合は空リスト。
    """
    # --- 層を検出 ---
    layers = separate_layers(
        atoms,
        return_type="indices",
        z_tolerance=z_tolerance,
        sort_by_z=True,
        use_substrate_mask=use_substrate_mask,
    )
    if not layers:
        return []

    # --- 中間層を返却（偶数の場合は上側寄り） ---
    middle_index = len(layers) // 2
    return list(layers[middle_index])


def _farthest_point_sampling(
    positions_xy: NDArray[np.floating], k: int
) -> list[int]:
    """2次元点群のfarthest-point samplingでk点のローカルindexを返す。

    Args:
        positions_xy (NDArray[np.floating]): 2次元点群（形状: (N, 2)）。
        k (int): 選択点数。

    Returns:
        list[int]: ローカルindexのリスト。

    Raises:
        ValueError: k が1未満、または点数を超える場合。
    """
    # --- 入力検証 ---
    n_points = positions_xy.shape[0]
    if k <= 0:
        raise ValueError("k は1以上の値を指定してください。")
    if k > n_points:
        raise ValueError("k は点数以下である必要があります。")

    # --- 初期点（重心から最遠）を選択 ---
    centroid = positions_xy.mean(axis=0)
    dist_sq = np.sum((positions_xy - centroid) ** 2, axis=1)
    selected = [int(np.argmax(dist_sq))]

    # --- 反復的に最遠点を追加 ---
    while len(selected) < k:
        selected_pos = positions_xy[selected]
        # 全点(N)×選択済み点(M)の差分ベクトルを、ブロードキャストで一括生成（形状: N×M×2）
        diff = positions_xy[:, None, :] - selected_pos[None, :, :]
        dist_sq = np.sum(diff**2, axis=2)
        min_dist_sq = dist_sq.min(axis=1)
        min_dist_sq[selected] = -1.0
        next_idx = int(np.argmax(min_dist_sq))
        selected.append(next_idx)

    return selected


def select_points_middle_layer_farthest_k(
    atoms: Atoms,
    *,
    k: int = 6,
    z_tolerance: float = 0.3,
    use_substrate_mask: Literal["auto", True, False] = False,
) -> list[int]:
    """真ん中あたりの層から互いに離れたk点を選んで返す。

    Args:
        atoms (Atoms): 対象構造。
        k (int): 選択する点数（3以上）。
        z_tolerance (float): 同一層とみなすz許容幅 [Å]。
        use_substrate_mask (Literal["auto", True, False]): 基板マスクの使用設定。

    Returns:
        list[int]: 中間層から選ばれた原子インデックス。

    Raises:
        ValueError: k が3未満、または内部サンプリング条件が満たせない場合。
    """
    # --- 入力検証 ---
    if k < 3:
        raise ValueError("k は3以上の値を指定してください。")

    # --- 中間層の取得 ---
    layer_indices = select_points_middle_layer(
        atoms, z_tolerance=z_tolerance, use_substrate_mask=use_substrate_mask
    )
    if len(layer_indices) <= k:
        return layer_indices

    # --- 中間層から遠い点を選択 ---
    positions_xy = atoms.positions[layer_indices, :2]
    local_indices = _farthest_point_sampling(positions_xy, k)
    return [layer_indices[i] for i in local_indices]


# ----------
# ---傾き補正のメイン関数
# ----------
def align_slab(
    slab: Atoms,
    *,
    selector: PlanePointSelector | Sequence[PlanePointSelector] | None = None,
    rotate_cell: bool = False,
    recenter_to_cell: bool = True,
    max_iters: int = 1,
    angle_tolerance_deg: float | None = None,
    inplace: bool = False,
    logger: ConditionalLogger | None = None,
    enable_logging: bool = False,
) -> Atoms:
    """指定した点群選択で平面をフィットし、傾きを補正する。

    selector は単体関数または関数リストを受け付けます。
    Noneの場合、[select_points_seed_neighbors, select_points_middle_layer_farthest_k] を使用します。

    Args:
        slab (Atoms): 対象スラブ構造。
        selector (PlanePointSelector | Sequence[PlanePointSelector] | None):
            フィット点の選択関数（単体またはリスト）。
        rotate_cell (bool): セルも回転させるか。デフォルトはFalse。
        recenter_to_cell (bool): 回転後にセル中心へ平行移動するか。
        max_iters (int): 指定したセレクタ列を繰り返す回数。
        angle_tolerance_deg (float | None): 収束角度 [deg]。0.1など。調整角度がこの値以下なら収束と判定する。Noneなら判定なし。
        inplace (bool): Trueの場合はslabを直接変更する。
        logger (ConditionalLogger | None): ロガー。
        enable_logging (bool): ログ出力を有効にするか。

    Returns:
        Atoms: 傾き補正後の構造。

    Notes:
        - 処理は基本的に_apply_alignment関数で行われる(align_slab()はラッパーのような存在)
            - 平面フィッティングに使用する原子を選択する(selector())
            - 平面フィッティングをする(fit_plane_normal_svd())
            - 回転を適用する
            - 平行移動を適用する
            - ログを出力する
    """
    # --- 入力検証 ---
    if max_iters <= 0:
        raise ValueError("max_iters は1以上の値を指定してください。")

    # --- ログと作業用構造の準備 ---
    logger = ensure_logger("align_slab", enable_logging, logger)
    work = slab if inplace else slab.copy()

    # --- セレクタ列の準備（デフォルトを含む） ---
    default_selectors = [
        select_points_seed_neighbors,
        select_points_middle_layer_farthest_k,
    ]
    selectors = _normalize_selectors(selector, default_selectors)

    # --- 収束判定の閾値を準備 ---
    angle_tolerance_rad = (
        None if angle_tolerance_deg is None else np.deg2rad(angle_tolerance_deg)
    )

    # --- 指定セレクタ列で傾き補正を反復 ---
    for i in range(max_iters):
        max_angle = 0.0
        has_angle = False

        for stage_idx, stage_selector in enumerate(selectors, start=1):
            angle = _apply_alignment(
                work,
                stage_selector,
                rotate_cell=rotate_cell,
                recenter_to_cell=recenter_to_cell,
                logger=logger,
                label=f"stage{stage_idx}[{i + 1}]",
            )
            if angle is not None:
                has_angle = True
                max_angle = max(max_angle, angle)

        # --- 収束判定（指定時のみ） ---
        if angle_tolerance_rad is not None and has_angle:
            if max_angle <= angle_tolerance_rad:
                break

    return work
