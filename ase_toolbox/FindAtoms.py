"""
# FindAtoms.py
ASEのAtomsオブジェクトから、特定の原子を探すための関数をまとめたファイル

## 関数一覧
- find_atom_by_index(): 原子を、インデックスでピンポイントに指定する
- find_indices_by_symbol(): 指定した原子のインデックスを調べる
- get_neighbors(): 隣接原子を探す
- separate_layers(): (平面用)層別に分ける・層ごとのlistにする (is_substrate マスク対応)
- classify_surface_atoms(): (クラスター用)表面・内側を探す
- find_central_atom(): 重心に最も近い原子を探す
- get_appended_atom_indices(): くっつけた後の構造の中で、くっつけた原子のインデックスを知る
- get_neighbors_with_coordination_condition(): 自身の隣接原子数と同じ数だけ、隣接原子を持つ隣接原子を探す
"""

from typing import Literal, Optional, Sequence
from ase import Atoms, Atom
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
from .CalcValue import coordination_number


# インデックスでピンポイントに指定する
def find_atom_by_index(atoms: Atoms, index: int) -> Atom:
    """
    指定したインデックスの原子を取得する。

    Args:
        atoms (ase.Atoms): 原子構造を保持するASEのAtomsオブジェクト。
        index (int): 探索対象の原子インデックス。

    Returns:
        ase.Atom: 指定したインデックスの原子。
    """
    return atoms[index]


# 指定した原子のインデックスを調べる
def find_indices_by_symbol(atoms: Atoms, symbol: str):
    """
    指定した元素シンボルに一致する原子インデックスを取得する。

    Args:
        atoms (ase.Atoms): 原子構造を保持するASEのAtomsオブジェクト。
        symbol (str): 探索対象の元素シンボル（例: 'O', 'Fe'）。

    Returns:
        list[int]: 該当する原子インデックスのリスト。該当なしの場合は空リスト。
    """
    # 大文字小文字の揺れを吸収（ASEは大文字始まりが基本）
    target_symbol = symbol.capitalize()

    # 条件一致するインデックスを抽出
    indices = [i for i, atom in enumerate(atoms) if atom.symbol == target_symbol]

    return indices


# 隣接原子を探す
def get_neighbors(atoms: Atoms, target_atom: int | Atom, return_type: str = "atoms"):
    """
    指定した原子に対する隣接原子を取得する。

    Args:
        atoms (ase.Atoms): 原子構造を保持するASEのAtomsオブジェクト。
        target_atom (int | ase.Atom): 対象原子。インデックスまたはAtomオブジェクトを指定可能。
        return_type (str, optional): 返却形式を指定する。
            "atoms"   -> list[ase.Atom] 形式で返す（デフォルト）
            "indices" -> 隣接原子のインデックス（list[int]）で返す

    Returns:
        list[ase.Atom] または list[int]:
            隣接原子のリスト。`return_type` に応じて形式が変わる。

    Raises:
        ValueError: return_type が "atoms" または "indices" 以外の場合。
        IndexError: 指定インデックスが範囲外の場合。
        TypeError: target_atom の型が int でも ase.Atom でもない場合。
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

    # --- インデックス範囲チェック ---
    if index < 0 or index >= len(atoms):
        raise IndexError(f"インデックス {index} は範囲外です。")

    # --- return_type チェック ---
    if return_type not in ("atoms", "indices"):
        raise ValueError("return_type は 'atoms' または 'indices' を指定してください。")

    # --- 隣接リスト構築 ---
    cutoffs = natural_cutoffs(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    # --- 隣接原子取得 ---
    neighbor_indices, _ = nl.get_neighbors(index)

    # --- 出力形式に応じて返す ---
    if return_type == "atoms":
        return [atoms[i] for i in neighbor_indices]
    else:
        return neighbor_indices.tolist()


# (平面用)層別に分ける・層ごとのlistにする
def separate_layers(
    atoms: Atoms,
    return_type: Literal["atoms", "indices"] = "atoms",
    *,
    decimals: int = 4,
    sort_by_z: bool = True,
    use_substrate_mask: Literal["auto", True, False] = "auto",
) -> list[list[Atom]] | list[list[int]]:
    """
    (平面用)スラブ構造を層別に分離し、各層の原子をリストとして返す。

    z座標に基づいて原子を層ごとに分類し、指定された形式で返します。
    層は z 座標の昇順（bottom -> top）または降順でソートできます。

    Args:
        atoms (ase.Atoms): 分離対象のスラブ構造。
        return_type (Literal["atoms", "indices"], optional): 返却形式。
            "atoms": 各層の原子オブジェクトのリスト。デフォルト。
            "indices": 各層の原子インデックスのリスト。
        decimals (int, optional): z座標の丸め精度（小数点以下の桁数）。
            デフォルトは4。層の判定精度に影響します。
        sort_by_z (bool, optional): z座標で層をソートするか。
            True: z座標昇順（下層から上層）、False: 検出順。デフォルトはTrue。
        use_substrate_mask (Literal["auto", True, False], optional): 基板マスクの使用設定。
            "auto": atoms.arrays に "is_substrate" が存在する場合、is_substrate==True の原子のみで層を検出。
            True: 基板マスクを使用（存在しない場合は全原子）。
            False: マスクを無視して全原子を対象。
            デフォルトは "auto"。

    Returns:
        list[list[ase.Atom]] | list[list[int]]:
            各層の原子リスト。return_type に応じて形式が変わります。
            - "atoms": [[layer0_atoms], [layer1_atoms], ...]
            - "indices": [[layer0_indices], [layer1_indices], ...]

            sort_by_z=True の場合、layered_atoms[0] が最下層、
            layered_atoms[-1] が最上層になります。

    Raises:
        ValueError: return_type が "atoms" または "indices" 以外の場合。

    Examples:
        >>> from ase.build import surface, bulk
        >>> slab = surface(bulk('Cu'), (1,1,1), layers=3)
        >>> layers = separate_layers(slab, return_type="indices")
        >>> print(f"層数: {len(layers)}, 最下層原子数: {len(layers[0])}")

        >>> # 最表面の原子を取得
        >>> top_layer = separate_layers(slab)[-1]
        >>> print(f"最表面原子数: {len(top_layer)}")
    """
    # --- return_type の検証 ---
    if return_type not in ("atoms", "indices"):
        raise ValueError("return_type は 'atoms' または 'indices' を指定してください。")

    # --- 基板マスクの適用判定 ---
    should_use_mask = False
    if use_substrate_mask == "auto":
        should_use_mask = "is_substrate" in atoms.arrays
    elif use_substrate_mask is True:
        should_use_mask = "is_substrate" in atoms.arrays
    # use_substrate_mask == False の場合は should_use_mask = False のまま

    # --- 層検出対象の原子インデックスを決定 ---
    if should_use_mask:
        substrate_mask = atoms.arrays["is_substrate"].astype(bool)
        target_indices = np.where(substrate_mask)[0]
    else:
        target_indices = np.arange(len(atoms))

    # --- z座標を丸めて一意な層を特定 ---
    z_coords = atoms.positions[target_indices, 2]
    rounded_z = np.round(z_coords, decimals=decimals)
    unique_z_values = np.unique(rounded_z)

    # --- 層をz座標でソート（昇順：下層から上層） ---
    if sort_by_z:
        unique_z_values.sort()

    # --- 各層に属する原子インデックスを収集 ---
    layers_indices: list[list[int]] = []

    for z_value in unique_z_values:
        # 該当するz座標を持つ原子のインデックスを取得
        layer_mask = np.isclose(rounded_z, z_value, atol=10 ** (-decimals - 1))
        # target_indices の中でのマスク位置を、元の atoms でのインデックスに変換
        layer_indices = target_indices[layer_mask].tolist()
        layers_indices.append(layer_indices)

    # --- 出力形式に応じて返却 ---
    if return_type == "atoms":
        layers_atoms: list[list[Atom]] = []
        for layer_indices in layers_indices:
            layer_atoms = [atoms[i] for i in layer_indices]
            layers_atoms.append(layer_atoms)
        return layers_atoms
    else:  # return_type == "indices"
        return layers_indices


# (クラスター用)表面・内側を探す
def classify_surface_atoms(
    atoms: Atoms,
    return_type: Literal["atoms", "indices"] = "atoms",
    upper_tolerance: int = 3,
) -> tuple[list[Atoms], list[Atoms]] | tuple[list[int], list[int]]:
    """
    (クラスター用)クラスターの原子を外表面原子と内側原子に分類する。

    配位数（隣接原子数）が最小配位数からの範囲内にある原子を外表面原子とみなす。

    Args:
        atoms (ase.Atoms): 対象の原子構造。
        return_type (Literal["atoms", "indices"], optional):
            出力形式。
            "atoms" なら原子オブジェクトリスト、
            "indices"なら原子インデックスリストを返す。
            デフォルトは "atoms"。
        upper_tolerance (int, optional):
            最小配位数から上方向への許容範囲。デフォルトは3。
            最小配位数 + upper_tolerance が表面判定の上限となる。
            配位数が多いもの(=頂点性が少し低いもの)も表面原子として判定するようになる

    Returns:
        tuple[list[Atoms], list[Atoms]] または tuple[list[int], list[int]]:
            (表面原子リスト, 内側原子リスト)
    """
    # カットオフ距離を自動設定し、隣接リストを作成
    cutoffs = natural_cutoffs(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    coordination_numbers = []

    # 各原子の配位数（隣接原子数）を計算
    for i in range(len(atoms)):
        neighbors, _ = nl.get_neighbors(i)
        coordination_numbers.append(len(neighbors))

    # 最小の配位数を取得
    min_coordination = min(coordination_numbers)

    # 表面原子判定の範囲を計算
    surface_min = min_coordination
    surface_max = min_coordination + upper_tolerance

    # 配位数が範囲内の原子は表面、それ以外は内側と分類
    surface_indices = [
        i
        for i, cnum in enumerate(coordination_numbers)
        if surface_min <= cnum <= surface_max
    ]
    inner_indices = [
        i for i, cnum in enumerate(coordination_numbers) if cnum > surface_max
    ]

    # --- 出力形式に応じて返す ---
    if return_type == "atoms":
        surface_atoms = [atoms[i] for i in surface_indices]
        inner_atoms = [atoms[i] for i in inner_indices]
        return surface_atoms, inner_atoms
    else:
        return surface_indices, inner_indices


# 重心に最も近い原子を探す
def find_central_atom(
    atoms: Atoms | list[Atom], return_type: Literal["atom", "index"] = "atom"
) -> Atom | int:
    """
    xy面の重心に最も近い原子を返す。

    指定された原子構造またはリストの中で、xy面での重心座標に最も近い位置にある
    原子を見つけて返します。

    Args:
        atoms (ase.Atoms | list[ase.Atom]): 対象の原子構造またはリスト。
        return_type (Literal["atom", "index"], optional): 返却形式。
            "atom": 原子オブジェクト（デフォルト）。
            "index": 原子インデックス。list[Atom]の場合はリスト内のインデックス。

    Returns:
        ase.Atom | int: xy重心に最も近い原子。return_type に応じて形式が変わります。

    Raises:
        ValueError: return_type が不正、または原子が存在しない場合。
        TypeError: atoms の型が不正な場合。

    Examples:
        >>> from ase.build import surface, bulk
        >>> slab = surface(bulk('Cu'), (1,1,1), layers=3)

        >>> # 全体での中央原子
        >>> central_atom = find_central_atom(slab)
        >>> print(f"中央原子: {central_atom}")

        >>> # 特定層での中央原子（separate_layersと組み合わせ）
        >>> layers = separate_layers(slab, return_type="atoms")
        >>> surface_central = find_central_atom(layers[-1])  # 最表面
        >>> print(f"表面中央原子: {surface_central}")
    """
    # --- return_type の検証 ---
    if return_type not in ("atom", "index"):
        raise ValueError("return_type は 'atom' または 'index' を指定してください。")

    # --- 入力型に応じて処理を分岐 ---
    if isinstance(atoms, Atoms):
        # Atoms オブジェクトの場合
        target_atoms = [atoms[i] for i in range(len(atoms))]
        target_positions = atoms.positions
    elif isinstance(atoms, list) and all(isinstance(atom, Atom) for atom in atoms):
        # list[Atom] の場合
        if not atoms:
            raise ValueError("原子リストが空です。")
        target_atoms = atoms
        target_positions = np.array([atom.position for atom in atoms])
    else:
        raise TypeError("atoms は ase.Atoms または list[ase.Atom] を指定してください。")

    if len(target_atoms) == 0:
        raise ValueError("対象原子が存在しません。")

    # --- xy面での重心を計算 ---
    xy_positions = target_positions[:, :2]  # x, y座標のみ
    centroid_xy = np.mean(xy_positions, axis=0)

    # --- 重心に最も近い原子を探索 ---
    distances_squared = np.sum((xy_positions - centroid_xy) ** 2, axis=1)
    closest_idx = np.argmin(distances_squared)

    # --- 出力形式に応じて返却 ---
    if return_type == "atom":
        return target_atoms[closest_idx]
    else:  # return_type == "index"
        return closest_idx


# くっつけた後の構造の中で、くっつけた原子のインデックスを知る
def get_appended_atom_indices(before_atoms: Atoms, after_atoms: Atoms) -> list[int]:
    """くっつけた後の構造から、追加したAtoms（原子）のインデックスを返す。

    Args:
        before_atoms (Atoms): 追加前のAtomsオブジェクト。
        after_atoms (Atoms): 追加後のAtomsオブジェクト。

    Returns:
        list[int]: 追加されたAtoms（原子）がafter_atoms中で占めるインデックスリスト。

    Raises:
        ValueError: before_atoms の原子数が after_atoms より多い場合。

    Note:
        原子数の多い方から少ない方を引いてインデックスを求めるシンプルな動作です。
    """
    # --- 入力チェック ---
    if len(before_atoms) > len(after_atoms):
        raise ValueError("after_atoms の方が before_atoms より原子数が少ないです。")

    # --- インデックス計算処理 ---
    n_before = len(before_atoms)
    n_after = len(after_atoms)

    # --- 追加された原子のインデックスリストを作成 ---
    appended_indices = list(range(n_before, n_after))

    return appended_indices


# 自身の隣接原子数と同じ数だけ、隣接原子を持つ隣接原子を探す
def get_neighbors_with_coordination_condition(
    atoms: Atoms,
    target_atom: int | Atom,
    return_type: str = "atoms",
    *,
    cutoffs: Optional[Sequence[float]] = None,
    cutoff_scaling: float = 1.0,
    upper_tolerance: int = 1,
    lower_tolerance: int = 1,
) -> list[Atom] | list[int]:
    """
    get_neighbors() で得た隣接原子のうち、各隣接原子の配位数がtarget_atom周辺の配位数に収まるものを返す。

    Args:
        atoms (ase.Atoms): 構造。
        target_atom (int | ase.Atom): 対象原子（インデックスまたは Atom）。
        return_type (str): "atoms"（Atom のリスト）または "indices"（インデックスのリスト）。
        cutoffs (Sequence[float] | None): 配位数計算に用いるカットオフ配列。
            None の場合は coordination_number() 内で natural_cutoffs(atoms) を用いる。
        cutoff_scaling (float): natural_cutoffs 使用時のスケール係数。
        upper_tolerance (int): 最小配位数から上方向への許容範囲。デフォルトは1。
        lower_tolerance (int): 最小配位数から下方向への許容範囲。デフォルトは1。
            これらにより、target_atomの配位数-lower_tolerance ~ 配位数+upper_tolerance の範囲内の原子を隣接原子として採用する。

    Returns:
        list[ase.Atom] | list[int]: フィルタ済みの隣接原子（形式は return_type に依存）。

    Raises:
        ValueError: return_type が不正な場合。
        IndexError, TypeError: 下位関数に準じる。
    """
    if return_type not in ("atoms", "indices"):
        raise ValueError("return_type は 'atoms' または 'indices' を指定してください。")

    # 対象原子の配位数
    target_cn, _ = coordination_number(
        atoms,
        target_atom,
        return_type="indices",
        cutoffs=cutoffs,
        cutoff_scaling=cutoff_scaling,
    )

    # まずは通常の隣接原子（インデックス）を取得
    neighbor_indices = get_neighbors(atoms, target_atom, return_type="indices")

    # 各隣接原子の配位数を計算し、target と同じもののみ残す
    filtered_indices: list[int] = []
    for nbr_idx in neighbor_indices:
        nbr_cn, _ = coordination_number(
            atoms,
            nbr_idx,
            return_type="indices",
            cutoffs=cutoffs,
            cutoff_scaling=cutoff_scaling,
        )
        # print(f"{target_cn - lower_tolerance} ~ {target_cn + upper_tolerance} | {nbr_cn}")
        if target_cn - lower_tolerance <= nbr_cn <= target_cn + upper_tolerance:
            filtered_indices.append(nbr_idx)

    if return_type == "atoms":
        return [atoms[i] for i in filtered_indices]
    return filtered_indices
