"""
# CalcValue.py
ASEのAtomsオブジェクトから、特定の値を計算するための関数をまとめたファイル

## 関数一覧
- coordination_number(): 配位数を計算する
- compute_surface_normal(): (クラスター用)法線ベクトルを計算する
"""

from typing import Sequence, Optional, Literal
from ase import Atoms, Atom
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
from numpy.typing import NDArray
from FindAtoms import get_neighbors


# 配位数を計算する
def coordination_number(
    atoms: Atoms,
    target_atom: int | Atom,
    return_type: Literal["atoms", "indices"] = "atoms",
    *,
    cutoffs: Optional[Sequence[float]] = None,
    cutoff_scaling: float = 1.0,
) -> tuple[int, list[Atom] | list[int]]:
    """
    指定した原子の配位数（coordination number）を計算し、同時に隣接原子リストを返す。

    ASE の隣接判定（NeighborList / natural_cutoffs）を用いて隣接原子を決定し、
    配位数はその隣接原子数として返します。

    Args:
        atoms (ase.Atoms): 対象の構造を保持する ASE Atoms オブジェクト。
        target_atom (int | ase.Atom): 配位数を調べたい原子。インデックス（int）または
            `ase.Atom` オブジェクトのいずれかを指定可能。
        return_type (str, optional): 隣接原子リストの返却形式。
            - "atoms": list[ase.Atom] を返す（デフォルト）
            - "indices": list[int] を返す
        cutoffs (Sequence[float] | None, optional): 原子ごとのカットオフ半径の配列を
            直接与える場合に使用。None の場合は `natural_cutoffs(atoms)` を使う。
            （長さは `len(atoms)` と一致する必要があります）
        cutoff_scaling (float, optional): `natural_cutoffs` を使う場合のスケーリング係数。
            デフォルトは 1.0（例えば 1.2 にすると少し広めに隣接を拾います）。

    Returns:
        Tuple[int, Union[List[ase.Atom], List[int]]]:
            (coordination_number, neighbors)
            - coordination_number (int): 隣接原子数
            - neighbors: return_type に応じて list[Atom] か list[int] を返す

    Raises:
        TypeError: target_atom が int でも ase.Atom でもない場合。
        IndexError: 指定インデックスが範囲外の場合。
        ValueError: return_type が "atoms" でも "indices" でもない場合、
                    または cutoffs を与えたときに長さが不適切な場合。
    """
    # --- target_atom をインデックスに変換 ---
    if isinstance(target_atom, int):
        idx = target_atom
    elif isinstance(target_atom, Atom):
        try:
            idx = target_atom.index
        except ValueError:
            raise ValueError(
                "指定された Atom オブジェクトは `atoms` 内に存在しません。"
            )
    else:
        raise TypeError(
            "`target_atom` は int（インデックス）または ase.Atom を指定してください。"
        )

    # --- インデックス範囲チェック ---
    if idx < 0 or idx >= len(atoms):
        raise IndexError(f"インデックス {idx} は範囲外です（0 .. {len(atoms)-1}）。")

    # --- return_type チェック ---
    if return_type not in ("atoms", "indices"):
        raise ValueError(
            "`return_type` は 'atoms' または 'indices' を指定してください。"
        )

    # --- カットオフ半径の用意 ---
    if cutoffs is None:
        # natural_cutoffs を用いて原子種に基づくデフォルトのカットオフを得る
        base_cutoffs = natural_cutoffs(atoms)
        # スケーリングを反映
        cutoffs_used = [c * cutoff_scaling for c in base_cutoffs]
    else:
        # ユーザー指定の cutoffs を使用（長さチェック）
        if len(cutoffs) != len(atoms):
            raise ValueError(
                "`cutoffs` の長さは atoms の長さと一致する必要があります。"
            )
        cutoffs_used = list(cutoffs)

    # --- 隣接リストの構築 ---
    # self_interaction=False: 自分自身を隣接に含めない
    nl = NeighborList(cutoffs_used, self_interaction=False, bothways=True)
    nl.update(atoms)

    # --- 指定原子の隣接情報取得 ---
    neighbor_result = nl.get_neighbors(idx)
    if neighbor_result is None:
        # 隣接がない場合は空リスト
        neighbor_indices = []
    else:
        neighbor_indices, _ = neighbor_result

    # --- 配位数 (int) と希望の形式で隣接リストを作成 ---
    coord_num = len(neighbor_indices)
    if return_type == "atoms":
        neighbors = [atoms[i] for i in neighbor_indices]
    else:  # "indices"
        neighbors = list(neighbor_indices)

    return coord_num, neighbors


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
