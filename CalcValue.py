"""
# CalcValue.py
ASEのAtomsオブジェクトから、特定の値を計算するための関数をまとめたファイル

## 関数一覧
- coordination_number(): 配位数を計算する
"""

from typing import Sequence, Optional, Literal
from ase import Atoms, Atom
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
from numpy.typing import NDArray


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
