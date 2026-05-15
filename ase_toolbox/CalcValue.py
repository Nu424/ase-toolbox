"""
# CalcValue.py
ASEのAtomsオブジェクトから、特定の値を計算するための関数をまとめたファイル

## 関数一覧
- coordination_number(): 配位数を計算する
"""

from itertools import combinations
import math
from typing import Iterator, Sequence, Optional, Literal
from ase import Atoms, Atom
from ase.neighborlist import NeighborList, natural_cutoffs


def _validate_step(
    step: float,
    *,
    max_total: float = 1.0,
    tol: float = 1e-12,
) -> tuple[int, float, float]:
    """
    刻み幅 step と最大合計値 max_total を検証し、分割数と正規化値を返す。
    """
    try:
        step_value = float(step)
    except (TypeError, ValueError):
        raise TypeError("step は float と互換性のある数値を指定してください。") from None

    try:
        max_total_value = float(max_total)
    except (TypeError, ValueError):
        raise TypeError("max_total は float と互換性のある数値を指定してください。") from None

    if not (0.0 < max_total_value <= 1.0):
        raise ValueError("max_total は 0 より大きく 1 以下の値を指定してください。")
    if not (0.0 < step_value <= max_total_value):
        raise ValueError("step は 0 より大きく max_total 以下の値を指定してください。")

    total_units = int(round(max_total_value / step_value))
    if total_units <= 0:
        raise ValueError("step から算出される分割数が不正です。別の刻み幅を指定してください。")

    if abs(max_total_value - total_units * step_value) > tol:
        raise ValueError(
            "step は max_total を分割できる値のみ指定できます。"
            "例: max_total=1.0 なら 0.5, 0.25, 0.1、"
            "max_total=0.15 なら 0.05 など。"
        )

    return total_units, step_value, max_total_value


def _normalize_elements(elements: Sequence[str]) -> list[str]:
    """
    元素リストを検証し、ユニークなリストに整形する。
    """
    if not elements:
        raise ValueError("elements には 1 つ以上の元素記号を指定してください。")

    normalized: list[str] = []
    seen: set[str] = set()

    for symbol in elements:
        if not isinstance(symbol, str):
            raise TypeError("elements には str のみ指定可能です。")
        if symbol in seen:
            raise ValueError(f"elements に重複が含まれています: {symbol}")
        normalized.append(symbol)
        seen.add(symbol)

    return normalized


def _validate_rounding_decimals(rounding_decimals: int | None) -> int | None:
    """
    丸め桁指定を検証し、正規化した値を返す。
    """
    if rounding_decimals is None:
        return None

    if not isinstance(rounding_decimals, int):
        raise TypeError("rounding_decimals は int または None を指定してください。")
    if rounding_decimals < 0:
        raise ValueError("rounding_decimals は 0 以上の整数を指定してください。")

    return rounding_decimals


def _iter_simplex_counts(num_elements: int, total_units: int) -> Iterator[list[int]]:
    """
    Bars and Stars を用いて、合計 total_units となる非負整数解を列挙する。

    Args:
        num_elements (int): 元素の数。1 以上の整数を指定。
        total_units (int): 合計値。0 以上の整数を指定。

    Yields:
        list[int]: 各元素の個数。["Cu", "Au", "Pd"] の場合は [0, 0, 1] など。

    Raises:
        ValueError: num_elements が 1 未満、または total_units が 0 未満の場合。
    """
    if num_elements <= 0:
        raise ValueError("num_elements は 1 以上を指定してください。")
    if total_units < 0:
        raise ValueError("total_units は 0 以上を指定してください。")

    if num_elements == 1:
        yield [total_units]
        return

    upper = total_units + num_elements - 1 # upper(=最大数)は、星の数(=total_units) + 棒の数(=num_elements - 1)
    for dividers in combinations(range(upper), num_elements - 1): # 最大の枠のうち、どこに棒を入れるか
        counts: list[int] = []
        prev = -1 # 前の棒の位置
        for divider in (*dividers, upper):
            counts.append(divider - prev - 1) # 前の棒の位置から、新しい棒の位置まで、星の数を計算
            prev = divider
        yield counts


def _counts_to_composition(
    elements: Sequence[str],
    counts: Sequence[int],
    total_units: int,
    step: float,
    max_total: float,
    rounding_decimals: int | None,
) -> dict[str, float]:
    """
    カウント列を組成辞書に変換する。

    Args:
        elements (Sequence[str]): 元素のリスト。
        counts (Sequence[int]): 各元素の個数。_iter_simplex_counts() で生成されたリスト。
        total_units (int): 合計値。
        step (float): 刻み幅。
        max_total (float): 組成比の合計値。
        rounding_decimals (int | None): 戻り値を丸める小数点以下桁数。

    Returns:
        dict[str, float]: 元素記号をキー、組成比を値とする辞書。

    Raises:
        ValueError: total_units が 0 未満の場合。
    """
    if total_units <= 0:
        raise ValueError("total_units は正の整数である必要があります。")

    fractions = [count * step for count in counts] # countは星の数。unitsの値まで取りうるので、step * 星数で、max_totalになりうる

    if rounding_decimals is None:
        return {sym: frac for sym, frac in zip(elements, fractions)} # sym: symbol(元素記号)

    rounded = [round(frac, rounding_decimals) for frac in fractions]
    if len(rounded) >= 2:
        partial_sum = sum(rounded[:-1])
        rounded[-1] = round(max_total - partial_sum, rounding_decimals)
    else:
        rounded[0] = round(max_total, rounding_decimals)

    return {sym: val for sym, val in zip(elements, rounded)}


def iter_simplex_compositions(
    elements: Sequence[str],
    step: float,
    *,
    max_total: float = 1.0,
    rounding_decimals: int | None = None,
) -> Iterator[dict[str, float]]:
    """
    シンプレックス格子点上の組成（合計 max_total）を逐次的に生成する。

    Args:
        elements (Sequence[str]): 構成元素のリスト。重複は不可。
        step (float): 刻み幅。max_total/step が整数（許容誤差内）となる値を指定する。
        max_total (float, optional): 組成比の合計値。0 より大きく 1 以下。
            デフォルトは 1.0。
        rounding_decimals (int | None, optional): 戻り値を丸める小数点以下桁数。
            None の場合は丸めずに返す。デフォルトは None。

    Yields:
        dict[str, float]: 元素記号をキー、組成比を値とする辞書。

    Raises:
        TypeError: elements に str 以外が含まれる場合、または rounding_decimals の型が不正な場合。
        ValueError: elements が空、重複を含む、step や max_total が条件を満たさない場合。

    Examples:
        >>> gen = iter_simplex_compositions(["Cu", "Au"], 0.5)
        >>> next(gen)
        {'Cu': 0.0, 'Au': 1.0}
        >>> next(gen)
        {'Cu': 0.5, 'Au': 0.5}
        >>> enumerate_simplex_compositions(["Cu", "Au"], 0.05, max_total=0.15, rounding_decimals=2)
        [{'Cu': 0.0, 'Au': 0.15}, {'Cu': 0.05, 'Au': 0.1}, {'Cu': 0.1, 'Au': 0.05}, {'Cu': 0.15, 'Au': 0.0}]

    Note:
        - 戻り値の辞書は `HandleAtoms.substitute_elements()` の new 引数に直接渡せる。
        - 大規模探索ではこの逐次版を使い、必要な組成のみ評価するのが効率的。
    """
    normalized = _normalize_elements(elements)
    rounding = _validate_rounding_decimals(rounding_decimals)
    total_units, step_value, max_total_value = _validate_step(
        step, max_total=max_total
    )

    for counts in _iter_simplex_counts(len(normalized), total_units):
        yield _counts_to_composition(
            normalized, counts, total_units, step_value, max_total_value, rounding
        )


def enumerate_simplex_compositions(
    elements: Sequence[str],
    step: float,
    *,
    max_total: float = 1.0,
    rounding_decimals: int | None = None,
) -> list[dict[str, float]]:
    """
    シンプレックス格子点上の全組成をリストとして取得する。

    Args:
        elements (Sequence[str]): 構成元素のリスト。順序は出力にも反映される。
        step (float): 刻み幅。max_total/step が整数（許容誤差内）となる値を指定する。
        max_total (float, optional): 組成比の合計値。0 より大きく 1 以下。
            デフォルトは 1.0。
        rounding_decimals (int | None, optional): 組成比を丸める桁数。None なら丸めない。

    Returns:
        list[dict[str, float]]: 各元素の組成辞書を要素とするリスト。

    Raises:
        TypeError: elements に str 以外が含まれる場合、または rounding_decimals の型が不正な場合。
        ValueError: elements が空、重複を含む、step や max_total が条件を満たさない場合。

    Examples:
        >>> enumerate_simplex_compositions(["Cu", "Au"], 0.5)
        [{'Cu': 0.0, 'Au': 1.0}, {'Cu': 0.5, 'Au': 0.5}, {'Cu': 1.0, 'Au': 0.0}]
    """
    return list(
        iter_simplex_compositions(
            elements,
            step,
            max_total=max_total,
            rounding_decimals=rounding_decimals,
        )
    )


def count_simplex_compositions(
    num_elements: int,
    step: float,
    *,
    max_total: float = 1.0,
) -> int:
    """
    刻み幅と元素数から、列挙される組成の総数を計算する。

    Args:
        num_elements (int): 元素の数。1 以上の整数を指定。
        step (float): 刻み幅。max_total/step が整数（許容誤差内）となる値を指定する。
        max_total (float, optional): 組成比の合計値。0 より大きく 1 以下。
            デフォルトは 1.0。

    Returns:
        int: 対応する組成の総数。

    Raises:
        TypeError: num_elements が int ではない場合。
        ValueError: num_elements が 1 未満、または step や max_total が条件を満たさない場合。

    Examples:
        >>> count_simplex_compositions(3, 0.25)
        15
    """
    if not isinstance(num_elements, int):
        raise TypeError("num_elements には int を指定してください。")
    if num_elements <= 0:
        raise ValueError("num_elements は 1 以上の整数を指定してください。")

    total_units, _, _ = _validate_step(step, max_total=max_total)
    return math.comb(total_units + num_elements - 1, num_elements - 1)



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
