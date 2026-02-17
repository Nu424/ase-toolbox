# ----------
# ---格子定数の算出
# ----------
import re
from typing import Literal, Mapping
import warnings
from dataclasses import dataclass
from typing import Iterator

import numpy as np

from ase import Atoms
from ase.build import bulk
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer
from ase.filters import UnitCellFilter
from matlantis_features.ase_ext.optimize import FIRELBFGS


from ..util import normalize_composition
from ..HandleAtoms import substitute_elements


class LatticeConstantNotFoundError(KeyError):
    """格子定数が見つからない場合に投げる例外。"""


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
_NON_CUBIC_ELEMENTS: set[str] = {
    "Mg",
    "Ti",
    "Zn",
    "Zr",
    "Co",
    "Cd",
    "Be",
    "Ru",
    "Os",
    "Re",
}
_BCC_ELEMENTS: set[str] = {"Fe", "Cr", "W", "Mo", "V", "Nb", "Ta", "Ba"}


@dataclass
class LatticeConstant:
    """格子定数を格納するデータクラス

    Attributes:
        a (float): a軸の長さ [Å]
        b (float): b軸の長さ [Å]
        c (float): c軸の長さ [Å]
        alpha (float): α角 [度]
        beta (float): β角 [度]
        gamma (float): γ角 [度]
    """

    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float

    def __iter__(self) -> Iterator[float]:
        """イテレータプロトコルの実装

        Yields:
            float: 格子定数の各値 (a, b, c, alpha, beta, gamma の順)
        """
        yield self.a
        yield self.b
        yield self.c
        yield self.alpha
        yield self.beta
        yield self.gamma

    def __repr__(self) -> str:
        """文字列表現の実装

        Returns:
            str: 読みやすい格子定数の文字列表現
        """
        return (
            f"LatticeConstant(a={self.a:.4f}Å, b={self.b:.4f}Å, c={self.c:.4f}Å, "
            f"α={self.alpha:.2f}°, β={self.beta:.2f}°, γ={self.gamma:.2f}°)"
        )


# 元素の格子定数を取得
def _get_element_lattice_constant(
    symbol: str, user_map: dict[str, float] | None = None
) -> float:
    """
    元素の格子定数[a]（Å）を取得する。優先順は user_map > INTERNAL_LATTICE_MAP。

    Args:
        symbol: 元素記号（例: 'Cu'）。大文字小文字は自動正規化。
        user_map: ユーザ提供の格子定数マップ。キーは元素記号、値は格子定数[Å]。

    Returns:
        float: 格子定数[Å]。

    Raises:
        TypeError: symbol の型が不正な場合。
        ValueError: 元素記号が不正な形式の場合。
        LatticeConstantNotFoundError: 対象元素の格子定数が見つからない場合。
    """
    if not isinstance(symbol, str):
        raise TypeError("元素記号 symbol は str を指定してください。")

    sym = symbol.capitalize()
    if re.fullmatch(r"^[A-Z][a-z]?$", sym) is None:
        raise ValueError(f"元素記号の形式が不正です: {symbol}")

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
def calculate_lattice_constant_from_vegard(
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
        ValueError: 比率合計が1から外れる、負値が含まれる等。
        TypeError: 入力の型が不正な場合。
        LatticeConstantNotFoundError: 格子定数未定義の元素を含む場合。
        NotImplementedError: 未対応の method を指定した場合。
    """
    normalized = normalize_composition(
        composition,
        tol=tol,
        keep_zero=True,
        validate_symbols=True,
        renormalize=False,
    )
    symbols = list(normalized.keys())
    fractions = list(normalized.values())

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
        if any((s in _NON_CUBIC_ELEMENTS) for s in symbols) or any(
            (s in _BCC_ELEMENTS) for s in symbols
        ):
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


def calculate_lattice_constant_from_bulk(
    composition: Mapping[str, float],
    calculator: Calculator,
    *,
    is_pre_vegard: bool = False,
    lattice_map: dict[str, float] | None = None,
    seed: int | None = None,
    bulk_size: tuple[int, int, int] = (4, 4, 4),
    opt_fmax: float = 0.005,
    opt_maxsteps: int | None = None,
    optimizer_cls: type[Optimizer] = FIRELBFGS,
    crystal_structure: Literal["fcc", "bcc"] = "fcc",
) -> tuple[float, Atoms, LatticeConstant]:
    """バルクから格子定数を取得する

    Args:
        composition (Mapping[str, float]): 組成。
        calculator (Calculator): エネルギー計算に使用する計算機。
        is_pre_vegard (bool): 仮のベガード格子定数を使用するかどうか。
        lattice_map (dict[str, float] | None): 格子定数マップ。
        seed (int | None): 元素置換に用いる乱数シード。
        bulk_size (tuple[int, int, int]): バルクのサイズ。
        opt_fmax (float): 最適化のfmax。
        opt_maxsteps (int | None): 最適化の最大ステップ数。
        optimizer_cls (type[Optimizer]): 最適化に使用するクラス。
        crystal_structure (Literal["fcc", "bcc"]): 結晶構造(hcpは未対応)。

    Returns:
        tuple[float, Atoms, LatticeConstant]: 計算された格子定数、バルク構造、バルクから求めたままの格子定数オブジェクト。
    """
    if crystal_structure not in ("fcc", "bcc"):
        raise ValueError(f"Invalid crystal structure: {crystal_structure}")

    major_element = max(composition.items(), key=lambda item: item[1])[0]
    # ---バルクを作る
    bulk_atoms: Atoms = Atoms()
    if is_pre_vegard:
        # ---仮のベガード格子定数でバルクを作る
        vegard_a = calculate_lattice_constant_from_vegard(composition, lattice_map=lattice_map)
        bulk_atoms = bulk(
            major_element,
            a=vegard_a,
            crystalstructure=crystal_structure,
            cubic=True,
        )
    else:
        # ---ASE標準のバルクを作る
        bulk_atoms = bulk(
            major_element,
            crystalstructure=crystal_structure,
            cubic=True,
        )
    bulk_atoms = bulk_atoms.repeat(bulk_size)
    bulk_atoms.pbc = True
    # ---元素置換をおこなう
    bulk_atoms = substitute_elements(
        bulk_atoms,
        target=bulk_atoms,
        new=composition,
        inplace=False,
        seed=seed,
    )
    # ---バルクを最適化し、格子定数を取得する
    bulk_atoms.calc = calculator
    unit_cell_filter = UnitCellFilter(bulk_atoms)
    opt_dyn = optimizer_cls(unit_cell_filter)
    if opt_maxsteps is not None:
        opt_dyn.run(fmax=opt_fmax, steps=opt_maxsteps)
    else:
        opt_dyn.run(fmax=opt_fmax)
    # ---格子定数を取得する
    lengths_and_angles = bulk_atoms.get_cell_lengths_and_angles()
    raw_lattice_constant = LatticeConstant(
        a=lengths_and_angles[0],
        b=lengths_and_angles[1],
        c=lengths_and_angles[2],
        alpha=lengths_and_angles[3],
        beta=lengths_and_angles[4],
        gamma=lengths_and_angles[5],
    )
    # ---格子定数を計算する
    calculated_a = float(
        np.mean(
            [
                float(raw_lattice_constant.a) / float(bulk_size[0]),
                float(raw_lattice_constant.b) / float(bulk_size[1]),
                float(raw_lattice_constant.c) / float(bulk_size[2]),
            ]
        )
    )
    return calculated_a, bulk_atoms, raw_lattice_constant
