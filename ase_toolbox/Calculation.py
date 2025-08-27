"""
# Calculation.py
よく計算するものをまとめたファイル

## 関数一覧
- calculate_adsorption_energy(): 吸着エネルギーを計算する
- calculate_formation_energy(): 生成エネルギーを計算する
    - 付随する関数
        - analyze_composition(): 元素組成のdictを作成する
        - generate_reference_structure(): 純元素参照構造を生成する
- run_neb(): NEB計算を実行する
    - 付随する関数
        - plot_energy_profile(): エネルギープロファイルをプロットする
        - compute_barriers(): バリアを計算する
- calculate_delta_g(): ギブス自由エネルギーの差を計算する
    - 付随する関数
        - calculate_gibbs_free_energy(): ギブス自由エネルギーを計算する
- optimize_lattice_constant(): 格子定数を最適化する
"""

from typing import Literal, Optional, Any, Iterator
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from ase import Atoms
from ase.build import bulk, molecule
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer
from ase.neb import NEB
from ase.build.rotate import minimize_rotation_and_translation
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
from ase.constraints import UnitCellFilter
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator
from matlantis_features.ase_ext.optimize import FIRELBFGS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from .util import ConditionalLogger, ensure_logger, optimize_and_get_energy


# ----------
# ---吸着エネルギーを計算する
# ----------
@dataclass
class CAEInput:
    # calculate_adsorption_energy()で使用する入力を、構造的に扱うためのクラス。
    structure: Atoms
    calc_mode: Literal["molecule", "solid"]


# 吸着エネルギーを計算する
def calculate_adsorption_energy(
    calculator_molecule: Calculator,
    calculator_solid: Calculator,
    adsorbed_structure_input: CAEInput,
    reactant_structures_input: list[CAEInput],
    *,
    optimizer_cls: type[Optimizer] = FIRELBFGS,
    opt_fmax: float = 0.05,
    opt_maxsteps: int = 3000,
    logger: Optional[ConditionalLogger] = None,
    enable_logging: bool = True,
    copy_atoms: bool = True,
) -> float:
    """
    吸着エネルギーを計算する。

    指定された吸着後構造と吸着前構造群から、構造最適化とエネルギー計算を行い、
    吸着エネルギー E_ads = E(吸着後) - Σ[E(吸着前)] を算出します。

    Args:
        calculator_molecule (ase.calculators.calculator.Calculator): 分子用計算機。
            一般的に EstimatorCalcMode.MOLECULE を使用。
        calculator_solid (ase.calculators.calculator.Calculator): 固体用計算機。
            一般的に EstimatorCalcMode.CRYSTAL_U0 を使用。
            なお、それぞれ計算機を指定できるが、実際の計算の際には、1種類の計算機を使用するのがよい。計算機ごとにバイアスが異なるため。
        adsorbed_structure_input (CAEInput): 吸着後の構造（例: Cu-CO複合体）。
        reactant_structures_input (list[CAEInput]): 吸着前の構造群。
            例: [Cu表面, CO分子] のリスト。各構造は独立に最適化される。
        optimizer_cls (type[Optimizer], optional): 構造最適化アルゴリズムのクラス。
            デフォルトは FIRELBFGS。
        opt_fmax (float, optional): 構造最適化の力の収束閾値[eV/Å]。デフォルトは 0.05。
        opt_maxsteps (int, optional): 構造最適化の最大ステップ数。デフォルトは 3000。
        logger (ConditionalLogger | None, optional): ログ出力制御。
            Noneの場合は新規作成。
        enable_logging (bool, optional): ログ出力の有効/無効。デフォルトは True。
        copy_atoms (bool, optional): 原子構造をコピーするかどうか。Falseの場合は、入力の原子構造をそのまま使用する。デフォルトは True。

    Returns:
        float: 吸着エネルギー[eV]。負の値は吸着が熱力学的に有利であることを示す。

    Raises:
        ValueError: reactant_structures_input が空の場合。
        TypeError: 引数の型が不正な場合。
    """

    # --- ログ設定 ---
    logger = ensure_logger("adsorption", enable_logging, logger)

    # --- 引数検証 ---
    n_reactants = len(reactant_structures_input)
    if n_reactants == 0:
        raise ValueError(
            "reactant_structures は少なくとも1つの構造を含む必要があります。"
        )

    # --- メイン計算開始ログ ---
    logger.info("=" * 80)
    logger.info("吸着エネルギー計算開始")
    logger.info(
        f"吸着後構造: {adsorbed_structure_input.structure.symbols} ({len(adsorbed_structure_input.structure)} 原子)"
    )
    logger.info(f"反応物構造数: {n_reactants}")
    for i, atoms in enumerate(reactant_structures_input):
        logger.info(
            f"  反応物{i+1}: {atoms.structure.symbols} ({len(atoms.structure)} 原子)"
        )
    logger.info("計算モード設定:")
    logger.info(f"  吸着後: {adsorbed_structure_input.calc_mode}")
    for i, reactant in enumerate(reactant_structures_input):
        logger.info(f"  反応物{i+1}: {reactant.calc_mode}")

    # --- 1. 吸着後構造のエネルギー計算 ---
    adsorbed_calc = (
        calculator_molecule
        if adsorbed_structure_input.calc_mode == "molecule"
        else calculator_solid
    )
    e_adsorbed = optimize_and_get_energy(
        atoms=adsorbed_structure_input.structure,
        calculator=adsorbed_calc,
        optimizer_cls=optimizer_cls,
        fmax=opt_fmax,
        maxsteps=opt_maxsteps,
        label="吸着後構造",
        logger=logger,
        copy_atoms=copy_atoms,
    )

    # --- 2. 各反応物構造のエネルギー計算 ---
    reactant_energies: list[float] = []
    for i, reactant_input in enumerate(reactant_structures_input):
        label = f"反応物{i+1}"
        reactant_calc = (
            calculator_molecule
            if reactant_input.calc_mode == "molecule"
            else calculator_solid
        )
        e_reactant = optimize_and_get_energy(
            atoms=reactant_input.structure,
            calculator=reactant_calc,
            optimizer_cls=optimizer_cls,
            fmax=opt_fmax,
            maxsteps=opt_maxsteps,
            label=label,
            logger=logger,
            copy_atoms=copy_atoms,
        )
        reactant_energies.append(e_reactant)

    # --- 3. 吸着エネルギー計算 ---
    e_reactants_total = sum(reactant_energies)
    e_adsorption = e_adsorbed - e_reactants_total

    # --- 結果ログ出力 ---
    logger.info("=" * 80)
    logger.info("吸着エネルギー計算結果")
    logger.info(f"吸着後構造エネルギー: {e_adsorbed:.6f} eV")
    logger.info("反応物エネルギー:")
    for i, e in enumerate(reactant_energies):
        logger.info(
            f"  反応物{i+1} ({reactant_structures_input[i].structure.symbols}): {e:.6f} eV"
        )
    logger.info(f"反応物合計エネルギー: {e_reactants_total:.6f} eV")
    logger.info(f"吸着エネルギー: {e_adsorption:.6f} eV")
    if e_adsorption < 0:
        logger.info("→ 吸着は熱力学的に有利")
    else:
        logger.info("→ 吸着は熱力学的に不利")
    logger.info("=" * 80)

    # サマリーをロガーに出力
    logger.info(f"吸着エネルギー計算完了: E_ads = {e_adsorption:.3f} eV")

    return e_adsorption


# ----------
# ---生成エネルギーを計算する
# ----------
def analyze_composition(atoms: Atoms) -> dict[str, int]:
    """
    元素組成のdictを作成する。

    指定された原子構造に含まれる各元素の個数をカウントし、辞書形式で返す。

    Args:
        atoms (Atoms): 解析対象の原子構造。

    Returns:
        dict[str, int]: 元素記号をキー、個数を値とする辞書。
            例: {'Cu': 3, 'Au': 1}

    Raises:
        ValueError: 原子構造が空の場合。

    Examples:
        >>> from ase.build import bulk
        >>> cu3au = bulk('Cu', 'fcc', a=3.6).repeat((2,2,1))  # 4原子
        >>> # Cu3Auの構造を手動作成後
        >>> composition = analyze_composition(cu3au)
        >>> print(composition)
        {'Cu': 3, 'Au': 1}
    """
    # --- 入力検証 ---
    if len(atoms) == 0:
        raise ValueError("原子構造が空です。")

    # --- 元素組成の解析 ---
    symbols = [atom.symbol for atom in atoms]
    composition = dict(Counter(symbols))

    return composition


def generate_reference_structure(
    element: str,
    *,
    crystal_structure: str = "auto",
    lattice_parameter: Optional[float] = None,
    logger: Optional[ConditionalLogger] = None,
    enable_logging: bool = True,
) -> Atoms:
    """
    指定元素の純元素参照構造を生成する。

    金属元素の標準的な結晶構造（fcc, bcc, hcp）を自動判別または手動指定し、
    単原子の参照構造を生成します。

    Args:
        element (str): 元素記号（例: 'Cu', 'Fe', 'Au'）。
        crystal_structure (str, optional): 結晶構造の指定。
            'auto': 自動判別（デフォルト）
            'fcc', 'bcc', 'hcp': 手動指定
        lattice_parameter (float | None, optional): 格子定数[Å]。
            Noneの場合は ASE の標準値を使用。
        logger (ConditionalLogger | None, optional): ロガー。
            Noneの場合はログ出力しない。

    Returns:
        Atoms: 純元素の参照構造。結晶構造に応じて1つまたは複数の原子を含む。

    Raises:
        ValueError: 対応していない元素または結晶構造の場合。

    Note:
        自動判別は一般的な金属の結晶構造に基づいています：
        - fcc: Cu, Au, Ag, Al, Ni, Pt, Pd など
        - bcc: Fe, Cr, W, Mo, V など
        - hcp: Zn, Mg, Ti, Zr など

        生成される構造の原子数は結晶構造とASEのbulk()実装に依存します。
        生成エネルギー計算では自動的に原子あたりのエネルギーに正規化されます。
    """
    # --- ログ設定 ---
    logger = ensure_logger("reference", enable_logging, logger)

    # --- 元素記号の正規化 ---
    element = element.capitalize()

    # --- 結晶構造の自動判別 ---
    if crystal_structure == "auto":
        # 一般的な金属の結晶構造データベース
        crystal_db = {
            # fcc構造
            "Cu": "fcc",
            "Au": "fcc",
            "Ag": "fcc",
            "Al": "fcc",
            "Ni": "fcc",
            "Pt": "fcc",
            "Pd": "fcc",
            "Rh": "fcc",
            "Ir": "fcc",
            "Pb": "fcc",
            "Ca": "fcc",
            "Sr": "fcc",
            # bcc構造
            "Fe": "bcc",
            "Cr": "bcc",
            "W": "bcc",
            "Mo": "bcc",
            "V": "bcc",
            "Nb": "bcc",
            "Ta": "bcc",
            "Ba": "bcc",
            # hcp構造
            "Zn": "hcp",
            "Mg": "hcp",
            "Ti": "hcp",
            "Zr": "hcp",
            "Co": "hcp",
            "Cd": "hcp",
            "Be": "hcp",
            "Ru": "hcp",
            "Os": "hcp",
            "Re": "hcp",
        }

        if element not in crystal_db:
            raise ValueError(
                f"元素 '{element}' の結晶構造を自動判別できません。"
                f"crystal_structure を 'fcc', 'bcc', 'hcp' から手動指定してください。"
            )

        determined_structure = crystal_db[element]
        logger.info(f"元素 {element} の結晶構造を自動判別: {determined_structure}")
    else:
        # 手動指定の場合
        valid_structures = ["fcc", "bcc", "hcp"]
        if crystal_structure not in valid_structures:
            raise ValueError(
                f"結晶構造 '{crystal_structure}' は対応していません。"
                f"'auto' または {valid_structures} から選択してください。"
            )
        determined_structure = crystal_structure
        logger.info(f"元素 {element} の結晶構造を手動指定: {determined_structure}")

    # --- 参照構造の生成 ---
    try:
        if lattice_parameter is not None:
            # 格子定数を手動指定
            ref_structure = bulk(element, determined_structure, a=lattice_parameter)
            logger.info(f"格子定数を手動指定: a={lattice_parameter:.4f} Å")
        else:
            # ASEの標準値を使用
            ref_structure = bulk(element, determined_structure)
            # 生成された格子定数をログ出力
            cell_params = ref_structure.get_cell_lengths_and_angles()
            logger.info(f"ASE標準格子定数を使用: a={cell_params[0]:.4f} Å")

        logger.info(
            f"参照構造生成完了: {element} ({determined_structure}, {len(ref_structure)} 原子)"
        )
        return ref_structure

    except Exception as e:
        raise ValueError(f"元素 '{element}' の参照構造生成に失敗しました: {e}")


def calculate_formation_energy(
    calculator: Calculator,
    compound_structure: Atoms,
    *,
    optimizer_cls: type[Optimizer] = FIRELBFGS,
    opt_fmax: float = 0.05,
    opt_maxsteps: int = 3000,
    reference_crystal_structures: Optional[dict[str, str]] = None,
    reference_lattice_parameters: Optional[dict[str, float]] = None,
    logger: Optional[ConditionalLogger] = None,
    enable_logging: bool = True,
    copy_atoms: bool = True,
) -> float:
    """
    金属の生成エネルギーを計算する。

    指定された化合物構造から元素組成を自動解析し、純元素参照構造との
    エネルギー差から生成エネルギー E_formation = E(化合物) - Σ[n_i × E(純元素_i)] を算出します。

    Args:
        calculator (Calculator): 固体用計算機。
            一般的に EstimatorCalcMode.CRYSTAL_U0 を使用。
        compound_structure (Atoms): 化合物の原子構造。
        optimizer_cls (type[Optimizer], optional): 構造最適化アルゴリズムのクラス。
            デフォルトは FIRELBFGS。
        opt_fmax (float, optional): 構造最適化の力の収束閾値[eV/Å]。デフォルトは 0.05。
        opt_maxsteps (int, optional): 構造最適化の最大ステップ数。デフォルトは 3000。
        reference_crystal_structures (dict[str, str] | None, optional):
            純元素の結晶構造を手動指定する辞書。キーは元素記号、値は 'fcc', 'bcc', 'hcp'。
            例: {'Cu': 'fcc', 'Fe': 'bcc'}
        reference_lattice_parameters (dict[str, float] | None, optional):
            純元素の格子定数を手動指定する辞書。キーは元素記号、値は格子定数[Å]。
            例: {'Cu': 3.615, 'Au': 4.078}
        logger (ConditionalLogger | None, optional): ログ出力制御。
            Noneの場合は新規作成。
        enable_logging (bool, optional): ログ出力の有効/無効。デフォルトは True。
        copy_atoms (bool, optional): 原子構造をコピーするかどうか。Falseの場合は、入力の原子構造をそのまま使用する。デフォルトは True。

    Returns:
        float: 生成エネルギー[eV]。負の値は化合物形成が熱力学的に有利であることを示す。

    Raises:
        ValueError: 化合物構造が空、または未対応の元素が含まれる場合。
        TypeError: 引数の型が不正な場合。

    Examples:
        >>> from ase.build import bulk
        >>> from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
        >>> from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

        >>> # 計算機の設定
        >>> calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.CRYSTAL_U0))

        >>> # Cu3Au合金の生成エネルギー計算
        >>> cu3au = bulk('Cu', 'fcc', a=3.6).repeat((2,2,1))  # 仮の構造
        >>> formation_energy = calculate_formation_energy(calc, cu3au)
        >>> print(f"Cu3Au生成エネルギー: {formation_energy:.3f} eV")
    """
    # --- ログ設定 ---
    logger = ensure_logger("formation", enable_logging, logger)

    # --- 引数検証 ---
    if len(compound_structure) == 0:
        raise ValueError("化合物構造が空です。")

    # --- メイン計算開始 ---
    logger.info("=" * 80)
    logger.info("金属生成エネルギー計算開始")
    logger.info(
        f"化合物構造: {compound_structure.symbols} ({len(compound_structure)} 原子)"
    )

    # --- 1. 化合物の元素組成解析 ---
    composition = analyze_composition(compound_structure)
    logger.info("元素組成解析結果:")
    for element, count in composition.items():
        logger.info(f"  {element}: {count} 原子")

    # --- 2. 化合物のエネルギー計算 ---
    e_compound = optimize_and_get_energy(
        atoms=compound_structure,
        calculator=calculator,
        optimizer_cls=optimizer_cls,
        fmax=opt_fmax,
        maxsteps=opt_maxsteps,
        label="化合物構造",
        logger=logger,
        copy_atoms=copy_atoms,
    )

    # --- 3. 各純元素のエネルギー計算 ---
    element_energies: dict[str, float] = {}

    for element in composition.keys():
        logger.info(f"\n純元素 {element} の参照構造準備")

        # 結晶構造の決定
        crystal_structure = "auto"
        if reference_crystal_structures and element in reference_crystal_structures:
            crystal_structure = reference_crystal_structures[element]
            logger.info(f"結晶構造を手動指定: {crystal_structure}")

        # 格子定数の決定
        lattice_parameter = None
        if reference_lattice_parameters and element in reference_lattice_parameters:
            lattice_parameter = reference_lattice_parameters[element]
            logger.info(f"格子定数を手動指定: {lattice_parameter} Å")

        # 純元素参照構造の生成
        try:
            ref_structure = generate_reference_structure(
                element,
                crystal_structure=crystal_structure,
                lattice_parameter=lattice_parameter,
                logger=logger,
                enable_logging=enable_logging,
            )

            # エネルギー計算実行
            e_element_total = optimize_and_get_energy(
                atoms=ref_structure,
                calculator=calculator,
                optimizer_cls=optimizer_cls,
                fmax=opt_fmax,
                maxsteps=opt_maxsteps,
                label=f"純元素 {element}",
                logger=logger,
                copy_atoms=copy_atoms,
            )

            # 原子あたりのエネルギーに正規化
            e_element_per_atom = e_element_total / len(ref_structure)
            logger.info(
                f"純元素 {element} の原子あたりエネルギー: {e_element_per_atom:.6f} eV/atom (構造中 {len(ref_structure)} 原子)"
            )

            element_energies[element] = e_element_per_atom

        except Exception as e:
            logger.error(f"純元素 {element} の計算に失敗: {e}")
            raise ValueError(
                f"純元素 {element} の参照エネルギー計算に失敗しました: {e}"
            )

    # --- 4. 生成エネルギー計算 ---
    # 純元素エネルギーの加重和を計算
    e_elements_total = sum(
        count * element_energies[element] for element, count in composition.items()
    )

    # 生成エネルギー = 化合物エネルギー - 純元素エネルギー合計
    e_formation = e_compound - e_elements_total

    # --- 結果ログ出力 ---
    logger.info("=" * 80)
    logger.info("生成エネルギー計算結果")
    logger.info(f"化合物エネルギー: {e_compound:.6f} eV")
    logger.info("純元素エネルギー (原子あたり):")
    for element, energy_per_atom in element_energies.items():
        count = composition[element]
        total = count * energy_per_atom
        logger.info(
            f"  {element}: {energy_per_atom:.6f} eV/atom × {count} = {total:.6f} eV"
        )
    logger.info(f"純元素合計エネルギー: {e_elements_total:.6f} eV")
    logger.info(f"生成エネルギー: {e_formation:.6f} eV")
    if e_formation < 0:
        logger.info("→ 化合物形成は熱力学的に有利")
    else:
        logger.info("→ 化合物形成は熱力学的に不利")
    logger.info("=" * 80)

    # ロガーサマリー出力
    logger.info(
        f"生成エネルギー計算完了: 組成 {composition}, E_formation = {e_formation:.3f} eV"
    )

    return e_formation


# ----------
# ---NEB計算を実行する
# ----------
# NEB計算を実行する
def run_neb(
    init_atoms: Atoms,
    final_atoms: Atoms,
    num_intermediate_images: int,
    optimizer_cls: type[Optimizer],
    estimator: Estimator,
    *,
    fmax: float = 0.05,
    steps: int = 500,
    trajectory_path: str | None = None,
    pre_align: bool = True,
    k: float = 0.1,
    climb: bool = True,
    parallel: bool = False,
    mic: bool | None = None,
    interpolate_kwargs: dict[str, Any] | None = None,
) -> tuple[list[Atoms], list[float]]:
    """NEB計算を実行し、構造とエネルギーのリストを返す。Matlantis環境用。

    初期構造と最終構造を指定して、指定された数の中間構造を生成し、
    NEB計算を実行します。計算後の全構造とそれぞれのエネルギーを返します。

    Args:
        init_atoms: 初期構造。
        final_atoms: 最終構造。
        num_intermediate_images: 中間構造の数（端点は含まない）。
        optimizer_cls: 使用するオプティマイザーのクラス（例: FIRE, BFGS）。
        estimator: 計算器作成に使用するEstimatorオブジェクト。
        fmax: 収束判定の力の閾値 [eV/Å]。デフォルトは0.05。
        steps: 最大最適化ステップ数。デフォルトは500。
        trajectory_path: 軌跡を保存するファイルパス。Noneの場合は保存しない。
        pre_align: 初期・最終構造のアライメント（回転・並進最適化）を行うか。デフォルトはTrue。
        k: NEBのばね定数。デフォルトは0.1。
        climb: Climbing Image NEBを使用するか。デフォルトはTrue。
        parallel: 並列計算を行うか。Trueの場合は各画像に個別の計算器を作成、
                 Falseの場合は計算器を共有。デフォルトはFalse。
        mic: 最短イメージ規約を使用するか。Noneの場合は自動判定。
        interpolate_kwargs: interpolateメソッドに渡す追加の引数。

    Returns:
        tuple[list[Atoms], list[float]]: 計算後の構造リストとエネルギーリスト [eV]。

    Raises:
        ValueError: 構造の原子数が一致しない、または中間構造数が負の場合。
    """
    # 入力検証
    if len(init_atoms) != len(final_atoms):
        raise ValueError(
            f"初期構造と最終構造の原子数が一致しません: "
            f"{len(init_atoms)} != {len(final_atoms)}"
        )

    if num_intermediate_images < 0:
        raise ValueError(
            f"中間構造数は0以上である必要があります: {num_intermediate_images}"
        )

    # 構造のアライメント（回転・並進最適化）
    if pre_align:
        minimize_rotation_and_translation(init_atoms, final_atoms)

    # 画像リストの生成（初期構造 + 中間構造 + 最終構造）
    images = [init_atoms.copy()]
    images += [init_atoms.copy() for _ in range(num_intermediate_images)]
    images += [final_atoms.copy()]

    # parallelフラグに基づいてallow_shared_calculatorを決定
    allow_shared_calculator = not parallel

    # 計算器の設定
    if allow_shared_calculator:
        # 計算器を共有する場合（parallel=False）
        calculator = ASECalculator(estimator)
        for image in images:
            image.calc = calculator
    else:
        # 各画像に個別の計算器を作成する場合（parallel=True）
        for image in images:
            calculator = ASECalculator(estimator)
            image.calc = calculator

    # NEB オブジェクトの構築
    neb = NEB(
        images,
        k=k,
        climb=climb,
        allow_shared_calculator=allow_shared_calculator,
        parallel=parallel,
    )

    # 構造の補間
    interpolate_args = interpolate_kwargs or {}
    if mic is not None:
        interpolate_args["mic"] = mic
    neb.interpolate(**interpolate_args)

    # 最適化の実行
    if trajectory_path is not None:
        optimizer = optimizer_cls(neb, trajectory=trajectory_path)
    else:
        optimizer = optimizer_cls(neb)

    optimizer.run(fmax=fmax, steps=steps)

    # エネルギーの計算
    energies = [image.get_potential_energy() for image in images]

    return images, energies


def plot_energy_profile(
    energies: Sequence[float],
    *,
    ax: Axes | None = None,
    xlabel: str = "replica",
    ylabel: str = "energy [eV]",
    title: str | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """エネルギープロファイルの折れ線グラフを描画する。

    NEB計算で得られたエネルギーのリストから、反応経路のエネルギープロファイルを
    折れ線グラフとして描画します。

    Args:
        energies: 各画像のエネルギーのリスト [eV]。
        ax: 描画に使用するAxesオブジェクト。Noneの場合は新しく作成。
        xlabel: x軸のラベル。デフォルトは"replica"。
        ylabel: y軸のラベル。デフォルトは"energy [eV]"。
        title: グラフのタイトル。Noneの場合はタイトルなし。
        show: グラフを表示するかどうか。デフォルトはTrue。

    Returns:
        tuple[Figure, Axes]: 描画に使用したFigureとAxesオブジェクト。

    Examples:
        >>> energies = [-10.5, -10.2, -9.8, -10.1, -10.6]
        >>> fig, ax = plot_energy_profile(energies, title="Reaction Pathway")
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    # エネルギープロファイルのプロット
    replica_indices = list(range(len(energies)))
    ax.plot(replica_indices, energies, "o-", linewidth=2, markersize=6)

    # 軸ラベルとタイトルの設定
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    # グリッドの追加
    ax.grid(True, alpha=0.3)

    # レイアウトの調整
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def compute_barriers(
    energies: Sequence[float],
) -> tuple[int, float, float]:
    """エネルギーリストから活性化エネルギーを計算する。

    NEB計算で得られたエネルギーリストから、遷移状態のインデックスと
    順方向・逆方向の活性化エネルギーを計算します。

    Args:
        energies: 各画像のエネルギーのリスト [eV]。

    Returns:
        tuple[int, float, float]:
            - 遷移状態のインデックス
            - 順方向活性化エネルギー [eV] (遷移状態 - 初期状態)
            - 逆方向活性化エネルギー [eV] (遷移状態 - 最終状態)

    Note:
        - 遷移状態は最もエネルギーが高い画像として定義されます。
        - 複数の極大値がある場合は、最もエネルギーが高い点を遷移状態とします。

    Examples:
        >>> energies = [-10.5, -10.2, -9.8, -10.1, -10.6]
        >>> ts_idx, e_forward, e_backward = compute_barriers(energies)
        >>> print(f"遷移状態: {ts_idx}, 順方向: {e_forward:.2f} eV, 逆方向: {e_backward:.2f} eV")
    """
    if len(energies) < 2:
        raise ValueError("エネルギーリストは少なくとも2つの要素が必要です")

    # NumPy配列に変換
    energy_array = np.array(energies)

    # 遷移状態（最大エネルギー）のインデックス
    ts_index = np.argmax(energy_array)

    # 活性化エネルギーの計算
    e_act_forward = energy_array[ts_index] - energy_array[0]
    e_act_backward = energy_array[ts_index] - energy_array[-1]

    return int(ts_index), float(e_act_forward), float(e_act_backward)


# ----------
# ---ギブス自由エネルギーを計算する
# ----------
@dataclass
class CGFEInput:
    """
    calculate_gibbs_free_energy()で使用する入力を、構造的に扱うためのクラス。
    構造、振動させる分子のインデックス、計算モードをまとめて管理する。

    Args:
        structure: 計算する構造
        calc_mode: 計算モード。IdealGasThermo: 気相分子向け、HarmonicThermo: 固体向け
        vibrate_indices: 振動させる分子のインデックス。Noneの場合は全ての分子を振動させる。
        geometry: 分子の幾何学的構造。'linear'（線形分子）または'nonlinear'（非線形分子）
        symmetry_number: 分子の対称数（回転対称性）
        spin_multiplicity: スピン多重度（2S+1、Sは総スピン量子数）
        do_opt: 最適化を行うかどうか。Trueの場合は最適化を行う。
    """

    structure: Atoms
    calc_mode: Literal["IdealGasThermo", "HarmonicThermo"]
    vibrate_indices: list[int] | None = None
    # IdealGasThermo用の引数
    geometry: Literal["linear", "nonlinear"] = "nonlinear"
    symmetry_number: int = 1
    spin_multiplicity: int = 1
    # opt用の引数
    do_opt: bool = True


# ギブス自由エネルギーを計算する
def calculate_gibbs_free_energy(
    calculator_molecule: Calculator,
    calculator_solid: Calculator,
    calc_input: CGFEInput,
    *,
    temperature: float = 298.15,
    pressure: float = 101325.0,
    optimizer_cls: type[Optimizer] = FIRELBFGS,
    opt_fmax: float = 0.05,
    opt_maxsteps: int = 3000,
    logger: Optional[ConditionalLogger] = None,
    enable_logging: bool = True,
    cleanup_vibrations: bool = True,
    copy_atoms: bool = True,
):
    """
    ギブス自由エネルギーを計算する。

    Args:
        calculator_molecule (ase.calculators.calculator.Calculator): 分子用計算機
        calculator_solid (ase.calculators.calculator.Calculator): 固体用計算機
        calc_input (CGFEInput): 計算する構造、振動させる分子のインデックス、計算モード
        temperature (float): 温度[K]
        pressure (float): 圧力[Pa]
        optimizer_cls (type[Optimizer]): 最適化エンジンのクラス
        opt_fmax (float): 最適化の閾値
        opt_maxsteps (int): 最適化の最大ステップ数
        logger (ConditionalLogger): ロガー。Noneの場合はログを出力しない
        cleanup_vibrations (bool): 振動計算後にファイルをクリーンアップするか
        copy_atoms (bool, optional): 原子構造をコピーするかどうか。Falseの場合は、入力の原子構造をそのまま使用する。デフォルトは True。

    Returns:
        float: ギブス自由エネルギー。Δではない。
    """
    # --- ログ設定 ---
    logger = ensure_logger("gibbs", enable_logging, logger)

    # ---構造に計算機を設定
    atoms = calc_input.structure
    # calc_mode に応じて適切な calculator を割り当て
    if calc_input.calc_mode == "IdealGasThermo":
        atoms.calc = calculator_molecule
    else:
        atoms.calc = calculator_solid

    # 構造情報をログ出力
    logger.info("=" * 60)
    logger.info(f"構造解析開始: {atoms.symbols}")
    logger.info(f"原子数: {len(atoms)}")
    logger.info(f"計算モード: {calc_input.calc_mode}")
    logger.info(f"振動対象インデックス: {calc_input.vibrate_indices}")
    logger.info(f"最適化実行: {calc_input.do_opt}")

    # 拘束条件の確認
    constraints = atoms.constraints
    if constraints:
        logger.info(f"拘束条件: {len(constraints)} 個")
        for i, constraint in enumerate(constraints):
            logger.info(f"  拘束{i}: {type(constraint).__name__}")
            if hasattr(constraint, "index"):
                logger.info(
                    f"    固定原子数: {len(constraint.index) if hasattr(constraint.index, '__len__') else 1}"
                )
    else:
        logger.info("拘束条件: なし")

    # ---　1. 最適化
    if calc_input.do_opt:
        e_opt = optimize_and_get_energy(
            atoms=atoms,
            calculator=atoms.calc,  # atomsに設定済みの計算機を使用
            optimizer_cls=optimizer_cls,
            fmax=opt_fmax,
            maxsteps=opt_maxsteps,
            label=f"構造最適化: {atoms.symbols}",
            logger=logger,
            copy_atoms=copy_atoms,
        )
    else:
        logger.info("構造最適化はスキップされました (do_opt=False)。")
        e_opt = atoms.get_potential_energy()
        logger.info(f"ポテンシャルエネルギー: {e_opt:.6f} eV")

    # ---　2. 振動計算
    vib_indices = calc_input.vibrate_indices
    if vib_indices is None:
        vib_indices = list(range(len(atoms)))
        logger.info("振動計算: 全原子を対象")
    elif len(vib_indices) == 0:
        logger.info("振動計算: スキップ (vibrate_indices = [])")
    else:
        logger.info(f"振動計算: 指定原子 {vib_indices} を対象")

    # 振動計算を実行するかどうか
    if len(vib_indices) == 0:
        vib_energies = np.array([])
    else:
        logger.info(f"振動対象原子数: {len(vib_indices)}")
        logger.info(f"振動モード数: {3 * len(vib_indices)} (理論値)")
        vib = Vibrations(atoms, indices=vib_indices)
        vib.clean()
        logger.info("振動計算開始")
        vib.run()
        logger.info("振動計算完了")
        vib_energies = vib.get_energies()

        # 振動計算後のクリーンアップ
        if cleanup_vibrations:
            logger.info("振動計算ファイルのクリーンアップ実行")
            vib.clean()
            logger.info("クリーンアップ完了")

    # 振動解析の詳細ログ
    logger.info(f"振動エネルギー数: {len(vib_energies)}")

    # 虚数モードのチェック
    imag_modes = []
    real_modes = []
    for i, energy in enumerate(vib_energies):
        if np.iscomplex(energy) or energy < 0:
            imag_modes.append((i, energy))
        else:
            real_modes.append((i, energy))

    logger.info(f"実数モード数: {len(real_modes)}")
    logger.info(f"虚数モード数: {len(imag_modes)}")

    if imag_modes:
        logger.warning("虚数モードが検出されました:")
        for i, energy in imag_modes:
            logger.warning(f"  モード{i}: {energy}")

    # 振動エネルギーの統計
    real_energies = [
        e.real if np.iscomplex(e) else e
        for e in vib_energies
        if (np.iscomplex(e) and e.real > 0) or (not np.iscomplex(e) and e > 0)
    ]
    if real_energies:
        logger.info(
            f"実振動エネルギー範囲: {min(real_energies):.6f} - {max(real_energies):.6f} eV"
        )
        logger.info(f"実振動エネルギー平均: {np.mean(real_energies):.6f} eV")

    # ---　3. 熱化学補正量計算
    if vib_energies.size == 0:
        logger.info("熱化学補正計算をスキップ (振動エネルギーなし)")
        g = e_opt
        thermal_correction = 0.0
    else:
        logger.info("熱化学補正計算開始")
        if calc_input.calc_mode == "IdealGasThermo":
            logger.info("IdealGasThermo計算:")
            logger.info(f"  幾何構造: {calc_input.geometry}")
            logger.info(f"  対称数: {calc_input.symmetry_number}")
            logger.info(f"  スピン多重度: {calc_input.spin_multiplicity}")

            thermo = IdealGasThermo(
                vib_energies=vib_energies,
                potentialenergy=e_opt,
                atoms=atoms,
                geometry=calc_input.geometry,
                symmetrynumber=calc_input.symmetry_number,
                spin=calc_input.spin_multiplicity,
            )
            zpe = thermo.get_ZPE_correction()
            enthalpy_correction = thermo.get_enthalpy(temperature)
            entropy = thermo.get_entropy(temperature, pressure)
            logger.info(f"  ZPE補正: {zpe:.6f} eV")
            logger.info(f"  エンタルピー補正: {enthalpy_correction:.6f} eV")
            logger.info(f"  エントロピー項 (-TS): {-temperature * entropy:.6f} eV")
            logger.info(f"  エントロピー: {entropy:.6f} eV/K")
            g = thermo.get_gibbs_energy(temperature=temperature, pressure=pressure)
        else:  # HarmonicThermo
            logger.info("HarmonicThermo計算:")
            thermo = HarmonicThermo(
                vib_energies=vib_energies,
                potentialenergy=e_opt,
                ignore_imag_modes=True,
            )
            zpe = thermo.get_ZPE_correction()
            internal_energy_correction = thermo.get_internal_energy(temperature)
            entropy = thermo.get_entropy(temperature)
            logger.info(f"  ZPE補正: {zpe:.6f} eV")
            logger.info(f"  内部エネルギー補正: {internal_energy_correction:.6f} eV")
            logger.info(f"  エントロピー項 (-TS): {-temperature * entropy:.6f} eV")
            logger.info(f"  エントロピー: {entropy:.6f} eV/K")
            g = thermo.get_helmholtz_energy(temperature=temperature)
        # 固体系では G ≈ F
        thermal_correction = g - e_opt

    # エネルギー成分の詳細分析
    logger.info(f"熱化学補正合計: {thermal_correction:.6f} eV")
    logger.info(f"最終ギブス自由エネルギー: {g:.6f} eV")
    logger.info(f"  内訳: E_pot({e_opt:.6f}) + 熱補正({thermal_correction:.6f})")
    logger.info("=" * 60)

    # サマリーをロガーに出力
    logger.info(f"ギブス自由エネルギー計算完了: G({atoms.symbols}) = {g:.3f} eV")
    return g


# ギブス自由エネルギーの差を計算する
def calculate_delta_g(
    calculator_molecule: Calculator,
    calculator_solid: Calculator,
    reactants: list[CGFEInput | Literal["CHE"]],
    products: list[CGFEInput | Literal["CHE"]],
    *,
    temperature: float = 298.15,
    pressure: float = 101325.0,
    electrode_potential: float = 0.0,
    pH: float = 7.0,
    optimizer_cls: type[Optimizer] = FIRELBFGS,
    opt_fmax: float = 0.05,
    opt_maxsteps: int = 3000,
    logger: Optional[ConditionalLogger] = None,
    enable_logging: bool = True,
    cleanup_vibrations: bool = True,
    copy_atoms: bool = True,
):
    """
    ギブス自由エネルギーを計算する。CHEモデルに対応している

    Args:
        calculator_molecule (ase.calculators.calculator.Calculator): 分子用計算機
        calculator_solid (ase.calculators.calculator.Calculator): 固体用計算機
        reactants (list[CGFEInput|Literal["CHE"]]): 反応物。"CHE"を指定すると、CHEモデルによるギブス自由エネルギーを計算する。
        products (list[CGFEInput|Literal["CHE"]]): 生成物。"CHE"を指定すると、CHEモデルによるギブス自由エネルギーを計算する。
        temperature (float): 温度[K]
        pressure (float): 圧力[Pa]
        electrode_potential (float): 電極電位（V vs SHE）
        pH (float): pH
        optimizer_cls (type[Optimizer]): 最適化エンジンのクラス
        opt_fmax (float): 最適化の閾値
        opt_maxsteps (int): 最適化の最大ステップ数
        logger (ConditionalLogger): ロガー。Noneの場合はログを出力しない
        cleanup_vibrations (bool): 振動計算後にファイルをクリーンアップするか
        copy_atoms (bool, optional): 原子構造をコピーするかどうか。Falseの場合は、入力の原子構造をそのまま使用する。デフォルトは True。

    Returns:
        float: 反応物と生成物のギブス自由エネルギーの差(ΔG)
    """
    # --- ログ設定 ---
    logger = ensure_logger("delta_g", enable_logging, logger)

    # ---(H+ + e-)のギブス自由エネルギーを計算する(CHEモデルで)
    # CHEの使用有無をチェック
    che_in_reactants = "CHE" in reactants
    che_in_products = "CHE" in products

    logger.info("CHEモデル計算開始")
    logger.info(f"反応物にCHE: {che_in_reactants}")
    logger.info(f"生成物にCHE: {che_in_products}")
    logger.info(f"電極電位: {electrode_potential:.6f} V vs SHE")
    logger.info(f"pH: {pH:.2f}")
    logger.info(f"温度: {temperature:.2f} K")

    g_che = 0.0
    if che_in_reactants or che_in_products:
        logger.info("CHEモデルによる(H+ + e-)の自由エネルギー計算")

        # CHEモデルの式: g_(H+ + e-) = 0.5 * g_(H2) - e * U + k_B * T * log(10) * pH
        logger.info("H2分子の自由エネルギー計算")
        g_h2 = calculate_gibbs_free_energy(
            calculator_molecule,  # 気相分子向けの計算機
            calculator_solid,  # 固体表面向けの計算機
            CGFEInput(
                structure=molecule("H2"),
                calc_mode="IdealGasThermo",
                vibrate_indices=None,
                geometry="linear",
                symmetry_number=2,
                spin_multiplicity=1,
            ),
            temperature=temperature,
            pressure=pressure,
            optimizer_cls=optimizer_cls,
            opt_fmax=opt_fmax,
            opt_maxsteps=opt_maxsteps,
            logger=logger,
            cleanup_vibrations=cleanup_vibrations,
            copy_atoms=copy_atoms,
        )

        # CHEモデル計算の各項
        e_charge = 1.0  # 素電荷（eV/V）
        kB = 8.617333262e-5  # ボルツマン定数（eV/K）

        term1 = 0.5 * g_h2
        term2 = -e_charge * electrode_potential
        term3 = kB * temperature * np.log(10) * pH

        g_che = term1 + term2 + term3

        logger.info("CHEモデル計算詳細:")
        logger.info(f"  G(H2): {g_h2:.6f} eV")
        logger.info(f"  0.5 * G(H2): {term1:.6f} eV")
        logger.info(f"  -e_charge * U: {term2:.6f} eV")
        logger.info(f"  kBT * ln(10) * pH: {term3:.6f} eV")
        logger.info(f"  G(H+ + e-): {g_che:.6f} eV")
    else:
        logger.info("CHEモデルは使用されません")

    # ---反応物のギブス自由エネルギーを計算する
    logger.info("反応物のギブス自由エネルギー計算開始")
    reactants_gs = []
    for i, reactant in enumerate(reactants):
        if reactant == "CHE":
            logger.info(f"反応物{i+1}: CHE (G = {g_che:.6f} eV)")
            reactants_gs.append(g_che)
        else:
            logger.info(f"反応物{i+1}: {reactant.structure.symbols}")
            g = calculate_gibbs_free_energy(
                calculator_molecule,  # 気相分子向けの計算機
                calculator_solid,  # 固体表面向けの計算機
                reactant,
                temperature=temperature,
                pressure=pressure,
                optimizer_cls=optimizer_cls,
                opt_fmax=opt_fmax,
                opt_maxsteps=opt_maxsteps,
                logger=logger,
                cleanup_vibrations=cleanup_vibrations,
                copy_atoms=copy_atoms,
            )
            reactants_gs.append(g)

    # ---生成物のギブス自由エネルギーを計算する
    logger.info("生成物のギブス自由エネルギー計算開始")
    products_gs = []
    for i, product in enumerate(products):
        if product == "CHE":
            logger.info(f"生成物{i+1}: CHE (G = {g_che:.6f} eV)")
            products_gs.append(g_che)
        else:
            logger.info(f"生成物{i+1}: {product.structure.symbols}")
            g = calculate_gibbs_free_energy(
                calculator_molecule,  # 気相分子向けの計算機
                calculator_solid,  # 固体表面向けの計算機
                product,
                temperature=temperature,
                pressure=pressure,
                optimizer_cls=optimizer_cls,
                opt_fmax=opt_fmax,
                opt_maxsteps=opt_maxsteps,
                logger=logger,
                cleanup_vibrations=cleanup_vibrations,
                copy_atoms=copy_atoms,
            )
            products_gs.append(g)

    # ---左辺右辺のギブス自由エネルギーを計算する
    g_reactants = sum(reactants_gs)
    g_products = sum(products_gs)
    g_delta = g_products - g_reactants

    logger.info("=" * 80)
    logger.info("反応ギブス自由エネルギー変化の最終計算")
    logger.info("反応物:")
    for i, g in enumerate(reactants_gs):
        species_name = (
            "CHE" if reactants[i] == "CHE" else str(reactants[i].structure.symbols)
        )
        logger.info(f"  {species_name}: {g:.6f} eV")
    logger.info(f"反応物合計: {g_reactants:.6f} eV")

    logger.info("生成物:")
    for i, g in enumerate(products_gs):
        species_name = (
            "CHE" if products[i] == "CHE" else str(products[i].structure.symbols)
        )
        logger.info(f"  {species_name}: {g:.6f} eV")
    logger.info(f"生成物合計: {g_products:.6f} eV")

    logger.info(
        f"ΔG = G(products) - G(reactants) = {g_products:.6f} - {g_reactants:.6f} = {g_delta:.6f} eV"
    )
    logger.info("=" * 80)

    return g_delta


# ----------
# ---格子定数を最適化する
# ----------
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


# 格子定数を最適化する
def optimize_lattice_constant(
    atoms: Atoms,
    calculator: Calculator | None = None,
    optimizer_cls: type[Optimizer] = FIRELBFGS,
    opt_fmax: float = 0.01,
    opt_maxsteps: int | None = None,
    copy_atoms: bool = True,
) -> LatticeConstant:
    """格子定数を最適化してLatticeConstantオブジェクトを返す

    この関数は、与えられた原子構造の格子定数を最適化し、結果をLatticeConstantデータクラス
    として返します。最適化にはUnitCellFilterを使用してセルの形状と体積の両方を最適化します。

    Args:
        atoms (Atoms): 最適化する原子構造
        calculator (Calculator | None): 使用する計算機。Noneの場合は既存の計算機を使用
        optimizer_cls (type[Optimizer]): 使用する最適化アルゴリズムクラス（デフォルト: FIRELBFGS）
        opt_fmax (float): 収束条件（原子にかかる力の最大値） [eV/Å]
        opt_maxsteps (int | None): 最大ステップ数。Noneの場合は制限なし
        copy_atoms (bool, optional): 原子構造をコピーするかどうか。Falseの場合は、入力の原子構造をそのまま使用する。デフォルトは True。
    Returns:
        LatticeConstant: 最適化後の格子定数

    Raises:
        RuntimeError: 計算機が設定されていない場合

    Example:
        >>> from ase.build import bulk
        >>> from ase.calculators.emt import EMT
        >>> atoms = bulk('Si', 'diamond', a=5.4)
        >>> atoms.calc = EMT()
        >>> result = optimize_lattice_constant(atoms, fmax=0.001)
        >>> print(result.a)  # 最適化後のa軸長
    """
    # --- 入力検証 ---
    if calculator is None and atoms.calc is None:
        raise RuntimeError(
            "計算機が設定されていません。calculatorを指定するか、atomsに計算機を設定してください。"
        )

    # --- 原子構造のコピーを作成（元の構造を保持） ---
    if copy_atoms:
        atoms_copy = atoms.copy()
    else:
        atoms_copy = atoms

    # --- 計算機の設定 ---
    if calculator is not None:
        atoms_copy.calc = calculator

    # --- UnitCellFilterでセル最適化を有効化 ---
    unit_cell_filter = UnitCellFilter(atoms_copy)

    # --- 最適化アルゴリズムの設定と実行 ---
    opt_dyn = optimizer_cls(unit_cell_filter)

    # ステップ数の制限がある場合は設定
    if opt_maxsteps is not None:
        opt_dyn.run(fmax=opt_fmax, steps=opt_maxsteps)
    else:
        opt_dyn.run(fmax=opt_fmax)

    # --- 最適化後の格子定数を取得 ---
    # get_cell_lengths_and_angles() は (a, b, c, alpha, beta, gamma) のタプルを返す
    lengths_and_angles = atoms_copy.get_cell_lengths_and_angles()

    # --- LatticeConstantオブジェクトを作成して返す ---
    return LatticeConstant(
        a=lengths_and_angles[0],
        b=lengths_and_angles[1],
        c=lengths_and_angles[2],
        alpha=lengths_and_angles[3],
        beta=lengths_and_angles[4],
        gamma=lengths_and_angles[5],
    )
