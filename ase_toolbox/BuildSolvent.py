import os
import subprocess as sub
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.units import _Nav
from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass
class ComponentSpec:
    """溶媒・溶質成分の仕様.
    
    Attributes:
        name: 成分名
        concentration_mol_L: 濃度 (mol/L)
        molecule: 分子構造 (SMILES文字列、ASE Atoms、またはファイルパス)
    """
    name: str
    concentration_mol_L: float
    molecule: Union[str, Atoms, Path]


@dataclass
class CellSpec:
    """セル寸法の仕様.
    
    Attributes:
        lx: X方向の長さ (Å)
        ly: Y方向の長さ (Å)
        lz: Z方向の長さ (Å)
    """
    lx: float
    ly: float
    lz: float

    @property
    def volume_A3(self) -> float:
        """セル体積をÅ³で取得.
        
        Returns:
            セル体積 (Å³)
        """
        return self.lx * self.ly * self.lz

    @property
    def volume_L(self) -> float:
        """セル体積をリットルで取得.
        
        Returns:
            セル体積 (L)
        """
        return self.volume_A3 * 1e-27


@dataclass
class SolvationBuilder:
    """固定構造周囲への任意濃度溶媒系構築器.
    
    固定構造（スラブ、クラスター等）の周囲に、指定濃度の溶媒・溶質を
    Packmolを用いて配置するビルダークラス。水密度固定モード（Aモード）を採用。
    
    Attributes:
        cell: セル寸法仕様
        structure: 固定構造（ASE Atoms）
        components: 溶媒・溶質成分のリスト
        water_density: 水密度 (g/cm³)
        structure_position: 構造配置位置。Noneなら中心化
        margin: Packmolボックス境界からの余白 (Å)
        gap: 固定構造周囲の除外領域幅 (Å)
        tolerance: Packmol収束判定値
        packmol_bin: Packmol実行ファイル名またはパス
        filetype: 出力ファイル形式
        pbc: 周期境界条件 (x, y, z)
        outputdir: 出力ディレクトリ
        verbose: 詳細ログ出力フラグ
    """
    cell: CellSpec
    structure: Atoms
    components: list[ComponentSpec]
    water_density: float = 0.997
    structure_position: tuple[float, float, float] | None = None
    margin: float = 1.0
    gap: float = 2.0
    tolerance: float = 2.0
    packmol_bin: str = "packmol"
    filetype: str = "xyz"
    pbc: tuple[bool, bool, bool] = (True, True, True)
    outputdir: Path = field(default_factory=lambda: Path("output"))
    verbose: bool = False

    def __post_init__(self) -> None:
        """初期化後処理：出力ディレクトリ作成と入力値検証."""
        self.outputdir = Path(self.outputdir)
        if not self.outputdir.exists():
            self.outputdir.mkdir(parents=True, exist_ok=True)

        # Basic validation
        if self.cell.lx <= 0 or self.cell.ly <= 0 or self.cell.lz <= 0:
            raise ValueError("Cell dimensions must be positive")

        if self.water_density <= 0:
            raise ValueError("Water density must be positive")

        if self.margin < 0 or self.gap < 0:
            raise ValueError("Margin and gap must be non-negative")

        if len(self.structure) == 0:
            raise ValueError("Structure must contain at least one atom")

        # Validate components
        for component in self.components:
            if component.concentration_mol_L < 0:
                raise ValueError(
                    f"Concentration must be non-negative for {component.name}"
                )

        if self.verbose:
            print("SolvationBuilder initialized:")
            print(
                f"  Cell: {self.cell.lx:.1f} x {self.cell.ly:.1f} x {self.cell.lz:.1f} Å"
            )
            print(
                f"  Volume: {self.cell.volume_A3:.1f} Å³ ({self.cell.volume_L:.2e} L)"
            )
            print(f"  Structure: {len(self.structure)} atoms")
            print(f"  Components: {len(self.components)}")
            for comp in self.components:
                print(f"    {comp.name}: {comp.concentration_mol_L:.3f} mol/L")

    # ==========================================================================
    # 分子数計算
    # ==========================================================================
    
    def _calculate_water_molecules(self) -> int:
        """水密度に基づく水分子数計算（Aモード）.
        
        セル体積と指定水密度から、含有すべき水分子数を算出。
        
        Returns:
            水分子数
        """
        # 体積換算: 1 Å³ = 1e-24 cm³
        volume_cm3 = self.cell.volume_A3 * 1e-24

        # 水の質量（グラム）
        water_mass_g = self.water_density * volume_cm3

        # 水の分子量: 18.015 g/mol
        water_mw = 18.015

        # 水分子数の算出
        water_moles = water_mass_g / water_mw
        water_molecules = round(water_moles * _Nav)

        return max(0, water_molecules)

    def _calculate_component_molecules(self, component: ComponentSpec) -> int:
        """成分濃度に基づく分子数計算.
        
        指定濃度とセル体積から、該当成分の分子数を算出。
        
        Args:
            component: 成分仕様
            
        Returns:
            成分分子数
        """
        # 体積（リットル）
        volume_L = self.cell.volume_L

        # 成分のモル数
        component_moles = component.concentration_mol_L * volume_L

        # 分子数の算出
        component_molecules = round(component_moles * _Nav)

        return max(0, component_molecules)

    def _get_molecular_counts(self) -> dict[str, int]:
        """全成分の分子数集計.
        
        水および全溶媒・溶質成分の分子数を算出して辞書で返す。
        
        Returns:
            成分名をキーとした分子数辞書
        """
        counts = {"water": self._calculate_water_molecules()}

        for component in self.components:
            count = self._calculate_component_molecules(component)
            counts[component.name] = count

        if self.verbose:
            print(f"Molecular counts: {counts}")

        return counts

    # ==========================================================================
    # 構造配置
    # ==========================================================================
    
    def _place_structure(self) -> Atoms:
        """固定構造の配置.
        
        指定位置への配置、またはセル中心への自動配置を行う。
        配置後、セル境界内に収まっているかを検証。
        
        Returns:
            配置済み固定構造
        """
        structure = self.structure.copy()

        if self.structure_position is not None:
            # 指定位置への配置
            translation = np.array(self.structure_position)
            structure.translate(translation)
        else:
            # セル中心への自動配置
            center = np.array([self.cell.lx, self.cell.ly, self.cell.lz]) / 2
            current_center = structure.get_center_of_mass()
            translation = center - current_center
            structure.translate(translation)

        # セル境界内チェック
        positions = structure.get_positions()
        min_pos = positions.min(axis=0)
        max_pos = positions.max(axis=0)

        cell_bounds = np.array([self.cell.lx, self.cell.ly, self.cell.lz])

        if (min_pos < 0).any() or (max_pos > cell_bounds).any():
            warnings.warn(
                f"Structure extends outside cell bounds. "
                f"Min: {min_pos}, Max: {max_pos}, Cell: {cell_bounds}"
            )

        return structure

    def _get_exclusion_box(
        self, placed_structure: Atoms
    ) -> tuple[float, float, float, float, float, float]:
        """固定構造周囲の除外ボックス計算.
        
        固定構造のAABB（軸平行境界ボックス）にギャップを加えた
        除外領域を算出。溶媒分子はこの領域に配置されない。
        
        Args:
            placed_structure: 配置済み固定構造
            
        Returns:
            除外ボックス座標 (x_min, y_min, z_min, x_max, y_max, z_max)
        """
        positions = placed_structure.get_positions()
        min_pos = positions.min(axis=0)
        max_pos = positions.max(axis=0)

        # 構造周囲にギャップを適用
        x_min = max(0, min_pos[0] - self.gap)
        y_min = max(0, min_pos[1] - self.gap)
        z_min = max(0, min_pos[2] - self.gap)
        x_max = min(self.cell.lx, max_pos[0] + self.gap)
        y_max = min(self.cell.ly, max_pos[1] + self.gap)
        z_max = min(self.cell.lz, max_pos[2] + self.gap)

        return (x_min, y_min, z_min, x_max, y_max, z_max)

    # ==========================================================================
    # 分子構造処理
    # ==========================================================================
    
    def _molecule_to_atoms(self, molecule: Union[str, Atoms, Path]) -> Atoms:
        """分子仕様をASE Atomsオブジェクトに変換.
        
        SMILES文字列、既存のAtoms、ファイルパスのいずれかを
        ASE Atomsオブジェクトに統一的に変換する。
        
        Args:
            molecule: 分子仕様
            
        Returns:
            ASE Atomsオブジェクト
            
        Raises:
            ValueError: 無効なSMILES文字列の場合
            TypeError: サポートされていない分子型の場合
        """
        if isinstance(molecule, Atoms):
            return molecule.copy()
        elif isinstance(molecule, (str, Path)):
            # ファイルパスかどうかの判定
            if isinstance(molecule, Path) or (
                isinstance(molecule, str) and os.path.exists(molecule)
            ):
                return read(str(molecule))
            else:
                # SMILES文字列として処理
                mol = Chem.MolFromSmiles(str(molecule))
                if mol is None:
                    raise ValueError(f"Invalid SMILES string: {molecule}")

                mol = Chem.AddHs(mol)
                # 基本的な3D構造生成（最適化なし）
                AllChem.EmbedMolecule(mol)

                atoms = Atoms(
                    positions=mol.GetConformer().GetPositions(),
                    numbers=np.array([a.GetAtomicNum() for a in mol.GetAtoms()]),
                )
                return atoms
        else:
            raise TypeError(f"Unsupported molecule type: {type(molecule)}")

    def _write_water_structure(self) -> str:
        """水分子構造ファイル書き出し.
        
        標準的な水分子（H-O-H）の構造ファイルを作成。
        
        Returns:
            水分子構造ファイルのパス
        """
        # 標準的な水分子構造を定義
        water_atoms = Atoms(
            symbols=["O", "H", "H"],
            positions=[[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
        )

        water_file = self.outputdir / f"water.{self.filetype}"
        write(water_file, water_atoms)
        return str(water_file)

    def _write_component_structures(self, counts: dict[str, int]) -> dict[str, str]:
        """各成分の構造ファイル書き出し.
        
        分子数が1以上の成分について、構造ファイルを作成し、
        ファイルパスのマッピングを返す。
        
        Args:
            counts: 成分別分子数辞書
            
        Returns:
            成分名→ファイルパスのマッピング
        """
        structure_files = {}

        for component in self.components:
            if counts[component.name] > 0:
                atoms = self._molecule_to_atoms(component.molecule)
                filename = self.outputdir / f"{component.name}.{self.filetype}"
                write(filename, atoms)
                structure_files[component.name] = str(filename)

        return structure_files

    # ==========================================================================
    # Packmol入力生成・実行
    # ==========================================================================
    
    def _generate_packmol_input(
        self,
        placed_structure: Atoms,
        counts: dict[str, int],
        structure_files: dict[str, str],
        water_file: str,
    ) -> str:
        """Packmol入力ファイル生成.
        
        固定構造、水、各成分を配置するためのPackmol入力ファイルを生成。
        固定構造は固定配置、溶媒類は除外ボックス制約付きで配置。
        
        Args:
            placed_structure: 配置済み固定構造
            counts: 成分別分子数辞書
            structure_files: 成分別構造ファイルパス辞書
            water_file: 水分子構造ファイルパス
            
        Returns:
            Packmol入力ファイルパス
        """
        input_file = self.outputdir / "packmol.inp"
        output_file = self.outputdir / f"solvated.{self.filetype}"

        # 固定構造ファイルの書き出し
        structure_file = self.outputdir / f"structure.{self.filetype}"
        write(structure_file, placed_structure)

        # 除外ボックスの取得
        exclusion_box = self._get_exclusion_box(placed_structure)

        with open(input_file, "w") as f:
            # ヘッダー設定
            print(f"tolerance {self.tolerance}", file=f)
            print(f"filetype {self.filetype}", file=f)
            print(f"output {output_file}", file=f)
            print("", file=f)

            # 固定構造の配置指定
            print(f"structure {structure_file}", file=f)
            print("  number 1", file=f)
            print("  fixed 0. 0. 0. 0. 0. 0.", file=f)
            print("end structure", file=f)
            print("", file=f)

            # 水分子の配置指定
            if counts["water"] > 0:
                print(f"structure {water_file}", file=f)
                print(f"  number {counts['water']}", file=f)
                print(
                    f"  inside box {self.margin} {self.margin} {self.margin} "
                    f"{self.cell.lx - self.margin} {self.cell.ly - self.margin} "
                    f"{self.cell.lz - self.margin}",
                    file=f,
                )
                print(
                    f"  outside box {exclusion_box[0]} {exclusion_box[1]} {exclusion_box[2]} "
                    f"{exclusion_box[3]} {exclusion_box[4]} {exclusion_box[5]}",
                    file=f,
                )
                print("end structure", file=f)
                print("", file=f)

            # その他成分の配置指定
            for component in self.components:
                if counts[component.name] > 0:
                    print(f"structure {structure_files[component.name]}", file=f)
                    print(f"  number {counts[component.name]}", file=f)
                    print(
                        f"  inside box {self.margin} {self.margin} {self.margin} "
                        f"{self.cell.lx - self.margin} {self.cell.ly - self.margin} "
                        f"{self.cell.lz - self.margin}",
                        file=f,
                    )
                    print(
                        f"  outside box {exclusion_box[0]} {exclusion_box[1]} {exclusion_box[2]} "
                        f"{exclusion_box[3]} {exclusion_box[4]} {exclusion_box[5]}",
                        file=f,
                    )
                    print("end structure", file=f)
                    print("", file=f)

        return str(input_file)

    def _run_packmol(self, input_file: str) -> str:
        """Packmol実行.
        
        生成した入力ファイルを使ってPackmolを実行し、
        溶媒化システムを構築する。
        
        Args:
            input_file: Packmol入力ファイルパス
            
        Returns:
            Packmol出力ファイルパス
            
        Raises:
            RuntimeError: Packmol実行失敗または出力ファイル生成失敗時
        """
        output_file = self.outputdir / f"solvated.{self.filetype}"

        cmd = [self.packmol_bin]
        stdout = None if self.verbose else sub.DEVNULL

        with open(input_file, "rt") as f:
            result = sub.run(
                cmd, text=True, stdout=stdout, stdin=f, capture_output=not self.verbose
            )

        if result.returncode != 0:
            raise RuntimeError(f"Packmol failed with return code {result.returncode}")

        if not output_file.exists():
            raise RuntimeError(f"Packmol output file not created: {output_file}")

        return str(output_file)

    # ==========================================================================
    # メイン構築処理
    # ==========================================================================
    
    def build(self) -> Atoms:
        """溶媒化システム構築.
        
        全処理を統合し、固定構造周囲に溶媒を配置した
        最終的なASE Atomsオブジェクトを生成する。
        
        Returns:
            溶媒化されたASE Atomsオブジェクト
        """
        # 1. 分子数計算
        if self.verbose:
            print("Calculating molecular counts...")
        counts = self._get_molecular_counts()

        # 2. 固定構造の配置
        if self.verbose:
            print("Placing structure...")
        placed_structure = self._place_structure()

        # 3. 構造ファイル書き出し
        if self.verbose:
            print("Writing structure files...")
        water_file = self._write_water_structure()
        structure_files = self._write_component_structures(counts)

        # 4. Packmol入力生成
        if self.verbose:
            print("Generating Packmol input...")
        input_file = self._generate_packmol_input(
            placed_structure, counts, structure_files, water_file
        )

        # 5. Packmol実行
        if self.verbose:
            print("Running Packmol...")
        output_file = self._run_packmol(input_file)

        # 6. 結果読み込みとセル設定
        if self.verbose:
            print("Reading output and setting cell...")
        solvated_atoms = read(output_file)
        solvated_atoms.set_cell([self.cell.lx, self.cell.ly, self.cell.lz])
        solvated_atoms.set_pbc(self.pbc)

        # 7. 最終ファイル書き出し
        final_output = self.outputdir / f"solvated_final.{self.filetype}"
        write(final_output, solvated_atoms)

        if self.verbose:
            print(f"Solvated system created with {len(solvated_atoms)} atoms")

        return solvated_atoms


# ==========================================================================
# パブリックAPI
# ==========================================================================


def build_solvated_system(
    cell: tuple[float, float, float],
    structure: Atoms,
    components: list[ComponentSpec],
    water_density: float = 0.997,
    structure_position: tuple[float, float, float] | None = None,
    margin: float = 1.0,
    gap: float = 2.0,
    tolerance: float = 2.0,
    packmol_bin: str = "packmol",
    filetype: str = "xyz",
    pbc: tuple[bool, bool, bool] = (True, True, True),
    outputdir: Union[str, Path] = "output",
    verbose: bool = False,
) -> Atoms:
    """固定構造周囲への任意濃度溶媒系構築.

    固定構造（スラブ、クラスター等）の周囲に、指定濃度の溶媒・溶質を
    Packmolを用いて配置する。水密度固定モード（Aモード）を採用し、
    水の密度から水分子数を決定し、その水の体積を基準に各成分の濃度から
    分子数を算出する。

    Args:
        cell: セル寸法 (lx, ly, lz) [Å]
        structure: 固定構造（スラブ、クラスター等）
        components: 溶媒・溶質成分のリスト
        water_density: 水密度 [g/cm³] (デフォルト: 0.997)
        structure_position: 構造配置位置。Noneなら中心化
        margin: Packmolボックス境界からの余白 [Å]
        gap: 固定構造周囲の除外領域幅 [Å]
        tolerance: Packmol収束判定値
        packmol_bin: Packmol実行ファイル名またはパス
        filetype: 出力ファイル形式 ("xyz", "pdb"等)
        pbc: 周期境界条件 (x, y, z)
        outputdir: 出力ディレクトリ
        verbose: 詳細ログ出力フラグ

    Returns:
        溶媒化されたASE Atomsオブジェクト

    Example:
        >>> from ase import Atoms
        >>> from solvent import build_solvated_system, ComponentSpec
        >>>
        >>> # 簡単なスラブを作成
        >>> slab = Atoms('Au4', positions=[[0,0,0], [2,0,0], [0,2,0], [2,2,0]])
        >>>
        >>> # 溶媒成分を定義
        >>> components = [
        >>>     ComponentSpec("NaCl", 0.1, "Cl[Na]"),  # 0.1 M NaCl
        >>>     ComponentSpec("ethanol", 0.05, "CCO")   # 0.05 M ethanol
        >>> ]
        >>>
        >>> # 溶媒化システムを構築
        >>> result = build_solvated_system(
        >>>     cell=(20, 20, 20),
        >>>     structure=slab,
        >>>     components=components,
        >>>     verbose=True
        >>> )
    """
    cell_spec = CellSpec(lx=cell[0], ly=cell[1], lz=cell[2])

    builder = SolvationBuilder(
        cell=cell_spec,
        structure=structure,
        components=components,
        water_density=water_density,
        structure_position=structure_position,
        margin=margin,
        gap=gap,
        tolerance=tolerance,
        packmol_bin=packmol_bin,
        filetype=filetype,
        pbc=pbc,
        outputdir=Path(outputdir),
        verbose=verbose,
    )

    return builder.build()
