# ase-toolbox
ASEを使った化学シミュレーションをサクッと進めるためのヘルパー関数集です。配位数や隣接探索、層の分離・固定、原子の置換、吸着エネルギーやNEB、溶媒化配置など、よくある作業を簡単に呼び出せます。

## 概要

### 想定環境とインストール
- 想定環境: Matlantis 環境での使用を想定しています（`pfp_api_client`, `matlantis_features` を利用する関数があります）。
- インストール
  1. このリポジトリ内の、`ase_toolbox`フォルダ、`try-ase-toolbox.ipynb`を、Matlantis環境にコピーする。
  2. Matlantisの各種ノートブックで、以下のようにインポートして使う。
    ```python
    from ase_toolbox.XXX import XXX # XXXはモジュール名・関数名
    ```
- 補足: 一部の機能では内部的に ASE/NumPy/Matplotlib に加え、RDKit（SMILES→3D生成）や Packmol（溶媒配置）を利用します。必要に応じて環境に合わせてご利用ください。

### 使い方（ノートブック）
使い方の流れや実行例は、同梱の `try-ase-toolbox.ipynb` をご覧ください。段階的に動かしながら理解できる構成になっています。

### このツールボックスの読み方（共通ポリシー）
- 🔁 **原子指定は柔軟に**: 多くの関数で、原子指定は `int`（インデックス）と `ase.Atom` の両方に対応しています。
- 🎛️ **出力形式は選べます**: 隣接原子などの返り値は、`return_type="atoms" | "indices"` で `list[ase.Atom]` または `list[int]` を選べます。
- 📝 **ログ**: いくつかの計算関数は簡易ロガー（`util.ConditionalLogger`）を使えます。必要に応じてON/OFFできます。

---

## API クイックガイド（目的・場面・引数・戻り値）
「何をする関数か」「どんな場面で使うか」「主な引数」「戻り値」を、できるだけコンパクトにまとめています。詳細は各関数のdocstringと `ase-toolbox.ipynb` をどうぞ。

### CalcValue.py（値の計算）
- **coordination_number(atoms, target_atom, return_type="atoms", cutoffs=None, cutoff_scaling=1.0)**
  - 🧩 何をする: 指定原子の配位数と、その隣接原子を返します。
  - 🗺️ 使う場面: クラスターや表面の局所配位の把握に。
  - 🔧 主な引数:
    - `atoms (ase.Atoms)`: 計算対象の構造。
    - `target_atom (int | ase.Atom)`: 配位数を調べる対象原子。
    - `return_type (Literal["atoms","indices"])`: 隣接リストの型。
    - `cutoffs (Sequence[float] | None)`: 各原子のカットオフ半径。`None` なら `natural_cutoffs(atoms)`。
    - `cutoff_scaling (float)`: `natural_cutoffs` に掛ける倍率（1.0でそのまま）。
  - ↩️ 戻り値: `(coord_num: int, neighbors: list[ase.Atom] | list[int])`。

### FindAtoms.py（原子の探索）
- **find_atom_by_index(atoms, index)**
  - 🧩 何をする: インデックス指定で原子を1つ返す。
  - 🗺️ 場面: 単純な参照をサッと。
  - 🔧 主な引数:
    - `atoms (ase.Atoms)`: 検索対象の構造。
    - `index (int)`: 取得したい原子のインデックス。
  - ↩️ 戻り値: `ase.Atom`。

- **find_indices_by_symbol(atoms, symbol)**
  - 🧩 何をする: 元素記号に一致する原子のインデックスを全部返す。
  - 🗺️ 場面: Oだけ拾う、などの前処理。
  - 🔧 主な引数: `atoms (ase.Atoms)`, `symbol (str)`（例: "O", "Fe"）。
  - ↩️ 戻り値: `list[int]`。

- **filter_by_symbols(atoms, symbols, return_type="atoms")**
  - 🧩 何をする: 指定した元素記号（単一/複数）に一致する原子をフィルタリングして返す。
  - 🗺️ 場面: 複数の元素を同時に抽出したい、柔軟な出力形式が欲しい。
  - 🔧 主な引数:
    - `atoms (ase.Atoms)`: フィルタリング対象の構造。
    - `symbols (str | Sequence[str])`: 元素記号。単一("O")または複数(["Cu", "O"])。
    - `return_type (Literal["atoms","indices"])`: 出力形式。
  - ↩️ 戻り値: `list[ase.Atom] | list[int]`。
  - 📝 メモ: 大文字小文字は自動正規化。元の構造の順序を保持。

- **get_neighbors(atoms, target_atom, return_type="atoms")**
  - 🧩 何をする: `natural_cutoffs` に基づく隣接原子を返す。
  - 🗺️ 場面: まずは近傍を知りたいとき。
  - 🔧 主な引数:
    - `atoms (ase.Atoms)`: 構造。
    - `target_atom (int | ase.Atom)`: 基準原子。
    - `return_type (Literal["atoms","indices"])`: 出力形式。
  - ↩️ 戻り値: `list[ase.Atom] | list[int]`。

- **separate_layers(atoms, return_type="atoms", decimals=4, sort_by_z=True, use_substrate_mask="auto")**
  - 🧩 何をする: z座標で層を検出し、層ごとに原子をグルーピング。
  - 🗺️ 場面: スラブで下層/上層に手を入れたいとき。
  - 🔧 主な引数:
    - `atoms (ase.Atoms)`: スラブ構造。
    - `return_type (Literal["atoms","indices"])`: 出力形式。
    - `decimals (int)`: z丸め桁。層の判定に影響。
    - `sort_by_z (bool)`: 下→上の順で並べるか。
    - `use_substrate_mask (Literal["auto",True,False])`: 基板マスクの使用設定。"auto" の場合、`is_substrate` マスクが存在すれば基板原子のみで層を検出。デフォルトは "auto"。
  - ↩️ 戻り値: `list[list[ase.Atom]] | list[list[int]]`（`layered[0]` が最下層）。
  - 📝 メモ: 吸着分子配置後も正しく基板の層を検出するには、事前に `set_substrate_mask_all()` でマスクを設定してください。

- **classify_surface_atoms(atoms, return_type="atoms", upper_tolerance=3)**
  - 🧩 何をする: 配位数の低い原子を「表面」、それ以外を「内側」に分類。
  - 🗺️ 場面: クラスターの表面原子をざっくり抽出。
  - 🔧 主な引数: `atoms (ase.Atoms)`, `return_type`, `upper_tolerance (int)`。
  - ↩️ 戻り値: `(surface, inner)` いずれも `list[ase.Atom] | list[int]`。

- **find_central_atom(atoms_or_list, return_type="atom")**
  - 🧩 何をする: xy面の重心に最も近い原子を返す（`Atoms` または `list[Atom]`対応）。
  - 🗺️ 場面: 「中央っぽい原子」を基準にしたいとき。
  - 🔧 主な引数: `atoms_or_list (ase.Atoms | list[ase.Atom])`, `return_type (Literal["atom","index"])`。
  - ↩️ 戻り値: `ase.Atom | int`（`return_type`で切替）。

- **get_appended_atom_indices(before_atoms, after_atoms)**
  - 🧩 何をする: 結合後の構造で「追加された原子」のインデックスを返す。
  - 🗺️ 場面: 2つの構造を足した後、追加部分だけ扱いたい。
  - 🔧 主な引数: `before_atoms (ase.Atoms)`, `after_atoms (ase.Atoms)`。
  - ↩️ 戻り値: `list[int]`。

- **get_neighbors_with_coordination_condition(atoms, target_atom, return_type="atoms", cutoffs=None, cutoff_scaling=1.0, upper_tolerance=1, lower_tolerance=1)**
  - 🧩 何をする: 対象原子の配位数±許容幅に収まる配位数をもつ隣接原子だけを抽出。
  - 🗺️ 場面: 局所環境が似た原子だけを拾いたい。
  - 🔧 主な引数: `atoms`, `target_atom`, `return_type`, `cutoffs`, `cutoff_scaling`, `upper_tolerance`, `lower_tolerance`。
  - ↩️ 戻り値: `list[ase.Atom] | list[int]`。

### HandleAtoms.py（原子操作）
- **smiles_to_atoms(smiles, optimize="UFF", random_seed=None)**
  - 🧩 何をする: SMILES文字列からASE Atomsオブジェクトを生成する。
  - 🗺️ 場面: 分子構造を文字列から手軽に作りたい。吸着分子の初期構造作成など。
  - 🔧 主な引数:
    - `smiles (str)`: SMILES文字列（例: "CCO", "c1ccccc1"）。
    - `optimize (Literal["UFF","MMFF"]|None)`: 力場最適化の種類。"UFF"/"MMFF"/None。デフォルトは"UFF"。
    - `random_seed (int|None)`: 3D座標埋め込み時の乱数シード。再現性確保用。
  - ↩️ 戻り値: `ase.Atoms`（3D構造、単位Å）。
  - 📝 メモ: RDKit依存（遅延import）。未導入時は明確なエラー。ETKDG法で3D座標生成。水素は自動付加。

- **set_substrate_mask_all(atoms, is_substrate=True, inplace=True)**
  - 🧩 何をする: 全原子に `is_substrate` マスクを設定する。
  - 🗺️ 場面: 吸着分子を複数配置する前に、基板原子をマークしておく。これにより、後続の層検出や高さ基準が基板のみに基づいて決定されます。
  - 🔧 主な引数:
    - `atoms (ase.Atoms)`: マスクを設定する原子構造。
    - `is_substrate (bool)`: 設定する値。True で基板、False で非基板。デフォルトは True。
    - `inplace (bool)`: True の場合は atoms を直接変更。False の場合はコピーを作成。デフォルトは True。
  - ↩️ 戻り値: マスクが設定された `ase.Atoms`。
  - 📝 メモ: 複数の吸着分子を配置する場合、最初に `set_substrate_mask_all(slab)` で基板をマークすることで、`separate_layers()` や `place_adsorbate_on_surface()` が既存の吸着分子を無視し、常に基板のみから層や高さを決定します。

- **move_atoms(base_structure, target, direction, distance, inplace=False)**
  - 🧩 何をする: 指定原子（複数指定OK）を、与えた方向へ距離だけ平行移動。
  - 🗺️ 場面: 手動でちょっと動かしたい・探索したい。
  - 🔧 主な引数:
    - `base_structure (ase.Atoms)`: 操作対象の構造。
    - `target (int|Atom|Atoms|list[int]|list[Atom])`: 移動対象。
    - `direction (array-like[3])`: 方向ベクトル（内部で正規化）。
    - `distance (float)`: 平行移動距離 [Å]。
    - `inplace (bool)`: 直接書換えるか（Falseでコピー返却）。
  - ↩️ 戻り値: 変更後 `ase.Atoms`（`inplace=True`なら引数のまま）。

- **fix_layers(atoms, fixed_layers, inplace=False, decimals=4, logger=None, enable_logging=True, use_substrate_mask="auto")**
  - 🧩 何をする: 下から `fixed_layers` 個の層に `FixAtoms` 制約を付与。
  - 🗺️ 場面: スラブ計算の下層固定。
  - 🔧 主な引数: `atoms`, `fixed_layers (int)`, `inplace`, `decimals`, `logger`, `enable_logging`, `use_substrate_mask ("auto"|True|False)`。
  - ↩️ 戻り値: 制約付き `ase.Atoms`。
  - 📝 メモ: `use_substrate_mask="auto"` の場合、`is_substrate` マスクが存在すれば基板原子のみで層を検出します。

- **substitute_elements(atoms, target, new, inplace=False, seed=None)**
  - 🧩 何をする: 指定原子を新しい元素に置換。`new` は単一記号 or 組成辞書（合計1）。
  - 🗺️ 場面: ドーピング、ランダム置換。
  - 🔧 主な引数: `atoms`, `target`, `new (str | Mapping[str,float])`, `inplace`, `seed`。
  - ↩️ 戻り値: 置換後 `ase.Atoms`。

- **compute_surface_normal(atoms, target_atom, include_target=True, reference_vector=None, normalize=True, return_plane=False)**
  - 🧩 何をする: 対象原子近傍をPCAで局所平面近似し、法線ベクトルを返す。
  - 🗺️ 場面: 「外向き」方向を知りたい（`reference_vector`で向きを安定化）。
  - 🔧 主な引数: `atoms`, `target_atom`, `include_target`, `reference_vector (ndarray|None)`, `normalize`, `return_plane`。
  - ↩️ 戻り値: `normal: ndarray(3,)` または `(normal, centroid, d)`。

- **place_adsorbate_along_normal(substrate, adsorbate, target_atom, distance, upper_tolerance=1, lower_tolerance=1)**
  - 🧩 何をする: 局所法線に +z を合わせるよう吸着分子を回転・配置し、基板と結合。
  - 🗺️ 場面: まずは「自然な初期配置」を素早く作る。
  - 🔧 主な引数: `substrate (ase.Atoms)`, `adsorbate (ase.Atoms)`, `target_atom`, `distance (float)`, `upper_tolerance`, `lower_tolerance`。
  - ↩️ 戻り値: 結合後 `ase.Atoms`。

- **place_adsorbate_on_surface(substrate, adsorbate, target_atom, height, position, rotation_deg=None, align_vector=None, rotate_about="com", separate_layers_decimals=4, allow_search_surface_atom=True, inplace=False, use_substrate_mask="auto")**
  - 🧩 何をする: 指定した構造表面に、吸着分子を配置する。add_adsorbate()の高性能なラッパー関数。回転機能付き。
  - 🗺️ 場面: 表面に吸着分子を配置したい。特定の方向や角度で配置したい。
  - 🔧 主な引数:
    - `substrate (ase.Atoms)`, `adsorbate (ase.Atoms)`, `target_atom`, `height (float)`, `position (Literal["top", "bridge", "hollow"])`
    - `rotation_deg (tuple[float,float,float]|None)`: オイラー角回転(rx,ry,rz)[度]。XYZ順に適用。Noneなら回転なし。
    - `align_vector (Sequence[float]|None)`: 吸着分子の整列方向ベクトル。指定時、このベクトルを+z軸に整列。Noneなら整列なし。
    - `rotate_about (Literal["com","cog"])`: 回転中心。"com"=質量中心、"cog"=幾何中心。デフォルトは"com"。
    - `separate_layers_decimals`, `allow_search_surface_atom`, `inplace`, `use_substrate_mask ("auto"|True|False)`
  - ↩️ 戻り値: 結合後 `ase.Atoms`。
  - 📝 メモ: 
    - 回転順序: align_vector整列 → rotation_deg回転 → 位置決定。併用可能。
    - 複数の吸着分子を配置する場合は、最初に `set_substrate_mask_all(substrate)` でマスクを設定してください。これにより、2つ目以降の吸着分子配置時も、層検出と高さ基準が常に基板のみに基づいて決定されます。`use_substrate_mask="auto"` がデフォルトで、マスクが存在すれば自動的に使用されます。

- **mix_lattice_constant(composition, lattice_map=None, method="vegard", tol=1e-6, return_detail=False)**
  - 🧩 何をする: 組成 `{symbol: fraction}` から混合格子定数 `a` を計算。
  - 🗺️ 場面: `surface()` でスラブを作る前に `a` を決めたいとき（Vegard/体積混合）。
  - 🔧 主な引数: `composition (dict[str,float])`, `lattice_map`, `method ("vegard"|"volume")`, `tol`, `return_detail`。
  - ↩️ 戻り値: `float`（Å）。`return_detail=True` で `(a, detail_dict)`。
  - 📝 メモ: 体積混合は立方晶想定。bcc/hcp を含む場合は注意喚起を出します。

### BuildSolvent.py（溶媒化の構築）
- **ComponentSpec(name, concentration_mol_L, molecule)**（dataclass）
  - 🧩 何をする: 成分名・濃度（mol/L）・分子指定（SMILES/Atoms/ファイルパス）を保持。
  - 🔧 主なフィールド: `name (str)`, `concentration_mol_L (float)`, `molecule (str|Atoms|Path)`。
- **CellSpec(lx, ly, lz)**（dataclass）
  - 🧩 何をする: セル寸法と、`volume_A3`・`volume_L` を提供。
  - 🔧 主なフィールド: `lx (float)`, `ly (float)`, `lz (float)`。
- **build_solvated_system(cell, structure, components, water_density=0.997, ..., filetype="xyz", pbc=(True,True,True), outputdir="output", verbose=False)**
  - 🧩 何をする: 固定構造の周囲に、水（密度基準）と任意成分（濃度基準）を Packmol で配置。
  - 🗺️ 場面: スラブ/クラスターを所望の溶媒中に初期配置したい。
  - 🔧 主な引数:
    - `cell (tuple[float,float,float])`: セル寸法（Å）。
    - `structure (ase.Atoms)`: 固定構造。
    - `components (list[ComponentSpec])`: 溶媒/溶質の仕様。
    - `water_density (float)`: g/cm³。
    - `structure_position (tuple[float,float,float] | None)`: 明示配置位置。
    - `margin/gap/tolerance (float)`: Packmol関連設定。
    - `packmol_bin (str)`: 実行コマンド名。
    - `filetype (str)`: 出力形式（xyz, pdbなど）。
    - `pbc (tuple[bool,bool,bool])`: 周期境界条件。
    - `outputdir (str|Path)`, `verbose (bool)`。
  - ↩️ 戻り値: 溶媒化済み `ase.Atoms`（セル・PBC設定済み）。
  - 📝 メモ: SMILES→3D には RDKit、分子配置には Packmol を利用。

### Calculation.py（エネルギー・NEB・熱化学）
- データクラス: **CAEInput(structure, calculator=None, energy_override=None, coefficient=1.0)**, **CAEOutput(...)**, **CGFEInput(...)**, **LatticeConstant(a,b,c,alpha,beta,gamma)**

- **calculate_adsorption_energy(adsorbed_structure_input, reactant_structures_input, optimizer_cls=FIRELBFGS, opt_fmax=0.05, opt_maxsteps=3000, logger=None, enable_logging=True, copy_atoms=True) -> CAEOutput**
  - 🧩 何をする: 吸着後構造と反応物群をそれぞれ最適化（またはエネルギーを直接指定）し、係数付きで吸着エネルギーを算出する。
  - 🗺️ 場面: 分子/固体の混在系での吸着評価（Matlantis計算を想定）。
  - 🔧 主な引数:
    - `adsorbed_structure_input (CAEInput)`: 吸着後構造。`calculator` または `energy_override` のいずれかを指定。
    - `reactant_structures_input (list[CAEInput])`: 反応物群。各項目で `calculator` / `energy_override` / `coefficient` を個別指定可能。
    - `optimizer_cls`, `opt_fmax`, `opt_maxsteps`, `logger`, `enable_logging`, `copy_atoms`。
  - ↩️ 戻り値: `CAEOutput`。吸着エネルギー（`adsorption_energy`）に加え、最適化後構造や個別エネルギー、係数情報を保持。
  - 📝 ヒント:
    - `energy_override`: 既知のエネルギーを再利用する場合に指定（最適化をスキップ）。
    - `coefficient`: 反応式係数を設定（例: 0.5×H2 を表現する場合は `coefficient=0.5`）。
    - `CAEOutput` には `optimized_adsorbed`, `optimized_reactants`, `reactant_weighted_energies` などが含まれるため、後処理で利用しやすい。
  - 💡 使用例:
    ```python
    from ase_toolbox.Calculation import CAEInput, calculate_adsorption_energy

    adsorbed_input = CAEInput(structure=cu_co_adsorbed, calculator=calc_solid)
    reactant_inputs = [
        CAEInput(structure=cu_surface, calculator=calc_solid),
        CAEInput(structure=co_molecule, calculator=calc_molecule, coefficient=0.5),
        CAEInput(
            structure=h2_reference,
            energy_override=precomputed_e_h2,
            coefficient=0.5,
        ),
    ]

    result = calculate_adsorption_energy(
        adsorbed_structure_input=adsorbed_input,
        reactant_structures_input=reactant_inputs,
        opt_fmax=0.05,
        opt_maxsteps=3000,
        enable_logging=True,
    )

    print(f"吸着エネルギー: {result.adsorption_energy:.3f} eV")
    ```

- **analyze_composition(atoms)** / **generate_reference_structure(element, crystal_structure="auto", lattice_parameter=None, ...)**
  - 🧩 何をする: 元素組成の辞書作成 / 純元素参照構造（fcc/bcc/hcp自動判別も可）の生成。
  - 🔧 引数の例: `element (str)`, `crystal_structure ("auto"|"fcc"|"bcc"|"hcp")`, `lattice_parameter (float|None)`。

- **calculate_formation_energy(calculator, compound_structure, optimizer_cls, opt_fmax, opt_maxsteps, reference_crystal_structures=None, reference_lattice_parameters=None, logger=None, enable_logging=True, copy_atoms=True)**
  - 🧩 何をする: 化合物のエネルギーと、純元素参照エネルギー（原子あたり）から生成エネルギー。
  - 🔧 主な引数: `calculator`, `compound_structure (ase.Atoms)`, 参照構造の上書き辞書など, `copy_atoms`。
  - ↩️ 戻り値: `float`（eV、負なら形成有利）。

- **run_neb(init_atoms, final_atoms, num_intermediate_images, optimizer_cls, estimator, fmax=0.05, steps=500, trajectory_path=None, pre_align=True, k=0.1, climb=True, parallel=False, mic=None, interpolate_kwargs=None)**
  - 🧩 何をする: NEBを実行して全画像とエネルギーを返す（Matlantisの `Estimator` を想定）。
  - 🔧 主な引数: `init_atoms`, `final_atoms`, `num_intermediate_images (int)`, `optimizer_cls`, `estimator`, `fmax`, `steps`, `trajectory_path`, `pre_align`, `k`, `climb`, `parallel`, `mic`, `interpolate_kwargs`。
  - ↩️ 戻り値: `(images: list[ase.Atoms], energies: list[float])`。

- **plot_energy_profile(energies, ax=None, xlabel="replica", ylabel="energy [eV]", title=None, show=True)**
  - 📈 何をする: エネルギープロファイルを簡単プロット。
  - 🔧 主な引数: `energies (Sequence[float])`, 軸ラベル、`show`。
  - ↩️ 戻り値: `(fig, ax)`。

- **compute_barriers(energies)**
  - 🧮 何をする: 遷移状態インデックスと順逆の活性化エネルギーを返す。
  - 🔧 主な引数: `energies (Sequence[float])`。
  - ↩️ 戻り値: `(ts_index: int, e_forward: float, e_backward: float)`。

- **calculate_gibbs_free_energy(calculator_molecule, calculator_solid, calc_input, temperature=298.15, pressure=101325.0, optimizer_cls, opt_fmax, opt_maxsteps, logger=None, enable_logging=True, cleanup_vibrations=True, copy_atoms=True)**
  - 🧩 何をする: 構造最適化＋振動解析→IdealGasThermo/HarmonicThermoで G（またはF）を評価。
  - 🔧 主な引数:
    - `calculator_molecule / calculator_solid (Calculator)`。
    - `calc_input (CGFEInput)`: 振動対象・モードなどを含む入力。
    - `temperature (K)`, `pressure (Pa)`, `optimizer_cls`, `opt_fmax`, `opt_maxsteps`, `cleanup_vibrations`, `copy_atoms` ほか。
  - ↩️ 戻り値: `float`（ギブス自由エネルギー。Δではなく個別G）。

- **calculate_delta_g(calculator_molecule, calculator_solid, reactants, products, temperature=298.15, pressure=101325.0, electrode_potential=0.0, pH=7.0, optimizer_cls, opt_fmax, opt_maxsteps, logger=None, enable_logging=True, cleanup_vibrations=True, copy_atoms=True)**
  - 🧩 何をする: 反応物と生成物の総G差（ΔG）を返す。`"CHE"` 指定でCHEモデル（0.5·G(H2) − e·U + kBT·ln10·pH）。
  - 🔧 主な引数: `reactants/products (list[CGFEInput | "CHE"])`, `electrode_potential (V vs SHE)`, `pH`, 温度・圧力など, `copy_atoms`。
  - ↩️ 戻り値: `float`（eV）。

- **optimize_lattice_constant(atoms, calculator=None, optimizer_cls=FIRELBFGS, opt_fmax=0.01, opt_maxsteps=None, copy_atoms=True)**
  - 🧩 何をする: `UnitCellFilter` でセル形状・体積を最適化し、格子定数を返す。
  - 🔧 主な引数: `atoms (ase.Atoms)`, `calculator (Calculator|None)`, `optimizer_cls`, `opt_fmax`, `opt_maxsteps`, `copy_atoms`。
  - ↩️ 戻り値: `LatticeConstant`（a,b,c,α,β,γ）。

### util.py（ユーティリティ）
- **ConditionalLogger / ensure_logger / setup_logger**
  - 🧩 何をする: ログ出力を簡単にON/OFFしつつ、ファイル/コンソールへ整形出力。
- **optimize_and_get_energy(atoms, calculator, optimizer_cls, fmax, maxsteps, label, logger, copy_atoms=True)**
  - 🧩 何をする: 構造最適化を実行し、最終エネルギーを返す（ログ込み）。
  - 🔧 主な引数: `atoms`, `calculator`, `optimizer_cls`, `fmax`, `maxsteps`, `label`, `logger`, `copy_atoms`。
  - ↩️ 戻り値: `float`（eV）。
  - 📝 メモ: 原子数が1の場合は最適化しても常に0 eVとなるため、警告ログが出力される。
- **resolve_target_indices(base_atoms, target)**
  - 🧩 何をする: 多様なターゲット指定（int/Atom/Atoms/list）をインデックス配列に正規化。
  - 🔧 主な引数: `base_atoms (ase.Atoms)`, `target (様々な型に対応)`。
  - ↩️ 戻り値: `list[int]`。
- **sanitize_atoms_for_xyz_write(atoms)**
  - 🧩 何をする: XYZ形式で安全に書き出すために、Atomsオブジェクトをクリーンアップする。
  - 🗺️ 場面: 配列長不一致エラーやシリアライズエラーを回避してファイル書き出ししたい。
  - 🔧 主な引数: `atoms (ase.Atoms)`: クリーンアップ対象の原子構造。
  - ↩️ 戻り値: `ase.Atoms`（元素記号・座標・セル・PBCのみ保持。tagsとinfo単純型は条件付き保存）。
  - 📝 メモ: 元のAtomsは変更しない。per-atom配列の長さ不一致を自動除去。

---

必要なときにサクッと呼べる、小回りの効くツールボックスを目指しています。まずは `ase-toolbox.ipynb` を開いて、手元の構造で試してみてください。


---
## 開発メモ
### 新しい関数を作ったときの流れ
1. 実装する
2. Matlantis環境で動作検証する
3. ドキュメントを更新する
   - モジュールの一番上のコメント部分(docstring?)
   - README.mdの関数一覧
   - try-ase-toolbox.ipynbのおためしコード

### ヘルパー関数を作ってもらうときのプロンプト例
```markdown
PythonのASEで使用できるヘルパー関数を用意し、シミュレーションを効率的に実装できるようにしようと考えています。

ヘルパー関数の実装では、以下の点に気をつけたいです
- 日本語かつGoogleスタイルのdocstringをつける
- 処理のまとまりごとに適度なコメントをつける
- 適切な型ヒントをつける(3.11以降のベストな型ヒントをつけてください)
- 引数で原子を指定する場合、ase.Atomと原子のインデックスの両方を受け付けられるようにする
- 私が実行するので、あなたがテストをする必要はない。実行結果を後でフィードバックします。
```　
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
```　
- 原子を出力する場合、return_type: str = Literal["atoms", "indices"]引数で、list[Atoms]とlist[int]を指定して出力できるようにする
```　
# --- 出力形式に応じて返す ---
if return_type == "atoms":
    return [atoms[i] for i in atom_indices]
else:
    return list(atom_indices)
```　
----------

```

### シミュレーション用のコードを作成してもらうときのプロンプト例
```
PythonのASEで化学シミュレーションをします。
「<概要を書く>」
シミュレーションの流れは以下のようなものを考えています。
1. <手順を書く>

## ポイント
- 使用可能な関数については、README.mdにまとまっています。使用できる場合はこれを優先して使い、再実装を避けてください。
- 実験条件や計算結果などは、可能な限りログとして保存するようにしてください。その際、「{実験名}_{indexなど}.json」のようなJSON形式でまとめてください。
- 構造も、xyz形式などで適度に保存してください。
- 実際の計算は、Jupyter Notebook環境で実行します。インタラクティブセル(#%%)を用い、実験のワークフローに合わせて、適度な粒度で実装してください。

まずは、実装計画を立ててください。

```