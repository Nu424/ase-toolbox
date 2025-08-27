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

- **get_neighbors(atoms, target_atom, return_type="atoms")**
  - 🧩 何をする: `natural_cutoffs` に基づく隣接原子を返す。
  - 🗺️ 場面: まずは近傍を知りたいとき。
  - 🔧 主な引数:
    - `atoms (ase.Atoms)`: 構造。
    - `target_atom (int | ase.Atom)`: 基準原子。
    - `return_type (Literal["atoms","indices"])`: 出力形式。
  - ↩️ 戻り値: `list[ase.Atom] | list[int]`。

- **separate_layers(atoms, return_type="atoms", decimals=4, sort_by_z=True)**
  - 🧩 何をする: z座標で層を検出し、層ごとに原子をグルーピング。
  - 🗺️ 場面: スラブで下層/上層に手を入れたいとき。
  - 🔧 主な引数:
    - `atoms (ase.Atoms)`: スラブ構造。
    - `return_type (Literal["atoms","indices"])`: 出力形式。
    - `decimals (int)`: z丸め桁。層の判定に影響。
    - `sort_by_z (bool)`: 下→上の順で並べるか。
  - ↩️ 戻り値: `list[list[ase.Atom]] | list[list[int]]`（`layered[0]` が最下層）。

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

- **fix_layers(atoms, fixed_layers, inplace=False, decimals=4, logger=None, enable_logging=True)**
  - 🧩 何をする: 下から `fixed_layers` 個の層に `FixAtoms` 制約を付与。
  - 🗺️ 場面: スラブ計算の下層固定。
  - 🔧 主な引数: `atoms`, `fixed_layers (int)`, `inplace`, `decimals`, `logger`, `enable_logging`。
  - ↩️ 戻り値: 制約付き `ase.Atoms`。

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
- データクラス: **CAEInput(structure, calc_mode)**, **CGFEInput(...)**, **LatticeConstant(a,b,c,alpha,beta,gamma)**

- **calculate_adsorption_energy(calculator_molecule, calculator_solid, adsorbed_structure_input, reactant_structures_input, optimizer_cls, opt_fmax, opt_maxsteps, logger=None, enable_logging=True, copy_atoms=True)**
  - 🧩 何をする: 吸着後構造と、反応物群をそれぞれ最適化→エネルギーから吸着エネルギーを返す。
  - 🗺️ 場面: 分子/固体の混在系での吸着評価（Matlantis計算を想定）。
  - 🔧 主な引数:
    - `calculator_molecule / calculator_solid (Calculator)`: 分子/固体用計算機。
    - `adsorbed_structure_input (CAEInput)`: 吸着後構造と計算モード。
    - `reactant_structures_input (list[CAEInput])`: 反応物群。
    - `optimizer_cls`, `opt_fmax`, `opt_maxsteps`, `logger`, `enable_logging`, `copy_atoms`。
  - ↩️ 戻り値: `float`（eV、負なら有利）。

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
- **resolve_target_indices(base_atoms, target)**
  - 🧩 何をする: 多様なターゲット指定（int/Atom/Atoms/list）をインデックス配列に正規化。
  - 🔧 主な引数: `base_atoms (ase.Atoms)`, `target (様々な型に対応)`。
  - ↩️ 戻り値: `list[int]`。

---

必要なときにサクッと呼べる、小回りの効くツールボックスを目指しています。まずは `ase-toolbox.ipynb` を開いて、手元の構造で試してみてください。


---
## 開発メモ
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
