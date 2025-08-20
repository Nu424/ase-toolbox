# %% [markdown]
# # お道具箱
# ASEによる化学シミュレーションでよく使うようなコードをまとめておくやつ

# %% [markdown]
# ## 🔴基本的なコード

# %% [markdown]
# ### 🟡Matlantis固有のコード

# %%
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

# バージョンによって計算結果が異なる場合があるため、毎回バージョンを確認する
print(f"pfp_api_client: {pfp_api_client.__version__}") 

# ASECalculatorを使用するための設定
# EstimatorCalcModeは、以下のように使い分ける
# - 一般の系： EstimatorCalcMode.CRYSTAL_U0 Uなしモード
# - 酸化物など： EstimatorCalcMode.CRYSTAL　Uありモード
# - 単体有機分子： EstimatorCalcMode.MOLECULE 分子モード
estimator = Estimator(calc_mode=EstimatorCalcMode.MOLECULE)
calculator = ASECalculator(estimator) # このcalculatorをatoms.calcに設定して使用する

# %%
calc_copy=calculator.copy()

# %% [markdown]
# ### 🟡matplotlibで日本語を表示する

# %%
# ---matplotlibで日本語を表示する
plt.rcParams['font.family'] = ['DejaVu Sans', 'Hiragino Sans', 'Yu Gothic', 'Meiryo', 'Takao', 'IPAexGothic', 'IPAPGothic', 'VL PGothic', 'Noto Sans CJK JP']
# Windows環境での代替設定
try:
    import matplotlib.font_manager as fm
    # 日本語フォントを探す
    font_list = [f.name for f in fm.fontManager.ttflist if 'gothic' in f.name.lower() or 'mincho' in f.name.lower() or 'meiryo' in f.name.lower()]
    if font_list:
        plt.rcParams['font.family'] = font_list[0]
    else:
        # フォールバック: Unicode対応
        plt.rcParams['font.family'] = 'DejaVu Sans'
except:
    pass

plt.rcParams['font.family']

# %% [markdown]
# ### 🟡可視化

# %%
from pfcc_extras.visualize.view import view_ngl

view_ngl(slab, representations=["ball+stick"], w=400, h=300)

# %% [markdown]
# ### 🟡最適化

# %% [markdown]
# #### 最適化アルゴリズムについて
# - 局所最適へのたどり着き方(=アルゴリズム)は、いろいろある
# - 選び方
#     - 基本的に`LBFGSLineSearch`か`FIRE`を使えば良いらしい
#         - `LBFGSLineSearch`がうまくいかなかったら`FIRE`を使う、みたいな
#     - 選ぶのが面倒なら、上の2つを組み合わせた`FIRELBFGS`を使うと良い
#         - Matlantisオリジナルの最適化アルゴリズム
#         - 使用方法は、2つ下のセルに示す
#     - 局所最適を探すアルゴリズムの違いなだけなので、シミュレーションに大きく差が生じることはないのかもしれない
# - 参考文献
#     - https://docs.matlantis.com/atomistic-simulation-tutorial/ja/2_3_opt-algorithm.html

# %%
from ase.optimize import FIRE, LBFGSLineSearch

opt=FIRE(atoms)
opt.run(fmax=0.001)

# %%
# FIRELBFGS 最適化アルゴリズムの使い方(Matlantisオリジナル)
from matlantis_features.ase_ext.optimize import FIRELBFGS

opt = FIRELBFGS(atoms)
opt.run(fmax=0.001)

# %% [markdown]
# ### 🟡書き出し

# %%
# 内蔵ライブラリからcifファイルを作成する

from ase.build import molecule
from ase.io import write

# 出力先ディレクトリ（必要に応じて変更）
output_dir = './'

# 分子リスト
species = ['CO2', 'CO', 'H2', 'H2O', 'O2']

for sym in species:
    mol = molecule(sym)           # ASE 内蔵の分子ライブラリから生成
    filename = f'{output_dir}{sym}.cif'
    write(filename, mol)          # CIF 形式で保存
    print(f'Written: {filename}')


# %% [markdown]
# ## 🔴構造を作る

# %% [markdown]
# ### 🟡基本的な構造

# %%
from ase.build import molecule, bulk, surface
from ase.io import read, write

# ---分子
co=molecule('CO')

# ---結晶
cu=bulk('Cu')
# 増やす
_cu=cu*(2,3,4)

# ---表面
slab=surface(
    lattice=cu,
    indices=(1,1,1),
    layers=2,
    vacuum=10.0
) # cuの、(111)面を2層で真空層10.0Åで切り出す

# ---ファイルから
# cif,xyz,vaspからいけるらしい
# cu2o=read('Cu2O.cif')

# ---構造のコピー
# Pythonのルールとして、cu_copy = cu とするだけだと、cu_copyもcuも同じものを指し示してしまう
# (「Python 参照渡し」とかを調べるとわかる)
# copy()を使えば、別のオブジェクトとして扱える


cu_copy = cu.copy()
cu_copy[0].symbol = 'O' # コピー先の原子を変えてみる
print("=== cu_copy = cu.copy()の場合 ===")
print(f"コピー元: {cu[0].symbol}") # コピー元: Cu
print(f"コピー先: {cu_copy[0].symbol}") # コピー先: O
# copy()することで、コピー元とコピー先が別のものを指し示すようになる

cu_ref = cu
cu_ref[0].symbol = 'O'
print("=== cu_ref = cuの場合 ===")
print(f"コピー元: {cu_ref[0].symbol}") # コピー元: O
print(f"コピー先: {cu[0].symbol}") # コピー先: O
# そのまま = で代入すると、コピー元とコピー先が同じものを指し示してしまうので、どちらも置き換わってしまう


# %%
slab

# %% [markdown]
# ### 🟡クラスター

# %%
from ase.cluster import Icosahedron, Octahedron, wulff_construction
from pfcc_extras.visualize.view import view_ngl

# ---八面体
# 正八面体
_cluster=Octahedron(
    "Cu",
    length=7, # 層の数
    cutoff=0 # 頂点をどれくらい切るか。0だと切らない。
)

# 切頂八面体: 頂点が切られている八面体
_cluster=Octahedron(
    "Cu",
    length=7,
    cutoff=2
)

# 正切頂八面体
_cluster=Octahedron(
    "Cu",
    length=7, # length=3*cutoff+1で、正切頂八面体となる
    cutoff=2
)

# 立方八面体
_cluster=Octahedron(
    "Cu",
    length=5, # length=2*cutoff+1で、立方八面体となる
    cutoff=2
)

# ---二十面体
_cluster=Icosahedron(
    "Cu",
    noshells=5, # 原子の数
)

# クラスターの原子数の確認
print(f"原子数: {len(_cluster)}")

# 構造を可視化してみる
view_ngl(_cluster, representations=["ball+stick"], w=400, h=300)

# %% [markdown]
# ### 🟡特定の原子を探す

# %%
from ase import Atoms
from ase.build import bulk

atoms=bulk("Cu",a=3.6,cubic=True)

# %% [markdown]
# #### 🟢インデックスでピンポイントに指定する

# %%
# `atoms[<目的のインデックス>]`で、目的の原子を取得できる
target_index=1
atom = atoms[target_index]
atom

# %% [markdown]
# #### 🟢指定した原子のインデックスを調べる

# %%
from ase import Atoms


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

find_indices_by_symbol(atoms,"Cu")

# %% [markdown]
# #### 🟢隣接原子を探す

# %%
from ase import Atoms, Atom
from ase.neighborlist import NeighborList, natural_cutoffs

def get_neighbors(atoms: Atoms, target_atom:int|Atom, return_type: str = "atoms"):
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

# %% [markdown]
# #### 🟢(平面用)層別に分ける・層ごとのlistにする

# %%
from typing import Literal
import numpy as np
from ase import Atoms, Atom
from ase.constraints import FixAtoms


def separate_layers(
    atoms: Atoms,
    return_type: Literal["atoms", "indices"] = "atoms",
    *,
    decimals: int = 4,
    sort_by_z: bool = True
) -> list[list[Atom]] | list[list[int]]:
    """
    スラブ構造を層別に分離し、各層の原子をリストとして返す。

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

    # --- z座標を丸めて一意な層を特定 ---
    z_coords = atoms.positions[:, 2]
    rounded_z = np.round(z_coords, decimals=decimals)
    unique_z_values = np.unique(rounded_z)
    
    # --- 層をz座標でソート（昇順：下層から上層） ---
    if sort_by_z:
        unique_z_values.sort()
    
    # --- 各層に属する原子インデックスを収集 ---
    layers_indices: list[list[int]] = []
    
    for z_value in unique_z_values:
        # 該当するz座標を持つ原子のインデックスを取得
        layer_mask = np.isclose(rounded_z, z_value, atol=10**(-decimals-1))
        layer_indices = np.where(layer_mask)[0].tolist()
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

# %% [markdown]
# #### 🟢(クラスター用)表面・内側を探す

# %%
from typing import Literal
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs

def classify_surface_atoms(
    atoms: Atoms, 
    return_type: Literal["atoms", "indices"] = "atoms",
    upper_tolerance: int = 3,
) -> tuple[list[Atoms], list[Atoms]] | tuple[list[int], list[int]]:
    """
    クラスターの原子を外表面原子と内側原子に分類する。

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
    surface_indices = [i for i, cnum in enumerate(coordination_numbers) 
                      if surface_min <= cnum <= surface_max]
    inner_indices = [i for i, cnum in enumerate(coordination_numbers) 
                    if cnum > surface_max]

    # --- 出力形式に応じて返す ---
    if return_type == "atoms":
        surface_atoms = [atoms[i] for i in surface_indices]
        inner_atoms = [atoms[i] for i in inner_indices]
        return surface_atoms, inner_atoms
    else:
        return surface_indices, inner_indices

# %% [markdown]
# #### 🟢重心に最も近い原子を探す

# %%
def find_central_atom(
    atoms: Atoms | list[Atom],
    return_type: Literal["atom", "index"] = "atom"
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

# %% [markdown]
# #### 🟢くっつけた後の構造の中で、くっつけた原子のインデックスを知る

# %%
from typing import List
from ase import Atoms, Atom

def get_appended_atom_indices(
    before_atoms: Atoms, 
    after_atoms: Atoms
) -> list[int]:
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
    append_n = n_after - n_before

    # --- 追加された原子のインデックスリストを作成 ---
    appended_indices = list(range(n_before, n_after))

    return appended_indices

# %% [markdown]
# ### 🟡構造を作るために計算する

# %% [markdown]
# #### 🟢配位数を計算する

# %%
from typing import Sequence, Optional
from ase import Atoms, Atom
from ase.neighborlist import NeighborList, natural_cutoffs

def coordination_number(
    atoms: Atoms,
    target_atom: int|Atom,
    return_type: Literal["atoms", "indices"] = "atoms",
    *,
    cutoffs: Optional[Sequence[float]] = None,
    cutoff_scaling: float = 1.0,
) -> tuple[int, list[Atom]|list[int]]:
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
            raise ValueError("指定された Atom オブジェクトは `atoms` 内に存在しません。")
    else:
        raise TypeError("`target_atom` は int（インデックス）または ase.Atom を指定してください。")

    # --- インデックス範囲チェック ---
    if idx < 0 or idx >= len(atoms):
        raise IndexError(f"インデックス {idx} は範囲外です（0 .. {len(atoms)-1}）。")

    # --- return_type チェック ---
    if return_type not in ("atoms", "indices"):
        raise ValueError("`return_type` は 'atoms' または 'indices' を指定してください。")

    # --- カットオフ半径の用意 ---
    if cutoffs is None:
        # natural_cutoffs を用いて原子種に基づくデフォルトのカットオフを得る
        base_cutoffs = natural_cutoffs(atoms)
        # スケーリングを反映
        cutoffs_used = [c * cutoff_scaling for c in base_cutoffs]
    else:
        # ユーザー指定の cutoffs を使用（長さチェック）
        if len(cutoffs) != len(atoms):
            raise ValueError("`cutoffs` の長さは atoms の長さと一致する必要があります。")
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

# %% [markdown]
# #### 🟢法線ベクトルを求める

# %%
from typing import Literal, overload
from numpy.typing import NDArray
import numpy as np
from ase import Atoms, Atom


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
        raise IndexError(f"指定されたインデックス {index} は範囲外です（0-{len(atoms)-1}）。")
    
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
    # 最小固有値が非常に小さい場合は共線と判定
    min_eigval = eigvals[min_eigval_idx]
    tolerance = 1e-10
    if min_eigval < tolerance:
        raise ValueError(
            "指定された点群が共線状態で、平面を一意に定義できません。"
            f"最小固有値: {min_eigval:.2e} < 許容値: {tolerance:.2e}"
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

# %% [markdown]
# ### 🟡処理する

# %% [markdown]
# #### 🟢手動で微妙に動かす

# %%
"""
指定した原子群を微妙に動かすヘルパー関数

ASEによる化学シミュレーションで、原子を指定方向に微小変位させるためのヘルパー関数を提供します。
複数の入力形式（Atom、int、Atoms、list）に対応し、型安全で使いやすいインターフェースを提供します。
"""

from typing import Sequence
import numpy as np
from numpy.typing import NDArray
from ase import Atoms, Atom


def move_atoms(
    base_structure: Atoms,
    target: int | Atom | Atoms | list[int] | list[Atom],
    direction: tuple[float, float, float] | Sequence[float] | NDArray[np.floating],
    distance: float,
    *,
    inplace: bool = False,
) -> Atoms:
    """
    指定した原子を指定方向に指定距離だけ移動する。

    複数の入力形式に対応し、原子を任意の方向に微小変位させることができます。
    方向ベクトルは自動的に正規化され、指定距離だけ移動されます。

    Args:
        base_structure (ase.Atoms): 操作対象の原子構造。
        target (int | Atom | Atoms | list[int] | list[Atom]): 移動させる原子の指定。以下の形式に対応：
            - int: 単一原子のインデックス
            - ase.Atom: 単一原子オブジェクト
            - ase.Atoms: 構造全体（全原子を移動）
            - list[int]: 複数原子のインデックスリスト
            - list[ase.Atom]: 複数原子オブジェクトのリスト
        direction (tuple[float, float, float] | Sequence[float] | NDArray[np.floating]): 移動方向を示すベクトル。
            tuple[float, float, float]、Sequence[float]、またはnumpy配列を指定可能。
            自動的に正規化されます。
        distance (float): 移動させる距離（Å）。正の値で指定方向、負の値で逆方向。
        inplace (bool, optional): Trueの場合、base_structureを直接変更。
            Falseの場合、コピーを作成して返します。デフォルトはFalse。

    Returns:
        ase.Atoms: 原子が移動された構造。inplace=Trueの場合は引数と同じオブジェクト。

    Raises:
        ValueError: 方向ベクトルがゼロベクトルの場合、または指定されたAtomが
                   base_structure内に存在しない場合。
        TypeError: target、directionの型が不正な場合。
        IndexError: 指定されたインデックスが範囲外の場合。

    Examples:
        >>> from ase.build import molecule
        >>> h2o = molecule('H2O')
        
        >>> # 単一原子をx方向に0.1Å移動
        >>> moved = move_atoms(h2o, 0, (1, 0, 0), 0.1)
        
        >>> # 複数原子をz方向に移動
        >>> moved = move_atoms(h2o, [0, 1], (0, 0, 1), 0.2)
        
        >>> # 構造全体をxy平面対角線方向に移動
        >>> moved = move_atoms(h2o, h2o, (1, 1, 0), 0.05)
        
        >>> # インプレース操作
        >>> move_atoms(h2o, [1, 2], (-1, 0, 0), 0.1, inplace=True)
    """
    # --- 入力検証：direction ---
    try:
        direction_vec = np.array(direction, dtype=float)
    except (ValueError, TypeError) as e:
        raise TypeError(f"directionは数値のsequenceまたはnumpy配列を指定してください: {e}")
    
    if direction_vec.shape != (3,):
        raise ValueError(f"directionは3次元ベクトルを指定してください。現在の形状: {direction_vec.shape}")
    
    # --- ゼロベクトルチェック ---
    norm = np.linalg.norm(direction_vec)
    if norm == 0:
        raise ValueError("移動方向ベクトルがゼロです。正しい方向を指定してください。")
    
    # --- 単位ベクトルと変位ベクトルを計算 ---
    unit_vector = direction_vec / norm
    displacement = unit_vector * distance
    
    # --- 作業用構造の準備 ---
    if inplace:
        structure = base_structure
    else:
        structure = base_structure.copy()
    
    # --- target から移動対象原子のインデックスリストを取得 ---
    target_indices: list[int] = []
    
    # --- int: 単一原子インデックス ---
    if isinstance(target, int):
        if not (0 <= target < len(base_structure)):
            raise IndexError(f"インデックス{target}は範囲外です（0-{len(base_structure)-1}）。")
        target_indices = [target]
    
    # --- Atom: 単一原子オブジェクト ---
    elif isinstance(target, Atom):
        try:
            index = target.index
            if not (0 <= index < len(base_structure)):
                raise IndexError(f"インデックス{index}は範囲外です（0-{len(base_structure)-1}）。")
            target_indices = [index]
        except (AttributeError, ValueError):
            raise ValueError("指定されたAtomはbase_structure内に存在しません。")
    
    # --- Atoms: 構造全体 ---
    elif isinstance(target, Atoms):
        target_indices = list(range(len(base_structure)))
    
    # --- list: 複数指定 ---
    elif isinstance(target, list):
        if not target:  # 空リスト
            target_indices = []
        elif all(isinstance(item, int) for item in target):
            # list[int]の場合
            for idx in target:
                if not (0 <= idx < len(base_structure)):
                    raise IndexError(f"インデックス{idx}は範囲外です（0-{len(base_structure)-1}）。")
            target_indices = list(target)
        elif all(isinstance(item, Atom) for item in target):
            # list[Atom]の場合
            for atom in target:
                try:
                    index = atom.index
                    if not (0 <= index < len(base_structure)):
                        raise IndexError(f"インデックス{index}は範囲外です（0-{len(base_structure)-1}）。")
                    target_indices.append(index)
                except (AttributeError, ValueError):
                    raise ValueError("指定されたAtomはbase_structure内に存在しません。")
        else:
            raise TypeError("リストの要素は全てintまたは全てAtomsである必要があります。")
    
    else:
        raise TypeError(
            "targetはint、Atom、Atoms、list[int]、またはlist[Atom]を指定してください。"
        )
    
    # --- 原子の移動を実行 ---
    for idx in target_indices:
        structure[idx].position += displacement
    
    return structure

# %% [markdown]
# #### 🟢(平面用)層を固定する

# %%
from typing import Literal
import numpy as np
from ase import Atoms, Atom
from ase.constraints import FixAtoms


def fix_layers(
    atoms: Atoms,
    fixed_layers: int,
    *,
    inplace: bool = False,
    decimals: int = 4
) -> Atoms:
    """
    指定された数の下層を固定する制約を追加する。

    z座標に基づいて層を特定し、下から指定された数の層に含まれる
    全ての原子に FixAtoms 制約を適用します。

    Args:
        atoms (ase.Atoms): 制約を適用する原子構造。
        fixed_layers (int): 固定する層数（下から数えて）。
            0 の場合は何も固定しません。
        inplace (bool, optional): True の場合は atoms を破壊的に更新。
            False の場合はコピーを作成して返します。デフォルトはFalse。
        decimals (int, optional): z座標の丸め精度（小数点以下の桁数）。
            デフォルトは4。層の判定精度に影響します。

    Returns:
        ase.Atoms: 制約が適用された原子構造。
            inplace=True の場合は引数 atoms 自身が返されます。

    Raises:
        ValueError: fixed_layers が負の値の場合。

    Examples:
        >>> from ase.build import surface, bulk
        >>> slab = surface(bulk('Cu'), (1,1,1), layers=4)
        
        >>> # 下から2層を固定
        >>> fixed_slab = fix_layers(slab, fixed_layers=2)
        >>> print(f"制約数: {len(fixed_slab.constraints)}")
        
        >>> # 元の構造を変更して固定
        >>> fix_layers(slab, fixed_layers=2, inplace=True)

    Note:
        この関数は separate_layers() を内部で使用します。
        固定される原子は z 座標の低い順から fixed_layers 層分です。
    """
    # --- 引数の検証 ---
    if fixed_layers < 0:
        raise ValueError("fixed_layers は0以上の値を指定してください。")

    # --- 作業用原子構造の準備 ---
    if inplace:
        result_atoms = atoms
    else:
        result_atoms = atoms.copy()

    # --- 固定層数が0の場合は何もしない ---
    if fixed_layers == 0:
        return result_atoms

    # --- 層別に分離して固定対象を特定 ---
    layers_indices = separate_layers(
        atoms, 
        return_type="indices", 
        decimals=decimals, 
        sort_by_z=True  # 下層から上層の順序
    )

    total_layers = len(layers_indices)
    
    # --- 固定層数が総層数以上の場合の警告 ---
    if fixed_layers >= total_layers:
        print(f"警告: 固定層数 ({fixed_layers}) が総層数 ({total_layers}) 以上です。")
        print("全ての原子が固定されます。")
        fixed_layers = total_layers

    # --- 固定対象原子のインデックスを収集 ---
    fixed_indices: list[int] = []
    for layer_idx in range(fixed_layers):
        fixed_indices.extend(layers_indices[layer_idx])

    # --- FixAtoms制約を作成・適用 ---
    if fixed_indices:
        # マスク配列を作成（True = 固定する原子）
        mask = np.zeros(len(result_atoms), dtype=bool)
        mask[fixed_indices] = True
        
        constraint = FixAtoms(mask=mask)
        result_atoms.set_constraint(constraint)
        
        # 最上固定層のz座標を取得（報告用）
        fixed_z_coords = result_atoms.positions[fixed_indices, 2]
        max_fixed_z = np.max(np.round(fixed_z_coords, decimals=decimals))
        
        print(f"Z≤{max_fixed_z:.{decimals}f} Å の原子 {len(fixed_indices)} 個を固定しました。")
    else:
        print("固定対象の原子が見つかりませんでした。")

    return result_atoms

# %% [markdown]
# #### 🟢置き換える

# %%
from typing import Mapping, Optional
import numpy as np

def substitute_elements(
    atoms: Atoms,
    target: int|Atom|Atoms|list[int]|list[Atom],
    new: str|Mapping[str, float],
    *,
    inplace: bool = False,
    seed: Optional[int] = None,
) -> Atoms:
    """
    指定した原子（または全原子）を指定した元素に置換します。

    置換後の元素は、単一の元素記号（例: "Pd"）または組成辞書
    （例: {"Cu": 0.8, "Pd": 0.2}）で指定できます。

    Args:
        atoms (ase.Atoms): 操作対象の構造。
        target (Union[int, ase.Atom, ase.Atoms, list[int], list[ase.Atom]]): 置換対象の指定。
            - int: 単一原子のインデックス。
            - list[int]: 複数原子のインデックス。
            - ase.Atom: 単一原子（atoms内に存在する必要あり）。
            - list[ase.Atom]: 複数原子（atoms内に存在する必要あり）。
            - ase.Atoms: 構造全体を対象（全原子を置換）。
        new (Union[str, Mapping[str, float]]): 置換後の指定。
            - str: 単一元素記号 (例: "Cu")。
            - Mapping[str, float]: 組成辞書 (例: {"Cu": 0.9, "Pd": 0.1})。合計は1。
        inplace (bool, optional): True の場合は atoms を破壊的に更新。
            False の場合はコピーを作成して返す。デフォルトは False。
        seed (Optional[int], optional): シャッフルの再現性確保のための乱数シード。

    Returns:
        ase.Atoms: 要求どおりに置換された Atoms オブジェクト。
                   inplace=True の場合は引数 atoms 自身が返ります。

    """
    # ----------
    # --- 置換対象インデックスの決定
    # ----------
    # --- target の型に応じてインデックスを取得 ---
    atom_indices: list[int] = []
    
    if isinstance(target, int):
        index = target
        # 範囲チェック
        if not (0 <= index < len(atoms)):
            raise ValueError(f"インデックスが範囲外です: {index}")
        atom_indices = [index]

    elif isinstance(target, list):
        # リストの場合、要素の型をチェック
        if not target:  # 空リストの場合
            atom_indices = []
        elif all(isinstance(item, int) for item in target):
            # list[int]の場合
            for index in target:
                if not (0 <= index < len(atoms)):
                    raise ValueError(f"インデックスが範囲外です: {index}")
            atom_indices = list(target)
        elif all(isinstance(item, Atom) for item in target):
            # list[Atom]の場合
            atom_indices = []
            for atom in target:
                try:
                    index = atom.index
                    if not (0 <= index < len(atoms)):
                        raise ValueError(f"インデックスが範囲外です: {index}")
                    atom_indices.append(index)
                except (AttributeError, ValueError):
                    raise ValueError("指定されたAtomがatoms内に存在しないか、不正です。")
        else:
            raise TypeError("リストの要素は全てintまたは全てase.Atomである必要があります。")

    elif isinstance(target, Atom):
        try:
            # Atomオブジェクトがatoms内に存在する場合、そのインデックスを取得
            index = target.index
        except (AttributeError, ValueError):
            raise ValueError("指定されたAtomはatoms内に存在しません。")
        # 念のため範囲チェック
        if not (0 <= index < len(atoms)):
            raise ValueError(f"インデックスが範囲外です: {index}")
        atom_indices = [index]

    elif isinstance(target, Atoms):
        # Atomsが与えられた場合は「全原子」を対象とする
        atom_indices = list(range(len(atoms)))

    else:
        raise TypeError("target は int、list[int]、ase.Atom、list[ase.Atom]、または ase.Atoms を指定してください。")

    n_targets = len(atom_indices)
    if n_targets == 0:
        # 置換対象なしの場合はそのまま返す
        return atoms if inplace else atoms.copy()

    # ----------
    # --- 置換後のシンボル列を生成
    # ----------
    replacement_symbols: list[str] = []
    # --- 単一シンボル指定の場合 ---
    if isinstance(new, str):
        replacement_symbols = [new] * n_targets

    # --- 組成辞書指定の場合 ---
    elif isinstance(new, Mapping):
        if n_targets < 0:
            raise ValueError("count は0以上である必要があります。")
        if n_targets == 0:
            replacement_symbols = []
        else:
            # 合計の検証（許容誤差内で1）
            total = float(sum(new.values()))
            if not np.isclose(total, 1.0):
                raise ValueError(f"組成の合計が1ではありません（現在: {total}）。")

            # 元素記号の検証
            symbols = list(new.keys())

            # --- 丸めベースで個数を算出し、最後の元素に差分を集約 ---
            counts = {s: int(round(n_targets * frac)) for s, frac in new.items()}

            # 丸め誤差の補正
            diff = n_targets - sum(counts.values())
            counts[symbols[-1]] += diff

            # --- シンボル列を生成し、シャッフル ---
            replacement_symbols = []
            for s, c in counts.items():
                replacement_symbols.extend([s] * c)

            # 万一長さが合わない場合は例外
            if len(replacement_symbols) != n_targets:
                raise RuntimeError(
                    f"置換シンボル列の長さが一致しません: {len(replacement_symbols)} != {n_targets}"
                )

            rng = np.random.default_rng(seed)
            rng.shuffle(replacement_symbols)
    else:
        raise TypeError("new は str または Mapping[str, float] を指定してください。")

    # ----------
    # --- 書き換え先の Atoms を準備
    # ----------
    if inplace:
        out = atoms
    else:
        out = atoms.copy()

    # --- インデックス順にシンボルを割り当て ---
    # 対象インデックスの順序に対応するようにする
    for idx, sym in zip(atom_indices, replacement_symbols):
        out[idx].symbol = sym

    return out

# %% [markdown]
# #### 平面にくっつける

# %%


# %% [markdown]
# #### 🟢クラスターにくっつける

# %%
from typing import Optional, Sequence

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

# %%
def place_adsorbate_along_normal(
    substrate: Atoms,
    adsorbate: Atoms,
    target_atom: int | Atom,
    distance: float,
    *,
    upper_tolerance: int = 1,
    lower_tolerance: int = 1,
) -> Atoms:
    """
    指定した基板の表面法線方向に、吸着分子を配置する。
    
    指定した基板の原子を基準として、その周辺の局所平面を主成分分析（PCA）で求め、
    法線方向に指定距離だけ離れた位置に吸着分子を配置する。吸着分子の+z軸が
    法線ベクトルに一致するように回転調整も行う。
    
    Args:
        substrate (ase.Atoms): 基板となる原子構造。
        adsorbate (ase.Atoms): 配置する吸着分子。
        target_atom (int | ase.Atom): 基準となる基板上の原子。
            インデックスまたはAtomオブジェクトを指定可能。
        distance (float): 法線方向に離す距離（Å）。正の値を指定。
        upper_tolerance (int): 最小配位数から上方向への許容範囲。デフォルトは1。
        lower_tolerance (int): 最小配位数から下方向への許容範囲。デフォルトは1。
    
    Returns:
        ase.Atoms: 基板と配置済み吸着分子を結合した構造。
    
    Raises:
        ValueError: distance が負の値の場合、または法線計算に必要な点数が不足している場合。
        TypeError: target_atom の型が int でも ase.Atom でもない場合。
        IndexError: 指定されたインデックスが範囲外の場合。
    
    Note:
        - 法線ベクトルの向きは「重心→target_atom」方向を参照として決定されるため、
          常に基板表面から外側を向く。
        - 吸着分子の内部+z軸が法線に整列するように回転される。
        - この関数は compute_surface_normal と get_neighbors を利用する。
        
    Example:
        >>> from ase.build import fcc111, molecule
        >>> slab = fcc111('Cu', size=(3, 3, 4))
        >>> co = molecule('CO')
        >>> combined = place_adsorbate_along_normal(slab, co, 0, 2.0)
        >>> len(combined) == len(slab) + len(co)
        True
    """
    # --- 入力検証 ---
    if distance < 0:
        raise ValueError(f"distance は正の値を指定してください。指定値: {distance}")
    
    # --- target_atom の型に応じてインデックスを取得 ---
    if isinstance(target_atom, int):
        target_index = target_atom
    elif isinstance(target_atom, Atom):
        try:
            target_index = target_atom.index
        except ValueError:
            raise ValueError("指定されたAtomはsubstrate内に存在しません。")
    else:
        raise TypeError("target_atom は int または ase.Atom を指定してください。")
    
    # --- インデックスの範囲チェック ---
    if not (0 <= target_index < len(substrate)):
        raise IndexError(f"指定されたインデックス {target_index} は範囲外です（0-{len(substrate)-1}）。")
    
    # --- 局所点群の重心を計算 ---
    neighbor_indices = get_neighbors_with_coordination_condition(substrate, target_index, return_type="indices", upper_tolerance=upper_tolerance, lower_tolerance=lower_tolerance)
    
    point_indices = list(neighbor_indices) + [target_index]  # include_target=True相当
    
    if len(point_indices) < 3:
        raise ValueError(
            f"法線計算に必要な点数が不足しています。"
            f"必要: 3点以上、取得: {len(point_indices)}点"
        )
    
    points = np.array([substrate[i].position for i in point_indices])
    centroid = np.mean(points, axis=0)
    
    # --- 参照ベクトルの計算（重心→target_atom方向） ---
    target_pos = substrate[target_index].position
    
    substrate_centroid = substrate.get_center_of_mass()
    reference_vector = target_pos - substrate_centroid
    
    # --- 法線ベクトルを計算 ---
    normal = compute_surface_normal(
        substrate, 
        target_index, 
        include_target=True,
        reference_vector=reference_vector,
        normalize=True,
        return_plane=False
    )
    
    # --- 吸着分子のコピーを作成 ---
    adsorbate_copy = adsorbate.copy()
    
    # --- 吸着分子の重心を計算 ---
    com = adsorbate_copy.get_center_of_mass()
    
    # --- z軸を法線に整列する回転行列を計算 ---
    z_axis = np.array([0.0, 0.0, 1.0])
    rotation_matrix = _compute_rotation_matrix(z_axis, normal)
    
    # --- 重心を中心とした回転を適用 ---
    for atom in adsorbate_copy:
        # 重心からの相対位置を計算
        rel_pos = atom.position - com
        # 回転を適用
        rotated_pos = rotation_matrix @ rel_pos
        # 重心に戻す
        atom.position = rotated_pos + com
    
    # --- 回転後の重心を再計算 ---
    com_rotated = adsorbate_copy.get_center_of_mass()
    
    # --- 目標位置への平行移動 ---
    target_position = target_pos + normal * distance
    translation = target_position - com_rotated
    
    for atom in adsorbate_copy:
        atom.position += translation
    
    # --- 基板と吸着分子を結合 ---
    combined = substrate + adsorbate_copy
    
    return combined


def _compute_rotation_matrix(v_from: NDArray[np.float64], v_to: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    ベクトル v_from を v_to に回転させる回転行列を Rodrigues の回転公式で計算する。
    
    Args:
        v_from (NDArray[np.float64]): 回転前のベクトル（正規化済みを想定）。
        v_to (NDArray[np.float64]): 回転後のベクトル（正規化済みを想定）。
    
    Returns:
        NDArray[np.float64]: 3x3 回転行列。
    
    Note:
        v_from と v_to が同方向または反対方向の場合は特別に処理される。
    """
    # --- ベクトルの正規化 ---
    v_from = v_from / np.linalg.norm(v_from)
    v_to = v_to / np.linalg.norm(v_to)
    
    # --- 内積で角度関係をチェック ---
    dot_product = np.dot(v_from, v_to)
    
    # --- 同方向の場合（ほぼ平行） ---
    if np.isclose(dot_product, 1.0, atol=1e-10):
        return np.eye(3)
    
    # --- 反対方向の場合（ほぼ反平行） ---
    if np.isclose(dot_product, -1.0, atol=1e-10):
        # 任意の垂直軸を見つけて180度回転
        # v_from に垂直なベクトルを作成
        if abs(v_from[0]) < 0.9:
            perpendicular = np.array([1.0, 0.0, 0.0])
        else:
            perpendicular = np.array([0.0, 1.0, 0.0])
        
        # グラム・シュミット法で正規直交化
        perpendicular = perpendicular - np.dot(perpendicular, v_from) * v_from
        perpendicular = perpendicular / np.linalg.norm(perpendicular)
        
        # 180度回転行列: R = 2 * n * n^T - I （nは回転軸）
        return 2.0 * np.outer(perpendicular, perpendicular) - np.eye(3)
    
    # --- 一般的な場合：Rodrigues の回転公式 ---
    # 回転軸は外積で求める
    k = np.cross(v_from, v_to)
    k = k / np.linalg.norm(k)
    
    # 回転角度
    theta = np.arccos(np.clip(dot_product, -1.0, 1.0))
    
    # Rodrigues の回転公式: R = I + sin(θ) * [k]× + (1 - cos(θ)) * [k]×^2
    # [k]× は k の外積行列（歪対称行列）
    K = np.array([
        [0, -k[2], k[1]],
        [k[2], 0, -k[0]],
        [-k[1], k[0], 0]
    ])
    
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
    
    return R



# %%
from ase.build import molecule, bulk, surface

# ---結晶
# _bulk=bulk('Cu')
# _bulk=_bulk*(6,6,3)
_cluster=Octahedron(
    "Cu",
    length=7,
    cutoff=2
)

_molecule = molecule("CO")

# compute_surface_normal(_bulk, 40)
TARGET_INDEX=200

# print(coordination_number(_cluster,TARGET_INDEX))
# for a in get_neighbors(_cluster, TARGET_INDEX):
#     print(coordination_number(_cluster,a)[0])

_bulk=substitute_elements(
    _cluster,
    get_neighbors_with_coordination_condition(_cluster, TARGET_INDEX, upper_tolerance=2, lower_tolerance=2),
    "Au"
)

# view_ngl(_bulk)
adsorbed = place_adsorbate_along_normal(_bulk, _molecule, TARGET_INDEX, 2, upper_tolerance=2)
view_ngl(adsorbed)

# %%
view_ngl(substitute_elements(
    _cluster, get_neighbors_with_equal_coordination(_cluster, 295), "Au"
))

# %%
view_ngl(_cluster)

# %% [markdown]
# ## 🔴計算する

# %% [markdown]
# ### 🟡吸着エネルギーを計算する

# %%
from typing import Literal, Optional, Sequence
from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray
from ase import Atoms, Atom
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer
from matlantis_features.ase_ext.optimize import FIRELBFGS
import logging
from datetime import datetime


class ConditionalLogger:
    """
    ログ出力を条件によって制御するラッパークラス。
    
    enabled=False の場合、すべてのログ出力を無効化する。
    try-delta-g.py から移植。
    """
    def __init__(self, base_logger, enabled: bool = True):
        self.base_logger = base_logger
        self.enabled = enabled
    
    def __getattr__(self, name):
        """
        logger.info(), logger.warning() などの呼び出しを透過的に処理する。
        enabled=False の場合は何もしない関数を返す。
        """
        if not self.enabled:
            # ログが無効な場合は何もしない関数を返す
            return lambda *args, **kwargs: None
        
        # ログが有効な場合は元のloggerのメソッドを返す
        return getattr(self.base_logger, name)

def setup_logger():
    """デバッグ用ログの設定"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"calc_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename, encoding='utf-8'),
            logging.StreamHandler()  # コンソールにも出力
        ]
    )
    
    base_logger = logging.getLogger(__name__)
    base_logger.info(f"デバッグログファイル: {log_filename}")
    return ConditionalLogger(base_logger, enabled=True)

@dataclass
class CAEInput:
    # calculate_adsorption_energy()で使用する入力を、構造的に扱うためのクラス。
    structure: Atoms
    calc_mode: Literal['molecule', 'solid']
    

def calculate_adsorption_energy(
    calculator_molecule: Calculator,
    calculator_solid: Calculator,
    adsorbed_structure_input: CAEInput,
    reactant_structures_input: list[CAEInput],
    *,
    optimizer: Optimizer = FIRELBFGS,
    opt_fmax: float = 0.05,
    opt_maxsteps: int = 3000,
    logger: Optional[ConditionalLogger] = None,
    enable_logging: bool = True,
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
        optimizer (Optimizer, optional): 構造最適化アルゴリズム。
            デフォルトは FIRELBFGS。
        opt_fmax (float, optional): 構造最適化の力の収束閾値（eV/Å）。デフォルトは 0.05。
        opt_maxsteps (int, optional): 構造最適化の最大ステップ数。デフォルトは 3000。
        logger (ConditionalLogger | None, optional): ログ出力制御。
            Noneの場合は新規作成。
        enable_logging (bool, optional): ログ出力の有効/無効。デフォルトは True。
    
    Returns:
        float: 吸着エネルギー（eV）。負の値は吸着が熱力学的に有利であることを示す。
    
    Raises:
        ValueError: calc_mode_reactants の長さが reactant_structures と不一致の場合。
        TypeError: 引数の型が不正な場合。
    """
    
    # --- ログ設定 ---
    if logger is None:
        logger = ConditionalLogger(None, enabled=False)
    
    # --- 引数検証 ---
    n_reactants = len(reactant_structures_input)
    if n_reactants == 0:
        raise ValueError("reactant_structures は少なくとも1つの構造を含む必要があります。")

    # --- 内部ヘルパー関数: 構造最適化とエネルギー計算 ---
    def _optimize_and_get_energy(
        calc_input: CAEInput,
        label: str
    ) -> float:
        """
        単一構造の最適化とエネルギー取得。
        
        Args:
            calc_input: 計算する構造、振動させる分子のインデックス、計算モード
            label: ログ用ラベル
            
        Returns:
            最適化後のポテンシャルエネルギー (eV)
        """
        # 作業用コピーを作成
        work_atoms = calc_input.structure.copy()
        
        # 計算機を設定
        if calc_input.calc_mode == 'molecule':
            work_atoms.calc = calculator_molecule
        else:  # 'solid'
            work_atoms.calc = calculator_solid 
        
        logger.info(f"--- {label} 処理開始 ---")
        logger.info(f"原子数: {len(work_atoms)}")
        logger.info(f"組成: {work_atoms.symbols}")
        logger.info(f"計算モード: {calc_input.calc_mode}")
        
        # 制約情報のログ出力
        constraints = work_atoms.constraints
        if constraints:
            logger.info(f"制約条件: {len(constraints)} 個")
            for i, constraint in enumerate(constraints):
                constraint_name = type(constraint).__name__
                logger.info(f"  制約{i}: {constraint_name}")
                # 固定原子数の情報（FixAtomsの場合）
                if hasattr(constraint, 'index'):
                    n_fixed = len(constraint.index) if hasattr(constraint.index, '__len__') else 1
                    logger.info(f"    固定原子数: {n_fixed}")
        else:
            logger.info("制約条件: なし")
        
        # 初期エネルギー取得
        e_initial = work_atoms.get_potential_energy()
        logger.info(f"初期エネルギー: {e_initial:.6f} eV")
        
        # 構造最適化実行
        logger.info("構造最適化開始")
        opt = optimizer(work_atoms)
        opt.run(fmax=opt_fmax, steps=opt_maxsteps)
        logger.info("構造最適化完了")
        
        # 最適化後エネルギー取得
        e_final = work_atoms.get_potential_energy()
        e_change = e_final - e_initial
        logger.info(f"最適化後エネルギー: {e_final:.6f} eV")
        logger.info(f"エネルギー変化: {e_change:.6f} eV")
        logger.info(f"--- {label} 処理完了 ---")
        
        return e_final
    
    # --- メイン計算開始ログ ---
    logger.info("="*80)
    logger.info("吸着エネルギー計算開始")
    logger.info(f"吸着後構造: {adsorbed_structure_input.structure.symbols} ({len(adsorbed_structure_input.structure)} 原子)")
    logger.info(f"反応物構造数: {n_reactants}")
    for i, atoms in enumerate(reactant_structures_input):
        logger.info(f"  反応物{i+1}: {atoms.structure.symbols} ({len(atoms.structure)} 原子)")
    logger.info(f"計算モード設定:")
    logger.info(f"  吸着後: {adsorbed_structure_input.calc_mode}")
    logger.info(f"  反応物: {reactant_structures_input[0].calc_mode}")
    
    # --- 1. 吸着後構造のエネルギー計算 ---
    e_adsorbed = _optimize_and_get_energy(
        adsorbed_structure_input, 
        "吸着後構造"
    )
    
    # --- 2. 各反応物構造のエネルギー計算 ---
    reactant_energies: list[float] = []
    for i, reactant_input in enumerate(reactant_structures_input):
        label = f"反応物{i+1}"
        e_reactant = _optimize_and_get_energy(reactant_input, label)
        reactant_energies.append(e_reactant)
    
    # --- 3. 吸着エネルギー計算 ---
    e_reactants_total = sum(reactant_energies)
    e_adsorption = e_adsorbed - e_reactants_total
    
    # --- 結果ログ出力 ---
    logger.info("="*80)
    logger.info("吸着エネルギー計算結果")
    logger.info(f"吸着後構造エネルギー: {e_adsorbed:.6f} eV")
    logger.info("反応物エネルギー:")
    for i, e in enumerate(reactant_energies):
        logger.info(f"  反応物{i+1} ({reactant_structures_input[i].structure.symbols}): {e:.6f} eV")
    logger.info(f"反応物合計エネルギー: {e_reactants_total:.6f} eV")
    logger.info(f"吸着エネルギー: {e_adsorption:.6f} eV")
    if e_adsorption < 0:
        logger.info("→ 吸着は熱力学的に有利")
    else:
        logger.info("→ 吸着は熱力学的に不利")
    logger.info("="*80)
    
    # コンソール出力用サマリー
    print(f"\n{'='*50}")
    print(f"吸着エネルギー計算完了")
    print(f"E_ads = {e_adsorption:.3f} eV")
    print(f"{'='*50}\n")
    
    return e_adsorption

# %% [markdown]
# ### 🟡生成エネルギーを計算する

# %%
from typing import Optional
from collections import Counter
import numpy as np
from ase import Atoms, Atom
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer
from ase.build import bulk
from matlantis_features.ase_ext.optimize import FIRELBFGS
import logging
from datetime import datetime


class ConditionalLogger:
    """
    ログ出力を条件によって制御するラッパークラス。
    
    enabled=False の場合、すべてのログ出力を無効化する。
    adsorption_energy_helper.py から移植。
    """
    def __init__(self, base_logger, enabled: bool = True):
        self.base_logger = base_logger
        self.enabled = enabled
    
    def __getattr__(self, name):
        """
        logger.info(), logger.warning() などの呼び出しを透過的に処理する。
        enabled=False の場合は何もしない関数を返す。
        """
        if not self.enabled:
            # ログが無効な場合は何もしない関数を返す
            return lambda *args, **kwargs: None
        
        # ログが有効な場合は元のloggerのメソッドを返す
        return getattr(self.base_logger, name)


def setup_formation_logger() -> ConditionalLogger:
    """生成エネルギー計算用ログの設定。
    
    Returns:
        ConditionalLogger: 設定済みのロガー。
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"formation_energy_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename, encoding='utf-8'),
            logging.StreamHandler()  # コンソールにも出力
        ]
    )
    
    base_logger = logging.getLogger(__name__)
    base_logger.info(f"生成エネルギー計算ログファイル: {log_filename}")
    return ConditionalLogger(base_logger, enabled=True)



def analyze_composition(atoms: Atoms) -> dict[str, int]:
    """
    原子構造から元素組成を解析する。
    
    指定された原子構造に含まれる各元素の個数をカウントし、
    辞書形式で返します。
    
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
    crystal_structure: str = 'auto',
    lattice_parameter: Optional[float] = None,
    logger: Optional[ConditionalLogger] = None
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
        lattice_parameter (float | None, optional): 格子定数（Å）。
            Noneの場合は ASE の標準値を使用。
        logger (ConditionalLogger | None, optional): ロガー。
            Noneの場合はログ出力しない。
    
    Returns:
        Atoms: 純元素の参照構造（単原子）。
    
    Raises:
        ValueError: 対応していない元素または結晶構造の場合。
    
    Note:
        自動判別は一般的な金属の結晶構造に基づいています：
        - fcc: Cu, Au, Ag, Al, Ni, Pt, Pd など
        - bcc: Fe, Cr, W, Mo, V など  
        - hcp: Zn, Mg, Ti, Zr など
    """
    # --- ログ設定 ---
    if logger is None:
        logger = ConditionalLogger(None, enabled=False)
    
    # --- 元素記号の正規化 ---
    element = element.capitalize()
    
    # --- 結晶構造の自動判別 ---
    if crystal_structure == 'auto':
        # 一般的な金属の結晶構造データベース
        crystal_db = {
            # fcc構造
            'Cu': 'fcc', 'Au': 'fcc', 'Ag': 'fcc', 'Al': 'fcc', 
            'Ni': 'fcc', 'Pt': 'fcc', 'Pd': 'fcc', 'Rh': 'fcc',
            'Ir': 'fcc', 'Pb': 'fcc', 'Ca': 'fcc', 'Sr': 'fcc',
            
            # bcc構造  
            'Fe': 'bcc', 'Cr': 'bcc', 'W': 'bcc', 'Mo': 'bcc',
            'V': 'bcc', 'Nb': 'bcc', 'Ta': 'bcc', 'Ba': 'bcc',
            
            # hcp構造
            'Zn': 'hcp', 'Mg': 'hcp', 'Ti': 'hcp', 'Zr': 'hcp',
            'Co': 'hcp', 'Cd': 'hcp', 'Be': 'hcp', 'Ru': 'hcp',
            'Os': 'hcp', 'Re': 'hcp'
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
        valid_structures = ['fcc', 'bcc', 'hcp']
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
        
        logger.info(f"参照構造生成完了: {element} ({determined_structure}, {len(ref_structure)} 原子)")
        return ref_structure
        
    except Exception as e:
        raise ValueError(f"元素 '{element}' の参照構造生成に失敗しました: {e}")


def calculate_formation_energy(
    calculator: Calculator,
    compound_structure: Atoms,
    *,
    optimizer: type[Optimizer] = FIRELBFGS,
    opt_fmax: float = 0.05,
    opt_maxsteps: int = 3000,
    reference_crystal_structures: Optional[dict[str, str]] = None,
    reference_lattice_parameters: Optional[dict[str, float]] = None,
    logger: Optional[ConditionalLogger] = None,
    enable_logging: bool = True,
) -> float:
    """
    金属の生成エネルギーを計算する。
    
    指定された化合物構造から元素組成を自動解析し、純元素参照構造との
    エネルギー差から生成エネルギー E_formation = E(化合物) - Σ[n_i × E(純元素_i)] を算出します。
    
    Args:
        calculator (Calculator): 固体用計算機。
            一般的に EstimatorCalcMode.CRYSTAL_U0 を使用。
        compound_structure (Atoms): 化合物の原子構造。
        optimizer (type[Optimizer], optional): 構造最適化アルゴリズム。
            デフォルトは FIRELBFGS。
        opt_fmax (float, optional): 構造最適化の力の収束閾値（eV/Å）。デフォルトは 0.05。
        opt_maxsteps (int, optional): 構造最適化の最大ステップ数。デフォルトは 3000。
        reference_crystal_structures (dict[str, str] | None, optional): 
            純元素の結晶構造を手動指定する辞書。キーは元素記号、値は 'fcc', 'bcc', 'hcp'。
            例: {'Cu': 'fcc', 'Fe': 'bcc'}
        reference_lattice_parameters (dict[str, float] | None, optional):
            純元素の格子定数を手動指定する辞書。キーは元素記号、値は格子定数（Å）。
            例: {'Cu': 3.615, 'Au': 4.078}
        logger (ConditionalLogger | None, optional): ログ出力制御。
            Noneの場合は新規作成。
        enable_logging (bool, optional): ログ出力の有効/無効。デフォルトは True。
    
    Returns:
        float: 生成エネルギー（eV）。負の値は化合物形成が熱力学的に有利であることを示す。
    
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
    if logger is None:
        if enable_logging:
            logger = setup_formation_logger()
        else:
            logger = ConditionalLogger(None, enabled=False)
    
    # --- 引数検証 ---
    if len(compound_structure) == 0:
        raise ValueError("化合物構造が空です。")
    
    # --- 内部ヘルパー関数: 構造最適化とエネルギー計算 ---
    def _optimize_and_get_energy(
        atoms: Atoms,
        label: str
    ) -> float:
        """
        単一構造の最適化とエネルギー取得。
        
        Args:
            atoms: 計算する原子構造
            label: ログ用ラベル
            
        Returns:
            最適化後のポテンシャルエネルギー (eV)
        """
        # 作業用コピーを作成
        work_atoms = atoms.copy()
        
        # 計算機を設定
        work_atoms.calc = calculator
        
        logger.info(f"--- {label} 処理開始 ---")
        logger.info(f"原子数: {len(work_atoms)}")
        logger.info(f"組成: {work_atoms.symbols}")
                
        # 初期エネルギー取得
        e_initial = work_atoms.get_potential_energy()
        logger.info(f"初期エネルギー: {e_initial:.6f} eV")
        
        # 構造最適化実行
        logger.info("構造最適化開始")
        opt = optimizer(work_atoms)
        opt.run(fmax=opt_fmax, steps=opt_maxsteps)
        logger.info("構造最適化完了")
        
        # 最適化後エネルギー取得
        e_final = work_atoms.get_potential_energy()
        e_change = e_final - e_initial
        logger.info(f"最適化後エネルギー: {e_final:.6f} eV")
        logger.info(f"エネルギー変化: {e_change:.6f} eV")
        logger.info(f"--- {label} 処理完了 ---")
        
        return e_final
    
    # --- メイン計算開始 ---
    logger.info("="*80)
    logger.info("金属生成エネルギー計算開始")
    logger.info(f"化合物構造: {compound_structure.symbols} ({len(compound_structure)} 原子)")
    
    # --- 1. 化合物の元素組成解析 ---
    composition = analyze_composition(compound_structure)
    logger.info("元素組成解析結果:")
    for element, count in composition.items():
        logger.info(f"  {element}: {count} 原子")
    
    # --- 2. 化合物のエネルギー計算 ---
    e_compound = _optimize_and_get_energy(
        compound_structure, 
        "化合物構造"
    )
    
    # --- 3. 各純元素のエネルギー計算 ---
    element_energies: dict[str, float] = {}
    
    for element in composition.keys():
        logger.info(f"\n純元素 {element} の参照構造準備")
        
        # 結晶構造の決定
        crystal_structure = 'auto'
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
                logger=logger
            )
            
            # エネルギー計算実行
            e_element = _optimize_and_get_energy(
                ref_structure,
                f"純元素 {element}"
            )
            
            element_energies[element] = e_element
            
        except Exception as e:
            logger.error(f"純元素 {element} の計算に失敗: {e}")
            raise ValueError(f"純元素 {element} の参照エネルギー計算に失敗しました: {e}")
    
    # --- 4. 生成エネルギー計算 ---
    # 純元素エネルギーの加重和を計算
    e_elements_total = sum(
        count * element_energies[element] 
        for element, count in composition.items()
    )
    
    # 生成エネルギー = 化合物エネルギー - 純元素エネルギー合計
    e_formation = e_compound - e_elements_total
    
    # --- 結果ログ出力 ---
    logger.info("="*80)
    logger.info("生成エネルギー計算結果")
    logger.info(f"化合物エネルギー: {e_compound:.6f} eV")
    logger.info("純元素エネルギー:")
    for element, energy in element_energies.items():
        count = composition[element]
        total = count * energy
        logger.info(f"  {element}: {energy:.6f} eV × {count} = {total:.6f} eV")
    logger.info(f"純元素合計エネルギー: {e_elements_total:.6f} eV")
    logger.info(f"生成エネルギー: {e_formation:.6f} eV")
    if e_formation < 0:
        logger.info("→ 化合物形成は熱力学的に有利")
    else:
        logger.info("→ 化合物形成は熱力学的に不利")
    logger.info("="*80)
    
    # コンソール出力用サマリー
    print(f"\n{'='*50}")
    print(f"生成エネルギー計算完了")
    print(f"組成: {composition}")
    print(f"E_formation = {e_formation:.3f} eV")
    print(f"{'='*50}\n")
    
    return e_formation

# %%
from ase.build import bulk
from ase.constraints import FixAtoms

# メイン実行部
print("="*60)
print("Cu3Au合金の生成エネルギー計算サンプル")
print("="*60)

# 1. 計算機の設定
print("計算機を設定中...")

import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

# バージョンによって計算結果が異なる場合があるため、毎回バージョンを確認する
print(f"pfp_api_client: {pfp_api_client.__version__}") 

# ASECalculatorを使用するための設定
estimator = Estimator(calc_mode=EstimatorCalcMode.CRYSTAL_U0)
calc_solid = ASECalculator(estimator) # 固体用計算機

estimator = Estimator(calc_mode=EstimatorCalcMode.MOLECULE)
calc_molecule = ASECalculator(estimator) # 分子用計算機

print("✓ Matlantis計算機を設定しました")

# 2. 構造の作成
print("\n構造を作成中...")

# Cu3Au合金構造の仮想的な作成（実際の構造は別途準備が必要）
# ここでは簡単のため、CuのfccベースにAuを置換した構造を作成
cu_base = bulk("Cu", "fcc", a=3.615)
cu3au_structure = cu_base.repeat((2, 2, 1))  # 2x2x1の拡張（8原子）

# 手動でCu3Au組成に調整（8原子中2個をAuに変更）
cu3au_structure[0].symbol = 'Au'
cu3au_structure[2].symbol = 'Au'

print(f"✓ Cu3Au構造: {len(cu3au_structure)} 原子")
print(f"  組成: {dict(Counter(atom.symbol for atom in cu3au_structure))}")

print("✓ 化合物構造を準備しました")

# 4. 生成エネルギー計算実行
print("\n生成エネルギー計算を開始...")
print("注意: 実際の計算には時間がかかります（数分〜数十分）")
print("-" * 60)

formation_energy = calculate_formation_energy(
    calculator=calc_solid,
    compound_structure=cu3au_structure,
    opt_fmax=0.05,
    opt_maxsteps=3000,
    logger=setup_formation_logger(),
    enable_logging=True
)

# 5. 結果の表示
print("\n" + "="*60)
print("最終結果")
print("="*60)
print(f"Cu3Au合金の生成エネルギー: {formation_energy:.4f} eV")

if formation_energy < 0:
    print("→ 合金形成は熱力学的に有利です（発熱的）")
else:
    print("→ 合金形成は熱力学的に不利です（吸熱的）")

print(f"\n参考:")
print(f"  負の値: 安定な化合物・合金")
print(f"  正の値: 不安定な化合物・合金")
print(f"  計算値: {formation_energy:.4f} eV")
print("="*60)


# %% [markdown]
# ### 🟡NEBする

# %%
from collections.abc import Sequence, Callable
from typing import Any
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ase.atoms import Atoms
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer
from ase.neb import NEB
from ase.build.rotate import minimize_rotation_and_translation

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode


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

    Note:
        - parallel=Trueを使用する場合は、MPI環境での実行が必要です。
        - parallel=Trueの場合、各画像に個別の計算器を作成します（allow_shared_calculator=False）。
        - parallel=Falseの場合、全画像で計算器を共有します（allow_shared_calculator=True）。

    Examples:
        >>> from ase.optimize import FIRE
        >>> from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
        >>> from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
        >>> 
        >>> # Estimatorの準備
        >>> estimator = Estimator(calc_mode=EstimatorCalcMode.PBE, model_version="v8.0.0")
        >>> 
        >>> # NEB計算の実行
        >>> images, energies = run_neb(
        ...     init_atoms, final_atoms, 5, FIRE, estimator,
        ...     fmax=0.1, steps=100, trajectory_path='neb.traj'
        ... )
    """
    # 入力検証
    if len(init_atoms) != len(final_atoms):
        raise ValueError(
            f"初期構造と最終構造の原子数が一致しません: "
            f"{len(init_atoms)} != {len(final_atoms)}"
        )
    
    if num_intermediate_images < 0:
        raise ValueError(f"中間構造数は0以上である必要があります: {num_intermediate_images}")

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
        parallel=parallel
    )

    # 構造の補間
    interpolate_args = interpolate_kwargs or {}
    if mic is not None:
        interpolate_args['mic'] = mic
    neb.interpolate(**interpolate_args)

    # 最適化の実行
    if trajectory_path is not None:
        optimizer = optimizer_cls(neb, trajectory=trajectory_path)
    else:
        optimizer = optimizer_cls(neb)
    
    optimizer.run(fmax=fmax, steps=steps)

    # エネルギーの計算
    energies = [image.get_total_energy() for image in images]

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
    ax.plot(replica_indices, energies, 'o-', linewidth=2, markersize=6)
    
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


# %% [markdown]
# ### 🟡ギブス自由エネルギーを計算する

# %%
"""
ギブス自由エネルギーを計算するヘルパー関数


基本的な部分: 構造・振動させる分子のインデックス・どっちで計算するか(IdealGasThermo/HarmonicThermo)→(opt)→最適構造→(振動計算・g計算)→自由エネルギー
CHEモデル対応: 電極電位・pH・温度→(CHEモデル)→(H+ + e-)の自由エネルギー
最終的な計算: 左辺右辺で自由エネルギーをまとめる

"""

from ase.atoms import Atoms
from ase.build import molecule
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer
from matlantis_features.ase_ext.optimize import FIRELBFGS
from typing import Literal, Optional
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
import numpy as np
import logging
from datetime import datetime
from dataclasses import dataclass

import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

# バージョンによって計算結果が異なる場合があるため、毎回バージョンを確認する
print(f"pfp_api_client: {pfp_api_client.__version__}") 

# ASECalculatorを使用するための設定
# EstimatorCalcModeは、以下のように使い分ける
# - 一般の系： EstimatorCalcMode.CRYSTAL_U0 Uなしモード
# - 酸化物など： EstimatorCalcMode.CRYSTAL　Uありモード
# - 単体有機分子： EstimatorCalcMode.MOLECULE 分子モード
estimator_mol = Estimator(calc_mode=EstimatorCalcMode.MOLECULE)
calc_molecule = ASECalculator(estimator_mol)
estimator_solid = Estimator(calc_mode=EstimatorCalcMode.CRYSTAL_U0)
calc_solid = ASECalculator(estimator_solid)

class ConditionalLogger:
    """
    ログ出力を条件によって制御するラッパークラス。
    enabled=False の場合、すべてのログ出力を無効化する。
    """
    def __init__(self, base_logger, enabled: bool = True):
        self.base_logger = base_logger
        self.enabled = enabled
    
    def __getattr__(self, name):
        """
        logger.info(), logger.warning() などの呼び出しを透過的に処理する。
        enabled=False の場合は何もしない関数を返す。
        """
        if not self.enabled:
            # ログが無効な場合は何もしない関数を返す
            return lambda *args, **kwargs: None
        
        # ログが有効な場合は元のloggerのメソッドを返す
        return getattr(self.base_logger, name)


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
    calc_mode: Literal['IdealGasThermo', 'HarmonicThermo']
    vibrate_indices: list[int] | None = None
    # IdealGasThermo用の引数
    geometry: Literal['linear', 'nonlinear'] = 'nonlinear'
    symmetry_number: int = 1
    spin_multiplicity: int = 1
    # opt用の引数
    do_opt: bool = True


def calculate_g(
    calculator_molecule: Calculator,
    calculator_solid: Calculator,
    calc_input: CGFEInput,
    *,
    temperature: float = 298.15,
    pressure: float = 101325.0,
    optimizer: Optimizer = None,
    opt_fmax: float = 0.05,
    opt_maxsteps: int = 3000,
    logger: Optional['ConditionalLogger'] = None,
):
    """
    ギブス自由エネルギーを計算する。

    Args:
        calculator_molecule (ase.calculators.calculator.Calculator): 分子用計算機
        calculator_solid (ase.calculators.calculator.Calculator): 固体用計算機
        calc_input (CGFEInput): 計算する構造、振動させる分子のインデックス、計算モード
        temperature (float): 温度（K）
        pressure (float): 圧力（Pa）
        optimizer (ase.optimize.optimize.Optimizer): 最適化エンジン
        opt_fmax (float): 最適化の閾値
        opt_maxsteps (int): 最適化の最大ステップ数
        logger (ConditionalLogger): ロガー。Noneの場合はログを出力しない

    Returns:
        float: ギブス自由エネルギー。Δではない。
    """
    # ロガーがNoneの場合は無効なロガーを作成
    if logger is None:
        logger = ConditionalLogger(None, enabled=False)
    
    # ---構造に計算機を設定
    atoms = calc_input.structure
    # calc_mode に応じて適切な calculator を割り当て
    if calc_input.calc_mode == 'IdealGasThermo':
        atoms.calc = calculator_molecule
    else:
        atoms.calc = calculator_solid
    
    # 構造情報をログ出力
    logger.info("="*60)
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
            if hasattr(constraint, 'index'):
                logger.info(f"    固定原子数: {len(constraint.index) if hasattr(constraint.index, '__len__') else 1}")
    else:
        logger.info("拘束条件: なし")

    # ---　1. 最適化
    e_initial = atoms.get_potential_energy()
    logger.info(f"初期ポテンシャルエネルギー: {e_initial:.6f} eV")
    
    if calc_input.do_opt:
        if optimizer is None:
            optimizer = FIRELBFGS

        logger.info("構造最適化開始")
        opt_dyn = optimizer(atoms)
        opt_dyn.run(fmax=opt_fmax, steps=opt_maxsteps)
        logger.info("構造最適化完了")
    
    e_opt = atoms.get_potential_energy()
    logger.info(f"最適化後ポテンシャルエネルギー: {e_opt:.6f} eV")
    logger.info(f"最適化によるエネルギー変化: {e_opt - e_initial:.6f} eV")

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
    real_energies = [e.real if np.iscomplex(e) else e for e in vib_energies if (np.iscomplex(e) and e.real > 0) or (not np.iscomplex(e) and e > 0)]
    if real_energies:
        logger.info(f"実振動エネルギー範囲: {min(real_energies):.6f} - {max(real_energies):.6f} eV")
        logger.info(f"実振動エネルギー平均: {np.mean(real_energies):.6f} eV")

    # ---　3. 熱化学補正量計算
    if vib_energies.size == 0:
        logger.info("熱化学補正計算をスキップ (振動エネルギーなし)")
        g = e_opt
        thermal_correction = 0.0
    else:
        logger.info("熱化学補正計算開始")
        if calc_input.calc_mode == 'IdealGasThermo':
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
                spin=calc_input.spin_multiplicity
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
    logger.info("="*60)

    print(f"==========\ng_{atoms.symbols} = {g:.3f} eV\n==========")
    return g


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
    optimizer: Optimizer = None,
    opt_fmax: float = 0.05,
    opt_maxsteps: int = 3000,
    logger: Optional['ConditionalLogger'] = None,
):
    """
    ギブス自由エネルギーを計算する。CHEモデルに対応している

    Args:
        calculator_molecule (ase.calculators.calculator.Calculator): 分子用計算機
        calculator_solid (ase.calculators.calculator.Calculator): 固体用計算機
        reactants (list[CGFEInput|Literal["CHE"]]): 反応物。"CHE"を指定すると、CHEモデルによるギブス自由エネルギーを計算する。
        products (list[CGFEInput|Literal["CHE"]]): 生成物。"CHE"を指定すると、CHEモデルによるギブス自由エネルギーを計算する。
        temperature (float): 温度（K）
        pressure (float): 圧力（Pa）
        electrode_potential (float): 電極電位（V vs SHE）
        pH (float): pH
        optimizer (ase.optimize.optimize.Optimizer): 最適化エンジン
        opt_fmax (float): 最適化の閾値
        opt_maxsteps (int): 最適化の最大ステップ数
        logger (ConditionalLogger): ロガー。Noneの場合はログを出力しない

    Returns:
        float: 反応物と生成物のギブス自由エネルギーの差(ΔG)
    """
    # ロガーがNoneの場合は無効なロガーを作成
    if logger is None:
        logger = ConditionalLogger(None, enabled=False)
    
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
        g_h2 = calculate_g(
            calculator_molecule, # 気相分子向けの計算機
            calculator_solid,    # 固体表面向けの計算機
            CGFEInput(
                structure=molecule('H2'),
                calc_mode='IdealGasThermo',
                vibrate_indices=None,
                geometry='linear',
                symmetry_number=2,
                spin_multiplicity=1,
            ),
            temperature=temperature,
            pressure=pressure,
            optimizer=optimizer,
            opt_fmax=opt_fmax,
            opt_maxsteps=opt_maxsteps,
            logger=logger,
        )
        
        # CHEモデル計算の各項
        e = 1.0  # 素電荷（eV/V）
        kB = 8.617e-5  # ボルツマン定数（eV/K）
        
        term1 = 0.5 * g_h2
        term2 = -e * electrode_potential
        term3 = kB * temperature * np.log(10) * pH
        
        g_che = term1 + term2 + term3
        
        logger.info("CHEモデル計算詳細:")
        logger.info(f"  G(H2): {g_h2:.6f} eV")
        logger.info(f"  0.5 * G(H2): {term1:.6f} eV")
        logger.info(f"  -e * U: {term2:.6f} eV")
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
            g = calculate_g(
                calculator_molecule, # 気相分子向けの計算機
                calculator_solid,    # 固体表面向けの計算機
                reactant,
                temperature=temperature,
                pressure=pressure,
                optimizer=optimizer,
                opt_fmax=opt_fmax,
                opt_maxsteps=opt_maxsteps,
                logger=logger,
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
            g = calculate_g(
                calculator_molecule, # 気相分子向けの計算機
                calculator_solid,    # 固体表面向けの計算機
                product,
                temperature=temperature,
                pressure=pressure,
                optimizer=optimizer,
                opt_fmax=opt_fmax,
                opt_maxsteps=opt_maxsteps,
                logger=logger,
            )
            products_gs.append(g)

    # ---左辺右辺のギブス自由エネルギーを計算する
    g_reactants = sum(reactants_gs)
    g_products = sum(products_gs)
    g_delta = g_products - g_reactants
    
    logger.info("="*80)
    logger.info("反応ギブス自由エネルギー変化の最終計算")
    logger.info("反応物:")
    for i, g in enumerate(reactants_gs):
        species_name = "CHE" if reactants[i] == "CHE" else str(reactants[i].structure.symbols)
        logger.info(f"  {species_name}: {g:.6f} eV")
    logger.info(f"反応物合計: {g_reactants:.6f} eV")
    
    logger.info("生成物:")
    for i, g in enumerate(products_gs):
        species_name = "CHE" if products[i] == "CHE" else str(products[i].structure.symbols)
        logger.info(f"  {species_name}: {g:.6f} eV")
    logger.info(f"生成物合計: {g_products:.6f} eV")
    
    logger.info(f"ΔG = G(products) - G(reactants) = {g_products:.6f} - {g_reactants:.6f} = {g_delta:.6f} eV")
    logger.info("="*80)
    
    return g_delta


# %%
%%time

# 例: HO* + H+ + e- -> H2O + *
# ---ロガーの設定

# ログ設定
def setup_logger():
    """デバッグ用ログの設定"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"gibbs_energy_debug_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename, encoding='utf-8'),
            logging.StreamHandler()  # コンソールにも出力
        ]
    )
    
    base_logger = logging.getLogger(__name__)
    base_logger.info(f"デバッグログファイル: {log_filename}")
    return ConditionalLogger(base_logger, enabled=True)

# グローバルロガーの初期化
logger = setup_logger()

# ---構造の用意

# 4x4x4のPt(111)表面を作成
from ase.build import surface, molecule, bulk, add_adsorbate
from ase.atoms import Atoms
from ase.constraints import FixAtoms

# ---Pt(111)表面（4x4の表面、4層）
pt_bulk = bulk('Pt', a=3.92)  # Ptの格子定数
pt_slab = surface(pt_bulk, indices=(1,1,1), layers=4, vacuum=10.0)
pt_slab = pt_slab * (4, 4, 1)  # 4x4に拡張

# ---清浄なPt表面（生成物用）
clean_pt_slab = pt_slab.copy()

# ---OH分子
oh_molecule = molecule('OH')
oh_molecule.rotate(180,v="y")

# ---Pt表面にOHを吸着させた構造
# 表面の中央付近の原子を選択（top site吸着を想定）
oh_on_pt = pt_slab.copy()
top_atoms = [atom for atom in oh_on_pt if atom.position[2] > oh_on_pt.positions[:,2].max() - 1.0]
center_position = top_atoms[len(top_atoms)//2].position
add_adsorbate(oh_on_pt, oh_molecule, height=2.0, position=center_position[0:2])
oh_on_pt.center(vacuum=10.0, axis=2)

# 下を固定

def fix_with_layer(atoms: Atoms, layers: int)-> Atoms:
    _atoms = atoms.copy()
    # Z座標に基づいて原子層を特定し、下層を固定する
    z_coords = np.unique(np.round(_atoms.positions[:, 2], decimals=4))
    z_coords.sort()
    if len(z_coords) > 2:
        fix_threshold = z_coords[layers-1]
        constraint = FixAtoms(mask=_atoms.positions[:, 2] <= fix_threshold)
        _atoms.set_constraint(constraint)
    return _atoms

oh_on_pt = fix_with_layer(oh_on_pt, 2)
clean_pt_slab = fix_with_layer(clean_pt_slab, 2)  # 拘束条件を統一

# H2O分子
h2o_molecule = molecule('H2O')

# ---計算の実行
logger.info("メイン計算開始")
logger.info("反応: HO* + H+ + e- -> H2O + *")
delta_g = calculate_delta_g(
    calculator_molecule=calc_molecule,
    calculator_solid=calc_solid,
    reactants=[
        CGFEInput(
            structure=oh_on_pt,  # Pt-OH構造
            calc_mode='HarmonicThermo',
            vibrate_indices=[64,65,19,35,51,23,39,55],
        ),
        "CHE",  # H+ + e-
    ],
    products=[
        CGFEInput(
            structure=h2o_molecule,  # H2O（気相）
            calc_mode='IdealGasThermo',
            vibrate_indices=None,
            geometry="nonlinear",
            symmetry_number=2,
            spin_multiplicity=1,
        ),
        CGFEInput(
            structure=clean_pt_slab,  # 清浄なPt表面
            calc_mode='HarmonicThermo',  # 固体表面なのでHarmonicThermo
            vibrate_indices=[],
        ),
    ],
    temperature=298.15,
    pressure=101325,
    electrode_potential=0.0,
    pH=0,
    opt_maxsteps=1500,
    logger=logger,
)


# ---結果の表示
print(f"ΔG for HO* + H+ + e- -> H2O + * = {delta_g:.3f} eV")


# %% [markdown]
# ### 🟡格子定数を計算する

# %%
from dataclasses import dataclass
from typing import Iterator
from ase import Atoms, Atom
from ase.calculators.calculator import Calculator
from ase.optimize.optimize import Optimizer
from matlantis_features.ase_ext.optimize import FIRELBFGS
from ase.constraints import UnitCellFilter


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
        return (f"LatticeConstant(a={self.a:.4f}Å, b={self.b:.4f}Å, c={self.c:.4f}Å, "
                f"α={self.alpha:.2f}°, β={self.beta:.2f}°, γ={self.gamma:.2f}°)")


def optimize_lattice_constant(
    atoms: Atoms,
    calculator: Calculator | None = None,
    optimizer: Optimizer = FIRELBFGS,
    fmax: float = 0.01,
    steps: int | None = None
) -> LatticeConstant:
    """格子定数を最適化してLatticeConstantオブジェクトを返す
    
    この関数は、与えられた原子構造の格子定数を最適化し、結果をLatticeConstantデータクラス
    として返します。最適化にはUnitCellFilterを使用してセルの形状と体積の両方を最適化します。
    
    Args:
        atoms (Atoms): 最適化する原子構造
        calculator (Calculator | None): 使用する計算機。Noneの場合は既存の計算機を使用
        optimizer (Optimizer): 使用する最適化アルゴリズムクラス（デフォルト: FIRELBFGS）
        fmax (float): 収束条件（原子にかかる力の最大値） [eV/Å]
        steps (int | None): 最大ステップ数。Noneの場合は制限なし
        
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
        raise RuntimeError("計算機が設定されていません。calculatorを指定するか、atomsに計算機を設定してください。")
    
    # --- 原子構造のコピーを作成（元の構造を保持） ---
    atoms_copy = atoms.copy()
    
    # --- 計算機の設定 ---
    if calculator is not None:
        atoms_copy.calc = calculator
    
    # --- UnitCellFilterでセル最適化を有効化 ---
    unit_cell_filter = UnitCellFilter(atoms_copy)
    
    # --- 最適化アルゴリズムの設定と実行 ---
    optimizer = optimizer(unit_cell_filter)
    
    # ステップ数の制限がある場合は設定
    if steps is not None:
        optimizer.run(fmax=fmax, steps=steps)
    else:
        optimizer.run(fmax=fmax)
    
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
        gamma=lengths_and_angles[5]
    )


# %% [markdown]
# ## 参考文献
# - [MatlantisによるASEチュートリアル](https://docs.matlantis.com/atomistic-simulation-tutorial/ja/)
#     - 一通り見ると流れがわかる
# - [Materials Project](https://next-gen.materialsproject.org/)
#     - cifファイルのデータベース的な

# %%
oh_on_pt.positions[:,2]

# %% [markdown]
# ## おあそびスペース

# %%
from ase.build import molecule, bulk, surface
from ase.io import read, write

# ---分子
co=molecule('CO')

# ---結晶
cu=bulk('Cu')

# ---表面
slab=surface(
    lattice=cu,
    indices=(3,2,1),
    layers=4,
    vacuum=10.0
)*(5,5,1)

view_ngl(slab)

# %%
layers = separate_layers(slab)
layers

# %%
new_slab = substitute_elements(slab, layers[-1], "Au")
view_ngl(new_slab)

# %%


new_slab = substitute_elements(slab, find_central_atom(layers[-1]), "Au")
adsorbed_slab = new_slab.copy()
view_ngl(new_slab)

# %%
from ase import Atoms
from ase.build import fcc111, add_adsorbate

add_adsorbate(adsorbed_slab, co, 3, position=find_central_atom(layers[-1]).position[:2])
view_ngl(adsorbed_slab, representations=["ball+stick"])

# %%
get_appended_atom_indices(new_slab, adsorbed_slab)

# %%
adsorbed_slab.constraints

# %%
new_adsorbed_slab = fix_layers(adsorbed_slab, 1)
new_adsorbed_slab.constraints

# %%
moved_slab = move_atoms(adsorbed_slab, [100,101], (0,1,0),2)
view_ngl(moved_slab)

# %%


# %%
randomized_slab = adsorbed_slab.copy()

for atom in randomized_slab:
    move_atoms(randomized_slab, atom, np.random.rand(3),np.random.rand()*100, inplace=True)
view_ngl(randomized_slab)

# %%
from pfcc_extras.structure.ase_rdkit_converter import smiles_to_atoms, atoms_to_smiles

atoms = smiles_to_atoms("Oc1ccccc1")
view_ngl(atoms, representations=["ball+stick"])


