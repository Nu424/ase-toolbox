# ase-toolbox
よく使うコードをまとめて、より楽に化学シミュレーションをできるようにするためのPythonのASE用ヘルパー関数集みたいな。

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
ーーーーーーーーーーーーーー

```