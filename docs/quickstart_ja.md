# 日本語クイックスタート — v1.0.0

## 1. 環境作成
```bash
mamba env create -f environment.yml
conda activate peerj-organelle-audit
make test
```

## 2. ローカル入力の設定
`config/species/silene.yaml` と `config/species/chicken.yaml` の `local_inputs` のみを実データの場所に合わせて変更します。公開版には解析機固有の絶対パスは含めません。

## 3. 実行前確認
```bash
make plan
make autotune
make doctor
make local-check
```
未解決runがある場合は `work/<species>/local_inputs/unresolved_runs.tsv` を確認します。

## 4. 参照配列と解析
```bash
make references
make pilot
make full
make readiness
```
参照配列準備時には `work/<species>/references/reports/` にマスク範囲、統計、ハッシュ等のrun固有provenanceが生成されます。再解析時にはこのディレクトリも保存してください。

## 5. 一時ファイルの整理
```bash
make clean-tmp
```
外部FASTQと参照FASTAは削除しません。
