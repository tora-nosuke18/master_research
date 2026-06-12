# Crater Profile Estimation

`methodology.md` と `implementation.md` に基づき、月面クレーターの直径と
任意の内壁傾斜観測から現在の軸対称形状プロファイルを推定する実装です。

処理は、重い拡散計算を行うオフライン事前計算と、ルックアップテーブルだけを
参照するオンボード推定に分離されています。

## Model

- 初期形状: `synthterrain` の `FTmod_Crater` と同じ Fassett-Thomson
  区分3次多項式と Stopar fresh `d/D`
- 初期 `d/D`: Stopar 段階値を平均、標準偏差 `0.02` の正規分布
- オフライン劣化計算: 固定外周境界を持つ2次元線形地形拡散
- 事前計算テーブル: `G[D,x0,s]` と `H[D,x0,s,r]`
- オンボード推定: 最近傍直径のテーブル参照とログ尤度計算
- 現在 `d/D`: Fassett et al. (2022) Figure 3 のデジタイズ近似
- 傾斜観測: `S = 151.377(d/D)` と正規観測誤差
- 推定量: 初期 `d/D` を周辺化した `kappa*t` の最尤値
- 出力: 推定後の半径方向高さの10%、50%、90%分位

事前計算では無次元量 `lambda = (kappa*t)/D^2` を使って安定に拡散を計算し、
テーブルには直径ごとの `s = kappa*t` を平方メートル単位で保存します。

## Installation

```bash
python -m pip install -e .
```

## Offline Precomputation

標準テーブルは `data/crater_lookup.npz` に置かれます。再生成する場合:

```bash
build-crater-lookup
```

範囲や解像度を変更する例:

```bash
build-crater-lookup \
  --diameter-min 1 \
  --diameter-max 177 \
  --diameter-points 24 \
  --x0-points 31 \
  --diffusion-points 101 \
  --output data/crater_lookup.npz
```

この処理だけが2次元拡散方程式を解き、次の配列をNPZへ保存します。

- `diameter_grid_m`: 対数間隔の直径グリッド
- `x0_grid`: 初期 `d/D` グリッド
- `s_grid_m2`: 直径ごとの `kappa*t` グリッド
- `current_dd`: `G[D,x0,s]`
- `profiles_m`: `H[D,x0,s,r]`

## Onboard Estimation

直径だけを使う場合:

```bash
estimate-crater-profile \
  --diameter 20
```

内壁傾斜を使う場合:

```bash
estimate-crater-profile \
  --diameter 20 \
  --slope 10 \
  --sigma-slope 2
```

実行結果は、既定でリポジトリ直下の `output/` に保存されます。

- `output/profile.csv`: 半径方向の10%、50%、90%分位プロファイル
- `output/profile.png`: 左右対称断面と10-90%不確実性帯の画像
- `output/likelihood.csv`: `kappa*t` ごとの対数尤度と正規化尤度
- `output/summary.json`: 最尤劣化量、推定 `d/D`、使用した直径ビン

保存先を変更するときだけ `--output`、`--plot`、`--summary` を指定します。

推定時には拡散方程式を解きません。`--lookup` で指定した事前計算テーブルから
最も近い直径を選び、`G` と `H` を参照します。形状の長さスケールのみ観測直径へ
合わせて補間します。

`--sigma-slope` はLiDAR誤差、回帰誤差、部分壁観測、局所地形差を合わせた
標準偏差です。方法論で数値が定められていないため、傾斜使用時は明示指定が
必要です。

## Figure 3 Digitization

Fassett et al. (2022), DOI `10.1029/2022JE007510`, Figure 3 の曲線をPDFから
3倍で描画し、Matplotlib標準色の画素をプロット軸座標へ変換しました。
各曲線のCDF 0.05刻みの `d/D` 値を
[`distributions.py`](src/crater_profile_estimation/distributions.py) に固定しています。
確率密度は、この等重み分位点に帯域幅 `0.003` のGaussian kernelを適用した
近似です。

対応する直径ビンは次の3区間です。

- `1 <= D < 10 m`
- `10 <= D < 39.8 m`
- `39.8 <= D <= 177 m`

Figure 3 の範囲外へ暗黙に外挿しないため、`D < 1 m` または `D > 177 m` は
エラーになります。ビン境界では曲線が離散的に切り替わります。

## Tests

```bash
pytest -q
```

## Important Limitations

- 出力は個別クレーターの真値ではなく、平衡集団統計と線形拡散に整合する推定です。
- Figure 3 は元数値ではなく図からの近似値です。
- 円形・軸対称、線形拡散、固定外周境界を仮定します。
- 分位プロファイルは最尤 `lambda` での初期 `d/D` の条件付き重みから計算します。
- `kappa*t` は形状劣化量であり、絶対年代ではありません。
- 標準テーブルとFigure 3の適用範囲は `1 <= D <= 177 m` です。
- 直径方向の初期実装は最近傍テーブル選択であり、`G` の線形補間は行いません。

## Quality Flags

`summary.json` には次の品質フラグが出力されます。

- `flag_no_slope_observation`
- `flag_low_likelihood`
- `flag_multimodal_likelihood`
- `flag_out_of_distribution_D`
- `flag_slope_outlier`
- `flag_occluded_inner_wall`

## synthterrain Comparison

`synthterrain` でクレーターを生成し、中央断面の半径 `0.2R-0.9R` を一次回帰して
得た内壁傾斜と生成時の直径を推定器へ入力し、真の断面と推定結果を比較できます。

```bash
python -m pip install -e '.[comparison]'
compare-synthterrain-profile \
  --diameter 20 \
  --lambda-true 0.015 \
  --sigma-slope 2
```

出力画像は、左側にsynthterrainの2次元高度場、右側に中央断面、推定中央値、
10-90%分位帯を表示します。
