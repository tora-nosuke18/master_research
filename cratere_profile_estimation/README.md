# Crater Profile Estimation

`methodology.md` に基づき、月面クレーターの直径と任意の内壁傾斜観測から、
現在の軸対称形状プロファイルを統計的に推定する Python 実装です。

## Model

- 初期形状: `synthterrain` の `FTmod_Crater` と同じ Fassett-Thomson
  区分3次多項式と Stopar fresh `d/D`
- 初期 `d/D`: Stopar 段階値を平均、標準偏差 `0.02` の正規分布
- 劣化: 固定外周境界を持つ2次元線形地形拡散
- 現在 `d/D`: Fassett et al. (2022) Figure 3 のデジタイズ近似
- 傾斜観測: `S = 151.377(d/D)` と正規観測誤差
- 推定量: 初期 `d/D` を周辺化した `kappa*t` の最尤値
- 出力: 推定後の半径方向高さの10%、50%、90%分位

計算では無次元量 `lambda = (kappa*t)/D^2` を探索し、結果を
`diffusion_amount_m2 = kappa*t` として平方メートル単位で返します。

## Installation

```bash
python -m pip install -e .
```

## Usage

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
- `output/summary.json`: 最尤劣化量、推定 `d/D`、使用した直径ビン

保存先を変更するときだけ `--output`、`--plot`、`--summary` を指定します。

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
