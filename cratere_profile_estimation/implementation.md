以下は、**最尤推定で (\widehat{\kappa t}) を求め、その (\widehat{\kappa t}) で初期プロファイルを拡散して現在形状を得る**ための実装手順です。
オンボード実装を想定し、重い拡散計算は原則として事前計算に寄せます。

---

# 0. 実装全体の構成

処理は大きく2段階に分ける。

## オフライン事前計算

地上側で行う。

[
(D,\ x_0,\ s)
\rightarrow
g_D(x_0,s)
]

および

[
(D,\ x_0,\ s)
\rightarrow
h(r;D,x_0,s)
]

のルックアップテーブルを作る。

ここで、

[
x_0=(d/D)_0
]

は初期深さ比、

[
s=\kappa t
]

は拡散量、

[
g_D(x_0,s)
]

は拡散後の現在深さ比、

[
h(r;D,x_0,s)
]

は拡散後の半径方向プロファイルである。

synthterrain 側でも、初期断面を作り、2次元線形拡散をかけ、その後に高度場の最大値と最小値から (d/D) を計算する流れになっている。 初期断面は Stopar 的な fresh (d/D) と Fassett & Thomson 系の区分多項式プロファイルから構成される。

## オンボード推定

ローバ上で行う。

LiDARから

[
D
]

および可能なら

[
\tilde{S}
]

を得る。ここで (\tilde{S}) は観測された内壁傾斜角である。

その後、事前計算テーブルを使って

[
\widehat{\kappa t}
]

を最尤推定し、対応するプロファイルを取り出す。

---

# 1. オフライン事前計算

## Step 1. 直径グリッドを作る

対象とするクレーター直径範囲を決める。

例：

[
D \in [1,\mathrm{m}, 200,\mathrm{m}]
]

または、ローバ運用上必要な範囲に限定する。

直径グリッドは対数刻みがよい。

[
D_i \in {D_1,D_2,\dots,D_N}
]

例：

```text
D_grid = logspace(log10(D_min), log10(D_max), N_D)
```

小クレーターではスケール依存が強いので、等間隔より対数間隔の方が扱いやすい。

---

## Step 2. 初期 (d/D) 分布を定義する

各直径 (D_i) に対して、

[
p_0(x_0\mid D_i)
]

を定義する。

初期値の平均は Stopar / Fassett Table 2 系の値を使う。

例：

```text
D >= 400 m : μ0 = 0.21
D >= 200 m : μ0 = 0.17
D >= 100 m : μ0 = 0.15
D >= 40 m  : μ0 = 0.13
D >= 10 m  : μ0 = 0.11
D < 10 m   : μ0 = 0.10
```

synthterrain でもこのような Stopar 段階値を用いて fresh (d/D) を決める整理になっている。

分布は例えば正規分布でよい。

[
x_0\sim \mathcal{N}(\mu_0(D),\sigma_0^2)
]

ただし、物理的に不自然な値を避けるため、範囲を切る。

[
x_0 > 0
]

[
x_0 \le x_{0,\max}
]

実装上は離散グリッドにする。

```text
x0_grid = linspace(x0_min, x0_max, N_x0)
p0[D_i, x0_j] = NormalPDF(x0_j; μ0(D_i), σ0)
normalize p0 over x0_grid
```

---

## Step 3. (\kappa t) グリッドを作る

拡散量を

[
s=\kappa t
]

と置く。

これを離散化する。

[
s_k \in {s_1,s_2,\dots,s_K}
]

例：

```text
s_grid = linspace(0, s_max(D_i), N_s)
```

ここで (s_\mathrm{max}(D_i)) は、クレーターがほぼ消失する程度までの拡散量にする。

重要なのは、ここでの (s=\kappa t) は絶対年代ではなく、形状劣化量として扱うこと。

---

## Step 4. 初期プロファイルを生成する

各

[
(D_i,x_{0j})
]

について、初期プロファイル

[
h_0(r;D_i,x_{0j})
]

を生成する。

synthterrain 型なら、処理は次の通り。

```text
D_i, x0_j
→ Fassett & Thomson 系の区分多項式プロファイル
→ depth = x0_j * D_i
→ floor clipping
→ 2D 高度場 u0(x,y) に展開
```

断面モデル自体は統計モデルではなく、直径 (D) と初期 (d/D) から決まる経験的・決定論的プロファイルである。

---

## Step 5. 拡散計算を実行する

各

[
(D_i,x_{0j},s_k)
]

について、初期高度場を拡散させる。

拡散方程式は

[
\frac{\partial u}{\partial t}
=============================

\kappa(D)\nabla^2u
]

だが、実装上は (s=\kappa t) を直接使う。

離散更新は概念的に

[
u^{n+1}
=======

u^n
+
\Delta s,\nabla^2u^n
]

でよい。

synthterrain では、陽解法で2次元高度場を更新し、格子間隔から安定なステップ数を決める実装になっている。

---

## Step 6. 拡散後の (d/D) を計算する

拡散後の高度場を

[
u_{i,j,k}(x,y)
]

とする。

そこから

[
g_D(x_{0j},s_k)
===============

\frac{\max(u_{i,j,k})-\min(u_{i,j,k})}{D_i}
]

を計算する。

synthterrain でも、各ステップ後に高度場全体の relief を直径で割って (d/D) を計算している。

保存するテーブルは最低限これ。

```text
G[i, j, k] = g_Di(x0_j, s_k)
```

さらにプロファイルも保存する。

```text
H[i, j, k, r] = h(r; D_i, x0_j, s_k)
```

2D全体を保存すると重いので、軸対称断面

[
h(r)
]

だけ保存すればよい。

---

## Step 7. Fassett Fig. 3 の現在 (d/D) 分布を用意する

Fassett 2022 Fig. 3 から得られる現在 (d/D) 分布を

[
p_c(x_c\mid D)
]

として用意する。

実装ではCDFからPDFを作る。

```text
CDF_Fig3[D_bin, x]
→ numerical derivative
→ PDF pc[D_bin, x]
→ smoothing
→ normalize
```

または、CDFのまま使い、区間確率で評価してもよい。

```text
pc_value = probability_density_or_bin_probability(x_c | D)
```

Fassett 2022 では、小クレーターの平衡状態、すなわち生成と消失の釣り合いを使ってクレーター寿命や劣化量を扱う。 Fig. 3 は、平衡状態のクレーターを合成的に拡散させたときの (d/D) 分布として使う。

---

# 2. オンボード推定

## Step 8. LiDARから直径 (D) を推定する

LiDAR点群からリム候補を抽出し、円フィッティングで直径を求める。

```text
input: point cloud
detect rim points
fit circle
D_obs = 2 * fitted_radius
```

直径の不確実性を無視するなら、

[
D=D_\mathrm{obs}
]

としてよい。

不確実性を入れるなら、

[
D\sim \mathcal{N}(D_\mathrm{obs},\sigma_D^2)
]

とする。ただし初期実装では (D) 固定で十分。

---

## Step 9. 内壁傾斜角 (\tilde{S}) を推定する

クレーター内部底面は見えなくても、手前側の内壁やリム近傍の一部は見える可能性がある。

その点群から局所平面または半径方向断面をフィットし、内壁傾斜を推定する。

```text
select visible inner-wall points
remove outliers
fit local plane or radial profile
S_obs = slope angle in degrees
```

この観測値を

[
\tilde{S}
]

とする。

観測誤差を

[
\sigma_S
]

とする。

[
\sigma_S^2
==========

\sigma_\mathrm{LiDAR}^2
+
\sigma_\mathrm{regression}^2
+
\sigma_\mathrm{partial\ wall}^2
+
\sigma_\mathrm{terrain}^2
]

のように大きめに置く。遠距離で壁面の一部しか見えないなら、(\sigma_S) は小さくしすぎない。

---

## Step 10. 対応する直径ビンを選ぶ

観測直径 (D_\mathrm{obs}) に最も近い事前計算テーブルを選ぶ。

```text
i = nearest_index(D_grid, D_obs)
```

より滑らかにするなら、隣接する2つの直径ビンで補間する。

```text
G_interp = interpolate over D
H_interp = interpolate over D
p0_interp = interpolate over D
pc_interp = interpolate over D
```

初期実装では最近傍でよい。

---

# 3. 最尤推定

## Step 11. 観測なしの場合の尤度を計算する

内壁傾斜 (S) が使えない場合は、各 (s_k) について

[
L(s_k)
======

\sum_j
p_c(G[i,j,k]\mid D_i)
p_0(x_{0j}\mid D_i)
\Delta x_0
]

を計算する。

実装では離散和でよい。

```python
for k in range(N_s):
    L[k] = 0
    for j in range(N_x0):
        x_pred = G[i, j, k]
        L[k] += pc(x_pred, D_i) * p0[i, j] * dx0
```

この場合、推定される (\kappa t) は、観測個体に強く条件付けされた値ではなく、

> その直径のクレーターとして、Fassett Fig. 3 の現在 (d/D) 分布に最も整合する代表的拡散量

である。

---

## Step 12. 内壁傾斜 (S) ありの場合の尤度を計算する

Li et al. 系の関係

[
S=a(d/D)
]

を使う。

ここで

[
a=151.377
]

とする。文献整理では、この関係は内壁傾斜から (d/D) を推定する回帰モデルとして整理されている。

各 (s_k) について、

[
L(s_k)
======

\sum_j
\mathcal{N}
\left(
\tilde{S};
aG[i,j,k],
\sigma_S^2
\right)
p_c(G[i,j,k]\mid D_i)
p_0(x_{0j}\mid D_i)
\Delta x_0
]

を計算する。

実装例：

```python
a = 151.377

for k in range(N_s):
    L[k] = 0.0
    for j in range(N_x0):
        x_pred = G[i, j, k]              # predicted current d/D
        S_pred = a * x_pred             # predicted inner-wall slope
        obs_like = normal_pdf(S_obs, mean=S_pred, std=sigma_S)
        prior_current = pc(x_pred, D_i)
        prior_initial = p0[i, j]
        L[k] += obs_like * prior_current * prior_initial * dx0
```

数値アンダーフローを避けるならログ空間で計算する。

```python
log_terms = []
for j in range(N_x0):
    x_pred = G[i, j, k)
    S_pred = a * x_pred
    log_terms.append(
        log_normal_pdf(S_obs, S_pred, sigma_S)
        + log_pc(x_pred, D_i)
        + log_p0[i, j]
        + log(dx0)
    )
logL[k] = logsumexp(log_terms)
```

---

## Step 13. 物理制約を入れる

通常の拡散では (d/D) は増えないので、

[
G[i,j,k]\le x_{0j}
]

を満たさないものは除外する。

実装：

```python
if x_pred > x0_grid[j]:
    continue
```

または低尤度にする。

```python
if x_pred > x0_grid[j]:
    log_term = -inf
```

さらに、

[
s_k\ge0
]

は必須。

クレーター消失後に相当する過大な (s_k) も除外してよい。

---

## Step 14. 最尤 (\widehat{\kappa t}) を選ぶ

尤度最大のインデックスを選ぶ。

```python
k_hat = argmax(L)
s_hat = s_grid[k_hat]
```

つまり、

[
\widehat{\kappa t}
==================

s_{k_\mathrm{max}}
]

である。

必要なら、尤度曲線を正規化して

[
p(s_k\mid \text{data})
\propto L(s_k)
]

のように扱い、信頼区間も出せる。

```python
posterior_like = L / sum(L)
s_mean = sum(s_grid * posterior_like)
s_ci = percentile_from_cdf(s_grid, posterior_like, [0.1, 0.5, 0.9])
```

厳密にはこれはMLの尤度正規化であり、MAPや完全ベイズではないが、オンボードの不確実性指標として有用。

---

# 4. 現在形状プロファイルの取得

## Step 15. 代表 (x_0) を選ぶ

(\widehat{\kappa t}=s_{\hat{k}}) が決まったら、次にプロファイル生成に使う (x_0) を決める。

方法は3つある。

### 方法A：初期分布の中央値を使う

```python
x0_hat = median_of_p0(D_i)
j_hat = nearest_index(x0_grid, x0_hat)
```

単純だが、観測との整合性はやや弱い。

### 方法B：同時に最も尤もらしい (x_0) を使う

尤度計算で最大寄与した (x_0) を使う。

[
\hat{x}_0
=========

\arg\max_{x_{0j}}
\left[
p_\mathrm{obs}(\tilde{S}\mid aG[i,j,\hat{k}])
p_c(G[i,j,\hat{k}]\mid D_i)
p_0(x_{0j}\mid D_i)
\right]
]

実装：

```python
j_hat = argmax_over_j(term[j, k_hat])
```

単一プロファイルを出すならこれがよい。

### 方法C：(x_0) を周辺化して分位プロファイルを出す

最も推奨。

(\hat{s}) 固定で、各 (x_{0j}) の重みを計算する。

[
w_j
\propto
p_\mathrm{obs}(\tilde{S}\mid aG[i,j,\hat{k}])
p_c(G[i,j,\hat{k}]\mid D_i)
p_0(x_{0j}\mid D_i)
]

正規化して、

[
\sum_j w_j=1
]

とする。

```python
w[j] = obs_like * pc_value * p0_value
w = w / sum(w)
```

その重みでプロファイル分布を作る。

---

## Step 16. プロファイルを取り出す

単一プロファイルなら、

```python
h_hat = H[i, j_hat, k_hat, :]
```

分位プロファイルなら、各半径 (r_m) ごとに、重み付き分位を計算する。

```python
for m in range(N_r):
    values = H[i, :, k_hat, m]     # all x0 profiles at fixed s_hat
    h10[m] = weighted_quantile(values, w, 0.10)
    h50[m] = weighted_quantile(values, w, 0.50)
    h90[m] = weighted_quantile(values, w, 0.90)
```

出力：

[
h_{10}(r),\quad h_{50}(r),\quad h_{90}(r)
]

この方が、初期 (d/D) の不確実性を反映できる。

---

# 5. 最終出力

最低限の出力は以下。

```text
D_obs
S_obs, if available
kappa_t_hat
x_current_pred = g_D(x0_hat, kappa_t_hat)
h50(r)
h10(r), h90(r)
likelihood curve L(s)
quality flags
```

品質フラグとしては、以下を持つとよい。

```text
flag_no_slope_observation
flag_low_likelihood
flag_multimodal_likelihood
flag_out_of_distribution_D
flag_slope_outlier
flag_occluded_inner_wall
```

---

# 6. 推奨する実装の擬似コード

```python
def estimate_crater_profile(D_obs, S_obs=None, sigma_S=None):
    # 1. choose/interpolate diameter bin
    i = nearest_index(D_grid, D_obs)

    logL = np.full(N_s, -np.inf)
    all_log_terms = np.full((N_x0, N_s), -np.inf)

    # 2. likelihood over kappa*t
    for k, s in enumerate(s_grid):
        log_terms = []

        for j, x0 in enumerate(x0_grid):
            x_pred = G[i, j, k]  # predicted current d/D

            # physical constraint
            if x_pred > x0:
                continue

            lp0 = log_p0[i, j]
            lpc = log_pc(x_pred, D_grid[i])

            if S_obs is not None:
                S_pred = 151.377 * x_pred
                lobs = log_normal_pdf(S_obs, S_pred, sigma_S)
            else:
                lobs = 0.0

            log_term = lobs + lpc + lp0 + np.log(dx0)
            log_terms.append(log_term)
            all_log_terms[j, k] = log_term

        logL[k] = logsumexp(log_terms)

    # 3. maximum likelihood estimate
    k_hat = np.argmax(logL)
    s_hat = s_grid[k_hat]

    # 4. weights over x0 at fixed s_hat
    logw = all_log_terms[:, k_hat]
    w = np.exp(logw - logsumexp(logw))

    # 5. representative current d/D
    x_current_values = G[i, :, k_hat]
    x_current_hat = weighted_quantile(x_current_values, w, 0.50)

    # 6. profile quantiles
    h10 = np.zeros(N_r)
    h50 = np.zeros(N_r)
    h90 = np.zeros(N_r)

    for m in range(N_r):
        vals = H[i, :, k_hat, m]
        h10[m] = weighted_quantile(vals, w, 0.10)
        h50[m] = weighted_quantile(vals, w, 0.50)
        h90[m] = weighted_quantile(vals, w, 0.90)

    return {
        "D": D_obs,
        "kappa_t_hat": s_hat,
        "x_current_hat": x_current_hat,
        "r": r_grid,
        "h10": h10,
        "h50": h50,
        "h90": h90,
        "logL": logL,
        "s_grid": s_grid,
    }
```

---

# 7. 実装上の優先順位

最初に作るべき最小構成はこれ。

```text
1. D のみ入力
2. p0(x0|D) と pc(xc|D) を用意
3. G[D,x0,s] テーブルを作る
4. L(s)=Σ pc(G)p0 を最大化
5. h50(r) を出す
```

次に追加する。

```text
6. LiDARから S を推定
7. S_pred = 151.377 * G を計算
8. Normal(S_obs; S_pred, sigma_S) を尤度に追加
9. h10/h50/h90 を出す
```

最後に追加する。

```text
10. D の不確実性
11. 尤度の多峰性チェック
12. MAP化する場合の p(s|D)
13. 楕円・重複・局所斜面外れ値の除外
```

---

# 8. 実装の要点

この方法では、オンボードで拡散方程式を毎回解く必要はない。

重い処理はオフラインで、

[
(D,x_0,\kappa t)
\rightarrow
(d/D)_\mathrm{current},\ h(r)
]

のテーブルにしておく。

オンボードでは、

[
L(\kappa t)
]

をテーブル参照と足し算で評価し、最大値を選ぶだけでよい。

したがって、実装の本質は次である。

[
\boxed{
\text{拡散方程式は forward table 作成に使う}
}
]

[
\boxed{
\kappa t \text{ はオンボードで最尤推定する}
}
]

[
\boxed{
\widehat{\kappa t} \text{ に対応する拡散済みプロファイルを出力する}
}
]
