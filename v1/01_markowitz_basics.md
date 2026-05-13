# 第1章 Markowitz の平均分散分析

> *"Diversification is both observed and sensible; a rule of behavior which does not imply the superiority of diversification must be rejected both as a hypothesis and as a maxim."*  — H. M. Markowitz (1952)

## 1.1 問題設定

投資家が $n$ 個の危険資産に資金を配分する。リターン $R$ について $\mu = \mathbb{E}[R]$, $\Sigma = \mathrm{Cov}(R) \succ 0$ を所与とする。ポートフォリオ $w$ の

- 期待リターン：$\mathbb{E}[R_w] = w^\top \mu$
- 分散（リスク）：$\mathrm{Var}(R_w) = w^\top \Sigma w$

Markowitz の中心思想は次の二つ：

1. **「リスク」を分散で測る**。リターンの確率分布は平均と分散の二つの値だけで（投資家にとって）十分である。
2. **平均と分散のトレードオフを明示的に最適化する**。これにより「分散投資の価値」が定量化される。

## 1.2 平均分散効率性

**定義 1.1**（平均分散効率的ポートフォリオ）
ポートフォリオ $w^\star$ が **効率的（efficient）** であるとは、これと同じか高い期待リターンを持ち、かつ厳密に小さい分散を持つ別のポートフォリオが存在しないことをいう。すなわち以下の問題の解として特徴付けられる：

$$
\min_{w \in \mathbb{R}^n} \tfrac{1}{2} w^\top \Sigma w \quad \text{s.t.}\quad w^\top \mu = \mu_P,\ \mathbf{1}^\top w = 1. \tag{P}
$$

## 1.3 効率的ポートフォリオの陽な表現

問題 (P) は第0章の二次計画の枠組みに当てはまる。$A = \begin{pmatrix} \mu^\top \\ \mathbf{1}^\top \end{pmatrix}$, $b = \begin{pmatrix} \mu_P \\ 1 \end{pmatrix}$ とおき、次の **基本スカラー量** を定義する：

$$
A_0 := \mathbf{1}^\top \Sigma^{-1} \mathbf{1}, \qquad
B_0 := \mathbf{1}^\top \Sigma^{-1} \mu, \qquad
C_0 := \mu^\top \Sigma^{-1} \mu, \qquad
D_0 := A_0 C_0 - B_0^2.
$$

**補題 1.2**
$\mu$ が $\mathbf{1}$ と平行でない（つまり全資産の期待リターンが同一でない）とき $D_0 > 0$。

*証明方針*：$\Sigma^{-1} \succ 0$ より $(x, y) \mapsto x^\top \Sigma^{-1} y$ は内積、コーシー・シュワルツより $B_0^2 = (\mathbf{1}^\top \Sigma^{-1} \mu)^2 \le A_0 C_0$。等号は $\mu \parallel \mathbf{1}$ のときのみ。$\square$

**定理 1.3**（効率的ポートフォリオの陽な解）
問題 (P) の最適解は
$$
\boxed{\;
w^\star(\mu_P) = \frac{C_0 - B_0 \mu_P}{D_0}\, \Sigma^{-1} \mathbf{1} + \frac{A_0 \mu_P - B_0}{D_0}\, \Sigma^{-1} \mu
\;}
$$
で与えられ、対応する最小分散は
$$
\sigma_P^2(\mu_P) = w^{\star\top} \Sigma w^\star = \frac{A_0 \mu_P^2 - 2 B_0 \mu_P + C_0}{D_0}. \tag{1.1}
$$

*証明方針*：Lagrangian $L(w, \lambda_1, \lambda_2) = \tfrac{1}{2} w^\top \Sigma w - \lambda_1 (w^\top \mu - \mu_P) - \lambda_2(\mathbf{1}^\top w - 1)$。
- $\partial_w L = 0 \Rightarrow w = \lambda_1 \Sigma^{-1}\mu + \lambda_2 \Sigma^{-1}\mathbf{1}$。
- 制約 $w^\top \mu = \mu_P, \mathbf{1}^\top w = 1$ に代入し $\lambda_1, \lambda_2$ について連立一次方程式を解く：
$$
\begin{pmatrix} C_0 & B_0 \\ B_0 & A_0 \end{pmatrix} \begin{pmatrix} \lambda_1 \\ \lambda_2 \end{pmatrix} = \begin{pmatrix} \mu_P \\ 1 \end{pmatrix}.
$$
- 行列式 $D_0$ より $\lambda_1 = (A_0 \mu_P - B_0)/D_0$、$\lambda_2 = (C_0 - B_0 \mu_P)/D_0$。
- 分散は $w^{\star\top}\Sigma w^\star = \lambda_1 \mu_P + \lambda_2$ より整理して (1.1) を得る。$\square$

**系 1.4**（大域最小分散ポートフォリオ；GMV）
(1.1) を $\mu_P$ について最小化すると $\mu_P^{\rm GMV} = B_0 / A_0$、対応する重みは
$$
w^{\rm GMV} = \frac{\Sigma^{-1} \mathbf{1}}{\mathbf{1}^\top \Sigma^{-1} \mathbf{1}}, \qquad \sigma_{\rm GMV}^2 = \frac{1}{A_0}.
$$

## 1.4 分散投資が機能する仕組み（直観）

二資産の例：$\mathrm{Var}(R_w) = w_1^2 \sigma_1^2 + w_2^2 \sigma_2^2 + 2 w_1 w_2 \rho \sigma_1 \sigma_2$。相関係数 $\rho < 1$ である限り、適切な $w$ により単一資産よりリスクを低減できる。極端例として $\rho = -1$ ならば完全ヘッジが可能（リスクゼロ）。一般 $n$ 資産では、共分散行列の **小さい固有値方向** に重みを配分するほど分散が減少する。

## 1.5 平均分散最適化の解釈

定理 1.3 の解は **二つの基底ポートフォリオの線形結合** として書ける：
$$
w^\star(\mu_P) = (1 - \alpha) w^{\rm GMV} + \alpha\, w^{\rm zero},
$$
ここで $w^{\rm zero}$ は第二の基底（後章の二基金分離と直結）。これにより「効率的ポートフォリオ全体が二次元アフィン空間を成す」という鋭い構造定理が得られる。詳細は第3章。

## 1.6 制約付き問題：ショート禁止など

実務では $w \ge 0$ や $w_i \le \bar w_i$ などの不等式制約が課される。この場合、解は KKT 条件を満たすが **陽な閉形では書けず**、二次計画ソルバ（active set, interior point）で数値的に解く。本書では理論を見通しよく示すため、以降 $w \in \mathbb{R}^n$（ショート可）を基本とし、不等式制約は実装上の注意として扱う。

## 1.7 推定誤差と頑健化

実務上、$\mu, \Sigma$ は標本平均・標本共分散で推定する。Markowitz の最適化は **入力誤差を増幅する**（"error maximization" 現象）。代表的な対処：

- **収縮推定（shrinkage）**：Ledoit–Wolf (2004) による $\hat\Sigma$ の収縮。
- **正則化目的関数**：$L^1$/$L^2$ ペナルティ追加（Brodie et al. 2009）。
- **Bayes 的事前情報の取込み**：第9章の **Black–Litterman**。
- **頑健最適化**：$\mu$ が楕円的不確実性集合に属するときの min–max。

これらは第8章・第9章で再訪する。

## 1.8 本章のまとめ

- 平均分散効率ポートフォリオは $\Sigma^{-1}\mathbf{1}$, $\Sigma^{-1}\mu$ の二つだけから合成される。
- 最小分散は目標リターンの二次関数（双曲線、第2章で詳述）。
- 推定誤差問題が実務上の核心的課題である。

---
[← 第0章](00_preliminaries.md) ｜ [次章 → 第2章 効率的フロンティア](02_efficient_frontier.md)
