# 第10章 連続時間最適化：Merton 問題

Markowitz は単期的な静的問題を扱う。実際には投資家は **時間を通じて連続的に再配分** し、しかも **消費** も行う。Merton (1969, 1971) はこの問題を確率制御として定式化し、CRRA 効用下では閉形解を得た。

## 10.1 市場モデル

無リスク資産：$dB_t = r B_t \, dt$

危険資産（幾何ブラウン運動）：
$$
dS_t = \mu S_t \, dt + \sigma S_t \, dW_t
$$
（一般化して $n$ 資産：$dS_t^i / S_t^i = \mu_i \, dt + \sum_j \sigma_{ij} dW_t^j$）。

投資家は時点 $t$ に富 $X_t$ のうち $\pi_t$ の割合を危険資産に投資、$c_t$ の率で消費する。

**富の動学**（自己融資条件）：
$$
dX_t = \bigl[(r + \pi_t (\mu - r)) X_t - c_t\bigr] dt + \pi_t \sigma X_t \, dW_t.
$$

## 10.2 最適化問題

投資家は無限期間の総割引期待効用を最大化：
$$
V(x) = \sup_{\pi, c} \mathbb{E}\left[\int_0^\infty e^{-\rho t} u(c_t) \, dt \;\middle|\; X_0 = x\right]
$$
（$\rho > 0$ は主観的割引率）。あるいは有限期間 $T$ の場合は
$$
V(t, x) = \sup_{\pi, c}\mathbb{E}\left[\int_t^T e^{-\rho(s-t)} u(c_s) ds + e^{-\rho(T-t)} U(X_T)\;\middle|\; X_t = x\right].
$$

## 10.3 動的計画原理と HJB 方程式

**動的計画原理**：
$$
V(t, x) = \sup_{\pi, c} \mathbb{E}\left[\int_t^{t+h} e^{-\rho(s-t)} u(c_s) ds + e^{-\rho h} V(t+h, X_{t+h})\;\middle|\; X_t = x\right].
$$

$h \to 0$ の極限（および $V$ の十分な滑らかさ）から、**Hamilton–Jacobi–Bellman (HJB) 方程式**：

$$
\boxed{\;
0 = \sup_{\pi, c}\Bigl\{ u(c) - \rho V + V_t + [(r + \pi(\mu - r)) x - c]\, V_x + \tfrac{1}{2} \pi^2 \sigma^2 x^2 V_{xx}\Bigr\}
\;}
$$

$\pi, c$ について一階条件：
- 消費：$u'(c) = V_x$（**包絡定理**）。
- ポートフォリオ：$(\mu - r) x V_x + \pi \sigma^2 x^2 V_{xx} = 0$、すなわち
$$
\pi^\star = -\frac{(\mu - r) V_x}{\sigma^2 x V_{xx}}.
$$

## 10.4 CRRA 効用での閉形解

$u(c) = c^{1-\gamma}/(1-\gamma)$（$\gamma > 0, \gamma \neq 1$）の場合、**値関数の同次性** から
$$
V(t, x) = h(t) \cdot \frac{x^{1-\gamma}}{1-\gamma}
$$
と推測する。

**定理 10.1**（Merton 1969, 無限期間）
$\rho > (1-\gamma)\bigl[r + (\mu-r)^2/(2\gamma\sigma^2)\bigr]$ の下で、最適ポートフォリオは
$$
\boxed{\;
\pi^\star = \frac{\mu - r}{\gamma \sigma^2}
\;}
$$
（時間不変な定数）、最適消費は富に比例
$$
c^\star_t = m \cdot X_t,\qquad m = \frac{1}{\gamma}\left[\rho - (1-\gamma)\left(r + \frac{(\mu-r)^2}{2\gamma\sigma^2}\right)\right].
$$

*証明方針*：仮定 $V = h \cdot x^{1-\gamma}/(1-\gamma)$ を HJB に代入し、$x^{1-\gamma}$ 部分を消去、$h(t)$ ないし定常解 $h$ が満たす ODE/代数方程式を解く。CRRA 同次性が「富比例」の最適消費・投資を導く。$\square$

**多資産化**：$\mu \in \mathbb{R}^n$、共分散 $\Sigma$ で $\pi^\star = \tfrac{1}{\gamma} \Sigma^{-1}(\mu - r\mathbf{1})$。これは **静的 Markowitz 接線ポートフォリオの動的版** であり、配分比率は時間不変（**myopic policy**）。

![Merton 問題の最適投資比率 $\pi^\star = (\mu-r)/(\gamma\sigma^2)$ を $(\mu, \sigma)$ 平面で等高線表示。ボラティリティが高いほど・期待リターンが低いほど危険資産投資を絞る。](figures/merton_policy.png)

![Merton 問題の最適制御下での富経路と消費経路のシミュレーション（8 軌跡、$\mu=0.08, \sigma=0.18, r=0.02, \gamma=3$）。富は対数線形成長、消費は富に比例する。](figures/merton_paths.png)

## 10.5 確率係数モデルと「ヘッジ需要」

市場係数（$\mu, \sigma, r$）が **状態変数** $Y_t$（金利、ボラティリティ、配当利回り等）に依存し $Y_t$ が独立な確率過程に従う場合、最適投資は二項に分解：

$$
\pi^\star_t = \underbrace{\frac{\mu(Y_t) - r(Y_t)}{\gamma \sigma(Y_t)^2}}_{\text{myopic}} + \underbrace{\frac{V_{xY}}{x V_{xx}} \cdot \frac{(\rho_{SY} \sigma_Y)}{\sigma(Y_t)}}_{\text{ヘッジ需要}}
$$

ヘッジ需要は **状態変動に対する効用の感度** を反映し、Merton (1973) の **異時点間 CAPM (ICAPM)** の起源となる。実証的には金利・ボラティリティのヘッジ需要が観察される。

## 10.6 マーチンゲール法（Cox–Huang 1989, Karatzas–Lehoczky–Shreve 1987）

完全市場では、HJB を解く代わりに **状態価格密度（state-price density）** $\xi_t$ を用い、富の動的予算制約を **単一の静的予算制約** に書き換える：
$$
\mathbb{E}\left[\int_0^T \xi_t c_t \, dt + \xi_T X_T\right] \le X_0.
$$
ラグランジュ乗数法で **凸双対** を取り、消費と最終富を $\xi$ の関数として陽に得る。これは確率制御を **凸最適化** に帰着する強力な技法である。

## 10.7 摩擦・制約のある拡張

- **取引費用**：Davis–Norman (1990), Shreve–Soner (1994)。最適政策は「無取引帯（no-trade region）」を持つ。
- **借入制約・空売り禁止**：自由境界問題。
- **不完全市場**：BSDE / 二重凸双対 (Cvitanic–Karatzas)。
- **ジャンプ過程**：Merton (1976) のジャンプ拡散モデル。

## 10.8 連続時間 CAPM・ICAPM

Merton (1973) は **異時点間 CAPM (ICAPM)** を提示：
$$
\mathbb{E}[R_i] - r = \beta_{iM}\, (\mathbb{E}[R_M] - r) + \sum_k \beta_{ik}^Y\, \lambda_k^Y.
$$
通常の市場 β に加え、**状態変数のヘッジリスクプレミアム** が現れる。これは Fama–French 等の経験的多因子モデルへの理論的橋渡しとなる。

## 10.9 確率割引因子（SDF）と一般枠組み

最適化の一階条件は、ある正の確率過程 $M_t$（SDF）について
$$
\mathbb{E}_t[M_{t+\Delta} R_{i,t+\Delta}] = 1
$$
と書ける。SDF は CCAPM（消費 CAPM、第4章末で言及）、Merton ICAPM、APT 等すべての資産価格モデルを統一する。

## 10.10 まとめ

- Merton 問題は HJB あるいはマーチンゲール法で解ける。
- CRRA 効用下では **「Markowitz 接線比率を持ち続け、富比例で消費」** が最適。
- 状態変数依存の係数では **ヘッジ需要** が現れ、ICAPM の多因子構造を生む。
- 取引費用・制約・ジャンプ等の拡張は活発な研究領域。

---
[← 第9章](09_black_litterman.md) ｜ [次章 → 第11章 パフォーマンス評価](11_performance_evaluation.md)
