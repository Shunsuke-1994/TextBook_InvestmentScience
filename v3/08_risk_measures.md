# 第8章 リスク尺度

## 8.1 リスク尺度とは何か

リスク尺度 $\rho$ は損失分布 $L$ に実数 $\rho(L)$ を対応させる写像。標準偏差はその一例だが、

- 損失の **左裾（テールリスク）** に十分敏感でない
- 対称な指標であり「損益の非対称性」を捉えない
- 規制目的では **必要資本量** の解釈が直接できない

これらの欠点を補うのが本章のテーマである。

## 8.2 Value-at-Risk (VaR)

**定義 8.1**
損失 $L$ の信頼水準 $\alpha \in (0, 1)$ における VaR は
$$
\mathrm{VaR}_\alpha(L) := \inf\{ x \in \mathbb{R} : P(L > x) \le 1 - \alpha\} = F_L^{-1}(\alpha)
$$
通常 $\alpha = 0.95$ や $0.99$。

**バーゼル規制等で長年標準だった**。しかし以下の欠点を持つ。

![VaR と CVaR の幾何的意味（損失が標準正規の場合）。VaR は分位点で「テールの境界」だけを見るのに対し、CVaR は分位点を超えた領域の **平均** を取るため、テール形状の情報を保持する。](figures/var_cvar.png)

## 8.3 VaR の問題点：非劣加法性

**事実 8.2**
VaR は一般に **劣加法性**（$\rho(L_1 + L_2) \le \rho(L_1) + \rho(L_2)$）を満たさない。すなわち分散投資により VaR が悪化することがある。

**反例**：独立同分布の二つのデフォルト債を考える。各々 $P(\text{損失}=100) = 0.04$、$P(0) = 0.96$。$\alpha = 0.95$ で個別 VaR はいずれも 0。しかし合算 $L_1 + L_2$ の損失分布は $P(L=200) = 0.0016$, $P(L=100) = 2 \cdot 0.04 \cdot 0.96 = 0.0768$, $P(L=0) = 0.9216$。よって $\mathrm{VaR}_{0.95}(L_1 + L_2) = 100 > 0 + 0 = \mathrm{VaR}(L_1) + \mathrm{VaR}(L_2)$。

これは「分散投資をペナルティする」異常な振る舞いで、合理的なリスク管理に反する。

![VaR の劣加法性違反例（独立 2 デフォルト債）。左：個別の損失分布、$\mathrm{VaR}_{0.95}=0$。右：合算した分布で $\mathrm{VaR}_{0.95}=100$ となり、$0+0$ を超える。分散投資が VaR を悪化させる典型例。](figures/var_subadditivity.png)

## 8.4 コヒーレントリスク尺度（Artzner et al. 1999）

**定義 8.3**（コヒーレントリスク尺度）
$\rho$ が次の四公理を満たすとき：
1. **単調性**：$L_1 \le L_2$ a.s. $\Rightarrow \rho(L_1) \le \rho(L_2)$
2. **平行移動不変性**：$\rho(L + c) = \rho(L) + c$（$c$ は確定損失）
3. **正同次性**：$\rho(\lambda L) = \lambda \rho(L)$（$\lambda \ge 0$）
4. **劣加法性**：$\rho(L_1 + L_2) \le \rho(L_1) + \rho(L_2)$

**定理 8.4**（双対表現）
$\rho$ がコヒーレントである必要十分条件は、ある確率測度の凸集合 $\mathcal{Q}$ が存在して
$$
\rho(L) = \sup_{Q \in \mathcal{Q}} \mathbb{E}_Q[L].
$$

*証明方針*：劣加法性と正同次性から $\rho$ は **凸関数（正同次凸 = 部分線形）**、Hahn–Banach の凸双対と Riesz の表現定理により、$\rho$ は確率測度の凸集合上の sup として書ける。$\square$

経済的解釈：コヒーレントリスク尺度は「最悪のシナリオ家族 $\mathcal{Q}$ における最大期待損失」。

## 8.5 Conditional VaR / Expected Shortfall

**定義 8.5**
$$
\mathrm{CVaR}_\alpha(L) = \mathbb{E}[L \mid L \ge \mathrm{VaR}_\alpha(L)]
$$
（厳密には離散的分布の場合の取扱いに注意、Rockafellar–Uryasev 2000 を参照）。同義：Expected Shortfall (ES)、Tail VaR、Average VaR。

**定理 8.6**
CVaR は **コヒーレント** リスク尺度。

*証明方針*：双対表現
$$
\mathrm{CVaR}_\alpha(L) = \sup\bigl\{ \mathbb{E}_Q[L] : Q \ll P,\ \tfrac{dQ}{dP} \le \tfrac{1}{1-\alpha}\bigr\}
$$
を直接確かめ、定理 8.4 を適用。$\square$

**実務的利点**：
- 劣加法性のため分散投資を促進する。
- 凸最適化と相性が良い（Rockafellar–Uryasev による線形計画化）。
- 損失の **平均的な大きさ** を測るため、テールリスクの情報を保つ。

バーゼル III から CVaR が市場リスク資本計算の標準指標となった。

## 8.6 CVaR 最適化（Rockafellar–Uryasev）

**定理 8.7**
ポートフォリオ $w$ の損失 $L_w(R) = -w^\top R$ の CVaR は補助関数
$$
F_\alpha(w, \zeta) := \zeta + \frac{1}{1-\alpha}\, \mathbb{E}[(L_w - \zeta)^+]
$$
について $\mathrm{CVaR}_\alpha(L_w) = \min_\zeta F_\alpha(w, \zeta)$。さらに $(w, \zeta)$ について同時最小化が可能で **凸最適化**（リターン分布が標本表現なら **線形計画問題**）に帰着する。

*証明方針*：$\zeta = \mathrm{VaR}_\alpha$ で最小達成、これは凸関数の劣勾配法で示せる。標本表現 $R \in \{R^{(1)}, \dots, R^{(S)}\}$ では $(L_w^{(s)} - \zeta)^+ = u_s$、$u_s \ge 0$、$u_s \ge L_w^{(s)} - \zeta$ の線形緩和で書ける。$\square$

これにより、**シナリオベース CVaR 最小化** は LP ソルバで効率的に解ける。

## 8.7 スペクトラルリスク尺度・歪曲リスク尺度

CVaR を一般化した **スペクトラルリスク尺度**：
$$
\rho_\phi(L) = \int_0^1 \phi(p) \mathrm{VaR}_p(L) \, dp
$$
（$\phi$ は確率密度型の重み関数）。

**事実 8.8**（Acerbi 2002）
$\phi$ が **非減少** ならば $\rho_\phi$ はコヒーレント。

これにより「テールに重みを大きく置く」リスク尺度を柔軟に設計できる。

## 8.8 凸リスク尺度

**定義 8.9**（Föllmer–Schied 2002）
劣加法性 + 正同次性を **凸性**：$\rho(\lambda L_1 + (1-\lambda) L_2) \le \lambda\rho(L_1) + (1-\lambda)\rho(L_2)$ に弱めたものを **凸リスク尺度** と呼ぶ。これにより流動性リスクや非線形なリスク（巨大ポジションの増大効果）を扱える。

双対表現は
$$
\rho(L) = \sup_{Q} \{\, \mathbb{E}_Q[L] - \alpha(Q) \,\}
$$
（$\alpha$ は罰金関数）。

## 8.9 動的リスク尺度（時間整合性）

複数期間ではリスクの **時間整合性（time consistency）** が問題となる。一般に「現在から見た 2 期間の CVaR」と「明日条件付き CVaR の現在 CVaR の合成」は一致しない。時間整合的リスク尺度は **入れ子型条件付期待値表現** を持つ。詳細は Cheridito–Delbaen–Kupper (2006), Ruszczyński (2010) を参照。

## 8.10 まとめ

- VaR は実務標準だったが劣加法性を欠き、分散投資を罰しうる。
- CVaR は劣加法性を満たし、凸最適化で扱える。
- Artzner 公理が「合理的なリスク尺度」を特徴付ける。
- 双対表現はリスク尺度を「最悪シナリオ族」として理解する強力な枠組み。

---
[← 第7章](07_stochastic_dominance.md) ｜ [次章 → 第9章 Black–Litterman](09_black_litterman.md)
