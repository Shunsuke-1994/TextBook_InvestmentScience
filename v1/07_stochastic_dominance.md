# 第7章 確率支配

平均分散だけでは捉えきれない「選好」のうち、ほぼ全ての合理的投資家が同意する **頑健な順序付け** が確率支配である。

## 7.1 一次確率支配（First-order Stochastic Dominance; FSD）

**定義 7.1**
確率変数 $X, Y$ について $X \succeq_{\rm FSD} Y$ とは、$F_X(t) \le F_Y(t)$ がすべての $t$ で成立すること。すなわち $X$ の累積分布関数が $Y$ より「下にある」。

**定理 7.2**
$X \succeq_{\rm FSD} Y \iff \mathbb{E}[u(X)] \ge \mathbb{E}[u(Y)]$ がすべての **単調非減少関数** $u$ に対して成立。

*証明方針*：
- ($\Leftarrow$) $u(t) = \mathbf{1}_{t > s}$ という指示関数（極限的に単調）を取ると $\mathbb{E}[u(X)] = 1 - F_X(s)$、これが全 $s$ で $\ge 1 - F_Y(s)$。
- ($\Rightarrow$) $u$ が単調非減少なら $u(t) = \int \mathbf{1}_{t > s} \mu(ds)$（Stieltjes 積分表示）で書け、各 indicator で成立する不等式を線形結合。$\square$

すなわち FSD は「より多くを好む」全ての投資家が同意する順序。

## 7.2 二次確率支配（Second-order Stochastic Dominance; SSD）

**定義 7.3**
$X \succeq_{\rm SSD} Y$ とは $\int_{-\infty}^t F_X(s) ds \le \int_{-\infty}^t F_Y(s) ds$ が全 $t$ で成立すること。

**定理 7.4**
$X \succeq_{\rm SSD} Y \iff \mathbb{E}[u(X)] \ge \mathbb{E}[u(Y)]$ がすべての **単調非減少凹関数** $u$ に対して成立。

*証明方針*：
- ($\Leftarrow$) $u(t) = -(s - t)^+$（凹かつ単調非減少）に対する期待値が CDF の積分で書けることを示す（部分積分）。
- ($\Rightarrow$) 任意の単調凹 $u$ を $-(s-t)^+$ 型関数の凸結合として近似（Riemann 和分解、Hardy–Littlewood–Pólya の補題）。$\square$

すなわち SSD は **「より多くを好み、かつリスク回避的な」全投資家** が同意する順序。

**重要な系**：$\mathbb{E}[X] = \mathbb{E}[Y]$ かつ $\mathrm{Var}(X) \le \mathrm{Var}(Y)$ でも、$X \succeq_{\rm SSD} Y$ は **一般には成立しない**（高次モーメントの効果）。逆に正規分布族の中では平均分散順序＝SSD 順序。

## 7.3 リスクのない/あるリスクの増大（Rothschild–Stiglitz 1970）

**定義 7.5**
$Y$ が $X$ の **平均保存スプレッド（mean-preserving spread; MPS）** であるとは、$\mathbb{E}[X] = \mathbb{E}[Y]$ かつ $Y \stackrel{d}{=} X + \varepsilon$（$\varepsilon$ は条件付き平均ゼロのノイズ）と表せること。

**定理 7.6**（Rothschild–Stiglitz）
以下は同値：
1. $Y$ は $X$ の MPS
2. $X \succeq_{\rm SSD} Y$（同一平均下で）
3. 全てのリスク回避効用 $u$ について $\mathbb{E}[u(X)] \ge \mathbb{E}[u(Y)]$

*証明方針*：MPS の構成性（独立ノイズ重ね合わせ）と Jensen の不等式の繰り返し適用、および「条件付き平均ゼロのノイズ追加」が Hardy–Littlewood–Pólya 型の majorization を生むこと。$\square$

## 7.4 平均分散効率性と SSD の関係

**事実 7.7**
- リターンが **多変量正規** な世界では：「平均分散効率」⇔「SSD 効率」⇔「ある CARA リスク回避者の最適化」。
- 非正規世界では：平均分散効率は SSD 効率の **必要条件ではない**。実際、平均分散効率なポートフォリオが SSD で支配されることがある（例：ポジティブ歪度資産を排除するケース）。

これが「Markowitz 最適化が常に賢いとは限らない」具体的な根拠となる。

## 7.5 高次の確率支配

**$n$ 次確率支配** は順次積分した分布関数の比較で定義され、$u, u'', u^{(4)}, \dots$ の符号についての条件を満たす全ての効用に対応する。3 次以上は **慎重度（prudence）** や **節制（temperance）** といった高次性質と関連する。

## 7.6 実用上の意義

- **退職基金管理**：SSD 効率ポートフォリオを選ぶことで、加入者の効用関数を特定せずに「合理的な誰にとっても劣らない」運用を保証できる。
- **オプション戦略**：プレミアム支払い・限定上限の戦略は平均分散では劣ることがあっても SSD では支配する/されない複雑な構造を持つ。

## 7.7 まとめ

- FSD：単調効用全体に対する順序（ほぼ普遍的）。
- SSD：単調凹効用全体に対する順序（リスク回避者全員）。
- MPS と SSD は同値（Rothschild–Stiglitz）。
- 平均分散効率 ≠ SSD 効率（非正規分布下）。

---
[← 第6章](06_utility_theory.md) ｜ [次章 → 第8章 リスク尺度](08_risk_measures.md)
