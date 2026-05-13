# Investment Science — マルコヴィッツから連続時間最適化まで

> **Version 3 — 現代理論拡張版**：v2 の可視化に加え、2010 年以降の標準的発展（縮小推定、ロバスト最適化、リスクパリティ・HRP、ファクター動物園・機械学習資産価格、バックテスト過剰適合）を 5 章追加。codex との 3 ラウンドの critical review を経て、各主張の caveat を明示。

## 本書の方針

本書は **平均分散分析（Markowitz, 1952）** の最小限の準備から出発し、現代投資理論の中核をなす結果を **証明（あるいは証明方針）付き** で展開する大学院初級レベルの教科書である。読者は

- 線形代数（固有値・二次形式・行列の正定値性）
- 多変量確率論（期待値・共分散行列・正規分布）
- 凸最適化と Lagrange 乗数法
- 確率過程と Itô 解析の入門（第10章以降）

を既習であることを想定する。Markowitz 理論の概略を知っている読者は第1章を流し読みし、第2章の解析的フロンティア導出から本格的に読み始めると良い。

## 目次

| 章 | 主題 | 主な結果 |
|---|------|---------|
| [00](00_preliminaries.md) | 記号と準備 | 線形代数・確率の道具 |
| [01](01_markowitz_basics.md) | Markowitz の平均分散分析 | 効率的ポートフォリオ、最小分散 |
| [02](02_efficient_frontier.md) | 効率的フロンティアの解析的導出 | 双曲線・接線ポートフォリオ |
| [03](03_two_fund_separation.md) | 二基金分離定理と無リスク資産 | Tobin の分離定理、シャープ最大化 |
| [04](04_capm.md) | CAPM | 証券市場線、β、Roll 批判 |
| [05](05_apt_factor_models.md) | APT とファクターモデル | Ross の裁定論証、Fama–French |
| [06](06_utility_theory.md) | 効用理論と期待効用 | von Neumann–Morgenstern、絶対/相対リスク回避 |
| [07](07_stochastic_dominance.md) | 確率支配 | FSD・SSD と効用の同値性 |
| [08](08_risk_measures.md) | リスク尺度 | VaR、CVaR、コヒーレント性 |
| [09](09_black_litterman.md) | Black–Litterman モデル | ベイズ的事後ポートフォリオ |
| [10](10_continuous_time_merton.md) | 連続時間最適化（Merton） | HJB 方程式、最適消費・投資 |
| [11](11_performance_evaluation.md) | 運用パフォーマンス評価 | Sharpe・Sortino・Jensen α・情報比 |
| [12](12_shrinkage_rmt.md) | 共分散の縮小推定と高次元漸近 | Marčenko–Pastur、Ledoit–Wolf 線形/非線形 |
| [13](13_robust_optimization.md) | ロバスト/分布ロバスト最適化 | Goldfarb–Iyengar、Garlappi–Uppal–Wang、DRO |
| [14](14_risk_parity_hrp.md) | リスクパリティと HRP | ERC、HRP、Maillard–Roncalli–Teiletche |
| [15](15_factor_zoo_ml.md) | ファクター動物園・機械学習資産価格 | FF5、q-factor、Harvey–Liu–Zhu、Gu–Kelly–Xiu |
| [16](16_backtest_overfit.md) | バックテスト過剰適合・DSR | Bailey–López de Prado、PBO、CSCV |
| [17](17_exercises.md) | 演習問題 | 各章の補強 |

## 表記の約束

- ベクトルは縦ベクトル、転置は $x^\top$。
- $\mathbf{1} = (1, \dots, 1)^\top \in \mathbb{R}^n$。
- 期待リターン $\mu \in \mathbb{R}^n$、共分散行列 $\Sigma \in \mathbb{R}^{n\times n}$（特に断らない限り正定値）。
- ポートフォリオ重み $w \in \mathbb{R}^n$、$\mathbf{1}^\top w = 1$（完全投資）。
- 無リスク利子率 $r_f$。
- 確率変数 $X$ の期待値 $\mathbb{E}[X]$、分散 $\mathrm{Var}(X)$。

## 引用慣行

本文中で挙げる古典文献は Markowitz (1952, *JF*), Tobin (1958), Sharpe (1964), Lintner (1965), Mossin (1966), Ross (1976), Merton (1969, 1971), Black–Litterman (1992), Artzner et al. (1999) などである。詳細書誌は各章末を参照のこと。
