# Investment Science — マルコヴィッツから連続時間最適化まで

> **Version 4 — 学部 4 年向けに難易度調整版**：v3 の現代理論を保ちながら、対象を学部 4 年（卒研着手レベル）に下げて書き直した。測度論・Itô 解析・凸双対などは本文から外し、各章末の「上級学習への道標」に整理。直観・数値例・図を増強し、4 種類のサイドボックス（直観／数学補足／実務注意／歴史 note）で読みやすさを改善。

## 本書の方針

本書は **平均分散分析（Markowitz, 1952）** から出発し、現代投資理論の中核をなす結果を **学部 4 年向けの直観と数値例** を中心に展開する教科書である。詳しい対象読者・前提知識・ボックスの使い方は [読者へ](00_to_reader.md) を参照。

想定する前提知識（最小限）：
- 線形代数（行列演算、固有値、対称行列の対角化）
- 微積分（多変数偏微分、Taylor 展開）
- 確率統計（期待値、分散、共分散、正規分布、大数の法則）
- ラグランジュ未定乗数法

仮定しない：測度論、Itô 解析、凸双対、SDP/SOCP の厳密理論（各章末で名前と道標だけ示す）。

## 目次

| 章 | 主題 | 主な結果 |
|---|------|---------|
| [読者へ](00_to_reader.md) | 対象読者・本書の使い方 | ボックス凡例、章ごとの読み進め方 |
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
