# Investment Science — マルコヴィッツから連続時間最適化、デリバティブまで

学部 4 年生（卒研着手レベル）を主な対象とした、日本語の投資理論教科書プロジェクト。Markowitz の平均分散分析から、CAPM・APT・効用理論・連続時間最適化・現代理論（縮小推定／ロバスト最適化／リスクパリティ／因子動物園／バックテスト過剰適合）・デリバティブ（無裁定／二項モデル／BSM・Greeks）まで、**直観・数値例・図** を中心に展開する。

最新版は **[v6/](v6/)**（全 20 章、100 頁、PDF [v6/book.pdf](v6/book.pdf)）。

## ハイライト

- **全章に証明スケッチ／段階的導出を併記**。学部生が手で再現できる粒度。
- **手計算による数値例** を多数収録（GMV、接線ポートフォリオ、β、BSM、Greeks、二項ツリー後ろ向き帰納など）。
- **40 枚以上の matplotlib 図**（payoff 図、双曲線フロンティア、SML、Marchenko–Pastur、IV smile、デルタヘッジ P&L 分解、etc.）。
- **4 種類のサイドボックス**：🎯 直観／📐 数学補足／⚠️ 実務注意／📜 歴史 note。
- **各章末「上級学習への道標」** で測度論・Itô 解析・凸双対などの先につながる文献を案内。
- **codex との critical review** によって、HRP の OOS 優位主張、FF5 の HML 冗長性、DSR の $\sqrt{2 \log N}$ 近似、リスク中立確率の誤読 など、入門書では省かれがちな caveat を明示。

## バージョン履歴

このリポジトリは **バージョンごとに `vN/` を作る** 方針で管理しており、各版を独立に再ビルドできる。

| ver | 頁 | 焦点 |
|---|---:|---|
| [v1](v1/) | 38 | 院級ベースライン — Markowitz, CAPM, APT, 効用, リスク尺度, BL, Merton, パフォーマンス |
| [v2](v2/) | 45 | + matplotlib 図 16 枚（Hiragino フォントで日本語対応） |
| [v3](v3/) | 67 | + 現代理論 5 章（Ledoit–Wolf 縮小、ロバスト/DRO、ERC/HRP、ファクター動物園・ML、バックテスト過剰適合）。codex critical review で caveat を強化 |
| [v4](v4/) | 72 | 学部 4 年向けに難易度調整 — 測度論・Itô・凸双対を本文から外し、各章末の「上級学習への道標」へ。サイドボックス導入 |
| [v5](v5/) | 88 | + デリバティブ 3 章（無裁定とパリティ、二項モデル、BSM・Greeks）。章番号ダブり問題を修正 |
| **[v6](v6/)** | **100** | **+ 主要章の大幅加筆**（手計算例・段階的導出・新図 9 枚） |

## 章構成（v6）

| 章 | 主題 |
|---|---|
| 00 to_reader | 対象読者・ボックス凡例・章ごとの読み進め方 |
| 00 preliminaries | 線形代数・確率・ラグランジュ法の復習 |
| 01 markowitz_basics | Markowitz 平均分散分析、3 資産手計算、推定誤差の実演 |
| 02 efficient_frontier | 双曲線フロンティア、ゼロ β、接線ポートフォリオ、CML |
| 03 two_fund_separation | Tobin 二基金分離 |
| 04 capm | CAPM、SML、β 手計算、Black ゼロ β、Roll 批判 |
| 05 apt_factor_models | APT、Fama–French、Carhart |
| 06 utility_theory | vNM 期待効用、Arrow–Pratt、CARA/CRRA |
| 07 stochastic_dominance | FSD/SSD、Rothschild–Stiglitz |
| 08 risk_measures | VaR・CVaR・コヒーレント性 |
| 09 black_litterman | ベイズ事後ポートフォリオ |
| 10 continuous_time_merton | Merton 連続時間（1 期間極限導出） |
| 11 performance_evaluation | Sharpe・Sortino・Jensen α・IR・MDD |
| 12 shrinkage_rmt | Marchenko–Pastur、Ledoit–Wolf |
| 13 robust_optimization | Goldfarb–Iyengar、Garlappi–Uppal–Wang、DRO |
| 14 risk_parity_hrp | ERC、HRP |
| 15 factor_zoo_ml | FF5、q-factor、Harvey–Liu–Zhu、Gu–Kelly–Xiu |
| 16 backtest_overfit | Bailey–López de Prado DSR、CSCV |
| 17 derivatives_intro | 無裁定、フォワード、プット・コール・パリティ |
| 18 binomial_model | CRR、複製ポートフォリオ、リスク中立確率、BSM 連続極限 |
| 19 black_scholes | BSM 公式の積分導出、Greeks、IV smile、デルタヘッジ P&L 分解 |
| 20 exercises | 3 カテゴリ × 3 難度の演習問題 |

## ビルド方法

各 `vN/` ディレクトリは自己完結している。図とPDFを再生成するには：

```sh
# 1. Python 仮想環境のセットアップ（リポジトリ直下）
python3 -m venv .venv
.venv/bin/pip install numpy scipy matplotlib

# 2. 図の生成（vN/ に入って）
cd v6
../.venv/bin/python figures/gen_figures.py

# 3. Markdown → LaTeX → PDF（lualatex 必須、目次のため 2 回回す）
python3 build.py
lualatex -interaction=nonstopmode book.tex
lualatex -interaction=nonstopmode book.tex
```

詳しい開発ガイド（pandoc を使わない理由、Markdown パーサの癖、matplotlib mathtext の落とし穴など）は [CLAUDE.md](CLAUDE.md) を参照。

## 演習問題の構造

各章末（および巻末の第17 / 20 章）の演習は 3 カテゴリで構成：

- **★ 計算問題**（手で解ける） — 必修
- **★★ シミュレーション問題**（Python） — 選択（プログラミング既習者）
- **★★★ 考察問題**（言葉で答える） — 卒研テーマ選定の練習、評価軸付き

## 編集上の方針

- **Markdown ファイルが正本**。`book.tex` / `book.pdf` は生成物。
- **章ごとの順序は `build.py` 内の `CHAPTERS` リスト** で管理。ファイル名のソート順ではない。
- **過去版は freeze**。改訂は `cp -r v6 v7` から始める。
- **作画は `gen_figures.py` の関数として実装**。Hiragino フォントが必要（macOS 想定）。

## 引用方針

本書は教育目的の二次資料。主要な一次文献（Markowitz 1952, Sharpe 1964, Ross 1976, Merton 1969/1971, Black–Scholes 1973, Artzner et al. 1999, Ledoit–Wolf 2004, ほか）は本文中に明示。詳細書誌は各章末のセクションを参照。

## ライセンス

未指定。個人学習・研究目的での参照を想定。

---

🤖 codex と Claude Code Opus 4.7 (1M context) の共同制作。codex には Round 1–3 の critical review で各章の caveat 発掘を担当してもらった（詳細は CLAUDE.md「Codex との協働プロトコル」参照）。
