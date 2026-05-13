# 第11章 運用パフォーマンスの評価

ポートフォリオ運用の **事後評価** は、リスク調整リターンを定量化することで「運用者は実力か運か」「ベンチマークを超えたか」を問う活動である。本章は古典的指標と統計的検出力を解説する。

## 11.1 Sharpe 比

**定義 11.1**（Sharpe 1966）
$$
\mathrm{SR} = \frac{\mathbb{E}[R_p] - r_f}{\sigma_p}.
$$

CML の傾きそのものであり、第2章定理 2.3 の最大化対象。

**統計的推定**：標本サイズ $T$、独立 i.i.d 仮定下で
$$
\hat{\mathrm{SR}} - \mathrm{SR} \sim \mathcal{N}\left(0,\ \frac{1 + \mathrm{SR}^2/2}{T}\right) \quad \text{(漸近)}.
$$
*証明方針*：デルタ法 + 標本平均・標本分散の漸近正規性。$\square$

**注意**：リターンが **正規** という仮定が標準。歪度・尖度がある場合は Sharpe 比が不適切なことがある（Mertens 2002 補正、Lo 2002）。

## 11.2 Treynor 比

$$
\mathrm{Treynor} = \frac{\mathbb{E}[R_p] - r_f}{\beta_p}.
$$
全分散の代わりに **β（市場リスク）** で割る。well-diversified なポートフォリオには適切だが、固有リスクが大きい場合に過大評価しがち。

## 11.3 Jensen の α

CAPM 回帰
$$
R_{pt} - r_{ft} = \alpha_p + \beta_p (R_{Mt} - r_{ft}) + \varepsilon_{pt}
$$
の切片 $\alpha_p$。CAPM が正しければ均衡で $\alpha = 0$ なので、$\alpha > 0$ は **「市場では説明できない超過リターン」**＝運用力（skill）の指標とされる。

**多因子 α**：Fama–French 3 因子回帰の切片で評価することで、サイズ・バリュー因子で説明される α を除外できる。

## 11.4 情報比率（Information Ratio; IR）

ベンチマーク $R_B$ に対する **アクティブリターン** $R_p - R_B$ について
$$
\mathrm{IR} = \frac{\mathbb{E}[R_p - R_B]}{\mathrm{SD}(R_p - R_B)} = \frac{\alpha}{\sigma_\varepsilon}.
$$
**トラッキングエラー** で標準化された超過リターン。アクティブ運用の評価指標として最も重視される。

## 11.5 Sortino 比

下方リスクのみを罰する指標：
$$
\mathrm{Sortino} = \frac{\mathbb{E}[R_p] - r_f}{\mathrm{SD}^-(R_p)},
$$
ここで $\mathrm{SD}^-$ は下方半分散の平方根。投資家がアップサイドのボラを嫌わないなら理論的に Sharpe より整合的。

## 11.6 最大ドローダウン・Calmar 比

**最大ドローダウン**（MDD）：累積価値のピークから谷底までの最大下落率
$$
\mathrm{MDD} = \sup_{0 \le s \le t \le T} \left(1 - \frac{V_t}{V_s}\right).
$$
**Calmar 比** $= \mathbb{E}[R] / \mathrm{MDD}$（典型的に年率）。

## 11.7 リスク調整パフォーマンスの数学的洞察

**事実 11.2**
任意の運用戦略について、Sharpe 比は **時間にスケールする**：
$$
\mathrm{SR}_{\text{年率}} = \sqrt{n} \cdot \mathrm{SR}_{\text{日次}}
$$
（独立同分布リターン仮定）。年率化は単純に $\sqrt{T}$ 倍。

**事実 11.3**（Bailey–López de Prado 2014）
バックテストにおいて、複数の戦略から最良を選んだ「事後選択 Sharpe」は **大幅にバイアスがかかる**：
$$
\mathbb{E}[\max_i \widehat{\mathrm{SR}}_i] \approx \sigma_{\mathrm{SR}} \cdot \sqrt{2\log N}
$$
（$N$ は戦略数）。**Deflated Sharpe Ratio** で補正が必要。

## 11.8 検出力：必要サンプル

「ある真の Sharpe 比が $\mathrm{SR}$ で正であるか有意性 5%・検出力 80% で確認するための最低期間」は概算
$$
T \gtrsim \frac{(1.96 + 0.84)^2 (1 + \mathrm{SR}^2/2)}{\mathrm{SR}^2}.
$$
たとえば $\mathrm{SR} = 0.5$ 年なら約 32 年必要。**多くのファンドが「実力か運か」統計的に判別不能**。

## 11.9 リスク帰属・パフォーマンス帰属

- **Brinson 分解**：超過リターンを「アロケーション効果・銘柄選択効果・相互作用効果」に分解。
- **リスク要因帰属**：ファクターモデルで分散を要因別に分解。
- **時系列回帰 α + 横断面回帰**：Fama–MacBeth 型評価。

## 11.10 リスク調整指標の関係

| 指標 | リスク尺度 | ベンチマーク | 仮定 |
|------|-----------|-------------|------|
| Sharpe | 全分散 | 無リスク | 正規分布 |
| Treynor | β | 無リスク | 完全分散 |
| Jensen α | 残差分散 | CAPM | CAPM 妥当性 |
| IR | TE（残差分散） | 指数 | — |
| Sortino | 下方分散 | 目標 | 下方リスク選好 |

## 11.11 まとめ

- Sharpe 比は最も普及した指標だが、正規分布・i.i.d 仮定に依存。
- Jensen α と IR は α 検出に有用。
- バックテスト Sharpe は **データスヌーピング** で著しくバイアス。
- 統計的に skill を実証するには **長期データ＋多重比較補正** が必須。

---
[← 第10章](10_continuous_time_merton.md) ｜ [次章 → 第12章 演習問題](12_exercises.md)
