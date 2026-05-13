"""本書 v2 用の図を一括生成。

実行: ../../.venv/bin/python gen_figures.py
出力: figures/*.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from scipy.stats import norm, lognorm

HERE = Path(__file__).parent

# Register Hiragino so Japanese chars render properly
for p in [
    "/System/Library/Fonts/ヒラギノ角ゴシック W4.ttc",
    "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
]:
    if Path(p).exists():
        try:
            font_manager.fontManager.addfont(p)
        except Exception:
            pass

plt.rcParams.update({
    "figure.figsize": (6.4, 4.2),
    "figure.dpi": 130,
    "font.size": 11,
    "font.family": ["Hiragino Sans", "Hiragino Maru Gothic Pro", "DejaVu Sans"],
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.unicode_minus": False,
})


def save(name: str):
    p = HERE / f"{name}.png"
    plt.tight_layout()
    plt.savefig(p, bbox_inches="tight")
    plt.close()
    print("wrote", p.name)


# -----------------------------------------------------------------------------
# 第1章: 二資産の分散投資効果 — 相関の値で軌跡が変わる
# -----------------------------------------------------------------------------
def fig_diversification():
    s1, s2 = 0.15, 0.25
    mu1, mu2 = 0.08, 0.12
    w = np.linspace(0, 1, 200)
    fig, ax = plt.subplots()
    for rho, color in zip([-1.0, -0.3, 0.0, 0.5, 1.0], plt.cm.viridis(np.linspace(0, 0.9, 5))):
        var = w**2 * s1**2 + (1 - w) ** 2 * s2**2 + 2 * w * (1 - w) * rho * s1 * s2
        sd = np.sqrt(np.maximum(var, 0))
        ret = w * mu1 + (1 - w) * mu2
        ax.plot(sd, ret, color=color, lw=2, label=fr"$\rho={rho:+.1f}$")
    ax.scatter([s1, s2], [mu1, mu2], color="black", zorder=5)
    ax.annotate("Asset 1", (s1, mu1), xytext=(7, -3), textcoords="offset points")
    ax.annotate("Asset 2", (s2, mu2), xytext=(7, -3), textcoords="offset points")
    ax.set_xlabel(r"標準偏差 $\sigma$")
    ax.set_ylabel(r"期待リターン $\mu$")
    ax.set_title("二資産の平均-分散軌跡（相関 $\\rho$ ごと）")
    ax.legend(loc="lower right", fontsize=9, frameon=True)
    save("diversification_two_asset")


# -----------------------------------------------------------------------------
# 第2章: 効率的フロンティア & 接線ポートフォリオ & CML
# -----------------------------------------------------------------------------
def _frontier_setup():
    mu = np.array([0.06, 0.10, 0.14])
    Sigma = np.array([
        [0.0100, 0.0018, 0.0011],
        [0.0018, 0.0400, 0.0026],
        [0.0011, 0.0026, 0.0900],
    ])
    rf = 0.02
    return mu, Sigma, rf


def fig_efficient_frontier():
    mu, Sigma, rf = _frontier_setup()
    ones = np.ones_like(mu)
    inv = np.linalg.inv(Sigma)
    A = ones @ inv @ ones
    B = ones @ inv @ mu
    C = mu @ inv @ mu
    D = A * C - B * B
    mu_grid = np.linspace(0.00, 0.20, 400)
    var = (A * mu_grid**2 - 2 * B * mu_grid + C) / D
    sd = np.sqrt(np.maximum(var, 0))
    # GMV
    gmv = inv @ ones / (ones @ inv @ ones)
    gmv_mu = gmv @ mu
    gmv_sd = np.sqrt(gmv @ Sigma @ gmv)
    # Tangency
    excess = mu - rf * ones
    tan = inv @ excess
    tan = tan / (ones @ tan)
    tan_mu = tan @ mu
    tan_sd = np.sqrt(tan @ Sigma @ tan)
    sr = (tan_mu - rf) / tan_sd

    fig, ax = plt.subplots()
    eff = mu_grid >= gmv_mu
    ax.plot(sd[eff], mu_grid[eff], color="C0", lw=2.4, label="効率的フロンティア")
    ax.plot(sd[~eff], mu_grid[~eff], color="C0", lw=1.4, ls="--", alpha=0.55, label="非効率部分")
    # CML
    s_line = np.linspace(0, 0.35, 100)
    ax.plot(s_line, rf + sr * s_line, color="C3", lw=2, label=f"CML (SR={sr:.2f})")
    # individual assets
    ax.scatter(np.sqrt(np.diag(Sigma)), mu, color="black", zorder=5, label="個別資産")
    for i in range(3):
        ax.annotate(f"$R_{i+1}$", (np.sqrt(Sigma[i, i]), mu[i]),
                    xytext=(6, 4), textcoords="offset points")
    # special points
    ax.scatter([gmv_sd], [gmv_mu], color="C1", s=70, zorder=6, label="GMV")
    ax.scatter([tan_sd], [tan_mu], color="C2", s=70, zorder=6, label="接線ポートフォリオ")
    ax.scatter([0], [rf], color="C3", marker="D", s=50, zorder=6)
    ax.annotate(r"$r_f$", (0, rf), xytext=(6, -2), textcoords="offset points", color="C3")
    ax.set_xlim(-0.005, 0.35)
    ax.set_ylim(-0.005, 0.20)
    ax.set_xlabel(r"$\sigma_P$")
    ax.set_ylabel(r"$\mu_P$")
    ax.set_title("効率的フロンティアと資本市場線")
    ax.legend(loc="lower right", fontsize=9)
    save("efficient_frontier")


def fig_two_fund():
    mu, Sigma, rf = _frontier_setup()
    ones = np.ones_like(mu)
    inv = np.linalg.inv(Sigma)
    A = ones @ inv @ ones
    B = ones @ inv @ mu
    C = mu @ inv @ mu
    D = A * C - B * B
    mu_grid = np.linspace(0.04, 0.18, 200)
    sd = np.sqrt((A * mu_grid**2 - 2 * B * mu_grid + C) / D)
    fig, ax = plt.subplots()
    ax.plot(sd, mu_grid, lw=2, color="C0", label="フロンティア")
    # Pick two base funds
    for mu_t, mark, name in [(0.06, "o", r"$w^{(1)}$"), (0.14, "s", r"$w^{(2)}$")]:
        s_t = np.sqrt((A * mu_t**2 - 2 * B * mu_t + C) / D)
        ax.scatter([s_t], [mu_t], color="C3", s=80, zorder=5, marker=mark)
        ax.annotate(name, (s_t, mu_t), xytext=(8, 0), textcoords="offset points", color="C3")
    # combinations
    ts = np.linspace(-0.4, 1.4, 60)
    mus = (1 - ts) * 0.06 + ts * 0.14
    sds = np.sqrt((A * mus**2 - 2 * B * mus + C) / D)
    ax.scatter(sds, mus, color="C2", s=14, zorder=4, label="2基金線形結合")
    ax.set_xlabel(r"$\sigma_P$")
    ax.set_ylabel(r"$\mu_P$")
    ax.set_title("二基金分離: 2点で全フロンティアを生成")
    ax.legend(loc="lower right")
    save("two_fund_separation")


# -----------------------------------------------------------------------------
# 第4章: SML (Security Market Line)
# -----------------------------------------------------------------------------
def fig_sml():
    rf = 0.02
    rm = 0.09
    beta = np.linspace(-0.5, 2.2, 100)
    er = rf + beta * (rm - rf)
    fig, ax = plt.subplots()
    ax.plot(beta, er, color="C0", lw=2.2, label="SML")
    # Some example assets
    np.random.seed(3)
    betas = np.array([0.3, 0.7, 1.0, 1.3, 1.8])
    noise = np.array([0.01, -0.005, 0.0, 0.012, -0.015])
    rets = rf + betas * (rm - rf) + noise
    colors = ["C2" if n < 0 else "C3" for n in noise]
    ax.scatter(betas, rets, color=colors, s=70, zorder=5, edgecolor="black")
    for b, r, n in zip(betas, rets, noise):
        label = r"$\alpha>0$" if n > 0 else (r"$\alpha<0$" if n < 0 else r"$\alpha=0$")
        ax.annotate(label, (b, r), xytext=(7, -3), textcoords="offset points", fontsize=9)
    ax.axhline(rf, color="gray", ls=":", lw=1)
    ax.axvline(0, color="gray", ls=":", lw=1)
    ax.annotate(r"$r_f$", (0, rf), xytext=(-25, -4), textcoords="offset points")
    ax.annotate("市場ポートフォリオ", (1.0, rm), xytext=(10, -18), textcoords="offset points",
                arrowprops=dict(arrowstyle="->"))
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\mathbb{E}[R]$")
    ax.set_title("証券市場線（SML）と Jensen の $\\alpha$")
    ax.legend(loc="upper left")
    save("sml")


# -----------------------------------------------------------------------------
# 第6章: 効用関数の比較
# -----------------------------------------------------------------------------
def fig_utility():
    W = np.linspace(0.2, 4.0, 300)
    fig, ax = plt.subplots()
    # Log
    ax.plot(W, np.log(W), label=r"$\log W$ (CRRA, $\gamma=1$)", lw=2)
    # CRRA γ=2
    ax.plot(W, -1.0 / W, label=r"CRRA $\gamma=2$ : $-1/W$", lw=2)
    # CRRA γ=0.5
    ax.plot(W, 2 * (np.sqrt(W) - 1), label=r"CRRA $\gamma=0.5$", lw=2)
    # CARA
    gamma = 1.5
    cara = -np.exp(-gamma * W) / gamma
    ax.plot(W, cara - cara[0] + np.log(W[0]), label=fr"CARA $\gamma={gamma}$ (shifted)", lw=2)
    ax.axvline(1, color="gray", ls=":", lw=1)
    ax.set_xlabel(r"富 $W$")
    ax.set_ylabel(r"効用 $u(W)$")
    ax.set_title("代表的なリスク回避効用関数")
    ax.legend(loc="lower right", fontsize=9)
    save("utility_functions")


def fig_jensen():
    fig, ax = plt.subplots()
    W = np.linspace(0.4, 3.6, 200)
    u = np.sqrt(W)
    ax.plot(W, u, color="C0", lw=2.2, label=r"凹な効用 $u(W)=\sqrt{W}$")
    W1, W2 = 0.8, 3.0
    p = 0.5
    EW = p * W1 + (1 - p) * W2
    Eu = p * np.sqrt(W1) + (1 - p) * np.sqrt(W2)
    ax.plot([W1, W2], [np.sqrt(W1), np.sqrt(W2)], color="C3", lw=2)
    ax.scatter([W1, W2], [np.sqrt(W1), np.sqrt(W2)], color="C3", zorder=5)
    ax.scatter([EW], [Eu], color="C3", s=80, zorder=6)
    ax.annotate(r"$\mathbb{E}[u(W)]$", (EW, Eu), xytext=(8, -16), textcoords="offset points", color="C3")
    ax.scatter([EW], [np.sqrt(EW)], color="C2", s=80, zorder=6)
    ax.annotate(r"$u(\mathbb{E}[W])$", (EW, np.sqrt(EW)), xytext=(8, 5), textcoords="offset points", color="C2")
    ax.vlines(EW, Eu, np.sqrt(EW), color="black", ls=":", lw=1.4)
    ax.set_xlabel("W")
    ax.set_ylabel("u(W)")
    ax.set_title("Jensen の不等式：リスク回避 $\\Leftrightarrow$ 凹性")
    ax.legend(loc="lower right", fontsize=9)
    save("jensen_inequality")


# -----------------------------------------------------------------------------
# 第7章: 確率支配
# -----------------------------------------------------------------------------
def fig_stochastic_dominance():
    x = np.linspace(-2, 6, 400)
    # FSD: X dominates Y
    fx = norm.cdf(x, loc=2.5, scale=1.0)
    fy = norm.cdf(x, loc=1.5, scale=1.0)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(x, fx, label="$F_X$ (X dominates)", color="C0", lw=2)
    axes[0].plot(x, fy, label="$F_Y$", color="C3", lw=2)
    axes[0].fill_between(x, fx, fy, where=(fy >= fx), alpha=0.15, color="C2")
    axes[0].set_title(r"FSD: $F_X \leq F_Y$ for all $t$")
    axes[0].set_xlabel("t"); axes[0].set_ylabel("CDF")
    axes[0].legend()
    # SSD: same mean, X less spread
    gx = norm.cdf(x, loc=2.0, scale=0.6)
    gy = norm.cdf(x, loc=2.0, scale=1.4)
    axes[1].plot(x, gx, label="$F_X$ (low variance)", color="C0", lw=2)
    axes[1].plot(x, gy, label="$F_Y$ (high variance)", color="C3", lw=2)
    axes[1].set_title("SSD: 同一平均、$X$ がリスク小")
    axes[1].set_xlabel("t"); axes[1].set_ylabel("CDF")
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(HERE / "stochastic_dominance.png", bbox_inches="tight")
    plt.close()
    print("wrote stochastic_dominance.png")


# -----------------------------------------------------------------------------
# 第8章: VaR vs CVaR
# -----------------------------------------------------------------------------
def fig_var_cvar():
    x = np.linspace(-4, 4, 800)
    pdf = norm.pdf(x)
    alpha = 0.95
    var = norm.ppf(alpha)
    cvar = norm.pdf(var) / (1 - alpha)
    fig, ax = plt.subplots()
    ax.plot(x, pdf, color="C0", lw=2, label="損失分布 $f_L$")
    ax.fill_between(x, 0, pdf, where=(x >= var), alpha=0.35, color="C3", label=f"確率 ${1-alpha:.2f}$ 領域")
    ax.axvline(var, color="C3", lw=2, label=fr"$\mathrm{{VaR}}_{{{alpha}}}={var:.2f}$")
    ax.axvline(cvar, color="C2", lw=2, ls="--", label=fr"$\mathrm{{CVaR}}_{{{alpha}}}={cvar:.2f}$")
    ax.set_xlabel("損失 L")
    ax.set_ylabel("密度")
    ax.set_title("VaR と CVaR の幾何 (損失が標準正規)")
    ax.legend(loc="upper left")
    save("var_cvar")


def fig_var_subadditivity_violation():
    """独立 2 デフォルト債の合算で VaR が劣加法性を破る図示。"""
    p = 0.04
    # individual loss distribution: 100 w.p. p, 0 w.p. 1-p
    # Combined loss for two independent assets
    outcomes = [(0, (1 - p) ** 2), (100, 2 * p * (1 - p)), (200, p**2)]
    losses, probs = zip(*outcomes)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].bar([0, 100], [1 - p, p], width=20, color=["C2", "C3"])
    axes[0].set_xticks([0, 100])
    axes[0].set_title(f"個別デフォルト債: $P(L=100)={p}$\n$\\mathrm{{VaR}}_{{0.95}}=0$")
    axes[0].set_xlabel("損失")
    axes[0].set_ylabel("確率")
    axes[1].bar(losses, probs, width=20, color=["C2", "C1", "C3"])
    axes[1].set_xticks(losses)
    axes[1].set_title("合算ポートフォリオ\n$\\mathrm{VaR}_{0.95}=100$  $>$  $0+0$")
    axes[1].set_xlabel("損失")
    axes[1].set_ylabel("確率")
    plt.tight_layout()
    plt.savefig(HERE / "var_subadditivity.png", bbox_inches="tight")
    plt.close()
    print("wrote var_subadditivity.png")


# -----------------------------------------------------------------------------
# 第9章: Black-Litterman の事前→事後シフト
# -----------------------------------------------------------------------------
def fig_black_litterman():
    pi = np.array([0.05, 0.07, 0.09])  # equilibrium
    tau = 0.05
    Sigma = np.array([
        [0.0225, 0.005, 0.002],
        [0.005, 0.04, 0.006],
        [0.002, 0.006, 0.0625],
    ])
    P = np.array([[0.0, -1.0, 1.0]])  # asset 3 - asset 2 = 0.03
    Q = np.array([0.03])
    Omega = np.array([[0.0005]])
    inv_tau = np.linalg.inv(tau * Sigma)
    post_cov = np.linalg.inv(inv_tau + P.T @ np.linalg.inv(Omega) @ P)
    post_mu = post_cov @ (inv_tau @ pi + P.T @ np.linalg.inv(Omega) @ Q)
    labels = ["資産1", "資産2", "資産3"]
    x = np.arange(3)
    w = 0.35
    fig, ax = plt.subplots()
    ax.bar(x - w / 2, pi * 100, w, label="均衡 $\\Pi$", color="C0")
    ax.bar(x + w / 2, post_mu * 100, w, label="BL 事後 $\\hat\\mu$", color="C3")
    for i, (a, b) in enumerate(zip(pi, post_mu)):
        ax.annotate(f"{a*100:.2f}", (i - w / 2, a * 100), xytext=(0, 3),
                    textcoords="offset points", ha="center", fontsize=9)
        ax.annotate(f"{b*100:.2f}", (i + w / 2, b * 100), xytext=(0, 3),
                    textcoords="offset points", ha="center", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("期待リターン (%)")
    ax.set_title("Black–Litterman: view 「資産3 − 資産2 = +3%」反映")
    ax.legend()
    save("black_litterman")


# -----------------------------------------------------------------------------
# 第10章: Merton の最適消費・富の経路
# -----------------------------------------------------------------------------
def fig_merton_paths():
    np.random.seed(0)
    T = 30
    dt = 1 / 252
    N = int(T / dt)
    t = np.linspace(0, T, N + 1)
    r = 0.02; mu = 0.08; sigma = 0.18; gamma = 3.0
    rho = 0.04
    pi_star = (mu - r) / (gamma * sigma**2)
    m = (rho - (1 - gamma) * (r + (mu - r) ** 2 / (2 * gamma * sigma**2))) / gamma
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for k in range(8):
        Z = np.random.randn(N)
        # Wealth dynamics under optimal control:
        # dX/X = (r + pi*(mu-r) - m) dt + pi*sigma dW
        drift = r + pi_star * (mu - r) - m
        vol = pi_star * sigma
        logX = np.cumsum((drift - 0.5 * vol**2) * dt + vol * np.sqrt(dt) * Z)
        X = np.exp(np.concatenate([[0], logX]))
        axes[0].plot(t, X, lw=1, alpha=0.7)
        axes[1].plot(t, m * X, lw=1, alpha=0.7)
    axes[0].set_title(f"富 $X_t$（$\\pi^*={pi_star:.2f}$、$m={m:.3f}$）")
    axes[0].set_xlabel("t (years)"); axes[0].set_ylabel(r"$X_t / X_0$")
    axes[0].set_yscale("log")
    axes[1].set_title("最適消費率 $c_t^* = m\\cdot X_t$")
    axes[1].set_xlabel("t (years)"); axes[1].set_ylabel(r"$c_t / X_0$")
    axes[1].set_yscale("log")
    plt.tight_layout()
    plt.savefig(HERE / "merton_paths.png", bbox_inches="tight")
    plt.close()
    print("wrote merton_paths.png")


def fig_merton_policy():
    mu_grid = np.linspace(0.03, 0.15, 60)
    sigma_grid = np.linspace(0.05, 0.4, 60)
    M, S = np.meshgrid(mu_grid, sigma_grid)
    r = 0.02; gamma = 3.0
    pi = (M - r) / (gamma * S**2)
    fig, ax = plt.subplots()
    cf = ax.contourf(M, S, np.clip(pi, 0, 3), levels=20, cmap="viridis")
    cs = ax.contour(M, S, pi, levels=[0.25, 0.5, 1.0, 1.5, 2.0], colors="white", linewidths=0.8)
    ax.clabel(cs, inline=True, fontsize=8)
    plt.colorbar(cf, ax=ax, label=r"$\pi^*$（危険資産比率）")
    ax.set_xlabel(r"$\mu$"); ax.set_ylabel(r"$\sigma$")
    ax.set_title(fr"Merton の最適投資比率 ($r={r}$, $\gamma={gamma}$)")
    save("merton_policy")


# -----------------------------------------------------------------------------
# 第11章: Sharpe 比のサンプリング分布と必要期間
# -----------------------------------------------------------------------------
def fig_sharpe_distribution():
    Ts = [12, 36, 60, 120, 240]
    SR = 0.5
    x = np.linspace(-0.6, 1.6, 400)
    fig, ax = plt.subplots()
    for T, color in zip(Ts, plt.cm.plasma(np.linspace(0.15, 0.85, len(Ts)))):
        std = np.sqrt((1 + SR**2 / 2) / T)
        ax.plot(x, norm.pdf(x, loc=SR, scale=std), color=color, lw=2, label=f"T={T} months")
    ax.axvline(0, color="gray", ls=":")
    ax.axvline(SR, color="black", ls="--", label=fr"真値 SR={SR}")
    ax.set_xlabel(r"$\hat{SR}$")
    ax.set_ylabel("密度")
    ax.set_title("標本 Sharpe 比の漸近分布")
    ax.legend(fontsize=9)
    save("sharpe_distribution")


def fig_drawdown():
    np.random.seed(2)
    T = 1000
    mu, sigma = 0.0006, 0.012
    r = np.random.randn(T) * sigma + mu
    V = np.exp(np.cumsum(r))
    peak = np.maximum.accumulate(V)
    dd = (V - peak) / peak
    fig, axes = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
    axes[0].plot(V, color="C0", lw=1.4, label="価値")
    axes[0].plot(peak, color="C3", lw=1, ls="--", label="ピーク")
    axes[0].set_ylabel("累積価値"); axes[0].legend(); axes[0].set_title("累積価値と最大ドローダウン")
    axes[1].fill_between(np.arange(T), dd, 0, color="C3", alpha=0.4)
    axes[1].plot(dd, color="C3", lw=1)
    axes[1].set_ylabel("Drawdown")
    axes[1].set_xlabel("t (日)")
    plt.tight_layout()
    plt.savefig(HERE / "drawdown.png", bbox_inches="tight")
    plt.close()
    print("wrote drawdown.png")


# -----------------------------------------------------------------------------
# 第5章: ファクターモデルのリスク分解
# -----------------------------------------------------------------------------
def fig_factor_decomposition():
    np.random.seed(1)
    n = 25
    # generate factor model
    B = np.random.randn(n, 3) * np.array([0.8, 0.4, 0.3])
    Omega = np.diag([0.04, 0.02, 0.015])
    D = np.diag(np.random.uniform(0.005, 0.04, n))
    Sigma = B @ Omega @ B.T + D
    w = np.ones(n) / n
    total = w @ Sigma @ w
    factor_var = w @ B @ Omega @ B.T @ w
    specific_var = w @ D @ w
    # By factor:
    per_factor = np.array([w @ (B[:, [k]] @ Omega[k:k+1, k:k+1] @ B[:, [k]].T) @ w for k in range(3)])
    labels = ["因子1\n(market)", "因子2\n(size)", "因子3\n(value)", "特異リスク"]
    vals = np.concatenate([per_factor, [specific_var]])
    pct = vals / total * 100
    fig, ax = plt.subplots()
    colors = ["C0", "C2", "C1", "C7"]
    bars = ax.bar(labels, pct, color=colors)
    for b, v in zip(bars, pct):
        ax.annotate(f"{v:.1f}%", (b.get_x() + b.get_width() / 2, v),
                    xytext=(0, 3), textcoords="offset points", ha="center")
    ax.set_ylabel("分散寄与率 (%)")
    ax.set_title("3 因子モデルによる分散分解 (等加重ポートフォリオ)")
    save("factor_decomposition")


# -----------------------------------------------------------------------------
# 第0章: 共分散行列の固有値スペクトル
# -----------------------------------------------------------------------------
def fig_eigenspectrum():
    np.random.seed(0)
    n = 50
    T = 200
    R = np.random.randn(T, n) * 0.02
    Sigma = np.cov(R, rowvar=False)
    eig = np.sort(np.linalg.eigvalsh(Sigma))[::-1]
    fig, ax = plt.subplots()
    ax.bar(np.arange(1, n + 1), eig, color="C0")
    ax.set_xlabel("インデックス k")
    ax.set_ylabel("固有値 $\\lambda_k$")
    ax.set_title("標本共分散行列の固有値（i.i.d.ノイズの場合 — Marchenko–Pastur 状）")
    save("eigenspectrum")


def main():
    fig_diversification()
    fig_efficient_frontier()
    fig_two_fund()
    fig_sml()
    fig_factor_decomposition()
    fig_utility()
    fig_jensen()
    fig_stochastic_dominance()
    fig_var_cvar()
    fig_var_subadditivity_violation()
    fig_black_litterman()
    fig_merton_paths()
    fig_merton_policy()
    fig_sharpe_distribution()
    fig_drawdown()
    fig_eigenspectrum()


if __name__ == "__main__":
    main()
