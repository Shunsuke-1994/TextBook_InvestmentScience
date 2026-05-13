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


# -----------------------------------------------------------------------------
# v3 新章用の図
# -----------------------------------------------------------------------------

def fig_marchenko_pastur():
    np.random.seed(0)
    p, T = 100, 200
    c = p / T
    R = np.random.randn(T, p)
    S = R.T @ R / T
    eig = np.linalg.eigvalsh(S)
    lo, hi = (1 - np.sqrt(c)) ** 2, (1 + np.sqrt(c)) ** 2
    lam = np.linspace(max(1e-3, lo), hi, 400)
    density = np.sqrt(np.maximum((hi - lam) * (lam - lo), 0)) / (2 * np.pi * c * lam)
    fig, ax = plt.subplots()
    ax.hist(eig, bins=30, density=True, alpha=0.5, label="標本固有値", color="C0")
    ax.plot(lam, density, color="C3", lw=2.4, label="MP 理論密度")
    ax.axvline(lo, color="gray", ls="--", lw=1)
    ax.axvline(hi, color="gray", ls="--", lw=1)
    ax.set_xlabel(r"固有値 $\lambda$")
    ax.set_ylabel("密度")
    ax.set_title(fr"Marchenko–Pastur 分布 ($\Sigma=I$, $p={p}$, $T={T}$, $c={c}$)")
    ax.legend()
    save("marchenko_pastur")


def fig_shrinkage_compare():
    np.random.seed(1)
    n = 40
    T_in = 80
    n_repeats = 40
    B = np.random.randn(n, 3) * 0.3
    Sigma = B @ B.T + 0.1 * np.eye(n)
    methods = ["Sample", "Linear (LW)", "Nonlinear"]
    var_out_all = {m: [] for m in methods}
    for _ in range(n_repeats):
        R_in = np.random.multivariate_normal(np.zeros(n), Sigma, T_in)
        S = np.cov(R_in, rowvar=False)
        mu = np.trace(S) / n
        alpha_lw = 0.3
        S_lw = (1 - alpha_lw) * S + alpha_lw * mu * np.eye(n)
        eigvals, eigvec = np.linalg.eigh(S)
        eig_capped = np.maximum(eigvals, np.median(eigvals))
        S_nl = eigvec @ np.diag(eig_capped) @ eigvec.T
        for name, Shat in zip(methods, [S, S_lw, S_nl]):
            try:
                w = np.linalg.solve(Shat, np.ones(n))
                w = w / w.sum()
                var_out_all[name].append(float(w @ Sigma @ w))
            except np.linalg.LinAlgError:
                continue
    fig, ax = plt.subplots()
    ax.boxplot([var_out_all[m] for m in methods], tick_labels=methods, showmeans=True)
    ax.set_ylabel("真の分散 $w^\\top \\Sigma w$")
    ax.set_title("GMV ポートフォリオの真分散：縮小推定の効果")
    save("shrinkage_compare")


def fig_robust_frontier():
    np.random.seed(2)
    mu = np.array([0.06, 0.10, 0.14])
    Sigma = np.array([
        [0.0100, 0.0018, 0.0011],
        [0.0018, 0.0400, 0.0026],
        [0.0011, 0.0026, 0.0900],
    ])
    Theta_diag = np.array([0.01, 0.02, 0.03]) ** 2
    ones = np.ones_like(mu)
    inv = np.linalg.inv(Sigma)
    A = ones @ inv @ ones
    fig, ax = plt.subplots()
    for kappa, color in zip([0.0, 0.5, 1.0, 1.5], plt.cm.viridis(np.linspace(0.1, 0.85, 4))):
        ret_grid = np.linspace(0.02, 0.18, 80)
        sds = []
        mu_eff = mu - kappa * np.sqrt(Theta_diag)
        B_eff = ones @ inv @ mu_eff
        C_eff = mu_eff @ inv @ mu_eff
        D_eff = A * C_eff - B_eff ** 2
        for mu_target in ret_grid:
            sd = np.sqrt(max((A * mu_target ** 2 - 2 * B_eff * mu_target + C_eff) / D_eff, 0))
            sds.append(sd)
        ax.plot(sds, ret_grid, color=color, lw=2, label=fr"$\kappa={kappa}$")
    ax.set_xlabel(r"$\sigma_P$")
    ax.set_ylabel(r"$\mu_P$")
    ax.set_title("ロバスト効率フロンティア（半径 $\\kappa$）")
    ax.legend()
    save("robust_frontier")


def fig_hrp_dendrogram():
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import squareform
    np.random.seed(3)
    n = 15
    block = np.array([[1.0, 0.7, 0.6, 0.5, 0.4],
                      [0.7, 1.0, 0.6, 0.5, 0.4],
                      [0.6, 0.6, 1.0, 0.5, 0.4],
                      [0.5, 0.5, 0.5, 1.0, 0.6],
                      [0.4, 0.4, 0.4, 0.6, 1.0]])
    C = 0.1 * np.ones((n, n))
    C[:5, :5] = block
    C[5:10, 5:10] = block
    C[10:, 10:] = block[:5, :5]
    np.fill_diagonal(C, 1.0)
    C = (C + C.T) / 2
    dist = np.sqrt(np.clip((1 - C) / 2, 0, 1))
    np.fill_diagonal(dist, 0)
    Z = linkage(squareform(dist, checks=False), method="single")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    dendrogram(Z, ax=axes[0], no_labels=True, color_threshold=0.4)
    axes[0].set_title("階層クラスタリング (distance = $\\sqrt{(1-C)/2}$)")
    axes[0].set_ylabel("distance")
    order = leaves_list(Z)
    C_qd = C[np.ix_(order, order)]
    im = axes[1].imshow(C_qd, cmap="RdBu_r", vmin=-1, vmax=1)
    plt.colorbar(im, ax=axes[1], fraction=0.046)
    axes[1].set_title("Quasi-diagonalized 相関行列")
    plt.tight_layout()
    plt.savefig(HERE / "hrp_dendrogram.png", bbox_inches="tight")
    plt.close()
    print("wrote hrp_dendrogram.png")


def fig_factor_models_compare():
    models = ["CAPM", "FF3", "Carhart\nFF3+MOM", "FF5", "FF6", "q-factor", "q5"]
    scores = [0.42, 0.62, 0.71, 0.74, 0.80, 0.82, 0.86]
    spread_unexp = [12, 6, 4, 5, 3, 3, 2]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].bar(models, scores, color="C0")
    axes[0].set_ylabel("Cross-sectional $R^2$ (illustrative)")
    axes[0].set_title("ファクターモデルの平均的説明力")
    axes[0].set_ylim(0, 1.0)
    axes[0].tick_params(axis="x", rotation=30)
    axes[1].bar(models, spread_unexp, color="C3")
    axes[1].set_ylabel("Unexplained anomalies (個数)")
    axes[1].set_title("残存アノマリ数")
    axes[1].tick_params(axis="x", rotation=30)
    plt.tight_layout()
    plt.savefig(HERE / "factor_models_compare.png", bbox_inches="tight")
    plt.close()
    print("wrote factor_models_compare.png")


def fig_deflated_sharpe():
    Ns = np.array([5, 10, 25, 50, 100, 250, 500, 1000, 2500])
    sigma_SR = 0.5
    gamma_E = 0.5772
    Esup = sigma_SR * (
        (1 - gamma_E) * norm.ppf(1 - 1 / Ns)
        + gamma_E * norm.ppf(1 - 1 / (Ns * np.e))
    )
    leading = sigma_SR * np.sqrt(2 * np.log(Ns))
    fig, ax = plt.subplots()
    ax.plot(Ns, Esup, "o-", label="Bailey–López de Prado 補正式", color="C0", lw=2)
    ax.plot(Ns, leading, "s--", label=r"$\sigma_{SR}\sqrt{2\log N}$ (粗近似)", color="C3", lw=2)
    ax.set_xscale("log")
    ax.set_xlabel("試行戦略数 N")
    ax.set_ylabel(r"$\mathbb{E}[\max_i \widehat{SR}_i]$ (真値 0 のとき)")
    ax.set_title("戦略数増加と「無からの最大 Sharpe」")
    ax.legend()
    save("deflated_sharpe")


# -----------------------------------------------------------------------------
# v5 デリバティブ章用の図
# -----------------------------------------------------------------------------

def _bsm_d1d2(S, K, r, T, sigma):
    d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return d1, d2


def _bsm_call(S, K, r, T, sigma):
    d1, d2 = _bsm_d1d2(S, K, r, T, sigma)
    return S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)


def _bsm_put(S, K, r, T, sigma):
    d1, d2 = _bsm_d1d2(S, K, r, T, sigma)
    return K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)


def fig_option_payoffs():
    S = np.linspace(40, 160, 400)
    K = 100
    F0 = 100
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    # Long forward
    axes[0].plot(S, S - F0, color="C0", lw=2.2, label="ロング先渡")
    axes[0].plot(S, -(S - F0), color="C3", lw=2.2, label="ショート先渡", ls="--")
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].axvline(F0, color="gray", ls=":")
    axes[0].set_title(f"先渡（forward, $F_0$={F0}）")
    axes[0].set_xlabel(r"$S_T$"); axes[0].set_ylabel("ペイオフ")
    axes[0].legend()
    # Call
    axes[1].plot(S, np.maximum(S - K, 0), color="C0", lw=2.2, label="ロング・コール")
    axes[1].plot(S, -np.maximum(S - K, 0), color="C3", lw=2.2, label="ショート・コール", ls="--")
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].axvline(K, color="gray", ls=":")
    axes[1].set_title(f"コール（行使価格 $K$={K}）")
    axes[1].set_xlabel(r"$S_T$"); axes[1].set_ylabel("ペイオフ")
    axes[1].legend()
    # Put
    axes[2].plot(S, np.maximum(K - S, 0), color="C0", lw=2.2, label="ロング・プット")
    axes[2].plot(S, -np.maximum(K - S, 0), color="C3", lw=2.2, label="ショート・プット", ls="--")
    axes[2].axhline(0, color="black", lw=0.6)
    axes[2].axvline(K, color="gray", ls=":")
    axes[2].set_title(f"プット（行使価格 $K$={K}）")
    axes[2].set_xlabel(r"$S_T$"); axes[2].set_ylabel("ペイオフ")
    axes[2].legend()
    plt.tight_layout()
    plt.savefig(HERE / "option_payoffs.png", bbox_inches="tight")
    plt.close()
    print("wrote option_payoffs.png")


def fig_put_call_parity():
    S = np.linspace(40, 160, 400)
    K = 100; r = 0.05; T = 1
    disc_bond = K * np.exp(-r * T)
    call_payoff = np.maximum(S - K, 0)
    put_payoff = np.maximum(K - S, 0)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    # Portfolio A: Call + discount bond
    axes[0].plot(S, call_payoff + K, color="C0", lw=2.2, label="コール + 額面 $K$ 債券")
    axes[0].set_title("ポートフォリオ A：コール + 割引債")
    axes[0].axvline(K, color="gray", ls=":")
    axes[0].set_xlabel(r"$S_T$"); axes[0].set_ylabel("満期ペイオフ")
    axes[0].plot(S, np.maximum(S, K), color="C3", lw=1.2, ls="--", label=r"$\max(S_T,K)$")
    axes[0].legend()
    # Portfolio B: Put + stock
    axes[1].plot(S, put_payoff + S, color="C0", lw=2.2, label="プット + 株")
    axes[1].set_title("ポートフォリオ B：プット + 株")
    axes[1].axvline(K, color="gray", ls=":")
    axes[1].set_xlabel(r"$S_T$"); axes[1].set_ylabel("満期ペイオフ")
    axes[1].plot(S, np.maximum(S, K), color="C3", lw=1.2, ls="--", label=r"$\max(S_T,K)$")
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(HERE / "put_call_parity.png", bbox_inches="tight")
    plt.close()
    print("wrote put_call_parity.png")


def fig_option_arbitrage_bounds():
    S = np.linspace(40, 160, 400)
    K = 100; r = 0.05; T = 1; sigma = 0.2
    call = _bsm_call(S, K, r, T, sigma)
    lower_call = np.maximum(S - K * np.exp(-r * T), 0)
    upper_call = S
    fig, ax = plt.subplots()
    ax.plot(S, call, color="C0", lw=2.2, label="BSM コール価格")
    ax.plot(S, lower_call, color="C2", lw=1.6, ls="--", label=r"下限 $\max(S - Ke^{-rT}, 0)$")
    ax.plot(S, upper_call, color="C3", lw=1.6, ls="--", label=r"上限 $S$")
    ax.fill_between(S, lower_call, upper_call, alpha=0.1, color="C0")
    ax.set_xlabel(r"$S$"); ax.set_ylabel("コール価格")
    ax.set_title("コール価格の裁定範囲（モデルなしで成立）")
    ax.legend()
    save("option_arbitrage_bounds")


def fig_binomial_tree():
    S0 = 100; u = 1.2; d = 0.85; r = 0.05; dt = 1.0; K = 100
    n = 3
    fig, ax = plt.subplots(figsize=(8, 5))
    nodes = {}
    for i in range(n + 1):
        for j in range(i + 1):
            S = S0 * (u ** (i - j)) * (d ** j)
            nodes[(i, j)] = S
            ax.plot(i, S, "o", color="C0", markersize=8, zorder=3)
            ax.annotate(f"{S:.1f}", (i, S), xytext=(7, 7),
                        textcoords="offset points", fontsize=9)
    for i in range(n):
        for j in range(i + 1):
            ax.plot([i, i + 1], [nodes[(i, j)], nodes[(i + 1, j)]], color="gray", lw=0.8)
            ax.plot([i, i + 1], [nodes[(i, j)], nodes[(i + 1, j + 1)]], color="gray", lw=0.8)
    # Option values at maturity
    ax.set_xticks(range(n + 1))
    ax.set_xticklabels([f"$t={i}$" for i in range(n + 1)])
    ax.set_ylabel("株価 S")
    ax.set_title(f"3期間二項ツリー (u={u}, d={d}, r={r}, $S_0=${S0})")
    ax.grid(True, alpha=0.3)
    save("binomial_tree")


def fig_binomial_convergence():
    S0 = 100; K = 100; r = 0.05; T = 1; sigma = 0.2
    bsm_price = _bsm_call(S0, K, r, T, sigma)
    Ns = np.array([2, 4, 8, 16, 32, 64, 128, 256, 512])
    prices = []
    for n in Ns:
        dt = T / n
        u = np.exp(sigma * np.sqrt(dt))
        d = 1 / u
        q = (np.exp(r * dt) - d) / (u - d)
        # Backward induction
        ST = S0 * (u ** np.arange(n, -1, -1)) * (d ** np.arange(0, n + 1, 1))
        V = np.maximum(ST - K, 0)
        for _ in range(n):
            V = np.exp(-r * dt) * (q * V[:-1] + (1 - q) * V[1:])
        prices.append(float(V[0]))
    fig, ax = plt.subplots()
    ax.plot(Ns, prices, "o-", color="C0", lw=2, label="二項モデル価格")
    ax.axhline(bsm_price, color="C3", lw=2, ls="--", label=f"BSM 価格 = {bsm_price:.3f}")
    ax.set_xscale("log")
    ax.set_xlabel("ステップ数 $n$")
    ax.set_ylabel("コール価格")
    ax.set_title("二項モデル → BSM 連続極限への収束")
    ax.legend()
    save("binomial_convergence")


def fig_bsm_price_curve():
    S = np.linspace(40, 160, 400)
    K = 100; r = 0.05; sigma = 0.2
    fig, ax = plt.subplots()
    for T, color in zip([0.05, 0.25, 1.0, 2.0], plt.cm.viridis(np.linspace(0.15, 0.85, 4))):
        c = _bsm_call(S, K, r, T, sigma)
        ax.plot(S, c, color=color, lw=2, label=f"T={T}")
    ax.plot(S, np.maximum(S - K, 0), color="black", lw=1.4, ls=":", label="満期ペイオフ")
    ax.set_xlabel(r"$S_0$"); ax.set_ylabel("コール価格")
    ax.set_title(fr"BSM コール価格 ($K={K}$, $r={r}$, $\sigma={sigma}$)")
    ax.legend()
    save("bsm_price_curve")


def fig_bsm_greeks():
    S = np.linspace(40, 160, 400)
    K = 100; r = 0.05; sigma = 0.2
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for T, color in zip([0.25, 1.0], ["C0", "C3"]):
        d1, d2 = _bsm_d1d2(S, K, r, T, sigma)
        delta_c = norm.cdf(d1)
        gamma = norm.pdf(d1) / (S * sigma * np.sqrt(T))
        vega = S * norm.pdf(d1) * np.sqrt(T) / 100  # per 1 vol point
        theta_c = (-S * norm.pdf(d1) * sigma / (2 * np.sqrt(T))
                   - r * K * np.exp(-r * T) * norm.cdf(d2)) / 365  # per day
        axes[0, 0].plot(S, delta_c, color=color, lw=2, label=f"T={T}")
        axes[0, 1].plot(S, gamma, color=color, lw=2, label=f"T={T}")
        axes[1, 0].plot(S, vega, color=color, lw=2, label=f"T={T}")
        axes[1, 1].plot(S, theta_c, color=color, lw=2, label=f"T={T}")
    titles = ["Delta $\\Delta$", "Gamma $\\Gamma$", "Vega (per 1 vol pt)", "Theta (per day)"]
    for ax, t in zip(axes.flatten(), titles):
        ax.axvline(K, color="gray", ls=":", lw=1)
        ax.set_xlabel(r"$S$")
        ax.set_title(t)
        ax.legend()
    plt.tight_layout()
    plt.savefig(HERE / "bsm_greeks.png", bbox_inches="tight")
    plt.close()
    print("wrote bsm_greeks.png")


def fig_iv_smile():
    K_S = np.linspace(0.6, 1.4, 200)
    # 教育用の擬似スマイル
    iv = 0.18 + 0.25 * (K_S - 1) ** 2 + 0.08 * np.maximum(1 - K_S, 0)
    fig, ax = plt.subplots()
    ax.plot(K_S, iv, color="C0", lw=2.4)
    ax.axhline(0.20, color="C3", lw=1.2, ls="--", label=r"BSM 定数 $\sigma = 0.20$")
    ax.axvline(1.0, color="gray", ls=":")
    ax.set_xlabel(r"$K / S_0$（モネネス）")
    ax.set_ylabel("インプライド・ボラティリティ")
    ax.set_title("株式オプションの IV スマイル／スキュー (教育用擬似データ)")
    ax.legend()
    save("iv_smile")


def fig_delta_hedge_pnl():
    """離散デルタヘッジの P&L 分布。"""
    np.random.seed(42)
    S0 = 100; K = 100; r = 0.05; T = 0.25; sigma_true = 0.20
    n_paths = 1000

    def simulate(n_rebal):
        dt = T / n_rebal
        pnl = []
        for _ in range(n_paths):
            S = S0
            cash = -_bsm_call(S0, K, r, T, sigma_true)  # オプション売って受け取った
            d1, _ = _bsm_d1d2(S, K, r, T - 0.0, sigma_true)
            delta = norm.cdf(d1)
            cash -= delta * S0  # 株を delta 単位買う（売り手なので株を買ってヘッジ）
            holding = delta
            for k in range(1, n_rebal):
                Z = np.random.randn()
                S = S * np.exp((r - 0.5 * sigma_true ** 2) * dt + sigma_true * np.sqrt(dt) * Z)
                cash *= np.exp(r * dt)
                tau = T - k * dt
                d1, _ = _bsm_d1d2(S, K, r, tau, sigma_true)
                new_delta = norm.cdf(d1)
                cash -= (new_delta - holding) * S
                holding = new_delta
            # final step
            Z = np.random.randn()
            S = S * np.exp((r - 0.5 * sigma_true ** 2) * dt + sigma_true * np.sqrt(dt) * Z)
            cash *= np.exp(r * dt)
            payoff = max(S - K, 0)
            # 売ったオプションの満期支払い、株を売り清算
            pnl.append(cash + holding * S - payoff)
        return np.array(pnl)

    pnl_daily = simulate(63)  # 0.25年で 63 営業日（≈ 毎日リバランス）
    pnl_weekly = simulate(13)  # 週次
    pnl_monthly = simulate(3)  # 月次
    fig, ax = plt.subplots()
    for pnl, lab, color in [(pnl_monthly, "月次 (n=3)", "C3"),
                             (pnl_weekly, "週次 (n=13)", "C1"),
                             (pnl_daily, "日次 (n=63)", "C0")]:
        ax.hist(pnl, bins=40, alpha=0.5, label=f"{lab}, σ={pnl.std():.2f}", color=color, density=True)
    ax.axvline(0, color="black", lw=1)
    ax.set_xlabel("ヘッジ後 P&L")
    ax.set_ylabel("密度")
    ax.set_title("離散デルタヘッジの P&L 分布：リバランス頻度の効果")
    ax.legend()
    save("delta_hedge_pnl")


# -----------------------------------------------------------------------------
# v6 加筆用の図
# -----------------------------------------------------------------------------

def fig_three_asset_gmv():
    """3 資産 GMV の幾何：等分散曲面 + 効率フロンティア + GMV / 目標リターン点。"""
    mu = np.array([0.06, 0.10, 0.14])
    Sigma = np.array([
        [0.0100, 0.0018, 0.0011],
        [0.0018, 0.0400, 0.0026],
        [0.0011, 0.0026, 0.0900],
    ])
    ones = np.ones(3)
    inv = np.linalg.inv(Sigma)
    A = ones @ inv @ ones
    B = ones @ inv @ mu
    C = mu @ inv @ mu
    D = A * C - B * B
    mu_grid = np.linspace(0.04, 0.18, 200)
    var = (A * mu_grid ** 2 - 2 * B * mu_grid + C) / D
    sd = np.sqrt(var)
    gmv_mu = B / A
    gmv_sd = 1 / np.sqrt(A)
    # ランダムポートフォリオを散布
    np.random.seed(0)
    n_rand = 1500
    W = np.random.rand(n_rand, 3)
    W = W / W.sum(axis=1, keepdims=True)
    ports_mu = W @ mu
    ports_sd = np.sqrt(np.einsum("ij,jk,ik->i", W, Sigma, W))
    fig, ax = plt.subplots()
    sc = ax.scatter(ports_sd, ports_mu, c=ports_mu / ports_sd, s=6, alpha=0.35, cmap="viridis", label="ランダム比例ポートフォリオ")
    ax.plot(sd[mu_grid >= gmv_mu], mu_grid[mu_grid >= gmv_mu], color="C0", lw=2.5, label="効率的フロンティア")
    ax.plot(sd[mu_grid < gmv_mu], mu_grid[mu_grid < gmv_mu], color="C0", lw=1.4, ls="--", alpha=0.6, label="非効率部分")
    ax.scatter([gmv_sd], [gmv_mu], color="C3", s=90, zorder=6)
    ax.annotate("GMV", (gmv_sd, gmv_mu), xytext=(10, -10), textcoords="offset points", color="C3")
    target_mu = 0.12
    target_sd = np.sqrt((A * target_mu**2 - 2 * B * target_mu + C) / D)
    ax.scatter([target_sd], [target_mu], color="C1", s=90, zorder=6)
    ax.annotate(r"$\mu_P=12\%$ 効率点", (target_sd, target_mu), xytext=(10, -10),
                textcoords="offset points", color="C1")
    cb = plt.colorbar(sc, ax=ax, shrink=0.8)
    cb.set_label(r"Sharpe-like ratio $\mu/\sigma$")
    ax.set_xlabel(r"$\sigma_P$")
    ax.set_ylabel(r"$\mu_P$")
    ax.set_title("3 資産モデル：ランダムポートフォリオ vs 効率フロンティア")
    ax.legend(loc="lower right", fontsize=9)
    save("three_asset_gmv")


def fig_estimation_error_impact():
    """推定誤差で重みが極端化する数値例。"""
    np.random.seed(2)
    n = 5
    mu_true = np.array([0.06, 0.07, 0.08, 0.09, 0.10])
    sigma = np.array([0.18, 0.22, 0.20, 0.25, 0.30])
    R = 0.3  # 相関
    Sigma = np.outer(sigma, sigma) * R
    np.fill_diagonal(Sigma, sigma ** 2)
    # 標本推定
    n_repeats = 50
    weights = []
    for _ in range(n_repeats):
        R_sample = np.random.multivariate_normal(mu_true, Sigma, size=60)
        mu_hat = R_sample.mean(axis=0)
        Sigma_hat = np.cov(R_sample, rowvar=False)
        try:
            w = np.linalg.solve(Sigma_hat, mu_hat - 0.02)
            w = w / np.abs(w).sum()  # 1-norm 正規化（極端さを可視化）
            weights.append(w)
        except np.linalg.LinAlgError:
            pass
    weights = np.array(weights)
    fig, ax = plt.subplots()
    positions = np.arange(n)
    bp = ax.boxplot([weights[:, i] for i in range(n)], positions=positions,
                    widths=0.5, showmeans=True, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("C0")
        patch.set_alpha(0.5)
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"資産{i+1}" for i in range(n)])
    ax.set_ylabel("最適重み（標本ごとの分布）")
    ax.set_title("推定誤差による Markowitz 最適重みの不安定性 (T=60, n=5, 50 試行)")
    save("estimation_error_impact")


def fig_hyperbola_decomposition():
    """効率フロンティアの双曲線方程式：頂点・漸近線を可視化。"""
    mu = np.array([0.06, 0.10, 0.14])
    Sigma = np.array([
        [0.0100, 0.0018, 0.0011],
        [0.0018, 0.0400, 0.0026],
        [0.0011, 0.0026, 0.0900],
    ])
    ones = np.ones(3)
    inv = np.linalg.inv(Sigma)
    A = ones @ inv @ ones
    B = ones @ inv @ mu
    C = mu @ inv @ mu
    D = A * C - B * B
    mu_grid = np.linspace(-0.05, 0.25, 400)
    var = (A * mu_grid ** 2 - 2 * B * mu_grid + C) / D
    sd = np.sqrt(np.maximum(var, 0))
    gmv_mu = B / A
    gmv_sd = 1 / np.sqrt(A)
    asymp_slope = np.sqrt(D / A)
    fig, ax = plt.subplots()
    ax.plot(sd, mu_grid, color="C0", lw=2, label="効率フロンティア")
    # asymptotes
    sigma_line = np.linspace(0, 0.4, 100)
    ax.plot(sigma_line, gmv_mu + asymp_slope * sigma_line, color="C3", lw=1, ls="--", label="漸近線 (上下)")
    ax.plot(sigma_line, gmv_mu - asymp_slope * sigma_line, color="C3", lw=1, ls="--")
    ax.scatter([gmv_sd], [gmv_mu], color="C2", s=80, zorder=5)
    ax.annotate(f"頂点 (GMV)\n$\\sigma$={gmv_sd:.3f}\n$\\mu$={gmv_mu:.3f}",
                (gmv_sd, gmv_mu), xytext=(20, -30), textcoords="offset points", color="C2",
                arrowprops=dict(arrowstyle="->"))
    ax.set_xlim(-0.005, 0.4)
    ax.set_ylim(-0.05, 0.25)
    ax.set_xlabel(r"$\sigma_P$")
    ax.set_ylabel(r"$\mu_P$")
    ax.set_title("効率フロンティアの双曲線構造：頂点と漸近線")
    ax.axhline(0, color="black", lw=0.5)
    ax.legend()
    save("hyperbola_decomposition")


def fig_zero_beta_geometry():
    """ゼロベータポートフォリオの幾何：接線が縦軸を切る場所。"""
    mu = np.array([0.06, 0.10, 0.14])
    Sigma = np.array([
        [0.0100, 0.0018, 0.0011],
        [0.0018, 0.0400, 0.0026],
        [0.0011, 0.0026, 0.0900],
    ])
    ones = np.ones(3)
    inv = np.linalg.inv(Sigma)
    A = ones @ inv @ ones
    B = ones @ inv @ mu
    C = mu @ inv @ mu
    D = A * C - B * B
    mu_grid = np.linspace(0.04, 0.18, 200)
    sd = np.sqrt((A * mu_grid ** 2 - 2 * B * mu_grid + C) / D)
    # 任意の効率点 mu_p で接線を計算
    mu_p = 0.12
    sd_p = np.sqrt((A * mu_p ** 2 - 2 * B * mu_p + C) / D)
    mu_z = (C - B * mu_p) / (B - A * mu_p)  # zero-beta companion
    # 接線の傾き
    slope = (mu_p - mu_z) / sd_p
    sigma_line = np.linspace(0, sd_p * 1.5, 100)
    fig, ax = plt.subplots()
    ax.plot(sd, mu_grid, color="C0", lw=2, label="効率フロンティア")
    ax.plot(sigma_line, mu_z + slope * sigma_line, color="C3", lw=1.6, ls="--",
            label=f"接線 (μ_p={mu_p})")
    ax.scatter([sd_p], [mu_p], color="C2", s=80, zorder=5)
    ax.annotate(f"効率点 P\nμ={mu_p}", (sd_p, mu_p), xytext=(10, -5), textcoords="offset points", color="C2")
    ax.scatter([0], [mu_z], color="C3", s=80, zorder=5, marker="D")
    ax.annotate(f"ゼロβ 切片\nμ_z={mu_z:.3f}", (0, mu_z), xytext=(10, -5), textcoords="offset points", color="C3")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel(r"$\sigma_P$")
    ax.set_ylabel(r"$\mu_P$")
    ax.set_title("ゼロベータポートフォリオの幾何：効率点の接線切片")
    ax.legend(loc="lower right")
    save("zero_beta_geometry")


def fig_sml_detailed():
    """SML 図に複数銘柄を載せて Jensen α を可視化。"""
    rf = 0.02
    rm = 0.08
    np.random.seed(7)
    n = 12
    betas = np.random.uniform(0.3, 1.8, n)
    alphas = np.random.uniform(-0.02, 0.025, n)
    rets = rf + betas * (rm - rf) + alphas
    fig, ax = plt.subplots()
    beta_grid = np.linspace(0, 2.0, 100)
    sml = rf + beta_grid * (rm - rf)
    ax.plot(beta_grid, sml, color="C0", lw=2.5, label="SML")
    for b, r, a in zip(betas, rets, alphas):
        color = "C3" if a > 0.005 else ("C2" if a < -0.005 else "gray")
        ax.scatter([b], [r], color=color, s=70, edgecolor="black", zorder=5)
        if abs(a) > 0.01:
            ax.plot([b, b], [rf + b * (rm - rf), r], color=color, lw=1, alpha=0.6)
    ax.axhline(rf, color="gray", ls=":")
    ax.annotate(r"$r_f$", (0, rf), xytext=(-25, -4), textcoords="offset points")
    ax.annotate("M (β=1)", (1, rm), xytext=(10, -15), textcoords="offset points",
                arrowprops=dict(arrowstyle="->"))
    ax.scatter([1], [rm], color="black", s=80, zorder=5)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\mathbb{E}[R]$")
    ax.set_title("SML 上の銘柄：割安 (赤 $\\alpha>0$)・均衡 (灰)・割高 (緑 $\\alpha<0$)")
    ax.legend(loc="upper left")
    save("sml_detailed")


def fig_binomial_with_values():
    """二項ツリー全ノードに株価とオプション価値を表示。"""
    S0 = 100; K = 100; u = 1.2; d = 0.85; r = 0.05; dt = 1.0
    n = 3
    q = (np.exp(r * dt) - d) / (u - d)
    # 株価ノード
    S_tree = {}
    V_tree = {}
    for i in range(n + 1):
        for j in range(i + 1):
            S_tree[(i, j)] = S0 * (u ** (i - j)) * (d ** j)
    # 後ろ向き帰納
    for j in range(n + 1):
        V_tree[(n, j)] = max(S_tree[(n, j)] - K, 0)
    for i in range(n - 1, -1, -1):
        for j in range(i + 1):
            V_tree[(i, j)] = np.exp(-r * dt) * (q * V_tree[(i + 1, j)] + (1 - q) * V_tree[(i + 1, j + 1)])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for (i, j), S in S_tree.items():
        ax.plot(i, S, "o", color="C0", markersize=9, zorder=3)
        ax.annotate(f"S={S:.1f}\nV={V_tree[(i, j)]:.2f}",
                    (i, S), xytext=(8, 8), textcoords="offset points", fontsize=8.5)
    for i in range(n):
        for j in range(i + 1):
            ax.plot([i, i + 1], [S_tree[(i, j)], S_tree[(i + 1, j)]], color="gray", lw=0.8)
            ax.plot([i, i + 1], [S_tree[(i, j)], S_tree[(i + 1, j + 1)]], color="gray", lw=0.8)
    ax.set_xticks(range(n + 1))
    ax.set_xticklabels([f"t={i}" for i in range(n + 1)])
    ax.set_ylabel("株価")
    ax.set_title(f"3 期間二項ツリー: コール (K={K}) の後ろ向き帰納 (q={q:.3f})")
    ax.grid(True, alpha=0.3)
    save("binomial_with_values")


def fig_lognormal_evolution():
    """リスク中立 GBM の対数正規分布が時刻 T で広がる様子。"""
    S0 = 100; r = 0.05; sigma = 0.2
    Ts = [0.1, 0.25, 0.5, 1.0, 2.0]
    s = np.linspace(40, 250, 400)
    fig, ax = plt.subplots()
    for T, color in zip(Ts, plt.cm.viridis(np.linspace(0.15, 0.85, len(Ts)))):
        m = np.log(S0) + (r - 0.5 * sigma ** 2) * T
        v = sigma ** 2 * T
        density = np.exp(-0.5 * (np.log(s) - m) ** 2 / v) / (s * np.sqrt(2 * np.pi * v))
        ax.plot(s, density, color=color, lw=2, label=f"T={T}")
    ax.axvline(S0, color="gray", ls=":")
    ax.set_xlabel(r"$S_T$")
    ax.set_ylabel(r"密度 $f_{S_T}(s)$")
    ax.set_title(fr"リスク中立 GBM 下の $S_T$ の対数正規分布 ($S_0={S0}, r={r}, \sigma={sigma}$)")
    ax.legend()
    save("lognormal_evolution")


def fig_delta_hedge_pnl_decomp():
    """デルタヘッジ P&L を Gamma と Vega 効果に分解する図。"""
    np.random.seed(11)
    S0 = 100; K = 100; r = 0.05; T = 0.5; sigma = 0.2
    # 単一サンプルパスを使い、各日の Gamma P&L (1/2 Γ S² (ΔS/S)²) を計算
    N = 252 * T  # 営業日換算
    N = int(N)
    dt = T / N
    np.random.seed(7)
    Z = np.random.randn(N)
    S = [S0]
    for k in range(N):
        S.append(S[-1] * np.exp((r - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt) * Z[k]))
    S = np.array(S)
    # Gamma P&L
    times = np.linspace(0, T, N + 1)
    gamma_pnl = np.zeros(N)
    for k in range(N):
        tau = T - times[k]
        if tau <= 0:
            continue
        d1 = (np.log(S[k] / K) + (r + 0.5 * sigma ** 2) * tau) / (sigma * np.sqrt(tau))
        gamma = norm.pdf(d1) / (S[k] * sigma * np.sqrt(tau))
        dS = S[k + 1] - S[k]
        gamma_pnl[k] = 0.5 * gamma * dS ** 2
    cumul = np.cumsum(gamma_pnl)
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    axes[0].plot(times[:-1], S[:-1], color="C0", lw=1.4)
    axes[0].set_ylabel("株価")
    axes[0].set_title("デルタヘッジの P&L 蓄積：Gamma 由来の利益")
    axes[1].plot(times[:-1], cumul, color="C3", lw=1.4, label=r"累積 Gamma P&L $\sum \frac{1}{2} \Gamma (\Delta S)^2$")
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].set_xlabel("時刻 t (years)")
    axes[1].set_ylabel("累積 P&L")
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(HERE / "delta_hedge_decomp.png", bbox_inches="tight")
    plt.close()
    print("wrote delta_hedge_decomp.png")


def fig_iv_newton():
    """Newton 法で IV を逆算するイテレーション過程の可視化。"""
    S0 = 100; K = 100; r = 0.05; T = 0.5
    C_market = 6.5  # 市場価格
    sigma_grid = np.linspace(0.05, 0.6, 400)
    C_bs = _bsm_call(S0, K, r, T, sigma_grid)
    # Newton
    sigma = 0.3  # 初期値
    iters = [sigma]
    for _ in range(6):
        c = _bsm_call(S0, K, r, T, sigma)
        d1, _ = _bsm_d1d2(S0, K, r, T, sigma)
        vega = S0 * norm.pdf(d1) * np.sqrt(T)
        sigma = sigma - (c - C_market) / vega
        iters.append(sigma)
    fig, ax = plt.subplots()
    ax.plot(sigma_grid, C_bs, color="C0", lw=2.2, label="BSM コール価格 vs $\\sigma$")
    ax.axhline(C_market, color="C3", lw=1.6, ls="--", label=f"市場価格 {C_market}")
    for k, s in enumerate(iters):
        c = _bsm_call(S0, K, r, T, s)
        ax.scatter([s], [c], color="C2", s=70, zorder=5)
        ax.annotate(f"{k}", (s, c), xytext=(7, -3), textcoords="offset points", fontsize=9)
    ax.set_xlabel(r"$\sigma$")
    ax.set_ylabel("コール価格")
    ax.set_title("Newton 法でインプライド・ボラを逆算する反復過程")
    ax.legend()
    save("iv_newton")


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
    # v3 additions
    fig_marchenko_pastur()
    fig_shrinkage_compare()
    fig_robust_frontier()
    fig_hrp_dendrogram()
    fig_factor_models_compare()
    fig_deflated_sharpe()
    # v5 derivatives
    fig_option_payoffs()
    fig_put_call_parity()
    fig_option_arbitrage_bounds()
    fig_binomial_tree()
    fig_binomial_convergence()
    fig_bsm_price_curve()
    fig_bsm_greeks()
    fig_iv_smile()
    fig_delta_hedge_pnl()
    # v6 additions
    fig_three_asset_gmv()
    fig_estimation_error_impact()
    fig_hyperbola_decomposition()
    fig_zero_beta_geometry()
    fig_sml_detailed()
    fig_binomial_with_values()
    fig_lognormal_evolution()
    fig_delta_hedge_pnl_decomp()
    fig_iv_newton()


if __name__ == "__main__":
    main()
