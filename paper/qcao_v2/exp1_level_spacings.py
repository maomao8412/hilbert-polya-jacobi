# -*- coding: utf-8 -*-
"""Physics experiment 1: energy-level spacing statistics of the Jacobi Hamiltonians.
The Jacobi matrices ARE tight-binding Hamiltonians (alpha_n = on-site potential,
b_n = hopping amplitude). Their eigenvalues = L-function zeros (gamma).
Quantum chaos predicts Wigner-Dyson (GUE) level repulsion for chaotic quantum
systems; integrable systems give Poisson (exponential) spacings.
Only Euler coefficients (primes / Gaussian primes / tau(n)) ever entered the matrices;
zeros here are the measured eigenvalues, not inputs.
"""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

OUT = "/Coze/Drive/黎曼猜想论文审核/所有对话/主对话/physics_experiments"

# ---- 1. extract energy levels (locked eigenvalues) from the three Hamiltonians ----
prod = json.load(open("/Coze/Drive/黎曼猜想论文审核/所有对话/主对话/GRH交互网页/乘积J100_20260830/product_J100_results.json"))
g_zeta = sorted(r[1] for r in prod["rows"] if r[2]=="ζ" and abs(r[4]) < 1e-4)
g_beta = sorted(r[1] for r in prod["rows"] if r[2]=="β" and abs(r[4]) < 1e-4)

dd = json.load(open("/tmp/delta_matrix/delta_J100_results.json"))
ev_inv = sorted(dd["zeros"]["J22"]["eigenvalues_inv"], reverse=True)  # 1/gamma^2, largest first
g_delta = [1/math.sqrt(x) for x in ev_inv[:15]]  # first 15 levels (first 11 verified vs roots)

families = {"ζ (Riemann)": g_zeta, "β (Dirichlet mod 4)": g_beta, "Δ (Ramanujan τ)": g_delta}

# ---- 2. unfold with smooth cubic fit of counting function (small-sample robust) ----
def unfold(g):
    g = np.array(g)
    n = np.arange(1, len(g)+1)
    p = np.polyfit(g, n, 3)
    xs = np.polyval(p, g)
    gaps = np.diff(xs)
    return gaps / gaps.mean()

all_gaps = {}
for name, g in families.items():
    s = unfold(g)
    all_gaps[name] = s
    print(f"{name:22s}: {len(g):2d} levels, {len(s):2d} gaps, "
          f"min gap = {s.min():.3f}, mean = {s.mean():.3f}, var = {s.var():.3f}, "
          f"#gaps<0.3 = {(s<0.3).sum()}")

pooled = np.concatenate(list(all_gaps.values()))
print(f"\nPOOLED: {len(pooled)} gaps, min = {pooled.min():.3f}, var = {pooled.var():.3f}, "
      f"#gaps<0.3 = {(pooled<0.3).sum()}")

# superposition control: union spectrum of two independent Hamiltonians (ζ+β)
g_union = sorted(g_zeta + g_beta)
s_union = unfold(g_union)
print(f"UNION ζ+β (superposition): {len(g_union)} levels, min gap = {s_union.min():.3f}, "
      f"var = {s_union.var():.3f}, #gaps<0.3 = {(s_union<0.3).sum()}  (Poisson expectation high)")

# ---- 3. Monte Carlo: what would Poisson (integrable) give? ----
rng = np.random.default_rng(42)
NMC = 20000
counts = np.array([len(s) for s in all_gaps.values()])
mc_min, mc_small, mc_var = [], [], []
for _ in range(NMC):
    sp = np.concatenate([-np.log(rng.random(c)) for c in counts])  # exponential gaps
    sp = sp/sp.mean()
    mc_min.append(sp.min()); mc_small.append((sp<0.3).sum()); mc_var.append(sp.var())
p_min   = (np.array(mc_min)   <= pooled.min()).mean()
p_small = (np.array(mc_small) <= (pooled<0.3).sum()).mean()
p_var   = (np.array(mc_var)   <= pooled.var()).mean()
print(f"\nPoisson Monte Carlo ({NMC} runs): p-value for observed min gap = {p_min:.4f}, "
      f"for #small gaps = {p_small:.5f}, for variance = {p_var:.5f}")

# ---- 4. figure ----
s_grid = np.linspace(0, 3, 300)
P_GUE = (32/math.pi**2)*s_grid**2*np.exp(-4*s_grid**2/math.pi)
P_POI = np.exp(-s_grid)

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axes[0]
ax.hist(pooled, bins=np.linspace(0, 2.6, 14), density=True, color="#4C72B0", alpha=0.75,
        label=f"Jacobi eigenvalues pooled (n={len(pooled)} gaps)")
ax.plot(s_grid, P_GUE, "r-", lw=2, label="Wigner–Dyson GUE (chaotic)")
ax.plot(s_grid, P_POI, "k--", lw=2, label="Poisson (integrable)")
ax.set_xlabel("unfolded level spacing s"); ax.set_ylabel("P(s)")
ax.set_title("Level spacing of the L-function Jacobi Hamiltonians")
ax.legend(fontsize=9); ax.set_xlim(0, 2.6)
ax.text(0.97, 0.95, f"observed min gap = {pooled.min():.2f}\nPoisson p = {p_min:.4f}\n"
        f"gaps<0.3: {(pooled<0.3).sum()} observed\n(Poisson expects ~{0.259*len(pooled):.0f})",
        transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round", fc="#FFF3CD", ec="gray"))

ax = axes[1]
ax.hist(s_union, bins=np.linspace(0, 2.6, 14), density=True, color="#8C8C8C", alpha=0.75,
        label=f"superposition ζ+β (n={len(s_union)} gaps)")
ax.plot(s_grid, P_GUE, "r-", lw=2, label="Wigner–Dyson GUE")
ax.plot(s_grid, P_POI, "k--", lw=2, label="Poisson")
ax.set_xlabel("unfolded level spacing s"); ax.set_ylabel("P(s)")
ax.set_title("Control: superposition of two independent spectra → Poisson")
ax.legend(fontsize=9); ax.set_xlim(0, 2.6)
plt.tight_layout()
plt.savefig(f"{OUT}/exp1_level_spacings.png", dpi=150)
print("saved exp1_level_spacings.png")

json.dump({"pooled_gaps": pooled.tolist(), "union_gaps": s_union.tolist(),
           "per_family": {k: v.tolist() for k, v in all_gaps.items()},
           "obs_min_gap": float(pooled.min()), "obs_var": float(pooled.var()),
           "obs_small_gaps": int((pooled<0.3).sum()),
           "poisson_p_min": float(p_min), "poisson_p_small": float(p_small),
           "poisson_p_var": float(p_var)},
          open(f"{OUT}/exp1_results.json", "w"), indent=1)
