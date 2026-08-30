# -*- coding: utf-8 -*-
"""组装「量子模拟基准数据集 v1」：arithmetic-origin 1D 紧束缚哈密顿。
四个哈密顿全部由 L 函数算术递推系数构造（Jacobi 三对角 alpha_n, b_n），
构造不输入任何零点/素数；零点仅事后对照。数据来源为 qcao 论文复现 bundle。
输出：每个哈密顿 hamiltonian.csv（三对角稠密矩阵）+ metadata.json，
外加 summary_table.csv。
用法：python3 build_dataset.py [SRC_BUNDLE_DIR] [OUT_DIR]
"""
import json, os, csv, sys
import numpy as np

SRC = sys.argv[1] if len(sys.argv) > 1 else \
    "/Coze/Drive/黎曼猜想论文审核/所有对话/主对话/zenodo_qcao_v1"
OUT = sys.argv[2] if len(sys.argv) > 2 else \
    "/Coze/Drive/黎曼猜想论文审核/所有对话/主对话/benchmark_v1/quantum_arithmetic_v1"
os.makedirs(OUT, exist_ok=True)

def write_hamiltonian(slug, name, alpha, b, meta_extra):
    N = len(alpha)
    H = np.diag(alpha) + np.diag(b, 1) + np.diag(b, -1)
    d = os.path.join(OUT, slug); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "hamiltonian.csv"), "w", newline="") as f:
        w = csv.writer(f)
        for row in H: w.writerow(["%.17g" % x for x in row])
    meta = {
        "slug": slug, "name": name,
        "matrix_size": N,
        "form": "real symmetric tridiagonal (Jacobi): H_nn=alpha_n, H_n,n+1=b_n",
        "units": "dimensionless; eigenvalues lambda ~ 1/gamma^2 where gamma = zero heights on the critical line",
        "alpha_diag": [float(x) for x in alpha],
        "b_offdiag": [float(x) for x in b],
        "construction": "entire L-function Xi -> contour sampling on circle R=2 -> Taylor/log moment recursion -> Gram-Schmidt orthonormal polynomials -> Jacobi coefficients; no zeros/primes enter the construction",
    }
    meta.update(meta_extra)
    json.dump(meta, open(os.path.join(d, "metadata.json"), "w"), indent=1, ensure_ascii=False)
    print(f"{slug}: N={N} written")
    return N

def read_matrix_csv(path):
    rows = []
    with open(path) as f:
        for line in f:
            vals = [float(x) for x in line.replace("\r", "").strip().split(",") if x.strip() != ""]
            if vals: rows.append(vals)
    H = np.array(rows)
    return np.diag(H).tolist(), np.diag(H, 1).tolist()

summary = []

# 1) zeta J50
a, b = read_matrix_csv(os.path.join(SRC, "coefficients/zeta_J50_matrix.csv"))
write_hamiltonian("zeta_J50", "Riemann zeta L-function, Jacobi truncation J=50", a, b, {
    "l_function": "Riemann zeta(s); completed xi(s)=1/2 s(s-1) pi^(-s/2) Gamma(s/2) zeta(s), FE xi(s)=xi(1-s)",
    "family": "Riemann zeta", "dirichlet_coefficients": "a_n=1 (zeta)",
    "pipeline": {"R": 2.0, "M": 1536, "dps": ">=230 (same pipeline as J100, lower precision)", "JMAX": 110},
    "locked_levels": 25, "locked_tol": 1e-4,
    "spectral_observation": "Wigner-Dyson level statistics; localized wavepackets; noise threshold eps_c~0.2-0.3",
    "reference_zeros": "zeta non-trivial zeros 1/2 + i*gamma_k",
})
summary.append(("zeta_J50", "Riemann zeta", len(a), 25))

# 2) beta J50
a, b = read_matrix_csv(os.path.join(SRC, "coefficients/beta_J50_matrix.csv"))
write_hamiltonian("beta_J50", "Dirichlet L-function mod 4 (beta), Jacobi truncation J=50", a, b, {
    "l_function": "beta(s)=4^-s(zeta(s,1/4)-zeta(s,3/4)); Dirichlet character chi_4 (Catalan G at s=2)",
    "family": "Dirichlet L mod 4", "dirichlet_coefficients": "a_n = chi_4(n) (1 for n=1 mod4, -1 for n=3 mod4, 0 even)",
    "pipeline": {"R": 2.0, "M": 1536, "dps": ">=230 (same pipeline as J100, lower precision)", "JMAX": 110},
    "locked_levels": 27, "locked_tol": 1e-4,
    "spectral_observation": "Wigner-Dyson level statistics; localized wavepackets; eps_c~0.2-0.3",
    "reference_zeros": "beta(s) non-trivial zeros on critical line",
})
summary.append(("beta_J50", "Dirichlet mod 4", len(a), 27))

# 3) Delta J22
dj = json.load(open(os.path.join(SRC, "coefficients/delta_jacobi_raw.json")))
alphas = [float(x) for x in dj["J22"]["alphas"]]
bvals = [float(np.sqrt(float(x))) for x in dj["J22"]["betas_sq"]]
write_hamiltonian("delta_J22", "Ramanujan tau modular L-function (weight-12 cusp form Delta), Jacobi J=22",
                  alphas, bvals, {
    "l_function": "L(Delta,s) = sum tau(n) n^-s; Ramanujan tau (weight-12 cusp form), FE s <-> 12-s",
    "family": "modular form (cusp)", "dirichlet_coefficients": "a_n = tau(n) (Ramanujan tau function)",
    "pipeline": {"R": 4.0, "M": 1024, "dps": 230,
                 "note": "Hankel conditioning floor at n~24 at 230dps; J22 is the reachable truncation"},
    "locked_levels": 15, "locked_tol": 1e-4,
    "spectral_observation": "Wigner-Dyson consistent (small sample: 15 locked levels, flagged throughout paper); eps_c~0.2-0.3",
    "reference_zeros": "L(Delta,s) non-trivial zeros",
})
summary.append(("delta_J22", "Ramanujan Delta", len(alphas), 15))

# 4) product J100
pj = json.load(open(os.path.join(SRC, "experiments/product_J100_results.json")))
write_hamiltonian("product_J100",
                  "Product xi_K = xi(zeta) * Lambda_beta (Gaussian-integer Dedekind zeta numerator), J=100",
                  pj["alpha"], pj["b"], {
    "l_function": "xi_K(s) = xi(s) * Lambda_beta(s); product of two primitive L-families with independent zero sets",
    "family": "product: zeta x Dirichlet-mod-4 (Gaussian field K=Q(i) numerator)",
    "dirichlet_coefficients": "convolution of zeta coefficients with chi_4",
    "pipeline": {"R": 2.0, "M": 1536, "dps": 650, "JMAX": 210},
    "locked_levels": 53, "locked_tol": 1e-4,
    "locked_composition": "35 beta-side + 18 zeta-side",
    "spectral_observation": "Wigner-Dyson across the union of two independently-arising families (Conjecture: superposition of independent families stays Wigner-Dyson)",
    "reference_zeros": "union of zeta and beta zeros",
})
summary.append(("product_J100", "zeta x beta product", len(pj["alpha"]), 53))

with open(os.path.join(OUT, "summary_table.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["slug", "family", "matrix_size_N", "locked_levels_relerr<1e-4"])
    for s in summary: w.writerow(s)
print("summary_table.csv written ->", OUT)
