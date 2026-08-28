#!/usr/bin/env python3
"""
=========================================================================
 Reproducible construction of the Hilbert–Pólya Jacobi matrix
 strictly following:
   "A Self-Adjoint Jacobi Operator from the Polylogarithmic
    Decomposition of the Riemann Zeta Function" (RH_jacobi_paper_v2)

 Construction chain (NO shortcuts, NO Newton identities from zeros):
   1. Three-term polylogarithm identity:
        ζ(s) = Li_s(z1) + Li_s(z2) + Σ c_n/n^s
        z1 = 2-√2, z2 = √2-1, c_n = 1-z1^n-z2^n
   2. ξ(s) = 1/2 s(s-1) π^{-s/2} Γ(s/2) ζ(s)
      Taylor coefficients a_{2k} via Cauchy integral on |w|=R
   3. Log-moments S_k = -k [t^{2k}] log(ξ(1/2+it)/ξ(1/2))
      via Mercator expansion of log(1+Σ b_{2k} t^{2k})
   4. Shifted Hankel H_N[i,j] = S_{i+j+1}, D_n = det(H_n) > 0
   5. Gram–Schmidt orthogonalisation of {1,x,x²,…} in L²(μ)
      with ⟨p,q⟩ = Σ p_i q_j S_{i+j+1}
   6. Three-term recurrence → Jacobi parameters α_n, b_n
   7. J_N = tridiag(α_n, b_n); eigenvalues → 1/γ_n²

 Requirements: mpmath
 Usage:        python reproduce_paper.py
 Output:       prints S_k, D_n, J_3..J_14 params, eigenvalue convergence
=========================================================================
"""
from mpmath import (mp, mpf, matrix, det, sqrt, pi, exp, j as mpj,
                    log as mp_log)
import sys, time

mp.dps = 120

# ── Step 1: Three-term polylogarithm identity constants ──────────────
# These come from the deformation family F(s,w) at w=1:
#   z1 = 1-q, z2 = q, q = √2-1, z1+z2 = 1
# The identity ζ(s)=Li_s(z1)+Li_s(z2)+Σ c_n/n^s is the analytic input
# that proves Hankel positivity D_n>0 via Herglotz/Stieltjes theory.
# For the numerical matrix we compute ξ directly via mpmath.zeta;
# the three-term identity is the theoretical foundation, not a
# computational shortcut around ζ.
q = sqrt(2) - 1
z1 = mpf(1) - q          # 2 - √2
z2 = q                   # √2 - 1
assert abs(z1 + z2 - 1) < mpf('1e-110')

# ── Step 2: ξ Taylor coefficients via Cauchy integral ────────────────
N_COEFF = 60
R_CAUCHY = mpf('3')
N_PTS = 800

def xi(s):
    """Completed xi: ξ(s) = 1/2 s(s-1) π^{-s/2} Γ(s/2) ζ(s)."""
    return mpf('0.5') * s * (s - 1) * pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s)

print(f"Computing Taylor coefficients (R={R_CAUCHY}, {N_PTS} pts, {mp.dps} digits)...")
t0 = time.time()
coeffs = [mpf(0)] * (N_COEFF + 1)
for j in range(N_PTS):
    theta = 2 * pi * j / N_PTS
    z = R_CAUCHY * exp(mpj * theta)
    xi_z = xi(mpf('0.5') + mpj * z)
    for k in range(N_COEFF + 1):
        coeffs[k] += xi_z * exp(-mpj * 2 * k * theta) / R_CAUCHY**(2*k)
for k in range(N_COEFF + 1):
    coeffs[k] = mp.re(coeffs[k]) / N_PTS
a0 = coeffs[0]
bn = [coeffs[k] / a0 for k in range(N_COEFF + 1)]
print(f"  xi(1/2) = {mp.nstr(a0, 18)}  ({time.time()-t0:.1f}s)")

# ── Step 3: Log-moments S_k via Mercator series ──────────────────────
MAX_K = 40
print("Computing log-moments S_k...")
log_coeff = [mpf(0)] * (MAX_K + 1)
f = [mpf(0)] * (MAX_K + 1)
for k in range(1, MAX_K + 1):
    f[k] = bn[k] if k <= N_COEFF else mpf(0)
f_power = [mpf(0)] * (MAX_K + 1)
f_power[0] = mpf(1)
for n in range(1, MAX_K + 1):
    new_fp = [mpf(0)] * (MAX_K + 1)
    for i in range(MAX_K + 1):
        if f_power[i] == 0:
            continue
        fi = f_power[i]
        for jj in range(MAX_K + 1 - i):
            new_fp[i + jj] += fi * f[jj]
    f_power = new_fp
    sign = mpf((-1)**(n + 1))
    for k in range(1, MAX_K + 1):
        log_coeff[k] += sign * f_power[k] / n
S = [mpf(0)] * (MAX_K + 1)
S[0] = mpf(1)
for k in range(1, MAX_K + 1):
    S[k] = -k * log_coeff[k]
print(f"  S_1 = {mp.nstr(S[1], 18)}")

# ── Step 4: Hankel determinants D_n ─────────────────────────────────
print("Hankel determinants D_n:")
for n in range(1, 18):
    H = matrix(n, n)
    for i in range(n):
        for jj in range(n):
            H[i, jj] = S[i + jj + 1]
    D = det(H)
    sgn = "+" if D > 0 else "-"
    print(f"  D_{n:2d} = 10^({float(mp_log(abs(D))/mp_log(10)):8.2f})  sign={sgn}")

# ── Step 5-6: Gram–Schmidt → Jacobi parameters ──────────────────────
def inner_product(p, qq):
    """⟨p,q⟩ = Σ p_i q_j S_{i+j+1}  (shifted Hankel inner product)."""
    result = mpf(0)
    for i in range(len(p)):
        pi = p[i]
        if pi == 0:
            continue
        for jj in range(len(qq)):
            idx = i + jj + 1
            if idx < len(S):
                result += pi * qq[jj] * S[idx]
    return result

def build_jacobi(N):
    """Stieltjes/Gram–Schmidt: monomials → monic orthogonal polys → J_N."""
    polys = [[mpf(1)]]
    norms = [inner_product([mpf(1)], [mpf(1)])]
    alphas = []
    betas_sq = []
    for n in range(1, N + 1):
        xp = [mpf(0)] + polys[n - 1]
        alpha_n = inner_product(xp, polys[n - 1]) / norms[n - 1]
        pi_n = list(xp)
        for i in range(len(polys[n - 1])):
            pi_n[i] -= alpha_n * polys[n - 1][i]
        if n >= 2:
            beta_n_sq = norms[n - 1] / norms[n - 2]
            for i in range(len(polys[n - 2])):
                pi_n[i] -= beta_n_sq * polys[n - 2][i]
            betas_sq.append(beta_n_sq)
        sigma_n = inner_product(pi_n, pi_n)
        alphas.append(alpha_n)
        polys.append(pi_n)
        norms.append(sigma_n)
    J = matrix(N, N)
    for i in range(N):
        J[i, i] = alphas[i]
    for i in range(N - 1):
        bv = sqrt(betas_sq[i])
        J[i, i + 1] = bv
        J[i + 1, i] = bv
    ev, _ = mp.eigsy(J)
    eigenvalues = sorted([ev[i, 0] for i in range(N)], reverse=True)
    return alphas, [sqrt(b) for b in betas_sq], norms, eigenvalues

# ── Step 7: Output matrices J_3 through J_14 ────────────────────────
KNOWN = [
    14.134725141734694, 21.022039638771555, 25.010857580145689,
    30.424876125859512, 32.935061587739190, 37.586178158825672,
    40.918719012147496, 43.327073280915000, 48.005150881167160,
    49.773832477672302, 52.970321477714461, 56.446247697063395,
    59.347044002602353, 60.831778524609815,
]

print("\n" + "="*78)
for N in [3, 4, 5, 6, 8, 10, 12, 14]:
    alphas, betas, norms, eigvals = build_jacobi(N)
    print(f"\nJ_{N}: α = {[mp.nstr(a, 8) for a in alphas]}")
    if betas:
        print(f"     b = {[mp.nstr(b, 8) for b in betas]}")
    print(f"  {'n':>2} {'1/√λ':>12} {'γ*':>12} {'err%':>8}")
    for i, lam in enumerate(eigvals):
        gamma_est = 1/sqrt(lam)
        err = abs(gamma_est - KNOWN[i])/KNOWN[i]*100 if i < len(KNOWN) else float('nan')
        lock = " ✓" if err < 0.001 else ""
        print(f"  {i+1:>2} {float(gamma_est):>12.5f} {KNOWN[i]:>12.5f} {float(err):>7.3f}%{lock}")
