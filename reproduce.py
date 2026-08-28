"""
Reproduce the Riemann-zero eigenvalues and prime-number reconstruction
from the explicit Hilbert–Polya Jacobi matrix J_50.

Data file: data/jacobi_N50_2000dps_result.json
  - eigenvalues_desc: list of 50 eigenvalues lambda_n (descending), each a
                      2000-digit-precision decimal string, satisfying
                      lambda_n -> 1/gamma_n^2, where gamma_n are the
                      imaginary parts of the nontrivial zeros of zeta(s).
  - alphas: diagonal entries alpha_0 .. alpha_49
  - betas : subdiagonal entries b_0 .. b_48   (b_n > 0)

The matrix is
    J = diag(alpha) + diag(beta, +1) + diag(beta, -1).

Usage:
    python reproduce.py
"""

import json
import math
import os
from sympy import primerange

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "jacobi_N50_2000dps_result.json")


def load():
    with open(DATA) as f:
        d = json.load(f)
    eigs = [float(x) for x in d["eigenvalues_desc"]]   # already descending
    return eigs


def gammas_from_eigenvalues(eigs):
    """lambda_n = 1/gamma_n^2, lambda sorted descending
    => gamma sorted ascending."""
    return sorted(1.0 / math.sqrt(lam) for lam in eigs)


def psi_explicit(x, gammas, K=None):
    """Riemann explicit formula truncated at K zeros."""
    if K is None:
        K = len(gammas)
    s = 0.0
    logx = math.log(x)
    sx = math.sqrt(x)
    for g in gammas[:K]:
        th = g * logx
        # pair rho, 1-rho: 2 Re[ x^{1/2+ig}/(1/2+ig) ]
        s += 2.0 * sx * (0.5 * math.cos(th) + g * math.sin(th)) / (0.25 + g * g)
    return x - math.log(2 * math.pi) - 0.5 * math.log1p(-1.0 / (x * x)) - s


def psi_true(x):
    """Exact Chebyshev psi(x) = sum_{p^k <= x} log p."""
    n = int(math.floor(x))
    total = 0.0
    for p in primerange(2, n + 1):
        pk = p
        while pk <= x:
            total += math.log(p)
            pk *= p
    return total


def main():
    eigs = load()
    gammas = gammas_from_eigenvalues(eigs)

    known = [14.134725141734695, 21.022039638771556, 25.010857580145688,
             30.424876125859513, 32.935061587739180, 37.586178158825670,
             40.918719012147500, 43.327073280915000]

    print("First 8 zeros from J_50:")
    print(f"{'n':>2}  {'gamma_matrix':>20}  {'gamma_known':>20}  {'|err|':>10}")
    for i, g in enumerate(gammas[:8]):
        err = abs(g - known[i])
        print(f"{i+1:>2}  {g:>20.15f}  {known[i]:>20.15f}  {err:.2e}")

    print()
    print("Chebyshev psi(x) from J_50 eigenvalues vs exact:")
    print(f"{'x':>5}  {'psi_exact':>12}  {'psi_J50':>12}  {'rel_err':>10}")
    for x in [10, 20, 50, 100]:
        exact = psi_true(x)
        approx = psi_explicit(x, gammas, K=50)
        print(f"{x:>5}  {exact:>12.4f}  {approx:>12.4f}  {abs(approx-exact)/abs(exact)*100:>9.2f}%")


if __name__ == "__main__":
    main()
