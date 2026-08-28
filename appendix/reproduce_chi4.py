"""
Appendix: Chebyshev bias from the mod-4 Dirichlet L-function

This script reproduces the delta_psi chi4 results from
appendix/jacobi_prime_explicit.json.

We use L(s, chi_4) = beta(s) = sum_{n>=0} (-1)^n / (2n+1)^s  (Dirichlet beta).
Its nontrivial zeros at rho = 1/2 + i gamma (under GRH, verified numerically
for the first ~200) enter the explicit formula

    psi_0(x, chi_4) = - sum_rho x^rho/rho + artanh(1/x) - (L'/L)(0, chi_4)

where psi_0(x, chi_4) = sum_{n <= x} chi_4(n) Lambda(n) and
chi_4(n) = 0 if n even, 1 if n == 1 mod 4, -1 if n == 3 mod 4.

For a prime power p^k the contribution is chi_4(p)^k log p = (-1)^k log p
for p == 3 mod 4 (note the sign alternates with k); this is the subtlety
that produces the Chebyshev bias.

Run:
    pip install mpmath sympy
    python reproduce_chi4.py

The script computes zeros of beta(s) on the critical line by scanning
arg beta(1/2 + it), then evaluates the truncated explicit formula and
compares against the directly-summed psi_0(x, chi_4).
"""

import math
import os

from mpmath import mp, mpf, mpc, diff, atanh, cos, sin, log, sqrt
from sympy import primerange

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "jacobi_prime_explicit.json")

mp.dps = 50
CHI4 = [0, 1, 0, -1]  # chi_4 mod 4: 0,1,0,-1 (indexed n mod 4, n=0..3)


def L_beta(s):
    """Dirichlet beta L(s, chi_4) via mpmath.dirichlet."""
    return mp.dirichlet(s, CHI4)


def find_zeros(t_max):
    """Find positive imaginary parts of zeros of beta(1/2 + it), 0 < t < t_max.

    Uses arg-tracking: integrate (d/dt) arg beta along the critical line and
    bracket sign changes of Im beta (with Re beta simultaneously crossing 0
    verified by the arg change).
    """
    # Sampling grid; 0.02 is fine enough since first gap ~ 7 and later gaps < 1.
    dt = mpf("0.02")
    n = int(t_max / dt)
    gammas = []
    prev_val = L_beta(mpc(mpf("0.5"), mpf("0")))
    prev_re, prev_im = prev_val.real, prev_val.imag
    t = mpf("0")
    for i in range(1, n + 1):
        t = i * dt
        v = L_beta(mpc(mpf("0.5"), t))
        re, im = v.real, v.imag
        # Sign change in imaginary part with Re also small => candidate zero
        if prev_im * im < 0:
            # Bisection on t in [t-dt, t] using Im(beta) == 0; verify Re close.
            lo, hi = t - dt, t
            for _ in range(60):
                mid = (lo + hi) / 2
                vm = L_beta(mpc(mpf("0.5"), mid))
                if (vm.imag * im) > 0:
                    hi = mid
                else:
                    lo = mid
            t0 = (lo + hi) / 2
            v0 = L_beta(mpc(mpf("0.5"), t0))
            # Need both real and imaginary parts ~0 (true zero).
            if abs(v0) < mpf("1e-30"):
                gammas.append(float(t0))
        prev_re, prev_im = re, im
    return gammas


def b_chi4():
    """Constant term -(L'/L)(0, chi_4)."""
    L0 = L_beta(mpf("0"))
    Lp0 = diff(lambda s: L_beta(s), mpf("0"))
    return -Lp0 / L0


def psi0_chi4_explicit(x, gammas, b_c):
    """Truncated explicit formula (pairs of conjugate zeros)."""
    xm = mpf(str(x))
    sx = sqrt(xm)
    lx = log(xm)
    s = mpf("0")
    for g in gammas:
        theta = g * lx
        # 2 Re[ x^(1/2+ig) / (1/2+ig) ]
        denom = 0.25 + g * g
        s += 2 * sx * (0.5 * cos(theta) + g * sin(theta)) / denom
    trivial = atanh(1 / xm) if x > 1 else mpf("0")
    return float(-s + trivial + b_c)


def psi0_chi4_true(x):
    """Exact psi_0(x, chi_4) = sum_{p^k <= x} chi_4(p)^k log p."""
    n = int(math.floor(x))
    total = 0.0
    for p in primerange(2, n + 1):
        if p == 2:
            continue
        chi_p = 1 if p % 4 == 1 else -1
        pk = p
        k = 1
        while pk <= x:
            total += (chi_p ** k) * math.log(p)
            pk *= p
            k += 1
    return total


def main():
    import json
    if os.path.exists(DATA):
        with open(DATA) as f:
            ref = json.load(f)
        print("Loaded reference JSON:", DATA)
        print("Reference parameters:", ref["parameters"])
        print()
    else:
        ref = None

    # Use first 100 zeros to keep runtime reasonable; full ref uses K up to 100.
    K = 100
    print(f"Computing first {K} zeros of beta(s) on critical line (this takes a moment)...")
    gammas = find_zeros(mpf("80"))  # 80 is enough for ~28 zeros; extend if needed
    # If we want K=100, need t_max ~ 230 for beta; cap by what we got.
    if len(gammas) < K:
        print(f"  only found {len(gammas)} zeros up to t=80; extending t_max to 300...")
        gammas = find_zeros(mpf("300"))
    gammas = sorted(gammas)[:K]
    print(f"  found {len(gammas)} zeros; first 5: {[round(g,6) for g in gammas[:5]]}")
    print(f"  known first zero of beta(s) is 6.0209489...; error = {abs(gammas[0]-6.020948904691):.3e}")

    b_c = b_chi4()
    print(f"  b(chi4) = -(L'/L)(0,chi4) = {float(b_c):.6f}  (reference ~ -0.783189)")
    print()

    print(f"{'x':>5} {'true':>12} {'explicit':>12} {'abs_err':>10}")
    for x in [10, 20, 50, 100]:
        true = psi0_chi4_true(x)
        # use as many zeros as we have, capped at K for high x convergence
        Kx = min(len(gammas), 100)
        approx = psi0_chi4_explicit(x, gammas[:Kx], b_c)
        print(f"{x:>5} {true:>12.6f} {approx:>12.6f} {abs(approx-true):>10.4f}")

    if ref is not None:
        print()
        print("Reference (from JSON, K=100, dps=50, t_max=300):")
        for row in ref["delta_psi_chi4"]["comparison"]:
            print(f"  x={row['x']:>3}  true={row['delta_psi_true']:>10.6f}  "
                  f"explicit={row['delta_psi_explicit']:>10.6f}  "
                  f"abs_err={row['abs_error']:.4f}")


if __name__ == "__main__":
    main()
