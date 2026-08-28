# Appendix — Chebyshev bias from the mod-4 Dirichlet L-function

This appendix extends the main construction to the completed Dirichlet
L-function

    Λ(s,χ₄) = (4/π)^((s+1)/2) Γ((s+1)/2) L(s,χ₄),

where χ₄ is the primitive character mod 4 (Dirichlet beta function
β(s) = L(s,χ₄)). The Newton–Hankel–Cholesky chain applied to
ξ_K(s) = ξ(s)Λ(s,χ₄) yields a real symmetric tridiagonal Jacobi matrix
whose eigenvalues carry both the Riemann zeros and the L(s,χ₄) zeros.

The L(s,χ₄) zeros enter the explicit formula

    ψ₀(x,χ₄) = − Σ_ρ x^ρ/ρ + artanh(1/x) − (L'/L)(0,χ₄),

where ψ₀(x,χ₄) = Σ_{n≤x} χ₄(n) Λ(n). For a prime power p^k the
contribution is χ₄(p)^k log p = (−1)^k log p for p ≡ 3 (mod 4); the sign
alternates with k. ψ₀(x,χ₄) is negative for most small x — this is the
**Chebyshev bias** (more primes ≡ 3 mod 4 than ≡ 1 mod 4).

## Files

| File | Description |
|------|-------------|
| `jacobi_prime_explicit.json` | Full numerical output (N=40, dps=50, t_max=300) |
| `jacobi_prime_explicit.png` | Two-panel figure: ψ(x) comparison and truncation error |
| `reproduce_chi4.py` | Standalone script that finds β(s) zeros and verifies the explicit formula |

## Reference results (N=40, dps=50, t_max=300)

Power-sum chain verification (m = 1..8):

| m | P_m direct | P_m from Jacobi | rel. error |
|---|-----------|-----------------|------------|
| 1 | 2.0516e-02 | 2.0516e-02 | 0 |
| 2 | 3.7164e-05 | 3.7164e-05 | 0 |
| 3 | 1.4417e-07 | 1.4417e-07 | 1.3e-50 |
| 6 | 1.5877e-14 | 1.5877e-14 | 0 |

Taylor-bridge verification: σ_{2m} recovered by Newton identities from
P_m agrees with direct DFT computation to 9.4e-17 relative. Chebyshev
b_n (b₁–b₆) from the Lanczos and Chebyshev algorithms agree to
10⁻⁴¹–10⁻⁵⁰.

ψ(x) from ξ_K zeros at K=137:

| x | exact | explicit | rel. error |
|---|-------|----------|-----------|
| 10 | 7.8320 | 7.7686 | 0.81% |
| 20 | 19.2657 | 19.3068 | 0.21% |
| 50 | 49.4854 | 49.3367 | 0.30% |
| 100 | 94.0453 | 93.7046 | 0.36% |
| 200 | 206.1459 | 205.7753 | 0.18% |

ψ₀(x,χ₄) at K=100 L-zeros (Chebyshev bias):

| x | true | explicit | abs. error |
|---|------|----------|-----------|
| 10 | −0.3365 | −0.3274 | 0.009 |
| 20 | −0.2806 | −0.3702 | 0.090 |
| 50 | −1.3130 | −1.1039 | 0.209 |
| 100 | −0.1126 | 0.1099 | 0.222 |
| 200 | 1.4480 | 3.0347 | 1.587 |

The error at x = 200 is larger because K = 100 L-zeros is not enough
to resolve oscillations at that range; it decreases monotonically with
more zeros (same Gibbs phenomenon as the ζ case). The negative sign of
ψ₀ for x ≤ 100 is the Chebyshev bias.

## Constant term

    b(χ₄) = −(L'/L)(0,χ₄) ≈ −0.783189

computed numerically by mpmath.diff at s=0. The trivial zeros of
L(s,χ₄) lie at s = −1, −3, −5, … and contribute
Σ x^{−(2k+1)}/(2k+1) = artanh(1/x).

## Scope

The β(s) zeros in `reproduce_chi4.py` are found by numerical
root-finding on the critical line. This is a numerical verification of
the explicit formula, not a proof that all β(s) zeros lie on the line
(that is the GRH for χ₄, addressed by the same Hankel positive
definiteness argument in the GRH paper, DOI 10.5281/zenodo.22143035).
