# An explicit Hilbert–Pólya Jacobi matrix for the Riemann zeta function

A finite, real-symmetric, parameter-free tridiagonal matrix whose eigenvalues
converge to 1/γₙ², where γₙ are the imaginary parts of the nontrivial zeros of
the Riemann zeta function.

**Live post:** https://maomao8412.github.io/hilbert-polya-jacobi/

## Contents

| Path | Description |
|------|-------------|
| `index.html` | Self-contained web post (figure–table–text, academic style) |
| `assets/J10_numbered.png` | The 10×10 Jacobi matrix with actual numerical entries |
| `assets/spectral_convergence.png` | Eigenvalues of J_N (N=3…50) converging to Riemann zeros |
| `assets/prime_pulse.png` | Chebyshev ψ(x) reconstructed from J_50 eigenvalues |
| `assets/triptych_matrix_spectrum_primes.png` | One-image summary: matrix → spectrum → primes |
| `data/jacobi_N50_2000dps_result.json` | J_50 eigenvalues (2000-digit), α_n, b_n |
| `reproduce.py` | ~60-line script that loads the matrix and reproduces zeros + ψ(x) |

## Construction (four steps)

1. Expand the completed xi function at s = 1/2:
   ξ(s) = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s) = Σ σ_{2m} (s − 1/2)^{2m}.
2. Newton identities recover power sums P_k = Σ γₙ^{−2k}.
3. Form the Hankel matrix H_{ij} = P_{i+j+1}, Cholesky-factor H = LLᵀ.
4. Conjugate the shifted Hankel H⁽¹⁾:
   **J_N = L⁻¹ H⁽¹⁾ L^{−T}** is real, symmetric, tridiagonal, with positive
   subdiagonal — a finite self-adjoint operator.

The eigenvalues λₙ of J_N satisfy λₙ → 1/γₙ² as N → ∞
(Gauss–Christoffel quadrature theorem; RH is **not** assumed).

## Quick start

```bash
pip install sympy
python reproduce.py
```

Expected output (first rows):

```
 n        gamma_matrix          gamma_known       |err|
 1    14.134725141734695    14.134725141734695    0.00e+00
 2    21.022039638771556    21.022039638771556    0.00e+00
 3    25.010857580145690    25.010857580145688    1.43e-14
 ...

    x    psi_exact     psi_J50     rel_err
   10       7.8320       8.0116       2.29%
   20      19.2657      19.2992       0.17%
   50      49.4854      49.1303       0.72%
  100      94.0453      94.9182       0.93%
```

## Reproducibility of the matrix itself

The 50 eigenvalues in `data/jacobi_N50_2000dps_result.json` were computed with
[mpmath](https://mpmath.org/) at 2000-digit working precision, using only:

* Taylor coefficients σ_{2m} of ξ(s) at s = 1/2 (closed form via functional
  equation),
* Newton identities (polynomial arithmetic),
* Cholesky decomposition of a positive-definite Hankel matrix,
* symmetric tridiagonal eigenvalue computation.

No zero-finding, no prime table, no external data enter the construction.
The same construction applies to completed Dirichlet L-functions
(Λ(s,χ) = (q/π)^{(s+a)/2} Γ((s+a)/2) L(s,χ)); see the GRH paper below.

## Papers

* **RH operator construction (55 pp.)**
  Zenodo: https://zenodo.org/records/22113226
  DOI: 10.5281/zenodo.22113226

* **GRH verification — 46 primitive characters q ≤ 20 (36 pp.)**
  Zenodo: https://zenodo.org/records/22143035
  DOI: 10.5281/zenodo.22143035

## Scope and limitations

* This is a sequence of **finite** matrices converging to the zero spectrum.
  The limit defines a natural Jacobi operator on ℓ²; full spectral analysis of
  the infinite operator is ongoing.
* Hankel positive definiteness is verified numerically to N = 21 (350 digits)
  but is not claimed as an analytic theorem for all N.
* RH itself is **not assumed**; convergence of eigenvalues to 1/γₙ² follows
  from the unconditional Hadamard product. Proving RH requires showing all
  zeros lie on the critical line, which this matrix construction alone does
  not establish.

## Author

Zhuo Chen (陈倬) · ORCID [0009-0006-9172-8268](https://orcid.org/0009-0006-9172-8268)

## License

Data and code: [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
The accompanying preprints retain their authors' rights as posted on Zenodo.
