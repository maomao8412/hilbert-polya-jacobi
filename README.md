# An explicit Hilbert–Pólya Jacobi matrix

**Constructed from the three-term polylogarithm decomposition of ζ(s)**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22113226.svg)](https://zenodo.org/records/22113226)

An explicit, finite, parameter-free real-symmetric tridiagonal Jacobi matrix whose eigenvalues converge monotonically to 1/γ_n², where γ_n are the imaginary parts of the non-trivial zeros of ζ(s).

- **English page**: https://maomao8412.github.io/hilbert-polya-jacobi/
- **中文页面**: https://maomao8412.github.io/hilbert-polya-jacobi/zh.html

## Construction chain (strictly per paper)

1. **Three-term polylogarithm identity**

   ζ(s) = Li_s(z₁) + Li_s(z₂) + Σ c_n / n^s,

   where z₁ = 2−√2, z₂ = √2−1, z₁+z₂ = 1, and c_n = 1−z₁^n−z₂^n satisfies a second-order recurrence. This is the analytic foundation; without it, Hankel positivity has no proof.

2. **ξ Taylor coefficients and log-moments** — Cauchy integral gives a_{2k}; define S_k = −k·d_{2k} from the Mercator expansion of log(ξ(½+it)/ξ(½)). S_k is unconditional.

3. **Hankel matrix and positivity** — H_N[i,j] = S_{i+j+1}; D_n = det(H_n) > 0 is **proved** from the three-term identity via Herglotz/Stieltjes theory, not assumed.

4. **Gram–Schmidt → Jacobi matrix** — Orthogonal polynomials in L²(μ) give α_n, b_n; J = tridiag(b_n, α_n, b_n). Entries are rational functions of S₁,…,S_{2n+2}; no zero-finding, no prime table.

## Reproducibility

```bash
pip install mpmath
python reproduce_paper.py
```

The script (~150 lines, no external data) reproduces:
- ξ(½), S_k, and Hankel determinants D_n
- Jacobi parameters α_n, b_n for J₃ through J₁₄
- Eigenvalues converging to 1/γ_n² (first five zeros locked at N=14)

All values match the paper tables to every printed digit.

## Files

| File | Description |
|------|-------------|
| `index.html` | English page (MathJax-rendered formulas) |
| `zh.html` | Chinese page (MathJax-rendered formulas) |
| `reproduce_paper.py` | Self-contained reproduction script (mpmath, 120-digit) |
| `data/jacobi_N50_2000dps_result.json` | J₅₀ at 2000-digit precision |
| `assets/` | Figures: matrix, convergence, primes |
| `appendix/` | GRH extension (χ₄) reproduction |

## Papers

- **RH operator construction** (55 pp.): https://zenodo.org/records/22113226
- **GRH verification** (46 primitive characters, 36 pp.): https://zenodo.org/records/22143035

## License

MIT
