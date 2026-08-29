# GRH companion: 50×50 Jacobi matrix for the Dirichlet beta function (L(s, χ₄))

Forward construction of the Hilbert–Pólya Jacobi matrix for the completed
Dirichlet beta function Λβ(s) = (π/4)^(-(s+1)/2) Γ((s+1)/2) β(s),
β(s) = L(s, χ₄) = 4^-s [ζ(s,¼) − ζ(s,¾)].

**No zero data enter the construction.** The only input is the completed
function; the zeros are the output (matrix spectrum → 1/√λₙ → β zeros).

## Pipeline

1. `sample600.py` — 768 samples of Λβ(½+z)/Λβ(½) on |z|=0.3, sympy 600-digit
   (this environment's mpmath is a double-only build; all high-precision work
   uses sympy). Produces `samples600.json` (not included; ~25 min to regenerate).
2. `full600.py` — even-index DFT (Cauchy integral) → Taylor coefficients d_{2j};
   log-series recurrence n c_n = n d_n − Σ_{k<n} k c_k d_{n−k} → moments
   T_m = Σₙ γₙ^(−2m); Stieltjes monic-orthogonal-polynomial recurrence →
   Jacobi αₙ (diagonal) and bₙ (off-diagonal); numpy eigenvalues; independent
   mpmath double-precision zero scan for comparison.
   Outputs: `beta_J50_results.json`, `beta_J50_matrix.csv`.

## Results (50×50, 600-digit forward computation)

- First 27 zeros locked: |1/√λₙ − γₙ|/γₙ < 10⁻⁴ for n = 1…27
  (first 24 to ~10⁻¹⁶, i.e. double precision); γ₁ = 6.020948904697597…
- All b²ₙ positive through n=50; smallest norm σ₅₀ ≈ 10⁻³⁶⁶, comfortably
  above the 10⁻⁶⁰⁰ numerical noise floor.
- `grh_web_data.json` — 46-primitive-character angular-monotonicity summary
  (10 representative characters: 798/95-point r–θ grids, all max T < 0;
  C(r)/D(r) envelopes; tight-region scan max T = −0.2057 over 790 points).

## Web page

`gen_svgs.py` + `assemble.py` build the interactive pages
(`grh_zh.html` / `grh.html` at repository root): lightbox-zoomable matrix and
spectrum SVGs with in-cell values, lock table, and an interactive per-character
angular-monotonicity heatmap.

Code: CC0 1.0. Paper: CC BY 4.0 (Zenodo 10.5281/zenodo.22143035).
