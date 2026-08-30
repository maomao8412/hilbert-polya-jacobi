# Quantum Arithmetic Benchmark Dataset — v1

**Arithmetic-origin one-dimensional tight-binding Hamiltonians for quantum-simulation
and spectral-statistics benchmarks.**

Four real-symmetric tridiagonal (Jacobi) matrices. Every matrix entry is derived,
by a deterministic numerical pipeline, purely from the arithmetic Dirichlet
coefficients of an L-function (Riemann zeta, Dirichlet character χ₄, Ramanujan τ
cusp form, and their product). **No zeros, primes, or any spectral information
enter the construction** — the matrices are built before and independently of any
knowledge of the L-function zeros. The zeros are used only afterwards, as an
independent check that the low-lying eigenvalues reproduce them.

Reference: Chen, Z. *Quantum Chaos of Arithmetic Origin*, Zenodo
[10.5281/zenodo.22180921](https://doi.org/10.5281/zenodo.22180921) (2026).
Companion construction papers: [10.5281/zenodo.22167882](https://doi.org/10.5281/zenodo.22167882),
[10.5281/zenodo.22167742](https://doi.org/10.5281/zenodo.22167742).

## What the dataset is for

- **Quantum simulation target / turnkey benchmark:** a single-chain
  nearest-neighbor Hamiltonian with no fine-tuned parameters — site energies and
  hoppings are arithmetic coefficients — whose finite-truncation spectrum shows
  Wigner–Dyson level statistics and localized wavepackets, i.e. quantum-chaotic
  signatures arising *without* any classical chaotic limit and *without* random
  disorder. Proposed as a calibration task for superconducting-qubit / ion-trap
  quantum simulators (see paper Appendix D for the predicted pass/fail protocol).
- **Spectral-statistics test suite:** Wigner–Dyson vs Poisson discrimination,
  level-rigidity, gap-ratio distributions, participation ratios, and noise
  threshold ε_c ≈ 0.2–0.3 (site-energy randomization destroys the signal at
  ~20–30% relative noise, three families consistently).

## Contents

| directory | L-function family | size N | locked levels (rel err < 1e-4) |
|---|---|---|---|
| `zeta_J50/` | Riemann ζ(s) | 50 | 25 |
| `beta_J50/` | Dirichlet L mod 4 (χ₄, Catalan) | 50 | 27 |
| `delta_J22/` | Ramanujan τ, weight-12 cusp form Δ | 22 | 15 (small sample) |
| `product_J100/` | ξ·Λ_β product (independent families) | 100 | 53 (35 β + 18 ζ) |

Each directory contains:

- `hamiltonian.csv` — dense N×N tridiagonal matrix, `H[i][i]=alpha_{i+1}`,
  `H[i][i+1]=H[i+1][i]=b_{i+1}`, all other entries zero. Dimensionless units.
  Eigenvalues satisfy λ_k ≈ 1/γ_k², where γ_k are the corresponding L-function
  zero heights on the critical line.
- `metadata.json` — family name, Dirichlet coefficients, pipeline parameters
  (contour radius R, sample count M, working precision), locked-level count,
  observed spectral statistics, and construction description.

`summary_table.csv` lists all four Hamiltonians; `loader_demo.py` loads each
matrix and prints basic spectral diagnostics.

## Construction pipeline (deterministic, reproducible)

1. Take the completed entire L-function Ξ(s) (functional-equation symmetric).
2. Sample Ξ on a circle |s − s₀| = R (R = 2; R = 4 for Δ) at M = 1024–1536
   points, high working precision (230–1300 decimal digits).
3. Extract Taylor coefficients via discrete Fourier; recursively obtain the
   moment sequence P_m from log-derivative coefficients.
4. Gram–Schmidt orthonormalization of monomials under the moment functional
   yields Jacobi coefficients α_n (diagonal) and b_n (off-diagonal).
5. Truncate at the first non-positive b² (Hankel conditioning floor); the
   locked subspace is the set of eigenvalues whose γ = 1/√λ matches an
   independent L-function zero to relative error < 1e-4.

Full pseudocode, parameter tables, seeds and reproduction scripts are in the
paper (Appendix A–D) and its reproducibility bundle on Zenodo.

## Scope and limitations (please read)

- All results are **finite-truncation numerical evidence**. Nothing here proves
  any statement about infinite-dimensional operators.
- The spectra are **computer-generated**, not measured; this is a simulation
  benchmark. Physical-simulation experiments are proposed but not yet performed.
- "Quantum chaos of arithmetic origin" is a **physical hypothesis** supported by
  the data, not an established fact.
- Δ reaches only 15 locked levels at 230-digit precision (Hankel conditioning
  floor); its statistics are reported with an explicit small-sample caveat.
- Gap-ratio estimates from 15–53 levels have a known upward small-sample bias
  (Atas et al. 2013); the dataset does not claim GUE vs GOE discrimination.

## License

CC BY 4.0. Cite as: Chen, Z. (2026), *Quantum Chaos of Arithmetic Origin*,
Zenodo 10.5281/zenodo.22180921; code and data:
github.com/maomao8412/hilbert-polya-jacobi.
