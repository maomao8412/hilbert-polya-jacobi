# GRH angular-monotonicity verification — 46 primitive Dirichlet characters

Numerical verification batch accompanying the GRH paper (Zenodo https://zenodo.org/records/22143035, CC BY 4.0).

Every result here is **forward-computed**: characters are generated programmatically
(Legendre symbols / roots of unity), completed L-functions are built from
`mpmath.dirichlet`, and the angular derivative `T_chi(r,theta)` is evaluated directly.
No zero data or fitted constants are used. Code in this repository is CC0 1.0.

## Result files

| File | Contents |
|------|----------|
| `grh_verify_output.txt` | Human-readable run log: 46 characters (q up to 101), c₂ constants, region checks |
| `grh_correct_results.json` | c₂ values, C/D envelope constants, grid results, tight q=3 point |
| `grh_Dr_verification.json` | T_chi negativity over r-grid per character (50-dps run, with timestamp) |
| `grh_xiK_angular.json` | ξ_K angular-monotonicity grid (full + reduced r ranges, theta samples) |
| `grh_matrix_primes_mod4.png` | Figure: matrix view of primes split mod 4 |

## Scripts (require `mpmath`)

| Script | Purpose |
|--------|---------|
| `grh_correct_verify.py` | Reference verification via Hurwitz zeta, 50-dps; characters generated programmatically |
| `grh_interval_verify.py` | Rigorous interval-arithmetic verification of `T_chi < 0` (regions A/B/C) |
| `grh_fast_verify.py` | Fast angular-monotonicity sweep, 30-dps |
| `grh_large_check.py` | Large-r and large-q spot checks |
| `grh_ultrafast.py` | Minimal point-density sanity sweep |

Run e.g. `python3 grh_correct_verify.py`. Outputs reproduced from scratch in minutes.
