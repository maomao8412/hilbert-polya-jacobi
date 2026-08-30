# Controls — discriminative checks for the benchmark pipeline

## e4_eisenstein/ — off-critical-line control (Eisenstein series E₄)

**Purpose.** Show that the L-function → moments → Jacobi pipeline *discriminates*
between zeros lying on the symmetry line and zeros lying off it. The four
benchmark Hamiltonians (ζ, β, Δ, product) all have their non-trivial zeros on
the critical line and yield Jacobi chains with b² > 0 through dozens–hundreds of
orders. This control deliberately uses an L-function whose zeros lie on **two
lines that are not the symmetry line**, and the chain collapses at low order —
exactly as the moment theory predicts.

**Object.** The Eisenstein series of weight 4,

- L(E₄, s) = ζ(s) ζ(s−3), Dirichlet coefficients σ₃(n) = Σ_{d|n} d³ (elementary, deterministic);
- completed Λ(s) = (2π)⁻ˢ Γ(s) ζ(s) ζ(s−3), functional equation Λ(s) = Λ(4−s), symmetry center s = 2;
- Λ has genuine poles at s = 0, 4; the entire function sampled is **Ξ_E₄(s) = s(s−4) Λ(s)** (s(s−4) invariant under s ↔ 4−s).

**Zeros.** Non-trivial zeros at s = ½ ± iγ_k and s = 7/2 ± iγ_k (γ_k = Riemann
zero heights) — two lines, both offset by 3/2 from the symmetry center. In the
t-variable s = 2 + it they are t = ±γ_k ∓ 3i/2: the moment measure is
**complex-supported**, not a positive measure on the real axis.

**Pipeline (identical to the benchmark families).** Circle sampling |s−2| = 2,
M = 1536 points, 650-digit precision; sampling grid shifted by half a step
(+π/M) so no sample sits at the canceled pole positions s = 0, 4. Taylor/log
moment recursion and Gram–Schmidt are line-for-line the same code.

**Result (650dps run, `e4_control_results.json`).**

- σ₀ = Ξ_E₄(2) = 1/72 = 0.013888…, reproduced to ~15 digits (imaginary part ~1e-654).
- Moments P[1]…P[7] positive; **P[8] negative**, P[50], P[100] negative —
  predicted: P[m] ≈ 2·Σ Re[(γ_k − 3i/2)^(−2m)], dominated at high m by k = 1
  with factor cos(2m·θ₁), θ₁ = arctan(3/(2γ₁)) ≈ 0.106, whose sign first
  changes at m ≈ π/(4θ₁) ≈ 7.4. Observed first negative at m = 8.
- Jacobi b² sequence: positive at order 2 only, then signs + − − + − − ++ − …,
  strictly real (imaginary part 0); **first non-positive b² at k+1 = 4**; the
  positive-prefix chain has size 2.

**Interpretation.** On-critical-line L-functions give a real positive moment
measure → Hankel moment matrices positive definite → long Jacobi chains
(b² > 0 through J = 22–100 in this benchmark). An off-line L-function gives a
complex-supported measure whose real projections cease to be positive definite
at low order, and the Gram–Schmidt chain breaks at J ≈ 4. The pipeline is
therefore not a generic number-producer: its central positive-definiteness
signal is tied to where the zeros sit. This is a discriminative control for the
numerical study, not a claim about any family's non-trivial zeros.

**Scope.** Finite-truncation numerical experiment at finite precision; the
collapse order is precision- and M-dependent and should be read as a
qualitative discriminator, not an invariant.
