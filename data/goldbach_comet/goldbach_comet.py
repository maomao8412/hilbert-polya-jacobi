#!/usr/bin/env python3
"""
Goldbach Comet — Three-Term Decomposition Visualization
R(N) = Σ_{i,j∈{1,2,3}} R_{ij}(N)

Correct formula: r(N) ~ N · (N) / (log N)²
where 𝔖(N) = C₂ · Π_{p|N, p>2} (p-1)/(p-2)
      C₂ = 0.66016... (twin prime constant)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ─── Configuration ───
MAX_N = 10000
OUT_DIR = '/Coze/Drive/黎曼猜想论文审核/所有对话/主对话/goldbach_comet_v2'

# ─── Font setup ───
plt.rcParams['font.sans-serif'] = ['WenQuanYi Micro Hei', 'Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ═══════════════════════════════════════════
# Step 1: Sieve
# ══════════════════════════════════════════
print("Step 1: Computing primes and Λ(n)...")
sieve = np.ones(MAX_N + 1, dtype=bool)
sieve[0] = sieve[1] = False
for i in range(2, int(MAX_N**0.5) + 1):
    if sieve[i]:
        sieve[i*i::i] = False
primes = np.where(sieve)[0]
is_prime = sieve.copy()

# von Mangoldt Λ(n)
Lambda = np.zeros(MAX_N + 1)
for p in primes:
    pk = p
    log_p = np.log(p)
    while pk <= MAX_N:
        Lambda[pk] += log_p
        pk *= p

# ═══════════════════════════════════════════
# Step 2: Goldbach prime pair count r(N)
# ═══════════════════════════════════════════
print("Step 2: Computing r(N) — Goldbach prime pair counts...")
even_nums = np.arange(4, MAX_N + 1, 2)
N_count = len(even_nums)

r_N = np.zeros(N_count, dtype=int)
prime_arr = primes[primes <= MAX_N // 2]

for p in prime_arr:
    start_N = max(4, 2 * p)
    start_idx = (start_N - 4) // 2
    q_vals = np.arange(start_N - p, MAX_N - p + 1, 2)
    r_N[start_idx:start_idx + len(q_vals)] += is_prime[q_vals]

print(f"  r(N) range: {r_N.min()} to {r_N.max()}, mean={r_N.mean():.1f}")

# ═══════════════════════════════════════════
# Step 3: Λ-convolution R(N) via FFT
# ═══════════════════════════════════════════
print("Step 3: Computing Λ-convolution R(N) via FFT...")
from numpy.fft import fft, ifft

fft_size = 1
while fft_size < 2 * MAX_N + 1:
    fft_size *= 2

L_fft = fft(Lambda[:MAX_N + 1], fft_size)
conv_result = np.real(ifft(L_fft * L_fft))
R_conv = np.array([conv_result[N] for N in even_nums])

# ═══════════════════════════════════════════
# Step 4: Singular series & Hardy-Littlewood main term
# ═══════════════════════════════════════════
print("Step 4: Computing singular series and HL prediction...")

C2 = 0.6601618158468696  # Twin prime constant

def singular_series(N):
    """𝔖(N) = C₂ · Π_{p|N, p>2} (p-1)/(p-2)"""
    result = C2
    temp = N
    for p in primes:
        if p == 2:
            continue
        if p * p > N:
            break
        if temp % p == 0:
            result *= (p - 1.0) / (p - 2.0)
            while temp % p == 0:
                temp //= p
    if temp > 2:
        result *= (temp - 1.0) / (temp - 2.0)
    return result

S_N = np.array([singular_series(N) for N in even_nums])
log_N = np.log(even_nums.astype(float))

# Hardy-Littlewood main term for prime pairs
HL_main = even_nums.astype(float) * S_N / (log_N ** 2)

print(f"  S(N) range: {S_N.min():.4f} to {S_N.max():.4f}")
print(f"  HL main range: {HL_main.min():.2f} to {HL_main.max():.2f}")

# ═══════════════════════════════════════════
# Step 5: Residuals
# ═══════════════════════════════════════════
residual = r_N.astype(float) - HL_main
sqrt_N = np.sqrt(even_nums.astype(float))

# Error bounds for visualization:
# Classical heuristic: ~C·N^{1/2}/log N (random model)
# GRH bound: ~C·N^{1/2}·log N
# Three-term decomposition (RH cancellation): ~C·log²N

bound_classical = 3.0 * sqrt_N / log_N
bound_GRH = 0.3 * sqrt_N * log_N
bound_3term = 0.012 * log_N ** 2

# ═══════════════════════════════════════════
# Step 6: Oscillatory cross-term from zeta zeros
# ═══════════════════════════════════════════
print("Step 5: Computing oscillatory correction from zeta zeros...")

known_zeros = [
    14.134725141734693790, 21.022039638771554993, 25.010857580145688763,
    30.424876125859513210, 32.935061587739189691, 37.586178158825671253,
    40.918719012147495186, 43.327073280914999520, 48.005150881167159728,
    49.773832477672302182, 52.970321477720275208, 56.446247697170741326,
    59.347044003098467894, 60.831778524617904771, 65.112544048081657133,
    67.079810744863853923, 69.546401147890683040, 72.067157543049867898,
    75.704695751674751940, 77.144840066965187048, 79.337375955849996960,
    82.910479248297098361, 84.735492980499525473, 87.425274615685067320,
    88.809111207633520973, 92.491899271382264933, 94.651344040519881060,
    95.870634228260375806, 98.831194382568689848, 101.31785100572819996,
    103.72553792788274808, 105.44662305409004984, 107.16861118367490228,
    111.02953554316994199, 111.87465917599859894, 114.32022091589620847,
]

K = 36
oscillatory = np.zeros(N_count)
for idx, N in enumerate(even_nums):
    logN = np.log(float(N))
    val = 0.0
    for g in known_zeros[:K]:
        val += np.cos(g * logN) / g
    oscillatory[idx] = val * even_nums[idx] / (logN ** 2)

# ═══════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════
print("Step 6: Generating plots...")

colors_map = {0: '#e74c3c', 2: '#3498db', 4: '#2ecc71'}
labels_map = {0: 'N ≡ 0 (mod 6)', 2: 'N ≡ 2 (mod 6)', 4: 'N ≡ 4 (mod 6)'}

def scatter_by_mod6(ax, x, y, s=3, alpha=0.6):
    for mod in [0, 2, 4]:
        mask = (x % 6) == mod
        ax.scatter(x[mask], y[mask], s=s, alpha=alpha, c=colors_map[mod],
                   label=labels_map[mod], edgecolors='none')

# ═══════════════════════════════════════════
# Plot 1: Classic Goldbach Comet
# ═══════════════════════════════════════════
print("  Plot 1: Classic Goldbach Comet...")
fig1, ax1 = plt.subplots(figsize=(10, 7))
scatter_by_mod6(ax1, even_nums, r_N, s=5, alpha=0.6)

ax1.plot(even_nums, HL_main, '-', color='#2c3e50', linewidth=1.5,
         alpha=0.7, label='Hardy-Littlewood N·𝔖(N)/(log N)²')

ax1.set_xlabel('N', fontsize=13)
ax1.set_ylabel('r(N) = #{p+q=N : p,q prime}', fontsize=13)
ax1.set_title('Classic Goldbach Comet\nThree-Arm Structure by N mod 6', fontsize=14, fontweight='bold')
ax1.legend(loc='upper left', fontsize=10, markerscale=3)
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig(f'{OUT_DIR}/comet_classic.png', dpi=150, bbox_inches='tight')
plt.close(fig1)
print("    → comet_classic.png saved")

# ═══════════════════════════════════════════
# Plot 2: Main term skeleton vs actual
# ═══════════════════════════════════════════
print("  Plot 2: Main term skeleton...")
fig2, ax2 = plt.subplots(figsize=(10, 7))

ax2.scatter(even_nums, r_N, s=4, alpha=0.4, c='#e74c3c',
            label='r(N) actual prime pairs', edgecolors='none')

ax2.plot(even_nums, HL_main, '-', color='#2980b9', linewidth=2.0,
         alpha=0.9, label='Main term: N·𝔖(N)/(log N)²')

ax2.plot(even_nums, HL_main + oscillatory, '--', color='#27ae60', linewidth=1.2,
         alpha=0.7, label='+ Oscillatory (Λ₁Λ₂ cross term)')

ax2.fill_between(even_nums, HL_main - bound_classical, HL_main + bound_classical,
                 alpha=0.12, color='#3498db', label='±3√N/log N heuristic band')

ax2.set_xlabel('N', fontsize=13)
ax2.set_ylabel('Count', fontsize=13)
ax2.set_title(u'Main Arc Skeleton: Λ₁Λ₁ vs Actual r(N)\nR₁₁ Main Term + Λ₁Λ₂ Oscillatory from Zeta Zeros',
              fontsize=13, fontweight='bold')
ax2.legend(loc='upper left', fontsize=9, markerscale=2)
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig(f'{OUT_DIR}/comet_main_term.png', dpi=150, bbox_inches='tight')
plt.close(fig2)
print("    → comet_main_term.png saved")

# ═══════════════════════════════════════════
# Plot 3: Residual plot with error bounds
# ═══════════════════════════════════════════
print("  Plot 3: Residual decomposition...")
fig3, ax3 = plt.subplots(figsize=(10, 7))

scatter_by_mod6(ax3, even_nums, residual, s=5, alpha=0.7)

ax3.plot(even_nums, bound_classical, '--', color='#e67e22', linewidth=2.0, alpha=0.85,
         label='Classical ~ C·√N/log N')
ax3.plot(even_nums, -bound_classical, '--', color='#e67e22', linewidth=2.0, alpha=0.85)

ax3.plot(even_nums, bound_GRH, '-.', color='#9b59b6', linewidth=1.8, alpha=0.8,
         label='GRH-type ~ C·√N·log N')
ax3.plot(even_nums, -bound_GRH, '-.', color='#9b59b6', linewidth=1.8, alpha=0.8)

ax3.plot(even_nums, bound_3term, ':', color='#27ae60', linewidth=2.5, alpha=0.9,
         label='Three-term ~ C·log²N (RH)')
ax3.plot(even_nums, -bound_3term, ':', color='#27ae60', linewidth=2.5, alpha=0.9)

ax3.axhline(y=0, color='black', linewidth=0.5, alpha=0.5)
ax3.set_xlabel('N', fontsize=13)
ax3.set_ylabel('r(N) − N·𝔖(N)/(log N)²', fontsize=12)
ax3.set_title(u'Residual Analysis: Classical vs Three-Term Decomposition Bounds\nΛ = Λ₁ + Λ₂ + Λ₃  →  R = Σᵢ Rᵢⱼ',
              fontsize=12, fontweight='bold')
ax3.legend(loc='upper left', fontsize=9, markerscale=2)
ax3.grid(True, alpha=0.3)
fig3.tight_layout()
fig3.savefig(f'{OUT_DIR}/comet_residual.png', dpi=150, bbox_inches='tight')
plt.close(fig3)
print("    → comet_residual.png saved")

# ═══════════════════════════════════════════
# Plot 4: Decay exponent log-log
# ═══════════════════════════════════════════
print("  Plot 4: Decay exponent...")
fig4, ax4 = plt.subplots(figsize=(10, 7))

abs_residual = np.abs(residual)
valid = abs_residual > 0.5
logN_vals = np.log(even_nums[valid].astype(float))
logR_vals = np.log(abs_residual[valid])

for mod in [0, 2, 4]:
    mask = (even_nums[valid] % 6) == mod
    if mask.any():
        ax4.scatter(logN_vals[mask], logR_vals[mask], s=5, alpha=0.5,
                    c=colors_map[mod], label=labels_map[mod], edgecolors='none')

coeffs = np.polyfit(logN_vals, logR_vals, 1)
fit_line = np.polyval(coeffs, logN_vals)
sort_idx = np.argsort(logN_vals)
ax4.plot(logN_vals[sort_idx], fit_line[sort_idx], 'k-', linewidth=2.5,
         label=f'Linear fit: slope α = {coeffs[0]:.3f}')

for exp, color, lbl in [(0.5, '#e67e22', 'α=0.5  (√N)'),
                         (0.0, '#27ae60', 'α=0.0  (const)'),
                         (1.0, '#9b59b6', 'α=1.0  (N)')]:
    mean_y = np.mean(logR_vals)
    mean_x = np.mean(logN_vals)
    ref = exp * (logN_vals - mean_x) + mean_y
    ax4.plot(logN_vals[sort_idx], ref[sort_idx], '--', color=color,
             linewidth=1.5, alpha=0.6, label=f'Ref: {lbl}')

ax4.set_xlabel('log(N)', fontsize=13)
ax4.set_ylabel('log |r(N) − N·𝔖(N)/(log N)²|', fontsize=12)
ax4.set_title(f'Residual Decay Rate (log-log)\nFitted exponent α = {coeffs[0]:.3f}  |  Expected under GRH: α ≈ 0.5',
              fontsize=13, fontweight='bold')
ax4.legend(loc='upper left', fontsize=9, markerscale=2)
ax4.grid(True, alpha=0.3)
fig4.tight_layout()
fig4.savefig(f'{OUT_DIR}/comet_decay.png', dpi=150, bbox_inches='tight')
plt.close(fig4)
print("    → comet_decay.png saved")

# ══════════════════════════════════════════
# Summary
# ══════════════════════════════════════════
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"  N range: 4 to {MAX_N}, step 2")
print(f"  Even N count: {N_count}")
print(f"  r(N) range: [{r_N.min()}, {r_N.max()}], mean={r_N.mean():.1f}")
print(f"  HL main term range: [{HL_main.min():.1f}, {HL_main.max():.1f}]")
print(f"  Residual r(N)−HL: [{residual.min():.1f}, {residual.max():.1f}]")
print(f"  Residual std: {residual.std():.2f}")
print(f"  Fitted decay exponent α = {coeffs[0]:.4f}")
print(f"  Theoretical (GRH): α ≈ 0.5")
print(f"\n  Three-term decomposition:")
print(f"    R₁₁(N) = N·𝔖(N)/(log N)² — main arc (Λ₁Λ₁)")
print(f"    R₁₂ + R₂₁ = cross terms (Λ₁Λ₂, oscillatory from zeros)")
print(f"    R₂₂ = zero-zero interaction")
print(f"    R₁₃ + R₃₁ + R₃₃ = trivial zeros (negligible O(1/N))")
print(f"\nAll plots saved to: {OUT_DIR}")
