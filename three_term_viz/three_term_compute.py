#!/usr/bin/env python3
"""
Three-Channel Decomposition of the Riemann Zeta Function
Precomputation script for interactive visualization.

Formula:
  ζ(s) = Li_s(z₁) + Li_s(z₂) − Γ(1−s)·L^{s−1} − δ(s)

where:
  z₁ = 2 − √2 ≈ 0.5858   (channel 1: polylog at algebraic point z₁)
  z₂ = √2 − 1 ≈ 0.4142   (channel 2: polylog at algebraic point z₂)
  L = ln(√2+1) = arsinh(1) ≈ 0.8814
  δ(s) = regular remainder / "residual potential"

Outputs three_term_data.json with all precomputed data for embedding in HTML.
"""

import json
import mpmath as mp

mp.mp.dps = 50  # high precision

# ── Constants ──────────────────────────────────────────────────────
z1 = 2 - mp.sqrt(2)
z2 = mp.sqrt(2) - 1
L = mp.log(mp.sqrt(2) + 1)  # = arsinh(1)


def compute_channels(s):
    """Compute all four channels at point s = σ + it.
    
    Returns dict with:
        zeta: ζ(s)
        li1:  Li_s(z₁)   [channel 1]
        li2:  Li_s(z₂)   [channel 2]
        gamma_c: Γ(1-s)·L^{s-1}  [gamma correction]
        delta: δ(s) = Li_s(z₁)+Li_s(z₂) − Γ(1-s)L^{s-1} − ζ(s)  [remainder]
    """
    zeta_val = mp.zeta(s)
    li1 = mp.polylog(s, z1)
    li2 = mp.polylog(s, z2)
    gamma_c = mp.gamma(1 - s) * L**(s - 1)
    delta = li1 + li2 - gamma_c - zeta_val
    
    return {
        'zeta': zeta_val,
        'li1': li1,
        'li2': li2,
        'gamma': gamma_c,
        'delta': delta,
    }


def complex_to_pair(c):
    """Convert mpmath complex to [real, imag] list."""
    return [float(mp.re(c)), float(mp.im(c))]


def mod_list(arr):
    """Convert list of complex pairs to list of moduli."""
    return [float((r**2 + i**2)**0.5) for r, i in arr]


# ── 1. Critical line data (Section 2) ─────────────────────────────
print("Computing critical line data...")
t_start = mp.mpf(0)
t_end = mp.mpf(50)
t_step = mp.mpf('0.05')
n_points = int((t_end - t_start) / t_step) + 1

critical_line = {
    't': [],
    'zeta_real': [],
    'zeta_imag': [],
    'zeta_mod': [],
    'li1_real': [],
    'li1_imag': [],
    'li1_mod': [],
    'li2_real': [],
    'li2_imag': [],
    'li2_mod': [],
    'gamma_real': [],
    'gamma_imag': [],
    'gamma_mod': [],
    'delta_real': [],
    'delta_imag': [],
    'delta_mod': [],
    'sum_real': [],
    'sum_imag': [],
    'sum_mod': [],
}

for i in range(n_points):
    t = t_start + i * t_step
    s = mp.mpc('0.5', t)
    ch = compute_channels(s)
    
    critical_line['t'].append(float(t))
    critical_line['zeta_real'].append(float(mp.re(ch['zeta'])))
    critical_line['zeta_imag'].append(float(mp.im(ch['zeta'])))
    critical_line['zeta_mod'].append(float(abs(ch['zeta'])))
    
    critical_line['li1_real'].append(float(mp.re(ch['li1'])))
    critical_line['li1_imag'].append(float(mp.im(ch['li1'])))
    critical_line['li1_mod'].append(float(abs(ch['li1'])))
    
    critical_line['li2_real'].append(float(mp.re(ch['li2'])))
    critical_line['li2_imag'].append(float(mp.im(ch['li2'])))
    critical_line['li2_mod'].append(float(abs(ch['li2'])))
    
    critical_line['gamma_real'].append(float(mp.re(ch['gamma'])))
    critical_line['gamma_imag'].append(float(mp.im(ch['gamma'])))
    critical_line['gamma_mod'].append(float(abs(ch['gamma'])))
    
    critical_line['delta_real'].append(float(mp.re(ch['delta'])))
    critical_line['delta_imag'].append(float(mp.im(ch['delta'])))
    critical_line['delta_mod'].append(float(abs(ch['delta'])))
    
    # li1 + li2 - gamma - delta should equal zeta
    s_total = ch['li1'] + ch['li2'] - ch['gamma'] - ch['delta']
    critical_line['sum_real'].append(float(mp.re(s_total)))
    critical_line['sum_imag'].append(float(mp.im(s_total)))
    critical_line['sum_mod'].append(float(abs(s_total)))
    
    if i % 100 == 0:
        print(f"  t={float(t):.1f} (point {i}/{n_points})")

print("  Done with critical line.")


# ── 2. 3D surface data (Section 3) ────────────────────────────────
print("Computing 3D surface data...")
sigma_start = 0.0
sigma_end = 0.99  # avoid σ=1 (zeta pole)
n_sigma = 30  # σ resolution
n_t_surf = 60  # t resolution
t_surf_end = 30.0

sigma_vals = [sigma_start + i * (sigma_end - sigma_start) / (n_sigma - 1)
              for i in range(n_sigma)]
t_surf_vals = [i * t_surf_end / (n_t_surf - 1) for i in range(n_t_surf)]

surface = {
    'sigma': sigma_vals,
    't': t_surf_vals,
    'zeta_mod': [],       # n_sigma × n_t_surf
    'li1_mod': [],
    'li2_mod': [],
    'gamma_mod': [],
    'delta_mod': [],
}

for si, sigma in enumerate(sigma_vals):
    zeta_row = []
    li1_row = []
    li2_row = []
    gamma_row = []
    delta_row = []
    for ti, t_val in enumerate(t_surf_vals):
        s = mp.mpc(sigma, t_val)
        ch = compute_channels(s)
        zeta_row.append(float(abs(ch['zeta'])))
        li1_row.append(float(abs(ch['li1'])))
        li2_row.append(float(abs(ch['li2'])))
        gamma_row.append(float(abs(ch['gamma'])))
        delta_row.append(float(abs(ch['delta'])))
    surface['zeta_mod'].append(zeta_row)
    surface['li1_mod'].append(li1_row)
    surface['li2_mod'].append(li2_row)
    surface['gamma_mod'].append(gamma_row)
    surface['delta_mod'].append(delta_row)
    print(f"  σ={sigma:.3f} (row {si+1}/{n_sigma})")

print("  Done with surface.")


# ── 3. Zero fingerprint data (Section 4) ──────────────────────────
print("Computing zero fingerprint data (first 20 zeros)...")
n_zeros = 20
zeros_data = {
    'zero_index': list(range(1, n_zeros + 1)),
    'gamma_n': [],            # imaginary part of nth zero
    'li1_mod': [],
    'li2_mod': [],
    'gamma_mod': [],
    'delta_mod': [],
    'zeta_mod': [],
    'li1_phase': [],
    'li2_phase': [],
    'gamma_phase': [],
    'delta_phase': [],
}

for n in range(1, n_zeros + 1):
    z = mp.zetazero(n)
    gamma_n = float(mp.im(z))
    s = z  # = 0.5 + i*gamma_n
    ch = compute_channels(s)
    
    zeros_data['gamma_n'].append(gamma_n)
    zeros_data['zeta_mod'].append(float(abs(ch['zeta'])))
    zeros_data['li1_mod'].append(float(abs(ch['li1'])))
    zeros_data['li2_mod'].append(float(abs(ch['li2'])))
    zeros_data['gamma_mod'].append(float(abs(ch['gamma'])))
    zeros_data['delta_mod'].append(float(abs(ch['delta'])))
    zeros_data['li1_phase'].append(float(mp.arg(ch['li1'])))
    zeros_data['li2_phase'].append(float(mp.arg(ch['li2'])))
    zeros_data['gamma_phase'].append(float(mp.arg(ch['gamma'])))
    zeros_data['delta_phase'].append(float(mp.arg(ch['delta'])))
    
    print(f"  zero {n}: γ_n = {gamma_n:.4f}, |ζ| = {float(abs(ch['zeta'])):.2e}")

print("  Done with zeros.")


# ── Constants for display ────────────────────────────────────────
constants = {
    'z1': float(z1),
    'z2': float(z2),
    'L': float(L),
}


# ── Assemble and write JSON ──────────────────────────────────────
output = {
    'constants': constants,
    'critical_line': critical_line,
    'surface': surface,
    'zeros': zeros_data,
}

output_path = './codeact/output/three_term_data.json'
with open(output_path, 'w') as f:
    json.dump(output, f)

print(f"\nData written to {output_path}")
print(f"Critical line points: {n_points}")
print(f"Surface grid: {n_sigma} × {n_t_surf}")
print(f"Zeros: {n_zeros}")

# Also compute file size for reference
import os
size_kb = os.path.getsize(output_path) / 1024
print(f"File size: {size_kb:.1f} KB")
