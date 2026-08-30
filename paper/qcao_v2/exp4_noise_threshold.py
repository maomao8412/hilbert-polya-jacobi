# -*- coding: utf-8 -*-
"""exp4: hardware-noise critical threshold eps_c.
Extends exp3 E6 to large noise (eps = 0.1, 0.2, 0.3, 0.5, 1.0, 2.0) to locate
the crossover from Wigner-Dyson repulsion (r~0.6) toward Poisson (r=0.3863).

Noise model (same as exp3 E6): multiplicative Gaussian on both diagonals,
  a_n -> a_n (1 + eps g_n),  b_n -> b_n |1 + eps h_n|,  g,h ~ N(0,1).
For eps >= ~0.3 the 'tracked levels' picture breaks down (spectrum moves far),
so we report TWO diagnostics:
  (i)  r_tracked  : mean gap ratio of the subset that still tracks the clean
                    zeros within 30% (the hardware 'can I still see my levels'
                    diagnostic; comparable to exp3 E6);
  (ii) r_full     : mean gap ratio over ALL unfolded eigenvalues (cubic
                    counting-function unfolding, consistent with exp1/exp3),
                    which remains valid for large noise.
40 independent noise realizations per (family, eps). RNG seed = 20260831.
Outputs exp4_results.json + exp4_noise_threshold.png.
"""
import numpy as np, json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

BASE = '/Coze/Drive/黎曼猜想论文审核/所有对话/主对话'
OUT = os.path.join(BASE, 'physics_experiments')
rng = np.random.default_rng(20260831)

def load_csv_matrix(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append([float(x) for x in line.split(',')])
    return np.array(rows)

def ab_from_H(H):
    a = np.diag(H).copy()
    b = np.diag(H, 1).copy()
    return a, b

GWIN = (5.0, 120.0)  # physical gamma window (zeros 14..~60 for J50; margin)

def gamma_all(H):
    """Return (gammas in physical window, frac_positive_eigenvalues).
    Large multiplicative noise destroys positive-definiteness (the moment /
    arithmetic structure): negative eigenvalues appear. They are excluded from
    the gamma mapping (1/sqrt(lambda) only defined on positive spectrum); the
    positive-eigenvalue fraction is itself reported as a structural diagnostic."""
    E = np.sort(np.linalg.eigvalsh(H))
    frac_pos = float(np.mean(E > 1e-18))
    Epos = E[E > 1e-18]
    g = 1.0/np.sqrt(Epos)
    g = np.sort(g)
    g = g[(g >= GWIN[0]) & (g <= GWIN[1])]
    return g, frac_pos

def unfold3(g):
    n = len(g)
    if n < 5:
        return None
    k = np.arange(n, dtype=float)
    try:
        c = np.polyfit(g, k, 3)
    except Exception:
        return None
    eu = np.polyval(c, g)
    eu = eu - eu[0]
    if eu[-1] == eu[0]:
        return None
    eu *= (n-1)/(eu[-1]-eu[0])
    return eu

def r_ratio_from_gaps(g):
    d = np.diff(g); d = d[d > 0]
    if len(d) < 3: return float('nan')
    r = np.minimum(d[:-1], d[1:])/np.maximum(d[:-1], d[1:])
    return float(np.mean(r))

def r_full(H):
    g, _ = gamma_all(H)
    eu = unfold3(g)
    if eu is None: return float('nan')
    return r_ratio_from_gaps(eu)

# ---- load matrices ----
H_zeta = load_csv_matrix(os.path.join(BASE, 'GRH交互网页/zeta_J50_matrix.csv'))
H_beta = load_csv_matrix(os.path.join(BASE, 'GRH交互网页/beta_J50_matrix.csv'))
djson = json.load(open(os.path.join(BASE, 'delta_matrix/results_20260831/delta_jacobi_raw.json')))
_j = djson['J22']
_ad = np.array([float(x) for x in _j['alphas']])
_bsq = np.array([float(x) for x in _j['betas_sq']])
_bd = np.sqrt(np.maximum(_bsq, 1e-300))
H_delta = np.diag(_ad)
for _i in range(len(_bd)):
    H_delta[_i, _i+1] = H_delta[_i+1, _i] = _bd[_i]

families = {'zeta': H_zeta, 'beta': H_beta, 'Delta': H_delta}
EPSS = [0.0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]
NREP = 40

results = {'meta': {'eps_list': EPSS, 'n_rep': NREP, 'seed': 20260831,
                    'noise_model': 'multiplicative Gaussian a*(1+eps g), b*|1+eps h|',
                    'theory': {'GUE': 0.5996, 'GOE': 0.5359, 'Poisson': 0.3863}}}
for name, H in families.items():
    a0, b0 = ab_from_H(H)
    g0, fp0 = gamma_all(H)
    r0 = r_full(H)
    row = {'0.0': {'r_full': r0, 'r_tracked': r0, 'n_tracked': int(len(g0)),
                   'frac_positive': fp0}}
    print(f'=== {name} (N={len(a0)}) clean r_full={r0:.3f}, '
          f'{len(g0)} levels in gamma window, frac_pos={fp0:.3f} ===', flush=True)
    for eps in EPSS[1:]:
        rf, rt, nt, fp = [], [], [], []
        for _ in range(NREP):
            a2 = a0*(1+eps*rng.standard_normal(len(a0)))
            b2 = b0*np.abs(1+eps*rng.standard_normal(len(b0)))
            Hn = np.diag(a2)
            for i in range(len(b2)):
                Hn[i, i+1] = Hn[i+1, i] = b2[i]
            rf.append(r_full(Hn))
            # tracked subset vs clean zeros
            gn, fpos = gamma_all(Hn)
            fp.append(fpos)
            tracked = []
            for gz in g0:
                if len(gn) == 0: break
                j = np.argmin(np.abs(gn-gz))
                if abs(gn[j]-gz)/gz < 0.30:
                    tracked.append(gn[j])
            tracked = np.sort(np.array(tracked))
            nt.append(len(tracked))
            if len(tracked) >= 4:
                eu = unfold3(tracked)
                if eu is not None:
                    rt.append(r_ratio_from_gaps(eu))
        row[str(eps)] = {'r_full_mean': float(np.nanmean(rf)),
                         'r_full_std': float(np.nanstd(rf)),
                         'r_tracked_mean': float(np.nanmean(rt)) if rt else float('nan'),
                         'r_tracked_std': float(np.nanstd(rt)) if rt else float('nan'),
                         'n_tracked_mean': float(np.mean(nt)),
                         'frac_positive_mean': float(np.mean(fp))}
        print(f"  eps={eps:g}: r_full={np.nanmean(rf):.3f}+-{np.nanstd(rf):.3f}  "
              f"r_tracked={np.nanmean(rt) if rt else float('nan'):.3f}  "
              f"ntrack={np.mean(nt):.1f}/{len(g0)}  frac_pos={np.mean(fp):.3f}", flush=True)
    results[name] = row

json.dump(results, open(os.path.join(OUT, 'exp4_results.json'), 'w'), indent=1)

# ---- plot ----
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
colors = {'zeta': 'tab:blue', 'beta': 'tab:orange', 'Delta': 'tab:green'}
for name in families:
    xs = [float(e) for e in results[name].keys()]
    yf = [results[name][str(x)].get('r_full', results[name][str(x)].get('r_full_mean')) for x in xs]
    yfs = [results[name][str(x)].get('r_full_std', 0) for x in xs]
    ax[0].errorbar(xs, yf, yerr=yfs, marker='o', capsize=3, color=colors[name], label=name)
ax[0].axhline(0.3863, color='k', ls='--', lw=1, label='Poisson 0.386')
ax[0].axhspan(0.5359-0.02, 0.5996+0.02, color='gray', alpha=0.15, label='GOE–GUE band')
ax[0].set_xscale('log'); ax[0].set_xlabel(r'relative coefficient noise $\epsilon$')
ax[0].set_ylabel(r'mean gap ratio $\langle r\rangle$ (full spectrum)')
ax[0].set_title('(a) Noise crossover: arithmetic GUE/GOE $\\to$ Poisson\n(40 realizations/point)')
ax[0].legend(fontsize=8)
for name in families:
    xs = [float(e) for e in results[name].keys() if float(e) > 0]
    yt = [results[name][str(x)]['r_tracked_mean'] for x in xs]
    ntr = [results[name][str(x)]['n_tracked_mean'] for x in xs]
    ax[1].plot(xs, ntr, marker='s', color=colors[name], label=name)
ax[1].set_xscale('log'); ax[1].set_xlabel(r'relative coefficient noise $\epsilon$')
ax[1].set_ylabel('mean number of trackable levels')
ax[1].set_title('(b) Hardware diagnostic: how many zero-levels remain\nidentifiable (within 30%) under noise')
ax[1].legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'exp4_noise_threshold.png'), dpi=150)
print('saved exp4_results.json + exp4_noise_threshold.png', flush=True)
