#!/usr/bin/env python3
"""
PRIMES DIRECTLY FROM THE MATRIX — no diagonalization, no zeros extracted.

Chain (forward only):
  seed q=sqrt2-1 -> three-term polylog identity -> xi Taylor coefficients
  -> log-moments S_k -> Hankel -> Gram-Schmidt -> Jacobi matrix J_N
Matrix function calculus (no eigenvalues ever computed):
  Hadamard (unconditional):  xi(s)/xi(1/2) = prod_gamma (1 + u^2/gamma^2), u=s-1/2
  spectrum lambda_n = 1/gamma_n^2  =>
  xi'/xi(s) = 2u Tr[ J (I + u^2 J)^{-1} ]
  zeta'/zeta(s) = xi'/xi(s) - 1/s - 1/(s-1) + (1/2)log pi - (1/2) digamma(s/2)
  Perron: psi(x) = -(x^c/pi) Re int_0^T [zeta'/zeta(c+it)/(c+it)] e^{it log x} dt
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpmath import mp, mpf, mpc, digamma as mpdigamma, zeta as mpzeta, log as mplog, pi as mppi, gamma as mpgamma, taylor as mptaylor

mp.dps = 60

def load_J(path):
    cp = json.load(open(path))
    alphas = np.array([float(mpf(s)) for s in cp['alphas']])
    betas = np.array([float(mpf(s)) for s in cp['betas_sq']]) ** 0.5
    n = len(alphas)
    J = np.zeros((n, n))
    J[np.arange(n), np.arange(n)] = alphas
    J[np.arange(n-1), np.arange(1, n)] = betas
    J[np.arange(1, n), np.arange(n-1)] = betas
    return J

J100 = load_J('jacobi100_checkpoint.json')
J50 = load_J('jacobi_checkpoint.json')
print(f"J_100 Tr = {np.trace(J100):.8f};  J_50 Tr = {np.trace(J50):.8f}")

LOG_PI_HALF = 0.5*float(mplog(mppi))

def make_funcs(J):
    n = J.shape[0]; I = np.eye(n)
    def xi_ld(u):
        # u complex; xi'/xi = 2u Tr[J(I+u^2 J)^{-1}]
        X = np.linalg.solve(I + u*u*J, J)
        return 2.0*u*np.trace(X)
    def zeta_ld(s):
        u = s - 0.5
        return xi_ld(u) - 1.0/s - 1.0/(s-1.0) + LOG_PI_HALF \
               - 0.5*float(mp.re(mpdigamma(mpc(s.real, s.imag)/2.0)))
    return xi_ld, zeta_ld

xi100, zl100 = make_funcs(J100)
xi50, zl50 = make_funcs(J50)

# ---- exact reference via mpmath (analytic derivative, no finite differences) ----
def xi_mp(s):
    return 0.5*s*(s-1)*mppi**(-s/2)*mpgamma(s/2)*mpzeta(s)
def zeta_ld_exact(s):
    return mpzeta(s, 1)/mpzeta(s)
def xi_ld_exact(s):
    return zeta_ld_exact(s) + 1/s + 1/(s-1) - 0.5*mplog(mppi) + 0.5*mpdigamma(s/2)

# ---------- validation panels ----------
print("\n=== zeta'/zeta(s) real axis: matrix vs exact ===")
print(f"{'s':>5} {'matrix':>16} {'exact':>16} {'abs err':>10}")
for s in [1.5, 2.0, 3.0, 4.0, 7.0]:
    vm = zl100(complex(s, 0)).real
    ve = float(zeta_ld_exact(mpf(s)))
    print(f"{s:5.1f} {vm:16.10f} {ve:16.10f} {abs(vm-ve):10.2e}")

print("\n=== Im xi'/xi(1/2+it): matrix vs exact ===")
print(f"{'t':>6} {'matrix':>14} {'exact':>14} {'rel err':>10}")
for t in [5.0, 10.0, 14.0, 20.0, 30.0, 50.0, 80.0, 120.0]:
    vm = xi100(1j*t).imag
    ve = float(mp.im(xi_ld_exact(mpc(0.5, t))))
    print(f"{t:6.1f} {vm:14.8f} {ve:14.8f} {abs(vm-ve)/abs(ve):10.2e}")

# ---------- true psi ----------
def sieve_primes(n):
    is_p = np.ones(n+1, dtype=bool); is_p[:2] = False
    for p in range(2, int(n**0.5)+1):
        if is_p[p]: is_p[p*p::p] = False
    return np.nonzero(is_p)[0]

xmax = 120.0
primes = sieve_primes(int(xmax)+2)
pos, wt = [], []
for p in primes:
    pk = int(p)
    while pk <= xmax:
        pos.append(pk); wt.append(np.log(p)); pk *= int(p)
idx = np.argsort(pos); pos = np.array(pos)[idx]; cum = np.cumsum(np.array(wt)[idx])

# ---------- Perron from matrix ----------
def perron_psi(zeta_ld, T, dt, c, xs):
    ts = np.arange(dt, T, dt)
    win = np.cos(np.pi*ts/(2*T))**2
    g = np.empty(len(ts), dtype=complex)
    for k, t in enumerate(ts):
        g[k] = zeta_ld(complex(c, t))/complex(c, t)*win[k]
    lx = np.log(xs)
    integ = np.exp(1j*np.outer(lx, ts)) @ g * dt
    return -(xs**c/np.pi)*integ.real

xs = np.linspace(6, xmax, 3000)
psi_t = cum[np.searchsorted(pos, xs, side='right')-1]
c = 1.5
psi_m100 = perron_psi(zl100, 140.0, 0.05, c, xs)
psi_m50 = perron_psi(zl50, 85.0, 0.05, c, xs)
m = xs > 40
rms100 = np.sqrt(np.mean((psi_m100-psi_t)[m]**2))
rms50 = np.sqrt(np.mean((psi_m50-psi_t)[m]**2))
print(f"\nPerron RMS (x>40): J_100/T=140 -> {rms100:.3f};   J_50/T=85 -> {rms50:.3f}")

# ---------- FIGURE 1: money shot ----------
fig, ax = plt.subplots(2, 1, figsize=(11.5, 8.5), gridspec_kw={'height_ratios':[2.3,1]})
ax[0].plot(xs, psi_t, color='black', lw=1.4, label=r'true $\psi(x)=\sum_{p^k\leq x}\log p$')
ax[0].plot(xs, psi_m100, color='#c0392b', lw=1.1, alpha=0.9,
           label=r'from $J_{100}$: invert $I+u^2J$, trace, Perron integral (no eigenvalues used)')
ax[0].plot(xs, xs-np.log(2*np.pi), color='#2c3e50', lw=0.7, ls='--', alpha=0.6, label=r'$x-\log 2\pi$')
ax[0].set_ylabel(r'$\psi(x)$', fontsize=12)
ax[0].set_title('Primes decoded straight from the matrix — no diagonalization, no zeros, no prime table',
                fontsize=12.5, fontweight='bold')
ax[0].legend(loc='upper left', fontsize=10); ax[0].grid(alpha=0.25)
ax[1].plot(xs, psi_m100-psi_t, color='#c0392b', lw=0.7, label='error: $J_{100}$ reconstruction')
ax[1].axhline(0, color='black', lw=0.5)
ax[1].set_ylabel('error', fontsize=11); ax[1].set_xlabel('x', fontsize=12)
ax[1].legend(fontsize=9.5); ax[1].grid(alpha=0.25)
ax[1].text(0.98, 0.9, f'RMS error (x>40): {rms100:.2f}\nfinite-matrix cutoff T=140',
           transform=ax[1].transAxes, ha='right', va='top', fontsize=10.5,
           bbox=dict(boxstyle='round', fc='#fdf2e9', ec='#c0392b', alpha=0.9))
plt.tight_layout(); plt.savefig('primes_direct_from_matrix.png', dpi=150); plt.close()
print('saved primes_direct_from_matrix.png')

# ---------- FIGURE 2: resolvent IS zeta ----------
tt = np.linspace(0.5, 140, 2000)
im_mat = np.array([xi100(1j*t).imag for t in tt])
tt_ref = np.linspace(0.5, 140, 280)
im_ref = np.array([float(mp.im(xi_ld_exact(mpc(0.5, t)))) for t in tt_ref])
fig, ax = plt.subplots(2, 1, figsize=(11.5, 8), gridspec_kw={'height_ratios':[2.3,1]}, sharex=True)
ax[0].plot(tt, im_mat, color='#c0392b', lw=1.5, label=r'matrix: $\mathrm{Im}\,2it\,\mathrm{Tr}[J(I-t^2J)^{-1}]$')
ax[0].plot(tt_ref, im_ref, color='black', lw=0.9, ls='--', label=r'exact $\mathrm{Im}\,\xi^{\prime}/\xi(\frac{1}{2}+it)$')
for g_ in [14.135, 21.022, 25.011, 30.425, 32.935, 37.586, 40.919, 43.327, 48.005, 49.774,
           52.970, 56.446, 59.347, 60.832, 65.113, 67.080, 69.546, 72.067, 75.705, 77.145]:
    ax[0].axvline(g_, color='#27ae60', lw=0.5, alpha=0.6)
ax[0].set_ylabel(r"Im $\xi'/\xi$", fontsize=12)
ax[0].set_title(r'The matrix resolvent is the zeta function: $\xi^{\prime}/\xi(\frac{1}{2}+it)=2it\,\mathrm{Tr}[J(I-t^2J)^{-1}]$'
                '\n(green: true zeros — poles of the resolvent; never supplied to the matrix)', fontsize=11.5)
ax[0].legend(loc='upper right', fontsize=10); ax[0].grid(alpha=0.25)
ax[1].plot(tt_ref, np.interp(tt_ref, tt, im_mat)-im_ref, color='#c0392b', lw=0.8)
ax[1].axhline(0, color='black', lw=0.5)
ax[1].set_ylabel('error', fontsize=11); ax[1].set_xlabel('t', fontsize=12)
ax[1].set_xlim(0, 140); ax[1].grid(alpha=0.25)
plt.tight_layout(); plt.savefig('matrix_resolvent_is_zeta.png', dpi=150); plt.close()
print('saved matrix_resolvent_is_zeta.png')

# ---------- FIGURE 3: matrix powers = moments ----------
ks = np.arange(1, 9)
trJk = [np.trace(np.linalg.matrix_power(J100, k)) for k in ks]
# exact S_k via Taylor coefficients of log xi(1/2+u)
f = lambda u: mplog(xi_mp(mpf('0.5')+u)/xi_mp(mpf('0.5')))
coefs = mptaylor(f, 0, 16)
S_exact = [float(k*abs(coefs[2*k])) for k in ks]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].bar(ks-0.18, trJk, width=0.36, color='#c0392b', label=r'$\mathrm{Tr}(J^k)$ from 199 matrix entries')
ax[0].bar(ks+0.18, S_exact, width=0.36, color='#2c3e50', alpha=0.55, label=r'$S_k=\sum_\gamma\gamma^{-2k}$ from $\xi$')
ax[0].set_yscale('log'); ax[0].set_xlabel('k'); ax[0].set_ylabel('moment')
ax[0].set_title(r'Matrix powers are zeta moments: $\mathrm{Tr}(J^k)=S_k$', fontsize=12)
ax[0].set_xticks(ks); ax[0].legend(fontsize=9.5); ax[0].grid(alpha=0.25, axis='y')
errs = [abs(trJk[k-1]-S_exact[k-1])/S_exact[k-1] for k in ks]
ax[1].semilogy(ks, errs, 'o-', color='#c0392b')
ax[1].set_xlabel('k'); ax[1].set_ylabel('relative error')
ax[1].set_title('Moment agreement (gap at k=1 = prime info beyond order 100)', fontsize=11)
ax[1].set_xticks(ks); ax[1].grid(alpha=0.25)
for k, e in zip(ks, errs):
    if e > 1e-12: ax[1].annotate(f'{e:.0%}', (k, e), textcoords='offset points', xytext=(0, 8), ha='center', fontsize=9)
plt.tight_layout(); plt.savefig('matrix_powers_are_moments.png', dpi=150); plt.close()
print('saved matrix_powers_are_moments.png')
print('moment rel errs:', ['%.2e' % e for e in errs])

# ---------- FIGURE 4: order 50 vs 100 ----------
fig, ax = plt.subplots(figsize=(11.5, 6))
ax.plot(xs, psi_t, color='black', lw=1.6, label=r'true $\psi(x)$')
ax.plot(xs, psi_m50, color='#2980b9', lw=1.0, alpha=0.85, label=r'$J_{50}$ (25 zeros locked), Perron $T=85$')
ax.plot(xs, psi_m100, color='#c0392b', lw=1.0, alpha=0.9, label=r'$J_{100}$ (50 zeros locked), Perron $T=140$')
ax.set_xlabel('x', fontsize=12); ax.set_ylabel(r'$\psi(x)$', fontsize=12)
ax.set_title('Bigger matrix = more primes decoded: the 199 entries of $J_{100}$ resolve finer prime structure',
             fontsize=12)
ax.legend(loc='upper left', fontsize=10.5); ax[0].grid(alpha=0.25) if False else ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig('matrix_order_primes.png', dpi=150); plt.close()
print('saved matrix_order_primes.png')
