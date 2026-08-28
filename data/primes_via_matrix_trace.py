#!/usr/bin/env python3
"""
PRIMES FROM PURE MATRIX ARITHMETIC — no zeros, no integrals, no diagonalization.

Encoding (unconditional):
  xi(s)/xi(1/2) = det(I + u^2 J),   u = s-1/2                      (Hadamard)
  xi'/xi(s)     = 2u Tr[ J (I + u^2 J)^{-1} ]                      (resolvent)
Zero-sum of the explicit formula, with A = J^{-1/2} (spectrum gamma_n):
  sum_rho x^rho/rho = 2 sqrt(x) Re Tr[ exp(i A log x) (1/2 I + i A)^{-1} ]
So psi(x) is a matrix trace:
  psi(x) = x - 2 sqrt(x) Re Tr[ e^{iA log x} (1/2 I + iA)^{-1} ] - log(2pi) - 1/2 log(1-x^-2)

All operations are finite matrix arithmetic:
  J = L L^T          Cholesky (J tridiagonal -> L bidiagonal)         [products]
  L^{-1}             triangular solve                                  [solve]
  polar iteration    Y <- 1/2 (Y + Y^{-T})  -> orthogonal U           [products + solves]
  A = U^T L^{-1}    = J^{-1/2}
  e^{iA t}           scaling-and-squaring with Taylor series          [products only]
  (1/2 I + iA)^{-1}  complex linear solve, done once                   [solve]
No eigenvalue is ever computed. No contour integral. No prime table in.
Primes "pop out": differences psi(n+1/2)-psi(n-1/2) equal von Mangoldt spikes.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpmath import mp, mpf, mpc, zeta as mpzeta, log as mplog, pi as mppi, gamma as mpgamma, digamma as mpdigamma

mp.dps = 50

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

def inv_sqrt_via_polar(J):
    """J^{-1/2} using Cholesky + polar iteration. Products and triangular solves only."""
    n = J.shape[0]
    L = np.linalg.cholesky(J)                 # J = L L^T, L lower bidiagonal
    Linv = np.linalg.solve(L, np.eye(n))      # L^{-1} by forward substitution
    Y = Linv.copy()
    for _ in range(40):
        YinvT = np.linalg.solve(Y.T, np.eye(n))      # (Y^T)^{-1} = (Y^{-1})^T
        Ynew = 0.5*(Y + YinvT)
        if np.linalg.norm(Ynew - Y) < 1e-13*np.linalg.norm(Y):
            Y = Ynew; break
        Y = Ynew
    U = Y                                     # orthogonal factor of L^{-1} = U H
    A = U.T @ Linv                            # H = J^{-1/2}
    return A

def expm_sq(M, order=12):
    """Matrix exponential by scaling-and-squaring + Taylor. Products only."""
    nrm = np.linalg.norm(M, ord=np.inf)
    q = max(0, int(np.ceil(np.log2(nrm + 1e-300))) + 1)
    E = M/(2.0**q)
    T = np.eye(M.shape[0], dtype=complex)
    term = np.eye(M.shape[0], dtype=complex)
    for k in range(1, order+1):
        term = term @ E / k
        T = T + term
    for _ in range(q):
        T = T @ T
    return T

# ---------------- build matrices ----------------
J100 = load_J('jacobi100_checkpoint.json')
J50  = load_J('jacobi_checkpoint.json')
A100 = inv_sqrt_via_polar(J100)
A50  = inv_sqrt_via_polar(J50)
n100, n50 = 100, 50
B100 = np.linalg.solve(0.5*np.eye(n100) + 1j*A100, np.eye(n100, dtype=complex))  # (1/2 I + iA)^{-1}
B50  = np.linalg.solve(0.5*np.eye(n50)  + 1j*A50,  np.eye(n50,  dtype=complex))
# sanity: A^2 J = I ?
chk100 = np.linalg.norm(A100 @ A100 @ J100 - np.eye(n100))/n100
chk50  = np.linalg.norm(A50 @ A50 @ J50 - np.eye(n50))/n50
print(f"A=J^-1/2 check  ||A^2 J - I||/n : J100 {chk100:.2e}, J50 {chk50:.2e}")

def psi_trace(x, A, B):
    t = np.log(x)
    E = expm_sq(1j*A*t)
    zsum = 2.0*np.sqrt(x)*np.trace(E @ B).real
    return x - zsum - np.log(2*np.pi) - 0.5*np.log(1.0 - x**-2)

# ---------------- exact references ----------------
def von_mangoldt(n):
    m = n
    for p in range(2, int(n**0.5)+1):
        if m % p == 0:
            while m % p == 0: m //= p
            return np.log(p) if m == 1 else 0.0
    return np.log(n) if n > 1 else 0.0
LAM = np.array([von_mangoldt(n) for n in range(0, 3000)])

def xi_ld_matrix(J, u):
    n = J.shape[0]
    X = np.linalg.solve(np.eye(n) + u*u*J, J)
    return 2.0*u*np.trace(X)
LOG_PI_HALF = 0.5*float(mplog(mppi))
def zeta_ld_matrix(J, s):
    u = s - 0.5
    return xi_ld_matrix(J, u) - 1.0/s - 1.0/(s-1.0) + LOG_PI_HALF \
           - 0.5*float(mp.re(mpdigamma(mpc(s.real, s.imag)/2.0)))
def zeta_ld_exact(s):
    return mpzeta(s, derivative=1)/mpzeta(s)
def xi_ld_exact(s):
    return zeta_ld_exact(s) + 1/s + 1/(s-1) - LOG_PI_HALF + 0.5*mpdigamma(s/2)

print("\n=== zeta'/zeta(s): matrix resolvent vs exact (mpmath) ===")
print(f"{'s':>5} {'matrix':>15} {'exact':>15} {'abs err':>10}")
for s in [1.5, 2.0, 3.0, 4.0, 7.0]:
    vm = zeta_ld_matrix(J100, complex(s, 0)).real
    ve = float(mp.re(zeta_ld_exact(mpf(s))))
    print(f"{s:5.1f} {vm:15.8f} {ve:15.8f} {abs(vm-ve):10.2e}")

print("\n=== Im xi'/xi(1/2+it): matrix vs exact ===")
print(f"{'t':>6} {'matrix':>14} {'exact':>14} {'rel err':>10}")
for t in [5.0, 10.0, 14.0, 20.0, 30.0, 50.0, 80.0]:
    vm = xi_ld_matrix(J100, 1j*t).imag
    ve = float(mp.im(xi_ld_exact(mpc(0.5, t))))
    print(f"{t:6.1f} {vm:14.8f} {ve:14.8f} {abs(vm-ve)/abs(ve):10.2e}")

# ---------------- psi reconstruction ----------------
xmax = 120.0
cum = np.cumsum(LAM[:int(xmax)+2])
xs = np.linspace(6, xmax, 300)
psi_t = np.interp(xs, np.arange(int(xmax)+2), cum)
psi_m100 = np.array([psi_trace(x, A100, B100) for x in xs])
psi_m50  = np.array([psi_trace(x, A50,  B50)  for x in xs])
m = xs > 40
rms100 = np.sqrt(np.mean((psi_m100-psi_t)[m]**2))
rms50  = np.sqrt(np.mean((psi_m50-psi_t)[m]**2))
print(f"\ntrace-formula psi RMS (x>40): J100 {rms100:.3f}, J50 {rms50:.3f}")

# ---------------- primes pop out at integers ----------------
Nmax = 52
spikes = np.zeros(Nmax+1)
for n in range(2, Nmax+1):
    spikes[n] = psi_trace(float(n)+0.5, A100, B100) - psi_trace(float(n)-0.5, A100, B100)
true_pw = [n for n in range(2, Nmax+1) if LAM[n] > 0]
recovered = [n for n in range(2, Nmax+1) if spikes[n] > 0.8]
hits = len(set(recovered) & set(true_pw))
fp = sorted(set(recovered) - set(true_pw))
missed = sorted(set(true_pw) - set(recovered))
print(f"\nprime powers recovered {hits}/{len(true_pw)}; false pos {fp}; missed {missed}")
print(f"spike decode RMS n=2..{Nmax}: {np.sqrt(np.mean((spikes[2:]-LAM[2:Nmax+1])**2)):.3f}")

# ---------------- FIGURE 1: psi from matrix arithmetic ----------------
fig, ax = plt.subplots(2, 1, figsize=(11.5, 8.5), gridspec_kw={'height_ratios':[2.3,1]})
ax[0].plot(xs, psi_t, color='black', lw=1.6, label=r'true $\psi(x)$')
ax[0].plot(xs, psi_m100, color='#c0392b', lw=1.1, alpha=0.9,
           label=r'$J_{100}$ trace: $x-2\sqrt{x}\,\mathrm{Re}\,\mathrm{Tr}[e^{iA\log x}(\frac{1}{2}I+iA)^{-1}]-\log 2\pi-\cdots$')
ax[0].plot(xs, xs-np.log(2*np.pi), color='#2c3e50', lw=0.7, ls='--', alpha=0.6)
ax[0].set_ylabel(r'$\psi(x)$', fontsize=12)
ax[0].set_title('Primes from matrix arithmetic alone — Cholesky, polar iteration, squaring. No zeros, no integrals',
                fontsize=11.5, fontweight='bold')
ax[0].legend(loc='upper left', fontsize=9.5); ax[0].grid(alpha=0.25)
ax[1].plot(xs, psi_m100-psi_t, color='#c0392b', lw=0.7, label='error')
ax[1].axhline(0, color='black', lw=0.5)
ax[1].set_ylabel('error', fontsize=11); ax[1].set_xlabel('x', fontsize=12)
ax[1].legend(fontsize=9.5); ax[1].grid(alpha=0.25)
ax[1].text(0.98, 0.9, f'RMS error (x>40): {rms100:.2f}', transform=ax[1].transAxes,
           ha='right', va='top', fontsize=10.5,
           bbox=dict(boxstyle='round', fc='#fdf2e9', ec='#c0392b', alpha=0.9))
plt.tight_layout(); plt.savefig('primes_matrix_trace_psi.png', dpi=150); plt.close()
print('saved primes_matrix_trace_psi.png')

# ---------------- FIGURE 2: spikes pop out ----------------
ns = np.arange(2, Nmax+1)
fig, ax = plt.subplots(2, 1, figsize=(12, 7.5), gridspec_kw={'height_ratios':[1,1]}, sharex=True)
ax[0].bar(ns-0.2, LAM[ns], width=0.4, color='black', alpha=0.75, label=r'true $\Lambda(n)$')
ax[0].bar(ns+0.2, np.maximum(spikes[ns], 0), width=0.4, color='#c0392b', alpha=0.85,
          label=r'matrix trace $\psi(n+\frac{1}{2})-\psi(n-\frac{1}{2})$')
ax[0].set_ylabel(r'$\Lambda(n)$', fontsize=12)
ax[0].set_title('Primes pop out at integers: von Mangoldt spikes from pure matrix arithmetic',
                fontsize=12, fontweight='bold')
ax[0].legend(fontsize=10); ax[0].grid(alpha=0.25, axis='y')
for n in [p for p in ns if LAM[p] > 0 and abs(LAM[p]-np.log(p)) < 1e-9]:
    ax[0].annotate(str(n), (n, LAM[n]), textcoords='offset points', xytext=(0, 3),
                   ha='center', fontsize=7.5, color='#1a4d8f')
err_sp = spikes[ns] - LAM[ns]
ax[1].bar(ns, err_sp, width=0.7, color=np.where(np.abs(err_sp) < 0.5, '#27ae60', '#c0392b'), alpha=0.85)
ax[1].axhline(0, color='black', lw=0.5)
ax[1].set_ylabel('decoding error', fontsize=11); ax[1].set_xlabel('n', fontsize=12)
ax[1].set_title(f'recovered prime powers {hits}/{len(true_pw)} up to {Nmax} (green: |err|<0.5)', fontsize=10.5)
ax[1].grid(alpha=0.25, axis='y')
plt.tight_layout(); plt.savefig('primes_pop_out_spikes.png', dpi=150); plt.close()
print('saved primes_pop_out_spikes.png')

# ---------------- FIGURE 3: J50 vs J100 ----------------
fig, ax = plt.subplots(figsize=(11.5, 6))
ax.plot(xs, psi_t, color='black', lw=1.7, label=r'true $\psi(x)$')
ax.plot(xs, psi_m50, color='#2980b9', lw=1.0, alpha=0.85, label=r'$J_{50}$ matrix arithmetic')
ax.plot(xs, psi_m100, color='#c0392b', lw=1.0, alpha=0.9, label=r'$J_{100}$ matrix arithmetic')
ax.set_xlabel('x', fontsize=12); ax.set_ylabel(r'$\psi(x)$', fontsize=12)
ax.set_title('More matrix entries = more primes: the trace sharpens with matrix order', fontsize=12)
ax.legend(loc='upper left', fontsize=10.5); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig('matrix_trace_order.png', dpi=150); plt.close()
print('saved matrix_trace_order.png')

# ---------------- FIGURE 4: determinant identity — zeta is a matrix determinant ----------------
us = np.linspace(-3, 3, 240)
logdet = np.array([np.linalg.slogdet(np.eye(n100) + u*u*J100)[1] for u in us])
# exact: 2 log |xi(1/2+u)/xi(1/2)| for real u
def xi_mp(s):
    return 0.5*s*(s-1)*mppi**(-s/2)*mpgamma(s/2)*mpzeta(s)
logdet_exact = np.array([2.0*float(mplog(abs(xi_mp(mpf('0.5')+u)/xi_mp(mpf('0.5'))))) for u in us[::4]])
fig, ax = plt.subplots(figsize=(11, 5.5))
ax.plot(us, logdet, color='#c0392b', lw=1.6, label=r'matrix: $\log\det(I+u^2 J_{100})$')
ax.plot(us[::4], logdet_exact, 'k--', lw=1.0, label=r'exact: $2\log|\,\xi(\frac{1}{2}+u)/\xi(\frac{1}{2})\,|$')
ax.set_xlabel(r'$u = s-\frac{1}{2}$ (real)', fontsize=12)
ax.set_ylabel(r'$\log\det$', fontsize=12)
ax.set_title(r'Converged band of $\xi(s)/\xi(\frac{1}{2})=\det(I+u^2J)$: the zeta function is a matrix determinant (finite $J_{100}$ resolves the band $|u|<3$ exactly; band widens with order)',
             fontsize=11.5, fontweight='bold')
ax.legend(fontsize=11); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig('zeta_is_determinant.png', dpi=150); plt.close()
print('saved zeta_is_determinant.png')
