# -*- coding: utf-8 -*-
"""E4 (Eisenstein E_k, k=4) 对照实验 — 低精度侦察版。

数学设定（标准教科书：自守 L 函数）：
  L(E4, s) = ζ(s) ζ(s-3)            Dirichlet 系数 σ_3(n)=Σ_{d|n} d^3（初等、算术确定）
  Λ(s) = (2π)^(-s) Γ(s) ζ(s) ζ(s-3) 函数方程 Λ(s)=Λ(4-s)，对称中心 s=2
  Λ 在 s=0,4 有真极点（Γ 极点 / ζ(1) 极点未被平凡零点抵消）；
  取整函数 Ξ_E4(s) = s(s-4) Λ(s)（s(s-4) 在 s↔4-s 下不变，FE 保持）。

Ξ_E4 零点（L 的非平凡零点）：s = 1/2 ± iγ_k 与 s = 7/2 ± iγ_k，γ_k = ζ 第 k 零点高
（两族高度完全重合 —— 相关叠加，与 ξ_K=ξ·Λ_β 的独立叠加不同）。
采样圆中心 s=2、半径 R=2；网格偏移半步避开 s=0/4 两个被消去极点的采样位置。

矩与 Jacobi 部分与 ζ/乘积管线逐行相同（R 圆 Cauchy 采样 → log 递推矩 → GS）。
输出：b^2 存活阶数、锁定能级、gap-ratio ⟨r⟩、IPR，并与 ζ 期望值对比。
"""
import mpmath as mp, json, time, os
import numpy as np
mp.mp.dps = 230
R = mp.mpf('2.0'); M = 512; NN = 50; JMAX = 120
t0 = time.time()
OUT = os.path.dirname(os.path.abspath(__file__))

def xi_e4(s):
    # 整函数：s(s-4) * (2π)^(-s) * Γ(s) * ζ(s) * ζ(s-3)
    return s*(s-4) * (2*mp.pi)**(-s) * mp.gamma(s) * mp.zeta(s) * mp.zeta(s-3)

print(f"scout: E4 Jacobi, {M} pts R=2, {mp.mp.dps}dps, NN={NN}", flush=True)
# 半格偏移：s_k = 2 + R exp(i(θ_k + π/M))，避开 s=0 与 s=4
fr = []
for k in range(M):
    th = 2*mp.pi*k/M + mp.pi/M
    s = mp.mpf(2) + R*mp.e**(1j*th)
    fr.append(xi_e4(s))
    if k % 128 == 0:
        print("  pt", k, round(time.time()-t0, 1), "s", flush=True)
print("sampled", round(time.time()-t0, 1), "s", flush=True)
# 数值健全性检查：采样值量级是否合理
mags = [float(abs(x)) for x in fr]
print("|f| min/median/max: %.3e / %.3e / %.3e" %
      (min(mags), float(np.median(mags)), max(mags)), flush=True)

sigma0 = sum(fr)/M
print("sigma0 =", mp.nstr(sigma0, 15), flush=True)

d = {}; c_log = {}
for n in range(1, JMAX+1):
    j = 2*n
    ss = mp.mpc(0)
    for k in range(M):
        ss += fr[k]*mp.e**(-1j*2*mp.pi*j*k/M)
    sig = ss/M/R**j
    d[n] = mp.re(sig/sigma0)   # FE 对称保证主系数为实；半格偏移仅给混叠项加符号
    sm = mp.mpf(0)
    for i in range(1, n):
        sm += i*c_log[i]*d[n-i]
    c_log[n] = d[n] - sm/n
P = [None] + [((-1)**(m+1))*m*c_log[m] for m in range(1, JMAX+1)]
print("P1 =", mp.nstr(P[1], 10), " P6 =", mp.nstr(P[6], 6),
      " P50 =", mp.nstr(P[50], 6), " P100 =", mp.nstr(P[100], 6), flush=True)

def Tm(m): return P[m]
def inner(p, q):
    r = mp.mpf(0)
    for i, pi_ in enumerate(p):
        if pi_ == 0: continue
        for j, qj in enumerate(q):
            if qj == 0: continue
            r += pi_*qj*Tm(i+j+1)
    return r

polys = [[mp.mpf(1)]]; norms = [inner(polys[0], polys[0])]
alphas = []; bsqs = []; bad_first = None
for k in range(1, NN+1):
    xp = [mp.mpf(0)] + polys[k-1]
    ak = inner(xp, polys[k-1])/norms[k-1]
    pk = list(xp)
    for i in range(len(polys[k-1])):
        pk[i] -= ak*polys[k-1][i]
    if k >= 2:
        bsq = norms[k-1]/norms[k-2]
        for i in range(len(polys[k-2])):
            pk[i] -= bsq*polys[k-2][i]
        bsqs.append(bsq)
    sig = inner(pk, pk)
    alphas.append(ak); polys.append(pk); norms.append(sig)
    if k >= 2 and bsqs[-1] <= 0 and bad_first is None:
        bad_first = k+1
    if k % 10 == 0:
        bval = mp.sqrt(bsqs[-1]) if (bsqs and bsqs[-1] > 0) else float('nan')
        print(f"  n={k:3d} a={mp.nstr(ak,6)} b={mp.nstr(bval,6) if bsqs and bsqs[-1]>0 else 'NEG'} {round(time.time()-t0,0)}s", flush=True)

nvalid = len(bsqs)
while nvalid > 0 and bsqs[nvalid-1] <= 0:
    nvalid -= 1
print("first non-positive b^2 at k+1 =", bad_first, "| valid b count =", nvalid, flush=True)
NNv = nvalid + 1
a = [float(x) for x in alphas[:NNv]]
b = [float(mp.sqrt(x)) for x in bsqs[:nvalid]]
Jmat = np.diag(a) + np.diag(b, 1) + np.diag(b, -1)
ev = np.sort(np.linalg.eigvalsh(Jmat))[::-1]
# t 变量：s=2+it，矩支撑 1/t^2；t_k = 1/sqrt(λ)；物理零点高 g=sqrt(t^2-9/4)
ts = 1.0/np.sqrt(np.maximum(ev, 1e-300))
ts = np.sort(ts)
gammas = np.sqrt(np.maximum(ts**2 - 2.25, 0.0))

# 事后对照 ζ 零点（不进构造）
zeta_zeros = [float(mp.im(mp.zetazero(k))) for k in range(1, 61)]
rows = []
for n, g in enumerate(gammas, 1):
    if g <= 0: continue
    dz = min((abs(g-z), z) for z in zeta_zeros)
    rows.append((n, float(g), float(dz[1]), float(dz[0]/dz[1])))
locked = [r for r in rows if r[3] < 1e-4]
print(f"J size={NNv}; gamma>0: {len(rows)}; locked(rel<1e-4): {len(locked)}", flush=True)
for r in rows[:15]:
    print(f"n={r[0]:3d} g={r[1]:.10f} ref={r[2]:.10f} relerr={r[3]:.2e}", flush=True)

# 谱统计：gap ratio（锁定能级，展开到 s 平面 t 间距不需要——比较 g 间距即可，两族同高）
def gap_ratio(xs):
    xs = np.sort(np.array(xs))
    d = np.diff(xs)
    r = np.minimum(d[:-1], d[1:])/np.maximum(d[:-1], d[1:])
    return float(np.mean(r)), r
if len(locked) >= 8:
    lg = sorted(r[1] for r in locked)
    rmean, rvals = gap_ratio(lg)
    print(f"locked gap-ratio <r> = {rmean:.4f}  (Poisson .386 / GUE .599 / GOE .536)", flush=True)
else:
    rmean = None; rvals = []
    print("locked too few for gap ratio", flush=True)

json.dump({"scout": True, "dps": 230, "M": M, "R": 2.0, "NN": NN, "JMAX": JMAX,
           "bad_first": bad_first, "n_valid": NNv,
           "alphas": a, "b": b,
           "gammas": [float(x) for x in gammas],
           "rows": rows, "n_locked": len(locked),
           "gap_ratio_locked": rmean,
           "P1": float(P[1]), "P6": float(P[6])},
          open(os.path.join(OUT, "e4_scout_results.json"), "w"), indent=1)
print("DONE -> e4_scout_results.json", round(time.time()-t0, 1), "s", flush=True)
