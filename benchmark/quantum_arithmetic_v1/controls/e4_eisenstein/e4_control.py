# -*- coding: utf-8 -*-
"""E4 (Eisenstein E_k, k=4) 对照实验 — 正式版（650dps / M=1536，与 ζ J100 同管线）。

数学设定（标准教科书：自守 L 函数）：
  L(E4, s) = ζ(s) ζ(s-3)             Dirichlet 系数 σ_3(n)=Σ_{d|n} d^3（初等、算术确定）
  Λ(s) = (2π)^(-s) Γ(s) ζ(s) ζ(s-3)  函数方程 Λ(s)=Λ(4-s)，对称中心 s=2
  Λ 在 s=0,4 有真极点；取整函数 Ξ_E4(s) = s(s-4) Λ(s)（s(s-4) 在 s↔4-s 下不变）。

Ξ_E4 非平凡零点：s = 1/2 ± iγ_k 与 s = 7/2 ± iγ_k（γ_k = ζ 第 k 零点高）。
  → 零点分布在偏离中心线 s=2 的两条平行线上；t 变量(2+it)下零点为
    t = ±γ_k ∓ 3i/2，即矩测度支撑在复平面（非实轴）。
  → 矩 P[m] = 2 Σ_k Re[(γ_k - 3i/2)^(-2m)]：高阶被首零点 k=1 主导，
    含因子 cos(2m·θ_1)，θ_1=arctan(3/(2γ_1))≈0.106，m≥8 起振荡变负
    → Hankel 低阶即失正定 → Jacobi b^2 在第 4 阶变负。
本实验作为「判别力对照」：ζ/β/Δ（零点全在临界线）b^2 全正、链存活至几十~上百阶；
E4（零点偏离临界线）链在第 4 阶崩溃 —— 流水线对「零点是否全在临界线上」有判别力。

采样：圆心 s=2、R=2，网格偏移半步(+π/M)避开 s=0/4 两个被消极点；
矩递推 + Gram-Schmidt 与 ζ/乘积管线逐行相同。
"""
import mpmath as mp, json, time, os
import numpy as np
mp.mp.dps = 650
R = mp.mpf('2.0'); M = 1536; NN = 60; JMAX = 130
t0 = time.time()
OUT = os.path.dirname(os.path.abspath(__file__))

def xi_e4(s):
    return s*(s-4) * (2*mp.pi)**(-s) * mp.gamma(s) * mp.zeta(s) * mp.zeta(s-3)

print(f"E4 control: {M} pts R=2, {mp.mp.dps}dps, NN={NN}", flush=True)
fr = []
for k in range(M):
    th = 2*mp.pi*k/M + mp.pi/M
    s = mp.mpf(2) + R*mp.e**(1j*th)
    fr.append(xi_e4(s))
    if k % 192 == 0:
        print("  pt", k, round(time.time()-t0, 1), "s", flush=True)
print("sampled", round(time.time()-t0, 1), "s", flush=True)
mags = [float(abs(x)) for x in fr]
print("|f| min/median/max: %.3e / %.3e / %.3e" %
      (min(mags), float(np.median(mags)), max(mags)), flush=True)

sigma0 = sum(fr)/M
print("sigma0 =", mp.nstr(sigma0, 20), flush=True)
print("theory  = 0.0138888... = 1/72", flush=True)

d = {}; c_log = {}
for n in range(1, JMAX+1):
    j = 2*n
    ss = mp.mpc(0)
    for k in range(M):
        ss += fr[k]*mp.e**(-1j*2*mp.pi*j*k/M)
    sig = ss/M/R**j
    d[n] = mp.re(sig/sigma0)
    sm = mp.mpf(0)
    for i in range(1, n):
        sm += i*c_log[i]*d[n-i]
    c_log[n] = d[n] - sm/n
P = [None] + [((-1)**(m+1))*m*c_log[m] for m in range(1, JMAX+1)]
for idx in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100):
    print(f"  P[{idx}] = {mp.nstr(P[idx], 8)}", flush=True)

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
    if k >= 2 and bad_first is None:
        bv = mp.re(bsqs[-1])
        if bv <= 0:
            bad_first = k+1
    if k <= 8:
        bv = mp.re(bsqs[-1]) if bsqs else mp.nan
        print(f"  k={k} alpha={mp.nstr(mp.re(ak),10)} b^2[{k-1}]={mp.nstr(bsqs[-1],10) if bsqs else '-'}", flush=True)

# b^2 序列诊断：实部符号（虚部应为数值噪声）
bsq_real = [float(mp.re(x)) for x in bsqs]
bsq_imag = [float(abs(mp.im(x))) for x in bsqs]
print("first non-positive b^2 at k+1 =", bad_first, flush=True)
print("b^2 real signs (first 20):",
      ["+" if x > 0 else "-" for x in bsq_real[:20]], flush=True)
print("max |Im b^2| / |Re b^2| (first 10):",
      [f"{bsq_imag[i]/max(abs(bsq_real[i]),1e-300):.1e}" for i in range(min(10,len(bsq_imag)))], flush=True)

# 有效链 = 从起点开始的最长 b^2>0 前缀
nvalid = 0
for x in bsq_real:
    if x > 0: nvalid += 1
    else: break
print("valid b prefix length =", nvalid, " (Jacobi size =", nvalid+1, ")", flush=True)

result = {"object": "Eisenstein E4 control",
          "L_function": "L(E4,s)=zeta(s)zeta(s-3), coefficients sigma_3(n)",
          "completed_Lambda": "(2pi)^(-s) Gamma(s) zeta(s) zeta(s-3), FE Lambda(s)=Lambda(4-s), center s=2",
          "entire_function": "Xi_E4(s) = s(s-4) Lambda(s)",
          "zeros": "s = 1/2 +/- i gamma_k and 7/2 +/- i gamma_k (gamma_k = zeta zeros); t-plane: t = +/-gamma_k -/+ 3i/2",
          "dps": 650, "M": M, "R": 2.0, "NN_target": NN, "JMAX": JMAX,
          "sigma0": float(mp.re(sigma0)), "sigma0_theory": 1.0/72.0,
          "P_first10": [float(mp.re(P[m])) for m in range(1, 11)],
          "bad_first_b2": bad_first,
          "valid_b_prefix": nvalid,
          "jacobi_size": nvalid+1,
          "b2_real_first20": bsq_real[:20],
          "b2_imag_max_ratio_first10": bsq_imag[:10],
          "interpretation": (
              "E4 zeros lie on TWO lines Re(s)=1/2 and Re(s)=7/2 (off the symmetry center s=2); "
              "the moment measure is complex-supported, Hankel loses positive definiteness at low order, "
              "Jacobi b^2 becomes non-positive at k+1=4. Contrast: zeta/beta/Delta (zeros ON critical line) "
              "give real positive measures with b^2>0 through dozens/hundreds of orders. "
              "The pipeline discriminates 'all zeros on the critical line' vs 'off-line zeros'.")}
json.dump(result, open(os.path.join(OUT, "e4_control_results.json"), "w"), indent=1)
print("DONE -> e4_control_results.json", round(time.time()-t0, 1), "s", flush=True)
