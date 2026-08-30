# -*- coding: utf-8 -*-
"""ζ 单矩阵 J50 vs ζ_K 乘积矩阵 J50：解码素数对比。有效区 n≤25（50阶矩阵截断点之后振荡）。"""
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpmath import mpf

def build_J(alphas, betas_b):
    a=np.array(alphas,float); bb=np.array(betas_b,float)
    n=len(a); J=np.zeros((n,n)); J[np.arange(n),np.arange(n)]=a
    J[np.arange(n-1),np.arange(1,n)]=bb; J[np.arange(1,n),np.arange(n-1)]=bb
    return J
def inv_sqrt_polar(J):
    n=J.shape[0]; L=np.linalg.cholesky(J)
    Linv=np.linalg.solve(L,np.eye(n)); Y=Linv.copy()
    for _ in range(60):
        Yit=np.linalg.solve(Y.T,np.eye(n)); Yn=0.5*(Y+Yit)
        if np.linalg.norm(Yn-Y)<1e-14*np.linalg.norm(Y): Y=Yn; break
        Y=Yn
    return Y.T@Linv
def expm_sq(M,order=14):
    nrm=np.linalg.norm(M,ord=np.inf); q=max(0,int(np.ceil(np.log2(nrm+1e-300)))+1)
    E=M/(2.0**q); T=np.eye(M.shape[0],dtype=complex); term=T.copy()
    for k in range(1,order+1): term=term@E/k; T=T+term
    for _ in range(q): T=T@T
    return T

def setup(alphas, betas_b):
    J=build_J(alphas,betas_b); A=inv_sqrt_polar(J)
    n=J.shape[0]
    B=np.linalg.solve(0.5*np.eye(n)+1j*A, np.eye(n,dtype=complex))
    return A,B
def psi_trace(x,A,B,C=0.0):
    t=np.log(x); E=expm_sq(1j*A*t)
    return x - 2.0*np.sqrt(x)*np.trace(E@B).real - C

# ζ 单矩阵
cp=json.load(open('/tmp/hpgit/data/J50_checkpoint.json'))
Az,Bz=setup([float(mpf(s)) for s in cp['alphas']],
            [float(mpf(s))**0.5 for s in cp['betas_sq']])
# ζ_K 乘积矩阵
R=json.load(open('/tmp/grh_product/product_J50_results.json'))
AK,BK=setup(R['alpha'], R['b'])

def lam(n):
    m=n
    for p in range(2,int(n**0.5)+1):
        if m%p==0:
            while m%p==0: m//=p
            return np.log(p) if m==1 else 0.0
    return np.log(n) if n>1 else 0.0
def chi4(n):
    n%=4; return 1 if n==1 else (-1 if n==3 else 0)
def lamK(n): return lam(n)*(1+chi4(n))

Nmax=25  # 有效区
LAM=np.array([lam(n) for n in range(Nmax+1)])
LAMK=np.array([lamK(n) for n in range(Nmax+1)])
sp_z=np.array([(psi_trace(n+0.5,Az,Bz)-psi_trace(n-0.5,Az,Bz)) if n>=2 else 0.0 for n in range(Nmax+1)])
sp_K=np.array([(psi_trace(n+0.5,AK,BK)-psi_trace(n-0.5,AK,BK)) if n>=2 else 0.0 for n in range(Nmax+1)])

print("=== 有效区 n=2..25 尖峰对比 ===")
print(f"{'n':>3} {'χ4':>3} {'Λ(有理)':>8} {'ζ矩阵':>8} {'Λ_K(高斯)':>9} {'ζ_K矩阵':>8}")
for n in range(2,Nmax+1):
    print(f"{n:>3} {chi4(n):>3} {LAM[n]:>8.3f} {sp_z[n]:>8.3f} {LAMK[n]:>9.3f} {sp_K[n]:>8.3f}")

# 相关/误差量化
mask=np.arange(2,Nmax+1)
rms_z=np.sqrt(np.mean((sp_z[mask]-LAM[mask])**2))
rms_K=np.sqrt(np.mean((sp_K[mask]-LAMK[mask])**2))
# 惰性素数 3,7,11,19,23 抑制比
inert=[3,7,11,19,23]
supp=np.mean([sp_z[n] for n in inert]), np.mean([max(sp_K[n],0) for n in inert])
print(f"\nζ矩阵解码有理Λ RMS={rms_z:.3f};  ζ_K矩阵解码高斯Λ_K RMS={rms_K:.3f}")
print(f"惰性素数{inert}平均尖峰: ζ矩阵={supp[0]:.3f} -> 乘积矩阵={supp[1]:.3f} (应被抑制到0)")

# ---------- 图：尖峰对比 ----------
ns=np.arange(2,Nmax+1)
fig,ax=plt.subplots(2,1,figsize=(12.5,8.5),gridspec_kw={'height_ratios':[1,1]},sharex=True)
ax[0].bar(ns-0.28,LAM[ns],width=0.28,color='black',alpha=.7,label=r'true von Mangoldt $\Lambda(n)$ (rational primes)')
ax[0].bar(ns,np.maximum(sp_z[ns],0),width=0.28,color='#c0392b',alpha=.85,label=r'ζ matrix $J_{50}$ trace spikes')
ax[0].bar(ns+0.28,LAMK[ns],width=0.28,color='#2c3e50',alpha=.0)  # 占位对齐
ax[0].set_ylabel('spike weight',fontsize=12)
ax[0].set_title(r'ζ matrix decodes RATIONAL primes: spikes at every prime (3,7,11 included)',fontsize=12,fontweight='bold')
ax[0].legend(fontsize=9.5,loc='upper left'); ax[0].grid(alpha=.25,axis='y')
for n in inert: ax[0].annotate(f'{n}',(n,LAM[n]),textcoords='offset points',xytext=(0,4),ha='center',fontsize=8,color='#c0392b')

ax[1].bar(ns-0.28,LAMK[ns],width=0.28,color='#1a4d8f',alpha=.75,label=r'true Gaussian-prime weight $\Lambda_K(n)=\Lambda(n)(1+\chi_4(n))$')
ax[1].bar(ns,np.maximum(sp_K[ns],0),width=0.28,color='#27ae60',alpha=.85,label=r'product matrix $\xi_K=\xi\cdot\Lambda_\beta$ trace spikes')
ax[1].set_ylabel('spike weight',fontsize=12); ax[1].set_xlabel('n (norm)',fontsize=12)
ax[1].set_title(r'Product matrix decodes GAUSSIAN primes: split p≡1 mod4 → 2log p at n=p; inert p≡3 mod4 → silent at n=p (3,7,11,19,23), spike at p²; ramified 2 → log2',
                fontsize=10.5,fontweight='bold')
ax[1].legend(fontsize=9.5,loc='upper left'); ax[1].grid(alpha=.25,axis='y')
# 标注惰性空缺
for n in inert:
    ax[1].annotate(f'{n}\n(silent)',(n,0.15),textcoords='offset points',xytext=(0,-2),ha='center',fontsize=7.5,color='#999')
for n in [5,13,17]: ax[1].annotate(f'{n}',(n,LAMK[n]),textcoords='offset points',xytext=(0,4),ha='center',fontsize=8,color='#1a4d8f')
ax[1].annotate('9=p² (3 inert)',(9,LAMK[9]),textcoords='offset points',xytext=(0,4),ha='center',fontsize=8,color='#1a4d8f')
plt.tight_layout(); plt.savefig('/tmp/grh_primes/gaussian_prime_decode.png',dpi=150); plt.close()
print("\nsaved gaussian_prime_decode.png")
