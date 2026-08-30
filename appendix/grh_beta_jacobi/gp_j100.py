# -*- coding: utf-8 -*-
"""J50 vs J100 乘积矩阵解码高斯素数：尖峰范围扩展对比。"""
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from mpmath import mpf
fp='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
fm.fontManager.addfont(fp)
plt.rcParams['font.family']=fm.FontProperties(fname=fp).get_name()
plt.rcParams['mathtext.fontset']='dejavusans'

def build_J(a,b):
    a=np.array(a,float); b=np.array(b,float); n=len(a)
    J=np.zeros((n,n)); J[np.arange(n),np.arange(n)]=a
    J[np.arange(n-1),np.arange(1,n)]=b; J[np.arange(1,n),np.arange(n-1)]=b
    return J
def inv_sqrt_polar(J):
    L=np.linalg.cholesky(J); Linv=np.linalg.solve(L,np.eye(J.shape[0])); Y=Linv.copy()
    for _ in range(80):
        Yit=np.linalg.solve(Y.T,np.eye(Y.shape[0])); Yn=0.5*(Y+Yit)
        if np.linalg.norm(Yn-Y)<1e-14*np.linalg.norm(Y): Y=Yn; break
        Y=Yn
    return Y.T@Linv
def expm_sq(M,order=16):
    nrm=np.linalg.norm(M,ord=np.inf); q=max(0,int(np.ceil(np.log2(nrm+1e-300)))+1)
    E=M/(2.0**q); T=np.eye(M.shape[0],dtype=complex); term=T.copy()
    for k in range(1,order+1): term=term@E/k; T=T+term
    for _ in range(q): T=T@T
    return T
def setup(a,b):
    J=build_J(a,b); A=inv_sqrt_polar(J); n=J.shape[0]
    B=np.linalg.solve(0.5*np.eye(n)+1j*A, np.eye(n,dtype=complex))
    return A,B
def psi_trace(x,A,B):
    t=np.log(x); E=expm_sq(1j*A*t)
    return x-2.0*np.sqrt(x)*np.trace(E@B).real

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

R5=json.load(open('/tmp/grh_product/product_J50_results.json'))
R10=json.load(open('/tmp/grh_product/product_J100_results.json'))
A5,B5=setup(R5['alpha'],R5['b'])
A10,B10=setup(R10['alpha'],R10['b'])

Nmax=60
sp5=np.array([(psi_trace(n+0.5,A5,B5)-psi_trace(n-0.5,A5,B5)) if n>=2 else 0.0 for n in range(Nmax+1)])
sp10=np.array([(psi_trace(n+0.5,A10,B10)-psi_trace(n-0.5,A10,B10)) if n>=2 else 0.0 for n in range(Nmax+1)])
th=np.array([lamK(n) for n in range(Nmax+1)])

print(" n  理论Λ_K   J50尖峰    J100尖峰   J50err   J100err   事件")
events=[]
for n in range(2,Nmax+1):
    t=th[n]; e5=abs(sp5[n]-t); e10=abs(sp10[n]-t)
    tag=''
    if t>0.1:
        p=n
        # 判断事件类型
        if n%2==0 and lam(n)>0.1: tag='ramified 2-power'
        else:
            # split prime p≡1 mod4 或其幂；inert 静默@p,尖峰@p²
            ispow=False; base=None
            for p0 in range(2,n+1):
                if p0==2: continue
                if n%(p0**2)==0 and lam(n)>0.1: ispow=True; base=p0
            if chi4(n)==1 and lam(n)>0.1: tag='split p≡1'
            elif ispow and base%4==3: tag='inert p² (p≡3)'
            else: tag='?'
        events.append((n,t,sp5[n],sp10[n],e5,e10,tag))
        print(f"{n:3d} {t:8.4f} {sp5[n]:9.4f} {sp10[n]:9.4f} {e5:8.4f} {e10:8.4f}  {tag}")

# RMS（尖峰位置）
mask=th>0.1
rms5=np.sqrt(np.mean((sp5[mask]-th[mask])**2)); rms10=np.sqrt(np.mean((sp10[mask]-th[mask])**2))
# 有效区：|err|<0.5
v5=max((n for n in range(2,Nmax+1) if abs(sp5[n]-th[n])<0.5 and th[n]>0.1), default=0)
v10=max((n for n in range(2,Nmax+1) if abs(sp10[n]-th[n])<0.5 and th[n]>0.1), default=0)
print(f"\nRMS 尖峰区: J50={rms5:.4f}  J100={rms10:.4f}")
print(f"有效尖峰最远 n: J50={v5}  J100={v10}")

ns=np.arange(2,Nmax+1)
fig,ax=plt.subplots(2,1,figsize=(12,7.2),sharex=True)
for axi,sp,title,rms in [(ax[0],sp5,'乘积矩阵 J50（有效区 n≤25）',rms5),
                          (ax[1],sp10,'乘积矩阵 J100',rms10)]:
    axi.bar(ns,sp[2:],width=0.6,color='#2471a3',alpha=.75,label='矩阵迹尖峰')
    # 理论尖峰（split/inert p²/ramified）位置标记
    for n in ns:
        if th[n]>0.1:
            axi.plot(n,th[n],'r^',ms=7)
    axi.axhline(0,color='k',lw=.5)
    axi.set_ylabel('尖峰高度')
    axi.set_title(f'{title}　尖峰 RMS={rms:.3f}')
    axi.legend(['理论 Λ_K 尖峰（▼）','矩阵迹尖峰'],fontsize=9,loc='upper left')
    axi.grid(alpha=.3)
ax[1].set_xlabel('n（整数）')
fig.suptitle('Q(i) 高斯素数解码：J50 vs J100 乘积矩阵\n'
             'split p≡1 mod4 → 尖峰@2log p；inert p≡3 → 尖峰@p²；ramified 2 → 尖峰@2^k',
             fontsize=12)
fig.tight_layout()
fig.savefig('/tmp/grh_primes/gaussian_primes_j100.png',dpi=150,bbox_inches='tight')
print("saved /tmp/grh_primes/gaussian_primes_j100.png")
