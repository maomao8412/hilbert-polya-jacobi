# -*- coding: utf-8 -*-
"""EN version: J50 vs J100 product matrix decoding Gaussian primes."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='DejaVu Sans'
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

mask=th>0.1
rms5=np.sqrt(np.mean((sp5[mask]-th[mask])**2)); rms10=np.sqrt(np.mean((sp10[mask]-th[mask])**2))

ns=np.arange(2,Nmax+1)
fig,ax=plt.subplots(2,1,figsize=(12,7.2),sharex=True)
for axi,sp,title,rms in [(ax[0],sp5,'Product matrix J50 (valid region n<=25)',rms5),
                          (ax[1],sp10,'Product matrix J100',rms10)]:
    hb=axi.bar(ns,sp[2:],width=0.6,color='#2471a3',alpha=.75)
    for n in ns:
        if th[n]>0.1:
            axi.plot(n,th[n],'r^',ms=7)
    ht,=axi.plot([],[],'r^',ms=7)
    axi.axhline(0,color='k',lw=.5)
    axi.set_ylabel('spike height')
    axi.set_title('%s   spike RMS=%.3f'%(title,rms))
    axi.legend([ht,hb],['theoretical $\\Lambda_K$ spikes','matrix-trace spikes'],fontsize=9,loc='upper left')
    axi.grid(alpha=.3)
ax[1].set_xlabel('n (integer)')
fig.suptitle('Q(i) Gaussian-prime decoding: J50 vs J100 product matrix\n'
             'split p=1 mod4 -> spike @ 2log p; inert p=3 mod4 -> spike @ p^2; ramified 2 -> spike @ 2^k',
             fontsize=12)
fig.tight_layout()
fig.savefig('/tmp/img_en/gaussian_primes_j100_en.png',dpi=150,bbox_inches='tight')
print("saved gaussian_primes_j100_en; RMS J50=%.3f J100=%.3f"%(rms5,rms10))
