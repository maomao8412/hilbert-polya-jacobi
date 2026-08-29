# -*- coding: utf-8 -*-
"""Full beta J50 pipeline at 600 digits: samples600 -> DFT(even) -> log recurrence
-> moments T1..T120 -> Stieltjes -> J50 -> eigenvalues vs beta zeros (mpmath double scan)."""
import json, csv, time
import sympy as sp
from sympy import pi, cos, sin, N, S, sqrt, Float

DPS=600; NN=50; M=768; R=sp.Rational(3,10)
t0=time.time()
dd=json.load(open("/tmp/grh_beta/samples600.json"))
fr=[Float(x,DPS) for x in dd["fr"]]; fi=[Float(x,DPS) for x in dd["fi"]]
JMAX=240  # T_120 needs c_240
ev_j=list(range(0,JMAX+1,2))
print("trig ...", flush=True)
cos_t={}; sin_t={}
for j in ev_j:
    cj=[]; sj=[]
    for k in range(M):
        ph=2*pi*((j*k)%M)/M
        cj.append(cos(ph)); sj.append(sin(ph))
    cos_t[j]=cj; sin_t[j]=sj
print(" trig", round(time.time()-t0,1), flush=True)
dre=[S.Zero]*(JMAX+1)
Rj=S.One
for j in ev_j:
    cr=S.Zero
    ct=cos_t[j]; st=sin_t[j]
    for k in range(M):
        cr+=fr[k]*ct[k]+fi[k]*st[k]
    dre[j]=N(cr/M/Rj,DPS)
    Rj*=R*R
    if j%60==0: print("  j",j, round(time.time()-t0,1), flush=True)
print(" d2 =", N(dre[2],14), flush=True)
c=[S.Zero]*(JMAX+1)
for n in ev_j[1:]:
    s=S.Zero
    for k in ev_j:
        if k>=n: break
        if k==0: continue
        s+=k*c[k]*dre[n-k]
    c[n]=dre[n]-s/n
T=[N(((-1)**(m+1))*m*c[2*m],DPS) for m in range(1,2*NN+21)]  # T1..T120
print(" T1 =",N(T[0],16)," T50 =",N(T[49],6)," T100 =",N(T[99],6)," T120 =",N(T[119],6), flush=True)

def Tm(m): return Float(str(T[m-1]), DPS) if not isinstance(T[m-1],Float) else T[m-1]
def inner(p,q):
    r=S.Zero
    for i,pi in enumerate(p):
        if pi==0: continue
        for j,qj in enumerate(q):
            if qj==0: continue
            r+=pi*qj*Tm(i+j+1)
    return r
n=NN
polys=[[S.One]]; norms=[inner([S.One],[S.One])]
alphas=[]; betas_sq=[]
for k in range(1,n+1):
    xp=[S.Zero]+polys[k-1]
    alpha_k=inner(xp,polys[k-1])/norms[k-1]
    pi=list(xp)
    for i in range(len(polys[k-1])):
        pi[i]-=alpha_k*polys[k-1][i]
    if k>=2:
        bsq=norms[k-1]/norms[k-2]
        for i in range(len(polys[k-2])):
            pi[i]-=bsq*polys[k-2][i]
        betas_sq.append(bsq)
    sig=inner(pi,pi)
    alphas.append(alpha_k); polys.append(pi); norms.append(sig)
    if k%5==0:
        print(f"  n={k:2d} alpha={N(alpha_k,6)} bsq={N(betas_sq[-1] if betas_sq else S.Zero,6)} sig={N(sig,4)} {round(time.time()-t0,1)}s", flush=True)
bad=[k+2 for k,x in enumerate(betas_sq) if N(x,40)<=0]
print("non-positive b_sq at n:", bad, flush=True)
a=[N(x,DPS) for x in alphas]
b=[N(sqrt(x),DPS) for x in betas_sq]
import numpy as np
ad=np.array([float(x) for x in a]); bd=np.array([float(x) for x in b])
evals=np.sort(np.linalg.eigvalsh(np.diag(ad)+np.diag(bd,1)+np.diag(bd,-1)))[::-1]
inv=np.sqrt(1.0/np.maximum(evals,1e-300))
import mpmath as mp
mp.dps=40
def Lam_t(t):
    s=mp.mpc(0.5,t)
    return float(mp.re((mp.pi/4)**(-(s+1)/2)*mp.gamma((s+1)/2)*mp.dirichlet(s,[0,1,0,-1])))
san=Lam_t(6.0209489046975975)
zs=[]; t=3.0; dt=0.02; prev=Lam_t(t)
while t<400 and len(zs)<NN+15:
    t+=dt; cur=Lam_t(t)
    if prev*cur<0:
        lo,hi=t-dt,t
        for _ in range(45):
            mid=(lo+hi)/2
            if Lam_t(lo)*Lam_t(mid)<0: hi=mid
            else: lo=mid
        zs.append((lo+hi)/2)
    prev=cur
rows=[]; nlock=0
for k in range(min(NN,len(zs))):
    err=abs(inv[k]-zs[k])/zs[k]
    if err<1e-4: nlock+=1
    rows.append({"k":k+1,"lambda":float(evals[k]),"inv_sqrt":float(inv[k]),"gamma":zs[k],"rel_err":float(err)})
    if k<6 or k in (9,14,19,24,29,34,39,44,49):
        print(f"  n={k+1:2d} inv={inv[k]:.10f} gamma={zs[k]:.10f} rel={err:.2e}", flush=True)
print("LOCKED:",nlock," sanity:",san," zeros:",len(zs), flush=True)
with open("/tmp/grh_beta/beta_J50_matrix.csv","w",newline="") as f:
    w=csv.writer(f)
    for i in range(NN):
        row=["0.0"]*NN; row[i]=str(a[i])
        if i>0: row[i-1]=str(b[i-1])
        if i<NN-1: row[i+1]=str(b[i])
        w.writerow(row)
json.dump({"N":NN,"dps":DPS,"M":M,"R":"0.3",
  "moments_T":[str(x) for x in T[:100]],
  "diag_a":[str(x) for x in a],"offdiag_b":[str(x) for x in b],
  "eigenvalues":[float(x) for x in evals],"inv_sqrt":[float(x) for x in inv],
  "beta_zeros_independent":zs,"comparison":rows,"n_locked":nlock,
  "sanity_Lam_at_gamma1":san,
  "method":"Cauchy DFT 1536pts/|z|=0.3, 600-digit sympy; log-series recurrence; Stieltjes monic orthogonal polynomial recurrence; no zero input",
  "elapsed":time.time()-t0},
  open("/tmp/grh_beta/beta_J50_results.json","w"),indent=1)
print("SAVED", round(time.time()-t0,1), flush=True)
