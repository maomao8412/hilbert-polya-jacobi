# -*- coding: utf-8 -*-
"""600-digit, 1536-point sampling of Lambda_beta(1/2+z)/Lambda_beta(1/2) on |z|=0.3."""
import json, time
import sympy as sp
from sympy import I, pi, gamma, cos, sin, exp, N, S, Rational, zeta

DPS=600; M=768; R=Rational(3,10)
def beta_s(s):
    return S(4)**(-s)*(zeta(s, Rational(1,4))-zeta(s, Rational(3,4)))
def Lam(s):
    return (pi/4)**(-(s+1)/2)*gamma((s+1)/2)*beta_s(s)
t0=time.time()
L0=N(Lam(S.Half),DPS)
print("L0 done", round(time.time()-t0,1), flush=True)
fr=[]; fi=[]
for k in range(M):
    th=2*pi*k/M
    z=R*exp(I*th)
    f=N(Lam(S.Half+z)/L0, DPS)
    fr.append(str(sp.re(f))); fi.append(str(sp.im(f)))
    if k%96==0: print(" sample",k, round(time.time()-t0,1), flush=True)
json.dump({"fr":fr,"fi":fi,"dps":DPS,"M":M,"R":"0.3"},
          open("/tmp/grh_beta/samples600.json","w"))
print("SAVED", round(time.time()-t0,1), flush=True)
