#!/usr/bin/env python3
"""Quick large-r and large-q verification."""
import mpmath as mp
mp.mp.dps = 15

def make_G(chars, q, a, bar=None):
    def L(s): return mp.dirichlet(s, chars)
    def Lam(s): return (q/mp.pi)**((s+a)/2)*mp.gamma((s+a)/2)*L(s)
    if bar:
        def Lb(s): return mp.dirichlet(s, bar)
        def Lamb(s): return (q/mp.pi)**((s+a)/2)*mp.gamma((s+a)/2)*Lb(s)
        return lambda s: Lam(s)*Lamb(s)
    return lambda s: Lam(s)**2

def T(F, r, th):
    s = mp.mpf('0.5')+r*mp.e**(1j*th)
    return mp.re(1j*r*mp.e**(1j*th)*mp.diff(lambda z: mp.log(F(z)), s, 1))

# Small q chars
Gs = {
    'q=3': make_G([0,1,-1],3,1),
    'q=4': make_G([0,1,0,-1],4,1),
    'q=5e': make_G([0,1,-1,-1,1],5,0),
    'q=5c': make_G([0,1,1j,-1j,-1],5,1,[0,1,-1j,1j,-1]),
    'q=7': make_G([0,1,1,-1,1,-1,-1],7,1),
    'q=8e': make_G([0,1,0,-1,0,-1,0,1],8,0),
    'q=8o': make_G([0,1,0,1,0,-1,0,-1],8,1),
}

print("=== Large r test (r=100,500,1000) ===")
for name, G in Gs.items():
    mx = -mp.inf
    for rv in [100, 500, 1000]:
        for tj in range(1, 16):
            th = mp.pi*tj/32
            try:
                t = T(G, mp.mpf(rv), th)
                if t > mx: mx = t
            except: pass
    print(f"  {name:6s}: maxT={float(mx):.6e} [{'PASS' if mx<0 else 'FAIL'}]")

print("\n=== q=101 Legendre (odd) ===")
chi101 = [0] + [1 if pow(n,50,101)==1 else -1 for n in range(1,101)]
G101 = make_G(chi101, 101, 1)
mx = -mp.inf; mpt = None
for rv in [0.5, 1, 2, 5, 10, 50, 100, 500]:
    for tj in range(1, 16):
        th = mp.pi*tj/32
        try:
            t = T(G101, mp.mpf(rv), th)
            if t > mx: mx = t; mpt = (rv, float(th))
        except: pass
print(f"  maxT={float(mx):.6e} at r={mpt[0]},th={mpt[1]:.3f} [{'PASS' if mx<0 else 'FAIL'}]")

# Also test q=11, q=13 (small odd primes)
print("\n=== q=11, q=13 ===")
for qq in [11, 13]:
    chi = [0] + [1 if pow(n,(qq-1)//2,qq)==1 else -1 for n in range(1,qq)]
    # Legendre symbol: chi(-1) = (-1)^{(q-1)/2}, so a=0 if q≡1 mod4, a=1 if q≡3 mod4
    aa = (qq-1)//2 % 2
    Gq = make_G(chi, qq, aa)
    mx = -mp.inf
    for rv in [0.5, 1, 2, 5, 10, 50, 100]:
        for tj in range(1, 16):
            th = mp.pi*tj/32
            try:
                t = T(Gq, mp.mpf(rv), th)
                if t > mx: mx = t
            except: pass
    print(f"  q={qq} a={aa}: maxT={float(mx):.6e} [{'PASS' if mx<0 else 'FAIL'}]")

# Test the Archimedean function f(theta; r, q) = sin(theta)*log(qr/(2pi)) + theta*cos(theta)
print("\n=== Archimedean f(theta;r,q) sign check ===")
import math
for qq in [1, 3, 4, 5, 7, 8, 11, 13, 101]:
    # Find minimum f over theta in (0, pi/2) for various r
    for rv in [0.5, 1, 2, 5, 7]:
        fmin = min(math.sin(th)*math.log(qq*rv/(2*math.pi)) + th*math.cos(th) 
                   for th in [i*math.pi/200 for i in range(1,100)])
        if fmin < 0:
            print(f"  q={qq:3d} r={rv}: min f = {fmin:.4f} < 0 (Archimedean alone insufficient)")

print("\nDONE")
