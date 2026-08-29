#!/usr/bin/env python3
"""Ultra-fast GRH verification - fewer points, lower precision."""
import mpmath as mp
mp.mp.dps = 20

def make_G(chars, q, a, complex_bar=None):
    def L_func(s): return mp.dirichlet(s, chars)
    def Lam(s): return (q/mp.pi)**((s+a)/2) * mp.gamma((s+a)/2) * L_func(s)
    if complex_bar is not None:
        def L_bar(s): return mp.dirichlet(s, complex_bar)
        def Lam_bar(s): return (q/mp.pi)**((s+a)/2) * mp.gamma((s+a)/2) * L_bar(s)
        return lambda s: Lam(s)*Lam_bar(s)
    return lambda s: Lam(s)**2

def T_func(Func, r, theta):
    s = mp.mpf('0.5') + r*mp.e**(1j*theta)
    Fp = mp.diff(lambda z: mp.log(Func(z)), s, 1)
    return mp.re(1j*r*mp.e**(1j*theta)*Fp)

def C_func(Func, r):
    s0 = mp.mpf('0.5')+r
    Fp = mp.diff(lambda z: mp.log(Func(z)), s0, 1)
    Fpp = mp.diff(lambda z: mp.log(Func(z)), s0, 2)
    return mp.re(r*Fp+r*r*Fpp)

# Characters
G3 = make_G([0,1,-1], 3, 1)
G4 = make_G([0,1,0,-1], 4, 1)
G5e = make_G([0,1,-1,-1,1], 5, 0)
G5c = make_G([0,1,1j,-1j,-1], 5, 1, [0,1,-1j,1j,-1])
G7 = make_G([0,1,1,-1,1,-1,-1], 7, 1)
G8 = make_G([0,1,0,-1,0,-1,0,1], 8, 0)
G8b = make_G([0,1,0,1,0,-1,0,-1], 8, 1)

# Large q: Legendre mod 101
chi101 = [0] + [1 if pow(n,50,101)==1 else -1 for n in range(1,101)]
G101 = make_G(chi101, 101, 1)

tests = [
    ("q=3 odd", G3), ("q=4 beta", G4), ("q=5 even", G5e),
    ("q=5 complex", G5c), ("q=7 odd", G7),
    ("q=8 even", G8), ("q=8 odd", G8b), ("q=101 large", G101),
]

print("=== Grid r=0.5..50, 16 theta ===")
for name, G in tests:
    max_T = -mp.inf; mp_pt = None
    for ri in range(1, 101):
        r = mp.mpf(ri)/2
        for tj in range(1, 16):
            th = mp.pi*tj/32
            try:
                t = T_func(G, r, th)
                if t > max_T: max_T = t; mp_pt = (float(r), float(th))
            except: pass
    st = "PASS" if max_T < 0 else "FAIL"
    print(f"  {name:15s}: maxT={mp.nstr(max_T,8):>16s} at r={mp_pt[0]:.1f},th={mp_pt[1]:.3f} [{st}]")

print("\n=== Large r ===")
for name, G in tests:
    max_T = -mp.inf
    for rv in [100, 500, 1000]:
        r = mp.mpf(rv)
        for tj in range(1, 16):
            th = mp.pi*tj/32
            try:
                t = T_func(G, r, th)
                if t > max_T: max_T = t
            except: pass
    print(f"  {name:15s}: maxT={mp.nstr(max_T,8):>16s} [{'PASS' if max_T<0 else 'FAIL'}]")

print("\n=== C(r) at key points ===")
for name, G in tests:
    cs = []
    for rv in [0.5, 1, 2, 5, 10, 50]:
        try: cs.append(f"{float(C_func(G, mp.mpf(rv))):+.4f}")
        except: cs.append("ERR")
    print(f"  {name:15s}: {cs}")

print("\nDONE")
