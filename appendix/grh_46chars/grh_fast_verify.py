#!/usr/bin/env python3
"""Fast verification of angular monotonicity for general Dirichlet L-functions."""
import mpmath as mp
mp.mp.dps = 30

def make_completed_L(chars, q, a):
    def L_func(s): return mp.dirichlet(s, chars)
    def Lambda(s): return (q/mp.pi)**((s+a)/2) * mp.gamma((s+a)/2) * L_func(s)
    return Lambda

def T_func(Func, r, theta):
    s = mp.mpf('0.5') + r * mp.e**(1j*theta)
    Fp = mp.diff(lambda z: mp.log(Func(z)), s, 1)
    return mp.re(1j * r * mp.e**(1j*theta) * Fp)

def C_func(Func, r):
    s0 = mp.mpf('0.5') + r
    Fp = mp.diff(lambda z: mp.log(Func(z)), s0, 1)
    Fpp = mp.diff(lambda z: mp.log(Func(z)), s0, 2)
    return mp.re(r * Fp + r*r * Fpp)

# Characters
chi3 = [0, 1, -1]; Lam3 = make_completed_L(chi3, 3, 1); G3 = lambda s: Lam3(s)**2
chi4 = [0, 1, 0, -1]; Lam4 = make_completed_L(chi4, 4, 1); G4 = lambda s: Lam4(s)**2
chi5q = [0, 1, -1, -1, 1]; Lam5q = make_completed_L(chi5q, 5, 0); G5q = lambda s: Lam5q(s)**2
chi5c = [0, 1, 1j, -1j, -1]; chi5cb = [0, 1, -1j, 1j, -1]
Lam5c = make_completed_L(chi5c, 5, 1); Lam5cb = make_completed_L(chi5cb, 5, 1)
G5c = lambda s: Lam5c(s)*Lam5cb(s)
chi7 = [0, 1, 1, -1, 1, -1, -1]; Lam7 = make_completed_L(chi7, 7, 1); G7 = lambda s: Lam7(s)**2
chi8 = [0, 1, 0, -1, 0, -1, 0, 1]; Lam8 = make_completed_L(chi8, 8, 0); G8 = lambda s: Lam8(s)**2
chi8b = [0, 1, 0, 1, 0, -1, 0, -1]; Lam8b = make_completed_L(chi8b, 8, 1); G8b = lambda s: Lam8b(s)**2

tests = [
    ("q=3 odd real", G3, 3), ("q=4 beta", G4, 4),
    ("q=5 even real", G5q, 5), ("q=5 odd complex", G5c, 5),
    ("q=7 odd real", G7, 7), ("q=8 even real", G8, 8),
    ("q=8 odd real", G8b, 8),
]

# Test 1: Medium grid r=0.5..50
print("=== Angular monotonicity (r=0.5..50, 32 theta) ===")
for name, G, q in tests:
    max_T = -mp.inf; max_pt = None; cnt = 0
    for ri in range(1, 101):
        r = mp.mpf(ri)/2
        for tj in range(1, 32):
            th = mp.pi*tj/64
            try:
                t = T_func(G, r, th); cnt += 1
                if t > max_T: max_T = t; max_pt = (float(r), float(th))
            except: pass
    st = "PASS" if max_T < 0 else "FAIL"
    print(f"  {name}: {cnt} pts, maxT={mp.nstr(max_T,8)} at r={max_pt[0]:.1f},th={max_pt[1]:.3f} [{st}]")

# Test 2: Large r
print("\n=== Large r (50,100,500,1000) ===")
for name, G, q in tests:
    max_T = -mp.inf; max_pt = None
    for rv in [50, 100, 500, 1000]:
        r = mp.mpf(rv)
        for tj in range(1, 64):
            th = mp.pi*tj/128
            try:
                t = T_func(G, r, th)
                if t > max_T: max_T = t; max_pt = (rv, float(th))
            except: pass
    st = "PASS" if max_T < 0 else "FAIL"
    print(f"  {name}: maxT={mp.nstr(max_T,8)} at r={max_pt[0]},th={max_pt[1]:.3f} [{st}]")

# Test 3: C(r) > 0
print("\n=== C(r) positivity ===")
for name, G, q in tests:
    cvs = []
    for rv in [0.5, 1, 2, 5, 10, 50]:
        try: cvs.append(mp.nstr(C_func(G, mp.mpf(rv)), 6))
        except: cvs.append("ERR")
    print(f"  {name}: {cvs}")

# Test 4: Functional equation
print("\n=== G(s)=G(1-s) check ===")
for name, G, q in tests:
    s = mp.mpf('0.3')+3j
    r = abs(G(s)/G(1-s)-1)
    print(f"  {name}: |G(s)/G(1-s)-1|={mp.nstr(r,4)}")

# Test 5: Very large q (q=101, a primitive character)
print("\n=== Large conductor q=101 (Legendre symbol, odd) ===")
chi101 = [0] + [1 if pow(n, 50, 101) == 1 else -1 for n in range(1, 101)]
Lam101 = make_completed_L(chi101, 101, 1)
G101 = lambda s: Lam101(s)**2
max_T = -mp.inf
for rv in [0.5, 1, 2, 5, 10, 50, 100]:
    r = mp.mpf(rv)
    for tj in range(1, 32):
        th = mp.pi*tj/64
        try:
            t = T_func(G101, r, th)
            if t > max_T: max_T = t
        except: pass
print(f"  q=101: maxT={mp.nstr(max_T,8)} [{'PASS' if max_T < 0 else 'FAIL'}]")
for rv in [0.5, 1, 2, 5]:
    print(f"  C({rv})={mp.nstr(C_func(G101, mp.mpf(rv)), 6)}")

print("\nDONE")
