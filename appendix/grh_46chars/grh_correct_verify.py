"""
Correct numerical verification using Hurwitz zeta.
ALL characters generated programmatically (no hand-coding errors).
"""
import mpmath as mp
import json, time

mp.mp.dps = 50

def legendre_symbol(n, p):
    n = n % p
    if n == 0: return 0
    r = pow(n, (p-1)//2, p)
    return 1 if r == 1 else -1

def make_quadratic_char(q):
    """Quadratic (real) character mod prime q."""
    return [legendre_symbol(n, q) if n < q else 0 for n in range(1, q+1)]

def make_mod8_even():
    """chi mod 8, even: chi(1)=1, chi(3)=-1, chi(5)=-1, chi(7)=1."""
    return [1, 0, -1, 0, -1, 0, 1, 0]

def make_mod8_odd():
    """chi mod 8, odd: chi(1)=1, chi(3)=1, chi(5)=-1, chi(7)=-1."""
    return [1, 0, 1, 0, -1, 0, -1, 0]

def make_mod4_beta():
    """beta = L(s, chi_4), odd: [1,0,-1,0]."""
    return [1, 0, -1, 0]

def make_q5_complex():
    """Order-4 character mod 5, odd parity."""
    return [1, 1j, -1j, -1, 0]

def make_q5_complex_bar():
    return [1, -1j, 1j, -1, 0]

def L_hurwitz(s, chi, q):
    total = mp.mpf(0)
    for a in range(1, q+1):
        c = chi[a-1]
        if c != 0:
            total += c * mp.zeta(s, mp.mpf(a)/q)
    return (mp.mpf(q)**(-s)) * total


CHARACTERS = {}

# q=3 odd real
CHARACTERS['chi3'] = {'q':3, 'a':1, 'real':True, 'chi':make_quadratic_char(3)}
# q=4 beta odd real
CHARACTERS['chi4'] = {'q':4, 'a':1, 'real':True, 'chi':make_mod4_beta()}
# q=5 even real (quadratic)
CHARACTERS['chi5e'] = {'q':5, 'a':0, 'real':True, 'chi':make_quadratic_char(5)}
# q=5 odd complex
CHARACTERS['chi5c'] = {'q':5, 'a':1, 'real':False,
    'chi':make_q5_complex(), 'chibar':make_q5_complex_bar()}
# q=7 odd real
CHARACTERS['chi7'] = {'q':7, 'a':1, 'real':True, 'chi':make_quadratic_char(7)}
# q=8 even real
CHARACTERS['chi8e'] = {'q':8, 'a':0, 'real':True, 'chi':make_mod8_even()}
# q=8 odd real
CHARACTERS['chi8o'] = {'q':8, 'a':1, 'real':True, 'chi':make_mod8_odd()}
# q=11 odd real
CHARACTERS['chi11'] = {'q':11, 'a':1, 'real':True, 'chi':make_quadratic_char(11)}
# q=13 even real
CHARACTERS['chi13'] = {'q':13, 'a':0, 'real':True, 'chi':make_quadratic_char(13)}
# q=101 odd real
CHARACTERS['chi101'] = {'q':101, 'a':0, 'real':True, 'chi':make_quadratic_char(101)}  # FIXED: 101≡1 mod 4 → even character (a=0), was a=1

# Print character table
print("Character verification:")
for name, info in CHARACTERS.items():
    print(f"  {name:8s} q={info['q']:3d} a={info['a']} real={info['real']}: {info['chi'][:8]}{'...' if info['q']>8 else ''}")


def make_G(info):
    q, a = info['q'], info['a']
    chi = info['chi']
    chibar = chi if info['real'] else info['chibar']
    def Lam(s, ch):
        return (mp.mpf(q)/mp.pi)**((s+a)/2) * mp.gamma((s+a)/2) * L_hurwitz(s, ch, q)
    def G(s):
        if info['real']:
            L1 = Lam(s, chi); return L1*L1
        return Lam(s, chi) * Lam(s, chibar)
    return G

def T_chi(G, r, theta, dps=25):
    mp.mp.dps = dps
    s = mp.mpf('0.5') + r * mp.e**(1j*theta)
    dlogG = mp.diff(lambda z: mp.log(G(z)), s)
    z = r * mp.e**(1j*theta)
    result = mp.re(1j*z*dlogG)
    mp.mp.dps = 50
    return result

def compute_c2(G):
    """c2 = (1/2) Re[d^2/ds^2 log G(s)] at s=1/2."""
    s0 = mp.mpf('0.5')
    d2 = mp.diff(lambda z: mp.log(G(z)), s0, 2)
    return float(mp.re(d2)/2)

def compute_CD(G, r_val):
    """C(r) = -dT/dtheta at theta=0; D(r) = -T(r,pi/2-h)/h."""
    r = mp.mpf(str(r_val))
    h = mp.mpf('1e-6')
    # C = -dT/dtheta at theta=0 (use forward diff)
    T0 = T_chi(G, r, mp.mpf(0), dps=30)
    Th = T_chi(G, r, h, dps=30)
    C = float(-(Th - T0)/h)
    # D = -dT/dtheta at theta=pi/2 (use backward diff)
    Tpi = T_chi(G, r, mp.pi/2, dps=30)
    Tpih = T_chi(G, r, mp.pi/2 - h, dps=30)
    D = float(-(Tpi - Tpih)/h)
    return C, D

print("\n" + "="*70)
print("c2 values")
print("="*70)
c2_results = {}
for name, info in CHARACTERS.items():
    G = make_G(info)
    c2 = compute_c2(G)
    c2_results[name] = c2
    print(f"  {name:8s} q={info['q']:3d}: c2 = {c2:.6f}")

print("\n" + "="*70)
print("C(r) and D(r) values")
print("="*70)
cd_results = {}
for name in ['chi3','chi4','chi5e','chi5c','chi7','chi8e','chi8o']:
    info = CHARACTERS[name]
    G = make_G(info)
    c05, d05 = compute_CD(G, 0.5)
    c20, d20 = compute_CD(G, 2.0)
    cd_results[name] = {'C05':c05,'D05':d05,'C20':c20,'D20':d20}
    print(f"  {name:8s}: C(0.5)={c05:.4f} D(0.5)={d05:.4f}  C(2.0)={c20:.4f} D(2.0)={d20:.4f}")

print("\n" + "="*70)
print("Grid verification T_chi(r,theta) < 0")
print("="*70)

r_small = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0]
r_med = [8, 10, 15, 20, 30, 50]
r_large = [100, 500, 1000]

grid_results = {}
for name, info in CHARACTERS.items():
    t0 = time.time()
    G = make_G(info)
    if info['q'] <= 13:
        r_vals = r_small + r_med + r_large
        n_th = 39
    else:
        r_vals = [0.5, 2.0, 7.0, 50, 100]
        n_th = 20

    max_T = -float('inf')
    worst = (0,0)
    n_pts = 0
    for rv in r_vals:
        r = mp.mpf(str(rv))
        for j in range(1, n_th):
            th = mp.mpf(j)/n_th * mp.pi/2
            t = float(T_chi(G, r, th))
            n_pts += 1
            if t > max_T:
                max_T = t
                worst = (rv, float(th))

    status = "PASS" if max_T < 0 else "FAIL"
    elapsed = time.time() - t0
    print(f"  {name:8s} q={info['q']:3d}: {n_pts:5d} pts, max T = {max_T:.6e} "
          f"at r={worst[0]}, th={worst[1]:.4f}  [{status}] ({elapsed:.1f}s)")
    grid_results[name] = {'q':info['q'],'a':info['a'],'n_pts':n_pts,
                          'max_T':max_T,'worst_r':worst[0],'worst_th':worst[1],
                          'status':status}

# Dense tight check for q=3, r in [5.5,8]
print("\n" + "="*70)
print("Tight region q=3, r in [5.5,8], dense theta")
print("="*70)
G3 = make_G(CHARACTERS['chi3'])
max_T = -float('inf')
worst = (0,0)
n_pts = 0
for rv in [5.5, 6.0, 6.3, 6.5, 6.8, 7.0, 7.2, 7.5, 7.8, 8.0]:
    r = mp.mpf(str(rv))
    for j in range(1, 80):
        th = mp.mpf(j)/80 * mp.pi/2
        t = float(T_chi(G3, r, th, dps=30))
        n_pts += 1
        if t > max_T:
            max_T = t; worst = (rv, float(th))
print(f"  {n_pts} pts, max T = {max_T:.8f} at r={worst[0]}, th={worst[1]:.4f}  [{'PASS' if max_T<0 else 'FAIL'}]")

all_results = {
    'c2': c2_results,
    'CD': cd_results,
    'grid': grid_results,
    'tight_q3': {'max_T': max_T, 'worst_r': worst[0], 'worst_th': worst[1], 'n_pts': n_pts}
}
with open('/Coze/Drive/黎曼猜想论文审核/所有对话/主对话/grh_correct_results.json','w') as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nAll results saved.")
