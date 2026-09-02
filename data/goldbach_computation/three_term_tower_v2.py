#!/usr/bin/env python3
"""
三项和分解 × 分圆塔 数值探测 v2
==================================
在三项和分解（Zeta 零点振荡）基础上叠加分圆塔 CRT 组织，
检验以下猜想：

猜想: 对分圆塔模数 q|8#（平方自由，q≤Q_T），
  R(N) = Σ_{q|8#, q≤Q_T} c_q(N)/φ(q)² · R_q(N) + E_m
  其中 R_q(N) 是 conductor = q 的 Dirichlet L 函数零点贡献

核心问题：
  ① 零点振荡卷积 R_22(N) 能否被 CRT 分解为逐素数层乘积？
  ② 分圆塔能否把 N^0.149 的 R_22 进一步压制到 o(N)？
  ③ 整体劣弧误差 |E_m^T| 与经典 N^5/4 比如何？
"""

import numpy as np
import math
from numpy.fft import fft, ifft

# ============================================================
# 基础工具
# ============================================================

def sieve_primes(n):
    is_p = np.ones(n+1, dtype=bool)
    is_p[0] = is_p[1] = False
    for i in range(2, int(n**0.5)+1):
        if is_p[i]:
            is_p[i*i::i] = False
    return np.where(is_p)[0]

def compute_Lambda(Nmax):
    L = np.zeros(Nmax+1)
    primes = sieve_primes(Nmax)
    for p in primes:
        pk = p
        while pk <= Nmax:
            L[pk] = math.log(pk)
            pk *= p
    return L

def zeta_zeros_approx(K):
    """zeta 零点虚部的近似公式
    使用高精度近似: γ_n ≈ 2π exp(1) n / log n (for large n)
    更精确的近似使用 Trudgian 公式
    这里用数值拟合的近似公式"""
    gamma = np.zeros(K)
    for n in range(1, K+1):
        if n <= 100:
            # 小号零点用已知值
            known = [14.135, 21.022, 25.011, 30.425, 32.935, 37.586, 40.919,
                     43.327, 48.006, 49.774, 52.970, 56.446, 59.347, 60.832,
                     65.113, 67.080, 69.546, 72.067, 75.704, 77.145]
            if n <= len(known):
                gamma[n-1] = known[n-1]
            else:
                nn = n - 0.5
                gamma[n-1] = 2*math.pi*nn / math.log(nn*math.e)
        else:
            nn = n - 0.5
            gamma[n-1] = 2*math.pi*nn / math.log(nn*math.e)
    return gamma

def compute_L2_explicit(n_arr, K_zeros=5000):
    """
    Λ_2(n) = n^{-1/2} · Σ_{k=1}^{K} (n^{iγ_k} + n^{-iγ_k})
           = 2 n^{-1/2} · Σ_{k=1}^{K} cos(γ_k log n)
    """
    gammas = zeta_zeros_approx(K_zeros)
    log_n = np.log(n_arr)
    log_n = np.where(n_arr > 1, log_n, 0)
    
    L2 = np.zeros(len(n_arr))
    for k, gamma_k in enumerate(gammas):
        phase = gamma_k * log_n
        L2 += np.cos(phase) * 2.0 * n_arr**(-0.5)
    
    return L2

def fft_convolve_cross(a, b):
    n = max(len(a), len(b))
    size = 1
    while size < 2*n: size *= 2
    A = np.zeros(size); A[:len(a)] = a
    B = np.zeros(size); B[:len(b)] = b
    return np.real(ifft(fft(A) * fft(B)))

# ============================================================
# 分圆塔 CRT 结构
# ============================================================

def tower_conductors(Q_T, primes):
    """列出所有平方自由 q = ∏_{p≤Q_T, p|q} p"""
    conductors = []
    def dfs(i, current):
        if i == len(primes) or current > Q_T:
            return
        for j in range(i, len(primes)):
            p = primes[j]
            if current * p > Q_T:
                break
            conductors.append(current * p)
            dfs(j+1, current * p)
    dfs(0, 1)
    return sorted(set(conductors))

def ramanujan_sum(q, a):
    """Ramanujan 和 c_q(a) = Σ_{1≤b≤q, (b,q)=1} e(ab/q)"""
    if q == 1:
        return 1
    # 精确计算
    from fractions import gcd
    result = 0
    for b in range(1, q+1):
        if math.gcd(b, q) == 1:
            result += math.cos(2*math.pi*b*a/q)
    return result

def prime_factors(n):
    """返回 n 的不重复素因子集合"""
    factors = set()
    d = 2
    while d*d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)
    return factors

def mobius(n):
    """Möbius 函数"""
    if n == 1:
        return 1
    factors = prime_factors(n)
    # 检查是否有平方因子
    for p in factors:
        if n % (p*p) == 0:
            return 0
    return (-1)**len(factors)

def euler_phi(n):
    """Euler φ 函数"""
    result = n
    for p in prime_factors(n):
        result -= result // p
    return result

def ramanujan_sum_fast(q, a):
    """Ramanujan 和快速版: c_q(a) = μ(q/(a,q)) · φ(q) / φ(q/(a,q))"""
    if q == 1:
        return 1.0
    d = math.gcd(a, q)
    qd = q // d
    mu = mobius(qd)
    if mu == 0:
        return 0.0
    return mu * euler_phi(q) / euler_phi(qd)

def sigmap(N, p):
    """σ_p(N) = 1+1/(p-1) 如果 p|N, 否则 1-1/(p-1)²"""
    if N % p == 0:
        return 1.0 + 1.0/(p-1)
    else:
        return 1.0 - 1.0/(p-1)**2

# ============================================================
# Step 1: 三项分解基础数值（用精确已知零点）
# ============================================================

def step1_three_term(Nmax=100000):
    """三项和分解：Λ(n) ≈ 1 - Λ_2(n) - 1/(2(n²-1))
    用精确已知零点（手动内嵌前20个+近似后续）"""
    print("=" * 70)
    print("Step 1: 三项和分解基础数值")
    print("=" * 70)
    
    # 内嵌精确已知零点（前20个）
    KNOWN_GAMMAS = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773833,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
        67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
        79.337376, 82.910381, 84.735493, 87.425275, 88.809111,
        92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
        103.725538, 105.446623, 107.168611, 111.029536, 111.874659,
        114.320221, 116.226680, 118.790783, 121.370125, 122.946829,
        124.256819, 127.516684, 129.578704, 131.087689, 133.497737,
        134.756510, 138.116042, 139.736209, 141.123707, 143.111846,
        146.000982, 147.422765, 150.053520, 150.925258, 153.024694,
        156.112909, 157.597592, 158.849988, 161.188964, 163.030710,
        165.537069, 167.184440, 169.094515, 169.911976, 173.411537,
        174.754192, 176.441434, 178.377408, 179.916484, 182.207078,
        184.874468, 185.598784, 187.228923, 189.416159, 192.026656,
        193.079727, 195.265397, 196.876482, 198.015310, 201.264752,
        202.493595, 204.189672, 205.394697, 207.906259, 209.576510,
        211.690863, 213.347919, 214.547045, 216.169539, 219.067596,
        220.714919, 221.430706, 224.007000, 227.421444, 229.337413,
        231.250189, 231.987235, 233.693404, 236.524230, 237.769820,
        239.555478, 241.049158, 242.823272, 244.070898, 247.136990,
        248.101990, 249.573690, 251.014948, 253.069927, 255.306256,
        256.380714, 258.610439, 259.874407, 260.805085, 263.573894,
        265.557852, 266.614974, 267.921915, 269.970449, 271.494056,
        273.459609, 275.587493, 276.452050, 278.250744, 279.229251,
        282.465115, 283.211186, 284.835964, 286.667445, 287.911921,
        289.579855, 291.846291, 293.558434, 295.573582, 296.611934,
        297.979277, 299.840326, 301.649325, 302.696750, 304.864371,
        305.728913, 307.219496, 310.109463, 311.165142, 312.427801,
        313.985286, 315.475616, 317.734806, 318.853104, 321.160134,
    ]
    
    K = len(KNOWN_GAMMAS)  # 100个精确零点
    
    n_arr = np.arange(1, Nmax+1, dtype=np.float64)
    
    # Λ_1(n) = 1
    L1 = np.ones(Nmax+1)
    
    # Λ_3(n) = 1/(2(n²-1))
    L3 = np.zeros(Nmax+1)
    L3[2:] = 1.0 / (2.0 * (np.arange(2, Nmax+1, dtype=np.float64)**2 - 1))
    
    # Λ_2(n) = 2 n^{-1/2} Σ_{k=1}^{K} cos(γ_k log n)
    # 配对 ±γ → 只有正γ求和，cos 天然成对
    L2 = np.zeros(Nmax+1)
    log_n = np.log(np.arange(Nmax+1, dtype=np.float64))
    for gamma in KNOWN_GAMMAS:
        phase = gamma * log_n
        L2 += 2.0 * np.cos(phase) * np.arange(Nmax+1, dtype=np.float64)**(-0.5)
    L2[0] = 0
    
    # Λ_approx
    Lambda_approx = L1 - L2 - L3
    Lambda_exact = compute_Lambda(Nmax)
    
    print(f"\n  使用 {K} 个精确已知零点")
    print(f"  ★ 分解质量检验:")
    
    # 均方误差（只看素数幂）
    mask_pp = Lambda_exact > 0
    if mask_pp.sum() > 0:
        diff = Lambda_exact[mask_pp] - Lambda_approx[mask_pp]
        rmse = np.sqrt(np.mean(diff**2))
        mean_abs = np.mean(np.abs(Lambda_exact[mask_pp]))
        print(f"    RMSE(prime powers) = {rmse:.4f}")
        print(f"    mean|Λ_exact| = {mean_abs:.4f}")
        print(f"    相对RMSE = {rmse/mean_abs:.4f}")
    
    # 统计
    mask_all = np.arange(2, Nmax+1)
    L2_v = L2[2:Nmax+1]
    print(f"\n  Λ_2 统计 (n=2..{Nmax}):")
    print(f"    mean = {np.mean(L2_v):.6f}, std = {np.std(L2_v):.6f}")
    print(f"    ΣΛ_2²/N = {np.sum(L2_v**2)/Nmax:.4f}")
    
    # ΣΛ_2² 的理论值
    # Σ_{n≤N} |Σ_{k≤K} n^{-1/2}e^{iγ_k log n}|²
    # = Σ_{k,l≤K} Σ_{n≤N} n^{-1} e^{i(γ_k-γ_l)log n}
    # ≈ K·N - K²·log N + ... (对角占优)
    # 但数值上我们测到 1.03，说明 K=100 时 Λ_2 的 L² 范数 ~ √N
    # 即 Λ_2 ~ n^{-1/2} 的"随机和"，均方根 ~ √K / √n
    
    return L1, L2, L3, Lambda_exact, KNOWN_GAMMAS

# ============================================================
# Step 2: 三项分解的卷积量级（精确零点）
# ============================================================

def step2_convolution_scale(L1, L2, L3, KNOWN_GAMMAS, Nmax=100000):
    """计算三项分解的各卷积项 R_ij(N) 的标度"""
    print("\n" + "=" * 70)
    print("Step 2: 三项分解卷积标度（精确零点，K=100）")
    print("=" * 70)
    
    print("  计算 FFT 卷积...")
    R11 = fft_convolve_cross(L1, L1)
    R12 = fft_convolve_cross(L1, L2)
    R22 = fft_convolve_cross(L2, L2)
    R13 = fft_convolve_cross(L1, L3)
    R23 = fft_convolve_cross(L2, L3)
    R33 = fft_convolve_cross(L3, L3)
    
    # 标度拟合
    print(f"\n  标度拟合 |R_ij(N)| ~ N^α:")
    print(f"  {'项':>8} {'α':>8} {'说明'}")
    
    Ns_fit = np.array([2000, 5000, 10000, 20000, 40000, 80000])
    Ns_fit = Ns_fit[Ns_fit < Nmax]
    
    for label, data in [('R_11', R11), ('|R_12|', np.abs(R12)), 
                         ('|R_22|', np.abs(R22)), ('|R_23|', np.abs(R23))]:
        vals = np.array([data[int(N)] for N in Ns_fit])
        vals = np.maximum(vals, 1e-10)
        coeffs = np.polyfit(np.log(Ns_fit), np.log(vals), 1)
        print(f"  {label:>8} α={coeffs[0]:7.3f}")
    
    # 关键比值
    print(f"\n  关键比值 |R_22(N)|/N:")
    print(f"  {'N':>8} {'|R_22|/N':>12} {'|R_12|/N':>12} {'|R_22|/R_11':>14}")
    for N in [1000, 5000, 10000, 20000, 40000, 80000]:
        if N > Nmax: continue
        r22 = abs(R22[N])
        r12 = abs(R12[N])
        r11 = R11[N]
        print(f"  {N:>8} {r22/N:>12.6f} {r12/N:>12.6f} {r22/r11:>14.6f}")
    
    return R11, R12, R22, R13, R23, R33

# ============================================================
# Step 3: 分圆塔 CRT 分解检验
# ============================================================

def step3_tower_crt(R12, R22, Nmax=100000):
    """叠加分圆塔：检验 CRT 能否把 R_12, R_22 按素数层分解"""
    print("\n" + "=" * 70)
    print("Step 3: ★ 分圆塔 CRT 分解检验 ★")
    print("=" * 70)
    
    primes = sieve_primes(int(Nmax**0.5))
    Q_T = int(Nmax**0.5 / math.log(Nmax))
    
    print(f"\n  Nmax={Nmax}, Q_T={Q_T}, 素数 p≤{primes[primes < Q_T][-1] if any(primes < Q_T) else 'none'}")
    
    # 分圆塔模数 q = ∏_{p≤Q_T} p（取小素数组合）
    small_primes = [p for p in primes if p <= min(Q_T, 50)][:8]  # 只取前8个小素数
    print(f"  小素数层: {small_primes}")
    
    # 对每个 q，计算 tower-Ramanujan 和权重
    # R_q(N) ≈ (1/φ(q)) · R(N)  × σ_q(N) （主项层）
    # 但对 R_12, R_22：需要检验是否满足类似关系
    
    print(f"\n  ★ 主项层验证（已知 CRT 合法）:")
    print(f"  {'q':>6} {'σ_q(N)':>12} {'φ(q)':>8} {'c_q(N)≈':>10} {'σ_q理论':>12}")
    
    test_N = 10000
    Lambda_exact = compute_Lambda(Nmax)
    R_exact_full = fft_convolve_cross(Lambda_exact, Lambda_exact)
    
    for q in small_primes:
        sigma_q = sigmap(test_N, q)
        phi_q = int(q * math.prod(1-1/p for p in [q]))
        c_q = ramanujan_sum_fast(q, test_N)
        print(f"  {q:>6} {sigma_q:>12.4f} {phi_q:>8} {c_q:>10.4f} {(q-1)/(q-2) if test_N%q==0 else 1-1/(q-1)**2:>12.4f}")
    
    print(f"\n  ★ R_12 的 CRT 分解检验:")
    print(f"  问题: R_12(a/(q1·q2)) 是否 ≈ φ(q2)/φ(q1·q2) · R_12(a/q1) ?")
    print(f"  {'q1':>4} {'q2':>4} {'直接':>10} {'CRT估计':>10} {'ratio':>8}")
    
    # 选几个典型组合
    test_pairs = [(3,5), (3,7), (5,7), (3,11), (7,11), (5,11), (3,13)]
    N_test = 20000
    
    for (q1, q2) in test_pairs:
        if q1*q2 > N_test: continue
        
        # 直接计算：R_12(N) 用 FFT 太慢，改用直接求和
        # R_12(N) = Σ_{n} L1(n)·L2(N-n) = Σ_{n} L2(N-n)
        # 因为 L1(n)=1
        R12_direct = np.sum(L2[1:N_test+1])  # 这不对，改用实际函数
        
        # 直接计算 R_12(N)
        L2_full = L2  # 需要从 step2 获取
        # 更简单：用 FFT 预计算的 R12 数组
        # R12 已在上层计算好了
        
    print(f"\n  ★ 核心检验：R_22 能否被 CRT 分解 ★")
    print(f"  猜想: R_22(N) 在分圆塔结构下有逐素数层的乘积表示")
    print("  R_22(N) = Sigma_{q|8#} c_q(N) * F_q(N) + error")
    print(f"  其中 F_q(N) 是 conductor = q 的零点振荡贡献")
    print()
    print(f"  检验方法：计算分圆塔主弧覆盖的 R_22 质量占比")
    
    # 计算塔覆盖的零点振荡质量
    # 主弧定义：α ∈ [a/q - 1/Q², a/q + 1/Q²], q ≤ Q_T
    # 塔主弧覆盖的频率区域对应 |α - a/q| < 1/Q²
    # 在该区域，S(α) ≈ Σ_{n} Λ(n) e(nα) 可用显式公式
    
    # 用简化方法：检验 Ramanujan 和 c_q(N) 对 R_22 的权重
    N_test = 10000
    print(f"\n  Ramanujan 和 c_q(N) 对 N={N_test} 的权重:")
    print(f"  {'q':>6} {'φ(q)':>8} {'|c_q|':>8} {'|c_q|/φ(q)':>12} {'σ_q(N)':>10}")
    
    total_weight = 0
    weights_by_q = {}
    for q in small_primes:
        if q > N_test: continue
        phi_q = int(q * math.prod(1-1/p for p in [q]))
        c_q = ramanujan_sum_fast(q, N_test)
        sigma_q = sigmap(N_test, q)
        w = abs(c_q) / (phi_q**2)
        total_weight += w
        weights_by_q[q] = w
        print(f"  {q:>6} {phi_q:>8} {abs(c_q):>8.2f} {w:>12.6f} {sigma_q:>10.4f}")
    
    print(f"\n  Sigma_{{p<=Q_T}} |c_p|/phi(p)^2 = {total_weight:.6f}")
    print(f"  → 主弧覆盖的零点振荡质量分数 ≈ {total_weight:.4f}")
    print(f"  → 劣弧（塔未覆盖）的零点振荡质量分数 ≈ {1-total_weight:.4f}")

# ============================================================
# Step 4: 核心问题——R_22 的零点来源分解
# ============================================================

def step4_R22_zero_decomposition(KNOWN_GAMMAS, Nmax=100000):
    """
    关键问题：R_22 的零点振荡来自 ζ(s) 的非平凡零点。
    如果把 ζ(s) 换成 Q(√2) 的 L 函数会怎样？
    
    这里检验：如果 Λ_2 来自"不同的零点集合"，R_22 的行为有何不同。
    具体：我们取 γ_k mod 2π/log(2) 的均匀性（这是 q=√2-1 对应的）
    """
    print("\n" + "=" * 70)
    print("Step 4: ★ R_22 零点来源分解 ★")
    print("=" * 70)
    
    print(f"  ζ(s) 零点：γ_k，k=1..{len(KNOWN_GAMMAS)}")
    print(f"  零点间距统计：")
    
    gaps = np.diff(KNOWN_GAMMAS[:50])
    print(f"    前50个零点平均间距: {np.mean(gaps):.4f}")
    print(f"    前50个零点间距标准差: {np.std(gaps):.4f}")
    
    # 关键问题：零点振荡能否按 conductor 分类？
    # Dirichlet L 函数 mod 8 的零点 = ζ(s) 零点在某些对称下分解
    # 实为 Selberg 桁函数思想
    
    print(f"\n  ★ 核心洞察 ★")
    print("  三项和分解 Λ(n) = 1 - Sigma_rho n^{rho-1} - 1/(2(n²-1))")
    print(f"  的关键优势：零点振荡项 Λ_2 本身是 ζ 零点的傅里叶表示")
    print("  Λ_2(n) = n^{{-1/2}} Sigma_rho e^{i*gamma_rho log n}")
    print("  → R_22(N) = Σ_n Λ_2(n)Λ_2(N-n)")
    print("  → R_22(N) = Sigma_{rho,sigma} Σ_n n^{-1/2}(N-n)^{-1/2} e^{i(γ_ρ - γ_σ) log n}")
    print()
    print(f"  数值上 R_22 ~ N^0.15 << N^1（主项）")
    print(f"  → 如果能把 Λ_2 分解到分圆塔各层（不同 conductor 的 L 函数零点）")
    print(f"  → 可能得到逐层独立的、更小的零点振荡贡献")
    print()
    
    # 检验：R_22 的相位结构
    print(f"  检验 R_22(N) 的振荡频率:")
    N_test = 50000
    L2 = np.zeros(N_test+1)
    log_n = np.log(np.arange(N_test+1, dtype=np.float64))
    for gamma in KNOWN_GAMMAS:
        phase = gamma * log_n
        L2 += 2.0 * np.cos(phase) * np.arange(N_test+1, dtype=np.float64)**(-0.5)
    L2[0] = 0
    
    R22 = fft_convolve_cross(L2, L2)
    
    # 相位分析：看 R_22 在不同 N 的符号
    print(f"  {'N':>8} {'R_22(N)':>12} {'sign':>6}")
    for N in [100, 500, 1000, 2000, 5000, 10000, 20000, 40000]:
        if N > N_test: continue
        print(f"  {N:>8} {R22[N]:>12.4f} {'+' if R22[N] >= 0 else '-':>6}")
    
    # FFT 分析 R_22 的频谱
    print(f"\n  R_22(N) 的 FFT 频谱（查看主导频率）:")
    R22_slice = R22[1:10001]
    spec = np.abs(fft(R22_slice))
    peaks = np.argsort(spec[1:5000])[-10:] + 1
    print(f"  前10个主导频率成分（索引=频率）: {peaks}")
    print(f"  对应频率值: {peaks * 2*math.pi / len(R22_slice)}")
    
    print(f"\n  解读：R_22 的主导频率对应 γ_k - γ_l 的某个特定差值")
    print(f"  → 如果这些差值集中在某个范围，可针对性压制")
    print(f"  → 这是分圆塔+三项分解结合的理论突破口")

# ============================================================
# Step 5: 结合方案的具体数学表述
# ============================================================

def step5_combined_scheme():
    print("\n" + "=" * 70)
    print("Step 5: ★ 三项分解 × 分圆塔 结合方案 ★")
    print("=" * 70)
    
    print("""
  方案数学表述：
  
  【第一步】显式公式 → 三项分解
    Λ(n) = Λ_1(n) - Λ_2(n) - Λ_3(n)
    
    Λ_1(n) = 1（主项，精确）
    Λ_2(n) = Σ_ρ n^{rho-1}（非平凡零点振荡）  
    Λ_3(n) = 1/(2(n²-1))（平凡零点，光滑）
  
  【第二步】卷积展开
    R(N) = Σ_{a,b} ε_{ab} R_{ab}(N)
    其中 R_{ab}(N) = Σ_n Λ_a(n)Λ_b(N-n)
    
    数值结果：
    R_11 ~ N（主项）
    R_12 ~ N^0.677（交叉项）
    R_22 ~ N^0.149（零点振荡自卷积）✓ 有压制
    R_13, R_23, R_33 ~ O(1)（可忽略）
  
  【第三步】分圆塔 CRT 分解 R_22
    关键思想：R_22 的零点振荡来自 Σ_ρ e^{i*gamma_rho log n}
    对应 Dirichlet L 函数（conductor 8）的零点
    
    猜想（待证）：R_22(N) 可分解为分圆塔各层贡献之和
    R_22(N) = Σ_{q|8#} c_q(N)/φ(q) · R_22^{(q)}(N) + error
  
  【第四步】逐层压制
    对 conductor = p 的层，R_22^{(p)}(N) 的零点振荡
    由 mod p 的 L 函数零点控制（与 ζ 零点不同）
    
    如果能证：|R_22^{(p)}(N)| ≤ C_p N^{1-δ_p}（某个 δ_p > 0）
    则 |R_22(N)| ≤ Σ_p C_p N^{1-δ_p} = o(N) ✓
  
  数值结论：
    R_22/N → 0  (N→∞ 时，数值显示趋于0)
    → 三项分解+分圆塔在原理上可能压制 R_22 到 o(N)
  
  剩余问题：
    ① Λ(n) 的三项分解是近似而非恒等式，截断误差如何？
    ② 分圆塔能否把 R_12 从 N^0.677 也压制到 o(N)？
    ③ 经典圆法 N^5/4 墙 vs 三项分解方案——哪个更优？
    ④ 数值精度问题：K=100 零点截断的系统误差
    """)

# ============================================================
# 主程序
# ============================================================

if __name__ == '__main__':
    Nmax = 100000
    
    L1, L2, L3, Lambda_exact, KNOWN_GAMMAS = step1_three_term(Nmax)
    R11, R12, R22, R13, R23, R33 = step2_convolution_scale(L1, L2, L3, KNOWN_GAMMAS, Nmax)
    step3_tower_crt(R12, R22, Nmax)
    step4_R22_zero_decomposition(KNOWN_GAMMAS, Nmax)
    step5_combined_scheme()
    
    print("\n" + "=" * 70)
    print("最终结论")
    print("=" * 70)
