#!/usr/bin/env python3
"""
三项和分解 × 分圆塔 数值探测 v1
=================================
Step 1: Λ(n) = Λ_1(n) - Λ_2(n) - Λ_3(n)
  Λ_1(n) = 1        (主项，来自 zeta 极点)
  Λ_2(n) = Σ_ρ n^{ρ-1}  (非平凡零点振荡)
  Λ_3(n) = 1/(2(n²-1))  (平凡零点光滑修正)

Step 2: R(N) = Σ_{a,b,c ∈ {1,2,3}} ± R_{ab}(N)
  其中 R_{ab}(N) = Σ_{n≤N} Λ_a(n) Λ_b(N-n)

Step 3: 看每项量级，哪些可压、哪些是墙

Step 4: 叠加分圆塔 CRT 结构
"""

import numpy as np
import math
from collections import defaultdict

# ============================================================
# 基础工具
# ============================================================

def sieve_primes(n):
    """埃氏筛"""
    is_p = np.ones(n+1, dtype=bool)
    is_p[0] = is_p[1] = False
    for i in range(2, int(n**0.5)+1):
        if is_p[i]:
            is_p[i*i::i] = False
    return np.where(is_p)[0]

def compute_Lambda(Nmax):
    """经典 von Mangoldt 函数"""
    L = np.zeros(Nmax+1)
    primes = sieve_primes(Nmax)
    for p in primes:
        pk = p
        while pk <= Nmax:
            L[pk] = math.log(pk)
            pk *= p
    return L

def zeta_zeros(K=100000):
    """取前 K 个 Riemann zeta 零点虚部 (Hardy 列表 / 近似)
    用近似公式: γ_n ≈ 2π n / W(n) 其中 W 修正
    但为精确起见，直接内嵌前几千个已知零点"""
    # 使用近似: Gram 点附近的零点
    # γ_n ~ 2π(n - 11/8) / W(2π(n-11/8))
    # 更好的近似: γ_n ≈ 2π*n / log(n) 对大 n
    # 最实用: 用 Riemann-Siegel theta 函数的反函数
    # 但为了数值实验，用经典近似公式足够
    gamma = np.zeros(K)
    for n in range(1, K+1):
        # 近似公式 (Trudgian 2011 简化版)
        # γ_n ≈ 2π(n - 11/8) / log(2π(n - 11/8))  (一阶近似)
        nn = n - 11.0/8.0
        if nn > 0:
            gamma[n-1] = 2*math.pi*nn / math.log(2*math.pi*nn)
        else:
            gamma[n-1] = 14.137  # γ_1
    return gamma

def compute_Lambda_decomp(Nmax, K_zeros=5000):
    """
    三项和分解:
    Λ(n) ≈ Λ_1(n) - Λ_2(n) - Λ_3(n)
    
    Λ_1(n) = 1
    Λ_2(n) = Σ_{k=1}^{K} n^{iγ_k - 1/2} = n^{-1/2} Σ_k n^{iγ_k}
             (GRH 下 ρ = 1/2 + iγ_k)
    Λ_3(n) = 1/(2(n²-1))
    
    注意: 实际 Λ(n) 是离散算术函数，三项分解是连续近似。
    差异 Λ(n) - (Λ_1 - Λ_2 - Λ_3) 是截断误差。
    """
    gammas = zeta_zeros(K_zeros)
    
    L1 = np.ones(Nmax+1)  # 主项
    L3 = np.zeros(Nmax+1)  # 平凡零点
    L3[2:] = 1.0 / (2.0 * (np.arange(2, Nmax+1, dtype=np.float64)**2 - 1))
    
    # Λ_2: 零点振荡项
    # Λ_2(n) = n^{-1/2} * Σ_k cos(γ_k * log n)  (取实部, 因为 ρ+ρ̄ 配对)
    # 实际上 ρ = 1/2+iγ, n^{ρ-1} = n^{-1/2+iγ} = n^{-1/2} * e^{iγ log n}
    # 配对 ρ,ρ̄ → 2 Re = 2 cos(γ log n)
    # 但原始公式是 Σ_ρ (不加配对), 所以直接取实部
    L2 = np.zeros(Nmax+1)
    for n in range(2, Nmax+1):
        logn = math.log(n)
        # Σ_k n^{iγ_k} = Σ_k cos(γ_k log n) + i Σ_k sin(γ_k log n)
        # 因为 γ, -γ 配对(或 ρ, ρ̄ 配对), sin 项相消
        # 所以 Λ_2(n) = n^{-1/2} * 2 * Σ_k cos(γ_k * log n) (因子2来自配对)
        # 但原始求和是 over ALL zeros (both +γ and -γ), 所以不需要因子2
        # 我们只取 γ>0 的零点, 所以确实需要因子 2
        phase = gammas * logn
        L2[n] = n**(-0.5) * 2.0 * np.sum(np.cos(phase))
    
    return L1, L2, L3

# ============================================================
# 卷积计算
# ============================================================

def fft_convolve_cross(a, b):
    """用 FFT 计算 Σ_{n} a(n) * b(N-n) for all N simultaneously"""
    n = max(len(a), len(b))
    size = 1
    while size < 2*n: size *= 2
    A = np.zeros(size); A[:len(a)] = a
    B = np.zeros(size); B[:len(b)] = b
    return np.real(ifft(fft(A) * fft(B)))

from numpy.fft import fft, ifft

# ============================================================
# 主程序
# ============================================================

def main():
    Nmax = 50000
    K_zeros = 3000  # 零点截断数
    
    print("=" * 70)
    print("三项和分解 × 分圆塔 数值探测 v1")
    print(f"Nmax = {Nmax}, K_zeros = {K_zeros}")
    print("=" * 70)
    
    # ---- Step 0: 经典 Λ(n) ----
    print("\n[0] 计算经典 Λ(n)...")
    Lambda_exact = compute_Lambda(Nmax)
    
    # ---- Step 1: 三项分解 ----
    print(f"\n[1] 三项和分解 Λ(n) ≈ Λ_1 - Λ_2 - Λ_3 ...")
    print(f"    Λ_1(n) = 1 (主项)")
    print(f"    Λ_2(n) = n^{{-1/2}} · 2·Σ_k cos(γ_k·log n) (零点振荡, K={K_zeros})")
    print(f"    Λ_3(n) = 1/(2(n²-1)) (平凡零点)")
    
    L1, L2, L3 = compute_Lambda_decomp(Nmax, K_zeros)
    Lambda_approx = L1 - L2 - L3
    
    # 验证分解质量
    print("\n  ★ 分解质量检验 (Λ_exact vs Λ_1-Λ_2-Λ_3):")
    print(f"  {'n':>6} {'Λ_exact':>10} {'Λ_approx':>12} {'差':>10} {'Λ_2':>12} {'Λ_3':>10}")
    test_ns = [2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 100, 1000, 10000]
    for n in test_ns:
        if n <= Nmax:
            print(f"  {n:6d} {Lambda_exact[n]:10.4f} {Lambda_approx[n]:12.4f} {Lambda_exact[n]-Lambda_approx[n]:10.4f} {L2[n]:12.4f} {L3[n]:10.6f}")
    
    # 均方误差
    mask = Lambda_exact > 0
    if mask.sum() > 0:
        rmse = np.sqrt(np.mean((Lambda_exact[mask] - Lambda_approx[mask])**2))
        mean_abs = np.mean(np.abs(Lambda_exact[mask]))
        print(f"\n  RMSE(Λ_exact, Λ_approx) = {rmse:.4f}")
        print(f"  mean|Λ_exact| = {mean_abs:.4f}")
        print(f"  相对RMSE = {rmse/mean_abs:.4f}")
    
    # ---- Step 2: 9项卷积 ----
    print(f"\n[2] ★★★ 9项卷积分解 ★★★")
    print("  R(N) = Σ_n Λ(n)Λ(N-n)")
    print("  Λ = Λ_1 - Λ_2 - Λ_3")
    print("  R = R_11 - R_12 - R_13 - R_21 + R_22 + R_23 - R_31 + R_32 + R_33")
    print("  由对称性 R_ij(N) = R_ji(N) (当 i,j 交换)")
    print("  → R = R_11 - 2R_12 - 2R_13 + R_22 + 2R_23 + R_33")
    
    # 计算各分量
    comps = {1: L1, 2: L2, 3: L3}
    
    print(f"\n  计算 FFT 卷积...")
    R_ij = {}
    for i in range(1, 4):
        for j in range(i, 4):
            conv = fft_convolve_cross(comps[i], comps[j])
            R_ij[(i,j)] = conv
    
    # 输出各 N 的分解
    print(f"\n  ★★★ R(N) 三项和分解结果 ★★★")
    print(f"  {'N':>7} {'R_exact':>10} {'R_approx':>10} {'R_11':>10} {'-2R_12':>10} {'-2R_13':>10} {'R_22':>10} {'2R_23':>10} {'R_33':>10}")
    
    test_Ns = [100, 500, 1000, 2000, 5000, 10000, 20000, 30000, 50000]
    for N in test_Ns:
        if N > Nmax: continue
        R_exact = sum(Lambda_exact[n] * Lambda_exact[N-n] for n in range(1, N))
        R_approx = conv[N]  # 用 Λ_approx 的 FFT 卷积
        
        r11 = R_ij[(1,1)][N]
        r12 = -2 * R_ij[(1,2)][N]
        r13 = -2 * R_ij[(1,3)][N]
        r22 = R_ij[(2,2)][N]
        r23 = 2 * R_ij[(2,3)][N]
        r33 = R_ij[(3,3)][N]
        
        print(f"  {N:7d} {R_exact:10.1f} {R_approx:10.1f} {r11:10.1f} {r12:10.1f} {r13:10.1f} {r22:10.1f} {r23:10.1f} {r33:10.1f}")
    
    # ---- Step 3: 量级分析 ----
    print(f"\n[3] ★★★ 各项量级标度律 ★★★")
    print("  分析各项 R_ij(N) 随 N 的增长率")
    
    Ns_scale = np.array([n for n in range(500, Nmax+1, 500)])
    print(f"\n  {'N':>7}", end="")
    labels = ['R_exact', 'R_11', '|R_12|', '|R_13|', 'R_22', '|R_23|', '|R_33|', 'R_approx']
    for lab in labels:
        print(f" {lab:>10}", end="")
    print()
    
    for N in Ns_scale[::5]:  # 每隔5个取一个
        R_ex = sum(Lambda_exact[n]*Lambda_exact[N-n] for n in range(1,N))
        r11 = abs(R_ij[(1,1)][N])
        r12 = abs(R_ij[(1,2)][N])
        r13 = abs(R_ij[(1,3)][N])
        r22 = abs(R_ij[(2,2)][N])
        r23 = abs(R_ij[(2,3)][N])
        r33 = abs(R_ij[(3,3)][N])
        r_ap = abs(conv[N])
        print(f"  {N:7d} {abs(R_ex):10.1f} {r11:10.1f} {r12:10.1f} {r13:10.1f} {r22:10.1f} {r23:10.1f} {r33:10.1f} {r_ap:10.1f}")
    
    # 标度指数拟合
    print("\n  标度拟合: |R_ij(N)| ~ N^alpha")
    for label, data_fn in [
        ('R_11', lambda N: abs(R_ij[(1,1)][N])),
        ('|R_12|', lambda N: abs(R_ij[(1,2)][N])),
        ('|R_22|', lambda N: abs(R_ij[(2,2)][N])),
        ('|R_23|', lambda N: abs(R_ij[(2,3)][N])),
        ('R_exact', lambda N: abs(sum(Lambda_exact[n]*Lambda_exact[N-n] for n in range(1,N)))),
    ]:
        Ns_fit = np.array([500, 1000, 2000, 5000, 10000, 20000, 40000])
        Ns_fit = Ns_fit[Ns_fit <= Nmax]
        vals = np.array([data_fn(int(N)) for N in Ns_fit])
        vals = np.maximum(vals, 1e-10)
        logN = np.log(Ns_fit)
        logV = np.log(vals)
        # 线性拟合
        coeffs = np.polyfit(logN, logV, 1)
        alpha = coeffs[0]
        print(f"    {label:>10}: α ≈ {alpha:.3f} (N^{alpha:.3f})")
    
    # ---- Step 4: 关键问题——零点振荡项 R_22 有多大？ ----
    print(f"\n[4] ★★★ 关键: 零点振荡卷积 R_22 的量级 ★★★")
    print(f"  R_22(N) = Σ_n Λ_2(n)·Λ_2(N-n)")
    print(f"  Λ_2(n) = n^{{-1/2}} · 2Σ_k cos(γ_k·log n)")
    print(f"  如果 R_22 ~ N^α 且 α < 1, 则零点振荡项本身是 sub-linear")
    print(f"  如果 R_22 ~ N, 则振荡项同阶主项，是核心困难")
    
    for N in [1000, 5000, 10000, 20000, 50000]:
        if N > Nmax: continue
        r22 = R_ij[(2,2)][N]
        r11 = R_ij[(1,1)][N]
        print(f"  N={N:6d}: R_22 = {r22:12.1f}, R_11 = {r11:12.1f}, R_22/R_11 = {r22/r11:.4f}, |R_22|/N = {abs(r22)/N:.4f}")
    
    # ---- Step 5: Λ_2 的统计性质 ----
    print(f"\n[5] ★★★ Λ_2(n) 的统计性质 ★★★")
    mask_nz = np.arange(2, Nmax+1)
    L2_vals = L2[2:Nmax+1]
    print(f"  mean(Λ_2) = {np.mean(L2_vals):.6f}")
    print(f"  std(Λ_2) = {np.std(L2_vals):.6f}")
    print(f"  max|Λ_2| = {np.max(np.abs(L2_vals)):.4f}")
    print(f"  mean|Λ_2| = {np.mean(np.abs(L2_vals)):.4f}")
    print(f"  Σ|Λ_2(n)|/N = {np.sum(np.abs(L2_vals))/Nmax:.4f}")
    print(f"  (Λ_2 的均方大小决定了 R_22 的量级)")
    
    # Parseval 检验
    L2_energy = np.sum(L2_vals**2)
    print(f"\n  Σ Λ_2(n)² = {L2_energy:.2f}")
    print(f"  Σ Λ_2(n)² / N = {L2_energy/Nmax:.4f}")
    print(f"  如果 Λ_2 ~ n^{-1/2}·K^{1/2} (K个零点随机相消), 则 ΣΛ_2² ~ K·Σ1/n ~ K·log N")
    print(f"  预期: ΣΛ_2²/N ≈ K·log(N)/N ≈ {K_zeros}*{math.log(Nmax):.2f}/{Nmax} = {K_zeros*math.log(Nmax)/Nmax:.4f}")
    
    print(f"\n[6] ★★★ 结论与下一步 ★★★")
    print(f"  如果 R_22 << R_11: 零点振荡项可控，三项分解有压制效果")
    print(f"  如果 R_22 ~ R_11: 振荡项同阶主项，三项分解不直接解决问题")
    print(f"  下一步: 叠加分圆塔 CRT，看塔结构能否进一步分解 R_22")

if __name__ == '__main__':
    main()
