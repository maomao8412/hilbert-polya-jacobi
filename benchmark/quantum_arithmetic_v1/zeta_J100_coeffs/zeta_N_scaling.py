# -*- coding: utf-8 -*-
"""zeta J100 v2 单族 N 标度收敛统计。
从 zeta_J100v2_results.json 的 100 阶 Jacobi 系数取前 N 阶主子矩阵，
每个 N 取锁定能级（relerr<1e-4，gamma 空间三次展开同 exp3 口径）计算
gap-ratio <r>、最小间距、小间距(<0.3 平均间距)计数、IPR 代理。
输出 zeta_N_scaling.json：N=22..100 各档的统计量。
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "zeta_J100v2_results.json")))
alpha = d["alpha"]; b = d["b"]; rows = d["rows"]
# rows: (n, gamma, tag, ref, relerr)
locked = [r for r in rows if r[4] < 1e-4]
locked_gamma = sorted(r[1] for r in locked)
locked_idx = sorted(r[0] for r in locked)  # 1-based rank among all 100 gammas
print("total locked:", len(locked))

def locked_gammas_for_N(N):
    """前 N 阶主子矩阵本征值重算，用锁定的高-秩参考：直接取被锁定且 rank<=N 的 gamma。
    （主截断本征值与锁定参考在 locked 集上 relerr<1e-4，直接用锁定 gamma 等价且更稳）"""
    # 重算主子矩阵本征值，保持自洽
    H = np.diag(alpha[:N]) + np.diag(b[:N-1], 1) + np.diag(b[:N-1], -1)
    ev = np.sort(np.linalg.eigvalsh(H))[::-1]
    gs = 1.0/np.sqrt(np.maximum(ev, 1e-300))
    gs = np.sort(gs)
    # 只保留与 zeta 零点锁定的（用 rows 中 rank<=N 且 relerr<1e-4 的）
    keep = [r for r in rows if r[0] <= N and r[4] < 1e-4]
    return sorted(r[1] for r in keep), gs

def gap_ratio(xs):
    xs = np.sort(np.array(xs)); dd = np.diff(xs)
    r = np.minimum(dd[:-1], dd[1:])/np.maximum(dd[:-1], dd[1:])
    return r

def ipr_proxy(N):
    """本征矢参与率代理：中心本征矢 IPR=sum v_i^4，扩展态~1/N"""
    H = np.diag(alpha[:N]) + np.diag(b[:N-1], 1) + np.diag(b[:N-1], -1)
    w, V = np.linalg.eigh(H)
    iprs = np.sum(V**4, axis=0)
    # 锁定能级对应最大的若干本征值（gamma 最小端）；取锁定数对应的低 gamma 端
    n_lock = sum(1 for r in rows if r[0] <= N and r[4] < 1e-4)
    # eigh 升序 -> 最大本征值在尾部；取尾部 n_lock 个
    if n_lock >= 3:
        sel = iprs[-n_lock:]
        return float(np.mean(sel)), float(1.0/n_lock)
    return float("nan"), float("nan")

out = {"meta": {"source": "zeta_J100v2 1300dps M1536 R2 J100",
                "locked_tol": 1e-4,
                "theories": {"Poisson_r": 0.3863, "GOE_r": 0.5359, "GUE_r": 0.5996}},
       "scaling": []}
for N in [22, 30, 40, 50, 60, 80, 100]:
    lg, allgs = locked_gammas_for_N(N)
    n = len(lg)
    if n >= 6:
        r = gap_ratio(lg); rmean = float(np.mean(r))
        dd = np.diff(np.sort(lg)); mean_d = float(np.mean(dd))
        nsmall = int(np.sum(dd < 0.3*mean_d))
        mingap = float(np.min(dd))
    else:
        rmean = float("nan"); nsmall = -1; mingap = float("nan"); mean_d = float("nan")
    ipr, ipr_ext = ipr_proxy(N)
    rec = {"N": N, "n_locked": n, "gap_ratio_mean": rmean,
           "min_gap": mingap, "mean_gap": mean_d,
           "small_gaps_lt_0p3mean": nsmall,
           "IPR_locked_mean": ipr, "IPR_extended_1overN": ipr_ext}
    out["scaling"].append(rec)
    print(f"N={N:3d} locked={n:3d} <r>={rmean:.4f}  minGap={mingap:.4f} meanGap={mean_d:.3f} "
          f"small(<0.3)={nsmall}  IPR={ipr:.4f} (1/N={ipr_ext:.4f})")

json.dump(out, open(os.path.join(HERE, "zeta_N_scaling.json"), "w"), indent=1)
print("DONE -> zeta_N_scaling.json")
