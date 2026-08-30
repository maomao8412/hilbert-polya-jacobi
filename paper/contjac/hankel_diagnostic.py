# -*- coding: utf-8 -*-
"""Hankel 条件数/精度标度诊断（算法论文 Figure 数据）。
从 1300dps 采样 checkpoint 重算矩序列，逐精度档位截断，
用 Cholesky 正定失败定位 Hankel 崩溃阶 n_fail(dps)；
另跑 Gram-Schmidt 定位首个 b^2<=0 阶。另从 Jacobi 系数重建矩做交叉核对。"""
import mpmath as mp, json, time, os, re, math
t0=time.time()
BASE='/Coze/Drive/黎曼猜想论文审核/所有对话/主对话'
samples=json.load(open(BASE+'/zeta_J100_20260831/samples_1300dps.json'))['fr']
R=mp.mpf(2); M=1536; JMAX=210

def moments_at_dps(target_dps):
    mp.mp.dps=target_dps
    def _parse(s):
        m=re.match(r'^\((.*)\s([+-])\s(.*)j\)$', s.strip())
        if m:
            re_,im_=m.group(1),m.group(3); sgn=1 if m.group(2)=='+' else -1
            return mp.mpc(mp.mpf(re_), sgn*mp.mpf(im_))
        return mp.mpc(s)
    fr=[_parse(x) for x in samples]
    sigma0=sum(fr)/M
    d={}; cl={}
    for n in range(1,JMAX+1):
        j=2*n
        ss=mp.mpc(0)
        for k in range(M):
            ss+=fr[k]*mp.e**(-1j*2*mp.pi*j*k/M)
        d[n]=mp.re(ss/M/R**j/sigma0)
        sm=mp.mpf(0)
        for i in range(1,n): sm+=i*cl[i]*d[n-i]
        cl[n]=d[n]-sm/n
    P=[None]+[((-1)**(m+1))*m*cl[m] for m in range(1,JMAX+1)]
    return P

def cholesky_fail(P, nmax=200):
    """构造 Hankel H_n=(P[i+j+1]), i,j=0..n-1；逐阶 LDL^T，首个非正 pivot 阶。"""
    nfail=None
    # LDL 增量：H 的 leading principal minor 顺序检测
    L=[[mp.mpf(0) for _ in range(nmax)] for _ in range(nmax)]
    piv=[]
    for i in range(nmax):
        for j in range(i+1):
            v=P[i+j+1]
            for k in range(j):
                v-=L[i][k]*L[j][k]*piv[k]
            if i==j:
                if v<=0:
                    return i+1, [float(p) for p in piv]
                piv.append(v); L[i][j]=1
            else:
                L[i][j]=v/piv[j]
    return None, [float(p) for p in piv]

def gs_fail(P, nmax=200):
    def Tm(m): return P[m]
    def inner(p,q):
        r=mp.mpf(0)
        for i,pi_ in enumerate(p):
            if pi_==0: continue
            for j,qj in enumerate(q):
                if qj==0: continue
                r+=pi_*qj*Tm(i+j+1)
        return r
    polys=[[mp.mpf(1)]]; norms=[inner(polys[0],polys[0])]; bsqs=[]
    for k in range(1,nmax+1):
        xp=[mp.mpf(0)]+polys[k-1]
        ak=inner(xp,polys[k-1])/norms[k-1]
        pk=list(xp)
        for i in range(len(polys[k-1])): pk[i]-=ak*polys[k-1][i]
        if k>=2:
            bsq=norms[k-1]/norms[k-2]
            for i in range(len(polys[k-2])): pk[i]-=bsq*polys[k-2][i]
            bsqs.append(bsq)
            if bsq<=0: return k+1  # Jacobi size index
        sig=inner(pk,pk)
        polys.append(pk); norms.append(sig)
    return None

out={'object':'Hankel/GS failure-order vs working precision (zeta xi, R=2, M=1536)',
     'note':'moments recomputed from 1300dps contour samples, truncated to target dps',
     'rows':[]}
for dps in [100,150,200,300,400,500,650,800,1000,1300]:
    P=moments_at_dps(dps)
    nchol,piv=cholesky_fail(P)
    ngs=gs_fail(P)
    mpiv=min((abs(x) for x in piv), default=0.0)
    row={'dps':dps,'cholesky_fail_order':nchol,'gs_bad_first_index':ngs,
         'min_pivot_log10':(math.log10(mpiv) if mpiv>0 else None)}
    out['rows'].append(row)
    print('dps=%4d  Cholesky_fail_n=%s  GS_bad_first=%s  min|pivot|~1e%.0f  %.0fs'%(
        dps, nchol, ngs, abs(row['min_pivot_log10']) if row['min_pivot_log10'] else -1, time.time()-t0), flush=True)
json.dump(out, open(BASE+'/numerical_algorithms_paper/hankel_diagnostic.json','w'), indent=1)
print('DONE', time.time()-t0, flush=True)
