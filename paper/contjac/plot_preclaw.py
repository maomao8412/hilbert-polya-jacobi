# -*- coding: utf-8 -*-
"""失败阶 vs 工作精度 诊断图（论文 Figure 1）。"""
import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE='/Coze/Drive/黎曼猜想论文审核/所有对话/主对话/numerical_algorithms_paper/'
lo=json.load(open(BASE+'hankel_diagnostic_lodeps.json'))['rows']
hi=json.load(open(BASE+'hankel_diagnostic.json'))['rows']
rows=lo[:3]+hi  # 16,30,50 + 100..1300(实际到800, 若hideps合并了再补)
import os
hp=BASE+'hankel_diagnostic_hideps.json'
if os.path.exists(hp):
    rows += json.load(open(hp))['rows']
rows=[r for r in rows if r['cholesky_fail_order']]
D=[r['dps'] for r in rows]; nf=[r['cholesky_fail_order'] for r in rows]
# 只对 n<100 的点做线性拟合（1300 若 >100 被 nmax=100 截断则排除）
fit=[(d,n) for d,n in zip(D,nf) if n<99]
A=np.polyfit([d for d,_ in fit],[n for _,n in fit],1)
xs=np.linspace(0,1350,100)
plt.figure(figsize=(6.2,4.2))
plt.scatter(D,nf,s=45,zorder=3,color='#1f77b4',label='Cholesky failure order $n_{\\rm fail}$')
trunc=[(d,n) for d,n in zip(D,nf) if n>=99]
if trunc: plt.scatter([d for d,_ in trunc],[99 for _ in trunc],marker='^',s=70,color='#d62728',zorder=4,label='no failure up to order 100')
plt.plot(xs,np.polyval(A,xs),'--',color='gray',label='linear fit: $n_{\\rm fail}\\approx %.2f+%.3f\\,D$'%(A[1],A[0]))
plt.xlabel('working precision (decimal digits)'); plt.ylabel('Jacobi / Hankel failure order')
plt.title('Riemann $\\xi$, $R=2$, $M=1536$: precision budget of the moment pipeline')
plt.grid(alpha=.3); plt.legend(fontsize=8,loc='upper left')
plt.tight_layout(); plt.savefig(BASE+'figure_preclaw.png',dpi=160)
print('rows:',sorted(zip(D,nf)))
print('fit slope=%.4f orders/digit (%.1f per 100 digits), intercept=%.1f'%(A[0],A[0]*100,A[1]))
print('saved figure_preclaw.png')
