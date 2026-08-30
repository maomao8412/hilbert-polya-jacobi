# -*- coding: utf-8 -*-
"""中文版：矩阵幂 = ζ 矩（matrix_powers_are_moments）。"""
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from mpmath import mp, mpf, log as mplog, pi as mppi, gamma as mpgamma, zeta as mpzeta, taylor as mptaylor
fp='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
fm.fontManager.addfont(fp)
plt.rcParams['font.family']=fm.FontProperties(fname=fp).get_name()
plt.rcParams['mathtext.fontset']='dejavusans'
plt.rcParams['axes.unicode_minus']=False
mp.dps=60

cp=json.load(open('/tmp/hpgit/data/jacobi100_checkpoint.json'))
alphas=np.array([float(mpf(s)) for s in cp['alphas']])
betas=np.array([float(mpf(s)) for s in cp['betas_sq']])**0.5
n=len(alphas); J100=np.zeros((n,n))
J100[np.arange(n),np.arange(n)]=alphas
J100[np.arange(n-1),np.arange(1,n)]=betas
J100[np.arange(1,n),np.arange(n-1)]=betas

def xi_mp(s):
    return 0.5*s*(s-1)*mppi**(-s/2)*mpgamma(s/2)*mpzeta(s)

ks=np.arange(1,9)
trJk=[np.trace(np.linalg.matrix_power(J100,k)) for k in ks]
f=lambda u: mplog(xi_mp(mpf('0.5')+u)/xi_mp(mpf('0.5')))
coefs=mptaylor(f,0,16)
S_exact=[float(k*abs(coefs[2*k])) for k in ks]

fig,ax=plt.subplots(1,2,figsize=(12,4.6))
ax[0].bar(ks-0.18,trJk,width=0.36,color='#c0392b',label=r'$\mathrm{Tr}(J^k)$，来自 199 个矩阵元素')
ax[0].bar(ks+0.18,S_exact,width=0.36,color='#2c3e50',alpha=0.55,label=r'$S_k=\sum_\gamma\gamma^{-2k}$，来自 $\xi$')
ax[0].set_yscale('log'); ax[0].set_xlabel('k'); ax[0].set_ylabel('矩值')
ax[0].set_title(r'矩阵幂就是 ζ 矩：$\mathrm{Tr}(J^k)=S_k$',fontsize=12)
ax[0].set_xticks(ks); ax[0].legend(fontsize=9.5); ax[0].grid(alpha=0.25,axis='y')
errs=[abs(trJk[k-1]-S_exact[k-1])/S_exact[k-1] for k in ks]
ax[1].semilogy(ks,errs,'o-',color='#c0392b')
ax[1].set_xlabel('k'); ax[1].set_ylabel('相对误差')
ax[1].set_title('矩的吻合（k=1 的缺口 = 超过 100 阶的素数信息）',fontsize=11)
ax[1].set_xticks(ks); ax[1].grid(alpha=0.25)
for k,e in zip(ks,errs):
    if e>1e-12: ax[1].annotate(f'{e:.0%}',(k,e),textcoords='offset points',xytext=(0,8),ha='center',fontsize=9)
plt.tight_layout(); plt.savefig('/tmp/img_zh/matrix_powers_are_moments_zh.png',dpi=150,bbox_inches='tight')
print('saved moments_zh; errs:', ['%.2e'%e for e in errs])
