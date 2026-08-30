# -*- coding: utf-8 -*-
"""中文版：纯矩阵算术解码素数（trace 路线），生成 4 张中文图。"""
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from mpmath import mp, mpf, mpc, zeta as mpzeta, log as mplog, pi as mppi, gamma as mpgamma, digamma as mpdigamma
fp='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
fm.fontManager.addfont(fp)
plt.rcParams['font.family']=fm.FontProperties(fname=fp).get_name()
plt.rcParams['mathtext.fontset']='dejavusans'
plt.rcParams['axes.unicode_minus']=False
mp.dps=50
DATA='/tmp/hpgit/data/'

def load_J(path):
    cp=json.load(open(path))
    alphas=np.array([float(mpf(s)) for s in cp['alphas']])
    betas=np.array([float(mpf(s)) for s in cp['betas_sq']])**0.5
    n=len(alphas); J=np.zeros((n,n))
    J[np.arange(n),np.arange(n)]=alphas
    J[np.arange(n-1),np.arange(1,n)]=betas
    J[np.arange(1,n),np.arange(n-1)]=betas
    return J

def inv_sqrt_via_polar(J):
    n=J.shape[0]
    L=np.linalg.cholesky(J)
    Linv=np.linalg.solve(L,np.eye(n)); Y=Linv.copy()
    for _ in range(40):
        YinvT=np.linalg.solve(Y.T,np.eye(n)); Ynew=0.5*(Y+YinvT)
        if np.linalg.norm(Ynew-Y)<1e-13*np.linalg.norm(Y): Y=Ynew; break
        Y=Ynew
    U=Y; A=U.T@Linv
    return A

def expm_sq(M,order=12):
    nrm=np.linalg.norm(M,ord=np.inf)
    q=max(0,int(np.ceil(np.log2(nrm+1e-300)))+1)
    E=M/(2.0**q); T=np.eye(M.shape[0],dtype=complex); term=T.copy()
    for k in range(1,order+1): term=term@E/k; T=T+term
    for _ in range(q): T=T@T
    return T

J100=load_J(DATA+'jacobi100_checkpoint.json')
J50 =load_J(DATA+'J50_checkpoint.json')
A100=inv_sqrt_via_polar(J100); A50=inv_sqrt_via_polar(J50)
n100,n50=100,50
B100=np.linalg.solve(0.5*np.eye(n100)+1j*A100,np.eye(n100,dtype=complex))
B50 =np.linalg.solve(0.5*np.eye(n50) +1j*A50, np.eye(n50, dtype=complex))

def psi_trace(x,A,B):
    t=np.log(x); E=expm_sq(1j*A*t)
    zsum=2.0*np.sqrt(x)*np.trace(E@B).real
    return x-zsum-np.log(2*np.pi)-0.5*np.log(1.0-x**-2)

def von_mangoldt(n):
    m=n
    for p in range(2,int(n**0.5)+1):
        if m%p==0:
            while m%p==0: m//=p
            return np.log(p) if m==1 else 0.0
    return np.log(n) if n>1 else 0.0
LAM=np.array([von_mangoldt(n) for n in range(0,3000)])

xmax=120.0
cum=np.cumsum(LAM[:int(xmax)+2])
xs=np.linspace(6,xmax,300)
psi_t=np.interp(xs,np.arange(int(xmax)+2),cum)
psi_m100=np.array([psi_trace(x,A100,B100) for x in xs])
psi_m50 =np.array([psi_trace(x,A50, B50)  for x in xs])
m=xs>40
rms100=np.sqrt(np.mean((psi_m100-psi_t)[m]**2)); rms50=np.sqrt(np.mean((psi_m50-psi_t)[m]**2))
print("trace psi RMS (x>40): J100 %.3f J50 %.3f"%(rms100,rms50))

Nmax=52
spikes=np.zeros(Nmax+1)
for n in range(2,Nmax+1):
    spikes[n]=psi_trace(float(n)+0.5,A100,B100)-psi_trace(float(n)-0.5,A100,B100)
true_pw=[n for n in range(2,Nmax+1) if LAM[n]>0]
recovered=[n for n in range(2,Nmax+1) if spikes[n]>0.8]
hits=len(set(recovered)&set(true_pw))
print("prime powers recovered %d/%d"%(hits,len(true_pw)))

OUT='/tmp/img_zh/'
# ---- 图1: psi from matrix arithmetic ----
fig,ax=plt.subplots(2,1,figsize=(11.5,8.5),gridspec_kw={'height_ratios':[2.3,1]})
ax[0].plot(xs,psi_t,color='black',lw=1.6,label=r'真值 $\psi(x)$')
ax[0].plot(xs,psi_m100,color='#c0392b',lw=1.1,alpha=0.9,
           label=r'$J_{100}$ 迹：$x-2\sqrt{x}\,\mathrm{Re}\,\mathrm{Tr}[e^{iA\log x}(\frac{1}{2}I+iA)^{-1}]-\log 2\pi-\cdots$')
ax[0].plot(xs,xs-np.log(2*np.pi),color='#2c3e50',lw=0.7,ls='--',alpha=0.6)
ax[0].set_ylabel(r'$\psi(x)$',fontsize=12)
ax[0].set_title('纯矩阵算术解码素数——Cholesky、极迭代、平方；不用零点、不用积分',fontsize=11.5,fontweight='bold')
ax[0].legend(loc='upper left',fontsize=9.5); ax[0].grid(alpha=0.25)
ax[1].plot(xs,psi_m100-psi_t,color='#c0392b',lw=0.7,label='误差')
ax[1].axhline(0,color='black',lw=0.5)
ax[1].set_ylabel('误差',fontsize=11); ax[1].set_xlabel('x',fontsize=12)
ax[1].legend(fontsize=9.5); ax[1].grid(alpha=0.25)
ax[1].text(0.98,0.9,f'RMS 误差（x>40）：{rms100:.2f}',transform=ax[1].transAxes,
           ha='right',va='top',fontsize=10.5,bbox=dict(boxstyle='round',fc='#fdf2e9',ec='#c0392b',alpha=0.9))
plt.tight_layout(); plt.savefig(OUT+'primes_matrix_trace_psi_zh.png',dpi=150); plt.close()
print('saved primes_matrix_trace_psi_zh.png')

# ---- 图2: spikes pop out ----
ns=np.arange(2,Nmax+1)
fig,ax=plt.subplots(2,1,figsize=(12,7.5),gridspec_kw={'height_ratios':[1,1]},sharex=True)
ax[0].bar(ns-0.2,LAM[ns],width=0.4,color='black',alpha=0.75,label=r'真值 $\Lambda(n)$')
ax[0].bar(ns+0.2,np.maximum(spikes[ns],0),width=0.4,color='#c0392b',alpha=0.85,
          label=r'矩阵迹 $\psi(n+\frac{1}{2})-\psi(n-\frac{1}{2})$')
ax[0].set_ylabel(r'$\Lambda(n)$',fontsize=12)
ax[0].set_title('素数在整数点跳出来：纯矩阵算术给出 von Mangoldt 尖峰',fontsize=12,fontweight='bold')
ax[0].legend(fontsize=10); ax[0].grid(alpha=0.25,axis='y')
for n in [p for p in ns if LAM[p]>0 and abs(LAM[p]-np.log(p))<1e-9]:
    ax[0].annotate(str(n),(n,LAM[n]),textcoords='offset points',xytext=(0,3),ha='center',fontsize=7.5,color='#1a4d8f')
err_sp=spikes[ns]-LAM[ns]
ax[1].bar(ns,err_sp,width=0.7,color=np.where(np.abs(err_sp)<0.5,'#27ae60','#c0392b'),alpha=0.85)
ax[1].axhline(0,color='black',lw=0.5)
ax[1].set_ylabel('解码误差',fontsize=11); ax[1].set_xlabel('n',fontsize=12)
ax[1].set_title(f'到 {Nmax} 为止恢复素数幂 {hits}/{len(true_pw)}（绿色：|误差|<0.5）',fontsize=10.5)
ax[1].grid(alpha=0.25,axis='y')
plt.tight_layout(); plt.savefig(OUT+'primes_pop_out_spikes_zh.png',dpi=150); plt.close()
print('saved primes_pop_out_spikes_zh.png')

# ---- 图3: J50 vs J100 ----
fig,ax=plt.subplots(figsize=(11.5,6))
ax.plot(xs,psi_t,color='black',lw=1.7,label=r'真值 $\psi(x)$')
ax.plot(xs,psi_m50,color='#2980b9',lw=1.0,alpha=0.85,label=r'$J_{50}$ 矩阵算术')
ax.plot(xs,psi_m100,color='#c0392b',lw=1.0,alpha=0.9,label=r'$J_{100}$ 矩阵算术')
ax.set_xlabel('x',fontsize=12); ax.set_ylabel(r'$\psi(x)$',fontsize=12)
ax.set_title('矩阵元素越多，素数越多：迹随矩阵阶数变锐利',fontsize=12)
ax.legend(loc='upper left',fontsize=10.5); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(OUT+'matrix_trace_order_zh.png',dpi=150); plt.close()
print('saved matrix_trace_order_zh.png')

# ---- 图4: zeta is a matrix determinant ----
us=np.linspace(-3,3,240)
logdet=np.array([np.linalg.slogdet(np.eye(n100)+u*u*J100)[1] for u in us])
def xi_mp(s):
    return 0.5*s*(s-1)*mppi**(-s/2)*mpgamma(s/2)*mpzeta(s)
logdet_exact=np.array([2.0*float(mplog(abs(xi_mp(mpf('0.5')+u)/xi_mp(mpf('0.5'))))) for u in us[::4]])
fig,ax=plt.subplots(figsize=(11,5.5))
ax.plot(us,logdet,color='#c0392b',lw=1.6,label=r'矩阵：$\log\det(I+u^2 J_{100})$')
ax.plot(us[::4],logdet_exact,'k--',lw=1.0,label=r'精确：$2\log|\,\xi(\frac{1}{2}+u)/\xi(\frac{1}{2})\,|$')
ax.set_xlabel(r'$u=s-\frac{1}{2}$（实数）',fontsize=12)
ax.set_ylabel(r'$\log\det$',fontsize=12)
ax.set_title(r'$\xi(s)/\xi(\frac{1}{2})=\det(I+u^2J)$ 的收敛带：zeta 函数就是矩阵行列式'
             r'（有限 $J_{100}$ 精确分辨 $|u|<3$ 带；带宽随阶数增大）',fontsize=11.5,fontweight='bold')
ax.legend(fontsize=11); ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(OUT+'zeta_is_determinant_zh.png',dpi=150); plt.close()
print('saved zeta_is_determinant_zh.png')
