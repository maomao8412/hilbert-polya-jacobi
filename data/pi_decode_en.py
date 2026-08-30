# -*- coding: utf-8 -*-
import json
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='DejaVu Sans'
plt.rcParams['mathtext.fontset']='dejavusans'
plt.rcParams['axes.unicode_minus']=False
mp.mp.dps=60

def logdet(alphas,betas_sq,u):
    n=len(alphas); J=mp.zeros(n,n)
    for i in range(n):
        J[i,i]=mp.mpf(alphas[i])
        if i+1<n: J[i,i+1]=J[i+1,i]=mp.sqrt(mp.mpf(betas_sq[i]))
    L=mp.cholesky(mp.eye(n)+u*u*J)
    return 2*sum(mp.log(L[i,i]) for i in range(n))

z50=json.load(open('/tmp/hpgit/data/J50_checkpoint.json'))
z100=json.load(open('/tmp/hpgit/data/jacobi100_checkpoint.json'))
def pi_est(z,N):
    al,bs=z['alphas'],z['betas_sq']
    r=mp.e**(logdet(al[:N],bs[:N-1],mp.mpf('1.5'))-logdet(al[:N],bs[:N-1],mp.mpf('0.5')))
    return float(3*r)
N50=[10,20,30,40,50]; v50=[pi_est(z50,N) for N in N50]
N100=[10,20,30,40,50,75,100]; v100=[pi_est(z100,N) for N in N100]
lam_star=1+2**0.5; r_star=2**0.5-1

fig,ax=plt.subplots(1,2,figsize=(13,5))
ax[0].plot(N50,v50,'o-',color='#c0392b',label=r'$J_{50}$ ($\zeta$ matrix)')
ax[0].plot(N100,v100,'s-',color='#1a4d8f',label=r'$J_{100}$ ($\zeta$ matrix)')
ax[0].axhline(3.14159265,color='k',ls='--',lw=1,label=r'$\pi$ = 3.14159...')
ax[0].set_xlabel('matrix order N')
ax[0].set_ylabel(r'$3\,\det(I+1.5^2J_N)/\det(I+0.5^2J_N)$')
ax[0].set_title(r'Decoding $\pi$ from pure matrix determinants: $\xi(2)/\xi(1)=\pi/3$'+'\n'+
                r'(matrix multiplication + Cholesky only; no zeros, integrals, primes)')
ax[0].legend(loc='lower right'); ax[0].grid(alpha=.3); ax[0].set_ylim(3.0,3.2)

lams=[1.05+(10-1.05)*i/399 for i in range(400)]
r1=[1/l for l in lams]; r2=[(l-1)/(l+1) for l in lams]; rmax=[max(a,b) for a,b in zip(r1,r2)]
ax[1].plot(lams,r1,color='#1a6b2a',label=r'$1/\lambda$')
ax[1].plot(lams,r2,color='#8a6d1a',label=r'$(\lambda-1)/(\lambda+1)$')
ax[1].plot(lams,rmax,color='#c0392b',lw=2.2,label=r'$r(\lambda)=\max$ (dominant per-step ratio)')
ax[1].axvline(lam_star,color='k',ls=':',lw=1.2)
ax[1].plot(lam_star,r_star,'k*',ms=16)
ax[1].annotate(r'$\lambda^*=1+\sqrt{2}$'+'\n'+r'$r_{\min}=\sqrt{2}-1=q$ (seed)',
    xy=(lam_star,r_star),xytext=(3.4,.55),fontsize=11,
    arrowprops=dict(arrowstyle='->',color='k'))
ax[1].set_xlabel(r'parameter $\lambda$'); ax[1].set_ylabel(r'per-step convergence ratio $r$')
ax[1].set_title(r'Minimax optimality of the parametrised $\pi$ series: seed $q=\sqrt{2}-1$'+'\n'+
                r'fastest-converging point of the double-arctangent parametrisation')
ax[1].legend(loc='upper right'); ax[1].grid(alpha=.3); ax[1].set_ylim(0,1)
plt.tight_layout()
plt.savefig('/tmp/img_en/pi_decode_en.png',dpi=130,bbox_inches='tight')
print('saved pi_decode_en; J100 pi=%.4f'%(v100[-1],))
