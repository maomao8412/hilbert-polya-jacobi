# -*- coding: utf-8 -*-
import json, numpy as np, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='DejaVu Sans'
plt.rcParams['mathtext.fontset']='dejavusans'
plt.rcParams['axes.unicode_minus']=False

pk=json.load(open('/tmp/grh_product/product_J100_results.json'))
aK=np.array(pk['alpha']); bK=np.array(pk['b'])
zk=json.load(open('/tmp/hpgit/data/jacobi100_checkpoint.json'))
aZ=np.array([float(x) for x in zk['alphas']])
bZ=np.sqrt(np.array([float(x) for x in zk['betas_sq']]))

def det_ratio(a,b,N,u1=1.5,u0=0.5):
    def det(u):
        u2=u*u; u4=u2*u2
        d0=1.0; d1=1.0+u2*a[0]
        for k in range(1,N):
            d2=(1.0+u2*a[k])*d1-u4*b[k-1]**2*d0
            d0,d1=d1,d2
        return d1
    return det(u1)/det(u0)

Ns=list(range(10,101,5))
pi_c=[3.0*det_ratio(aZ,bZ,N) for N in Ns]
G_c=[det_ratio(aK,bK,N) for N in Ns]
pi_true=math.pi
G_lim=(4.0/3.0)*0.915965594177219015054603514932384110774

fig,ax=plt.subplots(1,2,figsize=(12.5,4.8))
ax[0].plot(Ns,pi_c,'o-',color='#c0392b',ms=5,lw=1.6,
           label=r'matrix reading $3\,\det(I+\frac{9}{4}J)/\det(I+\frac{1}{4}J)$')
ax[0].axhline(pi_true,color='#2c3e50',ls='--',lw=1.4,
              label=r'true value $\pi=%.5f$'%pi_true)
ax[0].set_xlabel('Jacobi matrix order N'); ax[0].set_ylabel('determinant reading')
ax[0].set_title(r'$\zeta$ matrix (RH paper): determinants read out $\pi$'+'\n'+
                r'$\xi(2)/\xi(1)=\pi/3$, $\ \xi(1)=1/2$, $\ \xi(2)=\pi/6$')
ax[0].legend(fontsize=8.5,loc='lower right'); ax[0].grid(alpha=.3)
ax[0].annotate('N=100: %.5f'%pi_c[-1],(100,pi_c[-1]),xytext=(38,3.137),
               fontsize=9,arrowprops=dict(arrowstyle='->',color='gray'))

ax[1].plot(Ns,G_c,'s-',color='#2471a3',ms=5,lw=1.6,
           label=r'matrix reading $\det(I+\frac{9}{4}J_K)/\det(I+\frac{1}{4}J_K)$')
ax[1].axhline(G_lim,color='#2c3e50',ls='--',lw=1.4,
              label=r'true value $4G/3=%.5f$ ($G$ = Catalan)'%G_lim)
ax[1].set_xlabel('Jacobi matrix order N'); ax[1].set_ylabel('determinant reading')
ax[1].set_title(r'Product matrix (GRH paper, $\xi_K=\xi\cdot\Lambda_\beta$): reads out Catalan'+'\n'+
                r'$\Lambda_\beta(1)=1$, $\ \Lambda_\beta(2)=4G/\pi$; $\pi$ cancels, $G$ remains')
ax[1].legend(fontsize=8.5,loc='lower right'); ax[1].grid(alpha=.3)
ax[1].annotate('N=100: %.5f'%G_c[-1],(100,G_c[-1]),xytext=(40,1.145),
               fontsize=9,arrowprops=dict(arrowstyle='->',color='gray'))
fig.suptitle(r'Constants follow the number field: one determinant reading — the $\zeta$ matrix gives $\pi$, the Q(i) product matrix gives Catalan',
             fontsize=12,y=1.02)
fig.tight_layout()
fig.savefig('/tmp/img_en/catalan_decode_en.png',dpi=150,bbox_inches='tight')
print("N=100 pi=%.6f  G=%.6f (lim %.6f)"%(pi_c[-1],G_c[-1],G_lim))
print("saved catalan_decode_en")
