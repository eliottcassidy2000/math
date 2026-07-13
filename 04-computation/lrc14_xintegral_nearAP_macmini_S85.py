#!/usr/bin/env python3
"""mac-mini-S85: the POSITIVE-DEFINITE x-INTEGRAL on the near-AP residue (klein-S286 THM-729 device).
For |core|=1, the covering x-integral is eps_1 = INT g(x) 1_{G'}(x) dx (g=runner-1 recentred,
1_{G'}=smooth-safe set = PROD(1-1_{D_w})). Good-period exists <=> the x-integral of the good
indicator > 0 <=> coreCover<1. THM-729: use the POSITIVE-DEFINITE autocorrelation of 1_{G'},
NOT the divergent Fourier expansion. TEST: (1) eps_1 matches my coreCover (connect frameworks);
(2) the autocorrelation A(y)=INT 1_{G'}(x)1_{G'}(x+y) structure; (3) does the pos-def form give a
lower bound on meas(G' n core-safe) (the safe set) where the union bound (S73) FAILED?"""
import numpy as np
c=1.0/14; P=[2,3,5,7,11,13]
def cop(v): return all(v%p!=0 for p in P)
def indicators(S,N):
    x=(np.arange(N)+0.5)/N
    core=[v for v in S if cop(v)]; nc=[v for v in S if not cop(v)]
    def dang(v): r=(v*x)%1.0; return (np.minimum(r,1-r)<c).astype(float)
    Gp=np.ones(N)
    for w in nc: Gp*= (1-dang(w))     # 1_{G'} = smooth safe set
    return x,core,nc,Gp,dang
for nm,S in [("{1..11,13,84}",[*range(1,12),13,84]),("deep well {1..12,182}",[*range(1,13),182])]:
    S=sorted(set(S)); N=200000
    x,core,nc,Gp,dang=indicators(S,N)
    mGp=Gp.mean()
    v=core[0]  # =1
    D1=dang(v)
    safe = Gp*(1-D1)                       # G' AND runner-1 safe = full safe set
    L=safe.mean(); coreCover=1-L/mGp
    # x-integral eps_1 = INT g(x)1_{G'} with g = 1/6 - (7/6)1_{D1} (mean-0 recentred safe... use danger form)
    g = -D1 + (1-D1)/6                      # g_1 recentred: -1 on danger,+1/6 on safe (mean 0 over [0,1))
    eps1 = (g*Gp).mean()
    # autocorrelation of 1_{G'} at a few lags (positive-definite: A(0)=INT Gp^2)
    A0=(Gp*Gp).mean()
    # pos-def 2nd-moment (Paley-Zygmund style) lower bound on meas(safe): needs the good-indicator 2nd moment
    print(f"{nm}: core={core}, meas(G')={mGp:.5f}, L=safe={L:.5f}, coreCover={coreCover:.4f}")
    print(f"   x-integral eps_1 = INT g*1_G' = {eps1:.5f}  (mean(g)=0; eps_1<0 <=> runner-1 danger over-represented in G')")
    print(f"   autocorrelation A(0)=INT 1_G'^2 = {A0:.5f} (=meas(G') since indicator); A(0)/meas(G')^2 = {A0/mGp**2:.2f}")
    # KEY: safe = Gp*(1-D1); L = mGp - INT Gp*D1. INT Gp*D1 = the runner-1 danger mass in G'.
    dmass=(Gp*D1).mean()
    print(f"   runner-1 danger mass in G' = INT 1_G' 1_D1 = {dmass:.5f} = coreCover*meas(G'); L=meas(G')-this={mGp-dmass:.5f}")
print("\nCONNECTION: the covering x-integral eps_1 on the |core|=1 near-AP residue IS my coreCover")
print("object (INT g*1_G'). L(safe) = meas(G') - [runner-1 danger mass in G']. The pos-def THM-729")
print("device would lower-bound meas(safe) via the AUTOCORRELATION of 1_G' WITHOUT Fourier-expanding")
print("(S73 showed the union bound fails: danger mass ~ meas(G')). HONEST: the pos-def x-integral is")
print("the right METRIC form (>=0, no divergence), but bounding [danger mass < meas(G')] on the near-AP")
print("residue is the covering middle-order sum = the live open crux (opus/kps own it). I provide the")
print("concrete near-AP x-integral data + the coreCover<->eps_v connection for the fleet's device.")
