# A-PRIORI existence via AVERAGING (kps-S95). If avg_j maxgap(j) > Vmax/7, then max_j maxgap > Vmax/7
# => a good period EXISTS -- no exhaustion, no certificate. Lemniscate inspiration: bound maxgap by the
# SECOND MOMENT of gaps (maxgap >= sum gap^2 / Vmax, since maxgap*sum gap >= sum gap^2), and the
# 2nd moment averaged over dilations is an ADDITIVE-ENERGY object (dissociated => big => big gap).
# Test: A = avg_j (7 maxgap/Vmax);  B = avg_j (7 sum gap^2 / Vmax^2) [the 2nd-moment lower bound on A].
import numpy as np
from math import gcd, floor
from functools import reduce
import random
random.seed(95011)
def prim(E):
    E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best
def stats(E,Vmax):
    Ea=np.array(sorted(E)); j=np.arange(1,Vmax)
    ph=np.mod(np.outer(j,Ea),Vmax); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+Vmax-ph[:,-1])[:,None]],axis=1).astype(float)
    mg=g.max(axis=1)                 # maxgap per j
    sq=(g*g).sum(axis=1)             # sum of gap^2 per j
    A=7*mg.mean()/Vmax               # avg (7 maxgap/Vmax)
    B=7*sq.mean()/Vmax**2            # avg (7 sumgap^2/Vmax^2)  <= A (since maxgap>=sumgap^2/Vmax)
    return A,B
def rand_diss(s,Lcap):
    for _ in range(300):
        mid=sorted(random.sample(range(1,s),11)); E=[0]+mid+[s]
        if len(set(E))==13 and prim(E) and longest_ap(E)<=Lcap: return E
    return None
print("A = avg_j(7 maxgap/Vmax); if A>1 => good period EXISTS (max>=avg). B = 2nd-moment lower bound.")
print("At critical Vmax=floor(7s/6) (hardest). Adversarial: MINIMIZE A over dissociated clusters.")
print(f"{'s':>4}{'Vmax':>6}{'minA':>8}{'minB':>8}{'A>1?':>6}")
worstA=(9.9,None)
for s in (17,20,25,30,40,60,90,120,160,200):
    Vmax=floor(7*s/6); best=(9.9,None,None)
    for _r in range(50):
        E=rand_diss(s,7)
        if E is None: continue
        A,B=stats(E,Vmax); cur=(A,tuple(E),B)
        # hill-climb to minimize A
        improved=True
        while improved:
            improved=False
            for idx in range(1,12):
                for d in (-2,-1,1,2):
                    nv=E[idx]+d
                    if nv<=0 or nv>=s or nv in E: continue
                    E2=sorted(E[:idx]+[nv]+E[idx+1:])
                    if not prim(E2) or longest_ap(E2)>7: continue
                    A2,B2=stats(E2,Vmax)
                    if A2<cur[0]: cur=(A2,tuple(E2),B2); E=E2; improved=True; break
                if improved: break
        if cur[0]<best[0]: best=cur
    tag="YES" if best[0]>1 else "NO"
    if best[0]<worstA[0]: worstA=(best[0],(s,Vmax,best[1]))
    print(f"{s:>4}{Vmax:>6}{best[0]:>8.4f}{best[2]:>8.4f}{tag:>6}")
print(f"\nglobal min avg-maxgap ratio A = {worstA[0]:.4f} at s,Vmax,E={worstA[1]}")
print("If minA>1 for s>=S0 => averaging PROVES existence for spread>=S0 (elementary, no exhaustion).")
print("Then only s<S0 needs exhaustion. B<=A shows how much the 2nd-moment bound alone gives.")
