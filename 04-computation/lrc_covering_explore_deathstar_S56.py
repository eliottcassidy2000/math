from fractions import Fraction as F
from math import gcd
import random
random.seed(3)
def covers(fam):  # multiple of every q in 2..14
    return all(any(v%q==0 for v in fam) for q in range(2,15))
def M_exact(fam,Qmax=None):
    if Qmax is None: Qmax=2*max(fam)+2
    best=F(0); arg=None
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q); arg=(a,q)
    return best,arg

# 1) natural small covering families and their M
print("(1) natural covering families (must have mult of each q in 2..14):")
naturals={
 "8,9,10,11,12,13,14,+{2,3,4,5,6,7}":[8,9,10,11,12,13,14,2,3,4,5,6,7],
 "small greedy {2..14}":list(range(2,15)),
 "{1..13}+14->{2..14}":list(range(2,15)),
}
for nm,fam in naturals.items():
    fam=sorted(set(fam))
    if len(fam)!=13: print(f"   {nm}: len={len(fam)} skip"); continue
    M,arg=M_exact(fam); print(f"   {nm}: covers={covers(fam)} M={M}={float(M):.4f} @t={arg[0]}/{arg[1]}")

# 2) descent to MINIMIZE M over covering families (find the ~1/9 minimizer)
def rand_covering(Vmax=40):
    while True:
        fam=sorted(random.sample(range(2,Vmax+1),13))
        if covers(fam): return fam
def neighbors(fam,Vmax=40):
    for i in range(13):
        for d in (-2,-1,1,2):
            nv=fam[i]+d
            if 2<=nv<=Vmax and nv not in fam:
                nf=sorted(fam[:i]+[nv]+fam[i+1:])
                if covers(nf): yield nf
best_overall=(F(1),None)
for start in range(40):
    cur=rand_covering(); curM,_=M_exact(cur)
    for _ in range(60):
        improved=False
        for nb in neighbors(cur):
            m,_=M_exact(nb)
            if m<curM: cur,curM=nb,m; improved=True; break
        if not improved: break
    if curM<best_overall[0]: best_overall=(curM,cur)
M,fam=best_overall; Marg=M_exact(fam)
print(f"\n(2) descent MIN M over covering families: M={M}={float(M):.5f} @ {fam}")
print(f"    witness t={Marg[1]}; has mult of 14? {any(v%14==0 for v in fam)}; 1/9={1/9:.5f} 1/14={1/14:.5f}")
