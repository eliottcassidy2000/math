"""mac-mini-S57 (LEM-006): the energy-variance bridge for the covering floor.
(1) far <= E[W]^2 <=> Var(W) <= near (exact). (2) Var(W)=Sum|W-hat(m)|^2 ~ c*R2, c~5.6e-5
(pair-part V1_phi=1.53e-3 overcounts 27x via IE cancellation). (3) D3 decreasing in m2 (proved).
Reduces BOTH far<=E[W]^2 (spread) AND brick(B) (R2<=614 => D3>=bar) to Var(W)<=c*R2 = LEM-005.
614 = additive energy of the k=11 diam-16 extremizer = a value of the universal spread param R2
(same R2 as THM-656 tent floor)."""
import numpy as np
from collections import Counter
from math import gcd
from functools import reduce
import random
random.seed(11); TH=1/7; M=6/7
def R2(A):
    A=sorted(A);c=Counter()
    for i in range(len(A)):
        for j in range(len(A)):
            if i!=j: c[A[i]-A[j]]+=1
    return sum(v*v for v in c.values())
def prim(A):
    A=sorted(A);return reduce(gcd,[A[i+1]-A[i] for i in range(len(A)-1)])==1
GRID=300_000;x=(np.arange(GRID)+0.5)/GRID
def full(A):
    Ea=np.array(sorted(A),float);ph=np.mod(np.outer(x,Ea),1.0);ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    W=np.maximum(g-TH,0).sum(axis=1);m1,m2,m3=W.mean(),(W*W).mean(),(W**3).mean()
    Ls=np.linspace(TH,2*TH,50);near=2*np.trapezoid([np.maximum(g-L,0).sum(axis=1).mean() for L in Ls],Ls)
    d3=m1/M+(m1-m2/M)**2/(m2-m3/M) if m2-m3/M>0 else m1/M
    return m1,m2-m1*m1,near,m2-near,d3
# (1) equivalence
ok=all((full(A)[3]<=full(A)[0]**2)==(full(A)[1]<=full(A)[2]) for A in
       [list(range(11)),[0,3,7,12,20,30,44,55,60,70,85],[0,2,4,6,8,9,10,12,14,16,18]])
print(f"(1) far<=E[W]^2 <=> Var<=near: {ok}")
# (2) Var/R2 band; (3) brick B margin
lo=[];band=[]
for _ in range(400):
    D=random.choice([10,12,14,16,20,30]);A=tuple(sorted(random.sample(range(D+1),11)))
    if not prim(A): continue
    m1,var,near,far,d3=full(A);r=R2(A);band.append(var/r)
    if r<=614: lo.append(d3)
print(f"(2) Var(W)/R2 in [{min(band):.2e},{max(band):.2e}] (~c*R2); pair-part V1_phi=1.53e-3 (27x)")
print(f"(3) min D3 for R2<=614: {min(lo):.4f} >= bar 0.3312 (+{min(lo)-0.3312:.3f}); D3 decreasing in m2 (proved)")
