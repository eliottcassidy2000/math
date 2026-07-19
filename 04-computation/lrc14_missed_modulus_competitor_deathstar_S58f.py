"""The missed-modulus foreign-denominator competitor (death-star-S58f).
For any speed set V, let k' = smallest k>=2 with k | no v_i. Then t=1/k' is a competitor:
min_i ||v_i/k'|| >= 1/k' (any k'-nonmultiple is >=1/k' from Z), so M(V) >= 1/k'.
The two tied runners have residues j,k'-j mod k' => k' | (v_i+v_j)  [= the q'|(v_i+v_j) form].
CONSEQUENCE: M(V)<1/13 => V covers {2,...,13}; M(V)<1/14 => V covers {2,...,14} (Cover14).
=> every strict-interior (1/14<M<1/13) family covers 2..13; the AP-extraction KERNEL reduces to
COVERING-2..13 families (the Freiman/E3 inverse case). Non-covering families have M>=1/13 outright."""
from fractions import Fraction as F
from math import gcd
import random
def Mexact(V):
    mx=max(V);Qs=set()
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for d in range(2,2*mx+1):
                if s%d==0:Qs.add(d)
    best=F(0)
    for q in Qs:
        for a in range(1,q):
            if gcd(a,q)!=1:continue
            m=min(min((v*a)%q,q-((v*a)%q)) for v in V)
            if F(m,q)>best:best=F(m,q)
    return best
random.seed(1);bad=0;checked=0
for _ in range(400):
    V=sorted(random.sample(range(1,60),13))
    missed=[k for k in range(2,20) if all(v%k for v in V)]
    if not missed:continue
    kp=missed[0];M=Mexact(V);checked+=1
    if not (M>=F(1,kp)):bad+=1;print("VIOLATION",V,kp,M)
print("checked=%d missing-a-modulus families; M>=1/k' violations=%d"%(checked,bad))
for V in [list(range(1,13))+[182],list(range(1,12))+[13,156]]:
    M=Mexact(V);miss=[k for k in range(2,14) if all(v%k for v in V)]
    print("V max=%d M=%s strict-interior=%s covers2..13=%s"%(V[-1],M,F(1,14)<M<F(1,13),miss==[]))
