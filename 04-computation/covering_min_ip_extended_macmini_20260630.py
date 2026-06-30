import functools, math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
print=functools.partial(print,flush=True)
def Mexact(S):
    S=sorted(set(S)); cand=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j],S[i]+S[j]):
                if d>0:
                    for k in range(1,d): cand.add(F(k,d))
    best=F(0); at=None
    for t in cand:
        g=min(min((v*t)%1,1-((v*t)%1)) for v in S)
        if g>best: best=g; at=t
    return best,at
def norm(x):
    f=x-int(x); f=f+1 if f<0 else f
    return min(f,1-f)
def cmin(n,V):
    speeds=list(range(1,V+1)); P=len(speeds)
    univ=sorted({F(k,d) for d in range(2,2*V+1) for k in range(1,d) if math.gcd(k,d)==1}); U=len(univ)
    dist=[[norm(speeds[p]*univ[u]) for u in range(U)] for p in range(P)]
    rvals=sorted({dist[p][u] for p in range(P) for u in range(U)})
    base=[[1.0]*P]; lb=[n-1]; ub=[n-1]
    for q in range(2,n+1): base.append([1.0 if v%q==0 else 0.0 for v in speeds]); lb.append(1); ub.append(P)
    primes=[p for p in range(2,V+1) if all(p%d for d in range(2,int(p**0.5)+1))]
    for p in primes: base.append([1.0 if v%p!=0 else 0.0 for v in speeds]); lb.append(1); ub.append(P)
    def feas(r):
        rows=[row[:] for row in base]; l=lb[:]; u=ub[:]
        for uu in range(U):
            rows.append([1.0 if dist[p][uu]<=r else 0.0 for p in range(P)]); l.append(1); u.append(P)
        res=milp(c=np.zeros(P),constraints=LinearConstraint(np.array(rows),l,u),integrality=np.ones(P),bounds=Bounds(0,1),options={"time_limit":20})
        return [speeds[i] for i in range(P) if res.x is not None and res.x[i]>0.5] if (res.success and res.x is not None) else None
    lo,hi,best=0,len(rvals)-1,None
    while lo<=hi:
        mid=(lo+hi)//2; S=feas(rvals[mid])
        if S is not None: best=(rvals[mid],S); hi=mid-1
        else: lo=mid+1
    return best
for n,V in [(10,42),(12,50),(14,58),(15,62)]:
    res=cmin(n,V)
    if res is None: print(f"n={n}: none (V={V})"); continue
    r,S=res; tM,at=Mexact(S); g=0
    for v in S: g=math.gcd(g,v)
    marg=tM-F(1,n)
    print(f"n={n:>2}(V={V}): M_prim={tM}={float(tM):.5f} t*={at} margin={marg}={float(marg):.5f} faray_det={n*tM.numerator-tM.denominator} prim={g==1} S={S}")
