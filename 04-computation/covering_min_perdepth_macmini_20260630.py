import functools, math
from fractions import Fraction as F
from scipy.optimize import milp, LinearConstraint, Bounds
import numpy as np
print=functools.partial(print,flush=True)
def norm(x):
    f=x-int(x); f=f+1 if f<0 else f
    return min(f,1-f)
def min_M_at_V(n,V,tl=8):
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
        res=milp(c=np.zeros(P),constraints=LinearConstraint(np.array(rows),l,u),integrality=np.ones(P),bounds=Bounds(0,1),options={'time_limit':tl})
        return res.success and res.x is not None
    lo,hi,best=0,len(rvals)-1,None
    while lo<=hi:
        mid=(lo+hi)//2
        if feas(rvals[mid]): best=rvals[mid]; hi=mid-1
        else: lo=mid+1
    return best
for n in [12,13,14]:
    nb=F(1,n-1); constr=F(n,n*n-n+1)
    print(f'n={n}: near-block 1/(n-1)={float(nb):.5f}, construction n/Phi6={float(constr):.5f}')
    found_a=None
    for a in range(2,8):
        m=(n-1)*a+1; target=F(a,m)
        mm=min_M_at_V(n,m)
        ach = (mm==target)
        if ach and found_a is None: found_a=a
        print(f'   a={a} (m={m}, target {target}={float(target):.5f}): min-M@V={m} = {mm}={float(mm):.5f}  depth-{a}-achievable:{ach}  <construction:{mm<constr}')
    print(f'  => smallest achievable depth a({n}) = {found_a}; spread beats construction: {found_a is not None and found_a<n}')
