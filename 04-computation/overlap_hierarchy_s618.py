from fractions import Fraction as F
from itertools import combinations
import math

def forbidden_intervals(v, delta):
    # union of arcs {t: ||v t||<delta}, as sorted disjoint [lo,hi) on [0,1)
    r=delta/v; out=[]
    for k in range(v):
        c=F(k,v); L=(c-r)%1; H=(c+r)%1
        if L<H: out.append((L,H))
        else: out+=[(L,F(1)),(F(0),H)]
    return merge(out)
def merge(iv):
    iv=sorted(iv); res=[]
    for l,h in iv:
        if res and l<=res[-1][1]: res[-1]=(res[-1][0],max(res[-1][1],h))
        else: res.append((l,h))
    return res
def inter(a,b):
    res=[]; i=j=0
    while i<len(a) and j<len(b):
        lo=max(a[i][0],b[j][0]); hi=min(a[i][1],b[j][1])
        if lo<hi: res.append((lo,hi))
        if a[i][1]<b[j][1]: i+=1
        else: j+=1
    return res
def meas(iv): return sum(h-l for l,h in iv)
def interN(lst):
    if not lst: return [(F(0),F(1))]
    cur=lst[0]
    for x in lst[1:]: cur=inter(cur,x)
    return cur

def analyze(V, name):
    n=len(V); delta=F(1,n+1)
    A=[forbidden_intervals(v,delta) for v in V]
    # S_k = sum over |I|=k of meas(intersection of A_i, i in I)
    S=[F(0)]*(n+1); S[0]=F(1)
    for k in range(1,n+1):
        for I in combinations(range(n),k):
            S[k]+=meas(interN([A[i] for i in I]))
    # p0 by inclusion-exclusion
    p0_ie=sum((-1)**k*S[k] for k in range(n+1))
    # p0 direct: measure of complement of union
    U=A[0]
    for a in A[1:]: U=merge(U+a)
    p0_dir=F(1)-meas(U)
    # Bonferroni partial sums T_m = sum_{k<=m}(-1)^k S_k ; brackets p0 (odd m lower, even m upper)
    T=[sum((-1)**k*S[k] for k in range(m+1)) for m in range(n+1)]
    # order where bracket closes to p0
    print(f"\n{name}  n={n} delta={delta}  p0={float(p0_dir):.4f} {'COLLAPSE' if p0_dir==0 else ''}")
    print(f"  S_k overlap hierarchy: "+", ".join(f"S{k}={float(S[k]):.3f}" for k in range(n+1)))
    print(f"  incl-excl p0 = {float(p0_ie):.4f}  (matches direct: {p0_ie==p0_dir})")
    print(f"  Bonferroni T_m: "+", ".join(f"T{m}={float(T[m]):+.3f}" for m in range(n+1)))
    # Helly: triple overlaps S3 and which triples are sum-relations a+b=c (up to the arc structure)
    triples_pos=[]; sumrel=[]
    for I in combinations(range(n),3):
        m=meas(interN([A[i] for i in I]))
        a,b,c=sorted(V[i] for i in I)
        if m>0: triples_pos.append((a,b,c,float(m)))
        if a+b==c: sumrel.append((a,b,c))
    print(f"  Helly-3: triples with positive triple-overlap: {[(t[0],t[1],t[2]) for t in triples_pos]}")
    print(f"           sum-relations a+b=c among speeds: {sumrel}")
    return S, p0_dir

print("="*70); print("OVERLAP HIERARCHY / inclusion-exclusion / Helly-3 = sum-relations")
for V,nm in [([1,2,3],"AP{1,2,3}"),([1,2,3,4],"AP{1,2,3,4}"),([1,3,4,7],"chain{1,3,4,7}"),
             ([1,3,4,5,9],"chain{1,3,4,5,9}"),([1,4,6,9],"random{1,4,6,9}"),([2,3,7,8],"random{2,3,7,8}")]:
    analyze(V,nm)

# ---- Lens 4: independence polynomial of the arc-correlation dependency graph ----
print("\n"+"="*70); print("PARTITION-FUNCTION SIBLING: independence polynomial of the arc-dependency graph")
def indep_poly_smallest_root(V):
    n=len(V); delta=F(1,n+1); A=[forbidden_intervals(v,delta) for v in V]
    # dependency edge i~j iff arcs CORRELATED: meas(A_i cap A_j) != (2delta)^2
    indep_pair=(2*delta)**2
    E=set()
    for i,j in combinations(range(n),2):
        if meas(inter(A[i],A[j]))!=indep_pair: E.add((i,j))
    # independence polynomial I(lambda)=sum over independent sets lambda^|S|
    coef=[0]*(n+1)
    for r in range(n+1):
        for S in combinations(range(n),r):
            if all((min(a,b),max(a,b)) not in E for a,b in combinations(S,2)):
                coef[r]+=1
    # smallest positive real root of I(lambda)=0 (Shearer/LLL critical activity), search negative axis
    import numpy as np
    roots=np.roots(coef[::-1]) if any(coef[1:]) else []
    realneg=[r.real for r in roots if abs(r.imag)<1e-9 and r.real<0]
    crit=max(realneg) if realneg else None   # closest-to-0 negative root = LLL threshold -lambda*
    return coef, (abs(crit) if crit is not None else None), len(E)
for V,nm in [([1,2,3,4],"AP{1,2,3,4}"),([1,3,4,7],"chain{1,3,4,7}"),([1,4,6,9],"random{1,4,6,9}"),([2,3,7,8],"random{2,3,7,8}")]:
    coef,thr,ne=indep_poly_smallest_root(V)
    print(f"  {nm:18} dep-edges={ne}  indep-poly coeffs={coef}  LLL critical activity lambda*={thr if thr is None else round(thr,3)}")
