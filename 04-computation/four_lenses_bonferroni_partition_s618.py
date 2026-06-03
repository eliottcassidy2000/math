from fractions import Fraction as F
from itertools import combinations

def forbidden_intervals(v, delta):
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
    cur=lst[0]
    for x in lst[1:]: cur=inter(cur,x)
    return cur

def Sk(V):
    n=len(V); delta=F(1,n+1); A=[forbidden_intervals(v,delta) for v in V]
    S=[F(0)]*(n+1); S[0]=F(1)
    for k in range(1,n+1):
        for I in combinations(range(n),k):
            S[k]+=meas(interN([A[i] for i in I]))
    U=A[0]
    for a in A[1:]: U=merge(U+a)
    return S, F(1)-meas(U)

print("=== HELLY NUMBER 3: order-3 Bonferroni lower bound T3 = 1-S1+S2-S3 decides loneliness ===")
print("  (odd-order Bonferroni is a rigorous LOWER bound on p0; T3>0 PROVES p0>0)")
print(f"  {'config':20} {'p0':>8} {'T2(upper)':>10} {'T3(lower)':>10} {'verdict@order3':>16}")
for V,nm in [([1,2,3],"AP{1,2,3}"),([1,2,3,4],"AP{1,2,3,4}"),([1,3,4,7],"chain{1,3,4,7}"),
             ([1,3,4,5,9],"chain{1,3,4,5,9}"),([1,4,6,9],"rand{1,4,6,9}"),([2,3,7,8],"rand{2,3,7,8}"),
             ([1,5,8,11,13],"rand{1,5,8,11,13}"),([2,5,9,11],"rand{2,5,9,11}")]:
    S,p0=Sk(V)
    T2=S[0]-S[1]+S[2]; T3=T2-S[3]
    verdict = "LONELY (decided)" if T3>0 else ("COLLAPSE" if p0==0 else "undecided@3")
    print(f"  {nm:20} {float(p0):>8.4f} {float(T2):>10.3f} {float(T3):>10.3f} {verdict:>16}")

print("\n=== VITALI WALL: at collapse the Bonferroni truncations never close until the TOP order ===")
for V,nm in [([1,3,4,7],"chain{1,3,4,7}"),([1,4,6,9],"rand{1,4,6,9}")]:
    S,p0=Sk(V); n=len(V)
    T=[sum((-1)**k*S[k] for k in range(m+1)) for m in range(n+1)]
    print(f"  {nm}: p0={float(p0):.4f}  T_m="+", ".join(f"{float(T[m]):+.3f}" for m in range(n+1))
          +f"   (first finite order m<n reaching p0: {'NONE -> Vitali' if p0==0 else 'n/a'})")

print("\n=== PARTITION-FUNCTION SIBLING: depth GF P(z)=sum p_k z^k, and dependency-graph independence polynomial ===")
def depth_dist(V):
    n=len(V); delta=F(1,n+1)
    ev={}
    for v in V:
        for (lo,hi) in forbidden_intervals(v,delta):
            ev[lo]=ev.get(lo,0)+1; ev[hi]=ev.get(hi,0)-1
    pts=sorted(set(list(ev)+[F(0),F(1)])); pk={}; d=0; prev=F(0)
    for p in pts:
        if p>prev: pk[d]=pk.get(d,F(0))+(p-prev)
        d+=ev.get(p,0); prev=p
    return pk
def indep_poly(V):
    n=len(V); delta=F(1,n+1); A=[forbidden_intervals(v,delta) for v in V]
    indep_pair=(2*delta)**2; E=set()
    for i,j in combinations(range(n),2):
        if meas(inter(A[i],A[j]))!=indep_pair: E.add((i,j))
    coef=[0]*(n+1)
    for r in range(n+1):
        for Sset in combinations(range(n),r):
            if all((a,b) not in E for a,b in combinations(Sset,2)): coef[r]+=1
    # smallest-magnitude negative real root of sum coef[k] lam^k via bisection on I(lam), lam<0
    def I(x): return sum(coef[k]*x**k for k in range(len(coef)))
    root=None; x=0.0
    while x>-50:
        if I(x)*I(x-0.001)<0: 
            a,b=x-0.001,x
            for _ in range(60):
                m=(a+b)/2
                if I(a)*I(m)<=0: b=m
                else: a=m
            root=(a+b)/2; break
        x-=0.001
    return coef,E,root
for V,nm in [([1,2,3,4],"AP{1,2,3,4}"),([1,3,4,7],"chain{1,3,4,7}"),([1,4,6,9],"rand{1,4,6,9}"),([2,3,7,8],"rand{2,3,7,8}")]:
    pk=depth_dist(V); coef,E,root=indep_poly(V)
    P=" ".join(f"p{k}={float(pk.get(k,0)):.3f}" for k in range(len(V)+1))
    lamstar = None if root is None else round(abs(root),3)
    print(f"  {nm:16} depthGF coeffs[{P}] p0={float(pk.get(0,0)):.4f} | dep-edges={len(E)} indep-poly={coef} Shearer lambda*={lamstar}")
