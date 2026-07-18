from fractions import Fraction as F
from math import gcd
def norm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return min(r,1-r)
def minnorm(fam,t): return min(norm(F(v)*t) for v in fam)
def MP(P):
    best=F(0)
    for q in range(2,2*max(P)+2):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            v=minnorm(P,F(a,q))
            if v>best: best=v
    return best
def good_t0(P):
    out=[]
    for q in range(2,2*max(P)+2):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            v=minnorm(P,F(a,q))
            if v>=F(1,14): out.append((F(a,q),v))
    out.sort(key=lambda z:-z[1]); return out
def lift_witness(P,K):
    S=sorted(P+K); pmax=max(P)
    for t0,vP in good_t0(P)[:25]:
        smax=(vP-F(1,14))/pmax
        if smax<=0: continue
        for i in range(-600,601):
            s=smax*F(i,600); t=t0+s
            if minnorm(S,t)>=F(1,14): return (t,minnorm(S,t))
    return None

tests=[
 ("j=7 NEAR-EQUAL, core{1..6}", list(range(1,7)), [300,301,302,303,304,305,306]),
 ("j=8 NEAR-EQUAL, core{1..5}", list(range(1,6)), [300,301,302,303,304,305,306,307]),
 ("j=7 SPREAD, core{1..6}",     list(range(1,7)), [200,260,330,410,500,620,760]),
 ("j=6 SPREAD, core{1..7}",     list(range(1,8)), [200,270,360,470,600,760]),
 ("j=7 mid-spread core{1..6}",  list(range(1,7)), [300,305,311,318,326,335,345]),
]
print("Does the cluster-lift find a witness (=> M(S)>=1/14)?  apex-7 test:")
for name,P,K in tests:
    if len(P)+len(K)!=13: print(f"  {name}: size {len(P)+len(K)} skip"); continue
    w=lift_witness(P,K)
    # also the TRUE M(S) at core scale (lower bound) to compare
    print(f"  {name}: M(P)={float(MP(P)):.4f}  witness={'FOUND '+str(float(w[1]))+'>=1/14' if w else 'NONE (construction fails)'}")
