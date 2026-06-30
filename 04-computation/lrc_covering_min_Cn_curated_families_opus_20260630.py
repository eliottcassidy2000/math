"""Wide-n picture: evaluate curated covering families across n, find the winner (C(n)-estimate) per n."""
import math
from fractions import Fraction
def M_exact(S,Qmax):
    best=Fraction(0)
    for q in range(2,Qmax+1):
        bb=0
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
                if m<=bb:break
            if m>bb:bb=m
        v=Fraction(bb,q)
        if v>best:best=v
    return best
def is_cov(S,n): return len(S)==n-1 and all(any(x%q==0 for x in S) for q in range(2,n+1))
def best_family(n):
    Qmax=(n*n-n+1)+12
    cands={}
    # construction: {1..n-2,(n-1)n}
    S=list(range(1,n-1))+[(n-1)*n]
    if is_cov(S,n): cands["construction"]=M_exact(S,Qmax)
    # drop-2: {1,3..n-1} + best large (mult of n or lcm)
    base=[1]+list(range(3,n))  # skip 2; covers up to n-1
    bestd2=Fraction(1)
    for L in [n,2*n,3*n,4*n,5*n,6*n,(n-1)*n,n*(n+1),2*n*n//3 if (2*n*n)%3==0 else n]:
        S=sorted(set(base+[L]))
        if is_cov(S,n):
            m=M_exact(S,Qmax)
            if m<bestd2: bestd2=m
    if bestd2<Fraction(1): cands["drop-2"]=bestd2
    # drop-2,3: {1,4,5..n-1}+2 larges
    base=[1]+list(range(4,n))
    best23=Fraction(1)
    for L1 in [n,2*n,3*n]:
        for L2 in [(n-1)*n,2*(n-1),3*(n-1),6*(n-1)]:
            S=sorted(set(base+[L1,L2]))
            if is_cov(S,n):
                m=M_exact(S,Qmax)
                if m<best23: best23=m
    if best23<Fraction(1): cands["drop-2,3"]=best23
    # interval {1..n-1} won't cover n; {1..n-2,n} won't cover n-1; skip
    if not cands: return None,None,None
    win=min(cands,key=lambda k:cands[k])
    return cands[win],win,cands
print(f"{'n':>3} {'C(n)~':>9} {'dec':>8} {'1/n':>8} {'win family':>14}  candidates")
prev=None
for n in range(4,16):
    C,win,cands=best_family(n)
    if C is None: print(f"{n:>3} (no candidates)"); continue
    cs=" ".join(f"{k}={float(v):.4f}" for k,v in sorted(cands.items()))
    trans=" <== TRANSITION" if prev and prev!=win else ""
    print(f"{n:>3} {str(C):>9} {float(C):>8.5f} {1/n:>8.5f} {win:>14}  [{cs}]{trans}")
    prev=win
