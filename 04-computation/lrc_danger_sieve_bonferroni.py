from fractions import Fraction as F
from itertools import combinations
import math
delta=F(1,14)
def danger_measure(T):
    T=list(T)
    if not T: return F(1)
    bps=set([F(0),F(1)])
    for v in T:
        for k in range(0,v+1):
            for off in (delta,1-delta,F(0)):
                t=(k+off)/v
                if 0<=t<=1: bps.add(t)
    bps=sorted(bps);tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0:continue
        mid=(x0+x1)/2
        if all(((v*mid)%1<delta or (v*mid)%1>1-delta) for v in T):tot+=x1-x0
    return tot

def lonely_and_bonferroni(S, maxlevel=None):
    S=list(S); n=len(S)
    if maxlevel is None: maxlevel=n
    # partial sums by level
    level_sum={}
    for k in range(0, maxlevel+1):
        tot=F(0)
        for T in combinations(S,k):
            tot += danger_measure(T)
        level_sum[k]=tot*((-1)**k)
    return level_sum

# Sample clusters of 13 distinct positive ints
clusters = {
 "consec 1..13": list(range(1,14)),
 "consec 2..14": list(range(2,15)),
 "arith 1,3,..,25": list(range(1,26,2)),
 "geom-ish smalls": [1,2,3,4,5,6,7,8,9,10,11,12,13],
 "with 7-multiples":[7,14,21,28,35,42,49,56,63,70,77,84,91],
}
for name,S in clusters.items():
    # exact L via full sieve
    ls = lonely_and_bonferroni(S, maxlevel=13)
    full = sum(ls.values())
    # Bonferroni partials: truncate at level k
    partials=[]
    run=F(0)
    for k in range(0,14):
        run+=ls[k]
        partials.append(run)
    print(f"\n=== {name}: S={S}")
    print(f"  EXACT L(S) = {full} = {float(full):.6f}  (target >0, and L>=1/14? actually L can be tiny)")
    print(f"  level sums (signed): " + " ".join(f"{float(ls[k]):+.4f}" for k in range(0,7)) + " ...")
    print(f"  Bonferroni partial after level k (lower bd at even k, upper at odd k):")
    for k in [1,2,3,4,5,6,13]:
        bound = "LOWER" if k%2==0 else "upper"
        print(f"    S_{k} = {float(partials[k]):+.5f}  ({bound})")
