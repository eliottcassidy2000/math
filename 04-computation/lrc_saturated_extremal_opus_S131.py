from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
def distZ(num,den):
    r=num%den; return F(min(r,den-r),den)
def reach_M(V):
    V=[abs(v) for v in V if v!=0]; n=len(V); dens=set()
    for i in range(n):
        dens.add(2*V[i])
        for j in range(i+1,n):
            dens.add(V[i]+V[j])
            if V[i]!=V[j]: dens.add(abs(V[i]-V[j]))
    best=F(0); bt=None
    for d in dens:
        if d==0: continue
        for m in range(1,d):
            mn=min(distZ(v*m,d) for v in V)
            if mn>best: best=mn; bt=F(m,d)
    return best,bt
def saturated(V): return all(any(v%q==0 for v in V) for q in range(2,15))
thr=F(1,14)

print("=== hunting the EXTREMAL (min-M) saturated 13-family (opus-S131) ===\n")
minM=F(1); worst=None; belowthr=[]
def test(V,name=""):
    global minM,worst
    if len(set(V))!=13 or any(v<=0 for v in V): return
    if reduce(gcd,V)!=1 or not saturated(V): return
    M,bt=reach_M(V)
    if M<thr: belowthr.append((tuple(sorted(V)),M))
    if M<minM: minM=M; worst=(tuple(sorted(V)),M,bt,name)

# 1) {1..13}\{k} + {14*j}
for k in range(1,14):
    for j in [1,2,3]:
        base=[x for x in range(1,14) if x!=k]+[14*j]; test(base,f"AP\{k}+{14*j}")
# 2) {1..13}\{k1,k2} + {14, x}
for k1,k2 in combinations(range(1,14),2):
    rem=[x for x in range(1,14) if x not in (k1,k2)]
    for x in [14,28,15,16,26]:
        if x not in rem: test(rem+[14,x] if x!=14 else rem+[14, k1+13],"")
# 3) exhaustive small: 13-subsets of {1..18} that are saturated (near-AP range)
import itertools
cnt=0
for extra in combinations([x for x in range(1,19)], 13):
    cnt+=1
    if cnt>60000: break
    test(list(extra))
print(f"scanned near-AP + 13-subsets of small ranges")
print(f"MIN M over saturated = {minM} = {float(minM):.6f}   (1/14={1/14:.6f}, 1/13={1/13:.6f}, 1/11={1/11:.6f})")
if worst: print(f"   extremal saturated family: {worst[0]}  M={worst[1]} at t={worst[2]}")
print(f"saturated families with M < 1/14 (REFUTE LRC14): {len(belowthr)}  {belowthr[:4]}")
# what's the min over ALL 13-subsets of {1..14}? (the tightest near-AP saturated)
print("\n=== all saturated 13-subsets of {1..15} (the most-AP-like saturated) ===")
res=[]
for S in combinations(range(1,16),13):
    if saturated(list(S)) and reduce(gcd,S)==1:
        M,bt=reach_M(list(S)); res.append((M,S,bt))
res.sort()
for M,S,bt in res[:8]: print(f"   M={float(M):.5f}={M} t={bt}  {S}")
