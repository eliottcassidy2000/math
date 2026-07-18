import sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
n=14
def is_tight(fam):
    Q=2*max(fam)+2; hit=False
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if 14*m>q: return False
            if 14*m==q: hit=True
    return hit
def prim(fam): return reduce(gcd,fam)==1

# 1) Is bounded-spread even TRUE? Test candidate LARGE primitive families for tightness.
print("(1) Do large primitive tight families exist? (M=1/14 exactly?)")
cands={
 "translated AP {15..27}":list(range(15,28)),
 "translated AP {29..41}":list(range(29,42)),
 "AP*? {14*2+1..14*2+13}={29..41}":list(range(29,42)),
 "GW+14 {15..25,27,38}":[15,16,17,18,19,20,21,22,23,24,25,27,38],
 "scaled-ish {2,4,..,26}(non-prim)":[2*k for k in range(1,14)],
 "GW {1..11,13,24}":[1,2,3,4,5,6,7,8,9,10,11,13,24],
}
for nm,fam in cands.items():
    fam=sorted(set(fam))
    if len(fam)!=13: print(f"   {nm}: (not 13 distinct) skip"); continue
    print(f"   {nm}: primitive={prim(fam)} tight={is_tight(fam)}"); sys.stdout.flush()

# 2) Exhaustive: ALL primitive tight families with Vmax<=20 (expect only AP)
print("\n(2) EXHAUSTIVE primitive tight families, Vmax<=20:")
found=[]
V0=20
for w in range(13,V0+1):          # largest element = w
    for rest in combinations(range(1,w),12):
        fam=list(rest)+[w]
        if reduce(gcd,fam)!=1: continue
        if is_tight(fam): found.append(tuple(sorted(fam)))
    print(f"   ...Vmax={w} done, cumulative tight found={len(found)}"); sys.stdout.flush()
print(f"   TIGHT primitive families with Vmax<=20: {found}")
