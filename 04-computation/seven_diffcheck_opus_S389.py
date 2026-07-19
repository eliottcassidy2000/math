from fractions import Fraction as F
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
random.seed(389)
found=[]
for _ in range(25):
    V=sorted(random.sample(range(1,60),13))
    sums={a+b for a,b in combinations(V,2)}
    diffs={abs(a-b) for a,b in combinations(V,2) if a!=b}
    peaks={2*v for v in V}
    cand=sorted((sums|diffs|peaks)-{0})
    best=F(0); bq=None; bp=None
    for q in cand:
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,bq,bp=g,q,p
    if bq not in sums:
        found.append((V,bq,bp,best,bq in diffs,bq in peaks))
print(f"difference-only optima found: {len(found)}")
for V,q,p,g,isd,isp in found[:2]:
    print()
    print(f"  V = {V}")
    print(f"  optimum at t = {p}/{q}, gap = {g} = {float(g):.8f}")
    print(f"  q={q} is a difference: {isd}; a half-period 2v: {isp}; a pair SUM: {q in {a+b for a,b in combinations(V,2)}}")
    # exhaustive confirmation: is the max over ALL pair-sum denominators smaller?
    sums=sorted({a+b for a,b in combinations(V,2)})
    bs=F(0); bsq=None
    for qq in sums:
        for pp in range(1,qq):
            gg=gap_at(V,F(pp,qq))
            if gg>bs: bs,bsq=gg,qq
    print(f"  best achievable over PAIR-SUM denominators only: {bs} = {float(bs):.8f} (at q={bsq})")
    print(f"  => pair-sums alone {'MISS the optimum' if bs<g else 'attain it'}"
          f"  (deficit {float(g-bs):.8f})")
    print(f"  which pairs give the difference q={q}: "
          f"{[(a,b) for a,b in combinations(V,2) if abs(a-b)==q][:3]}")
