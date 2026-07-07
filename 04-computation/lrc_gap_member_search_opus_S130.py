from fractions import Fraction as F
from math import gcd
from itertools import combinations, product
from functools import reduce

def dist_int(num, den):
    r = num % den
    return F(min(r, den - r), den)

def reach_M(V):
    """Exact max-min reach M(V)=max_t min_i ||v_i t|| for integer speeds (t rational,
    denominators among 2v_i, v_i+v_j, |v_i-v_j| -- peaks & crossings)."""
    V = [abs(v) for v in V if v != 0]
    dens = set()
    n = len(V)
    for i in range(n):
        dens.add(2*V[i])
        for j in range(i+1, n):
            dens.add(V[i]+V[j])
            if V[i]!=V[j]: dens.add(abs(V[i]-V[j]))
    best = F(0)
    for d in dens:
        if d==0: continue
        for m in range(1, d):
            mn = None
            for v in V:
                dd = dist_int(v*m, d)
                if mn is None or dd < mn: mn = dd
                if mn == 0: break
            if mn is not None and mn > best: best = mn
    return best

lo, hi = F(1,13), F(2,25)
def in_gap(M): return lo < M < hi

AP=list(range(1,13))
print(f"M(AP)={reach_M(AP)} (expect 1/13)   gap=(1/13,2/25)=({float(lo):.6f},{float(hi):.6f})")
print("SEARCH: any 12-family with M in open gap (1/13,2/25) => REFUTES (C)\n")

found=[]; tested=0; near=[]  # near = M within [1/13, 1/13+0.004) i.e. just above AP
def consider(V):
    global tested
    if reduce(gcd,V)!=1: return
    M=reach_M(V); tested+=1
    if in_gap(M): found.append((tuple(V),M))
    elif lo <= M < lo + F(1,250): near.append((tuple(V),M))

# 1) {1..11, x}
for x in range(13,120): consider(list(range(1,12))+[x])
# 2) {1..10, x, y}
for x in range(11,45):
    for y in range(x+1,55): consider(list(range(1,11))+[x,y])
# 3) two-entry near-AP perturbations: replace entries i,j of AP by AP value + shift
for i,j in combinations(range(12),2):
    for si,sj in product([12,13,24,25,-1,1,11,14],[12,13,24,25,-1,1,11,14]):
        V=AP.copy(); V[i]=AP[i]+si; V[j]=AP[j]+sj
        if any(v<=0 for v in V) or len(set(map(abs,V)))<12: continue
        consider(V)
# 4) small dilations of AP with one defect: 2*AP with a swap, etc.
for base in ([2*k for k in range(1,13)], [3*k for k in range(1,13)]):
    for i in range(12):
        for d in [-1,1,-2,2,base[0]//2 if base[0]%2==0 else 1]:
            V=base.copy(); V[i]=base[i]+d
            if any(v<=0 for v in V) or len(set(map(abs,V)))<12: continue
            consider(V)

print(f"tested: {tested}")
print(f"GAP MEMBERS (refuting (C)): {len(found)}")
for V,M in found[:20]: print("   ", V, "M=",M, float(M))
print(f"\nnearest-above-AP families (M in [1/13, 1/13+0.004)): {len(near)} -- these hug the AP")
for V,M in sorted(near, key=lambda x:x[1])[:12]:
    print(f"    M={float(M):.6f} ({M})  {V}")
