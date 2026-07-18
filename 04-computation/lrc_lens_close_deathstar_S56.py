from fractions import Fraction as F
from math import gcd
from itertools import combinations
def M_exact(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
# The composed relation to CLOSE the "core misses 13,14" case:
#   stability: v_max <= max(W)/(13 delta);  covering: v_max>=182  =>  max(W) >= 2366*delta.
#   CLAIM (to test): every non-AP 12-set W (no mult of 13 or 14) has delta=M(W)-1/13 > max(W)/2366.
#   If true, 2366*delta > max(W) contradicts max(W)>=2366*delta  =>  no non-AP core  =>  core is AP.
print("Testing the closing relation: non-AP 12-set (no mult 13,14) => delta > max(W)/2366 ?")
C=2366
viol=[]; tested=0; minmargin=None
# enumerate near-AP cores: {1..12} with k elements replaced by values in [1,30]\{13,14,26,28}
base=list(range(1,13))
pool=[v for v in range(1,31) if v%13 and v%14]
import itertools, random
random.seed(1)
# exhaustive: replace 1 or 2 positions
for krep in (1,2):
    for positions in combinations(range(12),krep):
        for newvals in itertools.product(pool,repeat=krep):
            W=base[:]
            for p,nv in zip(positions,newvals): W[p]=nv
            W=sorted(set(W))
            if len(W)!=12: continue
            if any(v%13==0 or v%14==0 for v in W): continue
            if is_AP(W): continue
            tested+=1
            M=M_exact(W); delta=M-F(1,13)
            if delta<=0: continue
            m=max(W); bound=F(m,C)
            margin=delta-bound
            if minmargin is None or margin<minmargin: minmargin=margin; minW=(W,float(delta),m)
            if delta<=bound: viol.append((W,float(delta),m,float(bound)))
    print(f"  after k={krep} replacements: tested {tested}, violations (delta<=max/2366): {len(viol)}")
print(f"\n  RESULT: violations = {len(viol)}")
if viol:
    for w in viol[:6]: print(f"    VIOLATION: W={w[0]} delta={w[1]:.5f} max={w[2]} bound={w[3]:.5f}")
else:
    print(f"  => the relation HOLDS on all tested non-AP compact cores.")
    print(f"     tightest: W={minW[0]} delta={minW[1]:.5f} max={minW[2]} (max/2366={minW[2]/2366:.5f}), margin={float(minmargin):.5f}")
