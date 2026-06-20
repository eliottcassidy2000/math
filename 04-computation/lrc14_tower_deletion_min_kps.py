import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def lonely_measure(C, theta=Fraction(1,14)):
    segs=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0,d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for sh in (-1,0,1):
                a=max(lo+sh,Fraction(0)); b=min(hi+sh,Fraction(1))
                if a<b: segs.append((a,b))
    segs.sort(); union=[]; cur=Fraction(-1)
    for a,b in segs:
        if a>cur: union.append([a,b]); cur=b
        elif b>cur: union[-1][1]=b; cur=b
    return Fraction(1)-sum(b-a for a,b in union)
thr2=Fraction(426,35035)
print(f"thr2 = 426/35035 = {float(thr2):.6f}")
print("=== MIN meas(G_C) over AP-tail 12-cores MISSING each tower bit v in {1,2,4,8} ===")
print("    (if all mins >= thr2, the tower-deletion bound holds; find the binding case)\n")
for v in [1,2,4,8]:
    best=(Fraction(2),None)
    # AP-tail: remove holes (incl v), add tails, 12-core. scan 1-tail and 2-tail.
    for ktail,TMAX in [(1,50),(2,30)]:
        nh=ktail+1
        for holes in itertools.combinations(range(1,14), nh):
            if v not in holes: continue
            for tails in itertools.combinations(range(14,TMAX+1), ktail):
                C=tuple(sorted([d for d in range(1,14) if d not in holes]+list(tails)))
                if len(C)!=13-nh+ktail or reduce(gcd,C)!=1: continue
                L=lonely_measure(C)
                if L<best[0]: best=(L,(holes,tails))
    L,wit=best
    print(f"  v={v} missing: min meas = {L} = {float(L):.6f}  {'>=thr2 OK' if L>=thr2 else '*** < thr2 ***'}  margin={float(L-thr2):.6f}  at {wit}")
print("\n=== for comparison: shell-1-FULL cores that ACHIEVE the threshold ===")
for name,holes,tails in [("drop-12",{12},()),("drop-6",{6},())]:
    C=tuple(sorted([d for d in range(1,14) if d not in holes]+list(tails)))
    print(f"  {name}: meas={lonely_measure(C)}  (shell-1 full, 1,2,4,8 all present)")
