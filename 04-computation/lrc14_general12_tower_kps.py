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
thr2=Fraction(426,35035); TOWER={1,2,4,8}
print(f"=== GENERAL primitive positive 12-cores in [1,B]: meas<426/35035 ==> {{1,2,4,8}} subset C ? ===")
print(f"    (the OPEN-Q-108 generalization beyond AP-tail; HYP-2651 atlas: drop-6 unique min over [1,19])\n")
for B in [16,18,19]:
    sub=0; viol=[]; cnt=0; minrow=(Fraction(2),None)
    for C in itertools.combinations(range(1,B+1),12):
        if reduce(gcd,C)!=1: continue
        cnt+=1
        L=lonely_measure(C)
        if L<minrow[0]: minrow=(L,C)
        if L<thr2:
            sub+=1
            if not TOWER<=set(C): viol.append((float(L),C,sorted(TOWER-set(C))))
    print(f"  B={B}: {cnt} primitive 12-cores; {sub} below thr2; tower-violations(sub-thr2 w/o full tower)={len(viol)}")
    print(f"      global min = {minrow[0]} at {minrow[1]} (tower-full? {TOWER<=set(minrow[1])})")
    if viol:
        print("      *** COUNTEREXAMPLES (sub-thr2, tower NOT full): ***")
        for L,C,miss in sorted(viol)[:8]: print(f"         meas={L:.6f} C={C} missing={miss}")
    else:
        print("      ==> every sub-426/35035 general 12-core in this box contains the full tower {1,2,4,8}.")
