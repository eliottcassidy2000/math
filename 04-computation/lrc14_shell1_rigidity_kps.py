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
def shell1_intact(C):  # are speeds 1,2,4,8 all present?
    return all(v in set(C) for v in (1,2,4,8))
thr2=Fraction(426,35035)
SHELL1={1,2,4,8}
print(f"=== RIGIDITY: meas(G_C) < 426/35035  ==>  shell-1 tower {{1,2,4,8}} intact ? ===")
print(f"    (i.e. no sub-threshold AP-tail row may delete any of 1,2,4,8)\n")
viol=[]; subcnt=0; tot=0
# k-tail: remove (k+1) holes from {1..13}, add k tails in [14, TMAX]
for ktail, TMAX in [(1,60),(2,40),(3,30)]:
    nh=ktail+1
    rowcnt=0; below=0; below_shell1_damaged=0
    for holes in itertools.combinations(range(1,14), nh):
        for tails in itertools.combinations(range(14,TMAX+1), ktail):
            C=tuple(sorted([d for d in range(1,14) if d not in holes]+list(tails)))
            if len(C)!=(13-nh+ktail): continue
            if reduce(gcd,C)!=1: continue
            rowcnt+=1; tot+=1
            L=lonely_measure(C)
            if L<thr2:
                below+=1; subcnt+=1
                if not shell1_intact(C):
                    below_shell1_damaged+=1
                    viol.append((float(L),holes,tails,sorted(SHELL1-set(C))))
    print(f"  {ktail}-tail (holes={nh}, tails<= {TMAX}): {rowcnt} rows, {below} below thr2, "
          f"{below_shell1_damaged} of those with shell-1 DAMAGED")
print(f"\n  TOTAL sub-threshold rows: {subcnt}; with shell-1 damaged (VIOLATIONS): {len(viol)}")
if viol:
    print("  *** RIGIDITY FAILS — counterexamples (sub-thr2 but shell-1 damaged): ***")
    for L,h,t,miss in sorted(viol)[:10]: print(f"      meas={L:.6f} holes={h} tails={t} missing-from-shell1={miss}")
else:
    print("  ==> RIGIDITY HOLDS: every sub-426/35035 AP-tail row keeps the full shell-1 tower {1,2,4,8}.")
    print("      The mouth (drop-6 collar) is owned by the dyadic-1 tower; you cannot go below the AP")
    print("      2nd value without the complete {1,2,4,8} chain. This is the carry conservation law.")
