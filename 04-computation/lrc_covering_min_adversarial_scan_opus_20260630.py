"""Adversarial: try the DENSEST covering sets to beat 14/183. Covering 13&14 with one element (182=13*14
or 364) frees 12 slots for a small interval. Test these + perturbations."""
import math
from fractions import Fraction
from itertools import combinations
target=Fraction(14,183)
def M_exact(S,Qmax=400):
    best=Fraction(0)
    for q in range(2,Qmax+1):
        bb=0
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
                if m<=bb: break
            if m>bb: bb=m
        v=Fraction(bb,q)
        if v>best: best=v
    return best
def is_cov(S): return all(any(x%qq==0 for x in S) for qq in range(2,15))
# densest covering sets: 12 small + 1 element covering 13&14
adversarial=[
 ("{1..12,182}", list(range(1,13))+[182]),
 ("{1..12,364}", list(range(1,13))+[364]),
 ("{1..12,182*3=546}", list(range(1,13))+[546]),
 ("{1..11,13,28}", list(range(1,12))+[13,28]),   # 13 covers q13, 28 covers q14
 ("{1..11,13,14}", list(range(1,12))+[13,14]),
 ("{1..11,13,42}", list(range(1,12))+[13,42]),
 ("{1..11,13,56}", list(range(1,12))+[13,56]),
 ("{1..11,13,70}", list(range(1,12))+[13,70]),
 ("{1..11,14,26}", list(range(1,12))+[14,26]),   # 14 covers q14, 26 covers q13
 ("{2..13,182}", list(range(2,14))+[182]),
 ("{1..10,12,13,14}", list(range(1,11))+[12,13,14]),
 ("{1..12,26} (no q14?)", list(range(1,13))+[26]),
]
print(f"ADVERSARIAL dense covering sets vs target 14/183={float(target):.5f}:")
best=Fraction(1); bestname=None
for name,S in adversarial:
    if len(set(S))!=13: 
        print(f"   {name:>22}: SKIP (not 13 distinct: {len(set(S))})"); continue
    cov=is_cov(S)
    if not cov:
        print(f"   {name:>22}: NOT covering"); continue
    M=M_exact(S)
    flag=" <<< BEATS 14/183!" if M<target else (" = 14/183" if M==target else "")
    if M<best: best=M; bestname=name
    print(f"   {name:>22}: M={str(M):>8}={float(M):.5f}{flag}")
print(f"\n   densest covering-min found = {best}={float(best):.5f} at {bestname}")
print(f"   => 14/183 {'IS' if best>=target else 'is NOT'} the covering-min among these adversarial dense sets")
# also: does ANY 13-mult+14-mult split beat the single-182?
print("\n   the single-element 182=13*14 (frees 12 slots for {1..12}) is the densest packing => hardest set.")
