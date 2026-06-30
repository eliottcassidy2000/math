"""
TASK 1 - linearize at the cusp (razor): M(S) for single-element perturbations of the AP {1..13} (the
M-minimum). Which perturbations are FLAT (keep M=1/14, the razor manifold) vs RELEVANT (M increases)?
Classify odd (cusp-changing) vs even (deeper-shell) perturbations.
"""
import math
from fractions import Fraction
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
                if m<=bb:break
            if m>bb:bb=m
        v=Fraction(bb,q)
        if v>best:best=v
    return best
AP=list(range(1,14)); base=M_exact(AP)
print(f"M(AP {{1..13}}) = {base} = {float(base):.5f}  (the razor / global M-minimum)")
print("\nSingle-element replacements e -> e' (e in 1..13, e' in 14..27): M and FLAT(=1/14)/RELEVANT(>1/14)?")
flat=[]; rel=[]
for e in range(1,14):
    for ep in range(14,28):
        S=sorted(set(AP)-{e}|{ep})
        if len(S)!=13: continue
        M=M_exact(S)
        rec=(e,ep,M)
        (flat if M==Fraction(1,14) else rel).append(rec)
print(f"  FLAT perturbations (M stays 1/14): {len(flat)} of {len(flat)+len(rel)}")
print(f"  RELEVANT (M>1/14): {len(rel)}")
# which e (removed) keeps flat? odd vs even
from collections import Counter
flat_by_e=Counter(e for e,_,_ in flat); rel_by_e=Counter(e for e,_,_ in rel)
print("  removed element e -> # flat / # relevant (e odd changes the Z_7 cusp core):")
for e in range(1,14):
    par="odd " if e%2 else "even"
    print(f"     e={e:>2} ({par}): flat {flat_by_e.get(e,0):>2}, relevant {rel_by_e.get(e,0):>2}")
# the flat ones: e' mod 14 != 0 keeps witness t=1/14 valid (gap >=1/14)
print("  => FLAT iff the new element e' has e' mod 14 != 0 AND removing e doesn't break the comb:")
flatset=set((e%2,ep%14==0) for e,ep,_ in flat)
print(f"     flat perturbations all have e' mod 14 != 0: {all(ep%14!=0 for _,ep,_ in flat)}")
print(f"     min M among RELEVANT = {min((M for _,_,M in rel), default=None)} (the steepest still >= 1/14)")
