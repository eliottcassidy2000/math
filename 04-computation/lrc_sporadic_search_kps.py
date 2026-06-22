"""
Map the sporadic tiler family (kps-S31n): search AP single-replacements for tightness.
For AP {1..13}, remove one speed, add a larger integer v; tight (M=1/14) => sporadic tiler.
Tests finiteness of OPEN-Q-108 sporadic residue and the gap+collision mechanism.
"""
from fractions import Fraction as F
def nf(x):
    r=x%1; return min(r,1-r)
def M_exact(S):
    # max-min over rational candidates: breakpoints (14m+-1)/(14s) and j/14d
    S=[s for s in S if s!=0]; cset=set()
    for s in S:
        a=abs(s)
        for m in range(0,a+1):
            for sg in (-1,1):
                t=F(14*m+sg,14*a)
                if 0<t<1: cset.add(t)
    for d in range(1,29):
        for j in range(1,d): cset.add(F(j,d))
    best=F(0); arg=None
    for t in cset:
        mn=min(nf(s*t) for s in S)
        if mn>best: best=mn; arg=t
    return best,arg
AP=list(range(1,14))
tight=[]
for rem in AP:
    for v in range(14,61):
        if v in AP: continue
        S=sorted([x for x in AP if x!=rem]+[v])
        if len(set(S))<13: continue
        M,arg=M_exact(S)
        if M==F(1,14):
            tight.append((rem,v,arg,S))
print(f"AP single-replacements (remove rem, add v<=60) that stay TIGHT (M=1/14 exactly):")
print(f"  found {len(tight)} sporadic tilers:")
for rem,v,arg,S in tight:
    # residue collision check
    r_v = v % 14; collides = [s for s in S if s!=v and s%14==r_v]
    print(f"   remove {rem:2d}, add {v:2d} (v={r_v} mod14, collides w/ {collides}); witness t={arg}; S={S}")
print("\n=> If ALL are v ≡ (some speed) mod 14 with a single gap+collision, the mechanism (THM-560c)")
print("   is confirmed; the COUNT bounds the single-replacement sporadic family.")
