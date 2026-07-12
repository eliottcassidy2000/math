#!/usr/bin/env python3
"""Slice 2: the gcd-7 / mod-7 reduction. For E = {c + 7 m_i} (mod-7-aligned), relate
its sector-emptiness distribution N to structural quantities. HYPOTHESIS: mod-7 families
have HIGH mu (mean empty count), dispatching them from the J-frontier -- make the
mechanism exact. Also test the general gcd-g family E = g*E' vs E' (dilation invariance
of N-distribution, already THM known) and the mod-7 SPECIAL case."""
from fractions import Fraction as F
def pdist(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return p
def stats(E):
    p=pdist(E); mu=sum(j*p[j] for j in range(8)); m2=sum(j*(j-1)*p[j] for j in range(8))
    return mu, 6*mu-m2  # mu, J
print("MOD-7-ALIGNED families E={c+7m_i} vs threshold: mu high => J high (dispatched)")
print(f"{'family':32s} {'mu':>7s} {'J':>7s}  vs Jmin=5.06, thr=4.75")
for c in [1,2,3]:
    for span in [[0,1,2,3,4,5,6,7,8]]:
        E=[c+7*m for m in span]
        mu,J=stats(E)
        print(f"  c={c} {{{c}+7m: m in 0..8}}={E[:3]}...: mu={float(mu):7.4f} J={float(J):7.4f}  {'DISPATCHED (J>thr)' if J>F(432,91) else 'CHECK'}")
print()
# WHY mu is high: a mod-7 family at x=a/7 hits FEW sectors (aliasing) => many empty => high N
print("mod-7 mechanism: at x near j/7, all e_i=c+7m_i have frac(e_i x)~ frac(c j/7) SAME sector")
print("  => k phases collapse to ~1 sector => 6 empty => N=6 spike => high mu. Check p[6],p[7]:")
for c in [1,2]:
    E=[c+7*m for m in range(9)]; p=pdist(E)
    print(f"  c={c}: p[5]={float(p[5]):.4f} p[6]={float(p[6]):.4f} p[7]={float(p[7]):.4f}  (consec {{1..9}} p[6]={float(pdist(list(range(1,10)))[6]):.4f})")
print()
# the general gcd reduction: E = g*E' has SAME N-distribution as E' (dilation-invariance)
print("gcd-g dilation invariance (N-dist of g*E' = N-dist of E'), exact check:")
for E1 in [[1,2,3,4,5,6,7,8,9],[1,2,4,7,11,16,22,29,37]]:
    for g in [2,3,5]:
        p1=pdist(E1); pg=pdist([g*e for e in E1])
        same=all(p1[j]==pg[j] for j in range(8))
        print(f"  E'={E1[:4]}... g={g}: N-dist identical: {same}")
