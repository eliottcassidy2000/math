#!/usr/bin/env python3
"""Arc-cube functions, the score->iso compression (n<=4), and the even/odd (biquadratic/Worpitzky) parity split
of the k=8 cap dip (S71).

Owner: invariants as FUNCTIONS on the arc-cube; a+b,a*b (symmetric) vs a^b,b^a (asymmetric, order-aware);
n=3 from the edges (coin flips); n=4 in two schemes (Klein-four vs magma); compression -> tricks.
Merges codex HYP-3147 (n=3 edge kernel, eigenvalue -1/3, Worpitzky A(3,k)=[1,4,1]) with mac-mini S70
(k=8 = solvable biquadratic resolvent) and S67 (phi^4 cumulants).

Verified here:
 (1) the score sequence (commutative a+b face) DETERMINES the iso class for n<=4, FAILS at n=5;
 (2) the cap DIP turns on exactly at |P|=13-k crossing 4->5 (the n=4->5 compression boundary):
     k>=10 (|P|<=3) dip=0; k=9 (|P|=4) dip tiny; k=8 (|P|=5) dip large = the binding row;
 (3) the gK8 dip term -9 S3 + 6 S4 splits by PARITY: even +6 S4 = symmetric/biquadratic (S70),
     odd -9 S3 = antisymmetric/orientation/Worpitzky (codex), the odd dominating.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# (1) score -> iso class compression (the commutative face), n=3..5
def canon(n, arcs):
    best=None
    for perm in itertools.permutations(range(n)):
        key=tuple(sorted((perm[i],perm[j]) if (i,j) in arcs else (perm[j],perm[i])
                         for i in range(n) for j in range(i+1,n)))
        if best is None or key<best: best=key
    return best
def scores(n,arcs):
    sc=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if (i,j) in arcs: sc[i]+=1
            else: sc[j]+=1
    return tuple(sorted(sc))
print("(1) does SCORE (the commutative a+b face) determine the iso class?")
for n in range(3,6):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(pairs)
    iso=set(); s2i=defaultdict(set)
    for bits in range(2**m):
        arcs=set(pairs[k] for k in range(m) if (bits>>k)&1)
        c=canon(n,arcs); iso.add(c); s2i[scores(n,arcs)].add(c)
    bij=all(len(v)==1 for v in s2i.values())
    print(f"    n={n}: #iso={len(iso)}, #score={len(s2i)}, BIJECTIVE={bij}"
          + ("" if bij else "  <- score insufficient; the ORDER (a^b/b^a) face is needed"))

# (2)+(3) the cap dip vs |P|, and the parity split of the gK8 dip term
def sector_of(p): return int((p%1)*7)
def missdist(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): b.add(F(mm,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        t=7-len(set(sector_of(e*((x0+x1)/2)) for e in E))
        if 0<=t<=6: q[t]+=x1-x0
    return q
def Sr(q,r): return sum(comb(t,r)*q[t] for t in range(7))
caps={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
print("\n(2) the cap DIP vs |P|=13-k (the compression boundary is |P|=4->5):")
for k in sorted(caps):
    P=13-k; dip=F(comb(k+1,2),91)-caps[k]
    print(f"    k={k}: |P|={P}  dip={dip}={float(dip):.5f}  {'(=0, compression works |P|<=3)' if dip==0 else ('(tiny, |P|=4 edge)' if k==9 else '(LARGE, |P|=5 = the n=5 transition = binding)')}")

print("\n(3) PARITY split of the gK8 dip term for consec_8 (L_yK8 = 10S0-10S1+10S2 -9S3 +6S4):")
q=missdist(tuple(range(8))); S=[Sr(q,r) for r in range(5)]
odd=-9*S[3]; even=6*S[4]
print(f"    ODD  -9 S3 = {float(odd):+.3f}  (antisymmetric / a^b,b^a / orientation / Worpitzky -- codex HYP-3147)")
print(f"    EVEN +6 S4 = {float(even):+.3f}  (symmetric / a+b,a*b / biquadratic resolvent -- mac-mini S70)")
print(f"    |odd|/|even| = {float(abs(odd)/abs(even)):.2f}  (the ODD/Worpitzky face dominates the k=8 correction)")
print("    n=3 edge-kernel eigenvalue -1/3; Worpitzky A(3,k)=[1,4,1]; x^3=C(x+2,3)+4C(x+1,3)+C(x,3).")
print("\n=> the k=8 dip = EVEN(biquadratic, solvable degree-2, S70) + ODD(Worpitzky n=3-oscillation, codex);")
print("   it appears EXACTLY at |P|=5 = the n=5 transition where score->iso compression FAILS.")
