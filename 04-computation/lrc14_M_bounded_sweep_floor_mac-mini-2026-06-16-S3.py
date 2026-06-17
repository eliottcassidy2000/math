#!/usr/bin/env python3
"""
LRC(14) PROVE route, part 2: BOUNDED-MAX EXHAUSTIVE SWEEP + floor structure.

From part 1:
  * M is SCALE-INVARIANT: M(cS)=M(S) (verified exactly).  Proof sketch below.
  * Empirically M(S) >= 1/(2 max(S)) (the slowest-runner-at-half floor),
    with the EXTREMAL ratio M*2max = 13/7 hit by the tight AP {1..13}.

This script:
  1. Proves/illustrates the scale-invariance identity precisely.
  2. Tests the floor M(S) >= 1/(2 max(S)) much harder, and asks the SHARPER
     question: what is inf over primitive S of M(S)*max(S)?  (a scale-invariant
     quantity, since M*max is invariant under S->cS).  If that inf is >= 1/14 * c
     for the relevant max, we'd be done -- check whether M*max has a positive inf.
  3. Does an EXHAUSTIVE sweep of ALL primitive 13-subsets of {1..W} for small W,
     reporting the minimum M and any M<1/14.  This is rigorous for bounded max.
  4. Connects to kind-pasteur's L-quantization to bound the lcm of a minimizer.

Disprove-seed: a counterexample is M=a/d<1/14; we enumerate which a/d are
reachable as envelope-vertex VALUES and check the exhaustive bounded space.
"""

from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at
def lcm(a,b): return a*b//gcd(a,b)
def lcml(xs): return reduce(lcm, xs, 1)
def prim(S): return reduce(gcd, S, 0) == 1

# fast float M for screening
def M_fast(S):
    Ss = sorted(set(S))
    cs = set()
    for v in Ss:
        k=0
        while 2*k+1 <= v: cs.add((2*k+1)/(2.0*v)); k+=1
    n=len(Ss)
    for i in range(n):
        for j in range(i+1,n):
            for d in (Ss[i]+Ss[j], Ss[j]-Ss[i]):
                if d>0:
                    k=1
                    while 2*k<=d: cs.add(k/float(d)); k+=1
    cs.add(0.5)
    def gf(t):
        m=1.0
        for v in Ss:
            r=(v*t)%1.0; r=min(r,1.0-r)
            if r<m: m=r
        return m
    bt=max(cs,key=gf); return bt, gf(bt)

print("="*70)
print("PART 1: SCALE-INVARIANCE of M -- exact statement & proof check")
print("="*70)
print("""
CLAIM:  M(cS) = M(S) for every positive integer c.
PROOF:  M(cS) = max_t min_v ||c v t||.  As t ranges over R/Z, so does u=c t
        over R/Z (it wraps c times but covers the SAME image, including all of
        the circle).  Sub u=c t:  M(cS) = max_u min_v ||v u|| = M(S).  QED.
COROLLARY: M depends only on the PROJECTIVE class of S.  So WLOG gcd(S)=1
        (primitive) -- exactly the conjecture's hypothesis.  BUT scale-
        invariance does NOT bound max(S): {1,2,..,13} and {2,4,..,26} have the
        SAME M, yet we can also scale UP arbitrarily.  Among PRIMITIVE reps,
        max is the genuine size and is NOT bounded.  This is why quantization
        needs a SECOND ingredient (a max-bound) to close.
""")
for c in (1,2,3,6,10,1000):
    S=[1,2,3,4,5,6,7,8,9,10,11,12,13]
    m,_=M([c*v for v in S])
    print(f"  M({c}*[1..13]) = {m}  (== 1/14? {m==F(1,14)})")
print()

print("="*70)
print("PART 2: the SCALE-INVARIANT quantity  P(S) = M(S) * max(S)")
print("="*70)
print("""
M*max is NOT scale-invariant (max scales, M doesn't).  So that's the WRONG
normalization.  The scale-invariant content is M itself together with the
PROJECTIVE shape.  Re-pose the floor as: among primitive S, is inf M(S) > 0?
That IS the conjecture (inf >= 1/14).  The floor M >= 1/(2max) only helps when
max is bounded.  So compute inf M over EXHAUSTIVE primitive families of bounded
max -- that is a genuine, finite, rigorous certificate for those families.
""")

print("="*70)
print("PART 3: EXHAUSTIVE primitive 13-subsets of {1..W}")
print("="*70)
print()
print("For W just above 13 the count C(W,13) is small enough to fully enumerate.")
print(f"  C(14,13)={1}, C(15,13)={105}, C(16,13)={560}, C(17,13)={2380},")
print(f"  C(18,13)={8568}, C(19,13)={27132}, C(20,13)={77520}, C(21,13)={203490}")
print()
print(f"{'W':>4} {'#prim':>10} {'min M (exact)':>16} {'float':>11} {'<1/14?':>8} {'argmin':>40}")
THR = F(1,14)
overall_below = []
for W in range(13, 22):
    minM = None; argm=None; nprim=0; cnt_below=0
    for comb in itertools.combinations(range(1, W+1), 13):
        if reduce(gcd, comb, 0) != 1:  # must be primitive
            continue
        nprim += 1
        bt, mf = M_fast(comb)
        if mf < float(THR) - 1e-9:
            # exact confirm
            em,_=M(comb)
            if em < THR:
                cnt_below += 1; overall_below.append(comb)
        if minM is None or mf < minM:
            minM = mf; argm = comb
    # exact-confirm the float argmin
    em, eat = M(argm)
    flag = "YES!" if em < THR else "no"
    print(f"{W:>4} {nprim:>10} {str(em):>16} {float(em):>11.7f} {flag:>8} {str(list(argm)):>40}")
print()
print(f"TOTAL exact M<1/14 found across W=13..21: {len(overall_below)}")
if overall_below:
    print("  *** COUNTEREXAMPLES (triple-check below) ***")
    for c in overall_below[:20]:
        em,eat=M(c); print(f"    {list(c)}  M={em}  tau={eat}")
else:
    print("  NONE.  LRC(14) holds for ALL primitive 13-sets with max <= 21.")
    print("  (Every such S has M(S) >= 1/14, with equality only for the tight")
    print("   configs {1..13} and sporadics like {1..11,13,24}? -- note 24>21.)")
print()

print("="*70)
print("PART 4: lcm bound on a putative minimizer (kind-pasteur L-quantization)")
print("="*70)
print("""
kind-pasteur THM-522: L in (1/(14 lcm S)) Z, and L is scale-invariant.  M is the
gap; the tight locus L=0 contains every M-minimizer's projective class.  The M
VALUE itself has denominator d | (v_a +/- v_b) or 2 v_i, hence d | 2 lcm(S) and
d <= 2 max(S).  For a counterexample M=a/d<1/14 we need d>14 (since a>=1 =>
d=14a... no: a/d<1/14 with a>=1 forces d>14).  Smallest possible: a=1,d=15.
Check whether M can EVER equal 1/15 (or any a/d in (0,1/14)) for a 13-set.
""")
# enumerate small a/d in (0,1/14) and ask: does the envelope ever realize it as
# the MINIMUM value (not just as a candidate tau)?  Scan exhaustive small space.
seen_below = {}
for W in range(13, 20):
    for comb in itertools.combinations(range(1, W+1), 13):
        if reduce(gcd, comb, 0) != 1: continue
        bt, mf = M_fast(comb)
        if mf < float(THR)-1e-9:
            em,_=M(comb)
            if em < THR:
                seen_below[em]=seen_below.get(em,0)+1
print(f"distinct exact M-values in (0,1/14) seen over W<=19: {len(seen_below)}")
for v in sorted(seen_below): print(f"   {v} (count {seen_below[v]})")
if not seen_below:
    print("   NONE.  No primitive 13-set with max<=19 achieves M below 1/14.")
print()
print("DONE part 2.")
