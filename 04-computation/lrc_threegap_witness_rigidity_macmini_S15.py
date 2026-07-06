#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S15 -- the THREE-GAP WITNESS RIGIDITY: is the loneliness
spectrum quantized because near-tight families have a THREE-GAP ({k*alpha}) orbit
at their witness?

GOVERNING PATTERN (synthesis of the repo's threads):
 * LRC(AP) IS the three-distance (Sos-Steinhaus) theorem (opus 20260630).
 * M(S) = 1/(smallest surviving resonance); killing-vs-compactness tension; the AP
   is the unique least-spread killer (kps S31p).
 * The spectrum is the OSTROWSKI ladder k/(k(n-1)+1); three-gap is the rigidity;
   'deep well & AP have g=2 gaps, a generic covering family has g=5' (mac-mini S38).
 * The huge tail is the construction DILATED (Steinhaus scale-invariance, S73).

HYPOTHESIS (the mechanism behind gap-emptiness): at the witness t* achieving M(S),
the runner phases {v_i t* mod 1} together with 0 form a set whose GAP MULTISET has
few distinct values (g<=3 = a three-gap/{k alpha} signature) EXACTLY for near-tight
families, and MANY for loose ones.  If near-tight => g<=3 => the phases are an
arithmetic {k alpha} orbit => M is a continued-fraction convergent (an Ostrowski
rung) => M cannot lie strictly BETWEEN consecutive rungs => the gap (1/13,2/25) is
empty.  This reframes (G) as a THREE-GAP RIGIDITY, not a measure problem.

THIS SCRIPT: for many covering 12-families, find the witness t*, compute the number
g of distinct gaps of {0} u {v_i t*}, and correlate g with M(S).  Test: is g small
(<=3) precisely for the near-tight (small-M) families?
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random, sys
sys.path.insert(0,'04-computation')
from lonely_profile import profile

def log(m=""): print(m, flush=True)

def M_and_witness(S, grid=400000):
    """max_t min_i ||v_i t||; return (M_lb, t*) via fine grid."""
    best=-1.0; targ=0.0
    for j in range(1,grid):
        t=j/grid; mn=1.0
        for v in S:
            x=(v*t)%1.0; d=min(x,1-x)
            if d<mn: mn=d
            if mn<=best: break
        if mn>best: best=mn; targ=t
    return best, targ

def num_gaps(S, t, tol=1e-6):
    """number of DISTINCT gap lengths of {0} u {v_i t mod 1} on the circle (Steinhaus g)."""
    pts=sorted(set([0.0]+[ (v*t)%1.0 for v in S ]))
    gaps=[]
    for i in range(len(pts)):
        g=(pts[(i+1)%len(pts)]-pts[i])%1.0
        if g<=0: g+=1.0
        gaps.append(g)
    # cluster gap values
    gaps.sort(); distinct=[]
    for g in gaps:
        if not distinct or abs(g-distinct[-1])>tol:
            distinct.append(g)
    return len(distinct)

def covering(W,upto=12): return all(any(w%m==0 for w in W) for m in range(2,upto+1))

# curated families across the spectrum + random covering ones
curated={
 "AP {1..12} (tight 1/13)": list(range(1,13)),
 "doubled-apex {1..11,24} (2/25)": list(range(1,12))+[24],
 "block lift (2/25)": [1,2,3,5,7,8,9,10,11,12,17,19],
 "deep well {1..11,168}": list(range(1,12))+[168],
 "{1..11,23} (1/12)": list(range(1,12))+[23],
 "{1..11,13} (1/12)": list(range(1,12))+[13],
}
log("family                              M(S)        g(#gaps at witness)")
log("-"*72)
for name,S in curated.items():
    M,t=M_and_witness(S); g=num_gaps(S,t)
    log(f"{name:36s} {M:.5f}     g={g}")

log("\nrandom covering 12-families (M vs g):")
random.seed(15); rows=[]
tries=0
while len(rows)<14 and tries<20000:
    tries+=1
    W=sorted(random.sample(range(1,40),12))
    if reduce(gcd,W)!=1 or not covering(W): continue
    M,t=M_and_witness(W); g=num_gaps(W,t)
    rows.append((M,g,tuple(W)))
for M,g,W in sorted(rows):
    tag="near-tight" if M<0.095 else ("mid" if M<0.13 else "loose")
    log(f"  M={M:.5f}  g={g:<2} {tag:11s} {list(W)}")

log("\nREADING: if g is SMALL (<=3) exactly for near-tight families and grows for loose")
log("ones, the loneliness spectrum is three-gap-quantized -- the witness of a near-tight")
log("family is a {k*alpha} orbit, its M is a CF/Ostrowski rung, and no value can sit")
log("strictly inside the Farey cell (1/13,2/25).  That is the RIGIDITY behind (G).")
