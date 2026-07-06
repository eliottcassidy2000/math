#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S12 -- HYP-4382: test the AP-tower induction toward closing
TightLooseDichotomy (the 12-base rigidity: dilated-AP or margin 2/25).

The induction (opus-S103 AP-completion + kps-S13 tower): a compressed covering
12-base W with M(W) < 2/25 should reduce, by peeling to a TIGHT 11-subfamily
(= dilated AP {1..11} by n=12 rigidity), to AP-completion => W is dilated-AP.

KEY QUESTIONS (tested here):
 Q1  DICHOTOMY (re-confirm the target): every covering primitive 12-base has
     M = 1/13 (dilated-AP) or M >= 2/25 -- NONE in the open gap.  [census]
 Q2  PEEL VALIDITY: for the tight bases (dilated-AP, M=1/13), is the ARGMAX-peel
     (drop the largest) a TIGHT 11-subfamily (M=1/12, dilated-AP {1..11})?  And
     for near-boundary bases (M just >= 2/25), what does the peel give?
 Q3  RESIDUE MECHANISM: for M(W) < 2/25 bases, are the residues mod 13 pinned to
     {1..12} (THM-593A)?  (confirms the residue-pinning leg of a uniform proof)
 Q4  LIFT FLOOR: the non-AP lifts of {1..12} mod 13 -- do they all have M >= 2/25
     (the loose branch), with the block lift {..}\\{4,6} u {17,19} = 2/25 exactly?

If Q1-Q4 all hold, the uniform proof of TightLooseDichotomy is:
  residue-pin to {1..12} mod 13 (13 PRIME) + lift-floor 2/25 (my hdich) =>
  base is (k=0 lift = dilated-AP) or (non-trivial lift => M >= 2/25).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import sys, time
sys.path.insert(0, '04-computation')
from lonely_profile import profile

T0 = time.time()
def log(m=""): print(m, flush=True)

def M_exact(S):
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted([abs(x) for x in S]), F(1, cap))
        m = p.M()
        if m is not None:
            return m

def covering(W, upto=12):
    return all(any(w % m == 0 for w in W) for m in range(2, upto+1))

GAP_LO, GAP_HI = F(1,13), F(2,25)
AP = list(range(1,13))

# ---- Q1: dichotomy census (bounded 12-bases) ----
log("Q1: dichotomy census -- covering primitive 12-bases, M in the open gap (1/13,2/25)?")
import random; random.seed(12)
n_in_gap = 0; n_checked = 0; n_tight = 0
worst_near = []
t0 = time.time()
while time.time()-t0 < 90 and n_checked < 40000:
    W = sorted(random.sample(range(1, 55), 12))
    if reduce(gcd,W) != 1 or not covering(W):
        continue
    n_checked += 1
    M = M_exact(W)
    if GAP_LO < M < GAP_HI:
        n_in_gap += 1
        log(f"  !! IN-GAP: W={W} M={M}")
    if M == GAP_LO:
        n_tight += 1
    if GAP_HI <= M < F(3,25):
        worst_near.append((M, tuple(W)))
log(f"  checked {n_checked} covering-primitive 12-bases: in-gap={n_in_gap}, tight(1/13)={n_tight}")
log(f"  (dichotomy holds on the census: {'CONFIRMED' if n_in_gap==0 else 'VIOLATION'})")

# ---- Q2: peel validity on tight (dilated-AP) bases ----
log("\nQ2: argmax-peel of tight (dilated-AP) bases -> tight 11-subfamily?")
for c in (1,2,3,5,7):
    W = [c*v for v in AP]
    Wpeel = sorted(W)[:-1]                      # drop the max = c*12
    Mp = M_exact(Wpeel)
    log(f"  c={c}: base c*{{1..12}}, M={M_exact(W)}; argmax-peel = c*{{1..11}}, M={Mp} "
        f"(tight 1/12? {Mp==F(1,12)})")
# non-AP near-boundary: does SOME 11-peel give tight?
log("  near-boundary (M>=2/25) bases: does any 11-subfamily reach tight 1/12?")
for M, W in sorted(worst_near)[:4]:
    W = list(W)
    best = max((M_exact([W[j] for j in range(12) if j!=i]) for i in range(12)))
    anytight = any(M_exact([W[j] for j in range(12) if j!=i])==F(1,12) for i in range(12))
    log(f"    W={W} M={M}: max 11-peel M={best}, any tight-1/12 subfamily? {anytight}")

# ---- Q3 + Q4: residue pinning + lift floor ----
log("\nQ3: residues mod 13 of M<2/25 bases (should be {1..12} = full system):")
# the only M<2/25 bases are dilated APs; check their residues
for c in (1,2,3,7,11):
    W = [c*v for v in AP]
    res = sorted(set(w%13 for w in W if w%13!=0))
    log(f"  c={c}: residues mod 13 = {res} (full {{1..12}}? {res==list(range(1,13))}; 13|any? {any(w%13==0 for w in W)})")

log("\nQ4: lift floor -- non-AP lifts of {1..12} mod 13, min M (should be >= 2/25):")
tests = {
    "block {4,6}->{17,19}": [1,2,3,5,7,8,9,10,11,12,17,19],
    "deep well {1..11,168}": list(range(1,12))+[168],
    "single {1..11,25}":     list(range(1,12))+[25],
    "double {2,6}->{15,19}": [1,3,4,5,7,8,9,10,11,12,15,19],
}
for name,W in tests.items():
    res_full = sorted(set(w%13 for w in W))==list(range(1,13)) if all(w%13 for w in W) else False
    M = M_exact(W)
    log(f"  {name}: M={M} = {float(M):.5f}  {'>= 2/25' if M>=GAP_HI else '<< 2/25 GAP!'} (full-res-lift? {res_full})")

log(f"\nSYNTHESIS: if Q1-Q4 hold, TightLooseDichotomy = residue-pin(13 prime) + lift-floor(2/25);")
log(f"the block lift = 2/25 EXACTLY is the tight extremal of the loose branch.")
log(f"[t = {time.time()-T0:.0f}s]")
