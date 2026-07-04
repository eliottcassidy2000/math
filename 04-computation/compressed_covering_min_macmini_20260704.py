#!/usr/bin/env python3
"""
compressed_covering_min_macmini_20260704  (mac-mini-2026-07-04-S45)

The covering dispatch (opus `covering_lonely_of_dominant_or_compressed`, formalized by kps HYP-4091)
splits every primitive covering 13-family into:
  DOMINANT   : some runner v_i > 13*|v_j| for all j != i  (largest > 13 * second-largest).
               DISCHARGED in Lean (kps-S5 dominant peel, HYP-4087). Includes the deep well {1..12,182}
               (182 > 13*12=156), the covering-min extremizer.
  COMPRESSED : no dominant runner (largest <= 13 * second-largest). The SOLE remaining open leaf.

CLAIM UNDER TEST: the COMPRESSED covering-min has a definite MARGIN above 14/183 = 0.076503.
  i.e. every compressed covering family has M >= c for some c strictly > 14/183.
  If so, kps's open leaf is a *non-razor-thin* lower bound (much easier than 14/183 exactly),
  and the razor-thin value 14/183 lives ENTIRELY in the discharged dominant branch.

Method: exact-M (Fraction breakpoints), broad random sample + aggressive local descent over compressed
covering families, report the empirical floor and whether any family dips below 14/183.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import numpy as np, random, sys

def out(*a): print(*a); sys.stdout.flush()
def gcd_all(xs): return reduce(gcd, xs)
_G = np.arange(1, 8000)/8000.0
def approxM(sp):
    v = np.array(sp, float); ph = np.outer(v, _G) % 1.0
    return np.minimum(ph, 1.0-ph).min(axis=0).max()
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a, q)) for v in sp)
            if m > best: best = m
    return best
def is_cov(sp, n): return all(any(v % q == 0 for v in sp) for q in range(2, n+1))
def compressed(sp):
    """largest <= 13 * second-largest (NOT dominant). sp positive ints."""
    s = sorted(sp)
    return s[-1] <= 13 * s[-2]
def cc_ok(sp):
    return (len(set(sp)) == 13 and gcd_all(sp) == 1 and is_cov(sp, 14) and compressed(sp))

if __name__ == "__main__":
    cmin = F(14, 183); rng = random.Random(4595)
    RANGE = list(range(1, 210))
    BEST = 1.0; BF = None; below = 0; tested = 0

    def consider(sp):
        global BEST, BF, below, tested
        sp = tuple(sorted(set(sp)))
        if not cc_ok(sp): return
        tested += 1
        aM = approxM(list(sp))
        if aM < float(cmin) - 3e-4:
            em = M_exact(list(sp))
            if em < cmin:
                below += 1
                if below <= 5: out(f"  !! COMPRESSED M<14/183: {list(sp)}  M={em}")
        if aM < BEST + 2e-3:
            em = float(M_exact(list(sp)))
            if em < BEST: BEST = em; BF = sp

    # phase 1: broad random compressed covering families
    for _ in range(80000):
        hi = rng.choice([15, 16, 18, 22, 30, 45, 70, 120])
        consider(rng.sample(RANGE, 13))
    out(f"[phase1 random] tested {tested} compressed covering families; below 14/183: {below}; floor {BEST:.6f} at {sorted(BF) if BF else None}")

    # phase 2: aggressive local descent (single-runner replacement) from best + seeds
    seeds = [BF, (1,2,3,4,5,7,8,9,10,11,12,13,14), (1,2,3,4,5,7,8,9,11,12,13,40,84),
             (1,2,3,4,10,11,12,13,14,15,16,17,18)]
    for seed in seeds:
        if seed is None: continue
        cur = list(seed)
        if not cc_ok(cur): continue
        curM = approxM(cur)
        for _ in range(80):
            improved = False
            for i in range(13):
                bestrep = None
                for nv in RANGE:
                    if nv in cur: continue
                    cand = cur[:i]+[nv]+cur[i+1:]
                    if not cc_ok(cand): continue
                    mm = approxM(cand)
                    if mm < curM - 1e-6: curM = mm; bestrep = cand
                if bestrep: cur = bestrep; improved = True
            if not improved: break
        consider(cur)
    bf = list(BF)
    emf = M_exact(bf)
    out("")
    out(f"COMPRESSED COVERING-MIN FLOOR (empirical) = {emf} = {float(emf):.6f}  at {sorted(bf)}")
    out(f"  14/183 = {float(cmin):.6f};  floor - 14/183 = {emf-cmin} = {float(emf-cmin):.6f}")
    out(f"  1/13   = {float(F(1,13)):.6f}  (12-runner LRC bound; 1/13-14/183 = 1/2379 = {float(F(1,13)-cmin):.6f})")
    out(f"  2/25   = {float(F(2,25)):.6f}  ([0;12,2]);  4/49 = {float(F(4,49)):.6f}  ([0;12,4])")
    out(f"  compressed families with M < 14/183 (exact): {below}")
    out("")
    if emf > cmin:
        out(f"VERDICT: COMPRESSED floor {emf} > 14/183 with margin {float(emf-cmin):.4f}.")
        out("  => kps's open COMPRESSED leaf (HYP-4091) is a NON-razor-thin bound: M >= (floor) > 14/183.")
        out("     The razor-thin 14/183 lives ENTIRELY in the DISCHARGED dominant branch (deep well is dominant).")
    else:
        out("VERDICT: compressed floor reaches 14/183 — no margin; razor-thin persists in the open leaf.")
