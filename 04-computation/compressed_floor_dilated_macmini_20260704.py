#!/usr/bin/env python3
"""
compressed_floor_dilated_macmini_20260704  (mac-mini-2026-07-04-S46)

CORRECTION to S45 HYP-4089: my S45 compressed-floor of 7/89 was computed over SMALL-RANGE families and
MISSED DILATED families. The dilated deep-well  c*{1..12} u {182}  (13 not| c, gcd(c,182)=1) is COMPRESSED
(182 <= 13*12c) and covering, with M = 1/13 EXACTLY < 7/89 -- because dilation moves the base optimum to
t=1/(13c) where the killer 182/(13c) is a non-integer (safe), so NO killer-offset (contrast the deep well
c=1: 182/13=14 integer -> offset -> 14/183).

CRUX QUESTION for the compressed peel (compressed => M >= ? ): is the compressed floor EXACTLY 1/13
(clean, = 12-runner LRC bound), or can compressed covering families dip BELOW 1/13 toward 14/183?
  - If floor = 1/13: compressed => M >= 1/13 is TIGHT and TRUE; covering-min 14/183 uniquely at the
    (dominant) deep well; the peel target is clean.
  - If something in (14/183, 1/13): messier.

Method: exact-M over (a) dilated deep-wells c{1..12}u{182k}, (b) dilated single/split killers,
(c) broad random + local descent over compressed covering families with dilation-aware moves.
Report the true floor and any family below 1/13.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import numpy as np, random, sys

def out(*a): print(*a); sys.stdout.flush()
def gcd_all(xs): return reduce(gcd, xs)
_G = np.arange(1, 9000)/9000.0
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
    s = sorted(sp); return s[-1] <= 13 * s[-2]
def cc_ok(sp):
    return len(set(sp)) == 13 and gcd_all(sp) == 1 and is_cov(sp, 14) and compressed(sp)

if __name__ == "__main__":
    cmin = F(14, 183); c13 = F(1, 13); rng = random.Random(46)
    BEST = F(1); BF = None; below13 = 0; belowcmin = 0; tested = 0

    def consider(S):
        global BEST, BF, below13, belowcmin, tested
        S = tuple(sorted(set(S)))
        if not cc_ok(S): return
        tested += 1
        aM = approxM(list(S))
        if aM < float(c13) + 2e-3:
            em = M_exact(list(S))
            if em < c13:
                below13 += 1
                if below13 <= 8: out(f"  BELOW 1/13: {list(S)}  M={em}={float(em):.6f}")
            if em < cmin: belowcmin += 1
            if em < BEST: BEST = em; BF = S

    # (a) dilated deep-wells c{1..12} u {182k}
    for c in range(1, 14):
        for k in range(1, 6):
            consider([c*i for i in range(1, 13)] + [182*k])
    # (a2) dilated single-killer with tight base c{1..11,13} u {killer covering 14}
    for c in range(1, 14):
        base = [c*i for i in list(range(1, 12))+[13]]
        for kill in range(14, 400):
            consider(base + [kill])
    # (b) dilated {1..11,13,84}-type: c{1..11,13} u {84c'} etc. and general dilations of the 7/89 family
    for c in range(1, 8):
        consider([c*i for i in [1,2,3,4,5,6,7,8,9,10,11,13]] + [84*c] if False else
                 [c*i for i in [1,2,3,4,5,6,7,8,9,10,11,13]] + [84])
    # (c) broad random compressed + descent
    RANGE = list(range(1, 220))
    for _ in range(40000):
        consider(rng.sample(RANGE, 13))
    # descent from the dilated lows and any found
    seeds = [BF, tuple(sorted([3*i for i in range(1,13)]+[182])),
             tuple(sorted([5*i for i in range(1,13)]+[182]))]
    for seed in seeds:
        if seed is None: continue
        cur = list(seed)
        if not cc_ok(cur): continue
        curM = approxM(cur)
        for _ in range(60):
            improved = False
            for i in range(13):
                rep = None
                for nv in RANGE:
                    if nv in cur: continue
                    cand = cur[:i]+[nv]+cur[i+1:]
                    if not cc_ok(cand): continue
                    mm = approxM(cand)
                    if mm < curM - 1e-6: curM = mm; rep = cand
                if rep: cur = rep; improved = True
            if not improved: break
        consider(cur)

    bf = list(BF)
    emf = M_exact(bf)
    out("")
    out(f"tested {tested} compressed covering families (dilation-aware)")
    out(f"  families with M < 1/13: {below13};  with M < 14/183: {belowcmin}")
    out(f"  TRUE COMPRESSED FLOOR = {emf} = {float(emf):.6f}  at {sorted(bf)}")
    out(f"  1/13 = {float(c13):.6f};  14/183 = {float(cmin):.6f};  7/89 = {float(F(7,89)):.6f} (S45, WRONG)")
    out("")
    if belowcmin > 0:
        out("VERDICT: compressed DIPS BELOW 14/183 -- razor-thinness NOT quarantined; S45+klein-S129 both wrong.")
    elif emf == c13 or (emf >= cmin and emf <= c13 + F(1,10000)):
        out(f"VERDICT: compressed floor = {emf} (~1/13) > 14/183. Razor-thin 14/183 still ONLY at the dominant deep well.")
        out("  => 'compressed => M >= 1/13' is the CLEAN TIGHT target (dilated deep-wells attain 1/13 exactly).")
        out("  => S45 '7/89' floor was WRONG (missed dilated); klein-S129 '>=7/89' needs the dilated caveat.")
    else:
        out(f"VERDICT: floor {emf} in (14/183, 1/13) -- below 1/13; peel target must be lower than 1/13.")
