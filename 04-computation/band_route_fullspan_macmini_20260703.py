#!/usr/bin/env python3
"""
Does the band route cover the FULL regime C (mac-mini-2026-07-03-S21, HYP-3877)?
kps regime C = D*L < 2 with D = far span, L ~ 2/w1  =>  span D < w1 (NOT just span<=6!).
Test far clusters with span up to w1 (the whole regime-C width) + track where the band {15..Q} fails,
and whether the pair-floor (D*L>=2, i.e. span>=w1) takes over -- looking for a GAP.
"""
from fractions import Fraction as F
from math import gcd
import random

def danger_residues(q):
    return {r for r in range(q) if min(r, q - r) * 14 < q}

def lonely_at_q(speeds, q):
    dang = danger_residues(q)
    for a in range(1, q):
        if gcd(a, q) == 1 and all((v * a) % q not in dang for v in speeds):
            return a
    return None

def smallest_band_q(speeds, qmax=600):
    for q in range(15, qmax + 1):
        if lonely_at_q(speeds, q) is not None:
            return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def far_cluster(rng, w1, span, nfar):
    """nfar runners in [w1, w1+span], w1 included, spread across the span."""
    if span < nfar - 1:
        span = nfar - 1
    pts = {w1, w1 + span}
    while len(pts) < nfar:
        pts.add(w1 + rng.randint(0, span))
    return sorted(pts)

if __name__ == "__main__":
    rng = random.Random(7)
    print("Band route vs far SPAN (regime C = span < w1). worst band-q as span grows toward w1.")
    print("=" * 86)
    print(f"{'w1':>6} {'span frac':>10} {'#covering':>10} {'worst band-q':>13} {'#fail(>600)':>12} {'regime (D*L)':>14}")
    for w1 in [50, 200, 1000, 5000]:
        for frac_label, frac in [("6/const", None), ("w1/8", F(1,8)), ("w1/4", F(1,4)),
                                  ("w1/2", F(1,2)), ("w1 (~C/B edge)", F(1,1)), ("2*w1 (spread/B)", F(2,1))]:
            span = 6 if frac is None else int(frac * w1)
            worst = 0; nfail = 0; ncov = 0
            for _ in range(250):
                nfar = rng.choice([7, 8, 9])
                nnear = 13 - nfar
                near = rng.sample(range(1, 23), nnear)
                far = far_cluster(rng, w1, span, nfar)
                sp = near + far
                if len(set(sp)) != 13 or not all(v > 22 for v in far) or not is_covering(sp):
                    continue
                ncov += 1
                q = smallest_band_q(sp, qmax=600)
                if q is None: nfail += 1
                else: worst = max(worst, q)
            # D*L ~ span * (2/w1); regime C iff <2 iff span<w1
            DL = float(F(span * 2, w1)) if w1 else 0
            reg = "C (band)" if DL < 2 else "B (pair-floor)"
            print(f"{w1:>6} {frac_label:>10} {ncov:>10} {worst:>13} {nfail:>12}   {DL:>5.2f} {reg}")
    print("\n=> reads off the band route's TRUE scope within regime C, and whether a gap opens as span->w1.")
    print("   (worst band-q staying small across span<w1 => band covers ALL of regime C, not just tight clusters)")
