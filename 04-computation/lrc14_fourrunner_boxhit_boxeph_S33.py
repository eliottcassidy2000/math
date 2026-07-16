#!/usr/bin/env python3
"""
THE UNIFORM FOUR-RUNNER BOX-HIT BOUND (boxeph-2026-07-16-S33)
Turning codex-THM-891's negative residue-6 kernel (rays {1,5} + {2,4}, K6 = -12,
target A15+A24) into a uniform four-runner law.

Objects: reduced coprime runners v = (v1..v4), boxes B_i subseteq Z_7,
  hit(v, B) = meas{x in [0,1) : floor(7 {v_i x}) in B_i for all i}   (exact, Fractions).
Strata:
  product term  PROD = prod |B_i|/7;
  pair stratum  P2 = sum_{i<j} [hit(v_i,v_j; B_i,B_j) - (|B_i||B_j|/49)] * prod_{k != i,j} |B_k|/7
                 (codex pair law: each bracket = D-matrix/(v_i v_j), |D| <= 12/7);
  remainder     R = hit - PROD - P2   (the genuinely >= 3-linear content).
CLAIM to certify: |R| <= C / (v_a v_b v_c) over the minimal triple products --
uniform in the speeds, residue-dependent only. Battery on the kernel's residue
pattern (1,5,2,4) mod 7 and controls.
"""

import sys
from fractions import Fraction as Fr
from math import gcd

def hit_measure(vs, boxes):
    """exact measure of {x : sec_{v_i}(x) in B_i for all i}."""
    bps = sorted(set(Fr(k, 7 * v) for v in vs for k in range(7 * v + 1)))
    tot = Fr(0)
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        ok = True
        for v, B in zip(vs, boxes):
            if int((v * mid % 1) * 7) not in B:
                ok = False
                break
        if ok:
            tot += bps[i + 1] - bps[i]
    return tot

def strata(vs, boxes):
    n = len(vs)
    full = hit_measure(vs, boxes)
    prod = Fr(1)
    for B in boxes:
        prod *= Fr(len(B), 7)
    p2 = Fr(0)
    for i in range(n):
        for j in range(i + 1, n):
            pij = hit_measure([vs[i], vs[j]], [boxes[i], boxes[j]])
            dev = pij - Fr(len(boxes[i]) * len(boxes[j]), 49)
            rest = Fr(1)
            for k in range(n):
                if k not in (i, j):
                    rest *= Fr(len(boxes[k]), 7)
            p2 += dev * rest
    R = full - prod - p2
    return full, prod, p2, R

if __name__ == "__main__":
    print("=" * 76)
    print("FOUR-RUNNER BOX-HIT STRATA (kernel residue pattern (1,5,2,4) mod 7 + controls)")
    # quadruples with residues (1,5,2,4): the A15+A24 rays
    batteries = [
        ("kernel(1,5,2,4)", [(8, 12, 9, 11), (15, 19, 16, 18), (22, 26, 23, 25),
                              (29, 33, 30, 32), (36, 40, 37, 39), (50, 54, 44, 46),
                              (8, 19, 23, 32), (15, 26, 30, 39), (22, 40, 37, 46)]),
        ("control(1,2,3,4)", [(8, 9, 10, 11), (15, 16, 17, 18), (29, 30, 31, 32)]),
        ("control(1,1,1,1)", [(8, 15, 22, 29), (15, 22, 29, 36)]),
    ]
    # boxes: the R_s-relevant miss-box (sec != 0) x4, and a singleton box battery
    missbox = tuple(range(1, 7))
    for name, quads in batteries:
        print(f"-- {name}")
        for vs in quads:
            assert all(gcd(vs[i], vs[j]) >= 1 for i in range(4) for j in range(4))
            for boxes, bname in [((missbox,) * 4, "miss^4"),
                                 (({0}, {0}, {0}, {0}), "sing 0000"),
                                 (({1}, {5}, {2}, {4}), "sing resid")]:
                boxes = tuple(set(B) if not isinstance(B, set) else B for B in boxes)
                full, prod, p2, R = strata(vs, boxes)
                # normalizations: min triple product and min pair product
                import itertools
                trip = min(a * b * c for a, b, c in itertools.combinations(vs, 3))
                pairmin = min(a * b for a, b in itertools.combinations(vs, 2))
                print(f"   v={vs} box={bname}: dev={float(full-prod):+.3e} "
                      f"pair-part={float(p2):+.3e} R={float(R):+.3e} "
                      f"R*minTriple={float(R)*trip:+.3f}")
    print("done")
