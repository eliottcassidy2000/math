# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont63: attack the general core-runner-1 case via a MEASURE sufficient condition.
#
# A core-runner-1 covering family = {1} U B, B = 12 non-core runners (each divisible by a prime <=13).
# Runner 1's bad set D_1 = {t : ||t|| < 1/14} is a SINGLE arc of measure 2/14 = 1/7. Good set of the body:
#     G' = {t : ||b t|| >= 1/14  for all b in B}.
# KEY: if |G'| > 1/7 = |D_1|, then G' cannot fit inside D_1, so there is t in G'\D_1 where EVERY runner
# (including runner 1) is lonely => M({1}UB) >= 1/14 => LRC(14). So "|G'| > 1/7 for all smooth bodies B"
# would CLOSE the |core|=1 residual for LRC(14). This tests whether it holds (esp. at the deep-well extremal).
from math import gcd
from fractions import Fraction as F
from functools import reduce

P = 2*3*5*7*11*13
def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)
def is_cov(v, N=14): return all(any(x % d == 0 for x in v) for d in range(2, N+1))
def core(v): return [x for x in v if gcd(x, P) == 1]

def good_measure(B, level, grid=600000):
    """measure of {t in [0,1): ||b t|| >= level for all b in B}."""
    cnt = 0
    bs = list(B)
    for k in range(grid):
        t = k / grid
        ok = True
        for b in bs:
            r = (b * t) % 1.0
            if min(r, 1.0 - r) < level:
                ok = False; break
        if ok: cnt += 1
    return cnt / grid

def main():
    lvl = 1.0/14
    thr = 1.0/7
    print(f"level = 1/14 = {lvl:.5f}; runner-1 arc |D_1| = 1/7 = {thr:.5f}")
    print("CLAIM to test: |G'| > 1/7 for every smooth body B (=> M({1}UB) >= 1/14, LRC(14) for |core|=1).\n")
    print(f"{'body B (12 non-core)':<34} | {'1UB cover?':>10} | |core(1UB)| | {'|G_body|':>9} | > 1/7?")
    bodies = {
        "deep well body {2..12,182}": list(range(2,13))+[182],
        "ladder S_2 body {2..12,364}": list(range(2,13))+[364],
        "{2..12, 84} (lcm12,14)": list(range(2,13))+[84],
        "{2..12, 156}(lcm12,13)": list(range(2,13))+[156],
        "{2,3,4,6,8,9,10,12,14,7,11,13}": [2,3,4,6,8,9,10,12,14,7,11,13],
        "even-heavy {2,4,6,8,10,12,14,3,9,5,11,13}": [2,4,6,8,10,12,14,3,9,5,11,13],
        "{4,6,8,9,10,12,14,15,21,22,26,33}": [4,6,8,9,10,12,14,15,21,22,26,33],
    }
    mn = 1.0; mnb = None
    for name, B in bodies.items():
        full = [1] + B
        cov = is_cov(full); nc = len(core(full))
        gm = good_measure(B, lvl)
        if gm < mn: mn, mnb = gm, name
        print(f"{name:<34} | {str(cov):>10} | {nc:>11} | {gm:>9.5f} | {'YES' if gm > thr else 'NO <<<'}")
    print()
    print(f"min |G_body| over these = {mn:.5f} at [{mnb}]; threshold 1/7 = {thr:.5f}")
    print(f"independent model (6/7)^12 = {(6/7)**12:.5f} (> 1/7 by {(6/7)**12 - thr:.5f})")
    print()
    print("=> if min |G'| > 1/7 holds universally, the |core|=1 residual is CLOSED for LRC(14) by a clean")
    print("   MEASURE bound (no discrepancy/mollification needed). The margin is thin -- test the extremal hard.")

if __name__ == "__main__":
    main()
