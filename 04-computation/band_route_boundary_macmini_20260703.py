#!/usr/bin/env python3
"""
Band-route BOUNDARY (mac-mini-2026-07-03-S21, HYP-3877). The arithmetic band {15..Q} closes
NEAR-EQUAL far families (regime C AND c>=8 near-equal) because they are LOOSE (not tight {AP,GW}).
It should FAIL for the tight DILATED AP {d,2d,..,13d} (=regime B spread, needs q~89). Pin the boundary.

Lonely at a/q <=> va !in {0,1,q-1} (mod q) for all v.
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

def smallest_band_q(speeds, qmax=400):
    for q in range(15, qmax + 1):
        if lonely_at_q(speeds, q) is not None:
            return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def M_scan(speeds, N=90000):
    best = F(0)
    for k in range(1, N):
        t = F(k, N); m = min(min((v*t)%1, 1-(v*t)%1) for v in speeds)
        if m > best: best = m
    return best

if __name__ == "__main__":
    print("(1) CONTROL: dilated AP {d,2d,...,13d} -- tight (M=14/183), SPREAD. band should FAIL / need big q")
    print("=" * 84)
    for d in [1, 2, 3, 5, 14]:
        sp = [d*i for i in range(1, 14)]
        q = smallest_band_q(sp, qmax=250)
        cov = is_covering(sp)
        print(f"  d={d}: {sp[:5]}... covering={cov} smallest band-q={q}  (14/183={float(F(14,183)):.5f})")
    print("  => dilated AP needs a LARGE q (not in {15..33}) -- it's the tight/spread regime B, NOT regime C.\n")

    print("(2) NEAR-EQUAL far families (regime C + c>=8 near-equal): band closes at SMALL q. all nfar.")
    print("=" * 84)
    rng = random.Random(11)
    print(f"{'nfar':>5} {'far_scale':>10} {'#covering':>10} {'worst band-q':>13} {'min M(scan) sample':>18} {'#fail':>6}")
    overall = 0
    for nfar in [7, 8, 9, 10, 11, 12, 13]:
        for far_scale in [200, 2000, 7392]:
            worst = 0; nfail = 0; ncov = 0; minM = F(1)
            for _ in range(300):
                nnear = 13 - nfar
                for _try in range(60):
                    near = rng.sample(range(1, 23), nnear) if nnear > 0 else []
                    w1 = rng.randint(23, far_scale)
                    drifts = sorted(rng.sample(range(0, max(nfar+2, 10)), nfar))  # near-equal, small drifts
                    far = [w1 + d for d in drifts]
                    sp = near + far
                    if len(set(sp)) == 13 and all(v > 22 for v in far) and is_covering(sp):
                        break
                else:
                    continue
                ncov += 1
                q = smallest_band_q(sp, qmax=400)
                if q is None: nfail += 1
                else: worst = max(worst, q)
                if ncov <= 3:  # sample M for a few
                    m = M_scan(sp, N=20000)
                    minM = min(minM, m)
            overall = max(overall, worst)
            print(f"{nfar:>5} {far_scale:>10} {ncov:>10} {worst:>13} {float(minM):>18.4f} {nfail:>6}")
    print(f"\nOVERALL worst band-q over near-equal families (all nfar 7..13) = {overall}, failures above shown.")
    print("=> if bounded + 0 fail: band {15..%d} closes ALL near-equal far blocks (regime C AND c>=8 near-equal)." % overall)
