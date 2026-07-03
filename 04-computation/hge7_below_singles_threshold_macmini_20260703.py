#!/usr/bin/env python3
"""
POTENTIAL CLOSURE of the hge7 leg (mac-mini-2026-07-03-S21, HYP-3877).
kps's crude SINGLES bound discharges regime B only for max-speed (w1) > 22638. BELOW that, does the
arithmetic band {15..Q} close EVERY covering gcd=1 hge7 family?  If worst band-q <= Q0 for all such families
with max-speed < 22638, then:  band {15..Q0} (small speeds)  +  singles bound (large speeds)  =  hge7 CLOSED.

Adversarial: random + composite band-blockers + near-AP + 2-cluster, ALL with max speed < 22638, gcd=1,
covering, >=7 far. Report the worst smallest-lonely-q. (>22638 is handled by the singles bound, not here.)
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)

def smallest_lonely_q(speeds, qmax=1500):
    for q in range(15, qmax + 1):
        dang = {r for r in range(q) if min(r, q - r) * 14 < q}
        for a in range(1, q):
            if gcd(a, q) == 1 and all((v * a) % q not in dang for v in speeds):
                return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def far_count(speeds): return sum(1 for v in speeds if abs(v) > 22)

THRESH = 22638

if __name__ == "__main__":
    rng = random.Random(424242)
    print(f"hge7 band closure BELOW the singles threshold (max speed < {THRESH}). worst smallest-lonely-q.")
    print("=" * 88)
    worst = (0, None); nfail = 0; ntot = 0; over40 = 0

    def consider(speeds):
        global worst, nfail, ntot, over40
        speeds = sorted(set(speeds))
        if len(speeds) != 13: return
        if max(speeds) >= THRESH: return
        if gcd_all(speeds) != 1: return
        if far_count(speeds) < 7: return
        if not is_covering(speeds): return
        ntot += 1
        q = smallest_lonely_q(speeds, qmax=1500)
        if q is None: nfail += 1
        else:
            if q > worst[0]: worst = (q, speeds)
            if q > 40: over40 += 1

    # (A) composite band-blockers (the hard direction), max speed < THRESH
    bandmods = [17,19,23,29,31,25,27,32,16,21,22,33,34,35,37,38,39,40,41,43,15,18,20,24,26,28,30]
    for _ in range(300000):
        speeds = []
        for q in rng.sample(bandmods, rng.randint(7, 12)):
            m = q * rng.randint(1, max(1, (THRESH - 1) // q))
            if m < THRESH: speeds.append(m)
        while len(set(speeds)) < 13:
            speeds.append(rng.randint(2, 22) if rng.random() < 0.4 else rng.randint(23, THRESH - 1))
        consider(speeds)

    # (B) random spread + near-equal + 2-cluster, max speed < THRESH
    for _ in range(300000):
        nfar = rng.choice([7,8,9,10,11,12,13]); nnear = 13 - nfar
        near = rng.sample(range(1,23), nnear) if nnear else []
        kind = rng.choice(["spread","near-equal","2cluster"])
        if kind == "spread":
            far = rng.sample(range(23, THRESH), nfar)
        elif kind == "near-equal":
            w1 = rng.randint(23, THRESH - 60); far = [w1 + rng.randint(0, 55) for _ in range(nfar)]
        else:
            a1 = rng.randint(23, THRESH//2); a2 = rng.randint(23, THRESH-1)
            far = [rng.choice([a1,a2]) + rng.randint(0,30) for _ in range(nfar)]
        consider(near + far)

    # (C) covering gcd=1 APs with max < THRESH
    for start in range(1, 60):
        for step in range(1, 300):
            speeds = [start + step*i for i in range(13)]
            consider(speeds)

    print(f"tested {ntot} covering gcd=1 hge7 families with max speed < {THRESH}")
    print(f"band FAILURES (no lonely q <= 1500): {nfail}")
    print(f"WORST smallest-lonely-q = {worst[0]}")
    print(f"   worst family: {worst[1]}")
    print(f"families needing q > 40: {over40}")
    Q0 = worst[0]
    print(f"\n=> if 0 failures: band {{15..{Q0}}} closes ALL covering gcd=1 hge7 with max speed < {THRESH};")
    print(f"   combined with kps singles bound (max speed > {THRESH}) => the hge7 leg CLOSES with a FINITE band.")
