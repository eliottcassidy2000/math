#!/usr/bin/env python3
"""
CORRECTED band-below size (mac-mini-2026-07-03-S22): my S21 {15..33} for near-equal was too tight --
ALIGNED near-equal drifts are band-blockers (witness q up to 47). Find the TRUE worst witness q over ALL
covering gcd=1 hge7 families with max-speed < 22638 (the below-singles-threshold slice), NOW INCLUDING the
aligned band-blockers (near-equal AND spread) my S21 search missed. This sets the honest regime-C / band-below
denominator Q0 so band {15..Q0} + singles(>22638) is the corrected magnitude-split closure.
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def danger(q): return {r for r in range(q) if min(r, q - r) * 14 < q}
def witness_denominator(speeds, qmax=2000):
    for q in range(15, qmax + 1):
        d = danger(q)
        for a in range(1, q):
            if gcd(a, q) == 1 and all((v * a) % q not in d for v in speeds):
                return q
    return None
def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))
def far_count(speeds): return sum(1 for v in speeds if abs(v) > 22)

THRESH = 22638

if __name__ == "__main__":
    rng = random.Random(9001)
    worst = (0, None); nfail = 0; ntot = 0
    def consider(speeds):
        global worst, nfail, ntot
        speeds = sorted(set(int(v) for v in speeds if v > 0))
        if len(speeds) != 13 or max(speeds) >= THRESH: return
        if gcd_all(speeds) != 1 or far_count(speeds) < 7 or not is_covering(speeds): return
        ntot += 1
        q = witness_denominator(speeds, qmax=2000)
        if q is None: nfail += 1
        elif q > worst[0]: worst = (q, speeds)

    band = list(range(15, 70))
    # (A) ALIGNED near-equal band-blockers: runners ~N (N<THRESH) each divisible by band moduli
    for _ in range(150000):
        N = rng.randint(200, THRESH - 1)
        rng.shuffle(band)
        chosen = band[:rng.randint(7, 12)]
        far = sorted({q * round(N / q) for q in chosen if q * round(N / q) > 22 and q * round(N / q) < THRESH})
        speeds = far[:]
        for q in [8,9,5,7,11,13,2,3,4,6]:
            if len(speeds) >= 13: break
            if not any(s % q == 0 for s in speeds): speeds.append(q)
        while len(speeds) < 13: speeds.append(rng.randint(2, 22))
        consider(speeds)

    # (B) ALIGNED SPREAD band-blockers: runners at various magnitudes each divisible by band moduli
    for _ in range(150000):
        rng.shuffle(band)
        chosen = band[:rng.randint(7, 12)]
        far = sorted({q * rng.randint(1, (THRESH - 1) // q) for q in chosen})
        far = [f for f in far if 22 < f < THRESH]
        speeds = far[:]
        for q in [8,9,5,7,11,13,2,3,4,6]:
            if len(speeds) >= 13: break
            if not any(s % q == 0 for s in speeds): speeds.append(q)
        while len(speeds) < 13: speeds.append(rng.randint(2, 22))
        consider(speeds)

    # (C) generic random (near-equal + spread), for coverage
    for _ in range(100000):
        nfar = rng.choice([7,8,9,10,11,12,13]); nnear = 13 - nfar
        near = rng.sample(range(1,23), nnear) if nnear else []
        if rng.random() < 0.5:
            w1 = rng.randint(23, THRESH-60); far = [w1 + rng.randint(0,55) for _ in range(nfar)]
        else:
            far = rng.sample(range(23, THRESH), nfar)
        consider(near + far)

    print(f"CORRECTED band-below: worst witness q over covering gcd=1 hge7 with max-speed < {THRESH}")
    print("=" * 74)
    print(f"tested {ntot} families (incl. ALIGNED near-equal + spread band-blockers)")
    print(f"band FAILURES (no witness q<=2000): {nfail}")
    print(f"WORST witness q = {worst[0]}")
    if worst[1]:
        far = [v for v in worst[1] if v > 22]
        print(f"   worst family far part: {far}")
        print(f"   (max-speed {max(worst[1])}, far {far_count(worst[1])}, span-ratio {max(far)/min(far):.2f})")
    print(f"\n=> CORRECTED band-below = {{15..{worst[0]}}} (was over-optimistically {{15..33}} in S21 for near-equal).")
    print(f"   band {{15..{worst[0]}}} (max-speed<{THRESH}) + singles(>{THRESH}) = the honest magnitude-split closure.")
