#!/usr/bin/env python3
"""
STRESS my S21 regime-C claim (mac-mini-2026-07-03-S22): can a NEAR-EQUAL family be a band-blocker?
Construction: far runners ~ N, runner_q = q * round(N/q) (drift <= q <= ~50, so SMALL span = near-equal),
each runner divisible by a distinct band modulus q => at t=a/q that runner sits at residue 0 (danger).
If enough band moduli are each blocked by a near-equal runner, the witness denominator could exceed 33 --
which would REFUTE 'near-equal => band {15..33}'. My S21 tests used small RANDOM drifts and may have missed
these ALIGNED drifts. Test honestly.
"""
from math import gcd
from functools import reduce
from sympy import primerange
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
def is_compressed(speeds):
    return all(any(j != i and abs(speeds[i]) < 13 * abs(speeds[j]) for j in range(len(speeds)))
               for i in range(len(speeds)))
def span_ratio(far): return max(far) / min(far)

if __name__ == "__main__":
    rng = random.Random(77)
    print("NEAR-EQUAL band-blocker stress: runners ~N each divisible by a band modulus. witness q vs N.")
    print("=" * 90)
    print(f"{'N':>8} {'far runners (aligned to band moduli)':>44} {'span-ratio':>10} {'q(witness)':>11}")
    worst_overall = 0
    for N in [200, 1000, 5000, 30000, 200000, 2000000]:
        # block as many band moduli as possible with near-equal runners ~N
        band = list(range(15, 60))
        best = (0, None)
        for _try in range(3000):
            rng.shuffle(band)
            chosen = band[:rng.randint(7, 12)]
            far = sorted({q * round(N / q) for q in chosen})
            far = [f for f in far if f > 22]
            if len(far) < 7: continue
            # near-equal? span must be small relative to N
            if span_ratio(far) > 1.5: continue  # enforce genuine near-equal (ratio < 1.5)
            # add small covering runners to fix covering + reach 13
            speeds = far[:]
            for q in [8,9,5,7,11,13,2,3,4,6]:
                if len(speeds) >= 13: break
                if not any(s % q == 0 for s in speeds): speeds.append(q)
            while len(speeds) < 13: speeds.append(rng.randint(2, 22))
            speeds = sorted(set(speeds))[:13]
            if len(speeds) != 13 or gcd_all(speeds) != 1: continue
            if far_count(speeds) < 7 or not is_covering(speeds): continue
            q = witness_denominator(speeds, qmax=1500)
            if q and q > best[0]: best = (q, speeds)
        if best[1]:
            far = [v for v in best[1] if v > 22]
            sr = span_ratio(far)
            worst_overall = max(worst_overall, best[0])
            print(f"{N:>8} {str(far[:6])+'...':>44} {sr:>10.3f} {best[0]:>11}")
        else:
            print(f"{N:>8} {'(no near-equal (ratio<1.5) blocker found)':>44} {'-':>10} {'-':>11}")
    print(f"\nWORST witness q over near-equal (ratio<1.5) band-blockers = {worst_overall}")
    print(f"=> if > 33: my S21 regime-C claim 'near-equal => band {{15..33}}' is REFUTED for ALIGNED drifts;")
    print(f"   the honest scope is near-equal with GENERIC (small, unaligned) drifts, OR bounded span+magnitude.")
    print(f"   if <= 33: near-equal genuinely closes at 33 even with aligned drifts (claim robust).")
