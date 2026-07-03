#!/usr/bin/env python3
"""
RECONCILE kps-S28's 'q<=35 independent of magnitude' with my HYP-4040 (witness denom -> infinity).
mac-mini-2026-07-03-S23. kps tested ~407 covering families with clusters up to ~1000 and found all lonely at
p/q with q<=35. But my MISTAKE-095 aligned near-equal band-blockers reach q=47. Are those covering AND
COMPRESSED (kps's exact class: no ratio-13 dominant, forall i exists j!=i |v_i|<13|v_j|)? If yes and q>35,
kps's empirical bound is under-sampled (same MISTAKE-095), and the TRUE bound grows with magnitude.

Uses the rational-witness sieve = kps's lonely14_of_ratio: t=p/q lonely iff for all v, all mult of q, and
p, min-dist(v*p mod q) >= q/14, i.e. (v*p mod q) not in danger. (Search reduced fraction a/q, gcd(a,q)=1 =
p/q; smallest such q is what kps reports.)
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def danger(q): return {r for r in range(q) if min(r, q - r) * 14 < q}
def witness_denominator(speeds, qmax=4000):
    for q in range(2, qmax + 1):          # q from 2 (kps allows the covering sieve q<=14 too)
        d = danger(q)
        for a in range(1, q):
            if gcd(a, q) == 1 and all((v * a) % q not in d for v in speeds):
                return q
    return None
def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))
def is_compressed(speeds):
    return all(any(j != i and abs(speeds[i]) < 13 * abs(speeds[j]) for j in range(len(speeds)))
               for i in range(len(speeds)))
def far_count(speeds): return sum(1 for v in speeds if abs(v) > 22)

if __name__ == "__main__":
    rng = random.Random(2023)
    print("RECONCILIATION: max witness denominator over COVERING + COMPRESSED families vs magnitude N.")
    print("kps-S28 claim: q<=35 independent of magnitude.  HYP-4040: grows ~log(magnitude).")
    print("=" * 88)
    print(f"{'N (far scale)':>14} {'#cov+compressed tested':>24} {'MAX witness q':>14} {'#q>35':>7} {'example far (q>35)':>22}")
    band = list(range(15, 70))
    overall = []
    for N in [200, 1000, 5000, 30000, 300000, 3000000]:
        worst = (0, None); ncc = 0; n_over = 0; ex = None
        for _ in range(30000):
            rng.shuffle(band)
            chosen = band[:rng.randint(7, 12)]
            far = sorted({q * round(N / q) for q in chosen})
            far = [f for f in far if f > 22]
            if len(far) < 7: continue
            speeds = far[:]
            for q in [8,9,5,7,11,13,2,3,4,6]:   # small covering fillers
                if len(speeds) >= 13: break
                if not any(s % q == 0 for s in speeds): speeds.append(q)
            while len(speeds) < 13: speeds.append(rng.randint(2, 22))
            speeds = sorted(set(speeds))[:13]
            if len(speeds) != 13 or gcd_all(speeds) != 1: continue
            if far_count(speeds) < 7 or not is_covering(speeds) or not is_compressed(speeds): continue
            ncc += 1
            q = witness_denominator(speeds, qmax=3000)
            if q and q > worst[0]: worst = (q, speeds)
            if q and q > 35:
                n_over += 1
                if ex is None: ex = [v for v in speeds if v > 22]
        overall.append((N, worst[0]))
        print(f"{N:>14} {ncc:>24} {worst[0]:>14} {n_over:>7} {str(ex[:5])+'...' if ex else '-':>22}")
    print("\nMAX witness q vs magnitude N:", [(N, q) for N, q in overall])
    grows = overall[-1][1] > overall[0][1]
    print(f"=> {'GROWS with magnitude (kps q<=35 is under-sampled; TRUE bound grows ~log N per HYP-4040)' if grows else 'bounded'}")
    print("   IMPLICATION: the bounded-denominator route needs an a-priori MAGNITUDE BOUND on the hard case,")
    print("   or the magnitude split. spread13_lonely (ratio<=13) closes global-ratio<=13; the WIDE-CHAIN")
    print("   compressed families (small near + far cluster, global ratio >>13) are the residual with q>35.")
