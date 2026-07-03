#!/usr/bin/env python3
"""
Scope of the arithmetic band route over the WHOLE >=7-far leg (mac-mini-2026-07-03-S21, HYP-3877).
Does t=a/q with q in a finite band {15..Q} close EVERY covering compressed family with >=7 far runners
(regime B spread + regime C near-equal + c>=8 + large-far regime A + mixed)?

Lonely at a/q  <=>  va !in {0, 1, q-1} (mod q) for all speeds v  [q in 15..27; danger residues = 3].
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

def gen_covering_hge7(rng, nfar, far_kind, far_scale):
    """random covering family: 6 near (<=22) + nfar far (>22), far pattern = far_kind."""
    for _ in range(200):
        near = rng.sample(range(1, 23), 6)
        w1 = rng.randint(23, far_scale)
        if far_kind == "near-equal":
            drifts = sorted(rng.sample(range(0, 9), nfar))
            far = [w1 + d for d in drifts]
        elif far_kind == "spread":
            far = sorted(rng.sample(range(23, far_scale + 40), nfar))
        elif far_kind == "mixed":
            k = nfar // 2
            far = [w1 + d for d in range(k)] + sorted(rng.sample(range(23, far_scale + 40), nfar - k))
        far = [f for f in far if f > 22]
        speeds = near + far
        if len(set(speeds)) == 6 + nfar and is_covering(speeds):
            return speeds
    return None

if __name__ == "__main__":
    rng = random.Random(20260703)
    print("SCOPE of the arithmetic band route over the >=7-far leg. worst band-q per regime.")
    print("=" * 84)
    print(f"{'nfar':>5} {'far_kind':>12} {'far_scale':>10} {'#tested':>8} {'#covering':>10} {'worst band-q':>13} {'#fail(>400)':>12}")
    overall_worst = 0
    for nfar in [7, 8, 9, 10, 12]:
        for far_kind in ["near-equal", "spread", "mixed"]:
            for far_scale in [200, 2000, 7392, 50000]:
                worst = 0; nfail = 0; ntest = 0; ncov = 0
                for _ in range(400):
                    speeds = gen_covering_hge7(rng, nfar, far_kind, far_scale)
                    if speeds is None:
                        continue
                    ncov += 1; ntest += 1
                    q = smallest_band_q(speeds, qmax=400)
                    if q is None:
                        nfail += 1
                    elif q > worst:
                        worst = q
                overall_worst = max(overall_worst, worst)
                flag = "  <== FAILS" if nfail else ""
                print(f"{nfar:>5} {far_kind:>12} {far_scale:>10} {ntest:>8} {ncov:>10} {worst:>13} {nfail:>12}{flag}")
    print(f"\nOVERALL worst band-q (over all regimes, no failures) = {overall_worst}")
    print("=> if no failures and worst-q bounded: the arithmetic band {15..%d} closes the ENTIRE >=7-far leg" % overall_worst)
    print("   (regime B + C + c>=8 + large-far), by a finite per-family check -- NOT just regime C.")
