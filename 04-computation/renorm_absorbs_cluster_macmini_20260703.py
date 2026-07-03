#!/usr/bin/env python3
"""
HYP-4041 core test (mac-mini-2026-07-03-S23): the aligned band-blockers FAIL the small-q witness but are
LONELY at fine t via the fast phase (renormalization absorbs the cluster). If M(family) >> 1/14 (loose), the
bounded-denominator census's failure is NOT the family being tight -- it's the census using pinned rationals;
the renormalization/fine-t route closes them. This completes the architecture: census for bounded magnitude,
renormalization for large clusters.

MECHANISM: for far cluster {N+c_i}, ||(N+c_i)t|| = ||phi + c_i t|| with phi={Nt}. At a small-q rational t=a/q,
phi=(Na mod q)/q is PINNED (aligned families arrange it into danger). At a FINE t, phi is FREE, and the cluster
{c_i t} (small drifts c_i) is a short arc placeable in the safe region by choosing phi -> family lonely.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1.0; return min(x, 1 - x)
def danger(q): return {r for r in range(q) if min(r, q - r) * 14 < q}
def witness_denominator(speeds, qmax=3000):
    for q in range(2, qmax + 1):
        d = danger(q)
        for a in range(1, q):
            if gcd(a, q) == 1 and all((v * a) % q not in d for v in speeds):
                return q
    return None
def M_fine(speeds, N=4000000):
    """max_t min_i ||v_i t|| via fine scan over t in (0,1); returns (M, argmax t)."""
    vmax = max(speeds); best = 0.0; bestt = 0.0
    # scan at resolution ~ 1/(vmax*K): enough to resolve the fastest runner
    K = 40
    steps = min(N, vmax * K)
    for k in range(1, steps):
        t = k / steps
        m = min(nd(v * t) for v in speeds)
        if m > best: best, bestt = m, t
    return best, bestt
def is_covering(speeds): return all(any(v % q == 0 for v in speeds) for q in range(2, 15))
def is_compressed(speeds):
    return all(any(j != i and abs(speeds[i]) < 13*abs(speeds[j]) for j in range(len(speeds))) for i in range(len(speeds)))

def make_aligned(N, band_moduli, rng):
    far = sorted({q * round(N / q) for q in band_moduli})
    far = [f for f in far if f > 22]
    speeds = far[:]
    for q in [8,9,5,7,11,13,2,3,4,6]:
        if len(speeds) >= 13: break
        if not any(s % q == 0 for s in speeds): speeds.append(q)
    while len(speeds) < 13: speeds.append(rng.randint(2, 22))
    return sorted(set(speeds))[:13]

if __name__ == "__main__":
    rng = random.Random(4041)
    print("Are aligned band-blockers (small-q FAILS) LONELY at fine t (renormalization absorbs the cluster)?")
    print("=" * 90)
    print(f"{'far scale N':>12} {'witness q (small)':>17} {'M (fine scan)':>14} {'M vs 1/14':>12} {'lonely@fine t':>14}")
    found = 0
    for N in [1000, 5000, 30000]:
        # find an aligned covering+compressed family with LARGE small-q witness
        band = list(range(15, 60))
        target = None
        for _ in range(20000):
            rng.shuffle(band)
            sp = make_aligned(N, band[:rng.randint(8,11)], rng)
            if len(sp) != 13 or gcd_all(sp) != 1 or not is_covering(sp) or not is_compressed(sp): continue
            q = witness_denominator(sp, qmax=2000)
            if q and q > 38:
                target = (sp, q); break
        if not target:
            print(f"{N:>12} {'(no q>38 found)':>17}")
            continue
        sp, q = target
        M, tstar = M_fine(sp, N=6000000)
        found += 1
        print(f"{N:>12} {q:>17} {M:>14.5f} {M/(1/14):>11.2f}x {str(M > 1/14):>14}")
        far = sorted(v for v in sp if v > 22)
        print(f"             family far={far}")
        print(f"             lonely fine t*={tstar:.8f}  (1/q would be {1/q:.5f}); cluster spread={max(far)-min(far)}")
    print(f"\n=> if M >> 1/14 for these small-q-FAILING families: they are LOOSE, lonely at FINE t --")
    print("   the census's small-q failure is the pinned-phi artifact, NOT tightness. Renormalization")
    print("   (fast phase phi absorbs the near-equal cluster) closes them; census handles bounded magnitude.")
    print("   => HYP-4041 architecture: the two routes are COMPLEMENTARY, together covering all compressed families.")
