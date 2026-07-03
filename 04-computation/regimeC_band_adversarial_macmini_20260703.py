#!/usr/bin/env python3
"""
Regime-C band argument -- ADVERSARIAL stress (mac-mini-2026-07-03-S21).
Can any covering regime-C family (7 near-equal SMALL far + 6 near<=22) BLOCK a large denominator band?
The potential failure: highly-composite / large w1 far runners blocking many band moduli q via q|(w1+j),
plus near runners blocking {15..22}. If the max smallest-band-q stays SMALL and BOUNDED over all such
adversarial families, regime C closes by a finite band check {15..Q*}.

We also isolate the PURE 'q divides no speed' question (the dominant condition): find the smallest q>=15
dividing NONE of the 13 speeds -- and whether va∉{±1} is then always satisfiable.
"""
from fractions import Fraction as F
from math import gcd

def danger_residues(q):
    return {r for r in range(q) if min(r, q - r) * 14 < q}

def lonely_at_q(speeds, q):
    dang = danger_residues(q)
    for a in range(1, q):
        if gcd(a, q) == 1 and all((v * a) % q not in dang for v in speeds):
            return a
    return None

def smallest_band_q(speeds, qmax=200):
    for q in range(15, qmax + 1):
        if lonely_at_q(speeds, q) is not None:
            return q
    return None

def smallest_q_dividing_none(speeds, qmax=200):
    for q in range(15, qmax + 1):
        if all(v % q != 0 for v in speeds):
            return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

if __name__ == "__main__":
    print("ADVERSARIAL regime-C band stress. worst smallest-band-q over many families.")
    print("=" * 78)
    # highly-composite / near-composite w1 candidates (to maximize far divisors in the band)
    hc = [w for w in range(23, 7393) if len([d for d in range(15, 60) if w % d == 0]) >= 2]
    # near-sets that block band bottom {15..22}
    adversarial_nears = [
        [15,16,17,18,19,20],[16,17,18,19,20,21],[17,18,19,20,21,22],
        [15,17,19,21,22,20],[15,16,18,20,21,22],[16,18,19,20,21,22],
        [14,15,16,18,20,21],[12,15,16,18,20,22],
    ]
    drifts_list = [list(range(7)), [0,1,2,3,4,5,7], [0,1,2,3,5,6,8], [0,2,4,6,8,10,12]]
    worst = (0, None); worst_none = (0, None)
    n = 0; failed = []
    # scan w1 over hc (highly composite) + a stride over the whole range
    w1s = sorted(set(hc[:400] + list(range(23, 7393, 37))))
    for near in adversarial_nears:
        for drifts in drifts_list:
            for w1 in w1s:
                far = [w1 + d for d in drifts]
                speeds = near + far
                if len(set(speeds)) != 13 or not is_covering(speeds):
                    continue
                n += 1
                q = smallest_band_q(speeds, qmax=300)
                qn = smallest_q_dividing_none(speeds, qmax=300)
                if q is None:
                    failed.append((near, far))
                elif q > worst[0]:
                    worst = (q, (near, far))
                if qn and qn > worst_none[0]:
                    worst_none = (qn, (near, far))
    print(f"tested {n} covering adversarial regime-C families (highly-composite + strided w1)")
    print(f"WORST smallest-band-q (lonely a/q exists) = {worst[0]}")
    if worst[1]:
        print(f"   worst family: near={worst[1][0]} far={worst[1][1]}")
    print(f"WORST smallest-q-dividing-no-speed = {worst_none[0]}  (family {worst_none[1]})")
    print(f"families with NO band q <= 300: {len(failed)}")
    if failed:
        for f_ in failed[:5]:
            print("   BAND-FAILURE:", f_)
    print(f"\n=> regime C {'CLOSES by finite band {15..%d}' % worst[0] if not failed and worst[0] else 'has band failures -- deeper'}")
    print("   (band q that works => t=a/q is 1/14-lonely for the whole 13-family; a FINITE check per family)")
