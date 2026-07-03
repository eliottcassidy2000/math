#!/usr/bin/env python3
"""
Regime-C BAND ARGUMENT (mac-mini-2026-07-03-S21): can the finite denominator band {15..27} always
close a covering regime-C family?

KEY FACT: at t=a/q with q in {15,...,27}, gcd(a,q)=1, a runner v is 1/14-LONELY iff
  min(va mod q, q - va mod q) >= q/14  <=>  va mod q NOT in {0, 1, q-1}   (the 3 danger residues).
So the family is lonely at a/q iff  va ≢ 0, ±1 (mod q) for ALL speeds v.
Since covering only constrains q<=14, a band modulus q in {15..27} may divide NO speed -> then just
need `a` with va ∉ {±1} mod q for all v.

THE MAKE-OR-BREAK QUESTION: can an adversarial COVERING family (7 near-equal far + near-set chosen to
block band moduli) BLOCK THE ENTIRE band {15..27} (no q works)?  If not -- if some q<=Qmax always
works -- regime C closes by a FINITE band check.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def danger_residues(q):
    """residues r mod q that are NOT 1/14-safe: min(r,q-r) < q/14."""
    return {r for r in range(q) if min(r, q - r) * 14 < q}

def lonely_at_q(speeds, q):
    """is there a coprime a with all (v*a mod q) safe? return the a (or None)."""
    dang = danger_residues(q)
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        if all((v * a) % q not in dang for v in speeds):
            return a
    return None

def smallest_band_q(speeds, qmax=60):
    for q in range(15, qmax + 1):
        a = lonely_at_q(speeds, q)
        if a is not None:
            return q, a
    return None, None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def near_equal_far(w1, drifts):
    return [w1 + d for d in drifts]

if __name__ == "__main__":
    print("=" * 80)
    print("REGIME-C BAND ARGUMENT: smallest q in band {15..} closing covering families")
    print("=" * 80)
    # Adversarial near-sets designed to BLOCK band moduli (near runners = band moduli themselves)
    # plus far = 7 near-equal. We hunt for families that need a LARGE band q (or block the whole band).
    worst_q = 0
    worst_fam = None
    n_tested = 0
    n_covering = 0
    # try many near-sets (6 from a pool that includes band-blockers 15..22) + far consecutive at many w1
    pool = list(range(1, 23))  # near runners <= 22
    import random
    rng = random.Random(7)
    # systematic: near-sets that block as much band as possible + cover 2..14
    adversarial_nears = [
        [15, 16, 17, 18, 19, 20],   # blocks band q=15..20
        [16, 17, 18, 19, 20, 21],   # blocks 16..21
        [17, 18, 19, 20, 21, 22],   # blocks 17..22
        [15, 18, 20, 21, 22, 12],   # scattered blockers + 12
        [16, 18, 20, 21, 22, 15],
    ]
    # far drift patterns (near-equal): consecutive, and small-gap
    drift_pats = [list(range(7)), [0,1,2,3,4,5,7], [0,1,2,4,5,6,8], [0,2,3,4,5,6,7]]
    for near in adversarial_nears:
        for drifts in drift_pats:
            for w1 in range(23, 900):
                far = near_equal_far(w1, drifts)
                speeds = near + far
                if len(set(speeds)) != 13:
                    continue
                if not is_covering(speeds):
                    continue
                n_covering += 1
                q, a = smallest_band_q(speeds, qmax=80)
                n_tested += 1
                if q is None:
                    print(f"!!! BAND FAILS (no q<=80): near={near} far={far}")
                    worst_fam = (near, far, None)
                elif q > worst_q:
                    worst_q = q
                    worst_fam = (near, far, q)
    print(f"tested {n_tested} covering regime-C families (adversarial near-sets x drift patterns x w1)")
    print(f"WORST smallest-band-q = {worst_q}   (family: near={worst_fam[0]}, far={worst_fam[1]})" if worst_fam else "none")
    print(f"=> if worst_q is a small fixed number, regime C closes by checking the FINITE band q in 15..{worst_q}")

    print("\n" + "=" * 80)
    print("Stress: for a FIXED far w1, how much of the band {15..30} does the family block?")
    print("=" * 80)
    for near, drifts, w1 in [([15,16,17,18,19,20], list(range(7)), 231),
                             ([17,18,19,20,21,22], list(range(7)), 391),
                             ([15,16,17,18,19,20], [0,1,2,3,4,5,7], 154)]:
        far = near_equal_far(w1, drifts); speeds = near + far
        cov = is_covering(speeds)
        works = [q for q in range(15, 31) if lonely_at_q(speeds, q) is not None]
        print(f"near={near} far={far} covering={cov}: band-q that WORK in 15..30 = {works}")
