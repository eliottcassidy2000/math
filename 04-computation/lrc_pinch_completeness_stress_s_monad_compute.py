#!/usr/bin/env python3
"""
PINCH-sieve completeness stress test for LRC@14 (13 runners).
monad-compute-2026-06-02.

CONTEXT
-------
HYP-2075 (opus-S562) claims the PINCH sieve -- try t = a/(v_i+v_j) over all
pair-sums, gcd(a, v_i+v_j)=1 -- is COMPLETE for n=14: it finds a witness
t with ||v_i t|| >= 1/14 for EVERY config, measured on 329 random+loaded
configs (100%). The source script lrc_multi_sieve_recursive_s562.py only
exercises ~8 named hard configs for the pinch path. This script tries hard to
FALSIFY the 100% claim at scale.

We do two things per config V (a set of 13 positive speeds):
  (1) PINCH sieve: search t = a/(v_i+v_j). Report whether it finds a 1/14 witness.
  (2) GROUND TRUTH witness search: search ALL t = a/q for q up to QMAX,
      gcd(a,q)=1. Report whether ANY 1/14 witness exists at all.

A genuine REFUTATION of HYP-2075 = a config where (2) finds a witness but (1)
does not (pinch misses a witness that exists). We also flag configs where (2)
finds NOTHING up to QMAX (would be an LRC-relevant near-miss / needs larger q).

All arithmetic is INTEGER (no Fraction):  for t = a/q, runner v is safe iff
r = (v*a) mod q  has  n*min(r, q-r) >= q.

Config families (adversarial + random), n=14 so 13 speeds:
  - AP (1..13)               : the known tight extremal set
  - V* sporadic tight        : (1..11,13,24)
  - contains multiple of 14  : apex obstruction
  - LCM-loaded               : kills small-modulus division sieves
  - geometric / powers       : sparse multiplicative structure
  - random primitive, several speed ranges and many seeds
"""

from math import gcd
import random
from functools import reduce

N = 14
K = N - 1  # 13 runners
THR_NUM, THR_DEN = 1, N  # threshold 1/14


def safe_at(V, a, q):
    """Is t = a/q safe (all ||v t|| >= 1/N)?  integer-only."""
    for v in V:
        r = (v * a) % q
        # need min(r, q-r)/q >= 1/N  <=>  N*min(r,q-r) >= q
        if N * min(r, q - r) < q:
            return False
    return True


def pinch_witness(V):
    """PINCH sieve: t = a/(v_i+v_j) over distinct pair-sums."""
    sums = sorted({V[i] + V[j] for i in range(len(V)) for j in range(i + 1, len(V))})
    for s in sums:
        for a in range(1, s):
            if gcd(a, s) == 1 and safe_at(V, a, s):
                return (a, s)
    return None


def any_witness(V, qmax):
    """Ground truth: search ALL reduced fractions a/q, 2<=q<=qmax."""
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) == 1 and safe_at(V, a, q):
                return (a, q)
    return None


def primitive(V):
    g = reduce(gcd, V)
    return tuple(sorted(v // g for v in set(V)))


def gen_configs(seed):
    rng = random.Random(seed)
    cfgs = {}
    cfgs["AP 1..13"] = tuple(range(1, N))
    cfgs["V* sporadic"] = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)
    cfgs["apex mult-of-14"] = (1, 3, 5, 9, 11, 13, 2, 4, 6, 8, 10, 12, 14)
    cfgs["LCM-loaded"] = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 143, 2520)
    cfgs["two apexes"] = (1, 2, 3, 5, 7, 9, 11, 13, 15, 17, 14, 28, 42)
    cfgs["powers-of-2-ish"] = (1, 2, 4, 8, 16, 32, 64, 128, 3, 5, 7, 9, 11)
    cfgs["geometric-3"] = primitive((1, 3, 9, 27, 81, 243, 2, 4, 5, 7, 8, 10, 11))
    cfgs["near-AP+jump"] = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 100)
    # many random primitive configs across speed ranges
    for rng_hi in (20, 30, 60, 120, 300, 1000):
        for i in range(220):
            while True:
                V = primitive(tuple(rng.sample(range(1, rng_hi + 1), K)))
                if len(V) == K:
                    break
            cfgs[f"rand<= {rng_hi} #{i}"] = V
    return cfgs


def main():
    QMAX = 400  # ground-truth denominator bound
    total = 0
    pinch_caught = 0
    refutations = []        # pinch MISS but witness exists
    nowit = []              # NO witness up to QMAX (LRC near-miss flag)
    pinch_only_big = []     # pinch found, but ground-truth needed q>max pair-sum
    seeds = [1, 7, 42, 100, 2026, 31337]
    for seed in seeds:
        cfgs = gen_configs(seed)
        for name, Vraw in cfgs.items():
            V = tuple(sorted(set(Vraw)))
            if len(V) != K:
                continue
            total += 1
            pw = pinch_witness(V)
            if pw is not None:
                pinch_caught += 1
            else:
                aw = any_witness(V, QMAX)
                if aw is not None:
                    refutations.append((seed, name, V, aw))
                else:
                    nowit.append((seed, name, V))
    print(f"==== PINCH completeness stress test, n={N} ({K} runners) ====")
    print(f"seeds: {seeds}")
    print(f"total configs tested: {total}")
    print(f"PINCH caught (1/{N} witness via pair-sum): {pinch_caught}/{total} "
          f"({100*pinch_caught/total:.3f}%)")
    print(f"PINCH missed: {total - pinch_caught}")
    print()
    if refutations:
        print(f"*** {len(refutations)} REFUTATIONS: pinch MISS but witness exists "
              f"(q<={QMAX}) -> HYP-2075 'pinch complete' is FALSE ***")
        for seed, name, V, aw in refutations[:20]:
            sums = sorted({V[i]+V[j] for i in range(len(V)) for j in range(i+1, len(V))})
            print(f"   seed={seed} {name}: V={V}")
            print(f"       witness t={aw[0]}/{aw[1]} (q={aw[1]}); max pair-sum={sums[-1]}; "
                  f"q in pair-sums? {aw[1] in sums}")
    else:
        print("No refutations: PINCH caught every config that has any witness "
              f"with q<={QMAX}.")
    print()
    if nowit:
        print(f"!! {len(nowit)} configs with NO witness found up to q<={QMAX} "
              f"(LRC near-miss flag; pinch also missed these):")
        for seed, name, V in nowit[:20]:
            print(f"   seed={seed} {name}: V={V}")
    else:
        print(f"Every config has a 1/{N} witness with q<={QMAX}.")
    print()
    print("DONE.")


if __name__ == "__main__":
    main()
