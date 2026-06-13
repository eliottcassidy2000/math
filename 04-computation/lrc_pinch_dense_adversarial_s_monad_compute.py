#!/usr/bin/env python3
"""
PINCH-sieve adversarial test in the SCARCE-witness regime, LRC@14.
monad-compute-2026-06-02. Follow-up to lrc_pinch_completeness_stress.

The previous broad test (7968 configs, many large-range/random) found 0 pinch
misses. But witnesses are most SCARCE for dense near-AP configs (13 speeds drawn
from a small range), which is exactly the regime where the AP is tight. If pinch
is going to miss, it is here. We hammer that regime and also test the refined
claim: is the minimal-denominator witness always realizable as a pair-sum?

Per config V (13 distinct speeds):
  - pinch_witness: smallest pair-sum a/(v_i+v_j) that is a 1/14 witness
  - min_q_witness: smallest denominator q (any reduced a/q) that is a witness
  - flag if pinch misses but a witness exists  -> REFUTES HYP-2075 completeness
  - record whether min_q is itself a pair-sum   -> tests 'optimal witness is pinch'

Integer arithmetic only.
"""

from math import gcd
import random
from functools import reduce
from itertools import combinations

N = 14
K = N - 1


def safe_at(V, a, q):
    for v in V:
        r = (v * a) % q
        if N * min(r, q - r) < q:
            return False
    return True


def pinch_witness(V):
    sums = sorted({a + b for a, b in combinations(V, 2)})
    for s in sums:
        for a in range(1, s):
            if gcd(a, s) == 1 and safe_at(V, a, s):
                return (a, s)
    return None


def min_q_witness(V, qmax):
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) == 1 and safe_at(V, a, q):
                return (a, q)
    return None


def primitive(V):
    g = reduce(gcd, V)
    return tuple(sorted(v // g for v in set(V)))


def main():
    QMAX = 600
    rng = random.Random(20260602)
    total = 0
    pinch_caught = 0
    refutations = []
    nowit = []
    minq_not_pairsum = []   # min-q witness denominator is NOT a pair-sum
    biggest_minq = (0, None)

    # dense ranges: 13 speeds from {1..hi}, hi just above 13 -> scarcest witnesses
    plans = []
    for hi in (14, 15, 16, 17, 18, 20, 24):
        plans.append((hi, 4000))
    # also include the AP and its single-coordinate perturbations explicitly
    explicit = [tuple(range(1, N))]
    base = list(range(1, N))
    for i in range(K):
        for d in (1, 2, 3, -1):
            w = base[:]
            w[i] = base[i] + d
            if len(set(w)) == K and all(x >= 1 for x in w):
                explicit.append(tuple(sorted(set(w))))

    configs = []
    for V in explicit:
        configs.append(("explicit-AP-perturb", primitive(V)))
    for hi, cnt in plans:
        pool = list(range(1, hi + 1))
        if len(pool) < K:
            continue
        for _ in range(cnt):
            V = primitive(tuple(rng.sample(pool, K)))
            if len(V) == K:
                configs.append((f"dense<= {hi}", V))

    seen = set()
    for name, V in configs:
        if V in seen:
            continue
        seen.add(V)
        total += 1
        pw = pinch_witness(V)
        mq = min_q_witness(V, QMAX)
        if pw is not None:
            pinch_caught += 1
        if mq is not None and mq[1] > biggest_minq[0]:
            biggest_minq = (mq[1], V)
        if pw is None:
            if mq is not None:
                refutations.append((name, V, mq))
            else:
                nowit.append((name, V))
        if mq is not None:
            sums = {a + b for a, b in combinations(V, 2)}
            if mq[1] not in sums and len(minq_not_pairsum) < 25:
                minq_not_pairsum.append((name, V, mq, sorted(sums)[-1]))

    print(f"==== PINCH adversarial dense test, n={N} ({K} runners) ====")
    print(f"distinct configs tested: {total}  (QMAX ground-truth = {QMAX})")
    print(f"PINCH caught: {pinch_caught}/{total} ({100*pinch_caught/total:.4f}%)")
    print(f"PINCH missed: {total - pinch_caught}")
    print(f"largest min-q witness denominator seen: q={biggest_minq[0]}  "
          f"V={biggest_minq[1]}")
    print()
    if refutations:
        print(f"*** {len(refutations)} REFUTATIONS (pinch MISS, witness exists) -> "
              f"HYP-2075 completeness FALSE ***")
        for name, V, mq in refutations[:20]:
            print(f"   {name}: V={V} witness {mq[0]}/{mq[1]}")
    else:
        print("No refutations: pinch caught every config with a witness q<=QMAX.")
    if nowit:
        print(f"!! {len(nowit)} configs with NO witness up to q<={QMAX}:")
        for name, V in nowit[:20]:
            print(f"   {name}: V={V}")
    print()
    print(f"Configs whose SMALLEST-q witness denominator is NOT a pair-sum: "
          f"{'examples below' if minq_not_pairsum else 'NONE in sample'}")
    for name, V, mq, maxsum in minq_not_pairsum[:15]:
        print(f"   {name}: V={V} min-q={mq[0]}/{mq[1]} (q not a pair-sum; "
              f"pinch uses a larger pair-sum<= {maxsum})")
    print()
    print("DONE.")


if __name__ == "__main__":
    main()
