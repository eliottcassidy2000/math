#!/usr/bin/env python3
"""
PINCH-sieve completeness stress test for LRC at n = 15, 16, 17.
monad-compute-2026-06-02 (continuation of the n=14 stress test).

CONTEXT
-------
HYP-2075 (opus-S562) claims the PINCH sieve -- try t = a/(v_i+v_j) over all
pair-sums with gcd(a, v_i+v_j)=1 -- is COMPLETE for LRC@14: it finds a witness
t with ||v_i t|| >= 1/n for EVERY config. The prior monad-compute session
(lrc_pinch_completeness_stress_s..., lrc_pinch_dense_adversarial_s...) corroborated
this at n=14 across 20072 configs, 0 pinch misses. HANDOFF: "worth checking pinch
completeness at n=15,16,17 (untested here)." This script does exactly that.

LRC@n: n-1 runners (positive integer speeds), observer threshold 1/n. By Tao's
reduction it suffices to consider integer speeds. A config V (a set of n-1 speeds)
is "lonely" at time t iff every runner is >= 1/n from the origin, i.e. there is a
witness t = a/q (gcd(a,q)=1) with ||v t|| >= 1/n for all v in V.

Per config V we do:
  (1) PINCH sieve: search t = a/(v_i+v_j) over distinct pair-sums. Witness?
  (2) GROUND TRUTH: search ALL reduced t = a/q for q up to QMAX. Any witness?

REFUTATION of completeness = a config where (2) finds a witness but (1) does not.
We also flag configs where (2) finds NOTHING up to QMAX (true LRC near-miss /
needs larger q).

All arithmetic INTEGER: for t = a/q, runner v safe iff n*min((v*a)%q, q-(v*a)%q) >= q.
"""

from math import gcd
import random
from functools import reduce


def safe_at(V, a, q, N):
    """Is t = a/q safe (all ||v t|| >= 1/N)?  integer-only."""
    for v in V:
        r = (v * a) % q
        if N * min(r, q - r) < q:
            return False
    return True


def pinch_witness(V, N):
    """PINCH sieve: t = a/(v_i+v_j) over distinct pair-sums (ascending)."""
    L = len(V)
    sums = sorted({V[i] + V[j] for i in range(L) for j in range(i + 1, L)})
    for s in sums:
        for a in range(1, s):
            if gcd(a, s) == 1 and safe_at(V, a, s, N):
                return (a, s)
    return None


def any_witness(V, qmax, N):
    """Ground truth: search ALL reduced fractions a/q, 2<=q<=qmax."""
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) == 1 and safe_at(V, a, q, N):
                return (a, q)
    return None


def primitive(V):
    g = reduce(gcd, V)
    return tuple(sorted(v // g for v in set(V)))


def gen_configs(N, seed):
    """Generate adversarial + random configs of K = N-1 distinct speeds."""
    K = N - 1
    rng = random.Random(seed)
    cfgs = {}

    # --- canonical / structured ---
    cfgs["AP 1..K"] = tuple(range(1, N))                      # the tight extremal set
    # apex: force a multiple of N into the set (apex obstruction)
    base = list(range(1, N))
    base[-1] = N                                              # replace K with N (mult of N)
    cfgs["apex mult-of-N"] = tuple(sorted(set(base))) if len(set(base)) == K else tuple(range(1, N))
    cfgs["LCM-loaded"] = primitive(tuple(list(range(1, K)) + [reduce(lambda a, b: a * b // gcd(a, b), range(2, 8))]))
    cfgs["powers-of-2-ish"] = primitive(tuple([2 ** i for i in range(K // 2)] + list(range(3, 3 + (K - K // 2) * 2, 2))))
    cfgs["geometric-3"] = primitive(tuple([3 ** i for i in range(5)] + list(range(2, 2 + (K - 5)))))
    cfgs["near-AP+jump"] = tuple(sorted(set(list(range(1, K)) + [100])))

    # --- dense near-AP adversarial: single-coordinate perturbations of 1..K ---
    ap = list(range(1, N))
    for idx in range(K):
        for delta in (-1, 1, 2, K):
            pert = ap.copy()
            pert[idx] = ap[idx] + delta
            P = primitive(tuple(pert))
            if len(P) == K and all(x > 0 for x in P):
                cfgs[f"nearAP idx{idx} d{delta}"] = P

    # --- dense small-range: K speeds drawn from {1..hi}, hi just above K ---
    for hi in range(N, N + 6):
        for i in range(120):
            S = tuple(rng.sample(range(1, hi + 1), K))
            P = primitive(S)
            if len(P) == K:
                cfgs[f"dense<= {hi} #{i}"] = P

    # --- random primitive across wider ranges ---
    for rng_hi in (2 * N, 3 * N, 6 * N, 12 * N, 300, 1000):
        for i in range(150):
            while True:
                P = primitive(tuple(rng.sample(range(1, rng_hi + 1), K)))
                if len(P) == K:
                    break
            cfgs[f"rand<= {rng_hi} #{i}"] = P

    return cfgs


def run_for_N(N, seeds, QMAX):
    K = N - 1
    total = 0
    pinch_caught = 0
    refutations = []
    nowit = []
    seen = set()
    for seed in seeds:
        cfgs = gen_configs(N, seed)
        for name, Vraw in cfgs.items():
            V = tuple(sorted(set(Vraw)))
            if len(V) != K:
                continue
            if V in seen:
                continue
            seen.add(V)
            total += 1
            pw = pinch_witness(V, N)
            if pw is not None:
                pinch_caught += 1
            else:
                aw = any_witness(V, QMAX, N)
                if aw is not None:
                    refutations.append((seed, name, V, aw))
                else:
                    nowit.append((seed, name, V))
    print(f"==== PINCH completeness stress test, n={N} ({K} runners) ====")
    print(f"seeds: {seeds}   QMAX(ground truth)={QMAX}")
    print(f"distinct configs tested: {total}")
    print(f"PINCH caught (1/{N} witness via pair-sum): {pinch_caught}/{total} "
          f"({100*pinch_caught/total:.3f}%)")
    print(f"PINCH missed: {total - pinch_caught}")
    if refutations:
        print(f"*** {len(refutations)} REFUTATIONS: pinch MISS but witness exists "
              f"(q<={QMAX}) -> 'pinch complete' FALSE at n={N} ***")
        for seed, name, V, aw in refutations[:25]:
            sums = sorted({V[i]+V[j] for i in range(len(V)) for j in range(i+1, len(V))})
            print(f"   seed={seed} {name}: V={V}")
            print(f"       witness t={aw[0]}/{aw[1]} (q={aw[1]}); max pair-sum={sums[-1]}; "
                  f"q in pair-sums? {aw[1] in sums}")
    else:
        print(f"No refutations: PINCH caught every config that has any witness with q<={QMAX}.")
    if nowit:
        print(f"!! {len(nowit)} configs with NO witness up to q<={QMAX} "
              f"(LRC near-miss flag; pinch also missed):")
        for seed, name, V in nowit[:25]:
            print(f"   seed={seed} {name}: V={V}")
    else:
        print(f"Every config has a 1/{N} witness with q<={QMAX}.")
    print()
    return total, pinch_caught, len(refutations), len(nowit)


def main():
    seeds = [1, 7, 42, 2026, 31337]
    QMAX = 600
    grand = []
    for N in (15, 16, 17):
        grand.append((N,) + run_for_N(N, seeds, QMAX))
    print("==== SUMMARY ====")
    print(f"{'n':>3} {'configs':>9} {'pinch_caught':>13} {'refutations':>12} {'no_witness':>11}")
    for N, tot, caught, refs, now in grand:
        print(f"{N:>3} {tot:>9} {caught:>13} {refs:>12} {now:>11}")
    total_refs = sum(g[3] for g in grand)
    print()
    if total_refs == 0:
        print("VERDICT: PINCH sieve COMPLETE across all tested configs at n=15,16,17 "
              "(0 refutations).")
    else:
        print(f"VERDICT: {total_refs} pinch refutations found -> completeness FAILS at some n in 15,16,17.")
    print("DONE.")


if __name__ == "__main__":
    main()
