#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_adelic_nonrelation_probe_0614.py

Honest probe of the NON-RELATION baseline: does the single-prime limit beta_p of a
genuine Sidon-ish set (no small exact relations) actually converge to main (6/7)^13,
or sit below it?  This pins down the exact wording of the correction to HYP-2503:
  - if beta_p -> main for true Sidon sets: 'generic local factor = 1 (returns main)';
  - if beta_p < main: there are still SOME exact relations (or slow convergence), and
    the precise statement is 'beta_p = L(S) for every p' (full archimedean limit),
    with main only when L = main (truly relation-free).

Approach:
 1. Build a strict Sidon set (all pairwise sums distinct) -> NO 2-term or simple 3-term
    exact relations with small coefficients. Compare beta_p across primes and vs main.
 2. Push the prime ladder further (p=2 to high e) to separate slow-convergence from a
    genuine sub-main limit.
 3. Count small exact relations in each test set (the actual driver).
"""
import sys, math
from math import gcd
from itertools import combinations, product

sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout, 'reconfigure') else None
MAIN = (6/7) ** 13


def deficit(q, S):
    rad = q // 14
    cnt = 0
    for a in range(q):
        ok = True
        for v in S:
            r = (v * a) % q
            if r <= rad or r >= q - rad:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt


def small_exact_relations(S, tmax=3, maxk=3):
    """count exact relations sum t_v v = 0, |t|<=tmax, |T|<=maxk, t_v!=0 (canonical first>0)."""
    S = list(S); n = len(S); cc = [t for t in range(-tmax, tmax+1) if t != 0]
    cnt = 0; examples = []
    for k in range(2, maxk+1):
        for T in combinations(range(n), k):
            vals = [S[i] for i in T]
            for cs in product(cc, repeat=k):
                if cs[0] > 0 and sum(c*v for c, v in zip(cs, vals)) == 0:
                    cnt += 1
                    if len(examples) < 5:
                        examples.append((vals, cs))
    return cnt, examples


def is_sidon(S):
    sums = set()
    S = sorted(S)
    for i in range(len(S)):
        for j in range(i, len(S)):
            s = S[i] + S[j]
            if s in sums:
                return False
            sums.add(s)
    return True


def main():
    print("=" * 78)
    print(f"NON-RELATION baseline probe.  main (6/7)^13 = {MAIN:.6f}")
    print("=" * 78)

    # A genuine Sidon set of size 13 containing a multiple of 14 (Mian-Chowla-ish), primitive.
    # Mian-Chowla: 1,2,4,8,13,21,31,45,66,81,97,123,148 ... ensure a mult of 14 present.
    mc = [1, 2, 4, 8, 13, 21, 31, 45, 66, 81, 97, 123, 148]
    # replace one element to insert 14 while keeping Sidon-ish; test plain mc first
    sets = {
        "Mian-Chowla (strict Sidon, no mult-14)": sorted(mc),
        "Sidon+14 (swap 13->14)": sorted([14 if x == 13 else x for x in mc]),
        "geometric-ish 2^k+small": sorted(set([14] + [2**k + (k % 3) for k in range(1, 13)]))[:13],
    }
    for name, S in sets.items():
        nrel, ex = small_exact_relations(S)
        sid = is_sidon(S)
        print(f"\n{name}: S={S}")
        print(f"   Sidon? {sid}   small exact relations (|t|<=3,|T|<=3): {nrel}", end="")
        if ex:
            print("  e.g. " + ", ".join(f"{cs}*{vals}" for vals, cs in ex[:2]))
        else:
            print()
        # large-q L
        Lvals = [deficit(q, S)/q for q in (13999, 14000, 15013, 16007, 17011)]
        L = sum(Lvals)/len(Lvals)
        print(f"   L (large-q avg) ~ {L:.5f}   (main {MAIN:.5f}; gap {L-MAIN:+.5f})")
        # deep p=2 ladder to separate slow convergence from genuine sub-main
        print("   p=2 ladder D(2^e)/2^e:", end=" ")
        row = []
        for e in range(5, 17):
            q = 2**e
            row.append(f"{deficit(q,S)/q:.4f}")
        print(" ".join(row))
        # a couple other primes' deep tails
        for p in (3, 5):
            row = []
            e = 1
            while p**e <= 60000:
                if p**e >= 28:
                    row.append(deficit(p**e, S)/(p**e))
                e += 1
            print(f"   p={p} tail: {['%.4f'%x for x in row[-4:]]}")

    print("\n" + "=" * 78)
    print("READING: if a set with 0 small exact relations has beta_p -> main for every p,")
    print("then the generic local factor is trivially 1 (returns main). Any sub-main value")
    print("tracks the count of small exact relations -> confirms L is the exact-relation sum,")
    print("reproduced identically by every single prime (no genuine per-prime factorization).")


if __name__ == "__main__":
    main()
