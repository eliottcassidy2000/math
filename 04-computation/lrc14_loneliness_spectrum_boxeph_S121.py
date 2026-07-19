#!/usr/bin/env python3
"""
boxeph-2026-07-19-S121 — the loneliness SPECTRUM for 12 speeds: the near-minimizers of
M(C)=max_t min_i ||v_i t|| among bounded-diameter 12-sets, and the rigidity margin.

WHY BOUNDED DIAMETER SUFFICES.  S119: spread configurations are loose (large v_max/v_min =>
large safe arc => M bounded away from 1/13).  So the near-minimizers of loneliness are
COMPACT, and enumerating 12-subsets of {1,...,N} for modest N captures them.

TOOL.  By the PINCH LEMMA (HYP-2059 / THM-401, opus): M(C) is attained at t=m/(v_i+v_j) for a
pairwise SUM; we search all numerators m over all pairwise sums (reduction != representation,
MISTAKE-173 -- so search m/q for q a sum, NOT just reduced denominators).

OUTPUT.
 (A) the smallest ~25 distinct M-values over primitive 12-subsets of {1,...,N}, with an example
     set each -> the loneliness spectrum and the runner-up gap above 1/13;
 (B) characterize the near-minimizers: Hamming distance from {1,...,12}, whether residues mod 13
     are a complete nonzero system, the maximizing pair and its sum mod 13
     (testing the corollary M=1/13 => 13 | (v_i+v_j) at the maximizing pair);
 (C) confirm {1,...,12} is the UNIQUE primitive minimizer at 1/13 in the range.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations


def fd(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def fmin(sp, t: F) -> F:
    return min(fd(F(v) * t) for v in sp)


def M_and_pair(sp):
    """M via pinch lemma: best over t=m/q, q a pairwise sum, all numerators m.
    Return (M, q*, pair(v_i,v_j) summing to q*)."""
    sums = {}
    n = len(sp)
    for i in range(n):
        for j in range(i + 1, n):
            sums.setdefault(sp[i] + sp[j], (sp[i], sp[j]))
    best, bq, bm = F(0), None, None
    for q in sums:
        for m in range(1, q):
            v = fmin(sp, F(m, q))
            if v > best:
                best, bq, bm = v, q, m
    return best, bq, sums[bq]


def primitive(sp):
    g = 0
    for v in sp:
        g = gcd(g, v)
    return g == 1


THIRT = F(1, 13)
BASE = tuple(range(1, 13))


def run(N):
    print("=" * 76, flush=True)
    print(f"LONELINESS SPECTRUM over primitive 12-subsets of {{1,...,{N}}}", flush=True)
    print("=" * 76, flush=True)
    results = []   # (M, set, q, pair)
    total = 0
    for combo in combinations(range(1, N + 1), 12):
        if not primitive(combo):
            continue
        total += 1
        M, q, pair = M_and_pair(list(combo))
        results.append((M, combo, q, pair))
    results.sort(key=lambda r: (r[0], r[1]))
    print(f"  enumerated {total} primitive 12-subsets", flush=True)

    # (A) spectrum: smallest distinct M values
    print("\n  (A) smallest distinct M-values (the spectrum near the minimum):", flush=True)
    print(f"      {'M':>9} {'~dec':>7} {'count':>6}  example set", flush=True)
    seen = {}
    for M, combo, q, pair in results:
        seen.setdefault(M, [0, combo])
        seen[M][0] += 1
    distinct = sorted(seen.items(), key=lambda kv: kv[0])
    for M, (cnt, ex) in distinct[:25]:
        exs = "{" + ",".join(map(str, ex)) + "}"
        print(f"      {str(M):>9} {float(M):>7.4f} {cnt:>6}  {exs}", flush=True)

    # runner-up gap
    m0 = distinct[0][0]
    m1 = distinct[1][0] if len(distinct) > 1 else None
    print(f"\n  minimum M = {m0} ({'== 1/13' if m0 == THIRT else '!!'}) ; "
          f"runner-up M = {m1} ({float(m1):.4f}); gap = {m1 - m0} ~ {float(m1-m0):.5f}", flush=True)

    # (C) uniqueness at the minimum
    at_min = [combo for M, combo, q, pair in results if M == THIRT]
    print(f"\n  (C) primitive 12-subsets of {{1..{N}}} with M == 1/13: {len(at_min)}", flush=True)
    for c in at_min:
        print(f"      {c}  (== {{1..12}}: {c == BASE})", flush=True)

    # (B) characterize the near-minimizers (smallest 12 sets)
    print(f"\n  (B) near-minimizer structure (smallest 12 by M):", flush=True)
    print(f"      {'M':>9} {'Ham':>4} {'complete13':>10} {'pair':>10} {'q':>4} {'q%13':>5}  set", flush=True)
    for M, combo, q, pair in results[:12]:
        ham = len(set(combo) ^ set(BASE)) // 1  # symmetric-difference size
        ham = len(set(combo) - set(BASE))       # #elements changed from {1..12}
        res = sorted(v % 13 for v in combo)
        complete = (res == list(range(0, 13))[1:]) or (len(set(v % 13 for v in combo)) == 12
                                                       and 0 not in (v % 13 for v in combo))
        exs = "{" + ",".join(map(str, combo)) + "}"
        print(f"      {str(M):>9} {ham:>4} {str(complete):>10} {str(pair):>10} {q:>4} {q%13:>5}  {exs}",
              flush=True)
    return distinct, results


d16, r16 = run(16)
print("", flush=True)
d17, r17 = run(17)
print("\nDONE.", flush=True)
