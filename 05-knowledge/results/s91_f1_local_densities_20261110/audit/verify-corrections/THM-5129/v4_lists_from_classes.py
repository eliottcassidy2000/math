#!/usr/bin/env python3
"""Verifier of the corrections to THM-5129 (2026-11-11): the 35 slope lists for D <= 73 recomputed with own code.

Input: the session's class representatives reps_le8.npz (orders 3..8; the class counts 1, 1, 2, 4, 12, 56, 456, 6880
were confirmed by two auditors with their own generators).  For every strong class with H <= 73 the path count H,
the cover counts U at every root and Pi_2 are recomputed here by subset dynamic programming (own code, module
v2_small_cases), the product law (A2) is applied, and the result is compared with slope_lists_le73.txt.
Also printed: the facts behind binding decision D6 (9-grid, level 6) and the dyadic slopes of (B) consequence 3.
"""
import io
import sys

import numpy as np

import v2_small_cases as v2

ROOT = "C:/tmp/wt-s91/"
RS = ROOT + "04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/"
LIM = 73


def out_masks(code, n):
    out = [0] * n
    for j in range(1, n):
        base = j * (j - 1) // 2
        for i in range(j):
            if (code >> (base + i)) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    return out


def hcount(out, n):
    """number of Hamiltonian paths (ends-only DP, faster than the full tables)."""
    N = 1 << n
    inn = [0] * n
    for i in range(n):
        for j in range(n):
            if out[i] >> j & 1:
                inn[j] |= 1 << i
    ends = [None] * N
    for S in range(1, N):
        row = [0] * n
        for v in range(n):
            if S >> v & 1:
                R = S ^ (1 << v)
                if R == 0:
                    row[v] = 1
                else:
                    m = inn[v] & R
                    e = 0
                    er = ends[R]
                    while m:
                        u = (m & -m).bit_length() - 1
                        e += er[u]
                        m &= m - 1
                    row[v] = e
        ends[S] = row
    return sum(ends[N - 1])


def main():
    z = np.load(RS + "reps_le8.npz")
    counts = []
    root = {}   # H -> set of U_C(y)
    suf1 = {}   # H -> set of Pi_2(C)
    nstrong = []
    m = {}
    for n in range(3, 9):
        codes = [int(c) for c in z["n%d" % n]]
        counts.append(len(codes))
        ns = 0
        for code in codes:
            out = out_masks(code, n)
            if len(v2.strong_components(out, n)) != 1:
                continue
            ns += 1
            h = hcount(out, n)
            if h > LIM:
                continue
            tab = v2.tables(out, n)
            assert tab[0][(1 << n) - 1] == h
            m[h] = m.get(h, 0) + 1
            suf1.setdefault(h, set()).add(v2.pi2(out, n, tab))
            for y in range(n):
                D, E, a, b = v2.stats(out, n, y, tab)
                root.setdefault(h, set()).add(a + b)
        nstrong.append(ns)
    print("classes of orders 3..8:", counts, "; strong:", nstrong)
    strongH = sorted(root)
    att = {1}
    changed = True
    while changed:
        changed = False
        for a in list(att):
            for s in strongH:
                if a * s <= LIM and a * s not in att:
                    att.add(a * s)
                    changed = True
    print("odd numbers <= 73 not attained:", [d for d in range(1, LIM + 1, 2) if d not in att])
    suf = {1: {1}}
    for d in range(3, LIM + 1, 2):
        res = set()
        for s in strongH:
            if d % s == 0 and d // s in suf:
                for p in suf1[s]:
                    for w in suf[d // s]:
                        res.add(p * w)
        if res:
            suf[d] = res
    lists = {}
    for D in sorted(att):
        prim = set()
        for d1 in att:
            if D % d1:
                continue
            rem = D // d1
            for hh in [1] + strongH:
                if rem % hh or rem // hh not in suf:
                    continue
                for u in ([1] if hh == 1 else root[hh]):
                    for w in suf[rem // hh]:
                        prim.add(d1 * u * w)
        odd = {}
        for val in prim:
            j = 0
            while val % 2 == 0:
                val //= 2
                j += 1
            if val not in odd or j < odd[val]:
                odd[val] = j
        lists[D] = odd
    sess = {}
    with io.open(RS + "slope_lists_le73.txt", "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            head, tail = line.split(": ", 1)
            sess[int(head)] = {int(t.split(":")[0]): int(t.split(":")[1]) for t in tail.split()}
    same = lists == sess
    print("offsets:", len(lists), "; pairs (offset, odd part):", sum(len(v) for v in lists.values()),
          "; identical to slope_lists_le73.txt:", same)
    print("strong classes with H = h, h <= 73 with exactly two:", [h for h in sorted(m) if m[h] == 2],
          "; m(67) =", m.get(67), "; m(63) =", m.get(63))
    print("9-grid: least exponent of d = 9 at offsets 5, 9, 17, 33, 65:", [lists[D].get(9) for D in (5, 9, 17, 33, 65)],
          "; (65, 72) a slope:", lists[65].get(9, 99) <= 3, "; (65, 144):", lists[65].get(9, 99) <= 4,
          "; (65, 288):", lists[65].get(9, 99) <= 5)
    print("dyadic slopes (d = 1) with least exponent, attained offsets <= 73:",
          {D: lists[D][1] for D in sorted(lists) if 1 in lists[D]})
    print("offset 55: least exponent of d = 11:", lists[55].get(11))
    ok = same and counts == [2, 4, 12, 56, 456, 6880] and nstrong == [1, 1, 6, 35, 353, 6008]
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
