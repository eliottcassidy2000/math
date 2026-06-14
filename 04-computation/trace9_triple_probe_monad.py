#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tr(A^9): the first TRIPLE-cycle config  (monad-explorer-2026-06-14).

Partitions of 9 into parts >= 3:  {9}, {3,6}, {4,5}, {3,3,3}.  The first three are
the by-now-familiar simple-cycle + distinct-pair terms (coeff 9). The {3,3,3} part
is NEW: the census gains a genuine TRIPLE term for the first time at k=9.

We peel off the known pair-level terms and study the residual
   R = tr(A^9) - 9 c9 - 9*O(3,6) - 9*O(4,5)
to identify the triple-triangle structure. Candidate carriers:
   c3                         : a single triangle traversed 3x (period 3 -> coeff 3)
   pdt = ordered overlapping triangle pairs (X twice + Y once, X~Y), coeff ?
   tt3 = connected triples of DISTINCT triangles (pairwise/chain overlap), coeff ?
"""

import sys, itertools, random
import numpy as np
from collections import defaultdict
from fractions import Fraction

random.seed(99)


def random_tournament(n):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if random.getrandbits(1):
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A


def cycles(A, n, maxL):
    adj = [[j for j in range(n) if A[i, j]] for i in range(n)]
    out = defaultdict(list)
    for start in range(n):
        stack = [(start, 1 << start, (start,))]
        while stack:
            v, vis, path = stack.pop()
            for w in adj[v]:
                if w == start and len(path) >= 3:
                    out[len(path)].append(frozenset(path))
                elif w > start and not (vis >> w) & 1 and len(path) < maxL:
                    stack.append((w, vis | (1 << w), path + (w,)))
    return out


def opairs(la, lb, same):
    cnt = 0
    if same:
        for i in range(len(la)):
            for j in range(i + 1, len(la)):
                if la[i] & la[j]:
                    cnt += 1
    else:
        for x in la:
            for y in lb:
                if x & y:
                    cnt += 1
    return cnt


def triple_features(tris):
    """tris: list of frozenset triangles. Return triple-config counts."""
    nt = len(tris)
    # ordered overlapping pairs X!=Y with X&Y nonempty  (X traversed twice, Y once)
    pdt = 0
    for i in range(nt):
        for j in range(nt):
            if i != j and tris[i] & tris[j]:
                pdt += 1
    # connected triples of distinct triangles (overlap graph on the 3 is connected)
    tt3 = 0
    for c in itertools.combinations(range(nt), 3):
        a, b, d = c
        e_ab = bool(tris[a] & tris[b]); e_ac = bool(tris[a] & tris[d]); e_bc = bool(tris[b] & tris[d])
        deg = e_ab + e_ac + e_bc
        if deg >= 2 or (deg == 1 and False):  # connected on 3 nodes needs >=2 edges
            tt3 += 1
    return pdt, tt3


def solve_rational(Am, bv):
    n = len(Am)
    M = [[Am[i][j] for j in range(n)] + [bv[i]] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] != 0), None)
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                fac = M[r][col]
                M[r] = [M[r][j] - fac * M[col][j] for j in range(n + 1)]
    return [M[i][n] for i in range(n)]


def main():
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    NS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rows = []
    for _ in range(NS):
        A = random_tournament(n)
        cyc = cycles(A, n, 9)
        tris = cyc[3]
        c9 = len(cyc[9])
        O36 = opairs(tris, cyc[6], False)
        O45 = opairs(cyc[4], cyc[5], False)
        c3 = len(tris)
        pdt, tt3 = triple_features(tris)
        tr9 = int(np.trace(np.linalg.matrix_power(A.astype(object), 9)))
        R = tr9 - 9 * c9 - 9 * O36 - 9 * O45
        rows.append(dict(tr9=tr9, c9=c9, O36=O36, O45=O45, c3=c3, pdt=pdt, tt3=tt3, R=R))

    print(f"=== tr(A^9) triple-config probe  (n={n}, {NS} tournaments) ===", flush=True)
    print("  R := tr A^9 - 9 c9 - 9 O(3,6) - 9 O(4,5)   (should be pure {3,3,3} triple term)", flush=True)
    cols = ['c3', 'pdt', 'tt3']
    M = [[Fraction(int(d[c])) for c in cols] for d in rows]
    y = [Fraction(int(d['R'])) for d in rows]
    K = len(cols)
    Mt = [[sum(M[r][a] * M[r][b] for r in range(NS)) for b in range(K)] for a in range(K)]
    Mty = [sum(M[r][a] * y[r] for r in range(NS)) for a in range(K)]
    sol = solve_rational(Mt, Mty)
    if sol:
        ok = all(sum(M[r][a] * sol[a] for a in range(K)) == y[r] for r in range(NS))
        print(f"  exact LS fit  R = {dict(zip(cols, sol))}   exact_on_all={ok}", flush=True)
    else:
        print("  singular normal equations", flush=True)
    # clean candidate: R = 3 c3 + 9 pdt + 9 tt3 ?  (single tripled coeff 3, others coeff 9)
    for name, fn in [
        ("R = 3 c3 + 9 pdt + 9 tt3", lambda d: 3 * d['c3'] + 9 * d['pdt'] + 9 * d['tt3']),
        ("R = 3 c3 + 9 pdt + 18 tt3", lambda d: 3 * d['c3'] + 9 * d['pdt'] + 18 * d['tt3']),
    ]:
        bad = sum(1 for d in rows if d['R'] != fn(d))
        print(f"  test {name:30s}: exact {NS-bad}/{NS}", flush=True)
    print("\n  sample (R, c3, pdt, tt3):", flush=True)
    for d in rows[:8]:
        print(f"    R={d['R']} c3={d['c3']} pdt={d['pdt']} tt3={d['tt3']}  (c9={d['c9']} O36={d['O36']} O45={d['O45']})", flush=True)


if __name__ == "__main__":
    main()
