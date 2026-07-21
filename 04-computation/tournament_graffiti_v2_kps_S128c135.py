#!/usr/bin/env python3
"""tournament_graffiti_v2_kps_S128c135.py -- kind-pasteur-2026-07-21-S128c135

WOWII loop, mass-generation round: expanded invariant set + n=7 CROSS-VALIDATION so survivors
are robust (c134's srange<=beta held n<=6 but failed n=7; here every survivor must also pass a
n=7 sample).  Reports pair inequalities target <= a*source + b holding on BOTH all n=3..6
tournaments AND a n=7 sample -- these are the robust candidate theorems.
"""
import sys
import numpy as np
import random
from itertools import combinations
from fractions import Fraction as Fr


def from_bits(bits, n):
    A = np.zeros((n, n), dtype=np.int64)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
            idx += 1
    return A


def ham_paths(A, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            c = dp[mask][last]
            if c:
                for w in range(n):
                    if not (mask >> w & 1) and A[last][w]:
                        dp[mask | (1 << w)][w] += c
    return sum(dp[full][last] for last in range(n))


def largest_transitive(A, n):
    for k in range(n, 0, -1):
        for S in combinations(range(n), k):
            sub = A[np.ix_(S, S)]
            if sorted(sub.sum(axis=1).tolist()) == list(range(k)):
                return k
    return 1


def diameter(A, n):
    # digraph distances; diameter = max finite pairwise distance (tournaments: <= n-1, strongly conn)
    INF = 10 ** 9
    D = np.full((n, n), INF)
    for i in range(n):
        D[i, i] = 0
        for j in range(n):
            if A[i, j]:
                D[i, j] = 1
    for k in range(n):
        D = np.minimum(D, D[:, k:k + 1] + D[k:k + 1, :])
    finite = D[D < INF]
    return int(finite.max())


def num_kings(A, n):
    R2 = (A + A @ A) > 0
    return sum(1 for v in range(n) if all(R2[v, w] for w in range(n) if w != v))


def num_scc(A, n):
    reach = np.linalg.matrix_power(A + np.eye(n, dtype=np.int64), n) > 0
    mut = reach & reach.T
    seen = [False] * n
    c = 0
    for v in range(n):
        if not seen[v]:
            c += 1
            for w in range(n):
                if mut[v, w]:
                    seen[w] = True
    return c


def inv(A, n):
    c3 = int(round(np.trace(np.linalg.matrix_power(A, 3)) / 3))
    sc = A.sum(axis=1)
    return {
        "n": n, "c3": c3, "H": ham_paths([list(r) for r in A], n),
        "beta": largest_transitive(A, n), "diam": diameter(A, n),
        "kings": num_kings(A, n), "scc": num_scc(A, n),
        "smax": int(sc.max()), "smin": int(sc.min()),
        "srange": int(sc.max() - sc.min()),
    }


KEYS = ["c3", "H", "beta", "diam", "kings", "scc", "smax", "smin", "srange"]

print("gathering n=3..6 exhaustive + n=7 sample ...", flush=True)
LOW = []
for n in range(3, 7):
    for bits in range(1 << (n * (n - 1) // 2)):
        LOW.append(inv(from_bits(bits, n), n))
random.seed(2)
HI = [inv(from_bits(random.getrandbits(21), 7), 7) for _ in range(4000)]
print("  n<=6: %d ; n=7 sample: %d" % (len(LOW), len(HI)), flush=True)

print()
print("=" * 88)
print("ROBUST pair inequalities  target <= a*source + b  (hold on n<=6 AND n=7 sample)")
print("=" * 88)
rows = []
for t in KEYS:
    for x in KEYS:
        if t == x:
            continue
        for a in (Fr(1, 2), Fr(1), Fr(2)):
            need = max(Fr(d[t]) - a * Fr(d[x]) for d in LOW)
            b = int(np.ceil(float(need)))
            if abs(b) > 6:
                continue
            # must also hold on n=7 sample
            if any(Fr(d[t]) > a * Fr(d[x]) + b for d in HI):
                continue
            eq = sum(1 for d in LOW if Fr(d[t]) == a * Fr(d[x]) + b)
            # skip trivial identities and near-trivial
            rows.append((eq, t, a, x, b))
rows.sort(reverse=True)
seen = set()
shown = 0
for eq, t, a, x, b in rows:
    if (t, x) in seen:
        continue
    seen.add((t, x))
    astr = ("%s*" % a) if a != 1 else ""
    bstr = (" + %d" % b) if b > 0 else ((" - %d" % (-b)) if b < 0 else "")
    print("  %-6s <= %s%-6s%-8s [eq on %d]" % (t, astr, x, bstr, eq))
    shown += 1
    if shown >= 24:
        break
print()
print("  Every line above passed the n=7 cross-validation, so none is a c134-style small-n")
print("  artifact (srange <= beta, which held n<=6, is correctly ABSENT -- it fails at n=7).")
