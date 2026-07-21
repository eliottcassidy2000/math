#!/usr/bin/env python3
"""wowii_candidates_largen_kps_S128c135.py -- kind-pasteur-2026-07-21-S128c135

WOWII LOOP, refute/verify step: push the two candidate inequalities from THM-1845 past n=7.
  (C1) c3 <= H           (3-cycles <= Hamiltonian paths)
  (C2) H <= 2^{n-2}*c3+1 (tight at transitive: 0-cycles, H=1)

Heuristic: max H over tournaments ~ n!/2^{n-1} (grows super-exponentially) while max c3 ~
n^3/24, so 2^{n-2}*c3 ~ 2^{n-2} n^3 is eventually BEATEN by H -- (C2) should fail at some n.
(C1) should keep holding (regular tournaments have H huge, c3 = O(n^3)).  Compute both on
structured HIGH-H families (Paley p=3 mod 4, quadratic-residue rotational, random near-regular)
for n up to ~15, and report the first violation of each.
"""
import sys
import numpy as np
import random
from math import comb


def ham_paths(A, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c:
                Al = A[last]
                for w in range(n):
                    if not (mask >> w & 1) and Al[w]:
                        dp[mask | (1 << w)][w] += c
    return sum(dp[full][last] for last in range(n))


def c3_count(A, n):
    M = np.array(A)
    return int(round(np.trace(np.linalg.matrix_power(M, 3)) / 3))


def rotational(n, S):
    """circulant tournament: i->j iff (j-i) mod n in S. S a 'sign' set (Sidon-like)."""
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A


def paley(p):
    QR = set((k * k) % p for k in range(1, p))
    return rotational(p, QR)


def regular_random(n, seed):
    """random near-regular tournament (scores as equal as possible)."""
    random.seed(seed)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


print("=" * 78)
print("candidate inequalities at large n on high-H structured tournaments")
print("=" * 78)
print("  %-22s %-3s %-10s %-14s %-14s %-6s %-6s"
      % ("family", "n", "c3", "H", "2^{n-2}c3+1", "C1", "C2"))
first_fail_C1 = None
first_fail_C2 = None
fams = []
for p in (7, 11, 19, 23):
    fams.append(("Paley-%d" % p, p, paley(p)))
# rotational on odd n with S = first (n-1)/2 residues (a regular tournament)
for n in (9, 13, 15):
    S = set(range(1, (n) // 2 + 1))
    fams.append(("rotational-%d" % n, n, rotational(n, S)))
for n in (8, 10, 12, 14):
    fams.append(("near-regular-%d" % n, n, regular_random(n, n)))
for name, n, A in fams:
    c3 = c3_count(A, n)
    H = ham_paths(A, n)
    bound = (1 << (n - 2)) * c3 + 1
    C1 = c3 <= H
    C2 = H <= bound
    if not C1 and first_fail_C1 is None:
        first_fail_C1 = (name, n)
    if not C2 and first_fail_C2 is None:
        first_fail_C2 = (name, n)
    print("  %-22s %-3d %-10d %-14d %-14d %-6s %-6s"
          % (name, n, c3, H, bound, C1, C2))
    sys.stdout.flush()
print()
print("  C1 (c3 <= H)          first violation: %s" % (first_fail_C1 or "NONE (holds)"))
print("  C2 (H <= 2^{n-2}c3+1) first violation: %s" % (first_fail_C2 or "NONE (holds)"))
print()
print("  If C2 fails at moderate n, it is a WOWII off-by-scale bound: true small n, false large")
print("  n -- exactly the loop's 'the pattern breaks' outcome, found by pushing past the filter.")
