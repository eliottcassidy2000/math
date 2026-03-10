import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
h_from_trace_n8.py
kind-pasteur-2026-03-07-S39b

Complete H(T) computation at n=8 using trace-based formulas.

H(T) = 1 + 2*alpha_1 + 4*alpha_2  (alpha_3 = 0 at n <= 8)

alpha_1 = c_3 + c_5 + c_7
  - c_3 via Moon's formula: O(n^2)
  - c_5 via tr(A^5)/5: O(n^3)
  - c_7 via... what? At n=8, there are C(8,7)=8 seven-vertex subsets.
    For each, compute c_7(T[S]) via tr(A_S^7)/7. Total O(8 * 7^3) = O(2744).

alpha_2 = alpha_2(3-cycles) + alpha_2(cross 3x5)
  - alpha_2(3-cycles) via THM-097: O(n^3)
  - alpha_2(cross) via sum_{5-subsets} c_5(S)*c_3(S^c): O(n^5 * 5^3)

This gives a FULLY POLYNOMIAL algorithm for H(T) at n=8.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from tournament_fast import c3_from_score, c5_fast, alpha2_from_trace
from itertools import combinations
from math import comb
import random
import time


def c7_via_subtournaments(T):
    """Compute c_7(T) by summing c_7 over all 7-vertex subsets.

    For each 7-subset S, c_7(T[S]) = tr(A_S^7)/7.
    Total 7-cycles = sum over all subsets (each cycle counted once).
    """
    n = len(T)
    if n < 7:
        return 0

    total_c7 = 0
    for subset in combinations(range(n), 7):
        S = list(subset)
        # Build submatrix
        Ts = [[T[S[i]][S[j]] for j in range(7)] for i in range(7)]
        # Compute tr(A^7) via matrix multiplication
        Ak = [[int(i == j) for j in range(7)] for i in range(7)]
        for _ in range(7):
            new = [[0]*7 for _ in range(7)]
            for i in range(7):
                for j in range(7):
                    for l in range(7):
                        new[i][j] += Ak[i][l] * Ts[l][j]
            Ak = new
        tr7 = sum(Ak[i][i] for i in range(7))
        # At k=7 on 7 vertices: the only non-simple walks would be (3,4) compound.
        # But 3+4-1=6 <= 7, so these CAN exist! So tr(A^7) != 7*c_7 at n=7.
        # Hmm, this means we can't use tr(A^7)/7 directly. Need the correction.
        # Actually wait: THM-118 says tr(A^k) = k*c_k for k=3,4,5 ONLY.
        # At k=7 on 7 vertices, the correction is needed.

        # Use subset DP instead
        dp = [[0]*7 for _ in range(1 << 7)]
        dp[1][0] = 1
        for mask in range(1, 1 << 7):
            for last in range(7):
                if dp[mask][last] == 0 or not (mask & (1 << last)):
                    continue
                for nxt in range(1, 7):
                    if mask & (1 << nxt):
                        continue
                    if Ts[last][nxt]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
        full = (1 << 7) - 1
        c7_sub = 0
        for last in range(1, 7):
            if Ts[last][0]:
                c7_sub += dp[full][last]
        total_c7 += c7_sub

    return total_c7


def alpha2_cross_35(T):
    """Compute # vertex-disjoint (5-cycle, 3-cycle) pairs.

    For each 5-vertex subset S: compute c_5(T[S]) and c_3(T[S^c]).
    Cross pairs = sum c_5(S) * c_3(S^c).
    """
    n = len(T)
    if n < 8:
        return 0

    total = 0
    for subset in combinations(range(n), 5):
        S = list(subset)
        Sc = [v for v in range(n) if v not in subset]

        # c_5 on S via trace
        Ts = [[T[S[i]][S[j]] for j in range(5)] for i in range(5)]
        A2 = [[sum(Ts[i][k]*Ts[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        A4 = [[sum(A2[i][k]*A2[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        A5 = [[sum(A4[i][k]*Ts[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        c5_S = sum(A5[i][i] for i in range(5)) // 5

        if c5_S == 0:
            continue

        # c_3 on S^c
        if len(Sc) == 3:
            a, b, c = Sc
            c3_Sc = 1 if (T[a][b] and T[b][c] and T[c][a]) or (T[a][c] and T[c][b] and T[b][a]) else 0
        else:
            c3_Sc = c3_from_score([[T[Sc[i]][Sc[j]] for j in range(len(Sc))] for i in range(len(Sc))])

        total += c5_S * c3_Sc

    return total


def h_full_trace(T):
    """Compute H(T) using full trace-based formulas. Valid for n <= 8."""
    n = len(T)
    if n <= 1:
        return 1

    # alpha_1
    c3 = c3_from_score(T)
    c5 = c5_fast(T) if n >= 5 else 0
    c7 = c7_via_subtournaments(T) if n >= 7 else 0
    alpha_1 = c3 + c5 + c7

    # alpha_2
    a2_3 = alpha2_from_trace(T) if n >= 6 else 0
    a2_cross = alpha2_cross_35(T) if n >= 8 else 0
    alpha_2 = a2_3 + a2_cross

    return 1 + 2 * alpha_1 + 4 * alpha_2


# ============================================================
# Verification
# ============================================================
print("=" * 70)
print("H(T) via full trace formulas at n=8")
print("=" * 70)

n = 8
m = n * (n - 1) // 2
random.seed(42)
sample = [random.randint(0, (1 << m) - 1) for _ in range(100)]

matches = 0
tested = 0
t0 = time.time()

for bits in sample:
    T = tournament_from_bits(n, bits)
    h_true = hamiltonian_path_count(T)
    h_trace = h_full_trace(T)

    if h_true == h_trace:
        matches += 1
    elif tested < 5:
        print(f"  MISMATCH bits={bits}: true={h_true}, trace={h_trace}")
    tested += 1

t_trace = time.time() - t0

print(f"\nn=8 (100 sampled): H(T) via trace: {matches}/{tested} "
      f"({'PASS' if matches == tested else 'FAIL'})")
print(f"  Time: {t_trace:.2f}s ({t_trace/tested:.4f}s per tournament)")

# Timing comparison
random.seed(42)
sample2 = [random.randint(0, (1 << m) - 1) for _ in range(100)]

t0 = time.time()
for bits in sample2:
    T = tournament_from_bits(n, bits)
    h_dp = hamiltonian_path_count(T)
t_dp = time.time() - t0

print(f"  DP time: {t_dp:.2f}s ({t_dp/100:.4f}s per tournament)")
print(f"  Speedup: {t_dp/max(t_trace, 0.001):.1f}x")

print("\nDone.")
