import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
h_polynomial_n9.py
kind-pasteur-2026-03-07-S39b

Can we compute H(T) at n=9 in polynomial time using trace formulas?

The OCF gives H(T) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3.

Complexity analysis at n=9:
- alpha_1 = c3 + c5 + c7 + c9
  - c3: O(n^2) via Moon
  - c5: O(n^3) via tr(A^5)/5
  - c7: O(C(n,7) * 7^3) = O(n^7 * 7^3/7!) ~ O(n^7/14) by subtournament trace
    BUT: tr(A^7)/7 != c7 due to compound walks. Need subset DP on 7-vertex subs.
    O(C(9,7) * 2^7 * 7) = O(36 * 128 * 7) = O(32256) -- FAST!
  - c9: O(2^9 * 9) = O(4608) via full DP -- FAST!

- alpha_2 = alpha_2(3,3) + alpha_2(3,5)
  - alpha_2(3,3): O(n^3) via THM-097
  - alpha_2(3,5): O(C(n,5) * 5^3) = O(n^5/120 * 125) ~ O(n^5) at general n
    At n=9: C(9,5)*125 = 126*125 = 15750 -- FAST!

- alpha_3 = # vertex-disjoint 3-cycle triples covering all 9 vertices
  This needs checking all C(9,3)*C(6,3)/6 = 84*20/6 = 280 partitions
  For each: check if all three triples are 3-cycles. O(280) -- TRIVIAL!

Total: O(n^5) at n=9, dominated by alpha_2(3,5).
At general n: O(n^5) if alpha_max stays small (2 or 3).

But wait -- at n=10: alpha_2 gets (5,5) pairs too! And alpha_3 gets (3,3,3) on 9 of 10 verts.
Let's compute alpha_2(5,5) complexity at n=10: C(10,5)*5^3*C(5,5)*5^3/2 -- nope,
actually (5,5) pairs need iterating over disjoint 5-subset pairs: C(10,5)*C(5,5)/2 = 126,
with O(5^3) work each = 126*125 = 15750. STILL fast.

This script verifies H(T) via polynomial trace formulas at n=9 with TIMING.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from tournament_fast import c3_from_score, c5_fast, alpha2_from_trace
from itertools import combinations
from math import comb
import random
import time


def c7_subtournament(T):
    """Count c7 by summing Ham cycles over all 7-vertex subtournaments."""
    n = len(T)
    if n < 7:
        return 0
    total = 0
    for subset in combinations(range(n), 7):
        S = list(subset)
        # Subset DP for Hamiltonian cycles starting at S[0]
        dp = [[0]*7 for _ in range(1 << 7)]
        dp[1][0] = 1
        for mask in range(1, 1 << 7):
            for last in range(7):
                if dp[mask][last] == 0 or not (mask & (1 << last)):
                    continue
                for nxt in range(1, 7):
                    if mask & (1 << nxt):
                        continue
                    if T[S[last]][S[nxt]]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
        full = (1 << 7) - 1
        for last in range(1, 7):
            if T[S[last]][S[0]]:
                total += dp[full][last]
    return total


def c9_ham_cycles(T):
    """Count directed Hamiltonian 9-cycles via subset DP."""
    n = len(T)
    if n < 9:
        return 0
    dp = [[0]*n for _ in range(1 << n)]
    dp[1][0] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if dp[mask][last] == 0 or not (mask & (1 << last)):
                continue
            for nxt in range(1, n):
                if mask & (1 << nxt):
                    continue
                if T[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    full = (1 << n) - 1
    total = 0
    for last in range(1, n):
        if T[last][0]:
            total += dp[full][last]
    return total


def alpha2_cross_35_n9(T):
    """Compute alpha_2(3,5) cross terms at n=9."""
    n = len(T)
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

        # c_3 on S^c (4 vertices) via Moon
        T_sc = [[T[Sc[i]][Sc[j]] for j in range(len(Sc))] for i in range(len(Sc))]
        scores_sc = [sum(row) for row in T_sc]
        c3_Sc = comb(len(Sc), 3) - sum(comb(s, 2) for s in scores_sc)

        total += c5_S * c3_Sc

    return total


def alpha3_triples(T):
    """Count alpha_3 = triples of mutually vertex-disjoint 3-cycles.
    At n=9: must use all 9 vertices (3+3+3=9).
    """
    n = len(T)
    if n < 9:
        return 0

    # Enumerate all ways to partition [9] into three triples
    # Fix first element in first triple to avoid overcounting
    count = 0
    verts = list(range(n))

    # For each partition of {0,...,8} into three triples:
    # First triple contains 0, pick 2 more from {1,...,8}: C(8,2) = 28
    for pair1 in combinations(range(1, n), 2):
        t1 = (0,) + pair1
        remaining = [v for v in range(1, n) if v not in pair1]
        # Pick smallest remaining for second triple, then 2 more
        for pair2 in combinations(remaining[1:], 2):
            t2 = (remaining[0],) + pair2
            t3 = tuple(v for v in remaining if v != remaining[0] and v not in pair2)

            # Check all three are directed 3-cycles
            is_cycle = True
            for tri in [t1, t2, t3]:
                a, b, c = tri
                if not ((T[a][b] and T[b][c] and T[c][a]) or
                        (T[a][c] and T[c][b] and T[b][a])):
                    is_cycle = False
                    break

            if is_cycle:
                count += 1

    return count


def h_trace_n9(T):
    """Compute H(T) via trace formulas at n=9. All polynomial time."""
    n = len(T)
    assert n == 9

    c3 = c3_from_score(T)
    c5 = c5_fast(T)
    c7 = c7_subtournament(T)
    c9 = c9_ham_cycles(T)
    alpha_1 = c3 + c5 + c7 + c9

    a2_33 = alpha2_from_trace(T)
    a2_35 = alpha2_cross_35_n9(T)
    alpha_2 = a2_33 + a2_35

    alpha_3 = alpha3_triples(T)

    return 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3


# ============================================================
print("=" * 70)
print("H(T) via POLYNOMIAL trace formulas at n=9")
print("=" * 70)

n = 9
m = n * (n - 1) // 2
random.seed(42)
sample = [random.randint(0, (1 << m) - 1) for _ in range(200)]

matches = 0
t_trace_total = 0
t_dp_total = 0

for idx, bits in enumerate(sample):
    T = tournament_from_bits(n, bits)

    t0 = time.time()
    h_trace = h_trace_n9(T)
    t_trace = time.time() - t0

    t0 = time.time()
    h_true = hamiltonian_path_count(T)
    t_dp = time.time() - t0

    if h_trace == h_true:
        matches += 1
    elif idx < 5:
        print(f"  MISMATCH bits={bits}: trace={h_trace}, true={h_true}")

    t_trace_total += t_trace
    t_dp_total += t_dp

print(f"\nn=9 ({len(sample)} sampled): trace formula: {matches}/{len(sample)} "
      f"({'PASS' if matches == len(sample) else 'FAIL'})")
print(f"  Trace time: {t_trace_total:.2f}s ({t_trace_total/len(sample):.4f}s per)")
print(f"  DP time: {t_dp_total:.2f}s ({t_dp_total/len(sample):.4f}s per)")
print(f"  Ratio (trace/DP): {t_trace_total/max(t_dp_total, 0.001):.2f}x")

# ============================================================
# Complexity breakdown
# ============================================================
print(f"\n" + "=" * 70)
print(f"COMPLEXITY BREAKDOWN at n=9")
print(f"=" * 70)

T_test = tournament_from_bits(n, sample[0])

t0 = time.time()
for _ in range(10):
    c3_from_score(T_test)
t_c3 = (time.time() - t0) / 10

t0 = time.time()
for _ in range(10):
    c5_fast(T_test)
t_c5 = (time.time() - t0) / 10

t0 = time.time()
for _ in range(10):
    c7_subtournament(T_test)
t_c7 = (time.time() - t0) / 10

t0 = time.time()
for _ in range(10):
    c9_ham_cycles(T_test)
t_c9 = (time.time() - t0) / 10

t0 = time.time()
for _ in range(10):
    alpha2_from_trace(T_test)
t_a2_33 = (time.time() - t0) / 10

t0 = time.time()
for _ in range(10):
    alpha2_cross_35_n9(T_test)
t_a2_35 = (time.time() - t0) / 10

t0 = time.time()
for _ in range(10):
    alpha3_triples(T_test)
t_a3 = (time.time() - t0) / 10

print(f"  c3 (Moon):           {t_c3*1000:.3f} ms")
print(f"  c5 (trace):          {t_c5*1000:.3f} ms")
print(f"  c7 (subtournament):  {t_c7*1000:.3f} ms")
print(f"  c9 (full DP):        {t_c9*1000:.3f} ms")
print(f"  alpha2(3,3) (trace): {t_a2_33*1000:.3f} ms")
print(f"  alpha2(3,5) (cross): {t_a2_35*1000:.3f} ms")
print(f"  alpha3 (partition):  {t_a3*1000:.3f} ms")
total = t_c3 + t_c5 + t_c7 + t_c9 + t_a2_33 + t_a2_35 + t_a3
print(f"  TOTAL:               {total*1000:.3f} ms")
print(f"  DP comparison:       {t_dp_total/len(sample)*1000:.3f} ms")

# Bottleneck analysis
parts = [
    ("c3", t_c3), ("c5", t_c5), ("c7", t_c7), ("c9", t_c9),
    ("a2_33", t_a2_33), ("a2_35", t_a2_35), ("a3", t_a3)
]
parts.sort(key=lambda x: -x[1])
print(f"\n  Bottleneck ranking:")
for name, t in parts:
    print(f"    {name}: {100*t/total:.1f}%")

print("\nDone.")
