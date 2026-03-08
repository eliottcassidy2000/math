import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
trace_cycle_theorem.py
kind-pasteur-2026-03-07-S39b

THEOREM: For any tournament T on n vertices:
  tr(A^k) = k * c_k  for k = 3, 5

where A is the adjacency matrix and c_k counts directed k-cycles.

PROOF IDEA: In a tournament, A[i][j]*A[j][i] = 0 for all i != j
(no bidirectional edges). So any closed walk has length >= 3.
A closed walk of length k that revisits a vertex v splits into two
closed sub-walks at v, each of length >= 3. So k >= 6.
For k = 3 or 5: no vertex repetition possible => walk IS a simple cycle.
Each simple directed k-cycle contributes k to tr(A^k) (k starting points).
Hence tr(A^k) = k * c_k.

For k >= 7: the formula FAILS because compound walks exist
(e.g., a 3-cycle and 4-cycle sharing a vertex give a length-7 walk).

CONSEQUENCE: c_5(T) = tr(A^5) / 5 for ANY tournament.
Computing tr(A^5) via matrix multiplication is O(n^3), much faster than
direct cycle enumeration which is O(n^5).

This script provides exhaustive verification.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations
from collections import defaultdict


def count_directed_k_cycles_dp(T, k):
    """Count directed k-cycles using subset DP."""
    n = len(T)
    if k > n or k < 3:
        return 0
    count = 0
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = [[0] * k for _ in range(1 << k)]
        dp[1][0] = 1
        for mask in range(1, 1 << k):
            for last in range(k):
                if dp[mask][last] == 0 or not (mask & (1 << last)):
                    continue
                for nxt in range(1, k):
                    if mask & (1 << nxt):
                        continue
                    if T[v[last]][v[nxt]]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
        full = (1 << k) - 1
        for last in range(1, k):
            if T[v[last]][v[0]]:
                count += dp[full][last]
    return count


def matrix_power_trace(T, k):
    """Compute tr(A^k) via repeated matrix multiplication. O(n^3 * k)."""
    n = len(T)
    Ak = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(k):
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                s = 0
                for l in range(n):
                    s += Ak[i][l] * T[l][j]
                new[i][j] = s
        Ak = new
    return sum(Ak[i][i] for i in range(n))


# ============================================================
# Exhaustive verification
# ============================================================
print("=" * 70)
print("THEOREM: tr(A^k) = k * c_k for k = 3, 5 in tournaments")
print("=" * 70)

# Test for ALL tournaments at n = 3, 4, 5, 6
for n in range(3, 8):
    m = n * (n - 1) // 2
    total = 1 << m
    if n >= 7:
        # Sample for n=7
        import random
        random.seed(42)
        sample = [random.randint(0, total - 1) for _ in range(2000)]
        sample_type = "sampled (2000)"
    else:
        sample = range(total)
        sample_type = "exhaustive"

    k3_match = 0
    k5_match = 0
    k7_match = 0
    tested = 0

    for bits in sample:
        T = tournament_from_bits(n, bits)

        tr3 = matrix_power_trace(T, 3)
        c3 = count_directed_k_cycles_dp(T, 3)
        if tr3 == 3 * c3:
            k3_match += 1

        if n >= 5:
            tr5 = matrix_power_trace(T, 5)
            c5 = count_directed_k_cycles_dp(T, 5)
            if tr5 == 5 * c5:
                k5_match += 1

        if n >= 7 and tested < 200:  # Only check k=7 for small sample
            tr7 = matrix_power_trace(T, 7)
            c7 = count_directed_k_cycles_dp(T, 7)
            if tr7 == 7 * c7:
                k7_match += 1
            elif tested < 5:
                print(f"  k=7 mismatch: bits={bits}, tr(A^7)={tr7}, 7*c_7={7*c7}, "
                      f"diff={tr7-7*c7}")

        tested += 1

    print(f"\nn={n} ({sample_type}, {tested} tournaments):")
    print(f"  k=3: tr(A^3) = 3*c_3: {k3_match}/{tested} "
          f"({'PASS' if k3_match == tested else 'FAIL'})")
    if n >= 5:
        print(f"  k=5: tr(A^5) = 5*c_5: {k5_match}/{tested} "
              f"({'PASS' if k5_match == tested else 'FAIL'})")
    if n >= 7:
        k7_tested = min(200, tested)
        print(f"  k=7: tr(A^7) = 7*c_7: {k7_match}/{k7_tested} "
              f"({'PASS' if k7_match == k7_tested else 'FAIL -- compound walks exist'})")


# ============================================================
# Demonstrate the compound walk that breaks k=7
# ============================================================
print("\n" + "=" * 70)
print("COMPOUND WALKS AT k=7")
print("=" * 70)

# Find a tournament where tr(A^7) != 7*c_7
n = 7
m = n * (n - 1) // 2
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    tr7 = matrix_power_trace(T, 7)
    c7 = count_directed_k_cycles_dp(T, 7)
    if tr7 != 7 * c7:
        c3 = count_directed_k_cycles_dp(T, 3)
        c4 = count_directed_k_cycles_dp(T, 4)
        diff = tr7 - 7 * c7
        print(f"\nExample: bits={bits}")
        print(f"  c_3 = {c3}, c_4 = {c4}, c_7 = {c7}")
        print(f"  tr(A^7) = {tr7}, 7*c_7 = {7*c7}")
        print(f"  Excess = {diff}")
        print(f"  Excess / 7 = {diff/7:.2f}")
        # The excess comes from compound walks: 3-cycle + 4-cycle sharing vertex
        # A compound walk (3,4) contributes 3*4 = 12 walks per pair? No...
        # Actually: for each pair of a directed 3-cycle and directed 4-cycle
        # sharing a vertex v, the compound walk is:
        # v -> (3-cycle vertices) -> v -> (4-cycle vertices) -> v
        # This is one walk per starting rotation, giving 7 starting points.
        # Wait, the walk has 7 edges, so 7 starting points? No, 7 starting points
        # only if the walk visits 7 DISTINCT vertices. The compound walk visits 6.
        # For a compound walk visiting 6 vertices with vertex v repeated:
        # v -> a -> b -> v -> c -> d -> e -> v: 7 edges, 6 distinct vertices.
        # Starting points: can start at any of the 7 positions in the walk
        # (which are: v, a, b, v, c, d, e). But since it's a closed walk,
        # each of the 7 positions gives a rotation. The 7 rotations are:
        # (v,a,b,v,c,d,e), (a,b,v,c,d,e,v), (b,v,c,d,e,v,a), etc.
        # Each is a distinct closed walk contributing 1 to tr(A^7).
        # But v appears twice, so the contribution to tr(A^7) is 7 (not 6).

        # Actually, each starting position contributes to tr(A^7)[start][start].
        # The walk starting at v contributes to [v][v].
        # The walk starting at a contributes to [a][a].
        # So total contribution = 7 (one per starting position).

        # Number of compound (3,4)-walks sharing a vertex:
        # For each vertex v, count pairs (directed 3-cycle through v, directed 4-cycle through v)
        # Each such pair gives 2 compound walks (clockwise or counterclockwise arrangement
        # of the two sub-cycles).
        # Actually: the arrangement is determined. The walk is:
        # v -> 3-cycle (2 edges) -> v -> 4-cycle (3 edges) -> v
        # OR
        # v -> 4-cycle (3 edges) -> v -> 3-cycle (2 edges) -> v
        # These are 2 different walks. Each contributes 7 to tr(A^7).
        # So total excess = 2 * 7 * (#pairs of 3-cycle and 4-cycle sharing vertex).

        # Actually this is more nuanced. Let me just print the excess.
        print(f"  Interpretation: compound walks (3+4) sharing a vertex")
        break

# ============================================================
# Show O(n^3) c_5 computation
# ============================================================
print("\n" + "=" * 70)
print("PRACTICAL: c_5(T) = tr(A^5)/5 via O(n^3) matrix multiplication")
print("=" * 70)

def c5_via_trace(T):
    """Compute c_5(T) = tr(A^5)/5 using matrix multiplication. O(n^3)."""
    n = len(T)
    if n < 5:
        return 0
    # A^2
    A2 = [[sum(T[i][k]*T[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    # A^4 = A^2 * A^2
    A4 = [[sum(A2[i][k]*A2[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    # A^5 = A^4 * A
    A5 = [[sum(A4[i][k]*T[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    tr5 = sum(A5[i][i] for i in range(n))
    assert tr5 % 5 == 0, f"tr(A^5) = {tr5} not divisible by 5!"
    return tr5 // 5


# Verify at n=5,6
for n in [5, 6]:
    m = n * (n - 1) // 2
    all_match = True
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        c5_dp = count_directed_k_cycles_dp(T, 5)
        c5_tr = c5_via_trace(T)
        if c5_dp != c5_tr:
            print(f"  n={n}, bits={bits}: MISMATCH dp={c5_dp}, trace={c5_tr}")
            all_match = False
    if all_match:
        print(f"  n={n}: c5_via_trace matches for ALL {1 << m} tournaments")


# Benchmark: compare timing
import time
n = 7
m = n * (n - 1) // 2
sample_bits = [42, 1000, 5000, 10000, 50000]

print(f"\n  n=7 timing comparison (5 tournaments):")
t0 = time.time()
for bits in sample_bits:
    T = tournament_from_bits(n, bits)
    c5_dp = count_directed_k_cycles_dp(T, 5)
t_dp = time.time() - t0

t0 = time.time()
for bits in sample_bits:
    T = tournament_from_bits(n, bits)
    c5_tr = c5_via_trace(T)
t_tr = time.time() - t0

print(f"    DP method: {t_dp:.4f}s")
print(f"    Trace method: {t_tr:.4f}s")
print(f"    Speedup: {t_dp/t_tr:.1f}x")


# ============================================================
# Also verify Moon's c_3 = C(n,3) - sum C(s_v, 2) = tr(A^3)/3
# ============================================================
print("\n" + "=" * 70)
print("COMBINED: c_3 = tr(A^3)/3, c_5 = tr(A^5)/5 for all tournaments")
print("=" * 70)

from math import comb

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    all_match = True
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        scores = [sum(T[i]) for i in range(n)]

        c3_moon = comb(n, 3) - sum(comb(s, 2) for s in scores)
        tr3 = matrix_power_trace(T, 3)
        c3_trace = tr3 // 3

        if c3_moon != c3_trace:
            print(f"  n={n}, bits={bits}: MISMATCH Moon={c3_moon}, trace={c3_trace}")
            all_match = False
    if all_match:
        print(f"  n={n}: c_3 = tr(A^3)/3 = Moon's formula for ALL {1 << m} tournaments")


print("\nDone.")
