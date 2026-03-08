import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
trace_cycle_k4.py
kind-pasteur-2026-03-07-S39b

THEOREM (extension of THM-096): tr(A^k) = k * c_k for k = 3, 4, 5 in tournaments.

The k=4 case was MISSED in the original THM-096. The proof is the same:
- Any non-simple closed walk of length k splits at a repeated vertex into
  two sub-walks, each a closed walk of length >= 2.
- In tournaments, there are no 2-cycles (bidirectional edges), so each
  sub-walk has length >= 3.
- For k <= 5: sub-walks would need lengths summing to k with each >= 3,
  impossible since 3+3=6 > 5.
- For k = 6: first failure (3+3 compound possible).
- For k = 7: 3+4 compound.

Wait, this argument is actually WRONG for k=4!
A closed walk of length 4 that repeats a vertex w at positions i and j:
  - If j-i = 1: sub-walks of length 1 and 3. Length 1 = self-loop, impossible.
  - If j-i = 2: sub-walks of length 2 and 2. Both are 2-cycles, impossible.
  - If j-i = 3: sub-walks of length 3 and 1. Length 1 impossible.
So for k=4: NO vertex repetition possible. Every closed walk is a simple cycle.
Hence tr(A^4) = 4 * c_4. QED.

More generally: for ANY k < 6, non-simple closed walks are impossible.
The minimum compound walk length is 3+3=6 (two 3-cycles sharing a vertex).

So tr(A^k) = k * c_k for k = 3, 4, 5 in ANY tournament.

This also means c_4(T) = tr(A^4)/4, giving O(n^3) computation of 4-cycles too!

VERIFICATION BELOW.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations
import random

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
    """Compute tr(A^k) via matrix multiplication."""
    n = len(T)
    Ak = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(k):
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new[i][j] += Ak[i][l] * T[l][j]
        Ak = new
    return sum(Ak[i][i] for i in range(n))


def c4_via_trace(T):
    """Compute c_4(T) = tr(A^4)/4 using matrix multiplication. O(n^3)."""
    n = len(T)
    if n < 4:
        return 0
    # A^2
    A2 = [[sum(T[i][k]*T[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    # A^4 = A^2 * A^2
    A4 = [[sum(A2[i][k]*A2[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    tr4 = sum(A4[i][i] for i in range(n))
    assert tr4 % 4 == 0, f"tr(A^4) = {tr4} not divisible by 4!"
    return tr4 // 4


# ============================================================
# Exhaustive verification
# ============================================================
print("=" * 70)
print("EXTENDED THM-096: tr(A^k) = k * c_k for k = 3, 4, 5")
print("=" * 70)

for n in range(3, 8):
    m = n * (n - 1) // 2
    total = 1 << m
    if n >= 7:
        random.seed(42)
        sample = [random.randint(0, total - 1) for _ in range(2000)]
        sample_type = "sampled (2000)"
    else:
        sample = range(total)
        sample_type = "exhaustive"

    k3_match = 0
    k4_match = 0
    k5_match = 0
    k6_match = 0
    tested = 0

    for bits in sample:
        T = tournament_from_bits(n, bits)

        # k=3
        tr3 = matrix_power_trace(T, 3)
        c3 = count_directed_k_cycles_dp(T, 3)
        if tr3 == 3 * c3:
            k3_match += 1

        # k=4
        if n >= 4:
            tr4 = matrix_power_trace(T, 4)
            c4 = count_directed_k_cycles_dp(T, 4)
            if tr4 == 4 * c4:
                k4_match += 1

        # k=5
        if n >= 5:
            tr5 = matrix_power_trace(T, 5)
            c5 = count_directed_k_cycles_dp(T, 5)
            if tr5 == 5 * c5:
                k5_match += 1

        # k=6 (expect FAILURE at n>=5)
        if n >= 5 and tested < 500:
            tr6 = matrix_power_trace(T, 6)
            c6 = count_directed_k_cycles_dp(T, 6) if n >= 6 else 0
            if tr6 == 6 * c6:
                k6_match += 1
            elif tested < 3:
                print(f"  k=6 mismatch at n={n}: bits={bits}, tr6={tr6}, 6*c6={6*c6}, "
                      f"diff={tr6-6*c6}")

        tested += 1

    print(f"\nn={n} ({sample_type}, {tested} tournaments):")
    print(f"  k=3: tr(A^3) = 3*c_3: {k3_match}/{tested} "
          f"({'PASS' if k3_match == tested else 'FAIL'})")
    if n >= 4:
        print(f"  k=4: tr(A^4) = 4*c_4: {k4_match}/{tested} "
              f"({'PASS' if k4_match == tested else 'FAIL'})")
    if n >= 5:
        print(f"  k=5: tr(A^5) = 5*c_5: {k5_match}/{tested} "
              f"({'PASS' if k5_match == tested else 'FAIL'})")
    if n >= 5:
        k6_tested = min(500, tested)
        print(f"  k=6: tr(A^6) = 6*c_6: {k6_match}/{k6_tested} "
              f"({'PASS' if k6_match == k6_tested else 'FAIL -- compound walks exist'})")


# ============================================================
# c_4 via trace: verification and timing
# ============================================================
print("\n" + "=" * 70)
print("c_4 via trace: c_4(T) = tr(A^4)/4 in O(n^3)")
print("=" * 70)

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    all_match = True
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        c4_dp = count_directed_k_cycles_dp(T, 4)
        c4_tr = c4_via_trace(T)
        if c4_dp != c4_tr:
            print(f"  n={n}, bits={bits}: MISMATCH dp={c4_dp}, trace={c4_tr}")
            all_match = False
    if all_match:
        print(f"  n={n}: c4_via_trace matches for ALL {1 << m} tournaments")


# ============================================================
# Timing comparison
# ============================================================
import time
print("\n" + "=" * 70)
print("Timing: DP vs trace for c_4 and c_5 at n=8")
print("=" * 70)

n = 8
m = n * (n - 1) // 2
random.seed(42)
sample_bits = [random.randint(0, (1 << m) - 1) for _ in range(10)]

# c4 timing
t0 = time.time()
for bits in sample_bits:
    T = tournament_from_bits(n, bits)
    c4_dp = count_directed_k_cycles_dp(T, 4)
t_c4_dp = time.time() - t0

t0 = time.time()
for bits in sample_bits:
    T = tournament_from_bits(n, bits)
    c4_tr = c4_via_trace(T)
t_c4_tr = time.time() - t0

print(f"  c_4 DP: {t_c4_dp:.4f}s, trace: {t_c4_tr:.4f}s, "
      f"speedup: {t_c4_dp/max(t_c4_tr, 1e-6):.1f}x")

# c5 timing
from tournament_fast import c5_fast

t0 = time.time()
for bits in sample_bits:
    T = tournament_from_bits(n, bits)
    c5_dp = count_directed_k_cycles_dp(T, 5)
t_c5_dp = time.time() - t0

t0 = time.time()
for bits in sample_bits:
    T = tournament_from_bits(n, bits)
    c5_tr = c5_fast(T)
t_c5_tr = time.time() - t0

print(f"  c_5 DP: {t_c5_dp:.4f}s, trace: {t_c5_tr:.4f}s, "
      f"speedup: {t_c5_dp/max(t_c5_tr, 1e-6):.1f}x")


# ============================================================
# BONUS: Can we get c_6 via correction?
# ============================================================
print("\n" + "=" * 70)
print("c_6 correction from (3,3)-compound walks")
print("=" * 70)
print("At k=6, non-simple walks are (3,3)-compound: two directed 3-cycles")
print("sharing a vertex. Can we compute the correction from c_3 data?")

for n in [5, 6]:
    m = n * (n - 1) // 2
    total_tours = 1 << m

    # For each tournament, compute the correction and test formula
    matches_formula = 0
    tested = 0

    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)
        nn = len(T)
        tr6 = matrix_power_trace(T, 6)
        c6 = count_directed_k_cycles_dp(T, 6) if nn >= 6 else 0
        excess = tr6 - 6 * c6

        # t3(v) = # directed 3-cycles through v
        t3v = [0] * nn
        for i in range(nn):
            for j in range(nn):
                if i == j or not T[i][j]:
                    continue
                for k in range(nn):
                    if k == i or k == j or not T[j][k] or not T[k][i]:
                        continue
                    t3v[i] += 1

        # Compound (3,3) walk at v: choose two directed 3-cycles through v.
        # Walk: v -> a -> b -> v -> c -> d -> v (length 6).
        # The two 3-cycles share vertex v. They CAN share other vertices.
        # Number of such walks starting at v with "first-3-cycle" arrangement:
        # t3(v) * t3(v) (all ordered pairs, including same cycle twice!)
        # Wait, can we pick the SAME 3-cycle twice?
        # Walk: v->a->b->v->a->b->v. This requires edges v->a, a->b, b->v, v->a, a->b, b->v.
        # Yes, it's valid! It's a walk (not a path). It contributes to tr(A^6).
        # So compound walks include the same cycle traversed twice.

        sum_t3_sq = sum(t3v[v] ** 2 for v in range(nn))

        # Each compound walk v->(3cyc)->v->(3cyc)->v has the walk rotated by 3
        # starting at the second v, giving v->(3cyc)->v->(3cyc)->v again.
        # So among 6 rotations, position 0 and 3 start at v.
        # So sum_t3_sq counts 2/6 = 1/3 of the total contribution.
        # excess = 3 * sum_t3_sq?

        if excess == 3 * sum_t3_sq:
            matches_formula += 1
        elif tested < 3 and excess != 0:
            print(f"  n={n} bits={bits}: excess={excess}, 3*sum_t3_sq={3*sum_t3_sq}, "
                  f"sum_t3_sq={sum_t3_sq}")

        tested += 1

    print(f"  n={n}: excess = 3*sum_v t3(v)^2: {matches_formula}/{tested} "
          f"({'PASS' if matches_formula == tested else 'FAIL'})")

    if matches_formula != tested:
        # Try other multipliers
        for mult_num, mult_den in [(1,1), (2,1), (6,1), (3,2), (7,2)]:
            m_count = 0
            for bits in range(total_tours):
                T = tournament_from_bits(n, bits)
                nn = len(T)
                tr6 = matrix_power_trace(T, 6)
                c6 = count_directed_k_cycles_dp(T, 6) if nn >= 6 else 0
                excess = tr6 - 6 * c6
                t3v = [0] * nn
                for i in range(nn):
                    for j in range(nn):
                        if i == j or not T[i][j]:
                            continue
                    for k in range(nn):
                        if k == i or k == j:
                            continue
                        if T[j][k] and T[k][i]:
                            t3v[i] += 1
                sum_t3_sq = sum(t3v[v] ** 2 for v in range(nn))
                if excess * mult_den == mult_num * sum_t3_sq:
                    m_count += 1
            if m_count > matches_formula:
                print(f"    mult={mult_num}/{mult_den}: {m_count}/{tested}")

print("\nDone.")
