import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
alpha2_formula.py
kind-pasteur-2026-03-07-S39b

NEW RESULT: alpha_2 (disjoint odd cycle pairs) is computable in O(n^3)
via matrix trace data, at least for n <= 7.

FORMULA: alpha_2 = C(c_3, 2) - sum_v C(t_3(v), 2) + s_2
where:
  c_3 = # directed 3-cycles (via Moon's formula, O(n^2))
  t_3(v) = (A^3)[v][v] = # directed 3-cycles through v (O(n^3))
  s_2 = sum_{edges a->b} C((A^2)[b][a], 2) = # edge-sharing 3-cycle pairs (O(n^3))

PROOF:
  Consider the 3-cycle conflict graph Gamma_3 (edges = pairs sharing >= 1 vertex).
  alpha_2 = C(c_3, 2) - |E(Gamma_3)|.

  |E(Gamma_3)| = (pairs sharing >= 1 vertex)
               = (pairs sharing exactly 1 vertex) + (pairs sharing exactly 2 vertices)
               = s_1 + s_2

  s_1 + 2*s_2 = sum_v C(t_3(v), 2)  [each 1-overlap pair counted once at shared vertex;
                                       each 2-overlap pair counted twice at both vertices]

  So s_1 = sum_v C(t_3(v), 2) - 2*s_2, and
  |E| = s_1 + s_2 = sum_v C(t_3(v), 2) - s_2

  Therefore: alpha_2 = C(c_3, 2) - sum_v C(t_3(v), 2) + s_2. QED.

For s_2: two 3-cycles sharing exactly 2 vertices share a directed edge a->b.
The number of 3-cycles containing directed edge a->b is (A^2)[b][a] (counting
the paths b->?->a). So s_2 = sum_{a->b} C((A^2)[b][a], 2).

TOTAL COMPLEXITY: O(n^3) (dominated by A^2 and A^3 computation).

AT n <= 7: alpha_2(full Omega) = alpha_2(3-cycles only) since
5+3=8 > 7 means no 5-cycle can be disjoint from a 3-cycle.
So this formula gives the EXACT alpha_2 for OCF at n <= 7.

AT n >= 8: 5-cycle + 3-cycle disjoint pairs can exist (5+3=8).
Need additional terms for the 5-cycle contribution.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from math import comb
import random


def alpha2_from_trace(T):
    """Compute alpha_2 (disjoint 3-cycle pairs) in O(n^3) from matrix data.

    Returns: alpha_2 value.
    For n <= 7, this equals the full alpha_2 of Omega(T).
    """
    n = len(T)
    if n < 6:
        return 0  # Need 6 vertices for two disjoint 3-cycles

    # c_3 via Moon's formula: O(n^2)
    scores = [sum(T[i]) for i in range(n)]
    c3 = comb(n, 3) - sum(comb(s, 2) for s in scores)

    if c3 < 2:
        return 0

    # A^2 matrix: O(n^3)
    A2 = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = 0
            for k in range(n):
                s += T[i][k] * T[k][j]
            A2[i][j] = s

    # A^3 diagonal (t_3(v) = (A^3)[v][v]): O(n^3)
    t3v = [0] * n
    for v in range(n):
        for j in range(n):
            t3v[v] += A2[v][j] * T[j][v]

    # s_2 = sum of C(e(a,b), 2) over directed edges a->b: O(n^2)
    s2 = 0
    for a in range(n):
        for b in range(n):
            if a != b and T[a][b]:
                eab = A2[b][a]
                s2 += eab * (eab - 1) // 2

    # sum_v C(t_3(v), 2)
    sum_ct3 = sum(t * (t - 1) // 2 for t in t3v)

    return comb(c3, 2) - sum_ct3 + s2


def h_from_trace(T):
    """Compute H(T) using O(n^3) trace formulas.

    Valid for n <= 7 (where alpha_3 = 0 and only 3-cycles matter for alpha_2).
    For n >= 8, this gives a LOWER BOUND (missing alpha_3+ and 5-cycle pairs).
    """
    n = len(T)
    if n <= 1:
        return 1

    # alpha_1 = total odd directed cycles
    scores = [sum(T[i]) for i in range(n)]
    c3 = comb(n, 3) - sum(comb(s, 2) for s in scores)

    # c_5 via trace: O(n^3)
    if n >= 5:
        A2 = [[sum(T[i][k]*T[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        A4 = [[sum(A2[i][k]*A2[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        A5 = [[sum(A4[i][k]*T[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        c5 = sum(A5[i][i] for i in range(n)) // 5
    else:
        c5 = 0

    # c_7: only exists at n >= 7, and only one 7-vertex subset.
    # At n=7: need to compute directly (O(n!) worst case for this one subset).
    # This is still fast since there's only C(7,7)=1 subset.
    if n == 7:
        from itertools import combinations
        c7 = 0
        verts = list(range(7))
        dp = [[0]*7 for _ in range(1 << 7)]
        dp[1][0] = 1
        for mask in range(1, 1 << 7):
            for last in range(7):
                if dp[mask][last] == 0 or not (mask & (1 << last)):
                    continue
                for nxt in range(1, 7):
                    if mask & (1 << nxt):
                        continue
                    if T[verts[last]][verts[nxt]]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
        full = (1 << 7) - 1
        for last in range(1, 7):
            if T[verts[last]][verts[0]]:
                c7 += dp[full][last]
    else:
        c7 = 0

    alpha_1 = c3 + c5 + c7
    alpha_2 = alpha2_from_trace(T)

    return 1 + 2 * alpha_1 + 4 * alpha_2


# ============================================================
# Verification
# ============================================================
print("=" * 70)
print("alpha_2 and H(T) via O(n^3) trace formulas")
print("=" * 70)

from tournament_lib import hamiltonian_path_count

for n in [5, 6, 7]:
    m_bits = n * (n - 1) // 2
    total_tours = 1 << m_bits

    if n >= 7:
        random.seed(42)
        sample = [random.randint(0, total_tours - 1) for _ in range(500)]
        sample_type = f"sampled (500)"
    else:
        sample = range(total_tours)
        sample_type = "exhaustive"

    matches = 0
    tested = 0

    for bits in sample:
        T = tournament_from_bits(n, bits)
        h_true = hamiltonian_path_count(T)
        h_trace = h_from_trace(T)

        if h_true == h_trace:
            matches += 1
        elif tested < 3:
            print(f"  MISMATCH at n={n} bits={bits}: true={h_true}, trace={h_trace}")
        tested += 1

    print(f"\nn={n} ({sample_type}, {tested}): "
          f"H(T) via trace: {matches}/{tested} "
          f"({'PASS' if matches == tested else 'FAIL'})")


# ============================================================
# Timing comparison at n=7
# ============================================================
import time
print("\n" + "=" * 70)
print("Timing: DP vs trace-based H(T) at n=7")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
random.seed(42)
timing_sample = [random.randint(0, (1 << m) - 1) for _ in range(20)]

t0 = time.time()
for bits in timing_sample:
    T = tournament_from_bits(n, bits)
    h_dp = hamiltonian_path_count(T)
t_dp = time.time() - t0

t0 = time.time()
for bits in timing_sample:
    T = tournament_from_bits(n, bits)
    h_tr = h_from_trace(T)
t_tr = time.time() - t0

print(f"  DP: {t_dp:.4f}s")
print(f"  Trace: {t_tr:.4f}s")
print(f"  Speedup: {t_dp/max(t_tr, 1e-6):.1f}x")


print("\nDone.")
