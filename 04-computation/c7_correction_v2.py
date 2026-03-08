import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
c7_correction_v2.py
kind-pasteur-2026-03-07-S39b

CORRECTED analysis of tr(A^7) - 7*c_7.

The correction comes from non-simple closed walks of length 7.
In a tournament, any such walk splits at the first repeated vertex
into sub-walks of length (3, 4), since both sub-walks must be >= 3
(no 2-cycles in tournaments) and 3+4=7.

KEY CORRECTION: The two sub-cycles can share MORE than just the junction
vertex v. At n=5, compound (3,4) walks MUST share >= 2 vertices since
3+4-1=6 > 5. The previous script only counted vertex-disjoint pairs.

New approach: directly enumerate all compound (3,4) walks as sequences
and compare with tr(A^7) - 7*c_7.
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


def count_compound_34_walks_full(T):
    """Count ALL compound (3,4) closed walks of length 7.

    A compound (3,4) walk: start at v, traverse a directed 3-cycle
    returning to v, then traverse a directed 4-cycle returning to v.
    Total length = 3 + 4 = 7.

    The walk is: v -> a -> b -> v -> c -> d -> e -> v
    where v->a->b->v is a directed 3-cycle and v->c->d->e->v is a directed 4-cycle.
    The vertices {a,b} and {c,d,e} CAN overlap (share vertices).

    Also count the reverse arrangement: v -> c -> d -> e -> v -> a -> b -> v.
    But this is just a rotation of the first walk (shifted by 3), so it's
    the SAME closed walk viewed from a different starting point.

    For tr(A^7), what matters is: how many closed walks of length 7 are
    NON-SIMPLE? Each non-simple closed walk contributes 1 to tr(A^7) at
    its starting vertex.

    Approach: directly enumerate all closed walks w_0, w_1, ..., w_7
    (w_7 = w_0) of length 7 that are NOT simple (some vertex repeated).
    """
    n = len(T)
    # Direct enumeration of all closed walks of length 7
    # Walk: w[0] -> w[1] -> ... -> w[6] -> w[0]
    # Total = tr(A^7) by definition.
    # We want to count the NON-SIMPLE ones (vertex repeated).

    total_walks = 0
    simple_walks = 0

    # For small n, enumerate all walks starting at each vertex
    for start in range(n):
        # Use DFS/recursion to enumerate walks of length 7
        count = [0, 0]  # [total, simple]
        _enumerate_walks(T, n, start, [start], 7, count)
        total_walks += count[0]
        simple_walks += count[1]

    return total_walks, simple_walks, total_walks - simple_walks


def _enumerate_walks(T, n, start, path, remaining, count):
    """Recursively enumerate closed walks."""
    if remaining == 0:
        if path[-1] == start:
            count[0] += 1
            # Check simplicity: all vertices in path[0:7] are distinct
            # (path has 8 entries, path[0]=path[7]=start)
            inner = path[:-1]  # length 7
            if len(set(inner)) == len(inner):
                count[1] += 1
        return

    last = path[-1]
    for nxt in range(n):
        if T[last][nxt]:
            path.append(nxt)
            _enumerate_walks(T, n, start, path, remaining - 1, count)
            path.pop()


# ============================================================
# Verification
# ============================================================
print("=" * 70)
print("c_7 CORRECTION v2: Full compound walk enumeration")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    total_tours = 1 << m

    print(f"\nn={n} (exhaustive, {total_tours} tournaments):")

    correction_formula_checks = defaultdict(int)

    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)

        tr7 = matrix_power_trace(T, 7)
        c7 = count_directed_k_cycles_dp(T, 7) if n >= 7 else 0

        # Direct walk count
        total_w, simple_w, nonsimple_w = count_compound_34_walks_full(T)

        excess = tr7 - 7 * c7

        # tr(A^7) = total_w by definition
        # 7 * c_7 = simple_w (each 7-cycle contributes 7 starting points)
        # So excess = nonsimple_w

        check1 = (tr7 == total_w)
        check2 = (7 * c7 == simple_w) if n >= 7 else (simple_w == 0)
        check3 = (excess == nonsimple_w)

        if not check1 or not check2 or not check3:
            print(f"  FAIL bits={bits}: tr7={tr7}, total_w={total_w}, "
                  f"simple_w={simple_w}, nonsimple_w={nonsimple_w}, "
                  f"7*c7={7*c7}, excess={excess}")

        correction_formula_checks[(check1, check2, check3)] += 1

    for key, cnt in sorted(correction_formula_checks.items()):
        print(f"  checks {key}: {cnt}/{total_tours}")

# Now try to find a FORMULA for the non-simple walk count
print("\n" + "=" * 70)
print("FORMULA SEARCH: nonsimple walks from cycle counts")
print("=" * 70)

for n in range(4, 7):
    m = n * (n - 1) // 2
    total_tours = 1 << m

    print(f"\nn={n}:")

    # Collect data
    data = []
    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)
        tr7 = matrix_power_trace(T, 7)
        c3 = count_directed_k_cycles_dp(T, 3)
        c4 = count_directed_k_cycles_dp(T, 4)
        c5 = count_directed_k_cycles_dp(T, 5)
        c7 = count_directed_k_cycles_dp(T, 7) if n >= 7 else 0
        excess = tr7 - 7 * c7

        # Also compute tr(A^3), tr(A^4) for potential formula
        tr3 = matrix_power_trace(T, 3)
        tr4 = matrix_power_trace(T, 4)

        data.append({
            'bits': bits, 'c3': c3, 'c4': c4, 'c5': c5, 'c7': c7,
            'tr3': tr3, 'tr4': tr4, 'tr7': tr7, 'excess': excess
        })

    # Try formulas: excess = f(c3, c4, c5, tr3, tr4, ...)
    # Since non-simple walks are (3,4)-compound, the correction should
    # involve the number of (3-cycle, 4-cycle) pairs that can be composed.

    # For each vertex v, let t3(v) = # directed 3-cycles through v,
    # t4(v) = # directed 4-cycles through v.
    # A compound walk at v: choose a 3-cycle through v AND a 4-cycle through v.
    # Walk: v -> (3-cycle) -> v -> (4-cycle) -> v.
    # The 3-cycle and 4-cycle CAN share vertices (other than v).
    # Each such pair gives ONE walk starting at v (and its reversal is a rotation).
    # So the number of non-simple walks starting at v is:
    # 2 * t3(v) * t4(v)?  No, need to be more careful.

    # Actually: a directed 3-cycle through v has a SPECIFIC traversal from v:
    # if the cycle is v->a->b->v, the walk segment is v->a->b->v.
    # Similarly a directed 4-cycle through v is v->c->d->e->v.
    # Each pair (3-cycle through v, 4-cycle through v) gives:
    # - Walk A: v->a->b->v->c->d->e->v (3-first then 4)
    # But Walk A rotated by 3 gives: v->c->d->e->v->a->b->v (4-first then 3)
    # These are DIFFERENT starting points of the SAME closed walk.
    #
    # So for tr(A^7), each compound walk contributes to 7 diagonal entries.
    # The number of compound walks starting at v (in the "3-first" arrangement)
    # equals t3(v) * t4(v) (for DIRECTED 3-cycles and 4-cycles through v).
    #
    # Wait: t3(v) counts directed 3-cycles through v. There are 2 per vertex set
    # (two orientations), but in tournaments only 1 exists per {a,b,v} triple.
    # So t3(v) = # 3-cycles CONTAINING v (each is a single directed cycle).
    #
    # For directed 4-cycles through v: same, exactly 1 directed cycle per vertex
    # set (well, in a tournament on 4 vertices, there can be 0 or 1 directed 4-cycle
    # on that set, with 2 possible directions... actually up to 3 ham cycles on 4 vertices).
    # Hmm, on 4 vertices there are 3 distinct Hamiltonian cycles, each with 2 directions.
    # In a tournament exactly 1 directed Hamiltonian cycle = 4-cycle exists per 4-vertex set? No.
    # A tournament on 4 vertices can have 0, 1, or 3 directed Hamiltonian cycles.
    #
    # Let me just compute t3(v) * t4(v) and see if the sum predicts the excess.

    for d in data[:3]:
        bits = d['bits']
        T = tournament_from_bits(n, bits)
        nn = len(T)

        # Count directed cycles through each vertex
        # t3(v) = # directed 3-cycles containing v
        t3v = [0] * nn
        for i in range(nn):
            for j in range(nn):
                if i == j or not T[i][j]:
                    continue
                for k in range(nn):
                    if k == i or k == j or not T[j][k] or not T[k][i]:
                        continue
                    # i->j->k->i is a directed 3-cycle through i
                    t3v[i] += 1

        # t4(v) = # directed 4-cycles containing v
        t4v = [0] * nn
        for i in range(nn):
            for j in range(nn):
                if j == i or not T[i][j]:
                    continue
                for k in range(nn):
                    if k == i or k == j or not T[j][k]:
                        continue
                    for l in range(nn):
                        if l == i or l == j or l == k or not T[k][l] or not T[l][i]:
                            continue
                        t4v[i] += 1

        sum_t3_t4 = sum(t3v[v] * t4v[v] for v in range(nn))
        excess = d['excess']

        if bits < 5 or excess != 0:
            print(f"  bits={bits}: excess={excess}, sum_v t3(v)*t4(v)={sum_t3_t4}, "
                  f"c3={d['c3']}, c4={d['c4']}")

    # Now systematic test: does excess = sum_v t3(v)*t4(v)?
    matches = 0
    for d in data:
        bits = d['bits']
        T = tournament_from_bits(n, bits)
        nn = len(T)

        t3v = [0] * nn
        for i in range(nn):
            for j in range(nn):
                if i == j or not T[i][j]:
                    continue
                for k in range(nn):
                    if k == i or k == j or not T[j][k] or not T[k][i]:
                        continue
                    t3v[i] += 1

        t4v = [0] * nn
        for i in range(nn):
            for j in range(nn):
                if j == i or not T[i][j]:
                    continue
                for k in range(nn):
                    if k == i or k == j or not T[j][k]:
                        continue
                    for l in range(nn):
                        if l == i or l == j or l == k or not T[k][l] or not T[l][i]:
                            continue
                        t4v[i] += 1

        sum_t3_t4 = sum(t3v[v] * t4v[v] for v in range(nn))
        excess = d['excess']

        if excess == sum_t3_t4:
            matches += 1

    print(f"  excess = sum_v t3(v)*t4(v): {matches}/{len(data)} "
          f"({'PASS' if matches == len(data) else 'FAIL'})")

    if matches != len(data):
        # Try with coefficient
        ratios = set()
        for d in data:
            if d['excess'] != 0:
                bits = d['bits']
                T = tournament_from_bits(n, bits)
                nn = len(T)
                t3v = [0] * nn
                for i in range(nn):
                    for j in range(nn):
                        if i == j or not T[i][j]:
                            continue
                        for k in range(nn):
                            if k == i or k == j or not T[j][k] or not T[k][i]:
                                continue
                            t3v[i] += 1
                t4v = [0] * nn
                for i in range(nn):
                    for j in range(nn):
                        if j == i or not T[i][j]:
                            continue
                        for k in range(nn):
                            if k == i or k == j or not T[j][k]:
                                continue
                            for l in range(nn):
                                if l == i or l == j or l == k or not T[k][l] or not T[l][i]:
                                    continue
                                t4v[i] += 1
                s = sum(t3v[v] * t4v[v] for v in range(nn))
                if s != 0:
                    ratios.add((d['excess'], s))
        if len(ratios) <= 10:
            print(f"  (excess, sum_t3t4) pairs: {sorted(ratios)[:10]}")


# Test at n=7 (sample)
print("\n" + "=" * 70)
print("n=7 (sampled)")
print("=" * 70)

import random
random.seed(42)
n = 7
m = n * (n - 1) // 2
sample = [random.randint(0, (1 << m) - 1) for _ in range(200)]

matches = 0
for bits in sample:
    T = tournament_from_bits(n, bits)
    nn = len(T)

    tr7 = matrix_power_trace(T, 7)
    c7 = count_directed_k_cycles_dp(T, 7)
    excess = tr7 - 7 * c7

    t3v = [0] * nn
    for i in range(nn):
        for j in range(nn):
            if i == j or not T[i][j]:
                continue
            for k in range(nn):
                if k == i or k == j or not T[j][k] or not T[k][i]:
                    continue
                t3v[i] += 1

    t4v = [0] * nn
    for i in range(nn):
        for j in range(nn):
            if j == i or not T[i][j]:
                continue
            for k in range(nn):
                if k == i or k == j or not T[j][k]:
                    continue
                for l in range(nn):
                    if l == i or l == j or l == k or not T[k][l] or not T[l][i]:
                        continue
                    t4v[i] += 1

    sum_t3_t4 = sum(t3v[v] * t4v[v] for v in range(nn))
    if excess == sum_t3_t4:
        matches += 1
    elif bits in sample[:5]:
        print(f"  bits={bits}: excess={excess}, sum_t3t4={sum_t3_t4}")

print(f"  excess = sum_v t3(v)*t4(v): {matches}/{len(sample)} "
      f"({'PASS' if matches == len(sample) else 'FAIL'})")

print("\nDone.")
