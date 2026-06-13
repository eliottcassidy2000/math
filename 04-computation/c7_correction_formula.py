import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
c7_correction_formula.py
kind-pasteur-2026-03-07-S39b

Extend THM-118 to k=7: find the correction term for compound walks.

tr(A^7) = 7*c_7 + correction

The correction comes from compound walks of length 7:
- (3,4): 3-cycle + 4-cycle sharing a vertex
- (3,4) reverse arrangement
- No (5,2) since min walk length is 3

Question: is the correction computable from c_3, c_4, and local structure?

If correction = f(c_3, c_4, ...) then we get
c_7 = (tr(A^7) - f(c_3, c_4)) / 7
giving O(n^3) computation for c_7 as well.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations, permutations
from collections import defaultdict
from math import comb


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


def count_compound_34_walks(T):
    """Count compound (3,4) closed walks of length 7.

    A compound walk: v -> (traverse 3-cycle through v) -> v -> (traverse 4-cycle through v) -> v
    or: v -> (4-cycle) -> v -> (3-cycle) -> v

    For each vertex v, each directed 3-cycle C3 through v, and each directed
    4-cycle C4 through v with C3 and C4 sharing only v:
    - The walk C3 then C4 is one compound walk starting at v.
    - The walk C4 then C3 is another.
    But we also get compound walks starting at other vertices of the walk.

    Actually, the compound walk v -> a -> b -> v -> c -> d -> e -> v has length 7.
    It's a closed walk, so its rotations start at: v, a, b, v, c, d, e.
    Each rotation contributes 1 to tr(A^7) at the diagonal entry of its starting vertex.

    So the TOTAL contribution to tr(A^7) from this compound walk is 7
    (7 rotations, each contributing 1 to some diagonal entry).

    But how many distinct compound walks are there?
    For vertex v with directed 3-cycles C3_1, ..., C3_a and directed 4-cycles C4_1, ..., C4_b:
    Each pair (C3_i, C4_j) with V(C3_i) intersect V(C4_j) = {v} gives:
    - Walk: C3_i then C4_j (traversing C3 starting/ending at v, then C4 starting/ending at v)
    - Walk: C4_j then C3_i

    Wait, but there's another subtlety: HOW you traverse the 3-cycle from v.
    The 3-cycle a -> b -> c -> a goes through a,b,c. Starting from v=a:
    v -> b -> c -> v. That's the unique traversal (since directed cycle has unique direction).

    For 4-cycle a -> b -> c -> d -> a, starting from v=a:
    v -> b -> c -> d -> v. Unique traversal.

    So: for each vertex v, each pair (C3 through v, C4 through v with C3 cap C4 = {v}):
    - Exactly 2 compound walks (C3 first, or C4 first)
    - Each contributes 7 to tr(A^7)
    - Total per pair: 14

    But compound walks starting at non-v vertices are just rotations.
    For the walk v->a->b->v->c->d->e->v:
    Starting at a: a->b->v->c->d->e->v->a (uses edge v->a, which is the 3-cycle's edge)
    This is ALSO a compound walk, centered at v, rotated.

    So the total compound walk contribution to tr(A^7) is:
    sum over v: sum over pairs (C3, C4) disjoint except at v: 2 * 7 = 14.

    Wait, but if we count "distinct compound walks" (as sequences), each pair
    gives 2 walks (two orderings). Each walk has 7 rotations. But some rotations
    might coincide if there's symmetry. For tournaments, automorphisms are rare,
    so generically 2 * 7 = 14 per pair.

    Actually, the rotations of a compound walk are ALL different compound walks
    (different starting points). They're counted separately in tr(A^7).
    So the contribution per pair IS 2 * 7 = 14.

    Hmm wait, I need to be more careful. Let me just compute it directly.
    """
    n = len(T)
    # Find all directed 3-cycles and 4-cycles through each vertex
    cycles_3 = defaultdict(list)  # v -> list of (a, b) such that v->a->b->v
    cycles_4 = defaultdict(list)  # v -> list of (a, b, c) such that v->a->b->c->v

    # 3-cycles
    for a in range(n):
        for b in range(n):
            if a == b or not T[a][b]:
                continue
            for c in range(n):
                if c == a or c == b or not T[b][c] or not T[c][a]:
                    continue
                # a->b->c->a is a 3-cycle
                cycles_3[a].append((b, c))

    # 4-cycles
    for a in range(n):
        for b in range(n):
            if a == b or not T[a][b]:
                continue
            for c in range(n):
                if c == a or c == b or not T[b][c]:
                    continue
                for d in range(n):
                    if d == a or d == b or d == c or not T[c][d] or not T[d][a]:
                        continue
                    cycles_4[a].append((b, c, d))

    # Count compound walks
    total_compound = 0
    for v in range(n):
        for c3 in cycles_3[v]:
            c3_verts = {v, c3[0], c3[1]}
            for c4 in cycles_4[v]:
                c4_verts = {v, c4[0], c4[1], c4[2]}
                if c3_verts & c4_verts == {v}:
                    # Disjoint except at v: 2 arrangements, 7 contributions each
                    total_compound += 2  # Two walk orderings

    # But wait: each compound walk starting at v has 7 rotations.
    # When I count 2 per pair per v, I'm counting walks with STARTING POINT v.
    # Other rotations start at a, b, c, d, e, or the second v.
    # Since tr(A^7) sums over ALL starting points, I need to count
    # each compound walk ONCE (with all its rotations contributing to tr).
    # So total contribution to tr(A^7) = total_compound * 7.
    # Hmm, but compound walks starting at v vs at a are different walks in the sum.

    # Let me reconsider. tr(A^7) = sum_{all closed walks of length 7}.
    # A compound walk (C3 through v, then C4 through v):
    # v -> a -> b -> v -> c -> d -> e -> v
    # This is ONE closed walk starting at v.
    # Its rotation: a -> b -> v -> c -> d -> e -> v -> a is ANOTHER closed walk
    # starting at a.
    # Total: 7 distinct closed walks from 1 compound walk.
    # But I've already counted the compound walk as "starting at v",
    # "starting at a", etc. through different v choices?

    # No! When I fix v and pair (C3, C4), I get ONE walk starting at v.
    # The rotation starting at a is NOT generated by any v' pairing because
    # the walk a -> b -> v -> c -> d -> e -> v -> a is centered at v (repeating v),
    # not at a (a appears only once).

    # So: for EACH pair (v, C3, C4), there are exactly 2 walks starting at v
    # (C3-first and C4-first). Each such walk contributes 1 to tr(A^7)[v][v].
    # The rotations contribute to tr(A^7)[a][a], [b][b], etc.
    # Total contribution of all rotations of one walk = 7.
    # But I only counted the walk starting at v.
    # The rotations starting at a, b, c, d, e, and the second v are NOT counted
    # by my per-v enumeration because they're not "compound walks starting at those vertices"
    # in my pairing.

    # Hmm, actually the rotation starting at the SECOND v IS a compound walk
    # centered at v: v -> c -> d -> e -> v -> a -> b -> v (C4-first then C3-first).
    # And I DO count this as one of the 2 arrangements!

    # So for the 7 rotations of the walk v->a->b->v->c->d->e->v:
    # - Starting at v (first occurrence): counted as C3-then-C4 arrangement
    # - Starting at a: NOT counted by any pairing
    # - Starting at b: NOT counted
    # - Starting at v (second occurrence): counted as C4-then-C3 arrangement
    # - Starting at c: NOT counted
    # - Starting at d: NOT counted
    # - Starting at e: NOT counted

    # So my counting gives 2 out of 7 rotations. The total is:
    # total_compound_walks * 7 / 2 = total rotations?
    # No. total_compound = number of (v, C3, C4) pairs * 2 arrangements.
    # Each pair generates one closed walk with 7 rotations.
    # Of those 7 rotations, 2 are captured by my counting (the two starting at v).
    # The other 5 start at other vertices and are NOT compound walks from those vertices
    # (they're non-trivial rotations that have the junction vertex in the middle).

    # So total contribution to tr(A^7) = total_compound / 2 * 7 = total_compound * 7 / 2.

    # Wait, let me re-derive. For each pair (v, C3, C4):
    # Walk1: v->a->b->v->c->d->e->v (C3 first)
    # Walk2: v->c->d->e->v->a->b->v (C4 first)
    # Walk1 and Walk2 are two DIFFERENT closed walks.
    # Walk1 has 7 rotations: starts at v(1), a, b, v(2), c, d, e
    # Walk2 has 7 rotations: starts at v(1), c, d, e, v(2), a, b

    # But Walk1 rotated by 3 (starting at v(2)) gives:
    # v->c->d->e->v->a->b->v = Walk2!
    # So Walk1 and Walk2 are ROTATIONS of each other!

    # That means the 2 "arrangements" I counted are actually the SAME walk
    # (just different rotations). So total_compound/2 = number of distinct
    # compound walk orbits under rotation.

    # Total contribution to tr(A^7) = (total_compound/2) * 7

    # Let me verify this with a small example.
    return total_compound  # We'll adjust the formula after verification


# ============================================================
# Verification
# ============================================================
print("=" * 70)
print("c_7 CORRECTION FORMULA")
print("=" * 70)

for n in range(5, 8):
    m = n * (n - 1) // 2
    total_tours = 1 << m
    if n >= 7:
        import random
        random.seed(42)
        sample = [random.randint(0, total_tours - 1) for _ in range(500)]
        sample_type = "sampled"
    else:
        sample = range(total_tours)
        sample_type = "exhaustive"

    print(f"\nn={n} ({sample_type}):")

    correction_check = defaultdict(int)
    for bits in sample:
        T = tournament_from_bits(n, bits)
        if n < 7:
            # Skip if no 7-cycles possible
            pass

        tr7 = matrix_power_trace(T, 7)
        c7 = count_directed_k_cycles_dp(T, 7) if n >= 7 else 0
        c3 = count_directed_k_cycles_dp(T, 3)
        c4 = count_directed_k_cycles_dp(T, 4)

        excess = tr7 - 7 * c7

        # Count compound (3,4) pairs
        compound = count_compound_34_walks(T)

        # Test formula: excess = compound * 7 / 2?
        expected = compound * 7 // 2
        match = (excess == expected)
        correction_check[match] += 1

        if not match and bits < 100:
            print(f"  bits={bits}: excess={excess}, compound*7/2={expected}, "
                  f"compound={compound}")

    total = sum(correction_check.values())
    matches = correction_check.get(True, 0)
    print(f"  excess = compound * 7 / 2: {matches}/{total} "
          f"({'PASS' if matches == total else 'FAIL'})")

    # Try other formulas
    if matches != total:
        print("  Trying other formulas...")
        for bits in (list(sample) if n < 7 else sample[:20]):
            T = tournament_from_bits(n, bits)
            tr7 = matrix_power_trace(T, 7)
            c7 = count_directed_k_cycles_dp(T, 7) if n >= 7 else 0
            compound = count_compound_34_walks(T)
            excess = tr7 - 7 * c7
            if excess != compound * 7 // 2:
                print(f"    bits={bits}: excess={excess}, compound={compound}, "
                      f"7*cpd/2={compound*7//2}, 7*cpd={7*compound}")


# Simpler: test excess = 7 * compound / 2 for n=6
print("\n--- Testing specifically at n=6 ---")
n = 6
m = n * (n - 1) // 2
formula_checks = {}
for formula_name, formula_fn in [
    ("7*cpd/2", lambda e, c: 7*c//2),
    ("7*cpd", lambda e, c: 7*c),
    ("cpd", lambda e, c: c),
    ("cpd/2", lambda e, c: c//2),
]:
    matches = 0
    total = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        tr7 = matrix_power_trace(T, 7)
        c3 = count_directed_k_cycles_dp(T, 3)
        c4 = count_directed_k_cycles_dp(T, 4)
        compound = count_compound_34_walks(T)

        # c_7 at n=6: there are no 7-cycles (need 7 vertices).
        # So excess = tr(A^7) - 0 = tr(A^7).
        excess = tr7
        predicted = formula_fn(excess, compound)

        if excess == predicted:
            matches += 1
        total += 1

    print(f"  {formula_name}: {matches}/{total}")
    formula_checks[formula_name] = matches

# Since n=6 has no 7-cycles, excess = tr(A^7) entirely from compound walks.
# Check: does excess = 7 * (number of distinct compound walk orbits)?
print("\n--- Detailed examples at n=6 ---")
for bits in [0, 1, 2, 3, 4, 5, 10, 20]:
    T = tournament_from_bits(6, bits)
    tr7 = matrix_power_trace(T, 7)
    c3 = count_directed_k_cycles_dp(T, 3)
    c4 = count_directed_k_cycles_dp(T, 4)
    compound = count_compound_34_walks(T)

    print(f"  bits={bits}: tr7={tr7}, c3={c3}, c4={c4}, compound_pairs={compound}")
    if compound > 0:
        print(f"    tr7/7 = {tr7/7:.1f}, compound/2 = {compound/2}")


print("\nDone.")
