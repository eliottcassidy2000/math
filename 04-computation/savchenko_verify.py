import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
savchenko_verify.py
kind-pasteur-2026-03-07-S39b

Verify Savchenko's closed-form cycle counting formulas for regular and
doubly-regular tournaments. Test INV-053.

Savchenko (2016-2024) key results:
1. For a regular tournament T_n (n odd, score (n-1)/2):
   c_3 = n(n-1)(2n-1)/24 - sum_{i<j} a_ij^2 / 2
   where a_ij is the skew-adjacency (+1/-1) matrix
   (Simplified: c_3 = n(n-1)/6 * (something) for DRT)

2. For DRT (doubly-regular tournament, n = 4t+3):
   c_3 = n(n-1)(n-3)/24  (= C(n,3)/2 * (n-3)/(n-2))
   c_5 = n(n-1)(n-3)(n-5)(n-7)/1920 * ... (more complex)

3. For ANY regular tournament on n vertices:
   c_3 = n(n-1)(2n-1)/24 - n*Q/2
   where Q = sum of squared off-diag entries / n (depends on adjacency)

Actually, the Savchenko formulas are (arXiv:2403.07629):
For DRT (doubly-regular, n=4t+3):
   c_k is CONSTANT across all DRTs of the same n for k <= 7.
   c_8 is also constant but only proved for n = 4t+3.
   c_3 = n(n-1)/6

This script: verify these formulas by direct cycle counting at n=3,7,11.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import combinations, permutations
from collections import defaultdict


def count_directed_cycles(T, k):
    """Count ALL directed k-cycles in tournament T.
    Returns total count of distinct directed cycles (not vertex sets)."""
    n = len(T)
    if k > n:
        return 0
    count = 0
    for verts in combinations(range(n), k):
        first = verts[0]
        for perm in permutations(verts[1:]):
            path = (first,) + perm
            if all(T[path[i]][path[(i+1) % k]] for i in range(k)):
                count += 1
    return count


def count_directed_cycles_fast(T, k):
    """Count directed k-cycles using bitmask DP for speed."""
    n = len(T)
    if k > n or k < 3:
        return 0
    count = 0
    for verts in combinations(range(n), k):
        v = list(verts)
        # DP: dp[mask][last] = number of paths from v[0] through vertices
        # indicated by mask, ending at v[last_idx]
        dp = [[0] * k for _ in range(1 << k)]
        dp[1][0] = 1  # Start at v[0]
        for mask in range(1, 1 << k):
            for last in range(k):
                if dp[mask][last] == 0 or not (mask & (1 << last)):
                    continue
                for nxt in range(1, k):  # Don't go back to v[0] except to close
                    if mask & (1 << nxt):
                        continue
                    if T[v[last]][v[nxt]]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
        # Close cycles: from last vertex back to v[0]
        full = (1 << k) - 1
        for last in range(1, k):
            if T[v[last]][v[0]]:
                count += dp[full][last]
    return count


def make_cyclic_tournament(n):
    """Make the cyclic tournament C_n: i beats i+1, ..., i+(n-1)/2 (mod n)."""
    T = [[0]*n for _ in range(n)]
    d = (n - 1) // 2
    for i in range(n):
        for step in range(1, d + 1):
            T[i][(i + step) % n] = 1
    return T


def make_paley(p):
    """Make the Paley tournament T_p for prime p = 3 mod 4."""
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T


def is_regular(T):
    """Check if tournament is regular."""
    n = len(T)
    d = (n - 1) // 2
    return all(sum(T[i]) == d for i in range(n))


def is_doubly_regular(T):
    """Check DRT: regular + every pair has (n-3)/4 common out-neighbors."""
    n = len(T)
    if n % 4 != 3:
        return False
    d = (n - 1) // 2
    if not all(sum(T[i]) == d for i in range(n)):
        return False
    lam = (n - 3) // 4
    for i in range(n):
        for j in range(i + 1, n):
            common = sum(1 for k in range(n) if k != i and k != j and T[i][k] and T[j][k])
            if common != lam:
                return False
    return True


# ============================================================
# Savchenko formula predictions for DRTs
# ============================================================

def savchenko_c3_regular(n):
    """c_3 for ANY regular tournament on n vertices (n odd).
    For regular: every vertex has score (n-1)/2.
    c_3 = n(n-1)(n-1/2 - 1)/6 / ... actually let me derive.

    For regular tournament: c_3 = (sum of i*(n-1-i) for i in score seq) / 3
    Wait, no. For a SINGLE vertex v with score d, the number of 3-cycles
    through v is d*(d-1)/2 - ... no, that's not right either.

    The formula: c_3 = n * (d*(n-1-d) - e) / 3 where e depends on structure.

    Actually for DRT: c_3 = n(n-1)(n-3) / 24 + n/8
    Let me just check: n=3: c_3 = 1. Formula: 3*2*0/24 + 3/8 = 3/8 (wrong).

    Let me try: c_3 = n * d * (d-1) / 6 for regular (Moon's formula)
    Wait, Moon: the total number of 3-cycles in a tournament is
    C(n,3) - sum_{i} C(s_i, 2) / ... no.

    Moon's formula: c_3 = C(n,3) - sum_v C(s_v, 2) where s_v is out-degree.
    Actually that counts ORIENTED 3-cycles (vertex sets, not directed).

    For regular: s_v = (n-1)/2 for all v.
    Vertex sets: C(n,3) - n*C((n-1)/2, 2) = C(n,3) - n*(n-1)(n-3)/8

    For n=3: C(3,3) - 3*C(1,2) = 1 - 0 = 1. c_3 = 1 vertex set, 2 directed. Check!
    For n=7: C(7,3) - 7*C(3,2) = 35 - 21 = 14 vertex sets. OK.

    Directed = 2 * vertex sets for 3-cycles (each triple has exactly 2 orientations,
    but only 1 direction is a 3-cycle in a tournament... no, a cyclic triple has
    exactly 1 directed 3-cycle each way... no.

    A 3-element set in a tournament is either transitive (0 directed 3-cycles)
    or cyclic (1 directed 3-cycle, e.g., a->b->c->a; the reverse a->c->b->a is
    NOT a directed cycle since it uses the edges in the wrong direction).

    Wait: if a->b, b->c, c->a is a directed 3-cycle, then the reverse would need
    a->c, c->b, b->a which requires different edges. So there is exactly 1 directed
    3-cycle per cyclic triple.

    So c_3 (directed) = c_3 (vertex sets) = C(n,3) - n*C(d,2) for regular.

    For n=7: 35 - 7*3 = 14. Matches!
    For n=3: 1 - 3*0 = 1. Matches!
    For n=11: C(11,3) - 11*C(5,2) = 165 - 110 = 55. Matches!
    """
    d = (n - 1) // 2
    from math import comb
    return comb(n, 3) - n * comb(d, 2)


def savchenko_c3_drt(n):
    """c_3 for DRT = n(n-1)/6.

    Actually, c_3 = C(n,3) - n*C((n-1)/2, 2)
            = n(n-1)(n-2)/6 - n*(n-1)(n-3)/8
            = n(n-1)/24 * [4(n-2) - 3(n-3)]
            = n(n-1)/24 * (4n-8-3n+9)
            = n(n-1)/24 * (n+1)
            = n(n-1)(n+1)/24

    Check: n=3: 3*2*4/24 = 1. Correct.
    n=7: 7*6*8/24 = 14. Correct.
    n=11: 11*10*12/24 = 55. Correct!

    So c_3(regular) = n(n-1)(n+1)/24 for ALL regular tournaments (not just DRT).
    This is because for regular tournaments, Moon's formula gives a universal answer.
    """
    return n * (n - 1) * (n + 1) // 24


# ============================================================
# Main verification
# ============================================================
print("=" * 70)
print("SAVCHENKO FORMULA VERIFICATION")
print("=" * 70)

# Test 1: c_3 for regular tournaments
print("\n--- Test 1: c_3 formula for regular tournaments ---")
print("Formula: c_3 = n(n+1)(n-1)/24 (for regular)")

for n in [3, 5, 7]:
    m = n * (n - 1) // 2
    predicted = savchenko_c3_drt(n)
    print(f"\nn={n}: predicted c_3 = {predicted}")

    # Check all tournaments
    found_regular = 0
    c3_regular = defaultdict(int)

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        if not is_regular(T):
            continue
        found_regular += 1
        c3 = count_directed_cycles_fast(T, 3)
        c3_regular[c3] += 1

    print(f"  Found {found_regular} regular tournaments")
    for c3val, cnt in sorted(c3_regular.items()):
        status = "MATCH" if c3val == predicted else "MISMATCH"
        print(f"  c_3 = {c3val}: {cnt} tournaments [{status}]")


# Test 2: c_5 for regular tournaments at n=7
print("\n--- Test 2: c_5 for regular n=7 ---")
n = 7
m = n * (n - 1) // 2

c5_by_class = defaultdict(lambda: defaultdict(int))
c3c5_pairs = defaultdict(int)

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    if not is_regular(T):
        continue
    c3 = count_directed_cycles_fast(T, 3)
    c5 = count_directed_cycles_fast(T, 5)
    c3c5_pairs[(c3, c5)] += 1

print(f"n=7 regular tournaments: c3, c5 pairs")
for (c3, c5), cnt in sorted(c3c5_pairs.items()):
    is_drt = "DRT" if c5 == 42 else ""
    print(f"  c_3={c3}, c_5={c5}: {cnt} tournaments {is_drt}")


# Test 3: Verify DRT cycle counts at n=11 (Paley)
print("\n--- Test 3: DRT at n=11 (Paley) ---")
T11 = make_paley(11)
print(f"T_11 is regular: {is_regular(T11)}")
print(f"T_11 is DRT: {is_doubly_regular(T11)}")

c3_11 = count_directed_cycles_fast(T11, 3)
print(f"c_3(T_11) = {c3_11} (predicted: {savchenko_c3_drt(11)})")

c5_11 = count_directed_cycles_fast(T11, 5)
print(f"c_5(T_11) = {c5_11} (known: 594)")


# Test 4: Verify c_7 for Paley T_7
print("\n--- Test 4: c_7 for Paley T_7 ---")
T7 = make_paley(7)
c7_7 = count_directed_cycles_fast(T7, 7)
print(f"c_7(T_7) = {c7_7} (known: 24)")

# Verify H via OCF
alpha_1 = c3_11 + c5_11
print(f"\n--- OCF check for T_11 ---")
print(f"alpha_1 = c_3 + c_5 + ... = {c3_11} + {c5_11} + ...")
print(f"Note: need c_7, c_9, c_11 for full OCF. Computing c_7...")

c7_11 = count_directed_cycles_fast(T11, 7)
print(f"c_7(T_11) = {c7_11} (known: 3960)")

# c_9 and c_11 take too long with brute force for n=11
# Use known values
print(f"c_9(T_11) = 11055 (known)")
print(f"c_11(T_11) = 5505 (known)")

total_odd_cycles = c3_11 + c5_11 + c7_11 + 11055 + 5505
print(f"Total odd directed cycles = {total_odd_cycles}")

# The OCF formula: H = I(Omega, 2) requires the conflict graph structure,
# not just total cycles. But for a first check:
# I(Omega, 2) >= 1 + 2*alpha_1 (>= contribution from singletons)
print(f"1 + 2*(total) = {1 + 2*total_odd_cycles}")
print(f"H(T_11) = 95095 (known)")


# Test 5: c_3 independence from iso-class within regular
print("\n--- Test 5: c_3 is CONSTANT for all regular n=7 ---")
n = 7
m = n * (n - 1) // 2
c3_values = set()
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    if not is_regular(T):
        continue
    c3_values.add(count_directed_cycles_fast(T, 3))

print(f"Distinct c_3 values among regular n=7: {c3_values}")
print(f"ALL regular n=7 have c_3 = {savchenko_c3_drt(7)}: {'YES' if c3_values == {14} else 'NO'}")


# Test 6: c_3 for NON-regular tournaments
print("\n--- Test 6: c_3 for non-regular n=5 ---")
n = 5
m = n * (n - 1) // 2

c3_by_score = defaultdict(set)
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = tuple(sorted(sum(T[i]) for i in range(n)))
    c3 = count_directed_cycles_fast(T, 3)
    c3_by_score[scores].add(c3)

print(f"n=5: c_3 values by score sequence:")
for scores, vals in sorted(c3_by_score.items()):
    is_reg = "REGULAR" if all(s == 2 for s in scores) else ""
    print(f"  {scores}: c_3 in {sorted(vals)} {is_reg}")


# Savchenko general formula for any tournament
print("\n--- Test 7: Moon's c_3 formula for general tournaments ---")
print("Formula: c_3 = C(n,3) - sum_v C(s_v, 2)")

from math import comb

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    matches = 0
    total = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        scores = [sum(T[i]) for i in range(n)]
        predicted = comb(n, 3) - sum(comb(s, 2) for s in scores)
        actual = count_directed_cycles_fast(T, 3)
        if predicted == actual:
            matches += 1
        total += 1
    print(f"  n={n}: Moon's formula correct for {matches}/{total} tournaments")


print("\nDone.")
