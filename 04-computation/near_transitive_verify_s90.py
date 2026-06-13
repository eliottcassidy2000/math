#!/usr/bin/env python3
"""
near_transitive_verify_s90.py — Verify the near-transitive construction
opus-2026-03-15-S90

CORRECTED: Near-transitive = total order σ with the SINGLE arc
from min(σ) to max(σ) reversed.
"""

from itertools import permutations, combinations
from math import factorial, comb

def count_ham_paths(T, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for p in range(n-1):
            if T[perm[p]][perm[p+1]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

def transitive_core_check(T, n):
    """Check if transitive core is a DAG and total order."""
    core = [[False]*n for _ in range(n)]
    for triple in combinations(range(n), 3):
        for a, b, c in permutations(triple):
            if T[a][b] and T[b][c] and T[a][c]:
                core[a][b] = True
                core[b][c] = True
                core[a][c] = True
                break
    # Kahn DAG check
    in_deg = [sum(1 for i in range(n) if core[i][j]) for j in range(n)]
    queue = [v for v in range(n) if in_deg[v] == 0]
    count = 0
    in_d = in_deg[:]
    while queue:
        v = queue.pop(0)
        count += 1
        for u in range(n):
            if core[v][u]:
                in_d[u] -= 1
                if in_d[u] == 0:
                    queue.append(u)
    if count < n:
        return False, False
    # Total order check
    reach = [row[:] for row in core]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][k] and reach[k][j]:
                    reach[i][j] = True
    total = all(reach[i][j] or reach[j][i] for i in range(n) for j in range(i+1, n))
    return True, total

print("=" * 70)
print("CORRECTED NEAR-TRANSITIVE CONSTRUCTION")
print("=" * 70)
print("""
Construction: Start with total order σ = (v₁ > v₂ > ... > vₙ).
Reverse the SINGLE arc from vₙ (min) to v₁ (max).
Result: vₙ → v₁ (instead of v₁ → vₙ), all other arcs unchanged.

This creates exactly n-2 directed 3-cycles: {v₁, vₖ, vₙ} for k=2,...,n-1.
Each cycle: vₙ → v₁ → vₖ → vₙ (since v₁ beats vₖ and vₖ beats vₙ).
""")

for n in range(3, 9):
    # Canonical: 0 > 1 > 2 > ... > n-1 (so max=0, min=n-1)
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            T[i][j] = 1  # i beats j (i < j means i is "higher" in the order)

    # Reverse arc from min (n-1) to max (0): currently 0→(n-1), reverse to (n-1)→0
    T[0][n-1] = 0
    T[n-1][0] = 1

    c3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        s = T[i][j] + T[j][k] + T[k][i]
        if s == 0 or s == 3:
            c3 += 1

    if n <= 7:
        H = count_ham_paths(T, n)
    else:
        H = "?"

    is_dag, is_total = transitive_core_check(T, n)

    print(f"n={n}: H={H}, c3={c3}, is_DAG={is_dag}, is_total_order={is_total}, "
          f"expected H={2**(n-2)+1}")

# Now verify that this construction matches the exhaustive near-transitive set
print()
print("=" * 70)
print("MATCHING CONSTRUCTION TO EXHAUSTIVE RESULTS")
print("=" * 70)

def all_tournaments(n):
    m = n*(n-1)//2
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    T[i][j] = 1
                else:
                    T[j][i] = 1
                idx += 1
        yield T, bits

for n in [4, 5, 6]:
    print(f"\n--- n={n} ---")

    # Construct near-transitive for each total order (permutation)
    constructed = set()
    for sigma in permutations(range(n)):
        # sigma[0] > sigma[1] > ... > sigma[n-1]
        T = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                # sigma[i] beats sigma[j]
                T[sigma[i]][sigma[j]] = 1

        # Reverse the single arc from min to max
        max_v = sigma[0]
        min_v = sigma[n-1]
        T[max_v][min_v] = 0
        T[min_v][max_v] = 1

        key = tuple(tuple(row) for row in T)
        constructed.add(key)

    # Find exhaustive near-transitive
    exhaustive = set()
    for T, bits in all_tournaments(n):
        is_dag, is_total = transitive_core_check(T, n)
        if is_total:
            H = count_ham_paths(T, n)
            if H > 1:
                key = tuple(tuple(row) for row in T)
                exhaustive.add(key)

    print(f"  Constructed: {len(constructed)}, Exhaustive: {len(exhaustive)}")
    print(f"  MATCH: {constructed == exhaustive}")

# Path structure analysis for the correctly constructed tournament
print()
print("=" * 70)
print("PATH STRUCTURE IN SINGLE-REVERSED-ARC TOURNAMENT")
print("=" * 70)

for n in range(4, 8):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            T[i][j] = 1
    T[0][n-1] = 0
    T[n-1][0] = 1

    paths = []
    for perm in permutations(range(n)):
        ok = True
        for p in range(n-1):
            if T[perm[p]][perm[p+1]] != 1:
                ok = False
                break
        if ok:
            paths.append(perm)

    print(f"\nn={n}: H = {len(paths)} (= 2^{n-2}+1 = {2**(n-2)+1})")
    if n <= 6:
        # Classify paths by whether they use the reversed arc (n-1)→0
        uses_reversed = 0
        for p in paths:
            for k in range(n-1):
                if p[k] == n-1 and p[k+1] == 0:
                    uses_reversed += 1
                    break
        print(f"  Paths using reversed arc (n-1)→0: {uses_reversed}")
        print(f"  Paths NOT using reversed arc: {len(paths) - uses_reversed}")

    # The simplicial path is (0, 1, 2, ..., n-1) using arcs i→(i+1)
    # Check: is this a valid path?
    sim_ok = all(T[i][i+1] == 1 for i in range(n-1))
    print(f"  Path (0,1,...,{n-1}) valid? {sim_ok}")

    # The topo order path is the simplicial one
    # What about the "totally reversed" path?
    rev_ok = all(T[n-1-i][n-2-i] == 1 for i in range(n-1))
    print(f"  Path ({n-1},{n-2},...,0) valid? {rev_ok}")

    # Count by position of the "jump" using arc (n-1)→0
    if n <= 7:
        # Binary encoding: for each path, record which "internal" edges
        # use the reversed arc relationship
        for p in paths:
            arc_types = []
            for k in range(n-1):
                a, b = p[k], p[k+1]
                if a < b:
                    arc_types.append('F')  # forward (i→j with i<j, original)
                elif a > b and not (a == n-1 and b == 0):
                    arc_types.append('B')  # backward (but not the reversed arc)
                else:
                    arc_types.append('R')  # the reversed arc (n-1)→0
            if n <= 5:
                print(f"    {p}  {''.join(arc_types)}")

print()
print("=" * 70)
print("FORMULA DERIVATION")
print("=" * 70)

print("""
WHY H = 2^{n-2} + 1:

In the single-reversed-arc tournament (order 0>1>...>n-1, arc n-1→0):
  Forward arcs: i→j for i<j (except 0→n-1 which is reversed)
  Reversed arc: (n-1)→0

A Hamiltonian path is a sequence using these arcs.

TYPE 1: The path does NOT use arc (n-1)→0.
  Then it only uses forward arcs i→j (i<j, i≠0 or j≠n-1) and
  backward arcs j→i (j>i) that come from the original tournament.
  Wait — the only "backward" arc is (n-1)→0. If we don't use it,
  we can only use forward arcs i→j (i<j). But these form a DAG,
  and the only Hamiltonian path in a DAG with a total order is
  the total order itself: (0, 1, 2, ..., n-1).

  So there's exactly 1 such path (the simplicial/topo-order path).

TYPE 2: The path USES arc (n-1)→0 at some position k.
  Path: (..., n-1, 0, ...)
  Before (n-1): some subset of {1,...,n-2} visited in decreasing order
                (since they form a transitive sub-tournament)
  After 0: the remaining subset of {1,...,n-2} in increasing order

  If the prefix visits vertices S ⊂ {1,...,n-2} (in some valid order),
  then the suffix visits {1,...,n-2}\\S (in some valid order).

  CLAIM: The prefix must end at (n-1) using only forward arcs within S∪{(n-1)},
  and the suffix must start at 0 and visit {1,...,n-2}\\S using forward arcs.

  How many valid prefixes end at (n-1)?
  Any subset S of {1,...,n-2} can serve: the prefix is
  (elements of S in increasing order, then n-1).
  This works because: within S, i→(i') is a forward arc (i<i'),
  and the last element of S → (n-1) is forward (since max(S) < n-1).
  Then (n-1) → 0 (reversed arc).
  Then 0 → remaining elements in increasing order.

  The number of subsets S ⊂ {1,...,n-2} is 2^{n-2}.

  So H = 1 (simplicial) + 2^{n-2} (using reversed arc) - 1 (double count?)
  = 2^{n-2} + 1  ... but wait, we need to check for double-counting.

  Actually, the subset S = ∅ gives prefix = (n-1), suffix = (0, 1, ..., n-2).
  The full path: (n-1, 0, 1, ..., n-2). This is valid!
  And S = {1,...,n-2} gives prefix = (1, 2, ..., n-2, n-1), suffix = (0).
  Full path: (1, 2, ..., n-2, n-1, 0). Valid!

  All 2^{n-2} subsets give distinct valid paths.
  Plus the simplicial path (0, 1, ..., n-1) which doesn't use (n-1)→0.
  Total: 2^{n-2} + 1. ∎
""")

# Verify the subset interpretation
for n in [5, 6, 7]:
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            T[i][j] = 1
    T[0][n-1] = 0
    T[n-1][0] = 1

    # Generate all 2^{n-2} + 1 paths via the subset construction
    constructed_paths = set()

    # Type 1: simplicial path
    constructed_paths.add(tuple(range(n)))

    # Type 2: for each subset S of {1,...,n-2}
    for mask in range(2**(n-2)):
        S = [k+1 for k in range(n-2) if (mask >> k) & 1]
        complement = [k+1 for k in range(n-2) if not ((mask >> k) & 1)]
        path = sorted(S) + [n-1, 0] + sorted(complement)
        constructed_paths.add(tuple(path))

    # Verify all constructed paths are valid
    all_valid = True
    for p in constructed_paths:
        for k in range(n-1):
            if T[p[k]][p[k+1]] != 1:
                all_valid = False
                break

    # Count actual H
    actual_H = 0
    actual_paths = set()
    for perm in permutations(range(n)):
        ok = True
        for p in range(n-1):
            if T[perm[p]][perm[p+1]] != 1:
                ok = False
                break
        if ok:
            actual_H += 1
            actual_paths.add(perm)

    print(f"n={n}: constructed={len(constructed_paths)}, actual H={actual_H}, "
          f"all_valid={all_valid}, match={constructed_paths == actual_paths}")
