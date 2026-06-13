import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
bjorklund_cycle_cover.py
kind-pasteur-2026-03-07-S39b

INVESTIGATION INV-034: Björklund cycle cover reduction adapted for OCF

Björklund (2010) reduces Hamiltonian CYCLE counting to cycle cover counting
via inclusion-exclusion over vertex labelings. Can the same technique,
adapted for Hamiltonian PATHS in tournaments, naturally produce an
odd-cycle formula matching OCF?

Key ideas tested:
1. Cycle cover polynomial: permanent of adjacency matrix of T[S] for subsets S
2. Inclusion-exclusion: H(T) = sum_{S} (-1)^{n-|S|} * cyclecover(T[S])
   (This is for Ham cycles. For paths, we augment with source-sink edges.)
3. Odd vs even cycle covers: do even-cycle covers cancel out for tournaments?
4. Irving-Omar trace formula: log W(z) = 2 * sum_{k odd} tr(A^k) z^k / k

This is exploratory: we want to find if there's a "Björklund-like" formula
that naturally produces I(Omega(T), 2) for tournaments.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from collections import defaultdict
from itertools import permutations


def cycle_cover_count(adj, n):
    """Count directed cycle covers of the digraph adj on n vertices.
    A cycle cover is a collection of directed cycles that covers every vertex.
    This equals the permanent of the adjacency matrix."""
    if n == 0:
        return 1
    # permanent via Ryser's formula
    result = 0
    for mask in range(1 << n):
        # Column sums for subset indicated by mask
        col_sums = [0] * n
        bits_count = 0
        for j in range(n):
            if mask & (1 << j):
                bits_count += 1
                for i in range(n):
                    col_sums[i] += adj[i][j]
        prod = 1
        for i in range(n):
            prod *= col_sums[i]
        if (n - bits_count) % 2 == 0:
            result += prod
        else:
            result -= prod
    return result * ((-1) ** n)


def cycle_cover_by_parity(adj, n):
    """Count cycle covers split by parity of number of even-length cycles.
    Returns (all_odd_covers, has_even_covers) counts."""
    if n == 0:
        return 1, 0
    all_odd = 0
    has_even = 0
    for perm in permutations(range(n)):
        # Check if perm is a valid cycle cover
        valid = True
        for i in range(n):
            if not adj[i][perm[i]]:
                valid = False
                break
        if not valid:
            continue
        # Decompose into cycles, count even-length ones
        visited = [False] * n
        even_count = 0
        for start in range(n):
            if visited[start]:
                continue
            length = 0
            cur = start
            while not visited[cur]:
                visited[cur] = True
                cur = perm[cur]
                length += 1
            if length % 2 == 0:
                even_count += 1
        if even_count == 0:
            all_odd += 1
        else:
            has_even += 1
    return all_odd, has_even


def odd_cycle_cover_weighted(adj, n):
    """For each cycle cover consisting only of odd-length cycles,
    compute 2^{number of cycles} and sum them.
    Compare to I(Omega(T), 2)."""
    if n == 0:
        return 1  # empty cycle cover
    total = 0
    for perm in permutations(range(n)):
        valid = True
        for i in range(n):
            if not adj[i][perm[i]]:
                valid = False
                break
        if not valid:
            continue
        visited = [False] * n
        num_cycles = 0
        all_odd = True
        for start in range(n):
            if visited[start]:
                continue
            length = 0
            cur = start
            while not visited[cur]:
                visited[cur] = True
                cur = perm[cur]
                length += 1
            num_cycles += 1
            if length % 2 == 0:
                all_odd = False
                break
        if all_odd:
            total += 2 ** num_cycles
    return total


def partial_odd_cycle_cover_weighted(adj, n):
    """Sum over all subsets S of vertices: for each S, count all-odd cycle
    covers of T[S] weighted by 2^{num_cycles}.
    This is Sum_S (weighted all-odd cycle covers of T[S])."""
    total = 0
    for mask in range(1 << n):
        verts = [i for i in range(n) if mask & (1 << i)]
        k = len(verts)
        if k == 0:
            total += 1  # empty cover
            continue
        # Build sub-adjacency
        sub = [[0]*k for _ in range(k)]
        for a in range(k):
            for b in range(k):
                sub[a][b] = adj[verts[a]][verts[b]]
        total += odd_cycle_cover_weighted(sub, k)
    return total


def trace_powers(adj, n, max_k):
    """Compute tr(A^k) for k=1,...,max_k."""
    # Matrix power
    A = [row[:] for row in adj]
    traces = []
    Ak = [[int(i == j) for j in range(n)] for i in range(n)]
    for k in range(1, max_k + 1):
        # Ak = Ak * A
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new[i][j] += Ak[i][l] * A[l][j]
        Ak = new
        traces.append(sum(Ak[i][i] for i in range(n)))
    return traces


# ============================================================
# Main
# ============================================================
print("=" * 70)
print("BJÖRKLUND-STYLE CYCLE COVER ANALYSIS FOR TOURNAMENTS")
print("=" * 70)

# Test 1: Do all-odd cycle covers give I(Omega, 2)?
print("\n--- Test 1: All-odd cycle covers weighted by 2^k ---")
print("Compare: full-vertex all-odd CC weighted vs I(Omega, 2) vs H(T)")

for n in range(3, 7):
    m = n * (n - 1) // 2
    limit = min(1 << m, 256)
    matches_full = 0
    matches_partial = 0

    print(f"\nn={n} (testing {limit} tournaments):")

    for bits in range(limit):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)

        # Full vertex all-odd cycle cover weighted
        wocc = odd_cycle_cover_weighted(T, n)
        if wocc == H:
            matches_full += 1

    print(f"  Full-vertex all-odd CC weighted == H: {matches_full}/{limit}")

# Test 2: Partial odd cycle cover weighted (sum over subsets)
print("\n--- Test 2: Partial all-odd cycle cover (sum over subsets) ---")
print("Compare: Sum_S (weighted all-odd CC of T[S]) vs H(T)")

for n in range(3, 6):
    m = n * (n - 1) // 2
    limit = min(1 << m, 64)

    print(f"\nn={n} (testing {limit} tournaments):")
    matches = 0
    examples = []

    for bits in range(limit):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)
        pocc = partial_odd_cycle_cover_weighted(T, n)

        if pocc == H:
            matches += 1
        else:
            if len(examples) < 3:
                examples.append((bits, H, pocc))

    print(f"  Partial all-odd CC == H: {matches}/{limit}")
    for b, h, p in examples:
        print(f"    bits={b}: H={h}, partial_occ={p}, diff={h-p}")


# Test 3: Inclusion-exclusion with signed cycle covers
print("\n--- Test 3: Inclusion-exclusion sum ---")
print("Sum_{S subseteq V} (-1)^{n-|S|} * perm(T[S])")

for n in range(3, 7):
    m = n * (n - 1) // 2
    limit = min(1 << m, 128)

    print(f"\nn={n} (testing {limit} tournaments):")
    ie_matches = 0

    for bits in range(limit):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)

        # Inclusion-exclusion sum
        ie_sum = 0
        for mask in range(1 << n):
            verts = [i for i in range(n) if mask & (1 << i)]
            k = len(verts)
            sub = [[0]*k for _ in range(k)]
            for a in range(k):
                for b in range(k):
                    sub[a][b] = T[verts[a]][verts[b]]
            cc = cycle_cover_count(sub, k)
            sign = (-1) ** (n - k)
            ie_sum += sign * cc

        if ie_sum == H:
            ie_matches += 1
        elif bits < 5:
            print(f"    bits={bits}: H={H}, IE_sum={ie_sum}")

    print(f"  IE sum == H: {ie_matches}/{limit}")


# Test 4: Trace formula (Irving-Omar)
print("\n--- Test 4: Irving-Omar odd trace contributions ---")
print("log W(z) = 2 * sum_{k odd} tr(A^k) z^k / k")
print("Check: sum_{k=1,3,5,...} 2*tr(A^k)/k relates to OCF terms")

for n in range(3, 7):
    m = n * (n - 1) // 2
    limit = min(1 << m, 64)

    print(f"\nn={n}:")
    for bits in range(min(limit, 4)):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)
        traces = trace_powers(T, n, n)

        # tr(A^k) for odd k counts closed walks of length k
        # These include non-simple walks, but the simple ones are directed cycles
        odd_traces = {k+1: traces[k] for k in range(0, n, 2)}
        print(f"  bits={bits}: H={H}, odd_traces={odd_traces}")


# Test 5: Permanent-based formulas
print("\n--- Test 5: perm(I + x*A) at x=1 vs H(T) ---")
print("Does perm(I + A) relate to H(T)?")

for n in range(3, 7):
    m = n * (n - 1) // 2
    limit = min(1 << m, 128)

    print(f"\nn={n} (testing {limit} tournaments):")
    for x_val in [1, 2]:
        matches = 0
        for bits in range(limit):
            T = tournament_from_bits(n, bits)
            H = hamiltonian_path_count(T)

            # Compute perm(I + x*A)
            M = [[int(i == j) + x_val * T[i][j] for j in range(n)] for i in range(n)]
            p = cycle_cover_count(M, n)

            if p == H:
                matches += 1
        print(f"  x={x_val}: perm(I + x*A) == H: {matches}/{limit}")


# Test 6: "Odd permanent" — permanent restricting to odd-cycle permutations
print("\n--- Test 6: Odd permanent (only odd-cycle permutations) ---")

for n in range(3, 7):
    m = n * (n - 1) // 2
    limit = min(1 << m, 128)

    print(f"\nn={n}:")
    matches = 0
    for bits in range(limit):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)

        # Sum over permutations with all odd cycles
        odd_perm = 0
        for perm in permutations(range(n)):
            valid = True
            for i in range(n):
                if not T[i][perm[i]]:
                    valid = False
                    break
            if not valid:
                continue
            # Check all cycles are odd
            visited = [False] * n
            all_odd = True
            for start in range(n):
                if visited[start]:
                    continue
                length = 0
                cur = start
                while not visited[cur]:
                    visited[cur] = True
                    cur = perm[cur]
                    length += 1
                if length % 2 == 0:
                    all_odd = False
                    break
            if all_odd:
                odd_perm += 1

        if odd_perm == H:
            matches += 1
        elif bits < 5:
            print(f"  bits={bits}: H={H}, odd_perm={odd_perm}")

    print(f"  Odd permanent == H: {matches}/{limit}")


print("\nDone.")
