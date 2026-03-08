#!/usr/bin/env python3
"""
Tournament Fast Computation Module
====================================
kind-pasteur-2026-03-07-S39

Optimized computation routines that exploit proven theorems to avoid
unnecessary brute-force enumeration.

KEY SPEEDUPS:
  1. F_k(T) mod 2 = C(n-1, k) mod 2 for ALL T (THM-094) — O(1) lookup
  2. H(T) via OCF = I(Omega(T), 2) — avoids 2^n DP when cycle count is small
  3. T(n) counting via Davis/Burnside — avoids exhaustive isomorphism enumeration
  4. Mod-3 Taylor zeros: c_j(T) = 0 mod 3 for j < 2*floor((n-1)/2) (THM-086)
  5. Savchenko cycle formulas for doubly-regular tournaments

Usage:
    from tournament_fast import *

    # Instead of computing F(T,x) mod 2 via DP:
    f_mod2 = f_poly_mod2(n)  # [C(n-1,0)%2, C(n-1,1)%2, ..., C(n-1,n-1)%2]

    # Instead of hamiltonian_path_count(T) when T has few cycles:
    h = hamiltonian_paths_ocf(T)

    # Count non-iso tournaments without enumeration:
    t = tournament_count(n)
"""

from math import gcd, factorial, comb
from fractions import Fraction
from itertools import combinations, permutations


# ===================================================================
# 1. THM-094: F_k(T) mod 2 = C(n-1, k) mod 2
# ===================================================================

def f_poly_mod2(n):
    """Return F(T,x) mod 2 for ANY tournament T on n vertices.
    THM-094: F_k(T) = C(n-1, k) mod 2, equivalently F(T,x) = (1+x)^{n-1} mod 2.
    This is completely tournament-independent!
    Returns list [F_0 mod 2, F_1 mod 2, ..., F_{n-1} mod 2]."""
    return [comb(n - 1, k) % 2 for k in range(n)]


def f_k_is_odd(n, k):
    """Check if F_k(T) is odd for ALL tournaments T on n vertices.
    By THM-094 + Lucas' theorem: F_k is odd iff every binary digit of k
    is <= the corresponding digit of n-1.
    Equivalently: (k & (n-1)) == k."""
    return (k & (n - 1)) == k


def f_k_parity_positions(n):
    """Return (odd_positions, even_positions) where F_k is odd/even for all T.
    Uses Lucas' theorem on C(n-1, k) mod 2."""
    odd = [k for k in range(n) if (k & (n - 1)) == k]
    even = [k for k in range(n) if (k & (n - 1)) != k]
    return odd, even


# ===================================================================
# 2. THM-086: Universal Taylor zeros mod 3
# ===================================================================

def taylor_zero_bound_mod3(n):
    """Return the number of universal Taylor zeros mod 3.
    THM-086: c_j(T) = 0 mod 3 for all j < val(n) = 2*floor((n-1)/2).
    Returns val(n)."""
    return 2 * ((n - 1) // 2)


def taylor_zero_bound_mod_p(n, p):
    """Return the (x-1)-adic valuation of A_n(x) mod p.
    This matches the universal Taylor zero count for:
      p=2: all n (val = n-1)
      p=3: all n >= 3
      p>=5: n >= p+2 only
    Returns the Eulerian valuation."""
    A = eulerian_numbers(n)
    for j in range(n):
        cj = sum(comb(k, j) * A[k] for k in range(n))
        if cj % p != 0:
            return j
    return n


# ===================================================================
# 3. Eulerian number utilities
# ===================================================================

_eulerian_cache = {}

def eulerian_numbers(n):
    """Compute Eulerian numbers A(n, 0), A(n, 1), ..., A(n, n-1).
    Uses recurrence with memoization."""
    if n in _eulerian_cache:
        return _eulerian_cache[n]
    if n == 1:
        result = [1]
    else:
        prev = eulerian_numbers(n - 1)
        A = [0] * n
        for k in range(n):
            A[k] = (k + 1) * prev[k] if k < len(prev) else 0
            if k > 0:
                A[k] += (n - k) * prev[k - 1]
        result = A
    _eulerian_cache[n] = result
    return result


# ===================================================================
# 4. OCF-based Hamiltonian path counting
# ===================================================================

def find_3_cycles(T):
    """Find all directed 3-cycles in tournament T. O(n^3)."""
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                # Check both orientations
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i, j, k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i, k, j))
    return cycles


def find_directed_cycles_of_length(T, length):
    """Find all directed cycles of given odd length. O(C(n,length) * length!)."""
    n = len(T)
    cycles = []
    for verts in combinations(range(n), length):
        first = verts[0]
        for perm in permutations(verts[1:]):
            path = (first,) + perm
            valid = True
            for i in range(length):
                if not T[path[i]][path[(i + 1) % length]]:
                    valid = False
                    break
            if valid:
                cycles.append(path)
    return cycles


def count_3_cycles(T):
    """Count number of directed 3-cycles. O(n^3)."""
    n = len(T)
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    count += 1
                elif T[i][k] and T[k][j] and T[j][i]:
                    count += 1
    return count


def hamiltonian_paths_ocf(T, max_cycle_length=None):
    """Compute H(T) via OCF: H(T) = I(Omega(T), 2).

    For tournaments with few odd cycles, this is MUCH faster than bitmask DP.

    Strategy:
      1. Find all odd directed cycles (3-cycles first, then 5, 7, ...)
      2. Build conflict graph (vertex-disjoint check)
      3. Compute I(Omega, 2) = sum over vertex-disjoint cycle collections of 2^|collection|

    For n <= 11, max independent set size <= 3, so I evaluation is O(m^3).

    Args:
        T: tournament adjacency matrix
        max_cycle_length: if set, only consider cycles up to this length.
            Default None = all odd lengths up to n.
    """
    n = len(T)
    if n <= 1:
        return 1

    # Find all odd cycles
    cycles = []
    for length in range(3, (max_cycle_length or n) + 1, 2):
        if length > n:
            break
        cycles.extend(find_directed_cycles_of_length(T, length))

    if not cycles:
        return 1

    # Use the fast independence polynomial evaluation
    m = len(cycles)
    vsets = [frozenset(c) for c in cycles]

    # Build adjacency bitmasks for conflict graph
    adj_bits = [0] * m
    for a in range(m):
        for b in range(a + 1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    # Max independent set size = max vertex-disjoint odd cycles <= floor(n/3)
    max_indep = n // 3

    # Evaluate I(Omega, 2) by counting independent sets up to max_indep size
    total = 1  # empty set

    # Size 1: each cycle contributes 2
    total += 2 * m

    if max_indep >= 2:
        # Size 2: pairs of vertex-disjoint cycles
        pairs = []
        for a in range(m):
            for b in range(a + 1, m):
                if not (adj_bits[a] & (1 << b)):
                    pairs.append((a, b))
        total += 4 * len(pairs)

        if max_indep >= 3:
            # Size 3: triples of vertex-disjoint cycles
            triples = []
            for a, b in pairs:
                for c in range(b + 1, m):
                    if not (adj_bits[a] & (1 << c)) and not (adj_bits[b] & (1 << c)):
                        triples.append((a, b, c))
            total += 8 * len(triples)

            if max_indep >= 4:
                # Size 4
                quads = 0
                for a, b, c in triples:
                    for d in range(c + 1, m):
                        if (not (adj_bits[a] & (1 << d)) and
                            not (adj_bits[b] & (1 << d)) and
                            not (adj_bits[c] & (1 << d))):
                            quads += 1
                total += 16 * quads

                if max_indep >= 5:
                    # Size 5+: fall back to brute force enumeration
                    # This only happens for n >= 15
                    for size in range(5, max_indep + 1):
                        count = _count_indep_sets_of_size(adj_bits, m, size)
                        total += (2 ** size) * count

    return total


def _count_indep_sets_of_size(adj_bits, m, size):
    """Count independent sets of given size via recursive enumeration."""
    if size == 0:
        return 1
    if size == 1:
        return m
    count = 0

    def recurse(start, chosen_mask, remaining):
        nonlocal count
        if remaining == 0:
            count += 1
            return
        for v in range(start, m):
            if not (adj_bits[v] & chosen_mask):
                recurse(v + 1, chosen_mask | (1 << v), remaining - 1)

    recurse(0, 0, size)
    return count


# ===================================================================
# 4b. Fast cycle finding via bitmask DP
# ===================================================================

def find_odd_cycles_dp(T):
    """Find all directed odd cycles using bitmask DP.

    Uses dp[(mask, v)] = # directed paths from min(mask) through exactly
    the vertices in mask, ending at v. A cycle of length k exists when
    there's an arc from v back to min(mask) and popcount(mask) = k.

    Complexity: O(2^n * n) — much faster than permutation enumeration
    for cycle lengths >= 7.

    Returns: list of (vertex_tuple, count) where count is the number of
    distinct directed cycles on that vertex set. Each directed cycle on
    the same vertex set is a separate node in Omega(T).
    """
    n = len(T)
    results = []  # list of vertex tuples (one per directed cycle)

    for start in range(n):
        # dp[mask][v] = # paths from start through mask ending at v
        # Only consider masks containing start, where start is minimum
        dp = [[0] * n for _ in range(1 << n)]
        dp[1 << start][start] = 1

        for mask in range(1 << start, 1 << n):
            if not (mask & (1 << start)):
                continue
            # Ensure start is minimum in mask
            if mask & ((1 << start) - 1):
                continue
            popcount = bin(mask).count('1')
            for last in range(n):
                if not (mask & (1 << last)):
                    continue
                cnt = dp[mask][last]
                if cnt == 0:
                    continue

                # Try to close cycle (only for odd lengths >= 3)
                if popcount >= 3 and popcount % 2 == 1 and last != start:
                    if T[last][start]:
                        # Found cnt directed cycles on this vertex set
                        verts = tuple(i for i in range(n) if mask & (1 << i))
                        for _ in range(cnt):
                            results.append(verts)

                # Extend path
                if popcount < n:
                    for nxt in range(start + 1, n):
                        if mask & (1 << nxt):
                            continue
                        if T[last][nxt]:
                            dp[mask | (1 << nxt)][nxt] += cnt

        # Clear to save memory
        dp = None

    return results


def hamiltonian_paths_ocf_fast(T):
    """Compute H(T) via OCF using bitmask DP cycle finding.

    Finds all odd directed cycles via bitmask DP (O(2^n * n)),
    builds conflict graph, evaluates I(Omega, 2).

    NOTE: For random tournaments, standard bitmask DP for H(T) is typically
    faster because cycle-finding has the same O(2^n) cost plus independence
    polynomial overhead. OCF is faster for STRUCTURED tournaments with
    few odd cycles (e.g., near-transitive, Paley at small n, circulant).
    """
    n = len(T)
    if n <= 1:
        return 1

    cycles = find_odd_cycles_dp(T)
    if not cycles:
        return 1

    m = len(cycles)
    vsets = [frozenset(c) for c in cycles]

    # Build adjacency bitmasks for conflict graph
    adj_bits = [0] * m
    for a in range(m):
        for b in range(a + 1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    max_indep = n // 3
    total = 1 + 2 * m

    if max_indep >= 2:
        pairs = []
        for a in range(m):
            for b in range(a + 1, m):
                if not (adj_bits[a] & (1 << b)):
                    pairs.append((a, b))
        total += 4 * len(pairs)

        if max_indep >= 3:
            triples = 0
            for a, b in pairs:
                for c in range(b + 1, m):
                    if not (adj_bits[a] & (1 << c)) and not (adj_bits[b] & (1 << c)):
                        triples += 1
            total += 8 * triples

            if max_indep >= 4:
                quads = 0
                for a, b in [(a, b) for a, b in pairs]:
                    for c in range(b + 1, m):
                        if not (adj_bits[a] & (1 << c)) and not (adj_bits[b] & (1 << c)):
                            for d in range(c + 1, m):
                                if (not (adj_bits[a] & (1 << d)) and
                                    not (adj_bits[b] & (1 << d)) and
                                    not (adj_bits[c] & (1 << d))):
                                    quads += 1
                total += 16 * quads

                if max_indep >= 5:
                    for size in range(5, max_indep + 1):
                        count = _count_indep_sets_of_size(adj_bits, m, size)
                        total += (2 ** size) * count

    return total


# ===================================================================
# 5. Davis/Burnside formula for T(n)
# ===================================================================

def partitions_into_odd_parts(n, max_part=None):
    """Generate all partitions of n into odd parts as [(size, mult), ...]."""
    if max_part is None:
        max_part = n
    if max_part % 2 == 0:
        max_part -= 1
    if n == 0:
        yield []
        return
    if max_part <= 0:
        return
    for count in range(n // max_part, 0, -1):
        remainder = n - count * max_part
        for rest in partitions_into_odd_parts(remainder, max_part - 2):
            yield [(max_part, count)] + rest
    yield from partitions_into_odd_parts(n, max_part - 2)


def tournament_count(n):
    """T(n) = number of non-isomorphic tournaments on n vertices.
    OEIS A000568. Uses Davis/Burnside formula: O(p(n/2)) partitions.
    Exact for any n; practical up to n~200."""
    if n <= 1:
        return 1
    nfact = factorial(n)
    total = 0
    for partition in partitions_into_odd_parts(n):
        t = 0
        for size, mult in partition:
            t += mult * (size - 1) // 2
            t += mult * (mult - 1) // 2 * size
        for i in range(len(partition)):
            for j in range(i + 1, len(partition)):
                s1, m1 = partition[i]
                s2, m2 = partition[j]
                t += m1 * m2 * gcd(s1, s2)
        z = 1
        for size, mult in partition:
            z *= (size ** mult) * factorial(mult)
        total += (nfact // z) * (1 << t)
    return total // nfact


# ===================================================================
# 6. Savchenko cycle formulas for regular tournaments
# ===================================================================

def c3_from_score(T):
    """Number of directed 3-cycles from score sequence. O(n^2).
    Formula: c3 = C(n,3) - sum_i C(s_i, 2) where s_i are out-degrees."""
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    return comb(n, 3) - sum(comb(s, 2) for s in scores)


def c3_regular(n):
    """Number of directed 3-cycles in ANY regular tournament on n vertices.
    For regular (odd n), all out-degrees = (n-1)/2, so:
    c3 = C(n,3) - n*C((n-1)/2, 2) = n(n-1)(n+1)/24.
    Derivation: C(n,3) - n*C((n-1)/2,2) = n(n-1)(n-2)/6 - n(n-1)(n-3)/8
    = n(n-1)[4(n-2) - 3(n-3)]/24 = n(n-1)(n+1)/24."""
    if n % 2 == 0:
        raise ValueError("Regular tournaments require odd n")
    return n * (n - 1) * (n + 1) // 24


# ===================================================================
# 7. Forward-edge polynomial (when full F_k values needed)
# ===================================================================

def compute_F_poly(T):
    """Compute full F-polynomial coefficients [F_0, F_1, ..., F_{n-1}].
    Uses bitmask DP tracking forward edges. O(2^n * n^2)."""
    n = len(T)
    if n <= 1:
        return [1]

    # dp[mask][last][fwd] = # paths through mask ending at last with fwd forward edges
    dp = [[[0] * n for _ in range(n)] for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v][0] = 1

    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            for fwd in range(n):
                if dp[mask][last][fwd] == 0:
                    continue
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    new_mask = mask | (1 << nxt)
                    if T[last][nxt]:
                        dp[new_mask][nxt][fwd + 1] += dp[mask][last][fwd]
                    else:
                        dp[new_mask][nxt][fwd] += dp[mask][last][fwd]

    full = (1 << n) - 1
    F = [0] * n
    for last in range(n):
        for fwd in range(n):
            F[fwd] += dp[full][last][fwd]
    return F


# ===================================================================
# Self-test
# ===================================================================

def self_test():
    """Verify all fast routines against known values."""
    # THM-094: F mod 2
    assert f_poly_mod2(3) == [1, 0, 1]  # C(2,k) mod 2
    assert f_poly_mod2(4) == [1, 1, 1, 1]  # C(3,k) mod 2 = all 1
    assert f_poly_mod2(5) == [1, 0, 0, 0, 1]  # C(4,k) mod 2
    assert f_poly_mod2(8) == [1, 1, 1, 1, 1, 1, 1, 1]  # C(7,k) mod 2 = all 1

    # F_k parity positions
    odd5, even5 = f_k_parity_positions(5)
    assert odd5 == [0, 4]
    assert even5 == [1, 2, 3]

    # Taylor zero bound mod 3
    assert taylor_zero_bound_mod3(3) == 2
    assert taylor_zero_bound_mod3(5) == 4
    assert taylor_zero_bound_mod3(7) == 6
    assert taylor_zero_bound_mod3(8) == 6

    # Tournament counting
    assert tournament_count(1) == 1
    assert tournament_count(2) == 1
    assert tournament_count(3) == 2
    assert tournament_count(4) == 4
    assert tournament_count(5) == 12
    assert tournament_count(6) == 56
    assert tournament_count(7) == 456

    # OCF-based H(T) for small cases
    T3c = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]  # cyclic
    T3t = [[0, 1, 1], [0, 0, 1], [0, 0, 0]]  # transitive
    assert hamiltonian_paths_ocf(T3c) == 3
    assert hamiltonian_paths_ocf(T3t) == 1

    # c3 from score sequence
    assert c3_from_score(T3c) == 1
    assert c3_from_score(T3t) == 0

    # c3 for regular tournaments
    assert c3_regular(3) == 1   # 3*2*4/24 = 1
    assert c3_regular(5) == 5   # 5*4*6/24 = 5
    assert c3_regular(7) == 14  # 7*6*8/24 = 14

    # Cross-check OCF with F-polynomial for cyclic T_3
    # H(T) = F_{n-1} (directed ham paths = all-forward permutations)
    F3c = compute_F_poly(T3c)
    assert F3c[-1] == 3  # H(T3c) = 3
    assert sum(F3c) == 6  # sum(F_k) = n! always

    F3t = compute_F_poly(T3t)
    assert F3t[-1] == 1  # H(T3t) = 1
    assert sum(F3t) == 6

    # Cross-check F mod 2 with THM-094
    F3c_mod2 = [f % 2 for f in F3c]
    assert F3c_mod2 == f_poly_mod2(3)

    # Eulerian numbers spot checks
    assert eulerian_numbers(1) == [1]
    assert eulerian_numbers(3) == [1, 4, 1]
    assert eulerian_numbers(4) == [1, 11, 11, 1]

    # Mod-p Taylor zero bound: val_2(n) = n-1 for all n >= 2
    for nn in range(2, 10):
        assert taylor_zero_bound_mod_p(nn, 2) == nn - 1

    # Rotational T_5: i->i+1, i->i+2 mod 5 (unique regular tournament on 5 vertices)
    T5 = [[0]*5 for _ in range(5)]
    for i in range(5):
        for d in [1, 2]:
            T5[i][(i + d) % 5] = 1
    assert hamiltonian_paths_ocf(T5) == 15  # H(rotational T_5) = 15
    assert compute_F_poly(T5)[-1] == 15

    print("All self-tests PASSED.")


if __name__ == "__main__":
    self_test()
