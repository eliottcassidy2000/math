#!/usr/bin/env python3
"""periodicity_universality.py — The k-periodicity principle.

Session: kind-pasteur-2026-03-20-S8

THE ABSTRACT PRINCIPLE:
When counting (objects + basepoint) modulo symmetry:
  P(n) = n * T(n) - D(n)
The deficit D decomposes by DEPTH of cycle-type partitions.
Level d is exact for n <= (min_cycle)*(d+1) - 1.

So: min_cycle = 3 (tournaments) → 3-periodicity
    min_cycle = 2 (graphs)      → 2-periodicity
    min_cycle = 5 (???)         → 5-periodicity

TEST:
1. Compute P_graph(n) for simple graphs — expect 2-periodicity
2. Compute P_tournament(n) — known 3-periodicity
3. What structure has 5-periodicity? (hint: 5-regular hypergraphs?)
4. What has PRIME periodicity? (structures where Aut has prime-order elements only)
5. The "odd-only" filter: tournaments force odd cycles only.
   What forces "prime cycles only"? Or "p-cycles only"?
"""

from math import factorial, comb
from fractions import Fraction
from collections import defaultdict

# ================================================================
# BURNSIDE MACHINERY
# ================================================================

def partitions(n, max_part=None):
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest

def cycle_type_to_a(partition, n):
    a = [0] * (n + 1)
    for part in partition:
        a[part] += 1
    return a

def num_perms(a, n):
    result = factorial(n)
    for k in range(1, n + 1):
        if a[k] > 0:
            result //= (k ** a[k]) * factorial(a[k])
    return result

def edge_orbits(perm, n):
    """Count orbits on UNORDERED pairs under permutation perm."""
    visited = set()
    orbits = 0
    for i in range(n):
        for j in range(i+1, n):
            if (i,j) not in visited:
                orbits += 1
                ci, cj = i, j
                while (min(ci,cj), max(ci,cj)) not in visited:
                    visited.add((min(ci,cj), max(ci,cj)))
                    ci, cj = perm[ci], perm[cj]
    return orbits

def build_perm(partition, n):
    """Build a concrete permutation with given cycle type."""
    a = cycle_type_to_a(partition, n)
    perm = list(range(n))
    pos = 0
    for k in range(1, n + 1):
        for _ in range(a[k]):
            for i in range(k - 1):
                perm[pos + i] = pos + i + 1
            perm[pos + k - 1] = pos
            pos += k
    return perm

def ordered_pair_free_orbits(perm, n):
    """Count free ordered-pair orbits for tournaments.
    Returns -1 if self-reverse orbit exists (c_T = 0)."""
    visited = [[False]*n for _ in range(n)]
    free = 0
    for i in range(n):
        for j in range(n):
            if i == j or visited[i][j]:
                continue
            orbit = []
            ci, cj = i, j
            while not visited[ci][cj]:
                visited[ci][cj] = True
                orbit.append((ci, cj))
                ci, cj = perm[ci], perm[cj]
            reverse_in = any(a == j and b == i for a, b in orbit)
            if not reverse_in:
                ci, cj = j, i
                while not visited[ci][cj]:
                    visited[ci][cj] = True
                    ci, cj = perm[ci], perm[cj]
                free += 1
            else:
                return -1
    return free


# ================================================================
# GRAPH CASE: min_cycle = 2 → expected 2-periodicity
# ================================================================

def compute_graph_P(max_n):
    """Compute P_graph(n) = rooted simple graph count."""
    print(f"\n  {'='*60}")
    print(f"  GRAPHS: min_cycle = 2 → expect 2-periodicity")
    print(f"  {'='*60}")

    results = {}

    for n in range(1, max_n + 1):
        total_T = 0
        total_P = 0

        for partition in partitions(n):
            a = cycle_type_to_a(partition, n)
            N = num_perms(a, n)
            fix = a[1]
            perm = build_perm(partition, n)
            f = edge_orbits(perm, n)
            c_G = 2 ** f  # graphs: every edge orbit freely present/absent

            total_T += N * c_G
            total_P += N * c_G * fix

        T_n = total_T // factorial(n)
        P_n = total_P // factorial(n)
        D_n = n * T_n - P_n

        results[n] = {'T': T_n, 'P': P_n, 'D': D_n}

    # Known graph counts A000088
    known = {1:1, 2:2, 3:4, 4:11, 5:34, 6:156, 7:1044, 8:12346}

    print(f"\n  {'n':>3} {'T_G(n)':>8} {'P_G(n)':>10} {'D_G(n)':>10} {'D/2':>10}")
    for n in sorted(results.keys()):
        r = results[n]
        D_half = r['D'] // 2 if r['D'] % 2 == 0 else f"{r['D']}/2"
        print(f"  {n:3d} {r['T']:8d} {r['P']:10d} {r['D']:10d} {str(D_half):>10}")

    # Depth decomposition for graphs
    # ALL partitions contribute (not just odd-part ones)
    # Depth = number of parts > 1
    # Level d exact for n <= 2*(d+1) - 1 = 2d+1

    print(f"\n  Depth decomposition:")
    for n in range(2, max_n + 1):
        by_depth = defaultdict(lambda: Fraction(0))
        for partition in partitions(n):
            if partition == (1,)*n:
                continue
            a = cycle_type_to_a(partition, n)
            N = num_perms(a, n)
            fix = a[1]
            perm = build_perm(partition, n)
            f = edge_orbits(perm, n)
            c_G = 2 ** f

            non_one = tuple(p for p in partition if p > 1)
            depth = len(non_one)
            by_depth[depth] += Fraction(N * c_G * (n - fix), factorial(n))

        line = f"  n={n}: "
        for d in sorted(by_depth.keys()):
            line += f"D{d}={float(by_depth[d]):.2f}  "
        print(line)

    return results


# ================================================================
# TOURNAMENT CASE: min_cycle = 3 → known 3-periodicity
# ================================================================

def compute_tournament_P(max_n):
    """Compute P_tournament(n) for comparison."""
    print(f"\n  {'='*60}")
    print(f"  TOURNAMENTS: min_cycle = 3 → known 3-periodicity")
    print(f"  {'='*60}")

    results = {}

    for n in range(1, max_n + 1):
        total_T = 0
        total_P = 0

        for partition in partitions(n):
            a = cycle_type_to_a(partition, n)
            N = num_perms(a, n)
            fix = a[1]
            perm = build_perm(partition, n)
            f = ordered_pair_free_orbits(perm, n)
            if f < 0:
                continue
            c_T = 2 ** f

            total_T += N * c_T
            total_P += N * c_T * fix

        T_n = total_T // factorial(n)
        P_n = total_P // factorial(n)
        D_n = n * T_n - P_n

        results[n] = {'T': T_n, 'P': P_n, 'D': D_n}

    print(f"\n  {'n':>3} {'T(n)':>8} {'P(n)':>10} {'D(n)':>10}")
    for n in sorted(results.keys()):
        r = results[n]
        print(f"  {n:3d} {r['T']:8d} {r['P']:10d} {r['D']:10d}")

    return results


# ================================================================
# THE ABSTRACT k-PERIODICITY PRINCIPLE
# ================================================================

def analyze_periodicity(results, label, min_cycle):
    """Analyze the depth decomposition for k-periodicity."""
    print(f"\n  {'='*60}")
    print(f"  {label}: PERIODICITY ANALYSIS (min_cycle = {min_cycle})")
    print(f"  {'='*60}")

    print(f"\n  Tower predictions (level d exact for n <= {min_cycle}*(d+1)-1 = {min_cycle}d+{min_cycle-1}):")
    for d in range(5):
        exact_to = min_cycle * (d + 1) - 1
        print(f"    Level {d}: exact for n <= {exact_to}")

    # Test: is D(n) = 0 for n < min_cycle?
    for n in sorted(results.keys()):
        D = results[n]['D']
        if n < min_cycle:
            status = "YES (D=0)" if D == 0 else f"NO (D={D})"
            print(f"    n={n} < min_cycle={min_cycle}: D=0? {status}")


# ================================================================
# CREATIVE: WHAT HAS 5-PERIODICITY?
# ================================================================

def explore_5_periodicity():
    """What combinatorial structure has min_cycle = 5?

    Key insight: we need automorphism groups with NO elements of order 2, 3, or 4.
    Only elements of order >= 5 (or 1).

    Possibility 1: Structures where Aut is a p-group for p=5.
    Possibility 2: Structures with a "5-orientation" constraint.
    Possibility 3: Complete 5-uniform hypergraphs with specific orientations.
    """
    print(f"\n  {'='*60}")
    print(f"  CREATIVE: WHAT HAS 5-PERIODICITY?")
    print(f"  {'='*60}")

    print(f"""
  For k-periodicity, we need structures where the minimum order of
  a non-identity automorphism is k. In tournaments, Moon's theorem
  forces all automorphisms to have odd order, making min_order = 3.

  To get 5-periodicity, we need min_order = 5. Possibilities:

  1. PRIME-ORDER TOURNAMENTS (hypothetical):
     A "5-tournament" on vertex set Z_5: complete directed graph
     where the automorphism group is restricted to elements of order
     1 or 5 (i.e., Z_5 or trivial).
     This would require: no automorphism of order 3.
     Paley T_7 has |Aut| = 21 = 3*7, which includes order-3 elements.
     We need tournaments where |Aut| is a power of 5 or trivial.

  2. PENTAGRAPH ORIENTATIONS:
     Consider the complete graph K_n with edges grouped into pentagons
     (5-cycles). Orient each pentagon consistently. The automorphism
     group of such a structure might have min_order = 5.

  3. THE ODD-PRIME FILTER:
     In tournaments, even cycles are killed (c_T = 0).
     What if we also kill 3-cycles? Then min_order = 5.
     This would require: arcs reversed by any 3-cycle automorphism
     must be contradictory.
     A "3-free tournament" has no order-3 automorphism.
     These are tournaments where |Aut| is coprime to 3.
     The deficit for 3-free tournaments would have 5-periodicity!

  4. QUANDLE COLORINGS:
     A quandle is a self-distributive algebraic structure.
     The fundamental quandle of a knot has automorphisms related
     to the knot group. For torus knots T(5,q), the quandle
     automorphisms have order 5, giving 5-periodicity.

  5. DESIGNS / STEINER SYSTEMS:
     A Steiner system S(2, 5, n) has blocks of size 5.
     The automorphism group acts on blocks, and the minimum
     non-trivial cycle on blocks has length related to 5.
     The "rooted Steiner system" count might have 5-periodicity.

  THE GENERAL PRINCIPLE:
  k-periodicity in the tower of approximations occurs when the
  SMALLEST BUILDING BLOCK has size k. In tournaments, the building
  block is the 3-cycle. In graphs, the building block is the edge (2-swap).
  In 5-uniform structures, the building block has size 5.

  This is fundamentally about the CHEEGER CONSTANT:
  - k-periodicity means the expansion bottleneck has scale k
  - Larger k = fewer symmetries = better expansion = sparser correction
  - The "ideal" structure has k = infinity (no symmetries at all)
""")


# ================================================================
# THE PRIME HIERARCHY
# ================================================================

def prime_hierarchy():
    """The hierarchy of periodicities corresponds to primes."""
    print(f"\n  {'='*60}")
    print(f"  THE PRIME HIERARCHY OF PERIODICITIES")
    print(f"  {'='*60}")

    print(f"""
  p=2: GRAPHS           | Level d exact for n <= 2d+1
       Every partition contributes. Densest correction.
       The "easiest" symmetry — transpositions are cheap.

  p=3: TOURNAMENTS       | Level d exact for n <= 3d+2
       Only odd partitions contribute. Sparser correction.
       3-cycles are the "cheapest" tournament symmetry.
       The FILTER: even cycles kill tournaments.

  p=5: ???-STRUCTURES    | Level d exact for n <= 5d+4
       Only 5-smooth partitions contribute. Very sparse.
       The FILTER: 2-cycles AND 3-cycles are killed.

  p=7: ????-STRUCTURES   | Level d exact for n <= 7d+6
       Even sparser. Only 7-smooth odd partitions.

  GENERAL p: Level d exact for n <= pd + (p-1)
       Only partitions into parts divisible by p contribute.

  THE DREAM: p = infinity (no symmetries at all)
       D(n) = 0 for all n. P(n) = n*T(n) exactly.
       This is the "generic" or "rigid" case.

  CONNECTION TO NUMBER THEORY:
  The prime hierarchy is exactly the TOWER OF SIEVING:
  - Sieve by 2: removes even-order symmetries (graph → tournament filter)
  - Sieve by 3: removes 3-order symmetries
  - Sieve by 5: removes 5-order symmetries
  - ...
  Each sieve level removes a class of symmetries and increases
  the periodicity of the approximation tower.

  This is the CHEBYSHEV SIEVE applied to automorphism groups:
  just as Chebyshev sieves primes from integers,
  we sieve cycle types from permutation groups.

  THE EULER PRODUCT ANALOGY:
  The deficit generating function should factor as an Euler product
  over odd primes, each encoding a cycle-length contribution.
  This is a multiplicative structure for the automorphism correction!
""")


# ================================================================
# COMPUTE: WHAT IF WE FILTER TO ODD PRIMES ONLY?
# ================================================================

def filtered_periodicity(max_n):
    """Compute P(n) for "super-tournaments" where only prime-order
    automorphisms are allowed (filter out composite-order elements)."""
    print(f"\n  {'='*60}")
    print(f"  PRIME-FILTERED TOURNAMENT COUNT")
    print(f"  {'='*60}")

    # Standard tournaments: only odd-cycle partitions contribute
    # "Prime-filtered": only partitions into 1s and odd PRIMES contribute
    # This kills (9, ...), (15, ...) etc. but keeps (3, ...), (5, ...), (7, ...)

    def is_prime(k):
        if k < 2:
            return False
        for i in range(2, int(k**0.5)+1):
            if k % i == 0:
                return False
        return True

    def is_prime_partition(partition):
        """All parts are 1 or odd primes."""
        return all(p == 1 or (p % 2 == 1 and is_prime(p)) for p in partition)

    for n in range(1, max_n + 1):
        total_all_odd = 0
        total_prime_only = 0
        total_P_all = 0
        total_P_prime = 0

        for partition in partitions(n):
            a = cycle_type_to_a(partition, n)
            N = num_perms(a, n)
            fix = a[1]
            perm = build_perm(partition, n)

            # Check if all-odd
            all_odd = all(p % 2 == 1 for p in partition)
            if not all_odd:
                continue

            f = ordered_pair_free_orbits(perm, n)
            if f < 0:
                continue
            c_T = 2 ** f

            total_all_odd += N * c_T
            total_P_all += N * c_T * fix

            if is_prime_partition(partition):
                total_prime_only += N * c_T
                total_P_prime += N * c_T * fix

        T_all = total_all_odd // factorial(n)
        P_all = total_P_all // factorial(n)
        D_all = n * T_all - P_all

        # Prime-filtered: what's the difference?
        diff_T = (total_all_odd - total_prime_only) // factorial(n)
        diff_P = (total_P_all - total_P_prime) // factorial(n)

        if diff_T > 0 or diff_P > 0:
            print(f"  n={n}: T_all={T_all}, T difference from prime filter: {diff_T}, "
                  f"P difference: {diff_P}")
        else:
            print(f"  n={n}: prime filter makes no difference (all contributing parts are prime)")

    print(f"""
  KEY FINDING: For n <= 8, ALL contributing odd partitions have only
  prime parts (1, 3, 5, 7). The first composite odd part is 9,
  which first appears at n=9. So the prime filter has NO EFFECT
  for n <= 8 — the prime hierarchy IS the full hierarchy at small n.

  At n=9: partition (9) has a COMPOSITE (non-prime) odd part 9=3^2.
  But wait — 9 IS not prime. So the prime filter would REMOVE the
  (9) contribution, changing D(9) slightly.

  This means: for n <= 8, the tournament tower IS the prime tower.
  The hierarchy 3, 5, 7 of increasing periodicity IS the prime sequence.
""")


# ================================================================
# MAIN
# ================================================================

def main():
    print("=" * 70)
    print("THE UNIVERSALITY OF k-PERIODICITY")
    print("=" * 70)

    # Graph case
    graph_results = compute_graph_P(8)
    analyze_periodicity(graph_results, "GRAPHS", 2)

    # Tournament case
    tourn_results = compute_tournament_P(8)
    analyze_periodicity(tourn_results, "TOURNAMENTS", 3)

    # Creative exploration
    explore_5_periodicity()
    prime_hierarchy()
    filtered_periodicity(10)

    # COMPARISON TABLE
    print(f"\n  {'='*70}")
    print(f"  COMPARISON: GRAPH vs TOURNAMENT TOWERS")
    print(f"  {'='*70}")
    print(f"\n  {'n':>3} {'T_G':>8} {'P_G':>10} {'D_G':>8} | {'T_T':>8} {'P_T':>10} {'D_T':>8}")
    for n in range(1, 9):
        g = graph_results.get(n, {})
        t = tourn_results.get(n, {})
        print(f"  {n:3d} {g.get('T',0):8d} {g.get('P',0):10d} {g.get('D',0):8d} | "
              f"{t.get('T',0):8d} {t.get('P',0):10d} {t.get('D',0):8d}")


if __name__ == "__main__":
    main()
