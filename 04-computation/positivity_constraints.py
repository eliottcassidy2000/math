#!/usr/bin/env python3
"""
Positivity constraints from THM-062 forward-edge distribution.

THM-062: a_k(T) = A(n,k) + sum_I 2^{parts(I)} * c_k^{(f_I, n-1)} * I(T) >= 0

This gives linear inequalities on tournament invariants (t3, t5, t7, bc, etc.)
that define a FEASIBLE REGION. We check which actual (t3, t5, ...) values
occur and whether they hit the boundary.

opus-2026-03-07
"""

from math import comb, factorial
from itertools import combinations, permutations

###############################################################################
# Part 0: Eulerian numbers and inflated coefficients
###############################################################################

def eulerian(n, k):
    """Eulerian number A(n,k) = number of perms of {1,...,n} with exactly k ascents."""
    if n == 0:
        return 1 if k == 0 else 0
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def inflated_coeff(f, d, k):
    """
    c_k^{(f,d)} = sum_j A(f+1, j) * C(d-f, k-j) * (-1)^{d-f-k+j}

    where A(f+1, j) is the Eulerian number, C is binomial coefficient.
    """
    total = 0
    for j in range(f + 2):  # A(f+1, j) is 0 for j > f
        eul = eulerian(f + 1, j)
        if eul == 0:
            continue
        binom_arg_n = d - f
        binom_arg_k = k - j
        if binom_arg_k < 0 or binom_arg_k > binom_arg_n:
            continue
        sign = (-1) ** (d - f - k + j)
        total += eul * comb(binom_arg_n, binom_arg_k) * sign
    return total

###############################################################################
# Part 1: n=5, d=4
###############################################################################

def part1_n5():
    print("=" * 70)
    print("PART 1: n=5, d=4")
    print("=" * 70)

    n, d = 5, 4

    # Compute Eulerian numbers A(5, k) for k=0,...,4
    print("\nEulerian numbers A(5, k):")
    A5 = [eulerian(5, k) for k in range(5)]
    for k in range(5):
        print(f"  A(5,{k}) = {A5[k]}")
    print(f"  Sum = {sum(A5)} (should be {factorial(5)})")

    # Compute inflated coefficients
    # At n=5: t3 has f=2, t5 has f=0
    print("\nInflated coefficients c_k^{(f,4)}:")
    print(f"  {'k':>3} | {'c_k^(2,4)':>10} | {'c_k^(0,4)':>10}")
    print(f"  {'-'*3}-+-{'-'*10}-+-{'-'*10}")

    c2 = []
    c0 = []
    for k in range(5):
        val2 = inflated_coeff(2, 4, k)
        val0 = inflated_coeff(0, 4, k)
        c2.append(val2)
        c0.append(val0)
        print(f"  {k:>3} | {val2:>10} | {val0:>10}")

    # The formula: a_k(T) = A(5,k) + 2*c_k^{(2,4)}*t3 + 2*c_k^{(0,4)}*t5 >= 0
    print("\nConstraints: a_k >= 0 means:")
    print("  A(5,k) + 2*c_k^(2,4)*t3 + 2*c_k^(0,4)*t5 >= 0")
    print()
    for k in range(5):
        print(f"  k={k}: {A5[k]} + {2*c2[k]}*t3 + {2*c0[k]}*t5 >= 0")

    return A5, c2, c0

###############################################################################
# Part 2-3: Feasible region at n=5
###############################################################################

def part2_feasible_n5(A5, c2, c0):
    print("\n" + "=" * 70)
    print("PART 2-3: Feasible region in (t3, t5) space at n=5")
    print("=" * 70)

    # Each constraint: A5[k] + 2*c2[k]*t3 + 2*c0[k]*t5 >= 0
    # This is a half-plane in (t3, t5) space.

    print("\nHalf-plane constraints:")
    for k in range(5):
        a, b, c_const = 2*c2[k], 2*c0[k], A5[k]
        if a == 0 and b == 0:
            print(f"  k={k}: {c_const} >= 0 (trivially {'true' if c_const >= 0 else 'FALSE'})")
        else:
            print(f"  k={k}: {a}*t3 + {b}*t5 >= {-c_const}")
            # Express as boundary line
            if b != 0:
                print(f"        => t5 {'<=' if b < 0 else '>='} ({-c_const - a}*t3) / {b}")

    # Check at corners of plausible range
    # t3 in [0, C(5,3)] = [0, 10], t5 in [0, C(5,5)] = [0, 1]
    # But actual t3 range is 0..5 (not all C(5,3)=10 subsets can be 3-cycles simultaneously)

    print("\nFeasibility check over grid t3 in [0,10], t5 in [0,3]:")
    feasible_points = []
    for t3 in range(11):
        for t5 in range(4):
            feasible = True
            for k in range(5):
                val = A5[k] + 2*c2[k]*t3 + 2*c0[k]*t5
                if val < 0:
                    feasible = False
                    break
            if feasible:
                feasible_points.append((t3, t5))

    print(f"  Feasible integer points: {len(feasible_points)}")
    for p in feasible_points:
        vals = [A5[k] + 2*c2[k]*p[0] + 2*c0[k]*p[1] for k in range(5)]
        print(f"    (t3={p[0]}, t5={p[1]}): a_k = {vals}, min={min(vals)}")

###############################################################################
# Part 4-5: Actual tournament data at n=5
###############################################################################

def generate_all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency matrices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    tournaments = []
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
        tournaments.append(A)
    return tournaments

def count_directed_cycles(A, n, cycle_len):
    """Count directed cycles of given length in tournament A."""
    count = 0
    for subset in combinations(range(n), cycle_len):
        # Check all cyclic orderings of subset
        sub = list(subset)
        # Number of distinct directed cycles on this vertex set
        for perm in permutations(sub):
            if perm[0] != sub[0]:  # fix first vertex to avoid counting rotations
                continue
            is_cycle = True
            for i in range(cycle_len):
                if A[perm[i]][perm[(i+1) % cycle_len]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    return count

def count_forward_edges_dist(A, n):
    """Count a_k = #{permutations with exactly k forward edges}."""
    a = [0] * n
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]] == 1)
        a[fwd] += 1
    return a

def part45_actual_n5(A5, c2, c0):
    print("\n" + "=" * 70)
    print("PARTS 4-5: Actual tournament data at n=5")
    print("=" * 70)

    n = 5
    tournaments = generate_all_tournaments(n)
    print(f"\nTotal tournaments: {len(tournaments)}")

    # Collect (t3, t5) -> list of a_k distributions
    data = {}
    for A in tournaments:
        t3 = count_directed_cycles(A, n, 3)
        t5 = count_directed_cycles(A, n, 5)
        a_k = count_forward_edges_dist(A, n)
        key = (t3, t5)
        if key not in data:
            data[key] = {'count': 0, 'a_k_set': set()}
        data[key]['count'] += 1
        data[key]['a_k_set'].add(tuple(a_k))

    print(f"\nDistinct (t3, t5) pairs: {len(data)}")
    print(f"\n{'t3':>4} {'t5':>4} | {'count':>6} | {'a_k values':>40} | {'predicted a_k':>30} | boundary?")
    print("-" * 110)

    boundary_hits = []
    for key in sorted(data.keys()):
        t3, t5 = key
        predicted = [A5[k] + 2*c2[k]*t3 + 2*c0[k]*t5 for k in range(5)]
        for a_k_tuple in sorted(data[key]['a_k_set']):
            a_k_list = list(a_k_tuple)
            match = (a_k_list == predicted)
            has_zero = 0 in a_k_list
            boundary = "BOUNDARY" if has_zero else ""
            min_val = min(a_k_list)
            print(f"  {t3:>3} {t5:>3} | {data[key]['count']:>6} | a_k={a_k_list} | pred={predicted} | match={match} | min={min_val} {boundary}")
            if has_zero:
                boundary_hits.append((t3, t5, a_k_list))

    if boundary_hits:
        print(f"\n  BOUNDARY HITS (a_k = 0 for some k):")
        for t3, t5, a_k in boundary_hits:
            zeros = [k for k in range(len(a_k)) if a_k[k] == 0]
            print(f"    (t3={t3}, t5={t5}): a_k={a_k}, zero at k={zeros}")
    else:
        print(f"\n  No boundary hits: all a_k > 0 for all tournaments.")

###############################################################################
# Part 6: n=7, d=6
###############################################################################

def part6_n7():
    print("\n" + "=" * 70)
    print("PART 6: n=7, d=6")
    print("=" * 70)

    n, d = 7, 6

    # Eulerian numbers A(7, k)
    print("\nEulerian numbers A(7, k):")
    A7 = [eulerian(7, k) for k in range(7)]
    for k in range(7):
        print(f"  A(7,{k}) = {A7[k]}")
    print(f"  Sum = {sum(A7)} (should be {factorial(7)})")

    # Compute inflated coefficients
    # At n=7: t3 has f=4, t5 has f=2, t7 has f=0, bc has f=2 (weight 4)
    print("\nInflated coefficients c_k^{(f,6)}:")
    print(f"  {'k':>3} | {'A(7,k)':>8} | {'c_k^(4,6)':>10} | {'c_k^(2,6)':>10} | {'c_k^(0,6)':>10}")
    print(f"  {'-'*3}-+-{'-'*8}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}")

    c4 = []
    c2 = []
    c0_7 = []
    for k in range(7):
        val4 = inflated_coeff(4, 6, k)
        val2 = inflated_coeff(2, 6, k)
        val0 = inflated_coeff(0, 6, k)
        c4.append(val4)
        c2.append(val2)
        c0_7.append(val0)
        print(f"  {k:>3} | {A7[k]:>8} | {val4:>10} | {val2:>10} | {val0:>10}")

    # Check against provided table (k=0..3)
    expected_c4 = [1, 24, 15, -80]
    expected_c2 = [1, 0, -9, 16]
    expected_c0 = [1, -6, 15, -20]

    print("\nVerification against provided table (k=0..3):")
    for k in range(4):
        ok4 = c4[k] == expected_c4[k]
        ok2 = c2[k] == expected_c2[k]
        ok0 = c0_7[k] == expected_c0[k]
        status = "OK" if (ok4 and ok2 and ok0) else "MISMATCH"
        print(f"  k={k}: c4={c4[k]}({'OK' if ok4 else 'FAIL'}) c2={c2[k]}({'OK' if ok2 else 'FAIL'}) c0={c0_7[k]}({'OK' if ok0 else 'FAIL'})")

    # Constraints: a_k >= 0 means
    # A(7,k) + 2*c4[k]*t3 + 2*c2[k]*t5 + 4*c2[k]*bc + 2*c0[k]*t7 >= 0
    # (bc uses weight 4 = 2^2, and f=2 same as t5)
    print("\nConstraints a_k >= 0:")
    for k in range(7):
        print(f"  k={k}: {A7[k]} + {2*c4[k]}*t3 + {2*c2[k]}*t5 + {4*c2[k]}*bc + {2*c0_7[k]}*t7 >= 0")

    # Which k gives tightest constraint?
    # Test with known tournament data
    # Transitive tournament n=7: t3=C(7,3)=35, t5=0 (no 5-cycles), t7=0, bc=0...
    # Actually transitive has all edges i->j for i<j, so ALL 3-subsets are transitive triples,
    # meaning t3=0 for the TRANSITIVE tournament (no directed 3-cycles)

    print("\n--- Testing with specific tournaments ---")

    # Transitive tournament: t3=0, t5=0, t7=0, bc=0
    print("\nTransitive (t3=0, t5=0, t7=0, bc=0):")
    for k in range(7):
        val = A7[k]
        print(f"  a_{k} = {val}")

    # Regular tournament (Paley T_7): t3=7, t5=21, t7=?, bc=?
    # From memory: all regular n=7 have c3=14 (vertex sets), so t3=14 directed 3-cycles
    # Wait, t3 counts DIRECTED cycles. Each 3-vertex-set with a directed 3-cycle has exactly 2
    # directed cycles (clockwise and counter). Actually no - a 3-cycle vertex set has exactly
    # 2 directed Hamiltonian cycles on those 3 vertices.
    # Actually for directed 3-cycles: each set of 3 vertices forming a 3-cycle contributes
    # exactly 2 directed cycles (two orientations).
    # Hmm, but in a tournament on 3 vertices, either it's transitive (1 Ham path, 0 directed cycles)
    # or cyclic (2 directed Ham cycles = directed 3-cycles).
    # From memory data: c3=14 for Paley T_7 means 14 cyclic triples, so t3 = 14*2 = 28 directed 3-cycles.
    # But wait - the problem statement says t3 is an invariant. Let me just use c3 (number of 3-element
    # subsets forming a 3-cycle).

    # Actually, re-reading the problem: "invariants at n=5 are t3 (f=2) and t5 (f=0)"
    # At n=5: a 3-cycle uses f=2 forward edges (out of d=4), and a 5-cycle uses f=0.
    # So t3 = number of 3-cycles, t5 = number of 5-cycles.
    # For the Paley T_7 (BIBD): From memory data, c3=14 vertex sets => t3 = 2*14 = 28 directed 3-cycles?
    # Actually I think in this context, "number of 3-cycles" might mean the number of
    # 3-element subsets that form a cycle, not directed cycles. Let me be careful.

    # The independence polynomial formulation: I(Omega(T), 2) = sum over independent sets of 2^|S|
    # In OCF: H(T) = sum_I 2^{parts(I)} * I(T) where I ranges over cycle compositions
    # For the a_k formula, the coefficients come from the cycle types.
    # Let me just compute the invariants computationally for specific n=7 tournaments.

    # For now, test with Paley T_7 known values
    # H(Paley_7) = 189, from known data
    # alpha_1 = 80 (total directed 3-cycle vertex sets contributing... hmm)

    # Let me just report the constraint structure and tightest k
    print("\nAnalysis of which k gives tightest constraint:")
    print("The constraint at k=3 has the largest NEGATIVE coefficients:")
    print(f"  k=3: {A7[3]} + {2*c4[3]}*t3 + {2*c2[3]}*t5 + {4*c2[3]}*bc + {2*c0_7[3]}*t7 >= 0")
    print(f"  i.e., {A7[3]} - {-2*c4[3]}*t3 + {2*c2[3]}*t5 + {4*c2[3]}*bc - {-2*c0_7[3]}*t7 >= 0")
    print(f"  This gives: {-2*c4[3]}*t3 + {-2*c0_7[3]}*t7 <= {A7[3]} + {2*c2[3]}*t5 + {4*c2[3]}*bc")

    # By symmetry A(n,k) = A(n,n-1-k), so A(7,3) = A(7,3) is the largest Eulerian number
    # The tightest constraints come from k=0 and k=6 (smallest Eulerian numbers)
    print("\nNote: A(7,0) = A(7,6) = 1 (smallest Eulerian numbers)")
    print("  k=0: 1 + 2*t3 + 2*t5 + 4*bc + 2*t7 >= 0  (always true for non-negative invariants)")
    print(f"  k=6: {A7[6]} + {2*c4[6]}*t3 + {2*c2[6]}*t5 + {4*c2[6]}*bc + {2*c0_7[6]}*t7 >= 0")

    return A7, c4, c2, c0_7

###############################################################################
# Part 7: Does a_k = 0 ever occur?
###############################################################################

def part7_zero_check():
    print("\n" + "=" * 70)
    print("PART 7: Does a_k = 0 ever occur?")
    print("=" * 70)

    # Check n=3, 4, 5 exhaustively
    for n in [3, 4, 5]:
        print(f"\n--- n={n} ---")
        tournaments = generate_all_tournaments(n)
        zero_cases = []
        min_a = {k: float('inf') for k in range(n)}

        for A in tournaments:
            a_k = count_forward_edges_dist(A, n)
            for k in range(n):
                min_a[k] = min(min_a[k], a_k[k])
                if a_k[k] == 0:
                    zero_cases.append((A, k, a_k))

        print(f"  Min a_k across all tournaments: {[min_a[k] for k in range(n)]}")
        if zero_cases:
            print(f"  ZERO cases found: {len(zero_cases)}")
            # Show first few
            shown = set()
            for A, k, a_k in zero_cases[:5]:
                key = tuple(a_k)
                if key not in shown:
                    print(f"    a_k = {a_k}, zero at k={k}")
                    shown.add(key)
        else:
            print(f"  No zero cases: all a_k > 0 for all tournaments")

    # For transitive tournament, we know a_0 = a_{n-1} = 1 (the identity perm and reverse)
    print("\n--- Transitive tournament analysis ---")
    for n in [3, 4, 5, 6, 7]:
        # Build transitive tournament: edge i->j iff i < j
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1
        if n <= 6:
            a_k = count_forward_edges_dist(A, n)
            print(f"  n={n}: a_k = {a_k}")
        else:
            print(f"  n={n}: (skipping full enumeration, too slow)")
            # For transitive tournament, a_k = A(n,k) (Eulerian numbers)
            # because all edges are "forward" (i<j => A[i][j]=1),
            # so forward edges in a perm sigma are positions where sigma(i) < sigma(i+1)
            # which is exactly the ascent count = Eulerian numbers
            a_k = [eulerian(n, k) for k in range(n)]
            print(f"  n={n}: a_k = {a_k} (= Eulerian numbers, since all edges forward)")

###############################################################################
# Detailed n=5 analysis with correct invariant computation
###############################################################################

def detailed_n5_analysis():
    print("\n" + "=" * 70)
    print("DETAILED n=5 ANALYSIS: Verifying THM-062 formula")
    print("=" * 70)

    n = 5
    d = n - 1  # = 4

    A5_eul = [eulerian(5, k) for k in range(5)]
    c2 = [inflated_coeff(2, 4, k) for k in range(5)]
    c0 = [inflated_coeff(0, 4, k) for k in range(5)]

    tournaments = generate_all_tournaments(n)

    # For each tournament, compute:
    # 1. The actual a_k distribution
    # 2. The cycle counts (c3 = cyclic triples, c5 = 5-cycles)
    # 3. Check THM-062 formula

    # In the formula, t3 and t5 are the number of DIRECTED cycles,
    # or perhaps the number of vertex-subsets forming cycles?
    # From the OCF: H(T) = sum_S 2^{cycles(S)} over subsets S
    # The invariants I(T) in the formula are the independence-polynomial related invariants.
    #
    # Let me figure out what t3 and t5 are by checking the formula against actual data.
    #
    # For n=5: the cycle compositions of {1,...,5} include:
    # - The single 5-cycle: f=0, parts=1, coefficient = 2^1 = 2
    # - A 3-cycle + 2 transitive: f=2, parts=1, coefficient = 2^1 = 2
    # Wait, "parts(I)" is the number of parts in the composition I.
    # For a single 3-cycle from a 3-element subset, the remaining 2 vertices form
    # a transitive pair. So the composition has 2 parts: (3-cycle, 2-path)?
    #
    # Actually, I think the formula involves the independent sets of the conflict graph.
    # Let me just try: t3 = number of 3-element subsets forming a directed 3-cycle,
    # t5 = number of 5-element subsets forming a directed 5-cycle (i.e., 0 or 1 for n=5,
    # actually can be more since there can be multiple directed cycles on same 5 vertices).
    #
    # For n=5: t5 = number of directed Hamiltonian cycles on all 5 vertices.
    # H(T) = number of directed Hamiltonian PATHS.
    #
    # Let me just compute both c3 (cyclic triple count) and directed-5-cycle count,
    # and see which combination makes the formula work.

    results = {}
    for A in tournaments:
        # Count cyclic triples (c3)
        c3_count = 0
        for triple in combinations(range(n), 3):
            i, j, k = triple
            # Check if i->j->k->i or i->k->j->i
            if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                c3_count += 1

        # Count directed 5-cycles (Hamiltonian cycles)
        hc_count = 0
        for perm in permutations(range(n)):
            if perm[0] != 0:  # fix vertex 0 to avoid rotation counting
                continue
            is_cycle = True
            for i_idx in range(n):
                if A[perm[i_idx]][perm[(i_idx+1) % n]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                hc_count += 1
        # hc_count now counts directed Ham cycles with vertex 0 fixed
        # Total directed Ham cycles = hc_count (since we fixed one vertex in the cycle)

        a_k = count_forward_edges_dist(A, n)

        key = (c3_count, hc_count)
        if key not in results:
            results[key] = {'count': 0, 'a_k_set': set(), 'a_k_example': a_k}
        results[key]['count'] += 1
        results[key]['a_k_set'].add(tuple(a_k))

    print(f"\nTotal tournaments: {len(tournaments)}")
    print(f"Distinct (c3, hc5) pairs: {len(results)}")

    # Check if a_k is uniquely determined by (c3, hc5)
    print("\nIs a_k uniquely determined by (c3, hc5)?")
    all_unique = True
    for key, val in sorted(results.items()):
        if len(val['a_k_set']) > 1:
            all_unique = False
            print(f"  (c3={key[0]}, hc5={key[1]}): {len(val['a_k_set'])} distinct a_k distributions")
    if all_unique:
        print("  YES - a_k is uniquely determined by (c3, hc5)")

    # Now try to match the formula
    # a_k(T) = A(5,k) + 2*c_k^{(2,4)}*t3 + 2*c_k^{(0,4)}*t5
    # where t3 = c3_count and t5 = hc5_count
    print("\nFormula verification: a_k = A(5,k) + 2*c_k^(2,4)*c3 + 2*c_k^(0,4)*hc5")
    print(f"  A(5,k) = {A5_eul}")
    print(f"  c_k^(2,4) = {c2}")
    print(f"  c_k^(0,4) = {c0}")

    all_match = True
    for key in sorted(results.keys()):
        c3_val, hc5_val = key
        actual = list(list(results[key]['a_k_set'])[0]) if len(results[key]['a_k_set']) == 1 else None
        if actual is None:
            print(f"  (c3={c3_val}, hc5={hc5_val}): Multiple a_k, skipping formula check")
            all_match = False
            continue
        predicted = [A5_eul[k] + 2*c2[k]*c3_val + 2*c0[k]*hc5_val for k in range(5)]
        match = (actual == predicted)
        if not match:
            all_match = False
        print(f"  (c3={c3_val}, hc5={hc5_val}): actual={actual}, predicted={predicted}, {'MATCH' if match else 'MISMATCH'}")

    if all_match:
        print("\n  ALL MATCH: THM-062 formula verified at n=5!")
    else:
        print("\n  MISMATCHES found. Trying alternative invariant definitions...")

    return results

###############################################################################
# Alternative: try with number of directed 3-cycles (not vertex sets)
###############################################################################

def detailed_n5_analysis_v2():
    print("\n" + "=" * 70)
    print("DETAILED n=5 ANALYSIS v2: Trying directed cycle counts")
    print("=" * 70)

    n = 5
    d = n - 1  # = 4

    A5_eul = [eulerian(5, k) for k in range(5)]
    c2 = [inflated_coeff(2, 4, k) for k in range(5)]
    c0 = [inflated_coeff(0, 4, k) for k in range(5)]

    tournaments = generate_all_tournaments(n)

    # Try: t3 = number of directed 3-cycles (= 2 * number of cyclic triples)
    # t5 = number of directed 5-cycles (Hamiltonian cycles)

    results = {}
    for A in tournaments:
        # Count cyclic triples
        c3_subsets = 0
        for triple in combinations(range(n), 3):
            i, j, k = triple
            if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                c3_subsets += 1

        # Count directed 5-cycles
        hc5 = 0
        for perm in permutations(range(n)):
            if perm[0] != 0:
                continue
            is_cycle = True
            for idx in range(n):
                if A[perm[idx]][perm[(idx+1) % n]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                hc5 += 1

        a_k = count_forward_edges_dist(A, n)

        key = (c3_subsets, hc5)
        if key not in results:
            results[key] = {'count': 0, 'a_k_list': []}
        results[key]['count'] += 1
        results[key]['a_k_list'].append(tuple(a_k))

    # Check multiple formulas
    print("\nTrying different scalings:")

    for t3_scale_name, t3_scale in [("c3_subsets", 1), ("2*c3_subsets", 2)]:
        for t5_scale_name, t5_scale in [("hc5", 1), ("hc5/2", 0.5)]:
            all_ok = True
            for key in sorted(results.keys()):
                c3_val, hc5_val = key
                t3_val = c3_val * t3_scale
                t5_val = hc5_val * t5_scale
                predicted = [A5_eul[k] + 2*c2[k]*t3_val + 2*c0[k]*t5_val for k in range(5)]
                actual_set = set(results[key]['a_k_list'])
                for actual in actual_set:
                    if list(actual) != [int(p) for p in predicted]:
                        all_ok = False
                        break
                if not all_ok:
                    break
            status = "WORKS" if all_ok else "fails"
            print(f"  t3={t3_scale_name}, t5={t5_scale_name}: {status}")

    # Print actual data
    print("\nActual data (c3_subsets, hc5) -> a_k distributions:")
    for key in sorted(results.keys()):
        c3_val, hc5_val = key
        actual_set = set(results[key]['a_k_list'])
        for a in sorted(actual_set):
            predicted_v1 = [A5_eul[k] + 2*c2[k]*c3_val + 2*c0[k]*hc5_val for k in range(5)]
            print(f"  (c3={c3_val}, hc5={hc5_val}): actual={list(a)}, pred_v1={predicted_v1}, cnt={results[key]['count']}")

###############################################################################
# Main
###############################################################################

if __name__ == "__main__":
    print("POSITIVITY CONSTRAINTS FROM THM-062")
    print("Forward-Edge Distribution Analysis")
    print()

    # Part 1: Compute coefficients at n=5
    A5, c2_5, c0_5 = part1_n5()

    # Part 2-3: Feasible region
    part2_feasible_n5(A5, c2_5, c0_5)

    # Detailed n=5 analysis
    detailed_n5_analysis()
    detailed_n5_analysis_v2()

    # Parts 4-5: Actual data
    part45_actual_n5(A5, c2_5, c0_5)

    # Part 6: n=7
    part6_n7()

    # Part 7: Zero check
    part7_zero_check()
