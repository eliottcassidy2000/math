#!/usr/bin/env python3
"""
vitali_tournament_structure.py -- kind-pasteur-2026-03-13-S61

THE VITALI SET AND HIDDEN HIGHER-DIMENSIONAL STRUCTURE IN TOURNAMENTS

Core idea: The Vitali set is constructed by:
  1. Partitioning [0,1] into equivalence classes under Q-translation
  2. Choosing one representative per class (Axiom of Choice)
  3. The resulting set V is non-measurable

Tournaments have an ANALOGOUS structure:
  1. The space of labeled tournaments on n vertices is {0,1}^{C(n,2)}
  2. S_n acts by relabeling, partitioning into isomorphism orbits
  3. A "canonical form" picks one representative per orbit

The connection deepens when we look at n=0,1,2:
  - n=0: Empty tournament. H=1. The space is a POINT. Dimension 0.
  - n=1: Single vertex. H=1. The space is a POINT. Dimension 0.
  - n=2: One edge. H=1. The space is {0,1}^1 = 2 points. Dimension 1.
           Both labeled tournaments are isomorphic. One orbit.
  - n=3: Three edges. C(3,2)=3. Space is {0,1}^3 = cube. Dimension 3.
           8 labeled tournaments, 2 orbits (transitive:6, cycle:2).

The "hidden higher-dimensional structure" is: the orientation space {0,1}^{C(n,2)}
has dimension growing as n^2/2, but the observable invariants (H, spectra, etc.)
project this space onto low-dimensional manifolds with FRACTAL-LIKE level sets.

This script explores:
1. The binary encoding and its connection to the unit interval [0,1]
2. The orbit/quotient structure as a Vitali analogue
3. The {-1, 0, +1} skew-symmetric representation (Z/3Z grading)
4. How the degenerate cases n=0,1,2 encode axioms of the theory
5. The measure-theoretic structure: are H-level sets "measurable"?

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import permutations, combinations
from collections import defaultdict


def tournament_to_binary(A, n):
    """Encode tournament adjacency matrix as binary string.
    Convention: bit position (i,j) for i<j, ordered lexicographically.
    bit = 1 if i->j, bit = 0 if j->i."""
    bits = 0
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j]:
                bits |= (1 << pos)
            pos += 1
    return bits


def binary_to_tournament(bits, n):
    """Decode binary string to adjacency matrix."""
    m = n * (n - 1) // 2
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def apply_permutation(A, n, perm):
    """Apply vertex permutation to tournament. perm[i] = new label of vertex i."""
    B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                B[perm[i]][perm[j]] = 1
    return B


def canonical_form(bits, n):
    """Compute canonical form = minimum binary encoding over all relabelings."""
    A = binary_to_tournament(bits, n)
    min_bits = bits
    for perm in permutations(range(n)):
        B = apply_permutation(A, n, list(perm))
        b = tournament_to_binary(B, n)
        if b < min_bits:
            min_bits = b
    return min_bits


def count_ham_paths(A, n):
    """Count Hamiltonian paths using Held-Karp DP."""
    if n == 0:
        return 1  # empty permutation
    if n == 1:
        return 1  # trivial path

    # dp[mask][v] = number of paths visiting vertices in mask, ending at v
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]

    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def skew_symmetric_representation(A, n):
    """Convert tournament to skew-symmetric matrix sigma_{ij} in {-1, 0, +1}.
    sigma(i,j) = +1 if i->j, -1 if j->i, 0 if i=j.
    This lives in the Z/3Z-graded space."""
    sigma = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                sigma[i][j] = 0
            elif A[i][j]:
                sigma[i][j] = 1
            else:
                sigma[i][j] = -1
    return sigma


def sigma_to_vector(sigma, n):
    """Extract the upper triangle of sigma as a vector in {-1, +1}^{C(n,2)}.
    This is the tournament's embedding in the n(n-1)/2 dimensional hypercube."""
    vec = []
    for i in range(n):
        for j in range(i+1, n):
            vec.append(sigma[i][j])
    return tuple(vec)


# ========================================================================
# ANALYSIS 1: THE DEGENERATE CASES (n=0,1,2) AS AXIOMS
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: DEGENERATE CASES AS AXIOMS OF TOURNAMENT THEORY")
print("=" * 70)

print("""
The degenerate tournaments n=0,1,2 encode the AXIOMS of the theory:

n=0: The EMPTY tournament
  - Vertices: none. Edges: none.
  - H(T_0) = 1 (empty permutation = identity in S_0)
  - Omega(T_0) = empty graph
  - I(empty, 2) = 1 (empty product)
  - OCF: H = 1 + 0 = 1 (no cycles exist)
  - AXIOM: The empty sum of odd cycles contributes 1 to H.
  - This is the MULTIPLICATIVE IDENTITY for H under disjoint union.
  - Dimension of orientation space: C(0,2) = 0. A single point.

n=1: The TRIVIAL tournament
  - Vertices: {0}. Edges: none.
  - H(T_1) = 1 (unique permutation of 1 element)
  - Omega(T_1) = empty graph (no odd cycles on 1 vertex)
  - OCF: H = 1
  - AXIOM: A single vertex contributes 1 path (itself).
  - Dimension: C(1,2) = 0. Still a single point.
  - The tournament space at n=0 and n=1 are INDISTINGUISHABLE
    from the orientation-space perspective. Both are 0-dimensional.

n=2: The FIRST NONTRIVIAL case
  - Vertices: {0,1}. Edges: one directed edge.
  - H(T_2) = 1 (unique Hamiltonian path: winner -> loser)
  - Orientation space: {0,1}^1 = two points (0->1 or 1->0)
  - Both are isomorphic! S_2 swaps them.
  - Omega(T_2) = empty (no odd cycles on 2 vertices)
  - OCF: H = 1
  - AXIOM: The first edge introduces DIRECTION but no cycles.
  - This is where ANTISYMMETRY enters: A[i][j] + A[j][i] = 1.
  - Dimension: C(2,2) = 1. The orientation space is a LINE SEGMENT.
  - The S_2 action identifies the two endpoints: Z_2 quotient.
""")

# ========================================================================
# ANALYSIS 2: THE BINARY HYPERCUBE AND VITALI QUOTIENT
# ========================================================================
print("=" * 70)
print("ANALYSIS 2: BINARY HYPERCUBE AND VITALI QUOTIENT STRUCTURE")
print("=" * 70)

for n in range(2, 7):
    m = n * (n - 1) // 2  # dimension of orientation space
    total = 1 << m         # total labeled tournaments

    # Compute all orbits
    orbits = defaultdict(list)
    for bits in range(total):
        canon = canonical_form(bits, n)
        orbits[canon].append(bits)

    n_orbits = len(orbits)
    orbit_sizes = sorted([len(v) for v in orbits.values()], reverse=True)

    # Compute H for each orbit representative
    orbit_H = {}
    for canon, members in orbits.items():
        A = binary_to_tournament(canon, n)
        H = count_ham_paths(A, n)
        orbit_H[canon] = H

    H_values = sorted(set(orbit_H.values()))

    # H-level sets
    H_level_sizes = defaultdict(int)
    H_level_orbits = defaultdict(int)
    for canon, H in orbit_H.items():
        H_level_sizes[H] += len(orbits[canon])
        H_level_orbits[H] += 1

    print(f"\nn={n}: dimension={m}, |space|={total}")
    print(f"  S_{n} orbits: {n_orbits}")
    print(f"  Orbit sizes: {orbit_sizes}")
    print(f"  |S_n| = {math.factorial(n)}")

    # Vitali analogy
    print(f"\n  VITALI ANALOGY:")
    print(f"    Base space: {{0,1}}^{m}  ({total} points)")
    print(f"    Group: S_{n} acting by relabeling ({math.factorial(n)} elements)")
    print(f"    Equivalence classes: {n_orbits} orbits")
    print(f"    'Vitali set' (canonical reps): {n_orbits} points out of {total}")
    print(f"    Selection ratio: {n_orbits/total:.6f}")
    print(f"    Average orbit size: {total/n_orbits:.1f}")

    # Measure theory
    print(f"\n  H-LEVEL SET STRUCTURE:")
    for h in sorted(H_level_sizes.keys()):
        count = H_level_sizes[h]
        norb = H_level_orbits[h]
        pct = 100 * count / total
        print(f"    H={h}: {count} tournaments ({pct:.2f}%), {norb} orbit(s)")

    # Check: H-level sets are S_n-invariant (unions of orbits)
    # This is trivially true since H is a graph invariant
    print(f"    H-level sets are MEASURABLE: unions of full S_n-orbits")

    if n <= 5:
        # Sigma vectors: embed in {-1,+1}^m
        print(f"\n  SIGMA EMBEDDING in {{-1,+1}}^{m}:")
        sigma_vecs = {}
        for bits in range(total):
            A = binary_to_tournament(bits, n)
            sigma = skew_symmetric_representation(A, n)
            vec = sigma_to_vector(sigma, n)
            sigma_vecs[bits] = vec
            H = count_ham_paths(A, n)

        # The sigma vectors live on vertices of a hypercube centered at origin
        # Each coordinate is +/-1, so ||sigma||^2 = m for all tournaments
        # H is a function on this hypercube

        # Check: does H have a nice expansion in the Fourier basis of {-1,+1}^m?
        # The Walsh-Hadamard transform of H
        # H(sigma) = sum_S hat(H)(S) * prod_{(i,j) in S} sigma_{ij}

        # Compute H for all sigma vectors
        H_on_cube = {}
        for bits in range(total):
            A = binary_to_tournament(bits, n)
            H = count_ham_paths(A, n)
            H_on_cube[sigma_vecs[bits]] = H

        # Walsh-Hadamard coefficients
        # hat(H)(S) = (1/2^m) * sum_sigma H(sigma) * prod_{(i,j) in S} sigma_{ij}
        print(f"\n  WALSH-HADAMARD FOURIER ANALYSIS OF H:")
        print(f"    H as function on {{-1,+1}}^{m}")

        # Compute all Fourier coefficients
        all_subsets = []
        for k in range(m + 1):
            for subset in combinations(range(m), k):
                all_subsets.append(subset)

        fourier_coeffs = {}
        for S in all_subsets:
            coeff = 0
            for vec, H in H_on_cube.items():
                chi_S = 1
                for idx in S:
                    chi_S *= vec[idx]
                coeff += H * chi_S
            coeff /= total
            if abs(coeff) > 1e-10:
                fourier_coeffs[S] = coeff

        # Report non-zero coefficients by degree
        by_degree = defaultdict(list)
        for S, c in sorted(fourier_coeffs.items(), key=lambda x: len(x[0])):
            by_degree[len(S)].append((S, c))

        for deg in sorted(by_degree.keys()):
            items = by_degree[deg]
            total_energy = sum(c**2 for _, c in items)
            print(f"    Degree {deg}: {len(items)} non-zero terms, energy = {total_energy:.4f}")
            if deg <= 2 and n <= 4:
                for S, c in items:
                    # Convert index back to edge
                    edges = []
                    pos = 0
                    for i in range(n):
                        for j in range(i+1, n):
                            if pos in S:
                                edges.append(f"({i},{j})")
                            pos += 1
                    print(f"      S={','.join(edges) if edges else '{}'}: hat(H) = {c:.4f}")

        total_energy = sum(c**2 for c in fourier_coeffs.values())
        print(f"    Total Parseval energy: {total_energy:.4f}")
        print(f"    Expected: E[H^2] = {sum(H**2 for H in H_on_cube.values())/total:.4f}")


# ========================================================================
# ANALYSIS 3: THE Z/3Z GRADING AND {-1,0,+1} STRUCTURE
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: THE Z/3Z GRADING AND NON-MEASURABILITY")
print("=" * 70)

print("""
The skew-symmetric matrix sigma_{ij} takes values in {-1, 0, +1}.
These are the elements of Z/3Z (identifying -1 with 2).

In the Vitali construction:
  - [0,1] is partitioned by Q-cosets: x ~ y iff x-y in Q
  - Representatives form V, a non-measurable set
  - Key: Q is DENSE in R, so cosets are "everywhere"

In tournaments:
  - {-1,+1}^{C(n,2)} is partitioned by S_n-orbits: T ~ T' iff T ~ T'
  - Representatives form the "canonical tournament set"
  - Key: S_n acts TRANSITIVELY on some orbits, creating structure

The Z/3Z grading {-1,0,+1} reveals:
  - sigma_{ij} + sigma_{ji} = 0 (antisymmetry = "completeness")
  - sigma_{ii} = 0 (diagonal = "irreflexivity")
  - These TWO constraints define the tournament variety in Z^{n^2}

At n=2: sigma = [[0, s], [-s, 0]] where s in {-1,+1}
  The ONLY degree of freedom is the sign of s.
  This is a 0-dimensional sphere S^0 = {-1, +1}.

At n=3: sigma has 3 free entries (upper triangle).
  Constraint: each 3-element minor has det != 0 (strong tournament)
  Actually: a tournament IS a complete antisymmetric digraph.
  The space is {-1,+1}^3, a 0-skeleton of the 3-cube.
""")

# Now the deeper question: what "non-measurable" structure exists?
# In finite spaces, everything is measurable. But in the LIMIT n->inf:
# the tournament space becomes a topological/measure space where
# certain subsets (like "Paley-like" tournaments) become singular.

print("THE VITALI PHENOMENON IN TOURNAMENT LIMITS:")
print()

for n in range(2, 7):
    m = n * (n - 1) // 2
    total = 1 << m

    # Compute orbits and H
    orbits = defaultdict(list)
    H_map = {}
    for bits in range(total):
        canon = canonical_form(bits, n)
        orbits[canon].append(bits)
        A = binary_to_tournament(bits, n)
        H_map[bits] = count_ham_paths(A, n)

    # The "measure" question: what fraction of tournaments have H = max?
    max_H = max(H_map.values())
    max_H_count = sum(1 for H in H_map.values() if H == max_H)

    # What fraction are H-maximal?
    frac_max = max_H_count / total

    # The "concentration" question: what is Var(H)/E[H]^2?
    mean_H = sum(H_map.values()) / total
    var_H = sum((H - mean_H)**2 for H in H_map.values()) / total

    print(f"  n={n}: max_H={max_H}, #maximal={max_H_count}/{total} ({100*frac_max:.2f}%)")
    print(f"    E[H]={mean_H:.2f}, Var(H)={var_H:.2f}, CV={var_H**0.5/mean_H:.4f}")
    print(f"    #distinct H values: {len(set(H_map.values()))}")


# ========================================================================
# ANALYSIS 4: THE HIDDEN DIMENSION TRANSITION
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: HIDDEN DIMENSIONAL TRANSITIONS AT n=0,1,2,3")
print("=" * 70)

print("""
The "hidden higher-dimensional structure" emerges from the DIMENSION JUMPS:

  n=0 -> n=1: dim 0 -> dim 0 (no new edges possible)
    INVARIANT: H=1 in both cases.
    TRANSITION: Nothing changes. This is the GROUND STATE.

  n=1 -> n=2: dim 0 -> dim 1 (one edge appears)
    INVARIANT: H=1 still! The first edge doesn't create cycles.
    TRANSITION: DIRECTION is born. The tournament lives on S^0 = {-1,+1}.
    The Z/2Z symmetry (swap 0<->1) identifies the two orientations.
    ANALOGY: In the Vitali construction, this is like [0,1] being
    identified under x <-> 1-x. A single orbit.

  n=2 -> n=3: dim 1 -> dim 3 (two more edges)
    INVARIANT: H can now be 1 or 3!
    TRANSITION: The FIRST ODD CYCLE (3-cycle) appears.
    The 3-cube {-1,+1}^3 has 8 vertices. S_3 (6 elements) acts.
    Two orbits: transitive (H=1, 6 copies) and cyclic (H=3, 2 copies).
    ANALOGY: The Vitali quotient now has 2 classes. The cyclic
    tournament is the "special" element (higher H), analogous
    to the irrational representative in Vitali's construction.

    KEY: At n=3, the conflict graph Omega(T) is either:
      - empty (transitive tournament, H=1)
      - single vertex (3-cycle tournament, H=3=1+2*1)
    The cycle creates a new "dimension" in the H-decomposition.

  n=3 -> n=4: dim 3 -> dim 6 (three more edges)
    The dimension DOUBLES. The tournament space explodes from 8 to 64.
    S_4 (24 elements) acts, giving 4 orbits.
    H-values: {1, 3, 5} — now H=5 appears (the regular tournament).

  n=4 -> n=5: dim 6 -> dim 10 (four more edges)
    1024 tournaments, 12 orbits.
    H-values: {1, 3, 5, 9, 11, 13, 15} — full odd spectrum opens up.
""")

# Detailed transition analysis
print("DIMENSION JUMP TABLE:")
print(f"{'n':>3} {'dim':>5} {'|space|':>10} {'|S_n|':>8} {'orbits':>8} {'H-values':>20}")
print("-" * 60)

for n in range(0, 7):
    m = n * (n - 1) // 2
    total = 1 << m
    sn = math.factorial(n)

    if n <= 6:
        orbits = defaultdict(list)
        H_vals = set()
        for bits in range(total):
            canon = canonical_form(bits, n)
            orbits[canon].append(bits)
            A = binary_to_tournament(bits, n)
            H = count_ham_paths(A, n)
            H_vals.add(H)
        n_orb = len(orbits)
        h_str = str(sorted(H_vals))
    else:
        n_orb = "?"
        h_str = "?"

    print(f"{n:>3} {m:>5} {total:>10} {sn:>8} {str(n_orb):>8} {h_str:>20}")


# ========================================================================
# ANALYSIS 5: THE OVERLAP WEIGHT AND VITALI COSET STRUCTURE
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: OVERLAP WEIGHT AS VITALI COSET MEASURE")
print("=" * 70)

print("""
Connection between Vitali cosets and overlap weight:

In the Vitali construction:
  - mu(V + q) = mu(V) for all q in Q (translation invariance)
  - sum_{q in Q} mu(V + q) = mu([0,1]) = 1
  - If mu(V) = 0, the sum is 0 != 1. Contradiction.
  - If mu(V) > 0, the sum diverges != 1. Contradiction.

In tournaments:
  - The overlap weight w(C) = |edges shared with neighbors in Omega|
  - The "Vitali invariant" is: is there a UNIFORM measure on cycles
    such that the overlap structure is translation-invariant?

For CIRCULANT tournaments: YES! Z_p acts by translation, and:
  - All vertices are equivalent (vertex-transitivity)
  - Cycle vertex sets are cosets: C + i = C shifted by i
  - The co-occurrence matrix is CIRCULANT

For NON-circulant tournaments: The Vitali phenomenon appears!
  - Some vertices participate in more cycles than others
  - The "co-occurrence variance" measures non-uniformity
  - Paley: variance = 0 (perfect uniformity) [Vitali-MEASURABLE]
  - Interval: variance > 0 (non-uniform) [Vitali-like non-uniformity]
""")

# Compute co-occurrence structure for small n
for n in [3, 4, 5]:
    m = n * (n - 1) // 2
    total = 1 << m

    print(f"\n  --- n={n}: Co-occurrence analysis ---")

    orbits = defaultdict(list)
    for bits in range(total):
        canon = canonical_form(bits, n)
        orbits[canon].append(bits)

    for canon, members in sorted(orbits.items()):
        A = binary_to_tournament(canon, n)
        H = count_ham_paths(A, n)

        # Count 3-cycles
        c3_sets = []
        for a, b, c in combinations(range(n), 3):
            if A[a][b] and A[b][c] and A[c][a]:
                c3_sets.append(frozenset([a, b, c]))
            if A[a][c] and A[c][b] and A[b][a]:
                c3_sets.append(frozenset([a, b, c]))
        c3_sets = list(set(c3_sets))

        # Vertex participation in cycles
        v_count = [0] * n
        for cs in c3_sets:
            for v in cs:
                v_count[v] += 1

        # Score sequence
        scores = [sum(A[i]) for i in range(n)]

        if len(c3_sets) > 0:
            # Co-occurrence variance
            co_var = sum((c - sum(v_count)/n)**2 for c in v_count) / n
            print(f"    canon={canon:>3d} (orbit={len(members):>2d}), H={H:>3d}, "
                  f"scores={sorted(scores)}, c3={len(c3_sets)}, "
                  f"v_part={v_count}, co_var={co_var:.2f}")
        else:
            print(f"    canon={canon:>3d} (orbit={len(members):>2d}), H={H:>3d}, "
                  f"scores={sorted(scores)}, c3=0 (transitive)")


# ========================================================================
# ANALYSIS 6: THE {2,1,0} INTERPRETATION — ADJACENCY AS TERNARY
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: THE {2,1,0} INTERPRETATION")
print("=" * 70)

print("""
The user asks about "tournaments of n 2, 1 and 0."
Three possible interpretations, ALL illuminating:

INTERPRETATION A: Score values {0, 1, 2}
  At n=3, scores are permutations of (0, 1, 2):
    - Transitive: scores (0, 1, 2) — a dominance hierarchy
    - 3-cycle: scores (1, 1, 1) — perfect balance
  The score SUM is always C(n,2) = 3.
  Score = out-degree = "power" of each vertex.
  {0, 1, 2} = {loser, neutral, winner} at n=3.

INTERPRETATION B: Overlap multiplicities {0, 1, 2}
  Two cycles can share 0, 1, or 2 vertices (never 3, since
  3 shared vertices = same cycle for 3-cycles).
  The overlap weight W[i,j] in {0, 1, 2} determines whether
  cycles are disjoint (0), adjacent (1-2), or identical (3).
  This is the CONFLICT GRAPH'S edge weighting.

INTERPRETATION C: The ternary "tournament state" per vertex pair
  For (i,j): 0 = "no edge" (diagonal), 1 = i->j, 2 = j->i
  This gives a TERNARY encoding: each pair has state in Z/3Z.
  A tournament is a function f: C(n,2) -> {1, 2} with f(j,i) = 3 - f(i,j).
  Including the diagonal: f: [n]^2 -> {0, 1, 2}.

  THIS IS THE TRIT REPRESENTATION:
    n=0: 0 trits. 1 state.
    n=1: 1 trit (the self-loop = 0). 1 state.
    n=2: 4 trits, but only 1 free (others constrained). 2 states = 1 free trit.
    n=3: 9 trits, but only 3 free. 8 states = 2^3. (NOT 3^3!)

  The constraint sigma_{ij} + sigma_{ji} = 0, sigma_{ii} = 0
  kills the ternary freedom: the actual DOF is BINARY, not ternary.

  But the APPEARANCE is ternary: {beat, lose, self} = {2, 0, 1}
  or equivalently {+1, -1, 0} in the skew-symmetric representation.
""")

# Compute the ternary structure concretely
for n in range(0, 5):
    print(f"\n  n={n}:")
    print(f"    Full matrix: {n}x{n} = {n*n} entries")
    print(f"    Diagonal (=0): {n} entries")
    print(f"    Antisymmetric pairs: {n*(n-1)//2}")
    print(f"    Free binary DOF: {n*(n-1)//2}")
    print(f"    Ternary 'apparent' DOF: {n*n} entries in {{0,1,2}}")
    print(f"    Constraint ratio: {n*n - n*(n-1)//2} constrained out of {n*n}")
    if n > 0:
        print(f"    Information: {n*(n-1)//2} bits = {n*(n-1)//2 * math.log(2, 3):.2f} trits")


# ========================================================================
# ANALYSIS 7: VITALI NON-MEASURABILITY IN THE TOURNAMENT LIMIT
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: ASYMPTOTIC NON-MEASURABILITY")
print("=" * 70)

print("""
The Vitali set is non-measurable because:
  1. The group Q is infinite, creating infinitely many cosets
  2. Translation invariance forces mu(V) = const for all cosets
  3. But sum of infinitely many equal non-negative numbers is 0 or infinity

In tournaments, the "limiting" version:
  1. S_n grows as n! but |space| grows as 2^{C(n,2)} = 2^{n(n-1)/2}
  2. Average orbit size = |space|/|orbits| ~ n! (by Burnside)
  3. Number of orbits ~ 2^{n(n-1)/2} / n!

The ratio 2^{n(n-1)/2} / n! grows SUPER-EXPONENTIALLY:
""")

for n in range(1, 16):
    m = n * (n - 1) // 2
    log2_space = m
    log2_sn = sum(math.log2(k) for k in range(1, n+1))
    log2_orbits = log2_space - log2_sn  # approximate

    print(f"  n={n:>2d}: log2|space|={log2_space:>5d}, log2(n!)={log2_sn:>7.1f}, "
          f"log2(orbits)~{log2_orbits:>8.1f}")

print("""
At n=15: ~53 bits for S_n vs ~105 bits for space => ~2^52 orbits.
The "canonical representative set" has density ~2^{-53} in the space.
This is extremely sparse — like the Vitali set's representatives
being measure-zero if measurable.

DEEP ANALOGY: As n -> infinity:
  - The "Vitali density" = |orbits|/|space| -> 0 super-exponentially
  - But |orbits| -> infinity super-exponentially
  - This is EXACTLY the Vitali paradox: each representative has
    "measure zero" but there are "too many" of them.

  - For tournaments: the H-maximizer (Paley/interval) becomes
    exponentially rare in the total space (density -> 0)
  - But it carries an exponentially large fraction of total H
  - This is a "tournament Vitali phenomenon": the maximizer is
    non-measurable in the sense that its contribution to E[H]
    is disproportionate to its measure.
""")

# Compute the "H-concentration" phenomenon
print("H-CONCENTRATION (fraction of total H in maximizers):")
for n in range(2, 7):
    m = n * (n - 1) // 2
    total = 1 << m

    H_vals = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        H_vals.append(H)

    total_H = sum(H_vals)
    max_H = max(H_vals)
    max_count = H_vals.count(max_H)
    max_fraction = max_count * max_H / total_H

    print(f"  n={n}: max_H={max_H}, count={max_count}/{total} ({100*max_count/total:.1f}%), "
          f"carries {100*max_fraction:.1f}% of total H-sum")


print("\n" + "=" * 70)
print("DONE.")
