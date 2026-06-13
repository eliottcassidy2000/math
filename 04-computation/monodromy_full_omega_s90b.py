#!/usr/bin/env python3
"""
monodromy_full_omega_s90b.py — Monodromy structure + full Ω computation
opus-2026-03-15-S90b

1. Compute the FULL conflict graph Ω (all odd cycles, not just 3-cycles)
   for near-transitive tournaments to verify H = 2^{n-2}+1 via OCF.
2. Compute the monodromy representation of near-transitive tournaments.
3. Extend the simplicial Rédei to understand the SPECTRAL STRUCTURE
   of near-transitive tournaments (eigenvalues of their adjacency matrices).
4. Connect to k-nacci by computing the number of odd cycles by length.
"""

import numpy as np
from itertools import permutations, combinations
from math import factorial, comb, log, pi

# ======================================================================
# PART 1: FULL Ω FOR NEAR-TRANSITIVE — VERIFY VIA OCF
# ======================================================================
print("=" * 70)
print("PART 1: FULL CONFLICT GRAPH Ω FOR NEAR-TRANSITIVE TOURNAMENTS")
print("=" * 70)

def build_near_transitive(n):
    """Build the canonical near-transitive tournament: 0>1>...>n-1, reverse (n-1)→0."""
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            T[i][j] = 1
    T[0][n-1] = 0
    T[n-1][0] = 1
    return T

def find_all_odd_cycles(T, n):
    """Find ALL directed odd cycles in tournament T."""
    cycles = set()
    for length in range(3, n+1, 2):  # odd lengths 3, 5, 7, ...
        for subset in combinations(range(n), length):
            S = list(subset)
            # Check if there's a directed cycle of this length on these vertices
            # A directed cycle visits each vertex exactly once
            # We only need to find ONE such cycle per vertex set
            for start_perm in permutations(S):
                is_cycle = True
                for p in range(length):
                    if T[start_perm[p]][start_perm[(p+1) % length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(frozenset(S))
                    break  # One cycle per vertex set is enough
    return list(cycles)

def conflict_graph_full(cycles):
    """Build full conflict graph: vertices=cycles, edges=share vertex."""
    nc = len(cycles)
    adj = [[0]*nc for _ in range(nc)]
    for a in range(nc):
        for b in range(a+1, nc):
            if cycles[a] & cycles[b]:
                adj[a][b] = adj[b][a] = 1
    return adj

def independence_poly(adj, nc):
    """Independence polynomial by brute force."""
    alpha = [0] * (nc + 1)
    for mask in range(2**nc):
        sel = [i for i in range(nc) if (mask >> i) & 1]
        ok = True
        for a in range(len(sel)):
            for b in range(a+1, len(sel)):
                if adj[sel[a]][sel[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            alpha[len(sel)] += 1
    return alpha

print("""
For the near-transitive tournament (total order 0>1>...>n-1, arc (n-1)→0):

Odd cycles of length 2k+1 use vertices {0, i₁, ..., i_{2k-1}, n-1} where
{i₁,...,i_{2k-1}} ⊂ {1,...,n-2}. The cycle is:
  0 → i₁ → i₂ → ... → i_{2k-1} → (n-1) → 0

For this to be a valid cycle, we need:
  0→i₁ (i₁ < n-1, so this is a forward arc ✓)
  i_{j} → i_{j+1} for each j (forward arcs if i_j < i_{j+1} ✓)
  i_{2k-1} → (n-1) (forward arc ✓)
  (n-1) → 0 (reversed arc ✓)

So we need the i_j to be in INCREASING order! The number of such cycles
of length 2k+1 is C(n-2, 2k-1) (choose 2k-1 vertices from {1,...,n-2}).

Total odd cycles = Σ_{k=1}^{⌊(n-1)/2⌋} C(n-2, 2k-1) = 2^{n-3}
(sum of odd-indexed binomial coefficients of n-2).
""")

for n in range(4, 10):
    T = build_near_transitive(n)

    if n <= 8:
        cycles = find_all_odd_cycles(T, n)
    else:
        # For n=9, brute force is too slow for 7-cycles and 9-cycles
        cycles = None

    # Predicted cycle counts by length
    predicted_total = 0
    cycle_by_length = {}
    for k in range(1, n//2 + 1):
        length = 2*k + 1
        if 2*k - 1 <= n - 2:
            count = comb(n-2, 2*k-1)
            cycle_by_length[length] = count
            predicted_total += count

    print(f"\nn={n}: predicted total odd cycles = {predicted_total} = 2^{n-3} = {2**(n-3)}")
    print(f"  By length: {cycle_by_length}")

    if cycles is not None:
        # Verify cycle count
        actual_by_length = {}
        for c in cycles:
            l = len(c)
            actual_by_length[l] = actual_by_length.get(l, 0) + 1
        print(f"  Actual:    {actual_by_length}")
        print(f"  Total actual: {len(cycles)}, predicted: {predicted_total}, match: {len(cycles) == predicted_total}")

        # Verify all cycles share vertex 0 AND vertex n-1
        all_share_0 = all(0 in c for c in cycles)
        all_share_n1 = all((n-1) in c for c in cycles)
        print(f"  All share vertex 0: {all_share_0}")
        print(f"  All share vertex {n-1}: {all_share_n1}")

        # Build and analyze Ω
        adj = conflict_graph_full(cycles)
        nc = len(cycles)
        edges = sum(sum(row) for row in adj) // 2
        max_edges = nc * (nc-1) // 2
        print(f"  Ω: {nc} vertices, {edges}/{max_edges} edges (complete? {edges == max_edges})")

        # Compute I(Ω, 2)
        if nc <= 20:
            alpha = independence_poly(adj, nc)
            H_ocf = sum(alpha[k] * 2**k for k in range(nc+1))
            print(f"  I(Ω, 2) = {H_ocf} (expected H = {2**(n-2)+1})")
            print(f"  α vector: {alpha[:6]}")
        else:
            # For complete graph K_m: I(K_m, 2) = 1 + 2m
            H_predicted = 1 + 2 * nc
            print(f"  Ω = K_{nc}, I(K_{nc}, 2) = {H_predicted} (expected H = {2**(n-2)+1})")

# ======================================================================
# PART 2: MONODROMY GROUP STRUCTURE
# ======================================================================
print()
print("=" * 70)
print("PART 2: MONODROMY GROUP OF NEAR-TRANSITIVE TOURNAMENTS")
print("=" * 70)

print("""
The monodromy group of the near-transitive tournament is (Z/2)^{n-2}.

GEOMETRIC INTERPRETATION:
The transitive core defines a total order σ = (v₁ > ... > vₙ).
The unique simplicial path follows σ.
The single reversed arc (min→max) creates a "portal" between the
bottom and top of the order.

Each Hamiltonian path using the portal has the form:
  (S in ↑ order, n-1, 0, {1,...,n-2}\\S in ↑ order)
for some subset S ⊆ {1,...,n-2}.

The MONODROMY ACTION: each intermediate vertex k ∈ {1,...,n-2}
has a binary choice: "before the portal" (k ∈ S) or "after" (k ∉ S).
This gives a (Z/2)^{n-2} action on the set of portal-using paths.

The group structure:
  - Generator g_k: flip vertex k between S and complement(S)
  - Relations: g_k² = 1, g_k g_j = g_j g_k (abelian)
  - The group is the n-2 dimensional cube / Boolean lattice
""")

# Verify: the paths biject with elements of (Z/2)^{n-2}
for n in [5, 6, 7]:
    T = build_near_transitive(n)

    # Find all Hamiltonian paths
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for p in range(n-1):
            if T[perm[p]][perm[p+1]] != 1:
                ok = False
                break
        if ok:
            paths.append(perm)

    # Classify: simplicial (doesn't use portal) vs portal
    portal_paths = []
    simplicial = None
    for p in paths:
        uses_portal = False
        for k in range(n-1):
            if p[k] == n-1 and p[k+1] == 0:
                uses_portal = True
                break
        if uses_portal:
            # Extract S = vertices before the portal
            portal_pos = next(k for k in range(n-1) if p[k] == n-1 and p[k+1] == 0)
            S = frozenset(p[:portal_pos])  # vertices before n-1
            portal_paths.append((p, S))
        else:
            simplicial = p

    print(f"\nn={n}: {len(paths)} total paths, {len(portal_paths)} portal, 1 simplicial")
    print(f"  Simplicial path: {simplicial}")
    print(f"  Portal path subsets S (vertices before portal):")

    # Show the (Z/2)^{n-2} structure
    S_sets = [S for _, S in portal_paths]
    # Convert to binary vectors
    for p, S in sorted(portal_paths, key=lambda x: sum(1 << k for k in x[1])):
        binary = ''.join('1' if k in S else '0' for k in range(1, n-1))
        print(f"    S={str(sorted(S)):15s} binary={binary} path={p}")

# ======================================================================
# PART 3: ADJACENCY MATRIX EIGENVALUES OF NEAR-TRANSITIVE
# ======================================================================
print()
print("=" * 70)
print("PART 3: EIGENVALUES OF NEAR-TRANSITIVE ADJACENCY MATRIX")
print("=" * 70)

print("""
The adjacency matrix of the near-transitive tournament is:
  A = J_upper - e_0 e_{n-1}^T + e_{n-1} e_0^T
where J_upper is the strict upper triangular all-ones matrix
(the transitive tournament).

The transitive tournament has eigenvalues:
  λ_k = (n-1)/2 * exp(2πik/n) ?? NO — the transitive tournament
  is NOT circulant, so DFT doesn't apply.

Actually, the adjacency matrix of the transitive tournament on
{0,...,n-1} (with i→j iff i<j) is the strict upper triangular matrix:
  A_trans[i][j] = 1 if i < j, 0 otherwise.

Eigenvalues of strict upper triangular: all eigenvalues are 0!
(It's nilpotent: A^n = 0.)

So the near-transitive perturbation (rank-2 update) changes the
eigenvalues from all-zero to something nontrivial.

The perturbation: remove arc 0→(n-1), add arc (n-1)→0.
This is a rank-2 update: A_nt = A_trans - e_0 e_{n-1}^T + e_{n-1} e_0^T.
""")

for n in [4, 5, 6, 7, 8]:
    T = build_near_transitive(n)
    A = np.array(T, dtype=float)

    eigenvalues = np.linalg.eigvals(A)
    eigenvalues = sorted(eigenvalues, key=lambda z: -z.real)

    print(f"\nn={n}: Eigenvalues of near-transitive adjacency matrix:")
    for i, e in enumerate(eigenvalues):
        if abs(e.imag) < 1e-10:
            print(f"  λ_{i} = {e.real:.6f}")
        else:
            print(f"  λ_{i} = {e.real:.6f} + {e.imag:.6f}i  (|λ| = {abs(e):.6f})")

    # Characteristic polynomial
    charpoly = np.round(np.poly(A), 6)
    print(f"  Char poly: {list(charpoly)}")

    # Determinant and trace
    det = np.linalg.det(A)
    trace = np.trace(A)
    print(f"  det = {det:.0f}, trace = {trace:.0f}")

    # Permanent (= H) via Ryser formula for small n
    if n <= 7:
        H = sum(1 for perm in permutations(range(n))
                if all(T[perm[i]][perm[(i+1)]] for i in range(n-1)))
        print(f"  H = {H} (= 2^{n-2}+1 = {2**(n-2)+1})")

# ======================================================================
# PART 4: THE CYCLE INDEX — CONNECTING TO k-NACCI
# ======================================================================
print()
print("=" * 70)
print("PART 4: CYCLE INDEX OF NEAR-TRANSITIVE = BINARY EXPANSION")
print("=" * 70)

print("""
The odd cycle count by length for the near-transitive tournament is:
  c_{2k+1} = C(n-2, 2k-1)

This is the ODD-INDEXED row of Pascal's triangle!

    n  | c_3   c_5   c_7   c_9  | total = 2^{n-3}
  -----+------------------------+------------------
    4  |  2                      |  2
    5  |  3    1                 |  4
    6  |  4    4                 |  8
    7  |  5   10    1            |  16
    8  |  6   20    6            |  32
    9  |  7   35   21    1       |  64
   10  |  8   56   56    8       |  128

The generating function: Σ_k c_{2k+1} x^k = Σ_k C(n-2, 2k-1) x^k.

At x=1: sum = 2^{n-3} (half the binomial coefficients of n-2).

The OCF gives: H = I(K_m, 2) = 1 + 2m where m = 2^{n-3}.

The cycle INDEX polynomial is:
  Z(x) = Σ_k C(n-2, 2k-1) x^{2k+1}
       = x · ((1+x)^{n-2} - (1-x)^{n-2}) / 2

This IS the ODD PART of (1+x)^{n-2}!

In terms of the Cayley transform Q(x) = (1+x)/(1-x):
  The odd part of Q(x)^m is related to sinh via:
  (Q^m - Q^{-m})/2 = sinh(2m·arctanh(x)) / (1-x^2)^m

This connects the cycle index to the HYPERBOLIC SINE —
the imaginary part of the helix!

The real part (cosh) gives the EVEN part (transitive triples),
and the imaginary part (sinh) gives the ODD part (cycles).

TOURNAMENT HELICAL DECOMPOSITION:
  even part ↔ cosh ↔ transitive triples ↔ simplicial skeleton
  odd part  ↔ sinh ↔ odd cycles ↔ monodromy loops
""")

# Verify the cycle count formula
for n in range(4, 12):
    counts = {}
    total = 0
    for k in range(1, n):
        length = 2*k+1
        if length > n:
            break
        c = comb(n-2, 2*k-1)
        if c > 0:
            counts[length] = c
            total += c

    print(f"  n={n:2d}: {counts}  total = {total} = 2^{n-3} = {2**(n-3)}")

print("""
THE k-NACCI CONNECTION:
The cycle counts C(n-2, 1), C(n-2, 3), C(n-2, 5), ... are the
odd-indexed binomial coefficients. These satisfy:

  Σ C(n-2, 2k-1) = 2^{n-3}

The RATIO of consecutive terms:
  C(n-2, 2k+1) / C(n-2, 2k-1) = (n-2k-1)(n-2k-2) / ((2k+1)(2k))

As n → ∞ with k fixed: this ratio → (n/2k)² / 4 → ∞.
As k → ∞ with n fixed: ratio → 0.

This is NOT a k-nacci-type recurrence. Instead, the cycle counts
form a UNIMODAL sequence (rising then falling), peaking near k ≈ n/4.

But the TOTAL 2^{n-3} IS the tribonacci connection: it equals
the number of binary strings of length n-3 (no constraints),
which is the k→∞ limit of the k-nacci count!
""")

print("=" * 70)
print("COMPLETE — opus-2026-03-15-S90b")
print("=" * 70)
