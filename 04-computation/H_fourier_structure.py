#!/usr/bin/env python3
"""
H_fourier_structure.py -- kind-pasteur-2026-03-13-S61

DISCOVERY: H(T) has ONLY EVEN-DEGREE Walsh-Hadamard Fourier components!

This was discovered in vitali_tournament_structure.py:
  n=3: degrees {0, 2} only
  n=4: degrees {0, 2} only
  n=5: degrees {0, 2, 4} only

This means H is a POLYNOMIAL of degree 2 in the pairwise products
sigma_{ij} * sigma_{kl} — it's a QUADRATIC FORM on the tournament space!

Wait, actually at n=5 degree 4 appears. So H is not purely quadratic.
But the key is: ALL ODD DEGREES VANISH.

WHY? Because H is invariant under GLOBAL REVERSAL: H(T^op) = H(T)
(the complement tournament has the same number of Hamiltonian paths).
Global reversal flips every sigma_{ij} -> -sigma_{ij}, which multiplies
each degree-k term by (-1)^k. Invariance forces odd-degree terms to zero.

This is THM-030 (transfer matrix symmetry) in the Fourier domain!

Let's verify this at larger n and analyze the spectral structure more deeply.

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import permutations, combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
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


def tournament_to_sigma(A, n):
    """Extract upper triangle sigma vector in {-1,+1}^m."""
    vec = []
    for i in range(n):
        for j in range(i+1, n):
            vec.append(1 if A[i][j] else -1)
    return tuple(vec)


def count_ham_paths(A, n):
    if n <= 1:
        return 1
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


def complement_tournament(bits, n):
    """Flip all edges: T -> T^op."""
    m = n * (n - 1) // 2
    return ((1 << m) - 1) ^ bits


print("=" * 70)
print("H FOURIER STRUCTURE: EVEN DEGREE VANISHING")
print("=" * 70)

# Verify H(T) = H(T^op)
print("\n--- Verification: H(T) = H(T^op) ---")
for n in range(3, 7):
    m = n * (n - 1) // 2
    total = 1 << m
    mismatch = 0
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        A_comp = binary_to_tournament(complement_tournament(bits, n), n)
        if count_ham_paths(A, n) != count_ham_paths(A_comp, n):
            mismatch += 1
    print(f"  n={n}: {mismatch}/{total} mismatches (H(T) != H(T^op))")

print("""
PROOF SKETCH for H(T) = H(T^op):
  H(T) counts Hamiltonian paths a_1 -> a_2 -> ... -> a_n.
  H(T^op) counts paths where all edges are reversed.
  Reading a path backward: a_n -> a_{n-1} -> ... -> a_1 is a path in T^op.
  This gives a BIJECTION: paths in T <-> paths in T^op (reversed).
  Therefore H(T) = H(T^op) for ALL tournaments.
  (This is THM-030: transfer matrix symmetry.)

CONSEQUENCE for Fourier analysis:
  T^op corresponds to sigma -> -sigma (flip all signs).
  The Fourier coefficient hat(H)(S) transforms as:
    hat(H)(S) -> (-1)^{|S|} * hat(H)(S)
  Invariance: hat(H)(S) = (-1)^{|S|} * hat(H)(S)
  => hat(H)(S) = 0 for |S| odd.

This is an EXACT structural result. H is an EVEN function on {-1,+1}^m.
""")

# Detailed Fourier analysis at n=3,4,5
for n in range(3, 6):
    m = n * (n - 1) // 2
    total = 1 << m

    print(f"\n{'='*70}")
    print(f"  FOURIER ANALYSIS AT n={n} (m={m}, |space|={total})")
    print(f"{'='*70}")

    # Compute H for all tournaments
    H_vals = {}
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        sigma = tournament_to_sigma(A, n)
        H = count_ham_paths(A, n)
        H_vals[sigma] = H

    # Compute Fourier coefficients
    fourier = {}
    for k in range(m + 1):
        for S in combinations(range(m), k):
            coeff = 0
            for sigma, H in H_vals.items():
                chi = 1
                for idx in S:
                    chi *= sigma[idx]
                coeff += H * chi
            coeff /= total
            if abs(coeff) > 1e-10:
                fourier[S] = coeff

    # Report by degree
    by_degree = defaultdict(list)
    for S, c in fourier.items():
        by_degree[len(S)].append((S, c))

    total_energy = sum(c**2 for c in fourier.values())

    print(f"\n  Degree | #coeffs | Energy    | % of total")
    print(f"  -------+---------+-----------+-----------")
    for deg in range(m + 1):
        items = by_degree.get(deg, [])
        if items:
            energy = sum(c**2 for _, c in items)
            pct = 100 * energy / total_energy
            print(f"  {deg:>5d}  | {len(items):>7d} | {energy:>9.4f} | {pct:>8.2f}%")

    # Check: are all coefficients of the form k/2^a for integer k?
    # (This would show H is "rational" in a specific sense)
    print(f"\n  Non-zero Fourier coefficients (degree <= 4):")
    edge_labels = []
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            edge_labels.append(f"e_{{{i},{j}}}")
            pos += 1

    for S, c in sorted(fourier.items(), key=lambda x: (len(x[0]), x[0])):
        if len(S) <= 4:
            label = " * ".join(edge_labels[idx] for idx in S) if S else "1"
            # Check if coefficient is a nice fraction
            for denom in [1, 2, 4, 8, 16, 32]:
                if abs(c * denom - round(c * denom)) < 1e-8:
                    num = round(c * denom)
                    frac = f"{num}/{denom}" if denom > 1 else str(num)
                    print(f"    {label:>30s}: hat(H) = {c:.6f} = {frac}")
                    break
            else:
                print(f"    {label:>30s}: hat(H) = {c:.6f}")

    # The QUADRATIC FORM interpretation
    # H(sigma) = c_0 + sum_{i<j} c_{ij} * sigma_i * sigma_j + higher even terms
    # where sigma_i = sigma_{edge i}
    if n <= 5:
        print(f"\n  QUADRATIC FORM STRUCTURE:")
        c0 = fourier.get((), 0)
        print(f"    Constant term: {c0}")

        # Quadratic matrix Q[i,j] = hat(H)(i,j) for i<j
        Q = [[0.0]*m for _ in range(m)]
        for (S, c) in fourier.items():
            if len(S) == 2:
                Q[S[0]][S[1]] = c
                Q[S[1]][S[0]] = c

        # Eigenvalues of Q
        # Simple power method for eigenvalue analysis
        # Actually just compute tr(Q) and tr(Q^2) for now
        trQ = sum(Q[i][i] for i in range(m))
        trQ2 = sum(Q[i][j] * Q[j][i] for i in range(m) for j in range(m))

        print(f"    tr(Q) = {trQ:.4f}")
        print(f"    tr(Q^2) = {trQ2:.4f} => ||Q||_F = {trQ2**0.5:.4f}")
        print(f"    Q energy = {trQ2:.4f}")

        # Fraction of energy in degree 2
        deg2_energy = sum(c**2 for S, c in fourier.items() if len(S) == 2)
        print(f"    Degree-2 fraction: {100*deg2_energy/total_energy:.2f}%")

        # Symmetry analysis of Q
        # For tournaments, the edges are ordered as (0,1),(0,2),(0,3),...,(n-2,n-1)
        # The quadratic form Q encodes which edge-edge products contribute to H

        # Check: which edge pairs have non-zero coefficients?
        print(f"\n    Edge-pair interaction matrix (non-zero entries of Q):")
        non_zero_count = 0
        for i in range(m):
            for j in range(i+1, m):
                if abs(Q[i][j]) > 1e-10:
                    non_zero_count += 1
        print(f"    {non_zero_count} non-zero entries out of {m*(m-1)//2}")

        # Check if Q has a nice structure related to the tournament graph
        # Edge i corresponds to vertex pair (a_i, b_i)
        edge_pairs = []
        for i in range(n):
            for j in range(i+1, n):
                edge_pairs.append((i, j))

        # Classify edge-pair interactions by vertex overlap
        by_overlap = defaultdict(list)
        for i in range(m):
            for j in range(i+1, m):
                if abs(Q[i][j]) < 1e-10:
                    continue
                ov = len(set(edge_pairs[i]) & set(edge_pairs[j]))
                by_overlap[ov].append(((edge_pairs[i], edge_pairs[j]), Q[i][j]))

        print(f"\n    Q entries by vertex overlap of the two edges:")
        for ov in sorted(by_overlap):
            entries = by_overlap[ov]
            values = [v for _, v in entries]
            if len(set(round(v, 6) for v in values)) == 1:
                print(f"      overlap={ov}: {len(entries)} pairs, ALL coefficient = {values[0]:.4f}")
            else:
                print(f"      overlap={ov}: {len(entries)} pairs, coefficients = {sorted(set(round(v,4) for v in values))}")


# ========================================================================
# DEEPER: THE OVERLAP WEIGHT CONNECTION
# ========================================================================
print(f"\n{'='*70}")
print("CONNECTION TO OVERLAP WEIGHT ANALYSIS")
print("=" * 70)

print("""
The Fourier structure reveals WHY overlap weight matters:

1. H(sigma) = sum of even-degree monomials in sigma.
2. Each degree-2 term sigma_i * sigma_j represents a PAIR of edges.
   - If both edges point "forward" (sigma=+1,+1) or both "backward" (-1,-1):
     contribution = +|coeff|
   - If edges disagree: contribution = -|coeff|
3. This is precisely the CONFLICT/ALIGNMENT structure of edges!

For 3-cycles:
  A 3-cycle on {a,b,c} uses edges (a,b), (b,c), (a,c).
  The cycle exists iff sigma_{ab} * sigma_{bc} * sigma_{ac} = -1
  (one must go "against" the canonical ordering).

  But this is a DEGREE-3 (odd) monomial, which has zero Fourier coefficient!
  The cycle count c_3 is NOT a linear function of sigma.
  Instead, c_3 depends on the signs through:
    c_3 = C(n,3) - (sum of degree-2 terms) / normalization

The overlap weight W[i,j] between two 3-cycles C_i, C_j is:
  W = |V(C_i) ∩ V(C_j)|
  If W=0: disjoint, they contribute to alpha_2
  If W>0: conflicting, they share vertex resources

The Fourier degree-2 terms ENCODE the edge-level conflicts that
DETERMINE the cycle-level conflicts in Omega(T).
""")

# Compute the edge-edge interaction and its relation to cycle conflicts
n = 5
m = n * (n - 1) // 2
total = 1 << m

# Recompute Fourier at n=5
H_vals = {}
sigma_to_bits = {}
for bits in range(total):
    A = binary_to_tournament(bits, n)
    sigma = tournament_to_sigma(A, n)
    H = count_ham_paths(A, n)
    H_vals[sigma] = H
    sigma_to_bits[sigma] = bits

# Degree-2 Fourier
Q = {}
for S in combinations(range(m), 2):
    coeff = 0
    for sigma, H in H_vals.items():
        coeff += H * sigma[S[0]] * sigma[S[1]]
    coeff /= total
    if abs(coeff) > 1e-10:
        Q[S] = coeff

# Map edges back to vertex pairs
edge_pairs = []
for i in range(n):
    for j in range(i+1, n):
        edge_pairs.append((i, j))

# Analyze: how do degree-2 Fourier terms relate to cycle conflicts?
# For each tournament, compute c_3 and alpha_2
print(f"\n  n={n}: Relating Fourier structure to cycles")

# Get all 3-cycle vertex sets
all_c3 = set()
for bits in range(total):
    A = binary_to_tournament(bits, n)
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            all_c3.add(frozenset([a, b, c]))
        if A[a][c] and A[c][b] and A[b][a]:
            all_c3.add(frozenset([a, b, c]))
all_c3 = list(all_c3)
print(f"  Total possible 3-cycle vertex sets: {len(all_c3)}")

# For each tournament, which 3-cycles are present?
# and what's the alpha decomposition?
print(f"\n  Tournament | H | c3 | disj_pairs | deg2_energy")
print(f"  ----------+---+----+------------+-----------")

# Group by isomorphism class
from collections import Counter

H_to_stats = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    sigma = tournament_to_sigma(A, n)
    H = count_ham_paths(A, n)

    # Count 3-cycles
    c3 = 0
    present_c3 = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            present_c3.append(frozenset([a, b, c]))
            c3 += 1
        if A[a][c] and A[c][b] and A[b][a]:
            present_c3.append(frozenset([a, b, c]))
            c3 += 1
    present_c3 = list(set(present_c3))

    # Disjoint 3-cycle pairs
    disj = 0
    for i in range(len(present_c3)):
        for j in range(i+1, len(present_c3)):
            if not (present_c3[i] & present_c3[j]):
                disj += 1

    # Degree-2 energy contribution
    deg2_val = sum(c * sigma[S[0]] * sigma[S[1]] for S, c in Q.items())

    H_to_stats[H].append((c3, len(present_c3), disj, deg2_val))

for H in sorted(H_to_stats.keys()):
    stats = H_to_stats[H]
    c3_vals = set(s[0] for s in stats)
    c3v_vals = set(s[1] for s in stats)
    disj_vals = set(s[2] for s in stats)
    deg2_vals = set(round(s[3], 4) for s in stats)
    print(f"  H={H:>3d} ({len(stats):>3d} tours): c3={c3_vals}, c3_sets={c3v_vals}, "
          f"disj={disj_vals}, deg2={deg2_vals}")


# ========================================================================
# THE VITALI-FOURIER CONNECTION
# ========================================================================
print(f"\n{'='*70}")
print("THE VITALI-FOURIER BRIDGE")
print("=" * 70)

print("""
The DEEP connection between Vitali sets and tournament Fourier analysis:

1. VITALI: The quotient R/Q is non-measurable because Q is dense.
   The "characters" of R/Q are exactly the exponential functions
   e^{2 pi i q x} for q in Q-hat. Fourier analysis on R/Q is ill-defined.

2. TOURNAMENTS: The quotient {-1,+1}^m / S_n has orbits of varying size.
   The "characters" of {-1,+1}^m are the Walsh functions chi_S.
   Fourier analysis works perfectly (finite group), BUT:

   The S_n-invariant functions (like H) have a REDUCED Fourier expansion.
   Only the S_n-averaged characters survive:

   tilde(chi_S)(sigma) = (1/n!) sum_{pi in S_n} chi_{pi(S)}(pi(sigma))

3. These averaged characters are the ZONAL FUNCTIONS of the Gelfand pair
   (S_n, S_n acting on {-1,+1}^m).

4. The "hidden higher-dimensional structure" at n=0,1,2 is:
   - n=0,1: dim=0, NO characters exist. H=1 is forced.
   - n=2: dim=1, ONE character (sigma_{01}). But H(T)=H(T^op) kills it.
     The averaged character is: tilde(chi_1) = 0. H=1 is FORCED again.
   - n=3: dim=3, FIRST non-trivial degree-2 characters survive averaging.
     H can distinguish transitive from cyclic.
     The dimension jump 0->3 creates a PHASE TRANSITION in H-complexity.

This is exactly the Vitali phenomenon: at n=2, the "quotient" has only
one class (all tournaments are isomorphic). The "measure" is trivially
consistent. At n=3, TWO classes appear, and the measure question becomes
non-trivial. At n->infinity, the quotient becomes "wild" (Vitali-like).
""")

# Verify: the S_n-averaged Fourier transform
print("S_n-AVERAGED FOURIER CHARACTERS:")
for n in [3, 4]:
    m = n * (n - 1) // 2
    total = 1 << m

    edge_pairs = []
    for i in range(n):
        for j in range(i+1, n):
            edge_pairs.append((i, j))

    # Map permutation pi to its action on edge indices
    def perm_on_edges(perm, n, edge_pairs):
        """Given vertex permutation, compute induced edge permutation."""
        edge_to_idx = {e: i for i, e in enumerate(edge_pairs)}
        result = [0] * len(edge_pairs)
        for idx, (a, b) in enumerate(edge_pairs):
            pa, pb = perm[a], perm[b]
            new_edge = (min(pa, pb), max(pa, pb))
            result[idx] = edge_to_idx[new_edge]
        return result

    # Compute the averaged characters for degree 2
    print(f"\n  n={n}: Degree-2 zonal characters")

    # For each degree-2 subset S of edges, average over S_n
    all_deg2 = list(combinations(range(m), 2))
    # Two subsets are S_n-equivalent if they differ by a vertex relabeling
    # Group them into orbits
    orbit_reps = {}
    visited = set()

    for S in all_deg2:
        if S in visited:
            continue
        orbit = set()
        for perm in permutations(range(n)):
            pi_edges = perm_on_edges(list(perm), n, edge_pairs)
            # Apply to subset
            new_S = tuple(sorted([pi_edges[s] for s in S]))
            orbit.add(new_S)
        for s in orbit:
            visited.add(s)
        orbit_reps[S] = len(orbit)

    print(f"    Number of degree-2 orbits: {len(orbit_reps)}")
    for rep, size in sorted(orbit_reps.items()):
        e1 = edge_pairs[rep[0]]
        e2 = edge_pairs[rep[1]]
        ov = len(set(e1) & set(e2))
        print(f"    {e1}-{e2} (overlap={ov}): orbit size = {size}")


# ========================================================================
# ANALYSIS: H AS EVEN POLYNOMIAL IN SIGMA
# ========================================================================
print(f"\n{'='*70}")
print("H AS EVEN POLYNOMIAL: EXPLICIT FORMULA")
print("=" * 70)

for n in [3, 4]:
    m = n * (n - 1) // 2
    total = 1 << m

    edge_pairs = []
    for i in range(n):
        for j in range(i+1, n):
            edge_pairs.append((i, j))

    H_vals = {}
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        sigma = tournament_to_sigma(A, n)
        H = count_ham_paths(A, n)
        H_vals[sigma] = H

    # Full Fourier expansion (only even degrees)
    fourier = {}
    for k in range(0, m + 1, 2):  # only even degrees
        for S in combinations(range(m), k):
            coeff = 0
            for sigma, H in H_vals.items():
                chi = 1
                for idx in S:
                    chi *= sigma[idx]
                coeff += H * chi
            coeff /= total
            if abs(coeff) > 1e-10:
                fourier[S] = coeff

    print(f"\n  n={n}: H(sigma) = ", end="")
    terms = []
    for S, c in sorted(fourier.items(), key=lambda x: (len(x[0]), x[0])):
        if S == ():
            terms.append(f"{c:.0f}")
        else:
            edge_str = " ".join(f"s_{{{edge_pairs[i][0]}{edge_pairs[i][1]}}}" for i in S)
            sign = "+" if c > 0 else "-"
            terms.append(f"{sign} {abs(c):.0f}*{edge_str}")

    # Print nicely
    formula = " ".join(terms)
    print(formula)

    # Verify formula
    max_err = 0
    for sigma, H in H_vals.items():
        H_computed = 0
        for S, c in fourier.items():
            chi = 1
            for idx in S:
                chi *= sigma[idx]
            H_computed += c * chi
        max_err = max(max_err, abs(H - H_computed))
    print(f"  Max error: {max_err:.10f}")

    # Express in terms of "agreement" variables
    # a_{ij,kl} = sigma_{ij} * sigma_{kl} = +1 if both agree, -1 if disagree
    print(f"\n  H in terms of edge-pair agreements (a_{{ij,kl}} = sigma_ij * sigma_kl):")
    print(f"  H(T) = {fourier.get((), 0):.0f}", end="")
    for S, c in sorted(fourier.items(), key=lambda x: (len(x[0]), x[0])):
        if len(S) == 2:
            e1 = edge_pairs[S[0]]
            e2 = edge_pairs[S[1]]
            sign = "+" if c > 0 else "-"
            # a_{ij,kl} represents "do edges ij and kl agree?"
            print(f" {sign} {abs(c):.0f}*a_{{{e1[0]}{e1[1]},{e2[0]}{e2[1]}}}", end="")
    print()
    if any(len(S) > 2 for S in fourier):
        print(f"  + higher even-degree terms")


print("\n" + "=" * 70)
print("DONE.")
