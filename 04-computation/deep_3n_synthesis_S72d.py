"""
DEEP SYNTHESIS: The 3/n fraction as the bridge between THM-217 and THM-224.

This script systematically explores three interconnected structures:

A. THM-224: C^TC = n·P (simplicial up-Laplacian), independence fraction 3/n
B. THM-217: Transfer matrix M(x), char poly λ³-λ²-xλ-x factoring at x=n(n-3)
C. The C_T^TC_T spectrum for arbitrary tournaments and its eigenvalue anatomy

GOAL: Understand WHY the fraction 3/n appears in both the homological
structure (THM-224) and the HP-counting structure (THM-217), and what
the spectral properties of C_T^TC_T tell us about this connection.
"""
import numpy as np
from math import comb, factorial
from collections import Counter, defaultdict

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def constraint_matrix(adj, n):
    """Returns C_T (transitive triples only) and C_R (cyclic triples)."""
    E = n*(n-1)//2
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i,j))
            edge_idx[(i,j)] = len(edges)-1

    trans_rows = []
    cyclic_rows = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                wins = [adj[a][b]+adj[a][c], adj[b][a]+adj[b][c], adj[c][a]+adj[c][b]]
                row = np.zeros(E)
                # boundary: ∂(a,b,c) = (b,c) - (a,c) + (a,b) for the simplicial complex
                row[edge_idx[(a,c)]] = 1
                row[edge_idx[(a,b)]] = -1
                row[edge_idx[(b,c)]] = -1
                if max(wins) == 2:  # transitive
                    trans_rows.append(row)
                else:  # cyclic (all wins = 1)
                    cyclic_rows.append(row)

    C_T = np.array(trans_rows) if trans_rows else np.zeros((0, E))
    C_R = np.array(cyclic_rows) if cyclic_rows else np.zeros((0, E))
    return C_T, C_R, edges

def scores(adj, n):
    return tuple(sorted([sum(adj[i]) for i in range(n)]))

# ====================================================================
print("=" * 70)
print("SECTION A: EIGENVALUE ANATOMY OF C_T^TC_T AT n=5")
print("=" * 70)
# ====================================================================

n = 5
E = n*(n-1)//2  # 10
max_rank = (n-1)*(n-2)//2  # 6
Cn3 = comb(n, 3)  # 10

# Exhaustive scan
spectra_by_t3 = defaultdict(list)
spectra_by_beta1 = defaultdict(list)

for bits in range(1 << E):
    adj = tournament_from_bits(n, bits)
    C_T, C_R, edges = constraint_matrix(adj, n)
    t3 = C_T.shape[0]

    if t3 == 0:
        continue

    CTC = C_T.T @ C_T
    eigs = np.sort(np.linalg.eigvalsh(CTC))[::-1]
    rk = sum(1 for e in eigs if e > 1e-8)
    beta1 = max_rank - rk

    # Round eigenvalues for classification
    eigs_round = tuple(round(e, 3) for e in eigs)
    spectra_by_t3[t3].append((eigs_round, beta1, bits))
    spectra_by_beta1[beta1].append((t3, eigs_round))

print(f"\nn={n}: Exhaustive scan of {1<<E} tournaments")
print(f"max_rank = {max_rank}, C(n,3) = {Cn3}")

# For each t₃, show the distinct spectra
for t3 in sorted(spectra_by_t3.keys()):
    entries = spectra_by_t3[t3]
    spec_counter = Counter()
    for (eigs, b1, bits) in entries:
        nz = tuple(e for e in eigs if e > 0.01)
        spec_counter[(nz, b1)] += 1

    print(f"\n  t₃={t3} ({len(entries)} tournaments):")
    for (nz, b1), count in spec_counter.most_common():
        nz_str = ", ".join(f"{e:.3f}" for e in nz)
        print(f"    β₁={b1}, nonzero eigs=[{nz_str}] × {count}")
        # Check eigenvalue identities
        eig_sum = sum(nz)
        eig_prod = 1
        for e in nz:
            eig_prod *= e
        print(f"      sum={eig_sum:.3f} (=trace), product={eig_prod:.3f}")

# ====================================================================
print("\n" + "=" * 70)
print("SECTION B: EIGENVALUE PAIRS WITH PRODUCT = n")
print("=" * 70)
# ====================================================================

print(f"\nSearching for eigenvalue pairs (λ,μ) with λ·μ = {n}...")
print(f"Expected from THM-224: all eigenvalues = n for transitive tournament.")
print(f"For general tournaments: eigenvalues split into pairs with product=n?")

pair_data = defaultdict(list)
for t3 in sorted(spectra_by_t3.keys()):
    for (eigs, b1, bits) in spectra_by_t3[t3]:
        nz = [e for e in eigs if e > 0.01]
        pairs = []
        used = set()
        for i in range(len(nz)):
            if i in used:
                continue
            for j in range(i+1, len(nz)):
                if j in used:
                    continue
                if abs(nz[i] * nz[j] - n) < 0.01:
                    pairs.append((round(nz[i],4), round(nz[j],4)))
                    used.add(i)
                    used.add(j)
                    break
            # Eigenvalue equal to n pairs with 1? No, n*1=n
            if i not in used and abs(nz[i] - n) < 0.01:
                pairs.append((n, 'self'))
                used.add(i)
        unpaired = [nz[i] for i in range(len(nz)) if i not in used]
        pair_data[t3].append((pairs, unpaired, b1))

for t3 in sorted(pair_data.keys()):
    entries = pair_data[t3]
    # Show unique patterns
    patterns = Counter()
    for (pairs, unpaired, b1) in entries:
        p_str = str(pairs) + " | unpaired=" + str([round(u,3) for u in unpaired])
        patterns[(p_str, b1)] += 1

    print(f"\n  t₃={t3}:")
    for (p_str, b1), count in patterns.most_common(3):
        print(f"    β₁={b1}: {p_str} × {count}")

# ====================================================================
print("\n" + "=" * 70)
print("SECTION C: THE CONSTRAINT LAPLACIAN DECOMPOSITION")
print("=" * 70)
# ====================================================================

print("""
THM-224 says C^TC = n·P for the FULL (transitive) constraint matrix.
For a general tournament T with t₃ transitive triples:
  C_T^TC_T = n·P - C_R^TC_R

where C_R has the cyclic triple rows. This gives:
  eigenvalues of C_T^TC_T = n - eigenvalues of C_R^TC_R (on im(C^T))

So the spectrum of C_T^TC_T is determined by how the CYCLIC triples
interact with the constraint column space.

KEY: C_R^TC_R has rank ≤ min(|C_R|, dim(im(C^T)∩im(C_R^T)))
""")

# Verify the decomposition C_T^TC_T + C_R^TC_R = n·P
print("Verifying C_T^TC_T + C_R^TC_R = n·P at n=5:")
for bits in [0, 42, 127, 511, 1023]:
    adj = tournament_from_bits(n, bits)
    C_T, C_R, edges = constraint_matrix(adj, n)
    if C_T.shape[0] == 0 or C_R.shape[0] == 0:
        continue
    CTC_T = C_T.T @ C_T
    CTC_R = C_R.T @ C_R
    CTC_sum = CTC_T + CTC_R

    # This should be n·P where P is projection onto constraint space
    eigs_sum = np.sort(np.linalg.eigvalsh(CTC_sum))[::-1]
    nz_sum = [e for e in eigs_sum if e > 0.01]
    print(f"  bits={bits}: t₃={C_T.shape[0]}, c₃={C_R.shape[0]}, "
          f"sum eigenvalues = {[round(e,2) for e in nz_sum[:8]]}")

# ====================================================================
print("\n" + "=" * 70)
print("SECTION D: THE 3/n FRACTION IN EIGENVALUE STATISTICS")
print("=" * 70)
# ====================================================================

print("""
The fraction 3/n = max_rank/C(n,3) appears as:
  1. Independence fraction of the constraint system
  2. The ratio α/(−n) where α = 3-n is the THM-217 third eigenvalue
  3. The "efficiency" of each constraint (each of n edges contributes to
     n-2 triples, but only 3/n of triples are independent)

Question: Does 3/n appear DIRECTLY in the eigenvalue statistics?
""")

# Compute average eigenvalue statistics at each t₃
for t3_val in sorted(spectra_by_t3.keys()):
    entries = spectra_by_t3[t3_val]
    traces = []
    top_eigs = []
    min_nz_eigs = []
    for (eigs, b1, bits) in entries:
        nz = [e for e in eigs if e > 0.01]
        traces.append(sum(nz))
        top_eigs.append(nz[0])
        min_nz_eigs.append(nz[-1])

    avg_trace = np.mean(traces)
    avg_top = np.mean(top_eigs)
    avg_min = np.mean(min_nz_eigs)

    # Independence ratio at this t₃
    indep_ratio = max_rank / t3_val if t3_val > 0 else float('inf')
    # The "effective n" from trace: trace = sum of eigenvalues
    # For transitive: trace = n·rank = n·(n-1)(n-2)/2
    # For general: trace = # of +1/-1 entries in C_T^T C_T diagonal
    # Each edge in a transitive triple contributes 1 to diagonal
    # So trace = 3·t₃ (NO! trace of C_T^TC_T = sum ||row of C_T||² = 3·t₃??)
    # Actually, (C_T^TC_T)_{ee} = # triples containing edge e = degree of e in triple hypergraph
    # For transitive tournament: each edge in n-2 triples → trace = E(n-2) = n(n-1)(n-2)/2
    # For general: trace = 3·t₃ (each triple contributes 3 to the trace)
    trace_from_t3 = 3 * t3_val

    print(f"  t₃={t3_val:2d}: avg trace={avg_trace:.1f}, 3·t₃={trace_from_t3}, "
          f"top={avg_top:.3f}, min_nz={avg_min:.3f}, "
          f"indep_ratio={indep_ratio:.3f}")

print(f"\n  NOTE: trace(C_T^TC_T) = 3·t₃ ALWAYS (each triple contributes 3 to diagonal)")
print(f"  For transitive: 3·C(n,3) = 3·{Cn3} = {3*Cn3} = n·max_rank = {n}·{max_rank} = {n*max_rank}")
print(f"  So: 3·C(n,3)/max_rank = {3*Cn3/max_rank} = n, confirming 3/n identity")

# ====================================================================
print("\n" + "=" * 70)
print("SECTION E: THE TRACE IDENTITY AND 3/n")
print("=" * 70)
# ====================================================================

print("""
FUNDAMENTAL IDENTITY:
  trace(C^TC) = 3 · C(n,3) = n(n-1)(n-2)/2
  rank(C) = (n-1)(n-2)/2

  ⟹ trace/rank = n     (the uniform eigenvalue!)

  And: rank/C(n,3) = 3/n (the independence fraction)

  So: trace = rank · n = C(n,3) · 3

This triple identity is the HEART of the 3/n connection:

  trace = 3 · C(n,3)         [each triple has 3 edges]
  rank = C(n,3) · (3/n)      [only 3/n of triples are independent]
  eigenvalue = n              [trace/rank = n]

The 3/n fraction is the ratio that makes the uniform eigenvalue n work!

For a GENERAL tournament T with t₃ transitive triples:
  trace(C_T^TC_T) = 3·t₃
  rank(C_T) = min(t₃, max_rank) - β₁

  "effective eigenvalue" = 3·t₃ / rank(C_T)

  For transitive: 3·C(n,3) / max_rank = n
  For regular tournament at n=5: 3·5 / 5 = 3 (not n!)
""")

for t3_val in sorted(spectra_by_t3.keys()):
    entries = spectra_by_t3[t3_val]
    for (eigs, b1, bits) in entries[:1]:  # just first example
        nz = [e for e in eigs if e > 0.01]
        rk = len(nz)
        tr = sum(nz)
        eff_eig = tr / rk if rk > 0 else 0
        print(f"  t₃={t3_val}, β₁={b1}: rank={rk}, trace={tr:.1f}, "
              f"trace/rank={eff_eig:.3f}, 3·t₃/rank={3*t3_val/rk:.3f}")

# ====================================================================
print("\n" + "=" * 70)
print("SECTION F: THM-217 EIGENVALUES IN THE CONSTRAINT SPECTRUM")
print("=" * 70)
# ====================================================================

print("""
Does the metallic mean M_n - 1 = (n-2+√(n²+4))/2 appear in C_T^TC_T?

At n=5: M₅-1 = (3+√29)/2 ≈ 4.193
At n=6: M₆-1 = (4+√40)/2 = 2+√10 ≈ 5.162
""")

for n_val in [5, 6]:
    metallic = (n_val - 2 + np.sqrt(n_val**2 + 4)) / 2
    print(f"\nn={n_val}: M_{n_val}-1 = {metallic:.6f}")
    print(f"  (n±√n)/2 = {(n_val+np.sqrt(n_val))/2:.6f}, {(n_val-np.sqrt(n_val))/2:.6f}")
    print(f"  n = {n_val}")
    print(f"  These are DIFFERENT numbers:")
    print(f"    metallic M_n-1 = {metallic:.6f} (from transfer matrix)")
    print(f"    (n+√n)/2 = {(n_val+np.sqrt(n_val))/2:.6f} (from C_T^TC_T regular tournament)")
    print(f"    n = {n_val} (from C^TC = n·P)")

# ====================================================================
print("\n" + "=" * 70)
print("SECTION G: THE DEEPER BRIDGE — CV² AT x = n(n-3)")
print("=" * 70)
# ====================================================================

print("""
THM-217 says CV²(H) = Σ_k 2·g_k(n-2k) / (n)_{2k}

The transfer matrix M(x) generates g_k via [x^k] Tr(M^{n-2k}).

At x = n(n-3), the char poly factors as (λ-(3-n))(λ²-(n-2)λ-n).

What does evaluating g_k at x = n(n-3) mean?

The generating parameter x is NOT the same as the constraint count.
In THM-217, x tracks domino matching weight.
The value x = n(n-3) is the point where the char poly factors via
the metallic mean.

But n(n-3) = n² · (1-3/n) = n² · redundancy_fraction.

MEANING: x = n(n-3) is the "redundancy-weighted" version of n².
If we think of n² as the "total constraint capacity", then x captures
the portion attributable to redundant constraints.

The THM-224 eigenvalue n appears as:
  trace/rank = 3·C(n,3)/max_rank = n

The THM-217 factorization parameter x = n(n-3) = n² - 3n appears as:
  n² · (C(n,3) - max_rank)/C(n,3) = n² · (n-3)/n = n(n-3) ✓

Both encode the SAME underlying ratio 3/n, but in different
algebraic frameworks.
""")

# The connection: α = 3-n is the eigenvalue that "remembers" 3/n
# in the transfer matrix world.
# In the constraint world, 3/n = max_rank/C(n,3).
# In both cases, 3 = the size of a triple (simplex dimension + 1).

print("The number 3 in '3/n' has multiple meanings:")
print("  1. Each constraint involves 3 edges (simplex dimension)")
print("  2. Each edge appears in n-2 constraints → trace = 3·C(n,3)")
print("  3. THM-217 uses a 3×3 transfer matrix")
print("  4. The metallic quadratic is DEGREE 2 (one less than 3)")
print()

# ====================================================================
print("=" * 70)
print("SECTION H: UNIVERSAL TOP EIGENVALUE — PROOF COMPLETION")
print("=" * 70)
# ====================================================================

print("""
THEOREM: For ANY tournament T on n vertices with t₃ ≥ 1,
the top eigenvalue of C_T^TC_T is exactly n.

PROOF (building on THM-224):

1. C_T^TC_T = n·P - C_R^TC_R where P = (1/n)·C^TC (THM-224).

2. Need: ∃ v ∈ im(P) with C_R·v = 0.
   i.e., v ∈ im(C^T) ∩ ker(C_R).

3. dim(im(C^T)) = (n-1)(n-2)/2
   dim(ker(C_R)) ≥ E - |C_R| = n(n-1)/2 - (C(n,3)-t₃)

   Dimension argument: intersection nonempty if
   (n-1)(n-2)/2 + n(n-1)/2 - (C(n,3)-t₃) > E = n(n-1)/2
   ⟺ (n-1)(n-2)/2 > C(n,3) - t₃
   ⟺ t₃ > C(n,3) - (n-1)(n-2)/2 = C(n-1,3)

4. For n ≤ 5: every tournament has t₃ ≥ C(n-1,3)+1 (checkable).
   For n ≥ 6: need STRONGER argument.

   Actually, the minimum t₃ for tournaments on n vertices is
   C(n,3) - n(n-1)(n-3)/24 (for odd n, regular tournament).
   This is n(n-1)(n-5)/24, which is < C(n-1,3) for n ≤ ???
""")

# Check: min t₃ vs C(n-1,3)
for n_val in range(3, 12):
    Cn3 = comb(n_val, 3)
    max_3cycles = n_val * (n_val-1) * (n_val-3) // 24 if n_val % 2 == 1 else None
    min_t3 = Cn3 - max_3cycles if max_3cycles is not None else "N/A (even n)"
    threshold = comb(n_val-1, 3)
    if isinstance(min_t3, int):
        print(f"  n={n_val}: min_t₃={min_t3}, C(n-1,3)={threshold}, "
              f"min_t₃ > C(n-1,3)? {min_t3 > threshold}")
    else:
        print(f"  n={n_val}: min_t₃=? (even n), C(n-1,3)={threshold}")

print("""
For odd n:
  min_t₃ = n(n-1)(n-5)/24
  C(n-1,3) = (n-1)(n-2)(n-3)/6

  Need: n(n-1)(n-5)/24 > (n-1)(n-2)(n-3)/6
  ⟺ n(n-5)/4 > (n-2)(n-3)
  ⟺ n²-5n > 4(n²-5n+6) = 4n²-20n+24
  ⟺ -3n² + 15n - 24 > 0
  ⟺ 3n² - 15n + 24 < 0

  Discriminant: 225 - 288 = -63 < 0. NEVER true.

  So for odd n ≥ 7: min_t₃ < C(n-1,3). The dimension argument FAILS.
  Yet the top eigenvalue is still n. WHY?
""")

# The answer must be more subtle than dimension counting.
# Let me look at the SPECIFIC structure of C_R.

print("ALTERNATIVE PROOF APPROACH:")
print("""
The rows of C (full constraint matrix) are boundary vectors ∂₂(face).
The rows of C_R are boundaries of CYCLIC triples.

Key observation: Each cyclic triple {a,b,c} has boundary
  ∂₂(a,b,c)_cyclic = e_{bc} - e_{ac} + e_{ab}
which is the SAME vector (up to sign) as the transitive boundary!

Actually no — in the simplicial complex, the boundary only depends on
the UNDIRECTED triple, not orientation. So C has one row per {a,b,c},
and C_T/C_R partition these rows based on tournament structure.

The key is: im(C_T^T) + im(C_R^T) = im(C^T).

And C^TC = n·P means n·P = C_T^TC_T + C_R^TC_R.

For v to satisfy C_T^TC_T·v = n·v, we need C_R^TC_R·v = 0,
i.e., C_R·v = 0. But v must also be in im(P) = im(C^T).

CLAIM: The gradient subspace ker(P) = im(∂₁^T) is contained in
ker(C_R) (because ∂₂·∂₁^T = 0 regardless of which triples we use).

So ker(C_R) ⊇ ker(P), and we need v ∈ im(P) ∩ ker(C_R).

This is the INTERSECTION of im(P) (dimension (n-1)(n-2)/2) with
ker(C_R) (dimension ≥ E - rank(C_R)).

The question is whether C_R has ANY row that is NOT in the span of C_T.
If all rows of C_R are already in im(C_T^T), then ker(C_R) ⊇ ker(C_T)
and the intersection might be empty.
""")

# Check: does im(C_R^T) ⊆ im(C_T^T)?
print("\nChecking if im(C_R^T) ⊆ im(C_T^T) at n=5:")
n = 5
for bits in range(1 << 10):
    adj = tournament_from_bits(n, bits)
    C_T, C_R, edges = constraint_matrix(adj, n)
    if C_T.shape[0] == 0 or C_R.shape[0] == 0:
        continue

    # Check if each row of C_R is in the row space of C_T
    # Equivalently: rank([C_T; C_R]) == rank(C_T)?
    C_full = np.vstack([C_T, C_R])
    rank_T = np.linalg.matrix_rank(C_T)
    rank_full = np.linalg.matrix_rank(C_full)

    if rank_T != rank_full:
        # C_R adds NEW directions
        t3 = C_T.shape[0]
        c3 = C_R.shape[0]
        print(f"  bits={bits}: t₃={t3}, c₃={c3}, rank(C_T)={rank_T}, rank(C_full)={rank_full}")

        # Now check: is there still a v ∈ im(P) ∩ ker(C_R)?
        # im(P) = column space of C^T (= row space of C_full)
        # We need v in this space with C_R·v = 0

        # Build projection P from C_full
        U, S, Vt = np.linalg.svd(C_full, full_matrices=False)
        nz = sum(1 for s in S if s > 1e-8)
        P_basis = Vt[:nz]  # rows are basis of im(C^T) = col space of C_full.T

        # Project C_R onto im(P): C_R restricted to im(P)
        # Represent C_R in the P_basis
        C_R_proj = C_R @ P_basis.T  # shape: (c3, nz)

        # kernel of C_R_proj: vectors in im(P) with C_R·v = 0
        # dim = nz - rank(C_R_proj)
        rank_CR_proj = np.linalg.matrix_rank(C_R_proj)
        dim_intersection = nz - rank_CR_proj

        if dim_intersection > 0:
            pass  # good, v exists
        else:
            print(f"    *** NO v in im(P) ∩ ker(C_R)! dim_intersection = {dim_intersection}")
        break

# Actually, let me just verify exhaustively at n=5 that top eig = n ALWAYS
print("\nExhaustive verification at n=5: top eigenvalue of C_T^TC_T = 5?")
n = 5
E = 10
violations = 0
for bits in range(1 << E):
    adj = tournament_from_bits(n, bits)
    C_T, C_R, edges = constraint_matrix(adj, n)
    if C_T.shape[0] == 0:
        continue
    CTC = C_T.T @ C_T
    top = np.max(np.linalg.eigvalsh(CTC))
    if abs(top - n) > 0.01:
        violations += 1
        print(f"  VIOLATION at bits={bits}: top eig = {top:.6f}")

print(f"  Violations: {violations} / {1<<E}")

# ====================================================================
print("\n" + "=" * 70)
print("SECTION I: WHAT DETERMINES THE NON-TOP EIGENVALUES?")
print("=" * 70)
# ====================================================================

print("""
Since the top eigenvalue is ALWAYS n, the interesting spectral information
is in the REMAINING eigenvalues. For C_T^TC_T = n·P - C_R^TC_R:

  eigenvalue λ_i of C_T^TC_T = n - μ_i where μ_i = eigenvalue of C_R^TC_R
  (restricted to im(P))

So the non-top eigenvalues are n minus the CYCLIC-triple Gram matrix eigenvalues.

For β₁ = 0 (all (n-1)(n-2)/2 eigenvalues nonzero):
  ALL eigenvalues of C_R^TC_R on im(P) are < n.

For β₁ = 1 (one eigenvalue drops to 0):
  C_R^TC_R has eigenvalue n on im(P) → one constraint direction is
  COMPLETELY cyclic.
""")

# At n=5, study the C_R^TC_R spectrum
print("C_R^TC_R eigenvalues at n=5 by t₃:")
n = 5
E = 10

cr_spectra = defaultdict(list)
for bits in range(1 << E):
    adj = tournament_from_bits(n, bits)
    C_T, C_R, edges = constraint_matrix(adj, n)
    if C_R.shape[0] == 0:
        continue

    CTC_R = C_R.T @ C_R
    eigs_R = np.sort(np.linalg.eigvalsh(CTC_R))[::-1]

    # Restrict to im(P) = im(C^T) where C = full constraint matrix
    C_full = np.vstack([C_T, C_R]) if C_T.shape[0] > 0 else C_R
    U, S, Vt = np.linalg.svd(C_full, full_matrices=False)
    nz = sum(1 for s in S if s > 1e-8)
    P_basis = Vt[:nz]

    CTC_R_restricted = P_basis @ CTC_R @ P_basis.T
    eigs_R_restricted = np.sort(np.linalg.eigvalsh(CTC_R_restricted))[::-1]

    t3 = C_T.shape[0] if C_T.shape[0] > 0 else 0
    c3 = C_R.shape[0]

    rk_T = np.linalg.matrix_rank(C_T) if C_T.shape[0] > 0 else 0
    beta1 = max_rank - rk_T

    eigs_round = tuple(round(e, 3) for e in eigs_R_restricted if abs(e) > 0.01)
    cr_spectra[(t3, beta1)].append(eigs_round)

for (t3, b1) in sorted(cr_spectra.keys()):
    entries = cr_spectra[(t3, b1)]
    c = Counter(entries)
    print(f"\n  t₃={t3}, β₁={b1} ({len(entries)} tournaments):")
    for spec, count in c.most_common(3):
        nz_str = ", ".join(f"{e:.3f}" for e in spec)
        print(f"    C_R^TC_R|_{'{im(P)}'} nonzero eigs = [{nz_str}] × {count}")
        # Check: does any eigenvalue = n?
        if any(abs(e - n) < 0.01 for e in spec):
            print(f"      ^^^ Contains eigenvalue {n} → accounts for β₁=1")
        # Check: eigenvalues are n - (C_T^TC_T eigenvalues)?
        ct_eigs = tuple(round(n - e, 3) for e in spec)
        print(f"      → C_T^TC_T eigs = n - these = [{', '.join(f'{e:.3f}' for e in ct_eigs)}]")
