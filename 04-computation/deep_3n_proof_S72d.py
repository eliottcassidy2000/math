"""
PART 2: Structural proof that top eigenvalue = n ALWAYS,
and the eigenvalue arithmetic connecting THM-217 and THM-224.

Key insight from Part 1:
  C_T^TC_T = n·P - C_R^TC_R
  eigenvalues of C_T^TC_T on im(P) = n - eigenvalues of C_R^TC_R on im(P)
  β₁=1 ⟺ C_R^TC_R has eigenvalue n on im(P)

Need to prove: C_R^TC_R NEVER has eigenvalue > n on im(P), and
  always has eigenvalue 0 on im(P) (so C_T^TC_T always has eigenvalue n).
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

def full_constraint_matrix(n):
    """Full simplicial boundary matrix ∂₂ for K_n."""
    E = n*(n-1)//2
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i,j))
            edge_idx[(i,j)] = len(edges)-1

    rows = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                row = np.zeros(E)
                row[edge_idx[(b,c)]] = 1
                row[edge_idx[(a,c)]] = -1
                row[edge_idx[(a,b)]] = 1
                rows.append(row)

    return np.array(rows), edges, edge_idx

def constraint_matrices(adj, n, C_full, edge_idx):
    """Split full constraint matrix into transitive and cyclic parts."""
    Cn3 = comb(n, 3)
    row = 0
    trans_idx = []
    cyclic_idx = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                wins = [adj[a][b]+adj[a][c], adj[b][a]+adj[b][c], adj[c][a]+adj[c][b]]
                if max(wins) == 2:
                    trans_idx.append(row)
                else:
                    cyclic_idx.append(row)
                row += 1

    C_T = C_full[trans_idx] if trans_idx else np.zeros((0, C_full.shape[1]))
    C_R = C_full[cyclic_idx] if cyclic_idx else np.zeros((0, C_full.shape[1]))
    return C_T, C_R

# ====================================================================
print("=" * 70)
print("PROOF THAT TOP EIGENVALUE = n: THE POTENTIAL VECTOR ARGUMENT")
print("=" * 70)
# ====================================================================

print("""
CLAIM: For any tournament T with at least one transitive triple,
  max eigenvalue of C_T^TC_T = n.

PROOF STRATEGY: Construct an explicit eigenvector.

The eigenvector for eigenvalue n of C^TC (full matrix) lives in im(C^T) = im(∂₂^T).
This is the space of "co-exact 1-forms" — edge functions that arise from
assigning values to triangles.

For the TRANSITIVE tournament: the potential function f(v) = rank(v) gives
  h(i,j) = f(j) - f(i) = rank(j) - rank(i)
This is a gradient (1-coboundary), which is in ker(∂₂^T ∂₂) = ker(P).
NOT useful for eigenvalue n.

Instead, we need v ∈ im(∂₂^T) with C_R·v = 0.
Meaning: v is orthogonal to all CYCLIC triple boundaries.

KEY: The space ker(C_R) ∩ im(C^T) consists of edge functions that:
  1. Can be written as linear combination of triangle "charges"
  2. Have zero "charge" on every cyclic triple

This means: assign a charge q_{abc} to each TRANSITIVE triple {a,b,c},
and compute v = Σ q_{abc} · ∂₂(a,b,c).
Then C_R·v = 0 automatically (since v involves no cyclic triples).

But wait: v = C_T^T · q where q = (q_{abc}). So v ∈ im(C_T^T).

C_T^TC_T · v = C_T^T (C_T C_T^T q) ... this doesn't immediately give n·v.

Let me think differently.
""")

# Actually, we need v ∈ im(P) (the constraint column space for the FULL matrix)
# such that C_R v = 0.
#
# v ∈ im(C^T) means v = C^T w for some w.
# C_R v = C_R C^T w. We need C_R C^T w = 0.
#
# C_R C^T = (C_R C_T^T | C_R C_R^T) composed.
# Actually C^T = [C_T^T | C_R^T], so C_R C^T = C_R [C_T^T | C_R^T] ... no.
# C^T has columns indexed by triples. C_R C^T = C_R · [columns of C^T].
# Actually C^T is E × C(n,3), and C_R is c₃ × E.
# C_R C^T is c₃ × C(n,3). This is the "cyclic-full" cross Gram.
#
# Hmm, let me think of it in terms of the simplicial structure.

# The STRUCTURAL argument:
# C^T has columns that are boundary vectors of ALL triples.
# C_R selects cyclic triple rows.
# C_R · C^T has entry (cyclic triple i, triple j) = <row_i(C_R), col_j(C^T)>
# = <∂₂(triple_i), ∂₂(triple_j)> = inner product of boundary vectors.
#
# Two triples share 0 or 2 edges. If they share 0 edges, inner product = 0.
# If they share exactly 1 edge, inner product = ±1.
# If they share 2 edges (same triple), inner product = 3.

# The key insight might be simpler: let v be in im(C_T^T) directly.
# Then C_R v = C_R C_T^T q for some q ∈ ℝ^{t₃}.
# We need ker(C_R C_T^T) to be nonempty.

n = 5
E = n*(n-1)//2
C_full, edges, edge_idx = full_constraint_matrix(n)

print("Checking the cross-Gram C_R C_T^T at n=5:")
print(f"C_full shape: {C_full.shape}")

# For a specific tournament, check rank of C_R C_T^T
for bits in [0b0011101010]:  # just one example
    adj = tournament_from_bits(n, bits)
    C_T, C_R = constraint_matrices(adj, n, C_full, edge_idx)
    if C_T.shape[0] == 0 or C_R.shape[0] == 0:
        continue

    cross = C_R @ C_T.T
    print(f"\n  bits={bits}: t₃={C_T.shape[0]}, c₃={C_R.shape[0]}")
    print(f"  C_R C_T^T shape: {cross.shape}")
    print(f"  C_R C_T^T rank: {np.linalg.matrix_rank(cross)}")
    print(f"  C_R C_T^T:\n{cross}")

    # Also compute C_R C_R^T (the cyclic-cyclic Gram)
    gram_RR = C_R @ C_R.T
    print(f"\n  C_R C_R^T (cyclic self-Gram):\n{gram_RR}")
    print(f"  Eigenvalues of C_R C_R^T: {np.sort(np.linalg.eigvalsh(gram_RR))[::-1].round(4)}")

# ====================================================================
print("\n\n" + "=" * 70)
print("THE SIMPLICIAL INNER PRODUCT STRUCTURE")
print("=" * 70)
# ====================================================================

print("""
For two triples σ₁ = {a,b,c} and σ₂ = {d,e,f}:
  <∂₂(σ₁), ∂₂(σ₂)> depends only on |σ₁ ∩ σ₂|:

  |intersection| = 3 (same triple): inner product = 3
  |intersection| = 2 (share an edge): inner product = ±1
  |intersection| = 1 (share a vertex only): inner product = 0
  |intersection| = 0 (disjoint): inner product = 0

So C_full · C_full^T is a {0, ±1, 3}-valued matrix.
Its block structure (T×T, T×R, R×T, R×R) determines everything.
""")

# Compute the full inner product matrix at n=5
G = C_full @ C_full.T
print(f"Full Gram matrix C·C^T at n=5: {G.shape}")
print(f"Diagonal entries: {np.diag(G)}")
print(f"Off-diagonal values: {sorted(set(G.flatten()) - {3})}")

# The off-diagonal = ±1 where triples share exactly 2 vertices (an edge)
# How many triples share an edge? Each edge is in n-2 triples.
# So each triple shares an edge with 3*(n-3) other triples.
print(f"Each triple shares edge with 3*(n-3) = {3*(n-3)} others")
print(f"Row sum of |off-diag|: should be {3*(n-3)}")
print(f"Actual: {sum(abs(G[0,j]) for j in range(G.shape[0]) if j != 0)}")

# The eigenvalues of C·C^T:
eigs_G = np.sort(np.linalg.eigvalsh(G))[::-1]
print(f"\nEigenvalues of C·C^T (full, n=5): {eigs_G.round(4)}")
# This should have eigenvalue n with multiplicity (n-1)(n-2)/2 = 6
# and eigenvalue 0 with multiplicity C(n,3) - (n-1)(n-2)/2 = 4
nz = sum(1 for e in eigs_G if e > 0.01)
print(f"Nonzero count: {nz}, expected (n-1)(n-2)/2 = {(n-1)*(n-2)//2}")
print(f"Nonzero eigenvalues: {eigs_G[:nz].round(4)}")

# ====================================================================
print("\n\n" + "=" * 70)
print("KEY DISCOVERY: THE n-EIGENSPACE OF C_R^TC_R")
print("=" * 70)
# ====================================================================

print("""
β₁ = 1 ⟺ C_R^TC_R has eigenvalue n on im(P).

Since C_R^TC_R|_{im(P)} has eigenvalues μ_i, and C_T^TC_T|_{im(P)} has
eigenvalues n - μ_i, we have β₁ = 1 iff some μ_i = n.

For μ_i = n: there exists v ∈ im(P) with ||C_R v||² = n||v||².
Since C^TC = n·P, this means ||C_R v||² = ||C v||² (on im(P)),
which means ||C_T v||² = ||C v||² - ||C_R v||² = 0.
i.e., C_T v = 0.

So: β₁ = 1 ⟺ ∃ v ∈ im(P) ∩ ker(C_T) ⟺ im(P) ⊄ im(C_T^T)^⊥
Actually: v ∈ im(P) ∩ ker(C_T) means v is in the constraint column
space but is killed by all transitive triple constraints.

This is a NON-TRIVIAL condition: it says the homological "cycle"
that defines β₁ is visible in the simplicial structure.
""")

# For each β₁=1 tournament at n=5, find the eigenvector
print("Finding the β₁=1 eigenvector at n=5:")
n = 5
E = n*(n-1)//2
C_full, edges, edge_idx = full_constraint_matrix(n)

# The projection P onto im(C^T)
U, S, Vt = np.linalg.svd(C_full, full_matrices=False)
nz = sum(1 for s in S if s > 1e-8)
P_basis = Vt[:nz]  # (6, 10) — 6 basis vectors in R^10

count = 0
v_examples = []
for bits in range(1 << E):
    adj = tournament_from_bits(n, bits)
    C_T, C_R = constraint_matrices(adj, n, C_full, edge_idx)
    t3 = C_T.shape[0]

    rk = np.linalg.matrix_rank(C_T)
    if rk < (n-1)*(n-2)//2:  # β₁ = 1
        # Find v ∈ im(P) ∩ ker(C_T)
        # Project ker(C_T) onto im(P)
        if C_T.shape[0] > 0:
            # Null space of C_T
            _, _, Vt_T = np.linalg.svd(C_T, full_matrices=True)
            null_start = np.linalg.matrix_rank(C_T)
            null_T = Vt_T[null_start:]  # null space of C_T

            # Intersect with im(P): project null vectors onto P
            null_in_P = null_T @ P_basis.T  # project
            # Nonzero rows of null_in_P are null vectors that have component in im(P)
            rank_null_in_P = np.linalg.matrix_rank(null_in_P)

            if rank_null_in_P > 0 and count < 3:
                # Find actual vector
                U_n, S_n, Vt_n = np.linalg.svd(null_in_P)
                # The vector in im(P) that is also in ker(C_T)
                v_in_P = Vt_n[0] @ P_basis  # back to R^E
                v_in_P /= np.linalg.norm(v_in_P)

                # Verify
                check_CT = np.linalg.norm(C_T @ v_in_P)
                check_P = np.linalg.norm(v_in_P - P_basis.T @ (P_basis @ v_in_P))

                scores_T = tuple(sorted([sum(adj[i]) for i in range(n)]))
                print(f"\n  Tournament bits={bits}, scores={scores_T}, t₃={t3}:")
                print(f"    ||C_T v|| = {check_CT:.6f} (should be ~0)")
                print(f"    ||v - Pv|| = {check_P:.6f} (should be ~0)")
                print(f"    v = {v_in_P.round(4)}")

                # What does this vector look like on the edges?
                for i, (a,b) in enumerate(edges):
                    if abs(v_in_P[i]) > 0.01:
                        print(f"      edge ({a},{b}): {v_in_P[i]:.4f}")

                v_examples.append((bits, v_in_P, t3, scores_T))
                count += 1

# ====================================================================
print("\n\n" + "=" * 70)
print("THE SPECTRAL ARITHMETIC: EIGENVALUE PATTERNS")
print("=" * 70)
# ====================================================================

print("""
From Section A (n=5 exhaustive), the eigenvalue spectra have EXACT algebraic values:

t₃=5 (regular): [5, φn, φn, φ̄n, φ̄n] where φn=(5+√5)/2, φ̄n=(5-√5)/2
t₃=6 (β₁=0):   [5, 5, φn, 1+φ̄n, φ̄n, 1-φn] = [5,5,3.618,2.618,1.382,0.382]
t₃=6 (β₁=1,A): [5, 5, 4, 2, 2]
t₃=6 (β₁=1,B): [5, 5, 3, 3, 2]
t₃=7 (β₁=0):   [5, 5, 5, 2+√2, 2, 2-√2]
t₃=7 (β₁=1):   [5, 5, 5, 3, 3]
t₃=8:           [5, 5, 5, 5, 3, 1]
t₃=9:           [5, 5, 5, 5, 5, 2]
t₃=10:          [5, 5, 5, 5, 5, 5]

OBSERVATION 1: As t₃ increases, eigenvalues "collapse" toward n=5.
  The C_R^TC_R eigenvalues shrink to 0 as fewer cyclic triples remain.

OBSERVATION 2: The golden ratio pair (φ,φ̄) at t₃=5 is:
  φ = (5+√5)/2, φ̄ = (5-√5)/2
  These solve λ² - 5λ + 5 = 0, i.e., λ(λ-5) = -5.

  At t₃=6 (β₁=0), we get φ, φ̄ PLUS shifted versions:
  1+φ̄ = (7-√5)/2 = 2.618,  1-φ+5 = (3+√5)/2 ...
  Actually: 2.618 = 1 + 1.618 = 1 + (5-√5)/2? No, 1+(5-√5)/2 = (7-√5)/2 ≈ 2.382.
  Let me recompute.
""")

# Recompute exact eigenvalues
phi = (5 + np.sqrt(5)) / 2  # 3.618
phibar = (5 - np.sqrt(5)) / 2  # 1.382

print(f"φ = (5+√5)/2 = {phi:.6f}")
print(f"φ̄ = (5-√5)/2 = {phibar:.6f}")
print(f"φ + φ̄ = {phi+phibar:.6f} = n")
print(f"φ · φ̄ = {phi*phibar:.6f} = n")

# At t₃=6 β₁=0: eigenvalues are [5, 5, 3.618, 2.618, 1.382, 0.382]
# The non-5 eigenvalues: 3.618, 2.618, 1.382, 0.382
# = φ, φ+1-φ̄? No. φ = 3.618, 2.618 = φ-1, 1.382 = φ̄, 0.382 = φ̄-1
# Check: 3.618 - 1 = 2.618 ✓, 1.382 - 1 = 0.382 ✓
# So the eigenvalues are φ, φ-1, φ̄, φ̄-1!
# Sum: φ+(φ-1)+φ̄+(φ̄-1) = 2(φ+φ̄)-2 = 2n-2 = 8. With two 5's: 8+10=18=3*6 ✓

print(f"\nt₃=6, β₁=0: non-n eigenvalues are φ, φ-1, φ̄, φ̄-1:")
print(f"  {phi:.3f}, {phi-1:.3f}, {phibar:.3f}, {phibar-1:.3f}")
print(f"  Sum = {phi + (phi-1) + phibar + (phibar-1):.3f} = 2n-2 = {2*n-2}")

# At t₃=7, β₁=0: [5,5,5, 3.414, 2, 0.586]
# 3.414 ≈ 2+√2, 0.586 ≈ 2-√2
# Check: 2+√2 = 3.414, 2-√2 = 0.586
print(f"\nt₃=7, β₁=0: non-n eigenvalues: 2+√2={2+np.sqrt(2):.3f}, 2, 2-√2={2-np.sqrt(2):.3f}")
print(f"  Sum = {(2+np.sqrt(2))+2+(2-np.sqrt(2)):.3f} = 6 = 3*7-15 = 21-15")
print(f"  Products: (2+√2)(2-√2) = {(2+np.sqrt(2))*(2-np.sqrt(2)):.3f} = 2")

# At t₃=8: [5,5,5,5,3,1]. Non-n: 3,1. Sum=4. Product=3.
# At t₃=9: [5,5,5,5,5,2]. Non-n: 2.
# At t₃=10: all 5.

print(f"\n--- PATTERN: Non-n eigenvalue sums ---")
for t3_val, non_n in [(5, [phi,phi,phibar,phibar]),
                       (6, [phi, phi-1, phibar, phibar-1]),
                       (7, [2+np.sqrt(2), 2, 2-np.sqrt(2)]),
                       (8, [3, 1]),
                       (9, [2]),
                       (10, [])]:
    s = sum(non_n)
    p = 1
    for x in non_n:
        p *= x
    n_count = 6 - len(non_n)  # how many eigenvalues = n
    total = 5*n_count + s
    print(f"  t₃={t3_val}: {len(non_n)} non-n eigenvalues, sum={s:.3f}, "
          f"5×{n_count}+{s:.3f}={total:.3f} = 3·t₃={3*t3_val}")

# KEY: sum of non-n eigenvalues = 3·t₃ - n·(6 - len(non_n))
# At t₃=5: non-n sum = 15-5=10 (φ+φ+φ̄+φ̄ = 2(φ+φ̄)=2·5=10 ✓)
# At t₃=6: non-n sum = 18-10=8 (check: 3.618+2.618+1.382+0.382=8 ✓)
# At t₃=7: non-n sum = 21-15=6 (3.414+2+0.586=6 ✓)
# At t₃=8: non-n sum = 24-20=4 (3+1=4 ✓)
# At t₃=9: non-n sum = 27-25=2 ✓
# At t₃=10: non-n sum = 30-30=0 ✓

print(f"\nNon-n sum = 3·t₃ - n·(max_rank - #non_n_eigs)")

# ====================================================================
print("\n\n" + "=" * 70)
print("THE BRIDGE: C_R^TC_R EIGENVALUES AND THE THM-217 CONNECTION")
print("=" * 70)
# ====================================================================

print("""
The C_R^TC_R eigenvalues on im(P) are μ_i = n - λ_i where λ_i are
the C_T^TC_T eigenvalues.

For the REGULAR tournament at n=5 (t₃=5):
  C_T^TC_T eigenvalues: [5, φ, φ, φ̄, φ̄, 0]
  C_R^TC_R eigenvalues: [0, n-φ, n-φ, n-φ̄, n-φ̄, n]
                       = [0, φ̄, φ̄, φ, φ, 5]

REMARKABLE: C_R^TC_R has the SAME spectrum as C_T^TC_T!
This is because t₃ = c₃ = 5 (symmetric partition).

Actually this makes sense: if T is the regular tournament, T^op is also
regular but with t₃ and c₃ SWAPPED... no, complement reversal doesn't
swap transitive and cyclic triples.

Wait — for the regular tournament at n=5, t₃ = c₃ = 5, so
C_T and C_R have the SAME number of rows. And by the dihedral symmetry
of the regular tournament, they might even be isomorphic.

The golden ratio eigenvalues (φ, φ̄) are SHARED between C_T^TC_T and
C_R^TC_R precisely because of this T ↔ R symmetry at t₃ = c₃.

For n=5: φ·φ̄ = n = 5. So the pair (φ,φ̄) satisfies λ² - nλ + n = 0.

Q: Does the equation λ² - nλ + n = 0 have algebraic significance?
  Rewrite: λ(n-λ) = n, or equivalently (n-λ)(n-(n-λ)) = n.
  This says the eigenvalue and its "complement" n-μ multiply to n.
  i.e., the C_T and C_R eigenvalues are DUAL under the map μ ↦ n-μ,
  and when they're equal (as at the regular tournament), they satisfy
  μ(n-μ) = n.
""")

# Verify: at any t₃, do the C_R^TC_R eigenvalues on im(P) pair with C_T^TC_T?
print("Verifying C_T^TC_T and C_R^TC_R duality at n=5:")

n = 5
E = n*(n-1)//2
C_full, edges, edge_idx = full_constraint_matrix(n)
U_f, S_f, Vt_f = np.linalg.svd(C_full, full_matrices=False)
nz_f = sum(1 for s in S_f if s > 1e-8)
P_basis = Vt_f[:nz_f]

# Check one example per t₃
seen_t3 = set()
for bits in range(1 << E):
    adj = tournament_from_bits(n, bits)
    C_T, C_R = constraint_matrices(adj, n, C_full, edge_idx)
    t3 = C_T.shape[0]
    c3 = C_R.shape[0]
    if t3 in seen_t3 or t3 == 0 or c3 == 0:
        continue
    seen_t3.add(t3)

    CTC_T = C_T.T @ C_T
    CTC_R = C_R.T @ C_R

    # Restrict to im(P)
    CTC_T_P = P_basis @ CTC_T @ P_basis.T
    CTC_R_P = P_basis @ CTC_R @ P_basis.T

    eigs_T = np.sort(np.linalg.eigvalsh(CTC_T_P))[::-1]
    eigs_R = np.sort(np.linalg.eigvalsh(CTC_R_P))[::-1]

    print(f"\n  t₃={t3}, c₃={c3}:")
    print(f"    C_T eigs on im(P): {eigs_T.round(3)}")
    print(f"    C_R eigs on im(P): {eigs_R.round(3)}")
    print(f"    Sum check: {(eigs_T + eigs_R).round(3)}")
    # Should all be n

    # Check: sorted C_T eigenvalues + sorted C_R eigenvalues (reversed) = n?
    eigs_R_rev = eigs_R[::-1]
    print(f"    C_T + rev(C_R): {(eigs_T + eigs_R_rev).round(3)}")

# ====================================================================
print("\n\n" + "=" * 70)
print("THE 3/n BRIDGE: UNIFIED PERSPECTIVE")
print("=" * 70)
# ====================================================================

print("""
SYNTHESIS: The fraction 3/n connects THM-217 and THM-224 through
a DUALITY between transitive and cyclic triples.

THM-224 (Homological):
  C^TC = n·P, where rank(P) = (n-1)(n-2)/2 = C(n,3) · (3/n).
  The uniform eigenvalue n arises because:
    n = trace/rank = 3·C(n,3) / [(n-1)(n-2)/2] = n.

  The number 3 comes from dim(∂Δ²) = 3 (each triangle has 3 edges).
  The number (n-1)(n-2)/2 comes from H₁(K_n) = 0 (K_n is simply connected).

THM-217 (HP-counting):
  M(x) = [[1,2x,0],[0,0,1],[1,x,0]], char poly λ³-λ²-xλ-x.
  At x = n(n-3) = n²·(1-3/n), the char poly factors as:
    (λ-(3-n))(λ²-(n-2)λ-n) = 0

  Eigenvalues: α=3-n, λ₁=M_n-1, λ₂=-n/λ₁.

  α = 3-n encodes the independence fraction: α/(-n) = (n-3)/n = 1-3/n.

THE BRIDGE:
  Both structures encode the same combinatorial fact:
    "Each edge appears in n-2 triples, and each triple involves 3 edges."

  In THM-224: this gives trace = 3·C(n,3) = E·(n-2), eigenvalue = n.
  In THM-217: this gives x = n(n-3) as the "redundancy capacity" of n².

  The eigenvalue n in THM-224 and the parameter x=n(n-3) in THM-217 are
  related by:
    x = n(n-3) = n·(n-3) = n·(max_rank - 3)  ... hmm

  Actually: x = n(n-3) = n·E·(3/n)·[(n-3)/3] ... let me think.

  The cleaner relationship:
    3·C(n,3) = n · rank(C) = n · (n-1)(n-2)/2
    n(n-3) = 3·C(n,3) - 3·rank(C) = n·rank(C) - 3·rank(C) = (n-3)·rank(C)

  So: x = (n-3) · rank(C) = (n-3)(n-1)(n-2)/2 = C(n,3)·(n-3)·3/(n·1)...
  Actually:
    x = n(n-3) and rank = (n-1)(n-2)/2
    x/rank = 2n(n-3)/((n-1)(n-2)) = 2n(n-3)/((n-1)(n-2))

  At n=6: x/rank = 2·6·3/(5·4) = 36/20 = 1.8. Not obviously clean.

  Better: x = C(n,3) - rank = REDUNDANCY COUNT
    C(n,3) - rank = C(n,3) - (n-1)(n-2)/2 = C(n-1,3)
    And C(n-1,3) = (n-1)(n-2)(n-3)/6.
    But x = n(n-3). And C(n-1,3) = (n-1)(n-2)(n-3)/6.

    x/C(n-1,3) = n(n-3)/[(n-1)(n-2)(n-3)/6] = 6n/[(n-1)(n-2)]
    At n=5: 10/4 = 2.5. At n=6: 18/10 = 1.8.

  Hmm, x ≠ redundancy count exactly. Let me recheck.
""")

for n_val in range(3, 10):
    x = n_val*(n_val-3)
    rank = (n_val-1)*(n_val-2)//2
    Cn3 = comb(n_val, 3)
    redundancy = Cn3 - rank  # = C(n-1,3)
    print(f"  n={n_val}: x=n(n-3)={x}, rank={rank}, C(n,3)={Cn3}, "
          f"redundancy=C(n-1,3)={redundancy}, x/redundancy={x/redundancy:.3f}" if redundancy > 0 else
          f"  n={n_val}: x={x}, rank={rank}, redundancy={redundancy}")

print("""
x/redundancy = n(n-3)/C(n-1,3) = n(n-3)/[(n-1)(n-2)(n-3)/6] = 6n/[(n-1)(n-2)]

At n=6: 6·6/(5·4) = 36/20 = 9/5 = 1.8
At n=∞: → 6n/n² = 6/n → 0

So x and redundancy grow at DIFFERENT rates. The connection is NOT
x = redundancy_count, but x = n² · (redundancy_fraction).

The fraction 3/n appears in both:
  THM-224: independence_fraction = 3/n
  THM-217: x = n² · (1 - 3/n) = n² - 3n

This means:
  x + 3n = n²
  x = n² - 3n = n(n-3)

The "3n" that's subtracted from n² is 3 times the system size.
This is the "cost" of each triple involving 3 edges.
The remaining n(n-3) is the "free capacity" for non-trivial structure.
""")

# ====================================================================
print("\n" + "=" * 70)
print("FINAL: THE EIGENVALUE EQUATION λ²-nλ+n=0 AS A BRIDGE")
print("=" * 70)
# ====================================================================

print("""
The equation λ²-nλ+n=0 appears in BOTH frameworks:

1. THM-224 context: At the regular tournament (t₃=c₃), the eigenvalues
   of C_T^TC_T on the non-top subspace satisfy λ(n-λ) = n.
   This is the SELF-DUALITY condition: an eigenvalue λ of C_T^TC_T
   corresponds to eigenvalue n-λ of C_R^TC_R, and when T↔R symmetry
   holds, λ = n-λ' where λ' is another eigenvalue of C_T^TC_T,
   forming pairs with product n.

2. THM-217 context: The metallic quadratic λ²-(n-2)λ-n=0 gives
   eigenvalues of M(x) at the special point. Compare:
     Metallic: λ² = (n-2)λ + n  (product = -n, sum = n-2)
     THM-224:  λ² = nλ - n      (product = +n, sum = n)

   These are related by the substitution λ ↦ n-λ:
     If μ = n-λ satisfies μ²-nμ+n=0, then
     (n-λ)² - n(n-λ) + n = n²-2nλ+λ² - n²+nλ + n = λ²-nλ+n = 0.
   Same equation! So the substitution λ↦n-λ is an INVOLUTION.

   The metallic equation λ²-(n-2)λ-n=0 is DIFFERENT: it gives
   the THM-217 eigenvalues, not the THM-224 eigenvalues.

   But: metallic eigenvalue λ₁ satisfies λ₁²=(n-2)λ₁+n.
   THM-224 eigenvalue φ satisfies φ²=nφ-n.

   At n=5: λ₁ = (3+√29)/2 ≈ 4.193, φ = (5+√5)/2 ≈ 3.618.
   These are NOT the same.

   HOWEVER: the DISCRIMINANTS are related.
   Metallic: Δ = (n-2)²+4n = n²+4
   THM-224:  Δ = n²-4n = n(n-4)

   Ratio: (n²+4)/(n²-4n) = (n²+4)/(n(n-4))
   At n=5: 29/5. At n=6: 40/12 = 10/3.

   Sum of discriminants: (n²+4) + n(n-4) = 2n²-4n+4 = 2(n²-2n+2) = 2((n-1)²+1)
   Difference: (n²+4) - n(n-4) = 4n+4 = 4(n+1)

   So Δ_metallic - Δ_THM224 = 4(n+1). At n=5: 24. At n=6: 28.
""")

# Numerical comparison
for n_val in [5, 6, 7, 8]:
    disc_met = n_val**2 + 4
    disc_224 = n_val * (n_val - 4)

    metallic = (n_val - 2 + np.sqrt(disc_met)) / 2
    if disc_224 >= 0:
        phi_224 = (n_val + np.sqrt(disc_224)) / 2
    else:
        phi_224 = complex(n_val/2, np.sqrt(-disc_224)/2)

    print(f"  n={n_val}: metallic={metallic:.4f}, φ_224={phi_224 if isinstance(phi_224, complex) else f'{phi_224:.4f}'}")
    print(f"    Δ_met={disc_met}, Δ_224={disc_224}, diff={disc_met-disc_224}=4(n+1)={4*(n_val+1)}")

    if isinstance(phi_224, float):
        # Is there any transformation between them?
        # metallic: μ²-(n-2)μ-n=0 → μ = [(n-2)+√(n²+4)]/2
        # φ_224: φ²-nφ+n=0 → φ = [n+√(n²-4n)]/2
        #
        # What if we evaluate the metallic quadratic at φ_224?
        val = phi_224**2 - (n_val-2)*phi_224 - n_val
        print(f"    metallic eq at φ_224: {val:.4f}")
        # And THM-224 at metallic?
        val2 = metallic**2 - n_val*metallic + n_val
        print(f"    THM-224 eq at metallic: {val2:.4f}")
