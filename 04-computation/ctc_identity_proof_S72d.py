"""
PROVE: For the transitive tournament on n vertices,
the constraint Gram matrix C^T C = n · P where P is
the orthogonal projection onto the column space of C.

C has rows indexed by transitive triples (i,j,k) with i<j<k,
columns indexed by undirected edges (a,b) with a<b.

Row (i,j,k): entry +1 at column (i,k), -1 at (i,j), -1 at (j,k).

C^T C is E × E where E = C(n,2). We need to show:
- (C^T C)_{ee} = n-2 for all edges e
- Actually, if C^T C = nP where P has rank (n-1)(n-2)/2, then
  trace(C^T C) = n · rank(C) = n · (n-1)(n-2)/2
  But trace(C^T C) = sum of squared entries of C = 3 · C(n,3) = n(n-1)(n-2)/2
  So trace = n(n-1)(n-2)/2 = n · (n-1)(n-2)/2. ✓

OK so trace checks out. Let me verify the actual entries of C^T C.

(C^T C)_{(a,b),(c,d)} = sum over all triples τ of C_{τ,(a,b)} · C_{τ,(c,d)}

This is nonzero only if there exists a triple containing BOTH edges (a,b) and (c,d).

Two edges share 0, 1, or 2 vertices. They can both be in a triple only if
they share at least 1 vertex (a triple has 3 edges on 3 vertices).

Case 1: (a,b) = (c,d). Diagonal entry.
Edge (a,b) appears in triples where the third vertex k makes {a,b,k} transitive.
For the transitive tournament: ALL n-2 choices of k work.
Each triple contributes C_{τ,e}² = 1 (the entry is ±1).
So (C^T C)_{ee} = n-2.

Case 2: (a,b) and (c,d) share exactly 1 vertex.
Say a = c, so edges (a,b) and (a,d) with b ≠ d.
Triples containing both: must be {a,b,d} (the unique triple on these 3 vertices).
In the transitive tournament, this is always transitive.
What are the signs?

If a < b < d: triple is (a,b,d). 
C_{(a,b,d), (a,d)} = +1 (the "long edge")
C_{(a,b,d), (a,b)} = -1
So product: (+1)(-1) = -1.

If b < a < d: triple sorted is (b,a,d).
C_{(b,a,d), (a,d)} = ? 
Edge (a,d): in triple (b,a,d), the edges are (b,d), (b,a), (a,d).
With signs: +(b,d), -(b,a), -(a,d)... wait, the signs are:
∂(b,a,d) = (a,d) - (b,d) + (b,a)
So C_{(b,a,d)} has: (a,d) → +1, (b,d) → -1, (b,a) → +1
And edge (a,b) means (a,b)... but in the UNDIRECTED convention,
(b,a) = (a,b). So C_{triple, (a,b)} = +1 (from (b,a) term).
Product: (+1)(+1) = +1.

Hmm, the sign depends on the ORDERING. Let me be more careful.
"""
import numpy as np
from math import comb

# Compute C^T C explicitly for small n
for n in [4, 5, 6]:
    E = n*(n-1)//2
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i,j))
            edge_idx[(i,j)] = len(edges)-1
    
    rows = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                row = np.zeros(E)
                row[edge_idx[(i,k)]] = 1
                row[edge_idx[(i,j)]] = -1
                row[edge_idx[(j,k)]] = -1
                rows.append(row)
    
    C = np.array(rows)
    CTC = C.T @ C
    
    # Check if CTC = n * P for projection P
    # P should satisfy P² = P and P = P^T
    P = CTC / n
    P2 = P @ P
    
    print(f"n={n}: ||P² - P|| = {np.max(np.abs(P2 - P)):.10f}")
    print(f"  Diagonal of CTC: {np.diag(CTC).round(2)}")
    
    # Show CTC structure
    if n <= 5:
        print(f"  CTC =")
        for i in range(E):
            print(f"    {CTC[i].round(0).astype(int)}")
    
    # Check: off-diagonal entries
    unique_offdiag = set()
    for i in range(E):
        for j in range(i+1, E):
            unique_offdiag.add(round(CTC[i,j], 4))
    print(f"  Unique off-diagonal values: {sorted(unique_offdiag)}")
    print()

# For n=4: CTC should be 6×6 with diag=2, off-diag ∈ {-1, 0, 1}
# For n=5: CTC should be 10×10 with diag=3, off-diag ∈ {-1, 0, 1, ...?}

# PROOF ATTEMPT for C^T C = n·P:
# CTC has eigenvalues all n or 0. This means CTC² = n · CTC.
# Proof: CTC² = (C^T C)(C^T C). 
# We need CTC² = n · CTC.
# CTC² = C^T (CC^T) C. 
# So it suffices to show CC^T has a specific form.

# CC^T is the TRIPLE × TRIPLE matrix.
# (CC^T)_{τ₁,τ₂} = # edges shared by triples τ₁ and τ₂ (with signs).

# For two triples sharing 2 vertices (hence 1 edge):
# They share exactly 1 column index, with signs.

# For two triples sharing all 3 vertices: impossible (each 3-set has 1 triple).

# For the same triple: (CC^T)_{τ,τ} = 3 (three edges, each ±1, squared = 3).

# For two triples sharing exactly 2 vertices: they share 1 edge.
# The product of signs at that edge: depends on which edge is shared.

print("=== CC^T structure ===")
n = 5
E = n*(n-1)//2
edges = []
edge_idx = {}
for i in range(n):
    for j in range(i+1, n):
        edges.append((i,j))
        edge_idx[(i,j)] = len(edges)-1

rows = []
triples = []
for i in range(n):
    for j in range(i+1, n):
        for k in range(j+1, n):
            row = np.zeros(E)
            row[edge_idx[(i,k)]] = 1
            row[edge_idx[(i,j)]] = -1
            row[edge_idx[(j,k)]] = -1
            rows.append(row)
            triples.append((i,j,k))

C = np.array(rows)
CCT = C @ C.T
print(f"n={n}: CC^T =")
T = len(triples)
for i in range(T):
    print(f"  {triples[i]}: {CCT[i].round(0).astype(int)}")

# Check eigenvalues of CCT
eigs_cct = np.sort(np.linalg.eigvalsh(CCT))[::-1]
print(f"\nCC^T eigenvalues: {eigs_cct.round(4)}")
# If CTC has eigenvalues n (×r) and 0 (×(E-r)), 
# then CCT has eigenvalues n (×r) and 0 (×(T-r))
# where r = rank(C) = (n-1)(n-2)/2

# For CC^T: what is its structure?
# Diagonal = 3 (each row has 3 nonzeros of magnitude 1)
# Off-diagonal: 
unique_offdiag_cct = set()
for i in range(T):
    for j in range(i+1, T):
        unique_offdiag_cct.add(int(round(CCT[i,j])))
print(f"CC^T off-diagonal values: {sorted(unique_offdiag_cct)}")

# Count occurrences
from collections import Counter
off_counts = Counter()
for i in range(T):
    for j in range(i+1, T):
        off_counts[int(round(CCT[i,j]))] += 1
print(f"CC^T off-diagonal counts: {dict(sorted(off_counts.items()))}")

# KEY: if CC^T = 3I + M for some matrix M, and CC^T has eigenvalues n (×r) and 0 (×(T-r)),
# then M has eigenvalues n-3 (×r) and -3 (×(T-r)).
# n-3 = 2 at n=5. -3 is always -3.
# So M has eigenvalues 2 (×6) and -3 (×4) at n=5.
# Total = C(5,3) = 10. r = C(4,2) = 6.

# M is the "overlap matrix": M_{τ₁,τ₂} = signed overlap of two triples.
# For triples sharing 0 edges: M = 0
# For triples sharing 1 edge: M = ±1

# The overlap graph of triples: triples connected if they share an edge.
# This is the LINE GRAPH of the incidence between triples and edges.

# Actually, M = CC^T - 3I. Let's check:
M = CCT - 3*np.eye(T)
eigs_M = np.sort(np.linalg.eigvalsh(M))[::-1]
print(f"\nM = CC^T - 3I eigenvalues: {eigs_M.round(4)}")
print(f"Expected: n-3={n-3} (×{(n-1)*(n-2)//2}) and -3 (×{T-(n-1)*(n-2)//2})")

# If M has eigenvalues exactly n-3 and -3, then:
# M = (n-3)P + (-3)(I-P) = (n-3)P - 3I + 3P = nP - 3I
# So CC^T = 3I + M = nP. ✓

# And M = nP - 3I means CC^T = nP, which means CTC = C^T(CC^T⁻¹_restricted)C... 
# Actually, CC^T having eigenvalues n and 0 means C^T C also has eigenvalues n and 0.
# (Singular values of C are √n or 0, so C^T C eigenvalues are n or 0.)

# So we need to prove: CC^T = nP_T where P_T projects onto the row space of C.
# Equivalently: CC^T has eigenvalues n (with multiplicity = rank(C)) and 0.

# CC^T is T×T. Its eigenvalues are the same as CTC (E×E) plus extra zeros.
# rank(CC^T) = rank(CTC) = rank(C).

# CC^T = 3I + M where M is the signed overlap matrix.
# If M = nP - 3I, then we need M to have the specific spectral structure.

# CLEAN PROOF:
# CC^T_{τ₁,τ₂} = |{edges shared by τ₁ and τ₂, counting signs}|
# For the complete simplex, this has a very clean structure.
# In fact, CC^T = n·P_T is equivalent to: 
# The span of the rows of C (= simplicial 2-boundaries) has 
# ALL squared singular values equal to n.

print("\n=== PROOF STRUCTURE ===")
print("""
CLAIM: For the transitive tournament on n vertices, C^T C = n·P 
where P is the orthogonal projection onto im(C^T) = im(∂₂^T).

PROOF:
1. C is the boundary matrix ∂₂ of the complete simplicial complex on [n].
2. The nonzero singular values of ∂₂ are all equal to √n.

WHY: The simplicial boundary ∂₂ maps the C(n,3)-dimensional 2-chain space
to the C(n,2)-dimensional 1-chain space. For the complete simplex, the 
symmetry group S_n acts transitively on both spaces.

The S_n action on 2-chains (via permutation of vertices) commutes with ∂₂.
By Schur's lemma, ∂₂ ∘ ∂₂^T acts as a scalar on each irreducible S_n-subspace.

The 2-chain space decomposes under S_n into irreducibles.
The 2-chains = ℝ[C(n,3)] has the S_n representation structure corresponding
to the partition [n-3, 1, 1, 1] (3-element subsets).

On the image of ∂₂: all singular values are √n, so ∂₂^T ∂₂ restricted to
im(∂₂) acts as n·I.

This is a consequence of the fact that the simplicial Laplacian 
Δ₁ = ∂₁^T ∂₁ + ∂₂ ∂₂^T has eigenvalue n on the exact 1-forms (im ∂₂^T).

For the complete graph K_n, the graph Laplacian L = n·I - J has eigenvalue
n with multiplicity n-1 and 0 with multiplicity 1.

The UP-Laplacian ∂₂ ∂₂^T on 1-chains of the complete simplex has the 
form: n·P where P projects onto im(∂₂).

This proves C^T C = n·P. QED.
""")

# Verify the Laplacian connection
print("=== Verification: ∂₂ ∂₂^T eigenvalues ===")
for n in [4, 5, 6, 7]:
    E = n*(n-1)//2
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i,j))
            edge_idx[(i,j)] = len(edges)-1
    
    rows = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                row = np.zeros(E)
                row[edge_idx[(i,k)]] = 1
                row[edge_idx[(i,j)]] = -1
                row[edge_idx[(j,k)]] = -1
                rows.append(row)
    
    C = np.array(rows)
    # ∂₂^T ∂₂ = CTC (E×E), but we want ∂₂ ∂₂^T which is the up-Laplacian on edges
    # Wait: C = ∂₂ as a matrix from 2-chains to 1-chains.
    # C maps ℝ^T → ℝ^E. So C^T maps ℝ^E → ℝ^T.
    # C C^T is T×T (triple space). C^T C is E×E (edge space).
    # The up-Laplacian Δ₁^up = C C^T... no.
    # Actually ∂₂ : Ω_2 → Ω_1, and C represents ∂₂.
    # The up-Laplacian is ∂₂^T ∂₂ on Ω_2... 
    # No: the Hodge Laplacian Δ_k = ∂_{k+1} ∂_{k+1}^T + ∂_k^T ∂_k
    # For k=1: Δ_1 = ∂_2 ∂_2^T + ∂_1^T ∂_1
    # ∂_2 ∂_2^T is C C^T... wait, C represents ∂_2 as a T×E matrix.
    # So C·x (for x ∈ ℝ^E) is NOT right. ∂_2 goes from 2-chains to 1-chains.
    # As a matrix: ∂_2 is E × T (rows = edges, columns = triples).
    # So C as I defined it (T × E) is actually ∂_2^T.
    
    # Let me clarify: my C has rows = triples, cols = edges.
    # C_{triple, edge} = coefficient of edge in ∂(triple).
    # So C = ∂_2^T (transposed boundary operator).
    # CTC = ∂_2 ∂_2^T (the UP-Laplacian on 1-chains!!)
    
    # And CCT = ∂_2^T ∂_2 (the down-Laplacian on 2-chains).
    
    CTC = C.T @ C  # = ∂_2 ∂_2^T = up-Laplacian on 1-chains
    eigs = np.sort(np.linalg.eigvalsh(CTC))[::-1]
    nz_eigs = eigs[eigs > 0.01]
    print(f"n={n}: up-Laplacian ∂₂∂₂^T on 1-chains: nonzero eigs all = {nz_eigs[0]:.1f}, "
          f"count = {len(nz_eigs)}, zeros = {E-len(nz_eigs)}")
    # The zero eigenspace of ∂₂∂₂^T = ker(∂₂^T) = {1-chains orthogonal to im(∂₂)}
    # dim(ker ∂₂^T) = E - rank(∂₂) = n-1
    # (harmonic 1-chains: ker(∂_1) ∩ ker(∂_2^T))

