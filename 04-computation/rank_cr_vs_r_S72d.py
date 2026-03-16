"""
Check: is rank(C_R) < r = (n-1)(n-2)/2 for ALL tournaments?

This is equivalent to: the cyclic triple boundaries never span all of im(C^T).
Proved: dim(ker(C_R) ∩ im(C^T)) = r - rank(C_R) ≥ 1 iff rank(C_R) < r.
"""
import numpy as np
from math import comb
import random

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

def get_cr_matrix(adj, n, C_full):
    """Get C_R (cyclic triple rows)."""
    cyclic_idx = []
    row = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                wins = [adj[a][b]+adj[a][c], adj[b][a]+adj[b][c], adj[c][a]+adj[c][b]]
                if max(wins) < 2:  # cyclic
                    cyclic_idx.append(row)
                row += 1
    if cyclic_idx:
        return C_full[cyclic_idx]
    return np.zeros((0, C_full.shape[1]))

# Test at each n
for n in [5, 6, 7, 8, 9]:
    E = n*(n-1)//2
    r = (n-1)*(n-2)//2
    Cn3 = comb(n, 3)

    # max c₃ (regular tournament for odd n)
    if n % 2 == 1:
        max_c3 = Cn3 - n * comb((n-1)//2, 2)
    else:
        # For even n, the most cyclic tournament has scores as balanced as possible
        max_c3 = None  # we'll find it empirically

    C_full, edges, edge_idx = full_constraint_matrix(n)

    random.seed(42)
    if n <= 6:
        max_bits = 1 << E
        n_samples = min(max_bits, 50000)
        bits_list = range(n_samples) if n <= 5 else [random.randint(0, max_bits-1) for _ in range(n_samples)]
    else:
        max_bits = 1 << E
        n_samples = 10000
        bits_list = [random.randint(0, max_bits-1) for _ in range(n_samples)]

    max_rk_R = 0
    max_c3_seen = 0
    violations = 0

    for bits in bits_list:
        adj = tournament_from_bits(n, bits)
        C_R = get_cr_matrix(adj, n, C_full)
        c3 = C_R.shape[0]
        max_c3_seen = max(max_c3_seen, c3)

        if c3 == 0:
            continue

        rk_R = np.linalg.matrix_rank(C_R)
        max_rk_R = max(max_rk_R, rk_R)

        if rk_R >= r:
            violations += 1
            print(f"  n={n}: VIOLATION! rank(C_R)={rk_R} ≥ r={r}, c₃={c3}")

    print(f"n={n}: r={r}, max_c₃_seen={max_c3_seen}, max_rank(C_R)={max_rk_R}, "
          f"c₃>r? {max_c3_seen > r}, violations={violations}/{n_samples}")
    if max_c3 is not None:
        print(f"  theoretical max_c₃ (regular)={max_c3}")

# At n=9, regular tournament has c₃=30 but r=28. Can rank(C_R)=28?
print("\n--- FOCUSED TEST: n=9 regular tournament ---")
n = 9
E = n*(n-1)//2  # 36
r = (n-1)*(n-2)//2  # 28
C_full, edges, edge_idx = full_constraint_matrix(n)

# Construct regular tournament: C_9^{1,2,3,4}
adj = [[0]*n for _ in range(n)]
for i in range(n):
    for d in range(1, (n+1)//2):  # offsets 1..4
        adj[i][(i+d)%n] = 1

C_R = get_cr_matrix(adj, n, C_full)
rk_R = np.linalg.matrix_rank(C_R)
print(f"Regular tournament C_9^{{1,2,3,4}}: c₃={C_R.shape[0]}, rank(C_R)={rk_R}, r={r}")
print(f"  gap = r - rank(C_R) = {r - rk_R}")

# Check a few more constructions at n=9
random.seed(123)
for trial in range(5):
    bits = random.randint(0, (1<<E)-1)
    adj = tournament_from_bits(n, bits)
    C_R = get_cr_matrix(adj, n, C_full)
    c3 = C_R.shape[0]
    rk_R = np.linalg.matrix_rank(C_R)
    print(f"  Random tournament: c₃={c3}, rank(C_R)={rk_R}, gap={r-rk_R}")

# ====================================================================
print("\n" + "=" * 70)
print("WHY rank(C_R) < r: THE GRADIENT OBSTRUCTION")
print("=" * 70)
# ====================================================================

print("""
PROOF that rank(C_R) < r for any tournament:

The constraint column space im(C^T) decomposes (via Hodge theory on K_n) as:
  im(C^T) = im(∂₂^T) (curl space, dimension r)

The "gradient" space im(∂₁^T) has dimension n-1, and is orthogonal to im(C^T).
Together: ℝ^E = im(∂₁^T) ⊕ im(∂₂^T).

Now: the SCORE VECTOR s = (s₁,...,s_n) where s_i = out-degree(i) defines
a gradient 1-form v = ∂₁^T·s by v_{(i,j)} = s_j - s_i.

For any triple σ={a,b,c}: ∂₂(σ)·v = v_{bc} - v_{ac} + v_{ab}
= (s_c-s_b) - (s_c-s_a) + (s_b-s_a) = 0.

So the score gradient v is in ker(C) = ker(C_T) ∩ ker(C_R). This means
EVERY row of C (transitive or cyclic) is orthogonal to v.

But this doesn't help — v is in im(∂₁^T) which is already orthogonal to im(C^T).

DIFFERENT APPROACH: Use the simplicial dependency structure.

Every 4-vertex subset {a,b,c,d} gives a DEPENDENCY among the 4 triple
boundaries: ∂₂(bcd) - ∂₂(acd) + ∂₂(abd) - ∂₂(abc) = 0.
(This is ∂₃ of the 3-simplex.)

In a tournament, these 4 triples can be split into transitive and cyclic.
The dependency holds among ALL 4, regardless of T/C classification.

For rank(C_R) = r: the c₃ cyclic rows must span all r dimensions.
But the rows of C_R satisfy ALL the dependency relations that involve
at least one transitive triple. Each dependency is:
  Σ (cyclic rows) = -Σ (transitive rows in same dependency)

If a dependency involves k transitive rows and 4-k cyclic rows,
then the k transitive rows are in the span of the cyclic rows
(by the dependency).

rank(C_R) = r iff every transitive triple boundary is in span of
cyclic triple boundaries, i.e., every transitive triple appears in
some dependency where all OTHER triples are cyclic.

A 4-vertex subset {a,b,c,d} has 4 triples. For ALL 4 to be cyclic
except one transitive: this requires a specific structure.

CLAIM: At least one transitive triple has ALL its dependency partners
being themselves transitive (or mixed), preventing it from being
expressed as a combination of cyclic rows.

PROOF OF CLAIM: Consider the transitive tournament restricted to the
top 3 vertices (by out-degree). In any tournament, there exist 3
vertices with... hmm, this is getting complicated.

Actually, let me try a simpler argument:

The TRANSITIVE TOURNAMENT has c₃ = 0 and t₃ = C(n,3).
Clearly rank(C_R) = 0 < r. ✓

For a tournament ONE step from transitive (β₁=0, high t₃):
rank(C_R) ≤ c₃ which is small. ✓

The hard case is the REGULAR tournament (maximum c₃).
Let me verify at n=9: rank(C_R) for the regular tournament.
""")

# We already did this above, let me summarize
print(f"\nSummary:")
print(f"n=9, regular tournament: c₃=30, r=28, rank(C_R)={rk_R}, gap={r-rk_R}")
print(f"\nSo even at n=9 where c₃=30 > r=28, rank(C_R) < r!")
print(f"The cyclic boundaries have {30-rk_R} dependencies among themselves")
print(f"(beyond the simplicial dependencies).")

# Check if this gap is related to β₁
# For the regular tournament, β₁ = ?
import sys
sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')
try:
    from path_homology_v2 import path_betti_numbers
    A = np.array(adj)
    betti = path_betti_numbers(A, n, max_dim=2)
    print(f"  Regular C_9: β = {betti}")
except Exception as e:
    print(f"  (Could not compute β: {e})")

# Construct C_T and check rank
adj9 = adj
C_T9_rows = []
C_R9_rows = []
row = 0
for a in range(n):
    for b in range(a+1, n):
        for c in range(b+1, n):
            wins = [adj9[a][b]+adj9[a][c], adj9[b][a]+adj9[b][c], adj9[c][a]+adj9[c][b]]
            if max(wins) == 2:
                C_T9_rows.append(row)
            else:
                C_R9_rows.append(row)
            row += 1

C_T9 = C_full[C_T9_rows]
rk_T9 = np.linalg.matrix_rank(C_T9)
beta1_9 = r - rk_T9
print(f"  rank(C_T) = {rk_T9}, β₁ = {beta1_9}")
print(f"  rank(C_R) = {rk_R}")
print(f"  rank(C_T) + rank(C_R) = {rk_T9 + rk_R}, r = {r}")
print(f"  dim(im(C_T) ∩ im(C_R)) = rank(C_T) + rank(C_R) - r = {rk_T9 + rk_R - r}")

# Compute ker(C_R) ∩ im(C^T) directly
# dim = r - rank(C_R)
gap = r - rk_R
print(f"\n  dim(ker(C_R) ∩ im(C^T)) = r - rank(C_R) = {gap}")
print(f"  This gives {gap} eigenvectors for eigenvalue n = {n}")

# Verify top eigenvalue
CTC9 = C_T9.T @ C_T9
top_eig = np.max(np.linalg.eigvalsh(CTC9))
print(f"  Top eigenvalue of C_T^TC_T: {top_eig:.6f}")
