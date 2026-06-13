#!/usr/bin/env python3
"""
β_2 = 0: THE SIMPLICIAL ARGUMENT

The transitive triples of a tournament T form a simplicial complex Δ(T):
  - Vertices = vertices of T
  - Edges = directed edges a→b
  - 2-simplices = transitive triples (a,b,c) with a→b→c and a→c
  - 3-simplices = transitive 4-subsets (a,b,c,d) = totally ordered 4-cliques

The GLMY path homology of T at dimension 2 is:
  H_2^{path}(T) = ker(∂_2: Ω_2 → Ω_1) / im(∂_3: Ω_3 → Ω_2)

We showed: im(∂_3|DT) suffices (DT paths = individual 3-simplices in Δ).

QUESTION: Is Δ(T) actually the ORDER COMPLEX of the partial order
given by the tournament? No — a tournament isn't a partial order.
But the TRANSITIVE CLOSURE gives a partial order on the strongly
connected components.

ALTERNATIVE: The "transitive 4-clique" structure is related to the
FLAG COMPLEX of the tournament (edges = directed edges, higher
simplices = cliques). But this isn't standard simplicial homology
because of the directionality.

KEY OBSERVATION: For the transitive tournament on n vertices
(total order 0→1→2→...→n-1), the transitive triples are ALL
C(n,3) triples. The homology of this simplicial complex is that
of the (n-1)-simplex, which is trivial.

For a general tournament, the transitive triples form a SUBCOMPLEX
of the full simplex. Its homology might not be trivial.

BUT: the GLMY boundary is NOT the simplicial boundary!
∂^{GLMY}(a,b,c) = (b,c) - (a,c) + (a,b) depends on DIRECTED edges,
while simplicial ∂(a,b,c) = (b,c) - (a,c) + (a,b) looks the same
BUT the simplicial version has NO orientation constraint.

Actually wait — in GLMY, the vertex ordering matters!
∂(a,b,c) = (b,c) - (a,c) + (a,b) with DIRECTED meaning.

Hmm, for the simplicial complex of transitive triples, with the
INDUCED ordering from the tournament, the GLMY boundary IS the
simplicial boundary. The key question is whether H_2^{simp}(Δ(T)) = 0.

Let me just TEST this hypothesis directly.
"""
import numpy as np
from itertools import combinations
from collections import Counter
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def tournament_order(A, vertices):
    """Return vertices sorted by tournament domination on this subset.
    For a transitive subtournament, this gives the unique total order."""
    vlist = list(vertices)
    # Sort by number of wins within subset
    wins = {v: sum(A[v][u] for u in vlist if u != v) for v in vlist}
    return sorted(vlist, key=lambda v: -wins[v])

# ===== Is the simplicial complex of transitive triples acyclic? =====
print("=" * 70)
print("SIMPLICIAL COMPLEX OF TRANSITIVE TRIPLES: H_2 = 0?")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    h2_nonzero = 0
    total = 0

    for A in all_tournaments_gen(n):
        total += 1
        # Find all transitive triples — these are 2-simplices
        tt = []
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    if len({a,b,c}) == 3 and A[a][b] and A[b][c] and A[a][c]:
                        tt.append((a,b,c))

        # Find all transitive 4-cliques — these are 3-simplices
        t4 = []
        for quad in combinations(range(n), 4):
            # Check if transitive
            ordered = tournament_order(A, quad)
            is_trans = all(A[ordered[i]][ordered[j]] for i in range(4) for j in range(i+1,4))
            if is_trans:
                t4.append(tuple(ordered))

        if len(tt) == 0:
            continue

        # Build ∂_2: span(tt) → span(edges)
        # directed edges in tournament
        edges = [(i,j) for i in range(n) for j in range(n) if A[i][j]]
        edge_idx = {e: i for i, e in enumerate(edges)}

        bd2 = np.zeros((len(edges), len(tt)))
        for j, (a, b, c) in enumerate(tt):
            # ∂(a,b,c) = (b,c) - (a,c) + (a,b)
            bd2[edge_idx[(b,c)], j] += 1
            bd2[edge_idx[(a,c)], j] -= 1
            bd2[edge_idx[(a,b)], j] += 1

        rank_bd2 = np.linalg.matrix_rank(bd2, tol=1e-8)
        ker_dim = len(tt) - rank_bd2

        if ker_dim == 0:
            continue

        # Build ∂_3: span(t4) → span(tt)
        if len(t4) == 0:
            h2 = ker_dim
        else:
            tt_idx = {t: i for i, t in enumerate(tt)}
            bd3 = np.zeros((len(tt), len(t4)))
            for j, clique in enumerate(t4):
                a,b,c,d = clique
                faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
                signs = [1, -1, 1, -1]
                for face, sign in zip(faces, signs):
                    if face in tt_idx:
                        bd3[tt_idx[face], j] += sign

            im_rank = np.linalg.matrix_rank(bd3, tol=1e-8)
            h2 = ker_dim - im_rank

        if h2 != 0:
            h2_nonzero += 1

    print(f"  Total: {total}, H_2 ≠ 0: {h2_nonzero}")

# ===== Now compare: is this the SAME as GLMY path homology? =====
print(f"\n\n{'='*70}")
print("COMPARISON: SIMPLICIAL H_2 vs GLMY H_2")
print("="*70)

# In GLMY, Ω_2 = transitive triples (same as simplicial 2-faces).
# But Ω_3 includes MORE than just DT 4-paths (also cancellation chains).
# The simplicial ∂_3 uses ONLY 4-cliques (= DT paths).
#
# So: simplicial H_2 ≥ GLMY H_2 (GLMY has more in im(∂_3)).
# We showed GLMY H_2 = 0, and DT-sufficiency says simplicial H_2 = 0 too.
# So they agree at H_2!

# But wait — the simplicial complex here is the "transitive closure flag complex"
# or "acyclicity complex" of T. Its homology being trivial at dim 2
# is itself a nontrivial statement.

# Let me check what this simplicial complex IS for small tournaments.

n = 5
print(f"\nSimplicial complex structure at n={n}:")
face_stats = Counter()

for t_idx, A in enumerate(all_tournaments_gen(n)):
    # Count: transitive k-subsets
    t2 = sum(1 for a in range(n) for b in range(n) for c in range(n)
             if len({a,b,c})==3 and A[a][b] and A[b][c] and A[a][c])
    t3 = 0
    for quad in combinations(range(n), 4):
        ordered = tournament_order(A, quad)
        if all(A[ordered[i]][ordered[j]] for i in range(4) for j in range(i+1,4)):
            t3 += 1
    t4 = 0
    ordered = tournament_order(A, range(n))
    if all(A[ordered[i]][ordered[j]] for i in range(5) for j in range(i+1,5)):
        t4 += 1
    face_stats[(t2, t3, t4)] += 1

print("(|Ω_2|, |transitive 4-cliques|, |transitive 5-cliques|): count")
for key in sorted(face_stats.keys()):
    print(f"  {key}: {face_stats[key]}")

# ===== THE DEEPER QUESTION: Is this related to shellability? =====
print(f"\n\n{'='*70}")
print("SHELLABILITY / CONTRACTIBILITY OF THE TRANSITIVE COMPLEX")
print("="*70)

print("""
The complex Δ(T) of transitive subsets of a tournament T is known in
combinatorial topology as the "acyclicity complex" or "order complex"
of the tournament.

CLAIM: For any tournament T, the acyclicity complex Δ(T) has the
homology of a WEDGE OF SPHERES, all of dimension β_1(T).

If β_1(T) = 0 (no directed 3-cycles): Δ(T) is contractible (it's
a simplex minus nothing — the whole tournament is transitive).

If β_1(T) = k: Δ(T) has H_1 = Z^k and all higher H = 0.

The β_2 = 0 result then says: the transitive triple complex has
NO 2-dimensional holes, regardless of the tournament structure.

This could be related to COHEN-MACAULAY properties or shellability
of the acyclicity complex.

Let's verify: does the acyclicity complex always have H_k = 0 for k ≥ 2?
""")

# ===== Check H_1 of the simplicial complex =====
n = 5
h1_stats = Counter()
h2_stats = Counter()

for A in all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(n) if A[i][j]]
    edge_idx = {e: i for i, e in enumerate(edges)}

    # 0-simplices: vertices
    # 1-simplices: edges a→b
    # 2-simplices: transitive triples

    tt = []
    for a in range(n):
        for b in range(n):
            for c in range(n):
                if len({a,b,c})==3 and A[a][b] and A[b][c] and A[a][c]:
                    tt.append((a,b,c))

    t4 = []
    for quad in combinations(range(n), 4):
        ordered = tournament_order(A, quad)
        if all(A[ordered[i]][ordered[j]] for i in range(4) for j in range(i+1,4)):
            t4.append(tuple(ordered))

    # ∂_1: edges → vertices
    bd1 = np.zeros((n, len(edges)))
    for j, (a, b) in enumerate(edges):
        bd1[b, j] += 1
        bd1[a, j] -= 1
    rank_bd1 = np.linalg.matrix_rank(bd1, tol=1e-8)
    ker_bd1 = len(edges) - rank_bd1

    # ∂_2: triples → edges
    if len(tt) > 0:
        bd2 = np.zeros((len(edges), len(tt)))
        for j, (a, b, c) in enumerate(tt):
            bd2[edge_idx[(b,c)], j] += 1
            bd2[edge_idx[(a,c)], j] -= 1
            bd2[edge_idx[(a,b)], j] += 1
        rank_bd2 = np.linalg.matrix_rank(bd2, tol=1e-8)
        im_bd2 = rank_bd2
        ker_bd2 = len(tt) - rank_bd2
    else:
        im_bd2 = 0
        ker_bd2 = 0

    h1 = ker_bd1 - im_bd2

    # H_2
    if ker_bd2 == 0:
        h2 = 0
    elif len(t4) == 0:
        h2 = ker_bd2
    else:
        tt_idx = {t: i for i, t in enumerate(tt)}
        bd3 = np.zeros((len(tt), len(t4)))
        for j, clique in enumerate(t4):
            a,b,c,d = clique
            faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in tt_idx:
                    bd3[tt_idx[face], j] += sign
        im_rank = np.linalg.matrix_rank(bd3, tol=1e-8)
        h2 = ker_bd2 - im_rank

    h1_stats[h1] += 1
    h2_stats[h2] += 1

print(f"\nn={n}: Simplicial H_1 of acyclicity complex:")
for k in sorted(h1_stats.keys()):
    print(f"  H_1 = {k}: {h1_stats[k]} tournaments")

print(f"\nSimplicial H_2 of acyclicity complex:")
for k in sorted(h2_stats.keys()):
    print(f"  H_2 = {k}: {h2_stats[k]} tournaments")

# ===== Compare with GLMY β_1 =====
print(f"\nGLMY β_1 at n={n}: 720 with β_1=0, 304 with β_1=1")
print("Simplicial H_1 should match GLMY β_1 (since Ω_1 = A_1 for tournaments)")

print("\nDone.")
