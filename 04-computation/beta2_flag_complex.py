#!/usr/bin/env python3
"""
beta2_flag_complex.py — Does H₂(flag complex) = 0 for all tournaments?

The FLAG COMPLEX of a tournament T has:
- vertices = V(T)
- k-simplices = transitive (k+1)-subsets of V(T)

For tournaments, this is the ORDER COMPLEX of the tournament's dominance relation.

If H₂^simp(flag) = 0 for all tournaments, then every TT 2-cycle is fillable
by DT_f boundaries. Combined with DT_c + cancel for NT parts, this would
give β₂ = 0.

Author: opus-2026-03-08-S49
"""
import sys, time, random
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def is_transitive(A, verts):
    """Check if sub-tournament on verts is transitive (no 3-cycle)."""
    for i, a in enumerate(verts):
        for j, b in enumerate(verts):
            if j <= i: continue
            for k, c in enumerate(verts):
                if k <= j: continue
                total = A[a][b] + A[b][c] + A[c][a]
                if total == 0 or total == 3:
                    return False
    return True


def flag_homology_2(A, n):
    """Compute H₂ of the flag complex of tournament A.

    Flag complex:
    - 2-simplices = TT triples (a,b,c) with a→b→c, a→c
    - 3-simplices = transitive 4-subsets

    H₂ = ker(∂₂) / im(∂₃)
    """
    # Enumerate TT triples (2-simplices)
    tt_triples = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b) or not A[b][c]: continue
                if A[a][c]:  # TT: a→b→c, a→c
                    tt_triples.append((a, b, c))

    if not tt_triples:
        return 0, 0, 0  # H₂ = 0 trivially

    # Enumerate arcs (1-simplices)
    arcs = [(i, j) for i in range(n) for j in range(n) if i != j and A[i][j]]
    arc_idx = {a: i for i, a in enumerate(arcs)}

    # Build ∂₂: TT triples → arcs
    # ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
    d2 = np.zeros((len(arcs), len(tt_triples)))
    for j, (a, b, c) in enumerate(tt_triples):
        if (b, c) in arc_idx: d2[arc_idx[(b, c)], j] += 1
        if (a, c) in arc_idx: d2[arc_idx[(a, c)], j] -= 1
        if (a, b) in arc_idx: d2[arc_idx[(a, b)], j] += 1

    rk_d2 = np.linalg.matrix_rank(d2, tol=1e-8)
    z2_simp = len(tt_triples) - rk_d2

    if z2_simp == 0:
        return 0, len(tt_triples), rk_d2

    # Enumerate transitive 4-subsets (3-simplices)
    trans4 = []
    for quad in combinations(range(n), 4):
        if is_transitive(A, quad):
            # Find the unique DT_f path through this quad
            # In a transitive 4-tournament, there's a unique total order
            verts = list(quad)
            # Sort by score within the quad
            scores = [(sum(A[v][w] for w in quad if w != v), v) for v in quad]
            scores.sort(reverse=True)
            # Path from highest to lowest score
            ordered = [v for _, v in scores]
            trans4.append(tuple(ordered))

    if not trans4:
        return z2_simp, len(tt_triples), rk_d2

    # Build ∂₃: trans4 → TT triples
    tt_idx = {t: i for i, t in enumerate(tt_triples)}
    d3 = np.zeros((len(tt_triples), len(trans4)))
    for j, (a, b, c, d) in enumerate(trans4):
        # ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
        faces = [(1, (b,c,d)), (-1, (a,c,d)), (1, (a,b,d)), (-1, (a,b,c))]
        for sign, face in faces:
            if face in tt_idx:
                d3[tt_idx[face], j] += sign

    rk_d3 = np.linalg.matrix_rank(d3, tol=1e-8)
    h2 = z2_simp - rk_d3

    return h2, len(tt_triples), rk_d2


print("=" * 70)
print("FLAG COMPLEX H₂ FOR TOURNAMENTS")
print("=" * 70)

for n in [4, 5, 6, 7]:
    m = n*(n-1)//2
    total = 1 << m

    if n <= 6:
        samples = range(total)
    else:
        random.seed(42)
        samples = random.sample(range(total), min(2000, total))

    h2_nonzero = 0
    h2_dist = Counter()

    t0 = time.time()
    for idx, bits in enumerate(samples):
        A = build_adj(n, bits)
        h2, n_tt, rk_d2 = flag_homology_2(A, n)
        h2_dist[h2] += 1
        if h2 > 0:
            h2_nonzero += 1
            if h2_nonzero <= 3:
                scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
                c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                         if A[i][j]+A[j][k]+A[k][i] in (0,3))
                print(f"  H₂≠0: bits={bits}, scores={scores}, c₃={c3}, H₂={h2}")

        if isinstance(samples, range) and bits % 10000 == 0 and bits > 0:
            elapsed = time.time() - t0
            print(f"  ... {bits}/{total} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    total_tested = len(list(samples)) if isinstance(samples, range) else len(samples)
    print(f"\nn={n}: {total_tested} tournaments in {elapsed:.1f}s")
    print(f"  H₂ distribution: {dict(sorted(h2_dist.items()))}")
    if h2_nonzero == 0:
        print(f"  ✓ H₂(flag complex) = 0 for ALL tournaments at n={n}")
    else:
        print(f"  ✗ H₂(flag complex) ≠ 0 for {h2_nonzero} tournaments")


# ============================================================
# PART 2: Relationship between flag H₂ and path H₂
# ============================================================
print(f"\n{'='*70}")
print("PART 2: FLAG H₂ vs PATH β₂")
print("=" * 70)

# If flag H₂ = 0 always, then every combination of TT triples that forms
# a simplicial cycle is a simplicial boundary. Combined with the NT filling
# from DT_c + cancel, this would prove path β₂ = 0.

# But first: is the simplicial Z₂ (TT cycles) even a subset of path Z₂?
# In path homology, Z₂ ⊂ Ω₂. The TT triples are in Ω₂.
# A simplicial 2-cycle is a combination of TT triples with ∂₂ = 0 in arcs.
# This is the same as ∂₂ = 0 in Ω₁ (since all TT faces are arcs = A₁).
# So simplicial Z₂^simp ⊂ path Z₂.

# And simplicial im(∂₃^simp) ⊂ path im(∂₃) (since transitive 4-sets give DT_f paths).
# So H₂^simp = 0 implies simplicial Z₂ ⊂ path im(∂₃).

# The remaining question: is path Z₂ = simplicial Z₂ + (NT corrections)?
# If so, and if NT corrections are covered by DT_c + cancel, we're done.

print("\nComputing: dim(path Z₂) vs dim(simp Z₂) at n=5")

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


n = 5
m = n*(n-1)//2
z2_comparison = Counter()  # (z2_path, z2_simp)

for bits in range(1 << m):
    A = build_adj(n, bits)

    # Path Z₂
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0:
        z2_path = 0
    else:
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2_path = d2 - rk2

    # Simplicial Z₂
    h2_simp, n_tt, rk_d2_simp = flag_homology_2(A, n)
    z2_simp = n_tt - rk_d2_simp

    z2_comparison[(z2_path, z2_simp)] += 1

print(f"\n(path Z₂, simp Z₂) distribution at n=5:")
for key in sorted(z2_comparison.keys()):
    print(f"  path Z₂={key[0]}, simp Z₂={key[1]}: {z2_comparison[key]}")

# The difference z2_path - z2_simp tells us how many Z₂ directions
# come from NT components (not captured by the simplicial structure)
diff_dist = Counter()
for (zp, zs), count in z2_comparison.items():
    diff_dist[zp - zs] += count
print(f"\npath Z₂ - simp Z₂ distribution:")
for d, count in sorted(diff_dist.items()):
    print(f"  diff={d}: {count}")

print("\nDone.")
