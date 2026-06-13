#!/usr/bin/env python3
"""
β_2 = 0 VIA EXACT SEQUENCE ARGUMENT

The key insight: we can analyze the exact sequence
  0 → Ω_3 → A_3 → J_3 → 0
where J_3 is the "junk space" — the quotient A_3/Ω_3.

For tournaments:
- A_3 = all 3-paths (a,b,c,d) with a→b→c→d
- J_3 encodes the non-allowed face constraints
- Ω_3 = ker(π: A_3 → J_3)

The boundary map ∂_3: A_3 → A_2 restricts to Ω_3 → Ω_2.

Consider the commutative diagram:
  A_3 --∂_3--> A_2
   |             |
   v             v
  J_3          J_2

Here J_2 = A_2/Ω_2 encodes the non-transitive 2-paths.

The question is: does the induced map ∂_3|Ω_3 surject onto ker(∂_2|Ω_2)?

ALTERNATIVE APPROACH: Use the long exact sequence from the short exact
sequence of chain complexes.

Let A_• be the full chain complex and Ω_• the ∂-invariant subcomplex.
We have 0 → Ω_• → A_• → A_•/Ω_• → 0.
This gives a long exact sequence in homology:
  H_3(A) → H_3(A/Ω) → H_2(Ω) → H_2(A) → H_2(A/Ω) → H_1(Ω) → ...

For tournaments, A_p consists of all p-paths with consecutive edges.
A is NOT a chain complex in the usual sense because ∂(A_p) ⊄ A_{p-1}
in general. Actually wait — is it?

Let's check: for a 2-path (a,b,c) with a→b→c:
∂(a,b,c) = (b,c) - (a,c) + (a,b)
All three are 1-paths. But is (a,c) a 1-path? It's in A_1 iff a→c.
If c→a instead, then (a,c) is NOT in A_1 and ∂(a,b,c) ∉ A_1.

So A is NOT a chain complex! The boundary of an A_p element may not
be in A_{p-1}. The whole point of Ω is to restrict to the subspace
where boundaries stay allowed.

Actually wait: ∂ is defined on the full path space, and A_p ⊂ all p-paths.
The boundary map ∂ always outputs formal combinations of (p-1)-paths,
but these may not all be in A_{p-1}. The Ω subspace is defined as
{c ∈ A_p : ∂c ∈ A_{p-1}}, ensuring the chain complex works.

So the correct picture is:
- A_p = {p-paths with all consecutive edges}
- Ω_p = {c ∈ Span(A_p) : ∂c ∈ Span(A_{p-1})}
- ∂: Ω_p → Ω_{p-1} is well-defined (since ∂Ω_p ⊂ A_{p-1} and then
  ∂²=0 ensures ∂c is ∂-invariant too)

The beta_2=0 question is: H_2(Ω_•) = ker(∂_2:Ω_2→Ω_1) / im(∂_3:Ω_3→Ω_2) = 0.
"""
import numpy as np
from collections import Counter
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths, compute_omega_basis

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A_adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A_adj[i][j] = 1
            else: A_adj[j][i] = 1
        yield A_adj

# ===== Dimension counting: |A_p|, dim(Ω_p), |junk types| =====
print("=" * 70)
print("DIMENSION COUNTING FOR TOURNAMENTS")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    dim_stats = Counter()
    dim_detail = {}

    for A_idx, A in enumerate(all_tournaments_gen(n)):
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)

        # Ω_2 = transitive triples
        tt = [p for p in a2 if A[p[0]][p[2]] == 1]

        # Ω_3 via compute_omega_basis
        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        # DT 4-paths
        dt = [p for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1]

        # Non-DT 3-paths with face analysis
        non_dt = [p for p in a3 if not (A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1)]

        # Count 3-cycles
        t3 = 0
        for i in range(n):
            for j in range(i+1,n):
                for k in range(j+1,n):
                    if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                        t3 += 1

        key = (len(a2), len(tt), len(a3), len(dt), dim_om3, t3)
        dim_stats[key] += 1

    print(f"  Statistics (|A_2|, |Ω_2|, |A_3|, DT, dim(Ω_3), t3): count")
    for key in sorted(dim_stats.keys()):
        a2_ct, om2_ct, a3_ct, dt_ct, om3_ct, t3_val = key
        extra = om3_ct - dt_ct
        print(f"  ({a2_ct}, {om2_ct}, {a3_ct}, DT={dt_ct}, Ω_3={om3_ct}, t3={t3_val}): {dim_stats[key]}  "
              f"[extra Ω_3 = {extra}]")

# ===== THE KEY: Does dim(Ω_3) always ≥ dim(ker ∂_2|Ω_2)? =====
print(f"\n\n{'='*70}")
print("SURJECTIVITY CHECK: im(∂_3|Ω_3) ⊇ ker(∂_2|Ω_2)?")
print("="*70)

n = 5
all_pass = True
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    if len(tt) == 0: continue

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # Build ∂_2 on Ω_2
    edge_list = [(i,j) for i in range(n) for j in range(n) if A[i][j] == 1]
    edge_idx = {e: i for i, e in enumerate(edge_list)}
    bd2 = np.zeros((len(edge_list), len(tt)))
    for j, (a, b, c) in enumerate(tt):
        bd2[edge_idx[(b,c)], j] += 1
        bd2[edge_idx[(a,c)], j] -= 1
        bd2[edge_idx[(a,b)], j] += 1

    rank2 = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(tt) - rank2

    if ker_dim == 0: continue

    # Build ∂_3 on Ω_3 → Ω_2
    if dim_om3 == 0:
        im_dim = 0
    else:
        tt_idx = {t: i for i, t in enumerate(tt)}
        a3_arr = [tuple(p) for p in a3]

        # ∂_3 maps Ω_3 basis to A_2 coordinates, then project to Ω_2
        bd3_full = np.zeros((len(a2), len(a3_arr)))
        a2_idx = {tuple(p): i for i, p in enumerate(a2)}
        for j, path in enumerate(a3_arr):
            v0, v1, v2, v3 = path
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in a2_idx:
                    bd3_full[a2_idx[face], j] += sign

        # Restrict to Ω_3 input
        bd3_omega = bd3_full @ om3

        # Project to Ω_2 output
        # The columns of bd3_omega are in A_2 coordinates
        # Project to the transitive triple subspace
        tt_proj = np.zeros((len(tt), len(a2)))
        for i, t in enumerate(tt):
            tt_proj[i, a2_idx[t]] = 1.0
        bd3_tt = tt_proj @ bd3_omega

        im_dim = np.linalg.matrix_rank(bd3_tt, tol=1e-8)

    beta_2 = ker_dim - im_dim
    if beta_2 != 0:
        all_pass = False
        print(f"  FAIL: ker={ker_dim}, im={im_dim}, β_2={beta_2}")

if all_pass:
    print(f"  ALL n={n} tournaments: ker(∂_2|Ω_2) = im(∂_3|Ω_3) ✓")

# ===== The relationship between dim(Ω_3) and ker(∂_2|Ω_2) =====
print(f"\n\n{'='*70}")
print("DIMENSION RELATIONSHIP: dim(Ω_3) vs ker(∂_2|Ω_2)")
print("="*70)

n = 5
data = []
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    if len(tt) == 0: continue

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    edge_list = [(i,j) for i in range(n) for j in range(n) if A[i][j] == 1]
    edge_idx = {e: i for i, e in enumerate(edge_list)}
    bd2 = np.zeros((len(edge_list), len(tt)))
    for j, (a, b, c) in enumerate(tt):
        bd2[edge_idx[(b,c)], j] += 1
        bd2[edge_idx[(a,c)], j] -= 1
        bd2[edge_idx[(a,b)], j] += 1

    rank2 = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(tt) - rank2

    dt = sum(1 for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1)
    data.append((ker_dim, dim_om3, dt, len(tt), len(a3)))

# Summarize
print(f"  n={n}: {len(data)} tournaments with Ω_2 ≠ 0")
pairs = Counter()
for ker, om3, dt, tt, a3 in data:
    pairs[(ker, om3, dt)] += 1

print(f"\n  (ker(∂_2), dim(Ω_3), DT_count): tournaments")
for key in sorted(pairs.keys()):
    ker, om3, dt = key
    print(f"    ({ker}, {om3}, {dt}): {pairs[key]}")

# Critical question: is dim(Ω_3) always ≥ ker(∂_2)?
always_ge = all(om3 >= ker for ker, om3, _, _, _ in data)
print(f"\n  dim(Ω_3) ≥ ker(∂_2)? {always_ge}")

# Is dim(Ω_3) - ker(∂_2) = dim(im(∂_3|Ω_2)) - boundary from Ω_4?
# That is, β_3 analysis
print(f"\n  For β_3 analysis, need dim(ker(∂_3|Ω_3)) - dim(im(∂_4|Ω_4))")

print("\nDone.")
