#!/usr/bin/env python3
"""Search for formulas for dim(Ω_k) in tournaments.

At n=5, dim(Ω_2) ranges from 8 to 10.
The dim(Ω_2) = |TT| + gap formula holds.

Question: Is there a formula for dim(Ω_2) in terms of t3, t4, t5?
And for dim(Ω_3)?

Also: compute dim(Ω_k) for all k at each tournament type,
looking for a pattern.
"""
import numpy as np
from itertools import combinations
import sys, time
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

def count_transitive_k(A, n, k):
    """Count transitive k-subsets (all C(k,2) edges form transitive tournament)."""
    count = 0
    for subset in combinations(range(n), k):
        s = list(subset)
        # Check: is the induced subtournament transitive (acyclic)?
        # Transitive iff there are 0 directed 3-cycles
        is_trans = True
        for a, b, c in combinations(s, 3):
            if (A[a][b] and A[b][c] and A[c][a]) or \
               (A[b][a] and A[a][c] and A[c][b]):
                is_trans = False
                break
        if is_trans:
            count += 1
    return count

# ===== n=5: Full chain complex dimensions =====
print("=" * 70)
print("CHAIN COMPLEX DIMENSIONS vs TOURNAMENT INVARIANTS (n=5)")
print("=" * 70)

n = 5
data = defaultdict(list)
for A in all_tournaments_gen(n):
    t3 = sum(1 for a, b, c in combinations(range(n), 3)
             if (A[a][b] and A[b][c] and A[c][a]) or
                (A[b][a] and A[a][c] and A[c][b]))
    sc = score_seq(A, n)

    # TT count (transitive triples)
    tt = sum(1 for a, b, c in combinations(range(n), 3)
             if not ((A[a][b] and A[b][c] and A[c][a]) or
                     (A[b][a] and A[a][c] and A[c][b])))
    # = C(n,3) - t3

    # DT count (doubly-transitive 4-paths)
    t4 = count_transitive_k(A, n, 4)  # transitive 4-subsets

    dims = {}
    for d in range(n):
        ap_d = enumerate_allowed_paths(A, n, d)
        if d == 0:
            dims[d] = n
        elif d == 1:
            dims[d] = len(ap_d)  # = C(n,2) for tournaments
        else:
            ap_dm1 = enumerate_allowed_paths(A, n, d-1)
            om = compute_omega_basis(A, n, d, ap_d, ap_dm1)
            dims[d] = om.shape[1] if om.ndim == 2 else 0

    key = (t3, sc, tt, t4)
    data[key].append(tuple(dims[d] for d in range(n)))

print(f"(t3, score, |TT|, t4) → Ω dims (distinct)")
for key in sorted(data):
    t3, sc, tt, t4 = key
    dim_set = set(data[key])
    for dims in sorted(dim_set):
        count = data[key].count(dims)
        chi = sum((-1)**i * dims[i] for i in range(n))
        print(f"  t3={t3}, sc={sc}, TT={tt}, t4={t4}: Ω={list(dims)}, χ={chi} ({count}x)")

# ===== Check: is dim(Ω_2) = TT + gap, and is gap determined by t3? =====
print(f"\n\n{'='*70}")
print("dim(Ω_2) vs |TT| and gap")
print("="*70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    gap_by_t3 = defaultdict(list)
    for A in all_tournaments_gen(n):
        t3 = sum(1 for a, b, c in combinations(range(n), 3)
                 if (A[a][b] and A[b][c] and A[c][a]) or
                    (A[b][a] and A[a][c] and A[c][b]))
        tt = sum(1 for a, b, c in combinations(range(n), 3)
                 if not ((A[a][b] and A[b][c] and A[c][a]) or
                         (A[b][a] and A[a][c] and A[c][b])))

        a2 = enumerate_allowed_paths(A, n, 2)
        a1 = enumerate_allowed_paths(A, n, 1)
        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
        gap = dim_om2 - tt
        gap_by_t3[t3].append(gap)

    for t3 in sorted(gap_by_t3):
        gaps = Counter(gap_by_t3[t3])
        tt = n*(n-1)*(n-2)//6 - t3  # C(n,3) - t3
        print(f"  t3={t3} (TT={tt}): gap values = {dict(gaps)}")

# ===== At n=5, look for formula for dim(Ω_3) =====
print(f"\n\n{'='*70}")
print("dim(Ω_3) vs tournament invariants (n=5)")
print("="*70)

n = 5
dt_by_type = defaultdict(list)
for A in all_tournaments_gen(n):
    t3 = sum(1 for a, b, c in combinations(range(n), 3)
             if (A[a][b] and A[b][c] and A[c][a]) or
                (A[b][a] and A[a][c] and A[c][b]))
    sc = score_seq(A, n)

    # Count DT paths directly
    a3 = enumerate_allowed_paths(A, n, 3)
    a2 = enumerate_allowed_paths(A, n, 2)
    a2_set = set(tuple(p) for p in a2)
    dt_count = 0
    for p in a3:
        a, b, c, d = tuple(p)
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        if all(f in a2_set for f in faces):
            dt_count += 1

    a1 = enumerate_allowed_paths(A, n, 1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    dt_by_type[(t3, sc)].append((dt_count, dim_om3))

for key in sorted(dt_by_type):
    t3, sc = key
    vals = Counter(dt_by_type[key])
    for (dt, om3), cnt in sorted(vals.items()):
        print(f"  t3={t3}, sc={sc}: |DT|={dt}, Ω_3={om3} ({cnt}x)")

print("\nDone.")
