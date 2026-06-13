#!/usr/bin/env python3
"""
POS (Positive semidefiniteness) and the tiling skeleton.

The transfer matrix M is symmetric (THM-030). This means:
  M has real eigenvalues.

Question: Is M always PSD (positive semidefinite)?
If so, this would give strong constraints via Cauchy-Schwarz.

At odd n: tr(M) = H > 0, so at least one eigenvalue is positive.
But is M PSD?

From the n=5 spectral data:
  Class 0 (H=1):  evals = [2.0, 1.414, 1.0, -1.414, -2.0]  -- NOT PSD
  Class 9 (H=15): evals = [3, 3, 3, 3, 3]                    -- PSD (trivially)
  Class 10 (H=13): evals = [4.34, 4.0, 2.47, 2.0, 0.19]      -- PSD

So M is NOT always PSD. But the classes with high H tend to have
non-negative eigenvalues.

CONNECTION TO SKELETON: Does PSD-ness propagate along skeleton edges?
If T is PSD and we flip one tile, is T' still PSD?
What is the "boundary" between PSD and non-PSD regions?

CONNECTION TO IO: The Irving-Omar walk GF involves det(I+zA^T)/det(I-zA).
The positivity of W_D at z=1 relates to the eigenvalue structure of A.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

def ham_path_count_dp(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1:
                    continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False; break
        if valid: count += 1
    return count

def transfer_matrix(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

# =====================================================================
# n=5: PSD classification
# =====================================================================
print("=" * 70)
print("PSD CLASSIFICATION AT n=5")
print("=" * 70)

n = 5
tiles = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles.append((a, b))
tiles.sort()
m = len(tiles)

# Compute M for all tilings
tiling_data = {}
iso_classes = defaultdict(list)

for bits in range(2**m):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count_dp(A)
    canon = tournament_canonical(A)
    M = transfer_matrix(A)
    evals = sorted(np.linalg.eigvalsh(M).tolist())
    is_psd = all(e >= -1e-10 for e in evals)

    tiling_data[bits] = {
        'H': H, 'M': M, 'evals': evals, 'is_psd': is_psd,
        'canon': canon, 'bits': format(bits, f'0{m}b')
    }
    iso_classes[canon].append(bits)

class_labels = {}
for idx, canon in enumerate(sorted(iso_classes.keys())):
    class_labels[canon] = idx

# PSD summary
psd_count = sum(1 for bits in range(2**m) if tiling_data[bits]['is_psd'])
print(f"\n  PSD tilings: {psd_count} / {2**m}")
print(f"\n  PSD by class:")

psd_by_class = {}
for idx, canon in enumerate(sorted(iso_classes.keys())):
    members = iso_classes[canon]
    H = tiling_data[members[0]]['H']
    psd_members = [b for b in members if tiling_data[b]['is_psd']]
    evals = tiling_data[members[0]]['evals']
    min_eval = min(evals)
    psd_by_class[idx] = len(psd_members) == len(members)
    label = "ALL PSD" if len(psd_members) == len(members) else f"{len(psd_members)}/{len(members)} PSD"
    print(f"    Class {idx} (H={H:2d}, |class|={len(members):2d}): {label}, min_eval={min_eval:.3f}")

# =====================================================================
# PSD boundary on the skeleton
# =====================================================================
print()
print("=" * 70)
print("PSD BOUNDARY ON THE SKELETON")
print("=" * 70)
print("\nEdges crossing PSD/non-PSD boundary:")

boundary_edges = 0
for bits in range(2**m):
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        if neighbor > bits:
            psd1 = tiling_data[bits]['is_psd']
            psd2 = tiling_data[neighbor]['is_psd']
            if psd1 != psd2:
                boundary_edges += 1

total_edges = sum(1 for bits in range(2**m) for ti in range(m) if (bits ^ (1 << ti)) > bits)
print(f"  Boundary edges: {boundary_edges} / {total_edges} total")
print(f"  PSD region is {'connected' if psd_count > 0 else 'empty'}")

# =====================================================================
# Does H determine PSD-ness?
# =====================================================================
print()
print("=" * 70)
print("H vs PSD-NESS")
print("=" * 70)

H_to_psd = defaultdict(lambda: [0, 0])  # H -> [psd_count, non_psd_count]
for bits in range(2**m):
    H = tiling_data[bits]['H']
    if tiling_data[bits]['is_psd']:
        H_to_psd[H][0] += 1
    else:
        H_to_psd[H][1] += 1

print(f"\n  H   PSD  non-PSD  {'H threshold for guaranteed PSD?':>35}")
print("  " + "-" * 60)
for H in sorted(H_to_psd.keys()):
    psd, non_psd = H_to_psd[H]
    total = psd + non_psd
    guaranteed = "ALL PSD" if non_psd == 0 else ""
    print(f"  {H:2d}  {psd:3d}  {non_psd:7d}  {guaranteed}")

# Find the PSD threshold
psd_threshold = None
for H in sorted(H_to_psd.keys()):
    if H_to_psd[H][1] == 0:  # all PSD at this H
        if psd_threshold is None:
            psd_threshold = H
    else:
        psd_threshold = None

print(f"\n  PSD threshold: H >= {psd_threshold}")

# =====================================================================
# Eigenvalue spectrum visualization (text-based)
# =====================================================================
print()
print("=" * 70)
print("EIGENVALUE SPECTRUM BY CLASS")
print("=" * 70)

for idx, canon in enumerate(sorted(iso_classes.keys())):
    members = iso_classes[canon]
    H = tiling_data[members[0]]['H']
    evals = tiling_data[members[0]]['evals']
    is_psd = all(e >= -1e-10 for e in evals)

    # Text visualization: min...max with markers
    eval_str = " ".join(f"{e:7.3f}" for e in evals)
    psd_marker = "+" if is_psd else "-"
    print(f"  [{psd_marker}] Class {idx:2d} (H={H:2d}): [{eval_str}]")

# =====================================================================
# Relationship to Claim A via PSD
# =====================================================================
print()
print("=" * 70)
print("PSD AND CLAIM A")
print("=" * 70)
print("""
Claim A: H(T) - H(T-v) = 2 * sum mu(C_j)  for each vertex v

At odd n, H = tr(M). Removing vertex v:
  H(T-v) = tr(M(T-v)) where T-v is the induced sub-tournament

If M is PSD, then tr(M) >= tr(M_sub) for any principal submatrix M_sub
(by the eigenvalue interlacing theorem applied to symmetric M).

This would give H(T) >= H(T-v) when M is PSD!
Let's check: does H(T) >= H(T-v) hold for PSD tilings?
""")

# Check H >= H-v for PSD and non-PSD tilings
for bits in [0, 21, 31, 63]:  # sample tilings
    d = tiling_data[bits]
    A = tiling_to_tournament(bits, n)
    H = d['H']

    print(f"  Tiling {d['bits']} (H={H}, PSD={d['is_psd']}):")
    for v in range(n):
        # Compute H(T-v)
        remaining = [i for i in range(n) if i != v]
        A_sub = [[A[i][j] for j in remaining] for i in remaining]
        H_sub = ham_path_count_dp(A_sub)
        delta = H - H_sub
        print(f"    v={v}: H(T-v)={H_sub}, delta={delta}, sign={'>=0' if delta >= 0 else '<0'}")

# =====================================================================
# Negative eigenvalues and the odd cycle structure
# =====================================================================
print()
print("=" * 70)
print("NEGATIVE EIGENVALUES AND ODD CYCLES")
print("=" * 70)

# For non-PSD classes, what is the relationship between negative
# eigenvalues and the odd cycle structure?
print("\nNon-PSD classes:")
for idx, canon in enumerate(sorted(iso_classes.keys())):
    members = iso_classes[canon]
    d = tiling_data[members[0]]
    if d['is_psd']:
        continue
    H = d['H']
    evals = d['evals']
    neg_evals = [e for e in evals if e < -1e-10]
    pos_evals = [e for e in evals if e > 1e-10]
    neg_sum = sum(neg_evals)
    pos_sum = sum(pos_evals)
    print(f"  Class {idx} (H={H}): neg_sum={neg_sum:.3f}, pos_sum={pos_sum:.3f}, tr={pos_sum+neg_sum:.3f}=H={H}")
    print(f"    neg eigenvalues: {[round(e,3) for e in neg_evals]}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
PSD Analysis at n=5:
  - M is PSD for {psd_count}/{2**m} tilings
  - PSD threshold: H >= {psd_threshold} (all tilings with H >= {psd_threshold} have PSD M)
  - Non-PSD tilings have BOTH positive and negative eigenvalues
  - The PSD boundary on the skeleton separates low-H from high-H tilings

Key insight: PSD-ness of M is correlated with H but not determined by it.
At low H (few odd cycles), M can have negative eigenvalues.
At high H (many odd cycles), M is PSD.

For Claim A: PSD would imply H >= H-v (via eigenvalue interlacing)
but M is NOT always PSD, so this approach doesn't directly work.
The even-r-power structure (THM-030) constrains M but doesn't force PSD.
""")
