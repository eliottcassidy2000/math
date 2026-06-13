#!/usr/bin/env python3
"""
beta2_rank_structure.py - Analyze rank structure forcing beta_2 = 0

For each tournament T, compute:
  dim(Omega_p), rank(d_p), ker(d_p), im(d_{p+1}), beta_p
for all p.

Goal: Find explicit formulas for rank(d_2|Omega_2) and rank(d_3|Omega_3)
that explain why ker(d_2) = im(d_3) always.

Key identity: beta_2 = 0 iff ker(d_2|Omega_2) = im(d_3|Omega_3)
i.e., rank(d_3|Omega_3) = dim(Omega_2) - rank(d_2|Omega_2)

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def full_chain_analysis(A, n, max_p=None):
    """Complete chain complex analysis: dims, ranks, kernels, images, bettis."""
    if max_p is None:
        max_p = n - 1

    allowed = {}
    for p in range(-1, max_p + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    omega = {}
    for p in range(max_p + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    results = {}
    for p in range(max_p + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0

        # rank of d_p: Omega_p -> A_{p-1}
        if dim_omega_p > 0 and len(allowed[p-1]) > 0:
            D_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
            D_p_omega = D_p @ omega[p]
            S = np.linalg.svd(D_p_omega, compute_uv=False)
            rank_dp = sum(s > 1e-8 for s in S)
        else:
            rank_dp = 0

        ker_dp = dim_omega_p - rank_dp

        # im(d_{p+1}): rank of d_{p+1} restricted to Omega_{p+1}
        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0 and len(allowed[p]) > 0:
            D_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            D_p1_omega = D_p1 @ omega[p+1]
            S1 = np.linalg.svd(D_p1_omega, compute_uv=False)
            im_dp1 = sum(s > 1e-8 for s in S1)
        else:
            im_dp1 = 0

        beta_p = ker_dp - im_dp1

        results[p] = {
            'dim_omega': dim_omega_p,
            'rank_d': rank_dp,
            'ker_d': ker_dp,
            'im_next': im_dp1,
            'beta': beta_p,
        }

    return results


def count_cycles(A, n, k):
    """Count directed k-cycles (vertex sets with at least one directed k-cycle)."""
    if k < 3:
        return 0
    count = 0
    for verts in combinations(range(n), k):
        # Check if any directed cycle exists on these k vertices
        vlist = list(verts)
        from itertools import permutations
        found = False
        for perm in permutations(vlist):
            is_cycle = True
            for i in range(k):
                if A[perm[i]][perm[(i+1)%k]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                found = True
                break
        if found:
            count += 1
    return count


# ============================================================
# PART 1: Complete rank structure at n=4,5
# ============================================================
print("=" * 70)
print("COMPLETE RANK STRUCTURE BY SCORE TYPE")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m

    score_data = defaultdict(list)

    for bits in range(total):
        A = build_adj(n, bits)
        scores = tuple(sorted([sum(row) for row in A]))
        res = full_chain_analysis(A, n, max_p=n-1)
        c3 = count_cycles(A, n, 3)

        entry = {'bits': bits, 'c3': c3}
        for p in range(n):
            entry[f'dim_{p}'] = res[p]['dim_omega']
            entry[f'rank_{p}'] = res[p]['rank_d']
            entry[f'ker_{p}'] = res[p]['ker_d']
            entry[f'im_{p}'] = res[p]['im_next']
            entry[f'beta_{p}'] = res[p]['beta']
        score_data[scores].append(entry)

    print(f"\nn={n}: {total} tournaments")
    for scores in sorted(score_data.keys()):
        items = score_data[scores]
        print(f"\n  Score {scores} ({len(items)} tournaments):")

        for p in range(min(n, 5)):
            dims = [x[f'dim_{p}'] for x in items]
            ranks = [x[f'rank_{p}'] for x in items]
            kers = [x[f'ker_{p}'] for x in items]
            ims = [x[f'im_{p}'] for x in items]
            betas = [x[f'beta_{p}'] for x in items]

            dim_vals = sorted(set(dims))
            rank_vals = sorted(set(ranks))
            ker_vals = sorted(set(kers))
            im_vals = sorted(set(ims))

            # Compact display
            if len(dim_vals) == 1:
                print(f"    p={p}: dim(Omega)={dim_vals[0]}, rank(d)={rank_vals}, ker={ker_vals}, im(d+1)={im_vals}, beta={sorted(set(betas))}")
            else:
                print(f"    p={p}: dim(Omega)={dim_vals}, rank(d)={rank_vals}, ker={ker_vals}, im(d+1)={im_vals}, beta={sorted(set(betas))}")


# ============================================================
# PART 2: ker(d2) and im(d3) — what determines them?
# ============================================================
print(f"\n{'='*70}")
print("ker(d2) vs im(d3) BY TOURNAMENT INVARIANTS")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# Collect (c3, c5, dim_Omega2, ker_d2, im_d3, rank_d2, rank_d3)
data_rows = []
for bits in range(total):
    A = build_adj(n, bits)
    c3 = count_cycles(A, n, 3)
    res = full_chain_analysis(A, n, max_p=4)

    row = {
        'bits': bits,
        'c3': c3,
        'scores': tuple(sorted([sum(row) for row in A])),
    }
    for p in range(5):
        row[f'dim_{p}'] = res[p]['dim_omega']
        row[f'rank_{p}'] = res[p]['rank_d']
        row[f'ker_{p}'] = res[p]['ker_d']
        row[f'im_{p}'] = res[p]['im_next']
    data_rows.append(row)

# Group by c3
c3_groups = defaultdict(list)
for row in data_rows:
    c3_groups[row['c3']].append(row)

print(f"\nn=5: ker(d2) and im(d3) grouped by c3:")
for c3 in sorted(c3_groups.keys()):
    items = c3_groups[c3]
    ker2_vals = sorted(set(x['ker_2'] for x in items))
    im3_vals = sorted(set(x['im_2'] for x in items))
    dim2_vals = sorted(set(x['dim_2'] for x in items))
    rank2_vals = sorted(set(x['rank_2'] for x in items))
    dim3_vals = sorted(set(x['dim_3'] for x in items))
    rank3_vals = sorted(set(x['rank_3'] for x in items))

    print(f"\n  c3={c3}: {len(items)} tournaments")
    print(f"    dim(Omega_2)={dim2_vals}, rank(d2)={rank2_vals}, ker(d2)={ker2_vals}")
    print(f"    dim(Omega_3)={dim3_vals}, rank(d3)={rank3_vals}=im(d3)={im3_vals}")
    print(f"    ker(d2)==im(d3)? {all(x['ker_2']==x['im_2'] for x in items)}")


# ============================================================
# PART 3: Test explicit formulas
# ============================================================
print(f"\n{'='*70}")
print("FORMULA TESTS")
print("=" * 70)

# Hypothesis: ker(d2) = c3? Or ker(d2) = f(c3)?
print("\nn=5: ker(d2) vs c3:")
c3_ker2 = defaultdict(list)
for row in data_rows:
    c3_ker2[row['c3']].append(row['ker_2'])

for c3 in sorted(c3_ker2.keys()):
    vals = sorted(set(c3_ker2[c3]))
    print(f"  c3={c3}: ker(d2) in {vals}")

# Check: ker(d2) = c3?
match = all(row['ker_2'] == row['c3'] for row in data_rows)
print(f"\n  ker(d2) == c3 for all? {match}")

# Check at n=4 too
n = 4
data_4 = []
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    c3 = count_cycles(A, n, 3)
    res = full_chain_analysis(A, n, max_p=3)
    data_4.append({
        'c3': c3,
        'ker_2': res[2]['ker_d'],
        'im_2': res[2]['im_next'],
        'dim_2': res[2]['dim_omega'],
    })

print(f"\nn=4:")
for row in data_4:
    pass
c3_ker2_4 = defaultdict(list)
for row in data_4:
    c3_ker2_4[row['c3']].append(row['ker_2'])
for c3 in sorted(c3_ker2_4.keys()):
    vals = sorted(set(c3_ker2_4[c3]))
    print(f"  c3={c3}: ker(d2) in {vals}")

match4 = all(row['ker_2'] == row['c3'] for row in data_4)
print(f"  ker(d2) == c3 for all n=4? {match4}")


# ============================================================
# PART 4: Check at n=6 (exhaustive, may be slow)
# ============================================================
print(f"\n{'='*70}")
print("n=6: EXHAUSTIVE ker(d2) vs c3 CHECK")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

c3_ker2_6 = defaultdict(set)
all_match = True
counter = 0

for bits in range(total):
    A = build_adj(n, bits)
    c3 = count_cycles(A, n, 3)
    res = full_chain_analysis(A, n, max_p=4)
    ker2 = res[2]['ker_d']
    c3_ker2_6[c3].add(ker2)

    if ker2 != c3:
        all_match = False
        if counter < 3:
            print(f"  MISMATCH at T#{bits}: c3={c3}, ker(d2)={ker2}")
        counter += 1

    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s), match so far: {all_match}")

elapsed = time.time() - t0
print(f"\nn=6: ({elapsed:.0f}s)")
print(f"  ker(d2) == c3 for ALL? {all_match}")
if not all_match:
    print(f"  Mismatches: {counter}")
    for c3 in sorted(c3_ker2_6.keys()):
        print(f"  c3={c3}: ker(d2) in {sorted(c3_ker2_6[c3])}")


# ============================================================
# PART 5: If ker(d2) = c3, then beta_2 = 0 iff im(d3) = c3
# What determines im(d3)?
# ============================================================
print(f"\n{'='*70}")
print("im(d3) ANALYSIS")
print("=" * 70)

n = 5
print(f"\nn=5: im(d3) vs c3:")
c3_im3 = defaultdict(list)
for row in data_rows:
    c3_im3[row['c3']].append(row['im_2'])

for c3 in sorted(c3_im3.keys()):
    vals = sorted(set(c3_im3[c3]))
    print(f"  c3={c3}: im(d3) in {vals}")

match_im = all(row['im_2'] == row['c3'] for row in data_rows)
print(f"  im(d3) == c3 for all? {match_im}")

# If both ker(d2) = c3 and im(d3) = c3, then beta_2 = 0. QED.
print(f"\n  If ker(d2) = im(d3) = c3 universally, then beta_2 = 0 follows!")


# ============================================================
# PART 6: Detailed rank structure — where do the ranks come from?
# ============================================================
print(f"\n{'='*70}")
print("RANK FORMULA: rank(d2) and rank(d3)")
print("=" * 70)

n = 5
print(f"\nn=5: Rank(d2) = dim(Omega_2) - ker(d2) = dim(Omega_2) - c3")
for row in data_rows[:20]:
    rank2 = row['dim_2'] - row['ker_2']
    print(f"  T#{row['bits']}: dim(O2)={row['dim_2']}, c3={row['c3']}, rank(d2)={rank2}, check={row['rank_2']}")

# Also look at rank(d3) = im(d3) for Omega_3
print(f"\nn=5: rank(d3|Omega_3) = im(d3) = c3?")
for row in data_rows[:20]:
    rank3 = row['rank_3']
    dim3 = row['dim_3']
    im3 = row['im_2']
    print(f"  T#{row['bits']}: dim(O3)={dim3}, rank(d3)={rank3}, im(d3)={im3}, c3={row['c3']}")


print("\n\nDone.")
