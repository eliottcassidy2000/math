#!/usr/bin/env python3
"""
beta2_omega2_analysis.py - Direct analysis of Omega_2 and H_2 for tournaments

KEY QUESTION: Why is beta_2(T) = 0 for all tournaments?

APPROACH: Analyze the chain complex directly.
- Omega_2 condition: for each non-allowed pair (a,c) with c->a,
  sum_{b: a->b->c} z_{(a,b,c)} = 0
- ker(d_2) = Omega_2-cycles
- im(d_3) = boundaries from Omega_3
- beta_2 = dim(ker d_2) - dim(im d_3)

We compute and display the EXPLICIT structure of Omega_2 and its boundary.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from itertools import permutations, combinations
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


def analyze_omega2(A, n):
    """Detailed analysis of Omega_2 structure."""
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    # Non-allowed pairs: (a,c) with c->a
    non_allowed = []
    for a in range(n):
        for c in range(n):
            if a != c and A[c][a] == 1:
                non_allowed.append((a,c))

    # For each non-allowed pair, find mediators
    constraints = {}
    for (a,c) in non_allowed:
        mediators = [b for b in range(n) if b != a and b != c
                     and A[a][b] == 1 and A[b][c] == 1]
        constraints[(a,c)] = mediators

    # Build Omega_2 constraint matrix
    path2_idx = {p: i for i, p in enumerate(paths2)}
    C = np.zeros((len(non_allowed), len(paths2)))
    for row, (a,c) in enumerate(non_allowed):
        for b in constraints[(a,c)]:
            if (a,b,c) in path2_idx:
                C[row, path2_idx[(a,b,c)]] = 1

    # Omega_2 = ker(C)
    if C.shape[0] > 0 and C.shape[1] > 0:
        U, S, Vt = np.linalg.svd(C, full_matrices=True)
        rank_C = sum(s > 1e-10 for s in S)
        omega2_basis = Vt[rank_C:].T
    else:
        omega2_basis = np.eye(len(paths2))
        rank_C = 0

    dim_omega2 = omega2_basis.shape[1] if omega2_basis.ndim == 2 else 0

    # Boundary d2: A_2 -> A_1
    path1_idx = {p: i for i, p in enumerate(paths1)}
    D2 = np.zeros((len(paths1), len(paths2)))
    for j, (a,b,c_) in enumerate(paths2):
        if (b,c_) in path1_idx: D2[path1_idx[(b,c_)], j] += 1
        if (a,c_) in path1_idx: D2[path1_idx[(a,c_)], j] -= 1
        if (a,b) in path1_idx: D2[path1_idx[(a,b)], j] += 1

    # d2 restricted to Omega_2
    D2_omega = D2 @ omega2_basis if dim_omega2 > 0 else np.zeros((len(paths1), 0))

    if D2_omega.shape[1] > 0:
        S2 = np.linalg.svd(D2_omega, compute_uv=False)
        rank_d2 = sum(s > 1e-8 for s in S2)
    else:
        rank_d2 = 0
    ker_d2 = dim_omega2 - rank_d2

    # Omega_3 and d3
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_omega3 = omega3.shape[1] if omega3.ndim == 2 and omega3.shape[1] > 0 else 0

    D3 = build_full_boundary_matrix(paths3, paths2)
    if dim_omega3 > 0:
        D3_omega = D3 @ omega3
        S3 = np.linalg.svd(D3_omega, compute_uv=False)
        im_d3 = sum(s > 1e-8 for s in S3)
    else:
        im_d3 = 0

    beta2 = ker_d2 - im_d3

    return {
        'paths2': paths2, 'paths1': paths1,
        'non_allowed': non_allowed, 'constraints': constraints,
        'dim_A2': len(paths2), 'dim_omega2': dim_omega2,
        'num_constraints': len(non_allowed), 'rank_constraints': rank_C,
        'rank_d2': rank_d2, 'ker_d2': ker_d2,
        'dim_omega3': dim_omega3, 'im_d3': im_d3,
        'beta2': beta2,
    }


# ============================================================
# PART 1: Omega_2 dimensions by score type at n=5
# ============================================================
print("=" * 70)
print("OMEGA_2 DIMENSIONS BY SCORE TYPE (n=5)")
print("=" * 70)

n = 5
score_data = defaultdict(list)
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    info = analyze_omega2(A, n)
    score_data[scores].append(info)

for scores in sorted(score_data.keys()):
    items = score_data[scores]
    omega2_vals = [x['dim_omega2'] for x in items]
    ker_vals = [x['ker_d2'] for x in items]
    im_vals = [x['im_d3'] for x in items]
    beta2_vals = [x['beta2'] for x in items]
    print(f"\nScore {scores}: {len(items)} tournaments")
    print(f"  dim(A_2) = {items[0]['dim_A2']}")
    print(f"  dim(Omega_2): {dict(Counter(omega2_vals))}")
    print(f"  ker(d_2): {dict(Counter(ker_vals))}")
    print(f"  im(d_3): {dict(Counter(im_vals))}")
    print(f"  beta_2: {dict(Counter(beta2_vals))}")


# ============================================================
# PART 2: Exhaustive check at n=4,5,6 - ker(d2)=im(d3)?
# ============================================================
print(f"\n{'='*70}")
print("KEY: ker(d_2) = im(d_3) ALWAYS? (beta_2 = 0)")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m
    all_zero = True
    t0 = time.time()

    for bits in range(total):
        A = build_adj(n, bits)
        info = analyze_omega2(A, n)
        if info['beta2'] != 0:
            all_zero = False
            print(f"  COUNTEREXAMPLE at n={n}: T#{bits}, beta_2={info['beta2']}")
            break

    elapsed = time.time() - t0
    print(f"n={n}: beta_2 = 0 for ALL {total} tournaments ({elapsed:.0f}s)")


# n=6 exhaustive
n = 6
m = n*(n-1)//2
total = 1 << m
all_zero = True
t0 = time.time()

for bits in range(total):
    A = build_adj(n, bits)
    info = analyze_omega2(A, n)
    if info['beta2'] != 0:
        all_zero = False
        print(f"  COUNTEREXAMPLE at n={n}: T#{bits}, beta_2={info['beta2']}")
        break

    if (bits+1) % 10000 == 0:
        elapsed = time.time() - t0
        print(f"  n={n}: {bits+1}/{total} ({elapsed:.0f}s)")

elapsed = time.time() - t0
if all_zero:
    print(f"n={n}: beta_2 = 0 for ALL {total} tournaments ({elapsed:.0f}s)")


# ============================================================
# PART 3: FORMULA ANALYSIS
# ============================================================
print(f"\n{'='*70}")
print("FORMULA: dim(Omega_2) = dim(A_2) - rank(C)")
print("=" * 70)

# For a tournament T:
# dim(A_2) = sum_{a->b} |N+(b) \ {a}|
# Each non-allowed pair (a,c) with c->a gives constraint
# Mediators(a,c) = N+(a) intersect N-(c) \ {a,c}

n = 5
print(f"\nn={n}: Regular tournaments (score (2,2,2,2,2)):")
reg_count = 0
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = sorted([sum(row) for row in A])
    if scores != [2,2,2,2,2]:
        continue
    reg_count += 1
    if reg_count > 3:
        continue

    info = analyze_omega2(A, n)
    print(f"\n  T#{bits}:")
    print(f"    dim(A_2) = {info['dim_A2']}")
    print(f"    Non-allowed pairs: {info['num_constraints']}")
    print(f"    rank(C) = {info['rank_constraints']}")
    print(f"    dim(Omega_2) = {info['dim_omega2']}")

    # Show mediator structure
    med_counts = Counter(len(info['constraints'][ac]) for ac in info['non_allowed'])
    print(f"    Mediator count distribution: {dict(sorted(med_counts.items()))}")

    # Show explicit constraints
    for (a,c) in info['non_allowed']:
        meds = info['constraints'][(a,c)]
        med_str = ','.join(str(b) for b in meds)
        print(f"      ({a},{c}) c->a: mediators=[{med_str}] ({len(meds)} meds)")


# ============================================================
# PART 4: dim(Omega_2) formula test
# ============================================================
print(f"\n{'='*70}")
print("dim(Omega_2) vs TRANSITIVE TRIPLES and 3-CYCLES")
print("=" * 70)

n = 5
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))

    # Count transitive triples (a,b,c with a->b, b->c, a->c)
    trans3 = 0
    cyc3 = 0
    for triple in combinations(range(n), 3):
        i,j,k = triple
        # Check all 6 orderings
        for a,b,c in [(i,j,k),(i,k,j),(j,i,k),(j,k,i),(k,i,j),(k,j,i)]:
            if A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1:
                trans3 += 1
                break  # count each triple once
            if A[a][b] == 1 and A[b][c] == 1 and A[c][a] == 1:
                cyc3 += 1
                break

    info = analyze_omega2(A, n)
    if bits < 5 or scores == (2,2,2,2,2):
        if bits < 5 or (scores == (2,2,2,2,2) and bits == 76):
            print(f"  T#{bits} scores={scores}: trans3={trans3}, cyc3={cyc3}, dim(Omega_2)={info['dim_omega2']}, ker={info['ker_d2']}, im={info['im_d3']}")


# Actually, let's collect systematically
print(f"\nSystematic: dim(Omega_2) as function of c3 (directed 3-cycles):")
c3_to_omega = defaultdict(list)
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    c3 = 0
    for triple in combinations(range(n), 3):
        i,j,k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
            c3 += 1
    info = analyze_omega2(A, n)
    c3_to_omega[c3].append(info['dim_omega2'])

for c3 in sorted(c3_to_omega.keys()):
    vals = c3_to_omega[c3]
    print(f"  c3={c3}: dim(Omega_2) in {sorted(set(vals))}, count={len(vals)}")

# Now check: is dim(Omega_2) = dim(A_2) - c3?
# dim(A_2) at n=5 is always 20 (for regular d+=2) or varies?
print(f"\nCheck dim(Omega_2) = dim(A_2) - c3:")
for bits in range(min(20, 1 << (n*(n-1)//2))):
    A = build_adj(n, bits)
    c3 = 0
    for triple in combinations(range(n), 3):
        i,j,k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
            c3 += 1
    info = analyze_omega2(A, n)
    dim_A2 = info['dim_A2']
    dim_O2 = info['dim_omega2']
    match = "YES" if dim_O2 == dim_A2 - c3 else "NO"
    if match == "NO":
        print(f"  T#{bits}: dim(A_2)={dim_A2}, c3={c3}, dim(Omega_2)={dim_O2}, A2-c3={dim_A2-c3} [{match}]")

# Check all
all_match = True
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    c3 = 0
    for triple in combinations(range(n), 3):
        i,j,k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
            c3 += 1
    info = analyze_omega2(A, n)
    if info['dim_omega2'] != info['dim_A2'] - c3:
        all_match = False
        break

print(f"  dim(Omega_2) = dim(A_2) - c3 for ALL n=5 tournaments? {all_match}")

# Also check n=4
n = 4
all_match_4 = True
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    c3 = 0
    for triple in combinations(range(n), 3):
        i,j,k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
            c3 += 1
    info = analyze_omega2(A, n)
    if info['dim_omega2'] != info['dim_A2'] - c3:
        all_match_4 = False
        print(f"  n=4 FAIL: T#{bits}: dim(A_2)={info['dim_A2']}, c3={c3}, dim(Omega_2)={info['dim_omega2']}")
        break
print(f"  dim(Omega_2) = dim(A_2) - c3 for ALL n=4 tournaments? {all_match_4}")


print("\n\nDone.")
