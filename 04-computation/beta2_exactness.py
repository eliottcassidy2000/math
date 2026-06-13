#!/usr/bin/env python3
"""
beta2_exactness.py - Exactness analysis of path chain complex for tournaments

beta2=0 means EXACTNESS at dimension 2: ker(d2) = im(d3) in Om_2.

Key question: What structural property of tournaments forces this exactness?

Approach: Compute Euler characteristic chi = sum(-1)^p dim(Om_p)
and check if it equals sum(-1)^p beta_p with beta_2=0.

Also: analyze the "exactness gap" = dim(Z2) - rk(d3) as a function
of tournament structure.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


# ============================================================
# Euler characteristic at n=5
# ============================================================
print("=" * 70)
print("PATH COMPLEX EULER CHARACTERISTIC")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

results = []
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

    betti = path_betti_numbers(A, n, max_dim=4)

    # Omega dimensions
    a_prev = [(i,) for i in range(n)]
    om_dims = [n]
    for p in range(1, 5):
        ap = enumerate_allowed_paths(A, n, p)
        om = compute_omega_basis(A, n, p, ap, a_prev)
        d = om.shape[1] if om.ndim == 2 else 0
        om_dims.append(d)
        a_prev = ap

    chi = sum((-1)**p * om_dims[p] for p in range(len(om_dims)))

    results.append({
        'bits': bits, 'scores': scores, 'c3': c3,
        'om_dims': om_dims, 'betti': betti, 'chi': chi,
    })

# Chi distribution
chi_vals = Counter(d['chi'] for d in results)
print(f"chi distribution at n=5: {dict(sorted(chi_vals.items()))}")

# Check: is chi always = 1?
if all(d['chi'] == 1 for d in results):
    print("chi = 1 for ALL tournaments at n=5!")
else:
    print("chi varies!")

# By score
by_score = defaultdict(list)
for d in results:
    by_score[d['scores']].append(d)

print(f"\n{'Scores':<20} {'c3':>3} {'chi':>4} {'Om0-4':<30} {'betti':<20}")
for sc in sorted(by_score.keys()):
    entries = by_score[sc]
    om_sets = set(tuple(d['om_dims']) for d in entries)
    for om_val in sorted(om_sets):
        cnt = sum(1 for d in entries if tuple(d['om_dims']) == om_val)
        d0 = [d for d in entries if tuple(d['om_dims']) == om_val][0]
        om_str = str(list(om_val))
        print(f"{str(sc):<20} {d0['c3']:>3} {d0['chi']:>4} {om_str:<30} {d0['betti']} ({cnt}x)")


# ============================================================
# KEY: Exactness at each dimension
# ============================================================
print(f"\n{'='*70}")
print("EXACTNESS AT EACH DIMENSION")
print("=" * 70)

# beta_p = 0 means exact at dimension p
# Check which dimensions are exact
for p in range(5):
    exact = sum(1 for d in results if d['betti'][p] == 0)
    print(f"  dim {p}: exact in {exact}/{len(results)} tournaments")

# The pattern should be:
# dim 0: beta0=1 always (never exact)
# dim 1: beta1 = 0 or 1 (sometimes exact)
# dim 2: beta2 = 0 always (ALWAYS exact) <-- this is our conjecture
# dim 3: beta3 = 0 mostly (n=5: always exact since beta3=0 at n=5)
# dim 4: beta4 = 0 always at n=5


# ============================================================
# CRITICAL: What's special about dim 2?
# For exactness at dim p, we need rk(d_{p+1}) = dim(Z_p)
# Equivalently: rk(d_{p+1}) + rk(d_p) = dim(Om_p)
# ============================================================
print(f"\n{'='*70}")
print("RANK SUM = dim(Om) AT EACH DIMENSION")
print("=" * 70)

for bits in [0, 341, 10]:
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))

    a_lists = {}
    for p in range(6):
        a_lists[p] = enumerate_allowed_paths(A, n, p) if p > 0 else [(i,) for i in range(n)]

    om_dict = {}
    om_dict[0] = np.eye(n)
    for p in range(1, 5):
        om = compute_omega_basis(A, n, p, a_lists[p], a_lists[p-1])
        om_dict[p] = om if om.ndim == 2 else np.zeros((len(a_lists[p]), 0))

    bd_dict = {}
    for p in range(1, 5):
        bd_dict[p] = build_full_boundary_matrix(a_lists[p], a_lists[p-1])

    # Compute ranks of d_p|_Om_p
    rk = {}
    for p in range(1, 5):
        dp = om_dict[p].shape[1]
        if dp > 0:
            S = np.linalg.svd(bd_dict[p] @ om_dict[p], compute_uv=False)
            rk[p] = sum(s > 1e-8 for s in S)
        else:
            rk[p] = 0

    print(f"\nbits={bits}, scores={scores}")
    print(f"  dim(Om_p): {[om_dict[p].shape[1] if om_dict[p].ndim==2 else 0 for p in range(5)]}")
    print(f"  rk(d_p):   {[0] + [rk[p] for p in range(1, 5)]}")

    # Check: rk(d_{p+1}) + rk(d_p) = dim(Om_p) ?
    for p in range(1, 4):
        dp = om_dict[p].shape[1] if om_dict[p].ndim == 2 else 0
        rk_in = rk.get(p+1, 0)  # rank of d_{p+1} restricted to Om_{p+1}
        rk_out = rk[p]  # rank of d_p restricted to Om_p
        gap = dp - rk_in - rk_out  # = beta_p
        print(f"    dim {p}: rk(d_{p+1})={rk_in} + rk(d_p)={rk_out} = {rk_in+rk_out}, dim(Om_p)={dp}, gap(=beta_p)={gap}")


# ============================================================
# LOOK FOR FORMULA: Can dim(Om_p) be expressed combinatorially?
# ============================================================
print(f"\n{'='*70}")
print("COMBINATORIAL FORMULAS FOR dim(Om_p)")
print("=" * 70)

n = 5
print("dim(Om_0) = n = 5 (always)")
print("dim(Om_1) = C(n,2) - n = C(n-1,2) = 6 for tournaments (always)")

# Check Om_2
om2_data = []
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0

    # Count junk pairs
    j2 = sum(1 for p in a2 if A[p[2]][p[0]])  # (a,c) backward

    om2_data.append({
        'bits': bits, 'scores': scores,
        'A2': len(a2), 'Om2': d_om2, 'J2': j2,
        'diff': len(a2) - d_om2,
    })

# Is Om2 = A2 - J2 + something?
print(f"\nOm2 analysis:")
by_score2 = defaultdict(list)
for d in om2_data:
    by_score2[d['scores']].append(d)

for sc in sorted(by_score2.keys()):
    entries = by_score2[sc]
    vals = set((d['A2'], d['Om2'], d['J2'], d['diff']) for d in entries)
    if len(vals) == 1:
        d = entries[0]
        print(f"  {sc}: |A2|={d['A2']}, Om2={d['Om2']}, J2={d['J2']}, |A2|-Om2={d['diff']}")
    else:
        for v in sorted(vals):
            cnt = sum(1 for d in entries if (d['A2'], d['Om2'], d['J2'], d['diff']) == v)
            print(f"  {sc}: |A2|={v[0]}, Om2={v[1]}, J2={v[2]}, |A2|-Om2={v[3]} ({cnt}x)")


# ============================================================
# KEY TEST: Is dim(Om2) = |A2| - J2 + c3?
# Or dim(Om2) = |A2| - J2 + (some function of junk structure)?
# ============================================================
print(f"\n{'='*70}")
print("Om2 FORMULA SEARCH")
print("=" * 70)

# Compute more invariants
for d in om2_data:
    A = build_adj(n, d['bits'])
    # Count unique backward pairs (a,c) with c->a that have at least one 2-path a->b->c
    backward_pairs_with_path = set()
    for p in enumerate_allowed_paths(A, n, 2):
        if A[p[2]][p[0]]:  # c->a
            backward_pairs_with_path.add((p[0], p[2]))
    d['n_bp'] = len(backward_pairs_with_path)

    # For each backward pair, count multiplicity
    bp_mult = Counter()
    for p in enumerate_allowed_paths(A, n, 2):
        if A[p[2]][p[0]]:
            bp_mult[(p[0], p[2])] += 1
    d['bp_mult'] = bp_mult
    d['sum_mult_choose_2'] = sum(m*(m-1)//2 for m in bp_mult.values())

# Check: Om2 = A2 - n_bp (number of junk constraint rows)?
# Or: Om2 = A2 - rk(constraint matrix)?
# The constraint matrix has one row per backward pair, columns for junk 2-paths
# If constraints are independent: Om2 = A2 - n_bp
print(f"Is Om2 = |A2| - n_bp (# backward pairs with paths)?")
check1 = sum(1 for d in om2_data if d['Om2'] == d['A2'] - d['n_bp'])
print(f"  Matches: {check1}/{len(om2_data)}")

# Actually: J2 = total number of junk 2-paths
# n_bp = number of backward pairs (groups)
# The constraint matrix has n_bp rows and J2 columns
# Each row sums the coefficients of all junk paths for that pair
# Rank <= min(n_bp, J2)
# If all constraints are independent: Om2 = A2 - n_bp

if check1 == len(om2_data):
    print("  CONFIRMED: Om2 = |A2| - (# backward pairs with intermediaries)")
    # What is n_bp in terms of tournament invariants?
    print(f"\n  n_bp values by score:")
    for sc in sorted(by_score2.keys()):
        entries = by_score2[sc]
        bp_vals = set(d['n_bp'] for d in entries)
        if len(bp_vals) == 1:
            print(f"    {sc}: n_bp={entries[0]['n_bp']}")
        else:
            for v in sorted(bp_vals):
                cnt = sum(1 for d in entries if d['n_bp'] == v)
                print(f"    {sc}: n_bp={v} ({cnt}x)")


print("\n\nDone.")
