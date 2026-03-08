#!/usr/bin/env python3
"""
beta2_surplus_algebraic.py — Algebraic analysis of surplus under arc flips

GOAL: Understand why surplus = dim(Ω₃) - dim(Z₂) ≥ 0 always.

The key identity (rank-nullity):
  surplus = dim(Ω₃) - dim(Z₂)
          = dim(Ω₃) - dim(ker ∂₂|Ω₂)
          = dim(Ω₃) - (dim(Ω₂) - rk(∂₂|Ω₂))
          = dim(Ω₃) - dim(Ω₂) + rk(∂₂|Ω₂)

And β₂ = 0 iff surplus = rk(∂₃|Ω₃), i.e., rk(∂₃|Ω₃) = dim(Ω₃) - dim(Z₂).

APPROACH 1: Exact formulas for dim(Ω_p) in terms of tournament structure.

For Ω₂: A 2-path (a,b,c) ∈ Ω₂ iff ∂₂(a,b,c) ∈ Ω₁.
  ∂₂(a,b,c) = (b,c) - (a,c) + (a,b).
  For Ω₁ = A₁ (all edges): need (a,c) ∈ A₁, i.e., a→c.
  So INDIVIDUAL 2-paths in Ω₂ iff (a,b,c) is TT.

  But LINEAR COMBINATIONS can be in Ω₂ with NT components, if the
  "bad face" terms cancel. Specifically:
  dim(Ω₂) = |TT| + dim(NT_cancel)
  where NT_cancel = {combinations of NT 2-paths whose bad faces cancel}

  The NT cancellation space has dimension = rank of the "bad face incidence matrix":
  For each NT triple (a,b,c) with c→a, the bad face is (a,c).
  Since (a,c) is NOT an edge (c→a instead), (a,c) lives in A₁^c.
  Multiple NT triples can share the same bad face (a,c).
  If k NT triples share bad face (a,c), they contribute (k-1) to NT_cancel dim.

  So: dim(NT_cancel) = |NT| - #(distinct bad faces among NT triples)
                     = |NT| - #(non-edges used as middle faces)

For Ω₃: A 3-path (a,b,c,d) ∈ Ω₃ iff ∂₃(a,b,c,d) ∈ Ω₂.
  ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c).
  Need each face IN Ω₂ or their bad-face terms to cancel.
  Face (a,b,c): TT iff a→c (true if DT condition a→c)
  Face (a,b,d): TT iff a→d
  Face (a,c,d): TT iff a→d
  Face (b,c,d): TT iff b→d (true if DT condition b→d)

  Individual 3-paths in Ω₃: need all 4 faces in Ω₂.
  This means all 4 faces are TT, which requires:
  a→c (DT), b→d (DT), AND a→d.
  So individual Ω₃ elements are DT paths WITH a→d (what I'll call "super-DT").

  Wait — is this right? Actually Ω₃ elements need ∂₃(element) ∈ Ω₂,
  not each face individually in Ω₂. A linear combination could have
  non-Ω₂ face terms cancel.

  Let me compute this more carefully...

Author: opus-2026-03-08-S44
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

print("=" * 70)
print("ALGEBRAIC SURPLUS ANALYSIS")
print("=" * 70)

# ===== DIMENSION FORMULAS =====
print("\n--- Testing dim(Ω₂) formula: |TT| + NT_cancel ---")

for n in [4, 5]:
    print(f"\n  n = {n}:")
    nt_cancel_formula_works = 0
    total = 0

    for A in all_tournaments(n):
        total += 1
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        if not a2:
            continue

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

        tt_count = sum(1 for p in a2 if A[p[0]][p[2]] == 1)

        # NT triples and their bad faces
        nt_triples = [p for p in a2 if A[p[0]][p[2]] == 0]  # c→a, so (a,c) is bad
        bad_faces = set()
        for p in nt_triples:
            a, b, c = p
            bad_faces.add((a, c))  # or should it be (c, a)? The non-edge pair

        nt_cancel_dim = len(nt_triples) - len(bad_faces)
        predicted = tt_count + nt_cancel_dim

        if predicted == dim_om2:
            nt_cancel_formula_works += 1

    print(f"    Formula dim(Ω₂) = |TT| + |NT| - |bad_faces|: "
          f"{nt_cancel_formula_works}/{total}")

# ===== Ω₃ STRUCTURE: INDIVIDUAL vs CANCELLATION =====
print(f"\n--- Testing Ω₃ individual elements ---")

for n in [4, 5]:
    print(f"\n  n = {n}:")
    individual_in_om3 = Counter()  # (is_DT, a_to_d) → count
    om3_dim_formula = 0
    total_checked = 0

    for A in all_tournaments(n):
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        if not a3:
            continue

        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        # Check each individual 3-path
        bd3 = build_full_boundary_matrix(a3, a2)

        # For each 3-path, check if it's individually in Ω₃
        # by checking if ∂₃(path) is in Ω₂
        om2 = compute_omega_basis(A, n, 2, a2,
                                   enumerate_allowed_paths(A, n, 1))

        for j, p in enumerate(a3):
            a, b, c, d = p
            is_dt = A[a][c] == 1 and A[b][d] == 1
            a_to_d = A[a][d] == 1

            # Check if boundary is in Ω₂
            boundary = bd3[:, j]  # in A₂ coordinates
            # Project onto Ω₂
            if om2.shape[1] > 0:
                coords, res, _, _ = np.linalg.lstsq(om2, boundary, rcond=None)
                recon = om2 @ coords
                error = np.linalg.norm(boundary - recon)
                in_om2 = error < 1e-8
            else:
                in_om2 = np.linalg.norm(boundary) < 1e-8

            individual_in_om3[(is_dt, a_to_d)] += (1 if in_om2 else 0)

        total_checked += len(a3)

        # Count super-DT (DT + a→d)
        super_dt = sum(1 for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1 and A[p[0]][p[3]] == 1)

        # Count DT with d→a
        dt_backward = sum(1 for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1 and A[p[3]][p[0]] == 1)

        # Non-DT
        non_dt = len(a3) - super_dt - dt_backward

    print(f"    Individual 3-paths in Ω₃ by type (is_DT, a→d):")
    for key in sorted(individual_in_om3):
        print(f"      {key}: {individual_in_om3[key]}")

# ===== EXACT SURPLUS FORMULA SEARCH =====
print(f"\n{'='*70}")
print("SEARCHING FOR EXACT SURPLUS FORMULA")
print("="*70)

n = 5
print(f"\n  n = {n}: exhaustive")

data_points = []
for A in all_tournaments(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

    om1 = compute_omega_basis(A, n, 1, a1, enumerate_allowed_paths(A, n, 0))
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    S = np.linalg.svd(coords, compute_uv=False)
    rk2 = int(sum(s > 1e-8 for s in S))
    z2 = dim_om2 - rk2

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    surplus = dim_om3 - z2

    tt = sum(1 for p in a2 if A[p[0]][p[2]] == 1)
    nt = len(a2) - tt
    dt = sum(1 for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1)
    super_dt = sum(1 for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1 and A[p[0]][p[3]] == 1)
    dt_back = dt - super_dt

    # Bad faces
    nt_triples = [p for p in a2 if A[p[0]][p[2]] == 0]
    bad_faces = set()
    for p in nt_triples:
        bad_faces.add((p[0], p[2]))
    nt_cancel = len(nt_triples) - len(bad_faces)

    # 3-path bad face analysis
    one_bad_3 = 0
    two_bad_3 = 0
    for p in a3:
        a, b, c, d = p
        # Faces of 3-path:
        # (b,c,d) bad iff b not→d
        # (a,c,d) bad iff a not→d
        # (a,b,d) bad iff a not→d
        # (a,b,c) bad iff a not→c
        bad_count = 0
        if A[b][d] == 0: bad_count += 1  # (b,c,d) bad
        if A[a][d] == 0: bad_count += 2  # (a,c,d) and (a,b,d) both bad
        if A[a][c] == 0: bad_count += 1  # (a,b,c) bad
        if bad_count == 1:
            one_bad_3 += 1
        elif bad_count >= 2:
            two_bad_3 += 1

    scores = tuple(sorted(sum(A[i]) for i in range(n)))

    data_points.append({
        'A2': len(a2), 'A3': len(a3),
        'Om2': dim_om2, 'Om3': dim_om3,
        'Z2': z2, 'rk2': rk2, 'surplus': surplus,
        'TT': tt, 'NT': nt, 'DT': dt,
        'super_DT': super_dt, 'DT_back': dt_back,
        'NT_cancel': nt_cancel, 'bad_faces_2': len(bad_faces),
        'one_bad_3': one_bad_3, 'two_bad_3': two_bad_3,
        'scores': scores
    })

# Check: is surplus = super_DT - rk2?
# Or: surplus = super_DT + something - Z₂?
print(f"\n  Testing formulas for surplus:")

# Formula 1: surplus = dim(Ω₃) - Z₂
# We need to understand dim(Ω₃) and Z₂ separately.

# Check: dim(Ω₂) = TT + NT_cancel
om2_formula = all(d['Om2'] == d['TT'] + d['NT_cancel'] for d in data_points)
print(f"    dim(Ω₂) = TT + NT_cancel: {om2_formula}")

# Check: Z₂ = dim(Ω₂) - rk(∂₂)
# Already by definition, this is always true.

# Check: is rk(∂₂) a simple function of the tournament?
rk2_values = defaultdict(set)
for d in data_points:
    rk2_values[d['scores']].add(d['rk2'])
print(f"\n    rk(∂₂) variation by score sequence:")
for scores, vals in sorted(rk2_values.items()):
    if len(vals) > 1:
        print(f"      scores={scores}: rk2 ∈ {sorted(vals)}")

# Check: surplus vs (DT - rk2) or similar
surplus_formulas = defaultdict(int)
for d in data_points:
    # Try: surplus = DT - rk2
    if d['surplus'] == d['DT'] - d['rk2']:
        surplus_formulas['DT - rk2'] += 1
    # Try: surplus = super_DT + DT_back - rk2
    if d['surplus'] == d['super_DT'] + d['DT_back'] - d['rk2']:
        surplus_formulas['super_DT + DT_back - rk2'] += 1
    # Try: surplus = super_DT - Z2 + DT_back
    if d['surplus'] == d['super_DT'] - d['Z2'] + d['DT_back']:
        surplus_formulas['super_DT - Z2 + DT_back'] += 1
    # Since dim(Ω₃) is complex, try just correlations
    surplus_formulas['any'] += 1

print(f"\n    Formula matches (out of {len(data_points)}):")
for formula, count in sorted(surplus_formulas.items()):
    if count > 0:
        print(f"      {formula}: {count}")

# Direct check: what IS dim(Ω₃)?
print(f"\n  dim(Ω₃) breakdown:")
om3_by_dt_parts = defaultdict(list)
for d in data_points:
    om3_by_dt_parts[(d['super_DT'], d['DT_back'])].append(d['Om3'])

for key in sorted(om3_by_dt_parts):
    vals = om3_by_dt_parts[key]
    uniqs = sorted(set(vals))
    print(f"    (super_DT={key[0]}, DT_back={key[1]}): Ω₃ ∈ {uniqs}")

# ===== CRUCIAL INSIGHT SEARCH: δ(surplus) UNDER ARC FLIP =====
print(f"\n{'='*70}")
print("δ(SURPLUS) DECOMPOSITION UNDER ARC FLIP")
print("="*70)

n = 5
# For each arc flip u→v: decompose δ(surplus) = δ(Ω₃) - δ(Z₂)
# And δ(Z₂) = δ(Ω₂) - δ(rk₂)
# So δ(surplus) = δ(Ω₃) - δ(Ω₂) + δ(rk₂)

delta_decomp = defaultdict(list)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    om1 = compute_omega_basis(A, n, 1, a1, enumerate_allowed_paths(A, n, 0))
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    S = np.linalg.svd(coords, compute_uv=False)
    rk2 = int(sum(s > 1e-8 for s in S))
    z2 = dim_om2 - rk2

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    surplus = dim_om3 - z2

    # For each arc flip
    for idx2, (i,j) in enumerate(pairs):
        bits2 = bits ^ (1 << idx2)
        A2 = [[0]*n for _ in range(n)]
        for idx3, (ii,jj) in enumerate(pairs):
            if (bits2 >> idx3) & 1: A2[ii][jj] = 1
            else: A2[jj][ii] = 1

        a1_2 = enumerate_allowed_paths(A2, n, 1)
        a2_2 = enumerate_allowed_paths(A2, n, 2)
        a3_2 = enumerate_allowed_paths(A2, n, 3)

        om2_2 = compute_omega_basis(A2, n, 2, a2_2, a1_2)
        dim_om2_2 = om2_2.shape[1] if om2_2.ndim == 2 else 0
        om1_2 = compute_omega_basis(A2, n, 1, a1_2, enumerate_allowed_paths(A2, n, 0))
        bd2_2 = build_full_boundary_matrix(a2_2, a1_2)
        bd2_om_2 = bd2_2 @ om2_2
        coords2, _, _, _ = np.linalg.lstsq(om1_2, bd2_om_2, rcond=None)
        S2 = np.linalg.svd(coords2, compute_uv=False)
        rk2_2 = int(sum(s > 1e-8 for s in S2))
        z2_2 = dim_om2_2 - rk2_2

        om3_2 = compute_omega_basis(A2, n, 3, a3_2, a2_2)
        dim_om3_2 = om3_2.shape[1] if om3_2.ndim == 2 else 0
        surplus_2 = dim_om3_2 - z2_2

        # Determine flip direction
        if A[i][j] == 1:
            u, v = i, j
        else:
            u, v = j, i
        du = sum(A[u])
        dv = sum(A[v])

        delta_decomp[(du - dv)].append({
            'dOm2': dim_om2_2 - dim_om2,
            'dOm3': dim_om3_2 - dim_om3,
            'dZ2': z2_2 - z2,
            'drk2': rk2_2 - rk2,
            'dsurplus': surplus_2 - surplus,
            'du': du, 'dv': dv,
        })

print(f"\n  δ(surplus) = δ(Ω₃) - δ(Ω₂) + δ(rk₂) decomposition:")
print(f"  {'d_u-d_v':>7}  {'mean_dOm3':>9}  {'mean_dOm2':>9}  {'mean_drk2':>9}  {'mean_dsurp':>10}  {'min_dsurp':>9}  {'N':>5}")

for diff in sorted(delta_decomp):
    entries = delta_decomp[diff]
    dOm3 = np.mean([e['dOm3'] for e in entries])
    dOm2 = np.mean([e['dOm2'] for e in entries])
    drk2 = np.mean([e['drk2'] for e in entries])
    dsurp = np.mean([e['dsurplus'] for e in entries])
    min_dsurp = min(e['dsurplus'] for e in entries)
    print(f"  {diff:>7}  {dOm3:>9.2f}  {dOm2:>9.2f}  {drk2:>9.2f}  {dsurp:>10.2f}  {min_dsurp:>9}  {len(entries):>5}")

# ===== KEY CHECK: When surplus drops, is the drop bounded? =====
print(f"\n  When surplus drops maximally, what's the tournament structure?")

max_drop = 0
for diff in delta_decomp:
    for e in delta_decomp[diff]:
        if e['dsurplus'] < max_drop:
            max_drop = e['dsurplus']
print(f"  Maximum surplus drop in single flip: {max_drop}")

# Find the actual tournaments where max drop occurs
print(f"\n  Cases with largest surplus drops:")
drop_cases = []
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    om2 = compute_omega_basis(A, n, 2, a2, enumerate_allowed_paths(A, n, 1))
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    om1 = compute_omega_basis(A, n, 1, enumerate_allowed_paths(A, n, 1),
                               enumerate_allowed_paths(A, n, 0))
    bd2 = build_full_boundary_matrix(a2, enumerate_allowed_paths(A, n, 1))
    bd2_om = bd2 @ om2
    coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    S = np.linalg.svd(coords, compute_uv=False)
    rk2 = int(sum(s > 1e-8 for s in S))
    z2 = dim_om2 - rk2
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    surplus = dim_om3 - z2

    if surplus >= 4:  # Only care about cases with enough surplus to drop from
        for idx2, (i,j) in enumerate(pairs):
            bits2 = bits ^ (1 << idx2)
            # Quick surplus computation
            A2 = [[0]*n for _ in range(n)]
            for idx3, (ii,jj) in enumerate(pairs):
                if (bits2 >> idx3) & 1: A2[ii][jj] = 1
                else: A2[jj][ii] = 1

            a2_2 = enumerate_allowed_paths(A2, n, 2)
            a3_2 = enumerate_allowed_paths(A2, n, 3)
            om2_2 = compute_omega_basis(A2, n, 2, a2_2, enumerate_allowed_paths(A2, n, 1))
            dim_om2_2 = om2_2.shape[1] if om2_2.ndim == 2 else 0
            om1_2 = compute_omega_basis(A2, n, 1, enumerate_allowed_paths(A2, n, 1),
                                         enumerate_allowed_paths(A2, n, 0))
            bd2_2 = build_full_boundary_matrix(a2_2, enumerate_allowed_paths(A2, n, 1))
            bd2_om_2 = bd2_2 @ om2_2
            coords2, _, _, _ = np.linalg.lstsq(om1_2, bd2_om_2, rcond=None)
            S2 = np.linalg.svd(coords2, compute_uv=False)
            rk2_2 = int(sum(s > 1e-8 for s in S2))
            z2_2 = dim_om2_2 - rk2_2
            om3_2 = compute_omega_basis(A2, n, 3, a3_2, a2_2)
            dim_om3_2 = om3_2.shape[1] if om3_2.ndim == 2 else 0
            surplus2 = dim_om3_2 - z2_2

            drop = surplus2 - surplus
            if drop <= -4:
                scores = tuple(sorted(sum(A[i2]) for i2 in range(n)))
                if A[i][j] == 1:
                    u, v = i, j
                else:
                    u, v = j, i
                drop_cases.append({
                    'surplus_before': surplus,
                    'surplus_after': surplus2,
                    'drop': drop,
                    'scores': scores,
                    'du': sum(A[u]), 'dv': sum(A[v]),
                    'u': u, 'v': v
                })

drop_cases.sort(key=lambda x: x['drop'])
for case in drop_cases[:10]:
    print(f"    surplus {case['surplus_before']}→{case['surplus_after']} "
          f"(drop={case['drop']}): scores={case['scores']}, "
          f"flip {case['u']}→{case['v']} (d_u={case['du']}, d_v={case['dv']})")

# ===== SURPLUS LOWER BOUND AS FUNCTION OF d_u - d_v =====
print(f"\n  Surplus drop bounded by (surplus_before, d_u - d_v)?")
drop_by_surplus_and_dd = defaultdict(list)
for diff in delta_decomp:
    for e in delta_decomp[diff]:
        # We need the surplus_before too
        pass

# Better: directly compute
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    om2 = compute_omega_basis(A, n, 2, a2, enumerate_allowed_paths(A, n, 1))
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    om1 = compute_omega_basis(A, n, 1, enumerate_allowed_paths(A, n, 1),
                               enumerate_allowed_paths(A, n, 0))
    bd2 = build_full_boundary_matrix(a2, enumerate_allowed_paths(A, n, 1))
    bd2_om = bd2 @ om2
    coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    S = np.linalg.svd(coords, compute_uv=False)
    rk2 = int(sum(s > 1e-8 for s in S))
    z2 = dim_om2 - rk2
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    surplus = dim_om3 - z2

    for idx2, (i,j) in enumerate(pairs):
        bits2 = bits ^ (1 << idx2)
        A2 = [[0]*n for _ in range(n)]
        for idx3, (ii,jj) in enumerate(pairs):
            if (bits2 >> idx3) & 1: A2[ii][jj] = 1
            else: A2[jj][ii] = 1

        a2_2 = enumerate_allowed_paths(A2, n, 2)
        a3_2 = enumerate_allowed_paths(A2, n, 3)
        om2_2 = compute_omega_basis(A2, n, 2, a2_2, enumerate_allowed_paths(A2, n, 1))
        dim_om2_2 = om2_2.shape[1] if om2_2.ndim == 2 else 0
        om1_2 = compute_omega_basis(A2, n, 1, enumerate_allowed_paths(A2, n, 1),
                                     enumerate_allowed_paths(A2, n, 0))
        bd2_2 = build_full_boundary_matrix(a2_2, enumerate_allowed_paths(A2, n, 1))
        bd2_om_2 = bd2_2 @ om2_2
        coords2, _, _, _ = np.linalg.lstsq(om1_2, bd2_om_2, rcond=None)
        S2 = np.linalg.svd(coords2, compute_uv=False)
        rk2_2 = int(sum(s > 1e-8 for s in S2))
        z2_2 = dim_om2_2 - rk2_2
        om3_2 = compute_omega_basis(A2, n, 3, a3_2, a2_2)
        dim_om3_2 = om3_2.shape[1] if om3_2.ndim == 2 else 0
        surplus2 = dim_om3_2 - z2_2

        if A[i][j] == 1:
            u, v = i, j
        else:
            u, v = j, i
        dd = sum(A[u]) - sum(A[v])

        drop_by_surplus_and_dd[(surplus, dd)].append(surplus2 - surplus)

print(f"\n  Max drop by (surplus_before, d_u-d_v):")
print(f"  {'(surp,dd)':>12}  {'max_drop':>8}  {'min_after':>9}  {'count':>5}")
for key in sorted(drop_by_surplus_and_dd):
    vals = drop_by_surplus_and_dd[key]
    max_drop_val = min(vals)
    min_after = key[0] + max_drop_val
    print(f"  {str(key):>12}  {max_drop_val:>8}  {min_after:>9}  {len(vals):>5}")

print("\nDone.")
