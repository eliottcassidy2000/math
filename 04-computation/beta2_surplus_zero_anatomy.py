#!/usr/bin/env python3
"""
beta2_surplus_zero_anatomy.py — Deep analysis of surplus=0 tournaments

At n=5: 17/64 tilings have dim(Ω₃) = dim(Z₂), meaning EVERY Ω₃ element
is needed to fill Z₂ — there is NO room to spare.

In these cases, β₂=0 requires im(∂₃) = ker(∂₂) EXACTLY, meaning:
1. ∂₃ restricted to Ω₃ is injective (ker(∂₃|Ω₃) = 0, so β₃ = 0)
2. The image of ∂₃ exactly covers ker(∂₂)

Understanding these cases should reveal the algebraic mechanism for β₂=0.

At n=6: tight cases (surplus=1) correspond to β₃=1 tournaments (80/32768).
Surplus=0 should give β₃=0 tournaments where filling is maximally efficient.

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[j][i] and A[i][k] and A[k][j]: t3 += 1
    return t3

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def transitive_triples(A, n):
    """Count directed transitive triples (a→b→c with a→c)."""
    count = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    count += 1
    return count

# ======================================================================
print("="*70)
print("SURPLUS-ZERO TOURNAMENT ANATOMY")
print("="*70)

for n in [5, 6]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    surplus_groups = defaultdict(list)  # surplus -> list of tournament data
    total = 0

    for A in all_tournaments(n):
        total += 1
        if total % 5000 == 0:
            print(f"  ... {total}", flush=True)

        a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]
        a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
        a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        om3 = compute_omega_basis(A, n, 3, a3, a2)

        d_om2 = om2.shape[1] if om2.ndim == 2 else 0
        d_om3 = om3.shape[1] if om3.ndim == 2 else 0

        if d_om2 > 0:
            bd2 = build_full_boundary_matrix(a2, a1)
            bd2_om = bd2 @ om2
            S = np.linalg.svd(bd2_om, compute_uv=False)
            rk2 = sum(s > 1e-8 for s in S)
            z2 = d_om2 - rk2
        else:
            z2 = 0

        surplus = d_om3 - z2
        t3 = count_3cycles(A, n)
        sc = score_seq(A, n)
        tt = transitive_triples(A, n)

        if surplus <= 2:  # detailed analysis for tight cases
            # Compute ker(∂₃|Ω₃) = β₃ + im(∂₄ hitting ker ∂₃)
            if d_om3 > 0:
                bd3 = build_full_boundary_matrix(a3, a2)
                bd3_om = bd3 @ om3
                # Project onto Ω₂
                if om2.ndim == 2 and om2.shape[1] > 0:
                    coords, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
                    S3 = np.linalg.svd(coords, compute_uv=False)
                    rk3 = sum(s > 1e-8 for s in S3)
                else:
                    rk3 = 0
                ker3 = d_om3 - rk3
            else:
                ker3 = 0
                rk3 = 0

            # DT 4-paths
            dt = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])

            surplus_groups[surplus].append({
                't3': t3, 'sc': sc, 'tt': tt,
                'd_om2': d_om2, 'd_om3': d_om3,
                'z2': z2, 'rk3': rk3, 'ker3': ker3,
                'dt': dt, 'a2': len(a2), 'a3': len(a3),
            })

    print(f"\n  Total: {total} tournaments")

    # Surplus distribution
    print(f"\n  Surplus distribution:")
    all_surpluses = []
    for A in all_tournaments(n):
        a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]
        a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
        a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
        om2 = compute_omega_basis(A, n, 2, a2, a1)
        om3 = compute_omega_basis(A, n, 3, a3, a2)
        d2 = om2.shape[1] if om2.ndim == 2 else 0
        d3 = om3.shape[1] if om3.ndim == 2 else 0
        if d2 > 0:
            bd2 = build_full_boundary_matrix(a2, a1)
            S = np.linalg.svd(bd2 @ om2, compute_uv=False)
            z2 = d2 - sum(s > 1e-8 for s in S)
        else:
            z2 = 0
        all_surpluses.append(d3 - z2)

    for s, cnt in sorted(Counter(all_surpluses).items()):
        pct = 100 * cnt / total
        print(f"    surplus={s}: {cnt} ({pct:.1f}%)")

    # Analyze surplus=0 cases
    print(f"\n  === SURPLUS=0 CASES (dim(Ω₃) = dim(Z₂) exactly) ===")
    s0 = surplus_groups[0]
    if s0:
        print(f"  Count: {len(s0)}")
        print(f"\n  {'t3':>3} {'score':<18} {'TT':>4} {'|A₂|':>5} {'Ω₂':>4} {'|A₃|':>5} {'Ω₃':>4} {'Z₂':>4} {'DT':>4} {'ker∂₃':>6} {'rk∂₃':>5}")
        score_cnt = Counter()
        for d in s0:
            score_cnt[d['sc']] += 1
            if len([x for x in s0 if x['sc'] == d['sc']]) <= 5 or d == s0[0]:
                print(f"  {d['t3']:>3} {str(d['sc']):<18} {d['tt']:>4} {d['a2']:>5} {d['d_om2']:>4} {d['a3']:>5} {d['d_om3']:>4} {d['z2']:>4} {d['dt']:>4} {d['ker3']:>6} {d['rk3']:>5}")

        print(f"\n  Score sequence distribution in surplus=0:")
        for sc, cnt in sorted(score_cnt.items()):
            print(f"    {sc}: {cnt}")

        # KEY: ker(∂₃) = 0 for surplus=0?
        # If surplus=0 and β₂=0, then rk(∂₃) = z2 = d_om3, so ker(∂₃) = 0
        all_ker0 = all(d['ker3'] == 0 for d in s0)
        print(f"\n  All surplus=0 have ker(∂₃)=0: {all_ker0}")
        if all_ker0:
            print(f"    This means: ∂₃ is INJECTIVE on Ω₃ for surplus=0 tournaments")
            print(f"    Equivalently: β₃=0 for all surplus=0 tournaments")
            print(f"    And: dim(im ∂₃) = dim(Ω₃) = dim(ker ∂₂) = dim(Z₂)")

        # DT vs Ω₃: are all Ω₃ elements DT paths?
        dt_equals_omega3 = all(d['dt'] == d['d_om3'] for d in s0)
        print(f"  DT = Ω₃ for surplus=0: {dt_equals_omega3}")
        if not dt_equals_omega3:
            for d in s0:
                if d['dt'] != d['d_om3']:
                    print(f"    Counterexample: DT={d['dt']}, Ω₃={d['d_om3']}")
                    break

    # Surplus=1 cases
    print(f"\n  === SURPLUS=1 CASES ===")
    s1 = surplus_groups[1]
    if s1:
        print(f"  Count: {len(s1)}")
        ker3_dist = Counter(d['ker3'] for d in s1)
        print(f"  ker(∂₃) distribution: {dict(ker3_dist)}")
        if 1 in ker3_dist:
            print(f"    ker(∂₃)=1 corresponds to β₃=1 (one 3-cycle)")

    # Surplus=2 cases
    print(f"\n  === SURPLUS=2 CASES ===")
    s2 = surplus_groups[2]
    if s2:
        print(f"  Count: {len(s2)}")
        ker3_dist = Counter(d['ker3'] for d in s2)
        print(f"  ker(∂₃) distribution: {dict(ker3_dist)}")

    if n >= 6:
        break  # n=6 computation already done above

# ======================================================================
# DEEP DIVE: n=5 surplus=0 — what makes ∂₃ injective?
# ======================================================================
print(f"\n\n{'='*70}")
print("DEEP DIVE: WHY IS ∂₃ INJECTIVE FOR SURPLUS=0?")
print("="*70)

n = 5
count = 0
for A in all_tournaments(n):
    count += 1
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)

    d2 = om2.shape[1] if om2.ndim == 2 else 0
    d3 = om3.shape[1] if om3.ndim == 2 else 0

    if d2 > 0:
        bd2 = build_full_boundary_matrix(a2, a1)
        S = np.linalg.svd(bd2 @ om2, compute_uv=False)
        z2 = d2 - sum(s > 1e-8 for s in S)
    else:
        z2 = 0

    surplus = d3 - z2

    if surplus == 0 and count <= 1024:
        # Show the explicit Ω₃ basis
        t3 = count_3cycles(A, n)
        scores = score_seq(A, n)

        if count <= 5 or scores == (0,1,2,3,4):  # first few + transitive
            print(f"\n--- Tournament #{count}, t3={t3}, scores={scores} ---")
            print(f"    Ω₂={d2}, Ω₃={d3}, Z₂={z2}, surplus=0")

            # Show Ω₃ basis
            if d3 > 0:
                for j in range(d3):
                    col = om3[:, j]
                    nz = [(a3[i], col[i]) for i in range(len(col)) if abs(col[i]) > 1e-10]
                    if len(nz) == 1:
                        path, coeff = nz[0]
                        is_dt = A[path[0]][path[2]] and A[path[1]][path[3]]
                        print(f"    Ω₃ basis {j}: {path} (DT={is_dt})")
                    else:
                        print(f"    Ω₃ basis {j}: {len(nz)} terms: {nz[:5]}")

            # Show Z₂ basis (2-cycles)
            if z2 > 0:
                bd2_full = build_full_boundary_matrix(a2, a1)
                bd2_om = bd2_full @ om2
                U, S2, Vt = np.linalg.svd(bd2_om)
                rk = sum(s > 1e-8 for s in S2)
                # ker(∂₂) = last d2-rk rows of Vt
                ker_basis = Vt[rk:]  # rows are in Ω₂ coords
                for j in range(z2):
                    v = ker_basis[j]
                    # Convert to A₂ coords
                    z_vec = om2 @ v
                    nz = [(a2[i], z_vec[i]) for i in range(len(z_vec)) if abs(z_vec[i]) > 1e-10]
                    print(f"    Z₂ basis {j}: {nz[:6]}{'...' if len(nz)>6 else ''}")

            # Show ∂₃ matrix (should be injective)
            if d3 > 0 and d2 > 0:
                bd3 = build_full_boundary_matrix(a3, a2)
                bd3_om = bd3 @ om3
                coords, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
                print(f"    ∂₃ matrix ({z2}×{d3}):")
                # Project to Z₂ coords
                # coords is d_om2 × d_om3, we want to restrict to Z₂
                ker_proj = Vt[rk:] @ coords  # z2 × d_om3
                for i in range(min(z2, 5)):
                    row = [f"{ker_proj[i,j]:.2f}" for j in range(d3)]
                    print(f"      [{', '.join(row)}]")

print(f"\n{'='*70}")
print("DONE")
print("="*70)
