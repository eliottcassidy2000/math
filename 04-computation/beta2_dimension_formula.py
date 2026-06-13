#!/usr/bin/env python3
"""
beta2_dimension_formula.py — Search for formulas that make β₂=0 obvious

If β₂=0, then rk(∂₃|_Ω₃) = dim(Z₂) = dim(Ω₂) - rk(∂₂|_Ω₂).

Known formulas:
- dim(Ω₁) = C(n,2)
- |A₂| = (n-1)*C(n,2) - Σs_b²  [s_b = out-degree of b]
- dim(Ω₂) = |A₂| - J₂
- z1 = C(n-1,2) always (tournament strongly connected => rk(∂₁) = n-1)

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def c3_count(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                s = A[i][j] + A[j][k] + A[k][i]
                if s == 0 or s == 3:
                    c += 1
    return c


print("=" * 70)
print("DIMENSION FORMULA SEARCH FOR β₂=0")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    cn2 = n*(n-1)//2
    cn12 = (n-1)*(n-2)//2

    data = []
    t0 = time.time()

    for bits in range(total):
        A = build_adj(n, bits)

        scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]
        s2 = sum(s*s for s in scores)
        c3 = c3_count(A, n)

        # |A₂| = Σ_b s_b*(n-1-s_b) = (n-1)*C(n,2) - s2
        a2_formula = (n-1)*cn2 - s2

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)

        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
        om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

        d2 = dim_om(om2)
        d3 = dim_om(om3)
        a2 = len(ap2) if ap2 else 0
        a3 = len(ap3) if ap3 else 0
        j2 = a2 - d2
        j3 = a3 - d3

        # |A₂| check
        assert a2 == a2_formula, f"A₂ formula wrong: {a2} vs {a2_formula}"

        # Ranks
        bd1 = build_full_boundary_matrix(ap1, ap0)
        rk_d1 = np.linalg.matrix_rank(bd1, tol=1e-8)

        if d2 > 0:
            bd2 = build_full_boundary_matrix(ap2, ap1)
            bd2_om = bd2 @ om2
            coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
            rk_d2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        else:
            rk_d2 = 0

        if d3 > 0 and d2 > 0:
            bd3 = build_full_boundary_matrix(ap3, ap2)
            bd3_om = bd3 @ om3
            coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
            rk_d3 = np.linalg.matrix_rank(coords3, tol=1e-8)
        else:
            rk_d3 = 0

        z1 = cn2 - rk_d1
        z2 = d2 - rk_d2
        beta1 = z1 - rk_d2

        # Count transitive 4-subsets
        t4 = 0
        for quad in combinations(range(n), 4):
            is_trans = True
            for trip in combinations(quad, 3):
                s = A[trip[0]][trip[1]] + A[trip[1]][trip[2]] + A[trip[2]][trip[0]]
                if s == 0 or s == 3:
                    is_trans = False
                    break
            if is_trans:
                t4 += 1

        data.append({
            'bits': bits, 'c3': c3, 's2': s2,
            'a2': a2, 'a3': a3, 'j2': j2, 'j3': j3,
            'd2': d2, 'd3': d3,
            'rk_d1': rk_d1, 'rk_d2': rk_d2, 'rk_d3': rk_d3,
            'z1': z1, 'z2': z2, 'beta1': beta1,
            't4': t4
        })

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")

    # ===== BASIC FACTS =====
    print(f"\n  rk(∂₁) = {n-1} always: {all(d['rk_d1'] == n-1 for d in data)}")
    print(f"  z1 = C(n-1,2) = {cn12}: {all(d['z1'] == cn12 for d in data)}")
    print(f"  β₂ = 0 always: {all(d['z2'] == d['rk_d3'] for d in data)}")

    # ===== J₂ FORMULA =====
    print(f"\n  J₂ analysis:")
    j2_by_c3 = defaultdict(set)
    for d in data:
        j2_by_c3[d['c3']].add(d['j2'])

    for key in sorted(j2_by_c3.keys()):
        vals = sorted(j2_by_c3[key])
        if len(vals) == 1:
            print(f"    c3={key} → J₂={vals[0]}")
        else:
            print(f"    c3={key} → J₂ ∈ {vals}")

    # Test J₂ = 3*c3
    test_j2_3c3 = all(d['j2'] == 3*d['c3'] for d in data)
    print(f"  J₂ = 3*c₃: {test_j2_3c3}")

    if not test_j2_3c3:
        # More careful: J₂ formula involving c3 and something else
        # J₂ = #{(a,c): c→a and ∃b with a→b→c}
        # = #{(a,c): c→a and A²[a,c] > 0}
        # For each such (a,c), how many b? A²[a,c] = #{b: a→b, b→c}
        # So J₂ = Σ_{(a,c): c→a, A²[a,c]>0} 1
        # Actually J₂ = Σ_{(a,c): c→a} [A²[a,c] > 0]
        # But also J₂ = Σ_{(a,c): c→a} A²[a,c] - (#{(a,c): c→a, A²[a,c]>0} - #{...})
        # Hmm, it's the number of JUNK CONSTRAINTS, not junk paths.
        # Actually dim(Ω₂) = |A₂| - J₂ where J₂ = #{junk pairs}
        # where a junk pair (a,c) with c→a means all paths (a,*,c) are constrained to sum to 0

        # Test: J₂ = #{NT pairs (a,c) with c→a such that A²[a,c]>0}
        all_match = True
        for d in data:
            A = build_adj(n, d['bits'])
            j2_test = 0
            for a in range(n):
                for c in range(n):
                    if a == c or A[a][c]:  # skip if a→c (TT direction)
                        continue
                    # c→a, count paths a→b→c
                    count = sum(1 for b in range(n) if b != a and b != c and A[a][b] and A[b][c])
                    if count > 0:
                        j2_test += 1
            if j2_test != d['j2']:
                all_match = False
                break
        print(f"  J₂ = #{'{NT pairs (a,c) with common out-nbr}'}: {all_match}")

    # ===== J₃ FORMULA =====
    print(f"\n  J₃ analysis:")
    j3_by_c3 = defaultdict(set)
    for d in data:
        j3_by_c3[d['c3']].add(d['j3'])

    for key in sorted(j3_by_c3.keys()):
        vals = sorted(j3_by_c3[key])
        if len(vals) == 1:
            print(f"    c3={key} → J₃={vals[0]}")
        else:
            print(f"    c3={key} → J₃ ∈ {vals}")

    # ===== z2 FORMULA =====
    print(f"\n  z2 = dim(Ω₂) - rk(∂₂|_Ω₂)")
    print(f"     = (|A₂| - J₂) - (z1 - β₁)")
    print(f"     = |A₂| - J₂ - C(n-1,2) + β₁")

    # Group by (c3, s2, beta1) -> z2
    z2_by = defaultdict(set)
    for d in data:
        z2_by[(d['c3'], d['s2'], d['beta1'])].add(d['z2'])

    print(f"\n  z2 determined by (c3, s2, β₁): {all(len(v)==1 for v in z2_by.values())}")
    for key in sorted(z2_by.keys()):
        vals = sorted(z2_by[key])
        c3_v, s2_v, b1 = key
        a2_v = (n-1)*cn2 - s2_v
        print(f"    c3={c3_v}, s2={s2_v}, β₁={b1} → z2={vals}, |A₂|={a2_v}")

    # ===== rk_d3 FORMULA =====
    # rk_d3 = d3 - dim(Z₃)
    # dim(Z₃) = dim(ker ∂₃|_Ω₃)
    # For β₂=0: rk_d3 = z2
    print(f"\n  rk_d3 analysis:")
    rk_d3_by = defaultdict(set)
    for d in data:
        rk_d3_by[(d['c3'], d['s2'], d['beta1'])].add(d['rk_d3'])

    for key in sorted(rk_d3_by.keys()):
        vals = sorted(rk_d3_by[key])
        print(f"    c3={key[0]}, s2={key[1]}, β₁={key[2]} → rk_d3={vals}")


# =============================================================
# PART 2: The A₃ and J₃ structure
# =============================================================
print(f"\n{'='*70}")
print("PART 2: Detailed J₂, J₃ formulas")
print("-" * 50)

for n in [5, 6]:
    m = n*(n-1)//2
    cn2 = n*(n-1)//2
    cn12 = (n-1)*(n-2)//2

    j2_formulas = []
    for bits in range(1 << m):
        A = build_adj(n, bits)
        scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]
        c3 = c3_count(A, n)
        s2 = sum(s*s for s in scores)

        # Count junk pairs (a,c) with c→a and A²[a,c]>0
        j2_count = 0
        a2_sum = 0  # = sum of A²[a,c] over (a,c) with c→a
        for a in range(n):
            for c in range(n):
                if a == c or A[a][c]:
                    continue
                cnt = sum(1 for b in range(n) if b not in (a,c) and A[a][b] and A[b][c])
                if cnt > 0:
                    j2_count += 1
                a2_sum += cnt

        # NT pairs with A²>0 vs total NT pairs vs 3*c3
        nt_pairs = sum(1 for a in range(n) for c in range(n) if a != c and not A[a][c])
        # 3*c3 = number of ordered 3-cycle triples
        ordered_c3 = 6 * c3  # 6 orderings per 3-cycle, but...
        # Actually 3*c3 = number of (directed) 3-cycle arcs? No.
        # 3*c3: each 3-cycle has 3 arcs, each arc is an NT pair
        # So 3*c3 counts the number of NT ARCS that are IN a 3-cycle
        # But there could be NT arcs not in any 3-cycle... no.
        # If c→a, then (a,c) is NT. Is (a,c) always in a 3-cycle?
        # Need b with a→b→c→a, i.e., A[a][b]=1 and A[b][c]=1
        # This is A²[a,c] > 0. So if A²[a,c]=0, (a,c) is not in any 3-cycle.
        # So j2_count = #{NT arcs in at least one 3-cycle}

        j2_formulas.append((bits, j2_count, 3*c3, nt_pairs, a2_sum))

    # Check j2 = 3*c3?
    match_3c3 = all(j == 3*c for _, j, c, _, _ in j2_formulas)
    print(f"\nn={n}: J₂ = 3*c₃: {match_3c3}")

    if not match_3c3:
        # When does j2 ≠ 3*c3?
        mismatches = [(b, j, c) for b, j, c, _, _ in j2_formulas if j != c]
        print(f"  Mismatches: {len(mismatches)}")
        # Show a few
        for b, j, c in mismatches[:5]:
            A = build_adj(n, b)
            scores = tuple(sorted(sum(A[i][jj] for jj in range(n) if jj!=i) for i in range(n)))
            c3_val = c3_count(A, n)
            print(f"    bits={b}, scores={scores}, c3={c3_val}, J₂={j}, 3*c3={c}")

        # Better formula?
        # J₂ = #{NT arcs (a,c) with A²[a,c]>0}
        # 3*c3 = #{NT arcs in a 3-cycle} ≤ J₂ always?
        # Actually if A²[a,c]>0 then ∃b with a→b→c and c→a, so (a,b,c) is a 3-cycle
        # So J₂ = #{NT arcs in a 3-cycle}
        # And 3*c3 = #{NT arcs in a 3-cycle} since each 3-cycle has 3 NT arcs?
        # No! A 3-cycle (a,b,c) with a→b→c→a has arcs a→b, b→c, c→a
        # The NT arcs are: (a,c) with c→a, (b,a) with a→b, (c,b) with b→c
        # Wait, NT means the PAIR goes against the tournament direction
        # Actually (a,c) is an "NT pair" meaning the path (a,*,c) passes through b
        # where c→a is the "bad" direction
        # Each 3-cycle {a,b,c} contributes 3 NT pairs: (a,c), (b,a), (c,b)
        # where c→a, a→b, b→c
        # But the same NT pair (a,c) with c→a could be in MULTIPLE 3-cycles!
        # If ∃b1,b2 both with a→b→c, then (a,c) is in 3-cycles through b1 AND b2
        # The NUMBER of 3-cycles containing (a,c) = A²[a,c]
        # So 3*c3 = Σ_{NT arcs (a,c)} A²[a,c]
        # And J₂ = #{NT arcs with A²[a,c]>0}
        # These differ when A²[a,c] > 1 for some NT arcs

        # Verify: 3*c3 = Σ A²[a,c] over NT arcs
        test_3c3 = all(a == 3*c for _, _, c, _, a in j2_formulas)
        print(f"  3*c₃ = Σ A²[a,c] over NT pairs: {test_3c3}")

        # So J₂ < 3*c3 when some NT arcs are in multiple 3-cycles
        under = sum(1 for _, j, c, _, _ in j2_formulas if j < c)
        over = sum(1 for _, j, c, _, _ in j2_formulas if j > c)
        print(f"  J₂ < 3*c₃: {under}, J₂ > 3*c₃: {over}")


print("\nDone.")
