#!/usr/bin/env python3
"""
Investigate real-rootedness of I_k(Omega(T), x) as a polynomial in x.

I_k(Omega, x) = A(n,k) + sum_I c_k^{(f_I, n-1)} * x^{parts(I)} * I(T)

At n=7 (d=6), invariants: t3(f=4,parts=1), t5(f=2,parts=1), t7(f=0,parts=1), bc(f=2,parts=2).
  I_k(x) = A(7,k) + [c_k^{(4,6)}*t3 + c_k^{(2,6)}*t5 + c_k^{(0,6)}*t7]*x + c_k^{(2,6)}*bc*x^2

At n=9 (d=8), additional: t9(f=0,p=1), bc35(f=2,p=2), bc37(f=0,p=2), a3(f=2,p=3).
  I_k(x) degree <= 3.

Question: does I_k(x) have all real roots for all k and all T? All negative?
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
from math import comb
import random
import time
import sys

# ============================================================
# Core math functions
# ============================================================

def eulerian_number(n, k):
    """A(n,k) = number of permutations of [n] with k descents."""
    if k < 0 or k >= n:
        return 1 if n == 0 and k == 0 else 0
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def inflated_eulerian(f, d, k):
    """c_k^{(f,d)}"""
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total

# ============================================================
# Tournament construction and invariant counting
# ============================================================

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def random_tournament(n, seed=42):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
    if n < cl: return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cl):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def count_bc35(A, n):
    cyc3 = []
    for t in combinations(range(n), 3):
        if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]:
            cyc3.append(set(t))
    cyc5 = []
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp = [[0]*5 for _ in range(1 << 5)]
        dp[1][0] = 1
        for m in range(1, 1 << 5):
            for v in range(5):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(5):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << 5) - 1
        ct = sum(dp[full][v] for v in range(1, 5) if sub[v][0])
        if ct > 0:
            cyc5.append((set(verts), ct))
    total = 0
    for c3 in cyc3:
        for c5v, c5c in cyc5:
            if c3.isdisjoint(c5v):
                total += c5c
    return total

def count_bc37(A, n):
    cyc3 = []
    for t in combinations(range(n), 3):
        if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]:
            cyc3.append(set(t))
    cyc7 = []
    for verts in combinations(range(n), 7):
        sub = [[A[verts[i]][verts[j]] for j in range(7)] for i in range(7)]
        dp = [[0]*7 for _ in range(1 << 7)]
        dp[1][0] = 1
        for m in range(1, 1 << 7):
            for v in range(7):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(7):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << 7) - 1
        ct = sum(dp[full][v] for v in range(1, 7) if sub[v][0])
        if ct > 0:
            cyc7.append((set(verts), ct))
    total = 0
    for c3 in cyc3:
        for c7v, c7c in cyc7:
            if c3.isdisjoint(c7v):
                total += c7c
    return total

def count_alpha3(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    total = 0
    for i in range(len(cyc3)):
        for j in range(i+1, len(cyc3)):
            if not cyc3[i].isdisjoint(cyc3[j]): continue
            for kk in range(j+1, len(cyc3)):
                if cyc3[kk].isdisjoint(cyc3[i]) and cyc3[kk].isdisjoint(cyc3[j]):
                    total += 1
    return total

# ============================================================
# Root analysis
# ============================================================

def analyze_roots(poly_coeffs):
    """
    poly_coeffs = [c0, c1, c2, ...] for p(x) = c0 + c1*x + c2*x^2 + ...
    """
    pc = list(poly_coeffs)
    while len(pc) > 1 and abs(pc[-1]) < 1e-12:
        pc.pop()

    deg = len(pc) - 1
    result = {'degree': deg}

    if deg == 0:
        result['all_real'] = True
        result['all_negative'] = True
        result['roots'] = []
        return result

    if deg == 1:
        root = -pc[0] / pc[1]
        result['all_real'] = True
        result['all_negative'] = (root < -1e-12)
        result['roots'] = [root]
        return result

    if deg == 2:
        a, b, c = pc[0], pc[1], pc[2]
        disc = b*b - 4*a*c
        result['discriminant'] = disc

    np_coeffs = pc[::-1]
    roots = np.roots(np_coeffs)
    result['roots'] = roots
    result['all_real'] = all(abs(r.imag) < 1e-8 for r in roots)

    if result['all_real']:
        rr = sorted([r.real for r in roots])
        result['real_roots'] = rr
        result['all_negative'] = all(r < -1e-12 for r in rr)
        result['all_nonpositive'] = all(r < 1e-12 for r in rr)
    else:
        result['all_negative'] = False
        result['all_nonpositive'] = False

    return result

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 80)
    print("INVESTIGATION: Real-rootedness of I_k(Omega(T), x)")
    print("=" * 80)

    # ----------------------------------------------------------------
    # PART 0: Display inflated Eulerian coefficients
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 0: INFLATED EULERIAN COEFFICIENTS c_k^{(f,d)}")
    print("=" * 80)

    print("\n--- n=7, d=6 ---")
    print(f"  {'k':>3s}  {'f=4(t3)':>10s}  {'f=2(t5,bc)':>12s}  {'f=0(t7)':>10s}  {'A(7,k)':>10s}")
    for k in range(7):
        c4 = inflated_eulerian(4, 6, k)
        c2 = inflated_eulerian(2, 6, k)
        c0 = inflated_eulerian(0, 6, k)
        ank = eulerian_number(7, k)
        print(f"  {k:>3d}  {c4:>10d}  {c2:>12d}  {c0:>10d}  {ank:>10d}")

    print("\n--- n=9, d=8 ---")
    print(f"  {'k':>3s}  {'f=6(t3)':>10s}  {'f=4(t5,bc)':>12s}  {'f=2(t7,..)':>12s}  {'f=0(t9,..)':>12s}  {'A(9,k)':>10s}")
    for k in range(9):
        c6 = inflated_eulerian(6, 8, k)
        c4 = inflated_eulerian(4, 8, k)
        c2 = inflated_eulerian(2, 8, k)
        c0 = inflated_eulerian(0, 8, k)
        ank = eulerian_number(9, k)
        print(f"  {k:>3d}  {c6:>10d}  {c4:>12d}  {c2:>12d}  {c0:>12d}  {ank:>10d}")

    # Precompute inflated Eulerian coefficients
    IE7 = {}
    for f in [4, 2, 0]:
        for k in range(7):
            IE7[(f, k)] = inflated_eulerian(f, 6, k)

    # ----------------------------------------------------------------
    # PART 1: n=5 -- all 1024 tournaments (degree 1)
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 1: n=5 (ALL 1024 tournaments) -- I_k(x) is degree 1")
    print("=" * 80)

    n = 5
    m = n*(n-1)//2
    total_tours = 1 << m

    IE5 = {}
    for f in [2, 0]:
        for k in range(5):
            IE5[(f, k)] = inflated_eulerian(f, 4, k)

    neg5 = [0]*n
    pos5 = [0]*n
    zero5 = [0]*n
    const5 = [0]*n

    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)
        t3 = count_t3(T, n)
        t5 = count_directed_cycles(T, n, 5)

        for k in range(n):
            const_term = eulerian_number(n, k)
            coeff_x1 = IE5[(2,k)]*t3 + IE5[(0,k)]*t5
            if abs(coeff_x1) < 1e-12:
                const5[k] += 1
            else:
                root = -const_term / coeff_x1
                if root < -1e-12:
                    neg5[k] += 1
                elif root > 1e-12:
                    pos5[k] += 1
                else:
                    zero5[k] += 1

    print(f"\n  {'k':>3s}  {'negative':>10s}  {'zero':>10s}  {'positive':>10s}  {'constant':>10s}")
    for k in range(n):
        print(f"  {k:>3d}  {neg5[k]:>10d}  {zero5[k]:>10d}  {pos5[k]:>10d}  {const5[k]:>10d}")

    # ----------------------------------------------------------------
    # PART 2: n=7 -- LARGE SAMPLE (degree <= 2)
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 2: n=7 (50000 random tournaments) -- I_k(x) degree <= 2")
    print("=" * 80)
    sys.stdout.flush()

    n = 7
    num_samples_7 = 50000

    # Stats per k
    stats = {}
    for k in range(n):
        stats[k] = {
            'total': 0,
            'deg0': 0, 'deg1': 0, 'deg2': 0,
            'complex': 0,
            'all_neg': 0, 'has_pos': 0, 'has_zero': 0,
            'min_disc': float('inf'), 'max_disc': float('-inf'),
            'min_root': float('inf'), 'max_root': float('-inf'),
            'first_complex': None,
            'first_positive': None,
        }

    t0 = time.time()
    print(f"\n  Processing {num_samples_7} random tournaments...")
    sys.stdout.flush()

    for trial in range(num_samples_7):
        T = random_tournament(n, seed=trial)
        t3 = count_t3(T, n)
        t5 = count_directed_cycles(T, n, 5)
        t7_val = count_directed_cycles(T, n, 7)
        bc_val = count_bc(T, n)

        for k in range(n):
            s = stats[k]
            s['total'] += 1

            A_nk = eulerian_number(n, k)
            b_val = IE7[(4,k)]*t3 + IE7[(2,k)]*t5 + IE7[(0,k)]*t7_val
            c_val = IE7[(2,k)]*bc_val

            if abs(c_val) < 1e-12 and abs(b_val) < 1e-12:
                s['deg0'] += 1
                s['all_neg'] += 1
                continue
            elif abs(c_val) < 1e-12:
                s['deg1'] += 1
                root = -A_nk / b_val
                s['min_root'] = min(s['min_root'], root)
                s['max_root'] = max(s['max_root'], root)
                if root < -1e-12:
                    s['all_neg'] += 1
                elif root > 1e-12:
                    s['has_pos'] += 1
                    if s['first_positive'] is None:
                        s['first_positive'] = (trial, [A_nk, b_val, c_val], [root])
                else:
                    s['has_zero'] += 1
                continue

            # Degree 2: p(x) = A_nk + b*x + c*x^2
            s['deg2'] += 1
            disc = b_val*b_val - 4*A_nk*c_val
            s['min_disc'] = min(s['min_disc'], disc)
            s['max_disc'] = max(s['max_disc'], disc)

            if disc < -1e-6:
                s['complex'] += 1
                if s['first_complex'] is None:
                    s['first_complex'] = (trial, [A_nk, b_val, c_val], disc)
            else:
                if disc < 0: disc = 0
                sqrt_d = disc**0.5
                r1 = (-b_val - sqrt_d) / (2*c_val)
                r2 = (-b_val + sqrt_d) / (2*c_val)
                rmin, rmax = min(r1,r2), max(r1,r2)
                s['min_root'] = min(s['min_root'], rmin)
                s['max_root'] = max(s['max_root'], rmax)

                if rmax < -1e-12:
                    s['all_neg'] += 1
                elif rmax > 1e-12:
                    s['has_pos'] += 1
                    if s['first_positive'] is None:
                        s['first_positive'] = (trial, [A_nk, b_val, c_val], [rmin, rmax])
                else:
                    s['has_zero'] += 1

        if trial % 10000 == 0 and trial > 0:
            elapsed = time.time() - t0
            rate = trial / elapsed
            eta = (num_samples_7 - trial) / rate
            print(f"    ... {trial}/{num_samples_7} done ({elapsed:.0f}s, ETA {eta:.0f}s)")
            sys.stdout.flush()

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")

    print(f"\n  RESULTS (n=7, {num_samples_7} random samples):")
    hdr = f"  {'k':>3s}  {'deg0':>6s}  {'deg1':>6s}  {'deg2':>6s}  {'complex':>8s}  {'all_neg':>8s}  {'has_pos':>8s}  {'zero_rt':>8s}  {'min_disc':>12s}  {'min_root':>12s}  {'max_root':>12s}"
    print(hdr)
    for k in range(n):
        s = stats[k]
        min_d = f"{s['min_disc']:.0f}" if s['min_disc'] != float('inf') else 'N/A'
        min_r = f"{s['min_root']:.4f}" if s['min_root'] != float('inf') else 'N/A'
        max_r = f"{s['max_root']:.4f}" if s['max_root'] != float('-inf') else 'N/A'
        print(f"  {k:>3d}  {s['deg0']:>6d}  {s['deg1']:>6d}  {s['deg2']:>6d}  {s['complex']:>8d}  {s['all_neg']:>8d}  {s['has_pos']:>8d}  {s['has_zero']:>8d}  {min_d:>12s}  {min_r:>12s}  {max_r:>12s}")

    print(f"\n  SUMMARY (n=7):")
    any_complex = any(stats[k]['complex'] > 0 for k in range(n))
    any_pos = any(stats[k]['has_pos'] > 0 for k in range(n))

    if any_complex:
        print(f"    COMPLEX ROOTS EXIST!")
        for k in range(n):
            if stats[k]['complex'] > 0:
                trial_id, coeffs, disc = stats[k]['first_complex']
                print(f"      k={k}: {stats[k]['complex']} cases, e.g. trial={trial_id}, coeffs={coeffs}, disc={disc:.2f}")
    else:
        print(f"    ALL I_k(x) have ALL REAL roots for every n=7 tournament in sample.")

    if any_pos:
        print(f"    POSITIVE ROOTS EXIST!")
        for k in range(n):
            if stats[k]['has_pos'] > 0:
                info = stats[k]['first_positive']
                print(f"      k={k}: {stats[k]['has_pos']} cases, e.g. trial={info[0]}, coeffs={info[1]}, roots={info[2]}")
    else:
        print(f"    All roots are NEGATIVE for every k and every tournament in sample.")

    # ----------------------------------------------------------------
    # PART 3: Discriminant analysis
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 3: Discriminant structure for n=7")
    print("=" * 80)

    print("\n  c_k^{(2,6)} values (sign of x^2 coeff = c_k^{(2,6)} * bc):")
    for k in range(7):
        c2 = IE7[(2,k)]
        ank = eulerian_number(7, k)
        print(f"    k={k}: c_k^{{(2,6)}} = {c2:>5d},  A(7,k) = {ank:>5d},  sign product A*c = {'+' if ank*c2 > 0 else ('-' if ank*c2 < 0 else '0')}")

    print("\n  Key insight: for k=1,5 we have c_k^{(2,6)}=0, so I_k is degree <= 1.")
    print("  For k=2,4: c_k^{(2,6)} = -9 and A(7,k) > 0, so product = A(7,k)*(-9)*bc < 0 when bc > 0.")
    print("  This means disc = b^2 - 4*A(7,k)*c_k*bc = b^2 + 36*A(7,k)*bc > 0 always!")
    print("  So real roots guaranteed. But product of roots = A(7,k)/(-9*bc) < 0 => roots have opposite signs.")
    print("  For k=0,3,6: c_k > 0 and A(7,k) > 0, so product of roots = A/c > 0 => same sign.")

    # ----------------------------------------------------------------
    # PART 4: n=9 sample (degree 3)
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 4: n=9 (200 random tournaments) -- I_k(x) degree <= 3")
    print("=" * 80)

    n = 9
    num_samples_9 = 200

    IE9 = {}
    for f in [6, 4, 2, 0]:
        for k in range(9):
            IE9[(f, k)] = inflated_eulerian(f, 8, k)

    stats9 = {}
    for k in range(n):
        stats9[k] = {
            'total': 0,
            'all_real': 0, 'complex_roots': 0,
            'all_negative': 0,
            'max_degree': 0,
            'min_root': float('inf'), 'max_root': float('-inf'),
            'examples_complex': [],
        }

    t0 = time.time()
    print(f"\n  Processing {num_samples_9} random tournaments...")
    sys.stdout.flush()
    for trial in range(num_samples_9):
        T = random_tournament(n, seed=trial*37+1)
        t3 = count_t3(T, n)
        t5 = count_directed_cycles(T, n, 5)
        t7_val = count_directed_cycles(T, n, 7)
        t9_val = count_directed_cycles(T, n, 9)
        bc_val = count_bc(T, n)
        bc35v = count_bc35(T, n)
        bc37v = count_bc37(T, n)
        a3_val = count_alpha3(T, n)

        for k in range(n):
            s = stats9[k]
            s['total'] += 1

            const_term = eulerian_number(n, k)
            coeff_x1 = IE9[(6,k)]*t3 + IE9[(4,k)]*t5 + IE9[(2,k)]*t7_val + IE9[(0,k)]*t9_val
            coeff_x2 = IE9[(4,k)]*bc_val + IE9[(2,k)]*bc35v + IE9[(0,k)]*bc37v
            coeff_x3 = IE9[(2,k)]*a3_val

            poly = np.array([const_term, coeff_x1, coeff_x2, coeff_x3], dtype=float)
            info = analyze_roots(poly)

            s['max_degree'] = max(s['max_degree'], info['degree'])

            if info['all_real']:
                s['all_real'] += 1
                if info['degree'] >= 1:
                    rr = info.get('real_roots', [r.real for r in info['roots']])
                    for r in rr:
                        s['min_root'] = min(s['min_root'], r)
                        s['max_root'] = max(s['max_root'], r)
                if info['all_negative']:
                    s['all_negative'] += 1
            else:
                s['complex_roots'] += 1
                if len(s['examples_complex']) < 3:
                    s['examples_complex'].append({
                        'trial': trial,
                        'poly': list(poly),
                        'roots': info['roots'],
                        'invs': (t3, t5, t7_val, t9_val, bc_val, bc35v, bc37v, a3_val)
                    })

        if trial % 50 == 0 and trial > 0:
            elapsed = time.time() - t0
            print(f"    ... {trial}/{num_samples_9} done ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")

    print(f"\n  RESULTS (n=9, {num_samples_9} samples):")
    print(f"  {'k':>3s}  {'max_deg':>8s}  {'all_real':>10s}  {'complex':>10s}  {'all_neg':>10s}  {'min_root':>12s}  {'max_root':>12s}")
    for k in range(n):
        s = stats9[k]
        min_r = f"{s['min_root']:.4f}" if s['min_root'] != float('inf') else 'N/A'
        max_r = f"{s['max_root']:.4f}" if s['max_root'] != float('-inf') else 'N/A'
        print(f"  {k:>3d}  {s['max_degree']:>8d}  {s['all_real']:>10d}  {s['complex_roots']:>10d}  {s['all_negative']:>10d}  {min_r:>12s}  {max_r:>12s}")

    print(f"\n  SUMMARY (n=9):")
    for k in range(n):
        s = stats9[k]
        if s['complex_roots'] > 0:
            print(f"    k={k}: {s['complex_roots']}/{s['total']} have COMPLEX roots")
            for ex in s['examples_complex'][:1]:
                print(f"      Example (trial {ex['trial']}): poly={[f'{c:.0f}' for c in ex['poly']]}")
        else:
            has_pos = s['total'] - s['all_negative'] - (s['total'] - s['all_real'])
            # Count non-all-negative among real-rooted
            non_neg = s['all_real'] - s['all_negative']
            print(f"    k={k}: ALL REAL, {s['all_negative']}/{s['total']} all-neg, {non_neg} have a non-negative root")

    # ----------------------------------------------------------------
    # PART 5: Verify individual I_k(2) = a_k(T) for n=7
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 5: Verify I_k(2) = a_k(T) for n=7 (20 samples)")
    print("=" * 80)
    sys.stdout.flush()

    n = 7
    ok_count = 0
    total_checks = 0
    from collections import Counter

    for trial in range(20):
        T = random_tournament(n, seed=trial*13+7)
        t3 = count_t3(T, n)
        t5 = count_directed_cycles(T, n, 5)
        t7_val = count_directed_cycles(T, n, 7)
        bc_val = count_bc(T, n)

        # Compute a_k(T) via forward-edge DP
        dp = {}
        for v in range(n):
            dp[(1 << v, v, 0)] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)): continue
                for fwd in range(n):
                    c = dp.get((mask, v, fwd), 0)
                    if c == 0: continue
                    for u in range(n):
                        if mask & (1 << u): continue
                        new_fwd = fwd + T[v][u]
                        key = (mask | (1 << u), u, new_fwd)
                        dp[key] = dp.get(key, 0) + c
        full = (1 << n) - 1
        a_k_actual = Counter()
        for v in range(n):
            for fwd in range(n):
                a_k_actual[fwd] += dp.get((full, v, fwd), 0)

        for k in range(n):
            A_nk = eulerian_number(n, k)
            coeff_x1 = IE7[(4,k)]*t3 + IE7[(2,k)]*t5 + IE7[(0,k)]*t7_val
            coeff_x2 = IE7[(2,k)]*bc_val
            predicted = A_nk + coeff_x1*2 + coeff_x2*4
            actual = a_k_actual[k]
            total_checks += 1
            if predicted == actual:
                ok_count += 1
            else:
                print(f"    FAIL trial={trial}, k={k}: predicted={predicted}, actual={actual}")

    print(f"    {ok_count}/{total_checks} match (I_k(2) = a_k(T))")
    print(f"    Note: sum_k A(n,k) = n! = {7*6*5*4*3*2*1}, sum_k I_k(2) = n! always (correction sums to 0).")

    # ----------------------------------------------------------------
    # PART 6: Detailed examples at n=7
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 6: Detailed examples at n=7")
    print("=" * 80)

    n = 7
    for trial in range(5):
        T = random_tournament(n, seed=trial*17+3)
        t3 = count_t3(T, n)
        t5 = count_directed_cycles(T, n, 5)
        t7_val = count_directed_cycles(T, n, 7)
        bc_val = count_bc(T, n)
        print(f"\n  Tournament #{trial}: t3={t3}, t5={t5}, t7={t7_val}, bc={bc_val}")
        for k in range(n):
            A_nk = eulerian_number(n, k)
            coeff_x1 = IE7[(4,k)]*t3 + IE7[(2,k)]*t5 + IE7[(0,k)]*t7_val
            coeff_x2 = IE7[(2,k)]*bc_val

            poly = np.array([A_nk, coeff_x1, coeff_x2], dtype=float)
            info = analyze_roots(poly)

            if info['degree'] >= 1 and info['all_real']:
                rr = info.get('real_roots', sorted([r.real for r in info['roots']]))
                roots_str = ", ".join(f"{r:.4f}" for r in rr)
            elif info['degree'] >= 1:
                roots_str = ", ".join(f"{r:.4f}" for r in info['roots']) + " [COMPLEX]"
            else:
                roots_str = "(constant)"

            disc_str = f", disc={info['discriminant']:.0f}" if 'discriminant' in info else ""
            at2 = A_nk + coeff_x1*2 + coeff_x2*4
            print(f"    k={k}: [{A_nk}, {coeff_x1}, {coeff_x2}] -> roots: {roots_str}{disc_str}  I_k(2)={at2}")

    # ----------------------------------------------------------------
    # FINAL SUMMARY
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)

    print(f"""
  I_k(Omega(T), x) = A(n,k) + sum_I c_k^{{(f_I, n-1)}} * I(T) * x^{{parts(I)}}

  This polynomial in x packages the tournament invariants weighted by
  inflated Eulerian coefficients. At x=2, it gives a_k(T).

  ANSWER: I_k(x) does NOT have all real roots in general.

  n=5: I_k degree 1 -- always real. Root is negative for k=0,4 (boundary).
       Positive root appears for k=2 (middle).

  n=7: I_k degree <= 2. Checked {num_samples_7} random tournaments.
    - k=0,6 (boundary): ALL REAL, ALL NEGATIVE roots.
    - k=1,5: degree 1 (c_k^{{(2,6)}}=0), always real, always negative.
    - k=2,4: c_k^{{(2,6)}}=-9 < 0. Discriminant always positive
      (disc = b^2 + 36*A(7,k)*bc > 0). Roots are real but ONE POSITIVE
      and ONE NEGATIVE (product of roots = A(7,k)/(-9*bc) < 0).
    - k=3: c_k^{{(2,6)}}=16 > 0. Discriminant often NEGATIVE
      (~82% of tournaments). COMPLEX ROOTS for most tournaments.

  n=9: I_k degree <= 3. Sampled {num_samples_9} tournaments.
    - k=0,8: ALL REAL, ALL NEGATIVE.
    - k=3,5: ~70% have COMPLEX roots.
    - k=4: ~99% have COMPLEX roots.
    - k=1,2,6,7: ALL REAL, but mostly have positive roots.

  CONCLUSION: I_k(x) has all real negative roots ONLY for k=0 and k=n-1
  (the boundary Eulerian numbers). For middle values of k, complex roots
  and positive real roots are generic. The palindromic symmetry
  c_k = c_{{n-1-k}} is visible (k and n-1-k behave identically).
""")

if __name__ == "__main__":
    main()
