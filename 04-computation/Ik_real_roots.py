#!/usr/bin/env python3
"""
Investigate real-rootedness of I_k(Omega(T), x) as a polynomial in x.

BACKGROUND:
  a_k(T) = A(n,k) + sum_I 2^{parts(I)} * c_k^{(f_I, n-1)} * I(T)

  The independence polynomial evaluated at x gives:
    I(Omega(T), x) = sum_k a_k(T) * x^k   ... NO, this is the forward-edge GF.

  Actually, I(Omega(T), x) = 1 + (number of odd cycles)*x + (vertex-disjoint pairs)*x^2 + ...

  The question is: define I_k(Omega, x) as the polynomial in x obtained by
  replacing each invariant I(T) by x^{parts(I)} weighted by c_k^{(f_I, n-1)}:

    I_k(x) = A(n,k) + sum_I c_k^{(f_I, n-1)} * I(T) * x^{parts(I)}

  Wait -- re-reading the user's formula more carefully:
    I_k(Omega, x) = A(n,k) + sum_I c_k^{(f_I, n-1)} * x^{parts(I)} * I(T)

  At n=7: invariants are t3(f=4,parts=1), t5(f=2,parts=1), t7(f=0,parts=1), bc(f=2,parts=2).
  So:
    I_k(x) = A(7,k) + [c_k^{(4,6)}*t3 + c_k^{(2,6)}*t5 + c_k^{(0,6)}*t7]*x + c_k^{(2,6)}*bc*x^2

  This is degree 2 in x. Question: does it have all real roots? All negative?

  At n=9: invariants include a3 (parts=3), so degree 3 in x.

  Note: I_0 = I(Omega, x) since c_0^{(f,d)} * 2^parts gives the x=2 coefficient
  ... actually, let's be careful. The standard OCF is:
    H(T) = I(Omega, 2) = sum over vertex-disjoint odd-cycle collections S of 2^|S|

  So I(Omega, x) = 1 + (sum of all odd cycles)*x + (sum of VD pairs)*x^2 + ...

  And a_k(T) = A(n,k) + sum_I 2^{parts(I)} * c_k^{(f_I, n-1)} * I(T)

  The user defines I_k(x) by replacing 2^{parts} with x^{parts}:
    I_k(x) = A(n,k) + sum_I c_k^{(f_I, n-1)} * x^{parts(I)} * I(T)

  So a_k(T) = I_k(2).  And I_0(x) evaluated at x=2 gives H(T).

  But what IS I_0? We have c_0^{(f,d)} = inflated_eulerian(f, d, 0).
  For the a_k formula: a_0(T) = A(n,0) + sum_I 2^{parts} * c_0^{(f,d)} * I(T).
  So I_0(x) = A(n,0) + sum_I c_0^{(f,d)} * x^{parts} * I(T).

  And we know a_0 = 1 for all T (the identity permutation has 0 forward edges in
  a transitive tournament... wait, a_0 is number of permutations with 0 forward edges).
  Actually A(n,0)=1 and A(n,n-1)=1.

Computations:
  1. For ALL n=7 tournaments: compute I_k(x) for each k=0..6, check real-rootedness.
  2. For a sample of n=9 tournaments: compute I_k(x) for each k=0..8, check real-rootedness.

opus-2026-03-07-S... (investigation script)
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from math import comb
import random
import sys

# ============================================================
# Core functions
# ============================================================

def eulerian_number(n, k):
    """A(n,k) = number of permutations of [n] with k descents."""
    if k < 0 or k >= n:
        return 1 if (n == 0 and k == 0) else (0 if k != 0 else 1) if k == 0 and n >= 0 else 0
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def inflated_eulerian(f, d, k):
    """
    c_k^{(f,d)} = coefficient of p^k q^{d-k} in F_f(r) expressed in degree-d basis.
    """
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                T[i][j] = 1
            else:
                T[j][i] = 1
            k += 1
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
    """Count directed 3-cycles."""
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
    """Count directed cycles of length cl using Hamiltonian path DP on subsets."""
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
    """Pairs of vertex-disjoint 3-cycles."""
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def count_bc35(A, n):
    """Pairs: one 3-cycle and one disjoint 5-cycle."""
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
        for c5_verts, c5_ct in cyc5:
            if c3.isdisjoint(c5_verts):
                total += c5_ct
    return total

def count_bc37(A, n):
    """Pairs: one 3-cycle and one disjoint 7-cycle."""
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
        for c7_verts, c7_ct in cyc7:
            if c3.isdisjoint(c7_verts):
                total += c7_ct
    return total

def count_alpha3(A, n):
    """Triples of vertex-disjoint 3-cycles."""
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
# I_k(x) polynomial construction
# ============================================================

def compute_Ik_poly_n7(T):
    """
    For tournament T on 7 vertices, compute I_k(x) for k=0..6.

    I_k(x) = A(7,k) + [c_k^{(4,6)}*t3 + c_k^{(2,6)}*t5 + c_k^{(0,6)}*t7]*x + c_k^{(2,6)}*bc*x^2

    Wait -- the user's formula has:
      I_k(x) = A(n,k) + sum_I c_k^{(f_I, n-1)} * x^{parts(I)} * I(T)

    For n=7, d=6:
      parts=1 invariants: t3(f=4), t5(f=2), t7(f=0)
      parts=2 invariants: bc(f=2)

    So I_k(x) = A(7,k) + [c_k^{(4,6)}*t3 + c_k^{(2,6)}*t5 + c_k^{(0,6)}*t7]*x + c_k^{(2,6)}*bc*x^2

    Returns list of 7 polynomials (each as numpy array of coefficients, constant first).
    """
    n = 7
    d = 6
    t3 = count_t3(T, n)
    t5 = count_directed_cycles(T, n, 5)
    t7 = count_directed_cycles(T, n, 7)
    bc = count_bc(T, n)

    polys = []
    for k in range(n):
        c_t3 = inflated_eulerian(4, 6, k)
        c_t5 = inflated_eulerian(2, 6, k)
        c_t7 = inflated_eulerian(0, 6, k)
        c_bc = inflated_eulerian(2, 6, k)  # bc has f=2 too

        # Constant term
        const = eulerian_number(n, k)
        # x^1 coefficient
        coeff_x1 = c_t3 * t3 + c_t5 * t5 + c_t7 * t7
        # x^2 coefficient
        coeff_x2 = c_bc * bc

        polys.append(np.array([const, coeff_x1, coeff_x2], dtype=float))

    return polys, (t3, t5, t7, bc)

def compute_Ik_poly_n5(T):
    """
    For n=5, d=4. Invariants: t3(f=2, parts=1), t5(f=0, parts=1).
    I_k(x) = A(5,k) + [c_k^{(2,4)}*t3 + c_k^{(0,4)}*t5]*x
    Degree 1 in x -- always has a real root.
    """
    n = 5
    d = 4
    t3 = count_t3(T, n)
    t5 = count_directed_cycles(T, n, 5)

    polys = []
    for k in range(n):
        c_t3 = inflated_eulerian(2, 4, k)
        c_t5 = inflated_eulerian(0, 4, k)
        const = eulerian_number(n, k)
        coeff_x1 = c_t3 * t3 + c_t5 * t5
        polys.append(np.array([const, coeff_x1], dtype=float))

    return polys, (t3, t5)

def compute_Ik_poly_n9(T):
    """
    For n=9, d=8. Invariants:
      parts=1: t3(f=6), t5(f=4), t7(f=2), t9(f=0)
      parts=2: bc(f=4), bc35(f=2), bc37(f=0)
      parts=3: a3(f=2)

    I_k(x) = A(9,k) + [c_k^{(6,8)}*t3 + c_k^{(4,8)}*t5 + c_k^{(2,8)}*t7 + c_k^{(0,8)}*t9]*x
                     + [c_k^{(4,8)}*bc + c_k^{(2,8)}*bc35 + c_k^{(0,8)}*bc37]*x^2
                     + c_k^{(2,8)}*a3*x^3
    """
    n = 9
    d = 8
    t3 = count_t3(T, n)
    t5 = count_directed_cycles(T, n, 5)
    t7 = count_directed_cycles(T, n, 7)
    t9 = count_directed_cycles(T, n, 9)
    bc = count_bc(T, n)
    bc35_val = count_bc35(T, n)
    bc37_val = count_bc37(T, n)
    a3 = count_alpha3(T, n)

    polys = []
    for k in range(n):
        c_f6 = inflated_eulerian(6, 8, k)
        c_f4 = inflated_eulerian(4, 8, k)
        c_f2 = inflated_eulerian(2, 8, k)
        c_f0 = inflated_eulerian(0, 8, k)

        const = eulerian_number(n, k)
        coeff_x1 = c_f6*t3 + c_f4*t5 + c_f2*t7 + c_f0*t9
        coeff_x2 = c_f4*bc + c_f2*bc35_val + c_f0*bc37_val
        coeff_x3 = c_f2*a3

        polys.append(np.array([const, coeff_x1, coeff_x2, coeff_x3], dtype=float))

    return polys, (t3, t5, t7, t9, bc, bc35_val, bc37_val, a3)

def analyze_roots(poly_coeffs):
    """
    Analyze roots of polynomial given as [c0, c1, c2, ...] where p(x) = c0 + c1*x + c2*x^2 + ...
    Returns dict with: degree, roots, all_real, all_negative, discriminant (for degree 2).
    """
    # Strip trailing zeros
    while len(poly_coeffs) > 1 and abs(poly_coeffs[-1]) < 1e-12:
        poly_coeffs = poly_coeffs[:-1]

    deg = len(poly_coeffs) - 1
    result = {'degree': deg, 'coeffs': poly_coeffs.copy()}

    if deg == 0:
        result['roots'] = np.array([])
        result['all_real'] = True
        result['all_negative'] = True
        result['nonneg_for_x_ge_0'] = (poly_coeffs[0] >= 0)
        return result

    if deg == 1:
        root = -poly_coeffs[0] / poly_coeffs[1]
        result['roots'] = np.array([root])
        result['all_real'] = True
        result['all_negative'] = (root < 0)
        result['nonneg_for_x_ge_0'] = (poly_coeffs[0] >= 0 and poly_coeffs[1] >= 0) or \
                                       (poly_coeffs[0] >= 0 and root <= 0)
        return result

    if deg == 2:
        a, b, c = poly_coeffs[0], poly_coeffs[1], poly_coeffs[2]
        # p(x) = a + bx + cx^2.  In standard form: cx^2 + bx + a
        disc = b*b - 4*a*c
        result['discriminant'] = disc

    # Use numpy to find roots: np.roots expects [leading, ..., constant]
    np_coeffs = poly_coeffs[::-1]  # reverse to get [leading, ..., constant]
    roots = np.roots(np_coeffs)
    result['roots'] = roots

    # Check if all real
    result['all_real'] = all(abs(r.imag) < 1e-8 for r in roots)

    if result['all_real']:
        real_roots = sorted([r.real for r in roots])
        result['real_roots'] = real_roots
        result['all_negative'] = all(r < -1e-12 for r in real_roots)
        result['all_nonpositive'] = all(r < 1e-12 for r in real_roots)
        # Check if polynomial is nonneg for x >= 0
        # If leading coeff > 0 and all roots are negative, then p(x) > 0 for x > 0
        leading = poly_coeffs[-1]
        if result['all_negative']:
            result['nonneg_for_x_ge_0'] = (leading > 0 and poly_coeffs[0] > 0) or \
                                           (leading > 0)  # all neg roots, leading > 0 => p(x)>0 for x>0
        else:
            # Evaluate at a few points
            result['nonneg_for_x_ge_0'] = None  # can't easily determine
    else:
        result['all_negative'] = False
        result['all_nonpositive'] = False
        result['nonneg_for_x_ge_0'] = None

    return result

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 80)
    print("INVESTIGATION: Real-rootedness of I_k(Omega(T), x)")
    print("=" * 80)

    # ----------------------------------------------------------------
    # PART 0: Display the inflated Eulerian coefficients
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("INFLATED EULERIAN COEFFICIENTS c_k^{(f,d)}")
    print("=" * 80)

    print("\n--- n=7, d=6 ---")
    print(f"  {'k':>3s}", end="")
    for f_label in ['f=4(t3)', 'f=2(t5,bc)', 'f=0(t7)']:
        print(f"  {f_label:>12s}", end="")
    print()
    for k in range(7):
        c4 = inflated_eulerian(4, 6, k)
        c2 = inflated_eulerian(2, 6, k)
        c0 = inflated_eulerian(0, 6, k)
        print(f"  {k:>3d}  {c4:>12d}  {c2:>12d}  {c0:>12d}")

    # ----------------------------------------------------------------
    # PART 1: n=5 -- all 2^10 = 1024 tournaments (degree 1)
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 1: n=5 (ALL 1024 tournaments) — I_k(x) is degree 1")
    print("=" * 80)

    n = 5
    m = n*(n-1)//2
    total_tours = 1 << m

    neg_root_counts = defaultdict(int)  # k -> count of tournaments with negative root
    pos_root_counts = defaultdict(int)
    zero_root_counts = defaultdict(int)
    constant_counts = defaultdict(int)  # degree 0

    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)
        polys, _ = compute_Ik_poly_n5(T)
        for k in range(n):
            info = analyze_roots(polys[k])
            if info['degree'] == 0:
                constant_counts[k] += 1
            elif info['all_negative']:
                neg_root_counts[k] += 1
            elif info['roots'][0] > 1e-12:
                pos_root_counts[k] += 1
            else:
                zero_root_counts[k] += 1

    print(f"\n  {'k':>3s}  {'negative':>10s}  {'zero':>10s}  {'positive':>10s}  {'constant':>10s}")
    for k in range(n):
        print(f"  {k:>3d}  {neg_root_counts[k]:>10d}  {zero_root_counts[k]:>10d}  {pos_root_counts[k]:>10d}  {constant_counts[k]:>10d}")

    # ----------------------------------------------------------------
    # PART 2: n=7 -- all 2^21 = 2097152 tournaments (degree 2)
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 2: n=7 (ALL 2097152 tournaments) — I_k(x) is degree <= 2")
    print("=" * 80)

    n = 7
    m = n*(n-1)//2  # = 21
    total_tours = 1 << m

    # For each k, track:
    stats = {}
    for k in range(n):
        stats[k] = {
            'total': 0,
            'deg0': 0, 'deg1': 0, 'deg2': 0,
            'all_real': 0, 'complex_roots': 0,
            'all_negative': 0, 'has_positive_root': 0,
            'all_nonpositive': 0,
            'min_disc': float('inf'), 'max_disc': float('-inf'),
            'min_root': float('inf'), 'max_root': float('-inf'),
            'nonneg_for_x_ge_0': 0,
            'Ik_at_2_values': defaultdict(int),  # track I_k(2) distribution
        }

    print(f"\n  Processing {total_tours} tournaments...")
    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)
        polys, invs = compute_Ik_poly_n7(T)

        for k in range(n):
            s = stats[k]
            s['total'] += 1
            p = polys[k]
            info = analyze_roots(p)

            deg = info['degree']
            if deg == 0: s['deg0'] += 1
            elif deg == 1: s['deg1'] += 1
            else: s['deg2'] += 1

            if info['all_real']:
                s['all_real'] += 1
                if deg >= 1:
                    rr = info.get('real_roots', [r.real for r in info['roots']])
                    for r in rr:
                        s['min_root'] = min(s['min_root'], r)
                        s['max_root'] = max(s['max_root'], r)
                if info['all_negative']:
                    s['all_negative'] += 1
                if info.get('all_nonpositive', False):
                    s['all_nonpositive'] += 1
                if info.get('has_positive_root', False):
                    s['has_positive_root'] += 1
            else:
                s['complex_roots'] += 1

            if deg == 2 and 'discriminant' in info:
                s['min_disc'] = min(s['min_disc'], info['discriminant'])
                s['max_disc'] = max(s['max_disc'], info['discriminant'])

            # Evaluate I_k(2)
            val_at_2 = sum(p[j] * 2**j for j in range(len(p)))
            s['Ik_at_2_values'][int(val_at_2)] += 1

        if bits % 500000 == 0 and bits > 0:
            print(f"    ... {bits}/{total_tours} done")

    print(f"\n  RESULTS (n=7):")
    print(f"  {'k':>3s}  {'deg0':>6s}  {'deg1':>6s}  {'deg2':>6s}  {'all_real':>10s}  {'complex':>10s}  {'all_neg':>10s}  {'all_nonpos':>10s}  {'min_disc':>12s}  {'min_root':>12s}  {'max_root':>12s}")
    for k in range(n):
        s = stats[k]
        min_d = s['min_disc'] if s['min_disc'] != float('inf') else 'N/A'
        max_d = s['max_disc'] if s['max_disc'] != float('-inf') else 'N/A'
        min_r = f"{s['min_root']:.4f}" if s['min_root'] != float('inf') else 'N/A'
        max_r = f"{s['max_root']:.4f}" if s['max_root'] != float('-inf') else 'N/A'
        min_d_str = f"{min_d:.1f}" if isinstance(min_d, float) else min_d
        print(f"  {k:>3d}  {s['deg0']:>6d}  {s['deg1']:>6d}  {s['deg2']:>6d}  {s['all_real']:>10d}  {s['complex_roots']:>10d}  {s['all_negative']:>10d}  {s['all_nonpositive']:>10d}  {min_d_str:>12s}  {min_r:>12s}  {max_r:>12s}")

    # Summary
    print(f"\n  SUMMARY (n=7):")
    all_real_all_k = True
    all_neg_all_k = True
    for k in range(n):
        s = stats[k]
        has_nontrivial = s['deg1'] + s['deg2']
        if s['complex_roots'] > 0:
            all_real_all_k = False
            print(f"    k={k}: {s['complex_roots']} tournaments have COMPLEX roots!")
        if has_nontrivial > 0 and s['all_negative'] < has_nontrivial:
            not_neg = has_nontrivial - s['all_negative']
            if s['deg0'] + s['all_negative'] + s['all_nonpositive'] < s['total']:
                pass  # some have positive roots

    print(f"    All I_k(x) have all real roots for every T? {'YES' if all_real_all_k else 'NO'}")

    # Check: for every tournament, is I_k(2) = a_k(T)?
    # I_k(2) should equal a_k(T) by construction.
    print(f"\n  Verification: I_k(2) should give a_k(T).")
    print(f"    Checking a few distinct I_k(2) values per k:")
    for k in range(n):
        vals = stats[k]['Ik_at_2_values']
        sorted_vals = sorted(vals.items())
        print(f"    k={k}: {len(sorted_vals)} distinct values, range [{sorted_vals[0][0]}, {sorted_vals[-1][0]}]")

    # ----------------------------------------------------------------
    # PART 3: Detailed analysis of discriminant for n=7
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 3: Discriminant analysis for degree-2 cases at n=7")
    print("=" * 80)

    # For degree 2: p(x) = a + bx + cx^2
    # disc = b^2 - 4ac
    # Real roots iff disc >= 0
    # Both roots negative iff: disc >= 0, -b/2c < 0 (vertex < 0), and a/c > 0 (product of roots > 0)

    # The discriminant is:
    # b = c_k^{(4,6)}*t3 + c_k^{(2,6)}*t5 + c_k^{(0,6)}*t7
    # c = c_k^{(2,6)}*bc
    # a = A(7,k)
    # disc = b^2 - 4*A(7,k)*c_k^{(2,6)}*bc

    # Note: c (the x^2 coeff) = c_k^{(2,6)}*bc. Since bc >= 0, the sign of c
    # depends on c_k^{(2,6)}.

    print("\n  c_k^{(2,6)} values (determines sign of x^2 coefficient):")
    for k in range(7):
        c2 = inflated_eulerian(2, 6, k)
        ank = eulerian_number(7, k)
        print(f"    k={k}: c_k^{{(2,6)}} = {c2:>5d},  A(7,k) = {ank:>5d}")

    # ----------------------------------------------------------------
    # PART 4: n=9 sample (degree 3 in x)
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 4: n=9 (sample of 500 random tournaments) — I_k(x) is degree <= 3")
    print("=" * 80)

    n = 9
    num_samples = 500

    stats9 = {}
    for k in range(n):
        stats9[k] = {
            'total': 0,
            'all_real': 0, 'complex_roots': 0,
            'all_negative': 0,
            'max_degree': 0,
            'min_root': float('inf'), 'max_root': float('-inf'),
        }

    print(f"\n  Processing {num_samples} random tournaments...")
    for trial in range(num_samples):
        T = random_tournament(n, seed=trial*37+1)
        polys, invs = compute_Ik_poly_n9(T)

        for k in range(n):
            s = stats9[k]
            s['total'] += 1
            p = polys[k]
            info = analyze_roots(p)

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

        if trial % 100 == 0 and trial > 0:
            print(f"    ... {trial}/{num_samples} done")

    print(f"\n  RESULTS (n=9, {num_samples} samples):")
    print(f"  {'k':>3s}  {'max_deg':>8s}  {'all_real':>10s}  {'complex':>10s}  {'all_neg':>10s}  {'min_root':>12s}  {'max_root':>12s}")
    for k in range(n):
        s = stats9[k]
        min_r = f"{s['min_root']:.4f}" if s['min_root'] != float('inf') else 'N/A'
        max_r = f"{s['max_root']:.4f}" if s['max_root'] != float('-inf') else 'N/A'
        print(f"  {k:>3d}  {s['max_degree']:>8d}  {s['all_real']:>10d}  {s['complex_roots']:>10d}  {s['all_negative']:>10d}  {min_r:>12s}  {max_r:>12s}")

    # Summary for n=9
    print(f"\n  SUMMARY (n=9):")
    for k in range(n):
        s = stats9[k]
        if s['complex_roots'] > 0:
            print(f"    k={k}: {s['complex_roots']}/{s['total']} tournaments have COMPLEX roots")
        else:
            neg_frac = s['all_negative'] / max(1, s['all_real'])
            print(f"    k={k}: ALL REAL roots, {s['all_negative']}/{s['total']} have all-negative roots ({neg_frac:.1%})")

    # ----------------------------------------------------------------
    # PART 5: Check the standard I_0(x) = I(Omega, x) specifically
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 5: Standard I(Omega, x) = I_0(Omega, x) check")
    print("=" * 80)

    # I_0(x) should be the standard independence polynomial of Omega(T).
    # At n=7: I_0(x) = A(7,0) + [c_0^{(4,6)}*t3 + c_0^{(2,6)}*t5 + c_0^{(0,6)}*t7]*x + c_0^{(2,6)}*bc*x^2
    # But A(7,0) = 1, c_0^{(f,6)} = ... let's check

    print("\n  I_0 coefficients (n=7):")
    for f in [4, 2, 0]:
        print(f"    c_0^{{({f},6)}} = {inflated_eulerian(f, 6, 0)}")
    print(f"    A(7,0) = {eulerian_number(7, 0)}")

    # The standard independence polynomial is I(Omega, x) = 1 + |Omega|*x + (VD pairs)*x^2 + ...
    # In the deformed Eulerian framework, a_0(T) = I_0(2) should be 1 for transitive tournament.
    # But the OCF says H(T) = I(Omega, 2) = sum_k a_k(T).
    # Wait: a_0(T) for the identity permutation... A(n,0) = 1.
    # And a_0 = I_0(2) = 1 + 2*(c_0^{(4,6)}*t3 + ...) + 4*c_0^{(2,6)}*bc.

    # Actually, I_0(x) is NOT the standard I(Omega, x).
    # The standard I(Omega, x) has I(Omega, 2) = H(T) = sum_k a_k(T).
    # But I_0(2) = a_0(T), which is just one coefficient.
    #
    # The STANDARD I(Omega, x) = sum_k I_k(x) evaluated... no.
    # Actually, I(Omega, x) = 1 + t3*x + t5*x + t7*x + bc*x^2 + ...  (each cycle contributes x)
    # Hmm, that's also not right. I(Omega, x) = sum over independent sets S of x^|S|.

    # The user's I_k(x) is a DIFFERENT object from I(Omega, x).
    # It's the deformed Eulerian generating function with variable x replacing 2^{parts}.

    # So let's verify: sum_k I_k(x) should give a polynomial whose value at x=2 is H(T).
    print("\n  Checking: sum_k I_k(2) = H(T) for a sample of n=7 tournaments...")
    n = 7
    ok_count = 0
    for trial in range(100):
        T = random_tournament(n, seed=trial*13+7)
        polys, _ = compute_Ik_poly_n7(T)
        total = sum(sum(p[j] * 2**j for j in range(len(p))) for p in polys)

        # Compute H(T) directly
        from math import factorial
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        full = (1 << n) - 1
        for mask in range(1, 1 << n):
            for v in range(n):
                cnt = dp[mask][v]
                if not (mask & (1 << v)) or cnt == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    if T[v][u]: dp[mask | (1 << u)][u] += cnt
        H = sum(dp[full][v] for v in range(n))

        if abs(total - H) < 0.5:
            ok_count += 1
    print(f"    {ok_count}/100 match (sum_k I_k(2) = H(T))")

    # ----------------------------------------------------------------
    # PART 6: Specific examples with root details
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART 6: Detailed examples")
    print("=" * 80)

    n = 7
    for trial in range(5):
        T = random_tournament(n, seed=trial*17+3)
        polys, invs = compute_Ik_poly_n7(T)
        t3, t5, t7, bc = invs
        print(f"\n  Tournament #{trial}: t3={t3}, t5={t5}, t7={t7}, bc={bc}")
        for k in range(n):
            p = polys[k]
            info = analyze_roots(p)
            roots_str = ""
            if info['degree'] >= 1:
                if info['all_real']:
                    rr = info.get('real_roots', sorted([r.real for r in info['roots']]))
                    roots_str = ", ".join(f"{r:.4f}" for r in rr)
                else:
                    roots_str = ", ".join(f"{r:.4f}" for r in info['roots'])
                    roots_str += " [COMPLEX]"
            else:
                roots_str = "(constant)"

            disc_str = ""
            if 'discriminant' in info:
                disc_str = f", disc={info['discriminant']:.1f}"

            print(f"    k={k}: [{p[0]:.0f}, {p[1]:.0f}, {p[2]:.0f}] -> roots: {roots_str}{disc_str}")

    # ----------------------------------------------------------------
    # FINAL SUMMARY
    # ----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)

    print("""
  I_k(Omega(T), x) = A(n,k) + sum_I c_k^{(f_I, n-1)} * I(T) * x^{parts(I)}

  This polynomial in x packages the tournament invariants weighted by
  inflated Eulerian coefficients. At x=2, it gives a_k(T).

  Key findings from exhaustive/sampling analysis above.
""")

if __name__ == "__main__":
    main()
