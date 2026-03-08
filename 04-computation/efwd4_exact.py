#!/usr/bin/env python3
"""
efwd4_exact.py — Exact rational formula for E[fwd^4].

DISCOVERY: E[fwd^4] = A4(n) + B4(n)*t3 + C4(n)*t5 + D4(n)*alpha_2

Find exact rational coefficients.

Author: opus-2026-03-07-S46c
"""
from itertools import permutations, combinations
from math import comb, factorial
from fractions import Fraction

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def count_3cycles_list(adj, n):
    cycles = []
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            cycles.append((i, j, k))
    return cycles

def count_5cycles(adj, n):
    count = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                count += 1
    return count // 5

print("=" * 60)
print("EXACT E[fwd^4] FORMULA")
print("=" * 60)

for n in [5, 6]:
    print(f"\nn={n}:")
    m_vals = n*(n-1)//2
    seen = set()
    data = []

    for bits in range(1 << m_vals):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        cycles3 = count_3cycles_list(adj, n)
        t3 = len(cycles3)
        t5 = count_5cycles(adj, n)
        alpha_2 = sum(1 for a in range(len(cycles3))
                      for b in range(a+1, len(cycles3))
                      if set(cycles3[a]).isdisjoint(set(cycles3[b])))

        total = factorial(n)
        m4 = Fraction(sum(k**4 * F[k] for k in range(n)), total)

        data.append({'t3': t3, 't5': t5, 'alpha_2': alpha_2, 'E4': m4})

    print(f"  {len(data)} F-classes")

    # Solve exact system: E4 = a + b*t3 + c*t5 + d*alpha_2
    # Need 4 data points with linearly independent (1, t3, t5, alpha_2)
    # Use Fraction arithmetic

    # Find 4 data points with different (t3, t5, alpha_2)
    from itertools import combinations as combo_iter
    for indices in combo_iter(range(len(data)), 4):
        rows = [data[i] for i in indices]
        # Matrix: [1, t3, t5, alpha_2]
        mat = [[Fraction(1), Fraction(r['t3']), Fraction(r['t5']), Fraction(r['alpha_2'])]
               for r in rows]
        rhs = [r['E4'] for r in rows]

        # Check if matrix is non-singular (using determinant)
        def det4(M):
            # Expand 4x4 determinant
            from itertools import permutations as perm_iter
            d = Fraction(0)
            for p in perm_iter(range(4)):
                sgn = 1
                inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
                if inv % 2: sgn = -1
                term = Fraction(sgn)
                for i in range(4):
                    term *= M[i][p[i]]
                d += term
            return d

        d = det4(mat)
        if d == 0:
            continue

        # Solve by Cramer's rule
        coeffs = []
        for col in range(4):
            mod_mat = [row[:] for row in mat]
            for row_idx in range(4):
                mod_mat[row_idx][col] = rhs[row_idx]
            coeffs.append(det4(mod_mat) / d)

        # Verify against ALL data
        ok = True
        for d_item in data:
            pred = coeffs[0] + coeffs[1]*d_item['t3'] + coeffs[2]*d_item['t5'] + coeffs[3]*d_item['alpha_2']
            if pred != d_item['E4']:
                ok = False
                break

        if ok:
            print(f"  EXACT FORMULA FOUND:")
            print(f"  E[fwd^4] = {coeffs[0]} + {coeffs[1]}*t3 + {coeffs[2]}*t5 + {coeffs[3]}*alpha_2")
            print(f"           = {float(coeffs[0]):.6f} + {float(coeffs[1]):.6f}*t3 + {float(coeffs[2]):.6f}*t5 + {float(coeffs[3]):.6f}*alpha_2")
            break
    else:
        print("  No exact formula found with (t3, t5, alpha_2)")

    # Also try just (t3, t5) for n=5
    if n == 5:
        for i1, i2 in combo_iter(range(len(data)), 2):
            r1, r2 = data[i1], data[i2]
            mat2 = [[Fraction(1), Fraction(r1['t3']), Fraction(r1['t5'])],
                     [Fraction(1), Fraction(r2['t3']), Fraction(r2['t5'])]]
            # Need 3 equations for 3 unknowns. Use 3 data points.
            pass

        # Actually solve (1, t3, t5) system
        for indices in combo_iter(range(len(data)), 3):
            rows = [data[i] for i in indices]
            mat3 = [[Fraction(1), Fraction(r['t3']), Fraction(r['t5'])] for r in rows]
            rhs3 = [r['E4'] for r in rows]

            # 3x3 determinant
            def det3(M):
                return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                       -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                       +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

            d3 = det3(mat3)
            if d3 == 0:
                continue

            c3 = []
            for col in range(3):
                m3 = [row[:] for row in mat3]
                for r in range(3):
                    m3[r][col] = rhs3[r]
                c3.append(det3(m3) / d3)

            ok3 = all(c3[0] + c3[1]*d_item['t3'] + c3[2]*d_item['t5'] == d_item['E4'] for d_item in data)
            if ok3:
                print(f"\n  At n=5 (t3, t5 sufficient):")
                print(f"  E[fwd^4] = {c3[0]} + {c3[1]}*t3 + {c3[2]}*t5")
                print(f"           = {float(c3[0]):.6f} + {float(c3[1]):.6f}*t3 + {float(c3[2]):.6f}*t5")
                break
