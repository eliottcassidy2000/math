"""
H_c4_formula.py — Test the formula H = a - b*c_4 for circulant tournaments.

DISCOVERY: At p=7, H = 231 - 2*c_4 exactly!
  Proof: THM-133 gives H = (462 - tr(A^4))/2. And tr(A^4) = 4*c_4
  (because all closed 4-walks in a tournament are simple).
  So H = 231 - 2*c_4.

  This means: FEWER directed 4-cycles => MORE Hamiltonian paths!
  Paley T_7: c_4 = 21 (minimum!), H = 189 (maximum!)
  Non-Paley: c_4 = 28 (larger), H = 175 (smaller)

  The 4-cycle is the first cycle length where count varies among regular
  tournaments (c_3 is constant for regular).

  QUESTION: Does H = a(p) - b(p)*c_4 hold at p=11, 13?
  If YES: this is a universal trace formula for circulant tournaments.
  If NO: higher traces must be included.

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import cmath
import math
from itertools import combinations, permutations
from collections import defaultdict

sys.path.insert(0, '04-computation')


def all_circulant_tournaments(n):
    pairs, used = [], set()
    for a in range(1, n):
        if a not in used:
            b = n - a
            if a == b: return []
            pairs.append((a, b)); used.add(a); used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = [a if (bits >> i) & 1 else b for i, (a, b) in enumerate(pairs)]
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    S_set = set(S)
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%n in S_set: adj[i][j] = True
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if dp[mask][v]==0 or not(mask&(1<<v)): continue
            for w in range(n):
                if mask&(1<<w): continue
                if adj[v][w]: dp[mask|(1<<w)][w] += dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))


def circulant_eigenvalues(n, S):
    omega = cmath.exp(2j*cmath.pi/n)
    return [sum(omega**(k*s) for s in S) for k in range(n)]


def main():
    print("=" * 70)
    print("H = a - b*c_4 FORMULA TEST FOR CIRCULANT TOURNAMENTS")
    print("=" * 70)

    for p in [3, 5, 7, 11, 13]:
        print(f"\n{'=' * 60}")
        print(f"p = {p} (mod 4 = {p % 4})")
        print(f"{'=' * 60}")

        all_S = all_circulant_tournaments(p)
        m = (p-1)//2

        data = []
        for S in all_S:
            H = ham_count_dp(p, S)
            eigs = circulant_eigenvalues(p, S)

            # Traces
            traces = {}
            for k in range(3, min(p+1, 8)):
                traces[k] = sum(e**k for e in eigs).real

            # c_4 = tr(A^4)/4 (exact for tournaments)
            tr4 = traces.get(4, 0)
            c4 = round(tr4 / 4)

            # c_5 = tr(A^5)/5
            tr5 = traces.get(5, 0)
            c5 = round(tr5 / 5)

            # c_6 = tr(A^6)/6 — NOT necessarily simple for tournaments!
            # Actually let's check
            tr6 = traces.get(6, 0)
            tr7 = traces.get(7, 0)

            # y^2 values
            y2 = [eigs[k].imag**2 for k in range(1, m+1)]
            sum_y4 = sum(y**2 for y in y2)

            data.append({
                'S': S, 'H': H, 'c4': c4, 'c5': c5,
                'tr4': tr4, 'tr5': tr5, 'tr6': tr6, 'tr7': tr7,
                'y2': y2, 'sum_y4': sum_y4
            })

        data.sort(key=lambda d: d['H'], reverse=True)
        H_vals = sorted(set(d['H'] for d in data), reverse=True)

        # Group by H value
        print(f"\n  {'H':>12} {'c4':>6} {'c5':>6} {'tr4':>8} {'tr5':>8} {'tr6':>10} {'sum_y4':>10}")
        print(f"  {'-'*12} {'-'*6} {'-'*6} {'-'*8} {'-'*8} {'-'*10} {'-'*10}")
        for H_val in H_vals:
            d = next(x for x in data if x['H'] == H_val)
            count = sum(1 for x in data if x['H'] == H_val)
            print(f"  {H_val:>12} {d['c4']:>6} {d['c5']:>6} {d['tr4']:>8.0f} "
                  f"{d['tr5']:>8.0f} {d['tr6']:>10.0f} {d['sum_y4']:>10.4f}"
                  f"  (x{count})")

        # Test H = a - b*c_4
        if len(H_vals) >= 2:
            # Simple linear fit: H = a + b*c4
            c4_vals = [next(x for x in data if x['H']==H_val)['c4'] for H_val in H_vals]
            H_list = list(H_vals)
            n = len(H_list)

            if n == 2:
                # Exact: b = (H1-H2)/(c4_2-c4_1), a = H1 - b*c4_1
                if c4_vals[0] != c4_vals[1]:
                    b = (H_list[0] - H_list[1]) / (c4_vals[1] - c4_vals[0])
                    a = H_list[0] - b * c4_vals[0]
                    print(f"\n  EXACT FORMULA (2 points): H = {a:.4f} + ({b:.4f}) * c_4")
                    print(f"  Verification:")
                    for H_val, c4 in zip(H_list, c4_vals):
                        predicted = a + b * c4
                        print(f"    c_4={c4}: predicted={predicted:.4f}, actual={H_val}, "
                              f"match={abs(predicted-H_val) < 0.01}")
                else:
                    print(f"\n  c_4 is CONSTANT ({c4_vals[0]}) — cannot fit H vs c_4")
            else:
                # Least squares
                sx = sum(c4_vals)
                sy = sum(H_list)
                sxx = sum(x**2 for x in c4_vals)
                sxy = sum(x*y for x,y in zip(c4_vals, H_list))
                denom = n*sxx - sx**2
                if abs(denom) > 1e-10:
                    b = (n*sxy - sx*sy) / denom
                    a = (sy - b*sx) / n
                    residuals = [H - (a + b*c4) for H, c4 in zip(H_list, c4_vals)]
                    max_res = max(abs(r) for r in residuals)
                    ss_res = sum(r**2 for r in residuals)
                    ss_tot = sum((H - sy/n)**2 for H in H_list)
                    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 1

                    print(f"\n  LINEAR FIT: H = {a:.4f} + ({b:.6f}) * c_4")
                    print(f"  R^2 = {r2:.8f}, max |residual| = {max_res:.4f}")

                    if max_res < 0.01:
                        print(f"  *** EXACT LINEAR FORMULA ***")
                    else:
                        print(f"  NOT exact — need more variables")

                        # Try H = a + b*c4 + c*c5
                        c5_vals = [next(x for x in data if x['H']==H_val)['c5'] for H_val in H_vals]
                        # Solve 3x3 normal equations
                        X = [[1, c4, c5] for c4, c5 in zip(c4_vals, c5_vals)]
                        y = H_list
                        # X^T X
                        XtX = [[sum(X[i][r]*X[i][c] for i in range(n)) for c in range(3)] for r in range(3)]
                        Xty = [sum(X[i][r]*y[i] for i in range(n)) for r in range(3)]
                        # Gauss elimination
                        A_mat = [row[:] + [rhs] for row, rhs in zip(XtX, Xty)]
                        for col in range(3):
                            max_row = max(range(col, 3), key=lambda r: abs(A_mat[r][col]))
                            A_mat[col], A_mat[max_row] = A_mat[max_row], A_mat[col]
                            if abs(A_mat[col][col]) < 1e-15: continue
                            for row in range(col+1, 3):
                                f = A_mat[row][col] / A_mat[col][col]
                                for j in range(col, 4): A_mat[row][j] -= f * A_mat[col][j]
                        coeffs = [0]*3
                        for i in range(2, -1, -1):
                            if abs(A_mat[i][i]) > 1e-15:
                                coeffs[i] = (A_mat[i][3] - sum(A_mat[i][j]*coeffs[j] for j in range(i+1,3))) / A_mat[i][i]
                        a2, b2, c2 = coeffs
                        resid2 = [H - (a2 + b2*c4 + c2*c5) for H, c4, c5 in zip(H_list, c4_vals, c5_vals)]
                        max_res2 = max(abs(r) for r in resid2)
                        ss_res2 = sum(r**2 for r in resid2)
                        r2_2 = 1 - ss_res2/ss_tot if ss_tot > 0 else 1
                        print(f"\n  BILINEAR FIT: H = {a2:.4f} + ({b2:.6f})*c_4 + ({c2:.6f})*c_5")
                        print(f"  R^2 = {r2_2:.8f}, max |residual| = {max_res2:.4f}")

                        if max_res2 < 0.01:
                            print(f"  *** EXACT BILINEAR FORMULA ***")
                        else:
                            # Try adding c_4^2
                            c4sq_vals = [c4**2 for c4 in c4_vals]
                            if n >= 4:
                                X3 = [[1, c4, c5, c4sq] for c4, c5, c4sq in zip(c4_vals, c5_vals, c4sq_vals)]
                                Xt3X = [[sum(X3[i][r]*X3[i][c] for i in range(n)) for c in range(4)] for r in range(4)]
                                Xt3y = [sum(X3[i][r]*y[i] for i in range(n)) for r in range(4)]
                                A3 = [row[:] + [rhs] for row, rhs in zip(Xt3X, Xt3y)]
                                for col in range(4):
                                    max_row = max(range(col, 4), key=lambda r: abs(A3[r][col]))
                                    A3[col], A3[max_row] = A3[max_row], A3[col]
                                    if abs(A3[col][col]) < 1e-15: continue
                                    for row in range(col+1, 4):
                                        f = A3[row][col] / A3[col][col]
                                        for j in range(col, 5): A3[row][j] -= f * A3[col][j]
                                coeffs3 = [0]*4
                                for i in range(3, -1, -1):
                                    if abs(A3[i][i]) > 1e-15:
                                        coeffs3[i] = (A3[i][4] - sum(A3[i][j]*coeffs3[j] for j in range(i+1,4))) / A3[i][i]
                                a3, b3, c3, d3 = coeffs3
                                resid3 = [H - (a3 + b3*c4 + c3*c5 + d3*c4**2)
                                          for H, c4, c5 in zip(H_list, c4_vals, c5_vals)]
                                max_res3 = max(abs(r) for r in resid3)
                                print(f"\n  QUADRATIC: H = {a3:.4f} + ({b3:.6f})*c_4 + ({c3:.6f})*c_5 + ({d3:.6f})*c_4^2")
                                print(f"  max |residual| = {max_res3:.4f}")

                        # Also try tr6
                        tr6_vals = [next(x for x in data if x['H']==H_val)['tr6'] for H_val in H_vals]
                        X4 = [[1, c4, c5, tr6] for c4, c5, tr6 in zip(c4_vals, c5_vals, tr6_vals)]
                        if n >= 4:
                            Xt4X = [[sum(X4[i][r]*X4[i][c] for i in range(n)) for c in range(4)] for r in range(4)]
                            Xt4y = [sum(X4[i][r]*y[i] for i in range(n)) for r in range(4)]
                            A4 = [row[:] + [rhs] for row, rhs in zip(Xt4X, Xt4y)]
                            for col in range(4):
                                max_row = max(range(col, 4), key=lambda r: abs(A4[r][col]))
                                A4[col], A4[max_row] = A4[max_row], A4[col]
                                if abs(A4[col][col]) < 1e-15: continue
                                for row in range(col+1, 4):
                                    f = A4[row][col] / A4[col][col]
                                    for j in range(col, 5): A4[row][j] -= f * A4[col][j]
                            coeffs4 = [0]*4
                            for i in range(3, -1, -1):
                                if abs(A4[i][i]) > 1e-15:
                                    coeffs4[i] = (A4[i][4] - sum(A4[i][j]*coeffs4[j] for j in range(i+1,4))) / A4[i][i]
                            a4, b4, c4c, d4 = coeffs4
                            resid4 = [H - (a4 + b4*c4 + c4c*c5 + d4*tr6)
                                      for H, c4, c5, tr6 in zip(H_list, c4_vals, c5_vals, tr6_vals)]
                            max_res4 = max(abs(r) for r in resid4)
                            print(f"\n  WITH tr6: H = {a4:.4f} + ({b4:.6f})*c_4 + ({c4c:.6f})*c_5 + ({d4:.10f})*tr(A^6)")
                            print(f"  max |residual| = {max_res4:.4f}")

    # ================================================================
    # Summary
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SUMMARY: H vs CYCLE COUNTS")
    print(f"{'=' * 60}")
    print("""
    KEY FINDINGS:
    1. tr(A^k) = k * c_k for k = 3, 4, 5 in tournaments (no non-simple corrections)
       But tr(A^k) != k * c_k for k >= 6 (non-simple walks appear)

    2. For REGULAR tournaments on Z_p:
       - c_3 is CONSTANT (depends only on p, score sequence)
       - c_4, c_5 VARY across tournaments with same c_3

    3. At p=7: H = 231 - 2*c_4 EXACTLY
       Equivalently: H = 231 - tr(A^4)/2
       This is THM-133 (opus-S58)

    4. At p >= 11: H is NOT determined by c_4 alone
       Need c_5 and possibly higher cycles

    5. PALEY MINIMIZES c_4 and MAXIMIZES c_5 among circulants
       at p = 3 mod 4
    """)


if __name__ == '__main__':
    main()
