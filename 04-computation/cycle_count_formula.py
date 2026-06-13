"""
opus-2026-05-23-S5: Formula for odd cycle counts in all-0 staircase T_k.

Data from hyp1732_full_investigation.py:
  k=2: #3=2, #5=0, #7=0
  k=3: #3=6, #5=6, #7=0
  k=4: #3=12, #5=28, #7=28
  k=5: #3=20, #5=80, #7=220
  k=6: #3=30, #5=180, #7=906

Discovered: #5-cycles = C(k,3) * (k+3) = k(k-1)(k-2)(k+3)/6

Goals:
  1. Verify #5-cycle formula for k=7,8
  2. Find #7-cycle formula
  3. Find #9-cycle formula if k=7 is feasible
  4. Compute full IP of T_4 (all 68 cycles, find alpha_2)
"""
import sys
import time
from math import comb
sys.path.insert(0, '04-computation')

from hyp1732_full_investigation import build_Tk, find_odd_cycles, build_conflict_adj

def count_cycles_by_length(k, max_len=None):
    """Count odd cycles of each length in T_k."""
    n = 2*k
    A = build_Tk(k)
    if max_len is None:
        max_len = n

    cycles = find_odd_cycles(A, n, max_length=max_len)
    counts = {}
    for c in cycles:
        L = len(c)
        counts[L] = counts.get(L, 0) + 1
    return counts, cycles

def predict_5cycle(k):
    """Predicted #5-cycles = C(k,3) * (k+3)."""
    return comb(k, 3) * (k + 3)

def fit_polynomial(data, degree):
    """Fit polynomial of given degree to data [(k, val)]."""
    from fractions import Fraction
    n = len(data)
    # Basis: 1, k, k^2, ..., k^degree
    A = [[Fraction(k**j) for j in range(degree+1)] for k,_ in data]
    b = [Fraction(v) for _,v in data]

    # Gaussian elimination
    for col in range(min(n, degree+1)):
        pivot = None
        for row in range(col, n):
            if A[row][col] != 0:
                pivot = row; break
        if pivot is None: continue
        A[col], A[pivot] = A[pivot], A[col]
        b[col], b[pivot] = b[pivot], b[col]
        for row in range(n):
            if row != col and A[row][col] != 0:
                f = A[row][col] / A[col][col]
                for c2 in range(degree+1):
                    A[row][c2] -= f * A[col][c2]
                b[row] -= f * b[col]

    coeffs = [b[i]/A[i][i] if A[i][i] != 0 else Fraction(0)
              for i in range(min(n, degree+1))]
    return coeffs

def eval_poly(coeffs, k):
    """Evaluate polynomial sum_i coeffs[i] * k^i at k."""
    from fractions import Fraction
    return sum(c * Fraction(k)**i for i, c in enumerate(coeffs))

if __name__ == '__main__':
    print("=== Cycle Count Formulas for All-0 Staircase T_k ===\n")

    # Known data
    known_data = {
        2: {3: 2, 5: 0, 7: 0},
        3: {3: 6, 5: 6, 7: 0},
        4: {3: 12, 5: 28, 7: 28},
        5: {3: 20, 5: 80, 7: 220},
        6: {3: 30, 5: 180, 7: 906},
    }

    # Verify 5-cycle formula: C(k,3)*(k+3)
    print("5-CYCLE COUNT FORMULA: C(k,3)*(k+3)")
    for k, counts in known_data.items():
        pred = predict_5cycle(k)
        actual = counts.get(5, 0)
        ok = "✓" if pred == actual else f"✗ pred={pred}"
        print(f"  k={k}: predicted={pred}, actual={actual} {ok}")

    print()

    # Compute k=7 to verify
    print("Computing k=7 odd cycles (max_len=9 to limit 9-cycles)...")
    t0 = time.time()
    counts7, cycles7 = count_cycles_by_length(7, max_len=9)
    print(f"  k=7: {counts7} [{time.time()-t0:.1f}s]")
    print(f"  Predicted #5-cycles: {predict_5cycle(7)}")
    print()

    # Update known data with k=7
    known_data[7] = counts7

    # Analyze 5-cycle formula
    print("VERIFYING #5-cycles = C(k,3)*(k+3):")
    for k in range(2, 8):
        pred = predict_5cycle(k)
        actual = known_data.get(k, {}).get(5, 0)
        ok = "✓" if pred == actual else f"✗ actual={actual}"
        print(f"  k={k}: C({k},3)*({k}+3) = {comb(k,3)}*{k+3} = {pred} {ok}")

    print()

    # Fit 7-cycle count polynomial
    print("ANALYZING #7-cycles:")
    data_7 = [(k, known_data[k].get(7,0)) for k in sorted(known_data) if k >= 4]
    print("  Data:", data_7)

    # Try formula: C(k,4) * f(k)
    print("  C(k,4) = ", [(k, comb(k,4)) for k,_ in data_7])
    print("  #7/C(k,4) = ", [(k, v/comb(k,4) if comb(k,4)>0 else 'N/A') for k,v in data_7])

    # Also try: check if #7 = sum_j a_j * C(k,j)
    print("  Fitting polynomial in k:")
    for degree in range(3, 7):
        if len(data_7) >= degree + 1:
            coeffs = fit_polynomial(data_7[:degree+1], degree)
            pred_vals = [(k, int(eval_poly(coeffs, k))) for k,_ in data_7]
            ok = all(p == v for (k,p), (_, v) in zip(pred_vals, data_7))
            if ok:
                print(f"  Degree {degree} fit: {coeffs}")
                # Try to express as binomial combination
                # Convert to C(k,j) basis using finite differences
                print(f"  Predictions: {pred_vals}")
                break

    print()

    # Compute full IP of T_4 (alpha_2 from all 68 cycles)
    print("=== FULL IP OF T_4 ===")
    n4 = 8
    A4 = build_Tk(4)
    t0 = time.time()
    all_odd4, _ = count_cycles_by_length(4)
    print(f"T_4: n=8, cycle counts: {all_odd4}")

    # Actually get the cycles list
    _, odd_list4 = count_cycles_by_length(4)
    m4 = len(odd_list4)
    print(f"Total odd cycles: {m4}")

    adj4 = build_conflict_adj(odd_list4)

    # Count pairs (alpha_2)
    pairs4 = [(i,j) for i in range(m4) for j in range(i+1,m4)
              if not ((adj4[i]>>j)&1)]
    alpha2_full = len(pairs4)
    print(f"alpha_2 (full IP) = {alpha2_full}")

    # Categorize pairs by cycle length combination
    pair_cats = {}
    for i,j in pairs4:
        L_i = len(odd_list4[i])
        L_j = len(odd_list4[j])
        key = (min(L_i,L_j), max(L_i,L_j))
        pair_cats[key] = pair_cats.get(key, 0) + 1
    print(f"Pair breakdown by (len1, len2): {sorted(pair_cats.items())}")
    print(f"Full I(Omega(T_4), x) = 1 + {m4}x + {alpha2_full}x^2")

    # Check real-rootedness
    disc = m4**2 - 4*alpha2_full
    print(f"Discriminant = {m4}^2 - 4*{alpha2_full} = {disc}")
    print(f"Real-rooted: {disc >= 0}")
    if disc >= 0:
        import math
        r1 = (-m4 + math.sqrt(disc)) / 2
        r2 = (-m4 - math.sqrt(disc)) / 2
        print(f"Roots: {r1:.4f}, {r2:.4f}")
        print(f"Roots x where p=|p| satisfies -1/p in [r2,r1]: p in [{-1/r1:.4f}, {-1/r2:.4f}]")

    print()
    print("=== SUMMARY TABLE ===")
    print(f"{'k':>3} {'#3':>6} {'#5':>8} {'#7':>8} {'#9':>8} {'#total':>8} {'alpha_2':>8}")
    for k in sorted(known_data):
        counts = known_data[k]
        n3 = counts.get(3,0)
        n5 = counts.get(5,0)
        n7 = counts.get(7,0)
        n9 = counts.get(9,0)
        total = sum(counts.values())
        print(f"  k={k}: #3={n3:4d}, #5={n5:4d}, #7={n7:4d}, #9={n9:4d}, total={total:5d}")

    print()
    print("=== CREATIVE STRUCTURE ===")
    print(f"#3-cycles = k(k-1) = 2*C(k,2)")
    print(f"#5-cycles = k(k-1)(k-2)(k+3)/6 = C(k,3)*(k+3)")
    print()
    print("Interesting ratio: #5/#3 = C(k,3)*(k+3) / (k*(k-1))")
    for k in range(3, 8):
        n3 = k*(k-1)
        n5 = comb(k,3)*(k+3)
        print(f"  k={k}: #5/#3 = {n5}/{n3} = {n5/n3:.4f}")
