"""
var_ratio_formula.py -- kind-pasteur-2026-03-14-S105d
EXACT computation of Var(H)/Mean(H)^2 at n=3,4,5,6
to find the correction formula from 1/3.

KNOWN:
  n=3: 1/3 exactly
  n=4: 1/3 exactly
  n=5: 19/60 = 1/3 - 1/60 = 1/3 - 1/(5*4*3)

QUESTION: What is the EXACT value at n=6?
  If the pattern is 1/3 - c_n/n(n-1)(n-2), what is c_6?

Also: compute E_4/E_0 at n=5 exactly to understand the
level-4 energy contribution.
"""

import sys, math
import numpy as np
from fractions import Fraction
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

def count_ham_paths(adj, n):
    """Count Hamiltonian paths using Held-Karp DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def bits_to_tournament(bits, n):
    """Convert integer bits to adjacency matrix."""
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def all_H_values(n):
    """Compute H for all tournaments on n vertices."""
    m = C(n, 2)
    values = []
    for bits in range(1 << m):
        adj = bits_to_tournament(bits, n)
        h = count_ham_paths(adj, n)
        values.append(h)
    return values

def main():
    print("=" * 70)
    print("EXACT VAR/MEAN^2 COMPUTATION")
    print("kind-pasteur-2026-03-14-S105d")
    print("=" * 70)

    # n=3,4,5: exhaustive
    for n in [3, 4, 5]:
        vals = all_H_values(n)
        N = len(vals)
        m = C(n, 2)

        # Exact computation using Fraction
        sum_h = sum(vals)
        sum_h2 = sum(h*h for h in vals)

        mean = Fraction(sum_h, N)
        mean_h2 = Fraction(sum_h2, N)
        var = mean_h2 - mean * mean
        ratio = var / (mean * mean) if mean != 0 else Fraction(0)

        # Level-2 energy formula
        level2_ratio = Fraction(2*(n-2), n*(n-1))

        print(f"\n  n={n} (m={m}, N=2^m={N}):")
        print(f"    sum(H) = {sum_h}")
        print(f"    sum(H^2) = {sum_h2}")
        print(f"    Mean = {mean} = {float(mean):.6f}")
        print(f"    Mean(H^2) = {mean_h2} = {float(mean_h2):.6f}")
        print(f"    Var = {var} = {float(var):.6f}")
        print(f"    Var/Mean^2 = {ratio} = {float(ratio):.10f}")
        print(f"    Level-2 formula = {level2_ratio} = {float(level2_ratio):.10f}")
        diff = ratio - Fraction(1, 3)
        print(f"    Deviation from 1/3 = {diff} = {float(diff):.10f}")

        if diff != 0:
            nnn = n*(n-1)*(n-2)
            c_n = -diff * nnn
            print(f"    c_n = -(dev)*n(n-1)(n-2) = {c_n} = {float(c_n):.6f}")

    # n=6: exhaustive (2^15 = 32768 tournaments)
    print(f"\n\n  Computing n=6 (m=15, N=32768)...")
    n = 6
    vals6 = all_H_values(n)
    N = len(vals6)

    sum_h = sum(vals6)
    sum_h2 = sum(h*h for h in vals6)

    mean = Fraction(sum_h, N)
    mean_h2 = Fraction(sum_h2, N)
    var = mean_h2 - mean * mean
    ratio = var / (mean * mean)

    level2_ratio = Fraction(2*(n-2), n*(n-1))

    print(f"\n  n={n} (m=15, N=32768):")
    print(f"    sum(H) = {sum_h}")
    print(f"    sum(H^2) = {sum_h2}")
    print(f"    Mean = {mean} = {float(mean):.6f}")
    print(f"    Mean(H^2) = {mean_h2} = {float(mean_h2):.6f}")
    print(f"    Var = {var} = {float(var):.6f}")
    print(f"    Var/Mean^2 = {ratio} = {float(ratio):.10f}")
    print(f"    Level-2 formula = {level2_ratio} = {float(level2_ratio):.10f}")
    diff = ratio - Fraction(1, 3)
    print(f"    Deviation from 1/3 = {diff} = {float(diff):.10f}")

    nnn = n*(n-1)*(n-2)
    c_n = -diff * nnn
    print(f"    c_n = -(dev)*n(n-1)(n-2) = {c_n} = {float(c_n):.6f}")

    # Now compare
    print(f"\n  SUMMARY TABLE:")
    print(f"  {'n':>3} {'Var/Mean^2':>15} {'1/3-ratio':>15} {'c_n':>10}")
    print(f"  {'-'*50}")

    results = {}
    for n_val in [3, 4, 5, 6]:
        if n_val == 6:
            vals = vals6
        else:
            vals = all_H_values(n_val)
        N = len(vals)
        sum_h = sum(vals)
        sum_h2 = sum(h*h for h in vals)
        mean = Fraction(sum_h, N)
        var = Fraction(sum_h2, N) - mean * mean
        ratio = var / (mean * mean)
        diff = Fraction(1, 3) - ratio
        nnn = n_val*(n_val-1)*(n_val-2)
        c_n = diff * nnn
        results[n_val] = (ratio, diff, c_n)
        print(f"  {n_val:3d} {str(ratio):>15} {str(diff):>15} {str(c_n):>10}")

    # The key question: what is the pattern of c_n?
    print(f"\n  The c_n sequence: ", end="")
    for n_val in [3, 4, 5, 6]:
        print(f"c_{n_val}={results[n_val][2]}", end="  ")
    print()

    # Let's also compute the LEVEL-4 energy at n=5 and n=6
    print(f"\n{'='*70}")
    print("LEVEL-4 ENERGY DECOMPOSITION")
    print(f"{'='*70}")

    for n_val in [5, 6]:
        if n_val == 6:
            vals = vals6
        else:
            vals = all_H_values(n_val)
        N = len(vals)
        sum_h = sum(vals)
        sum_h2 = sum(h*h for h in vals)
        mean = Fraction(sum_h, N)
        var = Fraction(sum_h2, N) - mean * mean
        E0 = mean * mean
        E2 = Fraction(2*(n_val-2), n_val*(n_val-1)) * E0
        E_rest = var - E2

        print(f"\n  n={n_val}:")
        print(f"    E_0 = {E0} = {float(E0):.6f}")
        print(f"    E_2 = {E2} = {float(E2):.6f}")
        print(f"    E_4+ = {E_rest} = {float(E_rest):.6f}")
        print(f"    E_2/E_0 = {E2/E0} = {float(E2/E0):.10f}")
        print(f"    E_4+/E_0 = {E_rest/E0} = {float(E_rest/E0):.10f}")
        print(f"    Var/Mean^2 = {float(var/E0):.10f}")

    # At n=6, deg(H) = 2*floor(5/2) = 4, so still only levels 0,2,4
    # At n=7, deg(H) = 2*floor(6/2) = 6, so levels 0,2,4,6

    # Now let's look for a formula for E_4
    print(f"\n{'='*70}")
    print("LOOKING FOR THE E_4 FORMULA")
    print(f"{'='*70}")

    # At n=5: E_4 = 15/16, E_0 = 900/16 = 225/4, ratio = (15/16)/(225/4) = 15/(16*225/4) = 15*4/(16*225) = 60/3600 = 1/60
    # At n=6: what is E_4/E_0?

    for n_val in [5, 6]:
        if n_val == 6:
            vals = vals6
        else:
            vals = all_H_values(n_val)
        N = len(vals)
        mean = Fraction(sum(vals), N)
        var = Fraction(sum(h*h for h in vals), N) - mean * mean
        E0 = mean * mean
        E2_ratio = Fraction(2*(n_val-2), n_val*(n_val-1))
        E4_ratio = var/E0 - E2_ratio

        print(f"\n  n={n_val}:")
        print(f"    E_4/E_0 = {E4_ratio} = {float(E4_ratio):.10f}")

        # Try to find pattern
        # n=5: E_4/E_0 = 1/60
        # n=6: E_4/E_0 = ???
        # 60 = 5*4*3 = 5!/(5-3)!
        if n_val == 5:
            print(f"    = 1/60 = 1/(5*4*3) = 1/P(5,3)")
        if n_val == 6:
            denom = 1 / float(E4_ratio) if E4_ratio > 0 else 0
            print(f"    1/(E_4/E_0) = {denom:.4f}")
            print(f"    6*5*4 = {6*5*4}")
            print(f"    6*5*4*3 = {6*5*4*3}")

    # Distribution at n=6
    print(f"\n{'='*70}")
    print("H-DISTRIBUTION AT n=6")
    print(f"{'='*70}")
    counts6 = Counter(vals6)
    print(f"\n  H values at n=6 ({len(counts6)} distinct):")
    for h in sorted(counts6.keys()):
        count = counts6[h]
        bar = "#" * (count // 100)
        print(f"    H={h:4d}: {count:5d} {bar}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
