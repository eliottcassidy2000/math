#!/usr/bin/env python3
"""
fourier_cycle_bridge_n7.py -- kind-pasteur-2026-03-13-S61

Test the Fourier-Cycle Bridge (THM-164) at n=7:
  - Degree 6 should encode 7-cycle counts
  - Degree 4 should encode 5-cycle counts (beyond score)
  - At n=7, ALL regular tournaments have same score => H_2=0 residual
    So H = H_0 + H_4 + H_6 for regular tournaments

Key questions:
1. Does degree-6 Fourier capture 7-cycle information?
2. What is the energy distribution (degree 0/2/4/6)?
3. For the three rigid regular classes at n=7 (H=171/175/189):
   How do H_4 and H_6 split?

Since full Fourier at m=21 is C(21,6)=54264 coefficients for degree-6,
we use SAMPLING: compute H_4 and H_6 as residuals, then correlate with
cycle counts directly.

Author: kind-pasteur-2026-03-13-S61
"""

import math
import random
from itertools import combinations, permutations
from collections import defaultdict

random.seed(42)


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    if n <= 1:
        return 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_ham_cycles(A, n):
    """Count directed Hamiltonian cycles starting from vertex 0."""
    if n <= 2:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    total = 0
    for v in range(1, n):
        if (full, v) in dp and A[v][0]:
            total += dp[(full, v)]
    return total


def count_directed_cycles_on_subset(A, n, subset):
    """Count directed Hamiltonian cycles on a vertex subset."""
    verts = list(subset)
    k = len(verts)
    if k < 3:
        return 0
    sub_A = [[A[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
    return count_ham_cycles(sub_A, k)


def compute_H2(A, n):
    c2 = math.factorial(n - 2) / (2 ** (n - 2))
    scores = [sum(A[v]) for v in range(n)]
    half = (n - 1) / 2
    Z = [-2 * (s - half)**2 + half for s in scores]
    return c2 * sum(Z)


def score_variance(A, n):
    scores = [sum(A[v]) for v in range(n)]
    half = (n - 1) / 2
    return sum((s - half)**2 for s in scores) / n


# ========================================================================
# ANALYSIS 1: ENERGY DISTRIBUTION AT n=5,6 (EXHAUSTIVE)
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: FOURIER ENERGY DISTRIBUTION (EXHAUSTIVE n=5,6)")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m
    EH = math.factorial(n) / (2 ** (n - 1))

    H_sum = 0
    H2_sum = 0
    Hsq_sum = 0
    H2sq_sum = 0

    data = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        H2 = compute_H2(A, n)
        H4 = H - EH - H2

        data.append((H, H2, H4))
        H_sum += H
        Hsq_sum += H * H
        H2_sum += H2
        H2sq_sum += H2 * H2

    avg_H = H_sum / total
    var_H = Hsq_sum / total - avg_H**2
    var_H2 = H2sq_sum / total
    var_H4 = sum(h4**2 for _, _, h4 in data) / total

    # Cross term: should be zero by orthogonality
    cross = sum(h2 * h4 for _, h2, h4 in data) / total

    print(f"\nn={n} (m={m}, total={total}):")
    print(f"  E[H] = {avg_H:.4f}")
    print(f"  Var(H) = {var_H:.4f}")
    print(f"  Var(H_2) = {var_H2:.4f}  ({100*var_H2/var_H:.2f}% of Var(H))")
    print(f"  Var(H_4) = {var_H4:.4f}  ({100*var_H4/var_H:.2f}% of Var(H))")
    print(f"  Cross(H2,H4) = {cross:.6f} (should be ~0)")

    residual_var = var_H - var_H2 - var_H4
    if abs(residual_var) > 0.01:
        print(f"  Var(H_6+) = {residual_var:.4f}  ({100*residual_var/var_H:.2f}% of Var(H))")
        print(f"  => DEGREE-6 TERMS ARE PRESENT!")
    else:
        print(f"  Residual variance = {residual_var:.6f} (negligible)")
        print(f"  => H = H_0 + H_2 + H_4 exactly at n={n}")


# ========================================================================
# ANALYSIS 2: DEGREE-6 AT n=7 VIA REGULAR TOURNAMENT SAMPLING
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: n=7 REGULAR TOURNAMENTS — CYCLE STRUCTURE")
print("=" * 70)

n = 7
m = n * (n - 1) // 2  # 21
total = 1 << m  # 2097152

print(f"\nn={n}: m={m}, total tournaments = {total}")
print(f"Scanning all 2^{m} tournaments for regular ones (all scores = 3)...")

EH = math.factorial(n) / (2 ** (n - 1))
print(f"E[H] = {EH:.4f}")

# At n=7, regular tournaments have scores = (3,3,3,3,3,3,3)
# There are 2640 regular tournaments on 7 vertices

regular_data = []  # (bits, H, c3_dir, c5_dir, c7_dir, AA^T_var)
count = 0

for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = [sum(A[v]) for v in range(n)]
    if any(s != 3 for s in scores):
        continue

    H = count_ham_paths(A, n)

    # Count directed 3-cycles
    c3_dir = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            c3_dir += 1
        if A[a][c] and A[c][b] and A[b][a]:
            c3_dir += 1

    # Count directed 5-cycles on each 5-vertex subset
    c5_dir = 0
    for subset in combinations(range(n), 5):
        c5_dir += count_directed_cycles_on_subset(A, n, subset)

    # Count directed 7-cycles (Hamiltonian cycles)
    c7_dir = count_ham_cycles(A, n)

    # AA^T off-diagonal variance
    off_diag = []
    for i in range(n):
        for j in range(i+1, n):
            off_diag.append(sum(A[i][k] * A[j][k] for k in range(n)))
    aat_var = sum((x - sum(off_diag)/len(off_diag))**2 for x in off_diag) / len(off_diag)

    regular_data.append({
        'bits': bits, 'H': H,
        'c3_dir': c3_dir, 'c5_dir': c5_dir, 'c7_dir': c7_dir,
        'aat_var': round(aat_var, 6)
    })
    count += 1

print(f"Found {count} regular tournaments on 7 vertices")

# Group by H value
by_H = defaultdict(list)
for d in regular_data:
    by_H[d['H']].append(d)

print(f"\nRegular tournament classes at n=7:")
print(f"  {'H':>5s} | {'count':>6s} | {'c3_dir':>7s} | {'c5_dir':>7s} | {'c7_dir':>7s} | {'AA^T var':>10s}")
print(f"  {'':->5s}-+-{'':->6s}-+-{'':->7s}-+-{'':->7s}-+-{'':->7s}-+-{'':->10s}")
for H in sorted(by_H.keys()):
    group = by_H[H]
    c3s = set(d['c3_dir'] for d in group)
    c5s = set(d['c5_dir'] for d in group)
    c7s = set(d['c7_dir'] for d in group)
    avs = set(d['aat_var'] for d in group)
    print(f"  {H:>5d} | {len(group):>6d} | {str(sorted(c3s)):>7s} | {str(sorted(c5s)):>7s} | "
          f"{str(sorted(c7s)):>7s} | {str(sorted(avs))}")


# ========================================================================
# ANALYSIS 3: FOURIER DECOMPOSITION OF REGULAR n=7
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: FOURIER DECOMPOSITION FOR REGULAR n=7")
print("=" * 70)

# For regular tournaments, Var_s = 0, so H_2 = c2 * sum(Z) where all Z_v = (n-1)/2 = 3
# H_2 = c2 * n * (n-1)/2 = (5!/2^5) * 7 * 3
c2 = math.factorial(n - 2) / (2 ** (n - 2))
H2_regular = c2 * n * (n - 1) / 2
print(f"\nFor regular n=7: c2 = {c2:.6f}")
print(f"H_2 = c2 * 7 * 3 = {H2_regular:.4f}")
print(f"E[H] = {EH:.4f}")
print(f"H_0 + H_2 = {EH + H2_regular:.4f}")

# So for regular tournaments:
# H_4 + H_6 = H - H_0 - H_2
for H in sorted(by_H.keys()):
    H46 = H - EH - H2_regular
    group = by_H[H]
    c5 = list(set(d['c5_dir'] for d in group))[0]
    c7 = list(set(d['c7_dir'] for d in group))[0]
    print(f"\n  H={H}: H_4 + H_6 = {H46:.4f}")
    print(f"    c5_dir = {c5}, c7_dir = {c7}")
    print(f"    2*c5_dir/5 = {2*c5/5:.4f}")
    print(f"    2*c7_dir/7 = {2*c7/7:.4f}")


# ========================================================================
# ANALYSIS 4: SEPARATE H_4 AND H_6 AT n=6 (EXACT)
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: EXACT SEPARATION OF H_4, H_6 AT n=6")
print("=" * 70)

n = 6
m = n * (n - 1) // 2  # 15
total = 1 << m  # 32768
EH6 = math.factorial(n) / (2 ** (n - 1))

# We need to compute H_6 explicitly via Fourier
# H_6 = sum_{|S|=6} hat(H)(S) * prod_{i in S} sigma_i
# First compute all degree-6 coefficients

print(f"\nn={n}: Computing ALL degree-6 Fourier coefficients...")
print(f"  C({m},6) = {math.comb(m,6)} subsets to check")

edge_pairs = []
for i in range(n):
    for j in range(i+1, n):
        edge_pairs.append((i, j))

# Precompute all H values
H_vals = {}
sigma_list = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    sigma = []
    for i in range(n):
        for j in range(i+1, n):
            sigma.append(1 if A[i][j] else -1)
    sigma = tuple(sigma)
    H = count_ham_paths(A, n)
    H_vals[sigma] = H
    sigma_list.append(sigma)

# Compute degree-6 coefficients
deg6_coeffs = {}
deg6_count = 0
for S in combinations(range(m), 6):
    coeff = 0
    for sigma in sigma_list:
        chi = 1
        for idx in S:
            chi *= sigma[idx]
        coeff += H_vals[sigma] * chi
    coeff /= total
    if abs(coeff) > 1e-10:
        deg6_coeffs[S] = coeff
        deg6_count += 1

print(f"  Non-zero degree-6 coefficients: {deg6_count}")

if deg6_count > 0:
    mags = sorted(set(round(abs(c), 8) for c in deg6_coeffs.values()))
    print(f"  Distinct magnitudes: {mags[:10]}{'...' if len(mags) > 10 else ''}")

    # Analyze vertex support
    support_dist = defaultdict(int)
    for S, c in deg6_coeffs.items():
        edges = [edge_pairs[idx] for idx in S]
        verts = set()
        for e in edges:
            verts.update(e)
        support_dist[len(verts)] += 1

    print(f"  Support size distribution:")
    for sz in sorted(support_dist):
        print(f"    |support|={sz}: {support_dist[sz]} terms")

    # Compute H_6 for each tournament
    print(f"\n  H_6 values by score class:")
    score_to_H6 = defaultdict(set)
    score_to_H4 = defaultdict(set)
    for sigma in sigma_list:
        H = H_vals[sigma]
        H2 = 0
        # Need to compute H_2 properly
        bits = 0
        for idx in range(m):
            if sigma[idx] == 1:
                bits |= (1 << idx)
        A = binary_to_tournament(bits, n)
        H2 = compute_H2(A, n)

        H6_val = 0
        for S, c in deg6_coeffs.items():
            chi = 1
            for idx in S:
                chi *= sigma[idx]
            H6_val += c * chi
        H6_val = round(H6_val, 4)

        H4_val = round(H - EH6 - H2 - H6_val, 4)

        scores = tuple(sorted([sum(A[v]) for v in range(n)]))
        score_to_H6[scores].add(H6_val)
        score_to_H4[scores].add(H4_val)

    print(f"  {'Score class':>25s} | {'#H4':>6s} | {'#H6':>6s} | {'H4 range':>20s} | {'H6 range':>20s}")
    print(f"  {'':->25s}-+-{'':->6s}-+-{'':->6s}-+-{'':->20s}-+-{'':->20s}")
    for sc in sorted(score_to_H6.keys()):
        h4s = sorted(score_to_H4[sc])
        h6s = sorted(score_to_H6[sc])
        h4_str = f"{h4s[0]:.1f}" if len(h4s) == 1 else f"{h4s[0]:.1f} to {h4s[-1]:.1f}"
        h6_str = f"{h6s[0]:.1f}" if len(h6s) == 1 else f"{h6s[0]:.1f} to {h6s[-1]:.1f}"
        print(f"  {str(sc):>25s} | {len(h4s):>6d} | {len(h6s):>6d} | {h4_str:>20s} | {h6_str:>20s}")

else:
    print("  ALL degree-6 coefficients are ZERO at n=6!")
    print("  This means H = H_0 + H_2 + H_4 exactly.")


# ========================================================================
# ANALYSIS 5: CYCLE DECOMPOSITION BRIDGE AT n=6
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: CYCLE-FOURIER BRIDGE AT n=6")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m
EH6 = math.factorial(n) / (2 ** (n - 1))

# For each tournament, compute H, H_2, H_4, and cycle counts
print(f"\nn={n}: Computing H, cycle counts for all {total} tournaments...")

all_data = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    H2 = compute_H2(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    # 3-cycles (directed)
    c3_dir = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            c3_dir += 1
        if A[a][c] and A[c][b] and A[b][a]:
            c3_dir += 1

    # 5-cycles (directed, on each 5-subset)
    c5_dir = 0
    for subset in combinations(range(n), 5):
        c5_dir += count_directed_cycles_on_subset(A, n, subset)

    all_data.append({
        'bits': bits, 'H': H, 'H2': H2, 'scores': scores,
        'c3_dir': c3_dir, 'c5_dir': c5_dir
    })

# H_4 = H - EH - H_2 (since H_6 = 0 at n=6 if confirmed, otherwise approximate)
# Let's check: correlate H_4 with c5_dir
print(f"\nCorrelation analysis:")

# Within each score class
for sc in sorted(set(d['scores'] for d in all_data)):
    group = [d for d in all_data if d['scores'] == sc]
    if len(group) < 3:
        continue

    H4_vals = [d['H'] - EH6 - d['H2'] for d in group]
    c3_vals = [d['c3_dir'] for d in group]
    c5_vals = [d['c5_dir'] for d in group]

    # Check if c3 is constant (score-determined)
    c3_unique = set(c3_vals)
    c5_unique = set(c5_vals)

    if len(c5_unique) == 1:
        continue  # No variation to correlate

    # Compute correlation of H_4 with c5_dir
    n_g = len(group)
    mean_h4 = sum(H4_vals) / n_g
    mean_c5 = sum(c5_vals) / n_g
    cov = sum((h - mean_h4) * (c - mean_c5) for h, c in zip(H4_vals, c5_vals)) / n_g
    var_h4 = sum((h - mean_h4)**2 for h in H4_vals) / n_g
    var_c5 = sum((c - mean_c5)**2 for c in c5_vals) / n_g

    if var_h4 > 0 and var_c5 > 0:
        corr = cov / (var_h4**0.5 * var_c5**0.5)
        print(f"  Score {sc}: c3 {'CONSTANT' if len(c3_unique)==1 else 'VARIES'} ({sorted(c3_unique)}), "
              f"c5 varies ({min(c5_unique)}-{max(c5_unique)}), "
              f"corr(H_4, c5_dir) = {corr:.4f}")


# ========================================================================
# ANALYSIS 6: THE BRIDGE EQUATION — DOES H_4 = f(c5)?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: THE BRIDGE EQUATION H_4 = f(c5) WITHIN SCORE CLASSES")
print("=" * 70)

for sc in sorted(set(d['scores'] for d in all_data)):
    group = [d for d in all_data if d['scores'] == sc]
    H4_vals = [round(d['H'] - EH6 - d['H2'], 4) for d in group]
    c5_vals = [d['c5_dir'] for d in group]

    pairs = set(zip(c5_vals, H4_vals))
    c5_to_H4 = defaultdict(set)
    for c5, h4 in pairs:
        c5_to_H4[c5].add(h4)

    if len(set(c5_vals)) <= 1:
        continue

    # Is H_4 a function of c5_dir?
    is_func = all(len(v) == 1 for v in c5_to_H4.values())

    print(f"\n  Score {sc}:")
    print(f"    c5_dir determines H_4? {'YES' if is_func else 'NO'}")
    for c5 in sorted(c5_to_H4.keys()):
        h4s = sorted(c5_to_H4[c5])
        print(f"    c5_dir={c5}: H_4 = {h4s}")

    if is_func and len(c5_to_H4) >= 2:
        # Fit linear: H_4 = a * c5_dir + b
        pts = [(c5, list(h4s)[0]) for c5, h4s in c5_to_H4.items()]
        pts.sort()
        if len(pts) >= 2:
            x1, y1 = pts[0]
            x2, y2 = pts[-1]
            slope = (y2 - y1) / (x2 - x1) if x2 != x1 else 0
            intercept = y1 - slope * x1
            print(f"    Linear fit: H_4 = {slope:.4f} * c5_dir + {intercept:.4f}")

            # Check if all points fit
            all_fit = all(abs(y - (slope * x + intercept)) < 0.01 for x, y in pts)
            print(f"    All points on line? {all_fit}")


# ========================================================================
# ANALYSIS 7: H_6 AND 7-CYCLE BRIDGE AT n=7 REGULAR
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: H_6 AND 7-CYCLE BRIDGE FOR REGULAR n=7")
print("=" * 70)

n = 7
EH7 = math.factorial(n) / (2 ** (n - 1))
c2 = math.factorial(n - 2) / (2 ** (n - 2))
H2_reg = c2 * n * (n - 1) / 2

print(f"\nFor regular n=7: E[H] = {EH7:.4f}, H_2 = {H2_reg:.4f}")
print(f"H_0 + H_2 = {EH7 + H2_reg:.4f}")

# From the three rigid classes:
# H=171: c5_dir=72, c7_dir=?
# H=175: c5_dir=56, c7_dir=?
# H=189: c5_dir=84, c7_dir=?

# Actually, let's get exact values from our data
if regular_data:
    classes = {}
    for d in regular_data:
        H = d['H']
        if H not in classes:
            classes[H] = d

    print(f"\nRegular tournament classes (representatives):")
    for H in sorted(classes.keys()):
        d = classes[H]
        H46 = H - EH7 - H2_reg
        c5 = d['c5_dir']
        c7 = d['c7_dir']
        print(f"\n  H = {H}:")
        print(f"    H_4 + H_6 = {H46:.4f}")
        print(f"    c5_dir = {c5}, c7_dir = {c7}")
        print(f"    Total directed cycles = {d['c3_dir']} + {c5} + {c7} = {d['c3_dir'] + c5 + c7}")

    # Try to separate H_4 and H_6 using the n=5 pattern
    # At n=5: H_4 = 2*c5_dir - 3 within score class
    # At n=7 regular: if H_4 = a*c5_dir + b and H_6 = c*c7_dir + d,
    # then H_4 + H_6 = a*c5_dir + c*c7_dir + b + d
    # We have 3 equations and 4 unknowns (a, b, c, d) — underdetermined!
    # But we can try: does a 2-param model H_4+H_6 = a*c5_dir + c*c7_dir fit?

    H_list = sorted(classes.keys())
    if len(H_list) == 3:
        data_pts = []
        for H in H_list:
            d = classes[H]
            H46 = H - EH7 - H2_reg
            data_pts.append((d['c5_dir'], d['c7_dir'], H46))

        print(f"\n  Fitting H_4 + H_6 = a*c5_dir + b*c7_dir + c:")
        c5_list = [p[0] for p in data_pts]
        c7_list = [p[1] for p in data_pts]
        H46_list = [p[2] for p in data_pts]

        # 3 equations, 3 unknowns (a, b, c)
        # H46[i] = a * c5[i] + b * c7[i] + c
        # Solve by matrix inversion
        import numpy as np
        try:
            M = np.array([[c5_list[i], c7_list[i], 1] for i in range(3)])
            rhs = np.array(H46_list)
            coeffs = np.linalg.solve(M, rhs)
            print(f"    a (c5 coeff) = {coeffs[0]:.6f}")
            print(f"    b (c7 coeff) = {coeffs[1]:.6f}")
            print(f"    c (constant) = {coeffs[2]:.6f}")

            # Verify
            for i, H in enumerate(H_list):
                pred = coeffs[0] * c5_list[i] + coeffs[1] * c7_list[i] + coeffs[2]
                actual = H46_list[i]
                print(f"    H={H}: predicted={pred:.4f}, actual={actual:.4f}, diff={abs(pred-actual):.6f}")

        except ImportError:
            # Manual solve for 3x3 system
            print("    (numpy not available, solving manually)")
            # Eliminate c using differences
            dc5_01 = c5_list[1] - c5_list[0]
            dc7_01 = c7_list[1] - c7_list[0]
            dH_01 = H46_list[1] - H46_list[0]

            dc5_02 = c5_list[2] - c5_list[0]
            dc7_02 = c7_list[2] - c7_list[0]
            dH_02 = H46_list[2] - H46_list[0]

            # a * dc5_01 + b * dc7_01 = dH_01
            # a * dc5_02 + b * dc7_02 = dH_02
            det = dc5_01 * dc7_02 - dc5_02 * dc7_01
            if abs(det) > 1e-10:
                a = (dH_01 * dc7_02 - dH_02 * dc7_01) / det
                b = (dc5_01 * dH_02 - dc5_02 * dH_01) / det
                c_const = H46_list[0] - a * c5_list[0] - b * c7_list[0]
                print(f"    a (c5 coeff) = {a:.6f}")
                print(f"    b (c7 coeff) = {b:.6f}")
                print(f"    c (constant) = {c_const:.6f}")

                for i, H in enumerate(H_list):
                    pred = a * c5_list[i] + b * c7_list[i] + c_const
                    actual = H46_list[i]
                    print(f"    H={H}: predicted={pred:.4f}, actual={actual:.4f}")
            else:
                print(f"    System is degenerate (det={det})")


# ========================================================================
# ANALYSIS 8: THE VITALI HIERARCHY — NON-MEASURABILITY GROWTH
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 8: VITALI HIERARCHY — NON-MEASURABILITY GROWTH WITH n")
print("=" * 70)

print("""
The "non-measurable" fraction of tournament structure grows with n:

At each n, let:
  - N_score = number of distinct score sequences
  - N_H = number of distinct H values
  - N_total = number of distinct tournaments (up to isomorphism)

The "measurability ratio" = N_score / N_total tells us what fraction
of structure is score-determined. As n grows, this ratio -> 0.
""")

for n in range(3, 8):
    m = n * (n - 1) // 2
    total_raw = 1 << m

    if total_raw > 2**21:
        print(f"n={n}: (too large for exhaustive)")
        continue

    scores_seen = set()
    H_seen = set()
    # For isomorphism classes, we'd need canonical forms — skip for now
    score_to_H = defaultdict(set)

    for bits in range(total_raw):
        A = binary_to_tournament(bits, n)
        scores = tuple(sorted([sum(A[v]) for v in range(n)]))
        H = count_ham_paths(A, n)
        scores_seen.add(scores)
        H_seen.add(H)
        score_to_H[scores].add(H)

    # Count how many score classes have >1 H value
    multi_H = sum(1 for sc, hs in score_to_H.items() if len(hs) > 1)

    print(f"n={n}: {len(scores_seen)} score classes, {len(H_seen)} H values, "
          f"{multi_H} score classes with multiple H values "
          f"({100*multi_H/len(scores_seen):.1f}%)")

    # The non-measurable dimension
    max_H_range = max(len(hs) for hs in score_to_H.values())
    print(f"  Max H values in a single score class: {max_H_range}")


# ========================================================================
# ANALYSIS 9: n=7 FULL SAMPLING — ALL SCORE CLASSES
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 9: n=7 SCORE CLASS STRUCTURE (SAMPLED)")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m
EH7 = math.factorial(n) / (2 ** (n - 1))

print(f"\nn={n}: Sampling 200000 tournaments...")

score_to_data = defaultdict(list)
sample_size = 200000
for _ in range(sample_size):
    bits = random.randint(0, total - 1)
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    H = count_ham_paths(A, n)
    H2 = compute_H2(A, n)
    H4_approx = H - EH7 - H2  # This is H_4 + H_6 actually
    score_to_data[scores].append({'H': H, 'H2': H2, 'H46': H4_approx})

print(f"Found {len(score_to_data)} score classes in sample")
print(f"\n  {'Score class':>25s} | {'#sampled':>8s} | {'#H vals':>7s} | {'H range':>15s} | {'H46 range':>15s}")
print(f"  {'':->25s}-+-{'':->8s}-+-{'':->7s}-+-{'':->15s}-+-{'':->15s}")

multi_count = 0
for sc in sorted(score_to_data.keys()):
    group = score_to_data[sc]
    H_vals = sorted(set(d['H'] for d in group))
    H46_vals = sorted(set(round(d['H46'], 2) for d in group))

    if len(H_vals) > 1:
        multi_count += 1

    H_range = f"{H_vals[0]}-{H_vals[-1]}" if len(H_vals) > 1 else f"{H_vals[0]}"
    H46_range = f"{H46_vals[0]:.1f}-{H46_vals[-1]:.1f}" if len(H46_vals) > 1 else f"{H46_vals[0]:.1f}"
    print(f"  {str(sc):>25s} | {len(group):>8d} | {len(H_vals):>7d} | {H_range:>15s} | {H46_range:>15s}")

print(f"\n  Score classes with multiple H values: {multi_count}/{len(score_to_data)} "
      f"({100*multi_count/len(score_to_data):.1f}%)")


print("\n" + "=" * 70)
print("DONE.")
