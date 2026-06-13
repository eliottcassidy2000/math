#!/usr/bin/env python3
"""
worpitzky_deeper.py — What invariants determine each Worpitzky coefficient?

KNOWN (from S46):
- c_{n-1} = n (universal)
- c_{n-2} = C(n,2) (universal)
- c_{n-3} = C(n,3) + 2(n-2)*t3
- c_{n-4} = C(n,4) + (n-2)(n-3)*t3  (at n=4,5)
- c_{n-5} at n=6: NOT determined by t3 alone

QUESTION: What invariant enters at the c_{n-5} level?
Candidate invariants:
- t3 = number of directed 3-cycles
- t4 = number of directed 4-cycles
- t3^2 = quadratic in t3
- w = some other function of score sequence
- sigma_k = k-th elementary symmetric function of scores

Author: opus-2026-03-07-S46
"""
from itertools import permutations, combinations
from math import comb, factorial
import numpy as np

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

def count_kcycles(adj, n, k):
    """Count directed k-cycles (each counted once, not k times)."""
    count = 0
    for combo in combinations(range(n), k):
        # Check all cyclic orderings
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%k]] for i in range(k)):
                count += 1
    # Each k-cycle is counted k times (cyclic rotations) * 1 (direction is fixed)
    return count // k

def worpitzky_a(F, n, m):
    return sum(F[k] * comb(m + n - 1 - k, n - 1) for k in range(n))

def get_worpitzky_coeffs(F, n):
    m_points = list(range(n + 2))
    a_vals = [worpitzky_a(F, n, m) for m in m_points]
    return np.polyfit(m_points, a_vals, n - 1)

# ============================================================
# n=6: What determines c_1 and c_0?
# ============================================================
print("=" * 60)
print("n=6: FINDING THE ADDITIONAL INVARIANT")
print("=" * 60)

n = 6
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

    t3 = count_kcycles(adj, n, 3)
    t4 = count_kcycles(adj, n, 4)

    # Score sequence
    scores = tuple(sorted([sum(adj[i][j] for j in range(n) if j != i) for i in range(n)]))

    # Sum of squares of scores (= variance * n + constant)
    s2 = sum(s**2 for s in scores)

    # tr(A^3) = 3*t3 (directed 3-cycles counted with multiplicity 3)
    A = np.array(adj, dtype=float)
    tr3 = int(round(np.trace(A @ A @ A)))
    tr4 = int(round(np.trace(A @ A @ A @ A)))

    coeffs = get_worpitzky_coeffs(F, n)
    # delta from binomial
    delta = [round(coeffs[i] - comb(n, i+1), 4) for i in range(n)]

    data.append({
        'F': F, 't3': t3, 't4': t4, 'scores': scores,
        's2': s2, 'tr3': tr3, 'tr4': tr4,
        'coeffs': coeffs, 'delta': delta
    })

print(f"{len(data)} distinct F-vectors\n")

# First, verify the universal and t3-determined levels
print("DELTA from binomial C(n, n-j):")
print(f"{'F':<32} {'t3':>3} {'t4':>3} {'s2':>3} {'d0':>6} {'d1':>6} {'d2':>6} {'d3':>6} {'d4':>6} {'d5':>6}")
for d in sorted(data, key=lambda x: (x['t3'], x['t4'])):
    F_str = str(d['F'])[:30]
    dd = d['delta']
    print(f"{F_str:<32} {d['t3']:>3} {d['t4']:>3} {d['s2']:>3} {dd[0]:>6.0f} {dd[1]:>6.0f} {dd[2]:>6.0f} {dd[3]:>6.0f} {dd[4]:>6.0f} {dd[5]:>6.0f}")

# Check which invariants determine delta[3] (= coeff of m^2 deviation)
print("\n\nChecking delta[3] (c_2 deviation):")
# This should be 12*t3
for d in data:
    pred = 12 * d['t3']
    actual = round(d['delta'][3])
    if pred != actual:
        print(f"  MISMATCH: t3={d['t3']}, pred={pred}, actual={actual}")
print("  All match 12*t3? ", all(round(d['delta'][3]) == 12*d['t3'] for d in data))

# Check delta[4] (c_1 deviation) — what determines it?
print("\n\nChecking delta[4] (c_1 deviation):")
# Group by (t3, t4)
from collections import defaultdict
t3t4_to_d4 = defaultdict(set)
for d in data:
    t3t4_to_d4[(d['t3'], d['t4'])].add(round(d['delta'][4]))

print("  (t3, t4) → delta[4]:")
for key in sorted(t3t4_to_d4.keys()):
    vals = t3t4_to_d4[key]
    print(f"    {key}: {sorted(vals)}")

# Check if (t3, s2) determines it
t3s2_to_d4 = defaultdict(set)
for d in data:
    t3s2_to_d4[(d['t3'], d['s2'])].add(round(d['delta'][4]))

print("\n  (t3, s2) → delta[4]:")
for key in sorted(t3s2_to_d4.keys()):
    vals = t3s2_to_d4[key]
    print(f"    {key}: {sorted(vals)}")

# Check if (t3, tr4) determines it
t3tr4_to_d4 = defaultdict(set)
for d in data:
    t3tr4_to_d4[(d['t3'], d['tr4'])].add(round(d['delta'][4]))

print("\n  (t3, tr4) → delta[4]:")
for key in sorted(t3tr4_to_d4.keys()):
    vals = t3tr4_to_d4[key]
    print(f"    {key}: {sorted(vals)}")

# Try linear regression: delta[4] = alpha + beta*t3 + gamma*t4
print("\n\nLinear regression: delta[4] = a + b*t3 + c*t4")
X = np.array([[1, d['t3'], d['t4']] for d in data])
y = np.array([d['delta'][4] for d in data])
coeffs_fit, residuals, _, _ = np.linalg.lstsq(X, y, rcond=None)
print(f"  a={coeffs_fit[0]:.4f}, b={coeffs_fit[1]:.4f}, c={coeffs_fit[2]:.4f}")
preds = X @ coeffs_fit
max_err = max(abs(preds - y))
print(f"  max error: {max_err:.6f}")

# Try with more invariants
print("\n  Adding s2:")
X2 = np.array([[1, d['t3'], d['t4'], d['s2']] for d in data])
coeffs_fit2, residuals2, _, _ = np.linalg.lstsq(X2, y, rcond=None)
print(f"  a={coeffs_fit2[0]:.4f}, b={coeffs_fit2[1]:.4f}, c={coeffs_fit2[2]:.4f}, d={coeffs_fit2[3]:.4f}")
preds2 = X2 @ coeffs_fit2
max_err2 = max(abs(preds2 - y))
print(f"  max error: {max_err2:.6f}")

# Try with t3^2
print("\n  Adding t3^2:")
X3 = np.array([[1, d['t3'], d['t3']**2] for d in data])
coeffs_fit3, _, _, _ = np.linalg.lstsq(X3, y, rcond=None)
print(f"  a={coeffs_fit3[0]:.4f}, b={coeffs_fit3[1]:.4f}, c={coeffs_fit3[2]:.4f}")
preds3 = X3 @ coeffs_fit3
max_err3 = max(abs(preds3 - y))
print(f"  max error: {max_err3:.6f}")

# Try (t3, t4, t3^2)
print("\n  Adding t3, t4, t3^2:")
X4 = np.array([[1, d['t3'], d['t4'], d['t3']**2] for d in data])
coeffs_fit4, _, _, _ = np.linalg.lstsq(X4, y, rcond=None)
print(f"  coeffs: {[f'{c:.4f}' for c in coeffs_fit4]}")
preds4 = X4 @ coeffs_fit4
max_err4 = max(abs(preds4 - y))
print(f"  max error: {max_err4:.6f}")

# Try tr4
print("\n  Using tr4 alone:")
X5 = np.array([[1, d['tr4']] for d in data])
coeffs_fit5, _, _, _ = np.linalg.lstsq(X5, y, rcond=None)
print(f"  a={coeffs_fit5[0]:.4f}, b={coeffs_fit5[1]:.4f}")
preds5 = X5 @ coeffs_fit5
max_err5 = max(abs(preds5 - y))
print(f"  max error: {max_err5:.6f}")

# Try (t3, tr4)
print("\n  Using (t3, tr4):")
X6 = np.array([[1, d['t3'], d['tr4']] for d in data])
coeffs_fit6, _, _, _ = np.linalg.lstsq(X6, y, rcond=None)
print(f"  a={coeffs_fit6[0]:.4f}, b={coeffs_fit6[1]:.4f}, c={coeffs_fit6[2]:.4f}")
preds6 = X6 @ coeffs_fit6
max_err6 = max(abs(preds6 - y))
print(f"  max error: {max_err6:.6f}")

# ============================================================
# FORMULA VERIFICATION: c_j = C(n,n-j) + sum of corrections
# ============================================================
print("\n\n" + "=" * 60)
print("FORMULA FOR t3 COEFFICIENT AT LEVEL d")
print("=" * 60)

for n in [4, 5, 6]:
    # At transitive tournament, delta = 0 at all levels
    # At level d (coefficient of m^{n-1-d}):
    # delta_{n-1-d} = alpha_d * t3 + ...

    # Level d=2: alpha_2 = 2(n-2)
    # Level d=3: alpha_3 = (n-2)(n-3)

    expected_alpha = [
        0,           # d=0 (universal)
        0,           # d=1 (universal)
        2*(n-2),     # d=2
        (n-2)*(n-3), # d=3
    ]

    print(f"\nn={n}: Expected t3 coefficients: d=2: {expected_alpha[2]}, d=3: {expected_alpha[3]}")

    # Verify at level d=3 when it's determined by t3
    if n <= 5:
        print(f"  (all levels ≤ n-2 determined by t3 at n={n})")
    else:
        print(f"  (levels beyond d=3 need additional invariants)")

# ============================================================
# BEAUTIFUL: For transitive, a_m = (m+1)^n - m^n
# ============================================================
print("\n\n" + "=" * 60)
print("TRANSITIVE TOURNAMENT: a_m = (m+1)^n - m^n")
print("=" * 60)

for n in [3, 4, 5, 6, 7]:
    # Transitive tournament: adj[i][j] = 1 iff i < j
    adj = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
    if n <= 6:
        F = compute_F(adj, n)
    else:
        # For n=7, compute F directly from Eulerian numbers
        F = [0]*n
        for P in permutations(range(n)):
            fwd = sum(1 for i in range(n-1) if P[i] < P[i+1])
            F[fwd] += 1

    print(f"\nn={n}: F = {F}")

    # Check: a_m = (m+1)^n - m^n
    for m in range(5):
        a_actual = worpitzky_a(F, n, m)
        a_expected = (m+1)**n - m**n
        match = "✓" if a_actual == a_expected else f"✗ ({a_actual} vs {a_expected})"
        print(f"  a({m}) = {a_actual}, (m+1)^n - m^n = {(m+1)**n - m**n}  {match}")
