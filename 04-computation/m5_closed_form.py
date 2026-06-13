#!/usr/bin/env python3
"""
Derive and verify closed-form m5 coefficients.

From m5_general_n.py output:
  t5/(n-4)!: 720, 960, 1200, 1440, 1680, 1920 at n=7,9,11,13,15,17
  -> 120*(n-1) (verified: 120*6=720, 120*8=960, etc.)

  bc = 2*t5 -> 240*(n-1)*(n-4)!

  t3/(n-2)!: 1320, 2960, 5600, 9480, 14840, 21920
  Third differences = 240 -> cubic: 5n^3 - 10n^2 + 15n - 10 = 5(n-1)(n^2-n+2)

opus-2026-03-06-S28
"""
from itertools import combinations
from math import factorial, comb
from collections import Counter
from fractions import Fraction
from functools import reduce
import numpy as np

def count_position_patterns(n, k):
    if k == 0: return {(): 1}
    patterns = Counter()
    for S in combinations(range(n-1), k):
        pos = sorted(S)
        comps, comp = [], [pos[0]]
        for i in range(1, len(pos)):
            if pos[i] == comp[-1] + 1:
                comp.append(pos[i])
            else:
                comps.append(len(comp))
                comp = [pos[i]]
        comps.append(len(comp))
        patterns[tuple(sorted(comps, reverse=True))] += 1
    return dict(patterns)

def sigma_full(pattern, n):
    sizes = list(pattern)
    total_verts = sum(s + 1 for s in sizes)
    free = n - total_verts
    if free < 0: return (0, 0, 0, 0)
    num_size1 = sizes.count(1)
    big_sizes = sorted([s for s in sizes if s > 1], reverse=True)
    F = factorial(free)
    def pair_prod(start):
        p, r = 1, start
        for i in range(num_size1):
            p *= comb(r, 2); r -= 2
        return p
    if len(big_sizes) == 0: return (factorial(n) // (2**len(sizes)), 0, 0, 0)
    if big_sizes == [2]:
        pp = pair_prod(n-3); return (F*comb(n,3)*pp, F*2*pp, 0, 0)
    if big_sizes == [3]:
        pp = pair_prod(n-4); return (F*comb(n,4)*pp, F*2*(n-3)*pp, 0, 0)
    if big_sizes == [4]:
        pp = pair_prod(n-5); return (F*comb(n,5)*pp, F*2*comb(n-3,2)*pp, F*2*pp, 0)
    if big_sizes == [5]:
        pp = pair_prod(n-6)
        return (F*comb(n,6)*pp, F*2*comb(n-3,3)*pp, F*2*(n-5)*pp, F*4*pp)
    if big_sizes == [2, 2]:
        pp = pair_prod(n-6)
        return (F*comb(n,3)*comb(n-3,3)*pp, F*4*comb(n-3,3)*pp, 0, F*8*pp)
    if big_sizes == [3, 2]:
        pp = pair_prod(n-7)
        const = comb(n,4)*comb(n-4,3)
        t3_c = 2*(n-3)*comb(n-4,3) + 2*comb(n-3,4)
        bc_c = 8*(n-6)
        return (F*const*pp, F*t3_c*pp, 0, F*bc_c*pp)
    return None

def compute_moment(j, n):
    stirling_table = {
        0: {0: 1},
        1: {1: 1},
        2: {1: 1, 2: 2},
        3: {1: 1, 2: 6, 3: 6},
        4: {1: 1, 2: 14, 3: 36, 4: 24},
        5: {1: 1, 2: 30, 3: 150, 4: 240, 5: 120},
        6: {1: 1, 2: 62, 3: 540, 4: 1560, 5: 1800, 6: 720},
    }
    wts = stirling_table[j]
    coeffs = [0, 0, 0, 0]
    for k in range(j + 1):
        wt = wts.get(k, 0)
        if wt == 0: continue
        for pat, cnt in count_position_patterns(n, k).items():
            r = sigma_full(pat, n)
            if r is None: return None
            for i in range(4): coeffs[i] += wt * cnt * r[i]
    return coeffs

# =====================================================================
# Verify t3, t5, bc closed forms for m5
# =====================================================================
print("=" * 70)
print("Verify m5 coefficient closed forms")
print("=" * 70)

print("\nt3_coeff = 5*(n-1)*(n^2-n+2)*(n-2)!")
print("t5_coeff = 120*(n-1)*(n-4)!")
print("bc_coeff = 240*(n-1)*(n-4)!")

for n in range(7, 22, 2):
    c = compute_moment(5, n)
    if c is None: continue

    pred_t3 = 5 * (n-1) * (n*n - n + 2) * factorial(n-2)
    pred_t5 = 120 * (n-1) * factorial(n-4)
    pred_bc = 240 * (n-1) * factorial(n-4)

    ok = (c[1] == pred_t3 and c[2] == pred_t5 and c[3] == pred_bc)
    if not ok:
        print(f"  n={n}: FAIL")
        if c[1] != pred_t3: print(f"    t3: got {c[1]}, expected {pred_t3}")
        if c[2] != pred_t5: print(f"    t5: got {c[2]}, expected {pred_t5}")
        if c[3] != pred_bc: print(f"    bc: got {c[3]}, expected {pred_bc}")
    else:
        print(f"  n={n:2d}: OK")

# =====================================================================
# Fit m5 constant
# =====================================================================
print(f"\n{'='*70}")
print("Fit m5 constant")
print(f"{'='*70}")

ns = list(range(7, 22, 2))
consts = {}
for n in ns:
    c = compute_moment(5, n)
    if c is not None:
        consts[n] = c[0]

# Find LCD of const/n!
print("\n  const/n! as fractions:")
fracs = {}
for n in sorted(consts.keys())[:8]:
    f = Fraction(consts[n], factorial(n))
    fracs[n] = f
    print(f"  n={n:2d}: {f}")

def lcm(a, b):
    from math import gcd
    return a * b // gcd(a, b)
denoms = [f.denominator for f in fracs.values()]
common_denom = reduce(lcm, denoms)
print(f"\n  LCD of denominators: {common_denom}")

print(f"\n  {common_denom}*const/n!:")
int_vals = {}
for n in sorted(fracs.keys()):
    val = fracs[n] * common_denom
    int_vals[n] = int(val)
    print(f"  n={n:2d}: {val}")

# Fit polynomial
ns_fit = sorted(int_vals.keys())
vals_fit = [int_vals[n] for n in ns_fit]

# Try various degrees
for deg in range(4, 7):
    if len(ns_fit) < deg + 1: continue
    A = np.array([[n**k for k in range(deg+1)] for n in ns_fit[:deg+1]], dtype=float)
    b = np.array(vals_fit[:deg+1], dtype=float)
    coeffs = np.linalg.solve(A, b)

    # Verify against remaining points
    all_ok = True
    for i in range(deg+1, len(ns_fit)):
        pred = sum(coeffs[k] * ns_fit[i]**k for k in range(deg+1))
        if abs(pred - vals_fit[i]) > 0.5:
            all_ok = False
            break

    if all_ok:
        print(f"\n  Degree {deg} polynomial fits all {len(ns_fit)} points!")
        for k in range(deg+1):
            frac = Fraction(coeffs[k]).limit_denominator(10000)
            print(f"    n^{k}: {coeffs[k]:.6f} ≈ {frac}")
        break
    else:
        print(f"\n  Degree {deg}: doesn't fit all points")

# =====================================================================
# Summary of ALL closed forms
# =====================================================================
print(f"\n{'='*70}")
print("COMPLETE SUMMARY: Forward-Arc Moments at General n")
print(f"{'='*70}")
print("""
THEOREM (General n Moment Hierarchy):

Let T be a tournament on n vertices with:
  t3 = directed 3-cycle count
  t5 = directed 5-cycle count
  bc = pairs of vertex-disjoint directed 3-cycles

Then m_j = sum_{P in S_n} f_P^j satisfies:

  m0 = n!
  m1 = n!(n-1)/2
  m2 = n!(3n^2-5n+4)/12 + 4(n-2)! * t3
  m3 = n!(n-1)(n^2-n+2)/8 + 6(n-1)! * t3
  m4 = n!(15n^4-30n^3+65n^2-82n+48)/240
       + 2(n-2)!(3n^2-5n+4) * t3
       + 48(n-4)! * t5
       + 96(n-4)! * bc
  m5 = C5(n)
       + 5(n-1)(n^2-n+2)(n-2)! * t3
       + 120(n-1)(n-4)! * t5
       + 240(n-1)(n-4)! * bc

STRUCTURAL PROPERTIES:
  1. bc_coeff = 2 * t5_coeff in BOTH m4 and m5
  2. Polynomial (n-1)(n^2-n+2) appears in m3_const AND m5_t3
  3. Polynomial 3n^2-5n+4 appears in m2_const AND m4_t3
  4. Pattern: m_{j+2}_t3 ∝ m_j_const / (n(n-1))
""")

# =====================================================================
# Cross-moment relations
# =====================================================================
print(f"\n{'='*70}")
print("Cross-moment structural relations")
print(f"{'='*70}")

# m4_t3 = 24 * m2_const / (n(n-1))
# m5_t3 = 40 * m3_const / (n(n-1))
print("\n  m4_t3 vs m2_const:")
for n in [5, 7, 9, 11, 13]:
    m2_c = factorial(n) * (3*n*n - 5*n + 4) // 12
    m4_t3 = 2 * factorial(n-2) * (3*n*n - 5*n + 4)
    ratio = Fraction(m4_t3 * n * (n-1), m2_c)
    print(f"  n={n}: m4_t3*n(n-1)/m2_const = {ratio}")

print("\n  m5_t3 vs m3_const:")
for n in [7, 9, 11, 13]:
    m3_c = factorial(n) * (n-1) * (n*n - n + 2) // 8
    m5_t3 = 5 * (n-1) * (n*n - n + 2) * factorial(n-2)
    ratio = Fraction(m5_t3 * n * (n-1), m3_c)
    print(f"  n={n}: m5_t3*n(n-1)/m3_const = {ratio}")

# m4_t5 = 48*(n-4)!
# m5_t5 = 120*(n-1)*(n-4)!
print("\n  m5_t5 vs m4_t5:")
for n in [7, 9, 11, 13]:
    m4_t5 = 48 * factorial(n-4)
    m5_t5 = 120 * (n-1) * factorial(n-4)
    ratio = Fraction(m5_t5, m4_t5)
    print(f"  n={n}: m5_t5/m4_t5 = {ratio} = {float(ratio):.4f}")

print("\nDONE")
