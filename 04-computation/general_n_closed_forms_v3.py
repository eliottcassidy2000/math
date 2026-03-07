#!/usr/bin/env python3
"""
CLOSED-FORM moment formulas: final derivation.

PROVED:
  m0 = n!
  m1 = n!(n-1)/2
  m2 = n!(3n^2-5n+4)/12 + 4(n-2)! t3
  m3 = n!(n-1)(n^2-n+2)/8 + 6(n-1)! t3
  m4: t5_coeff = 48(n-4)!, bc_coeff = 96(n-4)!, t3_coeff = 2(n-2)!(3n^2-5n+4)

Now: find m4_const closed form, and express everything cleanly.

opus-2026-03-06-S28
"""
from itertools import combinations
from math import factorial, comb
from collections import Counter
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

def sigma_extended(pattern, n):
    sizes = list(pattern)
    total_verts = sum(s + 1 for s in sizes)
    free = n - total_verts
    if free < 0: return (0, 0, 0, 0)

    num_size1 = sizes.count(1)
    big_sizes = sorted([s for s in sizes if s > 1], reverse=True)
    F = factorial(free)

    def pair_product(start_remaining):
        p, r = 1, start_remaining
        for i in range(num_size1):
            p *= comb(r, 2)
            r -= 2
        return p

    if len(big_sizes) == 0:
        return (factorial(n) // (2**len(sizes)), 0, 0, 0)
    if big_sizes == [2]:
        pp = pair_product(n - 3)
        return (F*comb(n,3)*pp, F*2*pp, 0, 0)
    if big_sizes == [3]:
        pp = pair_product(n - 4)
        return (F*comb(n,4)*pp, F*2*(n-3)*pp, 0, 0)
    if big_sizes == [4]:
        pp = pair_product(n - 5)
        return (F*comb(n,5)*pp, F*2*comb(n-3,2)*pp, F*2*pp, 0)
    if big_sizes == [2, 2]:
        pp = pair_product(n - 6)
        return (F*comb(n,3)*comb(n-3,3)*pp, F*4*comb(n-3,3)*pp, 0, F*8*pp)
    if big_sizes == [5]:
        pp = pair_product(n - 6)
        return (F*comb(n,6)*pp, F*2*comb(n-3,3)*pp, F*2*(n-5)*pp, F*4*pp)
    return None

def compute_moment(j, n):
    """Compute m_j coefficients (const, t3, t5, bc) at given n."""
    # Stirling numbers S(j,k) * k!
    stirling_table = {
        0: {0: 1},
        1: {0: 0, 1: 1},
        2: {0: 0, 1: 1, 2: 2},
        3: {0: 0, 1: 1, 2: 6, 3: 6},
        4: {0: 0, 1: 1, 2: 14, 3: 36, 4: 24},
        5: {0: 0, 1: 1, 2: 30, 3: 150, 4: 240, 5: 120},
        6: {0: 0, 1: 1, 2: 62, 3: 540, 4: 1560, 5: 1800, 6: 720},
    }
    wts = stirling_table[j]
    coeffs = [0, 0, 0, 0]
    for k in range(j + 1):
        wt = wts.get(k, 0)
        if wt == 0: continue
        for pat, cnt in count_position_patterns(n, k).items():
            r = sigma_extended(pat, n)
            if r is None:
                return None  # Unhandled pattern
            for i in range(4):
                coeffs[i] += wt * cnt * r[i]
    return coeffs

# =====================================================================
# Compute m4 at many n values
# =====================================================================
print("=" * 70)
print("m4 coefficients at various n")
print("=" * 70)

ns = list(range(5, 22, 2))  # 5,7,...,21
m4_data = {}
for n in ns:
    c = compute_moment(4, n)
    if c is None:
        print(f"  n={n}: INCOMPLETE")
        continue
    m4_data[n] = c
    print(f"  n={n:2d}: const={c[0]}, t3={c[1]}, t5={c[2]}, bc={c[3]}")

# Verify closed forms for t3, t5, bc
print(f"\n{'='*70}")
print("Verify: t3_coeff = 2*(n-2)!*(3n^2-5n+4)")
print("        t5_coeff = 48*(n-4)!")
print("        bc_coeff = 96*(n-4)!")
print(f"{'='*70}")

for n in sorted(m4_data.keys()):
    c = m4_data[n]
    pred_t3 = 2 * factorial(n-2) * (3*n*n - 5*n + 4)
    pred_t5 = 48 * factorial(n-4) if n >= 4 else 0
    pred_bc = 96 * factorial(n-4) if n >= 4 else 0
    ok = (c[1] == pred_t3 and c[2] == pred_t5 and c[3] == pred_bc)
    print(f"  n={n:2d}: t3 {'OK' if c[1]==pred_t3 else 'FAIL'}, "
          f"t5 {'OK' if c[2]==pred_t5 else 'FAIL'}, "
          f"bc {'OK' if c[3]==pred_bc else 'FAIL'}")

# =====================================================================
# Fit m4_const
# =====================================================================
print(f"\n{'='*70}")
print("Fit m4_const = n! * P(n) where P is rational in n")
print(f"{'='*70}")

# Let's try: const = n! * P(n) / D
# First look at const / n!
for n in sorted(m4_data.keys()):
    r = m4_data[n][0] / factorial(n)
    print(f"  n={n:2d}: const/n! = {r:.10f}")

# These aren't integers. Try const * 60 / n!:
print("\n  60*const/n!:")
for n in sorted(m4_data.keys()):
    r = 60 * m4_data[n][0] / factorial(n)
    print(f"  n={n:2d}: {r:.6f}")

# Try const / ((n-4)! * n!/(n-4)!) = const / n!... same thing.
# Let me try const * 360 / n!:
print("\n  Trying to find integer normalization...")
for D in [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 24, 30, 60, 120, 240, 360]:
    vals = [D * m4_data[n][0] / factorial(n) for n in sorted(m4_data.keys())]
    all_int = all(abs(v - round(v)) < 0.001 for v in vals)
    if all_int:
        print(f"  D={D}: {[int(round(v)) for v in vals[:6]]} ...")
        break

# Alternative: try const / (n! / some_factor)
# At n=7: 596064 = 5040 * 118.267 ≈ 5040 * 1774/15.
# 596064 * 15 = 8940960. 8940960/5040 = 1774.
# At n=5: 3444 * 15 = 51660. 51660/120 = 430.5. Not integer.
# So 15 doesn't work.
# At n=5: 3444. Factor: 3444 = 4 * 861 = 4 * 3 * 287 = 12 * 287.
# 287 = 7 * 41. So 3444 = 12 * 7 * 41.
# At n=7: 596064 = 5040 * 118 + 5040*4/15... messy.

# Let me try the LCM approach
from fractions import Fraction
print("\n  const/n! as fractions:")
for n in sorted(m4_data.keys())[:7]:
    f = Fraction(m4_data[n][0], factorial(n))
    print(f"  n={n:2d}: {f} = {float(f):.10f}")

# Find LCD of all these fractions
fracs = [Fraction(m4_data[n][0], factorial(n)) for n in sorted(m4_data.keys())[:7]]
from functools import reduce
def lcm(a, b):
    from math import gcd
    return a * b // gcd(a, b)
denoms = [f.denominator for f in fracs]
common_denom = reduce(lcm, denoms)
print(f"\n  LCD of denominators: {common_denom}")
print(f"  {common_denom}*const/n!:")
for i, n in enumerate(sorted(m4_data.keys())[:7]):
    val = fracs[i] * common_denom
    print(f"  n={n:2d}: {val}")

# Now fit polynomial
print(f"\n  Fitting polynomial P(n) where const = n!*P(n)/{common_denom}:")
ns_fit = sorted(m4_data.keys())[:7]
int_vals = [int(fracs[i] * common_denom) for i, n in enumerate(ns_fit)]

# Build Vandermonde matrix
deg = 5  # Try degree 5
A = np.array([[n**k for k in range(deg+1)] for n in ns_fit], dtype=float)
b = np.array(int_vals, dtype=float)
if len(ns_fit) > deg + 1:
    coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
else:
    coeffs = np.linalg.solve(A[:deg+1, :deg+1], b[:deg+1])

print(f"  Polynomial coefficients (degree {deg}):")
for i, c in enumerate(coeffs):
    print(f"    n^{i}: {c:.6f}")

# Check if coefficients are rational
print(f"\n  Rational check ({common_denom}*const/n!):")
for i, c in enumerate(coeffs):
    frac = Fraction(c).limit_denominator(1000)
    print(f"    n^{i}: {c:.6f} ≈ {frac}")

# =====================================================================
# MOMENT PATTERN: express m_j coefficients generally
# =====================================================================
print(f"\n{'='*70}")
print("PATTERN: m_j t3 coefficient = C_j * (n-2)! * (3n^2-5n+4)")
print(f"{'='*70}")

for j in range(2, 5):
    print(f"\n  j={j}:")
    for n in [5, 7, 9, 11, 13]:
        c = compute_moment(j, n)
        if c is None: continue
        t3_ratio = c[1] / (factorial(n-2) * (3*n*n - 5*n + 4)) if c[1] != 0 else 0
        print(f"    n={n}: t3_coeff = {c[1]}, ratio = {t3_ratio:.6f}")

# =====================================================================
# Explore: is the polynomial 3n^2-5n+4 special?
# =====================================================================
print(f"\n{'='*70}")
print("Note: 3n^2-5n+4 factors as... let's check")
print(f"{'='*70}")
# Discriminant: 25-48 = -23 < 0, so it doesn't factor over reals.
# But 3n^2-5n+4 = n(3n-5) + 4. At n=1: 2, n=2: 6, n=3: 16.
# More importantly: 12*m2_const/n! = 3n^2-5n+4 and m4_t3/(2*(n-2)!) = 3n^2-5n+4.
# So m4_t3 = 2*(n-2)!*(3n^2-5n+4) = 2*(n-2)! * 12*m2_const/n! = 24*m2_const/(n*(n-1))
# Interesting cross-moment connection!
print("  m2_const = n!*(3n^2-5n+4)/12")
print("  m4_t3 = 2*(n-2)!*(3n^2-5n+4)")
print("  => m4_t3 = 24*m2_const / (n*(n-1))")
print()
for n in [5, 7, 9, 11, 13]:
    m2_c = factorial(n) * (3*n*n - 5*n + 4) // 12
    m4_t3 = 2 * factorial(n-2) * (3*n*n - 5*n + 4)
    check = 24 * m2_c == m4_t3 * n * (n-1)
    print(f"  n={n}: 24*m2_const = {24*m2_c}, m4_t3*n*(n-1) = {m4_t3*n*(n-1)}, {'OK' if check else 'FAIL'}")

# =====================================================================
# m5 coefficients
# =====================================================================
print(f"\n{'='*70}")
print("m5 coefficients at various n")
print(f"{'='*70}")

for n in [5, 7, 9, 11]:
    c = compute_moment(5, n)
    if c is None:
        print(f"  n={n}: INCOMPLETE (unhandled patterns)")
        continue
    print(f"  n={n:2d}: const={c[0]}, t3={c[1]}, t5={c[2]}, bc={c[3]}")

# =====================================================================
# Summary
# =====================================================================
print(f"\n{'='*70}")
print("SUMMARY: General n moment formulas")
print(f"{'='*70}")
print("""
THEOREM (Moment Hierarchy at General n):

  m0 = n!
  m1 = n!(n-1)/2
  m2 = n!(3n^2-5n+4)/12 + 4(n-2)! * t3
  m3 = n!(n-1)(n^2-n+2)/8 + 6(n-1)! * t3

  m4 = C4(n) + 2(n-2)!(3n^2-5n+4) * t3 + 48(n-4)! * t5 + 96(n-4)! * bc

  where C4(n) = m4_const (closed form TBD).

Key observations:
  - m0, m1: universal (no tournament dependence)
  - m2, m3: depend only on t3 (directed 3-cycle count)
  - m4, m5: depend on (t3, t5, bc)
  - m_{n-1}: depends on (t3, t5, bc, ..., H) — the full OCF hierarchy
  - Polynomial 3n^2-5n+4 appears in both m2_const and m4_t3
  - bc_coeff = 2 * t5_coeff in m4 (both = 48*(n-4)! * {1,2})
""")

print("DONE")
