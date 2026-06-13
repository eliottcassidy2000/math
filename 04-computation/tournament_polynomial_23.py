#!/usr/bin/env python3
"""
tournament_polynomial_23.py — opus-2026-03-14-S71e

THE POLYNOMIAL z^2 - 5z + 6 = (z-2)(z-3) = 0

This polynomial encodes "the universe" of tournament theory at n<=8.

Key properties:
  Roots: 2 and 3 (the evaluation points)
  Sum of roots: 5 (the first prime that enters the Vandermonde at k=3)
  Product of roots: 6 (the Vandermonde determinant V(2,3))
  Discriminant: 25-24 = 1 (PERFECT SQUARE → rational roots)

Connection to k-nacci:
  k-nacci characteristic: x^{k+1} - 2x^k + 1 = 0
  Doubled k-nacci:        x^{k+1} - 3x^k + 2 = 0

  These are x^{k+1} - r*x^k + (r-1) = 0 for r = 2, 3.
  Factor: (x-1)(x^k - (r-1)*x^{k-1} - ... - (r-1)) = 0.
  So x = 1 is always a root, and the dominant root → r.

  The polynomial P(r) = x^{k+1} - r*x^k + (r-1) vanishes at x=1
  for ALL r. So the k-nacci family interpolates between 2 and 3
  along the parameter r.

  For r = 2 + 3 = 5: x^{k+1} - 5x^k + 4 = 0 → dominant root → 5?
  No, (x-1)(x^k - 4x^{k-1} - ... - 4) = 0, root → 5.

  GENERAL: For parameter r, the r-nacci has x^{k+1}-rx^k+(r-1)=0,
  dominant root → r as k → ∞.

Exploration: What invariant does I(5) represent? I(r) for general r?
"""

import sys
import numpy as np
from itertools import combinations
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {(1 << 0, 0): 1}
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

# ======================================================================
# PART 1: The tournament polynomial z^2 - 5z + 6
# ======================================================================
print("=" * 70)
print("PART 1: THE TOURNAMENT POLYNOMIAL z^2 - 5z + 6 = (z-2)(z-3)")
print("=" * 70)

print("""
  For a tournament T on n <= 8 vertices:
    I(Omega(T), x) = 1 + alpha_1 * x + alpha_2 * x^2

  This is a polynomial in x of degree 2.

  The TWO roots of I(x) = 0 are:
    x = (-alpha_1 +/- sqrt(alpha_1^2 - 4*alpha_2)) / (2*alpha_2)

  But more naturally, we evaluate at INTEGERS x = 2, 3.
  The polynomial z^2 - 5z + 6 = (z-2)(z-3) encodes these evaluation points.

  KEY RELATIONSHIPS:
    z1 + z2 = 5     (sum of roots = 5)
    z1 * z2 = 6     (product = 6 = Vandermonde V(2,3))
    z1 - z2 = -1    (difference = -1)
    |z1 - z2| = 1   (gap = 1)

  The Vandermonde det = z1*z2*(z2-z1) = 6*1 = 6.

  I(z1) = H, I(z2) = I_3.
  So: I evaluated at the roots of z^2-5z+6 gives (H, I_3).
""")

# ======================================================================
# PART 2: I(x) for various x — what do they mean?
# ======================================================================
print("=" * 70)
print("PART 2: I(x) AT VARIOUS INTEGER POINTS")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)
print(f"\n  n={n}: 200 random tournaments")

data = []
for trial in range(200):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7
    a2 = (H - 1 - 2*a1) // 4

    Ix = {}
    for x in [-2, -1, 0, 1, 2, 3, 4, 5]:
        Ix[x] = 1 + a1 * x + a2 * x**2
    data.append({'a1': a1, 'a2': a2, 'H': H, 'Ix': Ix})

print(f"\n  I(x) statistics:")
for x in [-2, -1, 0, 1, 2, 3, 4, 5]:
    vals = [d['Ix'][x] for d in data]
    print(f"    I({x:2d}): mean={np.mean(vals):10.1f}, range=[{min(vals):6d}, {max(vals):6d}]")

# ======================================================================
# PART 3: Factored form of I(x)
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: FACTORED FORM OF I(x)")
print("=" * 70)

print("""
  I(x) = 1 + alpha_1*x + alpha_2*x^2 = alpha_2*(x - r1)*(x - r2)
  where r1, r2 are the roots of I(x) = 0.

  If alpha_2 > 0: I(x) has two real roots (since discriminant
  = alpha_1^2 - 4*alpha_2 can be positive, zero, or negative).

  If alpha_2 = 0: I(x) = 1 + alpha_1*x is LINEAR, root at x = -1/alpha_1.
""")

# Analyze the roots of I(x) for our sample
print(f"\n  Roots of I(x) = 0 at n=7:")
root_data = []
for d in data[:20]:
    a1, a2 = d['a1'], d['a2']
    if a2 == 0:
        if a1 > 0:
            r = -1/a1
            print(f"    alpha1={a1:3d}, alpha2=0: LINEAR, root x={r:.4f}")
            root_data.append(('linear', r, None))
        else:
            print(f"    alpha1=0, alpha2=0: CONSTANT I(x)=1, no roots")
            root_data.append(('constant', None, None))
    else:
        disc = a1**2 - 4*a2
        if disc >= 0:
            r1 = (-a1 + disc**0.5) / (2*a2)
            r2 = (-a1 - disc**0.5) / (2*a2)
            print(f"    alpha1={a1:3d}, alpha2={a2:3d}: disc={disc:5d}, roots x={r1:.4f}, {r2:.4f}")
            root_data.append(('real', r1, r2))
        else:
            re = -a1 / (2*a2)
            im = (-disc)**0.5 / (2*a2)
            print(f"    alpha1={a1:3d}, alpha2={a2:3d}: disc={disc:5d}, COMPLEX roots {re:.4f} +/- {im:.4f}i")
            root_data.append(('complex', re, im))

# Distribution of discriminant sign
discs = []
for d in data:
    a1, a2 = d['a1'], d['a2']
    if a2 > 0:
        disc = a1**2 - 4*a2
        discs.append(disc)

discs = np.array(discs)
print(f"\n  Discriminant statistics (alpha_2 > 0 cases: {len(discs)}/200):")
print(f"    disc > 0 (real roots): {np.sum(discs > 0)}/{len(discs)}")
print(f"    disc = 0 (double root): {np.sum(discs == 0)}/{len(discs)}")
print(f"    disc < 0 (complex roots): {np.sum(discs < 0)}/{len(discs)}")

# ======================================================================
# PART 4: The I(x) = z^2 - 5z + 6 connection
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: WHEN DOES I(x) = z^2 - 5z + 6?")
print("=" * 70)

print("""
  z^2 - 5z + 6 = 6 - 5z + z^2
  Comparing with I(x) = 1 + alpha_1*x + alpha_2*x^2:
    alpha_0 = 6 (not 1!)

  So I(x) ≠ z^2-5z+6 directly. But we can ask:
  Is there a NORMALIZED form where the evaluation points 2,3
  are the roots?

  The polynomial with roots 2,3 and constant term 1:
    (1 - x/2)(1 - x/3) = 1 - 5x/6 + x^2/6

  So: I(x) = 1 + alpha_1*x + alpha_2*x^2
  would have roots at 2 and 3 iff:
    alpha_1 = -5/6, alpha_2 = 1/6

  This corresponds to alpha_1 = -5/6 (NEGATIVE — impossible since alpha_1
  counts cycles!) and alpha_2 = 1/6 (not an integer).

  So I(x) never literally has roots at 2 and 3. But there's a deeper
  connection: the RATIO I(3)/I(2) = I_3/H is the key invariant.
""")

# Compute I(3)/I(2) = I_3/H ratio
ratios = [d['Ix'][3] / d['Ix'][2] for d in data if d['Ix'][2] > 0]
print(f"  I(3)/I(2) = I_3/H statistics (200 samples at n=7):")
print(f"    mean: {np.mean(ratios):.4f}")
print(f"    range: [{min(ratios):.4f}, {max(ratios):.4f}]")
print(f"    For transitive T: I(3)/I(2) = 1/1 = 1")
print(f"    For I(x) = 1+ax+bx^2: I(3)/I(2) = (1+3a+9b)/(1+2a+4b)")

# I(3)/I(2) ≈ 3/2 for "typical" tournaments?
print(f"\n    Fraction I(3)/I(2) vs 3/2 = {3/2}:")
close_to_32 = sum(1 for r in ratios if abs(r - 1.5) < 0.1) / len(ratios)
print(f"    |I_3/H - 3/2| < 0.1: {close_to_32*100:.1f}%")

# ======================================================================
# PART 5: The r-nacci family and tournament evaluations
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: THE r-NACCI FAMILY")
print("=" * 70)

print("""
  For each integer r >= 2, define the r-nacci:
    x_n = r * (x_{n-1} + x_{n-2} + ... + x_{n-k}) / (r-1)

  Wait, let me be more precise.

  The characteristic equation x^{k+1} - r*x^k + (r-1) = 0
  has dominant root phi_k(r) → r as k → ∞.

  Error: phi_k(r) ≈ r - (r-1)/r^k

  This gives a FAMILY of convergences, one for each r:
    r=2: phi_k(2) → 2 at rate 1/2^k (standard k-nacci)
    r=3: phi_k(3) → 3 at rate 2/3^k (doubled k-nacci)
    r=4: phi_k(4) → 4 at rate 3/4^k (tripled k-nacci)
    r=5: phi_k(5) → 5 at rate 4/5^k (quadrupled k-nacci)

  Error(r) = (r-1)/r^k.
  Error ratio: err(r+1)/err(r) = (r/(r-1)) * (r/(r+1))^k

  For k=n and large n: err(r+1)/err(r) → 0, so EACH subsequent
  evaluation point converges FASTER than the previous.
  3 converges faster than 2, 4 faster than 3, etc.

  BUT: the Vandermonde extraction DIVIDES by det ∝ prod(r_i - r_j),
  which grows with the number of evaluation points. So the
  EXTRACTION precision DEGRADES even as the CONVERGENCE improves.

  This is the fundamental tension:
    - More evaluation points → more information (higher alpha_k)
    - More points → higher Vandermonde denominator → less precision
    - Net effect: the first two points (2,3) are always THE most
      informative, with diminishing returns for each additional point.
""")

# Verify the error formula err(r) = (r-1)/r^k
print("  Verification of error formula phi_k(r) ≈ r - (r-1)/r^k:")
for r in [2, 3, 4, 5]:
    for k in [5, 10, 15]:
        # Compute phi_k(r) numerically
        a = [0]*k + [1]
        for _ in range(300):
            a.append((r-1) * sum(a[-k:]))  # This gives x^{k+1}=r*x^k-(r-1)? No.
        # Actually need to find what recurrence gives x^{k+1}-rx^k+(r-1)=0
        # Factor: (x-1)(x^k - (r-1)(x^{k-1}+...+1)) = 0
        # So x^k = (r-1)(x^{k-1}+...+1)
        # Recurrence: x_n = (r-1)*(x_{n-1}+...+x_{n-k})
        a = [0]*k + [1]
        for _ in range(300):
            a.append(int((r-1)) * sum(a[-k:]))
        ratio = a[-1] / a[-2] if a[-2] != 0 else 0
        err = r - ratio
        theoretical = (r-1) / r**k
        match = abs(err - theoretical) / theoretical if theoretical > 1e-15 else 0
        if k <= 10:
            print(f"    r={r}, k={k:2d}: err={err:.8f}, (r-1)/r^k={theoretical:.8f}, match={1-match:.4f}")

# ======================================================================
# PART 6: The product I(2)*I(3) and the sum I(2)+I(3)
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: SYMMETRIC FUNCTIONS OF I(2) AND I(3)")
print("=" * 70)

# I(2)*I(3) = H * I_3 and I(2)+I(3) = H + I_3
# By the theory of symmetric functions, if z1=2, z2=3 are roots of
# z^2-5z+6=0, then for polynomial p(z) evaluated at these roots:
# p(2) + p(3) and p(2)*p(3) encode the symmetric functions.

# For I(x) = 1 + a1*x + a2*x^2:
# I(2) + I(3) = 2 + 5*a1 + 13*a2
# I(2) * I(3) = (1+2a1+4a2)(1+3a1+9a2)
#             = 1 + 5a1 + 13a2 + 6a1^2 + 30a1*a2 + 36a2^2

print("""
  H + I_3 = 2 + 5*alpha_1 + 13*alpha_2
  H * I_3 = (1+2a1+4a2)(1+3a1+9a2)
          = 1 + 5a1 + 13a2 + 6a1^2 + 30a1*a2 + 36a2^2

  The SUM: H + I_3 = 2 + 5*a1 + 13*a2
    Note: 5 = 2+3 (sum of roots), 13 = 4+9 = 2^2+3^2 (sum of squares)

  The PRODUCT: H * I_3 = 1 + 5a1 + 13a2 + 6(a1^2+5a1*a2+6a2^2)
    Note: 6 = 2*3 (product of roots)
    The quadratic part: 6(a1+2a2)(a1+3a2) = 6*I'(2)*I'(3)/4... not quite.
""")

# Check: H + I_3 mod 5
sums = [d['Ix'][2] + d['Ix'][3] for d in data]
prods = [d['Ix'][2] * d['Ix'][3] for d in data]

sum_mod5 = Counter(s % 5 for s in sums)
prod_mod6 = Counter(p % 6 for p in prods)

print(f"  H + I_3 mod 5: {dict(sum_mod5)}")
print(f"  H * I_3 mod 6: {dict(prod_mod6)}")

# H + I_3 = 2 + 5a1 + 13a2 ≡ 2 + 3a2 (mod 5)
# So H+I_3 mod 5 depends on a2 mod 5.
# H*I_3 = H*I_3. H ≡ 1 (mod 2), I_3 ≡ 1 (mod 3).
# H*I_3 mod 6: since H is odd and I_3 ≡ 1 (mod 3),
# H*I_3 mod 2 = 1 (odd), H*I_3 mod 3 = H mod 3.
# By CRT: H*I_3 mod 6 ∈ {1, 3, 5} (all odd).

print(f"\n  H + I_3 ≡ 2 + 3*alpha_2 (mod 5)")
for d in data[:10]:
    a2 = d['a2']
    s = d['Ix'][2] + d['Ix'][3]
    expected = (2 + 3*a2) % 5
    actual = s % 5
    print(f"    alpha_2={a2:3d}, H+I_3={s:5d}, mod5={actual} (expected {expected})")

# ======================================================================
# PART 7: The Newton identity connection
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: NEWTON IDENTITIES FOR THE 2-3 SYSTEM")
print("=" * 70)

print("""
  If we view 2 and 3 as "eigenvalues" of some operator, then:
    e_1 = 2 + 3 = 5 (1st elementary symmetric)
    e_2 = 2 * 3 = 6 (2nd elementary symmetric)
    p_1 = 2 + 3 = 5 (1st power sum)
    p_2 = 4 + 9 = 13 (2nd power sum)

  Newton's identities: p_k = e_1*p_{k-1} - e_2*p_{k-2}
    p_1 = 5
    p_2 = 5*5 - 2*6 = 25-12 = 13 ✓
    p_3 = 5*13 - 6*5 = 65-30 = 35 = 8+27 = 2^3+3^3 ✓
    p_4 = 5*35 - 6*13 = 175-78 = 97 = 16+81 = 2^4+3^4 ✓

  The power sums p_k = 2^k + 3^k satisfy the recurrence p_k = 5p_{k-1} - 6p_{k-2}.
  This is EXACTLY the characteristic equation z^2 - 5z + 6 = 0!

  Now: I(x) evaluated at integer x gives I(x) = 1 + a1*x + a2*x^2.
  The "power sum" S_k = I(2^k) + I(3^k)... no, that's not right.
  But the VALUES p_k = 2^k + 3^k do follow z^2-5z+6=0.

  CONNECTION TO TOURNAMENT EVALUATION:
    Define T_k = I(k) = 1 + a1*k + a2*k^2 for the tournament.
    Then T_2 = H, T_3 = I_3, T_4 = I_4, etc.

  The "step" T_{k+1} - T_k = a1 + (2k+1)*a2 (HYP-968).
  The "second step" (T_{k+2}-T_{k+1}) - (T_{k+1}-T_k) = 2*a2.
  So the second differences are CONSTANT = 2*alpha_2.

  This means: the sequence I(0), I(1), I(2), I(3), I(4), ...
  is a QUADRATIC sequence with constant second differences.

  I(0) = 1
  I(1) = 1 + a1 + a2
  I(2) = H = 1 + 2a1 + 4a2
  I(3) = I_3 = 1 + 3a1 + 9a2
  I(4) = 1 + 4a1 + 16a2

  Second differences: I(k+2) - 2I(k+1) + I(k) = 2a2 for all k.

  In particular:
    I(4) - 2I(3) + I(2) = 2*alpha_2
    H - 2I_3 + I(4) = ... no, H = I(2).
    I(4) - 2I(3) + H = 2*alpha_2

  And: I_3 - 2H + I(1) = 2*alpha_2
  So: I(1) = 2H - I_3 + 2*alpha_2 = 2H - I_3 + (2I_3 - 3H + 1)/3
  Hmm, let me just compute.
""")

# Verify second differences = 2*alpha_2
print("  Verification: I(k+2)-2I(k+1)+I(k) = 2*alpha_2")
for d in data[:5]:
    a1, a2 = d['a1'], d['a2']
    for k in [-1, 0, 1, 2, 3]:
        Ik = 1 + a1*k + a2*k**2
        Ik1 = 1 + a1*(k+1) + a2*(k+1)**2
        Ik2 = 1 + a1*(k+2) + a2*(k+2)**2
        diff2 = Ik2 - 2*Ik1 + Ik
        assert diff2 == 2*a2, f"FAIL: {diff2} != {2*a2}"
    print(f"    alpha_2={a2}: all second differences = {2*a2} ✓")

# ======================================================================
# PART 8: I(1) = 1 + alpha_1 + alpha_2 = "chromatic" evaluation
# ======================================================================
print("\n" + "=" * 70)
print("PART 8: I(1) AND THE EDGE COUNT")
print("=" * 70)

# I(1) = 1 + a1 + a2 = number of independent sets of Omega(T)
# (counting the empty set).
# I(0) = 1 always (just the empty set).
# I(1) = total number of independent sets.

I1_vals = [d['Ix'][1] for d in data]
a1_vals = [d['a1'] for d in data]
a2_vals = [d['a2'] for d in data]

print(f"\n  I(1) = 1 + alpha_1 + alpha_2 (total independent sets of Omega):")
print(f"    range: [{min(I1_vals)}, {max(I1_vals)}]")
print(f"    mean: {np.mean(I1_vals):.1f}")

print(f"\n  I(0) = 1 always (empty set)")
print(f"  I(1) = total number of independent sets")
print(f"  I(2) = H = binary-weighted count of independent sets")
print(f"  I(3) = ternary-weighted count")
print(f"  I(-1) = alternating count (Euler char of independence complex)")
print(f"  I(-2) = 1 - 2*alpha_1 + 4*alpha_2")

I_neg2_vals = [d['Ix'][-2] for d in data]
print(f"\n  I(-2) statistics at n=7:")
print(f"    range: [{min(I_neg2_vals)}, {max(I_neg2_vals)}]")
print(f"    mean: {np.mean(I_neg2_vals):.1f}")

# I(-2) = 1 - 2a1 + 4a2 = H - 4a1 = (1 + 2a1 + 4a2) - 4a1 = 1 - 2a1 + 4a2
# Alternatively: I(-2) = I(2) - 4a1 = H - 4a1
# So I(-2) > 0 iff H > 4a1 iff 1 + 2a1 + 4a2 > 4a1 iff 1 + 4a2 > 2a1.
print(f"\n  I(-2) = H - 4*alpha_1 = 1 - 2*alpha_1 + 4*alpha_2")
print(f"  I(-2) > 0 iff 1 + 4*alpha_2 > 2*alpha_1")

# Check
pos_neg2 = sum(1 for v in I_neg2_vals if v > 0)
print(f"  I(-2) > 0: {pos_neg2}/{len(I_neg2_vals)} = {pos_neg2/len(I_neg2_vals)*100:.1f}%")

print("\nDone.")
