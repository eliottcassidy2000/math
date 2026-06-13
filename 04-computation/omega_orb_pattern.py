#!/usr/bin/env python3
"""
Pattern analysis for Omega_orb sequences.

P_7 (m=3): [1, 1, 2, 3, 3, 2, 1]
P_11 (m=5): [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6]
P_19 (m=9): [1, 1, 8, 60, 417, 2648, 15140, 76474, ...]

Full Omega = m * Omega_orb for d >= 1.
Omega_0 = 1 always.
Omega_1 = m always.

Questions:
1. What is Omega_orb_2 = (m-1)? No: 2, 4, 8 for m=3,5,9 → m-1.
   m-1 = 2, 4, 8. Yes! Omega_orb_2 = m - 1.
   Omega_2 = m(m-1). Already known.

2. Omega_orb_3: 3, 14, 60 for m=3,5,9.
   Check polynomial in m:
   m=3: 3, m=5: 14, m=9: 60
   Try quadratic: a*m^2 + b*m + c
   9a + 3b + c = 3
   25a + 5b + c = 14
   81a + 9b + c = 60
   16a + 2b = 11
   56a + 4b = 46 → 28a + 2b = 23
   12a = 12 → a = 1
   2b = 11 - 16 = -5 → b = -5/2
   c = 3 - 9 + 15/2 = -6 + 7.5 = 1.5

   Not integer coefficients. Try ratio:
   m=3: 3 = 3
   m=5: 14 = 14
   m=9: 60 = 60

   3/(1*2*3) = 1/2
   14/(3*4*5) = 14/60 = 7/30
   60/(7*8*9) = 60/504 = 5/42

   Try C(m,2): m=3: 3, m=5: 10, m=9: 36. No.
   Try C(m,3)/something: m=3: 1, m=5: 10, m=9: 84. No.

   Actually from formula: |A_3| = m(m-1)(2m+1)/2.
   Omega_3 = |A_3| - rank(C_3).
   |A_3|/m = (m-1)(2m+1)/2.

   For m=3: 2*7/2 = 7. But orbits = 7. Omega_orb = 3 ≠ 7.
   So there IS a constraint at d=3.

   Orbits of A_3: |A_3|/m = 21/3=7, 110/5=22, 684/9=76.
   Omega_orb_3: 3, 14, 60.
   Junk orbits (rank of C_orb_3): 7-3=4, 22-14=8, 76-60=16.

   Junk count: 4, 8, 16. Pattern: 4(m-1)? m=3: 8≠4, m=5: 16≠8, m=9: 32≠16. No.
   Powers of 2: 4=2^2, 8=2^3, 16=2^4. But m=3,5,9 not powers.
   Actually: junk = 2^{log2(m)+1}? m=3: 2^2.58=6. No.
   Just 4, 8, 16 = 4·1, 4·2, 4·4. Ratios: 2, 2. So doubling.
   m=3→5: ×2. m=5→9: ×2. So junk_3 = 4 · 2^{(m-3)/2}?
   m=3: 4·1=4✓, m=5: 4·2=8✓, m=9: 4·4=16✓.
   i.e., junk_3 = 2^{(m-3)/2 + 2} = 2^{(m+1)/2}.
   m=3: 2^2=4✓, m=5: 2^3=8✓, m=9: 2^5=32≠16✗.
   Hmm.

   Let me just look at the data.

opus-2026-03-13-S71b
"""

from fractions import Fraction

# Data
m_vals = [3, 5, 9]
omega_orb = {
    3: [1, 1, 2, 3, 3, 2, 1],
    5: [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6],
    9: [1, 1, 8, 60, 417, 2648, 15140, 76474],
}

# Full A_d orbit counts
A_orb = {
    3: [1, 1, 3, 7, 13, 15, 9],
    5: [1, 1, 5, 22, 86, 286, 794, 1747, 2879, 3149, 1729],
    9: [1, 1, 9, 76, 604, 4476, 30804, 195546],
}

print("=== OMEGA ORBIT ANALYSIS ===\n")

for m in m_vals:
    O = omega_orb[m]
    A = A_orb[m]
    p = 2*m + 1
    print(f"P_{p} (m={m}):")
    print(f"  A_orb:     {A}")
    print(f"  Omega_orb: {O}")
    junk = [A[d] - O[d] for d in range(len(O))]
    print(f"  Junk rank: {junk}")
    print()

# Compare Omega_orb values at same degree
print("=== CROSS-PRIME COMPARISON ===\n")
for d in range(8):
    vals = {}
    for m in m_vals:
        if d < len(omega_orb[m]):
            vals[m] = omega_orb[m][d]
    if len(vals) >= 2:
        print(f"d={d}: {vals}")
        # Try to fit polynomial in m
        ms = sorted(vals.keys())
        if len(ms) == 3:
            # Lagrange interpolation
            y0, y1, y2 = vals[ms[0]], vals[ms[1]], vals[ms[2]]
            m0, m1, m2 = ms[0], ms[1], ms[2]
            # Check polynomial degree
            print(f"  Ratios: {y1/y0:.4f}, {y2/y1:.4f}")

            # First differences
            if m1 - m0 == m2 - m1:
                d1 = y1 - y0
                d2 = y2 - y1
                print(f"  First diff: {d1}, {d2}")
                print(f"  Second diff: {d2 - d1}")

# Known formulas
print("\n=== KNOWN FORMULAS ===\n")
for m in m_vals:
    p = 2*m + 1
    A = A_orb[m]
    O = omega_orb[m]

    # A_1/m = 1, A_2/m = m (since |A_2| = m^2)
    print(f"P_{p} (m={m}):")
    print(f"  |A_1|/m = {A[1]} = 1 ✓")
    print(f"  |A_2|/m = {A[2]} = {m} {'✓' if A[2] == m else '✗'}")
    if len(A) > 3:
        # |A_3| = m(m-1)(2m+1)/2. A_3/m = (m-1)(2m+1)/2.
        expected_A3_orb = (m-1)*(2*m+1)//2
        print(f"  |A_3|/m = {A[3]}, expected (m-1)(2m+1)/2 = {expected_A3_orb} "
              f"{'✓' if A[3] == expected_A3_orb else '✗'}")

    # Omega_1 = m, Omega_1/m = 1 ✓
    # Omega_2 = m(m-1), Omega_2/m = m-1
    expected_O2 = m - 1
    print(f"  Omega_2/m = {O[2]}, expected m-1 = {expected_O2} "
          f"{'✓' if O[2] == expected_O2 else '✗'}")

    # Omega_3 = m(2m^2-3m+1)/2? Let me check.
    # For m=3: Omega_3 = 9. So Omega_3/m = 3. 2*9-9+1 = 10/2 = 5 ≠ 3.
    # Actually from our data: Omega_3 = [9, 70, 540] for m=3,5,9.
    # Omega_3/m = [3, 14, 60].
    # Known formula: Omega_3 = ? Let me check against |A_3| - rank(C_3).
    print(f"  Omega_3/m = {O[3]}")

    print()

# Look for the pattern in Omega_orb_3
print("=== OMEGA_ORB_3 PATTERN ===\n")
vals = [omega_orb[m][3] for m in m_vals]
ms = m_vals
print(f"m = {ms}")
print(f"Ω_3^orb = {vals}")

# Try: Ω_3^orb = (m-1)(2m-3)/2 ?
# m=3: 2*3/2=3 ✓, m=5: 4*7/2=14 ✓, m=9: 8*15/2=60 ✓
for m, v in zip(ms, vals):
    formula = (m-1)*(2*m-3)//2
    print(f"  m={m}: (m-1)(2m-3)/2 = {formula} {'✓' if formula == v else '✗'}")

print(f"\nΩ_3^orb = (m-1)(2m-3)/2")
print(f"Ω_3 = m*(m-1)(2m-3)/2")
print(f"Compare: |A_3| = m*(m-1)(2m+1)/2")
print(f"Junk_3 = |A_3| - Ω_3 = m*(m-1)*(2m+1-2m+3)/2 = m*(m-1)*4/2 = 2m(m-1)")
for m in ms:
    junk = A_orb[m][3] - omega_orb[m][3]
    expected = 2 * (m-1)  # junk orbits, not full junk
    print(f"  m={m}: junk_orb = {junk}, 2(m-1) = {expected} "
          f"{'✓' if junk == expected else '✗'}")

# Omega_orb_4
print(f"\n=== OMEGA_ORB_4 PATTERN ===\n")
vals4 = [omega_orb[m][4] for m in m_vals if len(omega_orb[m]) > 4]
ms4 = [m for m in m_vals if len(omega_orb[m]) > 4]
print(f"m = {ms4}")
print(f"Ω_4^orb = {vals4}")

# A_4 = m(m-1)(2m^2-m-2)/2. A_4/m = (m-1)(2m^2-m-2)/2.
for m, v in zip(ms4, vals4):
    A4_orb = (m-1)*(2*m**2-m-2)//2
    print(f"  m={m}: |A_4|/m = {A4_orb}, Ω_4^orb = {v}")
    print(f"    junk_4 = {A4_orb - v}")

# Try polynomial: Ω_4^orb in m.
# m=3: 3, m=5: 41, m=9: 417
# 41/3 ≈ 13.67. 417/41 ≈ 10.17.
# Try cubic: a*m^3 + b*m^2 + c*m + d
# 27a + 9b + 3c + d = 3
# 125a + 25b + 5c + d = 41
# 729a + 81b + 9c + d = 417
# Three equations, four unknowns. Try d=0:
# 27a + 9b + 3c = 3
# 125a + 25b + 5c = 41
# 729a + 81b + 9c = 417
# From first: 9a + 3b + c = 1
# From second: 25a + 5b + c = 41/5. Not integer. Try different.
# Actually we need Fraction:
from fractions import Fraction
ms_f = [Fraction(m) for m in ms4]
vs_f = [Fraction(v) for v in vals4]

if len(ms_f) == 3:
    m0, m1, m2 = ms_f
    y0, y1, y2 = vs_f

    # Lagrange
    L0 = y0 / ((m0-m1)*(m0-m2))
    L1 = y1 / ((m1-m0)*(m1-m2))
    L2 = y2 / ((m2-m0)*(m2-m1))

    # Quadratic: a*m^2 + b*m + c
    a = L0 + L1 + L2
    b = -(L0*(m1+m2) + L1*(m0+m2) + L2*(m0+m1))
    c = L0*m1*m2 + L1*m0*m2 + L2*m0*m1

    print(f"\n  Lagrange quadratic fit:")
    print(f"    a = {a} = {float(a):.6f}")
    print(f"    b = {b} = {float(b):.6f}")
    print(f"    c = {c} = {float(c):.6f}")

    # Verify
    for m, v in zip(ms4, vals4):
        pred = a*m**2 + b*m + c
        print(f"    m={m}: pred={pred} = {float(pred):.1f}, actual={v}")

# Try (m-1)(2m^2 - am - b)/c patterns
print(f"\n  Searching (m-1)*f(m) formulas:")
for m, v in zip(ms4, vals4):
    ratio = Fraction(v, m-1)
    print(f"    m={m}: Ω_4^orb/(m-1) = {ratio} = {float(ratio):.4f}")

# Ω_4^orb/(m-1): 3/2=1.5, 41/4=10.25, 417/8=52.125
# Not clean. Try Ω_4^orb/(m-1)(m-2):
print(f"\n  Ω_4^orb/((m-1)(m-2)):")
for m, v in zip(ms4, vals4):
    ratio = Fraction(v, (m-1)*(m-2))
    print(f"    m={m}: {ratio} = {float(ratio):.6f}")
# 3/2=1.5, 41/12=3.417, 417/56=7.446

# Try dividing by (m-1):
# m=3: 3/2, m=5: 41/4, m=9: 417/8
# These have denominator m-1. Numerators: 3, 41, 417.
# Try the numerator as polynomial in m:
# m=3: 3, m=5: 41, m=9: 417
# Second differences: 41-3=38, 417-41=376. Not constant.
# Third diff would need another point.

# Just report the Omega values
print(f"\n=== SUMMARY: P_19 Omega values ===\n")
print(f"P_19 (m=9):")
for d, v in enumerate(omega_orb[9]):
    print(f"  Ω_{d} = {v * 9 if d > 0 else 1}")
print(f"  (values for d=8..18 still unknown)")

# Let me try to find Omega_orb_5 pattern
print(f"\n=== OMEGA_ORB_5 ===")
vals5 = [omega_orb[m][5] for m in m_vals if len(omega_orb[m]) > 5]
ms5 = [m for m in m_vals if len(omega_orb[m]) > 5]
print(f"m = {ms5}")
print(f"Ω_5^orb = {vals5}")
# m=5: 92, m=9: 2648
if len(ms5) == 2:
    print(f"  Only 2 data points, can't fit polynomial.")
    print(f"  Ratio: {vals5[1]/vals5[0]:.4f}")
