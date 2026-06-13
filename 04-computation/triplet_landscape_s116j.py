#!/usr/bin/env python3
"""triplet_landscape_s116j.py — The function f(a,b,c) = 1/(1 - 1/a - 1/b - 1/c).

When f > 0 (hyperbolic): f = the reciprocal of the angle deficit = Hurwitz-like constant.
When f < 0 (spherical): |f| relates to the order of the finite symmetry group.
When f = infinity (flat): the Euclidean tiling.

Vary a, b, c. See everything.
"""
from fractions import Fraction

print()
print("  THE TRIPLET FUNCTION f(a,b,c) = 1/(1 - 1/a - 1/b - 1/c)")
print()
print("="*70)
print()

def f_triplet(a, b, c):
    """Returns (value as Fraction, type string)."""
    s = 1 - Fraction(1,a) - Fraction(1,b) - Fraction(1,c)
    if s > 0:
        return 1/s, "HYPER"
    elif s == 0:
        return None, "FLAT"
    else:
        return 1/s, "SPHERE"

# ============================================================
print("  I. THE FULL TABLE FOR SMALL a <= b <= c")
print("  " + "-"*40)
print()
print(f"  {'(a,b,c)':<12s} {'1/a+1/b+1/c':>12s} {'deficit':>10s} {'f':>8s}  {'type':>8s}  {'product':>8s}")
print("  " + "-"*70)

results = []
for a in range(2, 10):
    for b in range(a, 15):
        for c in range(b, 20):
            s = Fraction(1,a) + Fraction(1,b) + Fraction(1,c)
            deficit = 1 - s
            prod = a*b*c
            if deficit > 0:
                f_val = 1/deficit
                typ = "HYPER"
            elif deficit == 0:
                f_val = None
                typ = "FLAT"
            else:
                f_val = 1/deficit
                typ = "SPHERE"
            results.append((a,b,c,s,deficit,f_val,typ,prod))

# Sort by f value (descending for hyperbolic)
hyper = [(a,b,c,s,d,f,t,p) for a,b,c,s,d,f,t,p in results if t=="HYPER"]
flat = [(a,b,c,s,d,f,t,p) for a,b,c,s,d,f,t,p in results if t=="FLAT"]
sphere = [(a,b,c,s,d,f,t,p) for a,b,c,s,d,f,t,p in results if t=="SPHERE"]

print("  SPHERICAL (finite groups):")
for a,b,c,s,d,f,t,p in sorted(sphere, key=lambda x: -float(x[5])):
    print(f"  ({a},{b},{c})    {str(s):>12s}  {str(d):>10s}  {str(f):>8s}  {t:>8s}  {p:>8d}")

print()
print("  FLAT (Euclidean tilings):")
for a,b,c,s,d,f,t,p in flat:
    print(f"  ({a},{b},{c})    {str(s):>12s}  {'0':>10s}  {'inf':>8s}  {t:>8s}  {p:>8d}")

print()
print("  HYPERBOLIC (top 20 by f value):")
hyper_sorted = sorted(hyper, key=lambda x: -float(x[5]))
for a,b,c,s,d,f,t,p in hyper_sorted[:20]:
    print(f"  ({a},{b},{c})    {str(s):>12s}  {str(d):>10s}  {str(f):>8s}  {t:>8s}  {p:>8d}")

print()

# ============================================================
print()
print("  II. THE SPHERICAL TRIPLES (f < 0)")
print("  " + "-"*40)
print()
print("  These give FINITE groups. The Platonic solids and their relatives.")
print()
for a,b,c,s,d,f,t,p in sorted(sphere, key=lambda x: x[7]):
    # The group order is related to the triple
    # For (2,2,n): dihedral D_n, order 2n
    # For (2,3,3): tetrahedral A_4, order 12
    # For (2,3,4): octahedral S_4, order 24
    # For (2,3,5): icosahedral A_5, order 60
    group = ""
    if a == 2 and b == 2:
        group = f"Dihedral D_{c}, order {2*c}"
    elif (a,b,c) == (2,3,3):
        group = "Tetrahedral A_4, order 12"
    elif (a,b,c) == (2,3,4):
        group = "Octahedral S_4, order 24"
    elif (a,b,c) == (2,3,5):
        group = "Icosahedral A_5, order 60"
    # f = 1/(1-s) = 1/deficit. deficit < 0. so f < 0.
    # |f| = 1/(s-1) = 1/(1/a+1/b+1/c - 1)
    abs_f = -float(f)
    print(f"  ({a},{b},{c}): product = {p:4d}, |f| = {abs_f:8.2f}, {group}")

print()

# ============================================================
print()
print("  III. THE FLAT TRIPLES (f = infinity)")
print("  " + "-"*40)
print()
for a,b,c,s,d,f,t,p in flat:
    print(f"  ({a},{b},{c}): product = {p}. Euclidean tiling by (pi/{a}, pi/{b}, pi/{c}) triangles.")

print()
print("  These are the THREE Euclidean triangle tilings:")
print("  (2,3,6): 30-60-90 triangle. Product 36.")
print("  (2,4,4): 45-45-90 triangle. Product 32.")
print("  (3,3,3): 60-60-60 equilateral. Product 27.")
print()
print("  Products: 27, 32, 36.")
print(f"  27 = 3^3. 32 = 2^5. 36 = 6^2 = (2*3)^2.")
print()

# ============================================================
print()
print("  IV. VARYING ONE ELEMENT: FIX a=2, b=3, VARY c")
print("  " + "-"*40)
print()
print("  The (2,3,c) family:")
print(f"  {'c':>4s}  {'deficit':>12s}  {'f':>10s}  {'type':>8s}  {'product':>8s}")
print("  " + "-"*50)
for c in range(2, 25):
    s = Fraction(1,2) + Fraction(1,3) + Fraction(1,c)
    d = 1 - s
    p = 2*3*c
    if d > 0:
        fv = float(1/d)
        typ = "HYPER"
    elif d == 0:
        fv = float('inf')
        typ = "FLAT"
    else:
        fv = float(1/d)
        typ = "SPHERE"
    mark = ""
    if c == 5: mark = " <-- icosahedral (LAST spherical)"
    elif c == 6: mark = " <-- FLAT"
    elif c == 7: mark = " <-- Hurwitz (f=42, MAXIMUM hyper)"
    elif c == 8: mark = " <-- f=24"
    elif c == 12: mark = " <-- f=12"
    print(f"  {c:4d}  {str(d):>12s}  {fv:>10.2f}  {typ:>8s}  {p:>8d}{mark}")

print()
print("  As c increases past 7, f DECREASES: 42, 24, 18, 15, ...")
print("  f(2,3,c) = 1/(1-1/2-1/3-1/c) = 1/(1/6 - 1/c) = 6c/(c-6).")
print()
print("  f = 6c/(c-6):")
for c in range(7, 25):
    f_formula = 6*c/(c-6)
    print(f"    c={c:2d}: f = 6*{c}/({c}-6) = {f_formula:.1f}")

print()
print("  The MAXIMUM is at c=7: f=42.")
print("  This is because 6c/(c-6) is DECREASING for c > 6.")
print("  The FIRST integer past the flat boundary (c=7) gives the HIGHEST f.")
print()
print("  f(2,3,7) = 42 is the maximum because 7 is the SMALLEST integer")
print("  that makes 1/2 + 1/3 + 1/c < 1.")
print("  The smallest violation of the flat condition gives the largest constant.")
print("  THE FIRST STEP INTO HYPERBOLICITY IS THE BIGGEST.")
print()

# ============================================================
print()
print("  V. VARYING TWO ELEMENTS: FIX a=2, VARY b AND c")
print("  " + "-"*40)
print()
print("  Top f values for (2,b,c) triples:")
print(f"  {'(2,b,c)':<12s} {'f':>10s}  {'product':>8s}")
print("  " + "-"*35)
results_2bc = []
for b in range(3, 20):
    for c in range(b, 50):
        d = 1 - Fraction(1,2) - Fraction(1,b) - Fraction(1,c)
        if d > 0:
            fv = float(1/d)
            results_2bc.append((b, c, fv, 2*b*c))

for b,c,fv,p in sorted(results_2bc, key=lambda x: -x[2])[:15]:
    mark = ""
    if (b,c) == (3,7): mark = " <-- HURWITZ"
    elif (b,c) == (3,8): mark = ""
    elif (b,c) == (4,5): mark = " <-- note: 4*5=20"
    print(f"  (2,{b},{c})      {fv:>10.1f}  {p:>8d}{mark}")

print()
print("  The hierarchy: (2,3,7)=42 > (2,3,8)=24 > (2,3,9)=18 > (2,4,5)=20 > ...")
print("  (2,3,7) is the GLOBAL MAXIMUM among all (2,b,c) hyperbolic triples.")
print()

# ============================================================
print()
print("  VI. WHAT a=2 MEANS")
print("  " + "-"*40)
print()
print("  a=2 means the triangle has a RIGHT ANGLE (pi/2).")
print("  The (2,b,c) triangles are all RIGHT TRIANGLES.")
print()
print("  Without a=2:")
print("  (3,3,4): d = 1-1/3-1/3-1/4 = 1/12. f = 12. Product = 36.")
print("  (3,3,3): d = 0. FLAT. Product = 27.")
print("  (3,4,4): d = 1-1/3-1/4-1/4 = 1/6. f = 6. Product = 48.")
print()
print("  The TOP non-right-angle triple: (3,3,4) with f = 12.")
print("  Compare (2,3,7) with f = 42. The right angle TRIPLES the constant.")
print()
print("  WHY is a=2 so powerful?")
print("  Because 1/2 is the LARGEST unit fraction.")
print("  Starting with the biggest 'bite' (1/2) leaves the most room")
print("  for the other two angles to be small (closer to the flat boundary).")
print("  Small b and c = large 1/b + 1/c = close to the flat line 1-1/a.")
print("  Closeness to flat = large f.")
print()
print("  The right angle IS the doubler (a=2).")
print("  The doubler makes the constant largest because it takes")
print("  the biggest possible first step, leaving the deficit smallest.")
print()

# ============================================================
print()
print("  VII. THE PRODUCT abc AND THE CONSTANT f")
print("  " + "-"*40)
print()
print("  Is f always equal to abc? No. But close?")
print()
print(f"  {'(a,b,c)':<12s} {'f':>8s}  {'abc':>8s}  {'f/abc':>8s}")
print("  " + "-"*40)
for a,b,c,s,d,fv,t,p in hyper_sorted[:15]:
    ratio = float(fv)/p
    print(f"  ({a},{b},{c})      {float(fv):>8.1f}  {p:>8d}  {ratio:>8.4f}")

print()
print("  f/abc -> 0 as the triple grows.")
print("  For (2,3,7): f/abc = 42/42 = 1. EXACTLY!")
print("  For (2,3,8): f/abc = 24/48 = 0.5.")
print("  For (2,3,9): f/abc = 18/54 = 0.333...")
print()
print("  (2,3,7) is the UNIQUE triple where f = abc.")
print("  The Hurwitz constant EQUALS the product of the triple.")
print("  For every other triple, f < abc.")
print()

# Verify
d_237 = 1 - Fraction(1,2) - Fraction(1,3) - Fraction(1,7)
f_237 = 1/d_237
print(f"  Verify: f(2,3,7) = {f_237} = {2*3*7} = 2*3*7. YES.")
print()
print("  WHY? f = 1/(1-1/a-1/b-1/c) = abc/(abc - bc - ac - ab).")
print("  f = abc iff abc - bc - ac - ab = 1.")
print("  iff bc + ac + ab = abc - 1.")
print()
print("  For (2,3,7): bc+ac+ab = 21+14+6 = 41. abc-1 = 41. YES!")
print()
print("  When does bc + ac + ab = abc - 1?")
print("  Rearranging: abc - ab - ac - bc = 1.")
print("  (a-1)(b-1)(c-1) = abc - ab - ac - bc + a + b + c - 1.")
print("  So abc - ab - ac - bc = (a-1)(b-1)(c-1) - a - b - c + 1.")
print("  Need: (a-1)(b-1)(c-1) - a - b - c + 1 = 1.")
print("  i.e., (a-1)(b-1)(c-1) = a + b + c.")
print()
print("  For (2,3,7): (1)(2)(6) = 12. a+b+c = 12. 12 = 12. YES!")
print()
print("  THE EQUATION: (a-1)(b-1)(c-1) = a + b + c.")
print()
print("  Solutions with 2 <= a <= b <= c:")

solutions = []
for a in range(2, 20):
    for b in range(a, 50):
        for c in range(b, 200):
            lhs = (a-1)*(b-1)*(c-1)
            rhs = a+b+c
            if lhs == rhs:
                solutions.append((a,b,c))
            elif lhs > rhs + 100:
                break

for a,b,c in solutions:
    print(f"    ({a},{b},{c}): (a-1)(b-1)(c-1) = {(a-1)*(b-1)*(c-1)}, a+b+c = {a+b+c}. f = abc = {a*b*c}.")

print()
if len(solutions) == 1:
    print(f"  THERE IS ONLY ONE SOLUTION: (2, 3, 7).")
    print()
    print(f"  (2,3,7) is the UNIQUE triple where f(a,b,c) = a*b*c.")
    print(f"  42 is the UNIQUE number that is both the Hurwitz constant")
    print(f"  AND the product of its generating triple.")
    print(f"  This is not numerology. This is a theorem.")
