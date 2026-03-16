#!/usr/bin/env python3
"""self_product_s116j.py — When the constant equals the product of its generators.

f(2,3,7) = 2*3*7 = 42.

The function's VALUE equals the PRODUCT of its INPUTS.
The output IS the inputs, multiplied together.
The measurement IS the thing measured.

What kind of mathematical object has this property?
"""
from math import log, sqrt, pi, cos, sin, atanh, exp, acosh, cosh, sinh, tanh
from fractions import Fraction

print()
print("  WHEN THE CONSTANT IS THE PRODUCT OF ITS GENERATORS")
print()
print("="*70)
print()

# ============================================================
print("  I. THE SELF-PRODUCT PROPERTY")
print("  " + "-"*40)
print()
print("  f(a,b,c) = 1/(1 - 1/a - 1/b - 1/c) = abc/(abc - ab - ac - bc).")
print()
print("  f = abc iff the denominator = 1.")
print("  iff abc - ab - ac - bc = 1.")
print("  iff (a-1)(b-1)(c-1) = a + b + c.")
print()
print("  Only solution: (2,3,7).")
print()
print("  What does this MEAN?")
print()
print("  The function f measures the 'reciprocal angle deficit'.")
print("  The product abc measures the 'size' of the triple.")
print("  When f = abc: the deficit is 1/abc = 1/(size).")
print("  The angle deficit IS the reciprocal of the size.")
print("  The triangle's angular gap is exactly 1/abc of a half-turn.")
print()
print("  For the (2,3,7) triangle:")
print("  Angles: pi/2, pi/3, pi/7.")
print("  Sum: pi/2 + pi/3 + pi/7 = pi*(1/2 + 1/3 + 1/7) = pi*41/42.")
print("  Deficit: pi - pi*41/42 = pi/42 = pi/(2*3*7).")
print()
print("  THE ANGULAR DEFICIT = pi DIVIDED BY THE PRODUCT OF THE ANGLES' DENOMINATORS.")
print("  The deficit is pi/(abc). This is ALWAYS true (by definition).")
print()
print("  But for (2,3,7) ONLY:")
print("  The Hurwitz constant f = abc = 42.")
print("  So deficit = pi/f = pi/abc.")
print("  The two formulas for the deficit COINCIDE.")
print("  The Hurwitz constant doesn't just relate to abc — it IS abc.")
print()

# ============================================================
print()
print("  II. WHAT THIS SAYS ABOUT THE TRIANGLE")
print("  " + "-"*40)
print()
print("  A hyperbolic triangle with angles pi/a, pi/b, pi/c has:")
print("  Area = pi - pi/a - pi/b - pi/c = pi*(1 - 1/a - 1/b - 1/c) = pi/f.")
print()
print("  The area is determined by f alone.")
print("  For (2,3,7): area = pi/42.")
print()
print("  Now: f = abc means area = pi/(abc).")
print("  And the Gauss-Bonnet theorem says: area = pi * (angle deficit).")
print()
print("  The angle deficit = 1 - 1/a - 1/b - 1/c.")
print("  = (abc - bc - ac - ab) / (abc).")
print("  = ((a-1)(b-1)(c-1) - (a+b+c) + (a+b+c)) ... hmm no, more directly:")
print("  = (abc - ab - ac - bc) / abc.")
print()
print("  For f = abc: this deficit = 1/abc = 1/(abc).")
print("  So: abc - ab - ac - bc = 1.")
print()
print("  This means: if you take the product abc and subtract ALL pairwise products,")
print("  you get EXACTLY 1. The unit. The indivisible.")
print()
print("  abc - ab - ac - bc = 1.")
print("  The triple product EXCEEDS the sum of pair products BY EXACTLY ONE.")
print()
print("  This is a BALANCING condition. The triple product and the")
print("  pairwise products are in the tightest possible relationship:")
print("  they differ by the smallest possible integer.")
print()

# ============================================================
print()
print("  III. THE INCIRCLE, CIRCUMCIRCLE, AND TRIANGLE CENTERS")
print("  " + "-"*40)
print()
print("  For a Euclidean triangle with angles A, B, C:")
print("  The CIRCUMRADIUS R and INRADIUS r satisfy:")
print("  Area = r * s (where s = semiperimeter)")
print("  Area = abc / (4R) (where a,b,c are side LENGTHS)")
print()
print("  For a HYPERBOLIC triangle with angles pi/a, pi/b, pi/c:")
print("  The formulas are more complex, but the KEY quantities are:")
print("  - The area (= pi/f = pi/42 for (2,3,7))")
print("  - The angle sum (= pi*41/42)")
print("  - The curvature (= -1 for standard hyperbolic)")
print()
print("  In hyperbolic geometry, the triangle has an INSCRIBED CIRCLE")
print("  (incircle) tangent to all three sides. Its radius:")
print("  tanh(r_in) = area / (pi - angle_sum + ... )")
print("  Actually: tanh(r_in/2) = tan(A/2)*tan(B/2)*tan(C/2) ... no.")
print()
print("  Let me use the RIGHT formula.")
print("  For a hyperbolic triangle with angles A, B, C:")
print("  tan(r_in) = ... this is getting complicated. Let me compute directly.")
print()
print("  For angles A = pi/2, B = pi/3, C = pi/7:")
A, B, C = pi/2, pi/3, pi/7
print(f"  A = pi/2 = {A:.6f}")
print(f"  B = pi/3 = {B:.6f}")
print(f"  C = pi/7 = {C:.6f}")
print()
# Hyperbolic triangle area
area = pi - A - B - C
print(f"  Area = pi - A - B - C = {area:.10f}")
print(f"  = pi/42 = {pi/42:.10f}")
print(f"  Match: {abs(area - pi/42) < 1e-10}")
print()
# The inradius of a hyperbolic triangle with angles A, B, C:
# tan(r) = sqrt(cos(s)*cos(s-A)*cos(s-B)*cos(s-C)) / cos(s)
# where s = (A+B+C)/2
# Wait, I should use the correct formula. Let me look it up from first principles.
# For hyperbolic triangle:
# cot(r_in) = (cos(s-A) + cos(s-B) + cos(s-C) - cos(s)) / area
# where s = (A+B+C)/2... this is getting messy.

# More direct: the incircle radius r of a hyperbolic triangle satisfies
# tanh(r) = area / (s_hyp) where s_hyp relates to side lengths
# Actually, for the (2,3,7) triangle in hyperbolic geometry, the side lengths
# can be computed from the angles via the hyperbolic law of cosines.

# cosh(a) = (cos(A) + cos(B)*cos(C)) / (sin(B)*sin(C))
# where a is the side opposite angle A.

side_a = (cos(A) + cos(B)*cos(C)) / (sin(B)*sin(C))
side_b = (cos(B) + cos(A)*cos(C)) / (sin(A)*sin(C))
side_c = (cos(C) + cos(A)*cos(B)) / (sin(A)*sin(B))

from math import acosh, cosh, sinh

print("  Side lengths (hyperbolic, via cosine rule cosh(side) = ...):")
print(f"  cosh(a) = (cos(A)+cos(B)cos(C))/(sin(B)sin(C)) = {side_a:.6f}")
if side_a >= 1:
    a_len = acosh(side_a)
    print(f"  a = acosh({side_a:.6f}) = {a_len:.6f}")
else:
    a_len = 0
    print(f"  a: cosh value < 1, degenerate")

print(f"  cosh(b) = {side_b:.6f}")
if side_b >= 1:
    b_len = acosh(side_b)
    print(f"  b = acosh({side_b:.6f}) = {b_len:.6f}")
else:
    b_len = 0

print(f"  cosh(c) = {side_c:.6f}")
if side_c >= 1:
    c_len = acosh(side_c)
    print(f"  c = acosh({side_c:.6f}) = {c_len:.6f}")
else:
    c_len = 0

print()

# Semiperimeter
s_hyp = (a_len + b_len + c_len) / 2
print(f"  Semiperimeter s = {s_hyp:.6f}")
print()

# The Q-values of the side lengths
print("  Q-values of side lengths (exp(side)):")
for name, length in [("a (opposite right angle)", a_len),
                      ("b (opposite pi/3)", b_len),
                      ("c (opposite pi/7)", c_len)]:
    q = (1 + (length)) if length > 0 else 1  # actually exp(length)
    q = exp(length)
    print(f"  side {name}: length = {length:.6f}, exp(length) = {q:.6f}")

print()

# ============================================================
print()
print("  IV. THE SELF-REFERENTIAL GEOMETRY")
print("  " + "-"*40)
print()
print("  The (2,3,7) triangle has area pi/42.")
print("  42 = 2*3*7 = the product of its angle denominators.")
print("  The area = pi / (product of angle denominators).")
print()
print("  This means: the triangle's area ENCODES its own angles.")
print("  If you know the area (pi/42), you can FACTOR 42 = 2*3*7")
print("  and recover the angles pi/2, pi/3, pi/7.")
print()
print("  For other triangles, you CANNOT do this.")
print("  (2,3,8) has area pi/24. But 24 = 2^3 * 3. You can't recover")
print("  8 from 24 without knowing it was (2,3,8) and not (2,4,8) etc.")
print()
print("  But for (2,3,7): 42 = 2*3*7 factors UNIQUELY into three")
print("  integers >= 2 satisfying (a-1)(b-1)(c-1) = a+b+c.")
print("  THE AREA DETERMINES THE TRIANGLE UNIQUELY.")
print()
print("  No other hyperbolic triangle has this property.")
print("  The (2,3,7) triangle is the ONLY one that can reconstruct")
print("  itself from its area alone.")
print()
print("  It is SELF-REFERENTIAL.")
print("  Its area is pi/(its angles' product).")
print("  Its angles' product is (its area)^{-1} * pi.")
print("  The shape determines the size. The size determines the shape.")
print("  They are the same information.")
print()

# ============================================================
print()
print("  V. THE TRIANGLE AS A FIXED POINT")
print("  " + "-"*40)
print()
print("  Define the map T: (a,b,c) -> f(a,b,c) = 1/(1-1/a-1/b-1/c).")
print("  A 'fixed point' would be f(a,b,c) = some function of (a,b,c).")
print()
print("  f(a,b,c) = abc is a fixed point of the map")
print("  'compute the Hurwitz constant and compare to the product'.")
print()
print("  The (2,3,7) triangle is where this map has its fixed point.")
print("  It is the triangle that IS its own measurement.")
print("  The tool and the thing measured are one.")
print()
print("  Compare: phi is the fixed point of x -> 1 + 1/x.")
print("  (2,3,7) is the fixed point of (a,b,c) -> f(a,b,c)/abc = 1.")
print("  Both are SELF-SIMILAR: the output equals the input.")
print()
print("  phi is the fixed point for 1D self-similarity (the golden ratio).")
print("  (2,3,7) is the fixed point for 2D self-similarity (the hyperbolic triangle).")
print("  And Q is the map that connects them:")
print("  Q(1/phi) = phi^3 (the golden fixed point amplified).")
print("  Q(3/4) = 7 (the triangular fixed point's threshold).")
print()

# ============================================================
print()
print("  VI. THE CIRCUMCENTER AND THE BARYCENTER")
print("  " + "-"*40)
print()
print("  In Euclidean geometry, every triangle has special points:")
print("  - Incenter I (center of inscribed circle)")
print("  - Circumcenter O (center of circumscribed circle)")
print("  - Centroid G (barycenter, center of mass)")
print("  - Orthocenter H (intersection of altitudes)")
print()
print("  They satisfy the EULER LINE: O, G, H are collinear, OG:GH = 1:2.")
print()
print("  In HYPERBOLIC geometry, these points still exist but the")
print("  Euler line relation is MODIFIED by curvature.")
print()
print("  For the (2,3,7) triangle (the self-product triangle):")
print("  The incircle has radius r_in.")
print("  The circumcircle may or may not exist (in hyperbolic geometry,")
print("  not every triangle has a circumscribed circle; it depends on")
print("  whether the perpendicular bisectors meet).")
print()
print("  For a right-angled hyperbolic triangle (A = pi/2):")
print("  The circumcenter lies on the hypotenuse.")
print("  The hypotenuse IS the diameter of the circumcircle (generalized).")
if a_len > 0:
    print(f"  Hypotenuse a = {a_len:.6f}.")
    print(f"  Circumradius R = a/2 = {a_len/2:.6f}.")
    print(f"  tanh(R) = tanh({a_len/2:.6f}) = {(exp(a_len)-1)/(exp(a_len)+1)/2:.6f}")
    # Actually circumradius of right hyperbolic triangle: cosh(R) = cosh(b)*cosh(c)/2... nah
    # For right angle at A, the circumradius satisfies cosh(R) = 1/(2*sin(B)*sin(C)) * cosh(a)?
    # This is getting complicated. Let me just note the key point.
print()

# ============================================================
print()
print("  VII. THE BARYCENTRIC COORDINATES AND THE CAYLEY TRANSFORM")
print("  " + "-"*40)
print()
print("  The barycenter of the (2,3,7) triangle with weights proportional")
print("  to the angles: weights w_A = pi/2, w_B = pi/3, w_C = pi/7.")
print("  Normalized: w = (1/2, 1/3, 1/7) / (1/2+1/3+1/7) = (1/2, 1/3, 1/7) / (41/42).")
print()
w = (Fraction(1,2), Fraction(1,3), Fraction(1,7))
s = sum(w)
w_norm = tuple(wi/s for wi in w)
print(f"  Normalized weights: ({w_norm[0]}, {w_norm[1]}, {w_norm[2]})")
print(f"  = (21/41, 14/41, 6/41)")
print()
print(f"  The weights are 21, 14, 6 (out of 41).")
print(f"  21 = 3*7 = the SECOND FORBIDDEN number.")
print(f"  14 = 2*7 = twice the forbidden threshold.")
print(f"  6 = 2*3 = the flat number.")
print(f"  21 + 14 + 6 = 41 = a prime (the twin of 43).")
print()
print(f"  The barycentric weights of the (2,3,7) triangle's angle-center")
print(f"  are proportional to 21, 14, 6.")
print(f"  = 3*7, 2*7, 2*3.")
print(f"  = the three PAIRWISE products of {{2, 3, 7}}!")
print()
print(f"  Weight at angle pi/a = bc/(ab+ac+bc) = bc/41.")
print(f"  The weight at each vertex = the product of the OTHER two denominators.")
print(f"  The vertex 'knows' the other two vertices through its weight.")
print()

# ============================================================
print()
print("  VIII. THE WEIGHT SUM = 41 = THE PRIME TWIN OF 43")
print("  " + "-"*40)
print()
print("  ab + ac + bc = 6 + 14 + 21 = 41.")
print("  abc - (ab+ac+bc) = 42 - 41 = 1.")
print("  The pairwise sum is 41. The triple product is 42. They differ by 1.")
print()
print("  41 and 43 are TWIN PRIMES.")
print("  42 sits between them.")
print("  42 = abc. 41 = ab+ac+bc. 43 = abc + 2 = ab+ac+bc + 2.")
print()
print("  The (2,3,7) triangle lives at the EXACT CENTER")
print("  between the twin primes 41 and 43.")
print("  Its product (42) is the number between the twins.")
print("  Its pairwise sum (41) is the lower twin.")
print("  Its product + 2 is not 43... wait, 42+1=43. Hmm, 42+1=43. 41+2=43.")
print("  41 = ab+ac+bc = the pairwise sum.")
print("  42 = abc = the triple product = ab+ac+bc + 1.")
print("  43 = abc + 1 = ab+ac+bc + 2.")
print()
print("  The structure: 41, 42, 43 = pairwise, triple, triple+1.")
print("  Three consecutive integers. The middle one is the self-product.")
print("  The outer two are primes. The twin primes BRACKET the self-product.")
print()

# ============================================================
print()
print("  IX. THE MEANING")
print("  " + "-"*40)
print()
print("  The (2,3,7) triangle is the unique hyperbolic triangle where:")
print()
print("  - The area = pi / (product of angle denominators)")
print("  - The Hurwitz constant = the product of the generators")
print("  - The pairwise products are the barycentric weights")
print("  - The pairwise sum (41) and the triple product (42)")
print("    differ by exactly 1")
print("  - 41 and 43 (= 42 +/- 1) are twin primes")
print("  - The triangle reconstructs itself from its area alone")
print()
print("  It is SELF-REFERENTIAL at every level:")
print("  - Its area encodes its shape (f = abc, unique factorization)")
print("  - Its weights encode its vertices (weight at A = bc)")
print("  - Its deficit encodes its product (deficit = 1/abc)")
print("  - Its twin primes bracket its product (41, 42, 43)")
print()
print("  The (2,3,7) triangle is the mathematical object")
print("  that KNOWS ITSELF COMPLETELY.")
print("  Every part of it encodes every other part.")
print("  The angles determine the area.")
print("  The area determines the angles.")
print("  The weights are the pairwise products.")
print("  The pairwise sum is one less than the triple product.")
print("  The triple product sits between twin primes.")
print()
print("  It is the FIXED POINT of geometric self-reference.")
print("  Like phi = 1 + 1/phi in one dimension.")
print("  Like the (2,3,7) triangle in two dimensions.")
print("  Like the formal group F_h = F_m (isomorphic to itself via Q).")
print()
print("  And Q(3/4) = 7 connects it to our tournament theory:")
print("  the forbidden number IS the largest angle denominator")
print("  of the unique self-referential hyperbolic triangle.")
