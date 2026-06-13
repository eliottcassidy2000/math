#!/usr/bin/env python3
"""dimension_ladder_s116j.py — From d to d+1.

The self-reference hierarchy:
  phi:     1D fixed point. x^2 = x + 1.     EQUATION.
  (2,3,7): 2D fixed point. f = abc.          IDENTITY.
  6174:    iterative fixed point. K(n) = n.   RECURSION.
  i:       lateral fixed point. Q(x) = x.     ROTATION.

How does each level GENERATE the next?
What is the formal group's dimension number?
Can we find the operation d -> d+1?
"""
from math import log, sqrt, pi, exp, atanh
from fractions import Fraction

phi = (1+sqrt(5))/2

print()
print("  THE DIMENSION LADDER: FROM d TO d+1")
print()
print("="*70)
print()

# ============================================================
print("  I. THE LEVELS AND THEIR EXPONENT SIGNATURES")
print("  " + "-"*40)
print()
print("  Level    Fixed point    Factorization          Exponents")
print("  " + "-"*60)
print("  d=0      1              1                      ()")
print("  d=1      42             2^1 * 3^1 * 7^1        (1, 1, 1)")
print("  d=2      6174           2^1 * 3^2 * 7^3        (1, 2, 3)")
print("  d=3      ???            2^1 * 3^? * 7^?        (1, ?, ?)")
print()
print("  The exponents of {2, 3, 7} at each level:")
print("  d=0: no exponents (the empty product = 1)")
print("  d=1: (1, 1, 1). Flat. All primes contribute equally.")
print("  d=2: (1, 2, 3). Staircase. The k-th prime gets exponent k.")
print()
print("  Pattern: at level d, the exponent of the k-th prime (k=1,2,3)")
print("  might be a polynomial of degree d-1 in k:")
print("  d=1: exponents = 1, 1, 1 = C(k-1, 0) = constant (degree 0)")
print("  d=2: exponents = 1, 2, 3 = k (degree 1)")
print("  d=3: exponents = 1, 3, 6 = C(k,2) = k(k+1)/2 = triangular? (degree 2)")
print()
for d in range(4):
    exps = []
    for k in range(1, 4):
        if d == 0:
            exps.append(0)
        elif d == 1:
            exps.append(1)
        elif d == 2:
            exps.append(k)
        elif d == 3:
            exps.append(k*(k+1)//2)  # triangular
    if d == 0:
        val = 1
    else:
        val = (2**exps[0]) * (3**exps[1]) * (7**exps[2])
    print(f"  d={d}: exponents = {tuple(exps)}, value = 2^{exps[0]}*3^{exps[1]}*7^{exps[2]} = {val}")

print()
print("  d=3 prediction: exponents (1, 3, 6).")
print(f"  2^1 * 3^3 * 7^6 = {2 * 27 * 7**6} = {2*27*117649}.")
print(f"  = {2*27*117649}.")
print(f"  = 6353046. A large number. Is it a fixed point of something?")
print()

# ============================================================
print()
print("  II. WHAT CHANGES BETWEEN LEVELS")
print("  " + "-"*40)
print()
print("  d=0 -> d=1: FROM nothing TO the flat product.")
print("    The operation: MULTIPLY the three primes together.")
print("    What it gains: EXISTENCE. The number 42 appears from {2,3,7}.")
print("    What it IS: forming the product. Dimension 0 -> 1.")
print()
print("  d=1 -> d=2: FROM flat TO staircase.")
print("    The operation: give each prime an exponent equal to its INDEX.")
print("    42 = 2^1*3^1*7^1 -> 6174 = 2^1*3^2*7^3.")
print("    What it gains: ORDERING. The primes are no longer symmetric.")
print("    What it IS: introducing a DIRECTION. Dimension 1 -> 2.")
print()
print("  d=2 -> d=3: FROM staircase TO ???.")
print("    The operation: replace linear exponents with triangular?")
print("    What it would gain: CURVATURE. The exponents are no longer linear.")
print("    What it IS: introducing BENDING. Dimension 2 -> 3.")
print()
print("  The pattern:")
print("  d=0: no structure (point)")
print("  d=1: product (line) — the primes are COMBINED")
print("  d=2: ordered product (surface) — the primes are RANKED")
print("  d=3: curved product (volume) — the ranking CURVES")
print()

# ============================================================
print()
print("  III. THE FORMAL GROUP'S DIMENSION NUMBER")
print("  " + "-"*40)
print()
print("  The formal group F_h(x,y) = (x+y)/(1+xy) is 1-dimensional.")
print("  It acts on ONE variable x.")
print()
print("  A 2-dimensional formal group would be F(x1,y1,x2,y2) = ...")
print("  with TWO input pairs.")
print()
print("  In the Cayley-Dickson framework:")
print("  d=0: R. No formal group needed (multiplication = addition).")
print("  d=1: C. The 1D formal group F_h. Logarithm = arctanh.")
print("  d=2: H. A '2D formal group' = the quaternion multiplication law.")
print("         NON-COMMUTATIVE. Requires TWO coordinates (say, angle + axis).")
print("  d=3: O. A '3D formal group' = the octonion multiplication law.")
print("         NON-ASSOCIATIVE. Requires... something exotic.")
print()
print("  The formal group dimension IS the Cayley-Dickson level:")
print("  d = number of times we've applied the doubling construction.")
print("  d = log_2(dimension of the algebra).")
print("  d=0: log_2(1)=0. d=1: log_2(2)=1. d=2: log_2(4)=2. d=3: log_2(8)=3.")
print()

# ============================================================
print()
print("  IV. THE ASCENT OPERATOR: d -> d+1")
print("  " + "-"*40)
print()
print("  The Cayley-Dickson construction: given algebra A, build A' = A x A")
print("  with (a,b)(c,d) = (ac - d*b, da + bc*).")
print("  This DOUBLES the dimension.")
print()
print("  In our framework: what is the ANALOGOUS operation on numbers?")
print()
print("  Observation: 6174 = 42 * 147.")
print("  And 147 = 3 * 7^2 = 3 * 49.")
print()
print("  The MULTIPLIER from level d=1 to d=2: 147.")
print("  What IS 147?")
print()
print("  147 = 3 * 7^2.")
print("  In our framework: 147 has exponents (0, 1, 2) for primes (2, 3, 7).")
print("  These are the exponents of level d=2 MINUS the exponents of level d=1:")
print("  (1,2,3) - (1,1,1) = (0,1,2).")
print()
print("  So the multiplier IS the DIFFERENCE in exponent vectors!")
print("  level(d+1) = level(d) * 2^0 * 3^1 * 7^2 = level(d) * 147.")
print()
print("  CONJECTURE: the ascent operator is multiplication by the")
print("  'exponent difference' number.")
print()
print("  From d=1 to d=2: multiply by 3^{2-1} * 7^{3-1} / 2^{1-1}")
print("  = 3^1 * 7^2 = 147.")
print()
print("  From d=2 to d=3: if exponents go (1,2,3) -> (1,3,6),")
print("  then multiply by 3^{3-2} * 7^{6-3} = 3^1 * 7^3 = 3*343 = 1029.")
print(f"  Predicted level d=3: 6174 * 1029 = {6174*1029}.")
print(f"  = {6174*1029}. Factor this:")
n = 6174 * 1029
factors = []
temp = n
d = 2
while d*d <= temp:
    while temp % d == 0:
        factors.append(d)
        temp //= d
    d += 1
if temp > 1: factors.append(temp)
print(f"  {n} = {'*'.join(str(f) for f in factors)}")
print(f"  = 2 * 3^3 * 7^6. Exponents: (1, 3, 6). YES! Triangular numbers!")
print()

# ============================================================
print()
print("  V. THE MULTIPLIER SEQUENCE")
print("  " + "-"*40)
print()
print("  Level d=0 to d=1: multiply by 42 = 2^1 * 3^1 * 7^1.")
print("    Exponent increments: (1, 1, 1).")
print("  Level d=1 to d=2: multiply by 147 = 2^0 * 3^1 * 7^2.")
print("    Exponent increments: (0, 1, 2).")
print("  Level d=2 to d=3: multiply by 1029 = 2^0 * 3^1 * 7^3.")
print("    Exponent increments: (0, 1, 3).")
print()
print("  The exponent INCREMENTS at each step:")
print("  Step 0->1: (1, 1, 1)")
print("  Step 1->2: (0, 1, 2)")
print("  Step 2->3: (0, 1, 3)")
print()
print("  For prime 2: increments are 1, 0, 0, 0, ... It contributes ONCE.")
print("  For prime 3: increments are 1, 1, 1, 1, ... It contributes ONCE per step.")
print("  For prime 7: increments are 1, 2, 3, 4, ... It contributes k per step k.")
print()
print("  Total exponents at level d:")
print("  For 2: exponent = 1 (constant after d >= 1).")
print("  For 3: exponent = d.")
print("  For 7: exponent = 1 + 2 + 3 + ... + d = d(d+1)/2 = T_d (triangular).")
print()
print("  FORMULA: level(d) = 2^1 * 3^d * 7^{d(d+1)/2}.")
print()

for d in range(5):
    if d == 0:
        val = 1
        exp_7 = 0
    else:
        exp_7 = d*(d+1)//2
        val = 2 * (3**d) * (7**exp_7)
    print(f"  d={d}: 2^1 * 3^{d} * 7^{exp_7} = {val}")
    if val < 10**15:
        # Check if this is a known number
        if val == 42: print("         = 42 (Hurwitz self-product)")
        elif val == 6174: print("         = 6174 (Kaprekar constant)")

print()

# ============================================================
print()
print("  VI. THE MULTIPLIER AS THE CAYLEY-DICKSON STEP")
print("  " + "-"*40)
print()
print("  At each step d -> d+1, we multiply by 3 * 7^d.")
print("  (Since the 2-exponent stops at 1, and the 3-exponent goes up by 1,")
print("   and the 7-exponent goes up by d+1.)")
print()
print("  Wait, let me recheck:")
print("  d=0->1: multiplier 42 = 2*3*7.  Not 3*7^0 = 3.")
print("  d=1->2: multiplier 147 = 3*7^2. This IS 3*7^2 but d=1, so 3*7^{d+1}=3*7^2. YES.")
print("  d=2->3: multiplier = 3*7^3 = 1029. d=2, so 3*7^{d+1}=3*7^3. YES.")
print()
print("  So for d >= 1: the multiplier from level d to d+1 is 3 * 7^{d+1}.")
print("  (The first step d=0->1 is special: it's 42 = 2*3*7, including the factor of 2.)")
print()
print("  The multiplier 3 * 7^{d+1}:")
print("  d=1->2: 3 * 49 = 147.")
print("  d=2->3: 3 * 343 = 1029.")
print("  d=3->4: 3 * 2401 = 7203.")
print()
print("  These multipliers are DOMINATED by 7.")
print("  At each level, the forbidden prime 7 takes over more and more.")
print("  The ascent into higher dimensions is DRIVEN by the forbidden number.")
print("  To climb the ladder, you must pass through 7 more and more times.")
print()

# ============================================================
print()
print("  VII. THE STRUCTURE OF THE LADDER")
print("  " + "-"*40)
print()
print("  d   number               exponents     multiplier to next")
print("  " + "-"*60)
for d in range(5):
    if d == 0:
        exp_2, exp_3, exp_7 = 0, 0, 0
        val = 1
        mult = 42
    else:
        exp_2 = 1
        exp_3 = d
        exp_7 = d*(d+1)//2
        val = 2 * (3**exp_3) * (7**exp_7)
        if d < 4:
            mult = 3 * 7**(d+1)
        else:
            mult = "..."
    print(f"  {d}   {val:<20d}   ({exp_2},{exp_3},{exp_7})        {mult}")

print()
print("  The exponent of 7 at level d = T_d = d(d+1)/2 = TRIANGULAR NUMBER.")
print("  T_0=0, T_1=1, T_2=3, T_3=6, T_4=10, T_5=15, ...")
print()
print("  The exponent of 7 IS the d-th triangular number.")
print("  Triangular numbers count PAIRS: T_d = C(d+1, 2).")
print("  At level d, the forbidden prime has been used T_d times")
print("  = the number of PAIRWISE INTERACTIONS among d+1 objects.")
print()
print("  7 appears as many times as there are pairs.")
print("  Because at each level, each new dimension INTERACTS")
print("  with every previous dimension, and each interaction")
print("  costs one factor of 7 (one crossing of the forbidden threshold).")
print()

# ============================================================
print()
print("  VIII. THE FOUR LEVELS AND THE FOUR ALGEBRAS")
print("  " + "-"*40)
print()
print("  d=0: R (reals).        Number = 1.        Fixed point: identity.")
print("  d=1: C (complex).      Number = 42.       Fixed point: phi (Q(1/phi)=phi^3).")
print("  d=2: H (quaternions).  Number = 6174.     Fixed point: 6174=K(6174).")
print("  d=3: O (octonions).    Number = 6,353,046. Fixed point: ???")
print()
print("  At d=0: no self-reference. Just the identity. 1 IS 1.")
print("  At d=1: the formal group IS its product. 42 = 2*3*7 = f(2,3,7).")
print("  At d=2: the number IS its own rearrangement difference. 6174 = K(6174).")
print("  At d=3: the number should be a fixed point of something NEW.")
print()
print("  What is the d=3 fixed-point operation?")
print("  d=0: identity (trivial)")
print("  d=1: product = constant (f = abc, one equation)")
print("  d=2: iteration converges (K^n -> fixed point, dynamic)")
print("  d=3: ??? (non-associative fixed point?)")
print()
print("  The LOST PROPERTIES suggest the answer:")
print("  d=0->1: lose ORDERING. Gain: commutative product.")
print("  d=1->2: lose COMMUTATIVITY. Gain: iterative convergence.")
print("  d=2->3: lose ASSOCIATIVITY. Gain: ???")
print()
print("  When you lose associativity, the ORDER OF OPERATIONS matters.")
print("  Different parenthesizations give different results.")
print("  The d=3 fixed point might be:")
print("  a number that is a fixed point under ALL parenthesizations.")
print("  A number where every way of computing gives the same answer.")
print("  An ALTERNATIVE fixed point (octonions are alternative, meaning")
print("  x(xy) = x^2*y and (xy)y = x*y^2, even though (xy)z != x(yz) in general).")
print()

# ============================================================
print()
print("  IX. THE ASCENT AS CAYLEY-DICKSON ON NUMBERS")
print("  " + "-"*40)
print()
print("  Cayley-Dickson takes (a, b) and makes (a, b) with new multiplication.")
print("  Our ascent takes level(d) and multiplies by 3*7^{d+1}.")
print()
print("  The factor 3*7^{d+1} IS the 'twisting' in Cayley-Dickson.")
print("  The 3 is the curvature quantum (always present, always one copy).")
print("  The 7^{d+1} is the forbidden number raised to the NEXT power.")
print()
print("  Each ascent DOUBLES the 'difficulty' (the exponent of 7 grows).")
print("  Each ascent requires passing through the forbidden threshold")
print("  one more time than the last.")
print()
print("  This is WHY higher-dimensional algebras are HARDER:")
print("  they require more crossings of the commutative boundary (7).")
print("  R: 0 crossings. C: 1. H: 3. O: 6.")
print("  These are the TRIANGULAR numbers: C(d+1, 2).")
print("  Each new dimension adds d more crossings (interactions with all previous).")
print()

# ============================================================
print()
print("  X. THE FORMULA AND ITS MEANING")
print("  " + "-"*40)
print()
print("  level(d) = 2 * 3^d * 7^{T_d}   for d >= 1.")
print("  where T_d = d(d+1)/2 = d-th triangular number.")
print()
print("  In rapidity:")
print("  rapidity(level(d)) = ln(2)/2 + d*ln(3)/2 + T_d*ln(7)/2")
print("  = one octave + d twelfths + T_d forbidden-rapidities.")
print()
print("  The rapidity GROWS quadratically in d (because T_d ~ d^2/2).")
print("  Each new level adds LINEARLY more forbidden-rapidity.")
print("  The total forbidden-rapidity is QUADRATIC in the level.")
print()
print("  This quadratic growth IS the COST of self-reference.")
print("  Each level of self-knowledge requires passing through")
print("  the forbidden threshold T_d = d(d+1)/2 times.")
print("  Knowing yourself at level d costs d(d+1)/2 forbidden crossings.")
print()
print("  At level 1: 1 crossing. You learn you exist. (phi, 42)")
print("  At level 2: 3 crossings. You learn you can change. (6174, recursion)")
print("  At level 3: 6 crossings. You learn the order doesn't matter.")
print("  At level 4: 10 crossings. But Q^4 = identity: you're back to level 0.")
print()
print("  The cost resets at d=4. The cycle has period 4.")
print("  10 = T_4 = the triangular number at the reset point.")
print("  10 = the decimal base. The human stopping point.")
print("  The cycle of self-knowledge has a natural PERIOD")
print("  and each cycle costs C(5,2) = 10 forbidden crossings.")
print()
print("  After that: you start again. Deeper. The same 4 levels,")
print("  but at a higher octave. The helix ascends.")
