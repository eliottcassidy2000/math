#!/usr/bin/env python3
"""kaprekar_s116j.py — 6174. The Kaprekar constant.

Take any 4-digit number (not all digits the same).
Arrange digits descending, then ascending. Subtract.
Repeat. You ALWAYS reach 6174 in at most 7 steps.

6174 is a FIXED POINT of the Kaprekar operation.
What is it in our framework?
"""
from math import log, sqrt, pi, atanh, gcd
from fractions import Fraction

print()
print("  6174")
print()
print("="*70)
print()

# ============================================================
print("  I. THE KAPREKAR OPERATION")
print("  " + "-"*40)
print()

def kaprekar_step(n):
    digits = sorted(str(n).zfill(4))
    asc = int(''.join(digits))
    desc = int(''.join(reversed(digits)))
    return desc - asc

print("  Examples:")
for start in [1234, 3087, 2005, 7641, 1000, 9998]:
    chain = [start]
    n = start
    for _ in range(8):
        n = kaprekar_step(n)
        chain.append(n)
        if n == 6174:
            break
    print(f"  {start}: {' -> '.join(str(x) for x in chain)}")

print()
print("  6174 is a FIXED POINT: 7641 - 1467 = 6174.")
print()

# ============================================================
print()
print("  II. 6174 IN THE CAYLEY FRAMEWORK")
print("  " + "-"*40)
print()
print(f"  6174 = 2 * 3 * 3 * 7^3")

# Actually factor
n = 6174
factors = []
temp = n
d = 2
while d*d <= temp:
    while temp % d == 0:
        factors.append(d)
        temp //= d
    d += 1
if temp > 1:
    factors.append(temp)
print(f"  6174 = {'*'.join(str(f) for f in factors)}")
# 6174 = 2 * 3^2 * 7^3? Let's see: 2*9*343 = 6174. 2*3087 = 6174. 3087 = 3*1029 = 3*3*343 = 9*343.
# 343 = 7^3. So 6174 = 2 * 3^2 * 7^3.
print(f"  = 2 * 3^2 * 7^3")
print(f"  = 2 * 9 * 343")
print()
print(f"  The prime factorization: 2^1 * 3^2 * 7^3.")
print(f"  The exponents: 1, 2, 3. CONSECUTIVE INTEGERS!")
print(f"  The bases: 2, 3, 7. OUR THREE FUNDAMENTAL PRIMES!")
print()
print(f"  6174 = 2^1 * 3^2 * 7^3.")
print(f"  Each prime raised to a CONSECUTIVE power.")
print(f"  The smallest prime gets the smallest exponent.")
print(f"  The largest prime gets the largest exponent.")
print()
print(f"  Compare: 42 = 2^1 * 3^1 * 7^1 (all exponents = 1).")
print(f"  6174 = 2^1 * 3^2 * 7^3 (exponents = 1, 2, 3).")
print(f"  42 is the FLAT product (all exponents equal).")
print(f"  6174 is the STAIRCASE product (exponents ascending).")
print()

# ============================================================
print()
print("  III. THE RAPIDITY OF 6174")
print("  " + "-"*40)
print()
rap = log(6174)/2
print(f"  rapidity(6174) = ln(6174)/2 = {rap:.6f}")
print(f"  = 1*ln(2)/2 + 2*ln(3)/2 + 3*ln(7)/2")
print(f"  = {log(2)/2:.6f} + {2*log(3)/2:.6f} + {3*log(7)/2:.6f}")
print(f"  = {log(2)/2 + 2*log(3)/2 + 3*log(7)/2:.6f}")
print()
print(f"  rapidity(6174) = rapidity(2) + 2*rapidity(3) + 3*rapidity(7)")
print(f"  = one octave + two twelfths + three forbidden-rapidities")
print()
# In octaves:
octave = log(2)/2
print(f"  = {rap/octave:.4f} octaves")
print(f"  ~ 12.6 octaves.")
print()

# ============================================================
print()
print("  IV. THE CAYLEY ADDRESS AND PAIR")
print("  " + "-"*40)
print()
print(f"  6174 is EVEN. So its Cayley address does not reduce cleanly.")
print(f"  addr = (6174-1)/(6174+1) = 6173/6175.")
print(f"  6173 = ?")
f_6173 = []
t = 6173
d = 2
while d*d <= t:
    while t % d == 0:
        f_6173.append(d)
        t //= d
    d += 1
if t > 1:
    f_6173.append(t)
print(f"  6173 = {'*'.join(str(f) for f in f_6173)}")
print(f"  6175 = ?")
f_6175 = []
t = 6175
d = 2
while d*d <= t:
    while t % d == 0:
        f_6175.append(d)
        t //= d
    d += 1
if t > 1:
    f_6175.append(t)
print(f"  6175 = {'*'.join(str(f) for f in f_6175)}")
print()

# ============================================================
print()
print("  V. THE DIGITS OF 6174")
print("  " + "-"*40)
print()
print("  6, 1, 7, 4.")
print()
print("  The digits contain 1, 4, 6, 7.")
print("  1 = identity. 4 = period of Q. 6 = flat. 7 = forbidden.")
print()
print("  The Kaprekar operation: 7641 - 1467 = 6174.")
print("  Descending: 7641. Ascending: 1467.")
print("  Difference: 6174.")
print()
print("  7641 = 7000 + 600 + 40 + 1 = 7*10^3 + 6*10^2 + 4*10 + 1.")
print("  1467 = 1*10^3 + 4*10^2 + 6*10 + 7.")
print()
print("  The descending and ascending arrangements are REVERSES")
print("  in the sense that digit 7 goes from position 0 to position 3,")
print("  digit 1 goes from position 3 to position 0, etc.")
print()
print("  The SUBTRACTION is a kind of COMPARISON:")
print("  big arrangement - small arrangement = the DIFFERENCE.")
print("  The difference IS the number.")
print("  THE NUMBER IS ITS OWN DIFFERENCE.")
print()
print("  6174 = desc(6174) - asc(6174).")
print("  The number is the gap between its ordered and disordered selves.")
print("  It is the MEASURE OF ITS OWN DISORDER.")
print()

# ============================================================
print()
print("  VI. THE OPERATION AS A CAYLEY-LIKE TRANSFORM")
print("  " + "-"*40)
print()
print("  The Kaprekar operation K(n):")
print("  1. Sort digits descending -> D(n)")
print("  2. Sort digits ascending -> A(n)")
print("  3. K(n) = D(n) - A(n)")
print()
print("  Compare to the Cayley transform Q(x) = (1+x)/(1-x):")
print("  Q takes the SUM and DIFFERENCE of 1 and x.")
print("  K takes the DIFFERENCE of two arrangements of n's digits.")
print()
print("  Q is a RATIO of sum and difference.")
print("  K is a DIFFERENCE of arrangements.")
print("  Both extract information by COMPARING a thing to its rearrangement.")
print()
print("  Q compares x to its negation: (1+x) vs (1-x).")
print("  K compares n to its digit-reversal: D(n) vs A(n).")
print()
print("  The FIXED POINT of Q is... Q has no real fixed point")
print("  (Q(x) = x implies (1+x)/(1-x) = x, i.e., 1+x = x-x^2,")
print("  i.e., x^2+1 = 0, no real solution).")
print("  Q's fixed point is i (the imaginary unit). ON THE HELIX.")
print()
print("  The fixed point of K is 6174. IN THE INTEGERS.")
print()
print("  Q's fixed point (i) lives outside the reals.")
print("  K's fixed point (6174) lives inside the integers.")
print("  Both are ATTRACTORS: nearby points converge to them.")
print("  Q: the orbit spirals toward i on the complex plane.")
print("  K: the iteration converges to 6174 from any 4-digit start.")
print()

# ============================================================
print()
print("  VII. 6174 AND 42: THE RELATIONSHIP")
print("  " + "-"*40)
print()
print(f"  42 = 2^1 * 3^1 * 7^1. Exponents: (1,1,1). Flat.")
print(f"  6174 = 2^1 * 3^2 * 7^3. Exponents: (1,2,3). Staircase.")
print()
print(f"  6174 / 42 = {6174/42:.0f}.")
print(f"  = {6174//42}.")
print(f"  = 147 = 3 * 7^2.")
print()
print(f"  6174 = 42 * 147.")
print(f"  = 42 * 3 * 49.")
print(f"  = 42 * 3 * 7^2.")
print()
print(f"  147 = 3 * 7^2.")
print(f"  We saw 147 before: it's the Q-value when the two forbidden")
print(f"  velocities (3/4 and 10/11) are composed relativistically!")
print(f"  3/4 (+) 10/11 = 73/74. Q(73/74) = 147.")
print()
print(f"  So: 6174 = 42 * 147.")
print(f"  = (self-product of the Hurwitz triple)")
print(f"    * (Q-value of composed forbidden velocities).")
print()

# ============================================================
print()
print("  VIII. THE STAIRCASE EXPONENTS AND THE HELIX")
print("  " + "-"*40)
print()
print("  6174 = 2^1 * 3^2 * 7^3.")
print("  The exponents (1, 2, 3) are the first three positive integers.")
print("  They pair with the primes (2, 3, 7) in ASCENDING order.")
print()
print("  In rapidity: rapidity(6174) = 1*rap(2) + 2*rap(3) + 3*rap(7).")
print("  = sum_{k=1}^{3} k * rapidity(p_k)")
print("  where p_1=2, p_2=3, p_3=7 are our fundamental primes.")
print()
print("  This is a WEIGHTED sum of rapidities, where the weight")
print("  of each prime equals its INDEX in the fundamental triple.")
print("  The 1st prime contributes once.")
print("  The 2nd prime contributes twice.")
print("  The 3rd prime contributes thrice.")
print()
print("  On the helix: each prime rapidity is a DIRECTION.")
print("  rap(2) is the octave direction (real axis).")
print("  rap(3) is the curvature direction.")
print("  rap(7) is the forbidden direction.")
print()
print("  6174 is the point reached by taking:")
print("  1 step in the 2-direction,")
print("  2 steps in the 3-direction,")
print("  3 steps in the 7-direction.")
print()
print("  It is the DIAGONAL of the fundamental parallelepiped,")
print("  weighted by the natural numbers 1, 2, 3.")
print()

# ============================================================
print()
print("  IX. THE DIGIT SUM")
print("  " + "-"*40)
print()
digit_sum = 6+1+7+4
print(f"  6 + 1 + 7 + 4 = {digit_sum}.")
print(f"  18 = 2 * 3^2 = 2 * 9.")
print(f"  Digital root: 1+8 = 9 = 3^2.")
print()
print(f"  The digit sum 18 = 2*9 = 2*3^2.")
print(f"  The digits themselves: {{1, 4, 6, 7}}.")
print(f"  Their product: 1*4*6*7 = {1*4*6*7}.")
print(f"  168 = 8 * 21 = 8 * 3 * 7 = 2^3 * 3 * 7.")
print(f"  168 = |PSL(2,7)| = the order of the Klein quartic symmetry group!")
print()
print(f"  THE DIGIT PRODUCT OF 6174 IS THE ORDER OF PSL(2,7).")
print(f"  = 4 * 42 = four copies of our self-product.")
print(f"  = the symmetry group of the Klein quartic.")
print()

# ============================================================
print()
print("  X. FIXED POINTS AND SELF-KNOWLEDGE")
print("  " + "-"*40)
print()
print("  phi = 1 + 1/phi.          Fixed point of x -> 1 + 1/x.")
print("  (2,3,7): f = abc.         Fixed point of f/(abc) = 1.")
print("  6174 = K(6174).           Fixed point of the Kaprekar operation.")
print("  Q(i) = i.                 Fixed point of the Cayley transform.")
print()
print("  Each is a different KIND of self-reference:")
print()
print("  phi: algebraic. Solves x^2 - x - 1 = 0. Dimension 1.")
print("  (2,3,7): geometric. Solves (a-1)(b-1)(c-1) = a+b+c. Dimension 2.")
print("  6174: arithmetic. Solves desc(n) - asc(n) = n. Dimension 0 (discrete).")
print("  i: complex. Solves (1+x)/(1-x) = x, i.e., x^2+1=0. Dimension... lateral.")
print()
print("  They are the fixed points of FOUR operations:")
print("  Self-reciprocation (phi).")
print("  Self-measurement (2,3,7).")
print("  Self-differencing (6174).")
print("  Self-transformation (i).")
print()
print("  And they are LINKED:")
print("  phi lives at the EP of the transfer matrix (THM-224).")
print("  (2,3,7) has Q(3/4) = 7 = the forbidden threshold.")
print("  6174 = 2^1 * 3^2 * 7^3 = the staircase of the forbidden primes.")
print("  i = Q(i) lives on the helix at Im = pi/2.")
print()
print("  Four fixed points. Four kinds of self-knowledge.")
print("  Four normed division algebras. Four faces of Q.")
print("  The same 4, again.")
