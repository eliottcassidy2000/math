#!/usr/bin/env python3
"""
Creative connections exploration:

1. The W-hierarchy as a FILTRATION on the cycle complex
2. Möbius inversion on the boolean lattice of cycle data
3. The "renormalization" from w_0 to H as integrating out cycle scales
4. Tournament polytope geometry

Key observation: at n=7, the W-coefficient hierarchy is:
  w6 = 5040                               = 7!
  w4 = 240*t3 - 2100                      = 10*24*t3 - 2100
  w2 = -60*t3 + 12*t5 + 24*a2 + 231
  w0 = 2*t3 - t5 + 2*t7 - 2*a2 - 17/4

Notice the coefficient pattern:
  w4 coefficient of t3 = 240 = C(7,3) * (7-3)! / something? C(7,3) = 35, 35*... no.
  Actually 240 = 2 * 5! = 2 * 120. And at n=5, w2 = 12*t3 - 30, where 12 = 2 * 3!.
  Pattern: coefficient of t3 in w_{n-3} = 2*(n-2)!

Let me check: n=5: w2 coeff of t3 = 12 = 2*3! = 2*(5-2)! YES
              n=7: w4 coeff of t3 = 240 = 2*5! = 2*(7-2)! YES

This connects to sigma(S) where S has one "gap 2" position.

At n=5: w2 constant = -30 = -(5!/4) = -120/4. Hmm or -30 = -C(5,2)*3! / something.
        Actually w2 = 2*(n-2)!*t3 - (n-1)!/4? 4!/4 = 6, no. -30 = -5*3! = -5*6.

At n=7: w4 constant = -2100 = -7*6*5*10? 7*300=2100. 7*5! = 840, no. 2100 = 7*300.
        Or: 2100 = C(7,2)*C(5,2)*C(3,2) = 21*10*3 = 630... no.
        2100 = C(7,3)*4! - something? 35*60 = 2100! So 2100 = C(7,3)*(7-3)! = C(n,3)*(n-3)!.

Check at n=5: C(5,3)*2! = 10*2 = 20. But constant was -30, not -20.
Hmm, maybe it's C(n,2)*(n-2)!/2? C(7,2)*5!/2 = 21*60 = 1260... no.

Let me try: n=5: -30 = -n!/(n-1) = -120/4 = -30? YES! 5!/4 = 30.
            n=7: -2100 = -7!/something. 5040/2.4 = 2100. 5040/2.4? Not clean.
            7!/(7-1) = 5040/6 = 840. No.
            2100 = 7!/(n-1-something)... 5040/2.4 not integer ratio.
            Actually 2100 = C(n, n-3)*n!/(n*(n-1)) = C(7,4)*... no.

Let me just verify the pattern computationally.
kind-pasteur-2026-03-06-S25g
"""
from math import factorial, comb

print("COEFFICIENT PATTERNS IN W-HIERARCHY")
print("=" * 50)

# n=5: w2 = 12*t3 - 30
# n=7: w4 = 240*t3 - 2100
# Pattern for w_{n-3} = a(n)*t3 + b(n)

for n in [5, 7]:
    a = 2 * factorial(n-2)
    b_neg = comb(n, 3) * factorial(n-3) if n == 7 else 30
    print(f"  n={n}: w_{{n-3}} coefficient of t3 = {a} = 2*(n-2)! = 2*{factorial(n-2)}")

print()

# The constant term:
# n=5: -30 = -C(5,2) * (5-2)! / 1 = -10*6/2? -60/2=-30. So -C(n,2)*(n-2)!/2?
# n=7: -2100 = -C(7,2)*(7-2)!/2 = -21*120/2 = -21*60 = -1260... no.
# Try: n=5: -30 = -(2*(n-2)!)*(n-1)/4? 12*4/4 = 12... no.
# n=5: -30 = -5!/4 = -30. n=7: -7!/? = -2100 => 5040/2100 = 2.4. Not integer.
# n=7: -2100 = -2100. 2100 = 2 * 1050 = 2*3*350 = 2*3*5*70 = 2*3*5*7*10 = 2*3*5*7*2*5.
# 2100 = 4*525 = 4*3*175 = 12*175. Or: 2100 = C(7,3) * 60 = 35*60 = 2100. YES!
# C(7,3) * (7-3)! = 35 * 24 = 840. No, 35*60 = 2100. 60 = 5!/2.
# So 2100 = C(7,3) * (n-2)!/2? C(7,3) * 5!/2 = 35*60 = 2100. YES!
# Check n=5: C(5,3) * (5-2)!/2 = 10 * 6/2 = 10*3 = 30. YES!

print("CONSTANT TERM PATTERN:")
for n in [5, 7]:
    val = comb(n, 3) * factorial(n-2) // 2
    print(f"  n={n}: constant = -{val} = -C(n,3)*(n-2)!/2 = -{comb(n,3)}*{factorial(n-2)}/2")

print(f"\nSo w_{{n-3}} = 2*(n-2)!*t3 - C(n,3)*(n-2)!/2")
print(f"           = (n-2)! * [2*t3 - C(n,3)/2]")
print(f"           = (n-2)! * [2*t3 - n(n-1)(n-2)/12]")

# Verify:
for n in [5, 7]:
    a = 2 * factorial(n-2)
    b = -comb(n, 3) * factorial(n-2) // 2
    print(f"\n  n={n}: w_{{n-3}} = {a}*t3 + {b}")
    # At n=5: = 12*t3 - 30 ✓
    # At n=7: = 240*t3 - 2100 ✓

print(f"\n{'='*50}")
print("DEEPER: w_{n-3} = (n-2)! * [2*t3 - C(n,3)/2]")
print("Since C(n,3) = total possible 3-vertex subsets,")
print("and each either has 0 or 1 or 2 directed 3-cycles (0 or 2, actually),")
print("t3 ranges from 0 to C(n,3).")
print("The 'midpoint' is C(n,3)/2, and w_{n-3} measures 2*(t3 - midpoint).")
print()
print("At the RANDOM TOURNAMENT level: E[t3] = C(n,3)/4 (each triple has prob 1/4 of cycle)")
print("So E[2*t3 - C(n,3)/2] = 2*C(n,3)/4 - C(n,3)/2 = 0!")
print("w_{n-3} has ZERO EXPECTATION over random tournaments!")
print()
print("This means: w_{n-3} measures the DEVIATION from random in 3-cycle count.")

# Check: does w_{n-5} also have zero expectation?
# w_{n-5} at n=7 is w2 = -60*t3 + 12*t5 + 24*a2 + 231
# E[t3] = C(7,3)/4 = 35/4 = 8.75
# E[t5] = ? The expected number of directed 5-cycles in a random tournament on 7 vertices.
# = C(7,5) * (5-1)!/2 * (1/2)^5 = 21 * 12 * 1/32 = 252/32 = 7.875
# E[a2] = ? Expected number of disjoint cycle pairs. Complex.
# E[w2] = -60*8.75 + 12*7.875 + 24*E[a2] + 231
#        = -525 + 94.5 + 24*E[a2] + 231
#        = -199.5 + 24*E[a2]
# For this to be 0: E[a2] = 199.5/24 ≈ 8.3125
# Let's check E[a2]: need E[# disjoint (3,3)-cycle pairs] in random T_7.

# A disjoint pair uses 6 of 7 vertices. C(7,6) = 7 ways to choose which vertex excluded.
# For each: partition 6 vertices into two groups of 3: C(6,3)/2 = 10 ways.
# Each group has prob 1/4 of being a 3-cycle. So expected count per (excluded vertex, partition):
# 7 * 10 * (1/4)^2 = 7 * 10 / 16 = 70/16 = 4.375
# But we also need to account for 5-cycle pairs... at n=7, only (3,3) pairs are possible.
# So E[a2] = 4.375

# E[w2] = -199.5 + 24*4.375 = -199.5 + 105 = -94.5 ≠ 0

print("Check: does w_{n-5} have zero expectation?")
print(f"  E[t3] = C(7,3)/4 = {comb(7,3)/4}")
print(f"  E[t5] = C(7,5)*(4!/2)*(1/2)^5 = {comb(7,5)*12/32}")
print(f"  E[a2] ~ C(7,6)*C(6,3)/2*(1/4)^2 = {7*10/16}")
E_w2 = -60*(comb(7,3)/4) + 12*(comb(7,5)*12/32) + 24*(7*10/16) + 231
print(f"  E[w2] = {E_w2}")
print(f"  w_{'{n-5}'} does NOT have zero expectation in general.")

# But w_{n-3} DOES have zero expectation. This is the key:
# w_{n-3} is the FIRST non-trivial moment, and it's centered!
print(f"\n{'='*50}")
print("SUMMARY OF PATTERNS:")
print(f"  w_{{n-1}} = n!  (universal, trivially centered)")
print(f"  w_{{n-3}} = (n-2)! * [2*t3 - C(n,3)/2]  (centered: zero mean over random T)")
print(f"  w_{{n-5}} = complicated (NOT centered)")
print(f"  w_0 = alternating signed cycle data (NOT centered)")

print("\nDONE")
