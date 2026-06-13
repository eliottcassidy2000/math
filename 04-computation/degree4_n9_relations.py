#!/usr/bin/env python3
"""
DEGREE-4 LINEAR RELATIONS AT n=9

From the degree4_n9_deep.py output, we have the degree-4 Fourier coefficients
of each OCF invariant on the two basis types:

                        Type P    Type Q
  [d4 of t5]            1         0
  [d4 of t7]            3         3
  [d4 of t9]            1.5       3
  [d4 of a33]           0         1
  [d4 of a35]           1         3
  [d4 of a333]          0         0.25

This gives the exact linear relations in the 2D degree-4 space.
Let P = [d4 of t5], Q = [d4 of a33] be the basis.

Then:
  [d4 of t7] = 3*P + 3*Q
  [d4 of t9] = (3/2)*P + 3*Q
  [d4 of a35] = 1*P + 3*Q
  [d4 of a333] = 0*P + (1/4)*Q

Verification of OCF at degree 4:
  w_4 = 240*P + 480*Q
  OCF RHS at degree 4 = 32*(P + 3P+3Q + 1.5P+3Q) + 64*(Q + P+3Q) + 128*(0.25Q)
                       = 32*(5.5P + 6Q) + 64*(P + 4Q) + 128*(0.25Q)
                       = 176P + 192Q + 64P + 256Q + 32Q
                       = 240P + 480Q  CHECK!

opus-2026-03-07-S35
"""
import numpy as np

print("="*70)
print("DEGREE-4 LINEAR RELATIONS AT n=9")
print("="*70)

# Coefficients on (P, Q) basis
inv_coeffs = {
    't5':   (1,   0),
    't7':   (3,   3),
    't9':   (1.5, 3),
    'a33':  (0,   1),
    'a35':  (1,   3),
    'a333': (0,   0.25),
}

print("\nDegree-4 invariants expressed in (P, Q) basis:")
print(f"  P = [deg-4 of t5]  (5-vertex spanning paths)")
print(f"  Q = [deg-4 of a33] (disjoint P2 pairs)")
print()

for name, (cp, cq) in inv_coeffs.items():
    terms = []
    if cp != 0:
        if cp == 1:
            terms.append("P")
        elif cp == int(cp):
            terms.append(f"{int(cp)}*P")
        else:
            from fractions import Fraction
            terms.append(f"{Fraction(cp).limit_denominator(100)}*P")
    if cq != 0:
        if cq == 1:
            terms.append("Q")
        elif cq == int(cq):
            terms.append(f"{int(cq)}*Q")
        else:
            from fractions import Fraction
            terms.append(f"{Fraction(cq).limit_denominator(100)}*Q")
    if not terms:
        terms = ["0"]
    print(f"  [d4 of {name:5s}] = {' + '.join(terms)}")

# Verify OCF
print("\n" + "="*70)
print("OCF VERIFICATION AT DEGREE 4")
print("="*70)

# OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
# alpha_1 = t3 + t5 + t7 + t9
# alpha_2 = a33 + a35
# alpha_3 = a333
#
# At degree 4:
# H_d4 = w_4 / 2^4 = w_4/16
# OCF_d4 = 2*[d4 of alpha_1] + 4*[d4 of alpha_2] + 8*[d4 of alpha_3]
#         = 2*([d4 of t5]+[d4 of t7]+[d4 of t9])
#         + 4*([d4 of a33]+[d4 of a35])
#         + 8*[d4 of a333]

# LHS: w_4/16 at each monomial type
w4_P = 240 / 16  # = 15
w4_Q = 480 / 16  # = 30

print(f"\nLHS (w_4/16):  Type P = {w4_P},  Type Q = {w4_Q}")

# RHS
rhs_P = (2*(1 + 3 + 1.5) + 4*(0 + 1) + 8*(0))
rhs_Q = (2*(0 + 3 + 3) + 4*(1 + 3) + 8*(0.25))

print(f"RHS (OCF d4):  Type P = {rhs_P},  Type Q = {rhs_Q}")
print(f"Match P: {abs(w4_P - rhs_P) < 1e-10}")
print(f"Match Q: {abs(w4_Q - rhs_Q) < 1e-10}")

# Extract the linear dependence relations
print("\n" + "="*70)
print("EXPLICIT LINEAR DEPENDENCE RELATIONS")
print("="*70)

# Since space is 2D (spanned by P, Q), every invariant's d4 part
# is a linear combo of P and Q. The RELATIONS among invariants are:

# [d4 of t7] = 3*P + 3*Q
# Comparing to t5 and a33: [d4 of t7] = 3*[d4 of t5] + 3*[d4 of a33]

print("\nRelation (A): [d4 of t7] = 3*[d4 of t5] + 3*[d4 of a33]")
print("  (At n=7 this was: [d4 of t7] = (1/2)*[d4 of t5] + 1*[d4 of a33])")
print()

print("Relation (B): [d4 of t9] = (3/2)*[d4 of t5] + 3*[d4 of a33]")
print("  (New at n=9)")
print()

print("Relation (C): [d4 of a35] = 1*[d4 of t5] + 3*[d4 of a33]")
print("  (New at n=9)")
print()

print("Relation (D): [d4 of a333] = (1/4)*[d4 of a33]")
print("  (New at n=9)")

# Pattern analysis
print("\n" + "="*70)
print("PATTERN ANALYSIS: n=7 vs n=9")
print("="*70)

print("\n  n=7 coefficients (from degree4_identity_n7.py):")
print("    [d4 of t5] = 1*P + 0*Q")
print("    [d4 of t7] = (1/2)*P + 1*Q")
print("    [d4 of a33] = 0*P + 1*Q")
print("    w_4 coeff on P: 12")
print("    w_4 coeff on Q: 24")

print("\n  n=9 coefficients:")
print("    [d4 of t5] = 1*P + 0*Q")
print("    [d4 of t7] = 3*P + 3*Q")
print("    [d4 of t9] = (3/2)*P + 3*Q")
print("    [d4 of a33] = 0*P + 1*Q")
print("    [d4 of a35] = 1*P + 3*Q")
print("    [d4 of a333] = 0*P + (1/4)*Q")
print("    w_4 coeff on P: 240")
print("    w_4 coeff on Q: 480")

# General formulas:
# w_4 coeff on P = (n-4)! * 2 * (#block positions)
# At n=7: #block positions = n-4 = 3, extra vertices = 2, orderings = 2!
# w_4 on P = 3 * 2 * 2! = 12 CHECK
# At n=9: #block positions = 5, extra = 4, orderings = 4!
# w_4 on P = 5 * 2 * 24 = 240 CHECK
#
# General: w_4 on P = (n-4) * 2 * (n-5)!

print("\n  General w_4 coefficient on Type P:")
print("    = (n-4) * 2 * (n-5)!")
print(f"    n=7: {3*2*2} = 12  CHECK")
print(f"    n=9: {5*2*24} = 240  CHECK")

# w_4 on Q: two P2 blocks (3 each) + (n-6) singles
# Interleaving: (n-6+2)! = (n-4)! arrangements of (n-4) objects
# Block orientations: 2*2 = 4
# w_4 on Q = (n-4)! * 4
print("\n  General w_4 coefficient on Type Q:")
print("    = (n-4)! * 4")
print(f"    n=7: {6*4} = 24  CHECK")
print(f"    n=9: {120*4} = 480  CHECK")

# The ratio w_4(Q)/w_4(P) = (n-4)! * 4 / ((n-4) * 2 * (n-5)!)
# = 4 / (2 * (n-4)) * (n-4)! / (n-5)! ... wait
# = (n-4)! * 4 / ((n-4) * 2 * (n-5)!) = 4 * (n-5)! / (2 * (n-5)!) = 2
print("\n  Ratio w_4(Q)/w_4(P) = 2 for all n")
print(f"    n=7: {24/12}")
print(f"    n=9: {480/240}")

print("\n" + "="*70)
print("COEFFICIENT TABLE FOR GENERAL PATTERN")
print("="*70)

# The t_k coefficient on Type P (at a specific P4 monomial):
# t_k contributes (1/2)^{k-4} * (signed count of directed k-cycles through the P4)
# For P4 = path 0-1-2-3-4:
# k=5: 2 cycles, coeff = 2 * (1/2)^1 = 1
# k=7: 24 cycles, coeff = 24 * (1/2)^3 = 3
# k=9: 48 cycles, coeff = 48 * (1/2)^5 = 1.5

# For k=7: C(4,2)=6 ways to choose 2 extra vertices from {5,6,7,8}
# Each gives 2 orderings of the extras in the cycle, times 2 orientations = 4
# Total: 6*4 = 24
# For k=9: the 4 extra vertices must all be in the cycle.
# Cycle is v0-v1-v2-v3-v4-w1-w2-w3-w4 (in some order).
# Number of directed 9-cycles containing the P4 subpath: ?
# 48 = ... Let me compute: (n-5)! * 2 = 4! * 2 = 48? Yes!
# (The P4 must appear contiguously in the cycle, extra 4 vertices ordered in (n-5)! ways,
#  times 2 orientations, but also the block can be placed in different cycle positions...)

# Actually for k=n=9 (Hamiltonian cycle containing P4):
# P4 must be contiguous in the cycle (same argument as for paths).
# Cycle has 9 positions. Fix the P4 block: since cycle is symmetric,
# # of distinct placements = 1 (all rotations equivalent).
# But wait, we fix v0=0 in the enumeration. Let me re-derive from the count = 48.
# 48 directed cycles / (1/2)^5 would give 48/32 = 1.5 as the coefficient.

# For general k odd, the t_k coefficient on Type P should be:
# = (signed count of k-cycles through P4) * (1/2)^{k-4}
# = C(n-5, k-5) * 2 * (k-5)! * (1/2)^{k-4}
# = C(n-5, k-5) * (k-5)! * (1/2)^{k-5}
# = (n-5)! / (n-k)! * (1/2)^{k-5}

# n=9, k=5: (4!/4!) * (1/2)^0 = 1 CHECK
# n=9, k=7: (4!/2!) * (1/2)^2 = 12/4 = 3 CHECK
# n=9, k=9: (4!/0!) * (1/2)^4 = 24/16 = 1.5 CHECK

print("\nCoefficient of [d4 of t_k] on Type P:")
print("  = (n-5)! / (n-k)! * (1/2)^{k-5}")
for k in [5, 7, 9]:
    import math
    val = math.factorial(4) / math.factorial(9-k) * (0.5)**(k-5)
    print(f"  k={k}: {val:.4f}")

# For Type Q coefficient of t_k:
# k=5: 0 (can't fit 5-cycle through 6 disjoint vertices)
# k=7: 3
# k=9: 3
print("\nCoefficient of [d4 of t_k] on Type Q:")
# k=7: how many 7-cycles through the two P2s?
# P2 pair uses 6 vertices. 7-cycle needs 7, so 1 extra vertex from {6,7,8}=3 choices.
# Within the 7-cycle, the two P2 blocks must be contiguous.
# Count = 3 (extra vertex choices) * ... = 24 directed cycles (from our computation)
# 24 * (1/2)^3 = 3 CHECK

# k=9: all 9 vertices in cycle. Two P2 blocks contiguous, 3 singles.
# 48 directed cycles * (1/2)^5 = 1.5... but we got 3.
# Wait, the computation said 48 cycles but coeff = 3.
# Let me recheck: 48 cycles at k=9, factor = (1/2)^{9-4} = (1/2)^5 = 1/32
# 48/32 = 1.5. But we got 3.0 for Type Q at k=9.
# So the count of 9-cycles through Type Q monomial edges is different.
# Re-reading: for Type Q, the computation returned:
#   t_9: 48 directed cycles contain these edges, [deg-4 of t_9] = 3.0000... wait
# Let me check: 48 * (1/2)^5 = 48/32 = 1.5 for P, but for Q the count might be 96.
# 96/32 = 3. So 96 directed 9-cycles contain the Type Q edges.
# The script just showed the coefficient, not the raw count for Q.
# But actually looking back at the output for Type Q at t_9:
# "t_9: 48 directed cycles..." -- that was for Type P.
# For Type Q, it just showed the coefficient directly.

# Let me just verify the general formulas
print("\n  k=5: 0 (P2 pair uses 6 vertices, can't fit in 5-cycle)")
print("  k=7: signed_count * (1/2)^3")
print("  k=9: signed_count * (1/2)^5")
print("  (Exact counts need separate computation)")

# Summary of the full relation matrix
print("\n" + "="*70)
print("COMPLETE DEGREE-4 RELATION MATRIX AT n=9")
print("="*70)

print("""
         P-coeff   Q-coeff
t5:        1         0
t7:        3         3
t9:       3/2        3
a33:       0         1
a35:       1         3
a333:      0        1/4

where P = [d4 of t5], Q = [d4 of a33].

Key relations:
  [d4 of t7] = 3*[d4 of t5] + 3*[d4 of a33]
  [d4 of t9] = (3/2)*[d4 of t5] + 3*[d4 of a33]
  [d4 of a35] = [d4 of t5] + 3*[d4 of a33]
  [d4 of a333] = (1/4)*[d4 of a33]

OCF identity at degree 4:
  w_4/16 = 2*(t5 + t7 + t9)_d4 + 4*(a33 + a35)_d4 + 8*(a333)_d4
         = 2*(1+3+3/2)*P + 2*(0+3+3)*Q + 4*(0+1)*P + 4*(1+3)*Q + 8*(0)*P + 8*(1/4)*Q
         = 2*(11/2)*P + 2*(6)*Q + 4*P + 16*Q + 2*Q
         = 11P + 4P  + 12Q + 16Q + 2Q
         = 15P + 30Q
         = (240/16)*P + (480/16)*Q  CHECK!
""")
