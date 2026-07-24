"""
mac-mini-2026-07-23-S168 -- THE EXACT eq (27) of the owner's certified-log snippet.

CONTRIBUTION: klein-S402 and opus-S2 both decoded the snippet correctly (rapidity / Abel-Dini
telescoping: 2*arctanh(t)=log(S_n/S_{n-1}), the t's are (S_n-S_{n-1})/(S_n+S_{n-1})) but BOTH
concluded "no small-coefficient exact fit for RHS(27); the rational part is large." That conclusion
is WRONG -- it comes from searching INTEGER coefficients. The exact fit uses a RATIONAL coefficient:

    RHS(27) = (2457/6592) * log(8847357/2974400)  -  log(1285/896)          [c=2457/6592, d=1, NO constant]

Certified lower bound (the snippet's own bounds):
    (2457/6592)*L_B - U_A - 1/25  ==  391926968594914200867482400554891567498742649630277
                                       ---------------------------------------------------------
                                       82738859282193417287303438726081463937219800938169600
reproduced EXACTLY (all 51 digits). The coefficient 2457/6592 is the UNIQUE low-height solution;
every other integer pairing (a*L_B - b*U_A) needs height ~10^52.

STRUCTURE of the coefficient (the load-bearing new fact):
    2457 = 3^3 * 7 * 13 = 27 * 91 = 3 * S_2({1..13})        [ 91 = 7*13 = C(14,2) = S_1({1..13}) ]
    6592 = 2^6 * 103                                        [ 103 also divides 5872957 = 19*103*3001 ]
So  6592 * RHS(27) = 2457*log_B - 6592*log_A > 6592/25 = 263.68   (the integer linear form).

The partial-sum pairs (opus-S2, verified):
    A:  (S_{n-1}, S_n) = (896, 1285),      x_n = 389,      t_A = 389/2181
    B:  (S_{n-1}, S_n) = (2974400, 8847357), x_n = 5872957, t_B = 5872957/11821757
So alpha, beta are ratios of consecutive partial sums; eq (27) is a strict inequality between two
Abel-Dini log-combinations (THM-2000 sec 3.1) -- certified > 1/25.
"""
from fractions import Fraction as F

# the snippet's two-sided arctanh bounds
def logb(x):
    x = F(x); b, a = x.numerator, x.denominator
    t = F(b - a, b + a)
    lo = 2 * (t + t**3/3 + t**5/5)
    hi = 2 * (t + t**3/3 + t**5/(5*(1 - t**2)))
    return lo, hi

alpha = F(1285, 896)
beta  = F(8847357, 2974400)
loA, U_A = logb(alpha)      # upper bound on the SMALLER (subtracted) log
L_B, hiB = logb(beta)       # lower bound on the LARGER (added) log
c = F(2457, 6592)

given = F(391926968594914200867482400554891567498742649630277,
          82738859282193417287303438726081463937219800938169600)

cert = c * L_B - U_A - F(1, 25)
print("Abel-Dini telescoping check (opus-S2):")
print(f"  1285-896 = {1285-896} = 389 = x_n(A);  1285+896 = {1285+896} = 2181 = q_A  -> t_A=389/2181  OK")
print(f"  8847357-2974400 = {8847357-2974400} = 5872957 = x_n(B);  sum = {8847357+2974400} = 11821757 = q_B  OK")
print()
print("eq (27) reconstructed:  RHS(27) = (2457/6592) log(8847357/2974400) - log(1285/896)")
print(f"  certified lower bound - 1/25  = {cert}")
print(f"  given snippet rational        = {given}")
print(f"  EXACT MATCH: {cert == given}")
print()
print("coefficient anatomy:")
print(f"  2457 = 27*91 = 3*S_2(AP13) ;  91 = 7*13 = C(14,2) = S_1(AP13)")
print(f"  6592 = 2^6*103")
print(f"  integer form:  2457*log_B - 6592*log_A > 6592/25 = {F(6592,25)} = {float(F(6592,25)):.4f}")
import math
print(f"  true RHS(27) = {float(c)*math.log(float(beta)) - math.log(float(alpha)):.6f}  > 1/25={1/25}")
