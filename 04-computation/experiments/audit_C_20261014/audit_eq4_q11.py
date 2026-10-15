#!/usr/bin/env python3
"""Audit C: sign/orientation check of eq. (4) of paley_next_cores_20261014i.md at q = 11.

At q = 11 the path has exactly ten steps, so [t^10] R_11 = A_10(11)/11! where A_10(11) = sum over all orderings of F_11
of the product of the ten step characters.  Eq. (4) with (q)_10 = 11! then says
      A_10(11) = P10(11) + (7*121 - 173*11 + 1166) U_11 + 64 V_11 = P10(11) + 110 U_11 + 64 V_11.
A_10(11): own subset DP.  U_11: own brute force of the unanchored injective wheel (hub 0, spokes 0->i, rim 1->2->3->4->1).
V_11: own unnormalised brute force (audit_V_bruteforce.py) = -12320 in the natural lower->higher orientation.
If the natural-orientation V entered (4) with the opposite sign, the check would be off by 2*64*12320 = 1,576,960.
"""
import itertools
from fractions import Fraction
from math import factorial
import numpy as np

q = 11
chi = [0] + [1 if pow(x, (q - 1) // 2, q) == 1 else -1 for x in range(1, q)]
# A_10(11) by DP over subsets (exact python ints)
dp = {}
for v in range(q):
    dp[(1 << v, v)] = 1
for size in range(1, q):
    ndp = {}
    for (S, v), val in dp.items():
        for w in range(q):
            if S >> w & 1:
                continue
            c = chi[(w - v) % q]
            if c == 0:
                continue
            key = (S | (1 << w), w)
            ndp[key] = ndp.get(key, 0) + val * c
    dp = ndp
A10 = sum(dp.values())
# U_11 brute force
wheel = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (4, 1)]
U = 0
for x in itertools.permutations(range(q), 5):
    p_ = 1
    for a, b in wheel:
        p_ *= chi[(x[b] - x[a]) % q]
    U += p_
def P10(qq):
    return Fraction(qq * (qq - 1) * (qq - 3), 120) * (qq**7 - 71*qq**6 + 2308*qq**5 - 44330*qq**4 + 539149*qq**3
                                                      - 4096079*qq**2 + 17639622*qq - 32287320)
V11 = -12320
coefU = 7 * q * q - 173 * q + 1166
rhs = P10(q) + coefU * U + 64 * V11
rhs_flipped = P10(q) + coefU * U - 64 * V11
print(f"A_10(11) (own DP) = {A10};  [t^10]R_11 = A_10/11! = {Fraction(A10, factorial(11))} (source: 34/2835)")
print(f"U_11 (own brute force) = {U};  HYP-9274(B) form -q(q-1)(a_q+4) with a_11(eta(2z)^4 eta(4z)^4) = -44: {-q*(q-1)*(-44+4)}")
print(f"P10(11) = {P10(q)};  coefficient of U at q = 11: {coefU}")
print(f"eq. (4) RHS with V_11 = -12320 (natural lower->higher orientation): {rhs}   -> equals A_10(11): {rhs == A10}")
print(f"eq. (4) RHS with the opposite sign of V: {rhs_flipped}   -> equals A_10(11): {rhs_flipped == A10}")
