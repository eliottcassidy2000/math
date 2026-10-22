"""Exact arithmetic for the cited Kervaire–Milnor formula; no H-count claim."""
from fractions import Fraction
from math import comb

B = [Fraction(1)]
for n in range(1, 9):
    B.append(-sum(Fraction(comb(n+1, j))*B[j] for j in range(n))/(n+1))
assert B[4] == Fraction(-1,30)
assert B[6] == Fraction(1,42)
assert B[8] == Fraction(-1,30)
print('m | B_(2m) | Mersenne factor | abs(4 B_(2m)/m) | numerator | bP_(4m) order | constant-one E_(2m) q coefficient')
for m, expected in [(2,28),(3,992),(4,8128)]:
    mersenne = 2**(2*m-1)-1
    bernoulli_factor = abs(4*B[2*m]/m)
    order = 2**(2*m-2)*mersenne*bernoulli_factor.numerator
    assert order == expected
    print(m, B[2*m], mersenne, bernoulli_factor, bernoulli_factor.numerator, order, -4*m/B[2*m], sep=' | ')
assert 2**6-1 == 63
assert all(2*m-1 != 6 for m in range(2, 100))
print('PASS: all three exact orders; 63 has even exponent 6, not 2m-1 for integer m.')
print('Theta_15 order from Milnor Table 1: 2*8128 =',2*8128)
