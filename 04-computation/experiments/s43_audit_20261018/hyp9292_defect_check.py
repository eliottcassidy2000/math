# S43 referee: check the small factual claims in HYP-9292 (opus) 'Structure of the defect'
from math import comb
from fractions import Fraction as Fr
from sympy import factorint
def H(m): return sum(Fr(1, j) for j in range(1, m + 1))
def S(n): return sum((-1)**k * comb(n, k)**3 * H(n - k) for k in range(n + 1))
def D(n): return sum((-1)**k * comb(n, k)**3 for k in range(n + 1))
for n in (3, 5, 9, 23):
    s = S(n)
    print('S_%d = %s ; numerator %s ; denominator %s' % (n, s, factorint(s.numerator), factorint(s.denominator)))
E10 = S(10) - H(10) * D(10)
print('E_10 =', E10, 'numerator factors', factorint(E10.numerator))
