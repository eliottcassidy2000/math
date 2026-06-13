from itertools import product
from math import prod, factorial, gcd
from fractions import Fraction
from sympy.utilities.iterables import partitions
def A000568(n): 
    return int(sum(Fraction(1<<(sum(p[r]*p[s]*gcd(r, s) for r, s in product(p.keys(), repeat=2))-sum(p.values())>>1), prod(q**p[q]*factorial(p[q]) for q in p)) for p in partitions(n) if all(q&1 for q in p))) # Chai Wah Wu, Jul 01 2024

print(A000568(75))