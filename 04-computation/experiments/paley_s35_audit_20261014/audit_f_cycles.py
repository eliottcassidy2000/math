"""Audit (f2): parity data of the rational Collatz cycles quoted in the atlas note (standard map C(x) = 3x+1 / x/2)."""
from fractions import Fraction as F

def orbit(x, n=40):
    seen = []
    while x not in seen and len(seen) < n:
        seen.append(x)
        x = 3 * x + 1 if x.denominator % 2 == 1 and x.numerator % 2 != 0 else x / 2
    return seen[seen.index(x):]

for start in (F(1, 13), F(-5), F(-17)):
    cyc = orbit(start)
    odd = sum(1 for y in cyc if y.numerator % 2 != 0)
    even = len(cyc) - odd
    print(f'cycle of {start}: length {len(cyc)}, odd steps {odd}, halvings {even}, 2^halvings - 3^odd = {2 ** even - 3 ** odd}')
