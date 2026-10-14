#!/usr/bin/env python3
"""Known Giuga numbers (OEIS A007850, factor rows A236434, plus the 10-factor one in A007850's comments):
are g - 1 and g + 1 prime?  Exact checks of the Giuga condition, sum 1/p - 1/g, the arithmetic derivative g' = g + 1,
and probable primality of g +- 1 (sympy.isprime: deterministic below 2^64, BPSW above).
"""
from fractions import Fraction as F
from math import prod
import sympy

ROWS = [[2,3,5],[2,3,11,13],[2,3,7,41],[2,3,11,17,59],[2,3,11,23,31,47057],[2,3,7,43,3041,4447],
        [2,3,7,59,163,1381,775807],[2,3,7,71,103,67213,713863],[2,3,7,71,103,61559,29133437],
        [2,3,11,23,31,47137,28282147,3892535183],[2,3,11,23,31,47059,2259696349,110725121051],
        [2,3,7,43,1831,138683,2861051,1456230512169437],
        [2,3,11,23,31,47059,2217342227,1729101023519,8491659218261819498490029296021,58254480569119734123541298976556403]]
A007850 = [30,858,1722,66198,2214408306,24423128562,432749205173838,14737133470010574,550843391309130318,
           244197000982499715087866346,554079914617070801288578559178,1910667181420507984555759916338506]
checks = 0
for i, ps in enumerate(ROWS):
    g = prod(ps)
    if i < len(A007850):
        assert g == A007850[i]; checks += 1
    assert all(sympy.isprime(p) for p in ps); checks += 1
    assert all(((g // p) - 1) % p == 0 for p in ps); checks += 1          # Giuga condition
    s = sum(F(1, p) for p in ps) - F(1, g)
    d = sum(g // p for p in ps)                                           # arithmetic derivative of squarefree g
    lo, hi = sympy.isprime(g - 1), sympy.isprime(g + 1)
    print(f"G{i+1:2d}: {len(str(g)):3d} digits, {len(ps):2d} primes; sum 1/p - 1/g = {s}; g' = g + 1: {d == g + 1}; "
          f"g-1 prime: {lo}; g+1 prime: {hi}")
print(f"exact checks: {checks}")
