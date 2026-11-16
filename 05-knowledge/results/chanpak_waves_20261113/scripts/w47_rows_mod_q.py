"""Doubling rows mod an odd q (orbits of x -> 2x on Z/q minus 0).  They can be oriented into a circulant tournament
(choose one orbit from each pair {O, -O}) iff no orbit is closed under negation, iff no divisor d > 1 of q divides any
2^k + 1.  The same question for x -> 3x mod 2^k - 1 ... and the Frobenius of F_9 (x -> 3x on exponents mod 8)."""
from math import gcd
def orbits(q, m):
    seen = set(); orbs = []
    for x in range(1, q):
        if x in seen: continue
        o = []; y = x
        while y not in o:
            o.append(y); y = y * m % q
        seen |= set(o); orbs.append(sorted(o))
    return orbs
def orientable(q, m):
    return all(sorted((-y) % q for y in o) != o for o in orbits(q, m))
qs = [q for q in range(3, 400, 2) if orientable(q, 2)]
print("odd q <= 400 whose doubling rows orient into a circulant tournament:", qs)
print("doubling rows mod 7:", orbits(7, 2), " (QR mod 7 = connection set of the Paley tournament P_7)")
print("doubling rows mod 21:", orbits(21, 2), " orientable:", orientable(21, 2))
print("Frobenius of F_9 = x -> 3x on exponents mod 8:", orbits(8, 3), " negation = +4 on exponents;",
      "self-paired orbit present:", any(sorted((y + 4) % 8 for y in o) == o for o in orbits(8, 3)))
print("Frobenius of F_8 = x -> 2x on exponents mod 7:", orbits(7, 2))
# Collatz +-: in characteristic 2, 3x+1 = 3x-1 = x+1; in characteristic 3, 3x+-1 = +-1
print("char 2: 3x+1 == 3x-1 == x+1 for x in F_2:", all((3 * x + 1) % 2 == (3 * x - 1) % 2 == (x + 1) % 2 for x in range(2)))
