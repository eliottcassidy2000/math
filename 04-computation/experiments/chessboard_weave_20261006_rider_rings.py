#!/usr/bin/env python3
"""Which ring does the (a,b)-rider line t -> t(a,b) (mod 1) reach on the 2m x 2m subdivision
of the unit torus cell?  Exact: the deepest point has loneliness delta(a,b) = max_t min(||ta||,||tb||);
the line meets the CLOSED central 2x2 block iff delta >= (m-1)/(2m).  delta is computed exactly over
candidate times (crossings m/(a+b), m/(b-a), peaks (2k+1)/(2a), (2k+1)/(2b))."""
from fractions import Fraction as Fr
from math import gcd

def dist(x):
    x = x - (x.numerator // x.denominator)
    return min(x, 1 - x)

def delta2(a, b):
    cands = set()
    for den in (a + b, abs(b - a), 2 * a, 2 * b):
        if den == 0:
            continue
        for k in range(0, den + 1):
            cands.add(Fr(k, den))
    return max(min(dist(t * a), dist(t * b)) for t in cands if 0 <= t <= 1)

bad = []
for a in range(1, 41):
    for b in range(a + 1, 41):
        if gcd(a, b) != 1:
            continue
        d = delta2(a, b)
        s = a + b
        assert d == Fr(s // 2, s), (a, b, d)
        bad.append((a, b, d))
print("delta(a,b) = floor((a+b)/2)/(a+b) verified for all coprime 1<=a<b<=40:", len(bad), "pairs")

def deepest_ring(d, m):
    # loneliness lam lies in cell-depth floor(2m*lam) from the boundary (closed cells: a point at
    # exact depth k/(2m) touches cells of depth k and k-1); ring = (m-1) - depth, clipped at 0
    depth = (2 * m * d).numerator // (2 * m * d).denominator   # floor
    return max(0, (m - 1) - depth)

for m in (2, 3, 4, 5, 6, 8):
    fails = sorted({(a, b) for (a, b, d) in bad if deepest_ring(d, m) > 0 and a + b <= 40})
    print(f"2m={2*m:2d}: riders (a<b<=40 coprime) NOT reaching the closed central 2x2 block: "
          f"{[(a,b) for (a,b) in fails][:12]}{' ...' if len(fails)>12 else ''}  (predicted: odd a+b < {m})")
print("8x8: knight (1,2) deepest ring =", deepest_ring(Fr(1, 3), 4), "; zebra (2,3):", deepest_ring(Fr(2, 5), 4),
      "; giraffe (1,4):", deepest_ring(Fr(2, 5), 4), "; camel (1,3):", deepest_ring(Fr(1, 2), 4))
