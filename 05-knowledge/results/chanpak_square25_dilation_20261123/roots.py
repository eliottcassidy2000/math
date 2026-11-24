"""Integer obstructions and exact inverses for the twelve fixed order-eleven seeds."""
from bisect import bisect_left
from math import isqrt

from extensions import RAILS, height

MODULI = (4, 15, 4, 8, 4, 19, 4, 8, 4, 9, 4, 15)


def square_obstruction(coeff, modulus, start):
    """Follow the joint powers until their state repeats, including the transient."""
    two, three = pow(2, start, modulus), pow(3, start, modulus)
    c1, c2 = coeff[:2]
    c3 = coeff[2] if len(coeff) == 3 else 0
    constant, linear = 2 + c1 - 2*c2 + 3*c3, c2 - 3*c3
    seen = set()
    residues = set()
    while (two, three) not in seen:
        seen.add((two, three))
        residues.add((constant + linear*two + c3*three) % modulus)
        two, three = 2*two % modulus, 3*three % modulus
    squares = {x*x % modulus for x in range(modulus)}
    assert residues.isdisjoint(squares)
    return sorted(residues)


def closest(coeff, target):
    """Lower-bound search in a strictly increasing infinite sequence of counts."""
    lo, hi = 1, target.bit_length() + 1
    i = bisect_left(range(lo, hi + 1), target, key=lambda w: height(coeff, w))
    above = lo + i
    return min((abs(height(coeff, w)-target), w, height(coeff, w)-target)
               for w in (above-1, above) if w >= 1)


if __name__ == "__main__":
    perfect = 2**60*(2**61-1)
    distances = []
    for i, (rail, modulus) in enumerate(zip(RAILS, MODULI), 1):
        coeff = rail[-1]
        start = 2 if i in (4,8) else 1
        for w in range(1,start):
            v = height(coeff,w)+2
            assert isqrt(v)**2 != v
        print(i, modulus, square_obstruction(coeff,modulus,start))
        distances.append((*closest(coeff,perfect), i))
    _, w, delta, i = min(distances)
    assert (i,w,delta)==(6,111,5192296858534828781452000936066879)
    assert height(RAILS[5][-1],3) == 8015
    print("Closest bank member to exponent 61: rail, width, signed offset", i, w, delta)
