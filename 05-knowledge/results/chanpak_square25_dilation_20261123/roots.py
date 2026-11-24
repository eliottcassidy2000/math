"""Integer obstructions and exact inverses for the twelve fixed order-eleven seeds."""
from bisect import bisect_left
from math import isqrt

from extensions import RAILS, height

MODULI = (4, 15, 4, 8, 4, 19, 4, 8, 4, 9, 4, 15)


def square_obstruction(coeff, modulus, start, offset=2):
    """Follow the joint powers until their state repeats, including the transient."""
    two, three = pow(2, start, modulus), pow(3, start, modulus)
    c1, c2 = coeff[:2]
    c3 = coeff[2] if len(coeff) == 3 else 0
    constant, linear = offset + c1 - 2*c2 + 3*c3, c2 - 3*c3
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


def any_local_square(modulus):
    forms = []
    for *_, coeff in RAILS:
        c1, c2 = coeff[:2]
        c3 = coeff[2] if len(coeff) == 3 else 0
        forms.append((2+c1-2*c2+3*c3, c2-3*c3, c3))
    squares = {j*j % modulus for j in range(modulus//2 + 1)}
    u, v = 2 % modulus, 3 % modulus
    seen = set()
    while (u, v) not in seen:
        if any((a+b*u+c*v) % modulus in squares for a,b,c in forms):
            return True
        seen.add((u, v))
        u, v = 2*u % modulus, 3*v % modulus
    return False


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
    assert next(m for m in range(2, 1369) if not any_local_square(m)) == 1368
    print("Smallest modulus rejecting every shifted height as a square: 1368 = 8*9*19")
    certificates = {
        0: (65, 4, 3, 4, 8, 3, 7, 4, 3, 4, 315, 4),
        4: (8, 4, 8, 4, 1353, 4, 255, 4, 41, 3, 8, 4),
    }
    square_hits = []
    for offset, moduli in certificates.items():
        for (rail, modulus) in zip(RAILS, moduli):
            coeff = rail[-1]
            square_obstruction(coeff, modulus, 5, offset)
            for width in range(1, 5):
                shifted = height(coeff, width) + offset
                if isqrt(shifted)**2 == shifted:
                    square_hits.append((offset, coeff[0] % 24, width, shifted))
    assert square_hits == [(0, 9, 1, 5625)]
    print("All square collisions at offsets 0 and 4:", square_hits)
