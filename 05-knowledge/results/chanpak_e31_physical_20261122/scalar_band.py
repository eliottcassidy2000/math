"""Extend the exponent-31 band by ordinal products. Run from any directory."""

from math import gcd, prod
from check import p, perfect_shape
from e31_certificates import stored_pairs

# Offset, ordinary counts inside the certified prefix, optional two-apex address pair.
PRODUCTS = (
    (-41, (19, 3061), (7, 3814081016907)),
    (-39, (281, 3217, 363533, 7016629), None),
    (-37, (3, 3), (2, 124425867154225936)),
    (-35, (51343,), (2, 21806308967184)),
    (-33, (5, 10889, 888653, 47658407), None),
    (-31, (3, 107), (2, 2460064459210281)),
    (-29, (91, 269, 153871, 612180211), None),
    (-27, (191644001,), (9, 704382041)),
    (27, (35, 73387, 2035151, 441109), None),
    (29, (3, 59833), (2, 4568227477517)),
    (31, (23, 227), (11, 37610280057688)),
    (33, (16249, 20021), (2, 2389595409)),
    (35, (3, 3, 19, 19, 19, 1237), (2, 11360166379)),
    (37, (5, 1049, 53197), (35, 208587820)),
    (39, (197, 119293), (57, 919978849)),
)


def check():
    center = perfect_shape(31)
    endpoints = {n for lower, *_ in stored_pairs() for n in (lower, lower + 2)}
    for offset, factors, addresses in PRODUCTS:
        assert all(1 <= factor < 2**32 and factor % 2 and factor not in (7, 21)
                   for factor in factors)
        height = prod(factors)
        if addresses:
            b, e = addresses
            height *= b * e + (1 + p(b - 1)) * (1 + p(e - 1)) - (b & e & 1)
        assert height == center + offset
        endpoints.add(height)
    assert all(center + j in endpoints for j in range(-39, 40, 2))
    heights = list(range(center - 39, center + 40, 2))
    width = heights[-1] - heights[0]
    residuals = []
    for height in heights:
        residual = height
        for divisor in range(2, width + 1):
            while residual % divisor == 0:
                residual //= divisor
        assert residual > 1 and gcd(height // residual, residual) == 1
        assert all(gcd(other, residual) == 1 for other in heights if other != height)
        residuals.append(residual)
    assert min(zip(residuals, heights)) == (29861325211, center - 3)
    print('All 40 odd counts at offsets from -39 through 39 have private divisors and tournament witnesses.')


if __name__ == '__main__':
    check()
