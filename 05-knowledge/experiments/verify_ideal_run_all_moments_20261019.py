#!/usr/bin/env python3
"""Exact focused controls for the sharp shuffle formulas; standard library only."""
from fractions import Fraction
from itertools import product
from math import comb
import json

def compositions(k, n):
    if n == 1:
        yield (k,)
        return
    for a in range(k + 1):
        for tail in compositions(k - a, n - 1):
            yield ((a,) + tail)

def main():
    tested = 0
    for k, ell in product(range(9), range(6)):
        runs = [a[0] for a in compositions(k, ell + 1)]
        mean = Fraction(sum(runs), len(runs))
        if not mean == Fraction(k, ell + 1):
            raise AssertionError('exact certificate check failed')
        for r in range(7):
            moment = Fraction(sum((comb(x, r) if x >= r else 0 for x in runs)), len(runs))
            formula = Fraction(comb(k + ell, ell + r), comb(k + ell, ell)) if r <= k else Fraction(0)
            shifted = Fraction(sum((comb(x + r, r) for x in runs)), len(runs))
            if not (moment == formula and moment <= mean ** r):
                raise AssertionError('exact certificate check failed')
            if not shifted <= (1 + mean) ** r:
                raise AssertionError('exact certificate check failed')
            if k > 0 and r >= 2:
                if not moment < mean ** r:
                    raise AssertionError('exact certificate check failed')
            tested += 1
    return dict(status='PASS', scope='54 exact shuffle laws, moment orders 0 through 6', tests=tested, general_proof='Stanley rank log-concavity, tail submultiplicativity, compositions, Vandermonde')
if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
