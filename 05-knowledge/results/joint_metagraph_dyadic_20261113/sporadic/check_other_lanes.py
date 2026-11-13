#!/usr/bin/env python3
"""Independent arithmetic probes for pairings P1/P2; no producer imports."""
from math import comb
from pathlib import Path
from hashlib import sha256
import json

gates = 0


def require(condition, label):
    global gates
    gates += 1
    if not condition:
        raise RuntimeError(label)


def valuation_two(x):
    return (x & -x).bit_length()-1


def divide_by_one_plus_x(mask):
    quotient = 0
    remainder = mask
    while remainder.bit_length() > 1:
        shift = remainder.bit_length()-2
        quotient ^= 1 << shift
        remainder ^= 3 << shift
    return quotient, remainder


def main():
    for R in range(1, 49):
        for k in range(1, 6):
            L = k + (R.bit_length()-1)
            P = 2**L
            for r in range(R+1):
                for upper in range(r, r+41):
                    require((comb(upper+P, r)-comb(upper, r)) % (2**k) == 0,
                            'ordinary binomial period')
            for t in range(R//2+1):
                for upper in range(t, t+41):
                    require((comb(upper+P//2, t)-comb(upper, t)) % (2**k) == 0,
                            'converse half-order period')
    synthetic_cases = 0
    for d in range(1, 9):
        for numerator in range(1 << (d+2)):
            if numerator == 0:
                exponent = 0
            else:
                reduced = numerator
                cancelled = 0
                while reduced.bit_count() % 2 == 0:
                    reduced, rem = divide_by_one_plus_x(reduced)
                    require(rem == 0, 'polynomial factor division')
                    cancelled += 1
                exponent = max(0, d-cancelled)
            predicted = 1 << ((max(1, exponent)-1).bit_length())
            coefficients = []
            for n in range(100):
                a = sum(comb(n-j+d-1, d-1) for j in range(numerator.bit_length())
                        if ((numerator >> j) & 1) and j <= n) % 2
                coefficients.append(a)
            onset = d+3
            periods = [p for p in range(1, 17)
                       if all(coefficients[n] == coefficients[n+p]
                              for n in range(onset, 100-p))]
            require(min(periods) == predicted, 'least eventual parity period')
            synthetic_cases += 1
    w = list(range(25)) + [70]
    z = [1, 0] + list(range(2, 25)) + [70]
    dot = lambda a,b: sum(x*y for x,y in zip(a[:-1],b[:-1]))-a[-1]*b[-1]
    require(dot(w,w) == dot(z,z) == 0 and dot(w,z) == -1, 'Lorentzian pair')
    for i in range(26):
        x = [0]*26
        x[i] = 2
        projection = [a+dot(x,z)*b+dot(x,w)*c for a,b,c in zip(x,w,z)]
        require(dot(projection,w) == dot(projection,z) == 0, 'orthogonal projection')
    root = Path(__file__).resolve().parent.parent
    targets = ['audit/check_small.py','pairings/PROOF.md']
    print(json.dumps({'status':'PASS','active_checks':gates,
                      'P1_universe':'R=1..48, k=1..5, all r<=R and t<=floor(R/2), 41 nonnegative upper shifts each',
                      'P2_universe':'all binary numerators of degree<=d+1, denominator (1+x)^d, d=1..8; 100 coefficients; all tested periods1..16',
                      'P2_synthetic_cases':synthetic_cases,
                      'source_sha256':{p:sha256((root/p).read_bytes()).hexdigest() for p in targets},
                      'Lorentzian_Gram':[[dot(w,w),dot(w,z)],[dot(z,w),dot(z,z)]]},
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
