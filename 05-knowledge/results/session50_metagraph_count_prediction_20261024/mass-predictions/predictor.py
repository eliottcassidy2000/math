#!/usr/bin/env python3
"""Blind odd-fibre predictor: only a list of integer moments, never atlas access."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

CHECKS = 0

def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)

def rational(x):
    return str(F(x))

def ceil(x):
    x = F(x)
    return -(-x.numerator // x.denominator)

def floor(x):
    x = F(x)
    return x.numerator // x.denominator

def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out

def solve(matrix, vector):
    """Exact Gaussian elimination; singular systems are not silently inverted."""
    n = len(vector)
    a = [[F(x) for x in row] + [F(y)] for row, y in zip(matrix, vector)]
    require(len(a) == n and all(len(row) == n+1 for row in a), 'matrix shape')
    for i in range(n):
        pivot = next((j for j in range(i, n) if a[j][i]), None)
        require(pivot is not None, 'singular moment matrix')
        a[i], a[pivot] = a[pivot], a[i]
        t = a[i][i]
        a[i] = [x/t for x in a[i]]
        for j in range(i+1, n):
            t = a[j][i]
            if t:
                a[j] = [x-t*y for x, y in zip(a[j], a[i])]
    x = [F(0)]*n
    for i in range(n-1, -1, -1):
        x[i] = a[i][-1] - sum(a[i][j]*x[j] for j in range(i+1, n))
    for row, y in zip(matrix, vector):
        require(sum(t*v for t, v in zip(row, x)) == y, 'linear solve residual')
    return x

def odd_polynomial(k):
    q = [1]
    for j in range(k):
        q = mul(q, [-(2*j+1), 1])
    return q

def congruence(moments):
    """Only moments K_0,...,K_(k-1). No target, histogram, or support."""
    k = len(moments)
    q = odd_polynomial(k)
    exponent = 2*k-k.bit_count()
    modulus = 1 << exponent
    require(q[0] % 2 == 1, 'constant coefficient must be a unit')
    inverse = pow(q[0], -1, modulus)
    residue = (-inverse * sum(q[j]*moments[j-1] for j in range(1, k+1))) % modulus
    return {'k': k, 'exponent': exponent, 'modulus': modulus, 'residue': residue,
            'polynomial_coefficients': q, 'constant_inverse_mod_modulus': inverse}

def moment_bounds(moments, a=1):
    """Continuous positive-support bounds. NO upper endpoint is an input."""
    k = len(moments)
    total = moments[0]
    if not total:
        return F(0), F(0), [], []
    lower = F(1)  # a nonempty odd-fibre partition has at least one class
    upper = F(total, a)
    low_q, up_q = [], []
    # Singular Gram matrices can occur when support is already exhausted.
    # Use the largest invertible leading block; no determinant tolerance.
    for r in range(1, k//2+1):
        b = moments[:r]
        matrix = [[moments[i+j+1] for j in range(r)] for i in range(r)]
        try:
            q = solve(matrix, b)
        except ValueError as error:
            if str(error) == 'singular moment matrix':
                break
            raise
        value = sum(x*y for x, y in zip(b, q))
        require(value >= lower, 'lower-bound monotonicity')
        lower, low_q = value, q
    for r in range(1, (k-1)//2+1):
        b = [moments[i+1]-a*moments[i] for i in range(r)]
        matrix = [[moments[i+j+2]-a*moments[i+j+1] for j in range(r)] for i in range(r)]
        try:
            q = solve(matrix, b)
        except ValueError as error:
            if str(error) == 'singular moment matrix':
                break
            raise
        value = (F(total)-sum(x*y for x, y in zip(b, q)))/a
        require(value <= upper, 'upper-bound monotonicity')
        upper, up_q = value, q
    require(lower <= upper, 'inverted bounds')
    return lower, upper, low_q, up_q

def predict(moments, peeled_count=0):
    """Blind interface: K vector and optional independently known h=1 count only."""
    congr = congruence(moments)
    residual = [x-peeled_count for x in moments]
    lower, upper, low_q, up_q = moment_bounds(residual, a=3 if peeled_count else 1)
    lower += peeled_count
    upper += peeled_count
    first = ceil(lower)
    first += (congr['residue']-first) % congr['modulus']
    last = floor(upper)
    count = max(0, 1+(last-first)//congr['modulus'])
    return {'moment_count': len(moments), 'max_word_moment_degree': len(moments)-1,
            'peeled_h1_count': peeled_count, 'lower': rational(lower), 'upper': rational(upper),
            'lower_decimal': round(float(lower), 12), 'upper_decimal': round(float(upper), 12),
            'lower_q_coefficients': list(map(rational, low_q)),
            'upper_q_coefficients': list(map(rational, up_q)),
            'congruence': congr, 'first_candidate': first if count else None,
            'candidate_count': count, 'prediction': first if count == 1 else None}

def predict_modular(universe_size, residues):
    """Only N and K residues; universal bound 1 <= A <= N, no support atlas."""
    result = congruence(residues)
    modulus = result['modulus']
    require(type(universe_size) is int and universe_size > 0, 'positive universe size')
    require(all(type(x) is int and 0 <= x < modulus for x in residues), 'canonical moment residues')
    require(residues[0] == universe_size % modulus, 'K0 must equal universe cardinality modulo M')
    first = result['residue'] or modulus
    count = max(0, 1+(universe_size-first)//modulus)
    return {'universe_size': universe_size, 'congruence': result,
            'first_candidate': first if count else None, 'candidate_count': count,
            'prediction': first if count == 1 else None}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--moments', required=True, help='JSON list K_0,...,K_(k-1)')
    parser.add_argument('--peeled-h1-count', type=int, default=0)
    parser.add_argument('--output', default='blind_prediction.json')
    parser.add_argument('--modular-universe-size', type=int)
    args = parser.parse_args()
    moments = json.loads(Path(args.moments).read_text())
    require(isinstance(moments, list) and moments and all(type(x) is int and x >= 0 for x in moments), 'integer moment list required')
    answer = predict_modular(args.modular_universe_size, moments) if args.modular_universe_size else predict(moments, args.peeled_h1_count)
    Path(args.output).write_text(json.dumps(answer, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'moment_count': len(moments), 'prediction': answer['prediction'], 'candidate_count': answer['candidate_count'], 'modulus': answer['congruence']['modulus']}, sort_keys=True))
