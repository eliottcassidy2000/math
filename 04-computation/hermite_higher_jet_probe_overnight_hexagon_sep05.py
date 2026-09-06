#!/usr/bin/env python3
"""Exact exploratory higher-Hasse-jet inverse-denominator metric probes.

The universe is all zero-translated node sets of requested sizes/diameter.
For jet multiplicity k, reciprocal Taylor coefficients of Q_i(x_i+T)^-k
give the exact inverse denominator. Comparing every integer dilation e>=0
uses the eventual maximum and each finite breakpoint, not just e=0.
No general metric-only theorem is assumed or inferred from a passing head.
"""

import argparse
from fractions import Fraction
import hashlib
import itertools
import math


def valuation(a, p):
    if not a:
        return math.inf
    if isinstance(a, Fraction):
        return valuation(a.numerator, p) - valuation(a.denominator, p)
    value = 0
    while a % p == 0:
        a //= p
        value += 1
    return value


def tree(nodes, p):
    if len(nodes) == 1:
        return (-1, ())
    depth = min(valuation(x - nodes[0], p) for x in nodes[1:])
    groups = {}
    for x in nodes:
        groups.setdefault(x % p ** (depth + 1), []).append(x)
    return depth, tuple(sorted(tree(group, p) for group in groups.values()))


def inverse_constants(nodes, k, p):
    constants = [-math.inf] * k
    for x in nodes:
        coefficients = [Fraction(1)] + [Fraction(0)] * (k - 1)
        base = 0
        for y in nodes:
            if x == y:
                continue
            d = x - y
            base += k * valuation(d, p)
            factor = [Fraction((-1) ** j * math.comb(k+j-1, j), d ** j)
                      for j in range(k)]
            coefficients = [sum(coefficients[j] * factor[l-j]
                                for j in range(l+1)) for l in range(k)]
        for l, coefficient in enumerate(coefficients):
            constants[l] = max(constants[l], base - valuation(coefficient, p))
    return tuple(constants)


def envelope(constants):
    # The common k(n-1)e term is omitted. Once the last finite line wins,
    # it wins forever; testing every integer through all intersections suffices.
    finite = [(l, c) for l, c in enumerate(constants) if math.isfinite(c)]
    last_l, last_c = finite[-1]
    terminal = max([0] + [max(0, math.ceil((c-last_c)/(last_l-l)))
                         for l, c in finite[:-1]])
    values = tuple(max(c+l*e for l, c in finite) for e in range(terminal+1))
    # Canonical representation strips an already linear terminal suffix.
    while len(values) > 1 and values[-1] - values[-2] == last_l:
        values = values[:-1]
    return values, last_l


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--height', type=int, default=16)
    parser.add_argument('--max-nodes', type=int, default=5)
    parser.add_argument('--jets', type=int, nargs='+', default=[3, 4, 5])
    parser.add_argument('--primes', type=int, nargs='+', default=[2, 3, 5, 7])
    args = parser.parse_args()
    digest = hashlib.sha256()
    rows = 0
    hostiles = 0
    for k, p, n in itertools.product(args.jets, args.primes,
                                     range(2, args.max_nodes+1)):
        seen = {}
        for rest in itertools.combinations(range(1, args.height+1), n-1):
            nodes = (0,) + rest
            key = tree(nodes, p)
            constants = inverse_constants(nodes, k, p)
            observed = envelope(constants)
            rows += 1
            digest.update(repr((k, p, nodes, constants, observed)).encode())
            if key in seen and seen[key][0] != observed:
                old, old_nodes, old_constants = seen[key]
                print('METRIC HOSTILE', 'k', k, 'p', p, 'n', n,
                      old_nodes, nodes, 'constants', old_constants, constants,
                      'envelopes', old, observed, flush=True)
                hostiles += 1
                break
            seen.setdefault(key, (observed, nodes, constants))
        else:
            print('HEAD PASS', 'k', k, 'p', p, 'n', n,
                  'metric types', len(seen), flush=True)
    print('ROWS', rows, 'HOSTILE SHAPES', hostiles,
          'SEMANTIC', digest.hexdigest())
    print('FINITE-EXACT exploratory inverse formula; literal matrix audit pending.')


if __name__ == '__main__':
    main()
