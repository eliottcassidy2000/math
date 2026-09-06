#!/usr/bin/env python3
"""Exact exploratory minors; no sampled regularity is an all-node theorem.

Universe: two Hasse jets at three distinct integer nodes, degree <6.
After translating first node to zero its unit block leaves a 4x4 matrix.
This probe tests whether the dyadic Smith partition is determined by the
three pairwise valuations. Every one of the 69 symbolic residual minors is
retained, not only a predicted optimal list.
"""

from collections import defaultdict
from itertools import combinations
from math import gcd
import argparse
import hashlib
import json
import sympy as s


def require(test, message):
    if not test:
        raise RuntimeError(message)


def valuation(n, p):
    n = abs(int(n))
    if n == 0:
        return 10**9
    result = 0
    while n % p == 0:
        result += 1
        n //= p
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--height', type=int, default=48)
    parser.add_argument('--show-minors', action='store_true')
    args = parser.parse_args()
    a, b = s.symbols('a b')
    matrix = s.Matrix([[x**q if r == 0 else q*x**(q-1)
                        for q in range(2, 6)]
                       for x in (a, b) for r in (0, 1)])
    minors = {}
    for h in range(1, 5):
        minors[h] = []
        for rr in combinations(range(4), h):
            for cc in combinations(range(4), h):
                expression = s.factor(matrix.extract(rr, cc).det())
                evaluator = s.lambdify((a,b), expression, modules='math')
                minors[h].append((rr, cc, expression, evaluator))
                if args.show_minors:
                    print('MINOR', h, rr, cc, expression)
    require(sum(map(len, minors.values())) == 69, 'minor universe')
    require(s.expand(minors[4][0][2] - a**4*b**4*(a-b)**4) == 0,
            'confluent determinant')
    census = []
    gates = 2
    global_cases = 0
    for x in range(1, args.height):
        for y in range(x+1, args.height+1):
            g = gcd(x,y)
            product = x*y*(y-x)
            expected = (g*gcd(g,2), gcd(g**4,6*product),
                        2*g**4*product*gcd(gcd(3,g),(x+y)//g), product**4)
            actual = tuple(gcd(*(int(f(x,y)) for _,_,_,f in minors[h]))
                           for h in range(1,5))
            require(expected == actual, ('global divisor formula',x,y,expected,actual))
            global_cases += 1
            gates += 1
    for p in (2,3,5,7):
        profiles = defaultdict(dict)
        for x in range(1, args.height):
            for y in range(x+1, args.height+1):
                d = [0]
                for h in range(1,5):
                    vals = [valuation(f(x,y), p) for _,_,_,f in minors[h]]
                    d.append(min(vals))
                exponents = (0,0)+tuple(d[i]-d[i-1] for i in range(1,5))
                distances = tuple(sorted((valuation(x,p), valuation(y,p),
                                          valuation(y-x,p))))
                require(exponents == tuple(sorted(exponents)), 'Smith order')
                require(sum(exponents) == 4*sum(distances), 'determinant budget')
                profiles[distances].setdefault(exponents, (0,x,y))
                census.append((p,x,y,distances,exponents))
                gates += 2
        twins = [(dist, list(types.items())) for dist, types in profiles.items()
                 if len(types) > 1]
        print('PRIME', p, 'distance_types', len(profiles), 'metric_twins', twins[:5])
        for dist, types in sorted(profiles.items()):
            if p == 2 or dist[0] > 0:
                print('PROFILE',p,dist,list(types.items()))
    print('UNIVERSE', 'nodes=(0,a,b),1<=a<b<=', args.height,
          'primes=2,3,5,7','rows',len(census),'global_cases',global_cases,
          'explicit_gates',gates)
    print('SEMANTIC_SHA256', hashlib.sha256(json.dumps(census).encode()).hexdigest())


if __name__ == '__main__':
    main()
