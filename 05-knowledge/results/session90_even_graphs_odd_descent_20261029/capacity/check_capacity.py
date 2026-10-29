#!/usr/bin/env python3
"""Exact necessary capacity obstructions; assertions remain active under -O."""
import json
from math import comb, factorial

def require(p, message):
    if not p:
        raise AssertionError(message)

def main():
    checks = 0
    result = {}
    for k, expected in ((1, 13), (2, 28), (3, 43)):
        hits = []
        for n in range(2, 501):
            E = comb(n, 2)
            lhs = 1 << E
            rhs = factorial(n) * comb(E, n-1)**k
            if lhs > rhs:
                hits.append(n)
            checks += 1
        require(hits == list(range(expected, 501)), f'capacity range k={k}')
        checks += 1
        n = expected
        # A second integer route: iterative labelled sparse subset count.
        E = n*(n-1)//2
        b = 1
        for j in range(1, n):
            b = b*(E-j+1)//j
        require(b == comb(E, n-1), 'binomial route')
        require((1 << E) > factorial(n)*b**k, 'independent boundary')
        checks += 2
        result[str(k)] = {'first_bound_failure': expected,
                          'consecutive_checked_through': 500}
    # Finite cycle/component accounting across all possible component ranks
    # in a small box. A connected component is a tree iff beta_i=0.
    from itertools import product
    for c in range(1, 8):
        for ranks in product(range(5), repeat=c):
            lhs = sum(ranks) == c-1
            rhs = ranks.count(0) == 1+sum(r-1 for r in ranks if r)
            require(lhs == rhs, 'Euler component identity')
            checks += 1
            if lhs:
                for r in ranks:
                    require(ranks.count(0) >= r, 'core tree debt')
                    checks += 1
    print(json.dumps({'status': 'PASS', 'checks': checks,
                      'capacity': result, 'scope': 'loose necessary bounds, not exact class counts'},
                     sort_keys=True, indent=2))

if __name__ == '__main__':
    main()
