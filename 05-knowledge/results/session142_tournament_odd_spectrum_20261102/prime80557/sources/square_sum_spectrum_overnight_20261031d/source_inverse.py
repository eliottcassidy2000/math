#!/usr/bin/env python3
"""Exact inversion of a fixed singleton source, with no suffix-length cutoff."""
from math import gcd


def f(y):
    value = 1
    while y:
        a = y // 2
        if not y % 2:
            value += a + a % 2
        y = a
    return value


class Inverse:
    def __init__(self):
        self.nodes = 0
        self.direct = 0

    def solve(self, target, q, g, h):
        """Return y with target=q*y+g*(y mod2)+h*f(y), or None.

        q,g>=0,h>=1. In the q=0 case this consumer requires g=0;
        its repeated odd branch is removed by choosing an even representative.
        """
        if min(q, g) < 0 or h < 1 or (q == 0 and g != 0):
            raise ValueError('invalid source type')
        return self._solve(target, q, g, h)

    def _solve(self, t, q, g, h):
        # Explicit DFS preserves the binary address without a language-stack cap.
        pending = [(t, q, g, 0, 1)]
        while pending:
            t, q, g, address, place = pending.pop()
            self.nodes += 1
            if t < h:
                continue
            if t == h:
                return address
            if q == 0:
                pending.append((t, h, h, address, 2*place))
                continue
            if q > t-h:
                continue
            modulus = gcd(2*q, 2*h)
            if (t-h) % modulus and (t-h-q-g) % modulus:
                continue
            lo = max(1, (t-h-g+q+h-1)//(q+h))
            hi = (t-h)//q
            if hi-lo < 8:
                for y in range(lo, hi+1):
                    self.direct += 1
                    if q*y+g*(y%2)+h*f(y) == t:
                        return address+place*y
                continue
            pending.append((t, 2*q+h, h, address, 2*place))
            pending.append((t-q-g, 2*q, 0, address+place, 2*place))
        return None


if __name__ == '__main__':
    import json
    checks = 0
    solver = Inverse()
    for h in (1, 5, 25):
        for q in (1, 2, 5, 13, 32, 127):
            for g in (0, h//2, h):
                image = {q*y+g*(y%2)+h*f(y) for y in range(1000)}
                for t in range(h, min(1500, q*999+h)):
                    y = solver.solve(t, q, g, h)
                    # q*y+h<=t bounds every possible y by999 in this universe.
                    if (y is not None) != (t in image):
                        raise RuntimeError((t,q,g,h,y))
                    if y is not None and q*y+g*(y%2)+h*f(y) != t:
                        raise RuntimeError('bad witness')
                    checks += 1
    for h in (1,5,25):
        image = {h*f(y) for y in range(3001)}
        for t in range(h,1501):
            y = solver.solve(t,0,0,h)
            if (y is not None) != (t in image):
                raise RuntimeError(('zero axis',t,h,y))
            checks += 1
    print(json.dumps({'status':'passed','exhaustive_source_queries':checks,
                      'recursive_nodes':solver.nodes,'direct_evaluations':solver.direct},sort_keys=True))
