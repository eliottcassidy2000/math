#!/usr/bin/env python3
"""Exact symbolic shuffle certificate for the all-common-neighbor two-slave bound.
Standard library only. Polynomial identities are derived, not interpolated.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations
from math import factorial
from pathlib import Path
import json

class Poly:

    def __init__(self, value=0):
        self.d = dict(value) if isinstance(value, dict) else {(0, 0, 0): F(value)} if value else {}
        self.d = {k: F(v) for k, v in self.d.items() if v}

    def __add__(self, other):
        other = other if isinstance(other, Poly) else Poly(other)
        d = self.d.copy()
        for k, v in other.d.items():
            d[k] = d.get(k, F(0)) + v
        return Poly(d)
    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.d.items()})

    def __sub__(self, other):
        return self + -aspoly(other)

    def __rsub__(self, other):
        return aspoly(other) + -self

    def __mul__(self, other):
        other = aspoly(other)
        d = {}
        for a, u in self.d.items():
            for b, v in other.d.items():
                k = tuple((x + y for x, y in zip(a, b)))
                d[k] = d.get(k, F(0)) + u * v
        return Poly(d)
    __rmul__ = __mul__

    def __truediv__(self, n):
        return self * F(1, n)

    def __pow__(self, n):
        out = Poly(1)
        for _ in range(n):
            out = out * self
        return out

    def __eq__(self, other):
        return self.d == aspoly(other).d

    def at(self, p, q, r):
        return sum((v * p ** i * q ** j * r ** k for (i, j, k), v in self.d.items()))

    def record(self):
        return [[*k, str(v)] for k, v in sorted(self.d.items())]

def aspoly(v):
    return v if isinstance(v, Poly) else Poly(v)
p, q, r = [Poly({tuple((int(i == j) for i in range(3))): 1}) for j in range(3)]

def choose(n, k):
    ans = Poly(1)
    for i in range(k):
        ans = ans * (n - i)
    return ans / factorial(k)

def derive(region, ba):
    E = D = Poly(0)
    summands = 0
    for order in permutations('bxya'):
        pos = {v: order.index(v) for v in order}
        if pos['b'] > pos['y'] or pos['x'] > pos['a'] or (ba and pos['b'] > pos['a']):
            continue
        for h1 in range(5):
            for h2 in range(5 - h1):
                h3 = 4 - h1 - h2
                ranks = {}
                if region == 'gap':
                    if pos['a'] < h1:
                        continue
                    force = pos['b'] - (h1 + h2) + 1 if pos['b'] >= h1 + h2 else 0
                    f3 = h3 - force
                    weight = choose(p + h1 - 1, h1) * choose(q + h2 - 1, h2) * choose(r + f3, f3)
                    for j, v in enumerate(order):
                        if j < h1:
                            ranks[v] = (p - 1) * F(j + 1, h1 + 1)
                        elif j < h1 + h2:
                            ranks[v] = p + (q - 1) * F(j - h1 + 1, h2 + 1)
                        elif j < h1 + h2 + force:
                            ranks[v] = p + q
                        else:
                            ranks[v] = p + q + r * F(j - h1 - h2 - force + 1, f3 + 1)
                    c = p + q + r
                else:
                    if pos['a'] < h1 + h2 or pos['b'] >= h1 + h2:
                        continue
                    force = pos['b'] - h1 + 1 if pos['b'] >= h1 else 0
                    f2 = h2 - force
                    weight = choose(p + h1 - 1, h1) * choose(q + f2, f2) * choose(r + h3, h3)
                    for j, v in enumerate(order):
                        if j < h1:
                            ranks[v] = (p - 1) * F(j + 1, h1 + 1)
                        elif j < h1 + force:
                            ranks[v] = p
                        elif j < h1 + h2:
                            ranks[v] = p + q * F(j - h1 - force + 1, f2 + 1)
                        else:
                            ranks[v] = p + q + 1 + r * F(j - h1 - h2 + 1, h3 + 1)
                    c = p + q + r + 1
                E = E + weight
                D = D + weight * (c + 3 * (ranks['x'] - ranks['y']))
                summands += 1
    return (E, D, summands)

def printed_formulas():
    A0 = (p + 1) * (r + 1) * (p ** 3 + r ** 3 + 2 * p * r * (p + r) + 3 * (p ** 2 + r ** 2) + 8 * p * r)
    A1 = p ** 4 + r ** 4 + 4 * p * r * (p ** 2 + r ** 2) + 4 * p ** 2 * r ** 2 + 6 * (p ** 3 + r ** 3) + 17 * p * r * (p + r) + 10 * (p ** 2 + r ** 2) + 22 * p * r + 5 * (p + r)
    A2 = 2 * (p ** 3 + r ** 3) + 3 * p * r * (p + r) + 8 * (p ** 2 + r ** 2) + 6 * p * r + 6 * (p + r)
    A3 = p ** 2 + r ** 2 + p + r
    D0 = (A0 + q * A1 + q ** 2 * A2 + q ** 3 * A3) / 4
    H = 2 * q ** 3 + (10 * (p + r) + 14) * q ** 2 + (5 * (p ** 2 + r ** 2) + 40 * p * r + 35 * (p + r) + 24) * q + 15 * (p + 1) * (p + r) * (r + 1)
    D1 = D0 + (q + 1) * (q + 2) * H / 60
    J = p ** 3 + r ** 3 + 4 * q * (p ** 2 + r ** 2) + 2 * p * r * (p + r) + 7 * (p ** 2 + r ** 2) + 7 * q ** 2 * (p + r) + 8 * p * q * r + 25 * q * (p + r) + 16 * p * r + 20 * (p + r) + 4 * q ** 3 + 20 * q ** 2 + 24 * q + 8
    Do = (p + 1) * (r + 1) * J / 4
    return (D0, D1, Do)

def counts(c, K, L, ba):

    @lru_cache(None)
    def visit(j, mask):
        if j == c and mask == 15:
            return (1, 0)
        E = N = 0
        if j < c and (j < K or mask & 1):
            e, n = visit(j + 1, mask)
            E += e
            N += n + e * (int(bool(mask & 2)) + int(not bool(mask & 4)))
        for v in range(4):
            if mask >> v & 1:
                continue
            if v == 2 and (not mask & 1):
                continue
            if v == 3 and (not mask & 2 or j < L or (ba and (not mask & 1))):
                continue
            e, n = visit(j, mask | 1 << v)
            E += e
            N += n
        return (E, N)
    return visit(0, 0)

def main():
    formulas = printed_formulas()
    records = {}
    derived = {}
    for key, reg, ba, want in [('gap_no_ba', 'gap', False, formulas[0]), ('gap_with_ba', 'gap', True, formulas[1]), ('overlap', 'overlap', True, formulas[2])]:
        E, D, terms = derive(reg, ba)
        if not D == want:
            raise AssertionError('exact certificate check failed')
        if not all((v >= 0 for v in D.d.values())):
            raise AssertionError('exact certificate check failed')
        if not (max(map(sum, E.d)) <= 4 and max(map(sum, D.d)) <= 5):
            raise AssertionError('exact certificate check failed')
        derived[key] = (E, D)
        records[key] = {'finite_shuffle_summands': terms, 'E': E.record(), 'D': D.record(), 'nonzero_positive_D_coefficients': len(D.d)}
    checks = 0
    for c in range(31):
        for K in range(c + 1):
            for L in range(c + 1):
                for ba in (False, True):
                    E, N = counts(c, K, L, ba)
                    if K < L:
                        key = 'overlap'
                        args = (K, L - K - 1, c - L)
                    else:
                        key = 'gap_with_ba' if ba else 'gap_no_ba'
                        args = (L, K - L, c - K)
                    pe, pd = derived[key]
                    if not (pe.at(*args) == E and pd.at(*args) == 4 * c * E - 3 * N >= 0):
                        raise AssertionError('exact certificate check failed')
                    checks += 1
    missing = 0
    for c in range(61):
        for L in range(c + 1):
            E = (c + 3) * sum((i + 1 for i in range(L, c + 1)))
            N = F(c, 2) * sum(((i + 1) * (3 * c + 7 - i) for i in range(L, c + 1)))
            if not 7 * c * E - 6 * N == c * L * (L + 1) * (c - L + 1):
                raise AssertionError('exact certificate check failed')
            missing += 1
    for C in range(1, 20):
        for D in range(1, 20):
            for ba in (False, True):
                E = 2 * C * D + C + D + (1 if ba else 2)
                if not (2 * (2 * C * D + C) >= E and 2 * (2 * C * D + D) >= E):
                    raise AssertionError('exact certificate check failed')
    out = {'theorem': 'all-c mixed-two-slave aggregate, symbolic certificate', 'polynomials': records, 'independent_chain_DP_cases': checks, 'missing_slave_identity_cases': missing, 'status': 'PASS'}
    target = Path(__file__).with_name('two_slave_all_common_audit_20261019.json')
    print(json.dumps(out, indent=2, sort_keys=True))
if __name__ == '__main__':
    main()
