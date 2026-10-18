#!/usr/bin/env python3
"""Independent Virasoro PBW engine for the s44c audit (written from scratch; imports nothing from the session).

Vacuum module V^c = U(Vir_-)|0>, L_{-1}|0> = 0.  Basis: L_{-k_1} ... L_{-k_r}|0> with k_1 >= ... >= k_r >= 2, stored as
the decreasing tuple (k_1, ..., k_r).  Relations: [L_m, L_n] = (m-n) L_{m+n} + (c/12)(m^3-m) delta_{m+n,0}.

Coefficients: exact Fractions (P=None) or integers mod a prime P.  act(n, mono) returns a dict {mono: coeff} for
L_n applied to the basis vector mono, by left-to-right commutation (memoized).
"""
from fractions import Fraction


def partitions_ge2(W, maxpart=None):
    """All partitions of W into parts >= 2, as decreasing tuples."""
    if maxpart is None:
        maxpart = W
    if W == 0:
        yield ()
        return
    for k in range(min(W, maxpart), 1, -1):
        for rest in partitions_ge2(W - k, k):
            yield (k,) + rest


def partitions(n, maxpart=None, minpart=1):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), minpart - 1, -1):
        for rest in partitions(n - k, k, minpart):
            yield (k,) + rest


class VirVacuum:
    def __init__(self, c, P=None):
        self.P = P
        if P is None:
            self.c = Fraction(c)
            self.c12 = self.c / 12
        else:
            cc = Fraction(c)
            self.c = cc.numerator % P * pow(cc.denominator % P, -1, P) % P
            self.c12 = self.c * pow(12, -1, P) % P
        self.cache = {}

    def _red(self, x):
        return x % self.P if self.P else x

    def act(self, n, mono):
        key = (n, mono)
        r = self.cache.get(key)
        if r is None:
            r = self._act(n, mono)
            self.cache[key] = r
        return r

    def _act(self, n, mono):
        P = self.P
        if n == 0:
            w = sum(mono)
            return {mono: w} if w else {}
        if not mono:
            return {(-n,): 1} if n <= -2 else {}
        k = mono[0]
        if n <= -2 and -n >= k:
            return {(-n,) + mono: 1}
        rest = mono[1:]
        out = {}
        # L_n L_{-k} rest = L_{-k} (L_n rest) + (n+k) L_{n-k} rest + (c/12)(n^3-n) delta_{n,k} rest
        for m2, a in self.act(n, rest).items():
            for m3, b in self.act(-k, m2).items():
                out[m3] = out.get(m3, 0) + a * b
        coef = n + k
        if coef:
            for m2, a in self.act(n - k, rest).items():
                out[m2] = out.get(m2, 0) + coef * a
        if n == k:
            out[rest] = out.get(rest, 0) + self.c12 * (n ** 3 - n)
        if P:
            return {m: v % P for m, v in out.items() if v % P}
        return {m: v for m, v in out.items() if v}

    def apply(self, n, vec):
        out = {}
        for m, a in vec.items():
            for m2, b in self.act(n, m).items():
                out[m2] = out.get(m2, 0) + a * b
        if self.P:
            return {m: v % self.P for m, v in out.items() if v % self.P}
        return {m: v for m, v in out.items() if v}


def central_charge(p, q):
    return 1 - Fraction(6 * (p - q) ** 2, p * q)


if __name__ == '__main__':
    # smoke tests: hand-computed values
    V = VirVacuum(Fraction(7, 10))
    assert V.act(1, (3,)) == {(2,): 4}
    assert V.act(2, (4,)) == {(2,): 6}
    assert V.act(1, (2, 2)) == {(3,): 3}                         # L_1 L_-2 L_-2|0> = 3 L_-3|0>
    assert V.act(2, (2, 2)) == {(2,): 8 + Fraction(7, 10)}        # (8 + c) L_-2|0>
    assert V.act(2, (2,)) == {(): Fraction(7, 20)}                # c/2
    assert V.act(-1, (2,)) == {(3,): 1}
    assert V.act(-3, (2,)) == {(3, 2): 1}
    # L_{-2} L_{-3} = L_{-3} L_{-2} + L_{-5}
    assert V.act(-2, (3,)) == {(3, 2): 1, (5,): 1}
    r = V.act(-2, (3, 2))   # L_-2 L_-3 L_-2 = L_-3 L_-2 L_-2 + L_-5 L_-2
    assert r == {(3, 2, 2): 1, (5, 2): 1}, r
    print('vir_pbw smoke tests OK')
