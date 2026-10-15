#!/usr/bin/env python3
"""McCarthy's normalised 3F2(chi4, phi, phi; eps, eps | 1)^* over F_Q, Q = p^r = 1 (mod 4) (helper for k6p6_descent_check.py).

3F2^*(A0,A1,A2; eps,eps | x) = (1/(Q-1)) sum_chi prod_i g(A_i chi)/g(A_i) * g(chibar)^3 * chi(-1) chi(x)
(McCarthy, FFA 18 (2012); with B_j = eps and g(eps) = -1).  Gauss sums g(chi) = sum_x chi(x) theta(Tr x), complex floats.
McCarthy-Papanikolas Thm 1.8: for primes p = 1 mod 4 this equals b(p), g = 32.3.c.a."""
import cmath, itertools, sys
import numpy as np
from sympy.polys.galoistools import gf_irreducible_p
from sympy.polys.domains import ZZ
import sympy

B = {3: 4j, 5: 2, 7: -8j, 11: -4j, 13: -14, 17: 18, 19: -12j, 23: 40j, 29: -14, 31: -32j, 37: -30, 41: -14}


def field(p, r):
    Q = p ** r
    if r == 1:
        f = [1, 0]
    else:
        for tail in itertools.product(range(p), repeat=r):
            f = [1] + list(tail)
            if gf_irreducible_p(f, p, ZZ):
                break
    def mul(a, b):
        prod = [0] * (2 * r - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                prod[i + j] = (prod[i + j] + x * y) % p
        for k in range(2 * r - 2, r - 1, -1):
            c = prod[k]
            if c:
                prod[k] = 0
                for j, fc in enumerate(f[1:]):
                    prod[k - 1 - j] = (prod[k - 1 - j] - c * fc) % p
        return prod[:r]
    one = [1] + [0] * (r - 1)
    def pw(a, e):
        res, b = one[:], a[:]
        while e:
            if e & 1: res = mul(res, b)
            b = mul(b, b); e >>= 1
        return res
    fac = sympy.factorint(Q - 1)
    for cand in itertools.product(range(p), repeat=r):
        gen = list(cand)
        if any(gen) and all(pw(gen, (Q - 1) // l) != one for l in fac):
            break
    enc = lambda v: sum(c * p ** k for k, c in enumerate(v))
    log = {}
    cur = one[:]
    elems = []
    for k in range(Q - 1):
        log[enc(cur)] = k
        elems.append(cur[:])
        cur = mul(cur, gen)
    # trace to F_p: Tr(x) = x + x^p + ... ; compute via powering
    def tr(v):
        t = [0] * r
        y = v[:]
        for _ in range(r):
            t = [(a + b) % p for a, b in zip(t, y)]
            y = pw(y, p)
        assert all(c == 0 for c in t[1:]), t
        return t[0]
    trace = [tr(v) for v in elems]          # trace of gen^k
    minus_one_log = log[enc([(p - 1)] + [0] * (r - 1))]
    return Q, trace, minus_one_log


def F32_star(p, r):
    Q, trace, l_m1 = field(p, r)
    n = Q - 1
    zeta = cmath.exp(2j * cmath.pi / p)
    th = np.array([zeta ** t for t in trace])               # theta(gen^k)
    # Gauss sum of character chi_a (chi_a(gen^k) = e^{2 pi i a k / n}): g(a) = sum_k chi_a(gen^k) theta(gen^k)
    k = np.arange(n)
    G = np.array([np.sum(np.exp(2j * np.pi * a * k / n) * th) for a in range(n)])
    chi4 = n // 4                                          # a character of order 4
    phi = n // 2
    tot = 0
    for a in range(n):
        A0, A1 = (chi4 + a) % n, (phi + a) % n
        term = (G[A0] / G[chi4]) * (G[A1] / G[phi]) ** 2 * G[(-a) % n] ** 3 * cmath.exp(2j * cmath.pi * a * l_m1 / n)
        tot += term
    return tot / n


