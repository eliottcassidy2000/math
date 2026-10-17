#!/usr/bin/env python3
"""S39 audit, item 3 (statement 2 for odd prime powers, including characteristic 3) and a census check of the proof of
statement 3 (inert / principal / nonprincipal primes for Q(sqrt-15)).

G_q = Cay(F_(q^2), mu_(q+1)), q = p^m, F_(q^2) = F_p[x]/(f), deg f = 2m.  Exact: W_2, W_3, W_4 at the vertex 0 and the
codegree lemma R_2(d) = 1 - eta_q(N(d)(N(d) - 4)) for every d != 0 (eta_q the quadratic character of F_q, evaluated as
x^((q-1)/2) inside F_(q^2)).  q in {9, 25, 27, 49, 121, 125}.
Usage: python -B a10_statement2_prime_powers.py
"""
import numpy as np
import sympy
from sympy import primerange


def irreducible(p, n):
    X = sympy.symbols('X')
    import itertools
    for coeffs in itertools.product(range(p), repeat=n):
        if coeffs[0] == 0:
            continue
        f = sympy.Poly([1] + list(coeffs[::-1]), X, modulus=p)
        if f.is_irreducible:
            return [int(c) % p for c in f.all_coeffs()[::-1]]       # low -> high, monic
    raise RuntimeError


def make_field(p, n):
    f = irreducible(p, n)

    def mul(u, v):
        r = [0] * (2 * n - 1)
        for i in range(n):
            if u[i]:
                for j in range(n):
                    r[i + j] += u[i] * v[j]
        for k in range(2 * n - 2, n - 1, -1):
            if r[k]:
                co = r[k]
                for i in range(n + 1):
                    r[k - n + i] -= co * f[i]
        return tuple(x % p for x in r[:n])

    def power(u, e):
        res = tuple([1] + [0] * (n - 1)); base = u
        while e:
            if e & 1:
                res = mul(res, base)
            base = mul(base, base); e >>= 1
        return res
    return mul, power


for (p, m) in ((3, 2), (5, 2), (3, 3), (7, 2), (11, 2), (5, 3)):
    q = p ** m; n = 2 * m
    mul, power = make_field(p, n)
    N_el = p ** n
    els = [tuple((i // p ** k) % p for k in range(n)) for i in range(N_el)]
    one = tuple([1] + [0] * (n - 1)); mone = tuple([p - 1] + [0] * (n - 1)); zero = tuple([0] * n)
    norms = [power(e, q + 1) for e in els]
    mu = [i for i in range(N_el) if norms[i] == one]
    assert len(mu) == q + 1
    coords = np.array(els, dtype=np.int64)
    w = np.array([p ** k for k in range(n)], dtype=np.int64)
    shifts = [((coords - coords[u]) % p) @ w for u in mu]
    R = np.zeros(N_el, dtype=np.int64); R[0] = 1
    W = []
    for k in range(1, 5):
        Rn = np.zeros(N_el, dtype=np.int64)
        for s in shifts:
            Rn += R[s]
        if k == 2:
            R2 = Rn.copy()
        R = Rn
        W.append(int(R[0]))
    e3 = 0 if p == 3 else (1 if q % 3 == 1 else -1)
    okW = (W[1] == q + 1, W[2] == (q + 1) * (1 - e3), W[3] == 3 * q * (q + 1))

    def eta(xel):
        if xel == zero:
            return 0
        t = power(xel, (q - 1) // 2)
        return 1 if t == one else (-1 if t == mone else None)

    def sub(u, v):
        return tuple((a - b) % p for a, b in zip(u, v))
    four = tuple([4 % p] + [0] * (n - 1))
    bad = 0
    for i in range(1, N_el):
        c = norms[i]
        val = 1 - eta(mul(c, sub(c, four)))
        if R2[i] != val:
            bad += 1
    print('q = %3d (p = %d): W_2, W_3, W_4 as stated: %s; codegree lemma failures over all d != 0: %d' % (q, p, okW, bad))
    assert all(okW) and bad == 0

# census of the session's fitting claim: inert / principal / nonprincipal split primes 7 <= p <= 2999 for Q(sqrt-15)
inert = princ = nonp = 0
for p in primerange(7, 3000):
    if sympy.jacobi_symbol(-15 % p, p) == -1:
        inert += 1
        continue
    r = int((4 * p) ** 0.5) + 2
    if any(x * x + x * y + 4 * y * y == p for y in range(0, r) for x in range(-r, r)):
        princ += 1
    else:
        nonp += 1
print('primes 7 <= p <= 2999: inert %d, split principal %d, split nonprincipal %d (session: 218, 98, 111)' % (inert, princ, nonp))
