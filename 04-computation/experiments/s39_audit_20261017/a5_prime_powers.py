#!/usr/bin/env python3
"""S39 audit, item 5a: THM-4663 statement 6 (q = 49, 121, 169) by exact closed-walk counts.

G_q = Cay(F_(q^2), mu_(q+1)) with q = p^2: F_(q^2) = F_(p^4) = F_p[x]/(f), f an irreducible quartic; mu_(q+1) = {z : z^(q+1) = 1}.
W_k = closed k-walks at 0 (integer DP over the additive group (Z/p)^4), M_k(q) = (-1)^k (q^2 W_k - (q+1)^k)/(q+1).
Compared with the session's printed moments and with statement 3's forms at q, using the Frobenius-squared traces
a_p^2 - 2(-15/p) p^2, b_p^2 - 2p^3, d_p^2 - 2p^5 (forms from a2).  (The session used float Kloosterman sums over F_q.)
Usage: python -B a5_prime_powers.py
"""
import json
import numpy as np
import sympy

F = json.load(open('a2_forms_to_3000.json'))
A15, B46, D66 = F['f15'], F['b46'], F['f66']
SESSION = {49: [1, 2351, 2500, 227947, 718096, 28658090, 151801119, 4110000243],
           121: [1, 14519, 14884, 3498835, 10703056, 1090812458, 6241415104, 391460208171],
           169: [1, 28391, 28900, 9567427, 29104336, 4114382570, 23205956320, 2009849307387]}


def irreducible_quartic(p):
    x = sympy.symbols('x')
    for a in range(p):
        for b in range(p):
            for c in range(1, p):
                f = sympy.Poly(x ** 4 + a * x ** 2 + b * x + c, x, modulus=p)
                if f.is_irreducible:
                    return [c, b, a, 0, 1]      # coefficients of x^0..x^4
    raise RuntimeError


def polmul(u, v, f, p):
    r = [0] * 7
    for i in range(4):
        if u[i]:
            for j in range(4):
                r[i + j] += u[i] * v[j]
    for k in range(6, 3, -1):          # reduce x^k, monic f
        if r[k]:
            co = r[k]
            for i in range(5):
                r[k - 4 + i] -= co * f[i]
    return tuple(x % p for x in r[:4])


def polpow(u, e, f, p):
    res = (1, 0, 0, 0)
    base = u
    while e:
        if e & 1:
            res = polmul(res, base, f, p)
        base = polmul(base, base, f, p)
        e >>= 1
    return res


for p in (7, 11, 13):
    q = p * p
    f = irreducible_quartic(p)
    n = p ** 4
    els = [(i % p, (i // p) % p, (i // p ** 2) % p, i // p ** 3) for i in range(n)]
    mu = [i for i in range(1, n) if polpow(els[i], q + 1, f, p) == (1, 0, 0, 0)]
    assert len(mu) == q + 1
    coords = np.array(els, dtype=np.int64)
    idx = lambda C: C[:, 0] + p * C[:, 1] + p * p * C[:, 2] + p ** 3 * C[:, 3]
    shifts = [idx((coords - coords[u]) % p) for u in mu]       # v -> v - u
    R = np.zeros(n, dtype=np.int64); R[0] = 1
    W = []
    for k in range(1, 9):
        Rn = np.zeros(n, dtype=np.int64)
        for s in shifts:
            Rn += R[s]
        R = Rn
        W.append(int(R[0]))
    M = []
    for k in range(1, 9):
        num = q * q * W[k - 1] - (q + 1) ** k
        assert num % (q + 1) == 0
        M.append((-1) ** k * (num // (q + 1)))
    chi = sympy.jacobi_symbol(-15 % p, p)
    A2 = A15[p] ** 2 - 2 * chi * p * p
    B2 = B46[p] ** 2 - 2 * p ** 3
    D2 = D66[p] ** 2 - 2 * p ** 5
    pred = {1: 1, 2: q * q - q - 1, 3: q * q + 2 * q + 1, 4: 2 * q ** 3 - 3 * q * q - 3 * q - 1,
            5: 4 * q ** 3 + 5 * q * q + 4 * q + 1 + q * q * A2,
            6: 5 * q ** 4 - 10 * q ** 3 - 9 * q * q - 5 * q - 1 - q * q * B2,
            8: 14 * q ** 5 - 35 * q ** 4 - 28 * q ** 3 - 20 * q * q - 7 * q - 1 - 7 * q ** 3 * B2 - q * q * D2}
    ok = {k: M[k - 1] == v for k, v in pred.items()}
    print('q = %d: exact M_k(q) = %s' % (q, M))
    print('        = session values: %s; closed forms with Frobenius-squared traces: %s' % (M == SESSION[q], ok))
    assert M == SESSION[q] and all(ok.values())
print('ALL CHECKS PASSED')
