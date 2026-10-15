#!/usr/bin/env python3
"""Exact check of every intermediate identity in the proof of THM-4653 (statements 2-4), at all primes 5 <= q <= 101.

eps = chi(-1); D = {a, b not in {0,1}, a != b}; w(a,b) = chi(a)chi(b)chi(b-1); f4(x;a) = chi(x)chi(x-1)chi(x-a);
F4 = t; F5(a,b) = sum_x chi(x)chi(x-1)chi(x-a)chi(x-b).  Greene: 2F1(x) = eps(x)(phi(-1)/q) sum_y phi(y)phi(1-y)phi(1-xy),
3F2(x) = (phi(-1)/q) sum_y 2F1(xy) phi(y) phi(1-y)  (checked as exact rationals).
"""
from fractions import Fraction as Fr
import sympy


def run(q):
    c = [0] + [1 if pow(x, (q - 1) // 2, q) == 1 else -1 for x in range(1, q)]
    chi = lambda x: c[x % q]
    eps = chi(-1)
    R = range(q)
    t = [sum(chi(y) * chi(y - 1) * chi(y - a) for y in R) for a in R]
    F5 = [[sum(chi(x) * chi(x - 1) * chi(x - a) * chi(x - b) for x in R) for b in R] for a in R]
    w = lambda a, b: chi(a) * chi(b) * chi(b - 1)
    f4 = lambda x, a: chi(x) * chi(x - 1) * chi(x - a)
    D = [(a, b) for a in R for b in R if a not in (0, 1) and b not in (0, 1) and a != b]
    M = -sum(chi(a) * chi(x) * chi(x - 1) * chi(x - a) * t[a] * t[x] for a in R for x in R)
    main = sum(w(a, b) * F5[a][b] * t[a] for a in R for b in R)
    TA = sum(w(a, b) * F5[a][b] * t[a] for a, b in D)
    TB = -sum(w(a, b) * F5[a][b] * f4(b, a) for a, b in D)
    TC = sum(w(a, b) * (chi(-b) + chi(1 - b) + chi(a - b)) for a, b in D)
    ok = [
        t[0] == -eps and t[1] == -1,
        main == -eps * M,
        sum(w(1, b) * F5[1][b] * t[1] for b in R if b not in (0, 1)) == -1 - eps,
        all(F5[1][b] == -1 - chi(1 - b) for b in R if b not in (0, 1)),
        sum(w(a, a) * F5[a][a] * t[a] for a in R if a not in (0, 1)) == -2 - 2 * eps,
        all(F5[a][a] == -1 - chi(a) * chi(a - 1) for a in R if a not in (0, 1)),
        sum(chi(a - 1) * t[a] for a in R if a not in (0, 1)) == 1 + eps,
        sum(chi(a) * t[a] for a in R if a not in (0, 1)) == 1 + eps,
        TA - main == 3 + 3 * eps,
        TC == 2 + 5 * eps - eps * q,
        TB == 2 + 5 * eps - eps * q,
        all(sum(chi(b - a) * chi(x - b) for b in R if b not in (0, 1, a))
            == -eps - eps * chi(a) * chi(x) - chi(1 - a) * chi(x - 1)
            for a in R if a not in (0, 1) for x in R if x not in (0, 1, a)),
        all(sum(chi(a) * chi(x - a) for a in R if a not in (0, 1)) == -eps - chi(x - 1) for x in R if x not in (0, 1)),
        all(sum(chi(a) * chi(1 - a) * chi(x - a) for a in R) == t[x] for x in R),
    ]
    # statement 3 and 4
    h = [sum(chi(a) * chi(x - a) * t[a] for a in R) for x in R]
    ok.append(all(h[x] == sum(chi(y) * chi(1 - y) * t[x * y % q] for y in R) for x in range(1, q)))
    ok.append(M == -sum(chi(x) * chi(x - 1) * t[x] * h[x] for x in R))
    F21 = [Fr(0)] + [Fr(eps, q) * sum(chi(y) * chi(1 - y) * chi(1 - x * y) for y in R) for x in range(1, q)]
    F32 = [Fr(0)] + [Fr(eps, q) * sum(F21[x * y % q] * chi(y) * chi(1 - y) for y in R) for x in range(1, q)]
    ok.append(all(t[x] == eps * q * F21[x] for x in range(1, q)))
    ok.append(all(h[x] == q * q * F32[x] for x in range(1, q)))
    ok.append(M == -q ** 3 * sum(chi(x) * chi(1 - x) * F21[x] * F32[x] for x in R))
    W = TA + TB + TC + 0
    ok.append(W == -eps * M + 7 + 13 * eps - 2 * eps * q)
    return ok


bad = {}
P = list(sympy.primerange(5, 102))
for q in P:
    ok = run(q)
    if not all(ok):
        bad[q] = [i for i, v in enumerate(ok) if not v]
print(f"THM-4653 proof steps: {len(run(5))} identities at each of the {len(P)} primes 5..101; failures: {bad}")
