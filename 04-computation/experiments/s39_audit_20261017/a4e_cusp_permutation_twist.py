#!/usr/bin/env python3
"""S39 audit, item 4: is the linkage family L_c the 'same surface with the cusps permuted' as E'_c?

a4a: j(L_c) = j(E'_(sigma(c))) with sigma(c) = (c - 9)/(c - 1), and Phi_2(j(E'_c), j(L_c)) = 0.
Here: compare #L_c(F_p) with #E'_(sigma(c))(F_p) (isomorphism over Q(c), or a quadratic twist by some d(c)?), and record the
twist character if any.  Exact Legendre sums, primes 7 <= p < 100.
Usage: python -B a4e_cusp_permutation_twist.py
"""
import numpy as np
from sympy import primerange


def qr(p):
    t = -np.ones(p, dtype=np.int64); t[0] = 0
    t[(np.arange(1, p) ** 2) % p] = 1
    return t


def kr(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


same, twisted, other = 0, 0, 0
twist_vals = {}
for p in primerange(7, 100):
    chi = qr(p)
    tt = np.arange(p, dtype=np.int64)

    def nE(c):
        Q = ((tt * tt - 4 * c) % p) * ((tt - c - 1) % p) % p * ((tt - c + 3) % p) % p
        return p + 2 + int(chi[Q].sum())

    def nL(c):
        A = (tt * tt + (1 + c) * tt + 1) % p
        Q = (A * A - 4 * c * tt % p * ((tt + 1) * (tt + 1) % p)) % p
        return p + 2 + int(chi[Q].sum())

    for c in range(2, p):
        if c in (1, 9 % p) or (c - 1) % p == 0:
            continue
        s = (c - 9) * pow(c - 1, -1, p) % p
        if s in (0, 1, 9 % p):
            continue
        aL = p + 1 - nL(c); aE = p + 1 - nE(s)
        if aL == aE:
            same += 1
        elif aL == -aE:
            twisted += 1
            twist_vals.setdefault(p, []).append(c)
        else:
            other += 1
print('a_p(L_c) vs a_p(E\'_(sigma c)), sigma(c) = (c-9)/(c-1), all primes 7 <= p < 100: equal %d, opposite %d, neither %d'
      % (same, twisted, other))
# guess a twist character: test d(c) in a small list
cands = {'1': lambda c, p: 1, '-1': lambda c, p: kr(-1, p), '-3': lambda c, p: kr(-3, p), 'c': lambda c, p: kr(c, p),
         'c-1': lambda c, p: kr(c - 1, p), '-(c-1)': lambda c, p: kr(-(c - 1), p), 'c-9': lambda c, p: kr(c - 9, p),
         'c(c-1)': lambda c, p: kr(c * (c - 1), p), '-3(c-1)': lambda c, p: kr(-3 * (c - 1), p),
         '(c-1)(c-9)': lambda c, p: kr((c - 1) * (c - 9), p), '-(c-1)(c-9)': lambda c, p: kr(-(c - 1) * (c - 9), p),
         'c(c-9)': lambda c, p: kr(c * (c - 9), p), '-3c(c-9)': lambda c, p: kr(-3 * c * (c - 9), p)}
for name, f in cands.items():
    ok = True
    for p in primerange(7, 100):
        chi = qr(p); tt = np.arange(p, dtype=np.int64)
        for c in range(2, p):
            if c in (1, 9 % p) or (c - 1) % p == 0:
                continue
            s = (c - 9) * pow(c - 1, -1, p) % p
            if s in (0, 1, 9 % p):
                continue
            A = (tt * tt + (1 + c) * tt + 1) % p
            QL = (A * A - 4 * c * tt % p * ((tt + 1) * (tt + 1) % p)) % p
            QE = ((tt * tt - 4 * s) % p) * ((tt - s - 1) % p) % p * ((tt - s + 3) % p) % p
            aL = -1 - int(chi[QL].sum()); aE = -1 - int(chi[QE].sum())
            if aL != f(c, p) * aE:
                ok = False; break
        if not ok:
            break
    if ok:
        print('a_p(L_c) = (d/p) a_p(E\'_(sigma c)) holds for d =', name)
