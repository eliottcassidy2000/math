#!/usr/bin/env python3
"""Where a genuine Hasse-Weil term can first enter Paley TOURNAMENT pattern counts (q = 3 mod 4).

Edge-subset character sums over m distinct points, prod over edges of chi(v_head - v_tail):
  * even subdivisions of an edge act as 'no edge', odd ones as a single edge (path sums P_l(d) ~ chi(d)^(l mod 2));
  * series-parallel graphs reduce to Jacobi sums: polynomials in q;
  * K_m itself is chi(Vandermonde) and vanishes (a label transposition negates it, chi(-1) = -1); so do K_5 - e and
    K_5 minus two adjacent edges;
  * the first survivor is K_5 minus two disjoint edges = the wheel W_4 (hub + 4-cycle), which has a K_4 minor.
W(q) = sum over distinct nonzero a, b, c, d (hub at 0) of chi(a b c d) chi((c-a)(d-a)(c-b)(d-b)) (rim i->k, i->l, j->k,
j->l with i, j = a, b and k, l = c, d).  Test whether W(q) is a polynomial in q (interpolation leave-one-out) or carries
a fluctuating (Weil/Deligne) term; for comparison also q = 1 mod 4.  Session opus-2026-10-14-S34.
"""
import numpy as np
import numba as nb
from fractions import Fraction


@nb.njit(cache=True)
def W(q, chi):
    tot = 0
    for a in range(1, q):
        ca = chi[a]
        for b in range(1, q):
            if b == a:
                continue
            cab = ca * chi[b]
            for c in range(1, q):
                if c == a or c == b:
                    continue
                cabc = cab * chi[c] * chi[(c - a) % q] * chi[(c - b) % q]
                for d in range(1, q):
                    if d == a or d == b or d == c:
                        continue
                    tot += cabc * chi[d] * chi[(d - a) % q] * chi[(d - b) % q]
    return tot


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(3, n + 1) if s[i]]


def interp(xs, ys, x):
    tot = Fraction(0)
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        t = Fraction(yi)
        for j, xj in enumerate(xs):
            if j != i:
                t *= Fraction(x - xj, xi - xj)
        tot += t
    return tot


for cls in (3, 1):
    vals = {}
    for q in primes_upto(110):
        if q % 4 != cls or q < 7:
            continue
        chi = np.zeros(q, dtype=np.int64)
        sq = set((x * x) % q for x in range(1, q))
        for x in range(1, q):
            chi[x] = 1 if x in sq else -1
        vals[q] = int(W(q, chi))
    print(f'q = {cls} mod 4:')
    for q, w in vals.items():
        print(f'   q = {q:3d}: W = {w:12d},  W/(q-1) = {Fraction(w, q - 1)},  W/q^3 = {w / q ** 3:.4f}')
    qs = sorted(vals)
    # W/(q-1) is (if polynomial) of degree <= 3 (W ~ q^4 at most); interpolate through 4 primes, test the rest
    base = qs[:4]
    print(f'   cubic interpolation of W/(q-1) through {base}, residuals at the other primes:')
    for q in qs[4:]:
        pred = interp(base, [Fraction(vals[p], p - 1) for p in base], q)
        print(f'      q = {q}: predicted {float(pred):.2f}, exact {vals[q] / (q - 1):.2f}, residual {float(Fraction(vals[q], q - 1) - pred):+.2f}')
