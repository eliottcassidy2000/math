#!/usr/bin/env python3
"""HYP-9276 in its Frobenius form at larger prime-power fields, by an O(q^2) evaluation of M_q.

F_q = F_p[x]/(f) with a primitive element found by order testing; chi(y) = (-1)^log(y).  M_q = -(chi t) @ base @ t with
base[a, y] = chi(y) chi(y - 1) chi(y - a) and t = base.sum(1).  Prediction: -(alpha^(2r) + (alpha beta)^r + beta^(2r)) - 2q - 1
from a_p(g) in g_coeffs.out (g = newform of S_3^new(32, (-4/.))).  Usage: python3 k6p6_prime_powers_large.py > ..._large.out
"""
import sys
import time
import itertools
import numpy as np
import sympy
from sympy.polys.galoistools import gf_irreducible_p
from sympy.polys.domains import ZZ

sys.path.insert(0, '.')
from k6p6_prime_powers import predicted_M   # exact Gaussian-integer recursion (same folder)


def polymulmod(a, b, f, p):
    """a, b: coefficient lists low-degree first, length r; f monic high-first of degree r."""
    r = len(a)
    prod = [0] * (2 * r - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    prod[i + j] = (prod[i + j] + x * y) % p
    tail = f[1:]                                   # x^r = -(tail . (x^(r-1), ..., 1))
    for k in range(2 * r - 2, r - 1, -1):
        c = prod[k]
        if c:
            prod[k] = 0
            for j, fc in enumerate(tail):          # tail[j] multiplies x^(r-1-j)
                prod[k - r + (r - 1 - j)] = (prod[k - r + (r - 1 - j)] - c * fc) % p
    return prod[:r]


def field_chi_sub(p, r):
    q = p ** r
    for tailc in itertools.product(range(p), repeat=r):
        f = [1] + list(tailc)
        if gf_irreducible_p(f, p, ZZ):
            break
    one = [1] + [0] * (r - 1)
    fac = sympy.factorint(q - 1)

    def powm(a, e):
        res, base = one[:], a[:]
        while e:
            if e & 1:
                res = polymulmod(res, base, f, p)
            base = polymulmod(base, base, f, p)
            e >>= 1
        return res

    for cand in itertools.product(range(p), repeat=r):
        g = list(cand)
        if not any(g):
            continue
        if all(powm(g, (q - 1) // l) != one for l in fac):
            break
    pw = p ** np.arange(r)
    enc = lambda v: int(sum(c * p ** k for k, c in enumerate(v)))
    logt = np.full(q, -1, dtype=np.int64)
    cur = one[:]
    for k in range(q - 1):
        logt[enc(cur)] = k
        cur = polymulmod(cur, g, f, p)
    assert cur == one and (logt[1:] >= 0).all()
    chi = np.where(logt % 2 == 0, 1, -1).astype(np.int8)
    chi[0] = 0
    digs = np.array([[(e // p ** k) % p for k in range(r)] for e in range(q)], dtype=np.int16)
    sub = np.zeros((q, q), dtype=np.int32)
    for k in range(r):                               # accumulate digit-wise differences to save memory
        sub += (((digs[:, None, k] - digs[None, :, k]) % p) * int(pw[k])).astype(np.int32)
    return q, chi, sub


def main():
    ap = {}
    for line in open('g_coeffs.out'):
        parts = line.split()
        if len(parts) == 3 and parts[0].isdigit():
            ap[int(parts[0])] = (int(parts[1]), int(parts[2]))
    for p, r in [(5, 4), (3, 6), (11, 3), (3, 7), (7, 4)]:
        t0 = time.time()
        q, chi, sub = field_chi_sub(p, r)
        c = chi.astype(np.int64)
        base = c[None, :] * c[sub[:, 1]][None, :] * c[sub.T]          # base[a, y] = chi(y) chi(y-1) chi(y-a)
        t = base.sum(axis=1)
        M = -int((c * t) @ base @ t)
        del base, sub
        eps = int(c[(q - 1)]) if False else (1 if (q - 1) // 2 % 2 == 0 else -1)
        Mp = predicted_M(p, r, ap[p])
        V = q * (q - 1) * (-eps * M + 7 + 13 * eps - 2 * eps * q)
        print(f"q = {p}^{r} = {q:5d}  eps = {eps:+d}  M_q = {M}  predicted {Mp}  HYP-9276 {'ok' if M == Mp else 'FAIL'}  "
              f"V_q = {V}  [{time.time() - t0:.1f}s]", flush=True)


if __name__ == '__main__':
    main()
