#!/usr/bin/env python3
"""THM-4653 and HYP-9276 over prime-power fields F_q, q = p^r (r >= 2), and small primes.

For each q:
  (1) W_q from the normalised injective six-vertex sum (O(q^4) brute force for q <= QBRUTE),
  (2) W_q from the THM-4653 inclusion-exclusion formula (O(q^3)),
  (3) M_q = -sum_(a,x) chi(a) chi(x) chi(x-1) chi(x-a) t(a) t(x),  t(a) = sum_y chi(y) chi(y-1) chi(y-a)  (O(q^2)),
and checks  W_q = -eps M_q + 7 + 13 eps - 2 eps q  (eps = chi_q(-1))   [THM-4653, all odd q]
and the prime-power form of HYP-9276:
  M_q = -(alpha^(2r) + (alpha beta)^r + beta^(2r)) - 2q - 1,
alpha, beta the Satake roots at p of the weight-3 newform g in S_3^new(32, (-4/.)) (alpha + beta = a_p(g),
alpha beta = (-4/p) p^2), read from g_coeffs.out (gp).  Exact integer / Gaussian-integer arithmetic throughout.

Usage: python3 k6p6_prime_powers.py > k6p6_prime_powers.out
"""
import itertools
import numpy as np
from sympy.polys.galoistools import gf_irreducible_p
from sympy.polys.domains import ZZ

QBRUTE = 49


def field(p, r):
    """Elements 0..q-1 = base-p digit vectors; returns (add, sub, mul) tables as numpy int arrays."""
    q = p ** r
    if r == 1:
        x = np.arange(q)
        return (x[:, None] + x[None, :]) % q, (x[:, None] - x[None, :]) % q, (x[:, None] * x[None, :]) % q
    for tail in itertools.product(range(p), repeat=r):
        f = [1] + list(tail)                       # monic, high degree first
        if gf_irreducible_p(f, p, ZZ):
            break
    digs = np.array([[(e // p ** k) % p for k in range(r)] for e in range(q)])   # low degree first
    pw = p ** np.arange(r)
    add = ((digs[:, None, :] + digs[None, :, :]) % p) @ pw
    sub = ((digs[:, None, :] - digs[None, :, :]) % p) @ pw
    # multiplication: polynomial product reduced mod f
    red = {}                                     # x^k mod f as digit vector (low first), k < 2r - 1
    cur = [0] * r
    for k in range(2 * r - 1):
        if k < r:
            v = [0] * r; v[k] = 1
        else:
            v = [0] + red[k - 1][:-1]            # multiply previous by x
            top = red[k - 1][-1]
            # x^r = -(f[1] x^(r-1) + ... + f[r])
            for j in range(r):
                v[j] = (v[j] - top * f[r - j]) % p
        red[k] = v
    mul = np.zeros((q, q), dtype=np.int64)
    for a in range(q):
        da = digs[a]
        for b in range(q):
            db = digs[b]
            out = [0] * r
            for i in range(r):
                if da[i]:
                    for j in range(r):
                        if db[j]:
                            c = da[i] * db[j]
                            v = red[i + j]
                            for k in range(r):
                                out[k] += c * v[k]
            mul[a, b] = sum((out[k] % p) * p ** k for k in range(r))
    return add, sub, mul


def chi_table(mul, q):
    c = -np.ones(q, dtype=np.int64)
    c[np.unique(np.diag(mul)[1:])] = 1
    c[0] = 0
    return c


def tables(p, r):
    q = p ** r
    add, sub, mul = field(p, r)
    c = chi_table(mul, q)
    one = 1
    neg = sub[0]                                   # neg[x] = -x
    eps = int(c[neg[one]])
    return q, sub, c, eps


def W_brute(q, sub, c):
    tot = 0
    for x2 in range(2, q):                       # x2 != 0 (chi), x2 != 1 (non-edge 1-2)
        w2 = c[x2]
        for x3 in range(q):
            if x3 == x2:
                continue
            w3 = w2 * c[x3] * c[sub[x3, 1]]
            if w3 == 0:
                continue
            for x4 in range(q):
                if x4 == x3:
                    continue
                w4 = w3 * c[x4] * c[sub[x4, 1]] * c[sub[x4, x2]]
                if w4 == 0:
                    continue
                # x5 != x4: sum over all x5 of f5, minus the x5 = x4 term
                f5 = c * c[sub[:, 1]] * c[sub[:, x2]] * c[sub[:, x3]]
                tot += int(w4 * (f5.sum() - f5[x4]))
    return tot


def W_formula(q, sub, c):
    base = c[None, :] * c[sub[:, 1]][None, :] * c[sub.T]          # base[a, x] = chi(x) chi(x-1) chi(x-a)
    F4 = base.sum(axis=1)
    C = c[sub]                                                      # C[x, b] = chi(x - b)
    F5 = base @ C                                                   # F5[a, b]
    A = np.arange(q)[:, None]; B = np.arange(q)[None, :]
    mask = (A != 0) & (A != 1) & (B != 0) & (B != 1) & (A != B)
    w = c[A] * c[B] * c[sub[B, 1]]
    S = F5 * (F4[:, None] - base) + c[sub[0, B]] + c[sub[1, B]] + c[sub[A, B]]
    return int((mask * w * S).sum())


def M_value(q, sub, c):
    base = c[None, :] * c[sub[:, 1]][None, :] * c[sub.T]
    t = base.sum(axis=1)                                            # t(a)
    K = c[:, None] * base                                           # K[a, x] = chi(a) chi(x) chi(x-1) chi(x-a)
    return -int(t @ K @ t)


def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def predicted_M(p, r, ap):
    """-(alpha^(2r) + (alpha beta)^r + beta^(2r)) - 2 p^r - 1 with s_k = alpha^k + beta^k in Z[i]."""
    nb = (1 if p % 4 == 1 else -1) * p * p                          # alpha beta
    s = [(2, 0), ap]
    for k in range(2, 2 * r + 1):
        a1 = gmul(ap, s[k - 1])
        s.append((a1[0] - nb * s[k - 2][0], a1[1] - nb * s[k - 2][1]))
    s2r = s[2 * r]
    assert s2r[1] == 0, (p, r, s2r)
    return -(s2r[0] + nb ** r) - 2 * p ** r - 1


def main():
    ap = {}
    for line in open('g_coeffs.out'):
        parts = line.split()
        if len(parts) == 3 and parts[0].isdigit():
            ap[int(parts[0])] = (int(parts[1]), int(parts[2]))
    cases = [(3, 2), (3, 3), (5, 2), (7, 2), (3, 4), (11, 2), (5, 3), (13, 2), (3, 5), (17, 2), (7, 3), (19, 2)]
    cases += [(p, 1) for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)]
    nthm = nhyp = nbr = 0
    for p, r in cases:
        q, sub, c, eps = tables(p, r)
        Wf = W_formula(q, sub, c)
        M = M_value(q, sub, c)
        cq = 7 + 13 * eps - 2 * eps * q
        thm = (Wf == -eps * M + cq)
        nthm += thm
        br = ''
        if q <= QBRUTE:
            Wb = W_brute(q, sub, c)
            assert Wb == Wf, (q, Wb, Wf)
            nbr += 1
            br = ' = brute force'
        Mp = predicted_M(p, r, ap[p])
        hyp = (M == Mp)
        nhyp += hyp
        print(f"q = {p}^{r} = {q:4d}  eps = {eps:+d}  W_q = {Wf}{br}  V_q = {q * (q - 1) * Wf}  M_q = {M}  "
              f"THM-4653 {'ok' if thm else 'FAIL'}  predicted M_q = {Mp}  HYP-9276 {'ok' if hyp else 'FAIL'}")
        assert thm
    print(f"THM-4653 holds at all {nthm} fields (brute force at {nbr} with q <= {QBRUTE}); "
          f"HYP-9276 prime-power form holds at {nhyp} of {len(cases)}")


if __name__ == '__main__':
    main()
