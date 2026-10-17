#!/usr/bin/env python3
"""Exact power moments M_k(p) = sum_{a in F_p^*} K(a)^k of the Kloosterman sums K(a) = sum_{x in F_p^*} e((x + a/x)/p),
k = 1..KMAX, for every prime 5 <= p <= PMAX (S39).

Method (the Paley H(P_31) device of S34: exact integers from several word-size primes and the Chinese remainder theorem).
M_k is a rational integer (Galois permutes the K(a): sigma_c K(a) = K(c^2 a)), and |M_k| <= (p - 1)(2 sqrt p)^k (Weil).
For a prime l = 1 mod p with a primitive p-th root of unity w in F_l, the ring map Z[zeta_p] -> F_l, zeta_p -> w, sends
M_k to M_k mod l.  With g a primitive root mod p and f(n) = w^(g^n mod p), K(g^m) = sum_n f(n) f(m - n) is a cyclic
convolution of length p - 1 over F_l.  Four primes l < 2^31 give a modulus above 2^120, far beyond the Weil range.

Usage: python kloosterman_moments_exact_20261017.py PMAX KMAX out.json
Session opus-2026-10-17-S39.
"""
import sys, json, math
import numpy as np
import sympy
from numba import njit


@njit(cache=True)
def moments_mod(table, pw, l, kmax):
    n = pw.shape[0]
    f = np.empty(n, dtype=np.int64)
    for i in range(n):
        f[i] = table[pw[i]]
    mom = np.zeros(kmax + 1, dtype=np.int64)
    for m in range(n):
        s = 0
        for i in range(n):
            j = m - i
            if j < 0:
                j += n
            s = (s + f[i] * f[j]) % l
        x = 1
        for k in range(1, kmax + 1):
            x = (x * s) % l
            mom[k] = (mom[k] + x) % l
    return mom


def word_primes(p, count, start=2**31):
    """count primes l = 1 mod p below start, descending."""
    out = []
    l = start - 1 - ((start - 2) % p)        # l = 1 mod p
    while len(out) < count:
        if l % p == 1 and sympy.isprime(l):
            out.append(l)
        l -= p
    return out


def root_of_unity(p, l):
    e = (l - 1) // p
    for g in range(2, 10**6):
        w = pow(g, e, l)
        if w != 1:
            return w
    raise RuntimeError


def exact_moments(p, kmax, nprimes=4):
    g = sympy.primitive_root(p)
    pw = np.array([pow(g, i, p) for i in range(p - 1)], dtype=np.int64)
    residues, moduli = [], []
    for l in word_primes(p, nprimes):
        w = root_of_unity(p, l)
        table = np.array([pow(w, j, l) for j in range(p)], dtype=np.int64)
        residues.append(moments_mod(table, pw, l, kmax))
        moduli.append(l)
    Mod = 1
    for l in moduli:
        Mod *= l
    out = []
    for k in range(1, kmax + 1):
        r = int(sympy.ntheory.modular.crt(moduli, [int(res[k]) for res in residues])[0])
        if r > Mod // 2:
            r -= Mod
        bound = (p - 1) * (2 * math.sqrt(p)) ** k
        assert abs(r) <= bound + 1, (p, k, r, bound)
        assert 2 * bound < Mod
        out.append(r)
    return out


def float_check(p, kmax):
    g = sympy.primitive_root(p)
    pw = np.array([pow(g, i, p) for i in range(p - 1)], dtype=np.int64)
    f = np.exp(2j * np.pi * pw / p)
    K = np.fft.ifft(np.fft.fft(f) ** 2).real
    return [float(np.sum(K ** k)) for k in range(1, kmax + 1)]


if __name__ == '__main__':
    PMAX, KMAX, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    res = {}
    for p in sympy.primerange(5, PMAX + 1):
        M = exact_moments(p, KMAX)
        if p < 200:
            F = float_check(p, KMAX)
            for k in range(KMAX):
                assert abs(F[k] - M[k]) < 1e-6 * max(1.0, abs(M[k])) + 1e-3, (p, k + 1, F[k], M[k])
        assert M[0] == 1 and M[1] == p * p - p - 1, (p, M[:2])           # classical first and second moments
        res[str(p)] = M
    json.dump(dict(kmax=KMAX, pmax=PMAX, moments=res,
                   note='M_k(p) = sum over a in F_p^* of K(a)^k, K(a) = sum over x in F_p^* of e((x + a/x)/p); exact via CRT'),
              open(OUT, 'w'))
    print('primes:', len(res), 'largest', max(int(k) for k in res))
    for p in ['5', '7', '11', '13', '17', '19']:
        print(p, res[p])
