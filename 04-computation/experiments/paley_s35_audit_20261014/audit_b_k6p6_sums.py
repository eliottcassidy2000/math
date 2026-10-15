#!/usr/bin/env python3
"""Audit (b3) of S35: the K6 - P6 terminal by brute force and by an independent exact O(q^2) evaluator.

Graph (S35 labelling): E = (0,1),(0,3),(0,4),(0,5),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4), factor chi(x_v - x_u) on (u,v).
(Equal, as an oriented graph, to tiling-modular-atoms' K6 - P6 with path 0-1-2-3-4-5 removed and edges low -> high,
 under t0..t5 -> 0,2,1,3,5,4; checked numerically below.)
Y_q = U/(q(q-1)) (unrestricted), W_q = I/(q(q-1)) (injective).  AGL(1,q) normalisation x_0 = 0, x_1 = 1 (edge (0,1),
|E| = 10 even), validated against the un-normalised sum at q = 5, 7.

O(q^2) evaluator (derived for this audit): with L(l) = sum_t chi(t (t-1) (t-l)),
  sum_{x5} chi(x5) chi(x5-1) chi(x5-x2) = L(x2),
  sum_{x3} chi(x3) chi(x3-x2) chi(x4-x3) = chi(-1) chi(x2) L(x4/x2)  (x2 != 0),   = -chi(x4)  (x2 = 0),
  Y_q = sum_{x2, x4} chi(x4) chi(x4-1) chi(x4-x2) L(x2) h(x2, x4).
Exact int64, numba parallel.  Compared with: the session's O(q^3) table (124 primes), the newform side
|b(q)|^2 - q^2 + 2 chi_{-4}(q) q (audit_b_newform_bq.json), HYP-9276's W values, and the normalisation link
W_q - Y_q = 7 + 14 chi - 2 chi q.
"""
import json, re, sys, time
import numpy as np
import numba as nb
import sympy as sp

E = [(0, 1), (0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 4)]
EU = np.array([e[0] for e in E], dtype=np.int64)
EV = np.array([e[1] for e in E], dtype=np.int64)


def chi_table(q):
    c = np.full(q, -1, dtype=np.int64)
    c[0] = 0
    for x in range(1, q):
        c[(x * x) % q] = 1
    return c


@nb.njit(cache=True)
def brute(q, chi, eu, ev, injective, normalise):
    tot = 0
    v = np.zeros(6, dtype=np.int64)
    start = 2 if normalise else 0
    if normalise:
        v[0] = 0
        v[1] = 1
    m = 6 - start
    N = q ** m
    for code in range(N):
        c = code
        for i in range(start, 6):
            v[i] = c % q
            c //= q
        if injective:
            ok = True
            for i in range(6):
                for j in range(i + 1, 6):
                    if v[i] == v[j]:
                        ok = False
            if not ok:
                continue
        p = 1
        for k in range(eu.shape[0]):
            p *= chi[(v[ev[k]] - v[eu[k]]) % q]
            if p == 0:
                break
        tot += p
    return tot


@nb.njit(parallel=True, cache=True)
def legendre_table(q, chi):
    L = np.zeros(q, dtype=np.int64)
    for l in nb.prange(q):
        s = 0
        for t in range(q):
            s += chi[t] * chi[(t - 1) % q] * chi[(t - l) % q]
        L[l] = s
    return L


@nb.njit(parallel=True, cache=True)
def Y_fast(q, chi, L, inv):
    cm1 = chi[q - 1]
    part = np.zeros(q, dtype=np.int64)
    for x2 in nb.prange(q):
        if x2 == 0:
            s = 0
            for x4 in range(q):
                s += chi[x4] * chi[(x4 - 1) % q] * chi[x4] * (-chi[x4])
            part[0] = L[0] * s
            continue
        Lx2 = L[x2]
        if Lx2 == 0:
            continue
        ix2 = inv[x2]
        s = 0
        for x4 in range(q):
            w = chi[x4] * chi[(x4 - 1) % q] * chi[(x4 - x2) % q]
            if w != 0:
                s += w * L[(x4 * ix2) % q]
        part[x2] = Lx2 * cm1 * chi[x2] * s
    return part.sum()


def Y_of(q):
    chi = chi_table(q)
    L = legendre_table(q, chi)
    inv = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        inv[x] = pow(x, q - 2, q)
    return int(Y_fast(q, chi, L, inv))


def main():
    t0 = time.time()
    bq = {int(k): tuple(v) for k, v in json.load(open('audit_b_newform_bq.json')).items()}

    def pred(q):
        re_, im_ = bq[q]
        chi4 = 1 if q % 4 == 1 else -1
        return re_ * re_ + im_ * im_ - q * q + 2 * chi4 * q

    # normalisation check, un-normalised brute force
    for q in (5, 7):
        chi = chi_table(q)
        full = int(brute(q, chi, EU, EV, False, False))
        norm = int(brute(q, chi, EU, EV, False, True))
        print(f'q = {q}: un-normalised U = {full}, q(q-1) * normalised = {q * (q - 1) * norm}, equal: {full == q * (q - 1) * norm}')
    # tiling labelling check
    Et = [(i, j) for i in range(6) for j in range(i + 2, 6)]
    for q in (7, 11):
        chi = chi_table(q)
        a = int(brute(q, chi, EU, EV, False, False)) if q == 7 else None
        relab = {0: 0, 1: 2, 2: 1, 3: 3, 4: 5, 5: 4}
        Er = [(relab[i], relab[j]) for i, j in Et]
        assert sorted(Er) == sorted(E), (sorted(Er), sorted(E))
    print('tiling-modular-atoms K6-P6 (path 0-1-2-3-4-5 removed, edges low -> high) = S35 graph with identical orientation under t -> (0,2,1,3,5,4): True')
    # brute force vs fast vs newform
    HYP = {7: 8, 11: -112, 19: -224, 23: 1064, 31: 56, 43: -1072, 47: -1960}
    rows = []
    for q in sp.primerange(3, 48):
        chi = chi_table(q)
        Yb = int(brute(q, chi, EU, EV, False, True))
        Yf = Y_of(q)
        Wb = int(brute(q, chi, EU, EV, True, True)) if q >= 7 else 0
        chi4 = 1 if q % 4 == 1 else -1
        link = Wb - Yb == 7 + 14 * chi4 - 2 * chi4 * q
        rows.append((q, Yb, Yf, pred(q), Wb, link, HYP.get(q)))
        print(f'q = {q:3d}: Y brute = {Yb:7d}, Y O(q^2) = {Yf:7d}, newform |b|^2 - q^2 + 2chi q = {pred(q):7d};  '
              f'W (injective) brute = {Wb:6d}, W - Y = 7 + 14chi - 2chi q: {link}, HYP-9276 W: {HYP.get(q)}', flush=True)
        assert Yb == Yf
    # session table
    ses = {}
    for l in open('../../../05-knowledge/results/paley_k6p6_fast_table_20261014.out'):
        m = re.match(r'q =\s*(\d+)\s+q mod 16 =\s*\d+\s+Y_q =\s*(-?\d+)', l)
        if m:
            ses[int(m.group(1))] = int(m.group(2))
    bad_s, bad_p = [], []
    for q in sorted(ses):
        Y = Y_of(q)
        if Y != ses[q]:
            bad_s.append(q)
        if Y != pred(q):
            bad_p.append(q)
    print(f'O(q^2) evaluator vs session O(q^3) table at {len(ses)} primes: mismatches {bad_s};  vs newform square law: failures {bad_p}  [{time.time() - t0:.0f} s]', flush=True)
    # larger sample
    rng = np.random.default_rng(2026)
    big = sorted(set(int(sp.prevprime(int(x))) for x in rng.integers(1000, 30000, 25)))
    extra = [int(a) for a in sys.argv[1:]] if len(sys.argv) > 1 else [8191]
    fails = []
    for q in big + extra:
        t1 = time.time()
        Y = Y_of(q)
        ok = Y == pred(q)
        if not ok:
            fails.append(q)
        if q in extra:
            chi4 = 1 if q % 4 == 1 else -1
            print(f'q = {q}: Y_q = {Y}, Y_q - 1 = {Y - 1}, newform prediction {pred(q)}, agree: {ok}  [{time.time() - t1:.0f} s]', flush=True)
    print(f'random sample of {len(big)} primes in [1000, 30000] + {extra}: square-law failures {fails}  [{time.time() - t0:.0f} s]')


if __name__ == '__main__':
    main()
