#!/usr/bin/env python3
"""procgen_edim2_20261001_orchestrator_check.py -- the orchestrator's independent audit of the edim2 lane (2026-10-01).

Written from the statements in the lane note; the lane's code was not read (only the certificate DATA is parsed from the
CERT literal of procgen_edim2_20261001_run.py).
Checks:
  A. The explicit certificates for Q_6 .. Q_16 resolve every edge (own numpy checker); sizes 15, 19, 26, 38, 48, 65, 76,
     105, 125, 135, 171.
  B. Lemma A (Fourier atom bound) and Lemma A2 (atoms cannot be smaller) against exact atoms of Bin(n,q) - Bin(n',q).
  C. Proposition K: for parallel pairs par(h) in Q_d (d <= 8) the weighted level graph has weighted spanning-tree count
     2^n prod_{k=1}^{n} (C(n,k) - K_k(h)); it is connected iff h is odd and h < n (the note's 'iff h odd' fails at h = n).
  D. Constants: C* = (3 sqrt2 ln2)^(2/3) = 2.0526, c* = (3 ln2/(2 sqrt2))^(2/3) = 0.8146, C*/c* = 4^(2/3), and
     (sqrt2/3) C*^(3/2) = 2 ln 2.
"""
import ast
import math
import os
import random
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OKS = []


def ok(c, msg):
    OKS.append(bool(c))
    print(('[OK] ' if c else '[FAIL] ') + msg, flush=True)


PC = np.array([bin(i).count('1') for i in range(1 << 16)], dtype=np.int16)


def resolving(d, S):
    us = np.array([u for u in range(1 << d) for i in range(d) if not (u >> i) & 1], dtype=np.int64)
    vs = np.array([u | (1 << i) for u in range(1 << d) for i in range(d) if not (u >> i) & 1], dtype=np.int64)
    E = len(us)
    H = np.zeros((E, d), dtype=np.int16)
    rows = np.arange(E)
    for s in S:
        dm = np.minimum(PC[us ^ s], PC[vs ^ s])
        H[rows, dm] += 1
    return len(np.unique(H, axis=0)) == E, E


def check_A():
    src = open(os.path.join(HERE, 'procgen_edim2_20261001_run.py')).read()
    i = src.index('CERT = {')
    j = src.index('}', i)
    node = ast.parse(src[i + len('CERT = '):j + 1], mode='eval').body
    cert = {}
    for kn, vn in zip(node.keys, node.values):
        k = ast.literal_eval(kn)
        try:
            cert[k] = ast.literal_eval(vn)
        except ValueError:
            assert k == 6  # the paper's Q_6 set, written as a mask comprehension in the lane's runner
            cert[k] = [v for v in range(64) if (0x02283022a042a00a >> v) & 1]
    sizes = {6: 15, 7: 19, 8: 26, 9: 38, 10: 48, 11: 65, 12: 76, 13: 105, 14: 125, 15: 135, 16: 171}
    for d in sorted(cert):
        S = cert[d]
        good = len(set(S)) == len(S) == sizes[d] and all(0 <= s < 2 ** d for s in S)
        res, E = resolving(d, S)
        ok(good and res, f'A: Q_{d}: explicit set of size {len(S)} resolves all {E} edges  =>  edim_m(Q_{d}) <= {len(S)}')


def binom_pmf(n, q):
    p = [0.0] * (n + 1)
    lp = [math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * math.log(q) + (n - k) * math.log(1 - q)
          for k in range(n + 1)]
    return [math.exp(x) for x in lp]


def check_B():
    rng = random.Random(9169)
    bad_a = bad_a2 = 0
    for _ in range(300):
        n = rng.randint(1, 300)
        n2 = rng.randint(0, 300)
        q = rng.choice([0.5, 0.3, 0.25, 0.1, 0.03, 0.01])
        a = np.array(binom_pmf(n, q))
        b = np.array(binom_pmf(n2, q))
        conv = np.convolve(a, b[::-1])
        atom = conv.max()
        N = n + n2
        x = N * q * (1 - q)
        bound = (1 + 1 / (4 * x)) / math.sqrt(2 * math.pi * x) + math.exp(-x) / 2
        if atom > bound * (1 + 1e-9):
            bad_a += 1
        lower = (12 * x + 1) ** -0.5
        if atom < lower * (1 - 1e-9):
            bad_a2 += 1
    ok(bad_a == 0 and bad_a2 == 0, f'B: Lemma A and Lemma A2 on 300 random (n, n\', q): {bad_a} / {bad_a2} violations')


def krawtchouk(n, k, h):
    return sum((-1) ** j * math.comb(h, j) * math.comb(n - h, k - j) for j in range(0, k + 1))


def det_frac(M):
    M = [row[:] for row in M]
    n = len(M)
    det = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            det = -det
        det *= M[c][c]
        for r in range(c + 1, n):
            if M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return det


def check_C():
    good = True
    conn_ok = True
    for d in range(3, 9):
        n = d - 1
        for h in range(1, n + 1):
            # e = {0, e_0}; f = {x, x + e_0} with x = bits 1..h set
            x = sum(1 << i for i in range(1, h + 1))
            N = [[0] * (n + 1) for _ in range(n + 1)]
            for w in range(1 << d):
                a = min(bin(w).count('1'), bin(w ^ 1).count('1'))
                b = min(bin(w ^ x).count('1'), bin(w ^ x ^ 1).count('1'))
                if a != b:
                    N[a][b] += 1
            W = [[N[a][b] + N[b][a] for b in range(n + 1)] for a in range(n + 1)]
            L = [[Fraction(0)] * (n + 1) for _ in range(n + 1)]
            for a in range(n + 1):
                for b in range(n + 1):
                    if a != b and W[a][b]:
                        L[a][b] -= W[a][b]
                        L[a][a] += W[a][b]
            red = [row[1:] for row in L[1:]]
            kappa = det_frac(red)
            pred = 2 ** n
            for k in range(1, n + 1):
                pred *= (math.comb(n, k) - krawtchouk(n, k, h))
            if kappa != pred:
                good = False
            if (kappa != 0) != (h % 2 == 1 and h < n):
                conn_ok = False
    ok(good and conn_ok, 'C: Proposition K (weighted spanning-tree count = 2^n prod (C(n,k) - K_k(h))) for every par(h), 3 <= d <= 8; CORRECTION: connected iff h is odd AND h < n (at h = n the level graph is the matching {p, n-p}, disconnected even for odd n)')


def check_D():
    Cs = (3 * math.sqrt(2) * math.log(2)) ** (2 / 3)
    cs = (3 * math.log(2) / (2 * math.sqrt(2))) ** (2 / 3)
    ok(abs(Cs - 2.052616) < 1e-5 and abs(cs - 0.81458) < 1e-4 and abs(Cs / cs - 4 ** (2 / 3)) < 1e-12
       and abs(math.sqrt(2) / 3 * Cs ** 1.5 - 2 * math.log(2)) < 1e-12,
       f'D: C* = {Cs:.6f}, c* = {cs:.6f}, C*/c* = 4^(2/3), (sqrt2/3) C*^(3/2) = 2 ln 2')
    # Theorem B heuristic: integral of (lambda - 2t^2/d)/ln 2 over |t| <= sqrt(d lambda/2) equals d at lambda = c* d^(1/3)
    for d in (10 ** 6, 10 ** 9):
        lam = cs * d ** (1 / 3)
        a = math.sqrt(d * lam / 2)
        integral = (2 * a * lam - (4 / (3 * d)) * a ** 3) / math.log(2)
        ok(abs(integral / d - 1) < 1e-9, f'D: entropy-integral identity at d = {d}: (1/d) * integral = {integral / d:.12f}')


def main():
    t0 = time.time()
    print('==== A ====', flush=True); check_A()
    print('==== B ====', flush=True); check_B()
    print('==== C ====', flush=True); check_C()
    print('==== D ====', flush=True); check_D()
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')


if __name__ == '__main__':
    main()
