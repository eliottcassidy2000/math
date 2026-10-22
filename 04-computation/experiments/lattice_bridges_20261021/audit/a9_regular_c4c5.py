#!/usr/bin/env python3
"""AUDIT of proposed THM-4751 (6)(i): is c_5 + 2 c_4 constant only on circulant tournaments, or on ALL regular
tournaments of order n (which would make (i) a special case of a general identity)?
 * n = 5, 7: all labelled regular tournaments (exhaustive), c_4 = tr A^4 / 4, c_5 = tr A^5 / 5 (THM-118);
 * n = 9, 11, 13: random regular tournaments obtained from R_n by random cyclic-triangle reversals (preserve scores).
Also prints the values of c_4 at the circulant extremes and the closed forms c_4 = n(E - m)/4 and
2 c_5 + n E = (2m^5 + 5m^4 + 5m^2 + 3m)/5 derived in the audit report.
Output: a9_regular_c4c5.out
"""
import os, itertools, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a9_regular_c4c5.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()
def c45(A):
    A2 = A @ A; A4 = A2 @ A2; A5 = A4 @ A
    t4, t5 = int(np.trace(A4)), int(np.trace(A5))
    assert t4 % 4 == 0 and t5 % 5 == 0
    return t4 // 4, t5 // 5
def regular_all(n):
    m = (n - 1) // 2
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    res = []
    def rec(k, A, outdeg, indeg):
        if k == len(pairs):
            res.append(A.copy()); return
        i, j = pairs[k]
        for (u, v) in ((i, j), (j, i)):
            if outdeg[u] < m and indeg[v] < m:
                A[u, v] = 1; outdeg[u] += 1; indeg[v] += 1
                rec(k + 1, A, outdeg, indeg)
                A[u, v] = 0; outdeg[u] -= 1; indeg[v] -= 1
    rec(0, np.zeros((n, n), dtype=np.int64), [0] * n, [0] * n)
    return res
for n in (5, 7):
    m = (n - 1) // 2
    regs = regular_all(n)
    vals = {}
    for A in regs:
        c4, c5 = c45(A); vals.setdefault((c4, c5), 0); vals[(c4, c5)] += 1
    consts = {c5 + 2 * c4 for (c4, c5) in vals}
    say(f"n={n}: {len(regs)} labelled regular tournaments; (c4, c5): multiplicity = {vals}; c5 + 2 c4 values = {sorted(consts)}")
random.seed(20261022)
for n in (9, 11, 13):
    m = (n - 1) // 2
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for d in range(1, m + 1): A[i, (i + d) % n] = 1
    seen = {}
    for step in range(4000):
        # pick a random cyclic triangle and reverse it
        while True:
            a, b, c = random.sample(range(n), 3)
            if A[a, b] and A[b, c] and A[c, a]: break
        A[a, b] = A[b, c] = A[c, a] = 0; A[b, a] = A[c, b] = A[a, c] = 1
        if step % 10 == 0:
            c4, c5 = c45(A); seen[(c4, c5)] = seen.get((c4, c5), 0) + 1
    consts = sorted({c5 + 2 * c4 for (c4, c5) in seen})
    say(f"n={n}: random regular tournaments: {len(seen)} distinct (c4, c5) pairs; c5 + 2 c4 takes {len(consts)} values: {consts[:8]}{' ...' if len(consts) > 8 else ''}")
    circ_const = (2 * m ** 5 + 5 * m ** 4 + 5 * m ** 2 + 3 * m) // 5 - n * m  # 2c5+nE = K, 4c4 = nE - nm  ->  c5 + 2c4 = (K - nm)/2
    say(f"   circulant value of c5 + 2c4 (from the closed forms) = {((2*m**5 + 5*m**4 + 5*m**2 + 3*m)//5 - n*m)//2}")

# ---- symbolic proof of c5 + 2 c4 = n(n^2-1)(n^2-9)/160 for EVERY regular tournament (not only circulants) ----
# A = (J - I + B)/2 with B = A - A^T skew, BJ = JB = 0 (regularity).  With C = B - I: JC = CJ = -J, J^2 = nJ, so
# (J + C)^k = C^k + alpha_k J,  alpha_1 = 1, alpha_(k+1) = (-1)^k + (n-1) alpha_k.  tr B^odd = 0, tr B^2 = -n(n-1).
import sympy as sp
nn, t4 = sp.symbols('n t4')          # t4 = tr B^4 (the only free invariant)
alpha = {1: sp.Integer(1)}
for k in range(1, 5): alpha[k + 1] = sp.expand((-1) ** k + (nn - 1) * alpha[k])
trB = {0: nn, 1: 0, 2: -nn * (nn - 1), 3: 0, 4: t4, 5: 0}
def trC(k): return sum(sp.binomial(k, j) * (-1) ** (k - j) * trB[j] for j in range(k + 1))
trA4 = (trC(4) + nn * alpha[4]) / 16; trA5 = (trC(5) + nn * alpha[5]) / 32
c4s, c5s = trA4 / 4, trA5 / 5
say("symbolic: c5 + 2 c4 =", sp.factor(sp.simplify(c5s + 2 * c4s)), " (t4 = tr B^4 cancels:", sp.simplify(sp.diff(c5s + 2 * c4s, t4)) == 0, ")")
say("symbolic: c4 =", sp.factor(c4s))
out.close()
