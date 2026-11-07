"""Check (A') the three-block identity and (C) the balance lemma in its exact proved form:
   |A - S| <= 2 or max(A, S) <= 3, with A = a - dXi, S = s - dX (reduced words, non-converse... and converse pairs too)."""
import sys, random
import numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oalib import gtable

def C(x, w):
    r = 0
    for i in range(w):
        r |= (1 - ((x >> i) & 1)) << (w - 1 - i)
    return r

N = 22
g = gtable(N)
# (A') three-block identity on random words
random.seed(7)
bad = 0
for _ in range(200000):
    n = random.randint(3, N)
    v = random.randrange(1 << (n - 1), 1 << n)
    s = random.randint(0, n - 2)
    w = random.randint(1, n - 1 - s)
    a = n - w - s
    h = v >> (w + s)
    x = (v >> s) & ((1 << w) - 1)
    sig = v & ((1 << s) - 1)
    tau = C(sig, s) if s > 0 else 0
    eps = ((x & 1) * (tau & 1)) + ((h & 1) * (1 - ((x >> (w - 1)) & 1)))
    rhs = int(g[h]) + int(g[x]) + int(g[tau]) + h * tau * (1 << w) + h * C(x, w) + tau * x + eps - 2
    if rhs != int(g[v]):
        bad += 1
print("three-block identity violations:", bad)
# g(x) <= x + 1
print("max over x < 2^22 of g(x) - x - 1 =", int((g - np.arange(1 << N) - 1).max()))
# balance lemma, all collision pairs (converse included) of reduced words, n = 10..20
for n in range(10, 21):
    v = np.arange(1 << (n - 1), 1 << n, 2, dtype=np.int64)
    gv = g[v]
    order = np.argsort(gv, kind='stable')
    sv, gs = v[order], gv[order]
    st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
    en = np.r_[st[1:], len(gs)]
    viol = 0; pairs = 0; tight = 0
    for i0, i1 in zip(st, en):
        if i1 - i0 < 2:
            continue
        ws = sv[i0:i1].tolist()
        for p in ws:
            for q in ws:
                if p == q:
                    continue
                pairs += 1
                xx = p ^ q
                t = xx.bit_length() - 1
                s = (xx & -xx).bit_length() - 1
                a = n - 1 - t
                w = t - s + 1
                xp = (p >> s) & ((1 << w) - 1); xq = (q >> s) & ((1 << w) - 1)
                dX = w - abs(xp - xq).bit_length()
                dXi = w - abs(C(xp, w) - C(xq, w)).bit_length()
                A, S = a - dXi, s - dX
                if not (abs(A - S) <= 2 or max(A, S) <= 3):
                    viol += 1
                if abs(A - S) == 2:
                    tight += 1
    print("n=%d collision pairs (incl. converse) %d, balance-lemma violations %d, pairs with |A-S| = 2: %d" % (n, pairs, viol, tight))
