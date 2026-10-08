#!/usr/bin/env python3
"""Audit G: the rank-1 two-valued translation-only sketch (HYP-9244 update item 1 / results note section 4.1, working-tree
version that adds 'for mu > d^2 no balancing weight s exists').

Claim to test: for V = |f|^theta s^|k| (f = e / max(M,1), M = mu^k), the PER-STEP conditional expectation
E[V' | state] / V, computed exactly over the d fresh digits j (u's digit i = j + e mod d), equals kappa(theta)
= (1/d) sum (m_i/d)^theta at s = 1 for EVERY lag (up to the additive terms, made negligible by taking |f| huge),
so s slightly below 1 gives a per-step drift < 1 for every lag, whatever mu (including mu > d^2).
The 'pure flip' drift (movers only) is reported too: it exceeds 1 for mu > d^2, which is the author's obstruction --
but in a d-adic step the movers never occur alone (odd prime d: every lag cycle has stayers)."""
from fractions import Fraction as Fr
import math, random

def drift(d, m, r, k, e, theta, s):
    mu = max(m); M = Fr(mu) ** k
    f = e / max(M, 1)
    tot = 0.0; movers = [];
    for j in range(d):
        i = (j + (e.numerator * pow(e.denominator, -1, d))) % d   # M = 1 mod d (translation-only)
        Mp = M * Fr(m[i], m[j]); ep = (m[i] * e + r[i] - Fr(m[i], m[j]) * M * r[j]) / d
        fp = ep / max(Mp, 1); kp = k + (1 if m[i] > m[j] else -1 if m[i] < m[j] else 0)
        w = (abs(float(fp)) / abs(float(f))) ** theta * s ** (abs(kp) - abs(k))
        tot += w / d
        if kp != k: movers.append(w)
    return tot, (sum(movers) / len(movers) if movers else None)

rnd = random.Random(3)
for d, m in ((3, [1, 1, 16]), (5, [1, 1, 1, 1, 26]), (5, [1, 26, 1, 1, 1]), (3, [1, 1, 4]), (7, [1, 1, 1, 1, 1, 1, 64])):
    r = [(-m[i] * i) % d for i in range(d)]
    mu = max(m); theta = 0.1
    kappa = sum((x / d) ** theta for x in m) / d
    print(f"Z_{d} m = {m} (mu = {mu} {'>' if mu > d * d else '<='} d^2 = {d*d}), theta = {theta}: kappa = {kappa:.5f};"
          f"  pure-flip drift min_s = (mu/d^2)^(theta/2) = {(mu / d ** 2) ** (theta / 2):.5f}")
    for s in (1.0, 0.995, 0.98):
        worst = 0; worst_flip = 0
        for k in (-7, -3, -1, 0, 1, 2, 5, 9):
            for res in range(d):
                e = Fr(10 ** 40 * d + res + rnd.randint(0, 10 ** 6) * d)    # residue res mod d, huge
                if k > 0: e = e * Fr(mu) ** k                               # keep |f| huge at positive debt
                tdr, fl = drift(d, m, r, k, e, theta, s)
                worst = max(worst, tdr)
                if fl is not None: worst_flip = max(worst_flip, fl)
        print(f"   s = {s}: max over levels/lags of per-step drift = {worst:.5f};  max movers-only drift = {worst_flip:.5f}")
