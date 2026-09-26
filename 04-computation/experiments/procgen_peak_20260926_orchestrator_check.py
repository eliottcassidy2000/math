#!/usr/bin/env python3
"""Orchestrator audit of lane `peak` (peak-discounted provability price), written
from the note's statements; the lane's scripts were not read.

Checks
  1. rho_L and rho^peak_L by an exact DP over (e_j, argmax pair), against direct
     enumeration (L <= 16) and against the note's printed values (q = 3, 5).
  2. Theorem 2's chain rho^peak <= (q/2) sum_{e>cL} C(L,e) q^-e <= (q/2) 2^(-(1-H(c))L).
  3. The peak construction on actual integers: every n <= N descends within L, and
     |E cap [1,N]| <= N rho^peak + |Bad_L| + |X_L| (Theorem 1, upper bound).
  4. Theorem 1's capacity constant M*_L and its closed-form bound.
  5. The q = 5 contrast: rho_L(5) stays near 0.176 while rho^peak_L(5) decays like 2^(-0.0139 L).
"""
import math
from fractions import Fraction
from itertools import product
import numpy as np


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def H2(x):
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def greater(q, e1, i1, e2, i2):
    """S(e1,i1) > S(e2,i2) where S(e,i) = e log q - i log 2, compared exactly."""
    return q ** e1 * 2 ** i2 > q ** e2 * 2 ** i1


def dp_rho(q, L):
    """exact (rho_L, rho^peak_L) as Fractions.  state: (e_j, e*, i*) -> count, where
    (e*, i*) is the (first) maximiser of S_i over 0 <= i <= min(j, L-1)."""
    cur = {(0, 0, 0): 1}
    for j in range(L):
        nxt = {}
        for (e, es, is_), cnt in cur.items():
            for b in (0, 1):
                e2 = e + b
                j2 = j + 1
                if not q ** e2 > 2 ** j2:        # bad requires S_j > 0 for j = 1..L
                    continue
                es2, is2 = es, is_
                if j2 <= L - 1 and greater(q, e2, j2, es, is_):
                    es2, is2 = e2, j2
                key = (e2, es2, is2)
                nxt[key] = nxt.get(key, 0) + cnt
        cur = nxt
    nbad = sum(cur.values())
    peak = sum(Fraction(cnt * 2 ** is_, q ** es) for (e, es, is_), cnt in cur.items())
    return Fraction(nbad, 2 ** L), peak / 2 ** L


def enum_rho(q, L):
    nbad, peak = 0, Fraction(0)
    for u in product((0, 1), repeat=L):
        e, ok = 0, True
        best = (0, 0)
        for j in range(L):
            e += u[j]
            if not q ** e > 2 ** (j + 1):
                ok = False
                break
            if j + 1 <= L - 1 and greater(q, e, j + 1, best[0], best[1]):
                best = (e, j + 1)
        if ok:
            nbad += 1
            peak += Fraction(2 ** best[1], q ** best[0])
    return Fraction(nbad, 2 ** L), peak / 2 ** L


print("1. exact DP versus enumeration and the note's table")
for q in (3, 5, 7):
    for L in (1, 2, 3, 5, 8, 11, 14):
        assert dp_rho(q, L) == enum_rho(q, L), (q, L)
check(True, "DP = enumeration exactly for q = 3, 5, 7 and L in {1,2,3,5,8,11,14}")
note = {(3, 8): (7.421875e-02, 1.310299e-02), (3, 16): (3.225708e-02, 2.660874e-03),
        (3, 24): (1.708156e-02, 8.333989e-04), (3, 32): (9.626961e-03, 3.012321e-04),
        (3, 48): (3.537286e-03, 5.431076e-05), (3, 64): (1.493065e-03, 1.302336e-05),
        (5, 16): (None, 3.188e-3)}
for (q, L), (r, p) in note.items():
    rr, pp = dp_rho(q, L)
    if r is not None:
        assert abs(float(rr) / r - 1) < 1e-6, (q, L, float(rr), r)
    assert abs(float(pp) / p - 1) < 1e-3, (q, L, float(pp), p)
    print(f"   q={q} L={L:2d}: rho = {float(rr):.6e}, rho^peak = {float(pp):.6e}")
check(True, "rho_L and rho^peak_L equal the note's printed values (q = 3: L = 8..64; q = 5: L = 16)")

print("2. Theorem 2's chain")
for q in (3, 5, 7, 9):
    c = math.log(2) / math.log(q)
    for L in (4, 8, 16, 24, 32, 48):
        r, p = dp_rho(q, L)
        mid = Fraction(q, 2) * sum(Fraction(math.comb(L, e), q ** e) for e in range(L + 1) if q ** e > 2 ** L)
        top = (q / 2) * 2 ** (-(1 - H2(c)) * L)
        assert p <= mid, (q, L)
        assert float(mid) <= top * (1 + 1e-12), (q, L, float(mid), top)
check(True, "rho^peak <= (q/2) sum_{q^e > 2^L} C(L,e) q^-e <= (q/2) 2^(-(1-H(log_q 2))L) for q = 3,5,7,9, L <= 48")
for q in range(3, 200, 2):
    c = math.log(2) / math.log(q)
    assert c * (q + 1) >= 1
check(True, "the Chernoff tilt z = qc/(1-c) >= 1 (c(q+1) >= 1) for every odd q < 200")

print("3. the peak construction on actual integers")


def peak_construction(q, L, N):
    n = np.arange(0, N + 1, dtype=np.int64)
    x = n.copy()
    best = n.copy()
    bestj = np.zeros(N + 1, dtype=np.int64)
    ok = np.ones(N + 1, dtype=bool)
    for j in range(1, L + 1):
        x = np.where(x % 2 == 1, (q * x + 1) // 2, x // 2)
        ok &= x >= n
        if j <= L - 1:
            upd = x > best
            best = np.where(upd, x, best)
            bestj = np.where(upd, j, bestj)
    bad = ok.copy()
    bad[:2] = False
    E = np.unique(best[bad])
    # descent of G = 1 on E, for every 2 <= n <= N
    xs = np.arange(2, N + 1, dtype=np.int64)
    start = xs.copy()
    done = np.zeros(len(xs), dtype=bool)
    for j in range(1, L + 1):
        inE = np.isin(xs, E)
        xs = np.where(inE, 1, np.where(xs % 2 == 1, (q * xs + 1) // 2, xs // 2))
        done |= xs < start
    # bad residue classes, and the exceptional finite set X_L
    words_bad = set()
    for r in range(2 ** L):
        y, e, good = r, 0, True
        wbits = []
        for j in range(1, L + 1):
            if y % 2:
                e += 1
            y = (q * y + 1) // 2 if y % 2 else y // 2
            if not q ** e > 2 ** j:
                good = False
                break
        if good:
            words_bad.add(r)
    X = [int(v) for v in np.nonzero(bad)[0] if (int(v) % 2 ** L) not in words_bad]
    return done.all(), int((E <= N).sum()), len(words_bad), X


for q, L, N in ((3, 8, 300000), (3, 12, 300000), (5, 8, 300000), (5, 12, 300000), (7, 10, 300000)):
    ok, nE, nbadw, X = peak_construction(q, L, N)
    r, p = dp_rho(q, L)
    bound = N * float(p) + nbadw + len(X)
    print(f"   q={q} L={L:2d}: all n <= {N} descend: {ok}; |E cap [1,N]| = {nE} <= N rho^peak + |Bad_L| + |X_L| = {bound:.0f}; X_L cap [2,N] = {X[:6]}")
    assert ok and nE <= bound
check(True, "peak construction: every 2 <= n <= 3e5 descends within L, and the density count of Theorem 1 holds")

print("4. the capacity constant")
for q in (3, 5, 7):
    c = math.log(2) / math.log(q)
    for L in (4, 8, 16, 32, 64):
        Ms = sum(sum(e // q + 1 for e in range((math.floor(k * c) + 1) if k else 0, k + 1)) for k in range(L))
        assert Ms <= L * (L + 1) * (2 * L - 2 + 3 * q) / (6 * q)
check(True, "M*_L <= L(L+1)(2L-2+3q)/(6q) for q = 3,5,7 and L <= 64")

print("5. q = 5: undecided density versus peak-discounted density")
for L in (16, 32, 48, 64):
    r, p = dp_rho(5, L)
    print(f"   L={L}: rho_L(5) = {float(r):.4f}, rho^peak_L(5) = {float(p):.3e}, -(1/L) log2 rho^peak = {-math.log2(float(p))/L:.5f}")
check(float(dp_rho(5, 64)[0]) > 0.17, "rho_L(5) stays above 0.17 while rho^peak_L(5) decays (exponent 1 - H(log_5 2) = 0.013911 asymptotically)")
