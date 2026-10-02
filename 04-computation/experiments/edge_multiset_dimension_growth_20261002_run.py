#!/usr/bin/env python3
"""edge_multiset_dimension_growth_20261002_run.py -- edge multiset dimension of hypercubes II: ln edim_m(Q_d) = Theta(d^(1/3)).
Note: 05-knowledge/results/edge_multiset_dimension_growth_20261002.md

Sections
  S1  pair types and cell matrices against brute force (d = 3, 4, 5)
  S2  atom bounds (Lemma 2): exact atoms, sqrt(pi/(8x)), closed form (2), series enclosure, Taylor step
  S3  Lemma H, exactly, for every 2 <= L <= 60
  S4  Theorem A chain: toward-the-centre forest products dominate g(L, Lambda)
  S5  certified sparse table, 11 <= d <= 64 (stored certificates re-verified in interval arithmetic)
  S6  Theorem C: explicit union bound with lambda = 5 d^(1/3), every 64 <= d <= 10^5 (interval arithmetic)
  S7  Theorem C beyond 10^5: the hand estimates
  S8  constants kappa_+ and kappa_-
  S9  Theorem B: the entropy bound L5 evaluated exactly (sanity)
  S10 Q_7 (section 7): the two resolving 19-sets, orbit counts, and the deposited exhaustive-search records
      of 04-computation/edge_multiset_dimension_q7_20261002/ (k <= 12); with --q7 the C search is rebuilt
      and re-run on the fast subset (Q_6, k <= 13; Q_7, k <= 10) in a temporary directory

Run:  nice python3 -u 04-computation/experiments/edge_multiset_dimension_growth_20261002_run.py [--q7]
Prints results to stdout (deterministic) and timings to stderr; ends with ALL CHECKS PASSED.
"""
import itertools
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from fractions import Fraction
from math import comb

import numpy as np
from mpmath import iv, mp, mpf
from scipy.special import i0e
from scipy.stats import binom

iv.dps = 30
mp.dps = 30
NCHECK = 0
T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))


def check(cond, msg):
    global NCHECK
    if not cond:
        print('CHECK FAILED:', msg, flush=True)
        sys.exit(1)
    NCHECK += 1


def hdr(s):
    print('\n' + '=' * 78 + '\n' + s + '\n' + '=' * 78, flush=True)
    print('[%7.1fs] %s' % (time.time() - T0, s.split(':')[0]), file=sys.stderr, flush=True)


# ------------------------------------------------------------------------------- pair types and cells
def type_list(d):
    """Aut(Q_d)-types of unordered pairs of distinct edges: (name, h, count, {(a,b): cell size}), a != b."""
    T = []
    for h in range(1, d):
        nu = d - 1 - h
        M = {}
        for c in range(nu + 1):
            for j in range(h + 1):
                a, b = c + j, c + h - j
                if a != b:
                    M[(a, b)] = M.get((a, b), 0) + 2 * comb(nu, c) * comb(h, j)
        T.append(('P', h, d * 2 ** (d - 2) * comb(d - 1, h), M))
    for h in range(0, d - 1):
        nu = d - 2 - h
        M = {}
        for c in range(nu + 1):
            for j in range(h + 1):
                for x in (0, 1):
                    for y in (0, 1):
                        a, b = x + c + j, y + c + h - j
                        if a != b:
                            M[(a, b)] = M.get((a, b), 0) + comb(nu, c) * comb(h, j)
        T.append(('X', h, d * (d - 1) * 2 ** (d - 1) * comb(d - 2, h), M))
    return T


def edge_weights(M):
    W = {}
    for (a, b), n in M.items():
        k = (min(a, b), max(a, b))
        W[k] = W.get(k, 0) + n
    return W


def kruskal(W, d):
    par = list(range(d + 1))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    F = []
    for (a, b), n in sorted(W.items(), key=lambda kv: (-kv[1], kv[0])):
        ra, rb = f(a), f(b)
        if ra != rb:
            par[ra] = rb
            F.append(((a, b), n))
    return F


def s1():
    hdr('S1: pair types and cell matrices against brute force')
    for d in (3, 4, 5):
        V = range(1 << d)
        E = [(u, i) for u in V for i in range(d) if not (u >> i) & 1]
        pc = [bin(x).count('1') for x in range(1 << d)]
        D = [[pc[(u ^ s) & ~(1 << i)] for s in V] for (u, i) in E]       # projection lemma
        tally = Counter()
        for x in range(len(E)):
            for y in range(x + 1, len(E)):
                M = Counter((a, b) for a, b in zip(D[x], D[y]) if a != b)
                k1 = tuple(sorted(M.items()))
                k2 = tuple(sorted(((b, a), n) for (a, b), n in M.items()))
                tally[min(k1, k2)] += 1
        mine = Counter()
        for name, h, cnt, M in type_list(d):
            k1 = tuple(sorted(M.items()))
            k2 = tuple(sorted(((b, a), n) for (a, b), n in M.items()))
            mine[min(k1, k2)] += cnt
        check(tally == mine, 'S1 brute force d=%d' % d)
        Etot = d * 2 ** (d - 1)
        check(sum(cnt for _, _, cnt, _ in type_list(d)) == Etot * (Etot - 1) // 2, 'S1 counts d=%d' % d)
        print('d=%d: %d edge pairs, %d cell classes; formulas match brute force' % (d, Etot * (Etot - 1) // 2, len(tally)))
    for d in (6, 10, 20, 64):
        Etot = d * 2 ** (d - 1)
        check(sum(cnt for _, _, cnt, _ in type_list(d)) == Etot * (Etot - 1) // 2, 'S1 counts d=%d' % d)
    print('type counts add up to C(E,2) for d = 6, 10, 20, 64')


# ---------------------------------------------------------------------------------------- atom bounds
C4 = iv.mpf(75) / 1024 * iv.mpf(2) ** (iv.mpf(7) / 2)


def atom_closed_iv(x):
    return (1 + 1 / (8 * x) + iv.mpf(9) / (128 * x * x) + C4 / (x * x * x)) / iv.sqrt(2 * iv.pi * x) + iv.exp(-x) / 2


def atom_series_iv(x, K=None):
    """enclosure of e^{-x} I_0(x) from I_0(x) = sum (x^2/4)^k/(k!)^2 with a geometric tail bound."""
    K = int(float(x.b)) + 40 if K is None else K
    y = x * x / 4
    term = iv.mpf(1)
    s = iv.mpf(1)
    for k in range(1, K + 1):
        term = term * y / (k * k)
        s += term
    nxt = term * y / ((K + 1) * (K + 1))
    r = y / ((K + 2) * (K + 2))
    assert r.b < 1
    return iv.mpf([s.a, (s + nxt / (1 - r)).b]) * iv.exp(-x)


def log_atom_upper(x_iv):
    xlo = x_iv.a
    if xlo <= 0:
        return iv.mpf(0)
    X = iv.mpf(xlo)
    best = atom_series_iv(X).b if xlo < 60 else atom_closed_iv(X).b
    if best >= 1:
        return iv.mpf(0)
    return iv.log(iv.mpf(best))


def s2():
    hdr('S2: atom bounds (Lemma 2)')
    for N1, N2, q in [(3, 5, 0.3), (10, 10, 0.5), (1, 1, 0.5), (7, 0, 0.1), (50, 80, 0.02), (200, 100, 0.4), (400, 1, 0.01)]:
        p1 = binom.pmf(np.arange(N1 + 1), N1, q)
        p2 = binom.pmf(np.arange(N2 + 1), N2, q)
        a = np.convolve(p1, p2[::-1]).max()
        x = (N1 + N2) * q * (1 - q)
        check(a <= float(i0e(x)) + 1e-15, 'S2 exact atom (%d,%d,%g)' % (N1, N2, q))
    print('exact maximal atoms of Bin(N1,q) - Bin(N2,q) <= e^-x I_0(x) in 7 cases')
    worst_closed = 0.0
    for x in [0.01, 0.1, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 8, 13, 20, 30, 50, 59.9]:
        X = iv.mpf(x)
        enc = atom_series_iv(X)
        check(enc.a <= i0e(x) * (1 + 1e-12) and enc.b >= i0e(x) * (1 - 1e-12), 'S2 series enclosure %g' % x)
        check(math.sqrt(math.pi / (8 * x)) >= enc.b, 'S2 sqrt bound %g' % x)
        cl = atom_closed_iv(X)
        check(cl.a >= enc.b, 'S2 closed form (2) %g' % x)
        worst_closed = max(worst_closed, float(cl.b) / float(enc.a))
    for x in [60, 100, 1e3, 1e4, 1e5]:
        check(atom_closed_iv(iv.mpf(x)).a >= i0e(x), 'S2 closed form large %g' % x)
        check(math.sqrt(math.pi / (8 * x)) >= i0e(x), 'S2 sqrt bound large %g' % x)
    rel30 = float(atom_closed_iv(iv.mpf(30)).b) / float(atom_series_iv(iv.mpf(30)).a) - 1
    rel100 = float(atom_closed_iv(iv.mpf(100)).b) / float(i0e(100)) - 1
    print('series enclosure and both closed bounds checked on a grid x in [0.01, 1e5]; '
          'bound (2) exceeds the true value by %.1e (x=30) and %.1e (x=100)' % (rel30, rel100))
    # Taylor step: (1-y)^(-1/2) <= 1 + y/2 + 3y^2/8 + (5/16) 2^(7/2) y^3 on [0, 1/2]
    c3 = 5 / 16 * 2 ** 3.5
    ys = np.linspace(0, 0.5, 5001)
    check(np.all((1 - ys) ** -0.5 <= 1 + ys / 2 + 3 * ys ** 2 / 8 + c3 * ys ** 3 + 1e-12), 'S2 Taylor step')
    check(abs(float(C4.a) - 75 / 1024 * 2 ** 3.5) < 1e-12 and abs(float(C4.a) - 0.8286) < 1e-4, 'S2 C4')
    print('Taylor remainder step checked on [0, 1/2]; C_4 = (75/1024) 2^(7/2) = %.4f' % float(C4.a))


# ------------------------------------------------------------------------------------------- Lemma H
def s3():
    hdr('S3: Lemma H for every 2 <= L <= 60 (exact rational arithmetic)')
    worst = None
    cases = 0
    for L in range(2, 61):
        CL = [comb(L, a) for a in range(L + 1)]
        for h in range(1, L):
            nu = L - h
            for a in range(L + 1):
                t2 = 2 * a - L
                if abs(t2) <= 2:
                    continue
                aa = a if t2 > 0 else L - a
                t = Fraction(2 * aa - L, 2)
                best = Fraction(0)
                for k in range(max(0, aa - nu), min(h, aa) + 1):
                    if Fraction(h, 2) < k < Fraction(h, 2) + t:
                        p = Fraction(comb(h, k) * comb(nu, aa - k), CL[aa])
                        if p > best:
                            best = p
                bound = Fraction(1, 3 * (min(h, nu) + 1))
                check(best >= bound, 'S3 Lemma H (%d,%d,%d)' % (L, h, a))
                cases += 1
                r = best / bound
                if worst is None or r < worst[0]:
                    worst = (r, L, h, a)
    print('Lemma H holds in all %d cases; smallest ratio best/bound = %.4f at (L,h,a) = %s'
          % (cases, float(worst[0]), worst[1:]))


# ----------------------------------------------------------------------------------- Theorem A chain
def g_gen_f(L, Lam):
    tau = math.sqrt(L * Lam / 2)
    if tau < 3 or 2 * Lam >= L:
        return None
    return (2 / 3) * Lam * tau - 3 * Lam - (tau / 15 + 2 / 3) * Lam ** 2 / (L - 2 * Lam)


def s4():
    hdr('S4: Theorem A chain -- toward-the-centre forest sums dominate g(L, Lambda)')
    worst = None
    n = 0
    for d in [20, 30, 40, 64, 100, 150]:
        for C in [3, 4, 5, 6]:
            lam = C * d ** (1 / 3)
            if lam > (d - 2) * math.log(2):
                continue
            q = math.exp(lam - d * math.log(2))
            m = math.exp(lam)
            for typ in 'PX':
                L, kap = (d - 1, 2.0) if typ == 'P' else (d - 2, 0.5)
                for h in (range(1, d - 1) if typ == 'P' else range(1, d - 2)):
                    nu = L - h
                    tot = 0.0
                    used = set()
                    for a in range(L + 1):
                        t2 = 2 * a - L
                        if abs(t2) <= 2:
                            continue
                        aa = a if t2 > 0 else L - a
                        t = (2 * aa - L) / 2
                        best, bk = 0, None
                        for k in range(max(0, aa - nu), min(h, aa) + 1):
                            if h / 2 < k < h / 2 + t:
                                p = comb(h, k) * comb(nu, aa - k)
                                if p > best:
                                    best, bk = p, k
                        x = kap * (1 - q) * m * best / 2 ** L
                        b = int(round(aa - 2 * (bk - h / 2)))
                        e = tuple(sorted((a, b if t2 > 0 else L - b)))
                        check(e not in used, 'S4 distinct forest edges')
                        used.add(e)
                        tot += max(0.0, 0.5 * math.log(8 * x / math.pi))
                    Lam = lam - math.log(3 * math.pi * (L + 1) ** 2 / (8 * kap * (1 - q)))
                    gg = g_gen_f(L, Lam)
                    if gg is None:
                        continue
                    check(tot >= gg - 1e-9, 'S4 chain d=%d C=%d %s%d' % (d, C, typ, h))
                    n += 1
                    if worst is None or tot - gg < worst[0]:
                        worst = (tot - gg, d, C, typ, h)
    print('%d (d, C, type) cases; smallest slack %.1f nats at d=%d, C=%d, type %s(%d)' % ((n,) + worst))
    # the special types: X(0) (path forest, Lambda_X) and the antipodal types P(d-1), X(d-2) (matchings, g')
    worst2 = None
    n2 = 0
    for d in [20, 30, 40, 64, 100, 150]:
        for C in [3, 4, 5, 6]:
            lam = C * d ** (1 / 3)
            if lam > (d - 2) * math.log(2):
                continue
            q = math.exp(lam - d * math.log(2))
            m = math.exp(lam)
            lamX = lam - math.log(math.pi * (d - 1) / (4 * (1 - q)))
            lamP = lam - math.log(math.pi * d / (16 * (1 - q)))
            # X(0): levels 0..d-1, edge {c, c+1} with x = (1-q) m b_(d-2)(c)/2, each level takes the edge toward (d-1)/2
            tot = 0.0
            for a in range(d):
                s2 = 2 * a - (d - 1)
                if abs(s2) <= 2:
                    continue
                c = a - 1 if s2 > 0 else a
                x = (1 - q) * m * comb(d - 2, c) / 2 ** (d - 2) / 2
                tot += max(0.0, 0.5 * math.log(8 * x / math.pi))
            gg = g_gen_f(d - 2, lamX)
            if gg is not None:
                check(tot >= gg - 1e-9, 'S4 X(0) d=%d C=%d' % (d, C))
                n2 += 1
                if worst2 is None or tot - gg < worst2[0]:
                    worst2 = (tot - gg, d, C, 'X(0)')
            # antipodal matchings: P(d-1) with x = 2(1-q) m b_L(j), L = d-1; X(d-2) with x = (1-q) m b_L(j)/2, L = d-2
            for name, L, kap, lamM in (('P(d-1)', d - 1, 2.0, lamP), ('X(d-2)', d - 2, 0.5, lamX)):
                tot = 0.0
                for j in range(L + 1):
                    if 2 * j - L >= -2:
                        continue
                    x = kap * (1 - q) * m * comb(L, j) / 2 ** L
                    tot += max(0.0, 0.5 * math.log(8 * x / math.pi))
                tau = math.sqrt(L * lamM / 2)
                if tau < 3 or 2 * lamM >= L:
                    continue
                gm = (1 / 3) * lamM * tau - lamM - (tau / 15 + 2 / 3) * lamM ** 2 / (2 * (L - 2 * lamM))
                check(tot >= gm - 1e-9, 'S4 %s d=%d C=%d' % (name, d, C))
                n2 += 1
                if tot - gm < worst2[0]:
                    worst2 = (tot - gm, d, C, name)
    print('special types X(0), P(d-1), X(d-2): %d cases; smallest slack %.1f nats at d=%d, C=%d, %s' % ((n2,) + worst2))


# --------------------------------------------------------------------------------- certified table
CERTS = [
    (11, 222861, 1024000, 492), (12, 494559, 4096000, 543), (13, 5681, 81920, 622),
    (14, 10403, 256000, 724), (15, 155287, 6553600, 843), (16, 90673, 6553600, 978),
    (17, 1052691, 131072000, 1129), (18, 1223837, 262144000, 1306), (19, 1400923, 524288000, 1490),
    (20, 1614059, 1048576000, 1708), (21, 459133, 524288000, 1936), (22, 2082059, 4194304000, 2193),
    (23, 1176267, 4194304000, 2464), (24, 331207, 2097152000, 2774), (25, 1481233, 16777216000, 3095),
    (26, 1656007, 33554432000, 3454), (27, 1845761, 67108864000, 3833), (28, 1, 65536, 4254),
    (29, 22599, 2684354560, 4688), (30, 5004947, 1073741824000, 5179), (31, 5497383, 2147483648000, 5681),
    (32, 6046957, 4294967296000, 6242), (33, 52939, 68719476736, 6823), (34, 7246457, 17179869184000, 7463),
    (35, 7901081, 34359738368000, 8127), (36, 8621929, 68719476736000, 8857), (37, 9358381, 137438953472000, 9604),
    (38, 203499, 5497558138880, 10431), (39, 11017147, 549755813888000, 11283), (40, 5969909, 549755813888000, 12216),
    (41, 2576067, 439804651110400, 13172), (42, 6956827, 2199023255552000, 14221), (43, 14973433, 8796093022208000, 15294),
    (44, 197, 214748364800, 16469), (45, 17331657, 35184372088832000, 17674), (46, 18622661, 70368744177664000, 18978),
    (47, 4987297, 35184372088832000, 20317), (48, 1069591, 14073748835532800, 21777), (49, 5720211, 140737488355328000, 23268),
    (50, 24485831, 1125899906842624000, 24884), (51, 26096877, 2251799813685248000, 26536), (52, 27884361, 4503599627370496000, 28320),
    (53, 7429391, 2251799813685248000, 30150), (54, 3988127, 2251799813685248000, 32262), (55, 6754181, 7205759403792793600, 34186),
    (56, 35812997, 72057594037927936000, 36308), (57, 9586113, 36028797018963968000, 38726), (58, 2017107, 14411518807585587200, 40911),
    (59, 42805259, 576460752303423488000, 43347), (60, 45414131, 1152921504606846976000, 45966), (61, 48263011, 2305843009213693952000, 48740),
    (62, 639839, 57646075230342348800, 51658), (63, 13446593, 2305843009213693952000, 54415), (64, 28397833, 9223372036854775808000, 57738),
]


def U_iv(P, q):
    qq = iv.mpf(q.numerator) / q.denominator
    tot = iv.mpf(0)
    for name, h, cnt, Ns, antip in P:
        lp = iv.mpf(0)
        for n in Ns:
            lp += log_atom_upper(qq * (1 - qq) * n)
        tot += (iv.mpf(cnt) if antip else iv.mpf(cnt) / 2) * iv.exp(lp)
    return tot


def log_tail_iv(n, q, M):
    qq = iv.mpf(q.numerator) / q.denominator
    p = iv.mpf(M) / n
    return -n * (p * iv.log(p / qq) + (1 - p) * iv.log((1 - p) / (1 - qq)))


def s5():
    hdr('S5: certified sparse table, 11 <= d <= 64')
    print('  d  q (rational)                          M_d     U <=        tail <=     lnM/d^(1/3)')
    vals = []
    for d, qn, qd, M in CERTS:
        q = Fraction(qn, qd)
        P = []
        for name, h, cnt, Mc in type_list(d):
            F = kruskal(edge_weights(Mc), d)
            P.append((name, h, cnt, [n for _, n in F], name == 'P' and h == d - 1))
        U = U_iv(P, q)
        check(M + 1 > 2 ** d * q, 'S5 Chernoff regime d=%d' % d)
        tail = iv.exp(log_tail_iv(2 ** d, q, M + 1))
        check((U + tail).b < 1, 'S5 certificate d=%d' % d)
        r = math.log(M) / d ** (1 / 3)
        vals.append(r)
        print('%3d  %-36s %6d  %.8f  %.8f  %.4f' % (d, '%d/%d' % (qn, qd), M, float(U.b), float(tail.b), r))
    check(len(CERTS) == 54 and [c[0] for c in CERTS] == list(range(11, 65)), 'S5 range')
    check(2.73 <= min(vals) and max(vals) <= 2.79, 'S5 ratio range')
    print('every certificate verified: edim_m(Q_d) <= M_d; lnM/d^(1/3) in [%.4f, %.4f]' % (min(vals), max(vals)))


# --------------------------------------------------------------------------------------- Theorem C
LN2 = iv.log(2)


def g_gen(L, Lam):
    tau = iv.sqrt(L * Lam / 2)
    if tau.a < 3 or (2 * Lam).b >= L:
        return None
    return (iv.mpf(2) / 3) * Lam * tau - 3 * Lam - (tau / 15 + iv.mpf(2) / 3) * Lam ** 2 / (L - 2 * Lam)


def g_match(L, Lam):
    tau = iv.sqrt(L * Lam / 2)
    if tau.a < 3 or (2 * Lam).b >= L:
        return None
    return (iv.mpf(1) / 3) * Lam * tau - Lam - (tau / 15 + iv.mpf(2) / 3) * Lam ** 2 / (2 * (L - 2 * Lam))


def theorem_c_total(d, C=iv.mpf(5)):
    lam = C * iv.mpf(d) ** (iv.mpf(1) / 3)
    lq = lam - d * LN2
    if lq.b > -iv.log(4).a:
        return None
    q = iv.exp(lq)
    omq = 1 - q
    gs = []
    for L in (d - 1, d - 2):
        for kap in (iv.mpf(2), iv.mpf(1) / 2):
            g = g_gen(L, lam - iv.log(3 * iv.pi * (L + 1) ** 2 / (8 * kap * omq)))
            if g is None:
                return None
            gs.append(g)
    g = g_gen(d - 2, lam - iv.log(iv.pi * (d - 1) / (4 * omq)))
    if g is None:
        return None
    gs.append(g)
    gmin = min(gs, key=lambda z: z.a)
    gP = g_match(d - 1, lam - iv.log(iv.pi * d / (16 * omq)))
    gX = g_match(d - 2, lam - iv.log(iv.pi * (d - 1) / (4 * omq)))
    if gP is None or gX is None:
        return None
    lE = iv.log(d) + (d - 1) * LN2
    U = (iv.exp(2 * lE - iv.log(4) - gmin) + iv.exp(iv.log(d) + (d - 2) * LN2 - gP)
         + iv.exp(iv.log(d * (d - 1)) + (d - 2) * LN2 - gX))
    return U + iv.exp(-iv.exp(lam) / 3)


def s6():
    hdr('S6: Theorem C, lambda = 5 d^(1/3), every 64 <= d <= 10^5 (interval arithmetic)')
    worst = (0, None)
    for d in range(64, 100001):
        tot = theorem_c_total(d)
        check(tot is not None and tot.b < 1, 'S6 d=%d' % d)
        if tot.b > worst[0]:
            worst = (float(tot.b), d)
    check(worst[1] == 64, 'S6 worst at 64')
    print('U + e^(-m/3) < 1 for every d in [64, 100000]; largest value %.4f at d = %d' % worst)


def s7():
    hdr('S7: Theorem C beyond 10^5 -- the hand estimates')
    f = lambda d: 0.6 * d ** (1 / 3) - 2 * math.log(d) - 0.87
    check(f(1e5) > 0, 'S7 positivity at 1e5')
    check(all(0.2 * d ** (-2 / 3) - 2 / d > 0 for d in (1001, 1e4, 1e5, 1e8, 1e12)), 'S7 increasing')
    check(abs(math.log(3 * math.pi / 4) - 0.857) < 0.01, 'S7 constant 0.87')
    for d in [10 ** 5, 10 ** 6, 10 ** 8, 10 ** 12]:
        D = iv.mpf(d)
        lo, hi = 4.4 * D ** (iv.mpf(1) / 3), 5 * D ** (iv.mpf(1) / 3)
        tau_hi = iv.sqrt(D * hi / 2)
        check(tau_hi.b <= 1.6 * float(D.a) ** (2 / 3), 'S7 tau bound d=%g' % d)
        main = iv.sqrt(2) / 3 * iv.sqrt(D - 2) * lo ** iv.mpf(1.5)
        neg = 3 * hi + (tau_hi / 15 + iv.mpf(2) / 3) * hi ** 2 / (iv.mpf(0.99) * D)
        check((main - neg).a >= 4.2 * d, 'S7 g_* >= 4.2 d at d=%g' % d)
        mainp = iv.sqrt(2) / 6 * iv.sqrt(D - 2) * lo ** iv.mpf(1.5)
        negp = hi + (tau_hi / 15 + iv.mpf(2) / 3) * hi ** 2 / (2 * iv.mpf(0.99) * D)
        check((mainp - negp).a >= 2.1 * d, "S7 g' >= 2.1 d at d=%g" % d)
        tot = (iv.exp(2 * D * LN2 + 2 * iv.log(D) - 4.2 * D) + 2 * iv.exp(D * LN2 + 2 * iv.log(D) - 2.1 * D))
        check((iv.log(tot)).b < -1.4 * d, 'S7 U < e^(-1.4 d) at d=%g' % d)
    print('hand estimates of Theorem C checked at d = 1e5, 1e6, 1e8, 1e12 (all terms monotone in d)')


def s8():
    hdr('S8: constants')
    kp = (3 * mp.sqrt(2) * mp.log(2)) ** (mpf(2) / 3)
    km = kp / mpf(4) ** (mpf(2) / 3)
    check(abs((mp.sqrt(2) / 3) * kp ** 1.5 - 2 * mp.log(2)) < mpf(10) ** -25, 'S8 kappa_+')
    check(abs((2 * mp.sqrt(2) / 3) * km ** 1.5 - mp.log(2)) < mpf(10) ** -25, 'S8 kappa_-')
    check(abs(km - ((3 * mp.sqrt(2) / 4) * mp.log(2)) ** (mpf(2) / 3)) < mpf(10) ** -25, 'S8 kappa_- form')
    print('kappa_+ = (3 sqrt2 ln2)^(2/3) = %s; kappa_- = ((3 sqrt2/4) ln2)^(2/3) = %s; ratio 4^(2/3)'
          % (mp.nstr(kp, 12), mp.nstr(km, 12)))


def g_sum(ell):
    """sum of G(mu) = (mu+1) ln(mu+1) - mu ln mu over mu = exp(ell), computed without cancellation"""
    big = ell > 0
    x = np.exp(-ell[big])
    lp = np.log1p(x)
    ratio = np.where(x > 1e-12, lp / np.where(x > 1e-12, x, 1.0), 1 - x / 2)
    mu = np.exp(ell[~big])
    return (ell[big] + lp + ratio).sum() + ((1 + mu) * np.log1p(mu) - mu * ell[~big]).sum()


def s9():
    hdr('S9: Theorem B sanity -- the entropy bound L5 evaluated exactly')
    from scipy.special import gammaln
    for t in [-700.0, -30.0, -2.0, -0.1, 0.0, 0.1, 2.0, 30.0, 700.0]:
        mu = mp.exp(t)
        exact = mp.log1p(mu) + mu * mp.log1p(1 / mu)          # = G(mu), in 30-digit arithmetic
        check(abs(g_sum(np.array([t])) - float(exact)) <= 1e-13 * float(exact), 'S9 stable G at %g' % t)
    km = ((3 * math.sqrt(2) / 4) * math.log(2)) ** (2 / 3)
    rows = []
    for d in [50, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6]:
        L = d - 1
        r = np.arange(L + 1)
        lb = gammaln(L + 1) - gammaln(r + 1) - gammaln(L - r + 1) - L * math.log(2)
        target = (d - 1) * math.log(2) + math.log(d)
        lo, hi = 0.0, 5 * d ** (1 / 3) + 10
        check(g_sum(hi + lb) > target and g_sum(lo + lb) < target, 'S9 bracket d=%d' % d)
        for _ in range(100):
            mid = (lo + hi) / 2
            if g_sum(mid + lb) >= target:
                hi = mid
            else:
                lo = mid
        rows.append((d, hi, hi / d ** (1 / 3)))
    for d, lam, r in rows:
        print('d = %7d: L5 forces ln|S| >= %8.4f = %.4f d^(1/3)   (minus (1/2) ln(pi d/2): %.4f d^(1/3))'
              % (d, lam, r, (lam - 0.5 * math.log(math.pi * d / 2)) / d ** (1 / 3)))
    check(all(r > km for _, _, r in rows), 'S9 above kappa_-')
    check(all(rows[i][2] > rows[i + 1][2] for i in range(1, len(rows) - 1)), 'S9 decreasing from d = 100')
    print('the ratio decreases toward kappa_- = %.4f; most of the excess is the factor sqrt(2/(pi L)) in b_L' % km)


Q7_PKG = os.path.normpath(os.path.join(HERE, '..', 'edge_multiset_dimension_q7_20261002'))
Q7_RERUN = '--q7' in sys.argv
Q7_KNOWN19 = [5, 57, 54, 104, 35, 109, 115, 49, 6, 39, 55, 102, 15, 32, 85, 97, 21, 113, 41]
Q7_NEW19 = [2, 4, 21, 22, 32, 38, 42, 44, 47, 54, 70, 79, 84, 90, 110, 114, 116, 120, 126]
Q6_ORBITS = [1, 1, 6, 16, 103, 497, 3253, 19735, 120843, 681474, 3561696, 16938566]      # a = 0..11
Q7_LEAVES = [None, 1, 64, 384, 12154, 32283, 686390, 4300936, 68340374, 317832905, 4141968201,
             25086535792, 273963417196]                                                 # k = 1..12
Q6_LEAVES = [None, 1, 32, 160, 2506, 4971, 52535, 234240, 1808073, 4767589, 29955834, 96667005,
             491865822, 1234172511]                                                     # k = 1..13


def hist_defect(d, S):
    """#edges - #distinct histograms, straight from the definition (0 iff S is edge-multiset resolving)"""
    H = set()
    m = 0
    for u in range(1 << d):
        for i in range(d):
            if not (u >> i) & 1:
                v = u | (1 << i)
                h = [0] * d
                for s in S:
                    h[min(bin(u ^ s).count('1'), bin(v ^ s).count('1'))] += 1
                H.add(tuple(h))
                m += 1
    return m - len(H)


def q_burnside(n, amax):
    """orbits of Aut(Q_n) on the a-subsets of V(Q_n), a = 0..amax, by Burnside over all 2^n n! elements"""
    types = Counter()
    for pi in itertools.permutations(range(n)):
        P = [sum(((x >> i) & 1) << pi[i] for i in range(n)) for x in range(1 << n)]
        for t in range(1 << n):
            seen = 0
            cyc = []
            for x in range(1 << n):
                if not (seen >> x) & 1:
                    L, y = 0, x
                    while not (seen >> y) & 1:
                        seen |= 1 << y
                        y = P[y ^ t]
                        L += 1
                    cyc.append(L)
            types[tuple(sorted(cyc))] += 1
    total = [0] * (amax + 1)
    for cyc, mult in types.items():
        poly = [1] + [0] * amax
        for L in cyc:
            for a in range(amax, L - 1, -1):
                poly[a] += poly[a - L]
        for a in range(amax + 1):
            total[a] += mult * poly[a]
    G = math.factorial(n) << n
    check(all(x % G == 0 for x in total), 'S10 Burnside integrality')
    return [x // G for x in total]


def w_min(n, a):
    """least total weight of a distinct vertices of Q_n (lightest first)"""
    return sum(sorted(bin(x).count('1') for x in range(1 << n))[:a])


def q7_check_records(rundir, D, kmax, leaves, nrep_ok):
    """the driver/validate records of one run directory: tiling, C leaves = DP leaves, nothing found"""
    n = D - 1
    for k in range(1, kmax + 1):
        V = json.load(open(os.path.join(rundir, 'validation_k%d.json' % k)))
        check(V['D'] == D and V['k'] == k and V['ok'] and not V['problems'], 'S10 D=%d k=%d validation ok' % (D, k))
        check(V['resolving_leaves'] == 0, 'S10 D=%d k=%d nothing found' % (D, k))
        check(V['leaves_total_C'] == V['leaves_total_DP'] == leaves[k], 'S10 D=%d k=%d leaves' % (D, k))
        check(set(V['per_a']) == set(str(a) for a in range((k + 1) // 2, k + 1)), 'S10 D=%d k=%d a-range' % (D, k))
        for a in range((k + 1) // 2, k + 1):
            r, b = V['per_a'][str(a)], k - a
            if r.get('excluded_by_weight_lemma'):
                check(a - 2 * b > 0 and w_min(n, a) > n * b, 'S10 D=%d k=%d a=%d weight lemma' % (D, k, a))
            else:
                check(r['tiled'] and r['found'] == 0 and r['leaves_C'] == r['leaves_DP'], 'S10 D=%d k=%d a=%d' % (D, k, a))
                check(nrep_ok(a, r['nrep']), 'S10 D=%d k=%d a=%d representative count' % (D, k, a))


def q7_rerun():
    """rebuild the C search and re-run the fast subset in a temporary directory (about 4 min on 2 cores)"""
    work = tempfile.mkdtemp(prefix='edim_q7_')
    try:
        for sub in ('src', 'py'):
            shutil.copytree(os.path.join(Q7_PKG, sub), os.path.join(work, sub),
                            ignore=shutil.ignore_patterns('__pycache__', '*.o'))
        for sub in ('data/q5', 'data/q6', 'runs'):
            os.makedirs(os.path.join(work, sub))

        def sh(cmd):
            r = subprocess.run(cmd, cwd=work, shell=True, capture_output=True, text=True)
            check(r.returncode == 0, 'S10 --q7 command failed: %s\n%s' % (cmd, r.stdout[-1500:] + r.stderr[-1500:]))
            return r.stdout
        for prog in ('edimsearch', 'orbreps', 'canon'):
            sh('gcc -O3 -march=native -Wall -o src/%s src/%s.c' % (prog, prog))
        sh('./src/orbreps 5 32 data/q5')
        sh('./src/orbreps 6 10 data/q6')
        for D, kmax, leaves in ((6, 13, Q6_LEAVES), (7, 10, Q7_LEAVES)):
            rundir = os.path.join(work, 'runs', 'r%d' % D)
            for k in range(1, kmax + 1):
                sh('python3 py/driver.py %d %d --workers 2 --target 60 --engine 0 --run %s' % (D, k, rundir))
            out = sh('python3 py/validate.py %d %s %s' % (D, rundir, ' '.join(str(k) for k in range(1, kmax + 1))))
            good = [l for l in out.splitlines() if l.startswith('VALIDATE D=%d ' % D) and ' ok=True ' in l]
            check(len(good) == kmax, 'S10 --q7 validate D=%d' % D)
            q7_check_records(rundir, D, kmax, leaves, lambda a, nrep: True)
            print('--q7 rerun: Q_%d, k = 1..%d: %d leaves, C = DP, nothing found' % (D, kmax, sum(leaves[1:kmax + 1])))
    finally:
        shutil.rmtree(work, ignore_errors=True)


def s10():
    hdr('S10: Q_7 -- section 7 (pure-Python checks and the deposited exhaustive-search records)')
    # the two resolving 19-sets, from the definition; deleting any landmark of the known set breaks it
    check(hist_defect(7, Q7_KNOWN19) == 0, 'S10 known 19-set resolves Q_7')
    check(hist_defect(7, Q7_NEW19) == 0, 'S10 new 19-set resolves Q_7')
    dels = [hist_defect(7, Q7_KNOWN19[:j] + Q7_KNOWN19[j + 1:]) for j in range(19)]
    check(min(dels) > 0, 'S10 known 19-set minus a landmark')
    print('the known 19-set and the annealed 19-set both resolve Q_7; deleting one landmark of the known set')
    print('leaves defects', dels)
    # Aut(Q_6)-orbits of a-subsets (the representative counts), and the weight lemma (Lemma 2 of section 7)
    orb = q_burnside(6, 11)
    check(orb == Q6_ORBITS, 'S10 Aut(Q_6) orbit counts')
    print('Aut(Q_6)-orbits of a-subsets, a = 0..11:', orb)
    wm = [w_min(6, a) for a in (11, 12, 13, 14)]
    check(wm == [14, 16, 18, 20], 'S10 W_min')
    check(all(w_min(6, a) > 6 * (k - a) and a - 2 * (k - a) > 0 for k in range(1, 14) for a in range(11, k + 1)),
          'S10 weight lemma: every a >= 11 excluded for k <= 13')
    print('W_min(11..14) =', wm, '-> for k <= 13 every a >= 11 is excluded (Lemma 2)')
    # the deposited records of the exhaustive Q_7 search, k = 1..12
    rundir = os.path.join(Q7_PKG, 'runs', 'd7_e0')
    q7_check_records(rundir, 7, 12, Q7_LEAVES, lambda a, nrep: nrep == Q6_ORBITS[a])
    tot = sum(Q7_LEAVES[1:])
    check(tot == 303583126680, 'S10 total leaves')
    print('deposited records, Q_7, k = 1..12: every (k, a) tiled, C leaves = DP leaves (%d in total),' % tot)
    print('representative counts = Burnside, nothing found => edim_m(Q_7) >= 13')
    if Q7_RERUN:
        q7_rerun()
    else:
        print('(run with --q7 to rebuild the C search and re-run Q_6, k <= 13, and Q_7, k <= 10)')


if __name__ == '__main__':
    for sec in (s1, s2, s3, s4, s5, s6, s7, s8, s9, s10):
        sec()
    print('\n%d checks' % NCHECK)
    print('[%7.1fs] done' % (time.time() - T0), file=sys.stderr)
    print('ALL CHECKS PASSED')
