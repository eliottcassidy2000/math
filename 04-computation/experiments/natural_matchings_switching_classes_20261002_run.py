#!/usr/bin/env python3
"""natural_matchings_switching_classes_20261002_run.py -- OPEN-Q-060: natural matchings of tournament switching classes with
untwisted Euler graphs, the Royle et al. setting, and the Z/l switching count.
Note: 05-knowledge/results/natural_matchings_switching_classes_20261002.md

Sections
  S1  Theorem 1: odd n, unique member with all scores = (n-1)/2 mod 2
  S2  Theorem 2: level lemma (GF(2)), twisting witness, twisted Euler graph counts n <= 13
  S3  Theorem 3: graphs, even n
  S4  Theorem 4: census of natural matchings, classes vs untwisted Euler graphs (n <= 8; --full: 9)
  S5  Theorems 8, 9: Royle setting census (n <= 7; --full: 8) and the n = 5 partners
  S6  Lemma 5.3: rigid blocks
  S7  Theorems 5, 6: forced classes and class keys (5 <= n <= 60; --full: 100)
  S8  Theorem 6: the arithmetic condition A(n) >= 2 for n <= 10^6
  S9  Theorem 6, n = 0 mod 4: balanced-pair counts
  S10 prime blocks, doubly regular Paley tournaments
  S11 Theorem 10: Royle-forced families
  S12 Theorem 11: switching classes vs modular Eulerian matrices over Z/l

Run:  nice python3 -u 04-computation/experiments/natural_matchings_switching_classes_20261002_run.py [--full]
Prints results to stdout (deterministic) and timings to stderr; ends with ALL CHECKS PASSED.
"""
import itertools
import math
import os
import sys
import time
from collections import Counter

import numpy as np
from sympy import factorint, isprime
from sympy.utilities.iterables import partitions

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import natural_matchings_switching_classes_20261002_lib as L  # noqa: E402

FULL = '--full' in sys.argv
NCHECK = 0
T0 = time.time()


def check(cond, msg):
    global NCHECK
    if not cond:
        print('CHECK FAILED:', msg, flush=True)
        sys.exit(1)
    NCHECK += 1


def hdr(s):
    print('\n' + '=' * 78 + '\n' + s + '\n' + '=' * 78, flush=True)
    print('[%7.1fs] %s' % (time.time() - T0, s.split(':')[0]), file=sys.stderr, flush=True)


# ------------------------------------------------------------------------------------------------- S1
def s1():
    hdr('S1: Theorem 1 -- odd n, the member with all scores = (n-1)/2 mod 2')
    for n in (3, 5, 7) + ((9,) if FULL else ()):
        c = (n - 1) // 2 % 2
        k = sum(1 for A in L.gentourng(n) if all(s % 2 == c for s in L.scores(A)))
        check(k == L.A049313[n], 'S1 count n=%d' % n)
        print('n=%d: tournaments (up to isomorphism) with all scores = %d mod 2: %d = A049313(%d)' % (n, c, k, n))
    # labelled form: every labelled switching class has exactly one such member (n = 5, 7 exhaustive)
    for n in (5, 7):
        c = (n - 1) // 2 % 2
        P = [(i, j) for i in range(n) for j in range(i + 1, n)]
        m = len(P)
        codes = np.arange(1 << m, dtype=np.int64)
        bits = ((codes[:, None] >> np.arange(m)) & 1).astype(np.int8)    # bit t = 1: P[t][0] -> P[t][1]
        sc = np.zeros((1 << m, n), dtype=np.int16)
        for t, (i, j) in enumerate(P):
            sc[:, i] += bits[:, t]
            sc[:, j] += 1 - bits[:, t]
        good = np.all(sc % 2 == c, axis=1)
        # class normal form: switch the in-neighbours of vertex 0 (vertex 0 becomes a source)
        inn = np.zeros((1 << m, n), dtype=np.int8)
        for t, (i, j) in enumerate(P):
            if i == 0:
                inn[:, j] = 1 - bits[:, t]
        nf = np.zeros(1 << m, dtype=np.int64)
        for t, (i, j) in enumerate(P):
            if i == 0:
                continue
            b = bits[:, t] ^ (inn[:, i] ^ inn[:, j])
            nf |= b.astype(np.int64) << t
        per_class = np.bincount(nf[good], minlength=1 << m)
        classes = np.unique(nf)
        check(len(classes) == 1 << (m - n + 1), 'S1 number of labelled classes n=%d' % n)
        check(np.all(per_class[classes] == 1), 'S1 labelled uniqueness n=%d' % n)
        print('n=%d: all %d labelled switching classes contain exactly one such member' % (n, len(classes)))


# ------------------------------------------------------------------------------------------------- S2
def perm_of_type(lams):
    g, s = [], 0
    for l in lams:
        g += [s + (i + 1) % l for i in range(l)]
        s += l
    return g


def v2(x):
    return (x & -x).bit_length() - 1


def orbits_pairs_cyclic(g, n):
    seen, orbs = set(), []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in seen:
                continue
            o, a, b = [], i, j
            while (min(a, b), max(a, b)) not in seen:
                p = (min(a, b), max(a, b))
                seen.add(p)
                o.append(p)
                a, b = g[a], g[b]
            orbs.append(o)
    return orbs


def gf2_rank_member(rows, vec):
    basis = {}

    def red(x):
        while x:
            h = x.bit_length() - 1
            if h in basis:
                x ^= basis[h]
            else:
                return x
        return 0
    for r in rows:
        x = red(r)
        if x:
            basis[x.bit_length() - 1] = x
    return len(basis), red(vec) == 0


def class_size(lams, n):
    c = math.factorial(n)
    for l, mult in Counter(lams).items():
        c //= (l ** mult) * math.factorial(mult)
    return c


def s2():
    hdr('S2: Theorem 2 -- level lemma, twisting witness, twisted Euler graph counts')
    twisted = []
    for n in range(1, 14):
        Etot = Ctot = Ttot = 0
        for part in partitions(n):
            lams = sorted([l for l, mult in part.items() for _ in range(mult)], reverse=True)
            g = perm_of_type(lams)
            orbs = orbits_pairs_cyclic(g, n)
            k = len(orbs)
            rows = [0] * n
            for idx, o in enumerate(orbs):
                for (i, j) in o:
                    rows[i] ^= 1 << idx
                    rows[j] ^= 1 << idx
            cvec = 0
            for idx, o in enumerate(orbs):
                if sum(1 for (i, j) in o if g[i] > g[j]) & 1:
                    cvec |= 1 << idx
            r, inrow = gf2_rank_member(rows, cvec)
            dimK = k - r
            check(dimK == k - len(lams) + (1 if any(l & 1 for l in lams) else 0), 'S2 dim W_g %s' % lams)
            level = len({v2(l) for l in lams}) == 1
            check(inrow == level, 'S2 level lemma %s' % lams)       # eps vanishes on W_g iff g level
            if not level and n <= 10:
                # explicit witness F0 = antipodal orbit of a top-valuation cycle + one cross orbit
                starts = [sum(lams[:t]) for t in range(len(lams))]
                cyc = [list(range(starts[t], starts[t] + lams[t])) for t in range(len(lams))]
                top = max(range(len(lams)), key=lambda t: v2(lams[t]))
                low = min(range(len(lams)), key=lambda t: v2(lams[t]))
                x0, y0 = cyc[top][0], cyc[low][0]
                m = lams[top]
                anti = next(o for o in orbs if (min(x0, cyc[top][m // 2]), max(x0, cyc[top][m // 2])) in o)
                cross = next(o for o in orbs if (min(x0, y0), max(x0, y0)) in o)
                F0 = anti + cross
                check(L.is_euler(F0, n), 'S2 witness Euler %s' % lams)
                check(L.eps(F0, g) == -1, 'S2 witness twisted %s' % lams)
            cs = class_size(lams, n)
            Etot += cs * 2 ** dimK
            if level:
                Ctot += cs * 2 ** dimK
            else:
                Ttot += cs * 2 ** dimK
        f = math.factorial(n)
        check(Etot % f == 0 and Ctot % f == 0 and Ttot % f == 0, 'S2 integrality n=%d' % n)
        check(Etot // f == L.A002854[n] and Ctot // f == L.A049313[n], 'S2 Burnside totals n=%d' % n)
        twisted.append(Ttot // f)
    print('level lemma holds for every cycle type, n <= 13 (GF(2)); witness F0 checked for n <= 10')
    print('twisted Euler graphs = A002854 - A049313, n = 1..13:', twisted)
    check(twisted == [L.A002854[n] - L.A049313[n] for n in range(1, 14)], 'S2 twisted counts')


# ------------------------------------------------------------------------------------------------- S3
def s3():
    hdr('S3: Theorem 3 -- graphs, even n: [empty] and [K_n] are S_n-fixed, K_n is not Euler')
    for n in (4, 6, 8):
        Kn = [(i, j) for i in range(n) for j in range(i + 1, n)]
        check(not L.is_euler(Kn, n), 'S3 K_n not Euler')
        # K_n is not complete bipartite, so [K_n] != [empty]
        check(not any(set(Kn) == {(i, j) for (i, j) in Kn if ((i in U) != (j in U))}
                      for r in range(n + 1) for U in map(set, itertools.combinations(range(n), r))),
              'S3 classes distinct')
        print('n=%d: K_n has odd degrees; [K_n] != [empty]; both classes S_n-fixed -> both map to the empty graph' % n)


# ------------------------------------------------------------------------------------------------- S4
def untwisted_euler(n):
    graphs = [E for E in L.all_graphs(n) if L.is_euler(E, n)]
    tw = L.twisted_flags([(n, E) for E in graphs])
    can = L.canon_many([L.g6_of_edges(E, n) for E in graphs])
    return {c for c, t in zip(can, tw) if not t}, len(graphs)


def invariant_canons(H, n, euler_only):
    orbs = L.pair_orbits(n, H)
    g6s = []
    for mask in range(1 << len(orbs)):
        E = [p for i in range(len(orbs)) if mask >> i & 1 for p in orbs[i]]
        if euler_only and not L.is_euler(E, n):
            continue
        g6s.append(L.g6_of_edges(E, n))
    return set(L.canon_many(g6s))


def s4():
    hdr('S4: Theorem 4 -- natural matchings, switching classes -> untwisted Euler graphs')
    expect = {3: 0, 4: 0, 5: 1, 6: 1, 7: 2, 8: 3, 9: 4}
    for n in range(3, 10 if FULL else 9):
        unt, neuler = untwisted_euler(n)
        check(neuler == L.A002854[n] and len(unt) == L.A049313[n], 'S4 counts n=%d' % n)
        if n % 2:
            c = (n - 1) // 2 % 2
            cls = [A for A in L.gentourng(n) if all(s % 2 == c for s in L.scores(A))]
            stabs = [L.aut_tournament(A) for A in cls]          # Stab(C) = Aut(T*) (Corollary 1.1)
        else:
            tours = L.gentourng(n)
            can = L.canon_many([d for A in tours for d in L.descendant_d6s(A)])
            reps = {}
            for i, A in enumerate(tours):
                reps.setdefault(min(can[i * n:(i + 1) * n]), A)     # complete class invariant
            cls = list(reps.values())
            stabs = [L.class_stabilizer(A) for A in cls]
        check(len(cls) == L.A049313[n], 'S4 classes n=%d' % n)
        R = sorted(unt)
        adj = [R if len(H) == 1 else sorted(invariant_canons(H, n, True) & unt) for H in stabs]
        size, match_r = L.max_matching(adj, len(R))
        Ls, Rs = L.hall_violator(adj, match_r)
        check(len(cls) - size == expect[n], 'S4 deficiency n=%d' % n)
        check(len(Ls) - len(Rs) == len(cls) - size, 'S4 Konig n=%d' % n)
        print('n=%d: classes = untwisted Euler graphs = %d, Euler graphs %d; nontrivial stabilizer orders %s'
              % (n, len(cls), neuler, sorted((len(H) for H in stabs if len(H) > 1), reverse=True)))
        print('      maximum natural matching %d of %d (deficiency %d)%s' % (
            size, len(cls), len(cls) - size,
            '; Hall violator %d classes (stabilizer orders %s) vs %d Euler graphs'
            % (len(Ls), sorted((len(stabs[u]) for u in Ls), reverse=True), len(Rs)) if Ls else ''))
        if n == 5:
            names = sorted(len(E) for E in L.all_graphs(5) if L.is_euler(E, 5)
                           and L.canon_many([L.g6_of_edges(E, 5)])[0] in unt)
            check(names == [0, 4], 'S4 n=5 untwisted = empty, C4+K1')
            check(sorted(len(H) for H in stabs) == [3, 5], 'S4 n=5 stabilizers')
            print('      n=5: untwisted Euler graphs have 0 and 4 edges (empty, C4+K1); stabilizers Z3, Z5')


# ------------------------------------------------------------------------------------------------- S5
def even_graphs(n):
    graphs = L.all_graphs(n)
    tw = []
    B = 2000
    for i in range(0, len(graphs), B):
        tw += L.twisted_flags([(n, E) for E in graphs[i:i + B]])
    can = L.canon_many([L.g6_of_edges(E, n) for E in graphs])
    return {c for c, t in zip(can, tw) if not t}, dict(zip(can, graphs))


def s5():
    hdr('S5: Theorems 8, 9 -- Royle setting: tournaments -> even graphs')
    expect = {3: 0, 4: 0, 5: 1, 6: 2, 7: 5, 8: 12}
    for n in range(3, 9 if FULL else 8):
        even, byc = even_graphs(n)
        tours = L.gentourng(n)
        check(len(tours) == L.A000568[n] and len(even) == L.A000568[n], 'S5 RPGFD count n=%d' % n)
        R = sorted(even)
        auts = [L.aut_tournament(T) for T in tours]
        adj = [R if len(H) == 1 else sorted(invariant_canons(H, n, False) & even) for H in auts]
        size, match_r = L.max_matching(adj, len(R))
        Ls, Rs = L.hall_violator(adj, match_r)
        check(len(tours) - size == expect[n], 'S5 deficiency n=%d' % n)
        print('n=%d: tournaments = even graphs = %d; maximum natural matching %d (deficiency %d)%s' % (
            n, len(tours), size, len(tours) - size,
            '; Hall violator %d tournaments (|Aut| %s) vs %d even graphs'
            % (len(Ls), dict(sorted(Counter(len(auts[u]) for u in Ls).items())), len(Rs)) if Ls else ''))
        if n == 3:
            i3 = next(i for i in range(len(tours)) if len(auts[i]) == 3)
            check([len(byc[c]) for c in adj[i3]] == [0], 'S5 n=3: C3 -> empty graph')
            check(sorted(len(byc[c]) for c in R) == [0, 2], 'S5 n=3: even graphs empty, P3')
        if n == 5:
            cons = [i for i in range(len(tours)) if len(auts[i]) > 1]
            check(sorted(len(auts[i]) for i in cons) == [3, 3, 3, 3, 5], 'S5 n=5 constrained tournaments')
            partners = set()
            for i in cons:
                partners |= set(adj[i])
            check(sorted(len(byc[c]) for c in partners) == [0, 3, 4, 6], 'S5 n=5 partners')
            print('      n=5: 5 tournaments with |Aut| in {3,5}; their possible partners: %d even graphs with %s edges'
                  % (len(partners), sorted(len(byc[c]) for c in partners)))


# ------------------------------------------------------------------------------------------------- S6
def s6():
    hdr('S6: Lemma 5.3 -- rigid blocks')
    blocks = {'P3': L.paley(3), 'P7': L.paley(7), 'P11': L.paley(11), 'P19': L.paley(19),
              'R5': L.rquart(5), 'R13': L.rquart(13), 'R29': L.rquart(29),
              'P3[P3]': L.lex(L.paley(3), L.paley(3)), 'P3[R5]': L.lex(L.paley(3), L.rquart(5)),
              'R5[P3]': L.lex(L.rquart(5), L.paley(3)), 'P3[P7]': L.lex(L.paley(3), L.paley(7)),
              'P7[P3]': L.lex(L.paley(7), L.paley(3)), 'R5[R5]': L.lex(L.rquart(5), L.rquart(5))}
    for name, (A, gens) in blocks.items():
        n = len(A)
        check(L.is_tournament(A), 'S6 tournament ' + name)
        check(all(L.is_aut(A, g) for g in gens), 'S6 automorphisms ' + name)
        check(L.is_transitive_group(n, gens), 'S6 transitive ' + name)
        graphs, k = L.invariant_graphs(n, gens)
        tw = L.twisted_flags([(n, E) for E in graphs])
        check(all(tw), 'S6 rigid ' + name)
        print('%-7s order %2d: %2d orbits on pairs, all %d nonempty invariant graphs twisted' % (name, n, k, len(graphs)))


# ------------------------------------------------------------------------------------------------- S7
def s7():
    hdr('S7: Theorems 5, 6 -- forced classes, two non-isomorphic ones for every n in range')
    hi = 100 if FULL else 60
    for n in range(5, hi + 1):
        reps = L.representations(n)[:3]
        keys = {}
        for rep in reps:
            A, gens = L.build(rep)
            check(len(A) == n and L.is_tournament(A), 'S7 build %s' % (rep,))
            check(all(L.is_aut(A, g) for g in gens), 'S7 automorphisms %s' % (rep,))
            graphs, k = L.invariant_graphs(n, gens, euler_only=True)
            check(all(L.twisted_flags([(n, E) for E in graphs])), 'S7 forced %s' % (rep,))
            keys.setdefault(L.class_key(A), []).append(rep)
        check(len(keys) >= 2, 'S7 two classes n=%d' % n)
        if n <= 12 or n % 10 == 0:
            print('n=%3d: forced constructions %s -> %d non-isomorphic classes' % (n, reps, len(keys)))
    print('every n in [5, %d]: at least two non-isomorphic forced switching classes' % hi)


# ------------------------------------------------------------------------------------------------- S8
def arithmetic(N):
    idx = np.arange(N + 1)
    spf = np.zeros(N + 1, dtype=np.int64)
    for p in range(2, int(N ** 0.5) + 1):
        if spf[p] == 0:
            m = spf[p * p::p]
            m[m == 0] = p
            spf[p * p::p] = m
    isp = np.zeros(N + 1, dtype=bool)
    isp[2:] = spf[2:] == 0
    bad = np.zeros(N + 1, dtype=bool)
    for p in np.nonzero(isp & (idx % 8 == 1))[0]:
        bad[p::p] = True
    sig = np.zeros(N + 1, dtype=bool)
    sig[1::2] = True
    sig &= ~bad
    L2 = 1
    while L2 < 2 * (N + 1):
        L2 *= 2
    F = np.fft.rfft(sig.astype(float), L2)
    conv = np.rint(np.fft.irfft(F * F, L2)[:N + 1]).astype(np.int64)
    half = np.zeros(N + 1, dtype=np.int64)
    ev = np.arange(0, N + 1, 2)
    half[ev] = sig[ev // 2]
    unord = (conv + half) // 2                       # #{w1 <= w2 in Sigma, w1 + w2 = n}
    h = np.zeros(N + 1)
    w = np.arange(1, N // 2 + 1)
    h[2 * w] = sig[w]
    conv2 = np.rint(np.fft.irfft(np.fft.rfft(h, L2) * F, L2)[:N + 1]).astype(np.int64)   # 2w + w' = n
    big = np.zeros(N + 1, dtype=np.int64)
    sigI = sig.astype(np.int64)
    for q in np.nonzero(isp & (idx % 4 == 3))[0]:
        q = int(q)
        lo, hi = q + 1, min(2 * q - 1, N)
        if lo <= hi:
            big[lo:hi + 1] += sigI[1:hi - q + 1]          # n - q in Sigma, n/2 < q <= n - 1
    A = np.zeros(N + 1, dtype=np.int64)
    odd = np.arange(1, N + 1, 2)
    A[odd] = conv2[odd] + (isp[odd] & sig[odd])
    A[np.arange(2, N + 1, 4)] = unord[np.arange(2, N + 1, 4)]
    A[np.arange(4, N + 1, 4)] = big[np.arange(4, N + 1, 4)]
    return A, sig


def s8():
    hdr('S8: Theorem 6 -- the arithmetic condition A(n) >= 2, 5 <= n <= 10^6')
    N = 10 ** 6
    A, sig = arithmetic(N)
    # spot check against a direct count
    for n in (5, 6, 7, 8, 12, 16, 40, 41, 42, 44, 99, 100, 1001, 1002, 1004):
        if n % 2:
            d = sum(1 for w in range(1, n // 2 + 1) if L.in_sigma(w) and n - 2 * w >= 1 and L.in_sigma(n - 2 * w))
            d += 1 if (isprime(n) and L.in_sigma(n)) else 0
        elif n % 4 == 2:
            d = sum(1 for w1 in range(1, n // 2 + 1) if L.in_sigma(w1) and L.in_sigma(n - w1))
        else:
            d = sum(1 for q in range(n // 2 + 1, n) if isprime(q) and q % 4 == 3 and L.in_sigma(n - q))
        check(int(A[n]) == d, 'S8 direct count n=%d' % n)
    fails = [n for n in range(5, N + 1) if A[n] < 2]
    check(fails == [8, 16, 40], 'S8 exceptions')
    print('A(n) >= 2 for all 5 <= n <= 10^6 except n in %s (covered by S7 / S4)' % fails)
    for r in range(4):
        sel = np.arange(1000, N + 1)
        sel = sel[sel % 4 == r]
        print('  n = %d mod 4, 1000 <= n <= 10^6: min A(n) = %d' % (r, int(A[sel].min())))
    print('  Sigma up to 10^6: %d elements' % int(sig.sum()))


# ------------------------------------------------------------------------------------------------- S9
def s9():
    hdr('S9: Theorem 6, n = 0 mod 4 -- balanced pairs of B_w1 => P_q')
    for (w1, q) in [(1, 7), (5, 7), (1, 11), (5, 11), (1, 19), (9, 19), (5, 23), (13, 19), (1, 23), (9, 23)]:
        n = w1 + q
        if n % 4 or not (L.in_sigma(w1) and isprime(q) and q % 4 == 3 and 2 * q > n):
            continue
        A, _ = L.build((w1, q))
        s = L.scores(A)
        bal = 0
        for u in range(n):
            for v in range(u + 1, n):
                S = sum(1 for x in range(n) if x not in (u, v) and (A[x][u] + A[x][v]) == 1)
                par = [(s[x] - A[x][u] - A[x][v]) & 1 for x in range(n) if x not in (u, v)]
                k = sum(par)
                check(min(k, n - 2 - k) == min(S, n - 2 - S), 'S9 partition = S(u,v)')
                bal += (S == (n - 2) // 2)
        want = n * (n - 1) // 2 if w1 == 1 else w1 * q
        check(bal == want, 'S9 balanced (%d,%d)' % (w1, q))
        print('B_%d => P_%d (n=%d): %d balanced pairs = %s' % (w1, q, n, bal, 'C(n,2)' if w1 == 1 else 'w1*q'))


# ------------------------------------------------------------------------------------------------ S10
def s10():
    hdr('S10: prime blocks and doubly regular Paley tournaments')
    for p in (3, 7, 11, 19, 23, 31):
        A, _ = L.paley(p)
        check(L.prime_tournament(A), 'S10 P_%d prime' % p)
        vals = {sum(A[x][z] * A[y][z] for z in range(p)) for x in range(p) for y in range(p) if x != y}
        check(vals == {(p - 3) // 4}, 'S10 P_%d doubly regular' % p)
    for p in (5, 13, 29, 37):
        check(L.prime_tournament(L.rquart(p)[0]), 'S10 R_%d prime' % p)
    check(not L.prime_tournament(L.lex(L.paley(3), L.rquart(5))[0]), 'S10 lex product decomposable')
    print('P_p prime and doubly regular (p = 3..31), R_p prime (p = 5, 13, 29, 37), P_3[R_5] decomposable')


# ------------------------------------------------------------------------------------------------ S11
def s11():
    hdr('S11: Theorem 10 -- Royle-forced families')
    B1 = L.rigid(15, order=[3, 5])
    B2 = L.rigid(15, order=[5, 3])
    cases = {'P3[R5]': B1, 'R5[P3]': B2,
             'P3[R5] => P3[R5]': L.compose([B1, B1], L.transitive(2)),
             'R5[P3] => R5[P3]': L.compose([B2, B2], L.transitive(2)),
             'P3[R5] => R5[P3]': L.compose([B1, B2], L.transitive(2))}
    keys = {}
    for name, (A, gens) in cases.items():
        n = len(A)
        check(L.is_tournament(A) and all(L.is_aut(A, g) for g in gens), 'S11 build ' + name)
        graphs, k = L.invariant_graphs(n, gens)
        check(all(L.twisted_flags([(n, E) for E in graphs])), 'S11 Royle-forced ' + name)
        keys[name] = L.canon_many([L.d6(A)])[0]
        print('%-18s n=%d: %d orbits on pairs, every nonempty invariant graph odd -> Royle-forced' % (name, n, k))
    check(keys['P3[R5]'] != keys['R5[P3]'], 'S11 non-isomorphic n=15')
    check(len({keys[x] for x in cases if '=>' in x}) == 3, 'S11 non-isomorphic n=30')
    fam = sorted({m for m in range(3, 100, 2) if L.in_sigma(m) and len(factorint(m)) >= 2}
                 | {2 * m for m in range(3, 50, 2) if L.in_sigma(m) and len(factorint(m)) >= 2})
    check(fam[:7] == [15, 21, 30, 33, 35, 39, 42], 'S11 family list')
    print('no natural bijection (Royle setting) for n in', fam, '...')


# ------------------------------------------------------------------------------------------------ S12
def zl_counts(l, n):
    pairs = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    k = len(pairs)
    N = l ** k
    codes = np.arange(N, dtype=np.int64)
    dig = np.stack([(codes // l ** t) % l for t in range(k)]) if k else np.zeros((0, 1), dtype=np.int64)

    def full(block):
        M = np.zeros((n, n, block.shape[1]), dtype=np.int64)
        for t, (i, j) in enumerate(pairs):
            M[i, j] = block[t]
            M[j, i] = (-block[t]) % l
        return M
    Msw = full(dig)
    Meu = full(dig)
    for i in range(1, n):
        Meu[i, 0] = (-Meu[i].sum(axis=0)) % l
        Meu[0, i] = (-Meu[i, 0]) % l
    assert np.all(Meu.sum(axis=1) % l == 0)

    def encode(M):
        c = np.zeros(M.shape[2], dtype=np.int64)
        for t, (i, j) in enumerate(pairs):
            c += (M[i, j] % l) * l ** t
        return c
    fs = ft = 0
    for g in itertools.permutations(range(n)):
        gi = np.argsort(g)
        P = Msw[np.ix_(gi, gi)]
        s = (-P[0]) % l
        P = (P + s[None, :, :] - s[:, None, :]) % l
        for v in range(n):
            P[v, v] = 0
        fs += int(np.count_nonzero(encode(P) == codes))
        ft += int(np.count_nonzero(encode(Meu[np.ix_(gi, gi)]) == codes))
    f = math.factorial(n)
    assert fs % f == 0 and ft % f == 0
    return fs // f, ft // f


def s12():
    hdr('S12: Theorem 11 -- switching classes vs modular Eulerian matrices in Alt_n(Z/l)')
    plan = [(2, 6), (3, 6), (4, 6), (6, 5), (8, 5)] if FULL else [(2, 5), (3, 5), (4, 5), (6, 4), (8, 4)]
    for l, nmax in plan:
        row = []
        for n in range(2, nmax + 1):
            s, t = zl_counts(l, n)
            check(s == t, 'S12 l=%d n=%d' % (l, n))
            row.append(s)
        if l == 2:
            check(row == [L.A002854[n] for n in range(2, nmax + 1)], 'S12 l=2 is A002854')
        print('l=%d: s_{l,n} = t_{l,n} = %s (n = 2..%d)' % (l, row, nmax))


if __name__ == '__main__':
    for sec in (s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12):
        sec()
    print('\n%d checks; mode %s' % (NCHECK, 'full' if FULL else 'default'))
    print('[%7.1fs] done' % (time.time() - T0), file=sys.stderr)
    print('ALL CHECKS PASSED')
