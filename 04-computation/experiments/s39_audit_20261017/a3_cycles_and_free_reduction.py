#!/usr/bin/env python3
"""S39 audit, item 3: statements 2, 4, 5 of THM-4663 by independent graph computations.

Graph: G_p = Cay(F_(p^2), mu) with F_(p^2) = F_p[th]/(th^2 - th - r) (not the session's model).
  * W_2, W_3, W_4 and the codegree lemma 1 - eta(N(d)(N(d) - 4)) for every d != 0;
  * 5-cycles through 0 by a meet-in-the-middle count over ordered edges (b, c) with two 2-paths 0-a-b, 0-d-c;
  * 6-cycles through 0 as ordered pairs of internally disjoint simple 3-paths 0 -> v (each cycle twice);
  * networkx simple_cycles (length_bound) as a third method for p = 5, 7;
  * closed 5- and 6-walks from 0 enumerated (3-walk + reversed 3-walk), freely reduced, and classified by shape, to test
    every component count in the proofs of statements 4 and 5 (tree-like, 4-cycle + spur, lollipops, and for p = 2 mod 3
    the families A (two triangles at v), B (a triangle at a vertex u != v of a triangle through v), and their overlaps).
Closed forms use a_p, b_p from a2_forms_to_3000.json (built independently in a2).
Usage: python -B a3_cycles_and_free_reduction.py
"""
import json, itertools
from collections import defaultdict, Counter
import numpy as np
from numba import njit
import networkx as nx
from sympy import primerange, jacobi_symbol

F = json.load(open('a2_forms_to_3000.json'))
A15, B46 = F['f15'], F['b46']


def kr(d, n):
    return int(jacobi_symbol(d % n, n)) if d % n else 0


def graph(p):
    r = next(r for r in range(1, p) if pow((1 + 4 * r) % p, (p - 1) // 2, p) == p - 1)
    N = lambda a, b: (a * a + a * b - r * b * b) % p
    mu = [(a, b) for a in range(p) for b in range(p) if N(a, b) == 1]
    assert len(mu) == p + 1
    nbr = np.zeros((p * p, p + 1), dtype=np.int64)
    for a in range(p):
        for b in range(p):
            for j, (da, db) in enumerate(mu):
                nbr[a * p + b, j] = ((a + da) % p) * p + (b + db) % p
    norm = np.array([N(v // p, v % p) for v in range(p * p)], dtype=np.int64)
    return nbr, norm


@njit
def c5_through0(nbr):
    n, d = nbr.shape
    # L[x] = list of a in N(0) cap N(x)
    cnt2 = np.zeros(n, dtype=np.int64)
    L = -np.ones((n, d), dtype=np.int64)
    for i in range(d):
        a = nbr[0, i]
        for j in range(d):
            x = nbr[a, j]
            L[x, cnt2[x]] = a
            cnt2[x] += 1
    tot = 0
    for b in range(1, n):
        for j in range(d):
            c = nbr[b, j]
            if c == 0:
                continue
            for s in range(cnt2[b]):
                a = L[b, s]
                if a == c:
                    continue
                for t in range(cnt2[c]):
                    dd = L[c, t]
                    if dd == a or dd == b:
                        continue
                    tot += 1
    return tot // 2


@njit
def c6_through0(nbr):
    n, d = nbr.shape
    # simple 3-paths 0 - a - b - v grouped by v
    cnt = np.zeros(n, dtype=np.int64)
    for i in range(d):
        a = nbr[0, i]
        for j in range(d):
            b = nbr[a, j]
            if b == 0:
                continue
            for k in range(d):
                v = nbr[b, k]
                if v == 0 or v == a:
                    continue
                cnt[v] += 1
    mx = cnt.max()
    PA = -np.ones((n, mx), dtype=np.int64); PB = -np.ones((n, mx), dtype=np.int64)
    fill = np.zeros(n, dtype=np.int64)
    for i in range(d):
        a = nbr[0, i]
        for j in range(d):
            b = nbr[a, j]
            if b == 0:
                continue
            for k in range(d):
                v = nbr[b, k]
                if v == 0 or v == a:
                    continue
                PA[v, fill[v]] = a; PB[v, fill[v]] = b; fill[v] += 1
    pairs = 0
    for v in range(1, n):
        for s in range(fill[v]):
            a1 = PA[v, s]; b1 = PB[v, s]
            for t in range(fill[v]):
                a2 = PA[v, t]; b2 = PB[v, t]
                if a2 == a1 or a2 == b1 or b2 == a1 or b2 == b1:
                    continue
                pairs += 1
    return pairs // 2


def walks_from0(nbr, k):
    """all k-walks from 0 as tuples of vertices (0, w1, ..., wk)."""
    out = [(0,)]
    for _ in range(k):
        out = [w + (int(u),) for w in out for u in nbr[w[-1]]]
    return out


def free_reduce(w):
    st = [w[0]]
    for u in w[1:]:
        if len(st) >= 2 and st[-2] == u:
            st.pop()
        else:
            st.append(u)
    return tuple(st)


def closed_walks(nbr, k1, k2):
    """closed (k1 + k2)-walks from 0 as k1-walk to v followed by a reversed k2-walk from 0 to v."""
    W1 = walks_from0(nbr, k1); W2 = walks_from0(nbr, k2)
    by = defaultdict(list)
    for w in W2:
        by[w[-1]].append(w)
    out = []
    for w in W1:
        for x in by[w[-1]]:
            out.append(w + tuple(reversed(x[:-1])))
    return out


rep = []
print('Statement 2 and the cycle counts (model th^2 = th + r):')
print(' p  (p/3)  W2 W3 W4 ok  codeg-lemma ok   C5 (direct) = closed form   C6 (direct) = closed form   D(p) = derived')
for p in primerange(5, 60):
    nbr, norm = graph(p)
    d = p + 1
    e3 = kr(p, 3)
    # walk counts from 0
    R = np.zeros(p * p, dtype=np.int64); R[0] = 1
    Ws = []
    for k in range(1, 7):
        Rn = np.zeros_like(R)
        for j in range(d):
            np.add.at(Rn, nbr[:, j], R)
        R = Rn
        Ws.append(int(R[0]))
        if k == 2:
            R2 = R.copy()
    okW = Ws[1] == p + 1 and Ws[2] == (p + 1) * (1 - e3) and Ws[3] == 3 * p * (p + 1)
    eta = lambda x: kr(x, p)
    okC = all(R2[v] == 1 - eta(int(norm[v]) * (int(norm[v]) - 4)) for v in range(1, p * p))
    c5 = c5_through0(nbr)
    C5 = p * p * c5 // 5
    C5f = p * p * (p + 1) * (p * p + 1 - (1 - e3) * p - A15[p]) // 10
    c6 = c6_through0(nbr)
    C6 = p * p * c6 // 6
    C6f = p * p * (p ** 4 + p * p + 2 * p - (p + 1) * B46[p] - 2 * (1 - e3) * (p + 1) * (3 * p - 4)) // 12
    D = Ws[5] - 2 * c6
    Df = (p + 1) * (11 * p * p - 2 * p + 1) if e3 == 1 else (p + 1) * (11 * p * p + 10 * p - 15)
    print('%3d  %2d    %s       %s          %12d = %-12d   %14d = %-14d  %9d = %d' % (p, e3, okW, okC, C5, C5f, C6, C6f, D, Df))
    assert okW and okC and C5 == C5f and C6 == C6f and D == Df
    rep.append(p)

# networkx cross-check
for p in (5, 7):
    nbr, _ = graph(p)
    G = nx.Graph()
    for v in range(p * p):
        for u in nbr[v]:
            G.add_edge(v, int(u))
    cyc = Counter(len(c) for c in nx.simple_cycles(G, length_bound=6))
    print('networkx simple_cycles, p = %d: C_5 = %d, C_6 = %d' % (p, cyc[5], cyc[6]))

# free-reduction classification
print('\nFree-reduction classification of closed walks from a vertex (q = p, d = p + 1):')
for p in (5, 7, 11, 13, 17, 23):
    nbr, _ = graph(p)
    q, d = p, p + 1
    e3 = kr(p, 3)
    nb = [set(int(u) for u in nbr[v]) for v in range(p * p)]
    # 5-walks
    cl5 = closed_walks(nbr, 2, 3)
    cat5 = Counter()
    for w in cl5:
        r = free_reduce(w)
        L = len(r) - 1
        if L == 5:
            cat5['cycle' if len(set(r[:-1])) == 5 else ('lollipop' if r[1] == r[4] else 'other5')] += 1
        else:
            cat5['len%d' % L] += 1
    pred5 = {'len3': 2 * (q + 1) * (4 * d - 3) if e3 == -1 else 0, 'lollipop': (q + 1) * (2 * q - 2) if e3 == -1 else 0}
    print('p = %2d, closed 5-walks: %s; predicted triangle+spur %d, lollipops %d' % (p, dict(cat5), pred5['len3'], pred5['lollipop']))
    assert cat5.get('len3', 0) == pred5['len3'] and cat5.get('lollipop', 0) == pred5['lollipop'] and 'other5' not in cat5
    # 6-walks
    cl6 = closed_walks(nbr, 3, 3)
    cat = Counter()
    Aset, Bset, Lset = set(), set(), set()
    Binc = Counter()            # walk -> number of (T1, split point) representations
    for w in cl6:
        r = free_reduce(w)
        L = len(r) - 1
        if L == 0:
            cat['tree'] += 1; continue
        if L == 4:
            cat['4cyc+spur' if len(set(r[:-1])) == 4 else 'other4'] += 1; continue
        if L == 2:
            cat['len2'] += 1; continue
        assert L == 6
        if len(set(r[:-1])) == 6:
            cat['6cycle'] += 1; continue
        t = r
        inA = t[3] == 0
        rep_a = t[1] == t[4] and t[5] in nb[0]         # triangle (0, t1, t5) with a triangle at t1
        rep_b = t[2] == t[5] and t[1] in nb[0] and t[1] in nb[t[2]]   # triangle (0, t1, t2) with a triangle at t2
        inL = t[1] == t[5] and len(set([t[1], t[2], t[3], t[4]])) == 4 and t[2] != 0 and t[4] != 0
        if inA: Aset.add(t)
        if rep_a or rep_b: Bset.add(t); Binc[t] = int(rep_a) + int(rep_b)
        if inL: Lset.add(t)
        if not (inA or rep_a or rep_b or inL):
            cat['UNCLASSIFIED'] += 1
        cat['nonsimple6'] += 1
    dbl = {t for t in Aset if t[1] == t[4] and t[2] == t[5]}          # one directed triangle traversed twice
    print('p = %2d, closed 6-walks: %s' % (p, dict(cat)))
    print('        tree %d (pred d(5d^2-6d+2) = %d); 4-cycle+spur %d (pred (q^2-1)(5d-4) = %d); lollipops %d (pred (q+1)(q-1)^2 = %d)'
          % (cat['tree'], d * (5 * d * d - 6 * d + 2), cat['4cyc+spur'], (q * q - 1) * (5 * d - 4), len(Lset), (q + 1) * (q - 1) ** 2))
    assert cat['tree'] == d * (5 * d * d - 6 * d + 2) and cat['4cyc+spur'] == (q * q - 1) * (5 * d - 4)
    assert len(Lset) == (q + 1) * (q - 1) ** 2 and cat['UNCLASSIFIED'] == 0 and 'other4' not in cat
    if e3 == -1:
        incB_no_dbl = sum(v for t, v in Binc.items() if t not in dbl)
        both = sum(1 for t, v in Binc.items() if v == 2 and t not in dbl)
        BA = len((Bset & Aset) - dbl)
        LA = len(Lset & Aset)
        LB = len(Lset & Bset)
        per_split_adm = None
        print('        |A| = %d (pred 4q(q+1) = %d); doubled triangles %d (= 2(q+1) = %d)' % (len(Aset), 4 * q * (q + 1), len(dbl), 2 * (q + 1)))
        print('        B incidences incl. doubled triangle %d = 2(q+1)*2*(2q-1) = %d; excluding it %d = 8(q+1)(q-1) = %d'
              % (sum(Binc.values()), 2 * (q + 1) * 2 * (2 * q - 1), incB_no_dbl, 8 * (q + 1) * (q - 1)))
        print('        double counts: B walks at both split points (dbl excluded) %d = 2(q+1); B&A (dbl excluded) %d = 4(q+1); '
              'L&A %d = 2(q+1); L&B %d' % (both, BA, LA, LB))
        nonsimple = len(Aset | Bset | Lset)
        print('        |A u B u L| = %d; session sum 4q(q+1) + 8(q+1)(q-1) + (q+1)(q-1)^2 - 8(q+1) = %d'
              % (nonsimple, 4 * q * (q + 1) + 8 * (q + 1) * (q - 1) + (q + 1) * (q - 1) ** 2 - 8 * (q + 1)))
        assert len(Aset) == 4 * q * (q + 1) and len(dbl) == 2 * (q + 1)
        assert sum(Binc.values()) == 4 * (q + 1) * (2 * q - 1) and incB_no_dbl == 8 * (q + 1) * (q - 1)
        assert both == 2 * (q + 1) and BA == 4 * (q + 1) and LA == 2 * (q + 1) and LB == 0
        assert nonsimple == 4 * q * (q + 1) + 8 * (q + 1) * (q - 1) + (q + 1) * (q - 1) ** 2 - 8 * (q + 1)
    else:
        assert not Aset and not Bset
print('ALL CHECKS PASSED')
