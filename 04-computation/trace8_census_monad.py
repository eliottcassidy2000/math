#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The tr(A^8) closed-walk census  (monad-explorer-2026-06-14).

Extends the reducibility ladder one rung past THM-500's tr(A^7)=7(c7+TQ).

A closed k-walk in a tournament loop-erases into a multiset of simple directed
cycles whose lengths PARTITION k, each part >= 3 (no 1- or 2-cycles in a
tournament), and which are pairwise OVERLAPPING (a single closed walk is
connected). Partitions of 8 into parts >= 3:  {8}, {4,4}, {3,5}.

  {8}    : simple 8-cycle.                 8 rooted versions  -> 8 * c8
  {4,4}  : (a) ONE 4-cycle traversed twice -> 4 rooted        -> 4 * c4
           (b) TWO distinct 4-cycles glued at a shared vertex -> Q44_overlap
  {3,5}  : a triangle and a 5-cycle glued at a shared vertex  -> TF_overlap

We verify the lower rungs exactly and FIT the unknown gluing coefficients for the
Q44 / TF overlap terms (split by overlap size, as in trace7), solving for exact
integers. This pins the n=8 NON-SPECTRAL carriers: the spectrum fixes tr(A^8); the
simple count c8 then trades against {c4-doubled (spectral), Q44, TF}.
"""

import sys, itertools, random
import numpy as np
from fractions import Fraction
from collections import defaultdict

random.seed(8)


def random_tournament(n):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if random.getrandbits(1):
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A


def simple_cycles_by_len(A, n, lengths):
    """All simple directed cycles of the given lengths, each as a tuple (canonical
    rotation: min vertex first). Returns dict length -> list of vertex-tuples."""
    adj = [[j for j in range(n) if A[i, j]] for i in range(n)]
    Ls = set(lengths)
    maxL = max(Ls)
    out = defaultdict(list)
    for start in range(n):
        stack = [(start, 1 << start, (start,))]
        while stack:
            v, vis, path = stack.pop()
            if len(path) > maxL:
                continue
            for w in adj[v]:
                if w == start and len(path) in Ls:
                    out[len(path)].append(path)
                elif w > start and not (vis >> w) & 1 and len(path) < maxL:
                    stack.append((w, vis | (1 << w), path + (w,)))
    return out


def overlap_pairs(listA, listB, same=False):
    """count pairs (a in A, b in B) sharing >=1 vertex, split by |intersection|.
    If same=True, count unordered distinct pairs within one list."""
    f = defaultdict(int)
    setsA = [frozenset(x) for x in listA]
    if same:
        for i in range(len(setsA)):
            for j in range(i + 1, len(setsA)):
                inter = len(setsA[i] & setsA[j])
                if inter:
                    f[inter] += 1
    else:
        setsB = [frozenset(x) for x in listB]
        for sa in setsA:
            for sb in setsB:
                inter = len(sa & sb)
                if inter:
                    f[inter] += 1
    return f


def features8(A, n):
    cyc = simple_cycles_by_len(A, n, [3, 4, 5, 8])
    tri = cyc[3]; quad = cyc[4]; five = cyc[5]
    c3 = len(tri); c4 = len(quad); c5 = len(five); c8 = len(cyc[8])
    # Q44: distinct 4-cycle pairs sharing vertices, by overlap size 1..3
    q44 = overlap_pairs(quad, None, same=True)
    # TF: (triangle,5-cycle) pairs sharing vertices, by overlap size 1..3
    tf = overlap_pairs(tri, five, same=False)
    feat = dict(c3=c3, c4=c4, c5=c5, c8=c8)
    for s in (1, 2, 3):
        feat[f'q44_{s}'] = q44.get(s, 0)
        feat[f'tf_{s}'] = tf.get(s, 0)
    return feat


def trace(A, k):
    P = np.linalg.matrix_power(A.astype(object), k)
    return int(np.trace(P))


def solve_rational(Am, bv):
    n = len(Am)
    M = [[Am[i][j] for j in range(n)] + [bv[i]] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] != 0), None)
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                fac = M[r][col]
                M[r] = [M[r][j] - fac * M[col][j] for j in range(n + 1)]
    return [M[i][n] for i in range(n)]


def main():
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    NS = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    data = []
    for _ in range(NS):
        A = random_tournament(n)
        feat = features8(A, n)
        feat['tr8'] = trace(A, 8)
        data.append(feat)

    print(f"=== tr(A^8) census fit  (n={n}, {NS} random tournaments) ===", flush=True)
    # R8 = tr8 - 8*c8 - 4*c4   should be explained by overlap terms only
    cols = ['q44_1', 'q44_2', 'q44_3', 'tf_1', 'tf_2', 'tf_3']
    M = [[Fraction(int(d[c])) for c in cols] for d in data]
    y = [Fraction(int(d['tr8'] - 8 * d['c8'] - 4 * d['c4'])) for d in data]
    # exact least squares via normal equations
    K = len(cols)
    Mt = [[sum(M[r][a] * M[r][b] for r in range(NS)) for b in range(K)] for a in range(K)]
    Mty = [sum(M[r][a] * y[r] for r in range(NS)) for a in range(K)]
    sol = solve_rational(Mt, Mty)
    if sol is None:
        print("  normal equations singular (collinear features); trying float lstsq", flush=True)
        Mf = np.array([[float(d[c]) for c in cols] for d in data])
        yf = np.array([float(d['tr8'] - 8 * d['c8'] - 4 * d['c4']) for d in data])
        coef, *_ = np.linalg.lstsq(Mf, yf, rcond=None)
        print("  float coeffs:", dict(zip(cols, [round(c, 3) for c in coef])), flush=True)
    else:
        print("  exact LS coeffs for R8 = tr8 - 8*c8 - 4*c4:", flush=True)
        print("   ", dict(zip(cols, sol)), flush=True)
        ok = all(sum(M[r][a] * sol[a] for a in range(K)) == y[r] for r in range(NS))
        print(f"    exact fit on all {NS}: {ok}", flush=True)

    # test the clean uniform hypothesis: tr8 = 8*c8 + 4*c4 + 8*Q44 + 8*TF
    bad = 0; worst = None
    for d in data:
        Q44 = d['q44_1'] + d['q44_2'] + d['q44_3']
        TF = d['tf_1'] + d['tf_2'] + d['tf_3']
        lhs = d['tr8']
        rhs = 8 * d['c8'] + 4 * d['c4'] + 8 * Q44 + 8 * TF
        if lhs != rhs:
            bad += 1
            if worst is None:
                worst = (lhs, rhs, d)
    print(flush=True)
    print(f"=== CLEAN HYPOTHESIS:  tr(A^8) = 8*c8 + 4*c4 + 8*Q44 + 8*TF ===", flush=True)
    print(f"    Q44 = #distinct 4-cycle pairs sharing >=1 vertex", flush=True)
    print(f"    TF  = #(triangle, 5-cycle) pairs sharing >=1 vertex", flush=True)
    print(f"    exact on {NS-bad}/{NS} tournaments (bad={bad})", flush=True)
    if worst:
        lhs, rhs, d = worst
        print(f"    first mismatch lhs={lhs} rhs={rhs}: {d}", flush=True)

    # sample rows
    print("\n  sample (tr8, c8, c4, Q44=q44_*, TF=tf_*):", flush=True)
    for d in data[:6]:
        Q44 = d['q44_1'] + d['q44_2'] + d['q44_3']
        TF = d['tf_1'] + d['tf_2'] + d['tf_3']
        print(f"    tr8={d['tr8']} c8={d['c8']} c4={d['c4']} Q44={Q44} TF={TF}"
              f"  (q44={d['q44_1']},{d['q44_2']},{d['q44_3']} tf={d['tf_1']},{d['tf_2']},{d['tf_3']})", flush=True)


if __name__ == "__main__":
    main()
