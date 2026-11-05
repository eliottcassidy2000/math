#!/usr/bin/env python3
"""two_apex_drift_functionals_20261105.py -- opus-2026-11-05-S80

Test of HYP-9475 (B): for mixed blocks d of the two-apex grammar G(2,0), Session 194 proved
    H(u d^k) = A(u) + J(u) k + K k^2 + (D(u) + E(u) k) 2^(Lk) + G(u) 3^(Lk)   for k >= 3,
with L = |d|.  Here the six coefficients are solved exactly (fractions) from k = 3..8 for every prefix u of length
<= PMAX (the two-apex automaton dump of THM-4878, 200 cut states, letters 0..3 = (x1 beats v, x2 beats v) bits), the
form is checked on k = 9, 10, and the set of values of the degree-1 drift J(u) over strong prefixes (E of the final
apex placement... strongness is taken as: the tournament of u d^3 is strong, i.e. H(u d^3) > H of its ordered-join
factorisation is not checked; instead we record J over all prefixes and over prefixes whose word is "reduced") is
reported: parity pattern, integers attained in a window, and the least prefix for small targets.
Usage: python3 two_apex_drift_functionals_20261105.py dump_k2_w0.txt PMAX block1 block2 ...   (blocks as digit strings)
"""
import sys
from fractions import Fraction
from collections import defaultdict
import numpy as np

def parse(path):
    states = {}; trans = defaultdict(list); apex = defaultdict(list)
    for line in open(path):
        parts = line.split()
        if not parts: continue
        if parts[0] == 'S': states[int(parts[1])] = (int(parts[2][2:]), int(parts[3].split('=')[1]))
        elif parts[0] == 'T':
            sid = int(parts[1]); li = int(parts[2]); lst = [(int(x.split('*')[0]), int(x.split('*')[1])) for x in parts[7:]]
            trans[(sid, li)] = lst
        elif parts[0] == 'A':
            sid = int(parts[1]); a = int(parts[2]); lst = [(int(x.split('*')[0]), int(x.split('*')[1])) for x in parts[4:]]
            apex[(sid, a)] = lst
    return states, trans, apex

def main():
    path, PMAX = sys.argv[1], int(sys.argv[2]); blocks = sys.argv[3:]
    states, trans, apex = parse(path); n = max(states) + 1
    M = [np.zeros((n, n), dtype=object) for _ in range(4)]
    for (s, li), lst in trans.items():
        for (t, m) in lst: M[li][s, t] += m
    P = [np.zeros((n, n), dtype=object) for _ in range(2)]
    for (s, a), lst in apex.items():
        for (t, m) in lst: P[a][s, t] += m
    f = np.zeros(n, dtype=object)
    for s, (c, comp) in states.items():
        if comp: f[s] = 1
    fin = P[0].dot(P[1].dot(f))            # apex placement then complete indicator (column functional)
    e0 = np.zeros(n, dtype=object); e0[0] = 1
    for blk in blocks:
        L = len(blk); Md = np.identity(n, dtype=object)
        for ch in blk: Md = Md.dot(M[int(ch)])
        # powers of Md applied to fin (columns): fin_k = Md^k fin
        fins = [fin]
        for k in range(1, 11): fins.append(Md.dot(fins[-1]))
        # basis functions
        def basis(k): return [Fraction(1), Fraction(k), Fraction(k * k), Fraction(2 ** (L * k)), Fraction(k * 2 ** (L * k)), Fraction(3 ** (L * k))]
        ks = [3, 4, 5, 6, 7, 8]
        Bm = [basis(k) for k in ks]
        # solve coefficient vector c from values h_k = sum c_i basis_i(k): invert the 6x6 exactly once
        def solve(h):
            # Gaussian elimination over Fractions
            A = [row[:] + [Fraction(hv)] for row, hv in zip(Bm, h)]
            m = 6
            for col in range(m):
                piv = next(r for r in range(col, m) if A[r][col] != 0)
                A[col], A[piv] = A[piv], A[col]
                pv = A[col][col]; A[col] = [x / pv for x in A[col]]
                for r in range(m):
                    if r != col and A[r][col] != 0:
                        fac = A[r][col]; A[r] = [x - fac * y for x, y in zip(A[r], A[col])]
            return [A[r][m] for r in range(m)]
        Jvals = defaultdict(list); Kset = set(); bad = 0; count = 0
        # enumerate prefixes u of length 0..PMAX via row vectors
        frontier = [(e0, '')]
        for length in range(0, PMAX + 1):
            nxt = []
            for (v, w) in frontier:
                h = [int(v.dot(fins[k])) for k in range(0, 11)]
                c = solve(h[3:9])
                # check k = 9, 10
                for k in (9, 10):
                    pred = sum(ci * bi for ci, bi in zip(c, basis(k)))
                    if pred != h[k]: bad += 1
                count += 1
                Kset.add(c[2]); Jvals[c[1]].append(w)
                if length < PMAX:
                    for li in range(4): nxt.append((v.dot(M[li]), w + str(li)))
            frontier = nxt
        from math import lcm
        den = 1
        for j in Jvals: den = lcm(den, j.denominator)
        nums = sorted(int(j * den) for j in Jvals)
        kden = 1
        for kk in Kset: kden = lcm(kden, kk.denominator)
        residues = sorted(set(x % 4 for x in nums))
        window = [x for x in range(-60, 61) if Fraction(x, den) in Jvals]
        print(f"block {blk} (L={L}): prefixes up to length {PMAX}: {count}; form failures at k=9,10: {bad}; K values (times {kden}) {sorted(int(kk * kden) for kk in Kset)[:12]}{' ...' if len(Kset) > 12 else ''}")
        print(f"   J = numerator/{den}: {len(nums)} distinct values, range [{nums[0]}, {nums[-1]}], numerators mod 4 {residues}; numerators attained in [-60, 60]: {window}")
        print(f"   missing numerators in [-60, 60] with the attained residues: {[x for x in range(-60, 61) if x % 4 in residues and x not in window]}")
        for x in [x for x in (-8, -4, -2, 0, 2, 4, 8, 1, 3, -1, -3) if Fraction(x, den) in Jvals]:
            print(f"   J = {x}/{den}: shortest prefix {min(Jvals[Fraction(x, den)], key=len)!r} ({len(Jvals[Fraction(x, den)])} prefixes)")

main()
