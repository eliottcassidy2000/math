#!/usr/bin/env python3
"""AUDIT pass 2, item 1: HYP-474 ("Paley maximizes total odd cycles at EVERY length", circulants, p = 3 mod 4).
Reads a12_all{19,23}_hk.out (Held-Karp from 0, own code) and compares, for every odd k, P_q against every circulant
class; also checks the H column against the pass-1 class tables (a3_classes_n.txt), c_3, c_4, c_5 against traces,
and the session's C= column of circ_19.out / circ_23.out.  Output: a12b_analysis.out"""
import os, re, json, glob
import numpy as np
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a12b_analysis.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
def parse(fn):
    rows = {}
    for line in open(fn):
        left, right = line.split(' : ')
        S = tuple(int(x) for x in left.split()[1:])
        d = dict((k, int(v)) for k, v in re.findall(r'(H|c\d+)=(\d+)', right))
        rows[S] = d
    return rows
def theirC(n):
    d = {}
    for f in glob.glob(os.path.join(HERE, '..', '..', 'tournament_pseudoprimes_20261021', 'circ_*.out')):
        for line in open(f):
            m = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=\d+ H=\d+ H/n=\d+ C=(\d+)', line)
            if m and int(m.group(1)) == n: d[tuple(sorted(int(x) for x in m.group(2).split(',')))] = int(m.group(3))
    return d
for n in (19, 23):
    rows = parse(os.path.join(HERE, f'a12_all{n}_hk.out'))
    U = [a for a in range(1, n) if gcd(a, n) == 1]
    canon = lambda S: min(tuple(sorted((a * s) % n for s in S)) for a in U)
    m = (n - 1) // 2
    R = canon(range(1, m + 1)); P = canon(sorted({x * x % n for x in range(1, n)}))
    Hp1 = {tuple(json.loads(l)['S']): json.loads(l)['H'] for l in open(os.path.join(HERE, f'a3_classes_{n}.txt'))}
    say(f"n={n}: {len(rows)} classes; H (Held-Karp from 0) equals pass-1 class tables: {all(rows[S]['H'] == Hp1[S] for S in rows)}")
    th = theirC(n); thc = {canon(S): c for S, c in th.items()}
    say(f"   session C= column agrees with own c_{n} on {sum(1 for S in rows if thc.get(S) == rows[S]['c%d' % n])} of {len(rows)} classes")
    # traces for c3, c4, c5
    okt = True
    for S in rows:
        A = np.zeros((n, n), dtype=object)
        for i in range(n):
            for s in S: A[i, (i + s) % n] = 1
        A = A.astype(np.int64); A2 = A @ A; A4 = A2 @ A2
        t3, t4, t5 = np.trace(A2 @ A), np.trace(A4), np.trace(A4 @ A)
        if not (t3 == 3 * rows[S]['c3'] and t4 == 4 * rows[S]['c4'] and t5 == 5 * rows[S]['c5']): okt = False
    say(f"   c3, c4, c5 equal tr A^k / k (THM-118) on all classes: {okt}")
    for k in range(3, n + 1, 2):
        key = 'c%d' % k; best = max(rows[S][key] for S in rows); arg = [S for S in rows if rows[S][key] == best]
        tag = lambda S: 'P' if S == P else ('R' if S == R else str(list(S)))
        where = 'all classes (constant)' if len(arg) == len(rows) else [tag(S) for S in arg]
        say(f"   k={k:2d}: P_{n} = {rows[P][key]}, R_{n} = {rows[R][key]}, max over classes = {best} at {where}; P maximal: {rows[P][key] == best}")
out.close()
