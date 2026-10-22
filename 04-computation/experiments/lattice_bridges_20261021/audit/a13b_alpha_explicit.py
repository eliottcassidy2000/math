#!/usr/bin/env python3
"""AUDIT pass 2, item 2 (method E, explicit): alpha_k = number of k-sets of pairwise vertex-disjoint directed odd
cycles of a circulant tournament, by (i) explicit DFS enumeration of every directed odd cycle as a vertex sequence
(canonical start = smallest vertex), (ii) grouping by vertex set, (iii) explicit enumeration of all collections of
pairwise disjoint vertex sets (listed in increasing order of their index), weighted by the product of cycle counts.
Independent of alpha_decomp.c and of a13_alpha_ie.py.  Checks H = sum_k 2^k alpha_k against the pass-1 class tables.
Usage: python3 a13b_alpha_explicit.py n [all]   (default: R_n and the energy-minimum class(es))
Output: a13b_alpha_explicit_n.out"""
import sys, os, json
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
n = int(sys.argv[1]); allc = len(sys.argv) > 2
rows = [json.loads(l) for l in open(os.path.join(HERE, f'a3_classes_{n}.txt'))]
U = [a for a in range(1, n) if gcd(a, n) == 1]
canon = lambda S: min(tuple(sorted((a * s) % n for s in S)) for a in U)
m = (n - 1) // 2; R = canon(range(1, m + 1))
Emin = min(r['E'] for r in rows)
targets = rows if allc else [r for r in rows if tuple(r['S']) == R or r['E'] == Emin]
out = open(os.path.join(HERE, f'a13b_alpha_explicit_{n}.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()
for r in targets:
    S = r['S']; Sset = set(S)
    succ = [[(v + s) % n for s in S] for v in range(n)]
    cnt = {}                                   # vertex-set mask -> number of directed cycles on exactly that set
    def dfs(start, v, mask, length):
        for w in succ[v]:
            if w == start:
                if length % 2 == 1 and length >= 3: cnt[mask] = cnt.get(mask, 0) + 1
            elif w > start and not (mask >> w) & 1:
                dfs(start, w, mask | (1 << w), length + 1)
    for s0 in range(n): dfs(s0, s0, 1 << s0, 1)
    items = sorted(cnt.items())
    K = n // 3; alpha = [0] * (K + 1); alpha[0] = 1
    def rec(i0, used, k, prod):
        for i in range(i0, len(items)):
            mk, c = items[i]
            if mk & used: continue
            p = prod * c; alpha[k + 1] += p
            if k + 1 < K: rec(i + 1, used | mk, k + 1, p)
    sys.setrecursionlimit(10000)
    rec(0, 0, 0, 1)
    H = sum(a << k for k, a in enumerate(alpha))
    ck = {}
    for mk, c in items: L = bin(mk).count('1'); ck[L] = ck.get(L, 0) + c
    tag = 'R' if tuple(S) == R else ('Emin' if r['E'] == Emin else '')
    say(f"n={n} S={S} {tag}: alpha = {alpha}; sum 2^k alpha_k = {H}; equals class-table H: {H == r['H']}; "
        f"odd c_k = {[ck.get(L, 0) for L in range(3, n + 1, 2)]}; vertex sets with cycles: {len(items)}")
out.close()
