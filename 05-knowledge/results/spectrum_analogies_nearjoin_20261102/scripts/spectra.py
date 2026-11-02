#!/usr/bin/env python3
"""spectra.py (klein-2026-11-02) -- count spectra of tournaments from nauty's gentourng and scripts/hampath.c:
(1) directed Hamiltonian-cycle counts C(T), orders 3..9: attained values and gaps below the maximum;
(2) Hamiltonian paths from a fixed start vertex, orders 3..8: gaps; (3) 3-cycle counts, orders 3..9: gaps;
(4) strong realiser counts per Hamiltonian-path count H through order 9: the strong singletons, and the check that every
attained odd H in [15, 2881] has at least two strong realisers; (5) bipartite tournaments on <= 7 vertices: H values.
usage: python3 spectra.py [gentourng=/opt/homebrew/bin/gentourng]"""
import subprocess, sys, os, itertools
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
GEN = sys.argv[1] if len(sys.argv) > 1 else '/opt/homebrew/bin/gentourng'
HP = os.path.join(HERE, 'hampath')
if not os.path.exists(HP): subprocess.run(['cc', '-O2', '-o', HP, os.path.join(HERE, 'hampath.c')], check=True)
def gen(n): return subprocess.run([GEN, '-q', str(n)], capture_output=True, text=True).stdout
rows = {}
for n in range(3, 10):
    out = subprocess.run([HP, 'census2', str(n)], input=gen(n), capture_output=True, text=True).stdout
    open(os.path.join(DATA, 'census2_%d.txt' % n), 'w').write(out)
    rows[n] = [l.split() for l in out.split('\n') if l.strip()]
print('(1) directed Hamiltonian-cycle counts C(T)')
for n in range(3, 10):
    C = {}
    for h, c, st, code in rows[n]: C[int(c)] = C.get(int(c), 0) + 1
    mx = max(C); miss = [v for v in range(1, mx + 1) if v not in C]
    print('    order %d: max C = %d, distinct %d, values in [1, max] not attained: %s%s' % (n, mx, len(C), miss[:12], ' ...' if len(miss) > 12 else ''))
print('(2) Hamiltonian paths from a fixed start vertex')
allv = set()
for n in range(3, 9):
    out = subprocess.run([HP, 'starts', str(n)], input=gen(n), capture_output=True, text=True).stdout
    for l in out.split('\n'):
        if l.strip(): allv.update(int(x) for x in l.split()[1:1 + n])
    mx = max(allv); print('    orders <= %d: max %d, values in [1, max] not attained: %s' % (n, mx, [v for v in range(1, mx + 1) if v not in allv][:10]))
print('(3) 3-cycle counts')
for n in range(3, 10):
    vals = set()
    for h, c, st, code in rows[n]:
        s = [0] * n; idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if code[idx] == '1': s[i] += 1
                else: s[j] += 1
                idx += 1
        vals.add(comb(n, 3) - sum(comb(x, 2) for x in s))
    mx = max(vals); print('    order %d: max %d, values in [0, max] not attained: %s' % (n, mx, [v for v in range(mx + 1) if v not in vals]))
print('(4) strong realisers per H through order 9')
strong = {}; allc = {}
for n in range(3, 10):
    for h, c, st, code in rows[n]:
        h = int(h); allc[h] = allc.get(h, 0) + 1
        if st == '1': strong[h] = strong.get(h, 0) + 1
singles = [h for h in sorted(allc) if strong.get(h, 0) == 1]
print('    strong singletons:', singles)
print('    strong realisers of 11, 13, 47, 67:', {h: strong.get(h, 0) for h in (11, 13, 47, 67)})
print('    attained odd h in [15, 2881] with fewer than two strong realisers:', [h for h in sorted(allc) if 15 <= h <= 2881 and strong.get(h, 0) < 2])
print('(5) bipartite tournaments (orientations of K_{a,b}, a + b <= 7)')
def hcount(n, adj):
    f = [[0] * n for _ in range(1 << n)]
    for v in range(n): f[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = f[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in range(n):
                if adj[v][w] and not (S >> w & 1): f[S | 1 << w][w] += c
    return sum(f[(1 << n) - 1])
union = set()
for a in range(1, 5):
    for b in range(a, 8 - a):
        n = a + b; pairs = [(i, a + j) for i in range(a) for j in range(b)]; vals = set()
        for bits in range(1 << len(pairs)):
            adj = [[0] * n for _ in range(n)]
            for t, (i, j) in enumerate(pairs):
                if bits >> t & 1: adj[i][j] = 1
                else: adj[j][i] = 1
            vals.add(hcount(n, adj))
        union |= vals; print('    K_{%d,%d}: %s' % (a, b, sorted(vals)))
mx = max(union); print('    union: values in [1, %d] not attained: %s' % (mx, [v for v in range(1, mx + 1) if v not in union]))
