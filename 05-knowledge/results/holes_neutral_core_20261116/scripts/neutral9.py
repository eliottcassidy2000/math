#!/usr/bin/env python3
"""neutral9.py (N-7, N-8): trivial versus nontrivial balanced flips at order 9 (balanced9.txt, flips9.d6 from neutral9.c;
canonical forms by labelg), per fibre for h < 125 (complete below order 10 by Busch's floor 125), joined with the
orders <= 8 census of neutral.py; the neutral-richest strong classes of order 9."""
import subprocess
from collections import defaultdict, Counter
cn = subprocess.run('gentourng -q -z 9 | labelg -q', shell=True, capture_output=True, text=True).stdout.split(); index = {s: i for i, s in enumerate(cn)}
cf = subprocess.run('labelg -q flips9.d6', shell=True, capture_output=True, text=True).stdout.split()
rows = [tuple(map(int, l.split())) for l in open('balanced9.txt')]; assert len(rows) == len(cf)
h9 = [int(l.split()[1]) for l in open('counts.txt') if l.startswith('9 ')]
bal = Counter(); nontriv = defaultdict(set); edges = defaultdict(set)
for (i, a, b, h), s in zip(rows, cf):
    j = index[s]; assert h9[j] == h; bal[i] += 1
    if j != i: nontriv[h].add(i); edges[h].add((min(i, j), max(i, j)))
print(f"order 9: balanced arcs {len(rows)}, trivial {sum(1 for (i,a,b,h), s in zip(rows, cf) if index[s] == i)}, nontrivial {sum(1 for (i,a,b,h), s in zip(rows, cf) if index[s] != i)}")
prev_rigid = [1, 3, 5, 11, 13, 17, 19, 31, 35, 39, 47, 53, 55, 57, 63, 67, 75, 85, 95, 113, 115, 121, 125]
now = {h: len(edges[h]) for h in prev_rigid}
print("fibres rigid at orders <= 8 -> nontrivial neutral edges at order 9:", now)
print("fibres below 125 that remain rigid at order 9 (hence at every order, by the floor f(10) = 125):", [h for h in prev_rigid if h < 125 and now[h] == 0])
# strong classes
import numpy as np
def strongq(s):
    n = 9; bits = ''.join(format(ord(c) - 63, '06b') for c in s[2:]); A = [[int(bits[i*n+j]) for j in range(n)] for i in range(n)]
    def reach(v, A):
        r = {v}; f = [v]
        while f:
            x = f.pop()
            for y in range(n):
                if A[x][y] and y not in r: r.add(y); f.append(y)
        return r
    AT = [[A[j][i] for j in range(n)] for i in range(n)]
    return len(reach(0, A)) == n and len(reach(0, AT)) == n, sorted(sum(r) for r in A)
raw = subprocess.run(['gentourng', '-q', '-z', '9'], capture_output=True, text=True).stdout.split()
top = sorted(((c, i) for i, c in bal.items()), reverse=True)[:12]
best = []
for c, i in top:
    st, sc = strongq(raw[i])
    if st: best.append((c, i, h9[i], sc))
print("neutral-richest strong classes of order 9 (balanced arcs, index, h, scores):", best[:4])
print("largest fibres by nontrivial neutral edges at order 9 (h, edges):", sorted(((len(e), h) for h, e in edges.items()), reverse=True)[:8])
