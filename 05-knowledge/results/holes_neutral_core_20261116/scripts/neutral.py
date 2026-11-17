#!/usr/bin/env python3
"""neutral.py: balanced arcs (count-preserving flips) of every class of orders 3-8, from the pair matrices of
profiles8.npz (THM-5267 (A): the flip of a -> b changes the count by M(b, a) - M(a, b)); trivial (isomorphic result) versus
nontrivial (another class of the same count) via nauty labelg; per fibre F_h (odd h <= 125): classes, strong classes,
balanced arcs, nontrivial neutral edges, components of the neutral graph, rigidity, and the response set B(h); the
neutral-richest strong classes per order; the correlation between neutral edges and |B(h)|."""
import numpy as np, subprocess
from collections import defaultdict, Counter
P = np.load('profiles8.npz')
def d6_encode(n, M):
    bits = ''.join(str(int(M[i][j])) for i in range(n) for j in range(n)); bits += '0' * ((-len(bits)) % 6)
    return '&' + chr(n + 63) + ''.join(chr(int(bits[i:i+6], 2) + 63) for i in range(0, len(bits), 6))
def canon(lines):
    out = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split()
    assert len(out) == len(lines); return out
fib = defaultdict(list)          # h -> list of (n, i)
balanced = {}; nontriv = {}; edges = defaultdict(set); Bset = defaultdict(set); richest = {}
for n in range(3, 9):
    A, h, M, B, st = P[f'A{n}'], P[f'h{n}'], P[f'M{n}'], P[f'B{n}'], P[f'strong{n}']
    k = len(h); cn = canon([d6_encode(n, A[i]) for i in range(k)]); index = {s: i for i, s in enumerate(cn)}
    flips = []; owner = []
    for i in range(k):
        fib[int(h[i])].append((n, i))
        for b in B[i]: Bset[int(h[i])].add(int(b))
        bal = [(a, b) for a in range(n) for b in range(n) if A[i][a][b] and M[i][a][b] == M[i][b][a]]
        balanced[(n, i)] = bal
        for (a, b) in bal:
            T2 = A[i].copy(); T2[a][b] = 0; T2[b][a] = 1; flips.append(d6_encode(n, T2)); owner.append((i, a, b))
    cf = canon(flips) if flips else []
    for s, (i, a, b) in zip(cf, owner):
        j = index[s]; assert int(h[j]) == int(h[i])
        if j != i: nontriv.setdefault((n, i), []).append((a, b, j)); edges[int(h[i])].add((n, min(i, j), max(i, j)))
    # richest strong classes
    best = sorted(((len(balanced[(n, i)]), i) for i in range(k) if st[i]), reverse=True)[:3]
    richest[n] = [(cnt, i, int(h[i]), sorted(int(x) for x in A[i].sum(axis=1))) for cnt, i in best]
    print(f"order {n}: classes {k}; balanced arcs total {sum(len(balanced[(n,i)]) for i in range(k))}; classes with a balanced arc {sum(1 for i in range(k) if balanced[(n,i)])}; with a nontrivial neutral move {sum(1 for i in range(k) if (n,i) in nontriv)}; richest strong classes (balanced arcs, index, h, scores): {richest[n]}")
print("fibre h: classes | strong | balanced arcs | nontrivial neutral edges | neutral-graph components among classes with a move | |B(h)| (orders <= 8)")
rows = []
for hh in range(1, 126, 2):
    if hh not in fib: continue
    cl = fib[hh]; ns = sum(1 for (n, i) in cl if P[f'strong{n}'][i]); ba = sum(len(balanced[c]) for c in cl); ne = len(edges[hh])
    # components among the classes touched by nontrivial edges
    parent = {}
    def find(x):
        while parent.get(x, x) != x: x = parent[x]
        return x
    for (n, i, j) in edges[hh]:
        a, b = find((n, i)), find((n, j))
        if a != b: parent[a] = b
    touched = set(); [touched.update({(n, i), (n, j)}) for (n, i, j) in edges[hh]]
    comps = len(set(find(x) for x in touched)); rows.append((hh, len(cl), ns, ba, ne, comps, len(Bset[hh])))
    print(f"  {hh}: {len(cl)} | {ns} | {ba} | {ne} | {comps} | {len(Bset[hh])}")
rigid = [r[0] for r in rows if r[4] == 0]
print("fibres (odd h <= 125) with NO nontrivial neutral move:", rigid)
big = [r for r in rows if r[0] >= 23]; print(f"fibres with h >= 23: {len(big)}, with a nontrivial move: {sum(1 for r in big if r[4] > 0)}; fibre with the most nontrivial edges: h = {max(rows, key=lambda r: r[4])[0]} ({max(r[4] for r in rows)} edges)")
from scipy.stats import spearmanr
x = [r[4] for r in rows]; y = [r[6] for r in rows]; print(f"Spearman(nontrivial neutral edges, |B(h)|) over {len(rows)} fibres: {spearmanr(x, y).correlation:.3f}")
