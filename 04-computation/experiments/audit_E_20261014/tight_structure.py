#!/usr/bin/env python3
"""Audit E (2026-10-14), HYP-9279 'Sharpness' paragraph: structure of the vertex-wise tight strong tournaments
(min_v max_w min(P, 1-P) = 1/3 exactly) for n = 5..8, by own DFS path enumeration (digraph6 from gentourng -z -c).
For each tight vertex v: its exact partners (P = 1/3 or 2/3), whether v and two exact partners form a cyclic triangle,
how many arcs leave that triangle towards the rest, whether the rest is transitive, and whether T is transitive-plus-one."""
import subprocess
from fractions import Fraction as F
from itertools import combinations

def d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63; bits += [(v >> (5 - t)) & 1 for t in range(6)]
    return n, [[bits[i * n + j] for j in range(n)] for i in range(n)]

def tallies(n, A):
    H = 0; bef = [[0] * n for _ in range(n)]
    nb = [[j for j in range(n) if A[i][j]] for i in range(n)]
    def rec(p, used):
        nonlocal H
        if len(p) == n:
            H += 1
            for i in range(n):
                for j in range(i + 1, n): bef[p[i]][p[j]] += 1
            return
        for v in nb[p[-1]]:
            if not used >> v & 1:
                p.append(v); rec(p, used | 1 << v); p.pop()
    for s in range(n): rec([s], 1 << s)
    return H, bef

def transitive(A, S):
    sc = sorted(sum(A[v][w] for w in S) for v in S)
    return sc == list(range(len(S)))

for n in range(5, 9):
    lines = subprocess.run(['gentourng', '-q', '-z', '-c', str(n)], capture_output=True, text=True).stdout.split()
    tight = 0; summary = {}
    for s in lines:
        _, A = d6(s)
        H, bef = tallies(n, A)
        bal = lambda v, w: min(F(bef[v][w], H), 1 - F(bef[v][w], H))
        bv = [max(bal(v, w) for w in range(n) if w != v) for v in range(n)]
        if min(bv) != F(1, 3): continue
        tight += 1
        tpo = any(transitive(A, [u for u in range(n) if u != z]) for z in range(n))
        descr = []
        for v in range(n):
            if bv[v] != F(1, 3): continue
            ex = [w for w in range(n) if w != v and bal(v, w) == F(1, 3)]
            tri = None
            for a, b in combinations(ex, 2):
                if (A[v][a] and A[a][b] and A[b][v]) or (A[v][b] and A[b][a] and A[a][v]):
                    rest = [u for u in range(n) if u not in (v, a, b)]
                    back = sum(A[t][u] for t in (v, a, b) for u in rest)
                    tri = (back, transitive(A, rest))
            descr.append((v, len(ex), tri))
        key = (tpo, tuple(sorted((d[1], d[2]) for d in descr)))
        summary[key] = summary.get(key, 0) + 1
    print(f"n={n}: tight tournaments {tight}")
    for (tpo, d), c in sorted(summary.items(), key=lambda kv: -kv[1]):
        print(f"   x{c}: transitive-plus-one={tpo}; tight vertices as (#exact partners, (arcs from the C_3 to the rest, rest transitive) or None if no C_3 through v and two exact partners): {list(d)}")

# ---- refined: treat the triangle {v, a, b} as a contiguous block; over all orders of the rest and block positions,
# minimise the number of backward arcs outside the triangle; report (min back arcs, block positions attaining it).
from itertools import permutations
print("refined: triangle as a contiguous block; min #backward arcs outside the triangle and the block positions (0 = start) attaining it")
for n in range(5, 9):
    lines = subprocess.run(['gentourng', '-q', '-z', '-c', str(n)], capture_output=True, text=True).stdout.split()
    stats = {}
    for s in lines:
        _, A = d6(s)
        H, bef = tallies(n, A)
        bal = lambda v, w: min(F(bef[v][w], H), 1 - F(bef[v][w], H))
        bv = [max(bal(v, w) for w in range(n) if w != v) for v in range(n)]
        if min(bv) != F(1, 3): continue
        for v in range(n):
            if bv[v] != F(1, 3): continue
            ex = [w for w in range(n) if w != v and bal(v, w) == F(1, 3)]
            tris = [{v, a, b} for a, b in combinations(ex, 2) if (A[v][a] and A[a][b] and A[b][v]) or (A[v][b] and A[b][a] and A[a][v])]
            best = None; poss = set()
            for tri in tris:
              rest = [u for u in range(n) if u not in tri]
              for perm in permutations(rest):
                for pos in range(len(rest) + 1):
                    order = list(perm[:pos]) + ['T'] + list(perm[pos:])
                    idx = {}
                    for k, x in enumerate(order):
                        if x == 'T':
                            for t in tri: idx[t] = k
                        else: idx[x] = k
                    back = sum(1 for x in range(n) for y in range(n) if A[x][y] and not (x in tri and y in tri) and idx[x] > idx[y])
                    lab = 'start' if pos == 0 else 'end' if pos == len(rest) else 'middle'
                    if best is None or back < best: best, poss = back, {lab}
                    elif back == best: poss.add(lab)
            where = sorted(poss)
            key = (best, tuple(sorted(set(where))))
            stats[key] = stats.get(key, 0) + 1
    print(f"  n={n}: tight vertices by (min backward arcs outside the C_3 block, block positions): {stats}")
