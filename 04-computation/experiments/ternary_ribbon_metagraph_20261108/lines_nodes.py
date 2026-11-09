#!/usr/bin/env python3
"""lines_nodes.py n -- node-level census of the metagraph lines at order n (pynauty canonical forms).
For every line {t, t^c} (enumerated by its closed end t): the isomorphism classes of both ends, the blue flag.
Writes lines_nodes_n{n}.json: per class (certificate id): H, |Aut|, self-converse flag, strong-component count,
cyclic triangles; per line multiset: (class_closed, class_open, blue) -> count."""
import sys, os, json
from math import comb
from multiprocessing import Pool
import pynauty
n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
tiles = [(a, b) for a in range(3, n + 1) for b in range(1, a - 1) if not (a == n and b == 1)] + [(n, 1)]
m = len(tiles); tile_n1 = m - 1
cross = [sum(1 << (i - 1) for i in range(b, a)) for a, b in tiles]
sig = [tiles.index((n + 1 - b, n + 1 - a)) for a, b in tiles]
def adjacency(t, comp=False):
    adj = [0] * n
    for a in range(2, n + 1): adj[a - 1] |= 1 << (a - 2)
    for k, (a, b) in enumerate(tiles):
        bit = (t >> k) & 1
        if comp: bit ^= 1
        if bit: adj[b - 1] |= 1 << (a - 1)
        else: adj[a - 1] |= 1 << (b - 1)
    return adj
def cert(adj):
    g = pynauty.Graph(n, directed=True, adjacency_dict={v: [u for u in range(n) if (adj[v] >> u) & 1] for v in range(n)})
    return pynauty.certificate(g)
def hamcount(adj):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            out = adj[v] & ~mask
            while out:
                u = (out & -out).bit_length() - 1; out &= out - 1
                dp[mask | (1 << u)][u] += c
    return sum(dp[full])
def autsize(adj):
    g = pynauty.Graph(n, directed=True, adjacency_dict={v: [u for u in range(n) if (adj[v] >> u) & 1] for v in range(n)})
    gens, g1, g2, orbits, norb = pynauty.autgrp(g)
    return int(round(g1 * 10 ** g2))
def converse(adj):
    return [sum(1 << u for u in range(n) if (adj[u] >> v) & 1) for v in range(n)]
def components(adj):
    # strong components via reachability
    reach = []
    for v in range(n):
        seen = 1 << v; stack = [v]
        while stack:
            x = stack.pop(); out = adj[x] & ~seen
            while out:
                u = (out & -out).bit_length() - 1; out &= out - 1; seen |= 1 << u; stack.append(u)
        reach.append(seen)
    comps = set()
    for v in range(n):
        comps.add(frozenset(u for u in range(n) if (reach[v] >> u) & 1 and (reach[u] >> v) & 1))
    return len(comps)
def c3(adj):
    return comb(n, 3) - sum(comb(bin(a).count('1'), 2) for a in adj)
def work(rng):
    lo, hi = rng
    from collections import Counter
    lines = Counter(); info = {}
    for u in range(lo, hi):
        t = u | (1 << tile_n1)
        a1 = adjacency(t); a2 = adjacency(t, comp=True)
        c1 = cert(a1); c2 = cert(a2)
        blue = all(((t >> k) & 1) == ((t >> sig[k]) & 1) for k in range(m))
        lines[(c1, c2, blue)] += 1
        for c, adj in ((c1, a1), (c2, a2)):
            if c not in info:
                info[c] = (hamcount(adj), autsize(adj), cert(converse(adj)) == c, components(adj), c3(adj), adj)
    return lines, info
if __name__ == '__main__':
    total = 1 << (m - 1)
    parts = 64 if total > 4096 else 1
    ranges = [(total * p // parts, total * (p + 1) // parts) for p in range(parts)]
    from collections import Counter
    lines = Counter(); info = {}
    with Pool(8) as pool:
        for L, I in pool.imap_unordered(work, ranges):
            lines.update(L); info.update(I)
    ids = {c: i for i, c in enumerate(sorted(info))}
    conv = {c: ids.get(cert(converse(v[5]))) for c, v in info.items()}
    out = {"n": n, "m": m, "classes": {ids[c]: {"H": v[0], "aut": v[1], "sc": v[2], "k": v[3], "c3": v[4], "adj": v[5], "conv": conv[c]} for c, v in info.items()},
           "lines": [[ids[c1], ids[c2], int(b), cnt] for (c1, c2, b), cnt in lines.items()]}
    json.dump(out, open(f"lines_nodes_n{n}.json", "w"))
    print(f"n={n}: {total} lines, {len(info)} classes touched, {len(lines)} distinct (closed class, open class, blue) triples")
