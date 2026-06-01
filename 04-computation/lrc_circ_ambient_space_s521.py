#!/usr/bin/env python3
"""
lrc_circ_ambient_space_s521.py   claudebox-2026-06-01-S521

The ambient space of the LRC walk.

For any speed set, the runner sub-tournament at a generic time t is the half-turn
tournament of m=n-1 points on the FULL unit circle (positions v_i t mod 1).  Call
the set of iso-classes so realizable Circ(m).  The reduction (S521 seed):

    LRC(n)  <=>  every speed set's cyclic walk t |-> [T_v(t)] through Circ(m)
                 meets the arc-menu M(n) (the lonely classes).

with the nested structure   M(n) = A000016(m)  ⊆  Circ(m)  ⊆  A000568(m).

This script measures the middle term Circ(m) (sampling positions on the circle,
deduping raw adjacency matrices, then canonicalizing only the distinct ones) and
prints the nested counts.  Saturation is gauged by the number of new classes
found in the last 20% of samples.
"""
from __future__ import annotations
import random
from itertools import permutations

def circ_adj(xs):
    m = len(xs); adj = [[0]*m for _ in range(m)]
    for a in range(m):
        xa = xs[a]
        for b in range(m):
            if a == b: continue
            d = (xa - xs[b]) % 1.0
            if 0.0 < d < 0.5: adj[a][b] = 1
    return tuple(tuple(r) for r in adj)

def canon_brute(adj):
    m = len(adj); best = None
    for p in permutations(range(m)):
        flat = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or flat < best: best = flat
    return best

def refine_color(adj):
    m = len(adj); color = [sum(r) for r in adj]
    while True:
        sig = []
        for v in range(m):
            outc = tuple(sorted(color[w] for w in range(m) if adj[v][w]))
            inc = tuple(sorted(color[w] for w in range(m) if adj[w][v]))
            sig.append((color[v], outc, inc))
        order = sorted(set(sig)); rank = {s: i for i, s in enumerate(order)}
        nc = [rank[sig[v]] for v in range(m)]
        if nc == color: break
        color = nc
    return color

def canon_refine(adj):
    m = len(adj); color = refine_color(adj)
    cells = {}
    for v in range(m): cells.setdefault(color[v], []).append(v)
    cell_order = [cells[c] for c in sorted(cells)]
    best = None
    def gen(idx, perm):
        nonlocal best
        if idx == len(cell_order):
            p = perm
            flat = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
            if best is None or flat < best: best = flat
            return
        for q in permutations(cell_order[idx]): gen(idx + 1, perm + list(q))
    gen(0, [])
    return best

A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880]
A000016 = {3: 1, 4: 2, 5: 4, 6: 6, 7: 10, 8: 16}   # arc menu (m=3 is L=1/2 boundary)

def circ_count(m, trials, canon):
    random.seed(1000 + m)
    raw = set(); newcount_tail = 0; tail_start = int(trials * 0.8)
    for k in range(trials):
        xs = sorted(random.random() for _ in range(m))
        a = circ_adj(xs)
        if a not in raw:
            raw.add(a)
            if k >= tail_start: newcount_tail += 1
    classes = set(canon(a) for a in raw)
    return len(classes), len(raw), newcount_tail

def main():
    print("Ambient space of the LRC walk:  M(n)=A000016(m)  <=  Circ(m)  <=  A000568(m)\n")
    print(f" {'m':>2} {'menu M':>7} {'|Circ(m)|':>9} {'A000568':>8}   {'Circ/A000568':>12}  {'new-in-tail':>11}")
    for m in range(3, 9):
        trials = 200000 if m <= 7 else 400000
        canon = canon_brute if m <= 7 else canon_refine
        c, rawn, tail = circ_count(m, trials, canon)
        menu = A000016.get(m, '?')
        frac = f"{c}/{A000568[m]}"
        print(f" {m:>2} {str(menu):>7} {c:>9} {A000568[m]:>8}   {frac:>12}  {tail:>11}")
    print("\n(new-in-tail = distinct raw circular matrices first seen in the last 20% of")
    print(" samples; 0 indicates the class set is saturated / count is exact.)")

if __name__ == "__main__":
    main()
