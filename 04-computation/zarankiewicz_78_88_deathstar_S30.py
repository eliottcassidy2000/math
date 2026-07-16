#!/usr/bin/env python3
"""death-star-2026-07-16-S30 (HYP-7106): the cyclic bipartite book at K_{6,7}, K_{6,8},
K_{7,8}, K_{8,8} vs Zarankiewicz Z(m,n) = floor(m/2)floor((m-1)/2)floor(n/2)floor((n-1)/2).
Kleitman proves Z at min(m,n) <= 6 (controls); 7x8 and 8x8 are OPEN — the cyclic book's
minimum vs the conjectured values, honestly reported. Layout: parts = position parities on
Z_{m+n} (sizes split ceil/floor automatically); bipartite edges = chords between parities;
grouped by sum class mod (m+n) — each class's chords remain pairwise noncrossing (subsets
of parallel classes), so class-coloring is well-posed; enumerate ALL 2^{#classes-1}
colorings; ALSO free search (per-edge annealing) to probe below class-structure if the
class minimum exceeds Z."""
from itertools import combinations
import random, sys, time

def Z(m, n):
    return (m//2)*((m-1)//2)*(n//2)*((n-1)//2)

def cyclic_bipartite(mtot, parity_of):
    """edges between the two parity classes of positions on Z_mtot, grouped by sum class."""
    edges = []
    for a in range(mtot):
        for b in range(a+1, mtot):
            if parity_of(a) != parity_of(b):
                edges.append((a, b))
    classes = {}
    for e in edges:
        classes.setdefault(sum(e) % mtot, []).append(e)
    return edges, classes

def cross(e, f):
    a, b = e; c, d = f
    return (a < c < b < d) or (c < a < d < b)

def class_min(classes, mtot):
    keys = sorted(classes)
    k = len(keys)
    within = sum(1 for s in keys for e, f in combinations(classes[s], 2) if cross(e, f))
    X = {}
    for i in range(k):
        for j in range(i+1, k):
            X[(i, j)] = sum(1 for e in classes[keys[i]] for f in classes[keys[j]] if cross(e, f))
    best, bestpages = None, None
    for mask in range(1 << (k-1)):
        pages = [0] + [(mask >> i) & 1 for i in range(k-1)]
        tot = within + sum(X[(i, j)] for (i, j) in X if pages[i] == pages[j])
        if best is None or tot < best:
            best, bestpages = tot, pages
    return best, within, {keys[i]: len(classes[keys[i]]) for i in range(k)}, bestpages

def free_anneal(edges, rounds=40000, seed=30):
    """per-edge page assignment annealing (lower-bound probe below class structure)."""
    rnd = random.Random(seed)
    E = len(edges)
    crossing_pairs = [(i, j) for i in range(E) for j in range(i+1, E) if cross(edges[i], edges[j])]
    pages = [rnd.randint(0, 1) for _ in range(E)]
    def cost():
        return sum(1 for i, j in crossing_pairs if pages[i] == pages[j])
    cur = cost()
    best = cur
    # incremental annealing
    adj = [[] for _ in range(E)]
    for i, j in crossing_pairs:
        adj[i].append(j); adj[j].append(i)
    T0 = 2.0
    for it in range(rounds):
        T = T0 * (1 - it / rounds) + 0.01
        i = rnd.randrange(E)
        delta = 0
        for j in adj[i]:
            delta += -1 if pages[j] == pages[i] else 1
        if delta <= 0 or rnd.random() < pow(2.718, -delta / T):
            pages[i] ^= 1
            cur += delta
            if cur < best: best = cur
    return best

if __name__ == "__main__":
    t0 = time.time()
    for (m, n) in [(6, 7), (6, 8), (7, 8), (8, 8)]:
        mtot = m + n
        # layout: balanced-necklace interleaving — part A = positions {floor(i*mtot/m)}
        posA = set((i * mtot) // m for i in range(m))
        if len(posA) != m:
            posA = set(round(i * mtot / m) % mtot for i in range(m))
        par = lambda x: 0 if x in posA else 1
        edges, classes = cyclic_bipartite(mtot, par)
        assert len(edges) == m * n, (len(edges), m*n)
        cm, within, sizes, pages = class_min(classes, mtot)
        # free annealing probe (3 seeds)
        fa = min(free_anneal(edges, seed=s) for s in [1, 2, 3])
        Zc = Z(m, n)
        status = "PROVED-OPT range (Kleitman)" if min(m, n) <= 6 else "OPEN Zarankiewicz"
        verdict = "== Z" if cm == Zc else (f"gap +{cm-Zc}" if cm > Zc else f"*** BELOW Z by {Zc-cm}")
        print(f"K_{m},{n} on Z_{mtot}: edges={len(edges)} classes={len(classes)} "
              f"(sizes {sorted(set(sizes.values()))}) within={within}")
        print(f"  class-coloring min = {cm} vs Z({m},{n}) = {Zc}  [{verdict}]  "
              f"free-anneal best = {fa}  [{status}]  [{time.time()-t0:.0f}s]")
        sys.stdout.flush()
