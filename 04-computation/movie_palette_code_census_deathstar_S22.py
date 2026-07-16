#!/usr/bin/env python3
"""death-star-2026-07-16-S22 (HYP-7036): THE MOVIE PALETTE CODE [n,k,d] CENSUS.

STRUCTURE THEOREM (one paragraph, proved here): each wall is crossed EXACTLY ONCE per
period, so each wall lies on exactly one edge of the movie multigraph => the palette map
(cycle space -> F2^walls) is INJECTIVE and the palette code IS the cycle code of the movie
multigraph — a GRAPHIC code. Hence
    [n, k, d] = [#events, #events - #states + 1, girth(movie multigraph)].
(k by connectivity of the single closed walk; d = girth is the standard graphic-code
minimum distance; parallel edges give d = 2, the movie has no self-loops.)

CENSUS across the core taxonomy + the conjecture test: d vs the coherence meter kappa
(min L1 of a nonzero additive relation among the moving speeds, support <= 4, |k| <= 8).
Conjectured (S21): d is governed by coherence — coherent (small kappa / rich relations)
=> small girth; incoherent => large girth. The Moore bound (2607.14068) reads: the movie
graph's density forces girth <= O(log V) — checked against measured d.
"""
from fractions import Fraction as Fr
from math import gcd, log
from itertools import combinations, product as iproduct
from collections import defaultdict, deque
import sys, time

def movie(E):
    mov = [f for f in E if f > 0]
    evs = []
    for i, f in enumerate(mov):
        for w in range(7 * f):
            evs.append((Fr(w, 7 * f), i, w))
    evs.sort()
    pos, events = [], []
    i = 0
    while i < len(evs):
        p = evs[i][0]
        cols = set()
        while i < len(evs) and evs[i][0] == p:
            cols.add((evs[i][1], evs[i][2])); i += 1
        pos.append(p); events.append(cols)
    n = len(pos)
    cells = []
    for j in range(n):
        a = pos[j]; b = pos[j + 1] if j + 1 < n else Fr(1)
        mid = (a + b) / 2
        cells.append(tuple(int((f * mid % 1) * 7) for f in mov))
    return pos, events, cells, mov

def girth_multigraph(edges, V):
    """edges: list of (u, v) with u != v; V = #vertices. Returns girth (>= 2)."""
    seen_pairs = set()
    adj = defaultdict(list)
    best = None
    for idx, (u, v) in enumerate(edges):
        key = (min(u, v), max(u, v))
        if key in seen_pairs:
            return 2  # parallel edge => 2-cycle
        seen_pairs.add(key)
        adj[u].append(v); adj[v].append(u)
    # BFS girth (simple graph now)
    best = float('inf')
    for s in range(V):
        dist = {s: 0}
        par = {s: -1}
        dq = deque([s])
        while dq:
            u = dq.popleft()
            if dist[u] * 2 >= best - 1:
                continue
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1; par[w] = u; dq.append(w)
                elif w != par[u]:
                    c = dist[u] + dist[w] + 1
                    if c < best: best = c
        if best == 3: break
    return best if best < float('inf') else 0

def kappa(mov, supmax=4, kmax=8):
    """GLOBAL min L1 over all supports 2..supmax (bug fixed: no early return per support)."""
    found = []
    for sup in range(2, supmax + 1):
        for comb in combinations(range(len(mov)), sup):
            vals = [mov[i] for i in comb]
            for ks in iproduct(range(-kmax, kmax + 1), repeat=sup):
                if any(k == 0 for k in ks): continue
                if ks[0] < 0: continue  # sign symmetry
                if sum(k * v for k, v in zip(ks, vals)) == 0:
                    found.append(sum(abs(k) for k in ks))
    if found:
        return min(found), sum(1 for x in found if x == min(found))
    return None, 0

def census_row(E, label):
    t0 = time.time()
    pos, events, cells, mov = movie(E)
    n = len(cells)
    stid = {}
    seq = []
    for st in cells:
        if st not in stid: stid[st] = len(stid)
        seq.append(stid[st])
    V = len(stid)
    edges = [(seq[j], seq[(j + 1) % n]) for j in range(n)]
    k_dim = n - V + 1
    d = girth_multigraph(edges, V)
    km, kcnt = kappa(mov)
    # contiguous-return upper bound on d (shortest walk-segment loop)
    first = {}
    contig = None
    for j, s in enumerate(seq):
        if s in first:
            gap = j - first[s]
            if contig is None or gap < contig: contig = gap
        first[s] = j
    moore = 2 * log(V) / log(max(2, 2 * n / V - 1)) if n > V else float('nan')
    print(f"  {label}: [n,k,d] = [{n}, {k_dim}, {d}]  rate={k_dim/n:.3f}  "
          f"contig<= {contig}  kappa={km} (x{kcnt})  V={V}  moore~{moore:.0f}  [{time.time()-t0:.0f}s]")
    sys.stdout.flush()
    return dict(label=label, n=n, k=k_dim, d=d, kappa=km, V=V)

if __name__ == "__main__":
    t0 = time.time()
    print("THE MOVIE PALETTE CODE CENSUS  (palette code = cycle code of the movie graph)")
    rows = []
    print("\n-- species: EXACT DILATE (max coherence, non-primitive) --")
    rows.append(census_row([0] + [5 * f for f in [2, 3, 5, 7, 11, 13]], "5x[2,3,5,7,11,13]"))
    print("-- species: CONSECUTIVE compact (difference-coherent) --")
    for c in [10, 20, 30, 50]:
        rows.append(census_row([0] + list(range(c, c + 6)), f"[{c}..{c+5}]"))
    print("-- species: PLANTED single relation --")
    rows.append(census_row([0, 24, 48, 72, 97, 143, 201], "24+48=72 (gcd-24 sub)"))
    rows.append(census_row([0, 23, 48, 71, 97, 143, 201], "23+48=71 only"))
    print("-- species: GENERIC incoherent --")
    for c in [10, 20, 30, 50]:
        E = sorted(set([0, c, int(1.37*c)+1, int(1.91*c), int(2.83*c)+1, int(4.13*c), int(5.87*c)+1]))
        if len(E) == 7:
            rows.append(census_row(E, f"generic c={c}"))
    print("-- species: FAR BANK --")
    for t in [20, 50]:
        rows.append(census_row([0, 1, 2, 3, 4, 5, t], f"[0..5,{t}]"))
    print("-- species: THE TIGHT AP + GW-flavored --")
    rows.append(census_row([0, 1, 2, 3, 4, 5, 6], "AP [0..6]"))
    rows.append(census_row([0, 1, 2, 3, 4, 5, 12], "[0..5,12] GW-ish"))
    print("-- species: RANDOM baseline --")
    import random
    rnd = random.Random(20260716)
    for i in range(2):
        E = sorted(set([0] + rnd.sample(range(40, 300), 6)))
        if len(E) == 7:
            rows.append(census_row(E, f"random#{i}"))
    print("\nSUMMARY (d vs kappa):")
    for r in sorted(rows, key=lambda r: (r['d'], -(r['kappa'] or 99))):
        print(f"  d={r['d']:>3}  kappa={r['kappa']}  rate={r['k']/r['n']:.3f}  {r['label']}")
    print(f"[total {time.time()-t0:.1f}s]")
