#!/usr/bin/env python3
"""profile_caps.py: age-structured caps for the 119 profile coordinates of the literal span-3 transfer.
usage: profile_caps.py states3.txt [maxfrag [minfrag]]
For a span-3 ribbon T on old sites 0..n-1 and a profile state s (fragments (head age, tail age), -1 = dead head /
free tail), x_s(T) counts the path covers of T with profile s.  Iterated tail deletion (the argument of
SHARP_GROWTH.md, applied fragment by fragment) gives
    x_s(T) <= kappa_o(s) * H(T),
where o is the orientation of the twelve possible back arcs among the last six old sites and
    kappa_o(s) = min over (base fragment F_1 with the dead head if any, merge order F_2..F_k)
                 sum over tuples (S_2..S_k) of prod_i e_i(S_i)/H(T[S_i]),
S_i ranging over the subsets of N^+(b_i) (b_i = head of F_i, a window vertex) that are disjoint from each other and
from the heads and tails of fragments merged later, S_i = empty allowed with weight 1, and e_i(S) the number of
Hamiltonian paths of T[S] ending at the tail of the path merged so far (the tail of F_{i-1}; for a free tail, any
vertex of S among the three sites just left of the window).  Only the last six old sites matter: a window vertex of
age a beats at most 3 + a vertices, all among them, and every older site beats every one of the six.
The script prints kappa(s) = max over the 4096 orientations, with an argmax orientation, for every state with at most
`maxfrag` fragments (default 3; the eight four-fragment states are slower), as exact fractions."""
import sys, itertools
from fractions import Fraction
from functools import lru_cache
states = []
for line in open(sys.argv[1]):
    parts = line.split(); idx = int(parts[0]); k = int(parts[1])
    frags = [tuple(map(int, p.strip('()').split(','))) for p in parts[2:]]
    assert len(frags) == k; states.append((idx, frags))
maxfrag = int(sys.argv[2]) if len(sys.argv) > 2 else 3
minfrag = int(sys.argv[3]) if len(sys.argv) > 3 else 1
LOCAL = 6; WIN = {0: 5, 1: 4, 2: 3}                     # age -> local index (local 5 = newest old site)
pairs = [(p, q) for q in range(LOCAL) for p in range(q) if q - p <= 3]   # 12 possible back arcs
assert len(pairs) == 12
def build(o):
    adj = [[0] * LOCAL for _ in range(LOCAL)]
    for p in range(LOCAL):
        for q in range(p + 1, LOCAL): adj[p][q] = 1
    for bit, (p, q) in enumerate(pairs):
        if o >> bit & 1: adj[p][q] = 0; adj[q][p] = 1
    return adj
def path_table(adj):
    """paths[S][v] = number of Hamiltonian paths of T[S] ending at v (S as bitmask)."""
    paths = [[0] * LOCAL for _ in range(1 << LOCAL)]
    for v in range(LOCAL): paths[1 << v][v] = 1
    for S in range(1, 1 << LOCAL):
        for v in range(LOCAL):
            c = paths[S][v]
            if not c: continue
            for w in range(LOCAL):
                if not S >> w & 1 and adj[v][w]: paths[S | 1 << w][w] += c
    return paths
def bound_for(adj, paths, frags):
    k = len(frags)
    if k <= 1: return Fraction(1)
    dead = [i for i, (h, t) in enumerate(frags) if h == -1]
    bases = dead if dead else list(range(k))
    best = None
    for base in bases:
        others = [i for i in range(k) if i != base]
        for order in itertools.permutations(others):
            seq = [base] + list(order)
            total = Fraction(0)
            def rec(i, used, acc):
                nonlocal total
                if i == len(seq): total += acc; return
                fi = frags[seq[i]]; b = WIN[fi[0]]            # head of F_i is live (ages known)
                Nplus = [v for v in range(LOCAL) if adj[b][v]]
                # vertices of fragments merged later are not in the current path
                later = set()
                for j in seq[i + 1:]:
                    h, t = frags[j]
                    if h >= 0: later.add(WIN[h])
                    if t >= 0: later.add(WIN[t])
                cand = [v for v in Nplus if not used >> v & 1 and v not in later]
                tprev = frags[seq[i - 1]][1]                   # tail of the path merged so far
                rec(i + 1, used, acc)                          # S_i empty
                for r in range(1, len(cand) + 1):
                    for sub in itertools.combinations(cand, r):
                        S = 0
                        for v in sub: S |= 1 << v
                        HS = sum(paths[S])
                        if tprev >= 0:
                            v = WIN[tprev]
                            if not S >> v & 1: continue
                            e = paths[S][v]
                        else:
                            e = sum(paths[S][v] for v in sub if v <= 2)
                        if e == 0: continue
                        rec(i + 1, used | S, acc * Fraction(e, HS))
            rec(1, 0, Fraction(1))
            if best is None or total < best: best = total
    return best
results = {}
for idx, frags in states:
    if len(frags) > maxfrag or len(frags) < minfrag: continue
    kap = Fraction(0); arg = None
    for o in range(1 << 12):
        adj = build(o); paths = path_table(adj)
        bnd = bound_for(adj, paths, frags)
        if bnd > kap: kap, arg = bnd, o
    results[idx] = (kap, arg)
    print(f"state {idx} frags {frags}: kappa = {kap} ({float(kap):.4f}) at orientation {arg:012b}", flush=True)
