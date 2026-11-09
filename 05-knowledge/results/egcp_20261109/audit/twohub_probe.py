"""Probe cover.py's undocumented 'two-hub' family: d = 2^(n-1) + 2^f - 2^t + d(G), G in E(f), n = f + t + 2.
Lists small even graphs with that d, no universal vertex and a leaf, then tests a candidate construction."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_lib import *
def all_graphs(n):
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    for code in range(1 << len(pairs)):
        yield graph(n, [pairs[e] for e in range(len(pairs)) if (code >> e) & 1])
# candidate: hubs a, b NON-adjacent; a ~ all of V(G) and all of T; b ~ all of V(G) and a? ... test several
def cand(G, t, variant):
    f = G[0]; a, b = f, f + 1; T = list(range(f + 2, f + 2 + t))
    E = edges_of(G)
    if variant == 1:   # a ~ V(G) u {b}; b ~ T
        E += [(a, v) for v in range(f)] + [(a, b)] + [(b, x) for x in T]
    if variant == 2:   # a ~ V(G) u T; b ~ T  (T vertices have degree 2)
        E += [(a, v) for v in range(f)] + [(a, x) for x in T] + [(b, x) for x in T]
    if variant == 3:   # a ~ V(G) u {b} u T minus nothing; b ~ T ; a not adjacent to one leaf
        E += [(a, v) for v in range(f)] + [(a, b)] + [(b, x) for x in T] + [(a, x) for x in T[1:]]
    if variant == 4:   # a ~ V(G) u T; b ~ V(G) u {a}  (b has no leaves); T leaves at a
        E += [(a, v) for v in range(f)] + [(a, x) for x in T] + [(b, v) for v in range(f)] + [(a, b)]
    if variant == 5:   # a ~ everything except b; b ~ T only
        E += [(a, v) for v in range(f)] + [(a, x) for x in T] + [(b, x) for x in T]
    if variant == 6:   # a ~ V(G) u T; b ~ V(G); a !~ b
        E += [(a, v) for v in range(f)] + [(a, x) for x in T] + [(b, v) for v in range(f)]
    return graph(f + 2 + t, E)
tests = [(Pk(1), 1), (Pk(1), 2), (En(2), 1), (En(2), 2), (Pk(3), 1), (En(3), 1), (Pk(3), 2), (pad(Pk(3),4), 1)]
for variant in range(1, 7):
    ok = True
    for G, t in tests:
        f = G[0]; n = f + t + 2
        J = cand(G, t, variant)
        want = 2 ** (n - 1) + 2 ** f - 2 ** t + dom_count(G)
        if dom_count(J) != want: ok = False
    print("variant", variant, "matches formula on all tests:", ok)
