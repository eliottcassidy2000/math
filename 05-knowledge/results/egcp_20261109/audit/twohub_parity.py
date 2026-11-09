"""cover.py 'two-hub' family TH(G, t): hubs a, b non-adjacent, both joined to all of V(G); t pendant leaves at a.
d(TH) = 2^(n-1) + 2^f - 2^t + d(G) (derived by hand; checked below). Is TH(G, t) even for every even G, t >= 1?
cover.py uses it for every G in E(f) (values of even graphs on f vertices), 1 <= t <= n - 3, and adds the values to
E0, C and L."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_lib import *
def TH(G, t):
    f = G[0]; a, b = f, f + 1
    E = edges_of(G) + [(a, v) for v in range(f)] + [(b, v) for v in range(f)] + [(a, f + 2 + i) for i in range(t)]
    return graph(f + 2 + t, E)
def all_graphs(n):
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    for code in range(1 << len(pairs)):
        yield graph(n, [pairs[e] for e in range(len(pairs)) if (code >> e) & 1])
badv = 0; odd_cases = {}; tot = 0
for f in range(1, 6):
    for G in all_graphs(f):
        if not is_even(G): continue
        for t in range(1, 4):
            J = TH(G, t); n = J[0]; tot += 1
            if dom_count(J) != 2 ** (n - 1) + 2 ** f - 2 ** t + dom_count(G): badv += 1
            if not is_even(J):
                key = (f, t, dom_count(G))
                odd_cases.setdefault(key, edges_of(G))
print("TH(G,t), all labelled even G on f <= 5 vertices, t <= 3: %d cases, value failures %d" % (tot, badv))
print("odd TH(G,t) found for (f, t, d(G)) ->", len(odd_cases), "classes; examples:")
for k in sorted(odd_cases)[:12]:
    print("  f=%d t=%d d(G)=%d  G edges %s" % (k[0], k[1], k[2], odd_cases[k]))
# does any (f, t, d) have ONLY odd realisations among labelled even G? (then cover.py would add a value unsupported by TH)
only_odd = []
for f in range(1, 6):
    vals = {}
    for G in all_graphs(f):
        if not is_even(G): continue
        d = dom_count(G)
        for t in range(1, 4):
            vals.setdefault((t, d), []).append(is_even(TH(G, t)))
    for (t, d), lst in sorted(vals.items()):
        if not any(lst): only_odd.append((f, t, d))
print("(f, t, d(G)) for which EVERY even G gives an odd TH(G, t):", only_odd)
