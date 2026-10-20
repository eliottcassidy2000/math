#!/usr/bin/env python3
"""Even n = 6: tau_C(P) = # types of switching classes of tournaments fixed by (a conjugate of) P, versus
tau_EE(P) = # types of even Euler graphs with Aut containing a conjugate of P.  Classes enumerated directly:
a class is P-invariant iff for every g in P, gT - T is a cut.  Types via the canonical form min over switchings
of labelg(digraph6)."""
import sys, itertools, subprocess
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
sys.path.insert(0, '.')
from subgroups_fix import cyc, mul, closure, pair_orbits, graph6, labelg, to_d6, subgroup_reps
n = 6
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
idx = {p: t for t, p in enumerate(pairs)}
cuts = []
for A in range(1, 2 ** (n - 1)):          # subsets containing vertex 0 excluded by symmetry: use all A with bit0 = 0
    cut = 0
    for t, (i, j) in enumerate(pairs):
        if ((A >> i) & 1) != ((A >> j) & 1): cut |= 1 << t
    cuts.append(cut)
cutset = set(cuts) | {0}
def act(g, T):   # T as bitmask: bit t = 1 iff i -> j for pair (i,j)
    R = 0
    for t, (i, j) in enumerate(pairs):
        a, b = g[i], g[j]
        bit = (T >> t) & 1
        if a < b: R |= bit << idx[(a, b)]
        else: R |= (1 - bit) << idx[(b, a)]
    return R
def d6(T):
    arcs = [((i, j) if (T >> t) & 1 else (j, i)) for t, (i, j) in enumerate(pairs)]
    return to_d6(n, arcs)
# even/odd table of graphs
out = subprocess.run(['geng', '-q', str(n)], capture_output=True, text=True).stdout.split()
canon = labelg(out); even_of = {}
for s, c in zip(out, canon):
    G = nx.from_graph6_bytes(s.encode()); E = [(min(u, v), max(u, v)) for u, v in G.edges()]
    ev = True
    for m in GraphMatcher(G, G).isomorphisms_iter():
        if sum(1 for u, v in E if m[u] > m[v]) % 2: ev = False; break
    even_of[c] = ev
# all classes: representatives with pairs (0,j) oriented 0->j? no: class rep = T with bits for pairs (0,j) all 1 (switch at the vertices j)
reps_classes = []
for T in range(2 ** len(pairs)):
    if all((T >> idx[(0, j)]) & 1 for j in range(1, n)): reps_classes.append(T)
print("classes:", len(reps_classes))
for name, P in subgroup_reps(n).items():
    inv = [T for T in reps_classes if all(((act(g, T) ^ T) in cutset) for g in P)]
    # canonical form of each invariant class
    forms = set()
    for T in inv:
        lines = [d6(T ^ c) for c in cutset]
        forms.add(min(labelg(lines)))
    # even Euler graphs with Aut >= P
    orbs = pair_orbits(P, n); glines = []
    for choice in itertools.product([0, 1], repeat=len(orbs)):
        edges = [e for c, O in zip(choice, orbs) if c for e in O]
        deg = [0] * n
        for u, v in edges: deg[u] += 1; deg[v] += 1
        if all(d % 2 == 0 for d in deg): glines.append(graph6(n, edges))
    gc = set(labelg(glines)); ee = sum(1 for c in gc if even_of[c])
    print(f"  P = {name}: P-invariant classes (labelled) {len(inv)}, tau_C = {len(forms)}, tau_EE = {ee} {'' if len(forms) == ee else '<-- DIFFER'}")
