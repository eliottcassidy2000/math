"""AUDIT2: finite reduction of the HYP-9610 (A) statements about M0(n) and M0(n) - 2, for ALL n.

For G on n vertices with no universal vertex, H = comp(G) has no isolated vertex and (THM-5040 Lemma 0)
    d(G) = 2^n - D(H),  D(H) = 1 + n + e(H),
where the excess e(H) = #{sets of size >= 2 contained in some open neighbourhood of H} is additive over components
(Brown-Nowakowski union rule), K_2 has e = 0 and every other connected graph has e >= 1.  G is even iff H is co-even
(B2.5).  Write H = H' u sK_2 with H' having no K_2 component; then Aut(H) = Aut(H') x Aut(sK_2) and sK_2 is co-even, so
H is co-even iff H' is.  Hence, for every n,
    M0(n) - 2 in E0(n)  <=>  some co-even H' (components connected, >= 3 vertices) with e(H') = 2 (n even) or 5 (n odd),
                             |H'| = n (mod 2) and |H'| <= n.
Growth lemma: adding a vertex to a connected graph on >= 2 vertices raises e by >= 1, and deleting a non-cut vertex
does not raise e; so every connected graph with e <= 5 has <= 7 vertices and is reached by growing connected graphs
with e <= 5.  This script enumerates them (isomorphism classes by networkx), forms every multiset of components with
total excess E <= 7, and tests co-evenness of H' with a2_check.exe (evenness of comp(H')).
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a2_lib import *
import networkx as nx

def excess(g):
    n, adj = g
    seen = set()
    for v in range(n):
        N = adj[v]; sub = N
        while True:
            if bin(sub).count('1') >= 2: seen.add(sub)
            if sub == 0: break
            sub = (sub - 1) & N
    return len(seen)

def comp_graph(g):
    n, adj = g; full = (1 << n) - 1
    return (n, [(full ^ adj[v]) & ~(1 << v) for v in range(n)])

EMAX = 5
classes = []          # (graph, e)
def add_class(g):
    G = nx_graph(g)
    for h, e in classes:
        if h[0] == g[0] and len(edges(h)) == len(edges(g)) and nx.is_isomorphic(nx_graph(h), G):
            return False
    classes.append((g, excess(g))); return True

frontier = [from_edges(2, [(0, 1)])]
level = 2
while frontier:
    new = []
    for g in frontier:
        n, adj = g
        for S in range(1, 1 << n):
            h = (n + 1, [a | ((S >> v & 1) << n) for v, a in enumerate(adj)] + [S])
            if excess(h) <= EMAX and add_class(h):
                new.append(h)
    frontier = new; level += 1
conn = [(g, e) for g, e in classes]
print("connected classes with 3 <= order and excess <= %d: %d; max order %d" % (EMAX, len(conn), max(g[0] for g, e in conn)))
for E in range(1, 6):
    print("  e = %d: %s" % (E, sorted((g[0], len(edges(g))) for g, e in conn if e == E)))

# co-evenness of each connected class and of every multiset with total excess E
def union_all(gs):
    out = (0, [])
    for g in gs: out = union(out, g)
    return out

def coeven(h):
    return run_check([(comp_graph(h), -1)])[0][1] == 'EVEN'

results = {}
for E in range(1, EMAX + 1):
    combos = []
    def rec(start, remaining, chosen):
        if remaining == 0:
            combos.append(list(chosen)); return
        for i in range(start, len(conn)):
            g, e = conn[i]
            if e <= remaining:
                chosen.append(i); rec(i, remaining - e, chosen); chosen.pop()
    rec(0, E, [])
    items = [(comp_graph(union_all([conn[i][0] for i in c])), -1) for c in combos]
    res = run_check(items) if items else []
    co = [(c, sum(conn[i][0][0] for i in c)) for c, r in zip(combos, res) if r[1] == 'EVEN']
    results[E] = (combos, co)
    orders = sorted(set(o for c, o in co))
    print("E = %d: %d component multisets, %d co-even; orders of co-even H': %s"
          % (E, len(combos), len(co), orders))
    if E in (2, 5):
        par = 0 if E == 2 else 1
        bad = [(c, o) for c, o in co if o % 2 == par]
        print("   co-even H' with E = %d and order = %d (mod 2) [would put M0(n) - 2 in E0(n)]: %d %s"
              % (E, par, len(bad), [(o, [ (conn[i][0][0], len(edges(conn[i][0]))) for i in c]) for c, o in bad][:10]))
# smallest co-even H' per (E, parity) gives the least n with 2^n - n - 1 - E in E0(n)
for E in range(1, EMAX + 1):
    for par in (0, 1):
        o = [o for c, o in results[E][1] if o % 2 == par]
        print("E = %d, |H'| = %d mod 2: least co-even order %s" % (E, par, min(o) if o else None))
