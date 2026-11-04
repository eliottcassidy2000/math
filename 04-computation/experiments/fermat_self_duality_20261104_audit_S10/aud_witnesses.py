"""aud_witnesses.py -- independent audit check (fermat_self_duality_20261104_audit_S10) of the note's sections 6-7
claims about the four order-9 route-(b) witnesses v => C3 => x => y => C3' => v (S75, S74 audit):
two C3 modules, hc = 9, |Aut| = 9, forced arc x -> y, quotient types (AT_5 exactly when v beats both x, y or both
beat v), self-converse pattern, not {-3}-self-dual.  Own code (networkx only for isomorphism tests)."""
import itertools
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

def arc(out, u, v): return (out[u] >> v) & 1

def hc_count(n, out):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    dp[1][0] = 1
    for S in range(1, full + 1, 2):
        for v in range(n):
            c = dp[S][v]
            if c:
                for w in range(n):
                    if not (S >> w) & 1 and arc(out, v, w):
                        dp[S | 1 << w][w] += c
    return sum(dp[full][v] for v in range(1, n) if arc(out, v, 0))

def G(n, out):
    g = nx.DiGraph(); g.add_nodes_from(range(n))
    g.add_edges_from((u, v) for u in range(n) for v in range(n) if arc(out, u, v))
    return g

def converse(n, out):
    op = [0] * n
    for u in range(n):
        for v in range(n):
            if arc(out, u, v): op[v] |= 1 << u
    return op

def induced(out, S):
    S = sorted(S); k = len(S); sub = [0] * k
    for a in range(k):
        for b in range(k):
            if a != b and arc(out, S[a], S[b]): sub[a] |= 1 << b
    return sub

def modules(n, out):
    res = []
    for M in range(1, (1 << n) - 1):
        if not 2 <= bin(M).count("1") <= n - 1: continue
        ok = True
        for x in range(n):
            if (M >> x) & 1: continue
            o = out[x] & M; i = M & ~out[x]
            if o and i: ok = False; break
        if ok: res.append(M)
    return res

def at(n):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) == (0, n - 1): out[j] |= 1 << i
            else: out[i] |= 1 << j
    return out

def is_sd3(n, out):
    for X in itertools.combinations(range(n), 3):
        S = [v for v in range(n) if v not in X]
        sub = induced(out, S)
        if not nx.is_isomorphic(G(len(S), sub), G(len(S), converse(len(S), sub))): return False
    return True

def noncritical(n, out):
    return [v for v in range(n) if nx.is_strongly_connected(G(n - 1, induced(out, [u for u in range(n) if u != v])))]

# vertices: v = 0, C3 = 1,2,3, x = 4, y = 5, C3' = 6,7,8
V, C, X_, Y_, Cp = 0, [1, 2, 3], 4, 5, [6, 7, 8]
order = [C, [X_], [Y_], Cp]          # T - v = C3 => x => y => C3' (all arcs forward)
rows = []
outs = {}
for vx in (1, 0):
    for vy in (1, 0):
        out = [0] * 9
        def add(u, w): out[u] |= 1 << w
        for cyc in (C, Cp):
            for t in range(3): add(cyc[t], cyc[(t + 1) % 3])
        for a in range(4):
            for b in range(a + 1, 4):
                for u in order[a]:
                    for w in order[b]: add(u, w)
        for c in C: add(V, c)
        for c in Cp: add(c, V)
        add(V, X_) if vx else add(X_, V)
        add(V, Y_) if vy else add(Y_, V)
        name = ("v->x" if vx else "x->v") + "," + ("v->y" if vy else "y->v")
        outs[name] = out
        g = G(9, out)
        hc = hc_count(9, out)
        aut = sum(1 for _ in DiGraphMatcher(g, g).isomorphisms_iter())
        mods = modules(9, out)
        nontr = [M for M in mods if sorted(bin(out[v] & M).count("1") for v in range(9) if (M >> v) & 1) != list(range(bin(M).count("1")))]
        c3mods = [M for M in mods if bin(M).count("1") == 3 and M in nontr]
        # forced arc x -> y: deleting it kills every Hamiltonian cycle
        out2 = list(out); out2[X_] &= ~(1 << Y_)
        g2 = G(9, out2)
        forced = hc_count(9, out2) == 0
        # 2-arc-strong: removing any single arc leaves a strong digraph
        arc2 = all(nx.is_strongly_connected(nx.restricted_view(g, [], [e])) for e in g.edges())
        # quotient on (C3, x, y, C3', v) -> representatives 1, 4, 5, 6, 0
        reps = [1, 4, 5, 6, 0]
        q = [0] * 5
        for a in range(5):
            for b in range(5):
                if a != b and arc(out, reps[a], reps[b]): q[a] |= 1 << b
        q_at5 = nx.is_isomorphic(G(5, q), G(5, at(5)))
        q_hc = hc_count(5, q)
        q_nc = len(noncritical(5, q))
        sc = nx.is_isomorphic(g, G(9, converse(9, out)))
        rows.append((name, hc, aut, len(mods), len(c3mods), len(nontr), forced, arc2, q_at5, q_hc, q_nc, sc, is_sd3(9, out)))

print("witness      hc  |Aut|  #modules  #C3-modules  #non-transitive  x->y forced  2-arc-strong  quotient=AT_5  quotient hc  quotient #non-critical  self-converse  {-3}-self-dual")
for r in rows:
    print(f"{r[0]:10s}  {r[1]:3d}  {r[2]:5d}  {r[3]:8d}  {r[4]:11d}  {r[5]:15d}  {str(r[6]):11s}  {str(r[7]):12s}  {str(r[8]):13s}  {r[9]:11d}  {r[10]:22d}  {str(r[11]):13s}  {r[12]}")
a = outs["v->x,v->y"]; b = outs["x->v,y->v"]
print("v beats both  vs  both beat v: converses of each other:", nx.is_isomorphic(G(9, converse(9, a)), G(9, b)))
print("all four pairwise non-isomorphic:", all(not nx.is_isomorphic(G(9, outs[p]), G(9, outs[q])) for p, q in itertools.combinations(outs, 2)))
