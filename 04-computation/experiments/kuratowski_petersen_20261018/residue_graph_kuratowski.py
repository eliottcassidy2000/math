#!/usr/bin/env python3
"""Kuratowski graphs, the Petersen graph and the finite unit-circle graphs G_q = Cay(F_(q^2), mu_(q+1)) (mac-mini-2026-10-18-unitpath).

Checks (exact finite-field arithmetic, own implementation of F_(p^k) by an irreducible polynomial):
 (1) G_2 = K_4; G_4 is the Clebsch graph srg(16, 5, 0, 2); the ten non-neighbours of 0 in G_4 are the two-sums
     zeta^a + zeta^b of the residue pentagon mu_5 and induce the Petersen graph = Kneser K(5,2) (disjoint pairs adjacent);
 (2) for every q in a range: max common neighbours of two distinct vertices (<= 2, so K_(2,3)-free), the number of
     triangles through an edge, and whether K_4 occurs (only q = 2^odd);
 (3) the Petersen graph inside G_q for odd q (a homomorphism-injective copy), searched by backtracking.
Usage: python3 residue_graph_kuratowski.py
"""
import itertools, sys
from sympy import factorint, isprime
import networkx as nx

def field(p, k):
    """F_(p^k) as tuples via an irreducible monic polynomial found by search; returns (elements, add, mul, frob_q)."""
    from itertools import product
    if k == 1:
        mod = None
    for coeffs in product(range(p), repeat=k):
        poly = list(coeffs) + [1]                    # x^k + c_(k-1) x^(k-1) + ... + c_0
        # irreducible iff no roots and (for k <= 4) no quadratic factor: test via x^(p^k) = x and gcd conditions (brute)
        if k == 1:
            break
        if irreducible(poly, p):
            mod = poly; break
    return mod

def polymulmod(a, b, mod, p):
    k = len(mod) - 1
    res = [0] * (2 * k - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: res[i + j] = (res[i + j] + x * y) % p
    for d in range(len(res) - 1, k - 1, -1):
        c = res[d]
        if c:
            for i in range(k + 1):
                res[d - k + i] = (res[d - k + i] - c * mod[i]) % p
    return tuple(res[:k])

def irreducible(poly, p):
    k = len(poly) - 1
    # brute force: no factor of degree <= k/2 (monic polys of degree d)
    from itertools import product
    for d in range(1, k // 2 + 1):
        for cs in product(range(p), repeat=d):
            f = list(cs) + [1]
            # polynomial remainder of poly by f
            r = list(poly)
            for deg in range(len(r) - 1, d - 1, -1):
                c = r[deg]
                if c:
                    for i in range(d + 1):
                        r[deg - d + i] = (r[deg - d + i] - c * f[i]) % p
            if all(x == 0 for x in r[:d]): return False
    return True

def gq(q):
    """Return (elements, adjacency sets, mu) for G_q = Cay(F_(q^2), {x : x^(q+1) = 1})."""
    fac = factorint(q); (p, e), = fac.items()
    k = 2 * e
    mod = field(p, k)
    elems = list(itertools.product(range(p), repeat=k))
    one = tuple([1] + [0] * (k - 1)); zero = tuple([0] * k)
    def mul(a, b): return polymulmod(a, b, mod, p)
    def pw(a, n):
        r = one
        while n:
            if n & 1: r = mul(r, a)
            a = mul(a, a); n >>= 1
        return r
    mu = [x for x in elems if x != zero and pw(x, q + 1) == one]
    assert len(mu) == q + 1
    idx = {x: i for i, x in enumerate(elems)}
    add = lambda a, b: tuple((x + y) % p for x, y in zip(a, b))
    adj = [set() for _ in elems]
    for i, x in enumerate(elems):
        for m in mu:
            adj[i].add(idx[add(x, m)])
    return elems, adj, mu, idx, add, mul, zero

print('(1) characteristic two', flush=True)
for q in (2, 4):
    elems, adj, mu, idx, add, mul, zero = gq(q)
    G = nx.Graph(); G.add_nodes_from(range(len(elems)))
    for i in range(len(elems)):
        for j in adj[i]:
            if i < j: G.add_edge(i, j)
    if q == 2:
        print('  G_2: %d vertices, complete graph K_4: %s' % (G.number_of_nodes(), nx.is_isomorphic(G, nx.complete_graph(4))), flush=True)
    else:
        clebsch = nx.Graph()                      # folded 5-cube: F_2^4, connection set e1..e4, e1+e2+e3+e4
        vs = list(itertools.product(range(2), repeat=4)); conn = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,1,1,1)]
        for v in vs:
            for c in conn:
                w = tuple((a + b) % 2 for a, b in zip(v, c))
                clebsch.add_edge(v, w)
        z0 = idx[zero]
        non = [v for v in G.nodes if v != z0 and v not in adj[z0]]
        P = G.subgraph(non)
        # two-sums of the pentagon mu_5
        sums = {add(mu[a], mu[b]) for a in range(5) for b in range(a + 1, 5)}
        kneser = nx.Graph()
        pairs = list(itertools.combinations(range(5), 2))
        kneser.add_nodes_from(pairs)
        for s, t in itertools.combinations(pairs, 2):
            if not set(s) & set(t): kneser.add_edge(s, t)
        lab = {idx[add(mu[a], mu[b])]: (a, b) for a, b in pairs}
        same = all((lab[u], lab[v]) in kneser.edges or (lab[v], lab[u]) in kneser.edges for u, v in P.edges) and P.number_of_edges() == 15
        print('  G_4 isomorphic to the Clebsch graph: %s ; chromatic number of Clebsch = 4 is classical' % nx.is_isomorphic(G, clebsch))
        print('  non-neighbours of 0: %d = the ten two-sums of mu_5: %s ; they induce the Petersen graph: %s ; adjacency = disjoint pairs (Kneser K(5,2)): %s'
              % (len(non), set(elems[v] for v in non) == sums, nx.is_isomorphic(P, nx.petersen_graph()), same), flush=True)
print('(2) local Kuratowski data of G_q', flush=True)
rows = []
for q in [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32]:
    elems, adj, mu, idx, add, mul, zero = gq(q)
    z0 = idx[zero]
    # common neighbours of 0 and every other vertex
    mx = max(len(adj[z0] & adj[v]) for v in range(len(elems)) if v != z0)
    tri = len([v for v in adj[z0] if v in adj[idx[mu[0]]]]) if True else None
    nbr = sorted(adj[z0])
    k4 = any(all(b in adj[a] for a, b in itertools.combinations(c, 2)) for c in itertools.combinations(nbr, 3))
    rows.append((q, mx, tri, k4))
    print('  q = %2d: max common neighbours of distinct vertices %d (K_(2,3)-free: %s) ; common neighbours of an edge %d ; K_4 present: %s'
          % (q, mx, mx <= 2, tri, k4), flush=True)
print('(3) the Petersen graph inside G_q, symmetry-reduced backtracking (vertex 0 -> 0 and 1 -> 1 by transitivity)')
# Petersen: outer cycle 0-1-2-3-4-0, spokes i - (i+5), inner pentagram 5-7-9-6-8-5
def petersen_in(q):
    elems, adj, mu, idx, add, mul, zero = gq(q)
    z0 = idx[zero]; o1 = idx[mu[0]] if mu[0] != zero else None
    one = idx[tuple([1] + [0] * (len(zero) - 1))]
    N = adj
    count = 0
    for v2 in N[one] - {z0}:
        for v3 in N[v2] - {z0, one}:
            for v4 in (N[v3] & N[z0]) - {one, v2}:
                used = {z0, one, v2, v3, v4}
                for v5 in N[z0] - used:
                    for v7 in (N[v2] & N[v5]) - used - {v5}:
                        for v9 in (N[v4] & N[v7]) - used - {v5, v7}:
                            for v6 in (N[one] & N[v9]) - used - {v5, v7, v9}:
                                for v8 in (N[v3] & N[v6] & N[v5]) - used - {v5, v7, v9, v6}:
                                    return (z0, one, v2, v3, v4, v5, v6, v7, v8, v9), elems
    return None, elems
for q in [3, 5, 7, 9, 11, 13, 17, 19, 23, 25, 27, 29, 31]:
    emb, elems = petersen_in(q)
    print('  q = %2d: Petersen graph is a subgraph of G_q: %s%s' % (q, emb is not None, ('' if emb is None else ' ; e.g. images of vertices 0..9: ' + str([elems[v] for v in emb]))), flush=True)
