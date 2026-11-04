#!/usr/bin/env python3
"""audit_A (d): THM-4871 statement 4.
(a) K_{q,q} with parts Q and P u F (k = q - p) for (3,5) at n = 10 and (5,7) at n = 14: H-fixed? N <= Aut X? sign of
    every element of N (closure of the normaliser generators) on X; sign of the full automorphism group (nauty
    generators, and for (3,5) the full closure of Aut X, order 2 (5!)^2 = 28800); the side swap is odd.
(b) (3,7) at k = 3: the 2 t_7 = 4 configurations (K_{3,3} on P u F, multiplier-even C_Q, Q joined to all or none):
    H-fixed, N-even, sign-odd, pairwise non-isomorphic (certificates + VF2).
Sign convention: sgn_X(g) = (-1)^#{edges {u<v}: g(u) > g(v)}.
"""
import itertools, math
import pynauty, networkx as nx

def units(m): return [a for a in range(1, m) if math.gcd(a, m) == 1]

def normaliser(p, q, k):
    n = p + q + k
    h = list(range(n))
    for i in range(p): h[i] = (i + 1) % p
    for j in range(q): h[p + j] = p + (j + 1) % q
    gens = [tuple(h)]
    for (off, m) in ((0, p), (p, q)):
        for a in units(m):
            g = list(range(n))
            for i in range(m): g[off + i] = off + (a * i) % m
            gens.append(tuple(g))
    if k >= 2:
        g = list(range(n)); g[p + q], g[p + q + 1] = p + q + 1, p + q; gens.append(tuple(g))
        g = list(range(n))
        for i in range(k): g[p + q + i] = p + q + (i + 1) % k
        gens.append(tuple(g))
    return tuple(h), closure(gens, n)

def closure(gens, n, limit=10 ** 7):
    ident = tuple(range(n)); seen = {ident}; frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for s in gens:
                x = tuple(s[g[i]] for i in range(n))
                if x not in seen:
                    seen.add(x); nxt.append(x)
                    if len(seen) > limit: raise RuntimeError("group too large")
        frontier = nxt
    return seen

def sign(g, edges): return -1 if sum(1 for (u, v) in edges if g[u] > g[v]) & 1 else 1

def preserves(g, edgeset): return all((min(g[u], g[v]), max(g[u], g[v])) in edgeset for (u, v) in edgeset)

def check_Kqq(p, q):
    k = q - p; n = p + q + k
    A = list(range(p)) + list(range(p + q, n)); B = list(range(p, p + q))
    edges = sorted((min(a, b), max(a, b)) for a in A for b in B); E = set(edges)
    h, N = normaliser(p, q, k)
    print("(p,q)=(%d,%d) k=%d n=%d: K_{%d,%d} parts Q=%s and PuF=%s, %d edges" % (p, q, k, n, q, q, B, A, len(edges)))
    print("  H-fixed:", preserves(h, E), " |N| =", len(N), " N <= Aut X:", all(preserves(g, E) for g in N),
          " signs on N:", sorted(set(sign(g, edges) for g in N)))
    adj = {v: [] for v in range(n)}
    for (u, v) in edges: adj[u].append(v); adj[v].append(u)
    gens, grpsize1, grpsize2, orbits, norb = pynauty.autgrp(pynauty.Graph(n, directed=False, adjacency_dict=adj))
    print("  |Aut X| (nauty) = %d = 2 (q!)^2 = %d; generator signs: %s" % (int(round(grpsize1 * 10 ** grpsize2)), 2 * math.factorial(q) ** 2, [sign(tuple(g), edges) for g in gens]))
    swap = list(range(n))
    for a, b in zip(A, B): swap[a], swap[b] = b, a
    print("  side swap is an automorphism:", preserves(swap, E), " sign:", sign(swap, edges))
    if 2 * math.factorial(q) ** 2 <= 50000:
        G = closure([tuple(g) for g in gens], n)
        print("  full Aut X by closure: order %d, elements with odd sign: %d (half = %d)" % (len(G), sum(1 for g in G if sign(g, edges) < 0), len(G) // 2))
        print("  N cap Aut X by closure: %d elements, all even: %s" % (sum(1 for g in N if g in G), all(sign(g, edges) > 0 for g in N if g in G)))

def multiplier_even_classes(m):
    half = list(range(1, (m + 1) // 2)); seen, reps = set(), []
    for bits in itertools.product([0, 1], repeat=len(half)):
        C = frozenset(x for d, b in zip(half, bits) if b for x in (d, m - d))
        if C in seen: continue
        orb = {frozenset((a * c) % m for c in C) for a in units(m)}
        seen |= orb
        ed = sorted({(min(i, (i + c) % m), max(i, (i + c) % m)) for i in range(m) for c in C})
        even = all(sign([(a * i) % m for i in range(m)], ed) > 0 for a in units(m) if frozenset((a * c) % m for c in C) == C)
        if even: reps.append(C)
    return reps

def check_4b(q):
    p, k = 3, 3; n = p + q + k
    h, N = normaliser(p, q, k)
    confs = []
    for Cq in multiplier_even_classes(q):
        for cross in (0, 1):
            edges = set()
            for i in range(p):
                for t in range(k): edges.add((i, p + q + t))
            for i in range(q):
                for c in Cq:
                    j = (i + c) % q; edges.add((min(p + i, p + j), max(p + i, p + j)))
            if cross:
                for u in list(range(p)) + list(range(p + q, n)):
                    for j in range(q): edges.add((min(u, p + j), max(u, p + j)))
            edges = sorted(edges); confs.append((Cq, cross, edges))
    print("(3,%d) k=3 n=%d: %d configurations (2 e_q, e_q = %d)" % (q, n, len(confs), len(multiplier_even_classes(q))))
    certs = []
    for (Cq, cross, edges) in confs:
        E = set(edges)
        adj = {v: [] for v in range(n)}
        for (u, v) in edges: adj[u].append(v); adj[v].append(u)
        g = pynauty.Graph(n, directed=False, adjacency_dict=adj)
        gens = pynauty.autgrp(g)[0]
        Nstab = [x for x in N if preserves(x, E)]
        swap = list(range(n))
        for i in range(3): swap[i], swap[p + q + i] = p + q + i, i
        certs.append(pynauty.certificate(g))
        print("  C_Q=%s cross=%d: H-fixed %s; N-even %s (|N cap Aut| = %d); sign-odd %s (odd nauty generator exists: %s; swap P<->F sign %d, reverses %d edges)" % (
            sorted(Cq), cross, preserves(h, E), all(sign(x, edges) > 0 for x in Nstab), len(Nstab),
            any(sign(tuple(x), edges) < 0 for x in gens), any(sign(tuple(x), edges) < 0 for x in gens), sign(swap, edges), sum(1 for (u, v) in edges if swap[u] > swap[v])))
    print("  pairwise non-isomorphic (distinct certificates): %s" % (len(set(certs)) == len(certs)))
    Gs = []
    for (Cq, cross, edges) in confs:
        G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(edges); Gs.append(G)
    print("  VF2 pairwise non-isomorphic: %s" % all(not nx.is_isomorphic(Gs[i], Gs[j]) for i in range(len(Gs)) for j in range(i + 1, len(Gs))))

if __name__ == "__main__":
    check_Kqq(3, 5)
    check_Kqq(5, 7)
    check_4b(7)
