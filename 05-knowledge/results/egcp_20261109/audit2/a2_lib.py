"""AUDIT2 helpers, written from scratch (no import of egcp_lib.py).
Graphs are (n, adj) with adj[v] the neighbour bitmask of v.
"""
import os, subprocess, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
CHECK = os.path.join(HERE, 'a2_check.exe')

def from_edges(n, edges):
    adj = [0] * n
    for u, v in edges:
        assert u != v
        adj[u] |= 1 << v; adj[v] |= 1 << u
    return (n, adj)

def from_hex_pairs(n, s):
    return from_edges(n, [(int(s[i], 16), int(s[i + 1], 16)) for i in range(0, len(s), 2)])

def edges(g):
    n, adj = g
    return [(u, v) for u in range(n) for v in range(u + 1, n) if adj[u] >> v & 1]

def pad(g, m):
    n, adj = g
    assert n <= m
    return (m, adj + [0] * (m - n))

def cone(g):
    """K_1 v G: new vertex (index n) joined to all of G"""
    n, adj = g
    new = [a | (1 << n) for a in adj] + [(1 << n) - 1]
    return (n + 1, new)

def union(g, h):
    n, a = g; m, b = h
    return (n + m, a + [x << n for x in b])

def SR(g, q):
    """leaf stable recursion: hub A adjacent to V(G) - q, new pendant leaf U at A; returns (graph, U)"""
    n, adj = g
    A, U = n, n + 1
    new = adj[:] + [0, 0]
    for v in range(n):
        if v != q:
            new[A] |= 1 << v; new[v] |= 1 << A
    new[A] |= 1 << U; new[U] |= 1 << A
    return (n + 2, new), U

def threshold(y, m):
    """threshold graph on exactly m vertices with d = y (y odd, 1 <= y <= 2^m - 1): pad/cone induction"""
    assert y % 2 == 1 and 1 <= y <= 2**m - 1
    if m == 1:
        return (1, [0])
    if y <= 2**(m - 1) - 1:
        return pad(threshold(y, m - 1), m)
    return cone(threshold(y - 2**(m - 1), m - 1))

def d_brute(g):
    """dominating sets by direct check of every subset (independent of the C subset DP)"""
    n, adj = g
    full = (1 << n) - 1
    cnt = 0
    for S in range(1 << n):
        ok = True
        for v in range(n):
            if not (S >> v & 1) and not (adj[v] & S):
                ok = False; break
        if ok: cnt += 1
    return cnt

def run_check(items):
    """items: list of (graph, q or -1); returns list of (d, sign, qsign, univ, mindeg, maxdeg)"""
    inp = ''.join('%d %d %s\n' % (g[0], q, ' '.join(map(str, g[1]))) for g, q in items)
    out = subprocess.run([CHECK], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
    res = []
    for line in out:
        if not line.strip(): continue
        f = line.split()
        assert 'BADINPUT' not in f, line
        res.append((int(f[0]), f[1], f[2], int(f[3]), int(f[4]), int(f[5])))
    assert len(res) == len(items)
    return res

def nx_graph(g):
    import networkx as nx
    G = nx.Graph(); G.add_nodes_from(range(g[0])); G.add_edges_from(edges(g))
    return G

def nx_sign_odd(g, fix=None, limit=10**7):
    """independent parity engine: networkx VF2 automorphism enumeration; True iff some automorphism (fixing `fix`)
    reverses an odd number of edges.  Returns (odd?, number enumerated)."""
    from networkx.algorithms.isomorphism import GraphMatcher
    G = nx_graph(g)
    E = edges(g)
    cnt = 0
    for p in GraphMatcher(G, G).isomorphisms_iter():
        cnt += 1
        if fix is not None and p[fix] != fix:
            continue
        rev = sum(1 for u, v in E if p[u] > p[v])
        if rev % 2:
            return True, cnt
        if cnt >= limit:
            raise RuntimeError('limit')
    return False, cnt
