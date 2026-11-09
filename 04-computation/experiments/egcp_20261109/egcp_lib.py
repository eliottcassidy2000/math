"""Shared helpers for the S89 sharp even-graph Chan-Pak packet (THM-5040).

Graphs are (n, adj) with adj[v] a bitmask of neighbours.  Evenness is von Broemssen's: with the reference
orientation u -> v for u < v, every automorphism reverses an even number of edges.  Automorphisms are enumerated
exactly by backtracking (fine for the <= 13-vertex graphs used here).
"""

def M(n):
    """largest d of an even graph on n vertices"""
    return 2**n - n if n % 2 else 2**n - n - 1

def O(n):
    """the THM-5040 value set: odd k <= M(n), minus 3, 7, 13 and (n odd >= 5) X(n) = 2^n - n - 2"""
    s = set(range(1, M(n) + 1, 2)) - {3, 7, 13}
    if n % 2 and n >= 5:
        s.discard(2**n - n - 2)
    return s

def Tj(m):
    """j-set of the top band at order m: the values 2^(m-1) - j that need exactly m vertices"""
    N = m - 1
    if N % 2 == 0:
        return set(range(1, N, 2))
    return set(range(1, N - 1, 2)) | {N + 2}

def mk(n, edges):
    adj = [0] * n
    for u, v in edges:
        adj[u] |= 1 << v; adj[v] |= 1 << u
    return (n, adj)

def from_hex(n, e):
    adj = [0] * n
    for i in range(0, len(e), 2):
        u, v = int(e[i], 16), int(e[i + 1], 16); adj[u] |= 1 << v; adj[v] |= 1 << u
    return (n, adj)

def comp(g):
    n, adj = g; full = (1 << n) - 1
    return (n, [(full ^ adj[v]) & ~(1 << v) for v in range(n)])

def degs(g):
    return [bin(x).count('1') for x in g[1]]

def dcount(g):
    n, adj = g; full = (1 << n) - 1
    cl = [adj[v] | (1 << v) for v in range(n)]
    c = 0
    for S in range(1 << n):
        dom = 0; s = S
        while s:
            v = (s & -s).bit_length() - 1; s &= s - 1; dom |= cl[v]
        if dom == full:
            c += 1
    return c

def Dcount(g):
    """number of vertex sets contained in some open neighbourhood (the empty set included)"""
    n, adj = g; seen = {0} if n else set()
    for v in range(n):
        N = adj[v]; sub = N
        while True:
            seen.add(sub)
            if sub == 0: break
            sub = (sub - 1) & N
    return len(seen)

def automorphisms(g):
    n, adj = g; deg = degs(g)
    order = sorted(range(n), key=lambda v: -deg[v])
    res = []; perm = [-1] * n; used = [False] * n
    def bt(i):
        if i == n:
            res.append(perm[:]); return
        v = order[i]
        for w in range(n):
            if used[w] or deg[w] != deg[v]: continue
            if all(((adj[v] >> order[j]) & 1) == ((adj[w] >> perm[order[j]]) & 1) for j in range(i)):
                perm[v] = w; used[w] = True; bt(i + 1); used[w] = False; perm[v] = -1
    bt(0)
    return res

def edge_sign(g, p):
    n, adj = g; rev = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (adj[u] >> v) & 1 and p[u] > p[v]: rev += 1
    return -1 if rev % 2 else 1

def perm_sign(p):
    n = len(p); seen = [False] * n; s = 1
    for i in range(n):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]: seen[j] = True; j = p[j]; L += 1
            if L % 2 == 0: s = -s
    return s

def is_even(g):
    return all(edge_sign(g, p) == 1 for p in automorphisms(g))

def has_leaf(g):
    return 1 in degs(g)

def no_universal(g):
    n = g[0]
    return n >= 2 and max(degs(g)) < n - 1

def union(*gs):
    n = 0; adj = []
    for (m, a) in gs:
        adj += [x << n for x in a]; n += m
    return (n, adj)

def pad(g, n):
    m, a = g
    assert m <= n
    return (n, a[:] + [0] * (n - m))

# ---------------- the operations of THM-5040 ----------------
def cone(g):
    """K_1 v G (Lemma 4; needs |G| >= 2 and no universal vertex for evenness)"""
    n, a = g
    adj = a[:] + [(1 << n) - 1]
    for v in range(n): adj[v] |= 1 << n
    return (n + 1, adj)

def E2join(g, t):
    """E_2^(v t) v G (Lemma 5): d = 2^n - 2^g - 2t + d(G)"""
    n, a = g
    blocks = [list(range(n))] + [[n + 2 * i, n + 2 * i + 1] for i in range(t)]
    N = n + 2 * t; adj = a[:] + [0] * (2 * t)
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            for u in blocks[i]:
                for v in blocks[j]:
                    adj[u] |= 1 << v; adj[v] |= 1 << u
    return (N, adj)

def SR(g, q=None):
    """leaf stable recursion (Lemma 7): new hub a adjacent to V(G) - q plus a pendant leaf u; d -> 3*2^(n-1) + d"""
    n, a = g
    if q is None:
        q = next(v for v in range(n) if bin(a[v]).count('1') == 1)
    adj = a[:] + [0, 0]; A, U = n, n + 1
    for v in range(n):
        if v != q: adj[A] |= 1 << v; adj[v] |= 1 << A
    adj[A] |= 1 << U; adj[U] |= 1 << A
    return (n + 2, adj)

def jump(H, z):
    """jump family (Lemma 8): a + k leaves + comp(H), a adjacent to comp(H) - z, k = rho - deg_H(z);
    d = 2^(m-1) - D(H) with m = 2 rho + 1 - deg_H(z)"""
    rho, Ha = H
    h = bin(Ha[z]).count('1'); k = rho - h
    assert k >= 1
    m = 1 + k + rho
    GR = comp(H)[1]; adj = GR[:] + [0] * (1 + k); A = rho
    for v in range(rho):
        if v != z: adj[A] |= 1 << v; adj[v] |= 1 << A
    for i in range(k):
        u = rho + 1 + i; adj[A] |= 1 << u; adj[u] |= 1 << A
    return (m, adj)

def cocktail(t):
    """complement of a perfect matching on 2t vertices (Lemma 6): d = 2^(2t) - 2t - 1"""
    n = 2 * t; full = (1 << n) - 1
    return (n, [full & ~(1 << v) & ~(1 << (v ^ 1)) for v in range(n)])

def path_plus_K1(rho):
    """H = P_(rho-1) u K_1 on rho vertices (path 0-1-...-(rho-2), isolated vertex rho-1)"""
    return mk(rho, [(i, i + 1) for i in range(rho - 2)])

def spider_l12(rho):
    """H = spider with legs of lengths (rho-4, 1, 2), centre 0; returns (H, z) with z the end of the first leg"""
    l1 = rho - 4; edges = []; prev = 0
    for i in range(l1):
        edges.append((prev, 1 + i)); prev = 1 + i
    z = prev
    edges += [(0, 1 + l1), (0, 2 + l1), (2 + l1, 3 + l1)]
    return mk(rho, edges), z

def frontier(m, j):
    """Lemma 9: explicit leaf graphs for the new top-band j at order m"""
    if m % 2 == 0 and j == m - 3:
        return jump(path_plus_K1(m // 2), 0)
    if m % 2 == 0 and j == m + 1:
        H, z = spider_l12(m // 2); return jump(H, z)
    if m % 2 == 1 and j == m - 2:
        return jump(path_plus_K1((m + 1) // 2), 1)
    raise ValueError((m, j))
