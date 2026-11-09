"""Audit library for THM-5040 (written for the audit; shares NO code with egcp_lib.py).

Graphs are (n, adj) with adj a tuple of neighbour bitmasks.
- dom_count: brute force over all vertex subsets (subset DP of the dominated set).
- dom_count_ie: an independent second method, inclusion-exclusion over undominated sets.
- Automorphisms: individualisation-refinement search (1-WL colour refinement run in parallel on the two
  individualisations, so colour names are canonical across the pair). aut_generators returns one automorphism
  for every point of every basic orbit along a stabiliser chain; these generate Aut(G).
- is_even: every generator reverses an even number of edges (sign is a homomorphism).
- is_even_bruteforce: all n! permutations (small n only), used to validate is_even.
"""
import itertools


def popcount(x):
    return bin(x).count("1")


def graph(n, edges=()):
    adj = [0] * n
    for u, v in edges:
        assert u != v and 0 <= u < n and 0 <= v < n
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return (n, tuple(adj))


def edges_of(g):
    n, adj = g
    return [(u, v) for u in range(n) for v in range(u + 1, n) if (adj[u] >> v) & 1]


def complement(g):
    n, adj = g
    full = (1 << n) - 1
    return (n, tuple((full & ~adj[v]) & ~(1 << v) for v in range(n)))


def disjoint_union(*gs):
    edges = []
    off = 0
    for g in gs:
        edges += [(u + off, v + off) for u, v in edges_of(g)]
        off += g[0]
    return graph(off, edges)


def join(*gs):
    return complement(disjoint_union(*[complement(g) for g in gs]))


def pad(g, n):
    assert g[0] <= n
    return disjoint_union(g, graph(n - g[0]))


def Kn(n):
    return graph(n, [(u, v) for u in range(n) for v in range(u + 1, n)])


def En(n):
    return graph(n)


def Pk(k):
    return graph(k, [(i, i + 1) for i in range(k - 1)])


def Ck(k):
    return graph(k, [(i, (i + 1) % k) for i in range(k)])


def degrees(g):
    return [popcount(a) for a in g[1]]


def has_leaf(g):
    return 1 in degrees(g)


def has_universal(g):
    n = g[0]
    return any(d == n - 1 for d in degrees(g)) if n >= 1 else False


# ---------------- counting ----------------
def dom_count(g):
    n, adj = g
    if n == 0:
        return 1
    full = (1 << n) - 1
    closed = [adj[v] | (1 << v) for v in range(n)]
    dom = [0] * (1 << n)
    cnt = 0
    for S in range(1, 1 << n):
        low = S & -S
        dom[S] = dom[S ^ low] | closed[low.bit_length() - 1]
        if dom[S] == full:
            cnt += 1
    return cnt


def dom_count_ie(g):
    """d(G) = sum over U of (-1)^|U| 2^(n - |N[U]|)  (S must avoid N[U] for U to be undominated)."""
    n, adj = g
    closed = [adj[v] | (1 << v) for v in range(n)]
    NU = [0] * (1 << n)
    tot = 1 << n  # U = empty
    for U in range(1, 1 << n):
        low = U & -U
        NU[U] = NU[U ^ low] | closed[low.bit_length() - 1]
        term = 1 << (n - popcount(NU[U]))
        tot += -term if popcount(U) & 1 else term
    return tot


def D_count(h):
    """number of vertex sets contained in some open neighbourhood of h (empty set included if n >= 1)"""
    n, adj = h
    seen = set()
    for v in range(n):
        N = adj[v]
        sub = N
        while True:
            seen.add(sub)
            if sub == 0:
                break
            sub = (sub - 1) & N
    return len(seen)


# ---------------- automorphisms ----------------
def _nbr_lists(g):
    n, adj = g
    return [[u for u in range(n) if (adj[v] >> u) & 1] for v in range(n)]


def _refine_pair(nb, n, ind1, ind2):
    """parallel colour refinement of the same graph under two individualisations; returns (c1, c2) or None"""
    c1 = [0] * n
    c2 = [0] * n
    for i, v in enumerate(ind1):
        c1[v] = i + 1
    for i, v in enumerate(ind2):
        c2[v] = i + 1
    c1 = [(c1[v], len(nb[v])) for v in range(n)]
    c2 = [(c2[v], len(nb[v])) for v in range(n)]

    def relabel(a, b):
        if sorted(a) != sorted(b):
            return None
        keys = {k: i for i, k in enumerate(sorted(set(a)))}
        return [keys[x] for x in a], [keys[x] for x in b]

    r = relabel(c1, c2)
    if r is None:
        return None
    c1, c2 = r
    ncol = len(set(c1))
    while True:
        s1 = [(c1[v], tuple(sorted(c1[u] for u in nb[v]))) for v in range(n)]
        s2 = [(c2[v], tuple(sorted(c2[u] for u in nb[v]))) for v in range(n)]
        r = relabel(s1, s2)
        if r is None:
            return None
        c1, c2 = r
        k = len(set(c1))
        if k == ncol:
            return c1, c2
        ncol = k


def _search(g, nb, ind1, ind2):
    n, adj = g
    r = _refine_pair(nb, n, ind1, ind2)
    if r is None:
        return None
    c1, c2 = r
    if len(set(c1)) == n:
        inv2 = {c: y for y, c in enumerate(c2)}
        p = [inv2[c1[x]] for x in range(n)]
        for u in range(n):
            pu = p[u]
            a = adj[u]
            for v in range(u + 1, n):
                if ((a >> v) & 1) != ((adj[pu] >> p[v]) & 1):
                    return None
        return p
    # branch on the first vertex (in c1) of a smallest non-singleton cell
    cells = {}
    for x in range(n):
        cells.setdefault(c1[x], []).append(x)
    cell = min((c for c in cells.values() if len(c) > 1), key=len)
    x = cell[0]
    col = c1[x]
    for y in range(n):
        if c2[y] == col:
            p = _search(g, nb, ind1 + [x], ind2 + [y])
            if p is not None:
                return p
    return None


def find_automorphism(g, ind1, ind2):
    """an automorphism p with p[ind1[i]] = ind2[i] for all i, or None"""
    return _search(g, _nbr_lists(g), list(ind1), list(ind2))


def aut_generators(g, base_first=()):
    """generators of Aut(g): for a stabiliser chain along base = base_first + remaining vertices, one
    automorphism for each point of each basic orbit (a generating set by the Schreier argument)."""
    n = g[0]
    nb = _nbr_lists(g)
    order = list(base_first) + [v for v in range(n) if v not in base_first]
    gens = []
    base = []
    for v in order:
        r = _refine_pair(nb, n, base, base)
        if len(set(r[0])) == n:  # pointwise stabiliser of base is trivial
            break
        for w in range(n):
            if w == v or w in base or r[0][w] != r[0][v]:
                continue
            p = _search(g, nb, base + [v], base + [w])
            if p is not None:
                gens.append(p)
        base.append(v)
    return gens


def stabiliser_generators(g, z):
    """generators of Aut(g)_z"""
    n = g[0]
    nb = _nbr_lists(g)
    gens = []
    base = [z]
    for v in [u for u in range(n) if u != z]:
        r = _refine_pair(nb, n, base, base)
        if len(set(r[0])) == n:
            break
        for w in range(n):
            if w == v or w in base or r[0][w] != r[0][v]:
                continue
            p = _search(g, nb, base + [v], base + [w])
            if p is not None:
                gens.append(p)
        base.append(v)
    return gens


def edge_sign(g, p):
    rev = 0
    for u, v in edges_of(g):
        if p[u] > p[v]:
            rev ^= 1
    return -1 if rev else 1


def perm_sign(p):
    n = len(p)
    seen = [False] * n
    s = 1
    for i in range(n):
        if not seen[i]:
            j = i
            L = 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def is_even(g):
    return all(edge_sign(g, p) == 1 for p in aut_generators(g))


def is_automorphism(g, p):
    n, adj = g
    return all(((adj[u] >> v) & 1) == ((adj[p[u]] >> p[v]) & 1) for u in range(n) for v in range(n))


def all_automorphisms_bruteforce(g):
    n = g[0]
    return [list(p) for p in itertools.permutations(range(n)) if is_automorphism(g, p)]


def is_even_bruteforce(g):
    return all(edge_sign(g, p) == 1 for p in all_automorphisms_bruteforce(g))


# ---------------- the THM-5040 constructions, re-implemented from the text ----------------
def cone(g):
    return join(graph(1), g)


def e2_multijoin(g, t):
    return join(*([En(2)] * t + [g]))


def cocktail(n):
    assert n % 2 == 0
    return complement(graph(n, [(2 * i, 2 * i + 1) for i in range(n // 2)]))


def SR(g, q):
    """Lemma 7: hub a adjacent to V(G) - q, plus pendant leaf u at a"""
    n, adj = g
    assert popcount(adj[q]) == 1
    a, u = n, n + 1
    return graph(n + 2, edges_of(g) + [(a, v) for v in range(n) if v != q] + [(a, u)])


def jump(H, z):
    """Lemma 8: comp(H) on V(H), hub a adjacent to V(H) - z, k = rho - deg_H(z) leaves at a"""
    rho, adjH = H
    h = popcount(adjH[z])
    k = rho - h
    cH = complement(H)
    a = rho
    edges = edges_of(cH) + [(a, v) for v in range(rho) if v != z] + [(a, rho + 1 + i) for i in range(k)]
    return graph(rho + 1 + k, edges)


def path_plus_K1(rho):
    return pad(Pk(rho - 1), rho)


def spider(legs):
    """centre 0; legs listed in order; returns (H, list of leg ends)"""
    edges = []
    ends = []
    nxt = 1
    for L in legs:
        prev = 0
        for _ in range(L):
            edges.append((prev, nxt))
            prev = nxt
            nxt += 1
        ends.append(prev)
    return graph(nxt, edges), ends


def M(n):
    return 2 ** n - n if n % 2 else 2 ** n - n - 1


def X(n):
    return 2 ** n - n - 2 if (n % 2 == 1 and n >= 5) else None


def O(n):
    s = set(range(1, M(n) + 1, 2)) - {3, 7, 13}
    if X(n) is not None:
        s.discard(X(n))
    return s


def Tj(m):
    if m % 2 == 1:
        return set(range(1, m - 1, 2))
    return set(range(1, m - 2, 2)) | {m + 1}


def from_hex(n, e):
    edges = [(int(e[i], 16), int(e[i + 1], 16)) for i in range(0, len(e), 2)]
    return graph(n, edges)
