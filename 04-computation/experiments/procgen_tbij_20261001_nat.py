#!/usr/bin/env python3
"""procgen_tbij_20261001_nat.py -- orbit data and 'natural bijection' (equivariant map) analysis.

Species compared (all on vertex set [n] = {0..n-1}, S_n acting by relabelling):
  TOUR   tournaments                                  (orbits: A000568)
  EVENG  even graphs (no odd automorphism; DFGPR 2023) (orbits: A000568)
  CLASS  switching classes of tournaments             (orbits: A049313)
  EVENE  even Euler graphs ("untwisted"; THM-4524 E1)  (orbits: A049313)
  TWOG   switching classes of graphs = two-graphs     (orbits: A002854)
  EULER  Euler graphs                                 (orbits: A002854)

An S_n-equivariant map X -> Y inducing a bijection X/S_n -> Y/S_n exists iff the bipartite graph
  x ~ y  <=>  Aut(x) is conjugate in S_n to a subgroup of Aut(y)
on orbit representatives has a perfect matching (Lemma N1 of the note).  We test x ~ y by enumerating the
Aut(x)-invariant labelled members of Y and canonically labelling them (nauty labelg).

Only nauty binaries (geng, gentourng, labelg, dreadnaut) are called as subprocesses.
"""
import itertools
import math
import os
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_tbij_20261001_lib as L  # noqa: E402

# ----------------------------------------------------------------------------------------------
# small permutation-group utilities (groups given by generators; n <= 12)
# ----------------------------------------------------------------------------------------------


def perm_mul(p, q):
    """(p*q)(i) = p[q[i]]  (apply q first)"""
    return tuple(p[i] for i in q)


def group_elements(gens, n, limit=10 ** 6):
    """closure of the generators (BFS); raises if more than `limit` elements"""
    idt = tuple(range(n))
    seen = {idt}
    frontier = [idt]
    gens = [tuple(g) for g in gens]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = perm_mul(g, x)
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
                    if len(seen) > limit:
                        raise ValueError('group too large')
        frontier = nxt
    return seen


def sympy_order(gens, n):
    from sympy.combinatorics import Permutation, PermutationGroup
    if not gens:
        return 1
    G = PermutationGroup([Permutation(list(g)) for g in gens])
    return int(G.order())


def inverse(p):
    r = [0] * len(p)
    for i, a in enumerate(p):
        r[a] = i
    return tuple(r)


def group_order(gens, n):
    try:
        return len(group_elements(gens, n, limit=20000))
    except ValueError:
        return sympy_order(gens, n)


# ----------------------------------------------------------------------------------------------
# nauty wrappers
# ----------------------------------------------------------------------------------------------


def _dread_input(n, nbrs, digraph, cells=None):
    parts = ['n=%d %sg' % (n, 'd ' if digraph else '-d ')]
    rows = []
    for i in range(n):
        rows.append('%d:%s' % (i, ' '.join(str(j) for j in nbrs[i])))
    s = parts[0] + ' ' + '; '.join(rows) + '. '
    if cells:
        s += 'f=[' + '|'.join('%d:%d' % (a, b) for (a, b) in cells) + '] '
    else:
        s += '-f '
    return s + 'x\n'


def dreadnaut_generators(objs, digraph=False, chunk=2000, cells=None, nreal=None):
    """objs: list of (n, list-of-neighbour-lists). returns list of (generator list, grpsize-string).
    cells: optional list of (first, last) vertex ranges giving the initial colour partition (same for all
    objects); nreal: if given, generators are restricted to the first nreal vertices (which must be a
    union of cells)."""
    out_all = []
    for s in range(0, len(objs), chunk):
        part = objs[s:s + chunk]
        if cells is not None and len(cells) == len(objs) and isinstance(cells[0], list):
            cpart = cells[s:s + chunk]
        else:
            cpart = [cells] * len(part)
        inp = ''.join(_dread_input(n, nb, digraph, cc) for (n, nb), cc in zip(part, cpart)) + 'q\n'
        res = subprocess.run(['dreadnaut'], input=inp, capture_output=True, text=True, check=True).stdout
        cur = []
        curline = None
        results = []
        idx = 0
        for line in res.split('\n'):
            if line.startswith('('):
                if curline is not None:
                    cur.append(curline)
                curline = line.strip()
            elif line.startswith('   ') and curline is not None and line.strip().startswith(('(', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                curline += ' ' + line.strip()
            else:
                if curline is not None:
                    cur.append(curline)
                    curline = None
                if 'grpsize=' in line:
                    gs = line.split('grpsize=')[1].split(';')[0].strip()
                    n = part[idx][0]
                    gens = [_parse_cycles(c, n) for c in cur]
                    if nreal is not None:
                        gens = [tuple(g[:nreal]) for g in gens]
                        gens = [g for g in gens if g != tuple(range(nreal))]
                    results.append((gens, gs))
                    cur = []
                    idx += 1
        assert len(results) == len(part), (len(results), len(part))
        out_all.extend(results)
    return out_all


def _parse_cycles(s, n):
    p = list(range(n))
    for cyc in s.replace(')', ')\n').split('\n'):
        cyc = ' '.join(cyc.split())
        cyc = cyc.strip()
        if not cyc:
            continue
        assert cyc[0] == '(' and cyc[-1] == ')', s
        pts = [int(t) for t in cyc[1:-1].split()]
        for a, b in zip(pts, pts[1:] + pts[:1]):
            p[a] = b
    return tuple(p)


def graph_nbrs(adj):
    n = len(adj)
    return [[j for j in range(n) if (adj[i] >> j) & 1] for i in range(n)]


def canon_graphs(adjs):
    return L.labelg_batch([L.g6_encode(a) for a in adjs])


def canon_tours(outs):
    return L.labelg_batch([L.d6_encode(t) for t in outs])


# ----------------------------------------------------------------------------------------------
# signs / twists
# ----------------------------------------------------------------------------------------------


def sgn_graph(adj, g):
    """DFGPR sign: (-1)^{# edges {u<v} of the graph with g(u) > g(v)} (g an automorphism)"""
    return L.eps_orient(adj, g)


def is_even_graph(adj, gens):
    return all(sgn_graph(adj, g) == 1 for g in gens)


# ----------------------------------------------------------------------------------------------
# orbit data
# ----------------------------------------------------------------------------------------------


class Orbit:
    __slots__ = ('rep', 'canon', 'gens', 'order', 'info')

    def __init__(self, rep, canon, gens, order, info=None):
        self.rep = rep
        self.canon = canon
        self.gens = gens
        self.order = order
        self.info = info or {}


def euler_graph_orbits(n):
    """all Euler graphs on n vertices up to iso, built from all graphs on n-1 vertices (add a vertex joined
    to the odd-degree vertices), deduplicated by canonical form."""
    if n == 1:
        reps = [(0,)]
    else:
        reps = []
        for g in L.geng(n - 1):
            odd = [i for i in range(n - 1) if L.popcount(g[i]) % 2]
            adj = list(g) + [0]
            for i in odd:
                adj[i] |= 1 << (n - 1)
                adj[n - 1] |= 1 << i
            reps.append(tuple(adj))
    can = canon_graphs(reps)
    uniq = {}
    for c, a in zip(can, reps):
        if c not in uniq:
            uniq[c] = a
    keys = sorted(uniq)
    adjs = [uniq[k] for k in keys]
    gd = dreadnaut_generators([(n, graph_nbrs(a)) for a in adjs], digraph=False)
    res = []
    for k, a, (gens, gs) in zip(keys, adjs, gd):
        for g in gens:
            assert L.apply_perm_g(a, g) == a
        o = Orbit(a, k, gens, group_order(gens, n))
        o.info['even'] = is_even_graph(a, gens)
        o.info['edges'] = L.nedges(a)
        res.append(o)
    return res


def graph_orbits(n):
    reps = L.geng(n)
    can = canon_graphs(reps)
    gd = dreadnaut_generators([(n, graph_nbrs(a)) for a in reps], digraph=False)
    res = []
    for k, a, (gens, gs) in zip(can, reps, gd):
        o = Orbit(a, k, gens, group_order(gens, n))
        o.info['even'] = is_even_graph(a, gens)
        o.info['edges'] = L.nedges(a)
        res.append(o)
    return res


def tour_nbrs(out):
    n = len(out)
    return [[j for j in range(n) if (out[i] >> j) & 1] for i in range(n)]


def tournament_orbits(n):
    reps = L.gentourng(n)
    can = canon_tours(reps)
    gd = dreadnaut_generators([(n, tour_nbrs(t)) for t in reps], digraph=True)
    res = []
    for k, t, (gens, gs) in zip(can, reps, gd):
        for g in gens:
            assert L.apply_perm_t(t, g) == t
        res.append(Orbit(t, k, gens, group_order(gens, n)))
    return res


def source_member(out, w):
    """the member of the switching class of `out` in which w is a source"""
    n = len(out)
    full = (1 << n) - 1
    inn = full & ~out[w] & ~(1 << w)
    T = L.switch(out, inn)
    assert T[w] == full & ~(1 << w)
    return T


def delete_vertex_t(out, w):
    n = len(out)
    idx = [i for i in range(n) if i != w]
    pos = {v: k for k, v in enumerate(idx)}
    new = []
    for i in idx:
        m = 0
        for j in idx:
            if (out[i] >> j) & 1:
                m |= 1 << pos[j]
        new.append(m)
    return tuple(new)


def class_canon_batch(outs, chunk=20000):
    """complete invariant of the switching class: min over w of canon(T_w - w)"""
    if not outs:
        return []
    n = len(outs[0])
    res = []
    for s0 in range(0, len(outs), chunk):
        part = outs[s0:s0 + chunk]
        strs = []
        for T in part:
            for w in range(n):
                strs.append(L.d6_encode(delete_vertex_t(source_member(T, w), w)))
        can = L.labelg_batch(strs)
        res.extend(min(can[i * n:(i + 1) * n]) for i in range(len(part)))
    return res


def sdigraph(out):
    """Babai-Cameron S-digraph on 2n vertices: v -> 2v (+), 2v+1 (-); arc i->j of T gives
    i+ -> j+, i- -> j-, j+ -> i-, j- -> i+   (a directed 4-cycle on each pair of fibres)"""
    n = len(out)
    nb = [[] for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            if i != j and (out[i] >> j) & 1:
                nb[2 * i].append(2 * j)
                nb[2 * i + 1].append(2 * j + 1)
                nb[2 * j].append(2 * i + 1)
                nb[2 * j + 1].append(2 * i)
    return nb


def class_orbits(n):
    """switching classes of tournaments up to iso.  odd n: the mod-4 Eulerian members of gentourng(n);
    even n: classes of (tournaments on n-1 vertices + a source), deduplicated by class_canon."""
    if n == 1:
        return [Orbit((0,), '1', [], 1)]
    if n % 2 == 1:
        reps = [T for T in L.gentourng(n) if L.mod4_euler(T)]
        cans = class_canon_batch(reps)
        assert len(set(cans)) == len(cans)
    else:
        small = L.gentourng(n - 1)
        full = (1 << n) - 1
        cand = []
        for t in small:
            T = list(t) + [0]
            # new vertex n-1 is a source
            T[n - 1] = full & ~(1 << (n - 1))
            cand.append(tuple(T))
        cans0 = class_canon_batch(cand)
        uniq = {}
        for c, T in zip(cans0, cand):
            if c not in uniq:
                uniq[c] = T
        cans = sorted(uniq)
        reps = [uniq[c] for c in cans]
    res = []
    gd = dreadnaut_generators([(2 * n, sdigraph(T)) for T in reps], digraph=True)
    for c, T, (gens2, gs) in zip(cans, reps, gd):
        gens = []
        for h in gens2:
            # project: fibre {2v, 2v+1} -> v ; must preserve fibres
            g = [None] * n
            for v in range(n):
                a, b = h[2 * v] // 2, h[2 * v + 1] // 2
                assert a == b
                g[v] = a
            g = tuple(g)
            if g != tuple(range(n)) and g not in gens:
                gens.append(g)
        for g in gens:
            assert L.is_cut(L.flip_vector(T) ^ L.flip_vector(L.apply_perm_t(T, g)), n)
        o = Orbit(T, c, gens, group_order(gens, n))
        res.append(o)
    return res


def twograph_canon_batch(adjs):
    """complete invariant of a graph switching class: min over w of canon(G_w - w), G_w = member with w
    isolated"""
    n = len(adjs[0])
    strs = []
    for G in adjs:
        for w in range(n):
            Gw = graph_switch(G, G[w])
            assert Gw[w] == 0
            strs.append(L.g6_encode(delete_vertex_g(Gw, w)))
    can = L.labelg_batch(strs)
    return [min(can[i * n:(i + 1) * n]) for i in range(len(adjs))]


def graph_switch(adj, U):
    n = len(adj)
    full = (1 << n) - 1
    new = list(adj)
    for i in range(n):
        other = (full & ~U) if (U >> i) & 1 else U
        other &= ~(1 << i)
        new[i] = adj[i] ^ other
    return tuple(new)


def delete_vertex_g(adj, w):
    n = len(adj)
    idx = [i for i in range(n) if i != w]
    pos = {v: k for k, v in enumerate(idx)}
    new = []
    for i in idx:
        m = 0
        for j in idx:
            if (adj[i] >> j) & 1:
                m |= 1 << pos[j]
        new.append(m)
    return tuple(new)


def twograph_orbits(n):
    """two-graphs = switching classes of graphs up to iso: classes of (graphs on n-1 vertices + isolated
    vertex), deduplicated.  Aut computed by brute-force-free method: stabiliser of the class via the
    'double cover with fibre markers' digraph (see seidel_cover)."""
    if n == 1:
        return [Orbit((0,), '1', [], 1)]
    small = L.geng(n - 1)
    cand = [tuple(list(g) + [0]) for g in small]
    cans0 = twograph_canon_batch(cand)
    uniq = {}
    for c, G in zip(cans0, cand):
        if c not in uniq:
            uniq[c] = G
    cans = sorted(uniq)
    reps = [uniq[c] for c in cans]
    res = []
    # automorphisms of the two-graph = automorphisms of its 3-uniform hypergraph of odd triples:
    # point/triple incidence graph (every triple joined to its 3 points) with colour partition
    # [points | odd triples | even triples]; the action on points is faithful.
    trip = list(itertools.combinations(range(n), 3))
    objs = []
    cells = []
    for G in reps:
        odd = [t for t in trip if (((G[t[0]] >> t[1]) & 1) + ((G[t[1]] >> t[2]) & 1) + ((G[t[0]] >> t[2]) & 1)) % 2]
        even = [t for t in trip if t not in set(odd)]
        order = odd + even
        nb = [[] for _ in range(n + len(trip))]
        for k, t in enumerate(order):
            v = n + k
            for x in t:
                nb[v].append(x)
                nb[x].append(v)
        objs.append((n + len(trip), nb))
        cc = [(0, n - 1)]
        if odd:
            cc.append((n, n + len(odd) - 1))
        if even:
            cc.append((n + len(odd), n + len(trip) - 1))
        cells.append(cc)
    if n >= 3:
        gd = dreadnaut_generators(objs, digraph=False, cells=cells, nreal=n)
    else:
        gd = [([tuple(range(n))[::-1]] if n == 2 else [], '')] * len(reps)
    for c, G, (gens, gs) in zip(cans, reps, gd):
        gens = [g for g in gens if g != tuple(range(n))]
        for g in gens:
            assert L.is_cut(L.mask_from_adj(G) ^ L.mask_from_adj(L.apply_perm_g(G, g)), n)
        res.append(Orbit(G, c, gens, group_order(gens, n)))
    return res


def seidel_cover(G):
    """digraph on 3n vertices encoding the switching class of G with explicit fibres:
    fibre of v = {3v (hub), 3v+1 (+), 3v+2 (-)}; hub -> both sheet vertices (marks the fibre);
    edge uv: u+ - v+, u- - v-; non-edge: u+ - v-, u- - v+ (both directions = undirected)."""
    n = len(G)
    nb = [[] for _ in range(3 * n)]
    for v in range(n):
        nb[3 * v] += [3 * v + 1, 3 * v + 2]
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            if (G[u] >> v) & 1:
                nb[3 * u + 1].append(3 * v + 1)
                nb[3 * u + 2].append(3 * v + 2)
            else:
                nb[3 * u + 1].append(3 * v + 2)
                nb[3 * u + 2].append(3 * v + 1)
    return nb


# ----------------------------------------------------------------------------------------------
# invariant subspaces and neighbourhoods
# ----------------------------------------------------------------------------------------------


def pair_orbits_of_group(gens, n):
    """orbits of <gens> on unordered pairs, as lists of edge indices"""
    E = L.edge_list(n)
    idx = L.edge_index(n)
    parent = list(range(len(E)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for g in gens:
        for k, (i, j) in enumerate(E):
            a, b = g[i], g[j]
            if a > b:
                a, b = b, a
            r1, r2 = find(k), find(idx[(a, b)])
            if r1 != r2:
                parent[r1] = r2
    orb = defaultdict(list)
    for k in range(len(E)):
        orb[find(k)].append(k)
    return list(orb.values())


def invariant_graph_masks(gens, n, euler_only):
    """all <gens>-invariant edge sets (as int masks); if euler_only, only those with all degrees even."""
    orbs = pair_orbits_of_group(gens, n)
    E = L.edge_list(n)
    vecs = []
    for o in orbs:
        m = 0
        for k in o:
            m |= 1 << k
        vecs.append(m)
    if not euler_only:
        basis = vecs
        for bits in range(1 << len(basis)):
            m = 0
            b = bits
            i = 0
            while b:
                if b & 1:
                    m |= basis[i]
                b >>= 1
                i += 1
            yield m
        return
    # degree-parity map: each orbit-vector -> parity vector (int mask over vertices)
    def degpar(m):
        p = 0
        for k, (i, j) in enumerate(E):
            if (m >> k) & 1:
                p ^= (1 << i) | (1 << j)
        return p
    rows = [(degpar(v), v) for v in vecs]
    # kernel of the parity map over F_2: Gaussian elimination on (parity | combination)
    kernel = []
    piv = {}  # pivot bit -> (parity, combo-mask over edges)
    for par, v in rows:
        p, c = par, v
        while p:
            hb = p.bit_length() - 1
            if hb in piv:
                pp, cc = piv[hb]
                p ^= pp
                c ^= cc
            else:
                piv[hb] = (p, c)
                break
        if p == 0:
            kernel.append(c)
    # kernel spans the invariant Euler graphs
    d = len(kernel)
    for bits in range(1 << d):
        m = 0
        b = bits
        i = 0
        while b:
            if b & 1:
                m ^= kernel[i]
            b >>= 1
            i += 1
        yield m


def invariant_euler_dim(gens, n):
    cnt = 0
    for _ in invariant_graph_masks(gens, n, True):
        cnt += 1
    return cnt


def neighbourhood_graphs(gens, n, canon_to_id, euler_only, batch=200000):
    """set of orbit ids (via canon_to_id) of <gens>-invariant graphs (Euler only if euler_only) whose
    canonical form is a key of canon_to_id"""
    found = set()
    buf = []

    def flush():
        if not buf:
            return
        can = L.labelg_batch([L.g6_encode(L.adj_from_mask(m, n)) for m in buf])
        for c in can:
            if c in canon_to_id:
                found.add(canon_to_id[c])
        buf.clear()

    for m in invariant_graph_masks(gens, n, euler_only):
        buf.append(m)
        if len(buf) >= batch:
            flush()
    flush()
    return found


# ----------------------------------------------------------------------------------------------
# bipartite matching (Hopcroft-Karp)
# ----------------------------------------------------------------------------------------------


def max_matching(left, adjlist):
    """left: list of left ids; adjlist[l] = iterable of right ids. returns dict left->right of a maximum
    matching (simple augmenting paths; sizes here are small)"""
    match_r = {}
    match_l = {}

    def try_aug(u, seen):
        for v in adjlist[u]:
            if v in seen:
                continue
            seen.add(v)
            if v not in match_r or try_aug(match_r[v], seen):
                match_r[v] = u
                match_l[u] = v
                return True
        return False

    sys.setrecursionlimit(100000)
    for u in left:
        try_aug(u, set())
    return match_l


def hall_violator(left, adjlist, matching):
    """given a maximum matching that leaves some left vertex unmatched, return a Hall violator
    (S, N(S)) with |N(S)| < |S| via alternating-path reachability (Konig)."""
    match_r = {v: u for u, v in matching.items()}
    free = [u for u in left if u not in matching]
    if not free:
        return None
    u0 = free[0]
    S = {u0}
    N = set()
    stack = [u0]
    while stack:
        u = stack.pop()
        for v in adjlist[u]:
            if v not in N:
                N.add(v)
                w = match_r.get(v)
                assert w is not None, 'augmenting path exists: matching not maximum'
                if w not in S:
                    S.add(w)
                    stack.append(w)
    assert len(N) < len(S)
    return S, N


# ----------------------------------------------------------------------------------------------
# neighbourhoods for the matching analysis
# ----------------------------------------------------------------------------------------------


def perm_order(p):
    o = 1
    for c in L.cycles(p):
        o = o * len(c) // math.gcd(o, len(c))
    return o


def cyclic_generator(gens, n):
    """if <gens> is cyclic return a generator, else None"""
    elems = group_elements(gens, n, limit=100000)
    N = len(elems)
    for e in elems:
        if perm_order(e) == N:
            return e
    return None


def standard_perm(ct, n):
    """a fixed permutation of cycle type ct (tuple of cycle lengths, sum n)"""
    p = list(range(n))
    s = 0
    for l in ct:
        for i in range(l):
            p[s + i] = s + (i + 1) % l
        s += l
    return tuple(p)


LANDAU = [1, 1, 2, 3, 4, 6, 6, 12, 15, 20, 30, 30, 60, 60, 84, 105, 140]


def group_key(gens, n, order=None):
    """conjugacy-invariant cache key: ('cyc', cycle type of a generator) for cyclic groups, else None"""
    if not gens:
        return ('cyc', tuple([1] * n))
    if order is not None and order > LANDAU[n]:
        return None
    g = cyclic_generator(gens, n)
    if g is None:
        return None
    return ('cyc', L.cycle_type(g))


def class_space_action_matrix(h, n):
    """linear map on switching classes of graphs, in the normal form 'vertex 0 isolated' (edge sets on
    vertices 1..n-1, as masks over the edge indices of K_n), induced by the relabelling h"""
    E = L.edge_list(n)
    basis = [k for k, (i, j) in enumerate(E) if i != 0]   # edges not touching 0

    def normalise(mask):
        # switch at N(0) so that 0 becomes isolated
        adj = L.adj_from_mask(mask, n)
        G = graph_switch(adj, adj[0])
        return L.mask_from_adj(G)

    cols = []
    for k in basis:
        adj = L.adj_from_mask(1 << k, n)
        img = L.mask_from_adj(L.apply_perm_g(adj, h))
        cols.append(normalise(img))
    return basis, cols


def invariant_class_masks(gens, n):
    """all <gens>-invariant switching classes of graphs, each as its normal-form graph (vertex 0
    isolated), enumerated via the kernel of (A_h - I) over the generators"""
    E = L.edge_list(n)
    basis = [k for k, (i, j) in enumerate(E) if i != 0]
    pos = {k: t for t, k in enumerate(basis)}
    # build, for each basis vector b, the vector  sum_h ( A_h b + b )  stacked (one block per generator)
    mats = [class_space_action_matrix(h, n)[1] for h in gens]
    d = len(basis)
    rows = []  # for each basis vector: concatenated image under (A_h - I) as one big int
    for t, k in enumerate(basis):
        big = 0
        for gi, cols in enumerate(mats):
            v = cols[t] ^ (1 << k)
            # re-index v (mask over K_n edges, no edge at 0) into positions of basis
            w = 0
            mm = v
            while mm:
                b = mm & -mm
                e = b.bit_length() - 1
                w |= 1 << pos[e]
                mm ^= b
            big |= w << (gi * d)
        rows.append((big, 1 << k))
    kernel = []
    piv = {}
    for big, comb in rows:
        p, c = big, comb
        while p:
            hb = p.bit_length() - 1
            if hb in piv:
                pp, cc = piv[hb]
                p ^= pp
                c ^= cc
            else:
                piv[hb] = (p, c)
                break
        if p == 0:
            kernel.append(c)
    dim = len(kernel)
    for bits in range(1 << dim):
        m = 0
        b = bits
        i = 0
        while b:
            if b & 1:
                m ^= kernel[i]
            b >>= 1
            i += 1
        yield m


def neighbourhood_twographs(gens, n, canon_to_id, batch=100000):
    found = set()
    buf = []

    def flush():
        if not buf:
            return
        can = twograph_canon_batch([L.adj_from_mask(m, n) for m in buf])
        for c in can:
            if c in canon_to_id:
                found.add(canon_to_id[c])
        buf.clear()

    for m in invariant_class_masks(gens, n):
        buf.append(m)
        if len(buf) >= batch:
            flush()
    flush()
    return found


def analyse(Xorbs, Yorbs, n, kind, verbose=False):
    """kind in {'euler','graph','twog'}: Y's ambient space. Yorbs are the target orbits (their canonical
    strings index them; only these count). Returns dict with matching data."""
    canon_to_id = {o.canon: i for i, o in enumerate(Yorbs)}
    allY = set(range(len(Yorbs)))
    cache = {}
    adj = {}
    sym = []
    for xi, x in enumerate(Xorbs):
        if x.order == 1:
            adj[xi] = allY
            continue
        sym.append(xi)
        key = group_key(x.gens, n, x.order)
        if key is not None and key in cache:
            adj[xi] = cache[key]
            continue
        gens = [standard_perm(key[1], n)] if key is not None else x.gens
        if kind == 'euler':
            nb = neighbourhood_graphs(gens, n, canon_to_id, True)
        elif kind == 'graph':
            nb = neighbourhood_graphs(gens, n, canon_to_id, False)
        elif kind == 'twog':
            nb = neighbourhood_twographs(gens, n, canon_to_id)
        else:
            raise ValueError(kind)
        if key is not None:
            cache[key] = nb
        adj[xi] = nb
    M = max_matching(sym, adj)
    out = {'nX': len(Xorbs), 'nY': len(Yorbs), 'nsym': len(sym), 'matched_sym': len(M), 'adj': adj, 'sym': sym,
           'matching': M}
    out['perfect'] = (len(Xorbs) == len(Yorbs)) and len(M) == len(sym)
    if len(M) < len(sym):
        S, Nb = hall_violator(sym, adj, M)
        out['violator'] = (sorted(S), sorted(Nb))
    return out


# ----------------------------------------------------------------------------------------------
# lean exhaustive analysis TOUR -> EVENG for larger n (streaming; iterative Kuhn matching)
# ----------------------------------------------------------------------------------------------


def kuhn_matching(left, adjlist):
    """maximum bipartite matching, iterative DFS (left processed in order of increasing degree)"""
    match_r = {}
    match_l = {}
    for u0 in sorted(left, key=lambda u: len(adjlist[u])):
        # iterative DFS for an augmenting path from u0
        seen = set()
        stack = [(u0, iter(adjlist[u0]))]
        parent = {}
        found = None
        while stack:
            u, it = stack[-1]
            advanced = False
            for v in it:
                if v in seen:
                    continue
                seen.add(v)
                parent[v] = u
                if v not in match_r:
                    found = v
                    break
                stack.append((match_r[v], iter(adjlist[match_r[v]])))
                advanced = True
                break
            if found is not None:
                break
            if not advanced:
                stack.pop()
        if found is None:
            continue
        # augment along parents
        v = found
        while True:
            u = parent[v]
            prev = match_l.get(u)
            match_l[u] = v
            match_r[v] = u
            if u == u0:
                break
            v = prev
    return match_l


def tour_eveng_exhaustive(n, log=None):
    """returns dict with counts and the matching deficit for TOUR[n] -> EVENG[n]"""
    import time as _t
    t0 = _t.time()
    # even graphs: canonical strings
    graphs = L.geng(n)
    can_g = L.labelg_batch([L.g6_encode(a) for a in graphs])
    even_ids = {}
    nG = len(graphs)
    for s0 in range(0, nG, 5000):
        part = graphs[s0:s0 + 5000]
        gd = dreadnaut_generators([(n, graph_nbrs(a)) for a in part], digraph=False)
        for k, (a, (gens, gs)) in enumerate(zip(part, gd)):
            if is_even_graph(a, gens):
                even_ids[can_g[s0 + k]] = len(even_ids)
    del graphs, can_g
    if log:
        log('even graphs %d (%.1fs)' % (len(even_ids), _t.time() - t0))
    tours = L.gentourng(n)
    sym = []
    nT = len(tours)
    for s0 in range(0, nT, 5000):
        part = tours[s0:s0 + 5000]
        gd = dreadnaut_generators([(n, tour_nbrs(t)) for t in part], digraph=True)
        for t, (gens, gs) in zip(part, gd):
            gens = [g for g in gens if g != tuple(range(n))]
            if gens:
                sym.append((t, gens, group_order(gens, n)))
    del tours
    if log:
        log('tournaments %d, symmetric %d (%.1fs)' % (nT, len(sym), _t.time() - t0))
    cache = {}
    adj = {}
    for i, (t, gens, order) in enumerate(sym):
        key = group_key(gens, n, order)
        if key is not None and key in cache:
            adj[i] = cache[key]
            continue
        g2 = [standard_perm(key[1], n)] if key is not None else gens
        nb = neighbourhood_graphs(g2, n, even_ids, False)
        if key is not None:
            cache[key] = nb
            if log:
                log('  type %s: %d even graphs (%.1fs)' % (key[1], len(nb), _t.time() - t0))
        adj[i] = nb
    M = kuhn_matching(list(range(len(sym))), adj)
    res = {'n': n, 'ntour': nT, 'neven': len(even_ids), 'nsym': len(sym), 'matched': len(M),
           'perfect': nT == len(even_ids) and len(M) == len(sym)}
    if len(M) < len(sym):
        S, Nb = hall_violator(list(range(len(sym))), adj, M)
        res['violator_sizes'] = (len(S), len(Nb))
    return res
