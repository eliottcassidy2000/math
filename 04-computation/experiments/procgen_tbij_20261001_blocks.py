#!/usr/bin/env python3
"""procgen_tbij_20261001_blocks.py -- block-lemma certificates for 'no natural bijection CLASS -> EVENE'.

A block tournament is a list of blocks (B_1, ..., B_r) of odd sizes, each carrying a tournament invariant under
an odd-order transitive 'twist-rigid' group H_i, with all arcs between B_i and B_j (i < j) oriented B_i -> B_j.
Then H = H_1 x ... x H_r <= Aut(T) <= Aut([T]).  Block lemma (note, section 3.1): if r <= 2, or r = 3 with two
blocks of equal size, the only H-invariant even Euler graph is the empty graph, so by Lemma N1 any natural map
CLASS -> EVENE sends [T] to the empty graph.  Two non-isomorphic such classes on n vertices => no natural
bijection on types at that n.

Block kinds (all with odd-order groups):
  ('pt', 1)        a single vertex, trivial group
  ('qr', q)        Paley tournament QR_q, q prime = 3 (mod 4); H = {x -> a x + b : a a nonzero square}
                   (2-homogeneous, order q(q-1)/2)
  ('p58', q)       q prime = 5 (mod 8); H = {x -> a x + b : a in the odd part of the squares} (order q(q-1)/4);
                   tournament x -> y iff y - x in M u rM  (M = odd part of squares, r a primitive root)
  ('z3sq', 9)      Cayley tournament on Z_3^2 with S = {(1,0),(0,1),(1,1),(1,2)}; H = Z_3^2 (regular)
"""
import itertools
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_tbij_20261001_lib as L  # noqa: E402
import procgen_tbij_20261001_nat as N  # noqa: E402


def is_prime(p):
    if p < 2:
        return False
    for d in range(2, int(p ** 0.5) + 1):
        if p % d == 0:
            return False
    return True


def primitive_root(p):
    for r in range(2, p):
        if all(pow(r, (p - 1) // f, p) != 1 for f in prime_factors(p - 1)):
            return r
    return 1


def prime_factors(m):
    fs = []
    d = 2
    while d * d <= m:
        if m % d == 0:
            fs.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        fs.append(m)
    return fs


def block_data(kind, b):
    """returns (arcs on range(b) as set of (x,y), generators as tuples on range(b))"""
    if kind == 'pt':
        assert b == 1
        return set(), []
    if kind == 'qr':
        q = b
        assert is_prime(q) and q % 4 == 3
        sq = {(x * x) % q for x in range(1, q)}
        arcs = {(x, y) for x in range(q) for y in range(q) if x != y and (y - x) % q in sq}
        r = primitive_root(q)
        g1 = tuple((x + 1) % q for x in range(q))
        g2 = tuple((x * r * r) % q for x in range(q))
        return arcs, [g1, g2]
    if kind == 'p58':
        q = b
        assert is_prime(q) and q % 8 == 5
        r = primitive_root(q)
        M = {pow(r, 4 * k, q) for k in range((q - 1) // 4)}
        S = M | {(r * m) % q for m in M}
        assert len(S) == (q - 1) // 2 and all((-s) % q not in S for s in S)
        arcs = {(x, y) for x in range(q) for y in range(q) if x != y and (y - x) % q in S}
        g1 = tuple((x + 1) % q for x in range(q))
        g2 = tuple((x * pow(r, 4, q)) % q for x in range(q))
        return arcs, [g1, g2]
    if kind == 'z3sq':
        assert b == 9
        pts = [(i, j) for i in range(3) for j in range(3)]
        pos = {p: k for k, p in enumerate(pts)}
        S = {(1, 0), (0, 1), (1, 1), (1, 2)}
        arcs = set()
        for a in pts:
            for c in pts:
                if a != c and ((c[0] - a[0]) % 3, (c[1] - a[1]) % 3) in S:
                    arcs.add((pos[a], pos[c]))
        g1 = tuple(pos[((p[0] + 1) % 3, p[1])] for p in pts)
        g2 = tuple(pos[(p[0], (p[1] + 1) % 3)] for p in pts)
        return arcs, [g1, g2]
    raise ValueError(kind)


def block_kind_for(b):
    """a twist-rigid block kind of size b (odd), or None"""
    if b == 1:
        return 'pt'
    if b == 9:
        return 'z3sq'
    if is_prime(b) and b % 4 == 3:
        return 'qr'
    if is_prime(b) and b % 8 == 5:
        return 'p58'
    return None


RIGID_SIZES = [b for b in range(1, 200, 2) if block_kind_for(b) is not None]


def block_tournament(sizes):
    """tournament on n = sum(sizes) vertices and generators of H = prod H_i"""
    n = sum(sizes)
    out = [0] * n
    gens = []
    offs = []
    o = 0
    for b in sizes:
        offs.append(o)
        o += b
    for bi, b in enumerate(sizes):
        arcs, bg = block_data(block_kind_for(b), b)
        off = offs[bi]
        for (x, y) in arcs:
            out[off + x] |= 1 << (off + y)
        for g in bg:
            full = list(range(n))
            for x in range(b):
                full[off + x] = off + g[x]
            gens.append(tuple(full))
    for bi in range(len(sizes)):
        for bj in range(bi + 1, len(sizes)):
            for x in range(sizes[bi]):
                for y in range(sizes[bj]):
                    out[offs[bi] + x] |= 1 << (offs[bj] + y)
    T = tuple(out)
    # sanity: tournament, H <= Aut(T)
    for i in range(n):
        for j in range(i + 1, n):
            assert ((T[i] >> j) & 1) + ((T[j] >> i) & 1) == 1
    for g in gens:
        assert L.apply_perm_t(T, g) == T
    return T, gens


def is_odd_graph(adj):
    """True iff the graph has an automorphism reversing an odd number of edges (dreadnaut generators)"""
    n = len(adj)
    if L.nedges(adj) == 0:
        return False
    (gens, gs), = N.dreadnaut_generators([(n, N.graph_nbrs(adj))], digraph=False)
    return any(N.sgn_graph(adj, g) == -1 for g in gens)


def check_twist_rigid(kind, b):
    """every nonempty H-invariant graph on the block has an odd automorphism"""
    arcs, gens = block_data(kind, b)
    if b == 1:
        return True, 0
    cnt = 0
    for m in N.invariant_graph_masks(gens, b, False):
        if m == 0:
            continue
        cnt += 1
        if not is_odd_graph(L.adj_from_mask(m, b)):
            return False, cnt
    return True, cnt


def forced_to_empty(gens, n):
    """the only <gens>-invariant even Euler graph on n vertices is the empty graph"""
    k = 0
    for m in N.invariant_graph_masks(gens, n, True):
        if m == 0:
            continue
        k += 1
        if not is_odd_graph(L.adj_from_mask(m, n)):
            return False, k
    return True, k


def constructions(n):
    """block-size tuples of type (c) of the block lemma with twist-rigid sizes"""
    res = []
    rs = RIGID_SIZES
    if n in rs:
        res.append((n,))
    for b1 in rs:
        b2 = n - b1
        if b1 <= b2 and b2 in rs:
            res.append((b1, b2))
    for b in rs:
        b2 = n - 2 * b
        if b2 >= 1 and b2 in rs:
            res.append((b, b, b2))
    return res


# ----------------------------------------------------------------------------------------------
# DFGPR (TOUR -> EVENG) certificates
# ----------------------------------------------------------------------------------------------


def block_tournament_custom(sizes, inter):
    """tournament from blocks with explicit inter-block orientation: inter[(i,j)] = True means B_i => B_j"""
    n = sum(sizes)
    out = [0] * n
    gens = []
    offs = []
    o = 0
    for b in sizes:
        offs.append(o)
        o += b
    for bi, b in enumerate(sizes):
        arcs, bg = block_data(block_kind_for(b), b)
        off = offs[bi]
        for (x, y) in arcs:
            out[off + x] |= 1 << (off + y)
        for g in bg:
            full = list(range(n))
            for x in range(b):
                full[off + x] = off + g[x]
            gens.append(tuple(full))
    for bi in range(len(sizes)):
        for bj in range(bi + 1, len(sizes)):
            fwd = inter[(bi, bj)]
            for x in range(sizes[bi]):
                for y in range(sizes[bj]):
                    if fwd:
                        out[offs[bi] + x] |= 1 << (offs[bj] + y)
                    else:
                        out[offs[bj] + y] |= 1 << (offs[bi] + x)
    T = tuple(out)
    for g in gens:
        assert L.apply_perm_t(T, g) == T
    return T, gens


def invariant_even_graph_types(gens, n):
    """canonical forms of the <gens>-invariant graphs that have no odd automorphism"""
    masks = list(N.invariant_graph_masks(gens, n, False))
    adjs = [L.adj_from_mask(m, n) for m in masks]
    gd = N.dreadnaut_generators([(n, N.graph_nbrs(a)) for a in adjs], digraph=False)
    ev = [a for a, (g, s) in zip(adjs, gd) if all(N.sgn_graph(a, h) == 1 for h in g)]
    return set(N.canon_graphs(ev)) if ev else set()


def dfgpr_certificates(n):
    """list of Hall-violating families for TOUR[n] -> EVENG[n] built from blocks; each item:
    (label, [tournaments], union of neighbourhood types)"""
    rs = set(RIGID_SIZES)
    out = []
    if n in rs and (n - 2) in rs:
        fam = []
        T, g = block_tournament([n])
        fam.append((T, g))
        # (n-2, 1, 1): blocks B, a, b with the four inter-block patterns
        for (Ba, Bb, ab) in [(False, False, True), (False, True, True), (False, True, False), (True, True, True)]:
            # inter[(0,1)] = B => a ; inter[(0,2)] = B => b ; inter[(1,2)] = a => b
            T, g = block_tournament_custom([n - 2, 1, 1], {(0, 1): Ba, (0, 2): Bb, (1, 2): ab})
            fam.append((T, g))
        out.append(('D1 (n, n-2 rigid)', fam))
    if n % 2 == 0 and (n // 2) in rs and (n - 1) in rs:
        fam = []
        T, g = block_tournament([n // 2, n // 2])
        fam.append((T, g))
        for fwd in (True, False):
            T, g = block_tournament_custom([n - 1, 1], {(0, 1): fwd})
            fam.append((T, g))
        out.append(('D2 (n/2, n-1 rigid)', fam))
    res = []
    for label, fam in out:
        cans = N.canon_tours([T for T, g in fam])
        assert len(set(cans)) == len(cans), 'tournaments not distinct'
        nb = set()
        for T, g in fam:
            nb |= invariant_even_graph_types(g, n)
        res.append((label, len(fam), len(nb)))
    return res


def forced_tournaments(n):
    """tournaments T0 on n vertices whose block group admits only the empty even graph:
    single rigid block (n rigid) or two equal rigid blocks (n/2 rigid)"""
    rs = set(RIGID_SIZES)
    res = []
    if n in rs:
        res.append(('(%d)' % n,) + block_tournament([n]))
    if n % 2 == 0 and (n // 2) in rs:
        res.append(('(%d,%d)' % (n // 2, n // 2),) + block_tournament([n // 2, n // 2]))
    return res


def pattern_family(n, bprime, k):
    """all tournaments with blocks (bprime, 1^k): block B (rigid, its own tournament), k singletons with an
    arbitrary tournament among them and an arbitrary uniform orientation to B.  Returns (list of
    (T, gens) with distinct iso types, set of canonical forms of the H-invariant even graphs)."""
    assert bprime + k == n
    sizes = [bprime] + [1] * k
    pats = {}
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for bits in range(1 << (k + len(pairs))):
        inter = {}
        for i in range(k):
            inter[(0, 1 + i)] = bool((bits >> i) & 1)
        for t, (i, j) in enumerate(pairs):
            inter[(1 + i, 1 + j)] = bool((bits >> (k + t)) & 1)
        T, g = block_tournament_custom(sizes, inter)
        pats[T] = g
    Ts = list(pats)
    cans = N.canon_tours(Ts)
    uniq = {}
    for c, T in zip(cans, Ts):
        uniq.setdefault(c, (T, pats[T]))
    gens = next(iter(uniq.values()))[1]
    G = invariant_even_graph_types(gens, n)
    return uniq, G


def dfgpr_certificates_general(n, kmax=4):
    """Hall violators for TOUR[n] -> EVENG[n]: a forced T0 plus the pattern family of (n-k, 1^k).
    returns list of (T0 label, k, #patterns, #graph types, violation?)"""
    rs = set(RIGID_SIZES)
    out = []
    F = forced_tournaments(n)
    if not F:
        return out
    for (lab, T0, g0) in F:
        assert invariant_even_graph_types(g0, n) <= set(N.canon_graphs([tuple([0] * n)]))
    for k in range(1, kmax + 1):
        if (n - k) not in rs or n - k < 3:
            continue
        uniq, G = pattern_family(n, n - k, k)
        c0s = set(N.canon_tours([T0 for (lab, T0, g0) in F]))
        P = set(uniq) - c0s
        nforced = len(c0s)
        empty = N.canon_graphs([tuple([0] * n)])[0]
        assert empty in G
        viol = len(P) + nforced > len(G)
        out.append(([lab for (lab, T0, g0) in F], k, len(P), nforced, len(G), viol))
    return out


def dfgpr_block_hall(n, kmax=4):
    """Hall search for TOUR[n] -> EVENG[n] restricted to block-type tournaments: all patterns of (b', 1^k)
    (b' rigid, 1 <= k <= kmax) and the forced tournaments (n), (n/2, n/2).  Each tournament T gets the
    neighbourhood bound N_H(T) = types of H-invariant even graphs (intersection over its constructions);
    the true neighbourhood is contained in it.  Returns (#left, max matching, violator sizes or None)."""
    rs = set(RIGID_SIZES)
    left = {}     # canon -> list of neighbourhood bounds
    gcache = {}

    def add(T, gens, key):
        c = N.canon_tours([T])[0]
        if key not in gcache:
            gcache[key] = invariant_even_graph_types(gens, n)
        left.setdefault(c, []).append(gcache[key])

    for (lab, T0, g0) in forced_tournaments(n):
        add(T0, g0, ('forced', lab))
    for k in range(1, kmax + 1):
        bp = n - k
        if bp < 3 or bp not in rs:
            continue
        sizes = [bp] + [1] * k
        pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
        for bits in range(1 << (k + len(pairs))):
            inter = {}
            for i in range(k):
                inter[(0, 1 + i)] = bool((bits >> i) & 1)
            for t, (i, j) in enumerate(pairs):
                inter[(1 + i, 1 + j)] = bool((bits >> (k + t)) & 1)
            T, g = block_tournament_custom(sizes, inter)
            add(T, g, ('pat', bp, k))
    keys = sorted(left)
    adj = {}
    for i, c in enumerate(keys):
        nb = set(left[c][0])
        for x in left[c][1:]:
            nb &= x
        adj[i] = nb
    M = N.kuhn_matching(list(range(len(keys))), adj)
    if len(M) == len(keys):
        return len(keys), len(M), None
    S, Nb = N.hall_violator(list(range(len(keys))), adj, M)
    return len(keys), len(M), (len(S), len(Nb))
