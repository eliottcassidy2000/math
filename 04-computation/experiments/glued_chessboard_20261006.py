#!/usr/bin/env python3
"""Pac-man chessboards: the 8x8 board, its concentric rings and diagonal scaffolds, step versus slide
moves, and every way of gluing the edges into a closed (or partly closed) surface.

Session opus-2026-10-06-S15.  Note: 05-knowledge/results/glued_chessboard_rings_scaffolds_20261006.md

A board is the quotient of the infinite cell lattice Z^2 by a group G of lattice isometries (the edge
gluings), with walls where nothing is glued.  Cell maps are affine c -> A c + t with A a signed
permutation matrix.  Moves are computed by unfolding: step into the neighbouring plane cell, then pull it
back into the window F = {0..7}^2 with the group element that owns it (directions transform by A^-1).

Pieces: W wazir (orthogonal step), F ferz (diagonal step), K king (both); R rook, B bishop, Q queen (slide
along straight lines through seams until a wall; a closed line is reachable in full).

Reproduce:
    python3 04-computation/experiments/glued_chessboard_20261006.py            # all exact tables (~1 min)
    python3 04-computation/experiments/glued_chessboard_20261006.py --sat      # + SAT cross-checks,
                                                                                 # domination, chromatic numbers
Requires networkx, numpy; --sat needs python-sat.
"""
import sys, json, time
from collections import Counter, defaultdict, deque
from fractions import Fraction as Fr

import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism as iso

N = 8
F = [(x, y) for y in range(N) for x in range(N)]
FSET = set(F)
IDX = {c: i for i, c in enumerate(F)}
ORTH = [(1, 0), (0, 1), (-1, 0), (0, -1)]
DIAG = [(1, 1), (-1, 1), (-1, -1), (1, -1)]
ALLD = ORTH + DIAG
ORTHA, DIAGA = [(1, 0), (0, 1)], [(1, 1), (1, -1)]


# ----------------------------------------------------------------------------- affine cell maps
def mat_mul(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


def mat_vec(A, v):
    return (A[0][0] * v[0] + A[0][1] * v[1], A[1][0] * v[0] + A[1][1] * v[1])


def mat_inv(A):
    return ((A[0][0], A[1][0]), (A[0][1], A[1][1]))      # signed permutation: inverse = transpose


def det(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


I2 = ((1, 0), (0, 1))
ID = (I2, (0, 0))


def apply(g, c):
    v = mat_vec(g[0], c)
    return (v[0] + g[1][0], v[1] + g[1][1])


def compose(g, h):
    At = mat_vec(g[0], h[1])
    return (mat_mul(g[0], h[0]), (At[0] + g[1][0], At[1] + g[1][1]))


def inverse(g):
    Ai = mat_inv(g[0])
    v = mat_vec(Ai, g[1])
    return (Ai, (-v[0], -v[1]))


def point_map(g):
    """the same isometry on doubled coordinates X = 2p (corners even, cell centres odd)."""
    a1 = mat_vec(g[0], (1, 1))
    return (g[0], (2 * g[1][0] + 1 - a1[0], 2 * g[1][1] + 1 - a1[1]))


def T(a, b):
    return (I2, (a, b))


FLIP_Y = ((1, 0), (0, -1))
FLIP_X = ((-1, 0), (0, 1))
ROT90 = ((0, -1), (1, 0))
ANTI = ((0, -1), (-1, 0))
MINUS = ((-1, 0), (0, -1))
G_LR_FLIP = (FLIP_Y, (N, N - 1))        # off the right edge at row y -> in at the left edge at row 7-y
G_TB_FLIP = (FLIP_X, (N - 1, N))        # off the top at column x -> in at the bottom at column 7-x
R_A = (ROT90, (-1, 0))                   # quarter turn about the corner point (0,0): left edge <-> bottom edge
R_C = (ROT90, (2 * N - 1, 0))            # quarter turn about the corner point (8,8): right edge <-> top edge
PHI = (ANTI, (2 * N - 1, N - 1))         # diagonal glide: top edge -> right edge
PSI = (ANTI, (N - 1, -1))                # diagonal glide: left edge -> bottom edge


def half_turn_about(a, b):
    return (MINUS, (2 * a - 1, 2 * b - 1))


H_TOP, H_BOT = half_turn_about(N // 2, N), half_turn_about(N // 2, 0)
H_LEFT, H_RIGHT = half_turn_about(0, N // 2), half_turn_about(N, N // 2)

BOARDS = {
    'plane':       ([], 'no gluing (four walls)', 'disk'),
    'cylinder':    ([T(N, 0)], 'left<->right by translation; top, bottom walls', 'annulus'),
    'mobius':      ([G_LR_FLIP], 'left<->right with a flip; top, bottom walls', 'Mobius band'),
    'mobius_diag': ([PHI], 'top->right by a diagonal glide; left, bottom walls', 'Mobius band'),
    'torus':       ([T(N, 0), T(0, N)], 'both pairs by translation', 'torus (o, p1)'),
    'torus_k1':    ([T(N, 1), T(0, N)], 'helical: off the right edge one row down', 'torus (o, p1)'),
    'torus_k2':    ([T(N, 2), T(0, N)], 'off the right edge two rows down', 'torus (o, p1)'),
    'torus_k4':    ([T(N, 4), T(0, N)], 'off the right edge half a board down', 'torus (o, p1)'),
    'klein':       ([T(N, 0), G_TB_FLIP], 'left<->right translate, top<->bottom flipped (abab^-1)', 'Klein bottle (xx, pg)'),
    'klein_cc':    ([T(N, 0), (FLIP_X, (N - 2, N))], 'as klein, flip axis through cell centres', 'Klein bottle (xx, pg)'),
    'klein_diag':  ([PHI, PSI], 'adjacent edges by diagonal glides (abba)', 'Klein bottle (xx, pg)'),
    'rp2':         ([G_LR_FLIP, G_TB_FLIP], 'both pairs flipped (abab, antipodal boundary)', 'projective plane (22x, pgg)'),
    'rp2_fold':    ([G_LR_FLIP, H_TOP, H_BOT], 'left<->right flipped, top and bottom folded', 'projective plane (22x, pgg)'),
    'sphere442':   ([R_A, R_C], 'adjacent edges by quarter turns about two opposite corners', 'sphere (442, p4)'),
    'pillow':      ([H_TOP, H_BOT, H_LEFT, H_RIGHT], 'every edge folded at its midpoint', 'sphere (2222, p2)'),
    'pillow_cyl':  ([T(N, 0), H_TOP, H_BOT], 'left<->right translate, top and bottom folded', 'sphere (2222, p2)'),
}
ORDER = ['plane', 'cylinder', 'mobius', 'mobius_diag', 'torus', 'torus_k1', 'torus_k2', 'torus_k4',
         'klein', 'klein_cc', 'klein_diag', 'rp2', 'rp2_fold', 'sphere442', 'pillow', 'pillow_cyl']
DISTINCT = [n for n in ORDER if n != 'rp2_fold']     # rp2_fold is rp2 seen through a shifted window


class Board:
    def __init__(self, name, window=3 * N):
        self.name = name
        self.gens, self.desc, self.label = BOARDS[name]
        self._unfold(window)
        self.nbr = {(c, d): self.step(c, d) for c in F for d in ALLD}

    def _unfold(self, W):
        lo, hi = -W, N + W
        gens = []
        for g in self.gens:
            gens.append(g)
            if inverse(g) != g:
                gens.append(inverse(g))
        seen, q, self.elements, owner = {ID}, deque([ID]), [ID], {}

        def meets(g):
            pts = [apply(g, c) for c in [(0, 0), (N - 1, 0), (0, N - 1), (N - 1, N - 1)]]
            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            return not (max(xs) < lo or min(xs) > hi or max(ys) < lo or min(ys) > hi)
        while q:
            g = q.popleft()
            for c in F:
                p = apply(g, c)
                if lo <= p[0] <= hi and lo <= p[1] <= hi:
                    if p in owner:
                        raise ValueError(f'{self.name}: plane cell {p} covered twice (G does not act freely)')
                    owner[p] = (c, g)
            for s in gens:
                h = compose(g, s)
                if h not in seen and meets(h):
                    seen.add(h); q.append(h); self.elements.append(h)
        self.owner = owner

    def step(self, c, d):
        p = (c[0] + d[0], c[1] + d[1])
        if p in FSET:
            return p, d
        o = self.owner.get(p)
        if o is None:
            return None
        f, g = o
        return f, mat_vec(mat_inv(g[0]), d)

    def ray(self, c, d):
        out, state, start = [], (c, d), (c, d)
        while True:
            nxt = self.nbr[state]
            if nxt is None:
                break
            state = nxt
            if state == start:
                break
            out.append(state)
        return out

    def graph(self, kind):
        if kind in 'WFK':
            dirs = {'W': ORTH, 'F': DIAG, 'K': ALLD}[kind]
            adj = {c: set() for c in F}
            for c in F:
                for d in dirs:
                    t = self.nbr[(c, d)]
                    if t is not None and t[0] != c:
                        adj[c].add(t[0])
        else:
            dirs = {'R': ORTH, 'B': DIAG, 'Q': ALLD}[kind]
            adj = {}
            for c in F:
                s = set()
                for d in dirs:
                    s.update(st[0] for st in self.ray(c, d))
                s.discard(c)
                adj[c] = s
        for c in F:
            for e in adj[c]:
                assert c in adj[e]
        return adj

    def degeneracies(self, kind):
        dirs = {'W': ORTH, 'F': DIAG, 'K': ALLD}[kind]
        loops = sum(1 for c in F for d in dirs if self.nbr[(c, d)] and self.nbr[(c, d)][0] == c)
        rep = 0
        for c in F:
            ts = [self.nbr[(c, d)][0] for d in dirs if self.nbr[(c, d)] and self.nbr[(c, d)][0] != c]
            rep += len(ts) - len(set(ts))
        return loops, rep

    def colour_character(self):
        return [(g[1][0] + g[1][1]) % 2 for g in self.gens]

    def orientable(self):
        return all(det(g[0]) == 1 for g in self.gens)

    def lines(self, dirs):
        """undirected straight lines: list of (closed?, length, distinct cells, folds back?)."""
        states = [(c, d) for c in F for d in dirs]
        nxt = {s: self.nbr[s] for s in states if self.nbr[s] is not None}
        prv = {v: k for k, v in nxt.items()}
        done, orbits = set(), []
        for s in states:
            if s in done:
                continue
            cur, seen_local = s, {s}
            while cur in prv and prv[cur] not in seen_local:
                cur = prv[cur]; seen_local.add(cur)
            closed = cur in prv
            orbit, cur = [cur], cur
            done.add(cur)
            while cur in nxt:
                cur = nxt[cur]
                if cur == orbit[0]:
                    break
                orbit.append(cur); done.add(cur)
            orbits.append((closed, orbit))
        key = {st: i for i, (cl, o) in enumerate(orbits) for st in o}
        out, used = [], set()
        for i, (closed, o) in enumerate(orbits):
            if i in used:
                continue
            c0, d0 = o[0]
            r = key[(c0, (-d0[0], -d0[1]))]
            used.update((i, r))
            cells = [st[0] for st in o]
            fold = (r == i)
            length = len(cells) // 2 if (fold and closed) else len(cells)
            out.append(dict(closed=closed, length=length, distinct=len(set(cells)), fold=fold, cells=cells))
        return out

    def vertex_classes(self):
        pts = [(2 * i, 2 * j) for i in range(N + 1) for j in range(N + 1)]
        P, parent = set(pts), {p: p for p in pts}

        def find(p):
            while parent[p] != p:
                parent[p] = parent[parent[p]]; p = parent[p]
            return p
        for g in self.elements:
            A, Tt = point_map(g)
            for p in pts:
                v = mat_vec(A, p)
                q = (v[0] + Tt[0], v[1] + Tt[1])
                if q in P:
                    a, b = find(p), find(q)
                    if a != b:
                        parent[a] = b
        classes = defaultdict(list)
        for p in pts:
            classes[find(p)].append(p)
        return find, classes

    def euler(self, cells=None):
        cells = F if cells is None else cells
        find, classes = self.vertex_classes()
        allm = set()
        for (x, y) in F:
            X, Y = 2 * x, 2 * y
            allm.update([(X + 1, Y), (X + 1, Y + 2), (X, Y + 1), (X + 2, Y + 1)])
        mp_ = {m: m for m in allm}

        def mf(m):
            while mp_[m] != m:
                mp_[m] = mp_[mp_[m]]; m = mp_[m]
            return m
        for g in self.elements:
            A, Tt = point_map(g)
            for m in allm:
                v = mat_vec(A, m)
                q = (v[0] + Tt[0], v[1] + Tt[1])
                if q in mp_:
                    a, b = mf(m), mf(q)
                    if a != b:
                        mp_[a] = b
        verts, edges = set(), set()
        for (x, y) in cells:
            X, Y = 2 * x, 2 * y
            for p in [(X, Y), (X + 2, Y), (X, Y + 2), (X + 2, Y + 2)]:
                verts.add(find(p))
            for m in [(X + 1, Y), (X + 1, Y + 2), (X, Y + 1), (X + 2, Y + 1)]:
                edges.add(mf(m))
        return len(verts) - len(edges) + len(cells)

    def corners_at(self, pts):
        star, m = set(), 0
        for (X, Y) in pts:
            for (dx, dy) in [(1, 1), (-1, 1), (-1, -1), (1, -1)]:
                c = ((X + dx - 1) // 2, (Y + dy - 1) // 2)
                if c in FSET:
                    star.add(c); m += 1
        return m, star


def rings_plane():
    R = defaultdict(list)
    for (x, y) in F:
        k = (max(abs(2 * x + 1 - N), abs(2 * y + 1 - N)) + 1) // 2
        R[k].append((x, y))
    return dict(R)


def nxg(adj):
    G = nx.Graph(); G.add_nodes_from(F)
    for c in F:
        for e in adj[c]:
            G.add_edge(c, e)
    return G


# ----------------------------------------------------------------------------- analyses
def indep_poly(adj):
    """exact independence polynomial by a blocked-future-set DP (row-major order)."""
    later = []
    for i, v in enumerate(F):
        m = 0
        for e in adj[v]:
            if IDX[e] > i:
                m |= 1 << IDX[e]
        later.append(m)

    def padd(a, b):
        if len(a) < len(b):
            a, b = b, a
        return tuple(a[i] + (b[i] if i < len(b) else 0) for i in range(len(a)))
    states = {0: (1,)}
    for i in range(len(F)):
        bit, new = 1 << i, {}
        for B, p in states.items():
            B2 = B & ~bit
            new[B2] = padd(new[B2], p) if B2 in new else p
            if not B & bit:
                B3 = (B | later[i]) & ~bit
                sp = (0,) + p
                new[B3] = padd(new[B3], sp) if B3 in new else sp
        states = new
    tot = (0,)
    for p in states.values():
        tot = padd(tot, p)
    return list(tot)


def isometries(b):
    """automorphisms of the square-tiled surface via developing maps from cell (0,0)."""
    D4 = [((1, 0), (0, 1)), ((0, -1), (1, 0)), ((-1, 0), (0, -1)), ((0, 1), (-1, 0)),
          ((1, 0), (0, -1)), ((-1, 0), (0, 1)), ((0, 1), (1, 0)), ((0, -1), (-1, 0))]

    def lin(c, d):
        p = (c[0] + d[0], c[1] + d[1])
        if p in FSET or b.owner.get(p) is None:
            return I2
        return mat_inv(b.owner[p][1][0])
    found = []
    for target in F:
        for A in D4:
            phi, stack, ok = {(0, 0): (target, A)}, [(0, 0)], True
            while stack and ok:
                a = stack.pop()
                fa, La = phi[a]
                for d in ORTH:
                    s, t = b.nbr[(a, d)], b.nbr[(fa, mat_vec(La, d))]
                    if (s is None) != (t is None):
                        ok = False; break
                    if s is None:
                        continue
                    L2 = mat_mul(lin(fa, mat_vec(La, d)), mat_mul(La, mat_inv(lin(a, d))))
                    if mat_vec(L2, s[1]) != t[1]:
                        ok = False; break
                    if s[0] in phi:
                        if phi[s[0]] != (t[0], L2):
                            ok = False; break
                    else:
                        phi[s[0]] = (t[0], L2); stack.append(s[0])
            if not ok or len(phi) != len(F) or len({v[0] for v in phi.values()}) != len(F):
                continue
            good = True
            for a in F:
                fa, La = phi[a]
                for d in DIAG:
                    s, t = b.nbr[(a, d)], b.nbr[(fa, mat_vec(La, d))]
                    if (s is None) != (t is None) or (s is not None and phi[s[0]][0] != t[0]):
                        good = False; break
                if not good:
                    break
            if good:
                found.append(phi)
    return found


def lattice_quotient(L):
    e1, e2 = L
    D = e1[0] * e2[1] - e1[1] * e2[0]

    def key(v):
        return (Fr(v[0] * e2[1] - v[1] * e2[0], D) % 1, Fr(e1[0] * v[1] - e1[1] * v[0], D) % 1)
    cells, frontier = {key((0, 0)): (0, 0)}, [(0, 0)]
    while frontier:
        nf = []
        for v in frontier:
            for d in ORTH:
                w = (v[0] + d[0], v[1] + d[1])
                if key(w) not in cells:
                    cells[key(w)] = w; nf.append(w)
        frontier = nf
    assert len(cells) == abs(D)
    return cells, key


def lattice_graphs(L):
    cells, key = lattice_quotient(L)
    W = nx.Graph(); W.add_nodes_from(cells)
    for k, v in cells.items():
        for d in ORTH:
            k2 = key((v[0] + d[0], v[1] + d[1]))
            if k2 != k:
                W.add_edge(k, k2)
    return W


def rotate(v):
    assert (v[0] + v[1]) % 2 == 0
    return ((v[0] + v[1]) // 2, (v[0] - v[1]) // 2)


def line_multigraph(b, axis):
    dirs = []
    for d in axis:
        dirs += [d, (-d[0], -d[1])]
    states = [(c, d) for c in F for d in dirs]
    nxt = {s: b.nbr[s] for s in states if b.nbr[s] is not None}
    prv = {v: k for k, v in nxt.items()}
    lid, k = {}, 0
    for s in states:
        if s in lid:
            continue
        orbit, stack = {s}, [s]
        while stack:
            u = stack.pop()
            for v in (nxt.get(u), prv.get(u), (u[0], (-u[1][0], -u[1][1]))):
                if v is not None and v not in orbit:
                    orbit.add(v); stack.append(v)
        for u in orbit:
            lid[u] = k
        k += 1
    M = nx.MultiGraph(); M.add_nodes_from(range(k))
    for c in F:
        M.add_edge(lid[(c, axis[0])], lid[(c, axis[1])])
    return M


def lattice_multigraph(L):
    cells, key = lattice_quotient(L)
    lid, k = {}, 0
    for d in [(1, 0), (0, 1)]:
        for kk, v in cells.items():
            if (kk, d) in lid:
                continue
            cur = v
            while (key(cur), d) not in lid:
                lid[(key(cur), d)] = k
                cur = (cur[0] + d[0], cur[1] + d[1])
            k += 1
    M = nx.MultiGraph(); M.add_nodes_from(range(k))
    for kk in cells:
        M.add_edge(lid[(kk, (1, 0))], lid[(kk, (0, 1))])
    return M


def multi_iso(M1, M2):
    def w(M):
        G = nx.Graph(); G.add_nodes_from(M.nodes())
        for u, v in M.edges():
            if G.has_edge(u, v):
                G[u][v]['m'] += 1
            else:
                G.add_edge(u, v, m=1)
        return G
    G1, G2 = w(M1), w(M2)
    if G1.number_of_nodes() != G2.number_of_nodes():
        return False
    return iso.GraphMatcher(G1, G2, edge_match=iso.numerical_edge_match('m', 1)).is_isomorphic()


def mdesc(M):
    loops = sum(1 for u, v in M.edges() if u == v)
    return f'{M.number_of_nodes()} lines, {loops} loops, line lengths {sorted((d for _, d in M.degree()), reverse=True)}'


# ----------------------------------------------------------------------------- SAT parts (--sat)
def line_domination(b, axis):
    """gamma of a rook-type slider when every component's line multigraph is complete bipartite between
    two line families (each pair of lines from the two families meets): gamma = sum of min(|A|,|B|).
    Proof: one piece per line of the smaller family dominates; with fewer, an empty A-line and an empty
    B-line meet in an undominated square.  Returns None when the hypothesis fails."""
    M = line_multigraph(b, axis)
    total = 0
    for comp in nx.connected_components(M):
        H = nx.Graph(M.subgraph(comp))
        if any(u == v for u, v in M.subgraph(comp).edges()) or not nx.is_bipartite(H):
            return None
        A, B = nx.bipartite.sets(H)
        if H.number_of_edges() != len(A) * len(B):
            return None
        total += min(len(A), len(B))
    return total


def sat_parts(name, b, dp, budget=20_000_000):
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Solver
    out = {}
    var = {c: i + 1 for i, c in enumerate(F)}
    for kind in 'KRBQ':
        adj = b.graph(kind)
        # alpha cross-check: no independent set of size alpha+1; count size-alpha sets when <= 20000
        a, cnt = dp[kind]['alpha'], dp[kind]['count']
        cl = [[-var[c], -var[e]] for c in F for e in adj[c] if IDX[e] > IDX[c]]
        over = CardEnc.atleast(list(var.values()), bound=a + 1, top_id=len(F), encoding=EncType.seqcounter)
        with Solver(name='glucose4', bootstrap_with=cl + over.clauses) as s:
            assert not s.solve(), (name, kind, 'alpha too small')
        if cnt <= 20000:
            at = CardEnc.atleast(list(var.values()), bound=a, top_id=len(F), encoding=EncType.seqcounter)
            n = 0
            with Solver(name='glucose4', bootstrap_with=cl + at.clauses) as s:
                while s.solve():
                    m = [l for l in s.get_model() if 0 < l <= len(F)]
                    n += 1; s.add_clause([-l for l in m])
            assert n == cnt, (name, kind, n, cnt)
        # domination number: the line formula when it applies (proved), else SAT
        g = line_domination(b, ORTHA if kind == 'R' else DIAGA) if kind in 'RB' else None
        if g is not None:
            out['gamma_' + kind] = f'{g} (line formula)'
        else:
            base = [[var[v]] + [var[u] for u in adj[v]] for v in F]
            k = 1
            while True:
                card = CardEnc.atmost(list(var.values()), bound=k, top_id=len(F), encoding=EncType.seqcounter)
                with Solver(name='glucose4', bootstrap_with=base + card.clauses) as s:
                    if s.solve():
                        break
                k += 1
            out['gamma_' + kind] = k
    # chromatic numbers, starting at the proven lower bound max(omega, ceil(64/alpha))
    for kind in 'WFKRBQ':
        adj = b.graph(kind)
        G = nxg(adj)
        clique = max(nx.find_cliques(G), key=len)
        alpha = len(indep_poly(adj)) - 1
        lb = max(len(clique), -(-len(F) // alpha))
        trail = []
        for k in range(lb, len(F) + 1):
            def v_(v, c):
                return IDX[v] * k + c + 1
            cl = [[v_(v, c) for c in range(k)] for v in F]
            cl += [[-v_(v, c), -v_(e, c)] for v in F for e in adj[v] if IDX[e] > IDX[v] for c in range(k)]
            cl += [[v_(v, i)] for i, v in enumerate(clique)]
            with Solver(name='cadical153', bootstrap_with=cl) as s:
                s.conf_budget(budget)
                r = s.solve_limited()
            trail.append({True: 'sat', False: 'unsat', None: 'unknown'}[r])
            if r is True:
                break
        exact = all(t == 'unsat' for t in trail[:-1])
        out['chi_' + kind] = f'{k}' if exact else f'<= {k} (>= {lb + trail.index("unknown")})'
    return out


# ----------------------------------------------------------------------------- main report
def main():
    t0 = time.time()
    sat = '--sat' in sys.argv
    boards = {n: Board(n) for n in ORDER}
    R = rings_plane()
    ringof = {c: k for k in R for c in R[k]}
    J = {}
    print('PAC-MAN CHESSBOARDS (opus-2026-10-06-S15)\n')

    print('== 0. The owner\'s structures on the plane board')
    b = boards['plane']
    W = b.graph('W')
    print('   rings (Chebyshev radius k-1/2 about the centre point):', {k: len(v) for k, v in R.items()})
    print('   ring cycles: internal W-edges', {k: sum(1 for c in R[k] for e in W[c] if ringof[e] == k) // 2 for k in R},
          ' spokes', dict(Counter((ringof[c], ringof[e]) for c in F for e in W[c] if ringof[e] > ringof[c])))
    for col in (0, 1):
        dl = sorted(Counter(x - y for (x, y) in F if (x + y) % 2 == col).items())
        al = sorted(Counter(x + y for (x, y) in F if (x + y) % 2 == col).items())
        print(f'   colour {col}: diagonal (x-y) lengths {[v for k, v in dl]}  anti-diagonal (x+y) lengths {[v for k, v in al]}')
    Bm = {c: len(v) for c, v in b.graph('B').items()}
    Qm = {c: len(v) for c, v in b.graph('Q').items()}
    print('   bishop mobility by ring:', {k: sorted(set(Bm[c] for c in R[k])) for k in R},
          ' queen:', {k: sorted(set(Qm[c] for c in R[k])) for k in R})
    assert all(Bm[c] == 15 - 2 * ringof[c] and Qm[c] == 29 - 2 * ringof[c] for c in F)

    print('\n== 1. Boards: surface, Euler characteristic, orientability, colour character, isometries')
    J['boards'] = {}
    for n in ORDER:
        b = boards[n]
        isos = isometries(b)
        parent = {c: c for c in F}

        def fd(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for phi in isos:
            for a, (fa, L) in phi.items():
                ra, rb = fd(a), fd(fa)
                if ra != rb:
                    parent[ra] = rb
        orbits = len(set(fd(c) for c in F))
        find, classes = b.vertex_classes()
        cones = Counter()
        for r, pts in classes.items():
            m, star = b.corners_at(pts)
            cones[m] += 1
        J['boards'][n] = dict(label=b.label, chi=b.euler(), orientable=b.orientable(), colour=b.colour_character(),
                              isom=len(isos), cell_orbits=orbits)
        print(f'   {n:12s} {b.label:28s} chi={b.euler():2d} orientable={str(b.orientable()):5s} colour-char={b.colour_character()} '
              f'|Isom|={len(isos):3d} cell-orbits={orbits:2d}  vertex classes by #corners {dict(sorted(cones.items()))}')

    print('\n== 2. Seam localisation: rings 1-3 untouched; ring 4 carries the topology')
    for n in ORDER:
        b = boards[n]
        W = b.graph('W')
        internal = {k: sum(1 for c in R[k] for e in W[c] if ringof[e] == k) // 2 for k in R}
        sp = dict(Counter((ringof[c], ringof[e]) for c in F for e in W[c] if ringof[e] > ringof[c]))
        print(f'   {n:12s} ring-internal W-edges {internal}  spokes {sp}  chi(ring 4 region)={b.euler(R[4]):2d}  '
              f'chi(inner 6x6)={b.euler(R[1] + R[2] + R[3])}')

    print('\n== 3. King shells about every vertex class (cone points show their angle)')
    for n in ORDER:
        b = boards[n]
        K = b.graph('K')
        find, classes = b.vertex_classes()
        seqs = Counter()
        special = []
        for r, pts in classes.items():
            m, star = b.corners_at(pts)
            dist, fr = {c: 0 for c in star}, list(star)
            while fr:
                nf = []
                for c in fr:
                    for e in K[c]:
                        if e not in dist:
                            dist[e] = dist[c] + 1; nf.append(e)
                fr = nf
            sh = Counter(dist.values())
            seq = tuple(sh[i] for i in range(max(sh) + 1))
            if (N, N) in pts:
                centre = seq
            interior = all(0 < X < 2 * N and 0 < Y < 2 * N for (X, Y) in pts) or n not in ('plane', 'cylinder', 'mobius', 'mobius_diag')
            if m != 4 and interior:
                special.append((m, [(X // 2, Y // 2) for (X, Y) in pts], seq))
        print(f'   {n:12s} centre point: {centre}' + ''.join(f'\n{"":16s}cone point, {m} corner(s) at {p}: {s}' for m, p, s in special))

    print('\n== 4. Scaffold fusion: W bipartite <=> F has two components <=> colour character trivial')
    for n in ORDER:
        b = boards[n]
        GW, GF = nxg(b.graph('W')), nxg(b.graph('F'))
        fcomp = nx.number_connected_components(GF)
        triv = all(x == 0 for x in b.colour_character())
        assert nx.is_bipartite(GW) == (fcomp == 2) == triv
        Fbip = [nx.is_bipartite(GF.subgraph(c)) for c in nx.connected_components(GF)]
        print(f'   {n:12s} colour character {b.colour_character()}  W bipartite {nx.is_bipartite(GW)!s:5s}  '
              f'F components {fcomp}  (each F component bipartite: {Fbip})')

    print('\n== 5. Diagonal worlds are rotated colour covers')
    covers = {'torus': [(8, 0), (0, 8)], 'torus_k2': [(8, 2), (0, 8)], 'torus_k4': [(8, 4), (0, 8)],
              'torus_k1': [(16, 2), (0, 8)], 'klein': [(8, 0), (0, 16)]}
    for n, Lc in covers.items():
        b = boards[n]
        Lr = [rotate(v) for v in Lc]
        GF = nxg(b.graph('F'))
        Wl = lattice_graphs(Lr)
        okF = [nx.vf2pp_is_isomorphic(GF.subgraph(c).copy(), Wl) for c in nx.connected_components(GF)]
        Ml = lattice_multigraph(Lr)
        MB = line_multigraph(b, DIAGA)
        okB = [multi_iso(MB.subgraph(c).copy(), Ml) for c in nx.connected_components(MB)]
        print(f'   {n:9s} colour cover Z^2/<{Lc}> -> rotated <{Lr}>:  F ~ W(cover) {okF}   bishop lines ~ rook lines(cover) {okB}')
    MB = line_multigraph(boards['klein_cc'], DIAGA)
    Ml = lattice_multigraph([rotate(v) for v in covers['klein']])
    print('   hostile control: klein_cc (colour preserved) vs the klein cover:',
          [multi_iso(MB.subgraph(c).copy(), Ml) for c in nx.connected_components(MB)])
    print('   F(rp2) ~ F(sphere442):', nx.vf2pp_is_isomorphic(nxg(boards['rp2'].graph('F')), nxg(boards['sphere442'].graph('F'))))

    print('\n== 6. Straight lines (sliding): census (closed/open, length, distinct cells, folds back)')
    J['lines'] = {}
    for n in ORDER:
        b = boards[n]
        res = {}
        for lab, dirs in [('orth', ORTH), ('diag', DIAG)]:
            cnt = Counter((('closed' if l['closed'] else 'open'), l['length'], l['distinct'], 'fold' if l['fold'] else '')
                          for l in b.lines(dirs))
            res[lab] = {f'{k[0]} {k[1]}' + (f'({k[2]} cells)' if k[2] != k[1] else '') + (' fold' if k[3] else ''): v
                        for k, v in sorted(cnt.items(), key=lambda kv: (kv[0][0], kv[0][1]))}
        J['lines'][n] = res
        print(f'   {n:12s} orth {res["orth"]}\n{"":16s}diag {res["diag"]}')

    print('\n== 7. Line multigraphs (vertices = lines, edges = squares); slider graph = its line graph')
    M = line_multigraph(boards['plane'], DIAGA)
    for C in nx.connected_components(M):
        print('   plane bishop scaffold component:', mdesc(M.subgraph(C)))
    for kind, ax in [('R', ORTHA), ('B', DIAGA)]:
        Ms = {n: line_multigraph(boards[n], ax) for n in DISTINCT}
        classes = []
        for n in DISTINCT:
            for cl in classes:
                if multi_iso(Ms[n], Ms[cl[0]]):
                    cl.append(n); break
            else:
                classes.append([n])
        print(f'   {kind}: isomorphism classes of line multigraphs: {classes}')
    for kind in 'WFKRBQ':
        E = {n: frozenset(frozenset((c, e)) for c, s in boards[n].graph(kind).items() for e in s) for n in ORDER}
        groups = defaultdict(list)
        for n in ORDER:
            groups[E[n]].append(n)
        print(f'   {kind}: boards with literally identical attack sets: {[g for g in groups.values() if len(g) > 1]}')

    print('\n== 8. Basic move-graph data (edges, degree histogram, components, diameter)')
    J['graphs'] = {}
    for n in ORDER:
        b = boards[n]
        row = {}
        for kind in 'WFKRBQ':
            G = nxg(b.graph(kind))
            comps = [len(c) for c in nx.connected_components(G)]
            diam = [nx.diameter(G.subgraph(c)) for c in nx.connected_components(G)]
            row[kind] = dict(E=G.number_of_edges(), deg=dict(sorted(Counter(d for _, d in G.degree()).items())),
                             comps=comps, diam=sorted(set(diam)))
            if kind in 'WK' and len(comps) == 1:
                row[kind]['mean_dist'] = round(nx.average_shortest_path_length(G), 4)
        lw, rw = b.degeneracies('W'); lf, rf = b.degeneracies('F')
        J['graphs'][n] = row
        print(f'   {n:12s} ' + '  '.join(f'{k}:E{v["E"]} deg{v["deg"]} c{v["comps"]} d{v["diam"]}' + (f' m{v["mean_dist"]}' if 'mean_dist' in v else '')
                                         for k, v in row.items()) + f'  [W loops/repeats {lw}/{rw}, F loops {lf}]')

    print('\n== 9. Non-attacking pieces: maximum and number of maximum placements (exact DP)')
    J['alpha'] = {}
    for n in ORDER:
        b = boards[n]
        row = {}
        for kind in 'KRBQ':
            p = indep_poly(b.graph(kind))
            row[kind] = dict(alpha=len(p) - 1, count=p[-1])
        J['alpha'][n] = row
        print(f'   {n:12s} ' + '  '.join(f'{k}: {v["alpha"]:2d} in {v["count"]:7d} ways' for k, v in row.items()))

    print('\n== 10. Mobility maps (row 7 on top) where they are not constant')
    for n in ['plane', 'rp2', 'rp2_fold', 'sphere442', 'pillow', 'klein_cc']:
        b = boards[n]
        for kind in 'RBQ':
            m = {c: len(s) for c, s in b.graph(kind).items()}
            if len(set(m.values())) > 1:
                print(f'   {n} {kind}:')
                for y in range(N - 1, -1, -1):
                    print('      ' + ' '.join(f'{m[(x, y)]:2d}' for x in range(N)))

    if sat:
        print('\n== 11. SAT: alpha cross-check (no larger set; counts re-enumerated when <= 20000), domination and chromatic numbers')
        J['sat'] = {}
        for n in DISTINCT:
            r = sat_parts(n, boards[n], J['alpha'][n])
            J['sat'][n] = r
            print(f'   {n:12s} alpha confirmed; gamma K,R,B,Q = {r["gamma_K"]}, {r["gamma_R"]}, {r["gamma_B"]}, {r["gamma_Q"]};'
                  f'  chi W,F,K,R,B,Q = {r["chi_W"]}, {r["chi_F"]}, {r["chi_K"]}, {r["chi_R"]}, {r["chi_B"]}, {r["chi_Q"]}', flush=True)

    with open(sys.argv[sys.argv.index('--json') + 1] if '--json' in sys.argv else 'glued_chessboard_20261006.json', 'w') as f:
        json.dump(J, f, indent=1, default=str)
    print(f'\n[done in {time.time() - t0:.1f}s]')


if __name__ == '__main__':
    main()
