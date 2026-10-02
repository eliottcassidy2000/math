"""Camion, Busch and the two tournament gaps; tournaments as polyhedra; the Collatz comparison.

Companion to 05-knowledge/results/camion_busch_gaps_polyhedra_collatz_20261002.md.
The order-7 census is the C program camion_busch_strong7_20261002.c (+ .out).

Sections:
  A  forbidden independence polynomials: H = I(Omega, 2) equals 7 or 21 exactly for which conflict graphs
  B  the Camion-Moon bound B(m) against the strong floor f(m) (exhaustive m <= 6; m = 7 in the C census)
  C  the mod-4 locks around the gaps: n = 5 (gap 7) and strong n = 6 (gap 21)
  D  score vectors as weights: prod (x_i + x_j) = s_delta; the n = 4 shells (truncated octahedron, cube, octahedron),
     the arc-reversal roots (cuboctahedron), the score lattice (FCC, Voronoi cell the rhombic dodecahedron); n = 5
  E  planar duality on the tetrahedron: cyclic triangles of T = sources + sinks of the dual tournament T*
  F  P7 on the torus: the face-cyclic orientation of the 7-vertex torus triangulation; Aut(P7); dual (Heawood,
     Fano incidence orientation); medial graph; H(P7)
  G  tournaments with T - v transitive: H = 1 + 2 K(w) for a binary word w
  H  60 = p^2 q r = |A5|: the ten proper divisors as icosahedral stabilizer orders and orbit sizes
  I  the Collatz side: Syracuse data for 7 and 21, and the numerology of {1, 2, 4}
"""
import itertools
import math
from collections import Counter, defaultdict

import sympy as sp


def tour(n, mask):
    b = {}
    for t, (x, y) in enumerate(itertools.combinations(range(n), 2)):
        b[(x, y)] = bool(mask >> t & 1); b[(y, x)] = not b[(x, y)]
    return b


def hp(n, beats, verts=None):
    vs = list(range(n)) if verts is None else list(verts); k = len(vs)
    dp = [[0] * k for _ in range(1 << k)]
    for v in range(k): dp[1 << v][v] = 1
    for S in range(1 << k):
        for v in range(k):
            if dp[S][v]:
                for u in range(k):
                    if not S >> u & 1 and beats[(vs[v], vs[u])]:
                        dp[S | 1 << u][u] += dp[S][v]
    return sum(dp[-1])


def hc(beats, verts):
    vs = list(verts); k = len(vs); f = vs[0]
    return sum(1 for r in itertools.permutations(vs[1:])
               if all(beats[(c, d)] for c, d in zip((f,) + r, r + (f,))))


def strong(n, beats):
    for fwd in (True, False):
        seen, st = {0}, [0]
        while st:
            u = st.pop()
            for w in range(n):
                if w not in seen and (beats[(u, w)] if fwd else beats[(w, u)]):
                    seen.add(w); st.append(w)
        if len(seen) < n:
            return False
    return True


def cyclic3(beats, t):
    a, b, c = t
    return (beats[(a, b)] and beats[(b, c)] and beats[(c, a)]) or (beats[(a, c)] and beats[(c, b)] and beats[(b, a)])


# ============================================================================ A
COMPLEMENT_NAMES = {(): 'complete', (1, 1): 'K - e', (1, 1, 1, 1): 'K - 2K2', (1, 1, 2): 'K - P3',
                    (1, 1, 2, 2): 'P4 (self-complementary)', (1, 1, 1, 3): 'K3 + K1 (complement K_{1,3})'}


def section_A():
    print('== A. H = I(Omega, 2): which conflict graphs give 7 and 21')
    x = sp.symbols('x')
    for target in (7, 21):
        K = (target - 1) // 2      # alpha_1 + 2 alpha_2 + 4 alpha_3 + ... = K
        out = []
        for a1 in range(K + 1):
            for a2 in range((K - a1) // 2 + 1):
                for a3 in range((K - a1 - 2 * a2) // 4 + 1):
                    if a1 + 2 * a2 + 4 * a3 != K:
                        continue
                    # graphs on a1 vertices: the complement has a2 edges, a3 triangles and no K4
                    V = range(a1); found = set()
                    for E in itertools.combinations(list(itertools.combinations(V, 2)), a2):
                        tri = sum(1 for t in itertools.combinations(V, 3) if all(p in E for p in itertools.combinations(t, 2)))
                        k4 = any(all(p in E for p in itertools.combinations(t, 2)) for t in itertools.combinations(V, 4))
                        if tri == a3 and not k4:
                            found.add(tuple(sorted(Counter(v for e in E for v in e).values())))
                    out.append(((a1, a2, a3), sorted(found)))
        print(f'  I(Omega, 2) = {target}:')
        for a, found in out:
            if not found:
                continue
            poly = sp.expand(1 + a[0] * x + a[1] * x**2 + a[2] * x**3)
            names = [f'K{a[0]}' + COMPLEMENT_NAMES[f][1:] if COMPLEMENT_NAMES[f].startswith('K -') else
                     (f'K{a[0]}' if f == () else COMPLEMENT_NAMES[f]) for f in found]
            print(f'    alpha = {a}: {names};  I(Omega, x) = {sp.factor(poly)}')
        print(f'    (alpha-vectors with no graph: {[a for a, found in out if not found]})')
    assert sp.factor(1 + 4 * x + 3 * x**2) == (x + 1) * (3 * x + 1)
    print('  so H = 7 iff Omega = K3, and H = 21 iff Omega is K10, K8 - e, K6 - 2K2, K6 - P3, P4 or K3 + K1.')
    print('  I(P4, x) = I(K3 + K1, x) = (1 + x)(1 + 3x): the 7-polynomial 1 + 3x is a factor of one of the four')
    print('  21-polynomials, so the 21-problem contains the 7-problem')


# ============================================================================ B, C
def strong_data(n):
    rows = []
    for mask in range(1 << math.comb(n, 2)):
        b = tour(n, mask)
        st = strong(n, b)
        c3 = [t for t in itertools.combinations(range(n), 3) if cyclic3(b, t)]
        c5 = sum(hc(b, S) for S in itertools.combinations(range(n), 5)) if n >= 5 else 0
        d33 = sum(1 for A, B in itertools.combinations(c3, 2) if not set(A) & set(B))
        H = hp(n, b)
        assert H == 1 + 2 * (len(c3) + c5) + 4 * d33
        rows.append((st, len(c3), c5, d33, H))
    return rows


def camion_moon_bound(m):
    # Moon: c3 >= m - 2; vertex-pancyclicity: every vertex on a k-cycle, so c_k >= ceil(m / k); OCF with alpha_2 >= 0
    return 1 + 2 * ((m - 2) + sum(math.ceil(m / k) for k in range(5, m + 1, 2)))


def section_BC():
    print('== B. the Camion-Moon bound against the strong floor')
    floor = {}
    data = {}
    for n in range(3, 7):
        rows = strong_data(n)
        data[n] = rows
        sH = [H for st, *_, H in rows if st]
        floor[n] = min(sH)
        assert all(H >= n for st, *_, H in rows if st)      # Camion: a Hamiltonian cycle gives n Hamiltonian paths
        print(f'  m = {n}: strong labelled {len(sH)}, floor f(m) = {floor[n]}, Camion-Moon bound B(m) = '
              f'{camion_moon_bound(n)}, strong H values {sorted(set(sH))}')
    print('  m = 7 (C census, camion_busch_strong7_20261002.out): f(7) = 25, min alpha_1 = 9; B(7) =',
          camion_moon_bound(7), '; B(8) =', camion_moon_bound(8), '; B(9) =', camion_moon_bound(9))
    assert [camion_moon_bound(m) for m in range(3, 10)] == [3, 5, 9, 13, 17, 21, 25]
    assert floor == {3: 3, 4: 5, 5: 9, 6: 15}
    for n in (5, 6):
        at_floor = {(c3 + c5, d33) for st, c3, c5, d33, H in data[n] if st and H == floor[n]}
        print(f'  m = {n}: the floor is attained only with (alpha_1, alpha_2) in {sorted(at_floor)}')
        assert all(a2 == 0 for _, a2 in at_floor)
    print('  (m = 7: only (12, 0), one isomorphism class of 5040 labellings; see the C census)')
    print('  -> gap 7 lies in (f(4), f(5)) = (5, 9) and f(5) = 9 = B(5): Camion + Moon certify it exactly.')
    print('     gap 21 lies in (f(6), f(7)) = (15, 25), but B(7) = 17 and B(8) = 21 do not exceed 21: Camion + Moon')
    print('     cannot certify it; it needs the true floor f(7) = 25 (Busch 2006) and floor monotonicity (THM-1370).')

    print('== C. the mod-4 locks around the gaps (THM-466: H = 1 + 2 alpha_1 mod 4)')
    near7 = Counter((c3, c5, H) for st, c3, c5, d33, H in data[5] if 5 <= H <= 9)
    print('  n = 5, all tournaments with 5 <= H <= 9 (c3, c5, H):', dict(near7))
    assert set(H % 4 for (_, _, H) in near7) == {1}
    by_c3 = defaultdict(set); c5par = defaultdict(set); d33s = defaultdict(set)
    for st, c3, c5, d33, H in data[6]:
        if st:
            by_c3[c3].add(H); c5par[c3].add(c5); d33s[c3].add(d33)
    for c3 in sorted(by_c3):
        print(f'  strong n = 6, c3 = {c3}: c5 in {sorted(c5par[c3])}, d33 in {sorted(d33s[c3])}, H in {sorted(by_c3[c3])}')
    assert max(by_c3[4]) == 17 and min(by_c3[6]) == 23 and c5par[5] == {4, 6} and d33s[5] == {0, 1}
    assert {H % 4 for H in by_c3[5]} == {3} and by_c3[5] == {19, 23, 27}
    print('  -> the window around 7 is locked to H = 1 (mod 4) [7 = 3 mod 4]; the window around 21 is c3 = 5, where')
    print('     c5 is always even and d33 <= 1, so H = 11 + 2 c5 + 4 d33 is 3 (mod 4) [21 = 1 mod 4].')
    print('     c3 = 4 gives H <= 17 and c3 >= 6 gives H >= 23, so 21 falls through.')
    return data


# ============================================================================ D
def section_D():
    print('== D. score vectors are the weights of V_rho; the n = 4 shells')
    for n in (3, 4, 5):
        xs = sp.symbols(f'x0:{n}')
        lhs = sp.prod([xs[i] + xs[j] for i, j in itertools.combinations(range(n), 2)])
        num = sp.Matrix(n, n, lambda i, j: xs[i] ** (2 * (n - 1 - j))).det()
        den = sp.Matrix(n, n, lambda i, j: xs[i] ** (n - 1 - j)).det()
        assert sp.expand(sp.cancel(num / den) - lhs) == 0
        gen = Counter(tuple(sum(tour(n, m)[(v, w)] for w in range(n) if w != v) for v in range(n))
                      for m in range(1 << math.comb(n, 2)))
        poly = sp.Poly(sp.expand(lhs), *xs)
        assert all(poly.coeff_monomial(sp.prod([xs[i] ** e for i, e in enumerate(s)])) == c for s, c in gen.items())
        print(f'  n = {n}: prod_(i<j) (x_i + x_j) = s_delta (bialternant) = sum over tournaments of x^score: verified; '
              f'{len(gen)} score vectors')
    n = 4
    gen = Counter(tuple(sum(tour(n, m)[(v, w)] for w in range(n) if w != v) for v in range(n)) for m in range(64))
    c = 1.5
    shells = defaultdict(list)
    for s, mult in gen.items():
        d2 = sum((si - c) ** 2 for si in s)
        shells[d2].append(mult)
    for d2 in sorted(shells):
        print(f'  n = 4, squared distance {d2} from the centre: {len(shells[d2])} score vectors, multiplicity '
              f'{sorted(set(shells[d2]))}')
    assert {d2: (len(v), set(v)) for d2, v in shells.items()} == {1.0: (6, {4}), 3.0: (8, {2}), 5.0: (24, {1})}
    # geometry: project to the plane sum = 6 with an orthonormal basis; identify the three shells
    import numpy as np
    B = np.linalg.qr(np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float).T)[0]
    pts = {d2: np.array([np.array(s, float) - c for s, m in gen.items() if abs(sum((si - c) ** 2 for si in s) - d2) < 1e-9]) @ B
           for d2 in (1.0, 3.0, 5.0)}
    def edges_at_min(P):
        D = np.linalg.norm(P[:, None] - P[None], axis=2); m = D[D > 1e-9].min()
        return int(np.sum(np.abs(D - m) < 1e-9) // 2), round(float(m), 6)
    print('  shell 1 (6 strong score vectors): octahedron, edges', edges_at_min(pts[1.0]),
          '; shell 3 (8 "3-cycle + source/sink"): cube, edges', edges_at_min(pts[3.0]),
          '; shell 5 (24 transitive): truncated octahedron, edges', edges_at_min(pts[5.0]))
    assert edges_at_min(pts[1.0])[0] == 12 and edges_at_min(pts[3.0])[0] == 12 and edges_at_min(pts[5.0])[0] == 36
    # octahedron points sit at the centres of the cube's faces
    cube, octa = pts[3.0], pts[1.0]
    for o in octa:
        near = sorted(np.linalg.norm(cube - o, axis=1))[:4]
        assert np.allclose(near, near[0])
    print('  the 6 strong score vectors are the face centres of the cube of the 8 one-triangle score vectors')
    roots = np.array([np.eye(4)[i] - np.eye(4)[j] for i in range(4) for j in range(4) if i != j]) @ B
    print('  arc reversals change the score by a root e_j - e_i: 12 roots, edges at minimal distance',
          edges_at_min(roots), '= cuboctahedron (12 vertices, 24 edges), the medial graph of cube and octahedron')
    assert edges_at_min(roots)[0] == 24
    # the score lattice is A3 = FCC: the centre (1.5,...) is an octahedral hole (6 nearest lattice points at equal
    # distance); the 8 tetrahedral holes around a lattice point and its 6 octahedral holes are the 14 vertices of
    # the Voronoi cell, a rhombic dodecahedron
    lat = [np.array(v, float) for v in itertools.product(range(-3, 4), repeat=4) if sum(v) == 0]
    L = np.array(lat) @ B
    hole_o = (np.array([0.5, 0.5, -0.5, -0.5])) @ B          # an octahedral hole
    hole_t = (np.array([0.75, -0.25, -0.25, -0.25])) @ B     # a tetrahedral hole
    def neigh(h):
        D = np.linalg.norm(L - h, axis=1); m = D.min(); return int(np.sum(np.abs(D - m) < 1e-9))
    print('  score lattice A3 = FCC: an octahedral hole has', neigh(hole_o), 'nearest lattice points, a tetrahedral hole',
          neigh(hole_t), '; Voronoi cell vertices: 6 octahedral + 8 tetrahedral holes = rhombic dodecahedron')
    assert neigh(hole_o) == 6 and neigh(hole_t) == 4
    n = 5
    gen5 = Counter(tuple(sum(tour(n, m)[(v, w)] for w in range(n) if w != v) for v in range(n)) for m in range(1 << 10))
    sh5 = defaultdict(Counter)
    for s, mult in gen5.items():
        sh5[sum((si - 2) ** 2 for si in s)][mult] += 1
    print('  n = 5 shells around the centre (2,2,2,2,2): squared distance -> {multiplicity: number of score vectors}:',
          {d: dict(v) for d, v in sorted(sh5.items())})
    assert gen5[(2, 2, 2, 2, 2)] == 24


# ============================================================================ E
def section_E():
    print('== E. planar duality on the tetrahedron')
    import numpy as np
    P = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    faces = {}                    # face opposite vertex i, as its counterclockwise (outward) cyclic order
    for i in range(4):
        f = [v for v in range(4) if v != i]
        a, b, c = P[f[0]], P[f[1]], P[f[2]]
        if np.dot(np.cross(b - a, c - a), a + b + c - 3 * P[i]) < 0:
            f = [f[0], f[2], f[1]]
        faces[i] = f
    def left_face(x, y):          # the face whose ccw boundary contains x -> y
        for i, f in faces.items():
            for k in range(3):
                if f[k] == x and f[(k + 1) % 3] == y:
                    return i
    kinds = Counter(); invol = True
    def dual(b):
        d = {}
        for x, y in itertools.permutations(range(4), 2):
            if b[(x, y)]:
                lf, rf = left_face(x, y), left_face(y, x)
                d[(lf, rf)] = True; d[(rf, lf)] = False
        return d
    for mask in range(64):
        b = tour(4, mask)
        d = dual(b)
        c3 = sum(cyclic3(b, t) for t in itertools.combinations(range(4), 3))
        outd = [sum(d[(v, w)] for w in range(4) if w != v) for v in range(4)]
        ss = sum(1 for o in outd if o in (0, 3))
        assert c3 == ss
        cls = lambda bb: ('transitive' if not any(cyclic3(bb, t) for t in itertools.combinations(range(4), 3))
                          else 'strong' if strong(4, bb) else 'middle')
        kinds[(cls(b), cls(d))] += 1
        dd = dual(d)
        invol &= all(dd[(x, y)] == (not b[(x, y)]) for x, y in itertools.permutations(range(4), 2)) or \
            all(dd[(x, y)] == b[(x, y)] for x, y in itertools.permutations(range(4), 2))
    print('  c3(T) = #sources + #sinks of the dual tournament T*, for all 64 labelled 4-tournaments')
    print('  (class of T, class of T*):', dict(kinds), '; dual of dual = T or its converse:', invol)
    assert kinds == Counter({('transitive', 'strong'): 24, ('strong', 'transitive'): 24, ('middle', 'middle'): 16})
    # transitive tournaments = flags of the tetrahedron = vertices of the truncated octahedron
    import networkx as nx
    trans = [m for m in range(64) if not any(cyclic3(tour(4, m), t) for t in itertools.combinations(range(4), 3))]
    G1 = nx.Graph(); G1.add_nodes_from(trans)
    G1.add_edges_from((a, c) for a, c in itertools.combinations(trans, 2) if bin(a ^ c).count('1') == 1)
    flags = [(v, frozenset(e), frozenset(f)) for f in itertools.combinations(range(4), 3) for e in itertools.combinations(f, 2)
             for v in e]
    G2 = nx.Graph(); G2.add_nodes_from(flags)
    G2.add_edges_from((p, q) for p, q in itertools.combinations(flags, 2) if sum(p[i] != q[i] for i in range(3)) == 1)
    print(f'  transitive 4-tournaments joined by single arc reversals: {G1.number_of_nodes()} vertices, '
          f'{G1.number_of_edges()} edges; isomorphic to the flag graph of the tetrahedron: {nx.is_isomorphic(G1, G2)} '
          f'(the graph of the truncated octahedron)')
    assert nx.is_isomorphic(G1, G2) and G1.number_of_edges() == 36


# ============================================================================ F
def section_F():
    print('== F. P7 on the torus')
    n = 7; QR = {1, 2, 4}
    b = {(x, y): (y - x) % 7 in QR for x in range(7) for y in range(7) if x != y}
    A = [tuple(sorted({i, (i + 1) % 7, (i + 3) % 7})) for i in range(7)]
    Bf = [tuple(sorted({i, (i + 2) % 7, (i + 3) % 7})) for i in range(7)]
    faces = A + Bf
    cyc = [t for t in itertools.combinations(range(7), 3) if cyclic3(b, t)]
    assert sorted(cyc) == sorted(faces) and len(cyc) == 14
    edge_faces = Counter(p for f in faces for p in itertools.combinations(f, 2))
    assert len(edge_faces) == 21 and set(edge_faces.values()) == {2}
    assert all(sum(1 for f in A if set(p) <= set(f)) == 1 for p in itertools.combinations(range(7), 2))
    for v in range(7):            # vertex links are single 6-cycles
        link = nx_cycle([tuple(x for x in f if x != v) for f in faces if v in f])
        assert link == 6
    print('  the 14 cyclic triangles of P7 are the faces {i,i+1,i+3}, {i,i+2,i+3} of the 7-vertex torus triangulation')
    print('  (each edge in one face of each kind; links are hexagons; V - E + F = 7 - 21 + 14 = 0); A-faces and B-faces')
    print('  are the lines of the Fano planes {0,1,3} + i and {0,2,3} + i')
    # coherent orientation: A ccw (i, i+1, i+3), B cw-read (i, i+3, i+2); every edge traversed oppositely
    oA = [(i, (i + 1) % 7, (i + 3) % 7) for i in range(7)]
    oB = [(i, (i + 3) % 7, (i + 2) % 7) for i in range(7)]
    dirs = Counter()
    for f in oA + oB:
        for k in range(3):
            dirs[(f[k], f[(k + 1) % 3])] += 1
    assert all(dirs[(x, y)] == 1 for x, y in itertools.permutations(range(7), 2))
    agree_A = all(b[(f[k], f[(k + 1) % 3])] for f in oA for k in range(3))
    disagree_B = all(not b[(f[k], f[(k + 1) % 3])] for f in oB for k in range(3))
    print('  the surface is orientable with A-faces read (i,i+1,i+3) and B-faces read (i,i+3,i+2); every arc of P7 runs')
    print('  along its A-face and against its B-face:', agree_A and disagree_B,
          '-> in the dual every A-face is a source and every B-face a sink')
    assert agree_A and disagree_B
    # face-cyclic orientations of the map
    pairs = list(itertools.combinations(range(7), 2))
    count = 0
    for mask in range(1 << 21):
        bb = {}
        for t, (x, y) in enumerate(pairs):
            bb[(x, y)] = bool(mask >> t & 1); bb[(y, x)] = not bb[(x, y)]
        if all(cyclic3(bb, f) for f in faces):
            count += 1
    print('  orientations of the 21 edges making all 14 faces directed triangles:', count, '(P7 and its converse)')
    assert count == 2
    aut = [g for g in itertools.permutations(range(7)) if all(b[(g[x], g[y])] == b[(x, y)] for x, y in b)]
    assert len(aut) == 21
    arc_orb = {frozenset((g[x], g[y]) for g in aut) for (x, y) in b if b[(x, y)]}
    face_orb = {frozenset(tuple(sorted(g[x] for x in f)) for g in aut) for f in faces}
    print(f'  |Aut(P7)| = {len(aut)}: orbits on vertices 1 (7), on arcs {len(arc_orb)} (21, regular), on faces '
          f'{len(face_orb)} ({sorted(len(o) for o in face_orb)}); since every edge borders one face of each kind, the two '
          f'kinds alternate around each vertex')
    assert len(arc_orb) == 1 and sorted(len(o) for o in face_orb) == [7, 7]
    # dual graph: Heawood
    dual_adj = defaultdict(set)
    for i, f in enumerate(faces):
        for j, g in enumerate(faces):
            if i < j and len(set(f) & set(g)) == 2:
                dual_adj[i].add(j); dual_adj[j].add(i)
    bip = all((i < 7) != (j < 7) for i in dual_adj for j in dual_adj[i])
    girth = min(len(c) for c in short_cycles(dual_adj, 8))
    print(f'  dual graph: {len(dual_adj)} vertices, cubic {all(len(v) == 3 for v in dual_adj.values())}, bipartite '
          f'(A | B) {bip}, girth {girth}: the Heawood graph (Fano incidence graph); dual orientation = points -> lines')
    assert bip and girth == 6
    # medial graph on the 21 edges
    med = defaultdict(set)
    for f in faces:
        es = [frozenset(p) for p in itertools.combinations(f, 2)]
        for e1, e2 in itertools.combinations(es, 2):
            med[e1].add(e2); med[e2].add(e1)
    regular_on_edges = len({frozenset((g[x], g[y])) for g in aut for (x, y) in [(0, 1)]}) == 21
    print(f'  medial graph: {len(med)} vertices, 4-regular {all(len(v) == 4 for v in med.values())}; Aut(P7) acts '
          f'regularly on its vertices: {regular_on_edges} (a Cayley graph of the Frobenius group of order 21)')
    H = hp(7, b)
    cycles = set()
    for r in itertools.permutations(range(1, 7)):
        c = (0,) + r
        if all(b[(c[i], c[(i + 1) % 7])] for i in range(7)):
            cycles.add(frozenset((c[i], c[(i + 1) % 7]) for i in range(7)))
    cyc_orbits = {frozenset(frozenset((g[x], g[y]) for x, y in C) for g in aut) for C in cycles}
    print(f'  H(P7) = {H} = {H // 21} x |Aut(P7)|; Hamiltonian cycles {len(cycles)}, Aut-orbits of sizes '
          f'{sorted(len(o) for o in cyc_orbits)} (the size-3 orbit: the circulant cycles x -> x + d, d in {{1,2,4}})')
    assert H == 189 and len(cycles) == 24 and sorted(len(o) for o in cyc_orbits) == [3, 21]
    circ = frozenset(frozenset((x, (x + d) % 7) for x in range(7)) for d in (1, 2, 4))
    assert circ in cyc_orbits


def nx_cycle(edges):
    adj = defaultdict(set)
    for x, y in edges:
        adj[x].add(y); adj[y].add(x)
    if any(len(v) != 2 for v in adj.values()):
        return -1
    start = next(iter(adj)); prev, cur, L = None, start, 0
    while True:
        nxt = next(w for w in adj[cur] if w != prev)
        prev, cur, L = cur, nxt, L + 1
        if cur == start:
            return L if L == len(adj) else -1


def short_cycles(adj, maxlen):
    out = []
    for s in adj:
        stack = [(s, [s])]
        while stack:
            v, path = stack.pop()
            for w in adj[v]:
                if w == s and len(path) >= 3:
                    out.append(path)
                elif w not in path and len(path) < maxlen and w > s:
                    stack.append((w, path + [w]))
    return out


# ============================================================================ G
def K_of(w):
    return sum(1 << max(0, a - b - 2) for b in range(len(w)) if w[b] for a in range(b + 1, len(w)) if not w[a])


def section_G():
    print('== G. tournaments with T - v transitive: H = 1 + 2 K(w)')
    for N in range(1, 8):
        for w in itertools.product((0, 1), repeat=N):
            b = {}
            for x in range(N + 1):
                for y in range(N + 1):
                    if x != y:
                        b[(x, y)] = bool(w[y]) if x == N else (not w[x]) if y == N else x < y
            assert hp(N + 1, b) == 1 + 2 * K_of(w)
    print('  every odd cycle passes through v, so Omega is complete and H = 1 + 2 alpha_1, where alpha_1 =')
    print('  K(w) = sum over 10-inversions b < a of w (w_b = 1: v beats b; w_a = 0: a beats v) of 2^max(0, a-b-2):')
    print('  verified against the path DP for every word of length <= 7')
    vals = set()
    for N in range(15):
        for bits in range(1 << N):
            vals.add(K_of([(bits >> i) & 1 for i in range(N)]))
    miss = [k for k in range(256) if k not in vals]
    print(f'  attained K below 256 (words of length <= 14): {256 - len(miss)} of 256; first missing K {miss[:12]}')
    assert miss[:2] == [3, 6] and 10 in miss
    print('  K = 3 and K = 10 (H = 7, 21) are missing, as they must be; the family also misses about half of all values')


# ============================================================================ H
def section_H():
    print('== H. 60 = 2^2 3 5 = |A5|')
    N = 60
    divs = [d for d in sp.divisors(N) if d not in (1, N)]
    primes = [d for d in divs if sp.isprime(d)]
    sqfree = [d for d in divs if not sp.isprime(d) and sp.factorint(d) and max(sp.factorint(d).values()) == 1]
    sqcont = [d for d in divs if max(sp.factorint(d).values()) >= 2]
    print(f'  proper divisors {divs}: primes {primes}, squarefree composites {sqfree}, square-containing {sqcont} '
          f'({len(primes)}-{len(sqfree)}-{len(sqcont)})')
    assert (len(primes), len(sqfree), len(sqcont)) == (3, 4, 3)
    point = {2: 30, 3: 20, 5: 12}        # rotation order k: point orbit 60/k (edge midpoints, face centres, vertices)
    axis = {4: 15, 6: 10, 10: 6}         # stabilizer D_k of an axis: axis orbit 60/2k
    used = set(point) | set(point.values()) | set(axis) | set(axis.values())
    print('  icosahedral rotation group: point stabilizers C2, C3, C5 with orbits 30, 20, 12 (edges, dodecahedron')
    print('  vertices, icosahedron vertices); axis stabilizers D2, D3, D5 with 15, 10, 6 axes')
    print('  stabilizer orders and orbit sizes together are exactly the ten proper divisors:', sorted(used) == divs)
    assert sorted(used) == divs
    print('  V - E + F = 12 - 30 + 20 =', 12 - 30 + 20, '; passing from points to axes (the antipodal quotient)')
    print('  doubles each stabilizer (2 -> 4, 3 -> 6, 5 -> 10) and halves each orbit (30 -> 15, 20 -> 10, 12 -> 6)')
    for t in ((2, 2, 5), (2, 3, 3), (2, 3, 4), (2, 3, 5), (2, 3, 6), (2, 3, 7)):
        s = sum(sp.Rational(1, k) for k in t)
        order = 2 / (s - 1) if s > 1 else None
        print(f'  triangle group {t}: 1/p + 1/q + 1/r = {s} -> {"rotation group of order " + str(order) if order else ("flat (torus maps, e.g. K7 on the torus is {3,6})" if s == 1 else "hyperbolic (the Klein quartic is (2,3,7), automorphisms PSL(2,7))")}')


# ============================================================================ I
def section_I():
    print('== I. the Collatz side')
    def syr(m):
        m = 3 * m + 1
        while m % 2 == 0: m //= 2
        return m
    for m in (7, 21):
        orb = [m]
        while orb[-1] != 1: orb.append(syr(orb[-1]))
        print(f'  Syracuse orbit of {m} ({m} = {m % 4} mod 4): {orb}')
    print('  quadratic residues mod 7:', sorted({(x * x) % 7 for x in range(1, 7)}), '= the trivial Collatz cycle {1, 2, 4}')
    print('  1 + 2 + 4 =', 1 + 2 + 4, '; 1 + 4 + 16 =', 1 + 4 + 16, '; q^2 + q + 1 at q = 2, 4:', 2 * 2 + 2 + 1, 4 * 4 + 4 + 1)
    print('  verdict: NUMEROLOGY (no map between Syracuse dynamics and Hamiltonian path counts)')


if __name__ == '__main__':
    section_A()
    section_BC()
    section_D()
    section_E()
    section_F()
    section_G()
    section_H()
    section_I()
    print('ALL CHECKS PASSED')
