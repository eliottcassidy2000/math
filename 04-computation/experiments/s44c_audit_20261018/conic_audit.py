#!/usr/bin/env python3
"""Audit check (e)(ii)-(v) of the merge note pythagorean_form_conic_graphs_rp_margin_20261018.md (independent code).

Conic C: x^2 + y^2 - z^2 = 0 in PG(2,p).  Exterior = on two tangents (polar is a secant); interior = polar is a
passant.  G_int, G_ext = conjugacy graphs.  Checks:
  A. p = 5: G_int ~ Petersen (distance-regular {3,2;1,1}), G_ext ~ 5K3; exterior points <-> duads; interior points <->
     synthemes (3 secants through the point); the 5 all-exterior self-polar triangles form a synthematic total;
     self-polar triangle census; Petersen edges <-> exterior points; M_x (exterior points on the tangent at x) is a
     perfect matching whose contraction is K5 labelled by the conic points y != x; deletion of an interior point;
     the Stab(P) actions; equivariant bijections interior points <-> 3+3 splits (existence and uniqueness, for the
     Berggren image A5 and for the full SO(3,5) = PGL(2,5)); matchings / conic points / hemi-dodecahedron faces.
  B. p = 7: G_ext distance-regular with intersection array {3,2,2,1;1,1,1,2} (which characterizes the Coxeter graph,
     BCN) and isomorphic to an independent Coxeter construction; G_int ~ L(Heawood); exhaustive non-Hamiltonicity of
     Petersen and Coxeter.
  C. p <= 31: orders, degrees (p+1)/2, (p-1)/2, triangle pattern by p mod 4, connectivity; Hamiltonian cycles found
     and verified for the connected non-cubic graphs with p <= 13.
  D. Berggren image mod p via 2x2 matrices A' = [[2,-1],[1,0]], B' = [[2,1],[1,0]], C' = [[1,2],[0,1]]: exact
     Veronese intertwining with the 3x3 Berggren matrices, C'^-1 A' = S, S C' S^-1, image order in PGL(2,p)
     compared with |PSL(2,p)| / |PGL(2,p)| for odd primes p <= 31, and the projective 3x3 image on PG(2,p).
Output: conic_audit.out / conic_audit.json.
"""
import itertools, json, os, random, sys
from collections import Counter
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
R = {}


def pts_PG2(p):
    out = []
    for v in itertools.product(range(p), repeat=3):
        if v == (0, 0, 0):
            continue
        f = next(x for x in v if x)
        if f == 1:
            out.append(v)
    return out


def nrm(v, p):
    f = next(x % p for x in v if x % p)
    i = pow(f, -1, p)
    return tuple(x * i % p for x in v)


def B(u, v, p):
    return (u[0] * v[0] + u[1] * v[1] - u[2] * v[2]) % p


def classify(p):
    P = pts_PG2(p)
    con = [v for v in P if B(v, v, p) == 0]
    ext, inte = [], []
    for v in P:
        if B(v, v, p) == 0:
            continue
        # number of tangents through v = number of conic points x with v on the tangent at x, i.e. B(v,x) = 0
        k = sum(1 for x in con if B(v, x, p) == 0)
        assert k in (0, 2)
        (ext if k == 2 else inte).append(v)
    return P, con, ext, inte


def conj_graph(vs, p):
    G = nx.Graph()
    G.add_nodes_from(vs)
    for u, v in itertools.combinations(vs, 2):
        if B(u, v, p) == 0:
            G.add_edge(u, v)
    return G


def intersection_array(G):
    """Return (b, c) if G is distance-regular, else None."""
    D = dict(nx.all_pairs_shortest_path_length(G))
    diam = max(max(d.values()) for d in D.values())
    b = [None] * (diam + 1)
    c = [None] * (diam + 1)
    for u in G:
        for v in G:
            i = D[u][v]
            bi = sum(1 for w in G[v] if D[u][w] == i + 1)
            ci = sum(1 for w in G[v] if D[u][w] == i - 1)
            if b[i] is None:
                b[i], c[i] = bi, ci
            elif (b[i], c[i]) != (bi, ci):
                return None
    return tuple(b[:-1]), tuple(c[1:])


def line_pts(u, v, p):
    s = set()
    for a in range(p):
        for b in range(p):
            if (a, b) != (0, 0):
                s.add(nrm(tuple((a * x + b * y) % p for x, y in zip(u, v)), p))
    return s


def cross(u, v, p):
    return nrm(((u[1] * v[2] - u[2] * v[1]) % p, (u[2] * v[0] - u[0] * v[2]) % p, (u[0] * v[1] - u[1] * v[0]) % p), p)


def pole_of_line_through(u, v, p):
    """The pole of the line uv w.r.t. B: the point w with B(w, u) = B(w, v) = 0 (w = J (u x v), J = diag(1,1,-1))."""
    c = cross(u, v, p)
    return nrm((c[0], c[1], (-c[2]) % p), p)


# ------------------------------------------------------------------ group machinery (2x2 and 3x3)
def mat2_mul(X, Y, p):
    return ((X[0][0] * Y[0][0] + X[0][1] * Y[1][0]) % p, (X[0][0] * Y[0][1] + X[0][1] * Y[1][1]) % p), \
           ((X[1][0] * Y[0][0] + X[1][1] * Y[1][0]) % p, (X[1][0] * Y[0][1] + X[1][1] * Y[1][1]) % p)


def proj2(X, p):
    """Canonical representative of X in PGL(2,p): scale so that the first nonzero entry is 1."""
    flat = [X[0][0], X[0][1], X[1][0], X[1][1]]
    f = next(x for x in flat if x % p)
    i = pow(f, -1, p)
    return tuple(tuple(x * i % p for x in row) for row in X)


def generate(gens, mul, canon):
    e = canon(((1, 0), (0, 1))) if len(gens[0]) == 2 else None
    group = {canon(g) for g in gens}
    frontier = list(group)
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = canon(mul(s, g))
                if h not in group:
                    group.add(h)
                    new.append(h)
        frontier = new
    return group


def sym2(X):
    """3x3 matrix of the action of X on (m^2 - n^2, 2mn, m^2 + n^2) (exact over Z), via X acting on (m, n)."""
    import sympy as sp
    m, n = sp.symbols('m n')
    a, b = X[0]
    c, d = X[1]
    M2, N2 = a * m + b * n, c * m + d * n
    img = [sp.expand(M2 ** 2 - N2 ** 2), sp.expand(2 * M2 * N2), sp.expand(M2 ** 2 + N2 ** 2)]
    Xs, Ys, Zs = sp.symbols('X Y Z')
    sub = {m ** 2: (Xs + Zs) / 2, n ** 2: (Zs - Xs) / 2, m * n: Ys / 2}
    rows = []
    for e in img:
        poly = sp.Poly(e, m, n)
        expr = 0
        for (i, j), co in poly.terms():
            key = m ** i * n ** j
            expr += co * sub[key]
        expr = sp.expand(expr)
        rows.append([sp.Rational(expr.coeff(Xs)), sp.Rational(expr.coeff(Ys)), sp.Rational(expr.coeff(Zs))])
    if all(x.q == 1 for r in rows for x in r):
        return tuple(tuple(int(x) for x in r) for r in rows)
    return tuple(tuple(x for x in r) for r in rows)


A3 = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
B3 = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
C3 = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
A2 = ((2, -1), (1, 0))
B2 = ((2, 1), (1, 0))
C2 = ((1, 2), (0, 1))


def act3(M, v, p):
    return nrm(tuple(sum(M[i][k] * v[k] for k in range(3)) % p for i in range(3)), p)


def perm_group_on_points(p, P, mats):
    idx = {v: i for i, v in enumerate(P)}
    gens = [tuple(idx[act3(M, v, p)] for v in P) for M in mats]
    ident = tuple(range(len(P)))
    G = {ident}
    fr = [ident]
    while fr:
        new = []
        for g in fr:
            for s in gens:
                h = tuple(s[g[i]] for i in range(len(P)))
                if h not in G:
                    G.add(h)
                    new.append(h)
        fr = new
    return G, idx


def SO3_generators(p):
    """Generators of SO(3, F_p) for our form via Sym^2 of PGL(2,p) generators (all of GL(2,p) mod scalars):
    Sym^2(X)/det(X) for X in {[[1,1],[0,1]], [[1,0],[1,1]], [[g,0],[0,1]]} with g a primitive root."""
    g = next(x for x in range(2, p) if all(pow(x, (p - 1) // r, p) != 1 for r in prime_factors(p - 1)))
    mats = []
    for X in (((1, 1), (0, 1)), ((1, 0), (1, 1)), ((g, 0), (0, 1))):
        S = sym2(X)
        det = (X[0][0] * X[1][1] - X[0][1] * X[1][0]) % p
        inv = pow(det, -1, p)
        def redq(x):
            from fractions import Fraction
            x = Fraction(int(getattr(x, 'p', x)), int(getattr(x, 'q', 1)))
            return x.numerator % p * pow(x.denominator % p, -1, p) % p
        mats.append(tuple(tuple(redq(x) * inv % p for x in row) for row in S))
    return mats


def prime_factors(n):
    f, d = set(), 2
    while d * d <= n:
        while n % d == 0:
            f.add(d)
            n //= d
        d += 1
    if n > 1:
        f.add(n)
    return f


def exhaustive_hamiltonian(G):
    """Exact Hamiltonian-cycle test: DFS from a fixed vertex with the pruning rule that every unvisited vertex must
    keep >= 2 available neighbours (unvisited, or an end of the current path)."""
    nodes = list(G)
    n = len(nodes)
    s = nodes[0]
    sys.setrecursionlimit(100000)
    path, on = [s], {s}

    def feasible():
        ends = (path[0], path[-1])
        for w in nodes:
            if w in on:
                continue
            avail = sum(1 for x in G[w] if x not in on or x in ends)
            if avail < 2:
                return False
        return True

    def dfs():
        if len(path) == n:
            return s in G[path[-1]]
        if not feasible():
            return False
        last = path[-1]
        for w in G[last]:
            if w not in on:
                path.append(w); on.add(w)
                if dfs():
                    return True
                path.pop(); on.discard(w)
        return False
    return dfs()


def find_ham_cycle(G, tries=400, seed=7):
    rng = random.Random(seed)
    nodes = list(G)
    n = len(nodes)
    sys.setrecursionlimit(100000)
    for t in range(tries):
        s = rng.choice(nodes)
        path, on = [s], {s}
        budget = [50000]

        def dfs():
            budget[0] -= 1
            if budget[0] < 0:
                return False
            if len(path) == n:
                return s in G[path[-1]]
            cand = [w for w in G[path[-1]] if w not in on]
            rng.shuffle(cand)
            cand.sort(key=lambda w: sum(1 for x in G[w] if x not in on))
            for w in cand:
                path.append(w); on.add(w)
                if dfs():
                    return True
                path.pop(); on.discard(w)
            return False
        if dfs():
            cyc = list(path)
            assert len(set(cyc)) == n and all(G.has_edge(cyc[i], cyc[(i + 1) % n]) for i in range(n))
            return cyc
    return None


def part_A():
    p = 5
    P, con, ext, inte = classify(p)
    Gi, Ge = conj_graph(inte, p), conj_graph(ext, p)
    out = {'sizes': (len(con), len(ext), len(inte))}
    out['Gint_petersen'] = nx.is_isomorphic(Gi, nx.petersen_graph())
    out['Gint_intersection_array'] = intersection_array(Gi)
    comps = [Ge.subgraph(c) for c in nx.connected_components(Ge)]
    out['Gext_5K3'] = len(comps) == 5 and all(nx.is_isomorphic(c, nx.complete_graph(3)) for c in comps)
    duad = {E: frozenset(x for x in con if B(E, x, p) == 0) for E in ext}
    out['ext_to_duads_bijective'] = len(set(duad.values())) == 15
    # synthemes of interior points: the secants through P
    synth = {}
    for Pt in inte:
        pairs = set()
        for x in con:
            L = line_pts(Pt, x, p)
            others = [y for y in con if y in L and y != x]
            assert len(others) == 1
            pairs.add(frozenset((x, others[0])))
        synth[Pt] = frozenset(pairs)
        assert len(synth[Pt]) == 3 and len(set().union(*synth[Pt])) == 6
    # polar of P contains exactly the 3 exterior points whose duads are the secants through P
    polar_match = all(frozenset(duad[E] for E in ext if B(E, Pt, p) == 0) == synth[Pt] for Pt in inte)
    out['syntheme_of_P_equals_duads_of_exterior_points_on_polar'] = polar_match
    tris = [t for t in itertools.combinations(ext + inte, 3) if all(B(a, b, p) == 0 for a, b in itertools.combinations(t, 2))]
    out['selfpolar_triangles'] = len(tris)
    out['selfpolar_types_by_interior_count'] = dict(Counter(sum(1 for x in t if x in inte) for t in tris))
    eee = [t for t in tris if all(x in ext for x in t)]
    total = [frozenset(duad[x] for x in t) for t in eee]
    out['EEE_triangles_are_synthemes'] = all(len(set().union(*s)) == 6 for s in total)
    out['EEE_synthemes_form_a_total'] = len(set().union(*total)) == 15
    out['interior_and_EEE_synthemes_all_15'] = len(set(synth.values()) | set(total)) == 15
    # Petersen edges <-> exterior points
    e2x = {}
    for u, v in Gi.edges():
        w = pole_of_line_through(u, v, p)
        assert B(w, u, p) == 0 and B(w, v, p) == 0
        e2x[frozenset((u, v))] = w
    out['edge_label_exterior'] = all(w in ext for w in e2x.values()) and len(set(e2x.values())) == 15
    x2e = {w: e for e, w in e2x.items()}
    # M_x for each conic point x
    mres = []
    for x in con:
        tangent = [E for E in ext if B(E, x, p) == 0]
        Mx = [x2e[E] for E in tangent]
        vs = set().union(*Mx)
        perfect = len(Mx) == 5 and len(vs) == 10
        rep = {}
        for e in Mx:
            for v in e:
                rep[v] = e
        Kc = nx.Graph()
        Kc.add_nodes_from(Mx)
        for u, v in Gi.edges():
            if frozenset((u, v)) not in Mx:
                Kc.add_edge(rep[u], rep[v])
        isK5 = Kc.number_of_edges() == 10 and nx.is_isomorphic(Kc, nx.complete_graph(5))
        labels = [next(iter(duad[E] - {x})) for E in tangent]
        mres.append(perfect and isK5 and len(set(labels)) == 5 and x not in labels)
    out['M_x_perfect_matching_contracts_to_K5_labelled_by_other_conic_points'] = all(mres)
    pms = [frozenset(frozenset(e) for e in comb) for comb in itertools.combinations(Gi.edges(), 5)
           if len(set().union(*[set(e) for e in comb])) == 10]
    out['petersen_perfect_matchings'] = len(pms)
    # groups: Berggren image (A5) and full SO(3,5)
    Gb, idx = perm_group_on_points(p, P, (A3, B3, C3))
    Gso, idx2 = perm_group_on_points(p, P, SO3_generators(p))
    out['berggren_image_order'] = len(Gb)
    out['SO3_order'] = len(Gso)
    pts = P

    def act(g, v):
        return pts[g[idx[v]]]
    splits = set()
    for T in itertools.combinations(con, 3):
        T = frozenset(T)
        splits.add(frozenset((T, frozenset(con) - T)))
    splits = list(splits)

    def act_split(g, s):
        return frozenset(frozenset(act(g, x) for x in T) for T in s)

    def equivariant_bijections(G):
        G = list(G)
        cnt = []
        P0 = inte[0]
        for s0 in splits:
            phi = {}
            ok = True
            for g in G:
                a, b = act(g, P0), act_split(g, s0)
                if a in phi and phi[a] != b:
                    ok = False
                    break
                phi[a] = b
            if ok and len(phi) == 10 and len(set(phi.values())) == 10:
                # full equivariance check
                if all(phi[act(g, Q)] == act_split(g, phi[Q]) for g in G for Q in inte):
                    cnt.append(phi)
        return cnt
    for name, G in (('A5_berggren', Gb), ('SO3_full', Gso)):
        bij = equivariant_bijections(G)
        out[f'equivariant_bijections_interior_to_splits_{name}'] = len(bij)
        if bij:
            phi = bij[0]
            # orbit description: split of P = orbits of Stab_A5(P) on the conic points?
            if name == 'A5_berggren':
                desc = []
                for Q in inte:
                    st = [g for g in G if act(g, Q) == Q]
                    orbs = set()
                    for x in con:
                        orbs.add(frozenset(act(g, x) for g in st))
                    desc.append(frozenset(orbs) == phi[Q])
                out['split_equals_orbits_of_StabA5'] = all(desc)
                out['StabA5_interior_order'] = len([g for g in G if act(g, inte[0]) == inte[0]])
                # relation between split and syntheme of P
                rel = Counter()
                for Q in inte:
                    T, Tc = tuple(phi[Q])
                    rel[tuple(sorted((len(d & T), len(d & Tc)) for d in synth[Q]))] += 1
                out['split_vs_syntheme_pattern'] = {str(k): v for k, v in rel.items()}
            # Petersen adjacency on splits: intersection pattern of adjacent / nonadjacent splits
            adjpat, nonpat = Counter(), Counter()
            for Q1, Q2 in itertools.combinations(inte, 2):
                T1 = next(iter(phi[Q1]))
                pat = tuple(sorted(len(T1 & U) for U in phi[Q2]))
                (adjpat if Gi.has_edge(Q1, Q2) else nonpat)[pat] += 1
            out[f'split_intersection_patterns_{name}'] = {'adjacent': {str(k): v for k, v in adjpat.items()}, 'nonadjacent': {str(k): v for k, v in nonpat.items()}}
    # deletion of an interior point; Stab actions
    Pt = inte[0]
    dist = nx.single_source_shortest_path_length(Gi, Pt)
    far = [v for v in inte if dist[v] == 2]
    st = [g for g in Gb if act(g, Pt) == Pt]
    far_orbit_sizes = sorted(Counter(frozenset(act(g, v) for g in st) for v in far).keys(), key=len)
    con_orbits = set(frozenset(act(g, x) for g in st) for x in con)
    regular_on_far = len(st) == 6 and all(sum(1 for g in st if act(g, v) == v) == 1 for v in far)
    out['StabA5(P)_regular_on_far_vertices'] = regular_on_far
    out['StabA5(P)_orbits_on_conic'] = sorted(len(o) for o in con_orbits)
    # K3,3 bipartition preserved/swapped by involutions vs action on the two conic orbits
    K = nx.Graph()
    K.add_nodes_from(far)
    K.add_edges_from(Gi.subgraph(far).edges())
    for nb in Gi[Pt]:
        ends = [w for w in Gi[nb] if w != Pt]
        K.add_edge(*ends)
    out['deletion_K33'] = nx.is_isomorphic(K, nx.complete_bipartite_graph(3, 3))
    X, Y = nx.bipartite.sets(K)
    swap_bip = [frozenset(act(g, v) for v in X) != frozenset(X) for g in st]
    T0 = next(iter(con_orbits))
    swap_con = [frozenset(act(g, x) for x in T0) != T0 for g in st]
    out['some_Stab_element_swaps_K33_sides_but_fixes_conic_orbits'] = any(a and not b for a, b in zip(swap_bip, swap_con))
    R['p5'] = out


def coxeter_independent():
    """Coxeter graph from the 28 antiflags of the Fano plane: (point, line) with point not on line;
    (P,L) ~ (P',L') iff P' not on L... use the classical description: vertices = antiflags, adjacent iff
    P, P' and the point L cap L' are ... -- we instead use KG(7,3) minus Fano lines with a DIFFERENT Fano plane."""
    lines = [frozenset(l) for l in ((0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5))]
    verts = [frozenset(s) for s in itertools.combinations(range(7), 3) if frozenset(s) not in lines]
    G = nx.Graph()
    G.add_nodes_from(verts)
    for u, v in itertools.combinations(verts, 2):
        if not (u & v):
            G.add_edge(u, v)
    return G


def part_B():
    p = 7
    P, con, ext, inte = classify(p)
    Gi, Ge = conj_graph(inte, p), conj_graph(ext, p)
    out = {'sizes': (len(con), len(ext), len(inte))}
    out['Gext_intersection_array'] = intersection_array(Ge)
    out['Gext_is_Coxeter_array_3221_1112'] = out['Gext_intersection_array'] == ((3, 2, 2, 1), (1, 1, 1, 2))
    out['Gext_iso_independent_Coxeter'] = nx.is_isomorphic(Ge, coxeter_independent())
    out['Gint_iso_LHeawood'] = nx.is_isomorphic(Gi, nx.line_graph(nx.heawood_graph()))
    out['Gint_intersection_array'] = intersection_array(Gi)
    out['Coxeter_hamiltonian_exhaustive'] = exhaustive_hamiltonian(Ge)
    Gp = conj_graph(classify(5)[3], 5)
    out['Petersen_hamiltonian_exhaustive'] = exhaustive_hamiltonian(Gp)
    Gb, idx = perm_group_on_points(p, P, (A3, B3, C3))
    out['berggren_image_order_on_PG2'] = len(Gb)
    out['Aut_Gext_order'] = sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(Ge, Ge).isomorphisms_iter())
    R['p7'] = out


def part_C(primes=(3, 5, 7, 11, 13, 17, 19, 23, 29, 31)):
    rows = []
    for p in primes:
        P, con, ext, inte = classify(p)
        row = {'p': p, 'p_mod_4': p % 4, 'n_int': len(inte), 'n_ext': len(ext)}
        for name, vs in (('int', inte), ('ext', ext)):
            G = conj_graph(vs, p)
            degs = sorted(set(d for _, d in G.degree()))
            tri = sum(nx.triangles(G).values()) // 3
            row[name] = {'degree': degs, 'triangles': tri, 'connected': nx.is_connected(G)}
            if p <= 13 and nx.is_connected(G) and degs != [3]:
                row[name]['hamiltonian_cycle_found_and_verified'] = find_ham_cycle(G) is not None
        row['deg_ok'] = row['int']['degree'] == [(p + 1) // 2] and row['ext']['degree'] == [(p - 1) // 2]
        row['triangle_rule_ok'] = ((row['int']['triangles'] == 0) == (p % 4 == 1)) and ((row['ext']['triangles'] == 0) == (p % 4 == 3))
        rows.append(row)
        print(f"[C] p={p}: int deg {row['int']['degree']} tri {row['int']['triangles']} conn {row['int']['connected']}"
              f" ham {row['int'].get('hamiltonian_cycle_found_and_verified')}; ext deg {row['ext']['degree']} tri "
              f"{row['ext']['triangles']} conn {row['ext']['connected']} ham {row['ext'].get('hamiltonian_cycle_found_and_verified')};"
              f" degree rule {row['deg_ok']}, triangle rule {row['triangle_rule_ok']}", flush=True)
    R['survey'] = rows


def part_D(primes=(3, 5, 7, 11, 13, 17, 19, 23, 29, 31)):
    out = {}
    out['sym2_A'] = sym2(A2) == A3
    out['sym2_B'] = sym2(B2) == B3
    out['sym2_C'] = sym2(C2) == C3
    import sympy as sp
    Am, Cm = sp.Matrix(A2), sp.Matrix(C2)
    S = Cm.inv() * Am
    out['Cinv_A_equals_S'] = S == sp.Matrix([[0, -1], [1, 0]])
    out['S_C_Sinv'] = (S * Cm * S.inv()).tolist()
    out['dets'] = [sp.Matrix(X).det() for X in (A2, B2, C2)]
    rows = []
    for p in primes:
        G = generate([A2, B2, C2], lambda X, Y: mat2_mul(X, Y, p), lambda X: proj2(X, p))
        psl = p * (p * p - 1) // 2
        pgl = p * (p * p - 1)
        expect = psl if p % 4 == 1 else pgl
        Pp = pts_PG2(p)
        G3, _ = perm_group_on_points(p, Pp, (A3, B3, C3)) if p <= 13 else (None, None)
        rows.append({'p': p, 'order_PGL2_image': len(G), 'PSL': psl, 'PGL': pgl, 'matches_prop': len(G) == expect,
                     'order_on_PG2': len(G3) if G3 is not None else None})
        print(f"[D] p={p}: <A',B',C'> in PGL(2,p) has order {len(G)}; |PSL|={psl}, |PGL|={pgl}; proposition predicts "
              f"{expect}: {len(G) == expect}; 3x3 projective image on PG(2,p): {len(G3) if G3 is not None else '-'}", flush=True)
    out['rows'] = rows
    R['berggren_image'] = out


if __name__ == '__main__':
    part_A()
    print('[A] p=5:', json.dumps(R['p5'], default=str, indent=1), flush=True)
    part_B()
    print('[B] p=7:', json.dumps(R['p7'], default=str, indent=1), flush=True)
    part_C()
    part_D()
    print('[D]', {k: v for k, v in R['berggren_image'].items() if k != 'rows'}, flush=True)
    json.dump(R, open(os.path.join(HERE, 'conic_audit.json'), 'w'), indent=1, default=str)
