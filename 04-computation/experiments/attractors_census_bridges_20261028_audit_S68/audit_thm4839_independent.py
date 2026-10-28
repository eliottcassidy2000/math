"""Independent referee checks of THM-4839 (S68 audit, 2026-10-28).

Written from scratch; it does not import or copy the session's script.

  A  complete unit Cayley graphs Cay(O_K/P, image of O_K^x) over ALL imaginary quadratic fields with |D| <= DMAX
     (any class number) and all primes p <= PMAX; explicit graphs for N(P) <= 60
  B  toroidal maps built as genuine maps (darts, rotation, edge involution, faces) in Coxeter-Moser's convention:
     {4,4}_(b,c) = square tiling / <(b,c),(-c,b)>; {3,6}_(b,c) = triangular tiling / <b + c e^{i pi/3}>;
     {6,3}_(b,c) = dual of {3,6}_(b,c). Skeleton identified with networkx isomorphism; Euler characteristic,
     orientable regularity (monodromy acts regularly on darts) and chirality (mirror automorphism) computed.
  C  Petersen = hemi-dodecahedron: antipodal quotient of the dodecahedral map, its 6 faces, and the order of the map
     automorphism group (automorphisms of the Petersen graph preserving the face set)
  D  the Klein-four Belyi map lambda(t) = ((1-t^2)/(1+t^2))^2: fibres with multiplicities (including t = oo),
     numerical monodromy, Riemann-Hurwitz; identification with X(4) -> X(2) through theta functions
  E  Pythagorean Frey curves: order-4 point, full torsion (Z/2xZ/4 vs Z/2xZ/8), the Z/2xZ/4 family is strictly larger
  F  Berggren abc hits (independent (m,n)-tree generator)
  G  Collatz clocks (convergents and semiconvergents of log2 3): hit test and exact qualities
  H  1/p+1/q+1/r <= 41/42 on hyperbolic signatures; the S116 inequality at 1+8=9
"""
import itertools
import math
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import networkx as nx
import sympy as sp
from sympy import factorint, isprime

T0 = time.time()
DMAX = 2000
PMAX = 200


def say(*a):
    print(*a, flush=True)


def rad(n):
    r = 1
    for p in factorint(n):
        r *= p
    return r


# ====================================================================== A
say("=" * 100)
say("A. complete unit Cayley graphs over every imaginary quadratic field with |D| <= %d, primes p <= %d" % (DMAX, PMAX))
say("=" * 100)


def is_fund_disc(D):
    if D >= 0:
        return False
    if D % 4 == 1:
        return sp.factorint(-D) and all(e == 1 for e in sp.factorint(-D).values())
    if D % 4 == 0:
        m = D // 4
        if m % 4 in (2, 3):
            return all(e == 1 for e in sp.factorint(-m).values())
    return False


def class_number(D):
    """number of reduced primitive positive definite forms of discriminant D (D < 0)."""
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if c == a and b < 0:
                continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    return h


def residue_fields(D, p):
    """O_K = Z[w], w = (D + sqrt D)/2, min poly x^2 - D x + (D^2 - D)/4.
    Returns a list of (N, elements-of-unit-image as hashable residue representatives)."""
    c0 = (D * D - D) // 4
    # units as polynomials u0 + u1 w
    if D == -4:      # w = (-4 + 2i)/2 = -2 + i, i = w + 2
        units = [(1, 0), (-1, 0), (2, 1), (-2, -1)]
    elif D == -3:    # w = (-3 + sqrt-3)/2 = -1 + (1+sqrt-3)/2 = -1 - omega^2 ... 6th roots: zeta6 = w + 2
        z = (2, 1)   # zeta6 = (1 + sqrt-3)/2 = w + 2
        units = []
        cur = (1, 0)
        for _ in range(6):
            units.append(cur)
            # cur * (w + 2), with w^2 = D w - c0
            a0, a1 = cur
            # (a0 + a1 w)(2 + w) = 2a0 + (a0 + 2a1) w + a1 w^2 = 2a0 - a1 c0 + (a0 + 2 a1 + a1 D) w
            cur = (2 * a0 - a1 * c0, a0 + 2 * a1 + a1 * D)
        assert cur == (1, 0)
    else:
        units = [(1, 0), (-1, 0)]
    roots = [r for r in range(p) if (r * r - D * r + c0) % p == 0]
    out = []
    if roots:
        for r in sorted(set(roots)):
            img = frozenset(((u0 + u1 * r) % p,) for u0, u1 in units)
            out.append((p, img, "split" if len(set(roots)) == 2 else "ramified"))
    else:
        img = frozenset((u0 % p, u1 % p) for u0, u1 in units)
        out.append((p * p, img, "inert"))
    return out


def field_elements(N, p):
    if N == p:
        return [(x,) for x in range(p)]
    return [(x, y) for x in range(p) for y in range(p)]


def add_el(x, y, p):
    return tuple((a + b) % p for a, b in zip(x, y))


hits = []
nfields = 0
hist_h = {}
for D in range(-3, -DMAX - 1, -1):
    if not is_fund_disc(D):
        continue
    nfields += 1
    h = class_number(D)
    hist_h[h] = hist_h.get(h, 0) + 1
    for p in sp.primerange(2, PMAX + 1):
        for N, img, kind in residue_fields(D, p):
            nonzero_img = [x for x in img if any(x)]
            assert len(nonzero_img) == len(img)          # units never vanish mod P
            complete_by_count = (len(img) == N - 1)
            if N <= 60:                                  # explicit Cayley graph
                els = field_elements(N, p)
                G = nx.Graph()
                G.add_nodes_from(els)
                for x in els:
                    for s in img:
                        y = add_el(x, s, p)
                        if y != x:
                            G.add_edge(x, y)
                complete_graph = (G.number_of_edges() == N * (N - 1) // 2)
                assert complete_graph == complete_by_count
            if complete_by_count and N >= 4:
                hits.append((D, h, p, N, kind))
say("  fundamental discriminants scanned: %d (class numbers 1..%d; e.g. h histogram head %s)"
    % (nfields, max(hist_h), sorted(hist_h.items())[:6]))
say("  complete unit Cayley graphs with N(P) >= 4:")
for D, h, p, N, kind in hits:
    say("    D = %d (h = %d), p = %d (%s), N(P) = %d  -> K%d" % (D, h, p, kind, N, N))
assert sorted((D, N) for D, h, p, N, kind in hits) == [(-4, 5), (-4, 5), (-3, 4), (-3, 7), (-3, 7)]
say("  -> exactly (Q(sqrt-3), 4), (Q(i), 5) [both primes over 5], (Q(sqrt-3), 7) [both primes over 7]")
say("  small-norm cases excluded by N(P) >= 4: every field gives K2 at N = 2 and K3 at N = 3 (units +-1 suffice)")

# ====================================================================== B
say("")
say("=" * 100)
say("B. toroidal maps as genuine maps (Coxeter-Moser convention)")
say("=" * 100)


class TorusMap:
    """darts (v, k), k = direction index; sigma rotates counterclockwise; alpha reverses the dart."""

    def __init__(self, dirs, basis):
        (a, b), (c, d) = basis
        self.det = a * d - b * c
        self.basis = basis
        self.dirs = dirs
        self.n = len(dirs)
        Dt = self.det
        reps = set()
        R = abs(Dt) + 2
        for x in range(-R, R + 1):
            for y in range(-R, R + 1):
                reps.add(self.red((x, y)))
        self.V = sorted(reps)
        assert len(self.V) == abs(Dt)

    def red(self, v):
        (a, b), (c, d) = self.basis
        x, y = v
        D = self.det
        s = Fr(x * d - y * c, D)
        t = Fr(-x * b + y * a, D)
        s0, t0 = math.floor(s), math.floor(t)
        return (x - s0 * a - t0 * c, y - s0 * b - t0 * d)

    def darts(self):
        return [(v, k) for v in self.V for k in range(self.n)]

    def sigma(self, dart):
        v, k = dart
        return (v, (k + 1) % self.n)

    def alpha(self, dart):
        v, k = dart
        dx, dy = self.dirs[k]
        w = self.red((v[0] + dx, v[1] + dy))
        kk = (k + self.n // 2) % self.n           # opposite direction
        assert self.dirs[kk] == (-dx, -dy)
        return (w, kk)

    def faces(self):
        """faces = orbits of phi = sigma o alpha (go along dart, then turn)."""
        seen, faces = set(), []
        for d0 in self.darts():
            if d0 in seen:
                continue
            f, d = [], d0
            while d not in seen:
                seen.add(d)
                f.append(d)
                d = self.sigma(self.alpha(d))
            faces.append(f)
        return faces

    def skeleton(self):
        G = nx.MultiGraph()
        G.add_nodes_from(self.V)
        done = set()
        for d in self.darts():
            if d in done:
                continue
            e = self.alpha(d)
            done.add(d)
            done.add(e)
            G.add_edge(d[0], e[0])
        return G


def monodromy_order(darts, sig, alp):
    """size of the permutation group generated by sigma and alpha on darts (BFS on group elements)."""
    idx = {d: i for i, d in enumerate(darts)}
    s = tuple(idx[sig(d)] for d in darts)
    a = tuple(idx[alp(d)] for d in darts)
    ident = tuple(range(len(darts)))
    seen = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for h in (s, a):
                c = tuple(h[g[i]] for i in range(len(g)))
                if c not in seen:
                    seen.add(c)
                    nxt.append(c)
        frontier = nxt
        if len(seen) > 10 ** 6:
            raise RuntimeError("too big")
    return len(seen), s, a


def has_mirror_automorphism(s, a):
    """exists psi with psi s psi^-1 = s^-1 and psi a psi^-1 = a? (determined by the image of dart 0)."""
    n = len(s)
    sinv = [0] * n
    for i, j in enumerate(s):
        sinv[j] = i
    for target in range(n):
        psi = {0: target}
        stack = [0]
        ok = True
        while stack and ok:
            x = stack.pop()
            for (g, gimg) in ((s, sinv), (a, a)):
                y, yimg = g[x], gimg[psi[x]]
                if y in psi:
                    if psi[y] != yimg:
                        ok = False
                        break
                else:
                    psi[y] = yimg
                    stack.append(y)
        if ok and len(psi) == n and len(set(psi.values())) == n:
            return True
    return False


SQ_DIRS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
# triangular lattice in coordinates x + y*z6, z6 = e^{i pi/3}, z6^2 = z6 - 1; counterclockwise order of the 6 units
TRI_DIRS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]


def map44(b, c):
    return TorusMap(SQ_DIRS, ((b, c), (-c, b)))


def map36(b, c):
    # ideal generated by pi = b + c z6 ; pi*z6 = b z6 + c z6^2 = -c + (b + c) z6
    return TorusMap(TRI_DIRS, ((b, c), (-c, b + c)))


def dual_skeleton(M):
    faces = M.faces()
    fid = {}
    for i, f in enumerate(faces):
        for d in f:
            fid[d] = i
    G = nx.MultiGraph()
    G.add_nodes_from(range(len(faces)))
    done = set()
    for d in M.darts():
        if d in done:
            continue
        e = M.alpha(d)
        done.add(d)
        done.add(e)
        G.add_edge(fid[d], fid[e])
    return G, faces


def simple(G):
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    H.add_edges_from((u, v) for u, v in G.edges() if u != v)
    return H


def max_mult(G):
    cnt = {}
    for u, v in G.edges():
        key = frozenset((u, v))
        cnt[key] = cnt.get(key, 0) + 1
    return max(cnt.values()), any(u == v for u, v in G.edges())


heawood = nx.heawood_graph()
targets = {
    "K5": nx.complete_graph(5), "K4,4": nx.complete_bipartite_graph(4, 4), "K4": nx.complete_graph(4),
    "K7": nx.complete_graph(7), "K3,3": nx.complete_bipartite_graph(3, 3), "Heawood": heawood,
}
rows = []
for name, M, dual, tname in [
    ("{4,4}_(2,1)", map44(2, 1), False, "K5"),
    ("{4,4}_(2,2)", map44(2, 2), False, "K4,4"),
    ("{3,6}_(2,0)", map36(2, 0), False, "K4"),
    ("{3,6}_(2,1)", map36(2, 1), False, "K7"),
    ("{6,3}_(1,1)", map36(1, 1), True, "K3,3"),
    ("{6,3}_(2,1)", map36(2, 1), True, "Heawood"),
]:
    darts = M.darts()
    faces = M.faces()
    if not dual:
        G = M.skeleton()
        V, E, F = len(M.V), len(darts) // 2, len(faces)
        face_sizes = sorted(set(len(f) for f in faces))
        order, s, a = monodromy_order(darts, M.sigma, M.alpha)
    else:
        G, prim_faces = dual_skeleton(M)
        V, E, F = len(prim_faces), len(darts) // 2, len(M.V)
        face_sizes = [M.n]
        order, s, a = monodromy_order(darts, M.sigma, M.alpha)   # duality preserves the monodromy group
    mult, loops = max_mult(G)
    H = simple(G)
    iso = nx.is_isomorphic(H, targets[tname])
    planar, _ = nx.check_planarity(H)
    regular = (order == len(darts))
    chiral = regular and not has_mirror_automorphism(s, a)
    rows.append((name, V, E, F, V - E + F, face_sizes, mult, loops, tname, iso, planar, regular, chiral))
    say("  %-12s V=%-3d E=%-3d F=%-3d chi=%d faces%s max-mult=%d loops=%s  skeleton~%s: %s  planar=%s  "
        "orientably-regular=%s  chiral=%s"
        % (name, V, E, F, V - E + F, face_sizes, mult, loops, tname, iso, planar, regular, chiral))
    assert iso and V - E + F == 0 and regular
    assert planar == (tname == "K4")
chir = {r[0]: r[12] for r in rows}
assert chir == {"{4,4}_(2,1)": True, "{4,4}_(2,2)": False, "{3,6}_(2,0)": False, "{3,6}_(2,1)": True,
                "{6,3}_(1,1)": False, "{6,3}_(2,1)": True}
say("  -> all six identifications hold in Coxeter-Moser's convention (norm b^2+c^2, resp. b^2+bc+c^2);")
say("     {4,4}_(2,1) [K5], {3,6}_(2,1) [K7], {6,3}_(2,1) [Heawood] are CHIRAL (orientably regular, not reflexible)")
say("     {3,6}_(2,0): 4 vertices, 12 edges, every K4 edge doubled")

# girth forms of Euler's bound
say("  girth-Euler bound E <= g(V-2)/(g-2) for the simple skeletons:")
for nm, Gs in [("K5", nx.complete_graph(5)), ("K3,3", nx.complete_bipartite_graph(3, 3)), ("K7", nx.complete_graph(7)),
               ("Heawood", heawood), ("K4,4", nx.complete_bipartite_graph(4, 4)), ("Petersen", nx.petersen_graph()),
               ("K6", nx.complete_graph(6))]:
    g = nx.girth(Gs)
    V, E = Gs.number_of_nodes(), Gs.number_of_edges()
    bound = Fr(g * (V - 2), g - 2)
    say("    %-8s V=%-2d E=%-2d girth %d bound %-6s violated=%s  planar=%s"
        % (nm, V, E, g, str(bound), E > bound, nx.check_planarity(Gs)[0]))

# other small toroidal maps: which skeletons are simple and non-planar ("first non-planar" claim)
say("  small maps {4,4}_(b,c), {3,6}_(b,c), {6,3}_(b,c) (b >= c >= 0, b >= 1): simple? planar?")
for fam in ("{4,4}", "{3,6}", "{6,3}"):
    for b in range(1, 4):
        for c in range(0, b + 1):
            if fam == "{4,4}":
                M = map44(b, c)
                G = M.skeleton()
            elif fam == "{3,6}":
                M = map36(b, c)
                G = M.skeleton()
            else:
                M = map36(b, c)
                G, _ = dual_skeleton(M)
            mult, loops = max_mult(G)
            H = simple(G)
            say("    %s_(%d,%d): V=%-2d simple=%-5s planar(simple skeleton)=%s"
                % (fam, b, c, H.number_of_nodes(), mult == 1 and not loops, nx.check_planarity(H)[0]))

# ====================================================================== C
say("")
say("=" * 100)
say("C. Petersen = hemi-dodecahedron; map group")
say("=" * 100)
Dd = nx.dodecahedral_graph()
ok, emb = nx.check_planarity(Dd)
assert ok
faces = []
seenhe = set()
for u, v in emb.edges():
    if (u, v) in seenhe:
        continue
    f = emb.traverse_face(u, v, mark_half_edges=seenhe)
    faces.append(tuple(f))
assert len(faces) == 12 and all(len(f) == 5 for f in faces)
dist = dict(nx.all_pairs_shortest_path_length(Dd))
anti = {v: [w for w in Dd if dist[v][w] == 5] for v in Dd}
assert all(len(x) == 1 for x in anti.values())
anti = {v: x[0] for v, x in anti.items()}
assert all(anti[anti[v]] == v for v in Dd)
assert all(Dd.has_edge(anti[u], anti[v]) for u, v in Dd.edges())      # the antipodal map is a graph automorphism
cls = {v: min(v, anti[v]) for v in Dd}
Q = nx.Graph()
Q.add_edges_from((cls[u], cls[v]) for u, v in Dd.edges())
P = nx.petersen_graph()
iso_PQ = nx.is_isomorphic(Q, P)
qfaces = set(frozenset(cls[v] for v in f) for f in faces)
say("  dodecahedron: V=20 E=30 F=12; antipodal map = unique vertex at distance 5, an automorphism, maps faces to faces")
say("  quotient: V=%d E=%d, %d face pentagons, isomorphic to Petersen: %s; chi = %d (projective plane)"
    % (Q.number_of_nodes(), Q.number_of_edges(), len(qfaces), iso_PQ, Q.number_of_nodes() - Q.number_of_edges() + len(qfaces)))
assert iso_PQ and len(qfaces) == 6
GM = nx.algorithms.isomorphism.GraphMatcher(Q, Q)
auts = list(GM.isomorphisms_iter())
keep = [g for g in auts if set(frozenset(g[v] for v in f) for f in qfaces) == qfaces]
say("  |Aut(Petersen graph)| = %d; automorphisms preserving the 6 hemi-dodecahedral faces: %d" % (len(auts), len(keep)))
assert len(auts) == 120 and len(keep) == 60
# flags = (vertex, edge, face) incidences: 10 * 3 * 2 = 60; regular action?
flags = []
for f in qfaces:
    fl = list(f)
    for u, v in Q.edges():
        if u in f and v in f:
            # edge lies on the face boundary iff it is one of the 5 cycle edges
            sub = Q.subgraph(f)
            if sub.has_edge(u, v):
                flags.append((u, frozenset((u, v)), f))
                flags.append((v, frozenset((u, v)), f))
flags = set(flags)
orb = set()
f0 = next(iter(flags))
for g in keep:
    orb.add((g[f0[0]], frozenset(g[x] for x in f0[1]), frozenset(g[x] for x in f0[2])))
say("  flags: %d; orbit of one flag under the 60 map automorphisms: %d (regular action)" % (len(flags), len(orb)))
assert len(flags) == 60 and len(orb) == 60
# A_5 check: element orders 1:1, 2:15, 3:20, 5:24
def perm_order(g):
    k, h = 1, dict(g)
    while any(h[v] != v for v in h):
        h = {v: g[h[v]] for v in h}
        k += 1
    return k
ords = {}
for g in keep:
    o = perm_order(g)
    ords[o] = ords.get(o, 0) + 1
say("  element orders of the map group: %s  (A_5: {1: 1, 2: 15, 3: 20, 5: 24})" % dict(sorted(ords.items())))
assert ords == {1: 1, 2: 15, 3: 20, 5: 24}
say("  hemi-icosahedron (the other (2,3,5) quotient): K6 is non-planar too: planar(K6) = %s" % nx.check_planarity(nx.complete_graph(6))[0])

# ====================================================================== D
say("")
say("=" * 100)
say("D. the Klein-four Belyi map and X(4) -> X(2)")
say("=" * 100)
t, s_ = sp.symbols("t s")
lam = ((1 - t ** 2) / (1 + t ** 2)) ** 2
num, den = sp.fraction(sp.together(lam))
num1, den1 = sp.fraction(sp.together(lam - 1))
say("  lambda = %s / %s ; lambda - 1 = %s / %s" % (sp.factor(num), sp.factor(den), sp.factor(num1), sp.factor(den1)))
# behaviour at t = oo: substitute t = 1/s
lam_inf = sp.simplify(lam.subs(t, 1 / s_))
say("  at t = oo (t = 1/s): lambda = %s ; lambda - 1 = %s" % (sp.factor(lam_inf), sp.factor(sp.simplify(lam_inf - 1))))
fib0 = sp.roots(sp.Poly(num, t))
fib1 = sp.roots(sp.Poly(num1, t))
fibi = sp.roots(sp.Poly(den, t))
say("  fibre over 0: %s ; over 1: %s plus t = oo (order of vanishing of lambda-1 at s=0: %d) ; over oo: %s"
    % (fib0, fib1, sp.Poly(sp.fraction(sp.together(lam_inf - 1))[0], s_).as_expr().as_poly(s_).monoms()[-1][0], fibi))
assert fib0 == {1: 2, -1: 2} and fibi == {sp.I: 2, -sp.I: 2} and fib1 == {0: 2}
dl = sp.together(sp.diff(lam, t))
crit = sp.roots(sp.Poly(sp.fraction(dl)[0], t))
say("  finite critical points (zeros of lambda'): %s" % crit)
for g in (-t, 1 / t, -1 / t):
    assert sp.simplify(lam.subs(t, g) - lam) == 0
say("  deck group {t, -t, 1/t, -1/t} = Z/2 x Z/2 of order 4 = degree -> Galois; R-H: 2*0-2 = 4*(2*0-2) + 6 -> %s" % (-2 == 4 * (-2) + 6))
# numerical monodromy around 0 and 1, base point lambda = 1/2
mp.mp.dps = 30
lam_f = lambda z: ((1 - z * z) / (1 + z * z)) ** 2
dlam_f = lambda z: 2 * ((1 - z * z) / (1 + z * z)) * (-4 * z / (1 + z * z) ** 2)
base = mp.mpf(1) / 2
fibre = [complex(r) for r in sp.Poly(sp.fraction(sp.together(lam - sp.Rational(1, 2)))[0], t).nroots(n=30)]


def lift(z0, path):
    z = complex(z0)
    for k in range(1, len(path)):
        target = path[k]
        for _ in range(8):              # Newton to follow the branch
            z = z - (lam_f(z) - target) / dlam_f(z)
    return z


def loop(center, radius, n=3000):
    import cmath
    c, b0 = complex(center), complex(base)
    return [c + (b0 - c) * cmath.exp(2j * cmath.pi * k / n) for k in range(n + 1)]


perms = {}
for nm, cen in (("around 0", mp.mpf(0)), ("around 1", mp.mpf(1))):
    path = loop(cen, abs(base - cen))
    pm = []
    for z0 in fibre:
        z1 = lift(z0, path)
        j = min(range(4), key=lambda i: abs(fibre[i] - z1))
        pm.append(j)
    perms[nm] = pm
say("  numerical monodromy on the fibre over 1/2 (4 points): %s" % perms)
p0, p1 = perms["around 0"], perms["around 1"]
G4 = {tuple(range(4))}
fr = [tuple(range(4))]
while fr:
    nx_ = []
    for g in fr:
        for h in (p0, p1):
            c = tuple(h[g[i]] for i in range(4))
            if c not in G4:
                G4.add(c)
                nx_.append(c)
    fr = nx_
say("  monodromy group order %d, all non-identity elements involutions without fixed points: %s"
    % (len(G4), all(all(g[i] != i for i in range(4)) and all(g[g[i]] == i for i in range(4)) for g in G4 if g != tuple(range(4)))))
assert len(G4) == 4
# X(4) -> X(2): t(tau) = theta4^2/(theta3^2 + theta2^2), lambda(tau) = theta2^4/theta3^4


def thetas(tau):
    """theta_2, theta_3, theta_4 as functions of tau (direct sums; mpmath's jtheta(2, 0, q) uses the principal
    q^(1/4) and is not single-valued in tau, so it is not used here)."""
    tau = mp.mpc(tau)
    M = int(mp.ceil(mp.sqrt(mp.mpf(mp.mp.dps) * 2.4 / (mp.pi * mp.im(tau))))) + 5
    e = lambda x: mp.exp(1j * mp.pi * tau * x)
    th3 = 1 + 2 * mp.fsum(e(n * n) for n in range(1, M))
    th4 = 1 + 2 * mp.fsum((-1) ** n * e(n * n) for n in range(1, M))
    th2 = 2 * mp.fsum(e((n + mp.mpf(1) / 2) ** 2) for n in range(0, M))
    return th2, th3, th4


def t_of_tau(tau):
    th2, th3, th4 = thetas(tau)
    return th4 ** 2 / (th3 ** 2 + th2 ** 2)


def lam_of_tau(tau):
    th2, th3, th4 = thetas(tau)
    return (th2 / th3) ** 4


def act(g, tau):
    a, b, c, d = g
    return (a * tau + b) / (c * tau + d)


tau0 = mp.mpc("0.137", "1.21")
err_l = abs(lam_f(t_of_tau(tau0)) - lam_of_tau(tau0))
say("  theta check at tau = %s: |lambda(t(tau)) - lambda(tau)| = %.2e" % (tau0, float(err_l)))
assert err_l < 1e-20
mp.mp.dps = 30
gamma4 = [(1, 4, 0, 1), (1, 0, 4, 1), (5, -4, 4, -3), (-3, 4, -4, 5), (1, -4, 4, -15), (5, 4, 16, 13)]
gamma2_not4 = [(1, 2, 0, 1), (1, 0, 2, 1), (3, 2, 4, 3), (-1, 2, -2, 3)]
for g in gamma4 + gamma2_not4:
    a, b, c, d = g
    assert a * d - b * c == 1
    t1, t2 = t_of_tau(tau0), t_of_tau(act(g, tau0))
    rel = [abs(t2 - t1), abs(t2 + t1), abs(t2 - 1 / t1), abs(t2 + 1 / t1)]
    k = min(range(4), key=lambda i: rel[i])
    in4 = all(x % 4 == y for x, y in zip((a, b, c, d), (1, 0, 0, 1))) or all(x % 4 == y for x, y in zip((a, b, c, d), (3, 0, 0, 3)))
    say("    gamma = %-18s in +-Gamma(4): %-5s  t(gamma tau) = %s   (residual %.1e)"
        % (str(g), in4, ["t", "-t", "1/t", "-1/t"][k], float(rel[k])))
    assert rel[k] < 1e-15
    assert (k == 0) == in4
say("  -> t is a Hauptmodul for Gamma(4) and Gamma(2)/(+-Gamma(4)) acts as t -> +-t^(+-1): the Klein-four Belyi map")
say("     is the modular cover X(4) -> X(2) (also forced: F_2 = Gamma(2)-bar has a unique normal subgroup with quotient V_4)")

# ====================================================================== E
say("")
say("=" * 100)
say("E. Pythagorean Frey curves y^2 = x(x - a^2)(x + b^2): 4-torsion, full torsion, and the Kubert family")
say("=" * 100)


def ec_add(P1, P2, A2, A4):
    """y^2 = x^3 + A2 x^2 + A4 x (Fractions); None = O."""
    if P1 is None:
        return P2
    if P2 is None:
        return P1
    x1, y1 = P1
    x2, y2 = P2
    if x1 == x2 and y1 == -y2:
        return None
    if P1 == P2:
        l = (3 * x1 * x1 + 2 * A2 * x1 + A4) / (2 * y1)
    else:
        l = (y2 - y1) / (x2 - x1)
    x3 = l * l - A2 - x1 - x2
    return (x3, -(y1 + l * (x3 - x1)))


def is_sq(fr):
    fr = Fr(fr)
    if fr < 0:
        return False
    n, d = fr.numerator, fr.denominator
    return math.isqrt(n) ** 2 == n and math.isqrt(d) ** 2 == d


def mn_tree(depth):
    """Berggren tree in Euclid parameters: (m,n) -> (2m-n, m), (2m+n, m), (m+2n, n); root (2,1)."""
    level, out = [(2, 1)], [(2, 1)]
    for _ in range(depth):
        nxt = []
        for m, n in level:
            nxt += [(2 * m - n, m), (2 * m + n, m), (m + 2 * n, n)]
        out += nxt
        level = nxt
    return out


pairs = mn_tree(9)
trip = [(m * m - n * n, 2 * m * n, m * m + n * n) for m, n in pairs]
assert len(trip) == 29524 and len(set(trip)) == 29524
cnt_ord4 = 0
z8 = []
for (a, b, c) in trip[:3280]:            # depth <= 7
    for (A, B) in ((a, b), (b, a)):
        a2, b2 = Fr(A * A), Fr(B * B)
        A2, A4 = b2 - a2, -a2 * b2      # x(x - a2)(x + b2) = x^3 + (b2 - a2) x^2 - a2 b2 x
        Qp = (Fr(A * (A + c)), Fr(A * c * (A + c)))
        assert Qp[1] ** 2 == Qp[0] * (Qp[0] - a2) * (Qp[0] + b2)
        Q2 = ec_add(Qp, Qp, A2, A4)
        assert Q2 == (a2, 0) and ec_add(Q2, Q2, A2, A4) is None
        cnt_ord4 += 1
        # Z/8? a point R of order 4 is halvable iff x(R) - e_i are all squares (e_i = 0, a2, -b2)
        for T in (None, (Fr(0), Fr(0)), (-b2, Fr(0))):
            R = ec_add(Qp, T, A2, A4)
            x0 = R[0]
            if all(is_sq(x0 - e) for e in (0, a2, -b2)):
                z8.append((A, B, c))
say("  order-4 point (a(a+c), ac(a+c)) with 2Q = (a^2, 0): verified for %d (triple, leg order) pairs (depth <= 7)" % cnt_ord4)
say("  pairs whose curve has a halvable point of order 4 (torsion Z/2 x Z/8): %d %s" % (len(z8), z8[:5]))
# point counts for a few curves to confirm torsion order 8 (or 16)
def count_points(a2, b2, p):
    n = 1
    for x in range(p):
        r = (x * (x - a2) * (x + b2)) % p
        if r == 0:
            n += 1
        elif pow(r, (p - 1) // 2, p) == 1:
            n += 2
    return n
for (a, b, c) in trip[:6]:
    g = 0
    for p in sp.primerange(3, 400):
        if (a * b * c) % p == 0:
            continue
        g = math.gcd(g, count_points(a * a % p, b * b % p, p))
    say("    (%d,%d,%d): gcd of #E(F_p) over good p < 400 = %d (torsion order divides it)" % (a, b, c, g))
# the Z/2 x Z/4 family is strictly larger: y^2 = x(x+1)(x+4)
e = [0, -1, -4]
say("  y^2 = x(x+1)(x+4): (0,0) halvable (0-(-1)=1, 0-(-4)=4 squares) -> Z/2 x Z/4 subset of E(Q);")
lams = set()
for i, j, k in itertools.permutations(range(3)):
    lams.add(Fr(e[k] - e[i], e[j] - e[i]))
pyth = [l for l in lams if is_sq(l) and is_sq(1 - l)]
say("    its six lambda-values %s; values with lambda and 1-lambda both squares: %s -> NOT a Pythagorean Frey curve"
    % (sorted(lams), pyth))
assert not pyth
say("  So the Pythagorean Frey curves are the thin subfamily lambda, 1-lambda both squares = image of X(4)(Q),")
say("  i.e. a degree-2 image inside the Z/2 x Z/4 family (X_1(2,4) = X(Gamma(2) cap Gamma_1(4))); over Q(i) they have full 4-torsion.")

# ====================================================================== F
say("")
say("=" * 100)
say("F. Berggren abc hits (a^2, b^2, c^2), depth <= 9 (independent generator)")
say("=" * 100)
hitsF = []
for (a, b, c) in trip:
    R = rad(a * b * c)
    if R < c * c:
        hitsF.append((math.log(c * c) / math.log(R), a, b, c))
hitsF.sort(reverse=True)
say("  triples %d, hits %d" % (len(trip), len(hitsF)))
for q, a, b, c in hitsF[:8]:
    say("    q = %.6f  (%d, %d, %d)  rad = %d" % (q, a, b, c, rad(a * b * c)))
assert len(hitsF) == 60
assert (hitsF[0][1:] in ((36207, 18424, 40625), (18424, 36207, 40625)))
adu = [h for h in hitsF if sorted(h[1:3]) == [175, 288]]
say("  (175, 288, 337) present: %s (q = %.4f)" % (bool(adu), adu[0][0] if adu else float("nan")))
seven = [h for h in hitsF if sorted(h[1:3]) == [7, 24]]
say("  (7, 24, 25): q = %.4f" % seven[0][0])
# completeness for c <= 2000
eu = set()
for m in range(2, 50):
    for n in range(1, m):
        if (m - n) % 2 and math.gcd(m, n) == 1 and m * m + n * n <= 2000:
            eu.add((m * m + n * n, frozenset((m * m - n * n, 2 * m * n))))
say("  primitive triples with c <= 2000: %d" % len(eu))
assert len(eu) == 319

# ====================================================================== G
say("")
say("=" * 100)
say("G. Collatz clocks = convergents and semiconvergents K/L of log2 3")
say("=" * 100)
mp.mp.dps = 80
theta = mp.log(3) / mp.log(2)
cf = []
x = theta
for _ in range(14):
    a = int(mp.floor(x))
    cf.append(a)
    x = 1 / (x - a)
say("  log2 3 = [%s; %s]" % (cf[0], ", ".join(map(str, cf[1:]))))
convs = []
h1, h0, k1, k0 = 1, 0, 0, 1
for a in cf:
    h1, h0 = a * h1 + h0, h1
    k1, k0 = a * k1 + k0, k1
    convs.append((h1, k1))
clocks = set()
for i in range(len(convs)):
    clocks.add(convs[i])
    if i + 1 < len(cf):
        hm, km = convs[i - 1] if i >= 1 else (1, 0)
        for j in range(1, cf[i + 1]):
            clocks.add((hm + j * convs[i][0], km + j * convs[i][1]))
clocks = sorted(c for c in clocks if c[1] >= 1)
small = [c for c in clocks if c[1] <= 53]
say("  clocks with L <= 53: %s" % small)
assert small == [(1, 1), (2, 1), (3, 2), (5, 3), (8, 5), (11, 7), (19, 12), (27, 17), (46, 29), (65, 41), (84, 53)]
for K, L in clocks:
    if L > 20000:
        continue
    hi, lo = max(2 ** K, 3 ** L), min(2 ** K, 3 ** L)
    d = hi - lo
    sure = 6 * d < hi                    # then 6 rad(d) <= 6 d < hi: an abc hit without factoring
    if L <= 53:
        R = 6 * rad(d)
        q = math.log(hi) / math.log(R)
        say("    %3d/%-3d d = %-26d d < max/6: %-5s rad = %-26d q = %.4f %s"
            % (K, L, d, sure, R, q, "HIT" if R < hi else "no"))
        assert (R < hi) == (L >= 2)
    else:
        assert sure, (K, L)
say("  every clock with 2 <= L <= 20000 (%d clocks) satisfies |2^K - 3^L| < max/6, hence is an abc hit"
    % len([c for c in clocks if 2 <= c[1] <= 20000]))

# ====================================================================== H
say("")
say("=" * 100)
say("H. hyperbolic signatures and the S116 inequality")
say("=" * 100)
best, arg = Fr(0), []
for p in range(2, 60):
    for q in range(p, 60):
        for r in range(q, 400):
            sgm = Fr(1, p) + Fr(1, q) + Fr(1, r)
            if sgm < 1:
                if sgm > best:
                    best, arg = sgm, [(p, q, r)]
                elif sgm == best:
                    arg.append((p, q, r))
# beyond the box: if p >= 3 then sigma <= 1/3+1/3+1/4 = 11/12 or below; if p = 2, q >= 60 then sigma < 1/2+2/60; r >= 400 adds < 1/400
say("  max of 1/p+1/q+1/r < 1 over 2 <= p <= q <= r (p,q < 60, r < 400; the tail is trivially smaller): %s at %s" % (best, arg))
assert best == Fr(41, 42) and arg == [(2, 3, 7)]
eps_star = mp.log(9) / mp.log(6) - 1
say("  1 + 8 = 9, rad = 6: rad^{1+eps} > 9 fails exactly for eps <= log 9/log 6 - 1 = %s" % mp.nstr(eps_star, 8))
say("")
say("ALL ASSERTIONS PASSED (%.1f s)" % (time.time() - T0))
