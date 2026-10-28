"""opus-2026-10-28-S68: concrete bridges into the modular / triangle-group world (THM-4839).

  A  complete unit Cayley graphs Cay(O_K/P, O_K^x) of imaginary quadratic residue fields: exactly K4, K5, K7 (N(P) >= 4)
  B  toroidal regular maps {4,4}_(b,c), {3,6}_(b,c), {6,3}_(b,c) from Gaussian / Eisenstein ideals:
     K5 = {4,4}_(2,1), K4 = {3,6}_(2,0), K7 = {3,6}_(2,1), K3,3 = {6,3}_(1,1), Heawood = {6,3}_(2,1), K4,4 = {4,4}_(2,2)
  C  Petersen = hemi-dodecahedron (dodecahedron / antipodal map), checked by explicit isomorphism with Kneser K(5,2)
  D  the Klein-four Belyi map lambda(t) = ((1-t^2)/(1+t^2))^2 (Euclid's Pythagorean parametrisation, squared)
  E  Pythagorean Frey curves y^2 = x(x - a^2)(x + b^2): the rational point (a(a+c), ac(a+c)) has order 4
  F  abc qualities of the Berggren triples (a^2, b^2, c^2)
  G  the Collatz clocks K/L (convergents and semiconvergents of log2 3) as abc triples (min, |2^K - 3^L|, max)
Every claim is an assert.
"""
import itertools
import math
import time
from fractions import Fraction as Fr

import sympy as sp
from sympy import factorint, isprime

T0 = time.time()


def say(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- graph utilities
def norm_graph(V, E):
    adj = {v: set() for v in V}
    for u, v in E:
        if u == v:
            raise ValueError("loop")
        adj[u].add(v)
        adj[v].add(u)
    return adj


def simple_edges(V, pairs, allow_multi=False):
    """return (adj, multiplicity) from a list of directed neighbour pairs (both directions present);
    multiplicity = max edge multiplicity (1 = simple graph); loops are rejected."""
    multi = {}
    for u, v in pairs:
        if u == v:
            return None, 0
        key = (u, v)
        multi[key] = multi.get(key, 0) + 1
    mult = max(multi.values())
    if mult > 1 and not allow_multi:
        return None, mult
    adj = {v: set() for v in V}
    for (u, v) in multi:
        adj[u].add(v)
    for u in adj:
        for v in adj[u]:
            assert u in adj[v]
    return adj, mult


def girth(adj):
    best = math.inf
    for s in adj:
        dist = {s: 0}
        par = {s: None}
        q = [s]
        for u in q:
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    par[w] = u
                    q.append(w)
                elif par[u] != w:
                    best = min(best, dist[u] + dist[w] + 1)
    return best


def bipartite(adj):
    col = {}
    for s in adj:
        if s in col:
            continue
        col[s] = 0
        q = [s]
        for u in q:
            for w in adj[u]:
                if w not in col:
                    col[w] = 1 - col[u]
                    q.append(w)
                elif col[w] == col[u]:
                    return False
    return True


def isomorphic(A, B):
    """backtracking isomorphism test for small graphs (degree-refined)."""
    if len(A) != len(B):
        return False
    if sorted(len(A[v]) for v in A) != sorted(len(B[v]) for v in B):
        return False
    av = sorted(A, key=lambda v: -len(A[v]))
    bv = list(B)
    mp, used = {}, set()

    def ok(a, b):
        if len(A[a]) != len(B[b]):
            return False
        for x, y in mp.items():
            if (x in A[a]) != (y in B[b]):
                return False
        return True

    def rec(i):
        if i == len(av):
            return True
        a = av[i]
        for b in bv:
            if b in used or not ok(a, b):
                continue
            mp[a] = b
            used.add(b)
            if rec(i + 1):
                return True
            del mp[a]
            used.discard(b)
        return False

    return rec(0)


def complete(n):
    return {i: set(j for j in range(n) if j != i) for i in range(n)}


def kmn(m, n):
    adj = {("a", i): set(("b", j) for j in range(n)) for i in range(m)}
    adj.update({("b", j): set(("a", i) for i in range(m)) for j in range(n)})
    return adj


def kneser52():
    V = list(itertools.combinations(range(5), 2))
    return {v: set(w for w in V if not set(v) & set(w)) for v in V}


def heawood():
    lines = [frozenset({(i) % 7, (i + 1) % 7, (i + 3) % 7}) for i in range(7)]
    adj = {("p", i): set() for i in range(7)}
    adj.update({("l", j): set() for j in range(7)})
    for j, L in enumerate(lines):
        for i in L:
            adj[("p", i)].add(("l", j))
            adj[("l", j)].add(("p", i))
    return adj


# ---------------------------------------------------------------- A
say("=" * 78)
say("A. complete unit Cayley graphs of imaginary quadratic residue fields")
say("=" * 78)
# class-number-one fields are enough for the units; the argument in the note covers every imaginary quadratic field.
# residue fields are built explicitly: split/ramified p -> F_p via a root, inert p -> F_p[x]/(minimal polynomial).
FIELDS = {1: "Q(i)", 3: "Q(sqrt-3)", 2: "Q(sqrt-2)", 7: "Q(sqrt-7)", 11: "Q(sqrt-11)", 19: "Q(sqrt-19)", 43: "Q(sqrt-43)"}


def residue_images(d, p):
    """yield (N, set of unit images as tuples, field elements count) for each prime above p in O_K, K = Q(sqrt-d)."""
    # O_K = Z[theta], theta = sqrt(-d) (d = 1, 2) or (1 + sqrt(-d))/2 (d = 3 mod 4)
    if d in (1, 2):
        f = [d, 0, 1]                       # x^2 + d
        units_poly = {1: [[1], [p - 1], [0, 1], [0, p - 1]], 2: [[1], [p - 1]]}[d]
    else:
        f = [(1 + d) // 4, -1, 1]           # x^2 - x + (1+d)/4
        if d == 3:   # units +-1, +-theta, +-(theta - 1)   (theta = -omega^2 = 1 + omega)
            units_poly = [[1], [p - 1], [0, 1], [0, p - 1], [p - 1, 1], [1, p - 1]]
        else:
            units_poly = [[1], [p - 1]]
    # roots of f mod p
    roots = [r for r in range(p) if (f[0] + f[1] * r + f[2] * r * r) % p == 0]
    out = []
    if roots:
        for r in sorted(set(roots)):
            imgs = set()
            for u in units_poly:
                val = sum(c * pow(r, k, p) for k, c in enumerate(u)) % p
                imgs.add((val,))
            out.append((p, imgs))
    else:   # inert: F_{p^2} = F_p[x]/(f); elements as (c0, c1)
        imgs = set()
        for u in units_poly:
            c = [x % p for x in (u + [0])][:2]
            imgs.add(tuple(c))
        out.append((p * p, imgs))
    return out


hits = []
for d, name in FIELDS.items():
    for p in [q for q in range(2, 60) if isprime(q)]:
        for N, imgs in residue_images(d, p):
            nonzero = N - 1
            imgs_nz = set(x for x in imgs if any(x))
            if N >= 4 and len(imgs_nz) == nonzero:
                hits.append((name, N))
hits = sorted(set(hits), key=lambda t: t[1])
say("  complete unit Cayley graphs with N(P) >= 4 (fields d in %s, p < 60): %s" % (sorted(FIELDS), hits))
assert hits == [("Q(sqrt-3)", 4), ("Q(i)", 5), ("Q(sqrt-3)", 7)]

# ---------------------------------------------------------------- B
say("")
say("=" * 78)
say("B. toroidal regular maps from Gaussian and Eisenstein ideals")
say("=" * 78)


def lattice_quot(gens, N):
    """Z^2 modulo the sublattice spanned by the two integer vectors gens (index N): canonical reps via Smith-like reduction."""
    (a, b), (c, d) = gens
    det = abs(a * d - b * c)
    assert det == N
    # reduce a vector modulo the lattice using the inverse matrix (exact rationals)
    def red(v):
        x, y = v
        # solve s*(a,b) + t*(c,d) = (x,y)
        D = a * d - b * c
        s = Fr(x * d - y * c, D)
        t = Fr(-x * b + y * a, D)
        s0, t0 = math.floor(s), math.floor(t)
        return (x - s0 * a - t0 * c, y - s0 * b - t0 * d)
    reps = set()
    for x in range(-2 * N, 2 * N + 1):
        for y in range(-2 * N, 2 * N + 1):
            reps.add(red((x, y)))
    assert len(reps) == N, (len(reps), N)
    return red, sorted(reps)


def map44(b, c):
    # Gaussian ideal (b + c i): lattice spanned by (b, c) and (-c, b) in coordinates (re, im)
    red, V = lattice_quot(((b, c), (-c, b)), b * b + c * c)
    pairs = [(v, red((v[0] + dx, v[1] + dy))) for v in V for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))]
    return simple_edges(V, pairs)


def eis(x, y):
    """Eisenstein integer x + y*omega as a lattice vector (x, y); multiplication by omega: (x,y) -> (-y, x - y)."""
    return (x, y)


def emul_omega(v):
    x, y = v
    return (-y, x - y)


def ideal_lattice(pi):
    return (pi, emul_omega(pi))


def map36(b, c):
    # Eisenstein ideal pi = b + c(1 + omega) = (b + c) + c omega; triangular lattice, 6 unit neighbours
    pi = (b + c, c)
    N = b * b + b * c + c * c
    red, V = lattice_quot(ideal_lattice(pi), N)
    units = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]      # +-1, +-omega, +-(1+omega)
    pairs = [(v, red((v[0] + u[0], v[1] + u[1]))) for v in V for u in units]
    return simple_edges(V, pairs, allow_multi=True)


def map63(b, c):
    # honeycomb: A_x ~ B_x, B_(x+1), B_(x+1+omega)  (offsets 0, 1, 1+omega: a unit triangle), modulo pi
    pi = (b + c, c)
    N = b * b + b * c + c * c
    red, V = lattice_quot(ideal_lattice(pi), N)
    pairs = []
    for x in V:
        for off in ((0, 0), (1, 0), (1, 1)):
            y = red((x[0] + off[0], x[1] + off[1]))
            pairs.append((("A", x), ("B", y)))
            pairs.append((("B", y), ("A", x)))
    return simple_edges([("A", x) for x in V] + [("B", x) for x in V], pairs)


checks = [
    ("{4,4}_(2,1)", map44(2, 1), complete(5), "K5"),
    ("{4,4}_(2,2)", map44(2, 2), kmn(4, 4), "K4,4"),
    ("{3,6}_(2,0)", map36(2, 0), complete(4), "K4"),
    ("{3,6}_(2,1)", map36(2, 1), complete(7), "K7"),
    ("{6,3}_(1,1)", map63(1, 1), kmn(3, 3), "K3,3"),
    ("{6,3}_(2,1)", map63(2, 1), heawood(), "Heawood"),
]
for name, (adj, mult), target, tname in checks:
    assert mult == (2 if name == "{3,6}_(2,0)" else 1), (name, mult)
    E = sum(len(s) for s in adj.values()) // 2
    V = len(adj)
    g = girth(adj)
    planar_bound = (g * (V - 2)) // (g - 2)          # a planar graph of girth g has E <= g(V-2)/(g-2)
    iso = isomorphic(adj, target)
    say("  %-12s V = %-3d simple E = %-3d edge multiplicity %d girth %s bipartite %-5s girth-Euler bound %-3d -> %s: %s"
        % (name, V, E, mult, girth(adj), bipartite(adj), planar_bound, tname, iso))
    assert iso
    if tname in ("K5", "K3,3", "K7", "Heawood", "K4,4"):
        assert E > planar_bound          # non-planar by the girth form of Euler's bound
assert girth(heawood()) == 6

# ---------------------------------------------------------------- C
say("")
say("=" * 78)
say("C. Petersen = hemi-dodecahedron")
say("=" * 78)
phi = sp.GoldenRatio
pts = []
for s in itertools.product((1, -1), repeat=3):
    pts.append(tuple(sp.Integer(x) for x in s))
for s1, s2 in itertools.product((1, -1), repeat=2):
    pts.append((sp.Integer(0), s1 / phi, s2 * phi))
    pts.append((s1 / phi, s2 * phi, sp.Integer(0)))
    pts.append((s2 * phi, sp.Integer(0), s1 / phi))
assert len(pts) == 20
d2 = lambda u, v: sp.nsimplify(sp.expand(sum((x - y) ** 2 for x, y in zip(u, v))))
edge2 = sp.nsimplify(sp.expand((2 / phi) ** 2))
dod = {i: set() for i in range(20)}
for i, j in itertools.combinations(range(20), 2):
    if sp.simplify(d2(pts[i], pts[j]) - edge2) == 0:
        dod[i].add(j)
        dod[j].add(i)
assert all(len(s) == 3 for s in dod.values()) and sum(len(s) for s in dod.values()) == 60
anti = {}
for i in range(20):
    for j in range(20):
        if all(sp.simplify(a + b) == 0 for a, b in zip(pts[i], pts[j])):
            anti[i] = j
cls = {}
for i in range(20):
    cls[i] = min(i, anti[i])
Q = {c: set() for c in set(cls.values())}
for i in range(20):
    for j in dod[i]:
        assert cls[i] != cls[j]
        Q[cls[i]].add(cls[j])
Vq, Eq = len(Q), sum(len(s) for s in Q.values()) // 2
say("  dodecahedron 20 vertices, 30 edges; antipodal quotient: V = %d, E = %d, cubic %s, girth %d, Petersen: %s"
    % (Vq, Eq, all(len(s) == 3 for s in Q.values()), girth(Q), isomorphic(Q, kneser52())))
assert Vq == 10 and Eq == 15 and girth(Q) == 5 and isomorphic(Q, kneser52())

# ---------------------------------------------------------------- D
say("")
say("=" * 78)
say("D. Euclid's parametrisation squared is a Klein-four Belyi map")
say("=" * 78)
t = sp.symbols("t")
lam = ((1 - t ** 2) / (1 + t ** 2)) ** 2
for g in (-t, 1 / t, -1 / t):
    assert sp.simplify(lam.subs(t, g) - lam) == 0
num, den = sp.fraction(sp.together(sp.diff(lam, t)))
crit = sp.solve(sp.Eq(num, 0), t)
cvals = set(sp.simplify(lam.subs(t, c)) for c in crit)
assert cvals <= {0, 1}
assert sp.limit(lam, t, sp.oo) == 1          # t = infinity is a critical point over 1
one_minus = sp.factor(sp.simplify(1 - lam))
say("  lambda(t) = ((1-t^2)/(1+t^2))^2, invariant under t -> -t, 1/t, -1/t (Klein four-group)")
say("  finite critical points %s with critical values %s; 1 - lambda = %s" % (crit, cvals, one_minus))
say("  fibres: over 0 -> t = +-1 (e=2), over 1 -> t = 0, oo (e=2), over oo -> t = +-i (e=2): passport (2^2, 2^2, 2^2), genus 0")
assert sp.simplify(one_minus - 4 * t ** 2 / (1 + t ** 2) ** 2) == 0
# Riemann-Hurwitz: 2*0 - 2 = 4*(-2) + 6*(2-1)
assert 2 * 0 - 2 == 4 * (-2) + 6
# Euclid: t = n/m gives (a,b,c) = (m^2 - n^2, 2mn, m^2 + n^2) and lambda = (a/c)^2, 1 - lambda = (b/c)^2
for m_, n_ in [(2, 1), (3, 2), (4, 1), (16, 9)]:
    a_, b_, c_ = m_ * m_ - n_ * n_, 2 * m_ * n_, m_ * m_ + n_ * n_
    assert lam.subs(t, sp.Rational(n_, m_)) == sp.Rational(a_, c_) ** 2

# ---------------------------------------------------------------- E
say("")
say("=" * 78)
say("E. Berggren tree and Pythagorean Frey curves: a rational point of order 4")
say("=" * 78)
A = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
B = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
C = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))


def mv(M, v):
    return tuple(sum(M[i][j] * v[j] for j in range(3)) for i in range(3))


def berggren(depth):
    level = [(3, 4, 5)]
    out = [(3, 4, 5)]
    for _ in range(depth):
        nxt = []
        for v in level:
            for M in (A, B, C):
                w = mv(M, v)
                w = (abs(w[0]), abs(w[1]), abs(w[2]))
                nxt.append(w)
        out.extend(nxt)
        level = nxt
    return out


trip = berggren(9)
assert len(trip) == sum(3 ** k for k in range(10))
# completeness / uniqueness against Euclid's enumeration for c <= 2000
euclid = set()
for m_ in range(2, 50):
    for n_ in range(1, m_):
        if (m_ - n_) % 2 == 1 and math.gcd(m_, n_) == 1 and m_ * m_ + n_ * n_ <= 2000:
            a_, b_ = m_ * m_ - n_ * n_, 2 * m_ * n_
            euclid.add((min(a_, b_), max(a_, b_), m_ * m_ + n_ * n_))
def berggren_upto(cmax):
    out, stack = [], [(3, 4, 5)]
    while stack:
        v = stack.pop()
        out.append(v)
        for M in (A, B, C):
            w = mv(M, v)
            w = (abs(w[0]), abs(w[1]), abs(w[2]))
            if w[2] <= cmax:                 # hypotenuses strictly increase down the tree
                stack.append(w)
    return out


tree_small = [tuple(sorted(x[:2])) + (x[2],) for x in berggren_upto(2000)]
assert sorted(tree_small) == sorted(euclid) and len(tree_small) == len(set(tree_small))
say("  Berggren tree to depth 9: %d triples; every primitive triple with c <= 2000 appears exactly once (%d)" % (len(trip), len(euclid)))


def add(P, Q_, a2, b2):
    """group law on y^2 = x(x - a2)(x + b2) = x^3 + (b2 - a2) x^2 - a2 b2 x; points as (x, y) Fractions or None."""
    A2, A4 = b2 - a2, -a2 * b2
    if P is None:
        return Q_
    if Q_ is None:
        return P
    x1, y1 = P
    x2, y2 = Q_
    if x1 == x2 and y1 == -y2:
        return None
    if P == Q_:
        lmb = (3 * x1 * x1 + 2 * A2 * x1 + A4) / (2 * y1)
    else:
        lmb = (y2 - y1) / (x2 - x1)
    x3 = lmb * lmb - A2 - x1 - x2
    y3 = -(y1 + lmb * (x3 - x1))
    return (x3, y3)


count = 0
for (a_, b_, c_) in trip[:3000]:
    for (a1, b1) in ((a_, b_), (b_, a_)):
        a2, b2 = Fr(a1 * a1), Fr(b1 * b1)
        Qp = (Fr(a1 * (a1 + c_)), Fr(a1 * c_ * (a1 + c_)))
        x, y = Qp
        assert y * y == x * (x - a2) * (x + b2)
        twoQ = add(Qp, Qp, a2, b2)
        assert twoQ == (a2, Fr(0))
        assert add(twoQ, twoQ, a2, b2) is None
        count += 1
say("  for %d (triple, leg-order) pairs: Q = (a(a+c), ac(a+c)) lies on y^2 = x(x - a^2)(x + b^2), 2Q = (a^2, 0), 4Q = O" % count)

# ---------------------------------------------------------------- F
say("")
say("=" * 78)
say("F. abc qualities of the Berggren triples (a^2, b^2, c^2)")
say("=" * 78)


def rad(n):
    r = 1
    for p in factorint(n):
        r *= p
    return r


hitsF = []
for (a_, b_, c_) in trip:
    R = rad(a_) * rad(b_) * rad(c_)          # legs and hypotenuse are pairwise coprime
    if R < c_ * c_:
        hitsF.append((math.log(c_ * c_) / math.log(R), a_, b_, c_))
hitsF.sort(reverse=True)
say("  depth <= 9: %d triples, %d abc hits (rad(abc) < c^2); best qualities:" % (len(trip), len(hitsF)))
for q, a_, b_, c_ in hitsF[:6]:
    say("    q = %.4f  (%d, %d, %d)" % (q, a_, b_, c_))
assert any((a_, b_, c_) in ((175, 288, 337), (288, 175, 337)) for _, a_, b_, c_ in hitsF)

# ---------------------------------------------------------------- G
say("")
say("=" * 78)
say("G. Collatz clocks K/L (convergents and semiconvergents of log2 3) as abc triples")
say("=" * 78)
clocks = [(1, 1), (2, 1), (3, 2), (5, 3), (8, 5), (11, 7), (19, 12), (27, 17), (46, 29), (65, 41), (84, 53)]
rows = []
for K, L in clocks:
    x, y = 2 ** K, 3 ** L
    d = abs(x - y)
    lo, hi = min(x, y), max(x, y)
    if d == 0:
        continue
    t1 = time.time()
    R = 6 * rad(d) if d > 1 else 6
    q = math.log(hi) / math.log(R)
    rows.append((K, L, d, q))
    say("  %3d/%-3d |2^K - 3^L| = %-28d rad = 6*rad(d) = %-30d q = %.4f %s (%.1fs)"
        % (K, L, d, R, q, "HIT" if q > 1 else "", time.time() - t1))
assert [r[:2] for r in rows if r[3] > 1] == clocks[2:]       # every clock from 3/2 on is an abc hit
say("  every clock from 3/2 on is a hit: |2^K - 3^L| < max/6 forces q > 1; Baker's bound sends q -> 1 along convergents")
say("")
say("ALL ASSERTIONS PASSED (%.1fs)" % (time.time() - T0))
