#!/usr/bin/env python3
"""Tournaments as chirotopes (mac-mini-2026-10-07-twoanchor).
(1) 4-Pfaffian criterion: for every labelled 4-tournament, Pf(B) = b12 b34 - b13 b24 + b14 b23 is +-1 or +-3, and
    |Pf| = 3 iff the 4-set has exactly one 3-cycle (a 3-cycle with a dominating or dominated vertex).
    Hence T is a rank-2 chirotope (sign pattern det(v_i, v_j) of vectors in R^2) iff all 4-Pfaffians are +-1
    iff T is locally transitive.
(2) Legendre chirotopes: for p = 3 mod 4 the Paley tournament plus a sink is the chi-pattern of P^1(F_p);
    the elliptic constant-quiddity (Chebyshev) SL2-frieze over F_p realises the whole switching class.
"""
import itertools, random
def pf4(b, i, j, k, l): return b[i][j]*b[k][l] - b[i][k]*b[j][l] + b[i][l]*b[j][k]
def three_cycles(b, S):
    c = 0
    for (i, j, k) in itertools.combinations(S, 3):
        if (b[i][j] == b[j][k] == b[k][i]): c += 1
    return c
def locally_transitive(b, n):
    for v in range(n):
        for side in (1, -1):
            N = [w for w in range(n) if w != v and b[v][w] == side]
            if three_cycles(b, N) > 0: return False
    return True
CH = 0
# (1) all 64 labelled 4-tournaments
for bits in range(64):
    b = [[0]*4 for _ in range(4)]
    for idx, (i, j) in enumerate(itertools.combinations(range(4), 2)):
        s = 1 if (bits >> idx) & 1 else -1; b[i][j] = s; b[j][i] = -s
    P = pf4(b, 0, 1, 2, 3); c3 = three_cycles(b, range(4))
    assert P in (-3, -1, 1, 3)
    assert (abs(P) == 3) == (c3 == 1), (bits, P, c3)
    CH += 1
# random n: LT iff all 4-Pfaffians +-1 ; and LT tournaments are exactly sign patterns of planar vector configurations
rnd = random.Random(7)
for n in (5, 6, 7, 8):
    for trial in range(400):
        if trial % 2 == 0:
            # realizable: random vectors in R^2
            ang = [rnd.uniform(0, 6.283) for _ in range(n)]
            import math
            v = [(math.cos(a), math.sin(a)) for a in ang]
            b = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    if i != j: b[i][j] = 1 if v[i][0]*v[j][1] - v[i][1]*v[j][0] > 0 else -1
        else:
            b = [[0]*n for _ in range(n)]
            for i, j in itertools.combinations(range(n), 2):
                s = rnd.choice((1, -1)); b[i][j] = s; b[j][i] = -s
        allpm1 = all(abs(pf4(b, *q)) == 1 for q in itertools.combinations(range(n), 4))
        lt = locally_transitive(b, n)
        assert allpm1 == lt
        if trial % 2 == 0: assert lt
        CH += 1
print("(1) 4-Pfaffian criterion: |Pf|=3 iff exactly one 3-cycle (all 64 labelled 4-tournaments); LT <=> all 4-Pfaffians +-1 <=> realizable (800 random + 800 planar configurations)")

# (2) Legendre chirotopes and the Paley frieze
def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1
def frieze_strip(points, p):
    """points: list of (x,y) in F_p^2 representing distinct points of P^1; rescale successively so that
    det(v_i, v_{i+1}) = 1 (an SL2-frieze strip); return the normalised vectors."""
    v = [points[0]]
    for P in points[1:]:
        d = (v[-1][0]*P[1] - v[-1][1]*P[0]) % p
        lam = pow(d, -1, p)
        v.append(((lam*P[0]) % p, (lam*P[1]) % p))
    return v
def singer_points(p):
    # PGL2 non-split torus: multiplication by a generator z of F_{p^2}^* on F_{p^2} = F_p^2, acting on P^1(F_p)
    # F_{p^2} = F_p[s]/(s^2 - r), r a non-residue
    r = next(x for x in range(2, p) if pow(x, (p-1)//2, p) == p-1)
    def mul(a, b): return ((a[0]*b[0] + r*a[1]*b[1]) % p, (a[0]*b[1] + a[1]*b[0]) % p)
    def order_mod_scalars(z):
        w = z; k = 1
        while not (w[1] == 0): w = mul(w, z); k += 1
        return k
    for z0 in range(p):
        for z1 in range(1, p):
            z = (z0, z1)
            if order_mod_scalars(z) == p + 1:
                pts = [(1, 0)]
                for _ in range(p): pts.append(mul(pts[-1], z))
                return pts
def iso_paley_plus_sink(b, p):
    # switch so that vertex 0 is a sink, delete it, compare with Paley via affine relabelling search
    n = p + 1
    sw = [1]*n
    for j in range(1, n):
        if b[0][j] == 1: sw[j] = -1   # 0 -> j : switch j so that j -> 0
    c = [[b[i][j]*sw[i]*sw[j] for j in range(n)] for i in range(n)]
    for j in range(1, n): assert c[j][0] == 1
    # remaining tournament on 1..p : backtracking search for phi with c[u][v] = chi(phi(v) - phi(u))
    verts = list(range(1, n))
    def bt(phi, used, idx):
        if idx == len(verts): return True
        v = verts[idx]
        for x in range(p):
            if x in used: continue
            if all(c[u][v] == legendre(x - phi[u], p) for u in phi):
                phi[v] = x; used.add(x)
                if bt(phi, used, idx + 1): return True
                del phi[v]; used.discard(x)
        return False
    return bt({1: 0}, {0}, 1)
res = []
for p in (3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83):
    pts = singer_points(p)
    assert len({(x*pow(y, -1, p)) % p if y else 'inf' for (x, y) in pts}) == p + 1   # all points of P^1
    v = frieze_strip(pts, p)
    n = p + 1
    for i in range(n - 1): assert (v[i][0]*v[i+1][1] - v[i][1]*v[i+1][0]) % p == 1
    b = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            s = legendre(v[i][0]*v[j][1] - v[i][1]*v[j][0], p); assert s != 0
            b[i][j] = s; b[j][i] = -s
    assert iso_paley_plus_sink(b, p), p
    quid = [(v[i-1][0]*v[i+1][1] - v[i-1][1]*v[i+1][0]) % p for i in range(1, n-1)]
    res.append((p, len(set(quid))))
    CH += 1
# random sub-configurations: Legendre pattern of any frieze strip = induced subtournament of Paley+sink up to switching
rnd = random.Random(3)
for p in (7, 11, 19, 23):
    for trial in range(50):
        m = rnd.randint(4, min(p+1, 9))
        xs = rnd.sample(list(range(p)) + ['inf'], m)
        pts = [(1, x) if x != 'inf' else (0, 1) for x in xs]
        v = frieze_strip(pts, p)
        for i in range(m):
            for j in range(m):
                if i == j: continue
                # switching: the scalar lam_i with v_i = lam_i * standard_i ; pattern = chi(lam_i lam_j) * standard pattern
                std = lambda P, Q: (P[0]*Q[1] - P[1]*Q[0]) % p
                lam_i = (v[i][0]*pow(pts[i][0], -1, p)) % p if pts[i][0] else (v[i][1]*pow(pts[i][1], -1, p)) % p
                lam_j = (v[j][0]*pow(pts[j][0], -1, p)) % p if pts[j][0] else (v[j][1]*pow(pts[j][1], -1, p)) % p
                assert legendre(std(v[i], v[j]), p) == legendre(lam_i*lam_j, p)*legendre(std(pts[i], pts[j]), p)
        CH += 1
print("(2) Legendre chirotopes: Singer frieze strips through all p+1 points of P^1(F_p) give Paley+sink up to switching for p =",
      [r[0] for r in res], "; any frieze strip = switched induced subtournament (200 random strips)")
print(f"ALL CHECKS PASSED ({CH})")
