"""procgen_fencelim_20260926_audit.py -- Task A: audit of the junction inequality, the face angle
identity, the Corner Lemma, the per-face inequality and the finite bound, on exact configurations.

Session collatz-procgen-20260922, lane "fencelim" (2026-09-26).  Called by the runner.
"""
from fractions import Fraction as Fr
import math
import mpmath

import procgen_fencelim_20260926_geom as G
from procgen_fencelim_20260926_geom import QS, analyze, corner_report, angle_identity_check

mpmath.mp.dps = 50
G.MP_TOL = mpmath.mpf(10) ** -30

# ----------------------------------------------------------------------------- constants
PI = mpmath.pi


def P_reg(k):
    """perimeter of the regular k-gon of unit area."""
    return 2 * mpmath.sqrt(k * mpmath.tan(PI / k))


P5 = P_reg(5)
MU = 1 / (8 - P5)
RHO = (4 - P5) / (2 * (8 - P5))
LAM = (6 - P5) / (8 - P5)


def finite_bound(n):
    """largest A with A + 2 mu sqrt(pi) sqrt(A) <= 2(mu+rho) n - 6 rho."""
    c = 2 * (MU + RHO) * n - 6 * RHO
    if c <= 0:
        return mpmath.mpf(0)
    b = 2 * MU * mpmath.sqrt(PI)
    x = (-b + mpmath.sqrt(b * b + 4 * c)) / 2
    return x * x


def U_prev(n):
    """previous lane's Theorem F2 bound."""
    return ((mpmath.sqrt(1 + 4 * n / mpmath.sqrt(PI)) - 1) / 2) ** 2


# ----------------------------------------------------------------------------- builders
def F(x):
    return Fr(x)


def seg(x0, y0, x1, y1):
    return ((x0, y0), (x1, y1))


def grid(i, j):
    """i x j grid of unit squares, one fence per unit edge."""
    fs = []
    for y in range(j + 1):
        for x in range(i):
            fs.append(seg(F(x), F(y), F(x + 1), F(y)))
    for x in range(i + 1):
        for y in range(j):
            fs.append(seg(F(x), F(y), F(x), F(y + 1)))
    return fs


def polyomino(cells):
    """unit fences on the boundary edges of the union of unit cells (every edge of every cell)."""
    es = set()
    for (x, y) in cells:
        es.add((x, y, x + 1, y))
        es.add((x, y + 1, x + 1, y + 1))
        es.add((x, y, x, y + 1))
        es.add((x + 1, y, x + 1, y + 1))
    return [seg(F(a), F(b), F(c), F(d)) for (a, b, c, d) in sorted(es)]


def quasi_square(k):
    w = math.isqrt(k - 1) + 1 if k > 1 else 1
    cells = [(c % w, c // w) for c in range(k)]
    return polyomino(cells)


SQ = [seg(F(0), F(0), F(1), F(0)), seg(F(1), F(0), F(1), F(1)), seg(F(1), F(1), F(0), F(1)), seg(F(0), F(1), F(0), F(0))]


def shifted(fs, dx, dy):
    return [((a[0] + dx, a[1] + dy), (b[0] + dx, b[1] + dy)) for (a, b) in fs]


def brick_wall():
    """3 rows of running bond (half bricks at the ends of the middle row); straight joins and T's."""
    fs = []
    for y in range(4):
        for x in range(3):
            fs.append(seg(F(x), F(y), F(x + 1), F(y)))
    for y in range(3):
        fs.append(seg(F(0), F(y), F(0), F(y + 1)))
        fs.append(seg(F(3), F(y), F(3), F(y + 1)))
    for x in (1, 2):
        fs.append(seg(F(x), F(0), F(x), F(1)))
        fs.append(seg(F(x), F(2), F(x), F(3)))
    for x in (Fr(1, 2), Fr(3, 2), Fr(5, 2)):
        fs.append(seg(x, F(1), x, F(2)))
    return fs


def rect_with_X():
    """1 x 2 rectangle, middle chord y=1, two unit fences ending at (1/2,1) from opposite sides."""
    fs = [seg(F(0), F(0), F(0), F(1)), seg(F(0), F(1), F(0), F(2)), seg(F(1), F(0), F(1), F(1)),
          seg(F(1), F(1), F(1), F(2)), seg(F(0), F(0), F(1), F(0)), seg(F(0), F(2), F(1), F(2)),
          seg(F(0), F(1), F(1), F(1)), seg(Fr(1, 2), F(0), Fr(1, 2), F(1)), seg(Fr(1, 2), F(1), Fr(1, 2), F(2))]
    return fs


def q3(a, b=0):
    return QS(a, b, 3)


def hexagon_spokes():
    """regular hexagon of side 1 plus three spokes from the centre to alternate vertices (Y's)."""
    V = [(q3(1), q3(0)), (q3(Fr(1, 2)), q3(0, Fr(1, 2))), (q3(-Fr(1, 2)), q3(0, Fr(1, 2))),
         (q3(-1), q3(0)), (q3(-Fr(1, 2)), q3(0, -Fr(1, 2))), (q3(Fr(1, 2)), q3(0, -Fr(1, 2)))]
    O = (q3(0), q3(0))
    fs = [(V[i], V[(i + 1) % 6]) for i in range(6)]
    fs += [(O, V[0]), (O, V[2]), (O, V[4])]
    return fs


def nested_triangle():
    """unit square plus a disjoint equilateral unit triangle inside it (a nested component)."""
    c, s = Fr(24, 25), Fr(7, 25)          # rotation with cos=24/25, sin=7/25
    dx, dy = Fr(1, 50), Fr(1, 70)
    p0 = (q3(dx), q3(dy))
    p1 = (q3(c + dx), q3(s + dy))
    p2 = (q3(c / 2 + dx, -s / 2), q3(s / 2 + dy, c / 2))   # rotation of (1/2, sqrt3/2)
    sq = [((q3(a[0]), q3(a[1])), (q3(b[0]), q3(b[1]))) for (a, b) in SQ]
    return sq + [(p0, p1), (p1, p2), (p2, p0)]


def two_squares_bridge():
    fs = SQ + shifted(SQ, F(2), F(0))
    fs.append(seg(F(1), Fr(1, 2), F(2), Fr(1, 2)))
    return fs


def figure8_in_frame():
    """two unit squares sharing a corner, inside a 4 x 4 frame (pinched hole boundary)."""
    fr = []
    for x in range(4):
        fr.append(seg(F(x), F(0), F(x + 1), F(0)))
        fr.append(seg(F(x), F(4), F(x + 1), F(4)))
        fr.append(seg(F(0), F(x), F(0), F(x + 1)))
        fr.append(seg(F(4), F(x), F(4), F(x + 1)))
    return fr + shifted(SQ, F(1), F(1)) + shifted(SQ, F(2), F(2))


def pinched_outer():
    return SQ + shifted(SQ, F(1), F(1))


def square_two_chords():
    return SQ + [seg(F(0), Fr(3, 5), Fr(4, 5), F(0)), seg(F(1), Fr(2, 5), Fr(1, 5), F(1))]


def square_mid_chord():
    return SQ + [seg(F(0), Fr(1, 2), F(1), Fr(1, 2))]


def pythagorean_patch(a=Fr(3, 5), R2=Fr(9, 1), N=8):
    """fences (maximal segments, length a+b=1) of the Pythagorean tiling with squares a, b=1-a;
    keep fences inside the disc |p|^2<=R2 and prune fences with an unsupported end."""
    b = 1 - a
    v1, v2 = (a, b), (-b, a)
    hor, ver = {}, {}
    for i in range(-N, N + 1):
        for j in range(-N, N + 1):
            px, py = i * v1[0] + j * v2[0], i * v1[1] + j * v2[1]
            for (ox, oy, s) in ((px, py, a), (px + a - b, py + a, b)):
                hor.setdefault(oy, []).append((ox, ox + s))
                hor.setdefault(oy + s, []).append((ox, ox + s))
                ver.setdefault(ox, []).append((oy, oy + s))
                ver.setdefault(ox + s, []).append((oy, oy + s))

    def merge(iv):
        iv = sorted(iv)
        out = []
        for (l, r) in iv:
            if out and l <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], r))
            else:
                out.append((l, r))
        return out
    fs = []
    for y, iv in hor.items():
        for (l, r) in merge(iv):
            if r - l == 1:
                fs.append(((l, y), (r, y)))
    for x, iv in ver.items():
        for (l, r) in merge(iv):
            if r - l == 1:
                fs.append(((x, l), (x, r)))
    fs = [f for f in fs if all(p[0] ** 2 + p[1] ** 2 <= R2 for p in f)]
    while True:
        keep = []
        for f in fs:
            ok = True
            for e in f:
                sup = False
                for g in fs:
                    if g is f:
                        continue
                    if G.param_on(e, *g) is not None:
                        sup = True
                        break
                ok = ok and sup
            if ok:
                keep.append(f)
        if len(keep) == len(fs):
            break
        fs = keep
    return fs


def n6_record():
    """reproduce the n=6 record (unit-sided pentagon plus one unit chord), polished at 50 digits."""
    import numpy as np
    from scipy.optimize import minimize

    def build(v, sa, sb, lib):
        cos, sin = lib.cos, lib.sin
        ph = [0] + list(v[:4])
        V = [(0, 0)]
        for i in range(4):
            V.append((V[-1][0] + cos(ph[i]), V[-1][1] + sin(ph[i])))
        s, t = v[4], v[5]
        A = (V[sa][0] + s * (V[(sa + 1) % 5][0] - V[sa][0]), V[sa][1] + s * (V[(sa + 1) % 5][1] - V[sa][1]))
        B = (V[sb][0] + t * (V[(sb + 1) % 5][0] - V[sb][0]), V[sb][1] + t * (V[(sb + 1) % 5][1] - V[sb][1]))
        return V, A, B

    def area(P):
        return sum(P[i][0] * P[(i + 1) % len(P)][1] - P[(i + 1) % len(P)][0] * P[i][1] for i in range(len(P))) / 2

    def pieces(V, A, B, sa, sb):
        P1 = [A]
        i = (sa + 1) % 5
        while True:
            P1.append(V[i])
            if i == sb:
                break
            i = (i + 1) % 5
        P1.append(B)
        P2 = [B]
        i = (sb + 1) % 5
        while True:
            P2.append(V[i])
            if i == sa:
                break
            i = (i + 1) % 5
        P2.append(A)
        return P1, P2

    sa, sb = 0, 2
    rng = np.random.default_rng(7)
    best = None
    for trial in range(40):
        x0 = np.concatenate([np.cumsum(np.full(4, 2 * np.pi / 5)) + rng.normal(0, 0.2, 4), rng.uniform(0.2, 0.8, 2)])

        def cl(v):
            V, A, B = build(v, sa, sb, np)
            return (V[4][0] - 0) ** 2 + (V[4][1] - 0) ** 2 - 1

        def ch(v):
            V, A, B = build(v, sa, sb, np)
            return (A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 - 1

        def ar(v):
            V, A, B = build(v, sa, sb, np)
            P1, P2 = pieces(V, A, B, sa, sb)
            return area(P1), area(P2), area(V)
        cons = [{'type': 'eq', 'fun': cl}, {'type': 'eq', 'fun': ch},
                {'type': 'ineq', 'fun': lambda v: 1 - ar(v)[0]}, {'type': 'ineq', 'fun': lambda v: 1 - ar(v)[1]},
                {'type': 'ineq', 'fun': lambda v: np.array([v[4] - 0.01, 0.99 - v[4], v[5] - 0.01, 0.99 - v[5]])}]
        r = minimize(lambda v: -ar(v)[2], x0, constraints=cons, method='SLSQP', options={'maxiter': 600, 'ftol': 1e-14})
        if r.success and abs(cl(r.x)) < 1e-9 and abs(ch(r.x)) < 1e-9:
            a1, a2, at = ar(r.x)
            if a1 <= 1 + 1e-9 and a2 <= 1 + 1e-9 and a1 > 0 and a2 > 0 and (best is None or at > best[0]):
                best = (at, r.x.copy(), a1, a2)
    at, x, a1, a2 = best
    big = 0 if a1 > a2 else 1
    # polish: keep v[0], v[1] (two angles) fixed, solve closure(2), chord, big piece = 1 for v[2], v[3], s, t
    fixed = [mpmath.mpf(float(x[0]))]   # direction of side 2; sides 3,4 directions and s are solved for

    # closure in exact form: side 5 = V4 -> V0 must have unit length: use |V4|^2 = 1
    def system(p2, p3, s, t):
        v = fixed + [p2, p3, 0, s, t]   # v[3] is unused by build()
        V, A, B = build(v, sa, sb, mpmath)
        P1, P2 = pieces(V, A, B, sa, sb)
        return [V[4][0] ** 2 + V[4][1] ** 2 - 1, (A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2 - 1,
                (area(P1) if big == 0 else area(P2)) - 1]
    # 3 equations in 4 unknowns: also fix t
    tfix = mpmath.mpf(float(x[5]))
    sol = [mpmath.mpf(float(x[1])), mpmath.mpf(float(x[2])), mpmath.mpf(float(x[4]))]
    h = mpmath.mpf(10) ** -20
    for it in range(30):   # Newton with a central-difference Jacobian at 50 digits
        F0 = mpmath.matrix(system(sol[0], sol[1], sol[2], tfix))
        if mpmath.norm(F0) < mpmath.mpf(10) ** -45:
            break
        Jm = mpmath.matrix(3, 3)
        for c in range(3):
            sp = list(sol)
            sm = list(sol)
            sp[c] += h
            sm[c] -= h
            Fp = system(sp[0], sp[1], sp[2], tfix)
            Fm = system(sm[0], sm[1], sm[2], tfix)
            for r_ in range(3):
                Jm[r_, c] = (Fp[r_] - Fm[r_]) / (2 * h)
        dx = mpmath.lu_solve(Jm, -F0)
        sol = [sol[c] + dx[c] for c in range(3)]
    v = fixed + [sol[0], sol[1], 0, sol[2], tfix]
    V, A, B = build(v, sa, sb, mpmath)
    V = [(mpmath.mpf(p[0]), mpmath.mpf(p[1])) for p in V]
    fs = [(V[i], V[(i + 1) % 5]) for i in range(5)] + [(A, B)]
    return fs, float(at)


# ----------------------------------------------------------------------------- per-configuration audit
def audit_config(name, fs):
    R = analyze(fs)
    rep = corner_report(R)
    worst = angle_identity_check(R)
    areas = [f['A2'] / 2 for f in R['faces']]
    A = sum(areas, 0)
    Amp = G.tomp(A) if not isinstance(A, mpmath.mpf) else A
    valid = all(G.sgn(a - 1) <= 0 for a in areas)
    # per-face inequality (4): a <= mu P_kappa sqrt(a) + 2 rho (kappa-3) and a <= mu P + 2 rho (kappa - 3)
    ok4 = True
    for f in R['faces']:
        a = G.tomp(f['A2'] / 2) if not isinstance(f['A2'], mpmath.mpf) else f['A2'] / 2
        k = f['kappa']
        if a > 1:
            continue
        if not (a <= MU * f['per'] + 2 * RHO * (k - 3) + mpmath.mpf(10) ** -25):
            ok4 = False
        if not (f['per'] >= P_reg(k) * mpmath.sqrt(a) - mpmath.mpf(10) ** -25):
            ok4 = False
    n = len(fs)
    sumP = sum((f['per'] for f in R['faces']), mpmath.mpf(0))
    ok_perim = abs(sumP + R['P_o'] - 2 * n) < mpmath.mpf(10) ** -25
    fb = finite_bound(n)
    ok5 = (not valid) or (Amp <= fb)
    return dict(name=name, n=n, R=R, rep=rep, worst=worst, A=Amp, valid=valid, ok4=ok4, ok5=ok5, fb=fb, ok_perim=ok_perim)


def configurations():
    cs = []
    for (i, j) in ((1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (3, 4), (4, 4)):
        cs.append((f'grid {i}x{j}', grid(i, j)))
    for k in (3, 5, 7, 10, 13):
        cs.append((f'quasi-square {k}-omino', quasi_square(k)))
    cs.append(('L-tromino', polyomino([(0, 0), (1, 0), (0, 1)])))
    cs.append(('S-tetromino', polyomino([(0, 0), (1, 0), (1, 1), (2, 1)])))
    cs.append(('ring of 8 cells (enclosed centre)', polyomino([(x, y) for x in range(3) for y in range(3) if (x, y) != (1, 1)])))
    cs.append(('n=5 record: square + chord (0,3/5)-(4/5,0)', SQ + [seg(F(0), Fr(3, 5), Fr(4, 5), F(0))]))
    cs.append(('square + two chords', square_two_chords()))
    cs.append(('square + mid chord', square_mid_chord()))
    cs.append(('1x2 rectangle with an X junction', rect_with_X()))
    cs.append(('brick wall (running bond, 3 rows)', brick_wall()))
    cs.append(('hexagon + 3 spokes (Y junctions), Q(sqrt3)', hexagon_spokes()))
    cs.append(('nested: unit triangle inside unit square, Q(sqrt3)', nested_triangle()))
    cs.append(('two squares + bridge fence', two_squares_bridge()))
    cs.append(('two squares sharing a corner (pinched outer face)', pinched_outer()))
    cs.append(('figure-8 component inside a 4x4 frame (pinched hole)', figure8_in_frame()))
    cs.append(('two disjoint components', SQ + shifted(SQ, F(3), F(0))))
    return cs
