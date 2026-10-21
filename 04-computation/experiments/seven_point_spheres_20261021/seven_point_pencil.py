"""The seven labels of a tetrahedron lie on TWO spheres through a common circle (session opus-2026-10-21-S49).
Note: 05-knowledge/results/seven_point_pencil_fano_and_converse_gluck_20261021.md; THM-4742.

Labels (circumcentre O, radius R, centroid G, opposite face V1V2V3 with circumcentre O', centroid G', orthocentre H'):
  A = second point of the line V0G on the circumsphere;  B = pole of the face plane (tangent planes at V1, V2, V3);
  C = V0;  D = O + 4(G - O);  E = inverse of G in the circumsphere;  F = O + 4(O' - O);  H' = orthocentre of V1V2V3.
Normalise O = 0, R = 1 and write v = V0, s = V1 + V2 + V3, w = O' (so s.w = 3|w|^2 =: 3k), t = v + s = 4G, T = |t|^2.
Then C = v, D = t, E = 4t/T, F = 4w, B = w/k, H' = s - 2w, A = v + mu (s - 3v), mu = -2 v.(s - 3v)/|s - 3v|^2.

  (P1) PROVED (exact Gram-basis identities, sympy): with delta = v.(s x w) != 0, the spheres |x|^2 - 2u.x + 4 = 0 with
       S1: (2u.v, 2u.s, 2u.w) = (5, T - 1, 4k + 1)          [the concurrent session23b sphere through B, C, D, E, F]
       S2: (2u.v, 2u.s, 2u.w) = (2 v.s - 1, |s|^2 + 6, 4k + 1)
       satisfy: S1 contains B, C, D, E, F (and 4C); S2 contains A, B, D, E, F, H' (and 4A, 4H'/|H'|^2); both are invariant
       under the inversion x -> 4x/|x|^2 (radius 2R about O); Pow_S2(C) = Pow_S1(H') = 6 - 2 v.s = Sigma, where
       Sigma = |V0V1|^2 + |V0V2|^2 + |V0V3|^2 (in units of R^2); Pow_S1(A) = mu (16 - T); u1 - u2 is orthogonal to t and w,
       so the radical plane of S1, S2 is the plane O G O' and S1 = S2 iff Sigma = 0.
  (P1n) PROVED: the same identities for an n-simplex, symbolic in n and in R^2 (no normalisation, no square roots):
       with s = sum_(facet) (Vi - O), s.w = n|w|^2, constant term (n+1)R^2, H' = O' + n(G' - O'), the defect is
       2nR^2 - 2 v.s = sum |V0Vi|^2.
  (P2) FINITE-EXACT: random rational tetrahedra in actual coordinates (no normalisation): all of (P1), with
       Pow_S2(V0) = Sigma exactly, and the four points B, D, E, F concyclic on S1 cap S2.
       Also (V0 - O).(G - V0) = -Sigma/8 and mu = Sigma/|s - 3v|^2: Sigma = 0 iff the median V0G is tangent to the
       circumsphere at V0 iff A = C.
  (P3) FINITE-EXACT: over F_p (5 <= p <= 97), in every valid random configuration the six labels A, B, D, E, F, H' lie on
       one quadric, all seven do iff Sigma = 0, and Sigma = 0 iff A = C; such configurations exist for every such p.
  (P4) FINITE-EXACT: an explicit nondegenerate tetrahedron over Q(i) with Sigma = 0: all seven labels lie on one sphere,
       with A = C (six distinct points).  Over R (any ordered field) Sigma > 0 for a nondegenerate tetrahedron, so the
       seventh incidence never holds.
  (P5) FINITE-EXACT, every dimension n = 3..6 (corrected after the S49 audit: the first version tested spheres in R^n
       through n + 1 labels, but all labels lie in a 3-flat, so every system was singular and nothing was tested): with
       D = O + (n+1)(G - O), F = O + (n+1)(O' - O) and H' = O' + n(G' - O') (the facet's own D-point; the orthocentre for
       n = 3), the labels lie in the 3-flat O + span(V0 - O, s, O' - O), and there the whole pencil statement holds with
       power (n+1)R^2 and defect sum |V0Vi|^2.  In the plane (n = 2) v, s, w are always dependent: the statement is void.
  (P6) The classical Fano anchor: the seventh line of the real Fano picture has determinant -2 (char 2 only).
Prints ALL CHECKS PASSED.  Runtime about 10 seconds."""
import random, itertools, time
from fractions import Fraction as Fr
import sympy as sp

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


# ---------------------------------------------------------------- (P1) Gram-basis proof
def section_P1():
    print('(P1) identities in the basis (v, s, w)')
    a, b, p, k = sp.symbols('a b p k')          # a = v.s, b = |s|^2, p = v.w, k = |w|^2 ; s.w = 3k, |v| = 1
    Gm = sp.Matrix([[1, a, p], [a, b, 3 * k], [p, 3 * k, k]])
    sq = lambda X: sp.expand((sp.Matrix([X]) * Gm * sp.Matrix(X))[0, 0])
    T = sq([1, 1, 0])
    mu = -2 * (a - 3) / (b - 6 * a + 9)
    pts = {'A': [1 - 3 * mu, mu, 0], 'B': [0, 0, 1 / k], 'C': [1, 0, 0], 'D': [1, 1, 0],
           'E': [4 / T, 4 / T, 0], 'F': [0, 0, 4], 'H': [0, 1, -2]}
    S1 = (5, T - 1, 4 * k + 1)
    S2 = (2 * a - 1, b + 6, 4 * k + 1)
    Q = lambda S, X: sp.simplify(sq(X) - (X[0] * S[0] + X[1] * S[1] + X[2] * S[2]) + 4)
    q1 = {x: Q(S1, X) for x, X in pts.items()}
    q2 = {x: Q(S2, X) for x, X in pts.items()}
    print(f'     S1 powers: { {x: sp.factor(q) for x, q in q1.items()} }')
    print(f'     S2 powers: { {x: sp.factor(q) for x, q in q2.items()} }')
    ok = all(q1[x] == 0 for x in 'BCDEF') and all(q2[x] == 0 for x in 'ABDEFH')
    ok &= sp.simplify(q2['C'] - (6 - 2 * a)) == 0 and sp.simplify(q1['H'] - (6 - 2 * a)) == 0
    ok &= sp.simplify(q1['A'] - mu * (16 - T)) == 0
    # inversion invariance of any sphere with constant term 4: Q(4X/|X|^2) = 4 Q(X)/|X|^2
    X = [sp.Symbol('x1'), sp.Symbol('x2'), sp.Symbol('x3')]
    S = (sp.Symbol('l1'), sp.Symbol('l2'), sp.Symbol('l3'))
    nX = sq(X)
    Xi = [4 * xx / nX for xx in X]
    inv_ok = sp.simplify(Q(S, Xi) - 4 * Q(S, X) / nX) == 0
    # radical plane: u1 - u2 orthogonal to t = v + s and to w
    d = [S1[i] - S2[i] for i in range(3)]            # 2(u1 - u2).(v, s, w)
    rad_ok = sp.simplify(d[0] + d[1]) == 0 and d[2] == 0 and sp.simplify(d[0] - (6 - 2 * a)) == 0
    tang = sp.simplify((1 + a) / 4 - 1 - (-(6 - 2 * a) / 8)) == 0     # v.(t/4 - v) = -Sigma/8 in units R^2
    mu_ok = sp.simplify(mu - (6 - 2 * a) / (b - 6 * a + 9)) == 0
    check(ok and inv_ok and rad_ok and tang and mu_ok,
          'S1 contains B,C,D,E,F; S2 contains A,B,D,E,F,H\'; Pow_S2(C) = Pow_S1(H\') = 6 - 2 v.s; Pow_S1(A) = mu (16 - T); '
          'both spheres are invariant under x -> 4x/|x|^2; 2(u1 - u2) = (Sigma, -Sigma, 0) on (v, s, w), so the radical '
          'plane is the plane O G O\' and S1 = S2 iff Sigma = 0; (V0 - O).(G - V0) = -Sigma/8 and mu = Sigma/|s - 3v|^2, so '
          'Sigma = 0 iff V0G is tangent at V0 iff A = C (PROVED, exact identities)')


def section_P1n():
    print('(P1n) the same identities for an n-simplex, homogeneous in R^2 (no normalisation)')
    a, b, p, k, R2, n = sp.symbols('a b p k R2 n')     # a = v.s, b = |s|^2, p = v.w, k = |w|^2, |v|^2 = R2, s.w = n k
    Gm = sp.Matrix([[R2, a, p], [a, b, n * k], [p, n * k, k]])
    sq = lambda X: sp.expand((sp.Matrix([X]) * Gm * sp.Matrix(X))[0, 0])
    T = sq([1, 1, 0])
    q = [-n, 1, 0]                                    # q = t - (n+1) v = s - n v
    mu = -2 * (a - n * R2) / sq(q)                    # A = v + mu q (second point of V0G on |x|^2 = R2)
    pts = {'A': [1 - n * mu, mu, 0], 'B': [0, 0, R2 / k], 'C': [1, 0, 0], 'D': [1, 1, 0],
           'E': [(n + 1) * R2 / T, (n + 1) * R2 / T, 0], 'F': [0, 0, n + 1], 'H': [0, 1, -(n - 1)]}
    c0 = (n + 1) * R2
    S1 = ((n + 2) * R2, T - R2, (n + 1) * k + R2)
    S2 = (2 * a + (2 - n) * R2, b + 2 * n * R2, (n + 1) * k + R2)
    Q = lambda S, X: sp.simplify(sq(X) - (X[0] * S[0] + X[1] * S[1] + X[2] * S[2]) + c0)
    Sig = 2 * n * R2 - 2 * a
    ok = all(Q(S1, pts[x]) == 0 for x in 'BCDEF') and all(Q(S2, pts[x]) == 0 for x in 'ABDEFH')
    ok &= sp.simplify(Q(S2, pts['C']) - Sig) == 0 and sp.simplify(Q(S1, pts['H']) - Sig) == 0
    ok &= sp.simplify(sq(pts['A']) - R2) == 0
    d = [sp.simplify(S1[i] - S2[i]) for i in range(3)]
    ok &= sp.simplify(d[0] - Sig) == 0 and sp.simplify(d[1] + Sig) == 0 and d[2] == 0
    check(ok, 'for every n (symbolic) and every R^2: S1 contains B,C,D,E,F; S2 contains A,B,D,E,F,H\' (H\' = O\' + n(G\' - O\')); '
              'power (n+1)R^2 at O; Pow_S2(C) = Pow_S1(H\') = 2nR^2 - 2 v.s = sum |V0Vi|^2; radical plane O G O\' (PROVED, '
              'exact identities; needs only v, s, w independent)')


# ---------------------------------------------------------------- exact helpers (any field via Fractions)
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(a, s): return tuple(x * s for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def solve3(M, r, inv):
    d = det3(M)
    if d == 0:
        raise ZeroDivisionError
    di = inv(d); out = []
    for j in range(3):
        Mj = [list(row) for row in M]
        for i in range(3):
            Mj[i][j] = r[i]
        out.append(det3(Mj) * di)
    return tuple(out)


def configuration(V, inv, red=lambda x: x):
    """All labels for the tetrahedron V over a field with inverse inv; red reduces (mod p) when needed."""
    R = lambda X: tuple(red(x) for x in X)
    M = [mul(sub(V[i], V[0]), 2) for i in (1, 2, 3)]
    r = [dot(V[i], V[i]) - dot(V[0], V[0]) for i in (1, 2, 3)]
    O = R(solve3(M, r, inv))
    R2 = red(dot(sub(V[0], O), sub(V[0], O)))
    G = R(mul(add(add(V[0], V[1]), add(V[2], V[3])), inv(4)))
    n = cross(sub(V[2], V[1]), sub(V[3], V[1]))
    O1 = R(solve3([mul(sub(V[2], V[1]), 2), mul(sub(V[3], V[1]), 2), n],
                  [dot(V[2], V[2]) - dot(V[1], V[1]), dot(V[3], V[3]) - dot(V[1], V[1]), dot(n, V[1])], inv))
    g = sub(G, O); w = sub(O1, O); d = sub(G, V[0])
    gg, ww, dd = red(dot(g, g)), red(dot(w, w)), red(dot(d, d))
    if 0 in (R2, gg, ww, dd):
        raise ZeroDivisionError
    A = R(add(V[0], mul(d, red(-2 * dot(sub(V[0], O), d) * inv(dd)))))
    B = R(add(O, mul(w, red(R2 * inv(ww)))))
    E = R(add(O, mul(g, red(R2 * inv(gg)))))
    D = R(add(O, mul(g, 4)))
    F = R(add(O, mul(w, 4)))
    H = R(sub(add(add(V[1], V[2]), V[3]), mul(O1, 2)))
    Sig = red(sum(dot(sub(V[0], V[i]), sub(V[0], V[i])) for i in (1, 2, 3)))
    delta = red(dot(sub(V[0], O), cross(sub(add(add(V[1], V[2]), V[3]), mul(O, 3)), w)))
    return dict(A=A, B=B, C=tuple(V[0]), D=D, E=E, F=F, H=H, O=O, R2=R2, G=G, O1=O1, Sig=Sig, delta=delta)


def sphere_through(P, inv, red=lambda x: x):
    M = [mul(sub(P[i], P[0]), 2) for i in (1, 2, 3)]
    r = [dot(P[i], P[i]) - dot(P[0], P[0]) for i in (1, 2, 3)]
    c = tuple(red(x) for x in solve3(M, r, inv))
    kk = red(2 * dot(c, P[0]) - dot(P[0], P[0]))
    return c, kk


def power(c, kk, X, red=lambda x: x):
    return red(dot(X, X) - 2 * dot(c, X) + kk)


# ---------------------------------------------------------------- (P2) rational tetrahedra
def section_P2():
    print('(P2) random rational tetrahedra, actual coordinates')
    rng = random.Random(20261021)
    inv = lambda x: 1 / Fr(x)
    n = 0; ok = True; plane_ok = True
    while n < 150:
        V = [tuple(Fr(rng.randint(-25, 25), rng.randint(1, 6)) for _ in range(3)) for _ in range(4)]
        if dot(sub(V[1], V[0]), cross(sub(V[2], V[0]), sub(V[3], V[0]))) == 0:
            continue
        try:
            L = configuration(V, inv)
            if L['delta'] == 0 or L['D'] == L['E']:
                continue
            c1, k1 = sphere_through([L['B'], L['C'], L['D'], L['F']], inv)
            c2, k2 = sphere_through([L['D'], L['E'], L['F'], L['H']], inv)
        except ZeroDivisionError:
            continue
        n += 1
        O = L['O']; R2 = L['R2']
        iota = lambda X: add(O, mul(sub(X, O), 4 * R2 / dot(sub(X, O), sub(X, O))))
        on1 = all(power(c1, k1, L[x]) == 0 for x in 'BCDEF') and power(c1, k1, iota(L['C'])) == 0
        on2 = all(power(c2, k2, L[x]) == 0 for x in 'ABDEFH') and power(c2, k2, iota(L['A'])) == 0 \
            and power(c2, k2, iota(L['H'])) == 0
        ok &= on1 and on2
        ok &= power(c2, k2, L['C']) == L['Sig'] and power(c1, k1, L['H']) == L['Sig'] and L['Sig'] > 0
        ok &= power(c1, k1, O) == 4 * R2 and power(c2, k2, O) == 4 * R2
        # B, D, E, F coplanar (plane O G O') and C, A, D, E coplanar (plane O V0 G)
        cop = lambda P: dot(sub(P[1], P[0]), cross(sub(P[2], P[0]), sub(P[3], P[0]))) == 0
        plane_ok &= cop([L['B'], L['D'], L['E'], L['F']]) and cop([O, L['G'], L['O1'], L['B']])
    check(ok and plane_ok,
          f'{n} rational tetrahedra: S1 contains B,C,D,E,F,iota(C); S2 contains A,B,D,E,F,H\',iota(A),iota(H\'); both have power '
          f'4R^2 at O; Pow_S2(V0) = Pow_S1(H\') = Sigma = sum |V0Vi|^2 > 0; B,D,E,F lie in the plane O G O\' (FINITE-EXACT)')


# ---------------------------------------------------------------- (P3) finite fields
def section_P3():
    print('(P3) over F_p: seven labels on one quadric iff Sigma = 0')
    allok = True; rows = []
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97):
        rng = random.Random(p)
        def inv(x, p=p):
            if x % p == 0:
                raise ZeroDivisionError
            return pow(x % p, p - 2, p)
        red = lambda x, p=p: x % p
        valid = seven = sig0 = mism = 0
        for _ in range(3000):
            V = [tuple(rng.randrange(p) for _ in range(3)) for _ in range(4)]
            if dot(sub(V[1], V[0]), cross(sub(V[2], V[0]), sub(V[3], V[0]))) % p == 0:
                continue
            try:
                L = configuration(V, inv, red)
                if L['delta'] == 0 or L['D'] == L['E'] or L['B'] == L['F']:
                    continue
                c, kk = sphere_through([L['D'], L['E'], L['F'], L['H']], inv, red)
            except ZeroDivisionError:
                continue
            valid += 1
            all7 = all(power(c, kk, L[x], red) == 0 for x in 'ABCDEFH')
            six = all(power(c, kk, L[x], red) == 0 for x in 'ABDEFH')
            seven += all7; sig0 += L['Sig'] == 0
            mism += (all7 != (L['Sig'] == 0)) or not six or ((L['Sig'] == 0) != (L['A'] == L['C']))
        rows.append((p, valid, seven, sig0, mism))
        allok &= mism == 0 and seven > 0
    print(f'     (p, valid configurations, all seven on one quadric, Sigma = 0, mismatches): {rows}')
    check(allok, 'for every p in 5..97: in every valid configuration the six labels A,B,D,E,F,H\' share a quadric, all seven do '
          'iff Sigma = 0 iff A = C, and such configurations exist (FINITE-EXACT)')


# ---------------------------------------------------------------- (P4) Q(i) example and the real obstruction
def section_P4():
    print('(P4) an explicit seven-label sphere over Q(i)')
    I = sp.I
    inv = lambda x: 1 / sp.nsimplify(x)
    red = lambda x: sp.nsimplify(sp.expand(x))
    V0 = (sp.Integer(0), sp.Integer(0), sp.Integer(0))
    V = [V0, (sp.Integer(1), sp.Integer(0), sp.Integer(0)), (sp.Integer(0), sp.Integer(1), sp.Integer(0)),
         (sp.Integer(0), I, I)]
    L = configuration(V, inv, red)
    nondeg = (red(det3([sub(V[1], V[0]), sub(V[2], V[0]), sub(V[3], V[0])])) != 0 and L['delta'] != 0
              and L['D'] != L['E'] and L['B'] != L['F'] and L['R2'] != 0)
    c, kk = sphere_through([L['D'], L['E'], L['F'], L['H']], inv, red)
    all7 = all(red(power(c, kk, L[x])) == 0 for x in 'ABCDEFH')
    six_distinct = len({L[x] for x in 'ABCDEFH'}) == 6 and L['A'] == L['C']
    tangent = red(dot(sub(L['C'], L['O']), sub(L['G'], L['C']))) == 0
    print(f'     V0 = 0, V1 = e1, V2 = e2, V3 = (0, i, i): Sigma = {L["Sig"]}, R^2 = {L["R2"]}, delta = {L["delta"]}, A = C: {L["A"] == L["C"]}')
    check(nondeg and all7 and six_distinct and tangent and L['Sig'] == 0,
          'over Q(i) the tetrahedron V0 = 0, V1 = e1, V2 = e2, V3 = (0, i, i) is nondegenerate with Sigma = 0: all seven labels lie '
          'on one sphere, the median V0G is tangent at V0, and A = C, so six distinct points (FINITE-EXACT); over R, Sigma = sum of '
          'squared edge lengths > 0')


# ---------------------------------------------------------------- (P5) other dimensions
def solve_n(M, r):
    m = len(M); Aug = [list(M[i]) + [r[i]] for i in range(m)]
    for col in range(m):
        piv = next((i for i in range(col, m) if Aug[i][col] != 0), None)
        if piv is None:
            raise ZeroDivisionError
        Aug[col], Aug[piv] = Aug[piv], Aug[col]
        for i in range(m):
            if i != col and Aug[i][col] != 0:
                f = Aug[i][col] / Aug[col][col]
                Aug[i] = [x - f * y for x, y in zip(Aug[i], Aug[col])]
    return [Aug[i][m] / Aug[i][i] for i in range(m)]


def circ_affine(P):
    kd = len(P) - 1; base = P[0]
    Ev = [sub(P[i], base) for i in range(1, kd + 1)]
    cc = solve_n([[2 * dot(Ev[i], Ev[j]) for j in range(kd)] for i in range(kd)], [dot(Ev[i], Ev[i]) for i in range(kd)])
    X = tuple(base)
    for j in range(kd):
        X = add(X, mul(Ev[j], cc[j]))
    return X


def labels_nd(V):
    n = len(V) - 1
    O = circ_affine(V); R2 = dot(sub(V[0], O), sub(V[0], O))
    G = tuple(sum(v[i] for v in V) / Fr(n + 1) for i in range(n))
    d = sub(G, V[0])
    A = add(V[0], mul(d, -2 * dot(sub(V[0], O), d) / dot(d, d)))
    Of = circ_affine(V[1:]); Gf = tuple(sum(v[i] for v in V[1:]) / Fr(n) for i in range(n))
    w = sub(Of, O); g = sub(G, O)
    lab = dict(A=A, B=add(O, mul(w, R2 / dot(w, w))), C=tuple(V[0]), D=add(O, mul(g, n + 1)),
               E=add(O, mul(g, R2 / dot(g, g))), F=add(O, mul(w, n + 1)), H=add(Of, mul(sub(Gf, Of), n)))
    sv = (sum(V[i][j] for i in range(1, n + 1)) - n * O[j] for j in range(n))
    return lab, O, R2, sub(V[0], O), tuple(sv), w


def section_P5():
    print('(P5) other dimensions: the pencil inside the 3-flat O + span(v, s, w)')
    rng = random.Random(5)
    summary = {}
    for n in (2, 3, 4, 5, 6):
        good = degenerate = 0; ok = True; trials = 0
        while trials < 20:
            V = [tuple(Fr(rng.randint(-9, 9), rng.randint(1, 3)) for _ in range(n)) for _ in range(n + 1)]
            try:
                lab, O, R2, v, sv, w = labels_nd(V)
            except ZeroDivisionError:
                continue
            trials += 1
            basis = [v, sv, w]
            Gm = [[dot(a, b) for b in basis] for a in basis]
            try:
                coords = {}
                for x, X in lab.items():
                    rhs = [dot(e, sub(X, O)) for e in basis]
                    c = solve_n(Gm, rhs)
                    recon = tuple(O[j] + sum(c[i] * basis[i][j] for i in range(3)) for j in range(n))
                    if recon != tuple(X):
                        raise ValueError('label off the 3-flat')
                    coords[x] = c
            except ZeroDivisionError:
                degenerate += 1          # v, s, w dependent: delta = 0
                continue
            quad = lambda c: sum(c[i] * Gm[i][j] * c[j] for i in range(3) for j in range(3))
            def fit(names):
                M = [[-coords[x][0], -coords[x][1], -coords[x][2], 1] for x in names]
                return solve_n(M, [-quad(coords[x]) for x in names])
            def Q(sph, x):
                c = coords[x]
                return quad(c) - sum(c[i] * sph[i] for i in range(3)) + sph[3]
            Sig = sum(dot(sub(V[0], V[i]), sub(V[0], V[i])) for i in range(1, n + 1))
            S2 = fit(['D', 'E', 'F', 'H']); S1 = fit(['B', 'C', 'D', 'F'])
            ok &= all(Q(S2, x) == 0 for x in 'ABDEFH') and all(Q(S1, x) == 0 for x in 'BCDEF')
            ok &= S1[3] == (n + 1) * R2 and S2[3] == (n + 1) * R2
            ok &= Q(S2, 'C') == Sig and Q(S1, 'H') == Sig
            ok &= [S1[i] - S2[i] for i in range(3)] == [Sig, -Sig, 0]
            good += 1
        summary[n] = (good, degenerate, ok)
        print(f'     n = {n}: {good} generic simplices checked, {degenerate} with v, s, w dependent; identities hold: {ok}')
    check(summary[2][0] == 0 and summary[2][1] == 20 and all(summary[n][0] == 20 and summary[n][2] for n in (3, 4, 5, 6)),
          'for every n = 3..6 the labels (D = O + (n+1)(G - O), F = O + (n+1)(O\' - O), H\' = O\' + n(G\' - O\'), the facet\'s own '
          'D-point) lie in the 3-flat O + span(V0 - O, s, O\' - O); there S2 contains A,B,D,E,F,H\', S1 contains B,C,D,E,F, both have '
          'power (n+1)R^2 at O, Pow_S2(C) = Pow_S1(H\') = sum |V0Vi|^2 and the radical plane is O G O\'; in the plane (n = 2) v, s, w '
          'are always dependent and the statement is void (FINITE-EXACT, exact Gram-coordinate checks)')


# ---------------------------------------------------------------- (P6) the Fano anchor
def section_P6():
    print('(P6) the classical Fano anchor')
    pts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
    lines = [t for t in itertools.combinations(range(7), 3)
             if all((pts[t[0]][i] + pts[t[1]][i] + pts[t[2]][i]) % 2 == 0 for i in range(3))]
    dets = sorted(abs(det3([pts[t[0]], pts[t[1]], pts[t[2]]])) for t in lines)
    print(f'     |det| of the seven Fano lines in the 0/1 coordinates: {dets}')
    check(len(lines) == 7 and dets == [0, 0, 0, 0, 0, 0, 2],
          'six Fano lines are real collinearities; the seventh, {(1,1,0),(1,0,1),(0,1,1)}, has determinant -2: '
          'realizable only where 2 = 0 (CLASSICAL)')


if __name__ == '__main__':
    t0 = time.time()
    section_P1(); section_P1n(); section_P2(); section_P3(); section_P4(); section_P5(); section_P6()
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
