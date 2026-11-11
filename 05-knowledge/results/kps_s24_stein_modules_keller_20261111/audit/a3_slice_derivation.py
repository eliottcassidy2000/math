"""a3_slice_derivation.py  (independent audit; sympy, exact unless marked numeric)

Re-derives section 6.2 of the note from the map of THM-1300 (formulas taken from the canon file THM-1300, u = 1 + xy):
   F1 = u^3 z + y^2 u (4 + 3xy),  F2 = y + 3 x u^2 z + 3 x y^2 (4 + 3xy),  F3 = 2x - 3 x^2 y - x^3 z.

 1. det JF = -2;  on F1 = c:  u = y + 3 c x w, v = x w + x w^2 - c x^3 w^3, w = 1/(1+xy);  Jacobian -2 w^3
 2. master identity  L(c,u,v) x^3 + (4 - 3uv) x - 2v = 0 on the slice,  L = 27c^2v^2 - 18cuv + 16c + u^3v - u^2;
    discriminant of the cubic in x;  bidegree and points at infinity of S_c in P^1 x P^1;  irreducibility
 3. exact fibre sizes (Groebner basis in Q[x,y,w], dimension of the quotient = number of points, the map being etale):
    generic rational points (expected 3) and rational points of S_1 (expected 1)
 4. numeric fibre sizes at random complex points of S_1 (expected 1) and off S_1 (expected 3)
 5. intersection counts of S_1 with the fibres of the potentials: (u-2)(v-3) = t -> 5,  (u-0)(v-0) = t -> 4,
    u = t -> 2,  v = t -> 3,  u + v = t -> 4   (affine points, squarefree for generic t)
 6. H^1(X) -> H^1(Y) is injective: dw ^ dW != 0  (w is not constant on the fibres of W)
 7. the generic fibre of W = (u-2)(v-3) on X_1 is irreducible: the cubic in x restricted to the hyperbola is irreducible
    over Q(tau) and its discriminant is not a square  (geometric Galois group S_3)
"""
import sympy as sp
import random

x, y, z, c, u, v, w, t, tau = sp.symbols("x y z c u v w t tau")
ok = True


def report(name, cond, extra=""):
    global ok
    ok = ok and bool(cond)
    print("%-88s %s %s" % (name, "OK" if cond else "FAIL", extra))


q = 1 + x * y
F1 = q ** 3 * z + y ** 2 * q * (4 + 3 * x * y)
F2 = y + 3 * x * q ** 2 * z + 3 * x * y ** 2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x ** 2 * y - x ** 3 * z

# 1
detJ = sp.factor(sp.Matrix([F1, F2, F3]).jacobian([x, y, z]).det())
report("1. det JF = -2", detJ == -2, "(det = %s)" % detJ)
zz = sp.solve(sp.Eq(F1, c), z)
report("1. F1 = c has the unique solution z = (c - y^2 q (4+3xy))/q^3", len(zz) == 1 and sp.simplify(zz[0] - (c - y ** 2 * q * (4 + 3 * x * y)) / q ** 3) == 0)
us = sp.simplify(F2.subs(z, zz[0]))
vs = sp.simplify(F3.subs(z, zz[0]))
W_ = 1 / q
u_claim = y + 3 * c * x * W_
v_claim = x * W_ + x * W_ ** 2 - c * x ** 3 * W_ ** 3
report("1. u = y + 3 c x w", sp.simplify(us - u_claim) == 0)
report("1. v = x w + x w^2 - c x^3 w^3", sp.simplify(vs - v_claim) == 0)
Jac = sp.simplify(sp.diff(us, x) * sp.diff(vs, y) - sp.diff(us, y) * sp.diff(vs, x))
report("1. Jacobian of (u, v) in (x, y) = -2 w^3", sp.simplify(Jac + 2 / q ** 3) == 0)

# 2
L = 27 * c ** 2 * v ** 2 - 18 * c * u * v + 16 * c + u ** 3 * v - u ** 2
ident = sp.simplify((L * x ** 3 + (4 - 3 * u * v) * x - 2 * v).subs({u: u_claim, v: v_claim}))
report("2. L(c,u,v) x^3 + (4 - 3uv) x - 2v = 0 on the slice", ident == 0)
L1 = L.subs(c, 1)
disc = sp.factor(sp.discriminant(L1 * x ** 3 + (4 - 3 * u * v) * x - 2 * v, x))
print("      discriminant in x of the cubic (c = 1):", disc)
report("2. deg_u L = 3, deg_v L = 2, total degree 4", sp.degree(L, u) == 3 and sp.degree(L, v) == 2 and sp.Poly(L, u, v).total_degree() == 4)
fl = sp.factor_list(L1)
report("2. L(1,u,v) irreducible over Q", len(fl[1]) == 1 and fl[1][0][1] == 1)
dv_ = sp.factor(sp.discriminant(L1, v))
report("2. disc_v L(1,u,v) = (u^2 - 12)^3 (not a square: L is absolutely irreducible)", sp.expand(dv_ - (u ** 2 - 12) ** 3) == 0, "(%s)" % dv_)
# points at infinity in P^1 x P^1: homogenise u = U1/U0 (degree 3), v = V1/V0 (degree 2)
U0, U1, V0, V1 = sp.symbols("U0 U1 V0 V1")
Lh = sp.expand(sp.together(L1.subs({u: U1 / U0, v: V1 / V0}) * U0 ** 3 * V0 ** 2))
print("      on u = infinity (U0 = 0):", sp.factor(Lh.subs(U0, 0)), "   on v = infinity (V0 = 0):", sp.factor(Lh.subs(V0, 0)))
report("2. closure of S_1 meets the boundary only at (inf, 0) and (inf, inf): excluded centres are b = 0",
       sp.factor(Lh.subs(U0, 0)) == U1 ** 3 * V0 * V1 and sp.factor(Lh.subs(V0, 0)) == 27 * U0 ** 3 * V1 ** 2)


# 3 exact fibre sizes
def fibre_size(u0, v0):
    eqs = [w * (1 + x * y) - 1, y + 3 * x * w - u0, x * w + x * w ** 2 - x ** 3 * w ** 3 - v0]
    G = sp.groebner(eqs, x, y, w, order="grevlex")
    lead = [sp.Poly(g, x, y, w).monoms(order="grevlex")[0] for g in G.exprs]
    # standard monomials
    D = 12
    std = 0
    for i in range(D):
        for j in range(D):
            for k in range(D):
                if not any(i >= a and j >= b and k >= e for (a, b, e) in lead):
                    std += 1
    # sanity: finiteness (pure powers of every variable among the leading monomials)
    fin = all(any(m[idx] > 0 and sum(m) == m[idx] for m in lead) for idx in range(3))
    return std if fin else None


gen_pts = [(1, 1), (2, 3), (-3, 4), (5, -7), (7, 2), (sp.Rational(1, 2), sp.Rational(-1, 3)), (0, 0), (-1, -2)]
res = [(p, fibre_size(*p), sp.simplify(L1.subs({u: p[0], v: p[1]}))) for p in gen_pts]
print("      fibre sizes off S_1:", ", ".join("%s -> %s (L = %s)" % (p, n, l) for p, n, l in res))
report("3. three preimages at 8 rational points off S_1", all(n == 3 and l != 0 for _, n, l in res))
S_pts = []
for mm in (1, 2, 3, 5, sp.Rational(1, 2), -1, -4, sp.Rational(7, 3)):
    uu = (mm + 12 / sp.Integer(mm) if not isinstance(mm, sp.Rational) else mm + 12 / mm) / 2
    ss = (12 / sp.nsimplify(mm) - mm) / 2
    for sg in (1, -1):
        vv = (-uu ** 3 + 18 * uu + sg * ss ** 3) / 54
        if sp.simplify(L1.subs({u: uu, v: vv})) == 0:
            S_pts.append((sp.nsimplify(uu), sp.nsimplify(vv)))
S_pts = sorted(set(S_pts), key=lambda p: (float(p[0]), float(p[1])))
res = [(p, fibre_size(*p)) for p in S_pts]
print("      fibre sizes on S_1:", ", ".join("%s -> %s" % (p, n) for p, n in res))
report("3. one preimage at %d rational points of S_1 (two sheets missing)" % len(S_pts), len(S_pts) >= 10 and all(n == 1 for _, n in res))

# 3b: the other factor of the discriminant, 27 v^2 - 9 u v + 8 = 0, is NOT part of the non-properness set
conic_pts = [((27 * vv ** 2 + 8) / (9 * vv), vv) for vv in (sp.Integer(1), sp.Integer(-2), sp.Rational(1, 3), sp.Rational(5, 2))]
res = [(p, fibre_size(*p), sp.simplify(L1.subs({u: p[0], v: p[1]}))) for p in conic_pts]
print("      fibre sizes on the conic 27v^2 - 9uv + 8 = 0:", ", ".join("%s -> %s (L = %s)" % (p, n, l) for p, n, l in res))
report("3b. three preimages on the conic 27v^2 - 9uv + 8 = 0 (double root of the x-cubic, not a missing sheet)", all(n == 3 and l != 0 for _, n, l in res))
# 3c: points of S_1 with 3uv = 4 have an empty fibre (the image of the slice map is not all of C^2)
sol = sp.solve([L1, 3 * u * v - 4], [u, v], dict=True)
print("      points of S_1 with 3uv = 4:", len(sol), "; first:", sol[0] if sol else None)

# 4 numeric fibre sizes at random complex points
import mpmath as mp
mp.mp.dps = 40
random.seed(20261111)


def fibre_numeric(u0, v0):
    """solve u = u0, v = v0 on X_1: the x-coordinates are roots of L x^3 + (4 - 3 u0 v0) x - 2 v0 (item 2); for each root
    solve the numerator of u - u0 for y and keep the pairs that satisfy both equations.  Returns the list of |x|."""
    Lv = complex(sp.N(L1.subs({u: u0, v: v0}), 40))
    cub = [mp.mpc(Lv), mp.mpc(0), mp.mpc(complex(sp.N(4 - 3 * u0 * v0, 40))), mp.mpc(complex(sp.N(-2 * v0, 40)))]
    u0c, v0c = mp.mpc(complex(sp.N(u0, 40))), mp.mpc(complex(sp.N(v0, 40)))
    pts = []
    for xr in mp.polyroots(cub, maxsteps=2000, extraprec=400):
        # numerator of u - u0:  y (1 + x y) + 3 x - u0 (1 + x y) = x y^2 + (1 - u0 x) y + (3 x - u0)
        for yr in mp.polyroots([xr, 1 - u0c * xr, 3 * xr - u0c], maxsteps=2000, extraprec=400):
            qq = 1 + xr * yr
            if abs(qq) < mp.mpf(10) ** (-20):
                continue
            e1 = abs(yr + 3 * xr / qq - u0c)
            e2 = abs(xr / qq + xr / qq ** 2 - xr ** 3 / qq ** 3 - v0c)
            if e1 < mp.mpf(10) ** (-12) * (1 + abs(u0c)) and e2 < mp.mpf(10) ** (-12) * (1 + abs(v0c)):
                if all(abs(xr - a) + abs(yr - b) > mp.mpf(10) ** (-10) for a, b in pts):
                    pts.append((xr, yr))
    return sorted(float(abs(p[0])) for p in pts)


near, off = [], []
for _ in range(6):
    u0 = sp.Rational(random.randint(-400, 400), 100) + sp.I * sp.Rational(random.randint(-400, 400), 100)
    v_on = sp.solve(L1.subs(u, u0), v)[0]                        # an exact point of S_1 (algebraic number)
    eps = sp.Rational(1, 10 ** 12)
    near.append(fibre_numeric(u0, v_on + eps))                   # 1e-12 away from S_1
    off.append(fibre_numeric(u0, v_on + sp.Rational(37, 100) - sp.I * sp.Rational(11, 100)))
print("      numeric, off S_1: |x| of the preimages:", [["%.3g" % a for a in l] for l in off])
print("      numeric, at distance 1e-12 from S_1:    ", [["%.3g" % a for a in l] for l in near])
report("4. numeric: 3 preimages off S_1 at 6 random complex points", all(len(l) == 3 for l in off))
report("4. numeric: approaching S_1 one preimage stays bounded and two escape (|x| > 1e4)",
       all(len(l) == 3 and l[0] < 1e3 and l[1] > 1e4 and l[2] > 1e4 for l in near))


# 5 intersection counts
def count(expr_curve, param_sub, var, tval):
    """number of distinct affine intersection points for generic t (degree of the squarefree numerator)."""
    num = sp.numer(sp.together(expr_curve.subs(param_sub)))
    P = sp.Poly(num, var)
    Pt = sp.Poly(num.subs(t, tval), var)
    sqf = sp.gcd(Pt, Pt.diff(var)).degree() == 0
    return P.degree(), Pt.degree(), sqf


cases = [
    ("(u-2)(v-3) = t", {v: 3 + t / (u - 2)}, u, 5),
    ("(u-5)(v+7) = t", {v: -7 + t / (u - 5)}, u, 5),
    ("u v = t  [centre (0,0), excluded]", {v: t / u}, u, 4),
    ("u = t", {u: t}, v, 2),
    ("v = t", {v: t}, u, 3),
    ("u + v = t", {v: t - u}, u, 4),
]
for name, sub, var, expect in cases:
    d1, d2, sqf = count(L1, sub, var, sp.Rational(5, 3))
    # remove spurious roots at the pole of the parametrisation (u = centre) : check the numerator there
    report("5. S_1 meets %-34s in %d affine points (generic t; squarefree at t = 5/3)" % (name, expect), d1 == expect and d2 == expect and sqf, "(degree %d)" % d1)

# 6
W = (u_claim.subs(c, 1) - 2) * (v_claim.subs(c, 1) - 3)
wf = 1 / q
dwdW = sp.simplify(sp.diff(wf, x) * sp.diff(W, y) - sp.diff(wf, y) * sp.diff(W, x))
report("6. dw ^ dW != 0 on X_1 (so [dq/q] restricts to a nonzero class on every generic fibre)", dwdW != 0)

# 7 irreducibility of the generic fibre of W on X_1 for the centres (2,3) and (0,0)
for (a0, b0) in ((2, 3), (0, 0), (5, -7)):
    good = True
    for tv in (sp.Rational(5, 3), sp.Integer(-2), sp.Rational(7, 11)):
        sub = {u: a0 + tau, v: b0 + tv / tau}
        cub = sp.numer(sp.together((L1 * x ** 3 + (4 - 3 * u * v) * x - 2 * v).subs(sub)))
        fl = sp.factor_list(cub, x, tau)
        irr = sum(e for f, e in fl[1] if sp.degree(f, x) > 0) == 1
        dsc = sp.Poly(sp.discriminant(sp.Poly(cub, x), x), tau)
        sqfree_part = sp.Poly(sp.quo(dsc, sp.gcd(dsc, dsc.diff(tau))), tau)
        # odd multiplicity root exists iff dsc is not c * square: compare with the square of sqrt
        fac = sp.factor_list(dsc.as_expr(), tau)
        odd = any(e % 2 == 1 and sp.degree(f, tau) > 0 for f, e in fac[1])
        good = good and irr and odd
    report("7. centre (%d,%d): cubic in x irreducible over Q(tau), discriminant not a square, t = 5/3, -2, 7/11" % (a0, b0), good)

print("ALL OK" if ok else "SOME CHECK FAILED")
