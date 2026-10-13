"""Audit C of THM-4625 (independent auditor, 2026-10-12): the two named instances (THM-1300, d = 2; Gao's G, d = 3) at random rational base points.

  (a) both maps are polynomial Keller maps of the displayed shape (the session script's 'polynomial: False' for G3 is a denominator-8 artifact);
  (b) (F1,F2)-fibre: points with b x = 2 EXIST (p(w) = 0 there, F3 != 0): the chart of (3) with 'b x != 2' misses them -> the stated isomorphism
      must exclude them explicitly; the points with F3 = 0 are counted (finite);
  (c) (F1,F3)-fibre: with ac != 0 but c^2 a a critical value of q the puncture count drops below d+3 (THM-1300: c^2 a = 4/27);
  (d) chart <-> fibre bijection away from {F3 = 0} u {bx = 2} (random points), and K(x,y,z) = K(F)(w) identities;
  (e) special fibre: the (F2,F3)-fibre over (b, 0) contains the line {x = 0, y = b} (THM-1300) / {x = 0, y = b/2} (G): '(4) no coordinate-pair
      fibre ... union of lines' must be read for GENERIC fibres;
  (f) genus of the (F1,F2)-fibre of G is 0 (h quadratic), the session script's 'genus 2' is wrong.
"""
import sympy as sp, random
x, y, z, w, s = sp.symbols('x y z w s')
random.seed(20261012)

def family(p_expr, g0, al, be):
    q_expr = sp.expand(sp.integrate((s / 2) * sp.diff(p_expr, w).subs(w, s), (s, 0, w)))
    gam = g0 + al * x * y + be * x**2 * z
    u = 1 + x * y
    W = sp.expand(gam * u)
    F3 = sp.expand(gam * x)
    F2 = sp.expand(sp.cancel((p_expr.subs(w, W) + 2 * gam) / (gam * x)))
    F1 = sp.expand(sp.cancel((q_expr.subs(w, W) + gam * W) / (gam * x)**2))
    return (F1, F2, F3), p_expr, q_expr, gam, W

inst = {
    'THM-1300 (d=2)': family(-3 * w**2 + 4 * w, 2, -3, -1),
    "Gao G (d=3)": family(w**3 - 6 * w**2 + 6 * w, 2, -4, -1),
}
for name, (F, p, q, gam, W) in inst.items():
    F1, F2, F3 = F
    print("=" * 100); print(name, ": p =", p, " q =", q, " gamma =", gam)
    polys = [sp.Poly(sp.fraction(sp.together(Fi))[1], x, y, z).total_degree() == 0 for Fi in F]
    J = sp.Matrix([[sp.diff(Fi, v) for v in (x, y, z)] for Fi in F])
    print("(a) polynomial components:", polys, " det JF =", sp.factor(J.det()), " degrees", [sp.Poly(Fi, x, y, z).total_degree() for Fi in F])
    # K(x,y,z) = K(F)(w): x = 2F3/(F2F3 - p(w)), gamma = (F2F3 - p(w))/2
    print("    identity x = 2F3/(F2F3 - p(w)) in K(x,y,z):", sp.cancel(x - 2 * F3 / (F2 * F3 - p.subs(w, W))) == 0,
          "; gamma = (F2F3 - p(w))/2:", sp.expand(gam - (F2 * F3 - p.subs(w, W)) / 2) == 0)
    # (b) points of the (F1,F2)-fibre with b x = 2
    av, bv = sp.Rational(random.randint(2, 9), random.randint(1, 5)), sp.Rational(random.randint(2, 9), random.randint(1, 5))
    print("(b) base point (a, b) = (%s, %s)" % (av, bv))
    e1 = sp.expand(F1.subs(x, 2 / bv) - av); e2 = sp.expand(F2.subs(x, 2 / bv) - bv)
    res = sp.resultant(sp.Poly(e1, y, z), sp.Poly(e2, y, z), z)
    pts = []
    for yv in [sp.N(r, 40) for r in sp.Poly(res, y).all_roots()]:
        for zv in sp.Poly(sp.expand(e2.subs(y, yv)), z).nroots(n=30, maxsteps=500):
            if abs(sp.N(e1.subs({y: yv, z: zv}), 40)) < 1e-25:
                pts.append((2 / bv, yv, zv))
    print("    fibre points with b x = 2 (numeric, 40 digits):", len(pts))
    for pt in pts:
        sub = dict(zip((x, y, z), pt))
        vals = [sp.N(Fi.subs(sub), 40) for Fi in F]
        wv = sp.N(W.subs(sub), 40)
        print("      (x,y,z) ~", tuple(sp.N(v, 8) for v in pt), " F ~", tuple(sp.N(v, 8) for v in vals), " w ~", sp.N(wv, 8), " |p(w)| =", sp.N(abs(p.subs(w, wv)), 3), " |F3| =", sp.N(abs(vals[2]), 6))
    if 'THM' in name:   # exact rational instance: (a, b) = (2/3, 2) gives y in {0, 1}
        for yv in (0, 1):
            sub = {x: 1, y: yv}
            zz = sp.solve(F2.subs(sub) - 2, z)
            for zv in zz:
                sub2 = {x: 1, y: yv, z: zv}
                print("      exact (a,b) = (2/3, 2): point (1, %s, %s): F =" % (yv, zv), [sp.simplify(Fi.subs(sub2)) for Fi in F], " w =", W.subs(sub2), " p(w) =", p.subs(w, W.subs(sub2)))
    # points with F3 = 0 on the (F1,F2)-fibre: x = 0 branch and gamma = 0 branch
    s0 = sp.solve([F1.subs(x, 0) - av, F2.subs(x, 0) - bv], [y, z], dict=True)
    print("    F3 = 0 points on the fibre: x = 0 branch:", len(s0), s0)
    zg = sp.solve(gam, z)[0]                      # gamma = 0, x != 0
    g1 = sp.numer(sp.together(F1.subs(z, zg) - av)); g2 = sp.numer(sp.together(F2.subs(z, zg) - bv))
    rg = sp.factor(sp.resultant(sp.Poly(g1, x, y), sp.Poly(g2, x, y), y))
    print("    gamma = 0 branch: resultant in y =", rg)
    Rx = sp.Poly(sp.resultant(sp.Poly(g1, x, y), sp.Poly(g2, x, y), y), x)
    nx = len([r for r in Rx.all_roots() if r != 0])
    # bound: each x-root admits at most deg_y(g2) values of y; count exactly via gcd over the splitting field when cheap
    cnt = 0
    try:
        for r in Rx.all_roots():
            if r == 0: continue
            Kf = sp.QQ.algebraic_field(r)
            gg = sp.gcd(sp.Poly(g1.subs(x, r), y, domain=Kf), sp.Poly(g2.subs(x, r), y, domain=Kf))
            cnt += gg.degree()
        exact = True
    except Exception as ex:
        exact = False; cnt = nx * sp.Poly(g2, y).degree()
    print("    gamma = 0 branch: resultant nonzero => finitely many; x-values (x != 0):", nx, "; solutions counted exactly over Q(x-root):" if exact else "; crude upper bound:", cnt, " => F3 = 0 points on the fibre: finite")
    # (c) (F1,F3)-fibre puncture count: critical values of q are q(rho), p'(rho) = 0
    rhos = sp.Poly(sp.diff(p, w), w).all_roots()
    print("(c) critical points of q other than 0:", rhos, " critical values q(rho):", [sp.simplify(q.subs(w, r)) for r in rhos])
    cv = sp.Rational(3, 2)
    for aval in [sp.simplify(q.subs(w, rhos[0])) / cv**2, av]:
        f = sp.expand(q - cv**2 * aval)
        Kf = sp.QQ.algebraic_field(sp.sqrt(2)) if f.has(sp.sqrt(2)) else sp.QQ
        fp = sp.Poly(f, w, domain=Kf)
        ndist = fp.degree() - sp.gcd(fp, fp.diff(w)).degree()
        print("    c = %s, a = %s: q(w) = c^2 a has %d distinct roots (degree %d) -> punctures 0, infinity + roots = %d (d + 3 = %d)" % (cv, aval, ndist, fp.degree(), 2 + ndist, sp.Poly(p, w).degree() + 3))
    # (d) chart <-> fibre: random w on the chart, both quadratic roots x, mapped back
    d = sp.Poly(p, w).degree()
    P = sp.expand(sp.integrate(p.subs(w, s), (s, 0, w)))
    At = sp.expand((q * bv**2 - av * p**2) / w**2); Bt = sp.expand((w * p - 4 * q) * bv / w**2); Ct = sp.expand(-2 * P / w**2)
    ok = True; n = 0
    for wv in [sp.Rational(k, 3) for k in (1, 2, 5, 7, 11)]:
        for xv in sp.Poly((At * x**2 + Bt * x + Ct).subs(w, wv), x).all_roots():
            if xv == 0 or bv * xv == 2 or p.subs(w, wv) == 0:
                continue
            gv = p.subs(w, wv) / (bv * xv - 2); yv = (wv / gv - 1) / xv; zv = (gv - 2 - (-3 if 'THM' in name else -4) * xv * yv) / (-xv**2)
            sub = {x: xv, y: yv, z: zv}
            ok &= (sp.simplify(F1.subs(sub) - av) == 0 and sp.simplify(F2.subs(sub) - bv) == 0 and sp.simplify(F3.subs(sub) - gv * xv) == 0); n += 1
    print("(d) chart -> fibre map verified on %d points: %s" % (n, ok))
    h = sp.expand(bv**2 - 8 * av * P / w**2)
    print("(f) h = b^2 - 8aP/w^2 =", h, " degree", sp.Poly(h, w).degree(), " squarefree:", sp.gcd(sp.Poly(h, w), sp.Poly(sp.diff(h, w), w)).degree() == 0, " => genus", (sp.Poly(h, w).degree() - 1) // 2)
    # (e) special fibre over (b, c) = (b, 0) of (F2,F3): factor F2 restricted to x = 0, and the gamma = 0 part
    print("(e) (F2,F3)-fibre over (b, 0): F3 = x*(%s); on x = 0, F2 =" % sp.factor(F3 / x), sp.factor(F2.subs(x, 0)), " => contains the line {x = 0, %s = b} (z free)" % sp.factor(F2.subs(x, 0)))
    print("    on gamma = 0 (x != 0): F2 - b =", sp.factor(sp.numer(sp.together(F2.subs(z, zg) - sp.Symbol('b')))), " (not a line unless this factors linearly in a coordinate)")
