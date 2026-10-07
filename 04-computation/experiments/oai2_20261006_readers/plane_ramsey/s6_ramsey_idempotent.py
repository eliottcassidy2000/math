# For a spherical configuration with coordinates in a NUMBER FIELD F = Q(theta), the #172 tensor criterion
#   (p_i (x) 1)^T P (1 (x) p_i) = 0 in B = F (x)_Q F,   m_F(P_spatial) = I_d
# is solved by P = e * (H (x) 1), e = separability idempotent, H = sphere matrix.  Check exactly in
# B = Q[x,y]/(f(x), f(y)) for (i) the 18-point "three hexagons rotated by 0, +-theta_Moser" set
# (the M_L-type unit circle with ±w3-rotations, cos theta = 5/6), (ii) the algebraic twin of #172's
# 12-point three-squares example (Liouville angle -> Moser angle).
import sympy as sp
x, y, X = sp.symbols('x y X')

def reduce2(expr, f):
    e = sp.Poly(sp.expand(expr), x, y)
    # reduce modulo f(x) and f(y) (f monic)
    fx = sp.Poly(f.subs(X, x), x, y); fy = sp.Poly(f.subs(X, y), x, y)
    _, r = sp.reduced(e.as_expr(), [fx.as_expr(), fy.as_expr()], x, y)
    return sp.expand(r)

def check(name, f, pts_theta, H):
    """pts_theta: list of (X_i(theta), Y_i(theta)) as polynomials in X=theta; H 3x3 rational sphere matrix."""
    fprime = sp.diff(f, X)
    h = sp.cancel((f.subs(X, x) - f.subs(X, y)) / (x - y))
    # inverse of f'(x) in Q[x]/(f(x))
    inv = sp.invert(sp.Poly(fprime.subs(X, x), x), sp.Poly(f.subs(X, x), x)).as_expr()
    e = reduce2(h*inv, f)
    m_e = sp.simplify(sp.rem(sp.Poly(e.subs(y, x), x), sp.Poly(f.subs(X, x), x)).as_expr())
    ok_all = True; ok_without_e = True
    for (Xi, Yi) in pts_theta:
        px = [1, Xi.subs(X, x), Yi.subs(X, x)]; py = [1, Xi.subs(X, y), Yi.subs(X, y)]
        raw = sum(px[a]*H[a][b]*py[b] for a in range(3) for b in range(3))
        r_raw = reduce2(raw, f)
        r = reduce2(e*raw, f)
        ok_all &= (r == 0)
        ok_without_e &= (r_raw == 0)
    print(f"[{name}] m_F(e) = {m_e};  e*(tensor evaluation) == 0 for all {len(pts_theta)} points: {ok_all};"
          f"  evaluations vanish WITHOUT e: {ok_without_e}")

# F = Q(sqrt3, sqrt11), theta = sqrt3 + sqrt11, f = X^4 - 28 X^2 + 64
f = X**4 - 28*X**2 + 64
th = sp.sqrt(3) + sp.sqrt(11)
assert sp.simplify(f.subs(X, th)) == 0
s3 = (X**3 - 20*X)/16; s11 = (36*X - X**3)/16
assert sp.simplify(s3.subs(X, th) - sp.sqrt(3)) == 0 and sp.simplify(s11.subs(X, th) - sp.sqrt(11)) == 0
# unit vectors zeta6^j * w^eps, w = (5 + i sqrt11)/6 : real/imag parts as polys in theta
import itertools
c60 = [sp.Rational(1), sp.Rational(1, 2), sp.Rational(-1, 2), sp.Rational(-1), sp.Rational(-1, 2), sp.Rational(1, 2)]
s60 = [0, s3/2, s3/2, 0, -s3/2, -s3/2]
wc, ws = sp.Rational(5, 6), s11/6
pts = []
for j in range(6):
    for eps in (0, 1, -1):
        if eps == 0: re, im = c60[j], s60[j]
        else:
            re = c60[j]*wc - s60[j]*eps*ws
            im = s60[j]*wc + c60[j]*eps*ws
        pts.append((sp.expand(sp.sympify(re)), sp.expand(sp.sympify(im))))
H = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]
check("18 pts: three hexagons rotated by 0, +-arccos(5/6) on the unit circle, F=Q(sqrt3,sqrt11)", f, pts, H)
# three squares rotated by 0, +-theta_Moser  (algebraic twin of #172 Prop. 7.7, which uses a Liouville angle)
pts2 = []
for (qx, qy) in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
    for eps in (0, 1, -1):
        re = qx*wc - qy*eps*ws if eps else sp.Integer(qx)
        im = qy*wc + qx*eps*ws if eps else sp.Integer(qy)
        pts2.append((sp.expand(sp.sympify(re)), sp.expand(sp.sympify(im))))
check("12 pts: three squares rotated by 0, +-arccos(5/6)", f, pts2, H)
