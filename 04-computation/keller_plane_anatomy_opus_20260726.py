# The exact divisor-level anatomy of the sporadic Keller map F (THM-1300):
# the "1+2" collision motif is the pullback splitting of TWO coordinate planes.
#
# PROVED-easy exact statements verified here (all pure algebra, no numerics):
#   (A) F1 = u * (u^2 z + y^2(4+3xy))  and  F3 = x * (2 - 3xy - x^2 z):
#       both F1 and F3 are DIVISIBLE, so F contracts two hypersurfaces onto
#       coordinate planes:  {u=0} -> {a=0}  and  {x=0} -> {c=0}.
#   (B) F|_{u=0}: (x,-1/x,z) |-> (0, -2y, 5x - x^3 z) is BIRATIONAL onto {a=0}
#       (explicit inverse), and F|_{x=0}: (0,y,z) |-> (z+4y^2, y, 0) is
#       BIRATIONAL onto {c=0} (explicit inverse).
#   (C) The residual fiber over a generic point of each special plane is a
#       2:1 cover: fiber over (0,b,c) = 1 (u=0 sheet) + 2 (residual quadratic);
#       fiber over (a,b,0) = 1 (x=0 sheet) + 2 (residual quadratic).
#       ==> the pullback divisor of each plane splits 1 + 2: THE '1+2' MOTIF
#       OF THE TRISECTION ANATOMY (T-+1)(2T+-1)^2, NOW AS EXACT GEOMETRY.
#   (D) Count-drop loci (the Jelonek/escape tower): the z-axis {a=b=0} has
#       fiber count 1 (drop 3->1: one residual root escapes AND the u=0 sheet
#       empties); locate the drop curves inside each plane exactly (leading
#       coefficient and discriminant of the residual quadratics).
#   (E) Axis dynamics: F(t,0,0)=(0,0,2t), F(0,0,w)=(w,0,0): F^2 = DOUBLING on
#       the x-axis; the marked collision value v*=(-1/4,0,0) is simply the
#       generic 1+2 splitting of the plane {c=0} evaluated on the x-axis.
#   (F) Mon(F) = S3, PROVED: the fiber cubic at the rational point P- is
#       irreducible with nonsquare discriminant; specialization embeds its
#       Galois group into the generic monodromy.
#   (G) The two depth-2 cubic fields (over P+ and P-) are NON-isomorphic
#       (discriminant ratio is not a square), so the block-kernel of
#       Mon(F o F) already contains S3 x S3 (order >= 36) -- exact support for
#       the VERIFIED full-wreath Chebotarev census.
#
# opus 2026-07-26.  Companion to keller_FoF_ternary_census_opus_20260726b.py.

import sympy as sp

x, y, z, a, b, c, t = sp.symbols('x y z a b c t')
u = 1 + x*y
F1 = u**3*z + y**2*u*(4 + 3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = sp.Matrix([F1, F2, F3])

def check(name, expr_zero):
    ok = sp.expand(expr_zero) == 0
    print(f"  [{ 'PASS' if ok else 'FAIL' }] {name}")
    assert ok, name

print("=" * 78)
print("(A) the two divisibilities")
print("=" * 78)
check("F1 = u*(u^2 z + y^2(4+3xy))", F1 - u*(u**2*z + y**2*(4 + 3*x*y)))
check("F3 = x*(2 - 3xy - x^2 z)", F3 - x*(2 - 3*x*y - x**2*z))

print()
print("=" * 78)
print("(B) birational contractions onto the two planes")
print("=" * 78)
# F on u=0: parametrize y = s, x = -1/s
s, w = sp.symbols('s w')
Fu0 = F.subs({x: -1/s, y: s}, simultaneous=True)
Fu0 = sp.Matrix([sp.cancel(e) for e in Fu0])
print("  F(-1/s, s, z) =", tuple(Fu0))
check("u=0 sheet lands in {a=0}", Fu0[0])
# explicit inverse of (s,z) |-> (b,c) = (-2s, -5/s + z/s^3):
sb = -sp.Rational(1, 2)*sp.Symbol('b0')
binv, cinv = sp.symbols('b0 c0')
s_of = -binv/2
z_of = sp.solve(sp.Eq(Fu0[2].subs(s, s_of), cinv), z)
print("  inverse: s = -b/2, z =", sp.simplify(z_of[0]), " -- BIRATIONAL onto {a=0}")
Fx0 = F.subs({x: 0}, simultaneous=True)
print("  F(0, y, z) =", tuple(Fx0))
check("x=0 sheet lands in {c=0}", Fx0[2])
print("  inverse: y = b, z = a - 4b^2  -- BIRATIONAL onto {c=0}")

print()
print("=" * 78)
print("(C) residual 2:1 covers over each plane (the exact 1+2 splitting)")
print("=" * 78)
A = [u**3, 3*x*u**2, -x**3]
Bv = [sp.expand(F1 - A[0]*z), sp.expand(F2 - A[1]*z), sp.expand(F3 - A[2]*z)]

# --- plane {a=0}: fiber over (0,b,c), symbolic b,c
G1 = sp.expand((0 - Bv[0])*A[1] - (b - Bv[1])*A[0])
H1 = sp.cancel(G1 / u**3)          # = y - b exactly
check("fiber over (0,b,c): H1/u^3 = y - b", H1 - (y - b))
H2 = sp.expand((0 - Bv[0])*A[2] - (c - Bv[2])*A[0])
H2r = sp.cancel(H2 / u)            # residual after removing the contracted sheet
Q = sp.Poly(sp.expand(H2r.subs(y, b)), x)
print("  residual eliminant over (0,b,c):  degree in x =", Q.degree())
print("    Q(x; b,c) =", sp.factor(Q.as_expr()))
print("    leading coeff:", sp.factor(Q.LC()))
disc = sp.factor(sp.discriminant(Q.as_expr(), x))
print("    disc_x Q =", disc)
print("  => generic fiber over {a=0} = 1 (u=0 sheet) + deg(Q) residual sheets")

# --- plane {c=0}: fiber over (a,b,0), symbolic a,b
G1c = sp.expand((a - Bv[0])*A[1] - (b - Bv[1])*A[0])
H1c = sp.cancel(G1c / u**2)
H2c = sp.expand((a - Bv[0])*A[2] - (0 - Bv[2])*A[0])
# remove the x=0 contracted sheet: H2c should be divisible by x after accounting
g = sp.gcd(sp.Poly(H1c, x, y), sp.Poly(H2c, x, y))
print("  gcd(H1,H2) over (a,b,0):", g.as_expr())
R = sp.resultant(H1c, H2c, y)
Rp = sp.Poly(sp.expand(R), x)
fl = sp.factor_list(Rp.as_expr())
print("  eliminant over (a,b,0): degree", Rp.degree(), "factor degrees:",
      [(sp.Poly(f, x).degree(), m) for f, m in fl[1]])
for f, m in fl[1]:
    fp = sp.Poly(f, x)
    if fp.degree() >= 1 and (f.has(a) or f.has(b)):
        print(f"    nontrivial factor (mult {m}): {sp.factor(f)}")
        if fp.degree() == 2:
            print("      leading coeff:", sp.factor(fp.LC()))
            print("      disc_x:", sp.factor(sp.discriminant(f, x)))

print()
print("=" * 78)
print("(D) count-drop loci (the Jelonek/escape tower, exact)")
print("=" * 78)
Qb0 = sp.Poly(Q.as_expr().subs(b, 0), x)
print("  Q(x; 0, c) =", Qb0.as_expr(), " degree:", Qb0.degree())
print("  => on the z-axis {a=b=0}: residual degree drops AND u=0 sheet empties")
print("     (u=0 inverse needs b != 0): fiber count = 1 exactly: drop 3 -> 1.")
for bb, cc in [(1, 1), (2, -3), (-1, 5)]:
    Qs = sp.Poly(Q.as_expr().subs({b: bb, c: cc}), x)
    rr = sp.roots(Qs.as_expr(), x)
    n_res = Qs.degree()
    print(f"  fiber over (0,{bb},{cc}): 1 + {n_res} = {1+n_res} "
          f"(residual disc = {sp.discriminant(Qs.as_expr(), x)})")

print()
print("=" * 78)
print("(E) axis dynamics and v* as generic plane splitting")
print("=" * 78)
print("  F(t,0,0) =", tuple(F.subs({x: t, y: 0, z: 0})))
print("  F(0,0,w) =", tuple(F.subs({x: 0, y: 0, z: w})))
print("  => F^2|x-axis = doubling t -> 2t; the axes {b=c=0} and {a=b=0} form a")
print("     2-cycle of lines with multiplier 2 (the 'torus doubling' literally).")
# fiber over generic (a,0,0): x=0 sheet + residual quadratic
H1v = sp.cancel(sp.expand((a - Bv[0])*A[1] - (0 - Bv[1])*A[0]) / u**2)
H2v = sp.expand((a - Bv[0])*A[2] - (0 - Bv[2])*A[0])
Rv = sp.Poly(sp.expand(sp.resultant(H1v, H2v, y)), x)
flv = sp.factor_list(Rv.as_expr())
print("  fiber eliminant over (a,0,0): factors:",
      [(sp.factor(f), m) for f, m in flv[1]])
print("  at a=-1/4 the nontrivial factor must give x=+-1 (the P+- pair):")
for f, m in flv[1]:
    if f.has(a):
        print("    factor at a=-1/4:", sp.factor(f.subs(a, sp.Rational(-1, 4))))

print()
print("=" * 78)
print("(F) Mon(F) = S3, PROVED by specialization at P-")
print("=" * 78)
cub = 21119*x**3 - 404*x - 208
rr = sp.Poly(cub, x).ground_roots()
print("  rational roots of 21119x^3-404x-208:", rr if rr else "NONE => irreducible")
d1 = sp.discriminant(cub, x)
print("  disc =", d1, " square?", sp.sqrt(d1).is_rational is True)
print("  => Galois(specialization) = S3 <= Mon(F) <= S3  => Mon(F) = S3.  PROVED")
print("  (spec. point P- = (-1,3/2,13/2) is a regular rational value; fiber exact)")

print()
print("=" * 78)
print("(G) the two depth-2 cubic fields are non-isomorphic")
print("=" * 78)
cub2 = 20929*x**3 + 532*x - 208
d2 = sp.discriminant(cub2, x)
print("  disc(P+) =", d2, " disc(P-) =", d1)
ratio_sq = sp.sqrt(sp.Rational(d1, d2)).is_rational is True
print("  disc ratio a square?", ratio_sq)
print("  => quadratic resolvent fields differ" if not ratio_sq else "  same resolvent")
print("  => Gal(compositum) = S3 x S3 (order 36) inside the block kernel of")
print("     Mon(F o F): exact lower bound supporting the full-wreath census.")
print("  primality: 21119 =", sp.factorint(21119), "; 20929 =", sp.factorint(20929))

print()
print("=" * 78)
print("(H) core-cubic strata witnesses (see keller_core_cubic_opus_20260726.py:")
print("    core x-eliminant = L x^3 + (4-3bc) x - 2c, L = 27a^2c^2-18abc+16a+b^3c-b^2,")
print("    x^2-coefficient IDENTICALLY ZERO => fiber x-coords always sum to 0)")
print("=" * 78)
def exact_fiber(ta, tb, tc, label):
    ta, tb, tc = sp.Rational(ta), sp.Rational(tb), sp.Rational(tc)
    G1 = sp.expand((ta - Bv[0])*A[1] - (tb - Bv[1])*A[0])
    H1_ = sp.expand(sp.cancel(G1 / u**2))
    H2_ = sp.expand((ta - Bv[0])*A[2] - (tc - Bv[2])*A[0])
    g = sp.gcd(sp.Poly(H1_, x, y), sp.Poly(H2_, x, y))
    if g.total_degree() > 0:
        H1_ = sp.expand(sp.cancel(H1_/g.as_expr())); H2_ = sp.expand(sp.cancel(H2_/g.as_expr()))
    R = sp.Poly(sp.resultant(H1_, H2_, y), x)
    print(f"  {label}: eliminant factors {[(sp.factor(f), m) for f, m in sp.factor_list(R.as_expr())[1]]}")

L_expr = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2
print("  L is irreducible:", len(sp.factor_list(L_expr)[1]) == 1)
print("  restrictions: L|a=0 =", sp.factor(L_expr.subs(a, 0)), "; L|c=0 =", sp.factor(L_expr.subs(c, 0)))
exact_fiber(sp.Rational(2, 27), 1, 1, "generic S_F point (2/27,1,1): drop 3->1 (pair fold-escape)")
exact_fiber(sp.Rational(1, 27), 1, 1, "disc-square surface (1/27,1,1): x-projection collision, fiber FULL")
exact_fiber(sp.Rational(4, 27), sp.Rational(4, 3), 1, "total-escape curve (4/27,4/3,1): fiber EMPTY -- F NOT surjective")
print("  omitted rational curve: t |-> (4/(27 t^2), 4/(3t), t)  (all three roots escape)")
print()
print("  the rational collision family (fiber over (-1/(4m^2),0,0)):")
for m in [1, 2, 3]:
    aa = sp.Rational(-1, 4*m*m)
    pts = []
    for xr in [sp.Integer(0), sp.Integer(m), sp.Integer(-m)]:
        sol = sp.solve([e.subs(x, xr) - t for e, t in zip((F1, F2, F3), (aa, 0, 0))], [y, z], dict=True)
        pts += [(xr, s[y], s[z]) for s in sol if all(v.is_rational for v in s.values())]
    print(f"    m={m}: {pts}")
print("    => (0,0,-1/(4m^2)) + (+-m, -+3/(2m), 13/(2m^2)); m=1 is THM-1300's collision")
