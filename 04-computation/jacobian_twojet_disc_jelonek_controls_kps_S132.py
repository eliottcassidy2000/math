# P1 controls: (1) d=3 law reproduction (THM-1310/THM-1335 objects),
# (2) tame d=4 control on the P2-refutation Keller map.
import sympy as sp

a, b, g, x, y, z, s = sp.symbols('a b g x y z s')
fails = []

def check(name, expr_zero):
    ok = sp.simplify(sp.expand(expr_zero)) == 0
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        fails.append(name)
        print("   residual:", sp.expand(expr_zero))

print("=" * 72)
print("PART 1 : d=3 control (THM-1310 fiber cubic)")
print("=" * 72)

L = 27*a**2*g**2 - 18*a*b*g + 16*a + b**3*g - b**2
Q = 27*a*g**2 - 9*b*g + 8
E = 54*a**2*g - 18*a*b + b**3
N = L*x**3 + (4 - 3*b*g)*x - 2*g

# 1a. re-derive N from the conic pair (THM-1310 SS1-2)
u = 1 + s
C1 = 3*a*x**2 - u*(b*x - s)
C2 = a*x**3 + g*u**3 - x*u*(u + 1)
R = sp.expand(sp.resultant(sp.Poly(C1, s), sp.Poly(C2, s)))
check("Res_s(C1,C2) == a * x^3 * N(x)", R - sp.expand(a*x**3*N))

# 1b. discriminant law
Delta = sp.expand(sp.discriminant(sp.Poly(N, x)))
check("disc_x(N) == -4 * Q^2 * L", Delta - sp.expand(-4*Q**2*L))
print("   disc_x(N) =", sp.factor(Delta))

# equivalent depressed-cubic form of the same law
check("(4-3bg)^3 + 27 g^2 L == Q^2",
      sp.expand((4 - 3*b*g)**3 + 27*g**2*L - Q**2))

# 1c. master identity (THM-1335)
check("108 a^2 L == (12a - b^2)^3 + E^2",
      108*a**2*L - ((12*a - b**2)**3 + E**2))

# 1d. shape summary: disc = -4 * (square) * (Jelonek {L=0}), L to first power
check("L is not a factor of Q (L appears to FIRST power in disc)",
      0 if sp.rem(sp.Poly(Q**2, a, b, g), sp.Poly(L, a, b, g)).as_expr() != 0 else 1)

print()
print("=" * 72)
print("PART 2 : tame d=4 control (P2-refutation Keller map, THM-2446 S5)")
print("=" * 72)

F1 = x*z**2 + y*z + x*y - x**3 + x
F2 = z**2 + 3*x**2 - 2*y
F3 = x + z

# 2a. Keller
J = sp.Matrix([[sp.diff(Fi, v) for v in (x, y, z)] for Fi in (F1, F2, F3)])
detJ = sp.expand(J.det())
print("   det J =", detJ)
check("det J == -2 (constant, Keller)", detJ + 2)

# 2b. explicit triangular inverse => fiber is a SINGLE point
GX = a + (b*g - g**3)/2          # x = a + (bg - g^3)/2
GZ = g - GX                      # z = g - x
GY = sp.expand((GZ**2 + 3*GX**2 - b)/2)   # y from F2 = b
print("   inverse: x =", sp.expand(GX))
print("            y =", GY)
print("            z =", sp.expand(GZ))

sub = {x: GX, y: GY, z: GZ}
check("F1(G) == a", sp.expand(F1.subs(sub)) - a)
check("F2(G) == b", sp.expand(F2.subs(sub)) - b)
check("F3(G) == g", sp.expand(F3.subs(sub)) - g)

sub2 = {a: F1, b: F2, g: F3}
check("G1(F) == x  (x in C[F1,F2,F3])", sp.expand(GX.subs(sub2)) - x)
check("G2(F) == y  (y in C[F1,F2,F3])", sp.expand(GY.subs(sub2)) - y)
check("G3(F) == z  (z in C[F1,F2,F3])", sp.expand(GZ.subs(sub2)) - z)

# 2c. the "fiber quartic" degenerates: eliminate y (linear in F2), then
#     use F3 to eliminate one variable; the eliminant is LINEAR, not quartic
ysol = sp.solve(F2 - b, y)[0]
P = sp.expand((F1 - a).subs(y, ysol))     # polynomial in x, z over target
print("   after y-elimination:  P(x,z) =", P)
print("      deg_z P =", sp.degree(P, z), " deg_x P =", sp.degree(P, x))

elim_x = sp.expand(P.subs(z, g - x))      # Res_z(P, x+z-g) up to sign
elim_z = sp.expand(P.subs(x, g - z))
print("   x-eliminant:", elim_x, "  deg_x =", sp.degree(elim_x, x))
print("   z-eliminant:", elim_z, "  deg_z =", sp.degree(elim_z, z))
check("x-eliminant is LINEAR in x (deg 1, field degree 1)",
      sp.degree(elim_x, x) - 1)
check("z-eliminant is LINEAR in z (deg 1, field degree 1)",
      sp.degree(elim_z, z) - 1)
check("x-eliminant == x - GX (monic linear, unique fiber point)",
      elim_x - (x - GX))

# cross-check via resultants (no hand substitution)
Rz = sp.expand(sp.resultant(sp.Poly(P, z), sp.Poly(x + z - g, z)))
check("Res_z agrees with substitution eliminant (up to sign)",
      0 if sp.simplify(Rz - elim_x) == 0 or sp.simplify(Rz + elim_x) == 0 else 1)

# 2d. non-properness (Jelonek) set is EMPTY:
#     x, y, z are POLYNOMIALS in F1,F2,F3 (checked in 2b), i.e. each satisfies
#     a MONIC DEGREE-1 integral equation over C[F]; hence C[x,y,z] = C[F],
#     F is finite, proper at every target, A_F = emptyset.
mono_x = sp.expand(x - (F1 + (F2*F3 - F3**3)/2))
check("monic deg-1 integral equation for x over C[F]", mono_x)

print()
if fails:
    print("FAILED:", fails)
    raise SystemExit(1)
print("ALL CHECKS PASSED")
