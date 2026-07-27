"""THM-2546 z-integrality: from VERIFIED (5 sampled targets) to PROVED.

Object: the sporadic Keller map (THM-1300/THM-2473)

    F(x,y,z) = ( u^3 z + y^2 u (4+3xy),
                 y + 3x u^2 z + 3x y^2 (4+3xy),
                 2x - 3x^2 y - x^3 z ),          u = 1+xy,

target coordinates (a,b,c),

    L  = 27a^2c^2 - 18abc + 16a + b^3c - b^2,
    p  = 4 - 3bc,      q = b^2 - 12a,
    E_x(x) = L x^3 + p x - 2c                    (x-core, THM-2473, PROVED).

CLAIM to prove: the z-eliminant core (the cubic whose roots are the three
fibre z-values over a generic target) has CONSTANT leading coefficient
(a 2-power), i.e. z is integral over Q[a,b,c].

STRATEGY (two independent exact routes + one standalone certificate).

(I)  Rational fibre sections.  On every fibre, r = b-y satisfies the two
     PROVED identities P_r(r)=0 and x D_r(r)=r (THM-2546 (2)-(3); both
     re-certified below).  Reducing x^2 P_r(r) modulo the quadratic
     x D_r(r)-r yields a relation LINEAR in r:

         (q x^2 - b x - 2) r + 3ax (9acx - bx + 2) = 0   on every fibre,

     an exact polynomial identity in Q[x,y,z] after the substitution
     (a,b,c,r) -> (F1,F2,F3,F2-y) [certified below].  Hence, wherever
     Dq(x) := q x^2 - b x - 2 does not vanish,

         r(x) = -3ax (9acx - bx + 2) / Dq(x),   y(x) = b - r(x),

     and from the third component (z = (2x - c - 3x^2 y)/x^3),

         z(x) = N(x) / D(x),
         D(x) = x^3 Dq(x),
         N(x) = (2x - c - 3bx^2) Dq(x) - 9ax^3 (9acx - bx + 2),

     again an exact identity z*D - N = 0 on the fibre graph [certified].

(II) The z-eliminant as a single small resultant.  Define

         P(Z) := Res_x( E_x(x), D(x) Z - N(x) )  in  Z[a,b,c][Z].

     Poisson's product formula over K = Q(a,b,c)-bar gives, with
     x1,x2,x3 the roots of E_x,

         P(Z) = L^5 prod_i ( D(x_i) Z - N(x_i) )
              = [Res_x(E_x, D)] * prod_i ( Z - z(x_i) )  =: rho * chi(Z),

     valid wherever rho := Res_x(E_x,D) != 0 (rho is a NONZERO polynomial,
     computed below, so this holds generically, hence as an identity of
     rational functions).  Genuineness: off {L=0} the fibre has exactly 3
     points (THM-2473 (3)); each fibre x is a root of E_x (identity I_x);
     E_x is generically separable (disc_x = -4(27ac^2-9bc+8)^2 L != 0);
     and two fibre points cannot share an x where Dq(x) != 0 and x != 0
     (both y and z are then the rational functions of x above) -- the
     nonvanishing of rho = Res(E_x, x^3 Dq) excludes both degeneracies at
     all three roots simultaneously.  So generically the fibre is exactly
     {(x_i, y(x_i), z(x_i))} and chi(Z) = prod (Z - z_i) IS the monic
     z-eliminant.  The primitive integer core is P / cont(P).

(III) Degree-bounded exact interpolation (independent re-derivation of P,
     the method requested for the upgrade).  P is DEFINED by the 8x8
     Sylvester determinant of E_x (deg_x 3) and G := D Z - N (deg_x 5):
     5 rows of E_x-coefficients, 3 rows of G-coefficients.  By Leibniz,
     every term of the determinant is a product of FIVE E_x-coefficients
     and THREE G-coefficients.  The exact per-variable coefficient-degree
     tables (printed below) are

         E_x coeffs [L, 0, p, -2c]:      deg_a <= 2, deg_b <= 3, deg_c <= 2
         G   coeffs [qZ, -bZ-N4, ...]:   deg_a <= 2, deg_b <= 3, deg_c <= 1
                                          deg_Z <= 1

     hence the a-priori bounds

         deg_a P <= 5*2 + 3*2 = 16,
         deg_b P <= 5*3 + 3*3 = 24,
         deg_c P <= 5*2 + 3*1 = 13,
         deg_Z P <= 3*1      =  3.

     A determinant is a polynomial identity in its entries, so the integer
     value of the SAME 8x8 determinant at ANY integer point (a0,b0,c0,Z0)
     equals the value of the polynomial P there -- no genericity needed,
     no specialization caveat.  Interpolating on a product grid of sizes
     17 x 25 x 14 (strictly exceeding each degree bound) with exact
     rational arithmetic therefore RECOVERS P exactly; this is a proof.
     We then verify the interpolant against all 5950 grid nodes, against
     25 random off-grid integer targets, and against the independent
     symbolic resultant of route (II).

(IV) Standalone integrality certificate.  Whatever route produced it, the
     primitive core C(a,b,c;Z) is finally certified by the single
     polynomial identity  C(F1,F2,F3; z) = 0  in Q[x,y,z].  Since C has
     Z-degree 3 with CONSTANT lead, z is integral over Q[a,b,c]:
     monic cubic = C/lead.  This certificate alone upgrades THM-2546's
     z-claim to PROVED; routes (II) and (III) additionally prove that C
     is the genuine z-eliminant (its roots ARE the fibre z-values).

opus 2026-07-27.  Companion to THM-2473 and THM-2546.
"""

import random
import time
from fractions import Fraction

import sympy as sp

T0 = time.time()

x, y, z, a, b, c, r, Z, t = sp.symbols("x y z a b c r Z t")
u = 1 + x * y
F1 = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
F2 = sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
p = 4 - 3 * b * c
q = b**2 - 12 * a
Ex = L * x**3 + p * x - 2 * c
Sq = 27 * a * c**2 - 9 * b * c + 8  # the x-collision square factor

Pr = 27 * a**2 * c - 18 * a * r + 3 * b * r**2 - 2 * r**3
Dr = r**2 - b * r + 3 * a

SUB = {a: F1, b: F2, c: F3, r: F2 - y}


def on_fibre_zero(expr):
    """Assert expr vanishes identically after (a,b,c,r)->(F1,F2,F3,F2-y)."""
    assert sp.expand(expr.subs(SUB, simultaneous=True)) == 0


def banner(s):
    print("=" * 78)
    print(s)
    print("=" * 78)


banner("[1] Fibre identities (exact, in Q[x,y,z]) and the linear-in-r reduction")
on_fibre_zero(Ex)
print("  I_x : E_x = L x^3 + p x - 2c vanishes on every fibre        [PASS]")
on_fibre_zero(Pr)
on_fibre_zero(x * Dr - r)
print("  I_r : P_r(r)=0 and x D_r(r)=r on every fibre (THM-2546)     [PASS]")

# derivation of the linear relation: x^2 P_r reduced mod (x D_r - r) in r
quo, rem = sp.div(sp.expand(x**2 * Pr), sp.expand(x * Dr - r), r)
Gx = 9 * a * c * x - b * x + 2
Dq = q * x**2 - b * x - 2
lin = Dq * r + 3 * a * x * Gx
assert sp.expand(rem - lin) == 0
assert sp.expand(x**2 * Pr - quo * (x * Dr - r) - lin) == 0
on_fibre_zero(lin)
print("  I_lin: x^2 P_r = quo*(x D_r - r) + [Dq(x) r + 3ax Gx],")
print("         Dq = q x^2 - b x - 2,  Gx = 9acx - bx + 2,")
print("         so  Dq(x) r + 3ax Gx = 0 on every fibre              [PASS]")

N = sp.expand((2 * x - c - 3 * b * x**2) * Dq - 9 * a * x**3 * Gx)
D = sp.expand(x**3 * Dq)
on_fibre_zero(z * D - N)
print("  I_z : z * D(x) - N(x) = 0 on every fibre,")
print("        D = x^3 Dq,  N = (2x - c - 3bx^2) Dq - 9ax^3 Gx       [PASS]")
print("  => wherever Dq(x) != 0, x != 0:  y = b + 3ax Gx / Dq,  z = N/D.")

disc_x = sp.factor(sp.discriminant(Ex, x))
assert sp.expand(disc_x + 4 * Sq**2 * L) == 0
disc_r = sp.factor(sp.discriminant(Pr, r))
assert disc_r == -2916 * a**2 * L
print("  disc_x(E_x) = -4 (27ac^2-9bc+8)^2 L   (THM-2473)            [PASS]")
print("  disc_r(P_r) = -2916 a^2 L             (THM-2546)            [PASS]")
print()

banner("[2] Route II: symbolic z-eliminant P(Z) = Res_x(E_x, D Z - N)")
G = sp.expand(D * Z - N)
t1 = time.time()
P_sym = sp.expand(sp.resultant(Ex, G, x))
print(f"  resultant computed in {time.time() - t1:.2f}s")
P_poly = sp.Poly(P_sym, Z)
assert P_poly.degree() == 3
rho = sp.expand(sp.resultant(Ex, D, x))
assert sp.expand(P_poly.LC() - rho) == 0
rho_f = sp.factor(rho)
print("  deg_Z P = 3;   lead lc_Z(P) = rho = Res_x(E_x, D) =")
print("      rho =", rho_f)
R2 = sp.factor(sp.resultant(Ex, Dq, x))
print("  and Res_x(E_x, Dq) =", R2)
print("  (the x-collision square factor (27ac^2-9bc+8)^2 of disc_x reappears)")
assert sp.simplify(rho) != 0
print("  rho is a nonzero polynomial => Poisson: P = rho * prod(Z - z_i)")

# primitive core
cont = 0
for cf in P_poly.all_coeffs():
    cont = sp.gcd(cont, cf)
core = sp.expand(sp.cancel(P_sym / cont))
core_poly = sp.Poly(core, Z)
if core_poly.LC().is_number and core_poly.LC() < 0:
    core = sp.expand(-core)
    core_poly = sp.Poly(core, Z)
lead = core_poly.LC()
print("  content(P) =", sp.factor(cont), " ( = rho/8 )")
assert sp.expand(cont - rho / 8) == 0
print()
print("  PRIMITIVE z-ELIMINANT CORE  (the z-analogue of E_x and the y-cubic):")
cAz = sp.expand(core_poly.nth(2))
cBz6L = sp.expand(core_poly.nth(1))
cCzL = sp.expand(core_poly.nth(0))
print("    C(Z) = 8 Z^3 + A_z Z^2 + 6 L B_z Z + L C_z,  with")
Az = 324 * a**2 * c**2 - 216 * a * b * c + 408 * a - 15 * b**3 * c + 6 * b**2
Bz = 27 * a**2 * c**2 - 18 * a * b * c + 52 * a + b**3 * c + 14 * b**2
Cz = (729 * a**4 * c**4 - 972 * a**3 * b * c**3 + 2322 * a**3 * c**2
      + 54 * a**2 * b**3 * c**3 + 270 * a**2 * b**2 * c**2
      - 3735 * a**2 * b * c - 338 * a**2 - 36 * a * b**4 * c**2
      + 122 * a * b**3 * c + 1372 * a * b**2
      + b**6 * c**2 - 2 * b**5 * c - 80 * b**4)
assert lead == 8
assert sp.expand(cAz - Az) == 0
assert sp.expand(cBz6L - 6 * L * Bz) == 0
assert sp.expand(cCzL - L * Cz) == 0
print("    A_z =", Az)
print("    B_z =", Bz)
print("    C_z =", Cz)
print("  LEAD = 8, CONSTANT.  Matches THM-2546 eq.(6) exactly.       [PASS]")
print("  Tr(z) = -A_z/8 ; compare Tr(x) = 0 (depressed), Tr(y) = 3b/2.")
print("  On L=0 (A_z = 12L + 9(b^2 p - 2q)):  C = Z^2 (8Z + A_z),")
assert sp.expand(Az - 12 * L - 9 * (b**2 * p - 2 * q)) == 0
print("    the double escape shadow of THM-2546 eq.(14).             [PASS]")
print()

banner("[3] Route IV: standalone integrality certificate for the core")
t1 = time.time()
cert = sp.expand(core.subs({a: F1, b: F2, c: F3, Z: z}, simultaneous=True))
assert cert == 0
print(f"  C(F1,F2,F3; z) = 0 identically in Q[x,y,z]  ({time.time()-t1:.2f}s)"
      "  [PASS]")
print("  Since lead(C) = 8 is a nonzero CONSTANT, z satisfies the monic")
print("  cubic C/8 over Q[a,b,c]:  z IS INTEGRAL.  This certificate alone")
print("  upgrades the THM-2546 z-claim to PROVED.")
print()

banner("[4] Quasi-homogeneous structure (weights w(a,b,c)=(2,1,-1), w(z)=2)")
sc = {x: x / t, y: t * y, z: t**2 * z}
assert sp.simplify(F1.subs(sc) - t**2 * F1) == 0
assert sp.simplify(F2.subs(sc) - t * F2) == 0
assert sp.simplify(F3.subs(sc) - F3 / t) == 0
print("  F(x/t, ty, t^2 z) = (t^2 a, t b, c/t): exact equivariance   [PASS]")
assert sp.simplify(core.subs({a: t**2 * a, b: t * b, c: c / t, Z: t**2 * Z})
                   - t**6 * core) == 0
print("  C(t^2 a, tb, c/t; t^2 Z) = t^6 C(a,b,c;Z): C is quasi-")
print("  homogeneous of weight 6 (so are E_x [wt -1] and P_r [wt 3]).")
print()

banner("[5] disc_Z of the z-core, factored exactly")
disc_z = sp.factor(sp.discriminant(core, Z))
print("  disc_Z(C) =", disc_z)
dd = sp.factor(sp.cancel(disc_z / (-4 * L)))
print("  disc_Z(C) / (-4L) =", dd)
print("  All three integral-coordinate discriminants share the factor L:")
print("    disc_x = -4 (27ac^2-9bc+8)^2 L")
print("    disc_r = -2916 a^2 L = -4 (27a)^2 L")
print("    disc_z as printed above.")
assert sp.simplify(disc_z) != 0
print("  disc_Z(C) is not identically zero: generically the three fibre")
print("  z-values are DISTINCT (z separates the generic fibre).")
print()

banner("[6] Route III: degree-bounded exact interpolation (this is a proof)")
# --- exact coefficient-degree tables --------------------------------------
ex_coeffs = [L, sp.Integer(0), p, -2 * c]          # x^3..x^0
N_poly = sp.Poly(N, x)
assert N_poly.degree() == 4
Ncoef = N_poly.all_coeffs()                        # x^4..x^0 (5 entries)
G_coeffs = [q * Z,
            -b * Z - Ncoef[0],
            -2 * Z - Ncoef[1],
            -Ncoef[2],
            -Ncoef[3],
            -Ncoef[4]]                             # x^5..x^0


def degtab(polys, vs):
    return [max(sp.degree(cf, v) if cf != 0 else 0 for cf in polys) for v in vs]


ta = degtab(ex_coeffs, [a, b, c])
tg = degtab(G_coeffs, [a, b, c, Z])
print("  max coeff-degrees of E_x  in (a,b,c)  :", ta)
print("  max coeff-degrees of DZ-N in (a,b,c,Z):", tg)
DA = 5 * ta[0] + 3 * tg[0]
DB = 5 * ta[1] + 3 * tg[1]
DC = 5 * ta[2] + 3 * tg[2]
DZdeg = 3 * tg[3]
print(f"  Sylvester/Leibniz bounds: deg_a<= {DA}, deg_b<= {DB},"
      f" deg_c<= {DC}, deg_Z<= {DZdeg}")
assert (DA, DB, DC, DZdeg) == (16, 24, 13, 3)
for v, bd in ((a, DA), (b, DB), (c, DC), (Z, DZdeg)):
    assert sp.degree(P_sym, v) <= bd
print("  (actual degrees of the symbolic P respect the bounds; deg_c is")
print("   attained exactly: deg_c P = ", sp.degree(P_sym, c), ")", sep="")

# --- integer evaluators ----------------------------------------------------
VARS = (a, b, c)


def make_ev(expr):
    terms = [(tuple(int(e) for e in mono), int(cf))
             for mono, cf in sp.Poly(expr, *VARS).terms()]

    def ev(av, bv, cv):
        s = 0
        for (i, j, k), cf in terms:
            s += cf * av**i * bv**j * cv**k
        return s
    return ev


ev_L = make_ev(L)
ev_p = make_ev(p)
ev_q = make_ev(q)
ev_N = [make_ev(cf) for cf in Ncoef]               # N4..N0


def bareiss_det(M):
    n = len(M)
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
            M[i][k] = 0
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def P_det(av, bv, cv, Zv):
    """8x8 Sylvester determinant of (E_x, D*Z-N) at integers: value of P.

    Sign convention: with 5 rows of E_x-coefficients above 3 rows of
    G-coefficients, the determinant is (-1)^(3*5) = -1 times sympy's
    Res_x(E_x, G) (verified at 3 probes below); we return the negation so
    the interpolated polynomial is literally P = Res_x(E_x, G).  Either
    sign is a determinant formula, so the Leibniz degree bounds and the
    'polynomial identity in the entries' argument are unaffected.
    """
    f = [ev_L(av, bv, cv), 0, ev_p(av, bv, cv), -2 * cv]           # deg 3
    n4, n3, n2, n1, n0 = (e(av, bv, cv) for e in ev_N)
    g = [ev_q(av, bv, cv) * Zv, -bv * Zv - n4, -2 * Zv - n3,
         -n2, -n1, -n0]                                             # deg 5
    M = [[0] * 8 for _ in range(8)]
    for i in range(5):
        for j2, cf in enumerate(f):
            M[i][i + j2] = cf
    for i in range(3):
        for j2, cf in enumerate(g):
            M[5 + i][i + j2] = cf
    return -bareiss_det(M)


def cubic_from_values(v0, v1, v2, v3):
    """Exact coefficients (p3,p2,p1,p0) of the cubic with P(k)=vk, k=0..3."""
    p0 = Fraction(v0)
    p3 = Fraction(v3 - 3 * v2 + 3 * v1 - v0, 6)
    p2 = Fraction(v2 - 2 * v1 + v0 - 6 * p3, 2)
    p1 = Fraction(v1 - v0) - p2 - p3
    for f in (p3, p2, p1, p0):
        assert f.denominator == 1
    return (int(p3), int(p2), int(p1), int(p0))


# --- sanity: determinant convention matches sympy's resultant --------------
for (av, bv, cv, Zv) in [(1, 1, 1, 1), (2, -3, 1, 2), (-1, 4, -2, 0)]:
    lhs = P_det(av, bv, cv, Zv)
    rhs = int(P_sym.subs({a: av, b: bv, c: cv, Z: Zv}))
    assert lhs == rhs
print("  determinant convention agrees with sympy's Res_x at 3 probes [PASS]")

# --- the grid --------------------------------------------------------------
AV = list(range(-8, 9))       # 17 = DA+1 values
BV = list(range(-12, 13))     # 25 = DB+1 values
CV = list(range(-7, 7))       # 14 = DC+1 values
print(f"  grid: a in [-8,8] ({len(AV)} pts), b in [-12,12] ({len(BV)} pts),"
      f" c in [-7,6] ({len(CV)} pts): {len(AV)*len(BV)*len(CV)} nodes")
assert len(AV) == DA + 1 and len(BV) == DB + 1 and len(CV) == DC + 1

t1 = time.time()
node_vals = {}
for av in AV:
    for bv in BV:
        for cv in CV:
            vs = [P_det(av, bv, cv, k) for k in range(4)]
            node_vals[(av, bv, cv)] = cubic_from_values(*vs)
print(f"  {4*len(node_vals)} exact 8x8 integer determinants"
      f" in {time.time()-t1:.1f}s")


def interp1d(xs, ys):
    """Exact monomial coefficients of the interpolating polynomial."""
    n = len(xs)
    coef = [Fraction(v) for v in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    poly = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        newp = [Fraction(0)] * n
        for k2 in range(n - 1):
            if poly[k2]:
                newp[k2 + 1] += poly[k2]
                newp[k2] -= poly[k2] * xs[i]
        newp[0] += coef[i]
        poly = newp
    return poly


t1 = time.time()
interp = {}  # (kZ) -> dict {(da,db,dc): int}
for kZ in range(4):
    # stage 1: interpolate along c
    s1 = {}
    for av in AV:
        for bv in BV:
            s1[(av, bv)] = interp1d(CV, [node_vals[(av, bv, cv)][3 - kZ]
                                         for cv in CV])
    # stage 2: along b
    s2 = {}
    for av in AV:
        for kc in range(len(CV)):
            s2[(av, kc)] = interp1d(BV, [s1[(av, bv)][kc] for bv in BV])
    # stage 3: along a
    out = {}
    for kc in range(len(CV)):
        for kb in range(len(BV)):
            pa = interp1d(AV, [s2[(av, kc)][kb] for av in AV])
            for ka, cf in enumerate(pa):
                if cf:
                    assert cf.denominator == 1
                    out[(ka, kb, kc)] = int(cf)
    interp[kZ] = out
print(f"  tensor Lagrange interpolation (exact rationals)"
      f" in {time.time()-t1:.1f}s")

# --- interpolant == symbolic resultant (two independent routes agree) ------
sym_terms = {}
for mono, cf in sp.Poly(P_sym, a, b, c, Z).terms():
    da_, db_, dc_, dz_ = (int(e) for e in mono)
    sym_terms.setdefault(dz_, {})[(da_, db_, dc_)] = int(cf)
assert all(interp[kZ] == sym_terms.get(kZ, {}) for kZ in range(4))
print("  interpolated P == symbolic resultant P, coefficient-by-")
print("  coefficient (all four Z-coefficients)                       [PASS]")


def eval_interp(av, bv, cv, Zv):
    tot = 0
    for kZ in range(4):
        s = 0
        for (i, j, k2), cf in interp[kZ].items():
            s += cf * av**i * bv**j * cv**k2
        tot += s * Zv**kZ
    return tot


t1 = time.time()
for (av, bv, cv), (p3v, p2v, p1v, p0v) in node_vals.items():
    for Zv, want in ((0, p0v), (1, p3v + p2v + p1v + p0v)):
        assert eval_interp(av, bv, cv, Zv) == want
print(f"  interpolant reproduces ALL {len(node_vals)} grid nodes"
      f" (checked at Z=0,1) in {time.time()-t1:.1f}s               [PASS]")

rng = random.Random(20260727)
ev_core = [make_ev(sp.expand(core_poly.nth(k))) for k in range(4)]
prim_leads = {}
n_off = 0
while n_off < 25:
    av = rng.randint(-60, 60)
    bv = rng.randint(-60, 60)
    cv = rng.randint(-60, 60)
    if (av, bv, cv) in node_vals:
        continue
    n_off += 1
    vs = [P_det(av, bv, cv, k) for k in range(4)]
    cub = cubic_from_values(*vs)
    for Zv, v in zip(range(4), vs):
        assert eval_interp(av, bv, cv, Zv) == v
    # primitive per-node cubic and its lead (the earlier sampled statistic)
    if any(cub):
        g0 = 0
        for v in cub:
            g0 = sp.igcd(g0, abs(v))
        prim = [v // g0 for v in cub]
        if prim[0] < 0:
            prim = [-v for v in prim]
        # node cubic must be proportional to the specialized core
        cw = [e(av, bv, cv) for e in (ev_core[3], ev_core[2],
                                      ev_core[1], ev_core[0])]
        gc0 = 0
        for v in cw:
            gc0 = sp.igcd(gc0, abs(v))
        cwp = [v // gc0 for v in cw]
        if cwp[0] < 0:
            cwp = [-v for v in cwp]
        assert prim == cwp
        prim_leads[prim[0]] = prim_leads.get(prim[0], 0) + 1
print("  25 random off-grid integer targets in [-60,60]^3: interpolant")
print("  matches the determinant at Z=0..3, and each node's PRIMITIVE")
print("  cubic equals the primitivized specialized core              [PASS]")
print("  observed primitive per-node leads (counts):", dict(sorted(
    prim_leads.items())))
print("  -- exactly the {1,2,4,8} statistic of the old sampled run:")
print("     the symbolic lead is 8; per-node content can absorb 2-powers.")
print()

banner("[7] Exact fibre spot-checks (forward images, integer/rational)")
for (x0, y0, z0) in [(1, 1, 1), (2, -1, 3), (0, 2, -1), (1, 2, 3),
                     (-1, sp.Rational(3, 2), sp.Rational(13, 2))]:
    a0 = F1.subs({x: x0, y: y0, z: z0})
    b0 = F2.subs({x: x0, y: y0, z: z0})
    c0 = F3.subs({x: x0, y: y0, z: z0})
    val = core.subs({a: a0, b: b0, c: c0, Z: z0})
    assert sp.cancel(val) == 0
    exv = sp.cancel(Ex.subs({a: a0, b: b0, c: c0, x: x0}))
    assert exv == 0
    print(f"  (x,y,z)=({x0},{y0},{z0}) -> target ({a0},{b0},{c0}):"
          " C(...;z)=0 and E_x(x)=0  [PASS]")
print()

banner("[8] VERDICT")
print("  P(Z) = Res_x(E_x, D Z - N) = rho * chi(Z),")
print("     rho = 64 a^2 c^3 (27ac^2-9bc+8)^2,   cont(P) = rho/8,")
print("  PRIMITIVE z-ELIMINANT CORE (genuine cubic; roots = the three")
print("  fibre z-values, by Poisson + THM-2473 fibre law):")
print()
print("     C(Z) = 8 Z^3 + A_z Z^2 + 6 L B_z Z + L C_z")
print()
print("  with A_z, B_z, C_z as printed in [2].  LEADING COEFFICIENT = 8,")
print("  A CONSTANT 2-POWER: z is integral over Q[a,b,c].  PROVED by")
print("  (i) the degree-bounded interpolation (bounds 16/24/13, grid")
print("      17x25x14, all-node + 25 off-grid verification),")
print("  (ii) the independent symbolic resultant, and")
print("  (iii) the standalone certificate C(F;z) = 0 in Q[x,y,z].")
print("  This confirms and independently re-proves THM-2546 eq.(6)-(7):")
print("  the z-integrality claim is PROVED (no repair needed), and the")
print("  integral-coordinate dichotomy stands: x is the unique escape")
print("  coordinate; its eliminant lead is L; y-lead 2; z-lead 8 = |detJ|^3.")
print()
print(f"all exact checks passed   (total {time.time()-T0:.1f}s)")
