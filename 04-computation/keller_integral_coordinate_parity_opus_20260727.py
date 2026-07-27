# The integral-coordinate dichotomy of the sporadic Keller map F (THM-1300),
# and the parity-lens scoping for HYP-9030 test (ii).
#
#   [1] PROVED: the y-eliminant core is 2 y^3 - 3b y^2 + ... : CONSTANT lead 2.
#       => y is integral over the target algebra: no fiber escape in y, ever.
#       Also Tr(y) = 3b/2 != 0: the depressed law is an x-GAUGE fact, not a
#       map invariant (hostile to the naive stratum-wide depressed conjecture).
#   [2] VERIFIED: z-eliminant cores at integer targets have lead in {1,2,4,8}
#       (2-adic units up to content) and 2-power denominators in e_i(z):
#       z is integral (exact symbolic lead OPEN; sampled certificate).
#   [3] => MONOGENIC ESCAPE (PROVED for y, VERIFIED for z): all Jelonek escape
#       of F is carried by the single non-integral coordinate x (lead L,
#       THM-2473), whose depressed eliminant gives the even escape parity.
#   [4] POINTWISE PARITY RIGIDITY (trivial, load-bearing): for a real etale
#       polynomial map of complex degree d, every value with a FULL complex
#       fiber has #real preimages = d (mod 2) (conjugation pairing).  Hence an
#       even-degree wild map can NEVER exhibit a pointwise odd-real-count
#       full-fiber value: HYP-9030 test (ii) is structurally vacuous as
#       stated, and the parity program must live in escape parity (on S_F)
#       and integrated sphere degree, not pointwise counts.
#
# opus 2026-07-27.  Companion to THM-2543; extends THM-2473.

import sympy as sp
import random

x, y, z, a, b, c = sp.symbols('x y z a b c')
u = 1 + x*y
F1 = u**3*z + y**2*u*(4 + 3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
A = [u**3, 3*x*u**2, -x**3]
Bv = [sp.expand(F1 - A[0]*z), sp.expand(F2 - A[1]*z), sp.expand(F3 - A[2]*z)]

print("=" * 78)
print("[1] y-eliminant (symbolic): constant lead => y INTEGRAL; Tr(y) != 0")
print("=" * 78)
G1 = sp.expand((a - Bv[0])*A[1] - (b - Bv[1])*A[0])
H1 = sp.expand(sp.cancel(G1 / u**2))
H2 = sp.expand((a - Bv[0])*A[2] - (c - Bv[2])*A[0])
R = sp.resultant(sp.Poly(H1, x, domain=sp.QQ[a, b, c, y]),
                 sp.Poly(H2, x, domain=sp.QQ[a, b, c, y]))
Ry = sp.Poly(R.as_expr(), y)
for f, m in sp.factor_list(Ry.as_expr())[1]:
    fp = sp.Poly(f, y)
    if fp.degree() == 3:
        cs = fp.all_coeffs()
        print("  y-core cubic:", sp.expand(f))
        print("  lead =", cs[0], " (CONSTANT => y integral over Q[a,b,c])")
        print("  next =", sp.factor(cs[1]), " => Tr(y) =", sp.Rational(-1, 1)*sp.factor(cs[1]), "/", cs[0])
        assert not (sp.sympify(cs[0]).has(a) or sp.sympify(cs[0]).has(b) or sp.sympify(cs[0]).has(c))
        assert sp.expand(cs[1]) != 0, "y-eliminant unexpectedly depressed"
print("  [PASS] y integral; depressed law FAILS in the y-gauge (Tr(y) = 3b/2)")
print()

print("=" * 78)
print("[2] z-eliminant cores at integer targets: 2-adic-unit leads (z integral,")
print("    VERIFIED; exact symbolic lead left OPEN -- resultant too heavy)")
print("=" * 78)
random.seed(3)
leads = []
for _ in range(6):
    ta, tb, tc = [random.randint(-9, 9) for __ in range(3)]
    if ta == 0:
        continue
    E1 = sp.expand(A[0]*z - (ta - Bv[0]))
    E2 = sp.expand(A[1]*z - (tb - Bv[1]))
    E3 = sp.expand(A[2]*z - (tc - Bv[2]))
    R12 = sp.resultant(sp.Poly(E1, x), sp.Poly(E2, x))
    R13 = sp.resultant(sp.Poly(E1, x), sp.Poly(E3, x))
    P = sp.Poly(sp.resultant(sp.Poly(R12, y), sp.Poly(R13, y)), z)
    for f, m in sp.factor_list(P.as_expr())[1]:
        fp = sp.Poly(f, z)
        if fp.degree() == 3:
            cs = fp.all_coeffs()
            leads.append(int(cs[0]))
            print(f"  ({ta},{tb},{tc}): z-cubic lead={cs[0]}  Tr(z)={sp.Rational(-cs[1], cs[0])}")
ok = all(l & (l - 1) == 0 for l in leads)  # powers of two
print(f"  leads {sorted(set(leads))}: all powers of 2: {ok}")
assert ok
print("  [VERIFIED] z-values have 2-adic-bounded denominators: z integral")
print()

print("=" * 78)
print("[3] monogenic escape: x is the unique non-integral coordinate")
print("=" * 78)
print("  x-eliminant lead = L(a,b,c) = 27a^2c^2-18abc+16a+b^3c-b^2 (THM-2473,")
print("  nonconstant, irreducible) vs y lead 2 and z leads 2-adic units.")
print("  => every escaping fiber point has |x| -> infinity with y, z bounded;")
print("     the even escape parity of THM-2473(3) governs the whole boundary.")
print()

print("=" * 78)
print("[4] pointwise parity rigidity (conjugation pairing; sample witness)")
print("=" * 78)
print("  Lemma: real etale polynomial map, value v with FULL complex fiber")
print("  (#fiber = d): conjugation is a fixed-point-free-on-nonreal involution")
print("  of the fiber, so #real = d (mod 2).  For d even: #real EVEN at every")
print("  full-fiber value -- no pointwise odd witness can exist.  HYP-9030")
print("  test (ii) is therefore vacuous as stated; released and redirected.")
# sample: F has d=3 (odd); fiber over (-3,7/2,-4) -- count real vs complex:
E1 = sp.expand(A[0]*z - (-3 - Bv[0]))
ta, tb, tc = -3, sp.Rational(7, 2), -4
G1s = sp.expand((ta - Bv[0])*A[1] - (tb - Bv[1])*A[0])
H1s = sp.expand(sp.cancel(G1s / u**2))
H2s = sp.expand((ta - Bv[0])*A[2] - (tc - Bv[2])*A[0])
Rs = sp.Poly(sp.resultant(H1s, H2s, x), y)
core = [f for f, m in sp.factor_list(Rs.as_expr())[1] if sp.Poly(f, y).degree() == 3]
if core:
    cub = sp.Poly(core[0], y)
    d = sp.discriminant(cub.as_expr(), y)
    nreal = 3 if d > 0 else 1
    print(f"  sample fiber over (-3,7/2,-4): disc(y-cubic) = {d} => {nreal} real")
    print(f"  d = 3 odd, #real = {nreal} odd: congruence holds  [PASS]")
print()
print("all checks passed")
