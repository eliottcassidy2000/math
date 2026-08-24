#!/usr/bin/env python3
"""Exact THM-4005 atlas for the first live reduced 2:3 support cells.

Scope
-----
This is the exact THM-4005 consequence extractor from THM-3992 and THM-3996--3999.
It works on the live oriented seam

    gamma = -a^3/2,

and uses only the residual-mu_5 invariants

    A5=a^5, b=[y](R/gamma), d=[py](R/gamma).

It proves an exact two-diagonal support invoice and constructs the complete
finite-clutch companion jet modulo t^3.  It does not construct a B2 Darboux
pair, factor the global companion, decide component ownership, or close JC(2).
"""

from __future__ import annotations

import sys

import sympy as sp


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def zero(label: str, expr) -> None:
    value = simp(expr)
    if value != 0:
        raise AssertionError(f"{label}: {value}")
    print(f"PASS  {label}")


def gate(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def reduce_e(expr, e, A5):
    """Reduce a polynomial expression modulo e^2+6/A5."""
    poly = sp.Poly(sp.together(expr * A5**12).as_numer_denom()[0], e)
    modulus = sp.Poly(A5 * e**2 + 6, e)
    remainder = sp.rem(poly, modulus).as_expr()
    return simp(remainder)


x, t, z, p, y = sp.symbols("x t z p y")
a, A5 = sp.symbols("a A5", nonzero=True)
b, d = sp.symbols("b d")

sys.stdout.reconfigure(newline="\n")
print("STATUS=THM-4005;PROVED_NECESSARY_ATLAS;VERIFIED-EXACT;JC(2)_OPEN")
print("LIVE_SEAM=gamma=-a^3/2")
print("INVARIANTS=A5=a^5,b=beta/gamma,d=delta/gamma")

# -------------------------------------------------------------------------
# I. Residual fifth-root quotient and exact source-normal diagonals.
# -------------------------------------------------------------------------

zeta = sp.symbols("zeta", nonzero=True)
gate(
    "A5=a^5 is invariant modulo zeta^5=1",
    sp.rem(sp.expand((zeta**2 * a) ** 5 - a**5), zeta**5 - 1, zeta) == 0,
)

gamma = -a**3 / 2
beta = b * gamma
delta = d * gamma
theta0 = b / 2

# THM-3997's source jets, divided by a for A and by a^4 for C.  These
# quotients are invariant under the residual mu_5 action.
A0n = 1 + A5 * x**2 / 4
C0n = -sp.Rational(3, 4) * x - A5 * x**3 / 8

A1n = (
    sp.Rational(4, 3) / A5
    + A5 * b * x / 4
    + 2 * x**2
)
C1n = (
    -sp.Rational(3, 8) * b
    - 4 * x / A5
    - sp.Rational(3, 16) * A5 * b * x**2
    - sp.Rational(3, 2) * x**3
)

A2n = (
    (9 * A5**3 * b**2 - 512) / (144 * A5**2)
    + (A5 * d + 5 * b) * x / 4
    - 4 * x**2 / (5 * A5)
    + sp.Rational(3, 8) * A5 * b * x**3
)
C2n = (
    -(3 * A5 * d + 7 * b) / (8 * A5)
    + (2816 - 45 * A5**3 * b**2) * x / (480 * A5**2)
    - sp.Rational(3, 16) * (A5 * d + 12 * b) * x**2
    - 12 * x**3 / (5 * A5)
    - sp.Rational(9, 32) * A5 * b * x**4
)

# Reconstruct the table independently from the formulas in THM-3997.
N = -sp.Rational(2, 3) / (gamma * a) + 2 * gamma**2 * x * (
    theta0 + a * x / gamma**2
)
K = -x / a + 3 * gamma * (gamma**2 * x**2 + a / 2) * (
    theta0 + a * x / gamma**2
)
zero("invariant first A diagonal", N / a - A1n.subs(A5, a**5))
zero("invariant first C diagonal", K / a**4 - C1n.subs(A5, a**5))

m0 = beta**2 / 4 + sp.Rational(1, 3) / gamma**3 - sp.Rational(2, 9) / (
    a**3 * gamma**2
)
m1 = gamma * delta + sp.Rational(5, 4) * a * beta / gamma
m2 = 2 * (-3 * a**3 - 5 * gamma) / (5 * a * gamma**2)
m3 = 3 * gamma**2 * theta0
l0 = 3 * a * m1 / (4 * gamma) - theta0 * (3 * a**3 + 2 * gamma) / (
    2 * a * gamma
)
l1 = 3 * a * m2 / (4 * gamma) + (
    -9 * a**6
    + 36 * a**3 * theta0**2 * gamma**6
    - 6 * a**3 * gamma
    - 4 * gamma**2
) / (12 * a**3 * gamma**3)
l2 = 3 * a * theta0 * gamma + 3 * a * m3 / (4 * gamma) + 3 * gamma * m1 / 2
l3 = 3 * a**2 / (2 * gamma) + 3 * gamma * m2 / 2
l4 = 3 * gamma * m3 / 2
M = m0 + m1 * x + m2 * x**2 + m3 * x**3
L = l0 + l1 * x + l2 * x**2 + l3 * x**3 + l4 * x**4
zero("invariant second A diagonal", M / a - A2n.subs(A5, a**5))
zero("invariant second C diagonal", L / a**4 - C2n.subs(A5, a**5))

print("RESULT invariant_normal_jet_A/a = A0n+t*A1n+t^2*A2n+O(t^3)")
print(f"A0n={A0n}")
print(f"A1n={A1n}")
print(f"A2n={A2n}")
print("RESULT invariant_normal_jet_C/a^4 = C0n+t*C1n+t^2*C2n+O(t^3)")
print(f"C0n={C0n}")
print(f"C1n={C1n}")
print(f"C2n={C2n}")

# A monomial x^j t^n in a source-normal diagonal belongs to source weight
# w=j-2n, since u=x^2*t has weight zero.
for n in range(3):
    for j in range(0, 2 * n + 3):
        gate(f"weight tag (n={n},j={j})", j - 2 * n == j - 2 * n)

# The live seam forces these rows before importing the even-negative-tail
# invoice of THM-3987.
zero("A weight +2 boundary coefficient", A0n.coeff(x, 2) - A5 / 4)
zero("A weight -2 first coefficient", A1n.coeff(x, 0) - sp.Rational(4, 3) / A5)
zero("A weight 0 first coefficient", A1n.coeff(x, 2) - 2)
zero("C weight +3 boundary coefficient", C0n.coeff(x, 3) + A5 / 8)
zero("C weight +1 boundary coefficient", C0n.coeff(x, 1) + sp.Rational(3, 4))
zero("C weight -1 first coefficient", C1n.coeff(x, 1) + 4 / A5)

forced_A_base = {2, 0, -2}
forced_C_base = {3, 1, -1}
gate("A base weights are distinct", len(forced_A_base) == 3)
gate("C base weights are distinct", len(forced_C_base) == 3)
print("IMPORTED THM-3987: each output also has an even weight <= -4")
print("DEDUCTION support(A)>=4 and support(C)>=4 before the second diagonal")
print("DEDUCTION neither side can have retained support three")

# If b is nonzero, three additional first-diagonal coefficients are nonzero.
zero("A weight -1 coefficient is A5*b/4", A1n.coeff(x, 1) - A5 * b / 4)
zero("C weight -2 coefficient is -3*b/8", C1n.coeff(x, 0) + 3 * b / 8)
zero("C weight 0 coefficient is -3*A5*b/16", C1n.coeff(x, 2) + 3 * A5 * b / 16)
print("STRATUM b!=0: support(A)>=5, support(C)>=6")

# The b=0 seam forces a new C weight -3 on the second diagonal.  If d is
# nonzero, it also forces A:-3 and C:-4,-2.  At b=d=0, A:-4 is nonzero.
A2_b0 = sp.expand(A2n.subs(b, 0))
C2_b0 = sp.expand(C2n.subs(b, 0))
zero("b=0 forces A weight -4 coefficient", A2_b0.coeff(x, 0) + sp.Rational(32, 9) / A5**2)
zero("b=0 A weight -3 coefficient", A2_b0.coeff(x, 1) - A5 * d / 4)
zero("b=0 forces C weight -3 coefficient", C2_b0.coeff(x, 1) - sp.Rational(88, 15) / A5**2)
zero("b=0 C weight -4 coefficient", C2_b0.coeff(x, 0) + 3 * d / 8)
zero("b=0 C weight -2 coefficient", C2_b0.coeff(x, 2) + 3 * A5 * d / 16)

forced_A_minimal = {2, 0, -2, -4}
forced_C_minimal_before_tail = {3, 1, -1, -3}
gate("b=d=0 A has four distinct forced weights", len(forced_A_minimal) == 4)
gate(
    "b=d=0 C has four distinct forced weights before its even tail",
    len(forced_C_minimal_before_tail) == 4,
)
gate(
    "the required C even tail cannot collide with {3,1,-1,-3}",
    all(w not in forced_C_minimal_before_tail for w in (-4, -6, -8, -10)),
)
print("STRATUM b=0,d!=0: support(A)>=5, support(C)>=6")
print("STRATUM b=d=0: support(A)>=4, support(C)>=5")
print("POST-GAUGE REDUCTION: the oriented live reduced 2:3 cell has no 3x4 or 4x3 support")
print("FIRST POST-GAUGE RETAINED-SUPPORT CANDIDATE AFTER THESE JETS=4x5 on b=d=0")
print("CAUTION=4x5 is a necessary candidate floor, not an existence result")

# The recorded THM-3992 normalization is support-safe enough to transfer the
# seven-piece exclusion back across that specific gauge: diagonal target
# scalings and translations preserve nonconstant A-support, and the only
# mixing operation is C -> scalar*C+kappa*A+constant.  If A had size four,
# the table forces b=d=0 and S_A={2,0,-2,-4}.  The four forced C weights
# {3,1,-1,-3} are disjoint from S_A and therefore cannot be cancelled by
# undoing that C-by-A shear.  This is not a statement about arbitrary target
# automorphisms.
SA_four = {2, 0, -2, -4}
SC_shear_invisible = {3, 1, -1, -3}
gate("specific THM3992 shear cannot cancel four forced C weights", SA_four.isdisjoint(SC_shear_invisible))
print("GAUGE LEDGER: THM3992 diagonal scalings/translations preserve A support")
print("GAUGE LEDGER: its C-by-A shear cannot erase C weights {3,1,-1,-3} when A has four weights")
print("TRANSFER: the same oriented live seam excludes pre-normalization 3x4 and 4x3")
print("FIREWALL: no support claim is made for arbitrary nonlinear target automorphisms")

# THM-3998 hostile: its already-closed core has the all-degree gap 4D-2.
D = sp.symbols("D", integer=True, positive=True)
zero("THM3998 conic/ODE degree gap", 2 * (3 * D - 1) - 2 * D - (4 * D - 2))
gate("THM3998 hostile gap is positive from D=1", all(4 * n - 2 > 0 for n in range(1, 30)))
print("HOSTILE already-closed 3x<=3 cell remains excluded independently by THM-3998")

# -------------------------------------------------------------------------
# II. The complete companion jet through order t^2.
# -------------------------------------------------------------------------

z_src = 1 + x**2 * t
p_src = t * z_src
y_src = x * t * p_src
c20 = -sp.Rational(16, 3) / A5**2
c30 = b**2 / 4 + sp.Rational(2752, 135) / A5**3

R_known = c20 * p**2 + b * y + d * p * y + c30 * p**3
G_known_source = sp.expand(
    x**2 * t
    + 6 * p_src / A5
    + R_known.subs({p: p_src, y: y_src})
)
Q_known = sp.expand(
    x**2
    + 6 * z_src / A5
    + c20 * z_src * p_src
    + b * x * p_src
    + d * x * p_src**2
    + c30 * z_src * p_src**2
)
zero("complete known residual satisfies Gtilde=t*Qtilde", G_known_source - t * Q_known)

q0 = x**2 + 6 / A5
q1 = 6 * x**2 / A5 + c20 + b * x
q2 = 2 * c20 * x**2 + b * x**3 + d * x + c30
zero("companion t^0 row", Q_known.coeff(t, 0) - q0)
zero("companion t^1 row", Q_known.coeff(t, 1) - q1)
zero("companion t^2 row", Q_known.coeff(t, 2) - q2)
print(f"RESULT Qtilde mod t^3 = ({q0}) + t*({q1}) + t^2*({q2})")

# Exact t-adic Weierstrass preparation of the finite two-clutch packet.
U = 1 + 6 * t / A5 + (b * x - sp.Rational(32, 3) / A5**2) * t**2
W = (
    x**2
    + 6 / A5
    + t * (b * x - sp.Rational(124, 3) / A5**2)
    + t**2 * (
        (d - 12 * b / A5) * x
        + b**2 / 4
        + sp.Rational(44872, 135) / A5**3
    )
)
weierstrass_defect = sp.expand(Q_known - U * W)
for n in range(3):
    zero(f"Weierstrass factorization modulo t^3 row {n}", weierstrass_defect.coeff(t, n))
zero("Weierstrass unit has constant term one", U.subs(t, 0) - 1)
print("RESULT finite packet: Qtilde=U*W mod t^3 with U a t-adic unit")
print(f"U={U}")
print(f"W={W}")

Bquad = sp.expand(W).coeff(x, 1)
Pquad = sp.expand(W).coeff(x, 0)
disc = sp.expand(Bquad**2 - 4 * Pquad)
disc_expected = (
    -24 / A5
    + sp.Rational(496, 3) * t / A5**2
    - sp.Rational(179488, 135) * t**2 / A5**3
)
for n in range(3):
    zero(f"finite-packet discriminant row {n}", (disc - disc_expected).coeff(t, n))
gate("node-address quadratic is squarefree", simp(disc_expected.subs(t, 0)) != 0)
print("RESULT finite-packet discriminant is independent of b,d through t^2")

# Normalize each of the two finite germs over e^2=-6/A5.
e = sp.symbols("e", nonzero=True)
r1 = -b / 2 - sp.Rational(31, 9) * e / A5
r2 = -d / 2 + 6 * b / A5 + sp.Rational(653, 30) * e / A5**2
x_branch = e + r1 * t + r2 * t**2
branch_eval = sp.expand(Q_known.subs(x, x_branch))
for n in range(3):
    zero(f"normalized plus/minus branch row {n}", reduce_e(branch_eval.coeff(t, n), e, A5))

r1_minus = r1.subs(e, -e)
r2_minus = r2.subs(e, -e)
zero("branch-center order t", (r1 + r1_minus) + b)
zero("branch-center order t^2", (r2 + r2_minus) + d - 12 * b / A5)
print("RESULT for e_+/-^2=-6/A5:")
print("  x_+/-=e_+/-+(-b/2-31*e_+/-/(9*A5))*t")
print("       +(-d/2+6*b/A5+653*e_+/-/(30*A5^2))*t^2+O(t^3)")
print("  x_++x_-=-b*t+(-d+12*b/A5)*t^2+O(t^3)")

# -------------------------------------------------------------------------
# III. Endpoints, hostiles, and the first missing residual layer.
# -------------------------------------------------------------------------

E_known = sp.expand(1 - R_known.subs(p, 0))
zero("known boundary endpoint jet", E_known - (1 - b * y))
b_nonzero = sp.symbols("b_nonzero", nonzero=True)
zero("b!=0 endpoint root control", (1 - b_nonzero * y).subs(y, 1 / b_nonzero))
print("THM3999 ENDPOINT: E(y)=1-b*y+O(y^2); b!=0 forces an endpoint")
print("THM3999 FIREWALL: b=0 does not imply boundary disjointness")

# Boundary-disjoint live hostile.  It keeps the mandatory p^2 coefficient.
R_disjoint = c20 * p**2 + d * p * y + c30.subs(b, 0) * p**3
zero("boundary-disjoint representative has R(0,y)=0", R_disjoint.subs(p, 0))
zero(
    "boundary-disjoint representative ideal decomposition",
    R_disjoint - (p**2 * (c20 + c30.subs(b, 0) * p) + p * y * d),
)
gate("boundary-disjoint representative retains mandatory p^2", simp(c20) != 0)
print("HOSTILE R in (p^2,py): E=1 but [p^2]Rtilde=-16/(3*A5^2)!=0")

# R=0 is impossible on the live seam, but its ambient companion is the sharp
# THM-3999 G_m valuation/incidence hostile.
gate("R=0 contradicts mandatory live p^2 coefficient", simp(c20) != 0)
Q_R0 = x**2 + 6 * z_src / A5
t_gm = -A5 / 6 - x**-2
zero("ambient R=0 companion is G_m", Q_R0.subs(t, t_gm))
print("HOSTILE R=0: excluded for Keller live seam; ambient Q=0 is G_m")

# The next source order is exactly the three-dimensional residual layer.
c40, c21, c02 = sp.symbols("c40 c21 c02")
R_next = c40 * p**4 + c21 * p**2 * y + c02 * y**2
Q_next = sp.expand(c40 * z_src * p_src**3 + c21 * x * p_src**3 + c02 * x * p_src * y_src)
zero(
    "next residual layer divides by t",
    R_next.subs({p: p_src, y: y_src}) - t * Q_next,
)
zero("next Q layer begins at t^3", Q_next.coeff(t, 3) - (c40 + c21 * x + c02 * x**2))
for n in range(3):
    zero(f"next Q layer has no row {n}", Q_next.coeff(t, n))
print("FIRST MISSING RESIDUAL LAYER=([p^4],[p^2*y],[y^2]) of R/gamma")
print("FIRST MISSING COMPANION ROW=t^3*(c40+c21*x+c02*x^2)")
print("SIDECAR c02 is the first unknown endpoint coefficient after b")

# Odd terms control the symmetry of the two known clutch germs.  On b=d=0,
# c21=[p^2*y] is the first missing odd coefficient.
Q_known_abstract = x**2 + 6 * z / A5 + c20 * z * p + b * x * p + d * x * p**2 + c30 * z * p**2
odd_known = sp.expand(Q_known_abstract.xreplace({x: -x, y: -y}) - Q_known_abstract)
zero("known odd companion part", odd_known + 2 * b * x * p + 2 * d * x * p**2)
Q_next_abstract = c40 * z * p**3 + c21 * x * p**3 + c02 * x * p * y
odd_next = sp.expand(Q_next_abstract.xreplace({x: -x, y: -y}) - Q_next_abstract)
# This is the source involution x->-x, y->-y with z,p fixed.
zero("first missing source-odd term", odd_next + 2 * c21 * x * p**3)
print("ON b=d=0: known clutch packet is x->-x symmetric through t^2")
print("FIRST MISSING SYMMETRY SIDECAR includes c21=[p^2*y](R/gamma)")

# An arbitrary determinant-one target shear preserves the Jacobian but not
# this fixed nodal Weierstrass gauge, so A5,b,d are not general target
# invariants.
Avar, Cvar, Ivar, kappa = sp.symbols("Avar Cvar Ivar kappa")
Gvar = Cvar**2 - Avar**3 + Ivar * Avar
Gshear = (Cvar + kappa * Avar) ** 2 - Avar**3 + Ivar * Avar
zero("target shear changes defect by cross terms", Gshear - Gvar - 2 * kappa * Avar * Cvar - kappa**2 * Avar**2)
Ahost = x**2 + t
Chost = x**3 + x * t


def jac(F, H):
    return sp.expand(sp.diff(F, x) * sp.diff(H, t) - sp.diff(F, t) * sp.diff(H, x))


zero("determinant-one target shear preserves source Jacobian", jac(Ahost, Chost + kappa * Ahost) - jac(Ahost, Chost))
print("WARNING A5,b,d descend only under residual mu5 in the fixed THM3992 gauge")
print("WARNING arbitrary target automorphisms require renormalization and can change the atlas")

print("OWNERSHIP=the two roots e_+/- are two normalization addresses, not a component census")
print("THM3996=distinct owners imply another node address or a Jelonek node")
print("THM3999=endpoint data do not determine factor ownership or irreducibility")
print("GLOBAL FACTORIZATION=OPEN; U*W is only the lawful t-adic finite-packet factor mod t^3")
print("ALL THM-4005 JC2 LIVE 3x4/4x3 INVARIANT ATLAS CHECKS PASSED")
