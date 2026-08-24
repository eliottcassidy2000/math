#!/usr/bin/env python3
"""Exact companion for THM-3979 (two-color formal cusp lifting)."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


def trunc(expression: sp.Expr, variable: sp.Symbol, order: int) -> sp.Expr:
    return sp.series(expression, variable, 0, order).removeO().expand()


def reduce_l(expression: sp.Expr, L: sp.Symbol) -> sp.Expr:
    """Reduce a polynomial modulo L^2=1/6."""
    return sp.Poly(sp.expand(expression), L).rem(
        sp.Poly(L**2 - sp.Rational(1, 6), L)
    ).as_expr().expand()


x, t, y, p, z, Y, L = sp.symbols("x t y p z Y L")

# ---------------------------------------------------------------------------
# 1. The two components of V(x) and their formal coordinates.
# ---------------------------------------------------------------------------

# The height-two determinantal relations.
z_src = 1 + x**2 * t
p_src = sp.expand(z_src * t)
y_src = sp.expand(x * z_src * t**2)
zero(z_src * (z_src - 1) - x**2 * p_src, "Danielewski relation")
zero(p_src * (z_src - 1) - x * y_src, "modification relation")
zero(z_src * y_src - x * p_src**2, "saturated minor")
zero(z_src * (z_src - 1)**2 - x**3 * y_src, "D-chart equation")

# B/(x) is the product of the z=0 and z=1 colors.  The displayed maps
# check the defining quotient relations on each factor.
for zz, pp, yy, label in ((0, 0, y, "D"), (1, p, 0, "L1")):
    zero(zz * (zz - 1), f"{label}: idempotent color")
    zero(pp * (zz - 1), f"{label}: p relation")
    zero(zz * yy, f"{label}: y relation")

# Explicit Hensel idempotents.  q^2=W^2=1+4x^2p.
q, W = sp.symbols("q W", nonzero=True)
e_l = (1 + q / W) / 2
e_d = (1 - q / W) / 2
zero(e_l + e_d - 1, "Hensel idempotents sum to one")
zero(sp.together(4 * W**2 * (e_l**2 - e_l) - (q**2 - W**2)),
     "L1 idempotent modulo q^2-W^2")
zero(sp.together(4 * W**2 * (e_d**2 - e_d) - (q**2 - W**2)),
     "D idempotent modulo q^2-W^2")
zero(sp.together(4 * W**2 * e_l * e_d + (q**2 - W**2)),
     "orthogonality modulo q^2-W^2")
zero((2 * z_src - 1)**2 - (1 + 4 * x**2 * p_src),
     "source square-root relation")

# ---------------------------------------------------------------------------
# 2. D completion and the exact cusp Darboux normal form.
# ---------------------------------------------------------------------------

# Solve z(z-1)^2=x^3 y recursively.  The derivative in z is one at z=0,
# so every coefficient is uniquely determined.  We keep enough terms to
# verify the formal coordinate and the finite polynomial control below.
ORDER = 15
zeta = sp.Integer(0)
for degree in range(1, ORDER):
    aa = sp.Symbol(f"za{degree}")
    trial = zeta + aa * x**degree
    coefficient = sp.expand(trial * (trial - 1)**2 - x**3 * y).coeff(x, degree)
    solution = sp.solve(sp.Eq(coefficient, 0), aa)
    gate(len(solution) == 1, f"D implicit coefficient {degree} is unique")
    zeta = sp.expand(zeta + solution[0] * x**degree)

zero(trunc(zeta * (zeta - 1)**2 - x**3 * y, x, ORDER),
     "D implicit series")
gate(zeta.coeff(x, 3) == y, "D implicit leading coefficient")

# z_y=x^3/((z-1)(3z-1)); hence eta=dx wedge dt=x*w dx wedge dy.
w_d = trunc(1 / ((zeta - 1) * (3 * zeta - 1)), x, ORDER)
zero(trunc(sp.diff(zeta, y) - x**3 * w_d, x, ORDER),
     "D volume-unit derivative")
gate(w_d.subs(x, 0) == 1, "D volume unit")

# Center the cusp parameter exactly as in the corrected finite seed.
# With y=Y+x/2, integrate X dX=x*w dx at fixed Y.
w_shift = trunc(w_d.subs(y, Y + x / 2), x, ORDER)
H = trunc(2 * sp.integrate(x * w_shift, (x, 0, x)), x, ORDER)
gate(sp.expand(H).coeff(x, 2) == 1, "formal transverse square leading term")
phi = trunc(sp.sqrt(trunc(H / x**2, x, ORDER)), x, ORDER - 2)
Xser = trunc(x * phi, x, ORDER)
zero(trunc(Xser**2 - H, x, ORDER), "formal transverse square root")
zero(trunc(Xser * sp.diff(Xser, x) - x * w_shift, x, ORDER - 1),
     "formal transverse volume identity")
gate(Xser.coeff(x, 1) == 1, "formal transverse coordinate")

# The standard cusp pair has determinant X because 6L^2=1.
XX, YY = sp.symbols("XX YY")
Astd = YY**2 + 2 * L * XX
Cstd = YY**3 + 3 * L * XX * YY
det_std = sp.diff(Astd, XX) * sp.diff(Cstd, YY) - sp.diff(Astd, YY) * sp.diff(Cstd, XX)
zero(reduce_l(det_std - XX, L), "standard cusp determinant")

# Its initial D jets match the simultaneous corrected seed.
A_d = trunc(Astd.subs({XX: Xser, YY: Y}).subs(Y, y - x / 2), x, 5)
C_d = trunc(Cstd.subs({XX: Xser, YY: Y}).subs(Y, y - x / 2), x, 5)
A_d_expected = y**2 + (2 * L - y) * x + x**2 / 4
C_d_expected = (
    y**3 + (3 * L * y - sp.Rational(3, 2) * y**2) * x
    + (sp.Rational(3, 4) * y - sp.Rational(3, 2) * L) * x**2
    - x**3 / 8
)
zero(trunc(A_d - A_d_expected, x, 4), "corrected D A jets")
zero(trunc(C_d - C_d_expected, x, 4), "corrected D C jets")

# ---------------------------------------------------------------------------
# 3. L1 completion and its surjective jet operator.
# ---------------------------------------------------------------------------

# On z=1, z is a unit, t=p/z and (x,p) are formal coordinates.
# A=t+2Lx, C=-x is already an exact Darboux pair.
A_l = t + 2 * L * x
C_l = -x
j_l = sp.diff(A_l, x) * sp.diff(C_l, t) - sp.diff(A_l, t) * sp.diff(C_l, x)
zero(j_l - 1, "L1 exact Darboux pair")

# A correction x^r f(t) in C has leading Jacobian response -r f x^(r-1),
# so the L1 jet operator is onto with no cokernel.
f0, f1, f2 = sp.symbols("f0 f1 f2")
f = f0 + f1 * t + f2 * t**2
for r in range(2, 8):
    Ccorr = C_l + x**r * f
    response = sp.expand(
        sp.diff(A_l, x) * sp.diff(Ccorr, t)
        - sp.diff(A_l, t) * sp.diff(Ccorr, x) - 1
    )
    zero(response.coeff(x, r - 1) + r * f,
         f"L1 jet operator at order {r}")

# ---------------------------------------------------------------------------
# 4. D jet operator, its one-dimensional cokernel, and delayed payment.
# ---------------------------------------------------------------------------

# The new x^r jets (a,c) respond at order r-1 by
# r*y*(3*y*a-2*c).  Its image is (y), its cokernel is the constant row,
# and its kernel is c=(3/2)y*a.
a0, a1, a2, h0, h1, h2 = sp.symbols("a0 a1 a2 h0 h1 h2")
aa = a0 + a1 * y + a2 * y**2
hh = h0 + h1 * y + h2 * y**2
Afirst = y**2 + x * hh
Cfirst = y**3 + sp.Rational(3, 2) * x * y * hh
for r in range(2, 8):
    Ak = Afirst + x**r * aa
    Ck = Cfirst + sp.Rational(3, 2) * x**r * y * aa
    detk = sp.expand(
        sp.diff(Ak, x) * sp.diff(Ck, y)
        - sp.diff(Ak, y) * sp.diff(Ck, x)
    )
    det0 = sp.expand(
        sp.diff(Afirst, x) * sp.diff(Cfirst, y)
        - sp.diff(Afirst, y) * sp.diff(Cfirst, x)
    )
    response = sp.expand(detk - det0)
    zero(response.coeff(x, r - 1),
         f"D kernel leaves current row fixed at order {r}")
    next_constant = sp.expand(response.coeff(x, r).subs(y, 0))
    zero(next_constant - sp.Rational(3, 2) * (r + 1) * h0 * a0,
         f"D delayed constant payment at order {r}")

# h0=2L is nonzero, so the delayed scalar is invertible in characteristic 0.
for r in range(2, 8):
    gate(reduce_l(3 * (r + 1) * L, L) != 0,
         f"D delayed payment scalar nonzero at order {r}")

# ---------------------------------------------------------------------------
# 5. Exact finite polynomial control through simultaneous order five.
# ---------------------------------------------------------------------------

z0 = 1 + x**2 * t
p0 = sp.expand(z0 * t)
y0 = sp.expand(x * z0 * t**2)

A0 = y0**2 + 2 * L * x + p0 + x**2 * (1 - z0) / 4
C0 = (
    y0**3 + 3 * L * x * y0
    + sp.Rational(3, 2) * p0 * y0 * (1 - z0) - x * z0
    + x**2 * (-sp.Rational(3, 2) * L + sp.Rational(3, 8) * y0)
    + sp.Rational(3, 2) * L * x**2 * z0 * (1 - 2 * p0**2)
    + sp.Rational(3, 8) * x**2 * y0 * (1 - z0)
    - x**3 * (1 - z0) / 8
)

s2 = -p0 * (24 * p0**5 - 36 * p0**3 - 32 * p0**2 + 9 * p0 - 40) / 24
C0 += x**3 * z0 * s2 + x**4 * (1 - z0) * (y0 - sp.Rational(3, 2) * y0**3)

s3 = sp.Rational(2, 3) * L * (36 * p0**3 + 12 * p0**2 - 9 * p0 + 2)
C0 += x**4 * z0 * s3 / 4
A0 += sp.Rational(5, 6) * L * x**5 * (1 - z0)
C0 += x**5 * (1 - z0) * (
    -sp.Rational(2, 5) - sp.Rational(7, 4) * L * y0
    + sp.Rational(3, 4) * y0**2
)

s4 = (
    480 * p0**7 - 192 * p0**6 - 540 * p0**5 - 560 * p0**4
    + 135 * p0**3 - 464 * p0**2 + 17 * p0 - 6
) / 24
C0 += x**5 * z0 * s4 / 5 + x**6 * (1 - z0) * (y0 - 28 * L) / 16

J0 = reduce_l(
    sp.diff(A0, x) * sp.diff(C0, t) - sp.diff(A0, t) * sp.diff(C0, x) - 1,
    L,
)
poly_source = sp.Poly(J0, x)
source_order = min(monomial[0] for monomial, coefficient in poly_source.terms() if coefficient)
gate(source_order == 5, "finite control L1 order")
source_lead = sp.factor(poly_source.coeff_monomial(x**5))
source_expected = -sp.Rational(2, 15) * L * (
    324 * t**5 + 675 * t**4 + 330 * t**3 - 135 * t**2 + 51 * t - 4
)
zero(source_lead - source_expected, "finite control L1 leading row")

# Substitute the D-series z=zeta(x,Y), t=(z-1)/x^2.  Recompute enough
# implicit coefficients because Laurent substitution lowers x-order.
zeta_y = sp.Integer(0)
for degree in range(1, 28):
    aa_coeff = sp.Symbol(f"zb{degree}")
    trial = zeta_y + aa_coeff * x**degree
    coefficient = sp.expand(
        trial * (trial - 1)**2 - x**3 * Y
    ).coeff(x, degree)
    solution = sp.solve(sp.Eq(coefficient, 0), aa_coeff)
    gate(len(solution) == 1, f"finite D coefficient {degree} is unique")
    zeta_y = sp.expand(zeta_y + solution[0] * x**degree)

t_d = (zeta_y - 1) / x**2
J_d = reduce_l(trunc(J0.subs(t, t_d), x, 7), L)
d_orders = [degree for degree in range(-4, 7) if sp.expand(J_d).coeff(x, degree) != 0]
gate(d_orders[0] == 5, "finite control D order")
d_lead = sp.factor(sp.expand(J_d).coeff(x, 5))
d_expected = sp.Rational(7, 24) * (
    144 * L * Y**4 - 16 * L * Y**2 - 27 * L
    - 288 * Y**5 + 224 * Y**3
)
zero(d_lead - d_expected, "finite control D leading row")

summary = {
    "checks": CHECKS,
    "ring": "B2=k[x,1+x^2t,(1+x^2t)t,x(1+x^2t)t^2]",
    "completion": "Bhat_(x)=k[y][[x]] times k[p][[x]]",
    "D_volume": "eta=x*w dxdy; w=z_y/x^3 is a unit",
    "D_normal_form": "Y=y-x/2; X^2=2 integral x*w dx; A=Y^2+2LX; C=Y^3+3LXY",
    "L1_normal_form": "A=t+2Lx; C=-x",
    "jet_operators": "L1 surjective; D image (y), cokernel k, delayed kernel payment nonzero",
    "finite_control": "explicit polynomial lift has Jacobian error order exactly five on both colors",
    "scope": "formal/all-finite-order only; polynomial termination and Keller pair remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3979 two-color formal cusp lifting companion")
print(f"CHECKS={CHECKS}")
print("COMPLETION=K_Y_POWER_X_TIMES_K_P_POWER_X")
print("HENSEL_IDEMPOTENTS=EXPLICIT_SQRT_1_PLUS_4X2P")
print("D_VOLUME=X_TIMES_UNIT_DX_DY")
print("D_NORMAL_FORM=Y_YMINUSXOVER2;XDX=XWDX;CUSP_DETERMINANT_EXACT")
print("L1_NORMAL_FORM=ORDINARY_DARBOUX")
print("D_JET=IMAGE_Y;COKERNEL_K;DELAYED_KERNEL_PAYMENT_INVERTIBLE")
print("L1_JET=SURJECTIVE")
print(f"FINITE_SOURCE_ORDER={source_order};FINITE_SOURCE_LEAD={source_lead}")
print(f"FINITE_D_ORDER={d_orders[0]};FINITE_D_LEAD={d_lead}")
print("SCOPE=FORMAL_ALL_ORDER;POLYNOMIAL_TERMINATION_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
