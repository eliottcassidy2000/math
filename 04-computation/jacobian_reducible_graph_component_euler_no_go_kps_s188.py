#!/usr/bin/env python3
"""Exact identities and hostile controls for THM-3568."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


a, b, x, y, z, H = sp.symbols("a b x y z H")

phi = -2 * H**3 * a + b * H**2 - 2 * H
D = 3 * a * H**2 - b * H + 1
C = 12 * a**2 * H**2 - 4 * a * b * H + 16 * a - b**2
P = b**2 * H**2 - 2 * b * H + 4
T = b * H - 2
L = sp.expand(27 * a**2 * phi**2 + 18 * a * b * phi + 16 * a - b**3 * phi - b**2)

require(sp.expand(L - D**2 * C) == 0, "Jelonek D-square times C factorization")
require(
    sp.expand(H**2 * C + T**2 - 4 * (a * H**2 + 1) * D) == 0,
    "D/C reduced intersection support",
)

# On the omitted curve, parameterized by target b with a=b^2/12 and
# c=4/(3b), the graph equation is a perfect cube of T.
omitted_restriction = sp.factor(sp.Rational(4, 3) / b + phi.subs(a, b**2 / 12))
require(sp.cancel(omitted_restriction + T**3 / (6 * b)) == 0, "omitted support equals T=0")

# The C-section is a quadratic in a.  Its finite projection has two points
# except at roots of H or P, where it has one.
require(sp.factor(sp.discriminant(C, a) - 64 * P) == 0, "C projection discriminant")
require(sp.factor(C.subs(H, 0) - (16 * a - b**2)) == 0, "C linear fibre over H=0")

# Exact inverse for the degree-one component over D!=0.
xt = H / D
yt = b - 3 * a * H
zt = a * D**3 - yt**2 * D * (D + 3)
ut = sp.factor(1 + xt * yt)
F1t = sp.factor(ut**3 * zt + yt**2 * ut * (4 + 3 * xt * yt))
F2t = sp.factor(yt + 3 * xt * ut**2 * zt + 3 * xt * yt**2 * (4 + 3 * xt * yt))
F3t = sp.factor(2 * xt - 3 * xt**2 * yt - xt**3 * zt)
require(sp.factor(ut - 1 / D) == 0, "degree-one inverse unit")
require(sp.factor(F1t - a) == 0, "degree-one inverse F1")
require(sp.factor(F2t - b) == 0, "degree-one inverse F2")
require(sp.factor(F3t + phi) == 0, "degree-one inverse F3")

# The source equation itself forces D*(1+xy)=1.  This independently checks
# that the degree-one component lands in, and fills, exactly D!=0.
u_src = 1 + x * y
F1 = sp.expand(u_src**3 * z + y**2 * u_src * (4 + 3 * x * y))
F2 = sp.expand(y + 3 * x * u_src**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
D_src = sp.expand(3 * F1 * H**2 - F2 * H + 1)
R_src = x - u_src * H
require(sp.rem(sp.Poly(D_src * u_src - 1, x), sp.Poly(R_src, x)).as_expr() == 0, "source D-unit identity")
require(
    sp.rem(sp.Poly(y - (F2 - 3 * F1 * H), x), sp.Poly(R_src, x)).as_expr() == 0,
    "source y inverse identity",
)


def distinct_root_count(poly: sp.Expr) -> int:
    polynomial = sp.Poly(sp.expand(poly), b, domain=sp.QQ)
    if polynomial.degree() <= 0:
        return 0
    return int(polynomial.sqf_part().degree())


# Reduced-support controls, including a double omitted root for h=b^2-3.
test_h = [sp.Integer(1), sp.Integer(2), b, b**2, (b - 1) ** 2, b**2 - 3]
rows = []
for h_value in test_h:
    r_h = distinct_root_count(h_value)
    p_value = sp.expand(P.subs(H, h_value))
    t_value = sp.expand(T.subs(H, h_value))
    m_h = distinct_root_count(p_value)
    n_h = distinct_root_count(t_value)
    chi_D = 1 - r_h
    chi_C = 2 - r_h - m_h
    chi_linear = 1 - chi_D
    chi_quadratic = 2 - chi_D - 2 * chi_C + n_h
    require(chi_quadratic == -3 + 3 * r_h + 2 * m_h + n_h, f"Euler formula h={h_value}")
    if sp.degree(h_value, b) == 0:
        require((r_h, m_h, n_h, chi_quadratic) == (0, 2, 1, 2), f"constant row h={h_value}")
    else:
        require(r_h >= 1 and m_h >= 2 and n_h >= 1, f"nonconstant lower bounds h={h_value}")
        require(chi_quadratic >= 5, f"nonconstant Euler no-go h={h_value}")
    rows.append((str(h_value), r_h, m_h, n_h, chi_D, chi_C, chi_linear, chi_quadratic))

# The repeated-root hostile: b*h-2=(b+1)^2*(b-2), so Euler uses support two,
# not algebraic degree three.
repeated_T = sp.factor(T.subs(H, b**2 - 3))
require(repeated_T == (b - 2) * (b + 1) ** 2, "reduced omitted-support hostile")

print("THM-3568 reducible target-graph component Euler audit")
print("L_phi=(3*a*h^2-b*h+1)^2*(12*a^2*h^2-4*a*b*h+16*a-b^2)")
print("D intersect C and omitted support: b*h-2=0 (reduced)")
print("degree-one component: A2 minus V(D), with unit D")
print("h r m n chi(D) chi(C) chi(linear) chi(quadratic)")
for row in rows:
    print(" | ".join(str(entry) for entry in row))
print("formula: chi(quadratic)=-3+3*r+2*m+n")
print("constant h: chi=2; nonconstant h: chi>=5")
print("verdict: neither irreducible component is A2 for any nonzero polynomial h")
print("all active truth gates passed")
