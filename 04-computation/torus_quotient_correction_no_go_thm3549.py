#!/usr/bin/env python3
"""Exact planar quotient and low-degree correction audit for THM-3549.

All arithmetic is symbolic over QQ.  This file is the complete reproduction
artifact for the finite-exact claims in THM-3549.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def jac(left: sp.Expr, right: sp.Expr, a: sp.Symbol, b: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(left, a) * sp.diff(right, b)
        - sp.diff(left, b) * sp.diff(right, a)
    )


x, y, z = sp.symbols("x y z")
v, t, w = sp.symbols("v t w")
rho, sigma = sp.symbols("rho sigma")

unit = 1 + x * y
F1 = unit**3 * z + y**2 * unit * (4 + 3 * x * y)
F2 = y + 3 * x * unit**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x**2 * y - x**3 * z

# Source invariants are v=xy and t=x^2 z; w=2-3v-t is F3/x.
R = sp.expand(w * (4 * v + 6 - 3 * (v + 1) ** 2 * w))
S = sp.expand(w**2 * ((v + 1) * (v + 2) - (v + 1) ** 3 * w))
sub_invariants = {y: v / x, z: t / x**2}
R_from_F = sp.cancel((F2 * F3).subs(sub_invariants).subs(t, 2 - 3 * v - w))
S_from_F = sp.cancel((F1 * F3**2).subs(sub_invariants).subs(t, 2 - 3 * v - w))
require(sp.expand(R_from_F - R) == 0, "F2F3 quotient identity")
require(sp.expand(S_from_F - S) == 0, "F1F3^2 quotient identity")
J = sp.factor(jac(R, S, v, w))
require(J == -2 * w**2, "quotient Jacobian")

print("QUOTIENT")
print(f"R={sp.factor(R)}")
print(f"S={sp.factor(S)}")
print(f"Jac_vw(R,S)={J}")

# The two source-invariant collision points inherited from the triple collision.
p = {v: sp.Integer(0), w: sp.Integer(2)}
q = {v: -sp.Rational(3, 2), w: sp.Integer(0)}
require((R.subs(p), S.subs(p)) == (0, 0), "isolated collision point")
require((R.subs(q), S.subs(q)) == (0, 0), "contracted-line collision point")
print(f"collision_points={p},{q}; image=(0,0)")

# Exact generic fibre: resultant has an irrelevant w^6 factor and an essential cubic.
res_v = sp.factor(sp.resultant(R - rho, S - sigma, v))
essential = sp.expand(
    rho**3
    - rho**2
    - 3 * rho * w**2
    - 18 * rho * sigma
    - 2 * w**3
    + 4 * w**2
    + 27 * sigma**2
    + 16 * sigma
)
require(sp.expand(res_v + w**6 * essential) == 0, "generic fibre resultant")
require(sp.factor(essential) == essential, "generic cubic irreducible in QQ[rho,sigma,w]")
control_gb = sp.groebner([R - 1, S - 1], v, w, order="lex")
control_rows = tuple(sp.factor(row.as_expr()) for row in control_gb.polys)
expected_rows = (
    650 * v - 34 * w**2 - 133 * w + 675,
    (2 * w - 5) * (w**2 + 2 * w + 5),
)
require(control_rows == expected_rows, "three-point generic-fibre control")
require(sp.discriminant(expected_rows[1], w) == -67600, "control separability")
print(f"essential_fibre_cubic={sp.factor(essential)}")
print(f"control_(rho,sigma)=(1,1)_groebner={control_rows}")
print("generic_degree=3; control_discriminant=-67600")

# The exceptional zero fibre is exactly the contracted line plus one isolated point.
Rw2 = sp.factor(R.subs(w, 2))
Sw2 = sp.factor(S.subs(w, 2))
require(sp.gcd(sp.Poly(Rw2, v), sp.Poly(Sw2, v)).monic().as_expr() == v,
        "isolated point over zero")
sat_gb = sp.groebner([sp.cancel(R / w), sp.cancel(S / w**2)], v, w, order="lex")
sat_rows = tuple(sp.factor(row.as_expr()) for row in sat_gb.polys)
require(sat_rows == (v, w - 2), "zero-fibre residual saturation")
print("zero_fibre_set=V(w) union {(v,w)=(0,2)}")

# Keeping either quotient component is impossible at every correction degree:
# R has a critical point q, while grad(S) vanishes on the whole contracted line.
Rv, Rw = sp.diff(R, v), sp.diff(R, w)
Sv, Sw = sp.diff(S, v), sp.diff(S, w)
crit_R = sp.solve([Rv, Rw], [v, w], dict=True)
require(crit_R == [{v: -sp.Rational(3, 2), w: 0}], "R critical locus")
require(sp.expand(Sv.subs(w, 0)) == 0 and sp.expand(Sw.subs(w, 0)) == 0,
        "S critical line")
print("one_sided_correction_no_go=grad(R)(-3/2,0)=0; grad(S)|_(w=0)=0")


def correction_groebner(total_degree: int, preserve_collision: bool):
    monomials = [
        v**i * w**j
        for i in range(total_degree + 1)
        for j in range(total_degree + 1 - i)
        if i + j > 0
    ]
    aa = sp.symbols(f"a{total_degree}_0:{len(monomials)}")
    bb = sp.symbols(f"b{total_degree}_0:{len(monomials)}")
    kappa = sp.symbols(f"kappa_{total_degree}_{int(preserve_collision)}")
    A = sum(coefficient * monomial for coefficient, monomial in zip(aa, monomials))
    B = sum(coefficient * monomial for coefficient, monomial in zip(bb, monomials))
    equations = list(sp.Poly(jac(R + A, S + B, v, w) - kappa, v, w).coeffs())
    if preserve_collision:
        equations.extend([
            sp.expand(A.subs(p) - A.subs(q)),
            sp.expand(B.subs(p) - B.subs(q)),
        ])
    basis = sp.groebner(equations, *aa, *bb, kappa, order="grevlex")
    unit = any(row.as_expr() == 1 for row in basis.polys)
    kappa_remainder = sp.factor(basis.reduce(kappa)[1]) if not unit else sp.Integer(0)
    return len(monomials), len(equations), len(basis.polys), unit, kappa_remainder


print("LOW_TOTAL_DEGREE_CORRECTIONS")
for degree in (1, 2, 3):
    for collision in (False, True):
        row = correction_groebner(degree, collision)
        monomial_count, equation_count, basis_count, unit, kappa_remainder = row
        require(unit or kappa_remainder == 0, "constant Jacobian should be excluded")
        verdict = "EMPTY" if unit else "kappa_forced_0"
        print(
            f"D={degree};collision={collision};monomials_per_correction={monomial_count};"
            f"equations={equation_count};groebner_rows={basis_count};verdict={verdict}"
        )

# All-degree no-go for corrections affine in the transverse coordinate w.
# A=f(u)+a(u)w and B=g(u)+b(u)w, with u=v+1.
u = sp.symbols("u")
fu, au, gu, bu = (sp.Function(name)(u) for name in ("f", "a", "g", "b"))
P_aff = fu + w * (4 * u + 2 + au) - 3 * u**2 * w**2
Q_aff = gu + bu * w + u * (u + 1) * w**2 - u**3 * w**3
J_aff = sp.Poly(jac(P_aff, Q_aff, u, w), w)
c3 = sp.factor(J_aff.coeff_monomial(w**3))
c2 = sp.factor(J_aff.coeff_monomial(w**2))
expected_c3 = -3 * u**2 * (u * sp.diff(au, u) - au)
require(sp.expand(c3 - expected_c3) == 0, "affine-w highest coefficient")
c = sp.symbols("c")
c2_after = sp.factor(c2.subs({au: c * u, sp.diff(au, u): c}))
require(sp.expand(c2_after.subs(u, 0) + 2) == 0, "affine-w residual at u=0")
print("AFFINE_IN_w_CORRECTIONS")
print(f"coeff_w3={c3}")
print("coeff_w3=0 => a(u)=c*u for polynomial a in characteristic zero")
print(f"coeff_w2_after_a=cu={c2_after}")
print("coeff_w2_at_u=0=-2 => NO constant Jacobian (any degrees in u; collision not used)")

# Degree-floor identity for the first possible correction family.
A4 = sp.Function("A4")(u)
B4 = sp.Function("B4")(u)
top_equal_four = sp.factor(
    jac(A4 * w**4, B4 * w**4, u, w).coeff(w, 7)
)
require(
    sp.expand(top_equal_four - 4 * (sp.diff(A4, u) * B4 - A4 * sp.diff(B4, u))) == 0,
    "equal quartic top relation",
)
print(f"equal_w_degree_4_top_Jacobian_coefficient={top_equal_four}")
print("canonical_fibre_gate_consequence: collision-preserving Keller correction requires")
print("sorted final w-degrees >= (4,5); hence both outputs need new w^>=4 terms and one needs w^>=5")

# The first surviving (4,5) transverse top is itself a common-power skeleton.
a45 = sp.Function("a45")(u)
b45 = sp.Function("b45")(u)
top_four_five = sp.factor(jac(a45 * w**4, b45 * w**5, u, w).coeff(w, 8))
expected_four_five = 5 * sp.diff(a45, u) * b45 - 4 * a45 * sp.diff(b45, u)
require(
    sp.expand(top_four_five - expected_four_five) == 0,
    "quartic-quintic top relation",
)
print(f"w_degree_(4,5)_top_Jacobian_coefficient={top_four_five}")
print("coefficient=0 => a^5/b^4 is constant => a=c*h^4,b=d*h^5")
print("PASS")
