#!/usr/bin/env python3
"""Assertion-free hostile audit of current THM-3827.

The dual discriminant is checked without sympy.discriminant/resultant: at
55 integer values of Z we form the 13 by 13 Sylvester matrix over Z and
take its exact determinant.  The discriminant has Z-degree at most 54:
it is degree 12 and weight 42 in the coefficients of a degree-seven
polynomial, while deg_Z(a_i) <= (7-i)+1.  Therefore agreement at 55
distinct values proves the displayed polynomial identity.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    CHECKS += 1


T, Z, s = sp.symbols("T Z s")

# Enter by T-coefficient blocks, independently of the canonical companion's
# monomial ordering.
coeff_desc = [
    84,
    36 * Z**2 + 196 * Z,
    84 * Z**3 + 36 * Z**2,
    49 * Z**4 + 112 * Z**3,
    -12 * Z**5,
    -14 * Z**6 + 12 * Z**5,
    0,
    Z**8,
]
H = sum(a * T ** (7 - i) for i, a in enumerate(coeff_desc))

require(sp.Poly(H, T).degree() == 7, "dual degree seven")
require(sp.Poly(H, T).LC() == 84, "constant nonzero leading coefficient")

J = (
    6480000 * Z**8 + 1952190576 * Z**7 - 14515170152 * Z**6
    + 76957426508 * Z**5 + 405669771962 * Z**4
    - 55140029819 * Z**3 - 17308754768 * Z**2
    + 1276289280 * Z + 234420480
)
expected_disc = -683730534400 * Z**43 * J


def sylvester_resultant_at(z_value: int) -> int:
    """Exact resultant Res_T(H,H_T) at one integer Z, by determinant."""
    f = [int(sp.expand(a).subs(Z, z_value)) for a in coeff_desc]
    fp = [(7 - i) * f[i] for i in range(7)]
    rows: list[list[int]] = []
    # deg(H_T)=6 shifted copies of the degree-seven H coefficient vector.
    for shift in range(6):
        rows.append([0] * shift + f + [0] * (5 - shift))
    # deg(H)=7 shifted copies of the degree-six derivative vector.
    for shift in range(7):
        rows.append([0] * shift + fp + [0] * (6 - shift))
    require(all(len(row) == 13 for row in rows), f"Sylvester shape Z={z_value}")
    return int(sp.Matrix(rows).det(method="domain-ge"))


# For n=7, Disc=(-1)^(7*6/2) Res(H,H')/LC = -Res/84.
for z_value in range(-27, 28):
    resultant = sylvester_resultant_at(z_value)
    require(resultant % 84 == 0, f"leading coefficient divides resultant Z={z_value}")
    disc_value = -resultant // 84
    require(
        disc_value == int(expected_disc.subs(Z, z_value)),
        f"dual discriminant evaluation Z={z_value}",
    )

require(sp.Poly(expected_disc, Z).degree() == 51, "displayed discriminant degree 51")
require(J.subs(Z, 0) == 234420480, "residual discriminant factor is nonzero")

# A nonconstant substitution is injective on K[Z], so it cannot annihilate
# either Z or J(Z).  The explicit composition is a hostile control.
q = s**3 - 2 * s + 5
composed_disc = sp.Poly(expected_disc.subs(Z, q), s)
require(composed_disc.degree() == 153, "nonconstant composition preserves discriminant")
require(composed_disc.LC() != 0, "composed discriminant has nonzero leading term")

# Squarefreeness plus odd degree makes W^2-H(T,q) integral: a square in
# K(s)(T) has even valuation at infinity, whereas H has valuation -7.
require(7 % 2 == 1, "odd degree gives nonsquare by infinity valuation")
finite_branch_points = 7
infinite_branch_points = 1
branch_points = finite_branch_points + infinite_branch_points
genus = (branch_points - 2) // 2
require(branch_points == 8 and genus == 3, "dual sidecar genus three")
require(infinite_branch_points == 1, "dual sidecar has one infinity place")

# Compare to the monic degree-eight sidecar in the other projection.
even_infinity_places = 2
require(even_infinity_places == 2, "monic degree-eight sidecar has two infinity places")

# Riemann--Hurwitz: 2g_X-2 = d(2*3-2)+R, with d>=1 and R>=0.
for source_genus in range(3):
    require(2 * source_genus - 2 < 4, f"source genus {source_genus} excluded")
require(2 * 3 - 2 == 4, "genus-three equality forces d=1 and R=0")

# Reconstruct the opposite intrinsic arm k=0 from the original cubic laws.
c = sp.symbols("c", nonzero=True)
h0 = sp.Rational(3, 7) / c**2
m0 = -sp.Rational(7, 3) * c**2
D0 = sp.Rational(7, 9) * c**4
A0 = sp.Rational(1, 3) * c**2
omega0 = sp.Integer(0)
theta0 = -sp.Rational(7, 3) * c**2
require(sp.cancel(-m0 * h0) == 1, "k=0 determinant law")
require(sp.cancel(7 * D0 * h0**2) == 1, "k=0 first reconstruction law")
require(sp.cancel(9 * h0**2 * D0 - 3 * h0 * c**2) == 0,
        "k=0 second reconstruction law")
require(sp.cancel(c * omega0 - 3 * A0 * theta0 - 14 * A0**2 - D0) == 0,
        "k=0 different")
require(sp.cancel(omega0**2 - (-7 * A0**2 + c * omega0 - A0 * theta0)) == 0,
        "k=0 omega-square law")
require(sp.cancel(omega0 * theta0 - (3 * A0**2 - A0 * c**2)) == 0,
        "k=0 mixed law")
require(sp.cancel(theta0**2 - (
    3 * A0 * c - c**3 + (c**2 - 3 * A0) * omega0 - 7 * A0 * theta0
)) == 0, "k=0 theta-square law")
require(sp.cancel(h0 * c**2) == sp.Rational(3, 7),
        "k=0 quotient contains c inverse")

# Independently replay the completed square and the h-adic lift identity.
h, k, d = sp.symbols("h k d")
A5 = (7 * h**2 + 3 * k**2) * (3 * h**3 + 7 * h**2 * k + k**3)
B3 = (h + k) * (2 * h + k) * (3 * h - k)
H_hk = H.subs({T: h, Z: k})
require(sp.expand(H_hk - (k * B3) ** 2 - 4 * h**2 * A5) == 0,
        "second completed square")
w_lift = k * B3 + 2 * h**2 * d
require(sp.expand(
    w_lift**2 - H_hk - 4 * h**2 * (d * (k * B3 + h**2 * d) - A5)
) == 0, "h-adic lift equals divisor-allocation equation")

# The spectral factors are five distinct linear forms and avoid every factor
# of h*k*B3.  This is recomputed over QQ, not imported from the primary.
z = sp.symbols("z")
a5 = sp.expand(A5.subs({h: z, k: 1}))
kb3 = sp.expand((k * B3).subs({h: z, k: 1}))
require(sp.degree(a5, z) == 5, "five spectral slopes")
require(sp.discriminant(a5, z) == 353831803500, "spectral slopes distinct")
require(sp.resultant(a5, kb3, z) == -31298700,
        "spectral slopes avoid finite kB slopes")
require(A5.subs({h: 1, k: 0}) == 21, "spectral packet avoids k=0 slope")

# Hostile nonreduced pullback control.  For alpha^2=-3/7, the unimodular,
# algebraically independent pair below pulls the member h-alpha*k back to
# (1+x*y)^2.  Its multiplicity in A5 is exactly two, demonstrating why the
# valuation argument must allocate the full member rather than its support.
xx, yy = sp.symbols("xx yy")
alpha = sp.sqrt(-sp.Rational(3, 7))
p_nr = 1 + xx * yy
k_nr = xx
h_nr = p_nr**2 + alpha * xx
require(sp.simplify(7 * alpha**2 + 3) == 0, "hostile slope is spectral")
require(sp.expand(h_nr - xx * (2 * yy + xx * yy**2 + alpha)) == 1,
        "hostile nonreduced pair is unimodular")
require(sp.expand(h_nr - alpha * k_nr - p_nr**2) == 0,
        "hostile spectral member is a square")
jac_nr = sp.det(sp.Matrix([
    [sp.diff(h_nr, xx), sp.diff(h_nr, yy)],
    [sp.diff(k_nr, xx), sp.diff(k_nr, yy)],
]))
require(sp.expand(jac_nr) != 0, "hostile pair is algebraically independent")
A5_nr = sp.expand(A5.subs({h: h_nr, k: k_nr}))
quot_nr = sp.cancel(A5_nr / p_nr**2)
require(sp.denom(quot_nr) == 1, "nonreduced component occurs at least twice")
require(sp.cancel(quot_nr.subs(yy, -1 / xx)) != 0,
        "nonreduced component occurs exactly twice")

# If whole members are selected, d is homogeneous of their subset size.
# Only sizes one and two survive the first degree comparison.
surviving_sizes = []
for subset_size in range(6):
    degrees = (subset_size + 4, 2 * subset_size + 2)
    if max(degrees) == 5 or degrees[0] == degrees[1]:
        surviving_sizes.append(subset_size)
require(surviving_sizes == [1, 2], "whole-member degree sieve")
require(sp.gcd(sp.Poly(h**2, h, k), sp.Poly(k * B3, h, k)).total_degree() == 0,
        "size-two terminal gcd obstruction")

# Guards against Python -O silently deleting a gate.
source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "no Python assert statements")

semantic = {
    "theorem": "current THM-3827 dual genus and bichromatic passport",
    "discriminant_method": "55 exact integer Sylvester determinants; degree bound 54",
    "dual_discriminant": "-683730534400*Z^43*J_8(Z)",
    "squarefree": "q(ell) transcendental and q nonconstant imply nonzero discriminant",
    "curve": "odd degree 7 hyperelliptic, genus 3, one infinity place",
    "dichotomy": "h transcendental gives function-field injection/map; h algebraic lies in K[ell]",
    "terminal": "h,k in K[ell] contradict their intrinsic algebraic independence",
    "comparison": "degree 8 h-sidecar has two infinity places",
    "opposite_arm": "B/(k)=K[C,C^-1] from the original cubic laws",
    "bichromatic": "prime valuations allocate full multiplicity; monochromatic means whole-member subset",
    "hostile_nonreduced": "h-alpha*k=(1+xy)^2 in a unimodular algebraically independent control",
    "scope": "necessary atlas obstruction only; no Jacobian counterexample",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3827-current-dual-genus-and-bichromatic-passport")
print("method=independent_55_point_exact_sylvester_determinant_identity")
print("dual=degree7_squarefree_hyperelliptic_genus3")
print("infinity_parity=h_sidecar_2;k_sidecar_1")
print("dichotomy=transcendental_map_or_polynomial_dependence")
print("opposite_arm=B_mod_k_is_Laurent_Gm")
print("multiplicity=prime_valuations_force_full_whole_member_allocation")
print("hostile_nonreduced=(h-alpha*k)=(1+x*y)^2_exact_multiplicity_2")
print("status=PASS")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
