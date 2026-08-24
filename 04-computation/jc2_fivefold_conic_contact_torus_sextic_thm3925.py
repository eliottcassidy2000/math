#!/usr/bin/env python3
"""Exact companion for THM-3925's fivefold-contact torus sextic.

Reproduction:
  python3 04-computation/jc2_fivefold_conic_contact_torus_sextic_thm3925.py
  python3 -O 04-computation/jc2_fivefold_conic_contact_torus_sextic_thm3925.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


X, Y, Z, S, T = sp.symbols("X Y Z S T")
x, y, t, u = sp.symbols("x y t u")
c = sp.symbols("c", nonzero=True)
ell, m, n = sp.symbols("ell m n")
P, Q = sp.symbols("P Q")


# ---------------------------------------------------------------------------
# The complete linear ambiguity after fixing contact divisor 5P+Q.
# ---------------------------------------------------------------------------

Q2 = Y * Z - X**2
Q30 = Z**2 * (Z - X)
linear = ell * X + m * Y + n * Z
Q3 = c * Q30 + Q2 * linear
Delta_h = sp.expand(4 * Q2**3 - 27 * Q3**2)

conic_param = {X: S * T, Y: T**2, Z: S**2}
zero(Q2.subs(conic_param), "marked conic parametrization")
zero(
    Q3.subs(conic_param) - c * S**5 * (S - T),
    "five-plus-one conic intersection divisor",
)
gate(
    [sp.diff(Q2, variable) for variable in (X, Y, Z)] == [-2 * X, Z, Y],
    "smooth conic gradient",
)

infinity = sp.factor(Delta_h.subs(Z, 0))
zero(
    infinity + X**4 * (4 * X**2 + 27 * (ell * X + m * Y) ** 2),
    "complete infinity restriction",
)
gate(
    sp.Poly(-infinity / X**4, X, Y).coeff_monomial(Y**2) == 27 * m**2,
    "unique infinity support forces m zero",
)
zero(
    infinity.subs(m, 0) + (4 + 27 * ell**2) * X**6,
    "single-support infinity seam",
)


# ---------------------------------------------------------------------------
# Affine cusp pullback, its Jacobian, and the rational normalization.
# ---------------------------------------------------------------------------

p = y - x**2
q = c * (1 - x) + p * (ell * x + n)
Delta = sp.expand(4 * p**3 - 27 * q**2)
zero(
    Delta_h.subs({X: x, Y: y, Z: 1, m: 0}) - Delta,
    "affine torus sextic",
)
gate(sp.Poly(Delta, x, y).total_degree() == 6, "sextic degree on the noncancelled seam")

jacobian = sp.factor(sp.diff(p, x) * sp.diff(q, y) - sp.diff(p, y) * sp.diff(q, x))
zero(jacobian - (c - ell * p), "coefficient-map Jacobian")
gate(
    sp.factor(
        sp.diff(p, x) * sp.diff(q, y) - sp.diff(p, y) * sp.diff(q, x)
    ).subs({p: 0, x: 1, y: 1})
    == c,
    "finite cusp is transverse in coefficient coordinates",
)

D = c - 3 * ell * t**2
N = c + 3 * n * t**2 - 2 * t**3
x_t = N / D
y_t = x_t**2 + 3 * t**2
zero(p.subs({x: x_t, y: y_t}) - 3 * t**2, "normalization square coordinate")
zero(q.subs({x: x_t, y: y_t}) - 2 * t**3, "normalization cube coordinate")
zero(Delta.subs({x: x_t, y: y_t}), "normalization lies on sextic")
zero(jacobian.subs({x: x_t, y: y_t}) - D, "normalization denominator is pulled-back Jacobian")

t_inverse = sp.cancel(3 * q / (2 * p))
zero(p - 3 * t_inverse**2 - Delta / (4 * p**2), "birational inverse square identity")
zero(q - 2 * t_inverse**3 - q * Delta / (4 * p**3), "birational inverse cube identity")

resultant = sp.factor(sp.resultant(D, N, t))
zero(
    resultant - c**2 * (4 * c - 27 * ell * (ell + n) ** 2),
    "normalization basepoint resultant",
)
gate(sp.discriminant(D, t) == 12 * c * ell, "two distinct finite Jacobian zeros")


# ---------------------------------------------------------------------------
# Projective normalization and exact infinity-address count.
# ---------------------------------------------------------------------------

Dh = c * S**2 - 3 * ell * T**2
Nh = c * S**3 + 3 * n * T**2 * S - 2 * T**3
Xh = sp.expand(Nh * S * Dh)
Yh = sp.expand(Nh**2 + 3 * T**2 * Dh**2)
Zh = sp.expand(S**2 * Dh**2)

gate(all(sp.Poly(form, S, T).total_degree() == 6 for form in (Xh, Yh, Zh)), "degree-six projective normalization")
zero(Delta_h.subs({X: Xh, Y: Yh, Z: Zh, m: 0}), "projective normalization identity")
zero(Xh.subs(S, 0), "infinite address has X zero")
zero(Zh.subs(S, 0), "infinite address has Z zero")
zero(Yh.subs(S, 0) - (4 + 27 * ell**2) * T**6, "infinite address survives exactly off degree-drop seam")
gate(sp.rem(Xh, Dh, S) == 0, "finite Jacobian-zero addresses have X zero")
gate(sp.rem(Zh, Dh, S) == 0, "finite Jacobian-zero addresses have Z zero")
gate(sp.rem(Yh - Nh**2, Dh, S) == 0, "finite Jacobian-zero addresses have Y=N squared")
gate(sp.factor(Zh) == S**2 * Dh**2, "no other infinity addresses")

Xh_fold = sp.expand(Xh.subs(ell, 0))
Yh_fold = sp.expand(Yh.subs(ell, 0))
Zh_fold = sp.expand(Zh.subs(ell, 0))
gate(sp.factor(Zh_fold) == c**2 * S**6, "one-place fold has only S zero at infinity")
zero(Xh_fold.subs(S, 0), "one-place fold X at infinity")
zero(Yh_fold.subs(S, 0) - 4 * T**6, "one-place fold Y at infinity")


# ---------------------------------------------------------------------------
# The one-place seam is a polynomial coordinate fold and is monogenic.
# ---------------------------------------------------------------------------

q_fold = sp.expand(q.subs(ell, 0))
jac_fold = sp.factor(
    sp.diff(p, x) * sp.diff(q_fold, y) - sp.diff(p, y) * sp.diff(q_fold, x)
)
zero(jac_fold - c, "one-place coefficient map has constant Jacobian")

x_inverse = 1 + (n * P - Q) / c
y_inverse = P + x_inverse**2
zero(p.subs({x: x_inverse, y: y_inverse}) - P, "polynomial inverse first coordinate")
zero(q_fold.subs({x: x_inverse, y: y_inverse}) - Q, "polynomial inverse second coordinate")

cubic = u**3 - P * u - Q
zero(sp.discriminant(cubic, u) - (4 * P**3 - 27 * Q**2), "universal cubic discriminant")
different = 3 * u**2 - P
companion = 4 * P - 3 * u**2
zero(
    (4 * P**3 - 27 * Q**2).subs(Q, u**3 - P * u) - different**2 * companion,
    "different-companion factorization",
)
gate(sp.Poly(different, P, u).total_degree() == 2, "different is nonconstant")


semantic_payload = {
    "contact": "smooth_conic_cubic_intersection_divisor_5P_plus_Q",
    "infinity": "unique_support_forces_m_zero_then_irreducible_nonfold_has_three_places",
    "jacobian": "normalization_denominator_equals_coefficient_map_jacobian_c_minus_ell_p",
    "one_place": "ell_zero_only_and_coefficient_map_is_explicit_polynomial_automorphism",
    "cubic": "universal_monogenic_A2_cover_with_nonconstant_different_unit_on_etale_open",
    "scope": "fixed_fivefold_contact_torus_sextic_family_closed_JC2_open",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3925-fivefold-conic-contact-torus-sextic-one-place-fold")
print("contact=Q2|conic=0;Q3|conic=c*S^5*(S-T)")
print("infinity=unique_support_requires_m=0")
print("irreducible_seam=resultant=c^2*(4c-27ell(ell+n)^2)")
print("nonfold_places=ell_nonzero_gives_two_finite_Jacobian_poles_plus_infinity")
print("one_place=ell_zero;coefficient_map_polynomial_automorphism_Jacobian_c")
print("completion=universal_monogenic_cubic_A2;different_unit_obstruction")
print("scope=fivefold_conic_contact_torus_sextic_family_closed;JC2_open")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
