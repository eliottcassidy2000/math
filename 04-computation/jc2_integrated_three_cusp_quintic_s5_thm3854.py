#!/usr/bin/env python3
"""Exact companion for THM-3854's integrated three-cusp S5 quintic."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


def poly_pow_mod(
    base: sp.Poly, exponent: int, modulus: sp.Poly, prime: int
) -> sp.Poly:
    """Binary polynomial powering in F_prime[T]/(modulus)."""

    domain = sp.GF(prime)
    result = sp.Poly(1, T, domain=domain)
    power = sp.rem(base, modulus)
    while exponent:
        if exponent & 1:
            result = sp.rem(result * power, modulus)
        power = sp.rem(power * power, modulus)
        exponent //= 2
    return result


T, U, S, x, y, X, Y, Z = sp.symbols("T U S x y X Y Z")

# Integrate x'=4g and y'=15Tg for g=T(T^2-1).
g = T * (T**2 - 1)
x_T = T**4 - 2 * T**2
y_T = 3 * T**5 - 5 * T**3
zero(sp.diff(x_T, T) - 4 * g, "integrated x derivative")
zero(sp.diff(y_T, T) - 15 * T * g, "integrated y derivative")

delta = sp.expand(
    81 * x**5
    + 90 * x**4
    + 25 * x**3
    + 30 * x**2 * y**2
    + 30 * x * y**2
    - y**4
    + 8 * y**2
)
zero(delta.subs({x: x_T, y: y_T}), "normalization lies on quintic")
gate(sp.factor(delta) == delta, "quintic is irreducible")
gate(sp.Poly(delta, x, y).total_degree() == 5, "quintic degree")

# A rational inverse on a dense open set proves birationality.  Its
# denominator pulls back to precisely the finite exceptional address packet.
inverse_denominator = 27 * x**3 + 33 * x**2 + 10 * x + y**2
inverse_numerator = y * (9 * x**2 + 13 * x + 4)
exceptional_addresses = (
    T**2
    * (T - 1) ** 2
    * (T + 1) ** 2
    * (3 * T**2 - 5)
    * (9 * T**4 - 18 * T**2 + 4)
)
zero(
    inverse_denominator.subs({x: x_T, y: y_T}) - exceptional_addresses,
    "inverse denominator pullback",
)
zero(
    inverse_numerator.subs({x: x_T, y: y_T}) - T * exceptional_addresses,
    "inverse numerator pullback",
)

# The divided-difference resultant exhausts ramified and multiply-addressed
# finite points of the normalization.
x_U = U**4 - 2 * U**2
y_U = 3 * U**5 - 5 * U**3
x_difference = sp.cancel((x_T - x_U) / (T - U))
y_difference = sp.cancel((y_T - y_U) / (T - U))
collision_resultant = sp.factor(sp.resultant(x_difference, y_difference, U))
zero(collision_resultant - exceptional_addresses, "complete collision resultant")
zero(
    x_difference - (T + U) * (T**2 + U**2 - 2),
    "collision first factorization",
)
zero(
    y_difference.subs(U, -T) - T**2 * (3 * T**2 - 5),
    "opposite-address node branch",
)
second_collision_remainder = sp.rem(
    sp.Poly(y_difference, U), sp.Poly(U**2 - (2 - T**2), U)
).as_expr()
zero(
    second_collision_remainder - (3 * T**4 - 6 * T**2 + T * U + 2),
    "second node branch remainder",
)

# There are three A2 cusp addresses.  At T=+/-1 the displayed transverse
# coordinate cancels the common quadratic tangent.
v = sp.symbols("v")
for address, x_value, y_value, transverse_sign in (
    (0, 0, 0, None),
    (1, -1, -2, sp.Rational(-15, 4)),
    (-1, -1, 2, sp.Rational(15, 4)),
):
    local_x = sp.expand(x_T.subs(T, address + v) - x_value)
    local_y = sp.expand(y_T.subs(T, address + v) - y_value)
    gate(sp.Poly(local_x, v).terms()[-1][0] == (2,), f"cusp x order at {address}")
    if transverse_sign is None:
        transverse = local_y
    else:
        transverse = sp.expand(local_y + transverse_sign * local_x)
    gate(
        sp.Poly(transverse, v).terms()[-1][0] == (3,),
        f"cusp transverse order at {address}",
    )

# The three off-diagonal collision fibres are ordinary nodes.  The first is
# (-5/9,0); the other two have x=-4/9 and y^2=8/27.  At every non-cusp
# address the tangent slope is 15T/4, and paired addresses are distinct.
node_one_relation = 3 * T**2 - 5
zero(sp.rem(sp.Poly(x_T + sp.Rational(5, 9), T), sp.Poly(node_one_relation, T)),
     "first node x coordinate")
zero(sp.rem(sp.Poly(y_T, T), sp.Poly(node_one_relation, T)),
     "first node y coordinate")
node_pair_relation = 9 * T**4 - 18 * T**2 + 4
zero(sp.rem(sp.Poly(x_T + sp.Rational(4, 9), T), sp.Poly(node_pair_relation, T)),
     "paired nodes x coordinate")
zero(sp.rem(sp.Poly(y_T**2 - sp.Rational(8, 27), T), sp.Poly(node_pair_relation, T)),
     "paired nodes y coordinates")
node_partner = -sp.Rational(2, 3) / T
zero(
    sp.rem(
        sp.Poly(
            (x_U - x_T).subs(U, node_partner).together().as_numer_denom()[0],
            T,
        ),
        sp.Poly(node_pair_relation, T),
    ).as_expr(),
    "paired-node x equality",
)
zero(
    sp.rem(
        sp.Poly(
            (y_U - y_T).subs(U, node_partner).together().as_numer_denom()[0],
            T,
        ),
        sp.Poly(node_pair_relation, T),
    ).as_expr(),
    "paired-node y equality",
)
gate(sp.gcd(node_pair_relation, 3 * T**2 + 2) == 1,
     "paired node addresses are distinct")
zero(
    sp.cancel(sp.diff(y_T, T) / sp.diff(x_T, T)) - 15 * T / 4,
    "noncritical tangent slope",
)

# The projective curve has one smooth place at infinity, with fivefold line
# contact.  The finite 3A2+3A1 delta budget equals the quintic genus six.
delta_h = sp.expand(Z**5 * delta.subs({x: X / Z, y: Y / Z}))
zero(delta_h.subs(Z, 0) - 81 * X**5, "unique projective infinity point")
gate(
    sp.diff(delta_h, Z).subs({X: 0, Y: 1, Z: 0}) == -1,
    "projective infinity point is smooth",
)
projective_X = T**4 * S - 2 * T**2 * S**3
projective_Y = 3 * T**5 - 5 * T**3 * S**2
projective_Z = S**5
zero(
    delta_h.subs(
        {X: projective_X, Y: projective_Y, Z: projective_Z}, simultaneous=True
    ),
    "homogeneous normalization lies on projective quintic",
)
gate(3 + 3 == (5 - 1) * (5 - 2) // 2, "complete delta genus budget")

# The natural degree-five polynomial completion.  Its polynomial source is
# A2_(x,T), and its ramification curve maps birationally to delta.
F = 3 * T**5 - 10 * T**3 - 15 * x * T + 4 * y
F_T = sp.diff(F, T)
zero(F_T - 15 * (T**4 - 2 * T**2 - x), "quintic ramification equation")
zero(sp.discriminant(F, T) + 64800000 * delta, "exact quintic discriminant")
y_map = sp.expand((-3 * T**5 + 10 * T**3 + 15 * x * T) / 4)
zero(y_map.subs(x, x_T) - y_T, "ramification image parameter")
zero(sp.diff(y_map, T) + sp.Rational(15, 4) * (T**4 - 2 * T**2 - x),
     "polynomial map Jacobian")
gate(sp.Poly(F, y).degree() == 1 and sp.Poly(F, y).LC() == 4,
     "quintic total ring is a polynomial plane")

# The cuspidal normalization cannot instead be used as an interior smooth
# arm of a Keller map.  For arbitrary first-normal coefficients alpha,beta,
# the constant normal-bracket bucket remains divisible by g.
alpha, beta = sp.symbols("alpha beta")
arm_constant_bucket = alpha * sp.diff(y_T, T) - sp.diff(x_T, T) * beta
zero(
    arm_constant_bucket - g * (15 * T * alpha - 4 * beta),
    "all-degree arm constant-bucket factor",
)
gate(
    sp.gcd(sp.Poly(sp.diff(x_T, T), T), sp.Poly(sp.diff(y_T, T), T))
    == sp.Poly(g, T),
    "arm derivative ideal has nonunit gcd g",
)
gate(
    [arm_constant_bucket.subs(T, address) for address in (0, 1, -1)]
    == [0, 0, 0],
    "constant bracket vanishes at every cusp address",
)

# Arithmetic specialization proving the generic S5 group.  At (-3,1),
# reduction mod 29 is irreducible and reduction mod 67 has cycle pattern
# 1+1+1+2.  Both primes are unramified.
specialized = sp.Poly(F.subs({x: -3, y: 1}), T, domain=sp.ZZ)
specialized_discriminant = sp.discriminant(specialized.as_expr(), T)
gate(specialized_discriminant == 834688800000, "specialized discriminant")
for prime in (29, 67):
    gate(specialized_discriminant % prime != 0, f"unramified prime {prime}")

f29 = sp.Poly(specialized, T, modulus=29)
t29 = sp.Poly(T, T, modulus=29)
gate(
    sp.gcd(f29, poly_pow_mod(t29, 29, f29, 29) - t29).degree() == 0,
    "mod-29 degree-five irreducibility proper-divisor gate",
)
zero(
    (poly_pow_mod(t29, 29**5, f29, 29) - t29).as_expr(),
    "mod-29 degree-five irreducibility Frobenius gate",
)

factor67 = sp.Poly(
    3 * (T + 7) * (T + 21) * (T - 23) * (T**2 - 5 * T + 5),
    T,
    modulus=67,
)
gate(sp.Poly(specialized, T, modulus=67) == factor67,
     "mod-67 transposition factorization")
gate(pow(5, 33, 67) == 66, "mod-67 quadratic factor is irreducible")

# Two lowest seminormal functions do descend after squaring and cubing, but
# their depressed-cubic discriminants acquire nonsquare residual factors.
h_even = T**2 * (T**2 - 1) * (9 * T**2 - 14)
P_even = 81 * x**3 + 49 * x**2 + 8 * y**2
Q_even = (
    -243 * x**4
    - 143 * x**3
    + 81 * x**2 * y**2
    + 120 * x * y**2
    + 64 * y**2
)
zero(P_even.subs({x: x_T, y: y_T}) - h_even**2,
     "even seminormal square descent")
zero(Q_even.subs({x: x_T, y: y_T}) - h_even**3,
     "even seminormal cube descent")
even_residual = 6561 * x**4 + 3888 * x**3 - 512 * y**2
zero(P_even**3 - Q_even**2 - delta * even_residual,
     "even seminormal residual quotient")
gate(sp.factor(even_residual) == even_residual, "even residual nonsquare control")

h_odd = T * (T**2 - 1) * (9 * T**2 - 4) * (3 * T**2 - 5)
P_odd = -648 * x**3 - 720 * x**2 + 81 * x * y**2 - 200 * x + 49 * y**2
Q_odd = (
    6561 * x**4 * y
    + 17010 * x**3 * y
    + 18009 * x**2 * y
    + 243 * x * y**3
    + 8760 * x * y
    + 143 * y**3
    + 1600 * y
)
zero(P_odd.subs({x: x_T, y: y_T}) - h_odd**2,
     "odd seminormal square descent")
zero(Q_odd.subs({x: x_T, y: y_T}) - h_odd**3,
     "odd seminormal cube descent")
odd_residual_core = (
    41472 * x**2
    + 6561 * x * y**2
    + 46080 * x
    + 3888 * y**2
    + 12800
)
zero(
    P_odd**3
    - Q_odd**2
    + delta * (9 * x + 5) ** 2 * odd_residual_core,
    "odd seminormal residual quotient",
)
gate(sp.factor(odd_residual_core) == odd_residual_core,
     "odd residual core nonsquare control")

# The natural projective quadratic double-plane closure branches over
# Z*delta_h.  The two components have contact five at infinity, giving an
# A9 local packet.  This records why its 3-torsion is not inferred here.
gate(sp.Poly(delta_h.subs(Y, 1), X, Z).coeff_monomial(X**5) == 81,
     "line-quintic contact five leading term")
gate(sp.Poly(delta_h.subs(Y, 1), X, Z).coeff_monomial(Z) == -1,
     "line-quintic contact five transverse term")

semantic = {
    "front": "x=t4-2t2;y=3t5-5t3;one smooth infinity place",
    "singularities": "three A2 cusps plus three transverse A1 nodes",
    "completion": "F=3T5-10T3-15xT+4y;disc=-64800000*delta",
    "monodromy": "generic S5;over discriminant quadratic field A5;no C3 quotient",
    "arm_gate": "cuspidal normalization has derivative gcd g;no interior Keller arm",
    "seminormal_controls": "two square/cube descents have nonsquare residual quotients",
    "scope": "natural quintic splitting field only;Cl(W2=delta)[3] open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3854-integrated-three-cusp-quintic-s5-natural-completion-obstruction")
print("branch=integrated_three_cusp_quintic;finite_singularities=3A2_plus_3A1")
print("infinity_places=1;infinity_smooth=YES;line_contact=5")
print("natural_completion_degree=5;discriminant_constant=-64800000")
print("generic_group=S5;discriminant_kernel=A5;cyclic_cubic_quotient=NONE")
print("interior_smooth_arm=IMPOSSIBLE;derivative_gcd=t(t^2-1)")
print("seminormal_cubic_controls=2;nonsquare_residuals=2")
print("quadratic_resolvent_class_group_3_torsion=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
