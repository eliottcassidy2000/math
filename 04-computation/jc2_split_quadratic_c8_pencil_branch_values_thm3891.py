#!/usr/bin/env python3
"""Exact symbolic companion for THM-3891.

The branch-value lemma itself is proved by the all-degree Riemann--Hurwitz
argument in the theorem.  Here we freeze all algebraic identities, every
coordinate seam, representative gcd-degree controls, and a separately typed
finite hostile census.
"""

from __future__ import annotations

import hashlib
import itertools
import math

import sympy as sp


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def zero(label: str, expression: sp.Expr) -> None:
    check(sp.cancel(sp.expand(expression)) == 0, label)


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


def disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(b**2 * c**2 - 4 * a * c**3 - 4 * b**3 * d - 27 * a**2 * d**2 + 18 * a * b * c * d)


def coefficients(form: sp.Expr, U: sp.Symbol, V: sp.Symbol) -> tuple[sp.Expr, ...]:
    polynomial = sp.Poly(sp.expand(form), U, V)
    return tuple(
        polynomial.coeff_monomial(monomial)
        for monomial in (U**3, U**2 * V, U * V**2, V**3)
    )


def weighted_initial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, sp.Symbol],
    weights: tuple[int, int],
    degree: int,
    h: sp.Symbol,
) -> sp.Expr:
    scaled = sp.expand(
        expression.subs(
            {
                variables[0]: h ** weights[0] * variables[0],
                variables[1]: h ** weights[1] * variables[1],
            }
        )
    )
    polynomial = sp.Poly(scaled, h)
    check(
        all(exponent[0] >= degree for exponent, coefficient in polynomial.terms() if coefficient != 0),
        f"weighted order at least {degree}",
    )
    return sp.expand(polynomial.coeff_monomial(h**degree))


def support_inequality(
    expression: sp.Expr,
    first: sp.Symbol,
    second: sp.Symbol,
    first_weight: int,
    second_weight: int,
    bound: int,
    label: str,
) -> None:
    polynomial = sp.Poly(sp.expand(expression), first, second)
    check(
        all(
            first_weight * i + second_weight * j >= bound
            for (i, j), coefficient in polynomial.terms()
            if coefficient != 0
        ),
        label,
    )


# ---------------------------------------------------------------------------
# 1. Product discriminant and the two coordinate normal forms.
# ---------------------------------------------------------------------------

A, C, U, V = sp.symbols("A C U V")
L1, L2, M1, M2 = sp.symbols("L1 L2 M1 M2")
product_form = U * (L1 * U + M1 * V) * (L2 * U + M2 * V)
product_coefficients = coefficients(product_form, U, V)
equal(
    "product-of-determinants discriminant",
    disc(*product_coefficients),
    M1**2 * M2**2 * (L1 * M2 - L2 * M1) ** 2,
)

moving_form = U * (A * U + C * V) * ((A + C) * U + C * V)
moving_coefficients = coefficients(moving_form, U, V)
check(moving_coefficients == (A**2 + A * C, 2 * A * C + C**2, C**2, 0), "moving normal-form coefficients")
equal("moving normal-form discriminant", disc(*moving_coefficients), C**8)

f0a, f0b, f0c, f0d = sp.symbols("f0a f0b f0c f0d")
F0_discriminant = disc(f0a, f0b, f0c, f0d)
equal(
    "constant-carrier leading discriminant",
    disc(C**2 * f0a, C**2 * f0b, C**2 * f0c, C**2 * f0d),
    C**8 * F0_discriminant,
)


# ---------------------------------------------------------------------------
# 2. Self-contained replay of the moving-class Newton seams.
# ---------------------------------------------------------------------------

alpha, alpha1, beta, beta1, gamma, gamma1, delta, eta = sp.symbols(
    "alpha alpha1 beta beta1 gamma gamma1 delta eta"
)
aa = A**2 + A * C + alpha * A + alpha1 * C
bb = 2 * A * C + C**2 + beta * A + beta1 * C
cc = C**2 + gamma * A + gamma1 * C
dd = delta * A + eta * C
Delta_moving = disc(aa, bb, cc, dd)
x, z, w = sp.symbols("x z w")
H_moving = sp.expand(sp.cancel(z**8 * Delta_moving.subs({A: 1 / z, C: x / z})))
check(sp.denom(H_moving) == 1, "moving infinity chart polynomial")
equal("moving infinity restriction", H_moving.subs(z, 0), x**8)
support_inequality(H_moving, x, z, 1, 3, 6, "moving delta first support half-space")
support_inequality(H_moving, x, z, 1, 5, 8, "moving delta second support half-space")

def edge(expression: sp.Expr, order: int, degree: int) -> sp.Expr:
    substituted = sp.expand(expression.subs(z, w * x**order))
    return sp.expand(sp.Poly(substituted, x).coeff_monomial(x**degree))


equal("moving delta order-three edge", edge(H_moving, 3, 6), delta * w * (4 - 27 * delta * w))
equal("moving delta order-five edge", edge(H_moving, 5, 8), 1 + 4 * delta * w)

H_delta0 = sp.expand(H_moving.subs(delta, 0))
support_inequality(H_delta0, x, z, 1, 2, 6, "moving gamma first support half-space")
support_inequality(H_delta0, x, z, 1, 4, 8, "moving gamma unequal second support half-space")
equal(
    "moving gamma order-two edge",
    edge(H_delta0, 2, 6),
    w
    * (
        -4 * gamma**3 * w**2
        + (-27 * eta**2 + 36 * eta * gamma - 8 * gamma**2) * w
        + 4 * (eta - gamma)
    ),
)
equal("moving gamma order-four edge", edge(H_delta0, 4, 8), 1 + 4 * (eta - gamma) * w)
H_equal = sp.expand(H_delta0.subs(eta, gamma))
support_inequality(H_equal, x, z, 1, 2, 6, "moving equality first support half-space")
support_inequality(H_equal, x, z, 1, 3, 8, "moving equality second support half-space")
equal("moving equality order-two edge", edge(H_equal, 2, 6), gamma**2 * w**2 * (1 - 4 * gamma * w))
equal(
    "moving equality order-three edge",
    edge(H_equal, 3, 8),
    1 + 2 * (2 * beta + gamma - 2 * gamma1) * w + gamma**2 * w**2,
)
Delta_reducible = sp.expand(Delta_moving.subs({delta: 0, gamma: 0}))
check(sp.denom(sp.cancel(Delta_reducible / C)) == 1, "moving last seam C-divisible")


# ---------------------------------------------------------------------------
# 3. Both weighted initial forms in the constant-carrier class.
# ---------------------------------------------------------------------------

f1a, f1b, f1c, f1d = sp.symbols("f1a f1b f1c f1d")
f2a, f2b, f2c, f2d = sp.symbols("f2a f2b f2c f2d")
F0 = (f0a, f0b, f0c, f0d)
F1 = (f1a, f1b, f1c, f1d)
F2 = (f2a, f2b, f2c, f2d)
local_coefficients = tuple(x**2 * p + z * q + x * z * r for p, q, r in zip(F0, F1, F2))
H_constant = disc(*local_coefficients)
h = sp.symbols("h")
initial_12 = weighted_initial(H_constant, (x, z), (1, 2), 8, h)
expected_12 = disc(*(x**2 * p + z * q for p, q in zip(F0, F1)))
equal("constant-carrier weight-(1,2) initial", initial_12, expected_12)

lam, mu, s = sp.symbols("lam mu s", nonzero=True)
proportional_F1 = {q: lam * p for p, q in zip(F0, F1)}
local_proportional = tuple(sp.expand(value.subs(proportional_F1).subs(z, (s - x**2) / lam)) for value in local_coefficients)
expected_rows = tuple(s * p + x * s * r / lam - x**3 * r / lam for p, r in zip(F0, F2))
for index, (actual, expected) in enumerate(zip(local_proportional, expected_rows)):
    equal(f"proportional coordinate row {index}", actual, expected)
H_proportional = disc(*local_proportional)
initial_13 = weighted_initial(H_proportional, (x, s), (1, 3), 12, h)
expected_13 = disc(*(s * p - x**3 * r / lam for p, r in zip(F0, F2)))
equal("constant-carrier weight-(1,3) initial", initial_13, expected_13)

global_zero_F1 = tuple(C**2 * p + C * r for p, r in zip(F0, F2))
equal(
    "zero F1 reducible exit",
    disc(*global_zero_F1),
    C**4 * disc(*(C * p + r for p, r in zip(F0, F2))),
)
scalar_profile = C**2 + lam * A + mu * C
equal(
    "fully proportional fourth-power exit",
    disc(*(scalar_profile * p for p in F0)),
    scalar_profile**4 * F0_discriminant,
)


# ---------------------------------------------------------------------------
# 4. Exact controls for every gcd degree in the pencil lemma.
# ---------------------------------------------------------------------------

pencil_s, pencil_t = sp.symbols("pencil_s pencil_t")

# gcd degree two: two fixed roots give the two distinct collision values.
F_gcd2 = U * V * (U - V)
G_gcd2 = U * V * (U + V)
D_gcd2 = disc(*coefficients(pencil_s * F_gcd2 + pencil_t * G_gcd2, U, V))
equal("gcd-two collision divisor", D_gcd2, (pencil_t - pencil_s) ** 2 * (pencil_t + pencil_s) ** 2)

# gcd degree one: the residual quadratic has two branch values; the fixed
# root contributes the displayed additional collision factor.
F_gcd1 = U * (U - V) * (U + V)
G_gcd1 = U * (U**2 + U * V + 2 * V**2)
D_gcd1 = disc(*coefficients(pencil_s * F_gcd1 + pencil_t * G_gcd1, U, V))
equal(
    "gcd-one pencil divisor",
    D_gcd1,
    (-pencil_s + 2 * pencil_t) ** 2
    * (4 * pencil_s**2 - 4 * pencil_s * pencil_t - 7 * pencil_t**2),
)

# gcd degree zero: a squarefree part of projective degree at least two is a
# concrete hostile control for the Riemann--Hurwitz branch-value statement.
F_gcd0 = U * V * (U - V)
G_gcd0 = U**3 + V**3
D_gcd0 = sp.factor(disc(*coefficients(pencil_s * F_gcd0 + pencil_t * G_gcd0, U, V)))
D_gcd0_affine = sp.Poly(D_gcd0.subs(pencil_t, 1), pencil_s)
check(D_gcd0_affine.sqf_part().degree() >= 2, "gcd-zero at least two affine branch values control")

lambda_control = sp.symbols("lambda_control")
equal(
    "proportional one-support control",
    disc(*coefficients(pencil_s * F_gcd0 + pencil_t * lambda_control * F_gcd0, U, V)),
    (pencil_s + lambda_control * pencil_t) ** 4 * disc(*coefficients(F_gcd0, U, V)),
)


# ---------------------------------------------------------------------------
# 5. Declared 3^8 finite hostile census for F0=UV(U-V).
# ---------------------------------------------------------------------------

PolyDict = dict[tuple[int, int], int]

def pclean(poly: PolyDict) -> PolyDict:
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}

def padd(*polys: PolyDict) -> PolyDict:
    result: PolyDict = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            result[monomial] = result.get(monomial, 0) + coefficient
    return pclean(result)

def pscale(scale: int, poly: PolyDict) -> PolyDict:
    return pclean({monomial: scale * coefficient for monomial, coefficient in poly.items()})

def pmul(left: PolyDict, right: PolyDict) -> PolyDict:
    result: PolyDict = {}
    for (i, j), left_coefficient in left.items():
        for (r, t), right_coefficient in right.items():
            monomial = (i + r, j + t)
            result[monomial] = result.get(monomial, 0) + left_coefficient * right_coefficient
    return pclean(result)

def ppow(poly: PolyDict, exponent: int) -> PolyDict:
    result: PolyDict = {(0, 0): 1}
    for _ in range(exponent):
        result = pmul(result, poly)
    return result

def pdisc(pa: PolyDict, pb: PolyDict, pc: PolyDict, pd: PolyDict) -> PolyDict:
    return padd(
        pmul(ppow(pb, 2), ppow(pc, 2)),
        pscale(-4, pmul(pa, ppow(pc, 3))),
        pscale(-4, pmul(ppow(pb, 3), pd)),
        pscale(-27, pmul(ppow(pa, 2), ppow(pd, 2))),
        pscale(18, pmul(pmul(pmul(pa, pb), pc), pd)),
    )

def lower_hull(points: set[tuple[int, int]]) -> list[tuple[int, int]]:
    by_x: dict[int, int] = {}
    for i, j in points:
        by_x[i] = min(j, by_x.get(i, j))
    hull: list[tuple[int, int]] = []
    for point in sorted(by_x.items()):
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            cross = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1)
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    return hull

values = (-1, 0, 1)
pA: PolyDict = {(1, 0): 1}
pC: PolyDict = {(0, 1): 1}
pC2: PolyDict = {(0, 2): 1}
census_rows = 0
primitive_eligible = 0
for parameters in itertools.product(values, repeat=8):
    al, al1, be, be1, ga, ga1, de, et = parameters
    pa = padd(pscale(al, pA), pscale(al1, pC))
    pb = padd(pC2, pscale(be, pA), pscale(be1, pC))
    pc = padd(pscale(-1, pC2), pscale(ga, pA), pscale(ga1, pC))
    pd = padd(pscale(de, pA), pscale(et, pC))
    finite_delta = pdisc(pa, pb, pc, pd)
    census_rows += 1
    check(bool(finite_delta), "finite constant-carrier discriminant nonzero")
    degree = max(i + j for i, j in finite_delta)
    check(
        {(i, j): coefficient for (i, j), coefficient in finite_delta.items() if i + j == degree}
        == {(0, 8): 1},
        "finite constant-carrier pure C8 leading form",
    )
    local_points = {(j, degree - i - j) for i, j in finite_delta}
    hull = lower_hull(local_points)
    coordinate_free = min(j for i, j in finite_delta) == 0
    if len(hull) == 2:
        dx = hull[1][0] - hull[0][0]
        dy = hull[1][1] - hull[0][1]
        if coordinate_free and math.gcd(abs(dx), abs(dy)) == 1:
            primitive_eligible += 1

check(census_rows == 6561, "finite constant-carrier universe size")
check(primitive_eligible == 0, "finite constant-carrier primitive survivor count")


semantic_packet = "\n".join(
    (
        "split factor degree zero one one classification",
        "moving versus constant C squared carrier",
        "self contained moving Newton seams",
        "binary cubic pencil branch value lemma",
        "weight one two exceptional divisor",
        "weight one three proportional seam",
        "proportional reducible exits",
        "finite constant carrier hostile census",
    )
) + "\n"
semantic_sha256 = hashlib.sha256(semantic_packet.encode("utf-8")).hexdigest()

print("THM3891_CLASS moving_normal_form_or_constant_C2_carrier")
print("THM3891_MOVING delta_orders=3,5;gamma_orders=2,4_or_2,3;last_seam_reducible")
print("THM3891_PENCIL nonproportional_squarefree_cubic_pencil_has_at_least_two_branch_values")
print("THM3891_WEIGHTED first_weights=1,2;proportional_seam_weights=1,3")
print("THM3891_EXIT zero_first_profile_or_fully_proportional_profile_is_reducible")
print(f"THM3891_CENSUS rows={census_rows} primitive_irreducible_eligible={primitive_eligible}")
print("THM3891_SCOPE split_factor_degree_0_1_1_C8_only;nonsplit_quadratic_rows_open")
print(f"SEMANTIC_SHA256 {semantic_sha256}")
print(f"CHECKS {CHECKS}")
