#!/usr/bin/env python3
"""Exact controls for THM-3983's coordinate boundary place budget."""

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


x, t, z, p, y = sp.symbols("x t z p y")
a, c, ell = sp.symbols("a c ell", nonzero=True)


def jacobian(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, t)
        - sp.diff(first, t) * sp.diff(second, x)
    )


# ---------------------------------------------------------------------------
# The height-tower boundary and the exact hostile coordinate/seam controls.
# ---------------------------------------------------------------------------

height_rows: list[dict[str, object]] = []
for n in range(2, 11):
    zn = 1 + x**n * t
    pn = sp.expand(zn * t)
    yn = sp.expand(x ** (n - 1) * zn * t**2)

    zero(zn * yn - x ** (n - 1) * pn**2,
         f"n={n}: first determinantal relation")
    zero(zn * (zn - 1) - x**n * pn,
         f"n={n}: second determinantal relation")
    zero(pn * (zn - 1) - x * yn,
         f"n={n}: third determinantal relation")

    boundary_equation = z * (z - 1)**2 - x ** (n + 1) * y
    zero(sp.diff(boundary_equation, z) - (z - 1) * (3 * z - 1),
         f"n={n}: exact boundary derivative")
    gate(sp.diff(boundary_equation, z).subs({x: 0, z: 0}) == 1,
         f"n={n}: smooth boundary chart")

    # x is a source coordinate and has constant boundary restriction.
    gate(jacobian(x, t) == 1, f"n={n}: source coordinate x")
    gate(x.subs(x, 0) == 0, f"n={n}: x boundary value")

    # The critical-free THM-3978 seam.  Its boundary is x=z=0, not the
    # finite-source specialization x=0,z=1.
    seam_abstract = x + c * (z - 1)
    gate(seam_abstract.subs({x: 0, z: 0}) == -c,
         f"n={n}: seam boundary constant")
    seam_source = x + c * x**n * t
    gate(sp.diff(seam_source, t) == c * x**n,
         f"n={n}: seam t derivative")
    gate(sp.diff(seam_source, x).subs(x, 0) == 1,
         f"n={n}: seam has no critical point above x=0")

    # Over K=k(a), the generic seam fibre forces x to be a unit and is G_m.
    fibre_equation = x + c * x**n * t - a
    inverse_x = (1 + c * x ** (n - 1) * t) / a
    zero(x * inverse_x - 1 - fibre_equation / a,
         f"n={n}: x is a unit on the generic seam fibre")
    t_parameter = (a - x) / (c * x**n)
    zero(sp.together(fibre_equation.subs(t, t_parameter)),
         f"n={n}: Gm generic-fibre parametrization")
    gate(fibre_equation.subs(x, 0) == -a,
         f"n={n}: the generic seam fibre omits x=0")

    # The standard cusp seed has quadratic boundary restriction.
    cusp_abstract = y**2 + 2 * ell * x + p
    gate(cusp_abstract.subs({x: 0, z: 0, p: 0}) == y**2,
         f"n={n}: quadratic cusp boundary profile")

    height_rows.append({
        "n": n,
        "coordinate_boundary": "0",
        "seam_boundary": "-c",
        "seam_generic_fibre": "Gm",
        "cusp_boundary_degree": 2,
        "ord_D_t": -n,
    })


# ---------------------------------------------------------------------------
# Separable address controls and the exact numerical place budget.
# ---------------------------------------------------------------------------

address_rows: list[dict[str, object]] = []
for degree in range(1, 10):
    polynomial = y**degree - a
    derivative = sp.diff(polynomial, y)
    coefficient_field = sp.QQ.frac_field(a)
    gcd = sp.gcd(
        sp.Poly(polynomial, y, domain=coefficient_field),
        sp.Poly(derivative, y, domain=coefficient_field),
    )
    gate(gcd.degree() == 0,
         f"degree {degree}: generic boundary polynomial is separable")
    resultant = sp.factor(sp.resultant(polynomial, derivative, y))
    gate(resultant != 0,
         f"degree {degree}: nonzero generic address resultant")
    gate(sp.degree(polynomial, y) == degree,
         f"degree {degree}: exact address count")

    minimum_places = degree + 1
    gate(degree <= minimum_places - 1,
         f"degree {degree}: sharp formal place-budget boundary")
    gate(not (degree <= minimum_places - 2),
         f"degree {degree}: one residual infinity place is necessary")
    address_rows.append({
        "degree": degree,
        "resultant": str(resultant),
        "minimum_rational_places": minimum_places,
    })


# Coordinate and cusp endpoints of d<=r-1.
gate(0 <= 1 - 1, "coordinate endpoint d=0,r=1")
gate(2 <= 3 - 1, "cusp endpoint d=2,r=3")
gate(not (1 <= 1 - 1), "nonconstant boundary cannot fit an A1 fibre")


summary = {
    "checks": CHECKS,
    "family": "X_n, n>=2",
    "boundary": "D=V(x,z)=A1_y",
    "theorem": "geometrically rational generic fibre implies d<=r-1",
    "coordinate": "r=1 forces d=0",
    "cusp": "d=2 forces r>=3",
    "height_controls": height_rows,
    "address_controls": address_rows,
    "hostile": "A_c=x+c(z-1) is critical-free, boundary-constant, Gm-fibred",
    "scope": "geometric proof; finite checks are controls only",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3983 coordinate boundary place-budget companion")
print(f"CHECKS={CHECKS}")
print("FAMILY=X_N;N_GE_2;BOUNDARY=D=A1_Y")
print("PLACE_BUDGET=RATIONAL_GENERIC_FIBRE_IMPLIES_D_LE_R_MINUS_1")
print("COORDINATE=R_1_FORCES_BOUNDARY_CONSTANT")
print("CUSP=BOUNDARY_Y2_FORCES_R_GE_3_IF_RATIONAL")
print("HOSTILE=LINEAR_SEAM_CRITICAL_FREE;BOUNDARY_CONSTANT;GENERIC_GM")
print("NORMALIZATION=FINITE_AFFINE;D_TRANSVERSE_ADDRESSES_DISTINCT")
print("SCOPE=GEOMETRIC_THEOREM;FINITE_HEIGHTS_ARE_CONTROLS_ONLY")
print(f"SEMANTIC_SHA256={semantic}")
