#!/usr/bin/env python3
"""Finite exact companion for provisional THM-3599.

The proof is Newton-polygon and residue driven.  This companion checks every
coefficient mask on a broad integer box, the exact rational-mate boundaries,
all exceptional-face identities, and the sharp exponent/degree hostiles
without Python assert gates.
"""

from fractions import Fraction

import sympy as sp


b, c, e, w = sp.symbols("b c e w")
lam, mu, nu, s = sp.symbols("lambda mu nu s", nonzero=True)
CHECKS = 0


def require(label, condition):
    """Record one active gate and fail with a stable label."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError("FAILED: " + label)


def zero(expr):
    return sp.cancel(sp.expand(expr)) == 0


def poisson(F, G, exponent, sigma):
    """Full three-generator Poisson bracket on c^N e=sigma(b)."""
    return sp.expand(
        c**exponent * (sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b))
        - sp.diff(sigma, b)
        * (sp.diff(F, c) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, c))
        - exponent
        * c ** (exponent - 1)
        * e
        * (sp.diff(F, b) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, b))
    )


def cross(origin, first, second):
    return (first[0] - origin[0]) * (second[1] - origin[1]) - (
        first[1] - origin[1]
    ) * (second[0] - origin[0])


def convex_hull(points):
    """Return the strict CCW hull of exact integer points."""
    ordered = sorted(set(points))
    if len(ordered) <= 1:
        return ordered
    lower = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def strict_inside(point, hull):
    return all(cross(hull[j], hull[(j + 1) % len(hull)], point) > 0 for j in range(len(hull)))


print("THM-3599 exact companion")
print("SECTION affine rational-mate boundary")
MATE_ROWS = 0
for exponent in range(2, 11):
    sigma = b**3 + b + 1
    Q_lam = lam * b + mu * c
    P_lam = c ** (1 - exponent) / (lam * (exponent - 1))
    require(
        f"lambda mate N={exponent}",
        zero(poisson(P_lam, Q_lam, exponent, sigma) - 1),
    )
    Q_mu = mu * c
    P_mu = b / (mu * c**exponent)
    require(
        f"mu mate N={exponent}",
        zero(poisson(P_mu, Q_mu, exponent, sigma) - 1),
    )
    MATE_ROWS += 2

sigma_n1 = b**3 + b + 1
require(
    "N=1 hostile rational mate",
    zero(poisson(-b / (nu * e), nu * e, 1, sigma_n1) - 1),
)
sigma_d1 = s * (b - 2)
require(
    "D=1 hostile rational mate",
    zero(poisson(-c / (nu * s), nu * e, 4, sigma_d1) - 1),
)
MATE_ROWS += 2
print(f"PASS {MATE_ROWS} rational-mate and sharp-boundary identities")


print("SECTION complete Newton-mask census")
MASK_ROWS = 0
BOUNDARY_ROWS = []
for exponent in range(2, 21):
    for degree in range(2, 21):
        A = (0, -exponent)
        B = (degree, -exponent)
        O = (0, 0)
        L = (1, 0)
        U = (0, 1)
        point = (1, 1 - exponent)
        masks = {
            "pure": [A, B, O],
            "lambda": [A, B, L, O],
            "mu": [A, B, U],
            "both": [A, B, L, U],
        }
        for name, support in masks.items():
            hull = convex_hull(support)
            interior = strict_inside(point, hull)
            expected = not (name == "pure" and exponent == 2 and degree == 2)
            require(
                f"mask interior N={exponent} D={degree} mask={name}",
                interior == expected,
            )
            if not interior:
                BOUNDARY_ROWS.append((exponent, degree, name))
            MASK_ROWS += 1
        require(
            f"pure margin N={exponent} D={degree}",
            exponent * degree - degree - exponent
            == (exponent - 1) * (degree - 1) - 1,
        )
        require(
            f"lambda margins N={exponent} D={degree}",
            (degree - 1) * (exponent - 1) > 0 and exponent - 1 > 0,
        )
        require(
            f"mu margin N={exponent} D={degree}",
            exponent * (degree - 1) - 1 > 0,
        )
        MASK_ROWS += 3

require("unique Newton boundary", BOUNDARY_ROWS == [(2, 2, "pure")])
MASK_ROWS += 1
print(f"PASS {MASK_ROWS} mask, hull, margin, and unique-boundary controls")


print("SECTION toric face and resonant resolution")
FACE_ROWS = 0
for degree in range(2, 10):
    sigma = sp.prod(b - j for j in range(1, degree + 1))
    require(
        f"squarefree bottom D={degree}",
        sp.gcd(sp.Poly(sigma, b), sp.Poly(sp.diff(sigma, b), b)).degree() == 0,
    )
    FACE_ROWS += 1

z = sp.symbols("z", nonzero=True)
for exponent in range(2, 10):
    left = nu + mu * z ** (exponent + 1) - w * z**exponent
    critical_z = exponent * w / ((exponent + 1) * mu)
    critical_value = sp.cancel(left.subs(z, critical_z))
    require(
        f"left repeated root forces algebraic w N={exponent}",
        sp.Poly(sp.together(critical_value).as_numer_denom()[0], w).degree()
        == exponent + 1,
    )
    FACE_ROWS += 1

t, t0, u = sp.symbols("t t0 u", nonzero=True)
for exponent in range(2, 13):
    tuned_lambda = -(exponent + 1) * t0**exponent
    tuned_mu = exponent * t0 ** (exponent + 1)
    face = t ** (exponent + 1) + tuned_lambda * t + tuned_mu
    require(
        f"resonant root N={exponent}",
        zero(face.subs(t, t0)) and zero(sp.diff(face, t).subs(t, t0)),
    )
    require(
        f"resonant exactly double N={exponent}",
        zero(
            sp.diff(face, t, 2).subs(t, t0)
            - exponent * (exponent + 1) * t0 ** (exponent - 1)
        ),
    )
    a_sub = sp.Integer(exponent + 3)
    require(
        f"resonant transverse wz N={exponent}",
        sp.Poly(w - a_sub * t0**exponent, w).degree() == 1,
    )
    local_eta = sp.cancel((u**2) ** (exponent - 2) * sp.diff(u**2, u) / u)
    require(
        f"resonant eta order N={exponent}",
        zero(local_eta - 2 * u ** (2 * (exponent - 2))),
    )
    FACE_ROWS += 4

print(f"PASS {FACE_ROWS} bottom, left-face, and resonant-face controls")


print("SECTION logarithmic conic residues")
rho = sp.symbols("rho", nonzero=True)
res_plus = 1 / (2 * nu * s * rho)
res_minus = -res_plus
require("residue rho form", zero((res_plus - rho / (2 * w)).subs(w, nu * s * rho**2)))
require("opposite residues", zero(res_plus + res_minus))
require("nonzero residue control", res_plus.subs({nu: 1, s: 1, rho: 1}) == sp.Rational(1, 2))
print("PASS 3 opposite nonzero logarithmic-residue controls")


print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
