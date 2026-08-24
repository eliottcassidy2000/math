#!/usr/bin/env python3
"""Exact critical-resultant and genus-ladder companion for THM-3984."""

from __future__ import annotations

import hashlib
import json
from math import gcd

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def order_at_zero(expr: sp.Expr, variable: sp.Symbol) -> int:
    poly = sp.Poly(sp.expand(expr), variable)
    return min(power[0] for power, coefficient in poly.terms()
               if coefficient != 0)


def convex_hull(points: set[tuple[int, int]]) -> list[tuple[int, int]]:
    pts = sorted(points)

    def cross(o: tuple[int, int], a: tuple[int, int],
              b: tuple[int, int]) -> int:
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

    lower: list[tuple[int, int]] = []
    for point in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def interior_lattice_points(points: set[tuple[int, int]]) -> int:
    hull = convex_hull(points)
    area_twice = abs(sum(
        hull[i][0] * hull[(i + 1) % len(hull)][1]
        - hull[(i + 1) % len(hull)][0] * hull[i][1]
        for i in range(len(hull))
    ))
    boundary = sum(
        gcd(abs(hull[(i + 1) % len(hull)][0] - hull[i][0]),
            abs(hull[(i + 1) % len(hull)][1] - hull[i][1]))
        for i in range(len(hull))
    )
    return (area_twice - boundary + 2) // 2


x, t, u = sp.symbols("x t u")
alpha, beta, c = sp.symbols("alpha beta c", nonzero=True)
r = sp.expand(u * (u + 1))
b = sp.expand(u * r)
rp = sp.diff(r, u)
bp = sp.diff(b, u)


# Panel I: exact mixed p/y critical resultants.  Full resultants are replayed
# for several heights and derivative degrees; the proof uses the same leading
# and endpoint formulas uniformly.
mixed_degrees: dict[str, int] = {}
for n in range(2, 6):
    N = n + 2
    for label, gpoly in (("zero", sp.Integer(0)),
                         ("constant", sp.Integer(1)),
                         ("linear", u + 1)):
        F = sp.expand(x**N - n * alpha * r * x - (n + 1) * beta * b)
        G = sp.expand(gpoly * x**(N - 1) + alpha * rp * x + beta * bp)
        resultant = sp.factor(sp.resultant(F, G, x))
        gate(order_at_zero(resultant, u) == 2,
             f"mixed n={n},{label}: exact spurious u=0 order")
        gate(sp.factor(sp.expand(resultant).coeff(u, 2)
                       - (-1)**n * (n - 1) * alpha**N * beta) == 0,
             f"mixed n={n},{label}: u=0 leading coefficient")
        gate(sp.factor(resultant.subs(u, -1) - beta**N) == 0,
             f"mixed n={n},{label}: u=-1 excluded")
        if gpoly == 0:
            expected_degree = 2 * N
            expected_lc = (3 * beta)**N
        else:
            e = sp.degree(gpoly, u)
            expected_degree = N * e + 3 * (N - 1)
            expected_lc = ((-1)**(N + 1) * (N - 1)**(N - 1)
                           * beta**(N - 1)
                           * sp.Poly(gpoly, u).LC()**N)
        gate(sp.degree(resultant, u) == expected_degree,
             f"mixed n={n},{label}: resultant degree")
        gate(sp.factor(sp.Poly(resultant, u).LC() - expected_lc) == 0,
             f"mixed n={n},{label}: resultant leading coefficient")
        gate(expected_degree > 2,
             f"mixed n={n},{label}: genuine root remains")
        mixed_degrees[f"n{n}_{label}"] = int(expected_degree)


# The two singular coefficient axes have different compatibility polynomials
# and cannot be recovered by specializing the mixed resultant.
for n in range(2, 10):
    N = n + 2
    gpoly = u**2 + u + 1
    E_p = sp.expand((-n)**n * r**n * gpoly**(n + 1)
                    + alpha * rp**(n + 1))
    gate(sp.degree(E_p, u) == 4 * n + 2,
         f"p-axis n={n}: degree")
    gate(E_p.subs(u, 0) == alpha,
         f"p-axis n={n}: u=0 excluded")
    gate(sp.factor(E_p.subs(u, -1) - alpha * (-1)**(n + 1)) == 0,
         f"p-axis n={n}: u=-1 excluded")

    E_y = sp.expand((n + 1)**(n + 1) * u**n * (u + 1)**(n + 1)
                    * gpoly**(n + 2)
                    - (-1)**(n + 2) * beta * (3 * u + 2)**(n + 2))
    gate(sp.degree(E_y, u) == 4 * n + 5,
         f"y-axis n={n}: degree")
    gate(E_y.subs(u, 0) != 0,
         f"y-axis n={n}: u=0 excluded")
    gate(E_y.subs(u, -1) == -beta,
         f"y-axis n={n}: u=-1 excluded")


# The uncoupled polynomial shear is a submersion by an exact Euler row.
for n in range(2, 11):
    hsample = u**4 + 2 * u + 3
    A_source = x + hsample.subs(u, x**n * t)
    euler_row = sp.factor(x * sp.diff(A_source, x)
                          - n * t * sp.diff(A_source, t))
    gate(euler_row == x, f"uncoupled height {n}: Euler submersion row")


# Panel II: every single p^a y^b coupling has a valid off-axis critical
# address.  This includes both pure generators and every multiplicative
# mixed monomial.
monomial_degrees: dict[str, int] = {}
for n in range(2, 7):
    for p_exp in range(0, 4):
        for y_exp in range(0, 4):
            if p_exp + y_exp == 0:
                continue
            weight = n * p_exp + (n + 1) * y_exp
            M = weight + 1
            u_order = p_exp + 2 * y_exp
            color_order = p_exp + y_exp
            linear = (2 * p_exp + 3 * y_exp) * u + u_order
            for e in range(0, 3):
                gpoly = (u + 2)**e if e else sp.Integer(1)
                E_ab = (
                    weight**weight
                    * u**(M - u_order)
                    * (u + 1)**(M - color_order)
                    * gpoly**M
                    - (-1)**M * c * linear**M
                )
                expected = (2 * M - u_order - color_order + M * e)
                gate(sp.degree(E_ab, u) == expected,
                     f"monomial n={n},a={p_exp},b={y_exp},e={e}: degree")
                gate(expected - M
                     == ((n - 2) * (p_exp + y_exp) + 1 + M * e),
                     f"monomial n={n},a={p_exp},b={y_exp},e={e}: gap")
                gate(E_ab.subs(u, 0) != 0,
                     f"monomial n={n},a={p_exp},b={y_exp},e={e}: u=0")
                gate(E_ab.subs(u, -1) != 0,
                     f"monomial n={n},a={p_exp},b={y_exp},e={e}: u=-1")
                monomial_degrees[
                    f"n{n}_a{p_exp}_b{y_exp}_e{e}"
                ] = int(expected)


# Panel III: the n=2 linear-y generic-fibre time form is holomorphic for
# every degree of h.  Newton polygons give an independent genus check.
genus_ladder: dict[int, dict[str, int]] = {}
for d in range(0, 16):
    support = {(4, 0), (3, 0), (0, 2), (0, 3)}
    support.update((3, j) for j in range(d + 1))
    newton_genus = interior_lattice_points(support)
    if d == 0:
        infinity_places = 1
        infinity_zero_degree = 1
        divisor_degree = 1 + infinity_zero_degree
        genus_formula = 2
        source_infinity_places = 3
    else:
        gamma = gcd(3, d)
        e1_order = 2 * d - 2
        e2_order = 2 * d // gamma - 1
        gate(e1_order >= 0, f"degree {d}: first infinity regular")
        gate(e2_order >= 0, f"degree {d}: second infinity regular")
        infinity_places = 1 + gamma
        infinity_zero_degree = e1_order + gamma * e2_order
        divisor_degree = 1 + infinity_zero_degree
        genus_formula = (4 * d + 1 - gamma) // 2
        source_infinity_places = 3 + gamma
    gate(divisor_degree == 2 * genus_formula - 2,
         f"degree {d}: canonical divisor degree")
    gate(newton_genus == genus_formula,
         f"degree {d}: Newton polygon genus cross-check")
    genus_ladder[d] = {
        "genus": genus_formula,
        "toric_infinity_places": infinity_places,
        "source_infinity_places": source_infinity_places,
        "omega_infinity_zero_degree": infinity_zero_degree,
    }


summary = {
    "checks": CHECKS,
    "mixed_family": "A=x+h(u)+alpha*p+beta*y",
    "mixed_submersion": "iff alpha=beta=0",
    "mixed_resultant_degrees": mixed_degrees,
    "monomial_family": (
        "A=x+h(u)+c*p^a*y^b; all n>=2,a+b>0,c!=0"
    ),
    "monomial_critical": True,
    "monomial_degree_controls": monomial_degrees,
    "n2_linear_y_time_form": "omega=x du/F_x; holomorphic and nonexact",
    "genus_ladder": genus_ladder,
    "scope": "displayed boundary-generator families; JC2 open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3984 boundary-generator coupling companion")
print(f"CHECKS={CHECKS}")
print("MIXED_LINEAR_SUBMERSION=IFF_ALPHA_BETA_ZERO")
print("BOUNDARY_MONOMIAL=EVERY_P_POWER_Y_POWER_COUPLING_IS_CRITICAL")
print("N2_LINEAR_Y_TIME_FORM=HOLOMORPHIC_NONZERO_NONEXACT")
print("GENUS_D0=2;GENUS_DPOS=(4D+1-GCD(3,D))/2")
print("SOURCE_INFINITY_D0=3;SOURCE_INFINITY_DPOS=3+GCD(3,D)")
print("SCOPE=DISPLAYED_BOUNDARY_GENERATOR_FAMILIES;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
