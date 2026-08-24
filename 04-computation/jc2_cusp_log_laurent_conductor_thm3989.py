#!/usr/bin/env python3
"""Exact cusp-log Laurent conductor companion for THM-3989."""

from __future__ import annotations

import hashlib
import itertools
import json
from math import gcd

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    """Assertion-free exact gate retained under ``python -O``."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.cancel(expression)) == 0, message)


def jacobian(f: sp.Expr, g: sp.Expr,
             first: sp.Symbol, second: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(f, first) * sp.diff(g, second)
        - sp.diff(f, second) * sp.diff(g, first)
    )


s, tau, p, y = sp.symbols("s tau p y")
H = p**3 - y**2
p_log = s**2 + tau
y_log = s * p_log
x_log = s / tau
u_log = s**2 / tau
H_log = sp.expand(H.subs({p: p_log, y: y_log}))


# Panel I: exact cusp-log chart and logarithmic volume.
zero(H_log - tau * p_log**2, "H=tau*p^2")
zero(y_log / p_log - s, "s=y/p")
zero(H_log / p_log**2 - tau, "tau=H/p^2")
zero(x_log - y_log * p_log / H_log, "x=yp/H")
zero(u_log - y_log**2 / H_log, "u=y^2/H")
gate(sp.det(sp.Matrix([
    [sp.diff(x_log, s), sp.diff(x_log, tau)],
    [0, 1],
])) == 1 / tau, "dx wedge dtau=ds wedge dtau/tau")


# Panel II: exact negative-coefficient conductor law.  A generator monomial
# x^a*u^b*p^c*y^e has the displayed closed Laurent form.
def generator_monomial(a: int, b: int, c: int, e: int) -> sp.Expr:
    return sp.expand(x_log**a * u_log**b * p_log**c * y_log**e)


coefficient_controls: dict[str, int] = {}
for total in range(0, 8):
    tested = 0
    for a, b, c, e in itertools.product(range(total + 1), repeat=4):
        if a + b + c + e > total:
            continue
        expression = generator_monomial(a, b, c, e)
        closed = sp.expand(
            s**(a + 2 * b + e) * tau**(-a - b)
            * (s**2 + tau)**(c + e)
        )
        zero(expression - closed,
             f"monomial ({a},{b},{c},{e}): Laurent form")
        for j in range(1, a + b + 1):
            coefficient = sp.expand(expression).coeff(tau, -j)
            if coefficient != 0:
                quotient = sp.cancel(coefficient / s**j)
                gate(quotient.is_polynomial(s),
                     f"monomial ({a},{b},{c},{e}): s^{j} conductor")
        tested += 1
    coefficient_controls[str(total)] = tested


# Every s^(d+n), n>=0, occurs as a depth-d principal coefficient.  The
# semigroup <2,3> supplies all n>=2; n=1 uses x^(d-1)u.
def semigroup_23(n: int) -> tuple[int, int]:
    if n % 2 == 0:
        return n // 2, 0
    return (n - 3) // 2, 1


surjectivity_controls: dict[str, list[int]] = {}
for d in range(1, 10):
    offsets: list[int] = []
    for n in range(0, 19):
        if n == 0:
            witness = x_log**d
        elif n == 1:
            witness = x_log**(d - 1) * u_log
        else:
            c, e = semigroup_23(n)
            gate(c >= 0 and e >= 0 and 2 * c + 3 * e == n,
                 f"offset {n}: cusp semigroup witness")
            witness = x_log**d * p_log**c * y_log**e
        zero(sp.expand(witness).coeff(tau, -d) - s**(d + n),
             f"depth {d}, offset {n}: principal-symbol witness")
        offsets.append(n)
    surjectivity_controls[str(d)] = offsets


# Finite hostile for the depth-zero intersection: within generator degree
# <=4, cancel every negative Laurent coefficient and check that the tau=0
# residues span exactly the visible cusp exponents 0,2,3,... .
gens = (x_log, 1 + u_log, p_log, y_log)
intersection_controls: dict[int, dict[str, object]] = {}
for cutoff in range(1, 5):
    monomials: list[sp.Expr] = []
    for exponents in itertools.product(range(cutoff + 1), repeat=4):
        if sum(exponents) <= cutoff:
            monomials.append(sp.expand(sp.prod(
                generator**exponent
                for generator, exponent in zip(gens, exponents)
            )))
    negative_powers = sorted({
        int(term.as_powers_dict().get(tau, 0))
        for expression in monomials
        for term in sp.Add.make_args(expression)
        if term.as_powers_dict().get(tau, 0) < 0
    })
    rows: list[list[sp.Expr]] = []
    for power in negative_powers:
        polynomials = [expression.coeff(tau, power)
                       for expression in monomials]
        degree = max((sp.degree(poly, s) if poly else -1)
                     for poly in polynomials)
        for exponent in range(degree + 1):
            rows.append([poly.coeff(s, exponent) for poly in polynomials])
    matrix = sp.Matrix(rows) if rows else sp.zeros(0, len(monomials))
    nullspace = matrix.nullspace()
    residues = [sp.expand(sum(
        vector[index] * monomials[index].coeff(tau, 0)
        for index in range(len(monomials))
    )) for vector in nullspace]
    max_degree = max((sp.degree(poly, s) if poly else -1)
                     for poly in residues)
    residue_matrix = sp.Matrix([
        [poly.coeff(s, exponent) for poly in residues]
        for exponent in range(max_degree + 1)
    ])
    rank = residue_matrix.rank()
    visible: list[int] = []
    for exponent in range(max_degree + 1):
        unit = sp.zeros(max_degree + 1, 1)
        unit[exponent] = 1
        if residue_matrix.row_join(unit).rank() == rank:
            visible.append(exponent)
    gate(1 not in visible, f"cutoff {cutoff}: normalization parameter absent")
    gate(visible == [0] + list(range(2, max_degree + 1)),
         f"cutoff {cutoff}: exact cusp-semigroup residue span")
    intersection_controls[cutoff] = {
        "monomials": len(monomials),
        "residue_rank": rank,
        "visible_exponents": visible,
    }


# Panel III: complete first-depth cancellation syzygy.
l00, l10, l01, m00, m10, m01 = sp.symbols(
    "l00 l10 l01 m00 m10 m01"
)
L_poly = l00 + l10 * p + l01 * y
M_poly = m00 + m10 * p + m01 * y
G_poly = sp.expand(p**2 * L_poly + y * M_poly)
K_poly = sp.expand(-y * L_poly - p * M_poly)
zero(p * G_poly + y * K_poly - H * L_poly,
     "first-depth syzygy pG+yK=HL")
zero((y * p / H) * G_poly + (y**2 / H) * K_poly - y * L_poly,
     "first-depth cancellation descends to yL")


# Panel IV: exact Laurent bracket convolution and scalar moment.
amin, amax = -3, 4
cmin, cmax = -4, 3
a_coefficients = {
    i: (i + 5) * s**3 + (2 * i + 9) * s**2 + (i * i + 7) * s + 1
    for i in range(amin, amax + 1)
}
c_coefficients = {
    j: (j * j + 11) * s**3 + (3 * j + 13) * s**2 + (j + 8) * s + 2
    for j in range(cmin, cmax + 1)
}
A_generic = sp.Add(*(poly * tau**i for i, poly in a_coefficients.items()))
C_generic = sp.Add(*(poly * tau**j for j, poly in c_coefficients.items()))
log_bracket = sp.expand(tau * jacobian(A_generic, C_generic, s, tau))
for weight in range(amin + cmin, amax + cmax + 1):
    convolution = sp.expand(sum(
        j * sp.diff(a_coefficients[i], s) * c_coefficients[j]
        - i * a_coefficients[i] * sp.diff(c_coefficients[j], s)
        for i in a_coefficients for j in c_coefficients if i + j == weight
    ))
    zero(log_bracket.coeff(tau, weight) - convolution,
         f"Laurent bracket convolution weight {weight}")

moment = sp.expand(sum(
    i * a_coefficients[i] * c_coefficients.get(-i, 0)
    for i in a_coefficients
))
zero(log_bracket.coeff(tau, 0) + sp.diff(moment, s),
     "scalar coefficient is negative moment derivative")


# The tame automorphism hostile realizes the law but leaves B2 after the
# depth-lowering shear: A=x, C=tau+x^q.
shear_controls: dict[int, dict[str, object]] = {}
for q_value in range(2, 8):
    A_control = x_log
    C_control = sp.expand(tau + x_log**q_value)
    zero(tau * jacobian(A_control, C_control, s, tau) - 1,
         f"q={q_value}: log Keller hostile")
    control_moment = sp.expand(sum(
        i * sp.expand(A_control).coeff(tau, i)
        * sp.expand(C_control).coeff(tau, -i)
        for i in range(-q_value, q_value + 1)
    ))
    zero(control_moment + s, f"q={q_value}: exact scalar moment -s")
    zero(C_control - A_control**q_value - tau,
         f"q={q_value}: depth-lowering target shear")
    shear_controls[q_value] = {
        "depths": [1, q_value],
        "remainder": "tau",
    }


# Panel V: UFD common-power leading symbols and the s^g conductor.
alpha, beta = sp.symbols("alpha beta", nonzero=True)
leading_controls: dict[str, dict[str, int]] = {}
for d, e in ((1, 2), (2, 3), (2, 5), (4, 6), (6, 10), (6, 15)):
    common = gcd(d, e)
    h = s**common * (s**2 + 2 * s + 3)
    a_lead = alpha * h**(d // common)
    c_lead = beta * h**(e // common)
    zero(d * a_lead * sp.diff(c_lead, s)
         - e * sp.diff(a_lead, s) * c_lead,
         f"depths {d},{e}: leading common-power equation")
    gate(sp.rem(a_lead, s**d, s) == 0,
         f"depths {d},{e}: A principal conductor")
    gate(sp.rem(c_lead, s**e, s) == 0,
         f"depths {d},{e}: C principal conductor")
    leading_controls[f"{d}:{e}"] = {
        "gcd": common,
        "reduced_d": d // common,
        "reduced_e": e // common,
    }

bad_a = s**2 * (s + 1)
bad_c = s**3 * (s + 2)
gate(sp.factor(2 * bad_a * sp.diff(bad_c, s)
               - 3 * sp.diff(bad_a, s) * bad_c) != 0,
     "2:3 hostile without a common base fails")


summary = {
    "checks": CHECKS,
    "chart": "s=y/p=xt,tau=H/p^2=t; volume=ds^dtau/tau",
    "depth_zero": "B2 intersect k[s,tau]=k[p,y]",
    "negative_conductor": "[tau^-d]B2=s^d*k[s] for d>=1",
    "coefficient_controls": coefficient_controls,
    "surjectivity_controls": surjectivity_controls,
    "intersection_controls": intersection_controls,
    "first_depth_syzygy": "pole cancellation iff xG+uK=yL",
    "coefficient_law": (
        "sum_(i+j=k)(j*a_i'*c_j-i*a_i*c_j')=delta_(k,0)"
    ),
    "scalar_moment": "sum_i i*a_i*c_-i=-s",
    "leading_law": "a_-d=alpha*h^(d/g),c_-e=beta*h^(e/g),s^g|h",
    "reduction": "target shears leave positive nondividing depths; first 2:3",
    "shear_controls": shear_controls,
    "leading_controls": leading_controls,
    "scope": "height-two completion reduction; 2:3 and JC2 remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3989 cusp-log Laurent conductor companion")
print(f"CHECKS={CHECKS}")
print("CHART=S=Y/P=XT;TAU=H/P^2=T;VOLUME=DS_WEDGE_DTAU/TAU")
print("DEPTH_ZERO=B2_INTERSECT_K[S,TAU]=K[P,Y]")
print("NEGATIVE_COEFFICIENT=[TAU^-D]B2=S^D*K[S]")
print("FIRST_DEPTH_CANCELLATION=DESCENDS_TO_CUSP_PLANE")
print("SCALAR_MOMENT=SUM_I_I*A_I*C_-I=-S")
print("LEADING_SYMBOLS=COMMON_POWERS_WITH_S^G_BASE")
print("REDUCED_DEPTHS=POSITIVE_NONDIVIDING;FIRST_ARITHMETIC_TYPE_2_TO_3")
print("SCOPE=HEIGHT_TWO_COMPLETION_REDUCTION;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
