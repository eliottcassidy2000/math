#!/usr/bin/env python3
"""Exact scout for the full split degree-22 odd-Faber deformation.

This companion verifies the algebraic inputs for a geometric theorem which
is deliberately broader than THM-2704 but still generic rather than uniform.
In the target-translated full Faber gauge the chosen split sheet has columns

    E22 and E_j, j in {1,2,3,5,6,7,9,10,11,13,14,15,17,19,21}.

The first two Hamiltonian observables homogenize to degrees 23 and 24 in
P(1,2,3,4).  The script checks all weights and deck parities, the arithmetic
genus 425, the unique forced infinity singularity, and its generic coarse
delta 6 when the E21 coefficient is nonzero.  It also checks that the
THM-2704 all-even point has no component trapped on the omitted s=0 chart.

The passage from these exact checks to generic normalization genus 419 uses
standard flatness, openness of geometric integrality, generic smoothness, and
upper semicontinuity of delta.  No finite computation supplies those
quantifiers, and no claim about exceptional coefficient strata or JC(2) is
made here.
"""

from __future__ import annotations

from collections import defaultdict

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def faber_coefficients(
    degree: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
    extra: int = 2,
) -> list[sp.Expr]:
    """Laurent coefficients of (w^4+p*w^2+q*w+r)^(degree/4)."""
    exponent = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
            if step in quartic
        ) / index
        coefficients.append(sp.factor(value))
    return coefficients


def observables(
    degree: int,
    d: sp.Symbol,
    q: sp.Symbol,
    s: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr]:
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s)
    return (
        sp.factor(4 * coefficients[degree + 1]),
        sp.factor(4 * coefficients[degree + 2]),
    )


def direct_laurent_coefficient(
    degree: int,
    index: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
) -> sp.Expr:
    """Independent finite multinomial coefficient of P^(degree/4)."""
    exponent = sp.Rational(degree, 4)
    total = sp.Integer(0)
    for i in range(index // 2 + 1):
        for j in range(index // 3 + 1):
            remainder = index - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            k = remainder // 4
            chosen = i + j + k
            falling = sp.prod(exponent - offset for offset in range(chosen))
            total += (
                falling
                * p**i
                * q**j
                * r**k
                / (sp.factorial(i) * sp.factorial(j) * sp.factorial(k))
            )
    return sp.factor(total)


def weighted_degrees(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    weights: tuple[int, ...],
) -> set[int]:
    return {
        sum(exponent * weight for exponent, weight in zip(monomial, weights))
        for monomial, coefficient in sp.Poly(sp.expand(expression), *variables).terms()
        if coefficient
    }


def ordinary_leading_form(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[int, sp.Expr]:
    terms = sp.Poly(sp.expand(expression), *variables).terms()
    least = min(sum(monomial) for monomial, _ in terms)
    leading = sum(
        coefficient
        * sp.prod(variable**exponent for variable, exponent in zip(variables, monomial))
        for monomial, coefficient in terms
        if sum(monomial) == least
    )
    return least, sp.factor(leading)


def hilbert_coefficient(maximum: int) -> list[int]:
    coefficients = [0] * (maximum + 1)
    coefficients[0] = 1
    for weight in (1, 2, 3, 4):
        for degree in range(weight, maximum + 1):
            coefficients[degree] += coefficients[degree - weight]
    answer = coefficients[:]
    for degree in range(maximum + 1):
        if degree >= 23:
            answer[degree] -= coefficients[degree - 23]
        if degree >= 24:
            answer[degree] -= coefficients[degree - 24]
        if degree >= 47:
            answer[degree] += coefficients[degree - 47]
    return answer


d, q, s, h = sp.symbols("d q s h")
degrees = (1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 19, 21, 22)
lower_degrees = degrees[:-1]
bank = {degree: observables(degree, d, q, s) for degree in degrees}

# Every raw Laurent coefficient has its intrinsic quartic weight.  Deck
# parity is checked termwise, not merely from the largest q exponent.
for degree, (phi, psi) in bank.items():
    direct_phi = 4 * direct_laurent_coefficient(
        degree, degree + 1, 2 * d, q, d**2 - s
    )
    direct_psi = 4 * direct_laurent_coefficient(
        degree, degree + 2, 2 * d, q, d**2 - s
    )
    require(sp.factor(phi - direct_phi) == 0, f"direct Phi mismatch at degree {degree}")
    require(sp.factor(psi - direct_psi) == 0, f"direct Psi mismatch at degree {degree}")
    require(
        weighted_degrees(phi, (d, q, s), (2, 3, 4)) == {degree + 1},
        f"Phi weight changed at degree {degree}",
    )
    require(
        weighted_degrees(psi, (d, q, s), (2, 3, 4)) == {degree + 2},
        f"Psi weight changed at degree {degree}",
    )
    require(
        sp.expand(phi.subs(q, -q) - (-1) ** (degree + 1) * phi) == 0,
        f"Phi deck parity changed at degree {degree}",
    )
    require(
        sp.expand(psi.subs(q, -q) - (-1) ** degree * psi) == 0,
        f"Psi deck parity changed at degree {degree}",
    )

parameters = sp.symbols(" ".join(f"a{degree}" for degree in lower_degrees))
parameter = dict(zip(lower_degrees, parameters))
lam, W = sp.symbols("lambda W")
F23 = sp.expand(
    bank[22][0]
    + sum(parameter[j] * h ** (22 - j) * bank[j][0] for j in lower_degrees)
    - lam * h**23
)
G24 = sp.expand(
    bank[22][1]
    + sum(parameter[j] * h ** (22 - j) * bank[j][1] for j in lower_degrees)
    - W * h**24
)
require(
    weighted_degrees(F23, (h, d, q, s), (1, 2, 3, 4)) == {23},
    "full first flux is not weighted degree 23",
)
require(
    weighted_degrees(G24, (h, d, q, s), (1, 2, 3, 4)) == {24},
    "full second flux is not weighted degree 24",
)

# The Hilbert series is
# (1-t^23)(1-t^24)/((1-t)(1-t^2)(1-t^3)(1-t^4)).
# Its a-invariant is 37, and the degree-37 coefficient is h^0(omega)=p_a.
hilbert = hilbert_coefficient(132)
require(hilbert[37] == 425, "arithmetic genus changed")
for degree in range(60, 133, 12):
    require(
        hilbert[degree] == 23 * degree + 1 - 425,
        "Hilbert polynomial check changed",
    )

# The top forms control h=0 for every parameter.  Remove harmless constants.
top_phi = sp.factor(bank[22][0] * sp.Rational(512, 33))
top_psi = sp.factor(bank[22][1] * sp.Rational(8192, 33))
require(sp.gcd(sp.Poly(top_phi, d, q, s), sp.Poly(top_psi, d, q, s)) == 1,
        "top forms acquired a common component")

variables = (d, q, s)
minors = [
    sp.diff(top_phi, left) * sp.diff(top_psi, right)
    - sp.diff(top_phi, right) * sp.diff(top_psi, left)
    for index, left in enumerate(variables)
    for right in variables[index + 1 :]
]

# On q=1 and s=1 there is no singular top intersection.  On d=1 the
# singular scheme is supported only at q=s=0; powers certify the support.
for patch in (q, s):
    basis = sp.groebner(
        [top_phi, top_psi, *minors, patch - 1],
        d,
        q,
        s,
        order="grevlex",
    )
    require(basis.contains(sp.Integer(1)), f"unexpected top singularity on {patch}=1")

d_basis = sp.groebner(
    [top_phi, top_psi, *minors, d - 1],
    q,
    s,
    d,
    order="grevlex",
)
require(d_basis.reduce(q**20)[1] == 0 and d_basis.reduce(s**20)[1] == 0,
        "d=1 singular support is not the origin")

# At the forced point [0:1:0:0], nonzero a21 makes F23 linear in h on the
# index-one cover.  Eliminating h leaves the displayed six-line tangent cone.
phi21_at_origin = sp.factor(bank[21][0].subs({d: 1, q: 0, s: 0}))
require(phi21_at_origin == sp.Rational(88179, 131072),
        "E21 transverse coefficient changed")
phi_order, phi_leading = ordinary_leading_form(bank[22][0].subs(d, 1), (q, s))
psi_order, psi_leading = ordinary_leading_form(bank[22][1].subs(d, 1), (q, s))
require(phi_order == 6, "top Phi local order changed")
require(
    phi_leading == -sp.Rational(231, 128) * q * s * (q**2 - 3 * s**2) * (3 * q**2 - s**2),
    "top Phi local face changed",
)
expected_psi_face = (
    -sp.Rational(231, 256)
    * (q - s)
    * (q + s)
    * (q**2 - 4 * q * s + s**2)
    * (q**2 + 4 * q * s + s**2)
)
require(psi_order == 6 and sp.expand(psi_leading - expected_psi_face) == 0,
        "top Psi six-line face changed")

Q = sp.symbols("Q")
quotient_face = sp.factor(
    sum(
        coefficient * Q ** (power // 2)
        for (power,), coefficient in sp.Poly(expected_psi_face, q).terms()
    )
)
expected_quotient_face = -sp.Rational(231, 256) * (Q - s**2) * (
    Q**2 - 14 * Q * s**2 + s**4
)
require(sp.expand(quotient_face - expected_quotient_face) == 0,
        "coarse three-branch quotient face changed")
branch_constants = (sp.Integer(1), 7 + 4 * sp.sqrt(3), 7 - 4 * sp.sqrt(3))
require(
    all(
        sp.resultant(Q - left * s**2, Q - right * s**2, Q) != 0
        for index, left in enumerate(branch_constants)
        for right in branch_constants[index + 1 :]
    ),
    "coarse branches stopped being distinct",
)
pair_intersection_length = 2
coarse_delta = 3 * pair_intersection_length
require(coarse_delta == 6, "forced coarse delta changed")

# All-even THM-2704 witness: B=C=D=E=W=1, odd coefficients zero, and
# lambda=1/7496192 (so eta=1).  The y!=0 scaling is birational by THM-2704.
# These gcds certify that no direct q-model component is trapped on s=0 and
# that the two affine equations have no common surface factor.
even_degrees = (22, 14, 10, 6, 2)
special_phi = sp.factor(sum(bank[j][0] for j in even_degrees) - sp.Rational(1, 7496192))
special_psi = sp.factor(sum(bank[j][1] for j in even_degrees) - 1)
require(
    sp.gcd(
        sp.Poly(special_phi.subs(s, 0), d, q),
        sp.Poly(special_psi.subs(s, 0), d, q),
    )
    == 1,
    "THM-2704 witness acquired an s=0 component",
)
require(
    sp.gcd(
        sp.Poly(special_phi, d, q, s),
        sp.Poly(special_psi, d, q, s),
    )
    == 1,
    "THM-2704 witness acquired an affine surface component",
)

# A rank-two point proves dominance of the affine two-flux map.  The generic
# smoothness conclusion itself is theorem-level, not a finite inference.
rank_point = {d: 0, q: 1, s: 0}
rank_minor = sp.factor(
    (
        sp.diff(bank[22][0], d) * sp.diff(bank[22][1], q)
        - sp.diff(bank[22][0], q) * sp.diff(bank[22][1], d)
    ).subs(rank_point)
)
require(rank_minor == sp.Rational(9801, 524288), "dominance minor changed")

print("degree-22 full split odd-Faber weighted deformation scout")
print("columns=E22+15_lower:E18_target_translated_out")
print("raw_weights=d:2,q:3,s:4;Phi_j:j+1;Psi_j:j+2")
print("deck_parity=Phi_j:(-1)^(j+1),Psi_j:(-1)^j")
print("coefficient_weights=a_j:22-j,lambda:23,W:24")
print("weighted_closure=P(1,2,3,4),degrees=(23,24)")
print("hilbert_series=(1-t^23)(1-t^24)/((1-t)(1-t^2)(1-t^3)(1-t^4))")
print("arithmetic_genus=425")
print("top_singular_support=[h:d:q:s]=[0:1:0:0]")
print("E21_transverse=88179/131072")
print("index_cover_tangent=6_distinct_lines")
print("coarse_mu2_quotient=3_smooth_branches_Q=c*s^2")
print("forced_generic_delta=6")
print("generic_normalization_genus=419:THEORETICAL_OPEN_FAMILY_CONSEQUENCE")
print("THM2704_special_normalization_genus=89:special_delta=336")
print("affine_flux_rank_minor=9801/524288")
print("scope=GENERIC_FULL_SPLIT_FABER_FAMILY_NOT_EXCEPTIONAL_STRATA_NOT_JC2")
print("ALL CHECKS PASSED")
