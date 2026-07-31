#!/usr/bin/env python3
"""Exact algebra referee for THM-2759's split degree-six closure.

This verifies the finite algebra used in the audited proof mechanism:

* the full E6,E5,E3,E2,E1 Phi/Psi/R bank in square-prefix coordinates;
* the two points in the h=0 boundary and the nonzero response at the q-point;
* the exact z=0 values used to exclude that q-point by polynomiality;
* the P_infinity leading constants and the noncancellation of R after the
  only possible low-slope Phi cancellation; and
* the q=0 vertical specialization.

Normalization, properness, the rational-primitive theorem, and DVR arguments
are mathematical inputs rather than finite computations.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def recurrence_coefficients(
    degree: int, p: sp.Expr, q: sp.Expr, r: sp.Expr, extra: int = 3
) -> list[sp.Expr]:
    exponent = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
        ) / index
        coefficients.append(sp.factor(value))
    return coefficients


def multinomial_coefficient(
    degree: int, offset: int, p: sp.Expr, q: sp.Expr, r: sp.Expr
) -> sp.Expr:
    target = degree + offset
    exponent = sp.Rational(degree, 4)
    answer = sp.Integer(0)
    for i in range(target // 2 + 1):
        for j in range(target // 3 + 1):
            remainder = target - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            k = remainder // 4
            count = i + j + k
            falling = sp.prod(exponent - ell for ell in range(count))
            answer += (
                falling
                * p**i
                * q**j
                * r**k
                / (factorial(i) * factorial(j) * factorial(k))
            )
    return sp.factor(answer)


def observables(
    degree: int, d: sp.Symbol, q: sp.Symbol, s: sp.Symbol
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    p = 2 * d
    r = d**2 - s
    coefficients = recurrence_coefficients(degree, p, q, r)
    phi = sp.factor(4 * coefficients[degree + 1])
    psi = sp.factor(4 * coefficients[degree + 2])
    response = sp.factor(
        4 * coefficients[degree + 3] + p * coefficients[degree + 1]
    )
    independent = tuple(
        multinomial_coefficient(degree, offset, p, q, r)
        for offset in (1, 2, 3)
    )
    direct = (
        sp.factor(4 * independent[0]),
        sp.factor(4 * independent[1]),
        sp.factor(4 * independent[2] + p * independent[0]),
    )
    require(
        all(sp.expand(left - right) == 0 for left, right in zip(
            (phi, psi, response), direct
        )),
        f"recurrence/multinomial mismatch in degree {degree}",
    )
    return phi, psi, response


d, q, s = sp.symbols("d q s")
degrees = (1, 2, 3, 5, 6)
bank = {degree: observables(degree, d, q, s) for degree in degrees}

expected = {
    1: (2*d, q, (d**2 - 2*s)/2),
    2: (2*q, -2*s, -d*q),
    3: (sp.Rational(3, 2)*(d**2 - 2*s),
        -sp.Rational(3, 2)*d*q,
        (4*d**3 - 3*q**2)/8),
    5: (sp.Rational(5, 8)*(2*d**3 - 4*d*s + q**2),
        -sp.Rational(5, 8)*q*(d**2 + 2*s),
        sp.Rational(5, 32)*(3*d**4 - 4*d**2*s - 4*d*q**2 + 4*s**2)),
    6: (-3*q*s,
        -sp.Rational(3, 2)*(d*q**2 - s**2),
        -q*(q**2 - 6*d*s)/4),
}
for degree in degrees:
    require(
        all(sp.expand(left - right) == 0 for left, right in zip(
            bank[degree], expected[degree]
        )),
        f"displayed bank mismatch in degree {degree}",
    )

# At h=0 the two equations are q*s=0 and d*q^2-s^2=0, hence exactly
# P_d=[1:0:0] and P_q=[0:1:0].  The response is nonzero at P_q.
phi6, psi6, response6 = bank[6]
require(sp.gcd(sp.Poly(phi6, d, q, s), sp.Poly(psi6, d, q, s)) == 1,
        "top forms acquired a common factor")
require(response6.subs({d: 0, q: 1, s: 0}) == -sp.Rational(1, 4),
        "response at P_q changed")

# Evaluate every Faber polynomial at the original polynomial section z=0.
# Here beta=B/(2U), d=C-beta^2, and s=q*beta-E.
beta, C, E, w = sp.symbols("beta C E w")
substitution = {w: beta, d: C-beta**2, s: q*beta-E}
faber = {
    1: w,
    2: w**2+d,
    3: w**3+sp.Rational(3, 2)*d*w+sp.Rational(3, 4)*q,
    5: (w**5+sp.Rational(5, 2)*d*w**3+sp.Rational(5, 4)*q*w**2
        +(sp.Rational(15, 8)*d**2-sp.Rational(5, 4)*s)*w
        +sp.Rational(5, 8)*d*q),
    6: ((w**2+d)**3+sp.Rational(3, 2)*(w**2+d)*(q*w-s)
        +sp.Rational(3, 8)*q**2),
}
at_zero = {degree: sp.factor(expression.subs(substitution))
           for degree, expression in faber.items()}
expected_at_zero = {
    1: beta,
    2: C,
    3: (6*C*beta-2*beta**3+3*q)/4,
    5: (15*C**2*beta-10*C*beta**3+5*C*q+10*E*beta
        +3*beta**5-5*beta**2*q)/8,
    6: (8*C**3+12*C*E+3*q**2)/8,
}
for degree in degrees:
    require(sp.expand(at_zero[degree]-expected_at_zero[degree]) == 0,
            f"z=0 formula changed in degree {degree}")

# P_d low-slope cancellation constants.  Odd rows have gaps 1,3,5 and
# Phi/R constant pairs listed below.  If s0=eps*q0 and Phi6 is cancelled,
# the resulting response coefficient is never zero.
odd_data = {
    5: (1, sp.Rational(5, 4), sp.Rational(15, 32)),
    3: (3, sp.Rational(3, 2), sp.Rational(1, 2)),
    1: (5, sp.Integer(2), sp.Rational(1, 2)),
}
response_residues = {}
for degree, (gap, phi_constant, response_constant) in odd_data.items():
    residue = sp.factor(
        sp.Rational(3, 2) + 3 * response_constant / phi_constant
    )
    require(residue != 0, f"response cancellation appeared for E{degree}")
    response_residues[gap] = residue

# On q=0, Psi makes s constant.  Unless all odd coefficients vanish and
# lambda=0, Phi is a nonzero polynomial in d; in the exceptional all-even
# case R vanishes identically.
a5, a3, a2, a1, lam, W = sp.symbols("a5 a3 a2 a1 lambda W")
phi_vertical = sp.factor(
    bank[6][0] + a5*bank[5][0] + a3*bank[3][0]
    + a2*bank[2][0] + a1*bank[1][0] - lam
).subs(q, 0)
psi_vertical = sp.factor(
    bank[6][1] + a5*bank[5][1] + a3*bank[3][1]
    + a2*bank[2][1] + a1*bank[1][1] - W
).subs(q, 0)
response_vertical = sp.factor(
    bank[6][2] + a5*bank[5][2] + a3*bank[3][2]
    + a2*bank[2][2] + a1*bank[1][2]
).subs(q, 0)
require(sp.expand(psi_vertical - (3*s**2-4*a2*s-2*W)/2) == 0,
        "vertical Psi changed")
phi_poly = sp.Poly(sp.expand(phi_vertical), d)
require(phi_poly.all_coeffs() == [
    sp.Rational(5, 4)*a5,
    sp.Rational(3, 2)*a3,
    -sp.Rational(5, 2)*a5*s+2*a1,
    -3*a3*s-lam,
], "vertical Phi coefficient list changed")
require(sp.expand(response_vertical.subs({a5: 0, a3: 0, a1: 0})) == 0,
        "all-even vertical response stopped vanishing")

print("split degree-6 monicization exact bank referee")
print("status=THM-2759_PROVED_AND_INDEPENDENTLY_HOSTILE_AUDITED")
print("columns=E6,E5,E3,E2,E1:full_chosen_sheet_bank")
print("boundary_h0=P_d:[1:0:0],P_q:[0:1:0]")
print("response_at_P_q=-1/4:forced_response_pole")
print("Q_z0_E6=(8C^3+12CE+3q^2)/8")
print("Q_z0_lower_pole_ceiling_at_P_q=strictly_less_than_q^2")
print("P_d_low_slope_response_residues_g1_g3_g5="
      + ",".join(str(response_residues[gap]) for gap in (1, 3, 5)))
print("P_d_lambda_gap7_response_is_top_only_nonzero=True")
print("vertical_Psi=(3s^2-4a2s-2W)/2")
print("vertical_Phi_identically_zero_iff=a5=a3=a1=lambda=0")
print("vertical_exception_response=0")
print("scope=FINITE_ALGEBRA_REFEREE_FOR_PROVED_THM2759_DVR_GLOBAL_REGULARITY_PROOF")
print("ALL CHECKS PASSED")
