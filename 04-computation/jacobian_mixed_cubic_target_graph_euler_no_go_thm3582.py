#!/usr/bin/env python3
"""Exact symbolic audit for THM-3582's mixed-cubic Euler no-go.

The finite atlases below work over QQ after the weighted normalizations used
in the proof.  A truth gate is never an ``assert`` so that ``python -O`` runs
exactly the same audit.
"""

from __future__ import annotations

import math
import os
from functools import reduce

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


def progress(label: str) -> None:
    if os.environ.get("THM3582_PROGRESS") == "1":
        print(f"progress: {label}", flush=True)


def numerator(expr: sp.Expr) -> sp.Expr:
    """Expanded numerator after exact cancellation."""

    return sp.expand(sp.cancel(expr).as_numer_denom()[0])


def coeff(poly: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(poly, variable).coeff_monomial(variable**degree)


def rho(poly: sp.Expr, variable: sp.Symbol) -> int:
    """Number of distinct complex roots of a nonzero QQ polynomial."""

    P = sp.Poly(sp.expand(poly), variable, domain=sp.QQ)
    require(not P.is_zero, "rho called on zero polynomial")
    return sp.Poly(P.sqf_part(), variable).degree()


def gcd_many(polys: list[sp.Expr], variable: sp.Symbol) -> sp.Expr:
    nonzero = [sp.Poly(sp.expand(p), variable, domain=sp.QQ) for p in polys if p != 0]
    require(bool(nonzero), "nonempty polynomial gcd")
    return reduce(lambda x, y: sp.gcd(x, y), nonzero).as_expr()


def is_unit_ideal(
    equations: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    invert: sp.Expr = sp.Integer(1),
) -> bool:
    """Exact Rabinowitsch saturation followed by a QQ Groebner basis."""

    rab = sp.Dummy("rab")
    rows = [numerator(e) for e in equations]
    if sp.expand(invert) != 1:
        rows.append(sp.expand(rab * invert - 1))
        gens = (rab,) + variables
    else:
        gens = variables
    basis = sp.groebner(rows, *gens, order="grevlex", method="f5b")
    return any(p.total_degree() == 0 for p in basis.polys)


def strip_exact_factors(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    factors: tuple[sp.Expr, ...],
) -> sp.Expr:
    """Remove only explicitly listed exact polynomial factors."""

    P = sp.Poly(numerator(expr), *variables)
    for factor in factors:
        divisor = sp.Poly(factor, *variables)
        while True:
            quotient, remainder = sp.div(P, divisor)
            if remainder.is_zero:
                P = quotient
            else:
                break
    return P.as_expr()


def roots_are_allowed(poly: sp.Poly, allowed: sp.Expr, variable: sp.Symbol) -> bool:
    """Whether every root of ``poly`` is a root of the allowed product."""

    if poly.is_zero:
        return False
    radical = sp.Poly(poly.sqf_part(), variable).monic()
    return sp.rem(sp.Poly(allowed, variable), radical).is_zero


def binomial_target(
    constant: sp.Expr,
    slopes_and_multiplicities: tuple[tuple[sp.Expr, int], ...],
) -> sp.Expr:
    answer = constant
    for slope, multiplicity in slopes_and_multiplicities:
        answer *= (1 + slope * s) ** multiplicity
    return sp.expand(answer)


# ---------------------------------------------------------------------------
# Universal mixed-cubic family and strict Jelonek transform.
# ---------------------------------------------------------------------------

a, b, c, s = sp.symbols("a b c s")
A, B, C, D, E, F = sp.symbols("A B C D E F")

phi = a * (A * b**2 + B * b + C) + D * b**3 + E * b**2 + F * b + C / 4
require(sp.expand(phi.subs({a: -sp.Rational(1, 4), b: 0})) == 0, "collision value")

a_sheet = s * (2 * b - s) / 12
c_sheet = 4 * (-b + 2 * s) / (3 * s**2)
G = sp.expand(12 * s**2 * (c_sheet + phi.subs(a, a_sheet)))

g3 = 2 * s**2 * (A * s + 6 * D)
g2 = s**2 * (-A * s**2 + 2 * B * s + 12 * E)
g1 = -B * s**4 + 2 * C * s**3 + 12 * F * s**2 - 16
g0 = s * (-C * s**3 + 3 * C * s + 32)
require(sp.expand(G - (g3 * b**3 + g2 * b**2 + g1 * b + g0)) == 0, "strict transform")
require(sp.expand(G.subs(s, 0) + 16 * b) == 0, "sheet boundary G(0,b)=-16b")

Delta_formula = sp.expand(
    g2**2 * g1**2
    - 4 * g3 * g1**3
    - 4 * g2**3 * g0
    - 27 * g3**2 * g0**2
    + 18 * g3 * g2 * g1 * g0
)
Delta = sp.expand(sp.discriminant(G, b))
require(sp.expand(Delta - Delta_formula) == 0, "binary-cubic discriminant")
R, delta_remainder = sp.div(Delta, s**2, domain=sp.QQ[A, B, C, D, E, F])
require(delta_remainder == 0, "Delta=s^2 R")
R = sp.expand(R)
require(sp.Poly(R, s).degree() <= 14, "degree R<=14")
r = {j: coeff(R, s, j) for j in range(15)}

key_rows = {
    0: 196608 * D,
    1: 32768 * A,
    2: 36864 * (E**2 - 12 * D * F),
    3: -12288 * (6 * A * F - B * E + 6 * C * D + 108 * D * E),
    14: A**2 * (B**2 - 4 * A * C),
    13: 4 * A * B * (B**2 - 4 * A * C),
    12: 4
    * (
        3 * A**3 * C
        - 6 * A**2 * B * F
        - 8 * A**2 * C**2
        + 36 * A**2 * C * E
        - 2 * A * B**2 * C
        - 6 * A * B**2 * E
        - 54 * A * B * C * D
        + B**4
        + 12 * B**3 * D
    ),
    11: 4
    * (
        32 * A**3
        + 9 * A**2 * B * C
        + 120 * A**2 * C * F
        - 48 * A * B**2 * F
        - 16 * A * B * C**2
        - 12 * A * B * C * E
        - 216 * A * C**2 * D
        + 4 * B**3 * C
        + 12 * B**3 * E
        + 36 * B**2 * C * D
    ),
}
for j, expected in key_rows.items():
    require(sp.expand(r[j] - expected) == 0, f"R coefficient r{j}")

# A nonzero binary cubic is a cube of one linear form exactly when these three
# catalecticant rows vanish.  Whole zero fibres are removed from tau.
T1 = sp.expand(g2**2 - 3 * g3 * g1)
T2 = sp.expand(g1**2 - 3 * g2 * g0)
T3 = sp.expand(g2 * g1 - 9 * g3 * g0)

# Omitted curve a=b^2/12, c=4/(3b), after multiplying by 12b.
f_omit = sp.expand(12 * b * (c + phi).subs({a: b**2 / 12, c: 4 / (3 * b)}))
f_expected = A * b**5 + (B + 12 * D) * b**4 + (C + 12 * E) * b**3 + 12 * F * b**2 + 3 * C * b + 16
require(sp.expand(f_omit - f_expected) == 0, "omitted support polynomial")
require(coeff(f_omit, b, 0) == 16 and sp.Poly(f_omit, b).degree() <= 5, "rho(f)<=5")


# ---------------------------------------------------------------------------
# Case I: A*D != 0.  Here rho(g3)=2 and it suffices to exclude rho(R)=1.
# ---------------------------------------------------------------------------

d, e, fvar, u, v = sp.symbols("d e f u v")


def case_I_atlas() -> list[str]:
    rows: list[str] = []

    # The leading packet leaves only degrees 14, 12, 11.
    c_tangent = B**2 / (4 * A)
    r12_tangent = sp.factor(r[12].subs(C, c_tangent))
    r11_tangent = sp.factor(r[11].subs(C, c_tangent))
    expected_r12 = 3 * B * (A**2 * B - 8 * A**2 * F + 4 * A * B * E - 2 * B**2 * D)
    require(sp.expand(r12_tangent - expected_r12) == 0, "case I tangent r12")
    tangent_F = B * (A**2 + 4 * A * E - 2 * B * D) / (8 * A**2)
    require(sp.expand(r11_tangent.subs(F, tangent_F) - 128 * A**3) == 0, "case I r11 fallback")
    require(sp.expand(r11_tangent.subs(B, 0) - 128 * A**3) == 0, "case I B=0 r11 fallback")

    # If R has one root, R=196608*D*(1+p*s)^M.  Weighted homogeneity lets
    # us set p=1 and substitute A=6*M*d.

    # M=14.  The ratio of the top two coefficients gives B=294D; the
    # leading and s^12 rows then solve C and F.  The s^3 row is linear in E.
    # Its leading/constant coefficient gcd is one, and after solving it the
    # s^4/s^5 univariate gcd is one.
    M = 14
    subs14 = {A: 84 * d, D: d, B: u, C: v, E: e, F: fvar}
    target14 = sp.expand(196608 * d * (1 + s) ** M)
    difference14 = sp.expand(R.subs(subs14) - target14)
    eq14 = {j: sp.cancel(coeff(difference14, s, j)) for j in range(2, 15)}
    u14 = 294 * d
    v14 = sp.factor(sp.solve(eq14[14].subs(u, u14), v)[0])
    f14 = sp.factor(sp.solve(eq14[12].subs({u: u14, v: v14}), fvar)[0])
    low14 = {j: sp.cancel(eq14[j].subs({u: u14, v: v14, fvar: f14})) for j in (3, 4, 5)}
    linear14 = sp.Poly(numerator(low14[3]), e)
    e14_coefficient = sp.Poly(linear14.coeff_monomial(e), d)
    e14_constant = sp.Poly(linear14.coeff_monomial(1), d)
    require(sp.gcd(e14_coefficient, e14_constant).degree() == 0, "I M=14 E exceptional gcd")
    e14 = sp.solve(low14[3], e)[0]
    tail14 = [sp.Poly(numerator(low14[j].subs(e, e14)), d) for j in (4, 5)]
    require(sp.gcd(tail14[0], tail14[1]).degree() == 0, "I M=14 tail gcd")
    rows.append("I/M=14: two univariate gcds are units")

    # M=12.  The top packet is tangent, C=B^2/(4A).  The s^3 row solves F.
    # Put B=D*w.  Away from w=0,-54, the s^12 row solves E and resultants of
    # s^4 against s^5,s^6,s^7 have gcd supported only at w(w+54)=0.
    M = 12
    w = sp.symbols("w")
    A12 = 72 * d
    C12 = u**2 / (4 * A12)
    target12 = sp.expand(196608 * d * (1 + s) ** M)
    difference12 = sp.expand(
        R.subs({A: A12, D: d, B: u, C: C12, E: e, F: fvar}) - target12
    )
    eq12 = {j: sp.cancel(coeff(difference12, s, j)) for j in range(3, 13)}
    f12 = sp.solve(eq12[3], fvar)[0]
    P12 = {j: sp.cancel(eq12[j].subs(fvar, f12).subs(u, d * w)) for j in range(4, 13)}
    require(
        sp.cancel(P12[12].subs(w, 0) / d).is_Rational and P12[12].subs(w, 0) != 0,
        "I M=12 w=0 exception",
    )
    e12 = sp.solve(P12[12], e)[0]
    reduced12 = {j: numerator(P12[j].subs(e, e12)) for j in (4, 5, 6, 7)}
    res12 = [
        sp.Poly(sp.resultant(reduced12[4], reduced12[j], d), w)
        for j in (5, 6, 7)
    ]
    gcd12 = reduce(sp.gcd, res12)
    require(roots_are_allowed(gcd12, w * (w + 54), w), "I M=12 resultant roots")
    exceptional12 = [numerator(P12[j].subs(w, -54)) for j in (4, 5, 6, 7, 12)]
    require(is_unit_ideal(exceptional12, (d, e), d), "I M=12 w=-54 exception")
    rows.append("I/M=12: resultant support {0,-54}; both exceptions empty")

    # M=11.  Again C=B^2/(4A), and s^3 solves F.  Split B=0 from B=D*w.
    # In the latter branch s^12 solves E away from w=-54; three resultants
    # have gcd supported at w=0.  At w=-54, s^12 forces d=880/891 and the
    # remaining coefficient gcd in E is one.
    M = 11
    A11 = 66 * d
    C11 = u**2 / (4 * A11)
    target11 = sp.expand(196608 * d * (1 + s) ** M)
    difference11 = sp.expand(
        R.subs({A: A11, D: d, B: u, C: C11, E: e, F: fvar}) - target11
    )
    eq11 = {j: sp.cancel(coeff(difference11, s, j)) for j in range(3, 13)}
    f11 = sp.solve(eq11[3], fvar)[0]
    b0_rows = [eq11[j].subs({fvar: f11, u: 0}) for j in (4, 5, 11)]
    require(is_unit_ideal(b0_rows, (d, e), d), "I M=11 B=0 branch")

    P11 = {j: sp.cancel(eq11[j].subs(fvar, f11).subs(u, d * w)) for j in range(4, 13)}
    e11 = sp.solve(P11[12], e)[0]
    reduced11 = {j: numerator(P11[j].subs(e, e11)) for j in (4, 5, 6, 11)}
    res11 = [
        sp.Poly(sp.resultant(reduced11[4], reduced11[j], d), w)
        for j in (5, 6, 11)
    ]
    gcd11 = reduce(sp.gcd, res11)
    require(roots_are_allowed(gcd11, w, w), "I M=11 resultant roots")

    exceptional_equation = sp.factor(P11[12].subs(w, -54))
    exceptional_d = [root for root in sp.solve(exceptional_equation, d) if root != 0]
    require(exceptional_d == [sp.Rational(880, 891)], "I M=11 w=-54 forced d")
    d11 = exceptional_d[0]
    exceptional_rows = [
        sp.Poly(numerator(P11[j].subs({w: -54, d: d11})), e)
        for j in (4, 5, 6, 11)
    ]
    require(reduce(sp.gcd, exceptional_rows).degree() == 0, "I M=11 w=-54 tail gcd")
    rows.append("I/M=11: B=0 unit; resultant support {0}; w=-54 tail gcd unit")

    return rows


def multiplicity_pairs(total: int) -> list[tuple[int, int]]:
    return [(m, total - m) for m in range(1, total // 2 + 1)]


# ---------------------------------------------------------------------------
# Case II: A=0, D!=0.  Here rho(g3)=1 and R'(0)=0, so a one-root R is
# impossible.  The following finite atlas excludes every two-root R.
# ---------------------------------------------------------------------------


def case_II_problem(m: int, n: int) -> tuple[dict[int, sp.Expr], sp.Expr, sp.Expr]:
    subs = {A: 0, D: d, B: u, C: v, E: e, F: fvar}
    target = binomial_target(196608 * d, ((1, m), (-sp.Rational(m, n), n)))
    difference = sp.expand(R.subs(subs) - target)
    equations = {j: sp.cancel(coeff(difference, s, j)) for j in range(2, 15)}
    f_solution = sp.factor(sp.solve(equations[2], fvar)[0])
    v_solution = sp.factor(sp.solve(equations[3].subs(fvar, f_solution), v)[0])
    return equations, f_solution, v_solution


def case_II_atlas() -> list[str]:
    rows: list[str] = []

    # Exact high-degree stratification.  Together with (B,C)!=(0,0), these
    # rows leave degrees 8,9,10,11,12 and no others.
    high_subs = {A: 0}
    require(sp.factor(r[12].subs(high_subs) - 4 * B**3 * (B + 12 * D)) == 0, "II r12")
    require(
        sp.factor(r[11].subs({A: 0, B: -12 * D}) + 6912 * D**3 * (C + 12 * E)) == 0,
        "II tangent r11",
    )
    require(
        sp.factor(r[10].subs({A: 0, B: -12 * D, C: -12 * E}) + 82944 * D**3 * F) == 0,
        "II tangent r10",
    )
    require(
        sp.factor(r[9].subs({A: 0, B: -12 * D, C: -12 * E, F: 0}) - 248832 * D**3 * E) == 0,
        "II tangent r9",
    )
    require(
        sp.factor(r[8].subs({A: 0, B: -12 * D, C: 0, E: 0, F: 0}) + 110592 * D**3) == 0,
        "II terminal r8",
    )
    require(
        sp.factor(r[10].subs({A: 0, B: 0}) + 3888 * C**2 * D**2) == 0,
        "II B=0 r10",
    )

    # E=0: all 24 multiplicity pairs are directly empty after the low rows
    # solve F and C.  This includes the degree-eight terminal chamber.
    e0_count = 0
    for total in range(8, 13):
        for m, n in multiplicity_pairs(total):
            equations, f_solution, v_solution = case_II_problem(m, n)
            specialized = {
                fvar: f_solution.subs(e, 0),
                v: v_solution.subs(e, 0),
                e: 0,
            }
            tails = [equations[j].subs(specialized) for j in range(4, 15)]
            require(is_unit_ideal(tails, (d, u), d), f"II E=0 pair ({m},{n})")
            e0_count += 1
    require(e0_count == 24, "II E=0 atlas size")
    rows.append("II/E=0: 24/24 unit")

    # E!=0, degree nine: B=-12D, C=-12E, F=0.
    count = 0
    for m, n in multiplicity_pairs(9):
        equations, _, _ = case_II_problem(m, n)
        tails = [equations[j].subs({u: -12 * d, v: -12 * e, fvar: 0}) for j in range(2, 15)]
        require(is_unit_ideal(tails, (d, e), d * e), f"II M=9 pair ({m},{n})")
        count += 1
    require(count == 4, "II M=9 atlas size")
    rows.append("II/E!=0,M=9: 4/4 unit")

    # E!=0, degree ten has exactly two high-packet branches.
    b0_count = 0
    tangent_count = 0
    for m, n in multiplicity_pairs(10):
        equations, f_solution, v_solution = case_II_problem(m, n)
        b0_subs = {u: 0, fvar: f_solution.subs(u, 0), v: v_solution.subs(u, 0)}
        b0_tails = [equations[j].subs(b0_subs) for j in range(4, 15)]
        require(is_unit_ideal(b0_tails, (d, e), d * e), f"II M=10 B=0 pair ({m},{n})")
        b0_count += 1

        f_tangent = f_solution.subs(u, -12 * d)
        tangent_subs = {u: -12 * d, v: -12 * e, fvar: f_tangent}
        tangent_tails = [equations[j].subs(tangent_subs) for j in range(3, 15)]
        require(
            is_unit_ideal(tangent_tails, (d, e), d * e),
            f"II M=10 tangent pair ({m},{n})",
        )
        tangent_count += 1
    require((b0_count, tangent_count) == (5, 5), "II M=10 atlas sizes")
    rows.append("II/E!=0,M=10: B=0 5/5 unit; tangent 5/5 unit")

    # E!=0, degree eleven: B=-12D and C+12E!=0.
    count = 0
    for m, n in multiplicity_pairs(11):
        equations, f_solution, v_solution = case_II_problem(m, n)
        f11 = f_solution.subs(u, -12 * d)
        v11 = v_solution.subs({u: -12 * d, fvar: f11})
        subs11 = {u: -12 * d, fvar: f11, v: v11}
        tails = [equations[j].subs(subs11) for j in range(4, 15)]
        require(
            is_unit_ideal(tails, (d, e), d * e * numerator(v11 + 12 * e)),
            f"II M=11 pair ({m},{n})",
        )
        count += 1
    require(count == 5, "II M=11 atlas size")
    rows.append("II/E!=0,M=11: 5/5 unit")

    # E!=0, degree twelve.  Put B=D*w.  For m<n the s^11 row solves E
    # away from q(w)=w^2-81w-972.  Three exact resultants have gcd w^15;
    # w=0 is forbidden by the degree-twelve packet.  The q=0 exceptional
    # ideal is unit.  The symmetric (6,6) row forces q=0 directly and its
    # residual four-row ideal is unit.
    w = sp.symbols("w")
    q_exception = w**2 - 81 * w - 972
    generic_count = 0
    for m, n in multiplicity_pairs(12):
        equations, f_solution, v_solution = case_II_problem(m, n)
        P = {
            j: sp.cancel(
                equations[j]
                .subs(fvar, f_solution)
                .subs(v, v_solution)
                .subs(u, d * w)
            )
            for j in range(4, 13)
        }
        if m < n:
            e_coefficient = sp.factor(sp.Poly(P[11], e).coeff_monomial(e))
            require(
                sp.factor(e_coefficient / (d**3 * w**2 * q_exception)).is_Rational,
                f"II M=12 ({m},{n}) s11 coefficient",
            )
            e_solution = sp.solve(P[11], e)[0]
            reduced = {
                j: numerator(P[j].subs(e, e_solution)) for j in (4, 5, 10, 12)
            }
            resultants = [
                sp.Poly(sp.resultant(reduced[12], reduced[j], d), w)
                for j in (4, 5, 10)
            ]
            resultant_gcd = reduce(sp.gcd, resultants)
            gcd_expr = sp.factor(resultant_gcd.as_expr())
            quotient = sp.cancel(gcd_expr / w**15)
            require(quotient.is_Rational and quotient != 0, f"II M=12 ({m},{n}) resultant gcd")

            exceptional = [numerator(P[j]) for j in (11, 12)] + [q_exception]
            require(
                is_unit_ideal(exceptional, (d, e, w), d * w),
                f"II M=12 ({m},{n}) q-exception",
            )
            generic_count += 1
        else:
            require(
                sp.factor(P[11] / (d**3 * e * w**2 * q_exception)).is_Rational,
                "II M=12 (6,6) forces q-exception",
            )
            exceptional = [numerator(P[j]) for j in (9, 10, 12)] + [q_exception]
            require(
                is_unit_ideal(exceptional, (d, e, w), d * e * w),
                "II M=12 (6,6) exceptional ideal",
            )
    require(generic_count == 5, "II M=12 generic atlas size")
    rows.append("II/E!=0,M=12: 5 resultant cells + symmetric exceptional unit")

    return rows


# ---------------------------------------------------------------------------
# Case III: A!=0, D=0.  Now R=s*S, rho(g3)=1, and it remains to prove
# rho(S)>=3.  The possible degrees of S are 13,11,10.
# ---------------------------------------------------------------------------


def case_III_atlas() -> list[str]:
    rows: list[str] = []
    S, remainder = sp.div(R.subs(D, 0), s)
    require(remainder == 0, "III R=s*S")
    S = sp.expand(S)
    require(coeff(S, s, 0) == 32768 * A, "III S(0)=32768A")

    # The same leading packet as Case I gives deg(S) in {13,11,10}.
    tangent = {D: 0, C: B**2 / (4 * A)}
    r12_tangent = sp.factor(r[12].subs(tangent))
    r11_tangent = sp.factor(r[11].subs(tangent))
    tangent_relation = A * B - 8 * A * F + 4 * B * E
    require(
        sp.expand(r12_tangent - 3 * A * B * tangent_relation) == 0,
        "III tangent r12",
    )
    require(
        sp.expand(r11_tangent.subs(F, B * (A + 4 * E) / (8 * A)) - 128 * A**3) == 0,
        "III tangent r11 fallback",
    )
    require(sp.expand(r11_tangent.subs(B, 0) - 128 * A**3) == 0, "III B=0 r11 fallback")

    # One-root atlas.  S=32768*A*(1+s)^M and the first coefficient gives
    # A=9E^2/(8M); hence E!=0.  The second coefficient solves B.
    for M in (13, 11, 10):
        A1 = sp.Rational(9, 8 * M) * e**2
        target = sp.expand(32768 * A1 * (1 + s) ** M)
        difference = sp.expand(S.subs({A: A1, B: u, C: v, E: e, F: fvar}) - target)
        equations = {j: sp.cancel(coeff(difference, s, j)) for j in range(1, 14)}
        u_low = sp.factor(sp.solve(equations[2], u)[0])

        if M == 13:
            u_top = sp.factor(M * A1 / 4)
            f_top = sp.factor(sp.solve(sp.Eq(u_low, u_top), fvar)[0])
            v_top = sp.factor(sp.solve(equations[13].subs({u: u_top, fvar: f_top}), v)[0])
            residuals = [
                sp.Poly(numerator(equations[j].subs({u: u_top, fvar: f_top, v: v_top})), e)
                for j in range(3, 12)
            ]
            residual_gcd = reduce(sp.gcd, residuals)
            require(residual_gcd.degree() == 0, "III one-root M=13 coefficient gcd")
        else:
            v_tangent = sp.cancel(u_low**2 / (4 * A1))
            if M == 10:
                indices = (3, 10, 11)  # R rows 4,11,12.
            else:
                indices = (3, 4, 11)  # R rows 4,5,12.
            residuals = [
                equations[j].subs(u, u_low).subs(v, v_tangent) for j in indices
            ]
            require(
                is_unit_ideal(residuals, (e, fvar), e),
                f"III one-root M={M} tangent ideal",
            )
        rows.append(f"III/one-root M={M}: empty")
        progress(f"III one-root M={M}")

    rr, scale, V = sp.symbols("r scale V")

    def product_data(m: int, n: int) -> tuple[sp.Poly, callable]:
        product = sp.Poly(sp.expand((1 + s) ** m * (1 + rr * s) ** n), s)
        return product, lambda j: product.coeff_monomial(s**j)

    # Two-root atlas with E=0.  Then h1=m+n*r=0, so r=-m/n.  The second
    # coefficient solves F, and all sixteen exact residual ideals are unit.
    k = sp.symbols("k")
    e0_count = 0
    for total in (10, 11, 13):
        for m, n in multiplicity_pairs(total):
            r_value = -sp.Rational(m, n)
            product = sp.Poly(sp.expand((1 + s) ** m * (1 + r_value * s) ** n), s)
            h = lambda j: product.coeff_monomial(s**j)
            subs = {A: k, B: u, C: v, D: 0, E: 0, F: fvar}
            equations = {
                j: sp.cancel(r[j + 1].subs(subs) - 32768 * k * h(j))
                for j in range(1, 14)
            }
            f_zero = sp.solve(equations[2], fvar)[0]
            tails = [equations[j].subs(fvar, f_zero) for j in range(3, 14)]
            require(
                is_unit_ideal(tails, (k, u, v), k),
                f"III E=0 pair ({m},{n})",
            )
            e0_count += 1
    require(e0_count == 16, "III E=0 atlas size")
    rows.append("III/two-root E=0: 16/16 unit")
    progress("III E=0 atlas")

    # Two roots with E!=0.  Put h1=m+n*r and normalize
    # E=(8h1/9)*scale, A=(8h1/9)*scale^2.  Thus scale,r,h1 are nonzero.

    # Degree thirteen: high ratio and leading coefficient solve B,F,C.
    # Low-order resultants can vanish only at r=0, h1=0, or m*r+n=0;
    # the last apparent branch has coefficient gcd one.
    m13_count = 0
    for m, n in multiplicity_pairs(13):
        product, h = product_data(m, n)
        h1 = h(1)
        h2 = h(2)
        A13 = sp.Rational(8, 9) * h1 * scale**2
        E13 = sp.Rational(8, 9) * h1 * scale
        B13 = sp.cancel(A13 * (m * rr + n) / (4 * rr))
        F13 = sp.cancel((B13 / scale - sp.Rational(8, 3) * h2) / 6)
        C13 = sp.cancel((B13**2 - 32768 * rr**n / A13) / (4 * A13))
        subs13 = {A: A13, B: B13, C: C13, D: 0, E: E13, F: F13}
        equations = {
            j: sp.cancel(r[j + 1].subs(subs13) - 32768 * A13 * h(j))
            for j in range(3, 12)
        }
        reduced = {
            j: strip_exact_factors(
                equations[j],
                (rr, scale),
                (scale, rr, h1, m * rr + n),
            )
            for j in range(3, 8)
        }
        resultants = [
            sp.Poly(sp.resultant(reduced[3], reduced[j], scale), rr)
            for j in (4, 5, 6, 7)
        ]
        resultant_gcd = reduce(sp.gcd, resultants)
        require(
            roots_are_allowed(resultant_gcd, rr * h1 * (m * rr + n), rr),
            f"III M=13 ({m},{n}) resultant roots",
        )

        apparent = -sp.Rational(n, m)
        apparent_rows = [
            sp.Poly(numerator(equations[j].subs(rr, apparent)), scale)
            for j in range(3, 12)
        ]
        apparent_gcd = reduce(sp.gcd, apparent_rows)
        require(apparent_gcd.degree() == 0, f"III M=13 ({m},{n}) apparent branch")
        m13_count += 1
        progress(f"III M=13 ({m},{n})")
    require(m13_count == 6, "III M=13 atlas size")
    rows.append("III/E!=0,M=13: 6 resultant cells; 6 apparent branches empty")

    # Degree ten, B=0 branch.
    m10_b0 = 0
    for m, n in multiplicity_pairs(10):
        product, h = product_data(m, n)
        h1 = h(1)
        h2 = h(2)
        kappa = sp.Rational(8, 9) * h1
        A10 = kappa * scale**2
        E10 = kappa * scale
        F10 = -sp.Rational(4, 9) * h2
        subs_b0 = {A: A10, B: 0, C: 0, D: 0, E: E10, F: F10}
        tails = [
            sp.cancel(r[j + 1].subs(subs_b0) - 32768 * A10 * h(j))
            for j in range(3, 12)
        ]
        require(
            is_unit_ideal(tails, (scale, rr), scale * rr * (rr - 1) * h1),
            f"III M=10 B=0 pair ({m},{n})",
        )
        m10_b0 += 1
        progress(f"III M=10 B=0 ({m},{n})")
    require(m10_b0 == 5, "III M=10 B=0 atlas size")

    # Degree ten, generic tangent branch.  The leading row against four
    # other rows has resultant gcd supported only at r=0 or h1=0.
    m10_generic = 0
    for m, n in multiplicity_pairs(10):
        product, h = product_data(m, n)
        h1 = h(1)
        h2 = h(2)
        kappa = sp.Rational(8, 9) * h1
        A10 = kappa * scale**2
        E10 = kappa * scale
        hbar = sp.Rational(8, 3) * h2
        F10 = sp.cancel(-hbar * (4 + scale) / (16 + 6 * scale))
        B10 = sp.cancel(-8 * scale * hbar / (16 + 6 * scale))
        C10 = sp.cancel(B10**2 / (4 * A10))
        subs10 = {A: A10, B: B10, C: C10, D: 0, E: E10, F: F10}
        equations = {
            j: sp.cancel(r[j + 1].subs(subs10) - 32768 * A10 * h(j))
            for j in range(3, 11)
        }
        reduced = {
            j: strip_exact_factors(
                equations[j] / scale**2,
                (rr, scale),
                (scale, rr, h1),
            )
            for j in (3, 4, 8, 9, 10)
        }
        resultants = [
            sp.Poly(sp.resultant(reduced[10], reduced[j], scale), rr)
            for j in (3, 4, 8, 9)
        ]
        resultant_gcd = reduce(sp.gcd, resultants)
        require(
            roots_are_allowed(resultant_gcd, rr * h1, rr),
            f"III M=10 ({m},{n}) generic resultant roots",
        )

        # At the omitted denominator scale=-8/3, the tangent relation forces
        # h2=0.  The leading equation becomes 128*A^3=32768*A*r^n.
        A_exception = sp.factor(kappa * sp.Rational(64, 9))
        leading_exception = sp.factor(128 * A_exception**3 - 32768 * A_exception * rr**n)
        exceptional_gcd = sp.gcd(sp.Poly(h2, rr), sp.Poly(leading_exception, rr))
        require(exceptional_gcd.degree() == 0, f"III M=10 ({m},{n}) denominator exception")
        m10_generic += 1
        progress(f"III M=10 tangent ({m},{n})")
    require(m10_generic == 5, "III M=10 generic atlas size")
    rows.append("III/E!=0,M=10: B=0 5/5 unit; tangent 5/5 resultants")

    # Degree eleven tangent branch.  Use V=6F+(8/3)h2, B=scale*V,
    # C=V^2/(4kappa).  After removing scale^2, the s^3 row is linear in
    # scale with coefficient a nonzero rational multiple of -2368*h1^2.
    m11_count = 0
    for m, n in multiplicity_pairs(11):
        product, h = product_data(m, n)
        h1 = h(1)
        h2 = h(2)
        kappa = sp.Rational(8, 9) * h1
        subs11 = {
            A: kappa * scale**2,
            B: scale * V,
            C: V**2 / (4 * kappa),
            D: 0,
            E: kappa * scale,
            F: (V - sp.Rational(8, 3) * h2) / 6,
        }
        equations = {
            j: sp.cancel(
                (r[j + 1].subs(subs11) - 32768 * kappa * scale**2 * h(j))
                / scale**2
            )
            for j in range(3, 12)
        }
        scale_coefficient = sp.factor(sp.Poly(equations[3], scale).coeff_monomial(scale))
        require(
            sp.cancel(scale_coefficient / (-2368 * h1**2)).is_Rational,
            f"III M=11 ({m},{n}) linear scale coefficient",
        )
        scale_solution = sp.solve(equations[3], scale)[0]
        reduced = {
            j: numerator(equations[j].subs(scale, scale_solution))
            for j in (4, 5, 9, 10, 11)
        }
        resultants = [
            sp.Poly(sp.resultant(reduced[11], reduced[j], V), rr)
            for j in (4, 5, 9, 10)
        ]
        resultant_gcd = reduce(sp.gcd, resultants)
        require(
            roots_are_allowed(resultant_gcd, rr * h1, rr),
            f"III M=11 ({m},{n}) resultant roots",
        )
        m11_count += 1
        progress(f"III M=11 ({m},{n})")
    require(m11_count == 5, "III M=11 atlas size")
    rows.append("III/E!=0,M=11: 5/5 resultant cells empty")

    return rows


# ---------------------------------------------------------------------------
# Euler fibre table, exact hostiles, and the THM-3573 degree boundary.
# ---------------------------------------------------------------------------


def fibre_table_audit() -> None:
    representatives = (
        b * (b - 1) * (b - 2),
        b**2 * (b - 1),
        b**3,
        b**2 - 1,
        b**2,
        b,
        sp.Integer(1),
        sp.Integer(0),
    )
    for polynomial in representatives:
        P = sp.Poly(polynomial, b)
        coefficients = [P.coeff_monomial(b**j) for j in (3, 2, 1, 0)]
        aa, bb, cc, dd = coefficients
        discriminant = bb**2 * cc**2 - 4 * aa * cc**3 - 4 * bb**3 * dd - 27 * aa**2 * dd**2 + 18 * aa * bb * cc * dd
        cubes = (
            bb**2 - 3 * aa * cc == 0
            and cc**2 - 3 * bb * dd == 0
            and bb * cc - 9 * aa * dd == 0
            and polynomial != 0
        )
        if polynomial == 0:
            affine_euler = 1
        else:
            affine_euler = rho(polynomial, b)
        deficit = int(discriminant == 0) + int(aa == 0) + int(cubes)
        require(3 - affine_euler == deficit, f"binary-cubic fibre table {polynomial}")


def specialized_counts(values: tuple[sp.Expr, ...]) -> tuple[int, int, int, int, int, int]:
    substitution = dict(zip((A, B, C, D, E, F), values))
    gv = [sp.expand(g.subs(substitution)) for g in (g3, g2, g1, g0)]
    tv = [sp.expand(t.subs(substitution)) for t in (T1, T2, T3)]
    rho_delta = rho(Delta.subs(substitution), s)
    rho_g3 = rho(gv[0], s)
    cube_gcd = gcd_many(tv, s)
    zero_gcd = gcd_many(tv + gv, s)
    tau = rho(cube_gcd, s) - rho(zero_gcd, s)
    rho_f = rho(f_omit.subs(substitution), b)
    chi_D = 3 - rho_delta - rho_g3 - tau
    chi_X = -3 + 2 * (rho_delta + rho_g3 + tau) - rho_f
    return rho_delta, rho_g3, tau, rho_f, chi_D, chi_X


def hostile_audit() -> list[str]:
    rows: list[str] = []
    fibre_table_audit()

    hostile_1 = (0, -12, 0, 1, 0, 0)  # phi=b^3-12ab.
    delta_1 = sp.factor(Delta.subs(dict(zip((A, B, C, D, E, F), hostile_1))))
    require(
        sp.expand(delta_1 + 12288 * s**2 * (9 * s**8 + 132 * s**4 - 16)) == 0,
        "hostile b^3-12ab discriminant",
    )
    counts_1 = specialized_counts(hostile_1)
    require(counts_1 == (9, 1, 0, 0, -7, 17), "hostile b^3-12ab counts")
    rows.append("hostile b^3-12ab: (rhoDelta,rhoG3,tau,rhoF,chiD,chiX)=(9,1,0,0,-7,17)")

    hostile_2 = (1, 0, 0, 0, 0, 0)  # phi=ab^2.
    delta_2 = sp.factor(Delta.subs(dict(zip((A, B, C, D, E, F), hostile_2))))
    require(
        sp.expand(delta_2 - 128 * s**3 * (s**10 - 718 * s**5 + 256)) == 0,
        "hostile ab^2 discriminant",
    )
    counts_2 = specialized_counts(hostile_2)
    require(counts_2 == (11, 1, 0, 5, -9, 16), "hostile ab^2 counts")
    rows.append("hostile ab^2: (11,1,0,5,-9,16)")

    # The tau correction is load-bearing: at s=1 this strict fibre is b^3.
    hostile_3 = (
        sp.Rational(1, 2),
        sp.Rational(1, 4),
        -16,
        0,
        0,
        sp.Rational(193, 48),
    )
    substitution_3 = dict(zip((A, B, C, D, E, F), hostile_3))
    require(sp.expand(G.subs(substitution_3).subs(s, 1) - b**3) == 0, "triple fibre G(1,b)=b^3")
    counts_3 = specialized_counts(hostile_3)
    require(counts_3 == (12, 1, 1, 5, -11, 20), "triple-fibre hostile counts")
    rows.append("hostile triple fibre: G(1,b)=b^3, tau=1, chiX=20")

    # Scope boundaries: pure b^3 is the already closed deg_a=0 row, while
    # H=-b/2 in THM-3573 first produces reducibility in total degree four.
    pure_b3 = specialized_counts((0, 0, 0, 1, 0, 0))
    require(pure_b3 == (5, 1, 0, 4, -3, 5), "pure b^3 scope control")
    H_boundary = -b / 2
    phi_boundary = sp.expand(4 * H_boundary * (1 + b * H_boundary + 4 * a * H_boundary**2))
    require(phi_boundary == -2 * a * b**3 + b**3 - 2 * b, "THM-3573 quartic hostile")
    require(sp.Poly(phi_boundary, a, b).total_degree() == 4, "quartic hostile degree")
    rows.append("scope: pure b^3 has chiX=5; reducible H=-b/2 starts at total degree 4")

    return rows


if __name__ == "__main__":
    selected = os.environ.get("THM3582_SECTION", "all")
    print("THM-3582 mixed-cubic target-graph exact audit")
    print("universal: strict transform, discriminant Delta=s^2*R, omitted support verified")
    if selected in ("all", "I"):
        for line in case_I_atlas():
            print(line)
    if selected in ("all", "II"):
        for line in case_II_atlas():
            print(line)
    if selected in ("all", "III"):
        for line in case_III_atlas():
            print(line)
    if selected in ("all", "controls"):
        for line in hostile_audit():
            print(line)
    print("all active truth gates passed")
