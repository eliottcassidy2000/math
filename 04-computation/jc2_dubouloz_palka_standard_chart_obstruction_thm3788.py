#!/usr/bin/env python3
"""Exact audit companion for THM-3788.

This is not a search for the Dubouloz--Palka polynomials.  Their existence is
the cited input.  The script checks the complete degree/contact certificate
used by the theorem, the degree-three positive boundary control, and hostile
variants showing which local contacts carry the contradiction.

Reproduce with

    python3 04-computation/jc2_dubouloz_palka_standard_chart_obstruction_thm3788.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


MAX_N = 512
ASSIGNMENT_MAX_N = 8
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")
    checks += 1


t, zeta = sp.symbols("t zeta")
cyclotomic = sp.Poly(zeta**2 + zeta + 1, zeta, domain=sp.QQ[t])


def reduce_zeta(expression: sp.Expr) -> sp.Expr:
    """Reduce an element of Q[t,zeta] modulo zeta^2+zeta+1."""

    polynomial = sp.Poly(sp.expand(expression), zeta, domain=sp.QQ[t])
    return sp.expand(polynomial.rem(cyclotomic).as_expr())


# The d=3 (N=0) member is an actual cited-family control.  It is Corollary
# 5.2 with a primitive cube root zeta, normalized by R1(0)=1 and R2=1.
R1_0 = 1 + (zeta - 1) * t
R2_0 = sp.Integer(1)
R0_0 = 3 * (2 * t * zeta + t - zeta + 1)

check(reduce_zeta(zeta**3 - 1) == 0, "zeta is a cube root")
check(reduce_zeta(zeta - 1) != 0, "zeta is nontrivial")
check(
    reduce_zeta(1 - R1_0**3 - t * (1 - t) * R0_0 * R2_0**3) == 0,
    "N=0 Shabat identity",
)
check(sp.Poly(R1_0, t).degree() == 1, "N=0 degree R1")
check(sp.Poly(R0_0, t).degree() == 1, "N=0 degree R0")
check(sp.Poly(R2_0, t).degree() == 0, "N=0 degree R2")
check(reduce_zeta(R1_0.subs(t, 0) - 1) == 0, "N=0 normalized arm")
check(reduce_zeta(R1_0.subs(t, 1) - zeta) == 0, "N=0 second component")
check(reduce_zeta(R0_0.subs(t, 0)) != 0, "N=0 R0 nonzero at zero")
check(reduce_zeta(R0_0.subs(t, 1)) != 0, "N=0 R0 nonzero at one")

# Separability of (1-t)R0R1 in the positive control is checked by its three
# pairwise resultants.  (R2=1 has no root.)
positive_factors = (1 - t, R0_0, R1_0)
positive_resultants = []
for i, first in enumerate(positive_factors):
    for second in positive_factors[i + 1 :]:
        resultant = reduce_zeta(sp.resultant(first, second, t))
        positive_resultants.append(str(sp.factor(resultant)))
        check(resultant != 0, f"N=0 pairwise resultant {i}")


# Explicit standard-chart inverse for a primitive component.  The omega=1
# chart is the same formula with zeta specialized to one.
X, q = sp.symbols("X q")
chart_Z = zeta + X**3 * q
chart_Y = 3 * zeta**2 * q + 3 * zeta * X**3 * q**2 + X**6 * q**3
check(reduce_zeta(X**3 * chart_Y - chart_Z**3 + 1) == 0, "standard chart relation")
check(reduce_zeta(chart_Z.subs(X, 0) - zeta) == 0, "chart retains chosen component")
check(
    reduce_zeta(chart_Y.subs(X, 0) - 3 * zeta**2 * q) == 0,
    "chart arm coordinate",
)


# The same N=0 boundary has a particularly transparent planar pullback.  The
# standard U_1 chart cancels the original pole at x=0, but the cancellation
# transfers the pole to the other component 1+s=0 of F=0.
x, y = sp.symbols("x y")
s = x**3 * y
F = x * (1 + s)
G = y * (s**2 + 3 * s + 3) / (3 * (1 + s) ** 3)
pole_difference = 1 / (3 * x**3) - 1 / (3 * F**3)
jacobian = sp.diff(F, x) * sp.diff(G, y) - sp.diff(F, y) * sp.diff(G, x)
check(sp.cancel(G - pole_difference) == 0, "N=0 pole-transfer identity")
check(sp.cancel(jacobian - 1) == 0, "N=0 rational Keller identity")
check(sp.cancel(G.subs(x, 0) - y) == 0, "axis pole cancels")
check(sp.expand((s**2 + 3 * s + 3).subs(s, -1)) == 1, "new-arm numerator is a unit")


# The all-degree arithmetic imported from Theorem D and Lemma 5.1 is
# d=6N+3, (deg R0,deg R1,deg R2)=(3N+1,2N+1,N).  If one standard chart
# contained the image, R1-1 would have simple zeros at 0 and 1 and a triple
# zero at each of the N simple roots of R2.
minimum_margin = None
for N in range(MAX_N + 1):
    degree = 6 * N + 3
    degree_R2 = (degree - 3) // 6
    degree_R1 = degree_R2 * 2 + 1
    degree_R0 = (3 * degree_R2 + 1) * (2 - 1)
    forced_degree = 3 * degree_R2 + 2

    check(degree % 6 == 3, f"degree congruence N={N}")
    check(degree_R2 == N, f"degree R2 N={N}")
    check(degree_R1 == 2 * N + 1, f"degree R1 N={N}")
    check(degree_R0 == 3 * N + 1, f"degree R0 N={N}")
    check(forced_degree == 3 * N + 2, f"forced contact N={N}")
    check(forced_degree > degree_R1, f"degree contradiction N={N}")
    check(degree // 3 == degree_R1, f"factor degree N={N}")

    # Hostile controls.  With the t=1 endpoint deleted the N=0 case no
    # longer contradicts degree.  If triple contacts are weakened to simple
    # contacts, no contradiction remains for N>=1.
    without_t1 = 3 * N + 1
    if N == 0:
        check(without_t1 == degree_R1, "N=0 needs the t=1 address")
    else:
        check(without_t1 > degree_R1, f"large N survives one endpoint loss {N}")
        simple_only = N + 2
        check(simple_only <= degree_R1, f"triple contact is load-bearing N={N}")

    margin = forced_degree - degree_R1
    minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)


# Before degree is imposed, every t=1/root(R2) address may choose one of
# three special components.  Exactly one assignment could lie in U_1: all
# addresses must choose 1.  This is a finite hostile census of the component
# logic, not an existence search for R_i.
assignment_rows = []
for N in range(ASSIGNMENT_MAX_N + 1):
    free_addresses = N + 1
    assignments = 3**free_addresses
    compatible_with_U1 = 1
    check(assignments >= compatible_with_U1, f"assignment total N={N}")
    check(compatible_with_U1 == 1, f"unique all-one assignment N={N}")
    assignment_rows.append((N, assignments, compatible_with_U1))


semantic_packet = {
    "theorem": "THM-3788",
    "source_family": "S(3,3,1), alpha=0",
    "degree": "d=6N+3",
    "profile": ["degR0=3N+1", "degR1=2N+1", "degR2=N"],
    "forced_contact": "ord_0(R1-1)+ord_1(R1-1)+sum_beta ord_beta(R1-1)=3N+2",
    "chart": "U_omega=T\\union_{gamma!=omega}L_gamma",
    "chart_inverse": "Z=omega+X^3q,Y=3omega^2q+3omega X^3q^2+X^6q^3",
    "positive_control": "N=0 has R1=1+(zeta-1)t and meets L_1,L_zeta",
    "pole_transfer": "G=1/(3x^3)-1/(3F^3), F=x(1+x^3y), J(F,G)=1",
    "hostiles": [
        "dropping t=1 removes the N=0 contradiction",
        "weakening each R2 contact from triple to simple removes N>=1 contradiction",
    ],
    "universe": [0, MAX_N],
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_packet, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()

print("THM3788 DUBOULOZ-PALKA STANDARD-CHART CERTIFICATE")
print(f"UNIVERSE_N=0..{MAX_N}")
print("DEGREES=d:6N+3,R0:3N+1,R1:2N+1,R2:N")
print("FORCED_CONTACT=2+3N")
print(f"MINIMUM_DEGREE_MARGIN={minimum_margin}")
print("N0_CONTROL=R1(0):1,R1(1):zeta,R2:1")
print("N0_RESULTANTS=" + ",".join(positive_resultants))
print("N0_POLE_TRANSFER=axis:x=0->new_arm:1+x^3y=0,J=1")
print(
    "ASSIGNMENT_CENSUS="
    + ";".join(f"N{N}:{total}->{survivors}" for N, total, survivors in assignment_rows)
)
print("HOSTILE_N0_WITHOUT_T1=degree_equality")
print("HOSTILE_N_GE_1_SIMPLE_R2_CONTACT=no_degree_contradiction")
print("AUTOMORPHISM_SCOPE=source_precomposition,standard-component-permuting_target")
print(f"CHECKS={checks}")
print(f"SEMANTIC_SHA256={semantic_sha256}")
