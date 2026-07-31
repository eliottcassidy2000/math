#!/usr/bin/env python3
"""Exact companion for the split even-Faber ``y=0`` closure.

On the target-translated degree-22 even-Faber chosen split sheet, put
``Z=q^4`` and specialize the two exact flux numerators at ``y=0``.  The
nonzero first-flux constant gives

    G1 = q*N1 + 7496192*lambda = 0,

which is linear in ``u``; the second equation ``N2=0`` is cubic in ``u``.
This script computes ``Res_u(G1,N2)`` exactly and independently reconstructs
it by substituting the unique linear root.  Its primitive factor has degree
23 in ``q``, parameter-independent nonzero leading coefficient, and constant
term a nonzero multiple of ``lambda^3``.

The constant-field descent and third-flux contradiction are mathematical
steps in THM-2724.  No floating arithmetic or finite-field inference is used.
"""

from __future__ import annotations

import hashlib

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


q, u = s.symbols("q u")
B, C, D, E, W, lam = s.symbols("B C D E W lambda")
DEN = s.Integer(7496192)
Z = q**4

# Exact THM-2411 fluxes after y=0 and Z=q^4.
N1 = 1331 * (616 * B - 1089 * u) * Z + 4 * (2342560 * C * u - 3748096 * E)
N2 = (
    15944049 * Z**2
    - 206145280 * C * Z
    + 1443016960 * B * u**2
    - 1978994688 * D * u
    - 1319329792 * W
    - 1190488992 * u**3
)
G1 = s.expand(q * N1 + DEN * lam)

G1_poly = s.Poly(G1, u)
N2_poly = s.Poly(N2, u)
require(G1_poly.degree() == 1, "first split flux is not linear in u")
require(N2_poly.degree() == 3, "second split flux is not cubic in u")
require(
    s.expand(
        G1_poly.coeff_monomial(u)
        - (-14641 * q * (-640 * C + 99 * q**4))
    )
    == 0,
    "linear u coefficient changed",
)
require(
    s.expand(
        G1_poly.coeff_monomial(1)
        - 117128 * (7 * B * q**5 - 128 * E * q + 64 * lam)
    )
    == 0,
    "constant first-flux coefficient changed",
)
require(N2_poly.LC() == -1190488992, "second-flux cubic pivot changed")
require(G1.subs(q, 0) == DEN * lam, "q=0 did not retain nonzero lambda")


# Direct resultant and an independent linear-root substitution.
resultant = s.expand(s.resultant(G1, N2, u))
a = G1_poly.coeff_monomial(u)
b0 = G1_poly.coeff_monomial(1)
substitution_numerator = s.cancel(a**3 * N2.subs(u, -b0 / a))
require(s.denom(substitution_numerator) == 1, "linear substitution retained a denominator")
require(
    s.expand(resultant + substitution_numerator) == 0,
    "resultant and linear-root reconstruction disagree",
)

content, primitive = s.primitive(resultant)
if s.Poly(primitive, q).LC() < 0:
    content = -content
    primitive = -primitive
primitive = s.expand(primitive)
require(content == 505447028499293771, "resultant content changed")


# Freeze the complete 12-coefficient sparse eliminant, grouped by q-degree.
expected = s.expand(
    96059601 * q**23
    - 3104956800 * C * q**19
    + (
        1483603968 * B**3
        - 6744342528 * B * D
        + 36130406400 * C**2
        - 7948689408 * W
    )
    * q**15
    + (
        -17983078400 * B**3 * C
        - 30519853056 * B**2 * E
        + 87199580160 * B * C * D
        - 181665792000 * C**3
        + 154156400640 * C * W
        + 123325120512 * D * E
    )
    * q**11
    + (15259926528 * B**2 - 61662560256 * D) * lam * q**10
    + (
        657666867200 * B**2 * C * E
        - 281857228800 * B * C**2 * D
        - 372051542016 * B * E**2
        + 335544320000 * C**4
        - 996566630400 * C**2 * W
        - 1594506608640 * C * D * E
    )
    * q**7
    + (
        -328833433600 * B**2 * C
        + 372051542016 * B * E
        + 797253304320 * C * D
    )
    * lam
    * q**6
    - 93012885504 * B * lam**2 * q**5
    + (
        -6012954214400 * B * C * E**2
        + 2147483648000 * C**3 * W
        + 5153960755200 * C**2 * D * E
        + 7937099563008 * E**3
    )
    * q**3
    + (
        6012954214400 * B * C * E
        - 2576980377600 * C**2 * D
        - 11905649344512 * E**2
    )
    * lam
    * q**2
    + (-1503238553600 * B * C + 5952824672256 * E) * lam**2 * q
    - 992137445376 * lam**3
)
require(primitive == expected, "primitive degree-23 eliminant changed")

q_polynomial = s.Poly(primitive, q)
all_polynomial = s.Poly(primitive, q, B, C, D, E, W, lam)
q_support = sorted(power[0] for power, _ in q_polynomial.terms())
require(q_polynomial.degree() == 23, "q eliminant degree changed")
require(
    q_support == [0, 1, 2, 3, 5, 6, 7, 10, 11, 15, 19, 23],
    "q eliminant support changed",
)
require(len(all_polynomial.terms()) == 34, "q eliminant monomial count changed")
require(q_polynomial.LC() == 96059601, "q eliminant leading pivot changed")
require(
    q_polynomial.TC() == -992137445376 * lam**3,
    "q eliminant constant pivot changed",
)

digest = hashlib.sha256(s.srepr(all_polynomial.as_expr()).encode()).hexdigest()
require(
    digest == "e50e733ec85df8c1aff2bf805a6442424b82b3c3dbb4cb5563cea363c9f7b6e0",
    "q eliminant canonical digest changed",
)


print("split even-Faber y-zero exact resultant")
print("flux_specialization=y=0:Z=q^4")
print("first_flux=G1=q*N1+7496192*lambda:degree_u1")
print("second_flux=N2:degree_u3:leading_u=-1190488992")
print("independent_paths=direct_resultant,linear_root_substitution:match=True")
print(f"resultant_content={content}")
print(f"primitive_q_degree={q_polynomial.degree()}:q_coefficients={len(q_polynomial.terms())}:monomials={len(all_polynomial.terms())}")
print("q_support=0,1,2,3,5,6,7,10,11,15,19,23")
print("leading_q23=96059601:parameter_independent_nonzero=True")
print("constant_q0=-992137445376*lambda^3:nonzero_when_lambda_nonzero=True")
print(f"primitive_digest={digest}")
print("constant_field_consequence=q_constant_then_u_constant")
print("scope=EXACT_Y_ZERO_EVEN_FABER_BOUNDARY_NOT_ODD_BANK_NOT_JC2")
print("ALL CHECKS PASSED")
