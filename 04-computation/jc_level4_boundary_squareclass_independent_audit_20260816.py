#!/usr/bin/env python3
"""Independent audit of the old-L boundary and level-four square class.

The submitted companion extracts a pre-specified face directly.  This audit
instead reconstructs THM-3495's pinned ``J``, discovers the exposed face,
asks SymPy to factor that 16-term face, evaluates the finite branch by nested
Horner evaluation, and keeps a literal sign/exponent ledger for Delta_3 and
Delta_4.  It never constructs Norm(J) or assumes a newest-factor law.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import runpy
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
GLOBAL = ROOT / "04-computation/keller_level_three_global_norm_probe_20260816.py"

# J is already a proved and independently audited THM-3495 object.  Rebuild it
# from that theorem's raw H input, but perform every level-four operation below
# in a separate representation.
capture = io.StringIO()
with contextlib.redirect_stdout(capture):
    source = runpy.run_path(str(GLOBAL))
require(capture.getvalue().rstrip().endswith("all exact checks passed"), "J rebuild failed")

J = source["J"]
raw_terms = [
    (tuple(map(int, monomial)), int(coefficient))
    for monomial, coefficient in J.to_dict().items()
]
ledger = "\n".join(
    f"{monomial}:{coefficient}"
    for monomial, coefficient in sorted(raw_terms, reverse=True)
)
j_hash = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(
    j_hash == "9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe",
    "THM-3495 J ledger changed",
)

# Discover the exposed face for the valuation vector (-1,0,+1) of the inverse
# coordinates (x,y,z), i.e. maximize i-k.  Convert only that face to SymPy and
# factor it without supplying the claimed answer.
weight_histogram: dict[int, int] = defaultdict(int)
for (i, _j, k), _coefficient in raw_terms:
    weight_histogram[i - k] += 1
top_weight = max(weight_histogram)
face_terms = [item for item in raw_terms if item[0][0] - item[0][2] == top_weight]

x, y, z = sp.symbols("x y z")
face_expression = sum(
    coefficient * x**i * y**j * z**k
    for (i, j, k), coefficient in face_terms
)
face_unit, face_factors = sp.factor_list(face_expression)
require(top_weight == 43, "exposed J weight is not 43")
require(len(face_terms) == 16, "exposed J face does not have 16 terms")
require(
    face_factors == [(x, 43), (3 * x * z - 2 * y, 15)],
    f"unexpected independent face factorization: {face_factors}",
)
expected_unit = -(2**58) * (3**51) * (13**8) * (79**4) * (313**2)
require(face_unit == expected_unit, "exposed J face unit changed")

# On either divergent inverse branch u=1/w is a uniformizer after the tame
# quadratic extension.  The leading coordinate relation is
# (x,y,z)=(u^-1,d,-3*d*u)+higher terms, d=D/S a unit.  Direct substitution in
# the independently factored face leaves a nonzero coefficient.
d = sp.symbols("d", nonzero=True)
reduced = sp.expand(face_expression.subs({x: 1, y: d, z: -3 * d}))
expected_reduced = sp.expand(expected_unit * (-11 * d) ** 15)
require(reduced == expected_reduced and reduced != 0, "divergent face cancels")
reduced_integer = int(expected_reduced.subs(d, 1))

# Evaluate the finite sheet in a representation distinct from the submitted
# term-by-term Fraction sum: first Horner-evaluate in z, then y, then x.
def nested_horner_value(
    terms: list[tuple[tuple[int, int, int], int]],
    xv: Fraction,
    yv: Fraction,
    zv: Fraction,
) -> Fraction:
    by_i_j: dict[tuple[int, int], dict[int, int]] = defaultdict(dict)
    for (i, j, k), coefficient in terms:
        by_i_j[(i, j)][k] = coefficient

    by_i: dict[int, dict[int, Fraction]] = defaultdict(dict)
    for (i, j), z_coefficients in by_i_j.items():
        value = Fraction(0)
        for exponent in range(max(z_coefficients), -1, -1):
            value = value * zv + z_coefficients.get(exponent, 0)
        by_i[i][j] = value

    x_coefficients: dict[int, Fraction] = {}
    for i, y_coefficients in by_i.items():
        value = Fraction(0)
        for exponent in range(max(y_coefficients), -1, -1):
            value = value * yv + y_coefficients.get(exponent, 0)
        x_coefficients[i] = value

    value = Fraction(0)
    for exponent in range(max(x_coefficients), -1, -1):
        value = value * xv + x_coefficients.get(exponent, 0)
    return value


finite_point = (Fraction(2), Fraction(5, 6), Fraction(-7, 8))
finite_value = nested_horner_value(raw_terms, *finite_point)
require(finite_value != 0, "J vanishes on the finite inverse sheet")
finite_ledger = f"{finite_value.numerator}/{finite_value.denominator}"
finite_hash = hashlib.sha256(finite_ledger.encode("ascii")).hexdigest()

# Verify that the chosen point lies on L=0 and that all generic-DVR units used
# by the Puiseux expansion are indeed simultaneously nonzero there.
ta, tb, tc = Fraction(2, 27), Fraction(1), Fraction(1)
L0 = 27 * ta**2 * tc**2 - 18 * ta * tb * tc + 16 * ta + tb**3 * tc - tb**2
T0 = 4 - 3 * tb * tc
S0 = 27 * ta * tc**2 - 9 * tb * tc + 8
D0 = 18 * ta * tc - 3 * tb**2 * tc + 2 * tb
require((L0, tc, T0, S0, D0) == (0, 1, 1, 1, Fraction(1, 3)), "DVR control failed")

# A literal parity ledger catches both historical failure modes: dropping the
# nonmonic 2^35 normalization and dropping the odd-block outer sign.
# Entries are exponents modulo 2 for the generators -1, 2, L, J, G.
def parity(**entries: int) -> dict[str, int]:
    return {name: exponent & 1 for name, exponent in entries.items() if exponent & 1}


# -L*N(H) = -J/(2^35*L^6).
delta3_from_norm_h = parity(sign=1, two=-35, L=-6, J=1)
delta3_pinned = parity(sign=1, two=1, J=1)
require(delta3_from_norm_h == delta3_pinned, "Delta_3 normalization lost a 2 or sign")

# Norm has cubic degree: N(-2)=(-2)^3.  Multiplication by Delta_1=-L
# cancels the two minus signs and leaves 2^3*L*N(J).
delta4_before_g = parity(sign=3 + 1, two=3, L=1, NJ=1)
require(delta4_before_g == parity(two=1, L=1, NJ=1), "Delta_4 sign/2 ledger failed")

# G=L^43*N(J), so L*N(J)=G/L^42.  The old L exponent is even.
delta4_after_g = parity(two=3, L=1 - 43, G=1)
require(delta4_after_g == parity(two=1, G=1), "G substitution parity failed")

print("== independent level-four boundary and square-class audit ==")
print(f"J coefficient-ledger sha256={j_hash}")
print(f"discovered Newton support: max(i-k)={top_weight}; face terms={len(face_terms)}")
print("SymPy factorization: -2^58*3^51*13^8*79^4*313^2*x^43*(3*x*z-2*y)^15")
print(f"divergent leading coefficient at d=1: {reduced_integer} != 0")
print("two divergent branches: valuations -43/2 and -43/2; norm uses a product, so no cross-sheet cancellation")
print(f"finite Horner control: q={finite_point}; value sha256={finite_hash}; nonzero")
print("v_L(N(J))=-43; hence G=L^43*N(J) is polynomial and gcd(G,L)=1 in Q[a,b,c]")
print("Delta_3 ledger: [-L*J/(2^35*L^7)]=[-2J]")
print("Delta_4 ledger: [N(-2J)*(-L)]=[8*L*N(J)]=[8*G/L^42]=[2G]")
print("scope: no claim that G is primitive, squarefree, irreducible, or a new image component")
print("all independent exact checks passed")
