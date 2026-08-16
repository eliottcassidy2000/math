#!/usr/bin/env python3
"""Exact degree-nine primitivity probe for all coordinates of F composed F.

At the good target (a,b,c)=(1,1,1), build the outer cubic algebra in X and
the inner cubic algebra in xi.  Reconstruct the final source coordinates
(xi,eta,zeta), verify the sporadic Keller map equations inside the algebra,
and test that the power bases of xi, eta, and zeta all have rank nine.

One good specialization proves generic primitivity.  Together with THM-2582,
the basis-discriminant identity then gives the same square class [H] for all
three level-two coordinate eliminants.  This is a fixed-map statement only.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction

import sympy as sp


EXPECTED_CHARPOLY_HASHES = {
    "x": "af1933db4a9e9720c71ef56628ad9eb646662ba020929a726a77387fb8e04929",
    "y": "393eb575a2b74eb0662c74c6e8c390de656f7a3b84d7b0930fc3f4088e2a688a",
    "z": "6ac03f6767ef062f7b0ea76dd3f158fb4c225dd29608c74dfbf142bd908646c0",
}
EXPECTED_RATIO_HASH = "4ef93fd8a0384bf6b9bb44340194a7d8a159a8da4d13268173663f800b94c001"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


Q = Fraction
ZERO3 = (Q(0), Q(0), Q(0))
ONE3 = (Q(1), Q(0), Q(0))
X3 = (Q(0), Q(1), Q(0))


def kadd(left, right):
    return tuple(left[i] + right[i] for i in range(3))


def kneg(value):
    return tuple(-entry for entry in value)


def ksub(left, right):
    return kadd(left, kneg(right))


def kscale(value, scalar):
    scalar = Q(scalar)
    return tuple(scalar * entry for entry in value)


def kmul(left, right):
    # Q[X]/(25X^3+X-2), reduced from degree four downward.
    work = [Q(0)] * 5
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            work[i + j] += first * second
    for degree in range(4, 2, -1):
        coefficient = work[degree]
        work[degree] = 0
        work[degree - 2] -= coefficient / 25
        work[degree - 3] += 2 * coefficient / 25
    return tuple(work[:3])


def kpow(value, exponent):
    result = ONE3
    base = value
    while exponent:
        if exponent & 1:
            result = kmul(result, base)
        base = kmul(base, base)
        exponent //= 2
    return result


def to_sympy(value):
    return [sp.Rational(entry.numerator, entry.denominator) for entry in value]


def from_sympy(entries):
    return tuple(Q(int(entry.p), int(entry.q)) for entry in entries)


def kmultiplication_matrix(value):
    basis = (ONE3, X3, kmul(X3, X3))
    columns = [sp.Matrix(to_sympy(kmul(value, vector))) for vector in basis]
    return sp.Matrix.hstack(*columns)


def kinv(value):
    matrix = kmultiplication_matrix(value)
    require(matrix.det() != 0, "attempted to invert a zero outer-algebra element")
    return from_sympy(matrix.inv()[:, 0])


def kdiv(numerator, denominator):
    return kmul(numerator, kinv(denominator))


def kconst(value):
    return (Q(value), Q(0), Q(0))


# The first inverse stage over the fixed target (1,1,1).
outer_denominator = kadd(kadd(kscale(kpow(X3, 2), 11), X3), kconst(2))
outer_correction = kadd(kscale(kpow(X3, 2), 24), kscale(X3, 6))
middle_y = ksub(ONE3, kdiv(outer_correction, outer_denominator))
middle_z_numerator = ksub(
    ksub(kscale(X3, 2), kscale(kmul(kpow(X3, 2), middle_y), 3)),
    ONE3,
)
middle_z = kdiv(middle_z_numerator, kpow(X3, 3))

middle_a = X3
middle_b = middle_y
middle_c = middle_z

middle_L = kadd(
    kadd(
        kadd(
            kadd(
                kscale(kmul(kpow(middle_a, 2), kpow(middle_c, 2)), 27),
                kscale(kmul(kmul(middle_a, middle_b), middle_c), -18),
            ),
            kscale(middle_a, 16),
        ),
        kmul(kpow(middle_b, 3), middle_c),
    ),
    kscale(kpow(middle_b, 2), -1),
)
middle_T = ksub(kconst(4), kscale(kmul(middle_b, middle_c), 3))
inner_p = kdiv(middle_T, middle_L)
inner_q = kdiv(kscale(middle_c, -2), middle_L)

# Norm_{outer/Q}(L(q)) is the THM-2582 value H/(64L) at this target.
norm_middle_L = Q(int(kmultiplication_matrix(middle_L).det().p), int(kmultiplication_matrix(middle_L).det().q))
require(norm_middle_L == Q(951326441195, 1600), "specialized Norm(L(q)) changed")


def aadd(left, right):
    return tuple(kadd(left[i], right[i]) for i in range(3))


def aneg(value):
    return tuple(kneg(entry) for entry in value)


def asub(left, right):
    return aadd(left, aneg(right))


def ascale(value, scalar):
    return tuple(kscale(entry, scalar) for entry in value)


def amul(left, right):
    # K[xi]/(xi^3+inner_p*xi+inner_q).
    work = [ZERO3 for _ in range(5)]
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            work[i + j] = kadd(work[i + j], kmul(first, second))
    for degree in range(4, 2, -1):
        coefficient = work[degree]
        work[degree] = ZERO3
        work[degree - 2] = ksub(work[degree - 2], kmul(coefficient, inner_p))
        work[degree - 3] = ksub(work[degree - 3], kmul(coefficient, inner_q))
    return tuple(work[:3])


def apow(value, exponent):
    result = AONE
    base = value
    while exponent:
        if exponent & 1:
            result = amul(result, base)
        base = amul(base, base)
        exponent //= 2
    return result


def alift(value):
    return (value, ZERO3, ZERO3)


def aflatten(value):
    return [entry for xi_coefficient in value for entry in xi_coefficient]


def aunflatten(entries):
    fractions = [Q(int(entry.p), int(entry.q)) for entry in entries]
    return tuple(tuple(fractions[3 * j + i] for i in range(3)) for j in range(3))


AONE = alift(ONE3)
AXI = (ZERO3, ONE3, ZERO3)


def abasis():
    result = []
    for xi_degree in range(3):
        for x_degree in range(3):
            outer = [Q(0), Q(0), Q(0)]
            outer[x_degree] = Q(1)
            element = [ZERO3, ZERO3, ZERO3]
            element[xi_degree] = tuple(outer)
            result.append(tuple(element))
    return tuple(result)


ABASIS = abasis()


def amultiplication_matrix(value):
    columns = [
        sp.Matrix(
            [sp.Rational(entry.numerator, entry.denominator) for entry in aflatten(amul(value, vector))]
        )
        for vector in ABASIS
    ]
    return sp.Matrix.hstack(*columns)


def ainv(value):
    matrix = amultiplication_matrix(value)
    require(matrix.det() != 0, "attempted to invert a zero tower element")
    return aunflatten(matrix.inv()[:, 0])


def adiv(numerator, denominator):
    return amul(numerator, ainv(denominator))


# Reconstruct the final source y,z from the general inverse-section N,D.
A = alift(middle_a)
B = alift(middle_b)
C = alift(middle_c)
xi2 = apow(AXI, 2)
xi3 = apow(AXI, 3)

D = aadd(
    aadd(
        aadd(
            aadd(
                aadd(ascale(amul(amul(xi2, A), C), -3), amul(amul(xi2, apow(B, 2)), C)),
                aneg(amul(xi2, B)),
            ),
            ascale(amul(amul(AXI, B), C), 2),
        ),
        ascale(AXI, -2),
    ),
    C,
)
N = aadd(
    aadd(
        aadd(
            aadd(ascale(amul(amul(amul(xi3, A), B), C), -3), ascale(amul(xi3, A), 4)),
            ascale(amul(amul(xi2, A), C), -6),
        ),
        amul(amul(xi2, apow(B, 2)), C),
    ),
    aadd(
        aadd(aneg(amul(xi2, B)), ascale(amul(amul(AXI, B), C), 2)),
        aadd(ascale(AXI, -2), C),
    ),
)
s = aneg(adiv(N, D))
eta = adiv(s, AXI)
zeta = adiv(asub(asub(ascale(AXI, 2), ascale(amul(AXI, s), 3)), C), xi3)

# Verify all three defining map coordinates inside the nested algebra.
u = aadd(AONE, amul(AXI, eta))
four_plus_three_xy = aadd(ascale(AONE, 4), ascale(amul(AXI, eta), 3))
F1 = aadd(amul(apow(u, 3), zeta), amul(amul(apow(eta, 2), u), four_plus_three_xy))
F2 = aadd(
    aadd(eta, ascale(amul(amul(AXI, apow(u, 2)), zeta), 3)),
    ascale(amul(amul(AXI, apow(eta, 2)), four_plus_three_xy), 3),
)
F3 = asub(asub(ascale(AXI, 2), ascale(amul(xi2, eta), 3)), amul(xi3, zeta))
require(F1 == A and F2 == B and F3 == C, "inverse reconstruction failed the Keller equations")


def power_basis_matrix(value):
    columns = []
    power = AONE
    for _ in range(9):
        columns.append(
            sp.Matrix(
                [sp.Rational(entry.numerator, entry.denominator) for entry in aflatten(power)]
            )
        )
        power = amul(power, value)
    return sp.Matrix.hstack(*columns)


tvar = sp.symbols("t")
coordinate_rows = []
coordinate_polynomials = {}
for name, value in (("x", AXI), ("y", eta), ("z", zeta)):
    power_matrix = power_basis_matrix(value)
    basis_determinant = sp.factor(power_matrix.det())
    require(basis_determinant != 0, f"{name} failed to generate the nine-sheet algebra")
    characteristic = sp.Poly(amultiplication_matrix(value).charpoly(tvar).as_expr(), tvar, domain=sp.QQ)
    require(characteristic.degree() == 9, f"{name} characteristic degree changed")
    require(sp.gcd(characteristic, characteristic.diff()).degree() == 0, f"{name} specialized eliminant is inseparable")
    factor_degrees = tuple(sorted(factor.degree() for factor, exponent in sp.factor_list(characteristic)[1] for _ in range(exponent)))
    coefficient_ledger = ",".join(str(coefficient) for coefficient in characteristic.all_coeffs())
    coefficient_hash = hashlib.sha256(coefficient_ledger.encode("ascii")).hexdigest()
    require(factor_degrees == (3, 6), f"{name} specialized factor-degree ledger changed")
    require(coefficient_hash == EXPECTED_CHARPOLY_HASHES[name], f"{name} characteristic ledger changed")
    coordinate_rows.append((name, basis_determinant, factor_degrees, coefficient_hash))
    coordinate_polynomials[name] = characteristic

determinants = {name: determinant for name, determinant, _, _ in coordinate_rows}
ratio_y = sp.factor(determinants["y"] / determinants["x"])
ratio_z = sp.factor(determinants["z"] / determinants["x"])
discriminants = {
    name: sp.factor(sp.discriminant(polynomial.as_expr(), tvar))
    for name, polynomial in coordinate_polynomials.items()
}
require(
    sp.factor(discriminants["y"] / discriminants["x"] - ratio_y**2) == 0,
    "y/x discriminant ratio is not the basis-change square",
)
require(
    sp.factor(discriminants["z"] / discriminants["x"] - ratio_z**2) == 0,
    "z/x discriminant ratio is not the basis-change square",
)

# The middle-stage X coordinate is deliberately constant on each three-point
# inner block.  It is a degree-three hostile and must not pass the rank-nine
# primitive-element test.
require(power_basis_matrix(A).det() == 0, "intermediate degree-three hostile became primitive")
ratio_ledger = repr((ratio_y, ratio_z))
ratio_hash = hashlib.sha256(ratio_ledger.encode("ascii")).hexdigest()
require(ratio_hash == EXPECTED_RATIO_HASH, "basis-change ratio ledger changed")

print("== F^2 all-coordinate degree-nine primitivity gate ==")
print("target=(1,1,1); outer algebra=Q[X]/(25X^3+X-2); inner algebra rank=3")
print(f"Norm_outer(L(q))={norm_middle_L}=H/(64L): PASS")
print("inverse formulas satisfy all three sporadic Keller coordinates inside the rank-nine algebra: PASS")
for name, determinant, factor_degrees, coefficient_hash in coordinate_rows:
    numerator_digits = len(str(abs(int(determinant.p))))
    denominator_digits = len(str(int(determinant.q)))
    print(
        f"{name}: power-basis determinant nonzero; numerator/denominator digits="
        f"{numerator_digits}/{denominator_digits}; factor degrees={factor_degrees}; charpoly sha256={coefficient_hash}"
    )
print(f"basis-change ratio ledger sha256={ratio_hash}")
print("direct discriminant ratios equal the squared y/x and z/x basis-change determinants: PASS")
print("intermediate outer-root coordinate has singular rank-nine power basis: PASS_HOSTILE")
print("consequence: x,y,z are generically primitive for F^2, so all three discriminant square classes equal THM-2582 [H]")
print("scope: fixed sporadic map at level two; square class only, not multiplicities, higher levels, family classification, or JC(2)")
print("all exact checks passed")
