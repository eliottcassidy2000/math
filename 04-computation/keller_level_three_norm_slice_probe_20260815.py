#!/usr/bin/env python3
"""Finite-exact level-three norm probe for the fixed sporadic Keller map.

This does not construct the global next image divisor J and does not compute
the degree-27 discriminant.  Instead it evaluates the proved THM-2576
polynomial H in the generic cubic inverse algebra on three one-parameter
target slices and takes its 3-by-3 multiplication determinant.  Thus it
computes Norm(H) directly, with exact rational-function arithmetic.

The output tests the only new parity datum needed by THM-2582 at level three:
whether the old factor L occurs to odd order in the denominator of Norm(H).
All arithmetic is over QQ.  ``require`` remains active under ``python -O``.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


A, w = sp.symbols("A w")
a, b, c = sp.symbols("a b c")

ROOT = Path(__file__).resolve().parents[1]
H_ARTIFACT = ROOT / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
H_ARTIFACT_SHA256 = (
    "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
)
H_ARTIFACT_BYTES = H_ARTIFACT.read_bytes()
require(
    hashlib.sha256(H_ARTIFACT_BYTES).hexdigest() == H_ARTIFACT_SHA256,
    "transported H pickle artifact changed",
)
H = pickle.loads(H_ARTIFACT_BYTES)

H_poly = sp.Poly(H, a, b, c)
H_ledger = "\n".join(
    f"{monomial}:{coefficient}" for monomial, coefficient in H_poly.terms()
)
require(
    hashlib.sha256(H_ledger.encode("ascii")).hexdigest()
    == "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2",
    "transported H coefficient ledger changed",
)

EXPECTED_K_HASHES = {
    (1, 2): "de18f1f38b29b92cbbe46c913a1446fe017a4566615d60e37a96d82307ae84a7",
    (3, 1): "f55ce6cbdeadcb27f30a7aa0ca75094b5ded4a344f6d976cf31021b2eb6e6579",
    (1, 3): "c3c98f832b92316343ee24da5f3898c5e3465fc621edd3f063441686cc62febf",
}


def coefficient_hash(poly: sp.Poly) -> str:
    ledger = "\n".join(
        f"{monomial}:{coefficient}" for monomial, coefficient in poly.terms()
    )
    return hashlib.sha256(ledger.encode("ascii")).hexdigest()


def slice_norm(b_value: int, c_value: int, denominator_two_exponent: int) -> None:
    """Compute Norm(H) in QQ(A)[w]/(E) on one exact target slice."""

    b0 = sp.Integer(b_value)
    c0 = sp.Integer(c_value)
    coefficient_field = sp.QQ.frac_field(A)

    L = (
        27 * A**2 * c0**2
        - 18 * A * b0 * c0
        + 16 * A
        + b0**3 * c0
        - b0**2
    )
    T = 4 - 3 * b0 * c0
    E = sp.Poly(L * w**3 + T * w - 2 * c0, w, domain=coefficient_field).monic()

    def polynomial(expression) -> sp.Poly:
        return sp.Poly(expression, w, domain=coefficient_field).rem(E)

    def multiply(left: sp.Poly, right: sp.Poly) -> sp.Poly:
        return (left * right).rem(E)

    def divide(numerator: sp.Poly, denominator: sp.Poly) -> sp.Poly:
        return multiply(numerator, sp.invert(denominator, E))

    def power(base: sp.Poly, exponent: int) -> sp.Poly:
        value = polynomial(1)
        factor = base
        while exponent:
            if exponent & 1:
                value = multiply(value, factor)
            exponent //= 2
            if exponent:
                factor = multiply(factor, factor)
        return value

    # The inverse graph from the exact THM-2576/2582 chart.
    qx = polynomial(w)
    y_denominator = polynomial((12 * A - b0**2) * w**2 + b0 * w + 2)
    qy = polynomial(b0) - divide(
        polynomial(3 * A * w * ((9 * A * c0 - b0) * w + 2)),
        y_denominator,
    )
    qz = divide(
        polynomial(2 * w - c0) - multiply(polynomial(3 * w**2), qy),
        power(qx, 3),
    )

    # Hostile type check: this really is the inverse fibre graph on the slice.
    unit = polynomial(1) + multiply(qx, qy)
    four_plus = polynomial(4) + multiply(polynomial(3), multiply(qx, qy))
    image_a = multiply(power(unit, 3), qz) + multiply(
        multiply(power(qy, 2), unit), four_plus
    )
    image_b = qy + multiply(
        polynomial(3), multiply(multiply(qx, power(unit, 2)), qz)
    ) + multiply(
        polynomial(3), multiply(multiply(qx, power(qy, 2)), four_plus)
    )
    image_c = (
        multiply(polynomial(2), qx)
        - multiply(polynomial(3), multiply(power(qx, 2), qy))
        - multiply(power(qx, 3), qz)
    )
    require(image_a == polynomial(A), "inverse graph lost the first target row")
    require(image_b == polynomial(b0), "inverse graph lost the second target row")
    require(image_c == polynomial(c0), "inverse graph lost the third target row")

    qx_powers = [power(qx, index) for index in range(H_poly.degree(a) + 1)]
    qy_powers = [power(qy, index) for index in range(H_poly.degree(b) + 1)]
    qz_powers = [power(qz, index) for index in range(H_poly.degree(c) + 1)]

    H_at_q = polynomial(0)
    for (i, j, k), coefficient in H_poly.terms():
        term = multiply(multiply(qx_powers[i], qy_powers[j]), qz_powers[k])
        term = term.mul_ground(coefficient_field.convert(coefficient))
        H_at_q = (H_at_q + term).rem(E)

    def algebra_norm(element: sp.Poly):
        """Norm by the determinant of multiplication in the cubic basis."""

        columns = []
        for basis_element in (polynomial(1), polynomial(w), polynomial(w**2)):
            product = multiply(element, basis_element).as_dict()
            columns.append(
                [
                    product.get((degree,), coefficient_field.zero)
                    for degree in range(3)
                ]
            )
        multiplication_matrix = sp.Matrix(
            3, 3, lambda row, column: sp.sympify(columns[column][row])
        )
        return sp.cancel(multiplication_matrix.det())

    # Positive control: reconstruct THM-2582's proved level-two norm identity
    # in this independently implemented quotient-algebra route.
    L_at_q = (
        multiply(polynomial(27), multiply(power(qx, 2), power(qz, 2)))
        - multiply(polynomial(18), multiply(multiply(qx, qy), qz))
        + multiply(polynomial(16), qx)
        + multiply(power(qy, 3), qz)
        - power(qy, 2)
    )
    norm_L = algebra_norm(L_at_q)
    H_slice = sp.Poly(H.subs({a: A, b: b0, c: c0}), A)
    require(
        sp.cancel(norm_L - H_slice.as_expr() / (64 * L)) == 0,
        "sliced quotient algebra failed the proved Norm(L)=H/(64L) control",
    )

    # Norm in a cubic algebra is the determinant of multiplication by H(q).
    # A second route computes the same norm as a resultant against monic E.
    norm_H = algebra_norm(H_at_q)
    resultant_norm_H = sp.cancel(sp.resultant(E.as_expr(), H_at_q.as_expr(), w))
    require(norm_H == resultant_norm_H, "determinant and resultant norms disagree")
    numerator, denominator = norm_H.as_numer_denom()
    numerator_poly = sp.Poly(numerator, A)
    content, primitive_numerator = numerator_poly.primitive()
    primitive_numerator = sp.Poly(primitive_numerator, A)

    L_poly = sp.Poly(L, A)
    factor_content, factors = sp.factor_list(primitive_numerator)

    require(content == 1 and factor_content == 1, "new numerator lost primitivity")
    require(
        [(sp.degree(factor, A), exponent) for factor, exponent in factors]
        == [(86, 1)],
        "new numerator is not irreducible of degree 86",
    )
    require(
        sp.cancel(denominator / L**7) == 2**denominator_two_exponent,
        "Norm(H) denominator is not the expected constant times L^7",
    )
    require(sp.gcd(numerator_poly, L_poly).degree() == 0, "new numerator contains L")
    require(
        sp.gcd(numerator_poly, H_slice).degree() == 0,
        "new numerator contains the old image factor H",
    )
    require(
        sp.gcd(sp.Poly(denominator, A), H_slice).degree() == 0,
        "Norm(H) denominator contains H",
    )
    require(sp.discriminant(L, A) != 0, "slice makes L inseparable")

    # THM-2582 gives [Delta_3]=[-L Norm(H)].  Since the L exponent becomes
    # 1-7=-6 and the power of two is odd on every chosen slice, the resulting
    # square class is [-2 K], where K is this primitive degree-86 numerator.
    require(denominator_two_exponent % 2 == 1, "constant square class changed")
    numerator_hash = coefficient_hash(primitive_numerator)
    require(
        numerator_hash == EXPECTED_K_HASHES[(b_value, c_value)],
        "new numerator coefficient ledger changed",
    )

    print(f"slice (b,c)=({b_value},{c_value})")
    print(f"  L={sp.factor(L)}; disc_A(L)={sp.discriminant(L, A)}")
    print("  Norm(L)=H/(64L) rederived; determinant norm=resultant norm")
    print(
        "  Norm(H)=K/(2^%d L^7); deg(K)=86; K irreducible; gcd(K,LH)=1"
        % denominator_two_exponent
    )
    print(f"  K coefficient-ledger sha256={numerator_hash}")
    print("  [Delta_3]=[-L Norm(H)]=[-2 K]; old L exponent 1-7=-6 is even")


print("== finite-exact level-three Keller norm slices ==")
print(f"transported H artifact sha256={H_ARTIFACT_SHA256}")
print("transported H coefficient ledger: PASS")
slice_norm(1, 2, 21)
slice_norm(3, 1, 35)
slice_norm(1, 3, 35)
print("three independent exact slices: denominator L^7, no L/H numerator factor")
print("scope: finite slices only; no global J, degree-27 discriminant, or JC claim")
print("all exact checks passed")
