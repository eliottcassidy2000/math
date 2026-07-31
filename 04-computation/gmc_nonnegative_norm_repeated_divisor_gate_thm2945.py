#!/usr/bin/env python3
"""Exact companion for THM-2945.

The theorem is an abstract real-algebra gate:

* a nonzero polynomial which is nonnegative on a positive ray can only
  vanish in the ray interior at a repeated root;
* multiplying by a strictly positive clearing/flag factor preserves that
  interior zero set, although it can add irrelevant negative repeated roots;
* the endpoint is separate.

This companion checks an exact conjugate-pair complete-intersection model,
1,296 factor-structured ray controls, and the sharp endpoint, sign and
clearing-factor hostiles.  It intentionally does not claim a universal
factorial terminal-pole divisor formula.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import product
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


THM2843_SOURCE = Path(__file__).with_name(
    "gmc_four_slot_projective_resolvent_thm2843.py"
)
THM2843_SHA256 = (
    "4832a9e4cda1608473e3a4bcfcb880a7fa3f1c6db47b7e939b4b5183cb9549aa"
)
THM2843_BYTES = THM2843_SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
require(
    sha256(THM2843_BYTES).hexdigest() == THM2843_SHA256,
    "THM-2843 complete-intersection dependency hash changed",
)

THM2925_SOURCE = Path(__file__).with_name(
    "gmc_general_width_terminal_pole_cancellation_thm2925.py"
)
THM2925_SHA256 = (
    "83d70a95f0943992d0e4b7027eede431d4dc968b66655e37b43fd0acfc692e47"
)
THM2925_BYTES = THM2925_SOURCE.read_bytes().replace(
    b"\r\n", b"\n"
).replace(b"\r", b"\n")
require(
    sha256(THM2925_BYTES).hexdigest() == THM2925_SHA256,
    "THM-2925 clearing dependency hash changed",
)

t, u, v = sp.symbols("t u v", real=True)


def poly(expression) -> sp.Poly:
    return sp.Poly(sp.expand(expression), t, domain=sp.QQ)


def monic(polynomial: sp.Poly) -> sp.Poly:
    require(not polynomial.is_zero, "zero polynomial has no monic normalization")
    return polynomial.monic()


def derivative_gcd(polynomial: sp.Poly) -> sp.Poly:
    return monic(sp.gcd(polynomial, polynomial.diff()))


# A binary complete-intersection model.  The sextic roots are the three
# conjugate pairs +-i, +-2i, +-3i.  Homogenizing f by v^6 makes it an octic;
# c(1,0)=1, so the affine resultant loses no point at infinity.
c6 = (u**2 + 1) * (u**2 + 4) * (u**2 + 9)
f8_affine = t + u**2 - 2
model_resultant = poly(sp.resultant(c6, f8_affine, u))
expected_model_resultant = poly(
    (t - 3) ** 2 * (t - 6) ** 2 * (t - 11) ** 2
)
require(
    model_resultant == expected_model_resultant,
    "conjugate-pair resultant model changed",
)
model_gcd = derivative_gcd(model_resultant)
require(
    model_gcd == monic(poly((t - 3) * (t - 6) * (t - 11))),
    "conjugate-pair repeated divisor changed",
)

a, b = sp.symbols("a b", real=True)
pair_block = sp.Matrix([[a, -b], [b, a]])
require(
    sp.det(pair_block) == a**2 + b**2,
    "conjugate local norm stopped being a square modulus",
)
nonreduced_pair_block = sp.Matrix(
    [
        [0, 0, 0, -1],
        [1, 0, 0, 0],
        [0, 1, 0, -2],
        [0, 0, 1, 0],
    ]
)
require(
    sp.det(nonreduced_pair_block) == 1,
    "nonreduced conjugate-local norm sign changed",
)


# Exhaustive factor-structured controls.  Every polynomial is nonnegative on
# [0,infinity): positive-ray roots occur with even exponent; negative roots
# have arbitrary exponent; the two irreducible quadratics are positive there.
factor_t = poly(t)
positive_factors = (poly(t - 1), poly(t - 2))
negative_factors = (poly(t + 1), poly(t + 2))
quadratic_factors = (poly(t**2 + t + 1), poly(t**2 + 3 * t + 3))
all_factors = (
    (factor_t,)
    + positive_factors
    + negative_factors
    + quadratic_factors
)
for left_index, left in enumerate(all_factors):
    require(left.is_sqf, f"control factor {left_index} is not squarefree")
    for right in all_factors[left_index + 1 :]:
        require(
            sp.gcd(left, right).degree() == 0,
            "control factors stopped being pairwise coprime",
        )

records: list[str] = []
structured_controls = 0
simple_endpoint_controls = 0
negative_repeated_controls = 0
interior_zero_controls = 0
test_points = (
    sp.Rational(1, 2),
    sp.Rational(3, 2),
    sp.Rational(5, 2),
    sp.Rational(7, 2),
)
for (
    endpoint_exponent,
    positive_one_exponent,
    positive_two_exponent,
    negative_one_exponent,
    negative_two_exponent,
    quadratic_one_exponent,
    quadratic_two_exponent,
) in product(
    range(4),
    (0, 2, 4),
    (0, 2, 4),
    range(3),
    range(3),
    range(2),
    range(2),
):
    exponents = (
        endpoint_exponent,
        positive_one_exponent,
        positive_two_exponent,
        negative_one_exponent,
        negative_two_exponent,
        quadratic_one_exponent,
        quadratic_two_exponent,
    )
    expression = sp.Integer(1)
    expected_gcd_expression = sp.Integer(1)
    for factor, exponent in zip(all_factors, exponents):
        expression *= factor.as_expr() ** exponent
        expected_gcd_expression *= factor.as_expr() ** max(exponent - 1, 0)
    polynomial = poly(expression)
    common = derivative_gcd(polynomial)
    expected_common = monic(poly(expected_gcd_expression))
    require(common == expected_common, "structured derivative gcd changed")

    for root, exponent in (
        (sp.Integer(1), positive_one_exponent),
        (sp.Integer(2), positive_two_exponent),
    ):
        require(
            (polynomial.eval(root) == 0)
            == (common.eval(root) == 0)
            == (exponent > 0),
            "positive-ray zero/repeated-divisor equivalence failed",
        )
    for point in test_points:
        if point not in (1, 2):
            require(
                polynomial.eval(point) > 0,
                "structured ray control lost strict positivity off its roots",
            )
    require(
        (polynomial.eval(0) == 0) == (endpoint_exponent > 0),
        "endpoint zero typing changed",
    )
    require(
        (common.eval(0) == 0) == (endpoint_exponent >= 2),
        "endpoint multiplicity typing changed",
    )

    structured_controls += 1
    simple_endpoint_controls += int(endpoint_exponent == 1)
    negative_repeated_controls += int(
        negative_one_exponent >= 2 or negative_two_exponent >= 2
    )
    interior_zero_controls += int(
        positive_one_exponent > 0 or positive_two_exponent > 0
    )
    records.append(
        ":".join(
            (
                ",".join(map(str, exponents)),
                str(polynomial.degree()),
                str(common.degree()),
                sha256(
                    ",".join(map(str, polynomial.all_coeffs())).encode()
                ).hexdigest(),
                sha256(
                    ",".join(map(str, common.all_coeffs())).encode()
                ).hexdigest(),
            )
        )
    )

require(structured_controls == 1296, "structured control census changed")


# Sharp hostiles and clearing-factor boundary.
endpoint_hostile = poly(t)
require(
    endpoint_hostile.eval(0) == 0
    and derivative_gcd(endpoint_hostile).degree() == 0,
    "simple endpoint hostile changed",
)

sign_hostile = poly(t - 1)
require(
    sign_hostile.eval(1) == 0
    and derivative_gcd(sign_hostile).degree() == 0
    and sign_hostile.eval(0) < 0
    and sign_hostile.eval(2) > 0,
    "sign-indefinite interior hostile changed",
)

strict_positive_clearing = poly((t + 1) ** 2)
cleared_constant_norm = strict_positive_clearing
require(
    all(strict_positive_clearing.eval(point) > 0 for point in (0, 1, 2, 7))
    and derivative_gcd(cleared_constant_norm) == monic(poly(t + 1)),
    "strict positive clearing hostile changed",
)
require(
    derivative_gcd(cleared_constant_norm).count_roots(0, sp.oo) == 0,
    "negative clearing divisor entered the positive ray",
)

non_strict_clearing = poly((t - 3) ** 2)
require(
    non_strict_clearing.eval(3) == 0
    and derivative_gcd(non_strict_clearing).eval(3) == 0,
    "non-strict clearing hostile changed",
)

mixed_norm = poly((t - 2) ** 2 * (t**2 + 1))
mixed_clearing = poly((t + 1) ** 2 * (t**2 + 3 * t + 3))
mixed_numerator = mixed_norm * mixed_clearing
mixed_common = derivative_gcd(mixed_numerator)
require(
    mixed_numerator.eval(2) == 0
    and mixed_common.eval(2) == 0
    and mixed_common.eval(-1) == 0,
    "positive multiplier transport control changed",
)

control_digest = sha256(("\n".join(records) + "\n").encode()).hexdigest()

print("THM-2945 NONNEGATIVE NORM AND REPEATED-DIVISOR GATE")
print(f"thm2843_dependency_sha256={THM2843_SHA256}")
print(f"thm2925_dependency_sha256={THM2925_SHA256}")
print(
    "abstract_gate=nonzero_H>=0_on_[0,infinity);"
    "zeros_(0,infinity)=positive_roots_gcd(H,Hprime)"
)
print(
    "endpoint_gate=H(0)_separate;"
    "simple_endpoint_H=t_is_missed_by_derivative_gcd"
)
print(
    "sign_hostile=H=t-1_has_simple_positive_root_and_gcd_1;"
    "nonnegativity_is_load_bearing"
)
print(
    "clearing_hostile=N=1,D=(t+1)^2;"
    "gcd(DN,(DN)prime)=t+1_has_only_negative_root"
)
print(
    "strict_clearing_required=D=(t-3)^2_adds_false_positive_ray_divisor"
)
print(
    "binary_CI_model=c6=(u2+1)(u2+4)(u2+9),"
    "f8=v6*(u2+(t-2)v2)"
)
print(
    "model_resultant=(t-3)^2*(t-6)^2*(t-11)^2;"
    "model_gcd=(t-3)(t-6)(t-11)"
)
print("conjugate_local_norm=a^2+b^2;nonreduced_pair_control=PASS")
print(
    f"factor_structured_controls={structured_controls};"
    f"simple_endpoint_controls={simple_endpoint_controls};"
    f"negative_repeated_controls={negative_repeated_controls};"
    f"interior_zero_controls={interior_zero_controls}"
)
print(f"factor_structured_digest_sha256={control_digest}")
print(
    "factorial_scope=conditional_after_positive_clearing_and_positive_"
    "extraneous_multiplier;no_universal_terminal_pole_divisor_formula"
)
print("all_exact_checks=PASS")
