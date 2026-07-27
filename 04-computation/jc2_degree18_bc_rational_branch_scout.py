#!/usr/bin/env python3
"""Exact standard-library scout for the two rational B--C ratios.

The B--C weighted ratio is t=C^2/B^3.  After the constant weighted
normalization B=1 and the coordinate x=Cy, the spectral cubic is defined
over Q(t).  This script computes its cubic discriminant in u, its
squarefree decomposition in x, the three branches at infinity, and a
small modular irreducibility witness for each rational ratio in THM-2311.

No computer-algebra package is used.  Polynomials are coefficient lists
in ascending degree with Fraction coefficients.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import math


Poly = list[Fraction]


def trim(poly: Poly) -> Poly:
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def add(first: Poly, second: Poly) -> Poly:
    size = max(len(first), len(second))
    return trim(
        [
            (first[index] if index < len(first) else Fraction(0))
            + (second[index] if index < len(second) else Fraction(0))
            for index in range(size)
        ]
    )


def neg(poly: Poly) -> Poly:
    return [-coefficient for coefficient in poly]


def sub(first: Poly, second: Poly) -> Poly:
    return add(first, neg(second))


def scale(poly: Poly, scalar: Fraction) -> Poly:
    return trim([scalar * coefficient for coefficient in poly])


def mul(first: Poly, second: Poly) -> Poly:
    result = [Fraction(0)] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            result[i + j] += left * right
    return trim(result)


def power(poly: Poly, exponent: int) -> Poly:
    result = [Fraction(1)]
    base = list(poly)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = mul(result, base)
        base = mul(base, base)
        remaining >>= 1
    return result


def derivative(poly: Poly) -> Poly:
    if len(poly) == 1:
        return [Fraction(0)]
    return trim(
        [Fraction(index) * coefficient for index, coefficient in enumerate(poly)][1:]
    )


def divide(first: Poly, second: Poly) -> tuple[Poly, Poly]:
    numerator = trim(first)
    denominator = trim(second)
    if denominator == [0]:
        raise ZeroDivisionError
    if len(numerator) < len(denominator):
        return [Fraction(0)], numerator
    quotient = [Fraction(0)] * (len(numerator) - len(denominator) + 1)
    remainder = list(numerator)
    while remainder != [0] and len(remainder) >= len(denominator):
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] / denominator[-1]
        quotient[shift] = coefficient
        remainder = sub(
            remainder,
            [Fraction(0)] * shift + scale(denominator, coefficient),
        )
    return trim(quotient), trim(remainder)


def exact_quotient(first: Poly, second: Poly) -> Poly:
    quotient, remainder = divide(first, second)
    require(remainder == [0], "polynomial division was not exact")
    return quotient


def monic(poly: Poly) -> Poly:
    result = trim(poly)
    return scale(result, Fraction(1, 1) / result[-1])


def gcd(first: Poly, second: Poly) -> Poly:
    left, right = trim(first), trim(second)
    while right != [0]:
        _, remainder = divide(left, right)
        left, right = right, remainder
    return monic(left)


def discriminant_cubic(a: Poly, b: Poly, c: Poly, d: Poly) -> Poly:
    # b^2 c^2 - 4ac^3 - 4b^3d - 27a^2d^2 + 18abcd.
    return add(
        sub(
            sub(
                mul(power(b, 2), power(c, 2)),
                scale(mul(a, power(c, 3)), Fraction(4)),
            ),
            scale(mul(power(b, 3), d), Fraction(4)),
        ),
        add(
            scale(mul(power(a, 2), power(d, 2)), Fraction(-27)),
            scale(mul(mul(mul(a, b), c), d), Fraction(18)),
        ),
    )


def squarefree_layers(poly: Poly) -> list[tuple[int, Poly]]:
    """Return (multiplicity, monic factor) in characteristic zero."""

    f = monic(poly)
    repeated = gcd(f, derivative(f))
    layer = exact_quotient(f, repeated)
    multiplicity = 1
    answer: list[tuple[int, Poly]] = []
    while layer != [Fraction(1)]:
        overlap = gcd(layer, repeated)
        exact = exact_quotient(layer, overlap)
        if exact != [Fraction(1)]:
            answer.append((multiplicity, monic(exact)))
        layer = overlap
        repeated = exact_quotient(repeated, overlap)
        multiplicity += 1
    require(repeated == [Fraction(1)], "squarefree decomposition left a tail")
    return answer


def primitive_coefficients(poly: Poly) -> tuple[int, ...]:
    denominator = 1
    for coefficient in poly:
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in poly]
    content = 0
    for coefficient in integers:
        content = math.gcd(content, abs(coefficient))
    integers = [coefficient // content for coefficient in integers]
    if integers[-1] < 0:
        integers = [-coefficient for coefficient in integers]
    return tuple(integers)


def coefficient_hash(poly: Poly) -> str:
    payload = ",".join(str(value) for value in primitive_coefficients(poly))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def evaluate(poly: Poly, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def evaluate_u(poly: Poly, value: Fraction) -> Fraction:
    return evaluate(poly, value)


def cubic_discriminant_scalar(a: int, b: int, c: int, d: int) -> int:
    return b * b * c * c - 4 * a * c**3 - 4 * b**3 * d - 27 * a * a * d * d + 18 * a * b * c * d


def prime_numbers(limit: int) -> list[int]:
    primes: list[int] = []
    for value in range(2, limit + 1):
        if all(value % prime for prime in primes if prime * prime <= value):
            primes.append(value)
    return primes


def modular_cubic_witness(
    coefficients: tuple[Fraction, ...],
) -> tuple[int, tuple[int, int, int, int]]:
    """Find p where the specialized cubic has no F_p root."""

    for prime in prime_numbers(251):
        if any(coefficient.denominator % prime == 0 for coefficient in coefficients):
            continue
        residues = [
            (coefficient.numerator * pow(coefficient.denominator, -1, prime)) % prime
            for coefficient in coefficients
        ]
        if residues[0] == 0:
            continue
        if all(
            sum(residues[index] * value ** (3 - index) for index in range(4))
            % prime
            for value in range(prime)
        ):
            return prime, tuple(residues)
    raise RuntimeError("no modular irreducibility witness found")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def analyze(
    t: Fraction,
) -> tuple[
    list[tuple[int, Poly]],
    str,
    bool,
    Fraction,
    Fraction,
    int,
    int,
    tuple[int, int, int, int],
]:
    # H=a(x)u^3+b(x)u^2+c(x)u+d(x), after multiplying by t^3.
    a = [Fraction(-26040609) * t**3]
    b = [Fraction(49601160) * t**3, 0, Fraction(1607445) * t**2]
    c = [
        Fraction(-20995200) * t**3,
        0,
        Fraction(-2857680) * t**2,
        0,
        Fraction(-138915) * t,
    ]
    d = [
        0,
        Fraction(-5598720) * t**3,
        Fraction(777600) * t**2,
        Fraction(-435456) * t**2,
        Fraction(78120) * t,
        0,
        Fraction(1127),
    ]
    branch = discriminant_cubic(a, b, c, d)
    layers = squarefree_layers(branch)

    repeated_linear = [
        factor
        for multiplicity, factor in layers
        if multiplicity == 2 and len(factor) == 2
    ]
    require(
        len(repeated_linear) == 1,
        "expected one double discriminant value",
    )
    exceptional_x = -repeated_linear[0][0] / repeated_linear[0][1]
    fiber = [
        evaluate(d, exceptional_x),
        evaluate(c, exceptional_x),
        evaluate(b, exceptional_x),
        evaluate(a, exceptional_x),
    ]
    fiber_gcd = gcd(fiber, derivative(fiber))
    if len(fiber_gcd) == 2:
        repeated_root_factor = fiber_gcd
    else:
        repeated_root_factor = gcd(fiber_gcd, derivative(fiber_gcd))
    require(
        len(repeated_root_factor) == 2,
        "exceptional fibre repeated root was not rational",
    )
    exceptional_u = (
        -repeated_root_factor[0] / repeated_root_factor[1]
    )

    hx = [
        evaluate(derivative(d), exceptional_x),
        evaluate(derivative(c), exceptional_x),
        evaluate(derivative(b), exceptional_x),
        evaluate(derivative(a), exceptional_x),
    ]
    hx_value = evaluate_u(hx, exceptional_u)
    if hx_value != 0 and len(fiber_gcd) == 3:
        exceptional_type = "smooth_e3"
        exceptional_unramified = False
    elif hx_value == 0:
        huu = evaluate_u(derivative(derivative(fiber)), exceptional_u)
        hux = evaluate_u(derivative(hx), exceptional_u)
        hxx = evaluate_u(
            [
                evaluate(derivative(derivative(d)), exceptional_x),
                evaluate(derivative(derivative(c)), exceptional_x),
                evaluate(derivative(derivative(b)), exceptional_x),
                evaluate(derivative(derivative(a)), exceptional_x),
            ],
            exceptional_u,
        )
        tangent_determinant = huu * hxx - hux * hux
        require(
            tangent_determinant != 0,
            "exceptional singularity was not an ordinary node",
        )
        require(
            huu != 0,
            "an exceptional node branch was vertical over the x-line",
        )
        exceptional_type = "ordinary_node"
        exceptional_unramified = True
    else:
        raise RuntimeError("unclassified double discriminant value")

    # Leading equation in v=u/x^2 at x=infinity.
    infinity = (
        a[0],
        b[2],
        c[4],
        d[6],
    )
    infinity_discriminant = cubic_discriminant_scalar(
        *[
            coefficient.numerator
            for coefficient in infinity
        ]
    )
    # Clearing different denominators changes the scalar but not zero/nonzero.
    common_denominator = 1
    for coefficient in infinity:
        common_denominator *= coefficient.denominator
    cleared = tuple(
        int(coefficient * common_denominator) for coefficient in infinity
    )
    infinity_discriminant = cubic_discriminant_scalar(*cleared)
    require(infinity_discriminant != 0, "infinity cubic acquired a repeated root")

    # x=1 gives a cubic in u.  Irreducibility of this specialization,
    # together with constant nonzero leading u coefficient, proves
    # irreducibility of the bivariate cubic.
    specialized = tuple(
        evaluate(coefficient, Fraction(1)) for coefficient in (a, b, c, d)
    )
    prime, residues = modular_cubic_witness(specialized)
    return (
        layers,
        exceptional_type,
        exceptional_unramified,
        exceptional_x,
        exceptional_u,
        infinity_discriminant,
        prime,
        residues,
    )


def main() -> None:
    ratios = (Fraction(-2000, 15309), Fraction(-125, 1134))
    for ratio in ratios:
        (
            layers,
            exceptional_type,
            exceptional_unramified,
            exceptional_x,
            exceptional_u,
            infinity_discriminant,
            prime,
            residues,
        ) = analyze(ratio)
        print(f"ratio={ratio}")
        print(
            "branch_squarefree_layers="
            + ",".join(
                f"multiplicity{multiplicity}:degree{degree}"
                for multiplicity, factor in layers
                for degree in (len(factor) - 1,)
            )
        )
        for multiplicity, factor in layers:
            degree = len(factor) - 1
            if degree == 1:
                print(
                    f"multiplicity{multiplicity}_linear_coefficients="
                    + ",".join(
                        str(value) for value in primitive_coefficients(factor)
                    )
                )
            print(
                f"multiplicity{multiplicity}_factor_sha256="
                f"{coefficient_hash(factor)}"
            )
        print(f"exceptional_x={exceptional_x}")
        print(f"exceptional_u={exceptional_u}")
        print(f"exceptional_type={exceptional_type}")
        print(
            "exceptional_normalization_unramified="
            f"{int(exceptional_unramified)}"
        )
        print(f"infinity_discriminant_nonzero={int(infinity_discriminant != 0)}")
        print(f"x1_irreducible_mod_prime={prime}")
        print(
            "x1_mod_prime_coefficients="
            + ",".join(str(residue) for residue in residues)
        )
    print("status=BC_RATIONAL_BRANCH_SCOUT_EXACT")


if __name__ == "__main__":
    main()
