#!/usr/bin/env python3
"""Exact cyclotomic maximin certificate for THM-3667."""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENTS = (
    ROOT / "04-computation/lrc_minimal_three_twist_target_detector_thm3665.py",
    ROOT / "05-knowledge/results/lrc_minimal_three_twist_target_detector_thm3665.out",
)
EXPECTED_PARENT_HASHES = (
    "a5ef75f038d80c5d91308bb5379303970b44e9e538323eb49cf8779386356938",
    "172fe3e32fc27bb2abb21f4c7a7af59e71cdfa4c604586b2d4f8725a5ac6211a",
)
EXPECTED_SEMANTIC_SHA256 = "528591ae1b3448ba5a3e50b9c9e72a303a9d337ba4bbc46f8126006e319be12e"

P = 13
N = P * P
# Minimal polynomial of x=2*cos(pi/13), coefficients in ascending order.
MINPOLY = (-1, -3, 6, 4, -5, -1, 1)
X_INTERVAL = (Fraction(1_941_883, 1_000_000),
              Fraction(1_941_884, 1_000_000))


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def reduce_poly(poly: list[int]) -> tuple[int, ...]:
    values = trim(poly[:])
    while len(values) > 6:
        degree = len(values) - 1
        coefficient = values[-1]
        shift = degree - 6
        if coefficient:
            for index, minpoly_coefficient in enumerate(MINPOLY):
                values[shift + index] -= coefficient * minpoly_coefficient
        values = trim(values)
    values += [0] * (6 - len(values))
    return tuple(values)


def add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    size = max(len(left), len(right))
    return reduce_poly([
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ])


def neg(value: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-coefficient for coefficient in value)


def sub(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return add(left, neg(right))


def mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    product = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            product[i + j] += a * b
    return reduce_poly(product)


def constant(value: int) -> tuple[int, ...]:
    return (value, 0, 0, 0, 0, 0)


def interval_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def interval_mul(left, right):
    products = (
        left[0] * right[0], left[0] * right[1],
        left[1] * right[0], left[1] * right[1],
    )
    return min(products), max(products)


def interval_eval(poly: tuple[int, ...], interval):
    result = (Fraction(0), Fraction(0))
    for coefficient in reversed(poly):
        result = interval_add(interval_mul(result, interval),
                              (Fraction(coefficient), Fraction(coefficient)))
    return result


def rational_poly_eval(poly: tuple[int, ...], value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def canonical_cyclotomic(coefficients: list[int]) -> tuple[int, ...]:
    """Canonical basis 1,z,...,z^11 in Q[z]/Phi_13."""
    values = [0] * P
    for index, coefficient in enumerate(coefficients):
        values[index % P] += coefficient
    top = values[12]
    return tuple(values[index] - top for index in range(12))


def main() -> None:
    parent_hashes = tuple(lf_sha256(path) for path in PARENTS)
    require(parent_hashes == EXPECTED_PARENT_HASHES,
            ("parent hashes", parent_hashes))

    # Rigorous enclosure of x=2*cos(pi/13).  The classical rational pi
    # bounds and alternating cosine series place x inside X_INTERVAL.
    pi_low = Fraction(103_993, 33_102)
    pi_high = Fraction(104_348, 33_215)
    t_low, t_high = pi_low / P, pi_high / P
    cos_lower = sum(
        Fraction((-1) ** k, factorial(2 * k)) * t_high ** (2 * k)
        for k in range(6)
    )  # through the negative t^10 term
    cos_upper = sum(
        Fraction((-1) ** k, factorial(2 * k)) * t_low ** (2 * k)
        for k in range(5)
    )  # through the positive t^8 term
    require(2 * cos_lower > X_INTERVAL[0]
            and 2 * cos_upper < X_INTERVAL[1],
            ("cosine enclosure", 2 * cos_lower, 2 * cos_upper, X_INTERVAL))
    p_left = rational_poly_eval(MINPOLY, X_INTERVAL[0])
    p_right = rational_poly_eval(MINPOLY, X_INTERVAL[1])
    derivative = tuple(index * MINPOLY[index]
                       for index in range(1, len(MINPOLY)))
    derivative_interval = interval_eval(derivative, X_INTERVAL)
    require(p_left * p_right < 0 and derivative_interval[0] > 0,
            ("isolating interval", p_left, p_right, derivative_interval))

    one = constant(1)
    two = constant(2)
    xpoly = (0, 1, 0, 0, 0, 0)
    x2 = mul(xpoly, xpoly)
    ypoly = sub(x2, two)   # y=2*cos(2*pi/13)
    denominator = add(ypoly, one)
    # Optimal normalized positive weights are a=y/(y+1), b=1/(y+1).
    numerator_a = ypoly
    numerator_b = one
    require(interval_eval(numerator_a, X_INTERVAL)[0] > 0
            and interval_eval(denominator, X_INTERVAL)[0] > 0,
            "positive optimal weights")

    # C_n(x)=2*cos(n*pi/13), C_n=x*C_(n-1)-C_(n-2).
    chebyshev = [constant(2), xpoly]
    for _degree in range(2, 25):
        chebyshev.append(sub(mul(xpoly, chebyshev[-1]), chebyshev[-2]))
    sine4 = tuple(sub(two, chebyshev[2 * k]) for k in range(P))
    require(sine4[0] == constant(0), sine4[0])
    # sine4[k]=4*sin(k*pi/13)^2=2-C_(2k).

    gap_numerator = sine4[1]
    squared_rows = {}
    lower_differences = {}
    upper_differences = {}
    for u in range(P):
        for v in range(P):
            # D^2 |a+b*z^u-z^v|^2
            value = sub(add(
                mul(mul(numerator_a, denominator), sine4[v]),
                mul(denominator, sine4[(u - v) % P]),
            ), mul(numerator_a, sine4[u]))
            squared_rows[(u, v)] = value
            if (u, v) != (0, 0):
                lower_differences[(u, v)] = sub(value, gap_numerator)

    lower_equalities = tuple(
        frequency for frequency, value in lower_differences.items()
        if value == constant(0)
    )
    require(lower_equalities == ((1, 0), (2, 1), (11, 12), (12, 0)),
            ("lower equality set", lower_equalities))
    lower_intervals = {}
    for frequency, value in lower_differences.items():
        if value == constant(0):
            continue
        enclosure = interval_eval(value, X_INTERVAL)
        require(enclosure[0] > 0,
                ("lower spectral inequality", frequency, value, enclosure))
        lower_intervals[frequency] = enclosure

    upper_numerator = mul(mul(denominator, denominator), sine4[6])
    for frequency, value in squared_rows.items():
        upper_differences[frequency] = sub(upper_numerator, value)
    upper_equalities = tuple(
        frequency for frequency, value in upper_differences.items()
        if value == constant(0)
    )
    require(upper_equalities == ((0, 6), (0, 7)),
            ("upper equality set", upper_equalities))
    upper_intervals = {}
    for frequency, value in upper_differences.items():
        if value == constant(0):
            continue
        enclosure = interval_eval(value, X_INTERVAL)
        require(enclosure[0] > 0,
                ("upper spectral inequality", frequency, value, enclosure))
        upper_intervals[frequency] = enclosure

    # The swapped optimizer has the same multiplier magnitudes under
    # (u,v)->(-u,v-u), but equality of complex eigenvalues is not preserved
    # by the resulting frequency-dependent phase.  Record both equality
    # sets before doing the two separate cyclotomic collision censuses.
    swapped_squared_rows = {
        ((-u) % P, (v - u) % P): value
        for (u, v), value in squared_rows.items()
    }
    swapped_lower_equalities = tuple(
        frequency for frequency, value in swapped_squared_rows.items()
        if frequency != (0, 0) and sub(value, gap_numerator) == constant(0)
    )
    require(tuple(sorted(swapped_lower_equalities)) ==
            ((1, 1), (2, 1), (11, 12), (12, 12)),
            ("swapped lower equality set", swapped_lower_equalities))
    swapped_upper_equalities = tuple(
        frequency for frequency, value in swapped_squared_rows.items()
        if sub(upper_numerator, value) == constant(0)
    )
    require(tuple(sorted(swapped_upper_equalities)) == ((0, 6), (0, 7)),
            ("swapped upper equality set", swapped_upper_equalities))

    # Exact eigenvalue collision census at (a,b)=(y/D,1/D):
    # D*lambda(u,v)=y+z^u-(y+1)z^v.
    eigenclasses = {}
    for u in range(P):
        for v in range(P):
            coefficients = [0] * P
            coefficients[1] += 1
            coefficients[-1] += 1                         # y=z+z^-1
            coefficients[u] += 1
            coefficients[(v + 1) % P] -= 1               # -y*z^v
            coefficients[(v - 1) % P] -= 1
            coefficients[v] -= 1                         # -z^v
            key = canonical_cyclotomic(coefficients)
            eigenclasses.setdefault(key, []).append((u, v))
    multiplicities = tuple(sorted(len(value) for value in eigenclasses.values()))
    require(len(eigenclasses) == 156
            and multiplicities.count(2) == 13
            and multiplicities.count(1) == 143,
            ("collision multiplicities", len(eigenclasses), multiplicities))
    double_classes = tuple(sorted(
        tuple(sorted(value)) for value in eigenclasses.values() if len(value) == 2
    ))
    expected_doubles = tuple(sorted(
        tuple(sorted(((u, (u + 1) % P),
                      ((u + 3) % P, (u + 2) % P))))
        for u in range(P)
    ))
    require(double_classes == expected_doubles,
            ("double classes", double_classes, expected_doubles))
    zero_classes = tuple(value for key, value in eigenclasses.items()
                         if key == (0,) * 12)
    require(zero_classes == ([(0, 0)],), ("zero class", zero_classes))
    centralizer_dimension = 143 + 13 * 4
    require(centralizer_dimension == 195, centralizer_dimension)

    # At the equally optimal swapped orientation (a,b)=(1/D,y/D),
    # D*lambda(u,v)=1+y*z^u-(y+1)z^v.  Its spectrum is simple.
    swapped_eigenclasses = {}
    for u in range(P):
        for v in range(P):
            coefficients = [0] * P
            coefficients[0] += 1
            coefficients[(u + 1) % P] += 1             # +y*z^u
            coefficients[(u - 1) % P] += 1
            coefficients[(v + 1) % P] -= 1             # -(y+1)*z^v
            coefficients[(v - 1) % P] -= 1
            coefficients[v] -= 1
            key = canonical_cyclotomic(coefficients)
            swapped_eigenclasses.setdefault(key, []).append((u, v))
    swapped_multiplicities = tuple(
        sorted(len(value) for value in swapped_eigenclasses.values())
    )
    require(len(swapped_eigenclasses) == N
            and swapped_multiplicities == (1,) * N,
            ("swapped collision multiplicities",
             len(swapped_eigenclasses), swapped_multiplicities))
    swapped_zero_classes = tuple(
        value for key, value in swapped_eigenclasses.items() if key == (0,) * 12
    )
    require(swapped_zero_classes == ([(0, 0)],),
            ("swapped zero class", swapped_zero_classes))
    swapped_centralizer_dimension = N

    semantic = digest_json((
        parent_hashes, P, MINPOLY,
        tuple(str(value) for value in X_INTERVAL),
        tuple(str(value) for value in (2 * cos_lower, 2 * cos_upper)),
        ypoly, denominator,
        tuple(sine4),
        tuple((frequency, squared_rows[frequency]) for frequency in sorted(squared_rows)),
        lower_equalities, upper_equalities,
        tuple(sorted(swapped_lower_equalities)),
        tuple(sorted(swapped_upper_equalities)),
        tuple((key, tuple(value)) for key, value in sorted(eigenclasses.items())),
        double_classes, centralizer_dimension,
        tuple((key, tuple(value))
              for key, value in sorted(swapped_eigenclasses.items())),
        swapped_centralizer_dimension,
    ))
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3667 LRC optimal positive three-twist frame ==")
    print(f"parents_sha256_lf={parent_hashes}")
    print(f"x_minpoly={MINPOLY};isolating_interval={tuple(str(v) for v in X_INTERVAL)}")
    print("weights=a=y/(y+1),b=1/(y+1),negative=-1;y=2*cos(2*pi/13)")
    print("maximin_gap_squared=(4-x^2)/(x^2-1)^2")
    print(f"lower_equality_frequencies={lower_equalities};strict_others:{len(lower_intervals)}")
    print("upper_multiplier_squared=4*cos(pi/26)^2")
    print(f"upper_equality_frequencies={upper_equalities};strict_others:{len(upper_intervals)}")
    print("optimizers=(a,b)=(y/(y+1),1/(y+1)) and its swap")
    print(f"swapped_lower_equality_frequencies={tuple(sorted(swapped_lower_equalities))}")
    print(f"swapped_upper_equality_frequencies={tuple(sorted(swapped_upper_equalities))}")
    print(f"a_star_orientation_eigenvalue_classes={len(eigenclasses)};singletons:143;doubles:13")
    print(f"a_star_orientation_double_classes={double_classes}")
    print(f"a_star_orientation_centralizer_dimension={centralizer_dimension}")
    print(f"swapped_orientation_eigenvalue_classes={len(swapped_eigenclasses)};singletons:169;doubles:0")
    print(f"swapped_orientation_centralizer_dimension={swapped_centralizer_dimension}")
    print("orientation_asymmetry=optimal conditioning and simple spectrum coexist at swapped optimizer")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=positive three-site convolution frames on F13^2;not covering-row nonvanishing or LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
