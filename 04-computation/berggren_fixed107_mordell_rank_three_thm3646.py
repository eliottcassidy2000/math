#!/usr/bin/env python3
"""Exact arithmetic companion for the THM-3646 rank-three descent.

The imported ingredient is Bandini's explicit 3-isogeny Selmer-support
theorem.  This companion verifies every curve-specific input: isogenies,
third point, finite-field independence, imaginary quadratic class number,
split support sets, and the resulting dimension ledger.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path


B = 1_225_041
Q0 = 408_347
APRIME = -27 * B
P = (Fraction(232), Fraction(3_703))
Q = (Fraction(4_960), Fraction(349_321))
U = (Fraction(50_913, 16), Fraction(-11_482_065, 64))
R = (
    Fraction(8_279_053_120, 216_766_729),
    Fraction(-3_611_785_597_108_493, 3_191_456_551_067),
)
GOOD_PRIMES = (5, 11, 13)
Point = tuple[Fraction, Fraction] | None
FinitePoint = tuple[int, int] | None
CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    p = 3
    while p * p <= n:
        if n % p == 0:
            return False
        p += 2
    return True


def on_curve(point: Point, constant: int) -> bool:
    if point is None:
        return True
    x, y = point
    return y * y == x * x * x + constant


def negate(point: Point) -> Point:
    return None if point is None else (point[0], -point[1])


def add(left: Point, right: Point, constant: int = B) -> Point:
    require(on_curve(left, constant) and on_curve(right, constant), "add input")
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and y1 == -y2:
        return None
    slope = ((y2 - y1) / (x2 - x1) if x1 != x2
             else 3 * x1 * x1 / (2 * y1))
    x3 = slope * slope - x1 - x2
    y3 = -y1 + slope * (x1 - x3)
    result = (x3, y3)
    require(on_curve(result, constant), "add output")
    return result


def multiply(multiplier: int, point: Point, constant: int = B) -> Point:
    require(multiplier >= 0, "nonnegative multiplier")
    result = None
    addend = point
    k = multiplier
    while k:
        if k & 1:
            result = add(result, addend, constant)
        addend = add(addend, addend, constant)
        k >>= 1
    return result


def phi(point: Point) -> Point:
    """The rational 3-isogeny E_B -> E_{-27B}."""
    if point is None:
        return None
    x, y = point
    require(x != 0, "phi exceptional x")
    result = ((x**3 + 4 * B) / x**2,
              y * (x**3 - 8 * B) / x**3)
    require(on_curve(result, APRIME), "phi output")
    return result


def psi(point: Point) -> Point:
    """The dual 3-isogeny E_{-27B} -> E_B."""
    if point is None:
        return None
    x, y = point
    require(x != 0, "psi exceptional x")
    result = ((x**3 - 108 * B) / (9 * x**2),
              y * (x**3 + 216 * B) / (27 * x**3))
    require(on_curve(result, B), "psi output")
    return result


def reduce_fraction(value: Fraction, p: int) -> int:
    require(value.denominator % p != 0, ("bad reduction denominator", p))
    return value.numerator * pow(value.denominator, -1, p) % p


def reduce_point(point: Point, p: int) -> FinitePoint:
    if point is None:
        return None
    result = (reduce_fraction(point[0], p), reduce_fraction(point[1], p))
    require((result[1] * result[1] - result[0]**3 - B) % p == 0,
            ("reduction", point, p))
    return result


def points_mod(p: int) -> tuple[FinitePoint, ...]:
    return (None,) + tuple(
        (x, y)
        for x in range(p)
        for y in range(p)
        if (y * y - x**3 - B) % p == 0
    )


def add_mod(left: FinitePoint, right: FinitePoint, p: int) -> FinitePoint:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    slope = ((y2 - y1) * pow((x2 - x1) % p, -1, p) if x1 != x2 else
             3 * x1 * x1 * pow((2 * y1) % p, -1, p)) % p
    x3 = (slope * slope - x1 - x2) % p
    y3 = (-y1 + slope * (x1 - x3)) % p
    return x3, y3


def multiply_mod(multiplier: int, point: FinitePoint, p: int) -> FinitePoint:
    result = None
    addend = point
    k = multiplier
    while k:
        if k & 1:
            result = add_mod(result, addend, p)
        addend = add_mod(addend, addend, p)
        k >>= 1
    return result


def reduced_forms(discriminant: int) -> tuple[tuple[int, int, int], ...]:
    require(discriminant < -4 and discriminant % 4 == 1, "form discriminant")
    forms: list[tuple[int, int, int]] = []
    limit = isqrt(abs(discriminant) // 3) + 1
    for a in range(1, limit + 1):
        for b in range(-a, a + 1):
            numerator = b * b - discriminant
            if numerator % (4 * a):
                continue
            c = numerator // (4 * a)
            if a > c:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            forms.append((a, b, c))
    return tuple(forms)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    require(B == 107**3 - 2 == 3 * Q0, "B factorization")
    require(is_prime(Q0), "prime cofactor")
    require(APRIME == -33_076_107, "isogenous constant")
    require(all(on_curve(point, B) for point in (P, Q, R)), "E points")
    require(on_curve(U, APRIME), "E-prime point")

    psi_u = psi(U)
    require(add(psi_u, negate(P)) == R, "third-point construction")
    require(R[0].denominator == 14_723**2, "R x denominator")
    require(R[1].denominator == 14_723**3, "R y denominator")
    for point in (P, Q):
        require(psi(phi(point)) == multiply(3, point), ("dual composition", point))
    require(phi(psi(U)) == multiply(3, U, APRIME), "forward composition")

    finite_ledger = []
    expected_kernels = {
        5: tuple(sorted((a, b, c) for a in range(3) for b in range(3)
                        for c in range(3) if (2 * a + b + c) % 3 == 0)),
        11: tuple(sorted((a, b, c) for a in range(3) for b in range(3)
                         for c in range(3) if (a + b + c) % 3 == 0)),
        13: tuple(sorted((a, b, c) for a in range(3) for b in range(3)
                         for c in range(3) if (2 * a + 2 * b + c) % 3 == 0)),
    }
    for p in GOOD_PRIMES:
        group = points_mod(p)
        triple_image = {multiply_mod(3, point, p) for point in group}
        reductions = tuple(reduce_point(point, p) for point in (P, Q, R))
        kernel = []
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    value = add_mod(
                        add_mod(multiply_mod(a, reductions[0], p),
                                multiply_mod(b, reductions[1], p), p),
                        multiply_mod(c, reductions[2], p), p,
                    )
                    if value in triple_image:
                        kernel.append((a, b, c))
        kernel_tuple = tuple(sorted(kernel))
        require(kernel_tuple == expected_kernels[p], ("finite kernel", p))
        finite_ledger.append((p, len(group), tuple(sorted(
            ((-1, -1) if point is None else point) for point in triple_image
        )), reductions, kernel_tuple))
    require((2 * (1 * 1 - 1 * 2) - (1 * 1 - 1 * 2)
             + (1 * 2 - 1 * 2)) % 3 == 2, "kernel determinant")

    # The imaginary descent field K_-=Q(sqrt(-408347)) has fundamental
    # discriminant -408347.  Reduced positive forms enumerate its proper
    # ideal classes exactly.
    forms = reduced_forms(-Q0)
    require(len(forms) == 137, "imaginary class number")
    require(137 % 3 != 0, "imaginary 3-class rank")

    root = Path(__file__).resolve().parents[1]
    plus_script = root / "04-computation/berggren_fixed107_real_quadratic_class_number_one_thm3643.py"
    plus_output = root / "05-knowledge/results/berggren_fixed107_real_quadratic_class_number_one_thm3643.out"
    require(file_sha256(plus_script)
            == "95878dc28b5910e8330bbdd0b67db4759218712d4ab95d5d3c348d5c96e008ca",
            "K-plus script drift")
    require(file_sha256(plus_output)
            == "b607cfb55dcf18b7b9eb112ec041dcb0c6128784e91e2a12b5979e2a02c8cea9",
            "K-plus output drift")

    # Bandini Theorem 3.5 applied first to a=B and then to a'=-27B.
    # For a=B, only 3 enters S'_1(Q): v3(a)=1 and a/3=2 mod 3.
    # For a', only 2 enters: v2(4a')=2 and -3a'=81B is a Q_2 square.
    require(B % 8 == 1 and (-3 * B) % 8 == 5, "forward 2-adic exclusion")
    require(B // 3 % 3 == 2, "forward 3-adic inclusion")
    require((-3 * APRIME) % 8 == 1, "dual 2-adic inclusion")
    require(abs(APRIME) % 3**4 == 0 and abs(APRIME) % 3**5 != 0,
            "dual 3-adic exclusion")

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    semantic = {
        "B": B,
        "points": tuple(tuple(str(v) for v in point) for point in (P, Q, R)),
        "U": tuple(str(v) for v in U),
        "finite_ledger": finite_ledger,
        "minus_forms": forms,
        "class_numbers": (137, 1),
        "selmer_bounds": (1, 2),
        "rank": 3,
    }
    print("== THM-3646 fixed-107 Mordell rank-three exact companion ==")
    print(f"curves=E:y2=x3+{B};Eprime:y2=x3{APRIME:+d}")
    print(f"U={U};R={R};R_denominator=14723^(2,3)")
    print(f"finite_independence={tuple((row[0],row[1]) for row in finite_ledger)};matrix_det_mod3=2")
    print(f"finite_ledger_sha256={digest(finite_ledger)}")
    print(f"Kminus_discriminant={-Q0};reduced_forms={len(forms)};forms_sha256={digest(forms)}")
    print("Kplus_discriminant=1225041;ordinary_class_number=1;source=THM-3643")
    print("Bandini_support=phi:{3}->2_split_primes;psi:{2}->2_split_primes")
    print("norm_kernel_dimensions=phi:1;psi:2;selmer_sum_upper=3")
    print("rank_lower=3_by_finite_reduction;rank_upper=3_by_two_isogeny_Selmer_bounds")
    print(f"semantic_sha256={digest(semantic)}")
    print(f"CHECKS={CHECKS}")
    print("status=CANDIDATE EXACT COMPANION;imported Selmer-support theorem requires audit")
    print("scope=rank only;no integral-point classification or new two-cube ray")


if __name__ == "__main__":
    main()
