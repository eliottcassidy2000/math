#!/usr/bin/env python3
"""Exact certificate for the fixed-107 Mordell rank-two theorem target.

Standard-library integer and rational arithmetic only.  Every load-bearing
check uses ``require`` so ``python`` and ``python -O`` execute the same gates.
The bounded integral-point scan is explicitly finite evidence and is not used
to infer a complete integral-point classification.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = (
    ROOT
    / "01-canon/theorems/"
    / "THM-3370-berggren-two-cube-biquadratic-norm-collision.md"
)
EXPECTED_PARENT_SHA256_LF = (
    "9abf73a45d789fe5a804977235a3bce3e415c3ad6da865a628d8f338159fa53a"
)
EXPECTED_SEMANTIC_SHA256 = (
    "95622596e5eababb0bbfeb8ddf28d0587aabe0ad04408c0a44d350f20a8955bb"
)

B = 1_225_041
PRIME_FACTOR = 408_347
P = (232, 3_703)
Q = (4_960, 349_321)
GOOD_PRIMES = (7, 13, 19)
SCAN_LO = -107
SCAN_HI = 10_000_000
EXPECTED_SCAN_CANDIDATES = 4_422
EXPECTED_SCAN_POINTS = ((232, 3_703), (4_960, 349_321))

CHECKS = 0


def require(condition: bool, payload: object) -> None:
    """Optimization-safe exact gate."""
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return sha256(payload).hexdigest()


def is_prime_trial(number: int) -> bool:
    """Deterministic trial-division primality for the one pinned factor."""
    if number < 2:
        return False
    if number % 2 == 0:
        return number == 2
    divisor = 3
    limit = isqrt(number)
    while divisor <= limit:
        if number % divisor == 0:
            return False
        divisor += 2
    return True


Point = tuple[int, int] | None
RationalPoint = tuple[Fraction, Fraction] | None


def point_key(point: Point) -> tuple[int, int]:
    return (-1, -1) if point is None else point


def curve_points_mod(prime: int) -> tuple[Point, ...]:
    points: list[Point] = [None]
    for x_value in range(prime):
        rhs = (x_value**3 + B) % prime
        for y_value in range(prime):
            if y_value * y_value % prime == rhs:
                points.append((x_value, y_value))
    return tuple(points)


def point_neg_mod(point: Point, prime: int) -> Point:
    if point is None:
        return None
    return point[0], (-point[1]) % prime


def point_add_mod(left: Point, right: Point, prime: int) -> Point:
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right and (y_left + y_right) % prime == 0:
        return None
    if left == right:
        require(y_left % prime != 0, (prime, "invalid doubling denominator", left))
        slope = 3 * x_left * x_left * pow(2 * y_left, -1, prime) % prime
    else:
        require((x_right - x_left) % prime != 0,
                (prime, "invalid addition denominator", left, right))
        slope = (y_right - y_left) * pow(x_right - x_left, -1, prime) % prime
    x_sum = (slope * slope - x_left - x_right) % prime
    y_sum = (-y_left - slope * (x_sum - x_left)) % prime
    return x_sum, y_sum


def point_mul_mod(multiplier: int, point: Point, prime: int) -> Point:
    require(multiplier >= 0, (prime, "negative multiplier", multiplier))
    result: Point = None
    addend = point
    exponent = multiplier
    while exponent:
        if exponent & 1:
            result = point_add_mod(result, addend, prime)
        addend = point_add_mod(addend, addend, prime)
        exponent >>= 1
    return result


def rational_neg(point: RationalPoint) -> RationalPoint:
    if point is None:
        return None
    return point[0], -point[1]


def rational_add(left: RationalPoint, right: RationalPoint) -> RationalPoint:
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right and y_left == -y_right:
        return None
    if left == right:
        require(y_left != 0, ("rational doubling denominator", left))
        slope = 3 * x_left * x_left / (2 * y_left)
    else:
        require(x_left != x_right, ("rational addition denominator", left, right))
        slope = (y_right - y_left) / (x_right - x_left)
    x_sum = slope * slope - x_left - x_right
    y_sum = -y_left - slope * (x_sum - x_left)
    return x_sum, y_sum


def poly_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    size = max(len(left), len(right))
    values = [0] * size
    for index, value in enumerate(left):
        values[index] += value
    for index, value in enumerate(right):
        values[index] += value
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def poly_scale(polynomial: tuple[int, ...], scalar: int) -> tuple[int, ...]:
    return tuple(scalar * value for value in polynomial)


def poly_mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    values = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            values[i + j] += left_value * right_value
    return tuple(values)


def chord_polynomial(line_x: int, line_constant: int, y_scale: int) -> tuple[int, ...]:
    line = (line_constant, line_x)
    curve = (B, 0, 0, 1)
    return poly_add(poly_mul(line, line), poly_scale(curve, -(y_scale**2)))


def factored_chord(third_constant: int, third_x: int) -> tuple[int, ...]:
    factors = poly_mul((-P[0], 1), (-Q[0], 1))
    factors = poly_mul(factors, (third_constant, third_x))
    return poly_scale(factors, -1)


def fraction_payload(point: RationalPoint) -> object:
    if point is None:
        return None
    return tuple((value.numerator, value.denominator) for value in point)


def residue_support(modulus: int) -> tuple[int, ...]:
    return tuple(
        k_value
        for k_value in range(modulus)
        if any(
            ((2 * k_value + 1) ** 2 - x_value**3 - B) % modulus == 0
            for x_value in range(modulus)
        )
    )


def square_x_residues(modulus: int) -> frozenset[int]:
    squares = {value * value % modulus for value in range(modulus)}
    return frozenset(
        x_value
        for x_value in range(modulus)
        if (x_value**3 + B) % modulus in squares
    )


def bounded_integral_scan() -> tuple[int, tuple[tuple[int, int], ...], object]:
    """Finite-exact scan; modular wheel only prunes impossible residues."""
    wheel_moduli = (64, 9, 5, 7, 11)
    extra_moduli = (13, 17, 19, 23, 29, 31)
    allowed = {
        modulus: square_x_residues(modulus)
        for modulus in wheel_moduli + extra_moduli
    }
    wheel = 1
    for modulus in wheel_moduli:
        wheel *= modulus
    require(wheel == 221_760, ("wheel modulus", wheel))
    wheel_residues = tuple(
        residue
        for residue in range(wheel)
        if all(residue % modulus in allowed[modulus] for modulus in wheel_moduli)
    )

    points: list[tuple[int, int]] = []
    candidate_count = 0
    first_block = SCAN_LO // wheel
    last_block = SCAN_HI // wheel
    for block in range(first_block, last_block + 1):
        base = block * wheel
        for residue in wheel_residues:
            x_value = base + residue
            if x_value < SCAN_LO or x_value > SCAN_HI:
                continue
            if not all(x_value % modulus in allowed[modulus]
                       for modulus in extra_moduli):
                continue
            candidate_count += 1
            rhs = x_value**3 + B
            if rhs < 0:
                continue
            y_value = isqrt(rhs)
            if y_value * y_value == rhs:
                points.append((x_value, y_value))
    residue_ledger = tuple(
        (modulus, len(allowed[modulus]))
        for modulus in wheel_moduli + extra_moduli
    )
    return candidate_count, tuple(points), (wheel, len(wheel_residues), residue_ledger)


def main() -> None:
    parent_hash = lf_sha256(PARENT_PATH)
    require(parent_hash == EXPECTED_PARENT_SHA256_LF,
            ("THM-3370 parent hash", parent_hash, EXPECTED_PARENT_SHA256_LF))

    require(B == 107**3 - 2 == 3 * PRIME_FACTOR, ("B factorization", B))
    require(is_prime_trial(PRIME_FACTOR), ("factor primality", PRIME_FACTOR))
    discriminant = -432 * B * B
    require(discriminant == -(2**4) * (3**5) * (PRIME_FACTOR**2),
            ("displayed discriminant factorization", discriminant))
    good_reduction = tuple(discriminant % prime != 0 for prime in GOOD_PRIMES)
    require(good_reduction == (True, True, True),
            ("good reduction", GOOD_PRIMES, good_reduction))

    for label, point in (("P", P), ("Q", Q)):
        require(point[1] ** 2 == point[0] ** 3 + B,
                (label, "curve identity", point))
    depth_p = (P[1] - 1) // 2
    depth_q = (Q[1] - 1) // 2
    require((depth_p, depth_q) == (1_851, 174_660),
            ("Berggren depths", depth_p, depth_q))
    for point, depth in ((P, depth_p), (Q, depth_q)):
        require(point[0] ** 3 + 107**3 == point[1] ** 2 + 2,
                ("fixed-summand identity", point))
        require(point[1] ** 2 + 2 == 4 * depth * (depth + 1) + 3,
                ("U-spine identity", point, depth))

    finite_points = {prime: curve_points_mod(prime) for prime in GOOD_PRIMES}
    finite_counts = tuple((prime, len(finite_points[prime])) for prime in GOOD_PRIMES)
    require(finite_counts == ((7, 4), (13, 12), (19, 27)),
            ("finite field counts", finite_counts))
    torsion_bound = gcd(*(count for _prime, count in finite_counts))
    require(torsion_bound == 1, ("torsion gcd bound", torsion_bound))

    expected_triples = {
        13: (None, (1, 0), (3, 0), (9, 0)),
        19: (None, (0, 4), (0, 15)),
    }
    expected_kernels = {
        13: ((0, 0), (1, 2), (2, 1)),
        19: ((0, 0), (1, 1), (2, 2)),
    }
    quotient_ledger = []
    kernel_sets: dict[int, set[tuple[int, int]]] = {}
    for prime in (13, 19):
        points = finite_points[prime]
        triple_image = tuple(sorted(
            {point_mul_mod(3, point, prime) for point in points},
            key=point_key,
        ))
        require(triple_image == expected_triples[prime],
                (prime, "3E image", triple_image))
        p_mod = (P[0] % prime, P[1] % prime)
        q_mod = (Q[0] % prime, Q[1] % prime)
        rows = []
        kernel = set()
        for a_value in range(3):
            for b_value in range(3):
                combination = point_add_mod(
                    point_mul_mod(a_value, p_mod, prime),
                    point_mul_mod(b_value, q_mod, prime),
                    prime,
                )
                belongs = combination in triple_image
                expected = (
                    (a_value + b_value) % 3 == 0
                    if prime == 13
                    else (a_value - b_value) % 3 == 0
                )
                require(belongs == expected,
                        (prime, "coefficient condition", a_value, b_value,
                         combination, belongs, expected))
                if belongs:
                    kernel.add((a_value, b_value))
                rows.append((a_value, b_value, combination, belongs))
        kernel_tuple = tuple(sorted(kernel))
        require(kernel_tuple == expected_kernels[prime],
                (prime, "coefficient kernel", kernel_tuple))
        require(any(pair != (0, 0) for pair in kernel),
                (prime, "hostile one-prime kernel unexpectedly trivial"))
        kernel_sets[prime] = kernel
        quotient_ledger.append((prime, triple_image, p_mod, q_mod,
                                tuple(rows), kernel_tuple))

    kernel_intersection = tuple(sorted(kernel_sets[13] & kernel_sets[19]))
    require(kernel_intersection == ((0, 0),),
            ("two-prime coefficient intersection", kernel_intersection))
    rank_lower_bound = 2
    require(torsion_bound == 1 and kernel_intersection == ((0, 0),),
            "rank and saturation logic inputs")
    three_saturated = True

    p_rat = (Fraction(P[0]), Fraction(P[1]))
    q_rat = (Fraction(Q[0]), Fraction(Q[1]))
    difference = rational_add(q_rat, rational_neg(p_rat))
    point_sum = rational_add(p_rat, q_rat)
    expected_difference = (Fraction(3_448, 9), Fraction(-204_659, 27))
    expected_sum = (
        Fraction(94_164_361, 788**2),
        Fraction(1_062_189_113_125, 788**3),
    )
    require(difference == expected_difference, ("Q-P", difference))
    require(point_sum == expected_sum, ("P+Q", point_sum))
    require(difference[0].denominator == 3**2
            and difference[1].denominator == 3**3,
            ("Q-P denominators", difference))
    require(point_sum[0].denominator == (4 * 197) ** 2
            and point_sum[1].denominator == (4 * 197) ** 3,
            ("P+Q denominators", point_sum))

    difference_chord = chord_polynomial(224, -63_077, 3)
    difference_factor = factored_chord(-3_448, 9)
    sum_chord = chord_polynomial(57_603, -10_445_932, 788)
    sum_factor = factored_chord(-94_164_361, 620_944)
    require(difference_chord == difference_factor,
            ("Q-P chord factorization", difference_chord, difference_factor))
    require(sum_chord == sum_factor,
            ("P+Q chord factorization", sum_chord, sum_factor))

    p_mod_3 = (P[0] % 3, P[1] % 3)
    q_mod_3 = (Q[0] % 3, Q[1] % 3)
    require(p_mod_3 == q_mod_3 == (1, 1),
            ("collision modulo 3", p_mod_3, q_mod_3))
    p_mod_197 = (P[0] % 197, P[1] % 197)
    q_mod_197 = (Q[0] % 197, Q[1] % 197)
    require(q_mod_197 == point_neg_mod(p_mod_197, 197) == (35, 40),
            ("collision modulo 197", p_mod_197, q_mod_197))
    require(difference[0].denominator != 1 and point_sum[0].denominator != 1,
            "hostile: integral points closed under plus or minus")

    support_7 = residue_support(7)
    support_9 = residue_support(9)
    support_63 = tuple(
        k_value
        for k_value in range(63)
        if k_value % 7 in support_7 and k_value % 9 in support_9
    )
    require(support_7 == (3,), ("fixed-107 support mod 7", support_7))
    require(support_9 == (2, 6), ("fixed-107 support mod 9", support_9))
    require(support_63 == (24, 38), ("fixed-107 support mod 63", support_63))
    require((depth_p % 63, depth_q % 63) == (24, 24),
            ("known depth residues", depth_p % 63, depth_q % 63))
    require((-24 - 1) % 63 == 38 and (-38 - 1) % 63 == 24,
            "sign mirror exchanges the two depth classes")

    scan_candidates, scan_points, scan_ledger = bounded_integral_scan()
    require(scan_candidates == EXPECTED_SCAN_CANDIDATES,
            ("finite scan candidate count", scan_candidates))
    require(scan_points == EXPECTED_SCAN_POINTS,
            ("finite scan points", scan_points))

    semantic_payload = {
        "parent_sha256_lf": parent_hash,
        "B": B,
        "prime_factor": PRIME_FACTOR,
        "discriminant": discriminant,
        "good_primes": GOOD_PRIMES,
        "points": (P, Q),
        "depths": (depth_p, depth_q),
        "finite_counts": finite_counts,
        "torsion_bound": torsion_bound,
        "quotient_ledger": quotient_ledger,
        "kernel_intersection": kernel_intersection,
        "rank_lower_bound": rank_lower_bound,
        "three_saturated": three_saturated,
        "difference": fraction_payload(difference),
        "sum": fraction_payload(point_sum),
        "difference_chord": difference_chord,
        "sum_chord": sum_chord,
        "collisions": ((3, p_mod_3, q_mod_3),
                       (197, p_mod_197, q_mod_197)),
        "depth_support": ((7, support_7), (9, support_9), (63, support_63)),
        "scan": (SCAN_LO, SCAN_HI, scan_candidates, scan_points, scan_ledger),
    }
    semantic_hash = canonical_sha256(semantic_payload)
    if EXPECTED_SEMANTIC_SHA256 == "TO_BE_PINNED":
        raise RuntimeError(("UNPINNED semantic_sha256", semantic_hash))
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic_hash, EXPECTED_SEMANTIC_SHA256))

    script_hash = lf_sha256(Path(__file__).resolve())
    report = [
        "== THM-3620 fixed-107 Mordell rank-two exact companion ==",
        f"parent_thm3370_sha256_lf={parent_hash}",
        f"B={B}=3*{PRIME_FACTOR};prime_factor_is_prime=1",
        f"displayed_discriminant={discriminant};good_reduction={GOOD_PRIMES}",
        f"points=P{P},Q{Q};depths=({depth_p},{depth_q});curve_identities=PASS",
        f"finite_field_counts={finite_counts};torsion_gcd_bound={torsion_bound}",
        "three_image_13=" + repr(expected_triples[13]).replace(" ", ""),
        "coefficient_kernel_13=" + repr(expected_kernels[13]).replace(" ", "")
        + ";condition=a+b=0_mod_3",
        "three_image_19=" + repr(expected_triples[19]).replace(" ", ""),
        "coefficient_kernel_19=" + repr(expected_kernels[19]).replace(" ", "")
        + ";condition=a-b=0_mod_3",
        f"kernel_intersection={kernel_intersection};rank_lower_bound={rank_lower_bound}",
        "three_saturation=PASS;proof=3T=aP+bQ=>a=3a0,b=3b0=>"
        "3(T-a0P-b0Q)=O=>T_in_span_by_torsion_zero",
        "hostile_one_prime_only=13_leaves_(1,2),(2,1);"
        "19_leaves_(1,1),(2,2)",
        f"Q_minus_P={fraction_payload(difference)};denominators=(3^2,3^3)",
        f"P_plus_Q={fraction_payload(point_sum)};denominators=((4*197)^2,(4*197)^3)",
        "chord_Q_minus_P=(224X-63077)^2-9(X^3+B)="
        "-(X-232)(X-4960)(9X-3448)",
        "chord_P_plus_Q=(57603X-10445932)^2-788^2(X^3+B)="
        "-(X-232)(X-4960)(620944X-94164361)",
        f"collision_mod_3=P{p_mod_3}=Q{q_mod_3}",
        f"collision_mod_197=P{p_mod_197}=-Q{q_mod_197}",
        "hostile_integrality=Q-P_and_P+Q_are_nonintegral",
        f"fixed_107_depth_support=mod7:{support_7};mod9:{support_9};mod63:{support_63}",
        f"finite_exact_scan=X[{SCAN_LO},{SCAN_HI}];"
        f"modular_candidates={scan_candidates};points={scan_points}",
        "finite_scan_scope=no_complete_integral_point_classification",
        f"semantic_sha256={semantic_hash}",
        f"script_sha256_lf={script_hash}",
        "hash_basis=LF-normalized bytes;semantic=canonical sorted compact JSON",
        f"PASS explicit_require_gates={CHECKS}",
    ]
    report_payload_hash = sha256(("\n".join(report) + "\n").encode("utf-8")).hexdigest()
    report.append(f"report_payload_sha256={report_payload_hash}")
    report.append("status=EXACT CERTIFICATE;rank>=2 and 3-saturated subgroup")
    report.append("scope=no rank equality;no complete E_107(Z);finite scan only")
    print("\n".join(report))


if __name__ == "__main__":
    main()
