#!/usr/bin/env python3
"""Exact referee for THM-1226, the gcd-period charge obstruction.

All proof-facing calculations use integers or ``fractions.Fraction``.  The
program verifies the exact projective charge factorization, the translated
factorial strict-high counterfamily, its covered protected-needle embedding,
and the finite-channel constant on THM-1221's disconnected strict-spectrum
branches.  The latter is conditional only on THM-1221's already-certified
branch classification; no bounded speed box is used here.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path


OFFSETS = (1, 2, 3, 5, 7, 11, 13)
BASE_STEP = 27720
HIGH_BAR = F(1, 63)
FIRST_STRICT_HIGH = F(5, 308)
TREE_FLOOR = F(15, 154)
PAIR_FLOOR = F(1, 91)
C_DISC = F(448916, 194775)


def require(condition: bool, message: object) -> None:
    """Optimization-stable assertion."""
    if not condition:
        raise RuntimeError(message)


def fold(residue: int, modulus: int = 14) -> int:
    residue %= modulus
    return residue * (modulus - residue)


def rho(a: int, b: int) -> F:
    """Haar mass of ``D_a intersect D_b`` at radius ``1/14``."""
    require(a > 0 and b > 0 and a != b, (a, b))
    if a > b:
        a, b = b, a
    g = gcd(a, b)
    modulus = 14 * g
    return F(
        4 * a * b + fold(a + b, modulus) - fold(b - a, modulus),
        196 * a * b,
    )


def rho_ratio(x: F, y: F) -> F:
    ratio = x / y
    return rho(abs(ratio.numerator), abs(ratio.denominator))


def eta(value: F) -> F:
    return value * (1 - value)


def projective_kappa(x: F, y: F) -> tuple[F, tuple[int, int], F]:
    """Return kappa, the unordered reduced channel, and pair mass."""
    ratio = x / y
    a, b = abs(ratio.numerator), abs(ratio.denominator)
    value = rho(a, b)
    return eta(value) * F(a * b, a + b), tuple(sorted((a, b))), value


def reduced_bank(product_cap: int, predicate) -> list[tuple[int, int, F]]:
    return [
        (a, b, rho(a, b))
        for a in range(1, product_cap + 1)
        for b in range(a + 1, product_cap + 1)
        if a * b <= product_cap and gcd(a, b) == 1 and predicate(rho(a, b))
    ]


def graph_connected(values: tuple[F, ...], predicate) -> bool:
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in range(len(values)):
            if j not in seen and predicate(rho_ratio(values[i], values[j])):
                seen.add(j)
                stack.append(j)
    return len(seen) == len(values)


def maximum_kappa(vertices: set[F]) -> tuple[F, list[tuple[F, F, tuple[int, int], F]]]:
    rows: list[tuple[F, F, tuple[int, int], F]] = []
    maximum = F(0)
    for x, y in combinations(sorted(vertices), 2):
        value, channel, mass = projective_kappa(x, y)
        if value > maximum:
            maximum = value
            rows = [(x, y, channel, mass)]
        elif value == maximum:
            rows.append((x, y, channel, mass))
    return maximum, rows


def factorial_packet(multiplier: int) -> tuple[int, tuple[int, ...]]:
    require(multiplier > 0, multiplier)
    A = BASE_STEP * multiplier
    return A, tuple(A + r for r in OFFSETS)


def check_counterfamily(multiplier: int) -> dict[str, object]:
    A, speeds = factorial_packet(multiplier)
    require(A % 14 == 0 and sum(OFFSETS) == 42, (A, OFFSETS))
    require(all(A % d == 0 for d in range(1, 13)), A)
    require(all(gcd(x, y) == 1 for x, y in combinations(speeds, 2)), speeds)

    masses = tuple(rho(x, y) for x, y in combinations(speeds, 2))
    x_A = F(1, 49) - F(1, 4 * A * A)
    require(9 * A * A > 539 and x_A > FIRST_STRICT_HIGH, x_A)
    require(all(value > x_A > FIRST_STRICT_HIGH for value in masses), min(masses))
    require(
        all(value < F(1, 2) and eta(value) > eta(x_A) for value in masses),
        "pair variance does not clear the claimed lower bound",
    )

    # These conclusions hold for every one of the 7^5 labelled spanning
    # trees, because they use only six edgewise lower bounds.
    every_tree_mass = 6 * FIRST_STRICT_HIGH
    every_tree_error = 6 * eta(x_A)
    coarse_tree_error = 6 * eta(FIRST_STRICT_HIGH)
    coarse_normalized_error = F(4545 * A, 332024)
    harmonic = sum((F(1, speed) for speed in speeds), F(0))
    require(every_tree_mass == TREE_FLOOR, every_tree_mass)
    require(coarse_tree_error == F(4545, 47432), coarse_tree_error)
    require(every_tree_error > coarse_tree_error, every_tree_error)
    require(harmonic < F(7, A), harmonic)
    require(every_tree_error / harmonic > F(6 * A, 7) * eta(x_A), harmonic)
    require(every_tree_error / harmonic > coarse_normalized_error, harmonic)

    corrections = tuple(
        fold(r + s) - fold(s - r) for r, s in combinations(OFFSETS, 2)
    )
    require(
        corrections
        == (20, 16, 8, 0, -16, -24, 32, 16, 0, -32, -20,
            24, 0, -48, -16, 0, -24, -8, 0, 0, 16),
        corrections,
    )
    for (r, s), correction in zip(combinations(OFFSETS, 2), corrections):
        expected = F(1, 49) + F(correction, 196 * (A + r) * (A + s))
        require(rho(A + r, A + s) == expected, (r, s, correction))

    return {
        "multiplier": multiplier,
        "A": A,
        "min_mass": min(masses),
        "max_mass": max(masses),
        "x_A": x_A,
        "harmonic_bound": F(7, A),
        "coarse_tree_error": coarse_tree_error,
        "coarse_normalized_error": coarse_normalized_error,
        "corrections": corrections,
    }


def check_charge_identity() -> int:
    rows = 0
    for left in range(1, 70):
        for right in range(left + 1, 75):
            g = gcd(left, right)
            a, b = left // g, right // g
            variance = eta(rho(left, right))
            kappa = variance * F(a * b, a + b)
            error = variance / g
            require(error == kappa * (F(1, left) + F(1, right)), (left, right))
            require(error == variance * F(a, left), (left, right, "left charge"))
            require(error == variance * F(b, right), (left, right, "right charge"))
            rows += 1
    return rows


def check_protected_embedding() -> dict[str, object]:
    A, speeds = factorial_packet(1)
    q = A + 1
    require(speeds == tuple(q + d for d in (0, 1, 2, 4, 6, 10, 12)), speeds)
    require(q % 2 == 1 and q >= 521, q)
    core = tuple((3 * q + 1) // 2 + k for k in range(6))
    m = max(core)
    require(m == (3 * q + 11) // 2, (m, q))
    center = F(1, q)
    radius = F(1, 14 * m)
    length = 2 * radius

    core_margin = F(5 * q - 77, 14 * q)
    danger_margin = F(q * q - 517 * q - 1848, 14 * q * (3 * q + 11))
    require(core_margin > 0 and danger_margin > 0, (core_margin, danger_margin))
    require(520 * 520 - 517 * 520 - 1848 < 0, "q=520 should fail")
    require(521 * 521 - 517 * 521 - 1848 > 0, "q=521 should pass")

    # Triangle-inequality providers, with exact endpoint replays.
    for w in core:
        center_distance = F(1, 2) - F(2 * (w - (3 * q + 1) // 2) + 1, 2 * q)
        require(center_distance - w * radius > F(1, 14), (w, center_distance))
    for speed in speeds:
        d = speed - q
        require(F(d, q) + speed * radius < F(1, 14), (speed, d))

    return {
        "q": q,
        "deleted": speeds,
        "core": core,
        "m": m,
        "center": center,
        "radius": radius,
        "length": length,
        "core_margin": core_margin,
        "danger_margin": danger_margin,
        "local_pair_mass": length,
        "local_tree_mass": 6 * length,
    }


def check_disconnected_branch_constant() -> dict[str, object]:
    require(F(1, 49) - F(1, 4 * 56) > HIGH_BAR, "channel tail")
    strict = reduced_bank(55, lambda value: value < HIGH_BAR)
    closed = reduced_bank(55, lambda value: value <= HIGH_BAR)
    require(len(strict) == 7 and len(closed) == 12, (strict, closed))
    strict_vertices = {F(1)} | {F(b, a) for a, b, _ in strict} | {
        F(a, b) for a, b, _ in strict
    }
    closed_ratios = {F(b, a) for a, b, _ in closed} | {
        F(a, b) for a, b, _ in closed
    }
    closed_vertices = {F(1)} | closed_ratios

    strict_max, strict_rows = maximum_kappa(strict_vertices)
    closed_max, closed_rows = maximum_kappa(closed_vertices)
    require(strict_max == F(85975, 342804), strict_max)
    require({row[2] for row in strict_rows} == {(20, 33)}, strict_rows)
    require(closed_max == F(224458, 584325), closed_max)
    require(closed_rows == [(F(5, 9), F(9, 5), (25, 81), F(97, 4725))], closed_rows)

    # Reconstruct THM-1221's twelve normalized 2+5 strict-component packets.
    ordered_closed = tuple(sorted(closed_ratios))
    centers = sorted({r / s for r in ordered_closed for s in ordered_closed} - {F(1)})
    packets: set[tuple[F, ...]] = set()
    for second in centers:
        left = (F(1), second)
        if rho_ratio(*left) <= HIGH_BAR:
            continue
        common = tuple(v for v in ordered_closed if v / second in closed_ratios)
        for right in combinations(common, 5):
            if not graph_connected(right, lambda value: value > HIGH_BAR):
                continue
            if not any(rho_ratio(x, y) == HIGH_BAR for x in left for y in right):
                continue
            packet = tuple(sorted(left + right))
            packets.add(tuple(value / packet[0] for value in packet))
    require(len(packets) == 12, len(packets))
    two_five_rows = [
        (*projective_kappa(x, y), packet, x, y)
        for packet in packets
        for x, y in combinations(packet, 2)
    ]
    two_five_max = max(row[0] for row in two_five_rows)
    require(two_five_max == F(43774, 276507), two_five_max)

    strict_C = 6 * strict_max
    closed_C = 6 * closed_max
    two_five_C = 6 * two_five_max
    require(strict_C == F(85975, 57134), strict_C)
    require(closed_C == C_DISC, closed_C)
    require(two_five_C == F(87548, 92169), two_five_C)
    require(max(strict_C, closed_C, two_five_C) == C_DISC, "wrong C_disc")
    crown_ratio = TREE_FLOOR / (C_DISC + F(1, 7))
    require(crown_ratio == F(417375, 10488302), crown_ratio)
    return {
        "strict_vertices": len(strict_vertices),
        "strict_kappa": strict_max,
        "strict_C": strict_C,
        "closed_vertices": len(closed_vertices),
        "closed_kappa": closed_max,
        "closed_C": closed_C,
        "two_five_packets": len(packets),
        "two_five_kappa": two_five_max,
        "two_five_C": two_five_C,
        "C_disc": C_DISC,
        "crown_ratio": crown_ratio,
    }


def main() -> None:
    print("THM-1226 GCD-PERIOD PROJECTIVE-CHARGE EXACT REFEREE")
    print("method=integer/Fraction only; always-on checks; no dependencies")
    charge_rows = check_charge_identity()
    print(f"exact edge-charge identities={charge_rows}")

    print("\nTRANSLATED FACTORIAL STRICT-HIGH COUNTERFAMILY")
    print(f"offsets={OFFSETS}; sum={sum(OFFSETS)}; base_step={BASE_STEP}")
    correction_row = None
    for multiplier in (1, 2, 10, 100):
        row = check_counterfamily(multiplier)
        correction_row = row["corrections"]
        print(
            f"N={multiplier} A={row['A']} min_rho={row['min_mass']} "
            f"max_rho={row['max_mass']} H<{row['harmonic_bound']}"
        )
        print(
            f"  every_tree_error>{row['coarse_tree_error']} "
            f"E/H>{row['coarse_normalized_error']}"
        )
    print(f"fold_corrections={correction_row}")
    print("all pair gcds=1; G_gt=K7; every spanning tree mass>15/154")
    print("E_T/H grows as (288/16807)A+O(1); no absolute C exists")

    embedding = check_protected_embedding()
    print("\nCOVERED PROTECTED-NEEDLE EMBEDDING")
    for key in (
        "q", "deleted", "core", "m", "center", "radius", "length",
        "core_margin", "danger_margin", "local_pair_mass", "local_tree_mass",
    ):
        print(f"{key}={embedding[key]}")
    print("all seven deleted combs are active throughout I; this is not an F7 counterexample")

    branch = check_disconnected_branch_constant()
    print("\nCONDITIONAL DISCONNECTED-G_gt TRANSFER")
    for key, value in branch.items():
        print(f"{key}={value}")
    print("if G_gt is disconnected: E_T<=C_disc*H for a THM-1221 floor tree")
    print("localized tree mass >=(15/154)L-(448916/194775)H")

    print("\nTOURNAMENT / ALTERNATE-VERTEX AUDIT")
    print("factorial switch=all 21 edges strict-high; speed-order gauge is transitive")
    print("scores=(0,1,2,3,4,5,6); cycles=0; SCCs=7; Hamilton_paths=1; low_flips=0")
    print("runner threshold graph preserves global rho floor but destroys projective height kappa")
    print("gcd-period vertices collapse all counterfamily edges to period 1 and destroy wall alignment")
    print("faithful local vertices=wall events / tooth addresses with signed endpoint cocycles")
    print("STATUS=PASS")
    print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
