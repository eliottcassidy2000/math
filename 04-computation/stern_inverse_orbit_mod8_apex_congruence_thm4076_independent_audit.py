#!/usr/bin/env python3
"""Independent exact audit for THM-4076.

This path reconstructs Stern signs from Euclidean continued-fraction depth,
computes the apex imbalance directly over every lower vertex, and counts the
even-even modular hyperbola in half-box coordinates.  It does not import or
use the primary inverse-parity implementation.
"""

from __future__ import annotations

import hashlib
import json
from math import gcd


Q_MAX = 2999
CONTROLS = (3, 5, 9, 15, 25, 29, 35, 37, 45, 65, 105)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def stern_depth_sign(numerator: int, denominator: int) -> int:
    require(0 < numerator < denominator, "depth input outside open unit interval")
    require(gcd(numerator, denominator) == 1, "depth input is not primitive")
    digit_sum = 0
    top = numerator
    bottom = denominator
    while top:
        digit_sum += bottom // top
        bottom, top = top, bottom % top
    depth = digit_sum - 1
    return 1 if depth % 2 == 0 else -1


def factor_data(n: int) -> tuple[list[tuple[int, int]], int, int]:
    factors: list[tuple[int, int]] = []
    remaining = n
    prime = 3
    while prime * prime <= remaining:
        if remaining % prime == 0:
            exponent = 0
            while remaining % prime == 0:
                remaining //= prime
                exponent += 1
            factors.append((prime, exponent))
        prime += 2
    if remaining > 1:
        factors.append((remaining, 1))
    return factors, len(factors), sum(exponent for _, exponent in factors)


def divisors_from_factors(factors: list[tuple[int, int]]) -> list[int]:
    divisors = [1]
    for prime, exponent in factors:
        old = tuple(divisors)
        power = 1
        for _ in range(exponent):
            power *= prime
            divisors.extend(divisor * power for divisor in old)
    return sorted(divisors)


def packet_from_depth(q: int) -> tuple[int, int]:
    phi = 0
    packet_sum = 0
    for a in range(1, q):
        if gcd(a, q) == 1:
            phi += 1
            packet_sum += stern_depth_sign(a, q)
    return phi, packet_sum


def direct_apex_imbalance(q: int) -> int:
    total = 0
    for a in range(1, q):
        common = gcd(a, q)
        total += stern_depth_sign(a // common, q // common)
    return total


def half_box_hyperbola_count(q: int) -> int:
    half = (q - 1) // 2
    count = 0
    for r in range(1, half + 1):
        if gcd(r, q) != 1:
            continue
        s = pow(4 * r, -1, q)
        if s <= half:
            count += 1
    return count


def direct_root_counts(q: int) -> tuple[int, int]:
    root_list = [x for x in range(1, q) if x * x % q == 1]
    return len(root_list), sum(x % 2 == 0 for x in root_list)


def main() -> None:
    packets: dict[int, int] = {}
    packet_data: dict[int, tuple[int, int, int, int, int]] = {}
    rows: list[dict[str, int]] = []
    depth_evaluations = 0
    direct_apex_evaluations = 0
    zero_packets: list[int] = []

    for q in range(3, Q_MAX + 1, 2):
        factors, omega, big_omega = factor_data(q)
        phi, packet_sum = packet_from_depth(q)
        depth_evaluations += phi
        half_box = half_box_hyperbola_count(q)
        roots, even_roots = direct_root_counts(q)

        require(packet_sum == 4 * half_box - phi, f"depth/hyperbola identity failed at q={q}")
        require(roots == 2**omega, f"root count failed at q={q}")
        require(even_roots * 2 == roots, f"root parity split failed at q={q}")
        require((half_box - even_roots) % 2 == 0, f"orbit parity failed at q={q}")
        require(
            packet_sum % 8 == (-phi + (4 if omega == 1 else 0)) % 8,
            f"packet mod-eight failed at q={q}",
        )

        if packet_sum == 0:
            zero_packets.append(q)
            if omega == 1:
                require(factors[0][0] % 8 == 5, f"prime-power zero gate failed at q={q}")
            elif omega == 2:
                require(
                    any(prime % 4 == 1 for prime, _ in factors),
                    f"two-prime zero gate failed at q={q}",
                )

        packets[q] = packet_sum
        packet_data[q] = (phi, omega, big_omega, half_box, roots)

    for q in range(3, Q_MAX + 1, 2):
        factors, omega, big_omega = factor_data(q)
        phi, _, _, half_box, roots = packet_data[q]
        packet_sum = packets[q]
        direct_imbalance = direct_apex_imbalance(q)
        direct_apex_evaluations += q - 1
        divisor_imbalance = sum(
            packets[divisor]
            for divisor in divisors_from_factors(factors)
            if divisor > 1
        )
        require(direct_imbalance == divisor_imbalance, f"direct divisor-star failed at q={q}")
        require(
            direct_imbalance % 8 == (1 - q + 4 * big_omega) % 8,
            f"B mod-eight failed at q={q}",
        )

        direct_plus = 0
        direct_minus = 0
        for a in range(1, q):
            common = gcd(a, q)
            sign = stern_depth_sign(a // common, q // common)
            if sign == 1:
                direct_plus += 1
            else:
                direct_minus += 1
        direct_apex_evaluations += q - 1
        require(direct_plus - direct_minus == direct_imbalance, f"degree balance failed at q={q}")
        require(direct_plus + direct_minus == q - 1, f"degree total failed at q={q}")
        require(direct_plus % 4 == (2 * big_omega) % 4, f"indegree mod-four failed at q={q}")
        require(
            direct_minus % 4 == (q - 1 - 2 * big_omega) % 4,
            f"outdegree mod-four failed at q={q}",
        )
        if direct_imbalance == 0:
            require(q % 8 == (1 + 4 * big_omega) % 8, f"B-zero gate failed at q={q}")

        rows.append(
            {
                "q": q,
                "phi": phi,
                "omega": omega,
                "Omega": big_omega,
                "S": packet_sum,
                "N_half_box": half_box,
                "roots": roots,
                "B_direct": direct_imbalance,
                "indegree_direct": direct_plus,
                "outdegree_direct": direct_minus,
            }
        )

    row_by_q = {row["q"]: row for row in rows}
    require(len(zero_packets) == 67, "zero-packet aggregate changed")
    require(
        row_by_q[13]["S"] % 16 == 0 and row_by_q[13]["N_half_box"] % 4 == 3,
        "q=13 mod-sixteen control failed",
    )
    require(
        row_by_q[29]["S"] % 16 == 8 and row_by_q[29]["N_half_box"] % 4 == 1,
        "q=29 mod-sixteen hostile failed",
    )
    require(row_by_q[25]["S"] == -8, "q=25 prime-power zero-gate hostile changed")
    require(row_by_q[15]["S"] == 8, "q=15 two-prime zero-gate hostile changed")

    digest = hashlib.sha256(
        json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    controls = [row for row in rows if row["q"] in CONTROLS]
    print("THM-4076 independent Euclidean-depth audit")
    print(f"universe=odd q in [3,{Q_MAX}] rows={len(rows)}")
    print(f"primitive_depth_evaluations={depth_evaluations}")
    print(f"direct_apex_depth_evaluations={direct_apex_evaluations}")
    print(f"zero_S_count={len(zero_packets)}")
    print(f"controls={json.dumps(controls, sort_keys=True, separators=(',', ':'))}")
    print("mod16_hostile=q13:S=0,Nmod4=3 versus q29:S=8,Nmod4=1")
    print(f"semantic_sha256={digest}")
    print("PASS: Euclidean depth, half-box hyperbola, direct stars, and degree counts")


if __name__ == "__main__":
    main()
