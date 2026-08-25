#!/usr/bin/env python3
"""Exact primary audit for THM-4076.

This path evaluates the odd-denominator packet through the inverse-parity
formula of THM-4059, audits the inversion orbits on the even-even box, and
checks the divisor-star and apex congruences for every odd q <= Q_MAX.
"""

from __future__ import annotations

import hashlib
import json
from math import gcd


Q_MAX = 4999
CONTROLS = (3, 5, 9, 15, 25, 29, 35, 37, 45, 65, 105)
EXPECTED_PRIME_ZEROS = (5, 13, 61, 181, 293, 1109, 1429, 1453, 3637)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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


def packet(q: int) -> tuple[int, int, int]:
    phi = 0
    packet_sum = 0
    even_even = 0
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        inverse = pow(a, -1, q)
        phi += 1
        packet_sum += 1 if (a + inverse) % 2 == 0 else -1
        if a % 2 == 0 and inverse % 2 == 0:
            even_even += 1
    return phi, packet_sum, even_even


def root_counts(q: int) -> tuple[int, int]:
    roots = 0
    even_roots = 0
    for x in range(1, q):
        if x * x % q == 1:
            roots += 1
            if x % 2 == 0:
                even_roots += 1
    return roots, even_roots


def main() -> None:
    rows: list[dict[str, int]] = []
    packets: dict[int, int] = {}
    zero_packets: list[int] = []
    prime_zero_packets: list[int] = []
    orbit_pairs = 0
    root_scans = 0
    divisor_terms = 0

    for q in range(3, Q_MAX + 1, 2):
        factors, omega, big_omega = factor_data(q)
        phi, packet_sum, even_even = packet(q)
        roots, even_roots = root_counts(q)
        root_scans += q - 1

        require(roots == 2**omega, f"root count failed at q={q}")
        require(even_roots * 2 == roots, f"root parity split failed at q={q}")
        fixed_even = even_roots
        require(
            (even_even - fixed_even) % 2 == 0,
            f"nonfixed inversion orbit parity failed at q={q}",
        )
        orbit_pairs += (even_even - fixed_even) // 2
        require(packet_sum == 4 * even_even - phi, f"hyperbola identity failed at q={q}")
        expected_mod8 = (-phi + (4 if omega == 1 else 0)) % 8
        require(packet_sum % 8 == expected_mod8, f"packet mod-eight failed at q={q}")

        packets[q] = packet_sum
        if packet_sum == 0:
            zero_packets.append(q)
            if omega == 1 and factors[0][1] == 1:
                prime_zero_packets.append(q)
            if omega == 1:
                require(factors[0][0] % 8 == 5, f"prime-power zero gate failed at q={q}")
            elif omega == 2:
                require(
                    any(prime % 4 == 1 for prime, _ in factors),
                    f"two-prime zero gate failed at q={q}",
                )

        rows.append(
            {
                "q": q,
                "phi": phi,
                "omega": omega,
                "Omega": big_omega,
                "S": packet_sum,
                "N": even_even,
                "roots": roots,
                "even_roots": even_roots,
            }
        )

    require(tuple(prime_zero_packets) == EXPECTED_PRIME_ZEROS, "prime-zero control changed")
    require(len(zero_packets) == 90, "zero-packet aggregate changed")

    augmented_rows: list[dict[str, int]] = []
    row_by_q = {row["q"]: row for row in rows}
    for q in range(3, Q_MAX + 1, 2):
        factors, _, big_omega = factor_data(q)
        divisors = divisors_from_factors(factors)
        nontrivial = [divisor for divisor in divisors if divisor > 1]
        imbalance = sum(packets[divisor] for divisor in nontrivial)
        divisor_terms += len(nontrivial)
        indegree = (q - 1 + imbalance) // 2
        outdegree = (q - 1 - imbalance) // 2

        require(imbalance % 8 == (1 - q + 4 * big_omega) % 8, f"B mod-eight failed at q={q}")
        require(indegree % 4 == (2 * big_omega) % 4, f"indegree mod-four failed at q={q}")
        require(
            outdegree % 4 == (q - 1 - 2 * big_omega) % 4,
            f"outdegree mod-four failed at q={q}",
        )
        if imbalance == 0:
            require(q % 8 == (1 + 4 * big_omega) % 8, f"B-zero gate failed at q={q}")

        row = dict(row_by_q[q])
        row.update({"B": imbalance, "indegree": indegree, "outdegree": outdegree})
        augmented_rows.append(row)

    require(row_by_q[13]["S"] % 16 == 0 and row_by_q[13]["N"] % 4 == 3, "q=13 mod-sixteen control failed")
    require(row_by_q[29]["S"] % 16 == 8 and row_by_q[29]["N"] % 4 == 1, "q=29 mod-sixteen hostile failed")
    require(row_by_q[25]["S"] == -8, "q=25 prime-power zero-gate hostile changed")
    require(row_by_q[15]["S"] == 8, "q=15 two-prime zero-gate hostile changed")

    digest = hashlib.sha256(
        json.dumps(augmented_rows, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    controls = [row for row in augmented_rows if row["q"] in CONTROLS]

    print("THM-4076 primary exact audit")
    print(f"universe=odd q in [3,{Q_MAX}] rows={len(augmented_rows)}")
    print(f"root_scan_candidates={root_scans} nonfixed_even_even_orbit_pairs={orbit_pairs}")
    print(f"divisor_star_terms={divisor_terms}")
    print(f"zero_S_count={len(zero_packets)} prime_zero_S={prime_zero_packets}")
    print(f"controls={json.dumps(controls, sort_keys=True, separators=(',', ':'))}")
    print("mod16_hostile=q13:S=0,Nmod4=3 versus q29:S=8,Nmod4=1")
    print(f"semantic_sha256={digest}")
    print("PASS: inverse orbits, mod-eight laws, degree residues, and zero gates")


if __name__ == "__main__":
    main()
