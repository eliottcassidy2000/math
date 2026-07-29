#!/usr/bin/env python3
"""Exact companion for THM-2863 endpoint Prony / carry-character splitting."""

from __future__ import annotations

from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

N = 2366
P = 13
SIGMA = 1275
PRIME = 4733
XI = 25


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    ROOT / "05-knowledge/results/lrc14_q3_q11_transverse_endpoint_horn_thm2847.out":
        "155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d",
    COMP / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    ROOT / "05-knowledge/results/lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out":
        "424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac",
    COMP / "lrc14_endpoint_galois_carry_torsor_thm2857.py":
        "0bae59c9b1460f37e1879a81746154593cb0699ee13b3e5e800ba0af95ea5e4c",
    ROOT / "05-knowledge/results/lrc14_endpoint_galois_carry_torsor_thm2857.out":
        "ac1194c46db2cdf43c807ece781b63971c081cc5f9070964007fdecdc20f1583",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_q3_q11_transverse_endpoint_horn_thm2847 as horn


def transform(values, omega):
    return tuple(
        sum(
            values[r] * pow(omega, (-frequency * r) % P, PRIME)
            for r in range(P)
        )
        % PRIME
        for frequency in range(P)
    )


def order(value):
    product = 1
    for exponent in range(1, N + 1):
        product = product * value % PRIME
        if product == 1:
            return exponent
    raise RuntimeError("multiplicative order exceeded 2366")


def main():
    # Replay the literal endpoint interval at four 91-unit multiplier samples.
    allocation = horn.allocation
    (_module, full, _details, _e3, _clocks, _q_pairs) = (
        allocation.build_geometry()
    )
    period = full.T
    unit = period // P
    target = tuple(
        endpoint + allocation.physical.SHIFT
        for endpoint in allocation.ATOM_INTERVAL
    )
    target_atom = ((*target, 1),)
    endpoint_base = allocation.endpoint_base
    endpoint = endpoint_base.endpoint
    origins = (
        allocation.RIGHT_ORIGIN,
        allocation.add(allocation.RIGHT_ORIGIN, allocation.TARGET_STEP),
    )
    multipliers = tuple(12 + 26 * m for m in range(1, 5))
    require(
        multipliers == (38, 64, 90, 116)
        and all(gcd(multiplier, 91) == 1 for multiplier in multipliers),
        "eligible multiplier progression changed",
    )

    bank = {}
    exponent_bank = {}
    reduction_divisor = endpoint.NN // N
    require(reduction_divisor * N == endpoint.NN, "2366 reduction failed")
    for address in origins:
        ell = endpoint_base.REPS[address]
        present = tuple(endpoint.build_set(endpoint_base.PAT_E3, ell))
        starts = tuple(left for left, _right in present)
        for q in (3, 11):
            shifted = allocation.physical.overlap.shift_weighted(
                target_atom, q * unit
            )
            restricted = allocation.indexed_weighted_intersection(
                shifted, present, starts
            )
            values = []
            exponents = []
            for multiplier in multipliers:
                values.append(
                    tuple(
                        allocation.endpoint_sum(
                            restricted, -multiplier, embedding
                        )
                        for embedding in endpoint.MODS
                    )
                )
                if restricted:
                    raw = tuple(
                        (
                            multiplier * endpoint.RDIL * point
                        ) % endpoint.NN
                        for point in restricted[0][:2]
                    )
                    require(
                        all(value % reduction_divisor == 0 for value in raw),
                        "endpoint exponent lost its 2366 reduction",
                    )
                    exponents.append(
                        tuple(value // reduction_divisor for value in raw)
                    )
            bank[(address, q)] = tuple(values)
            exponent_bank[(address, q)] = tuple(exponents)

    expected_values = (
        (231164267889491750, 630230755085920022),
        (209019417557558827, 920317591576844127),
        (282053890632156827, 349031617258410408),
        (268435097921701104, 322681190093619272),
    )
    expected_exponents = tuple(
        ((2190 + 13 * m) % N, (2262 + 169 * m) % N)
        for m in range(1, 5)
    )
    require(
        expected_exponents
        == ((2203, 65), (2216, 234), (2229, 403), (2242, 572)),
        "displayed endpoint exponent progression changed",
    )
    require(
        bank[(origins[0], 3)]
        == bank[(origins[0], 11)]
        == bank[(origins[1], 11)]
        == expected_values
        and bank[(origins[1], 3)] == ((0, 0),) * 4,
        "four-sample physical endpoint bank changed",
    )
    require(
        exponent_bank[(origins[0], 3)]
        == exponent_bank[(origins[0], 11)]
        == exponent_bank[(origins[1], 11)]
        == expected_exponents
        and exponent_bank[(origins[1], 3)] == (),
        "four-sample endpoint exponent bank changed",
    )
    require(
        all(left and right for left, right in expected_values),
        "an endpoint factor vanished in a certified field",
    )

    # Independent exact finite-field certificate for the two-node sequence.
    require(
        pow(XI, N, PRIME) == 1
        and all(pow(XI, N // prime_factor, PRIME) != 1
                for prime_factor in (2, 7, 13)),
        "XI lacks exact order 2366",
    )
    omega = pow(XI, 182, PRIME)
    require(
        pow(omega, P, PRIME) == 1 and omega != 1,
        "omega lacks exact order 13",
    )
    alpha_l = pow(XI, 2190, PRIME)
    alpha_r = pow(XI, 2262, PRIME)
    lambda_l = pow(XI, 13, PRIME)
    lambda_r = pow(XI, 169, PRIME)
    require(
        lambda_l != lambda_r
        and order(lambda_l) == 182
        and order(lambda_r) == 14,
        "oriented Prony nodes changed",
    )

    def coefficient(m):
        return (
            alpha_l * pow(lambda_l, m, PRIME)
            - alpha_r * pow(lambda_r, m, PRIME)
        ) % PRIME

    coefficients = tuple(coefficient(m) for m in range(6))
    node_sum = (lambda_l + lambda_r) % PRIME
    node_product = lambda_l * lambda_r % PRIME
    require(
        all(
            coefficients[m + 2]
            == (
                node_sum * coefficients[m + 1]
                - node_product * coefficients[m]
            ) % PRIME
            for m in range(4)
        ),
        "two-node recurrence failed",
    )

    c1, c2, c3, c4 = coefficients[1:5]
    determinant = (c1 * c3 - c2 * c2) % PRIME
    expected_determinant = (
        -alpha_l
        * alpha_r
        * lambda_l
        * lambda_r
        * pow((lambda_l - lambda_r) % PRIME, 2, PRIME)
    ) % PRIME
    require(
        determinant == expected_determinant != 0,
        "Prony determinant identity failed",
    )
    determinant_inverse = pow(determinant, -1, PRIME)
    recovered_sum = (c1 * c4 - c2 * c3) * determinant_inverse % PRIME
    recovered_product = (c2 * c4 - c3 * c3) * determinant_inverse % PRIME
    require(
        (recovered_sum, recovered_product) == (node_sum, node_product),
        "four-sample Prony recovery failed",
    )
    require(
        {
            value
            for value in range(PRIME)
            if (
                value * value
                - recovered_sum * value
                + recovered_product
            ) % PRIME == 0
        }
        == {lambda_l, lambda_r},
        "recovered node polynomial changed",
    )

    node_difference_inverse = pow(
        (lambda_l - lambda_r) % PRIME, -1, PRIME
    )
    for j in range(4):
        recovered_l = (
            coefficients[j + 1] - lambda_r * coefficients[j]
        ) * node_difference_inverse % PRIME
        recovered_minus_r = (
            lambda_l * coefficients[j] - coefficients[j + 1]
        ) * node_difference_inverse % PRIME
        require(
            recovered_l
            == alpha_l * pow(lambda_l, j, PRIME) % PRIME
            and recovered_minus_r
            == -alpha_r * pow(lambda_r, j, PRIME) % PRIME,
            "oriented endpoint splitting failed",
        )

    # Relative Galois action: fixed right node/weight plus character-three left.
    require(
        SIGMA == 1 + 14 * 91
        and SIGMA * 13 % N == 13
        and SIGMA * 169 % N == 169
        and SIGMA * 2262 % N == 2262
        and (SIGMA * 2190 - 2190) % N == 3 * 182,
        "relative Galois character changed",
    )
    l1 = alpha_l * lambda_l % PRIME
    r1 = alpha_r * lambda_r % PRIME
    orbit = tuple(
        (pow(omega, 3 * r, PRIME) * l1 - r1) % PRIME
        for r in range(P)
    )
    orbit_transform = transform(orbit, omega)
    require(
        tuple(k for k, value in enumerate(orbit_transform) if value)
        == (0, 3)
        and orbit_transform[0] == -P * r1 % PRIME
        and orbit_transform[3] == P * l1 % PRIME,
        "trivial plus character-three endpoint splitting failed",
    )

    # Match the canonical character-three projection of the carry derivative.
    carry = tuple(
        449 * ((1 if r == 1 else 0) - (1 if r == 0 else 0))
        for r in range(P)
    )
    carry_transform = transform(carry, omega)
    carry_character = 449 * (pow(omega, -3, PRIME) - 1) % PRIME
    endpoint_character = P * l1 % PRIME
    require(
        carry_transform[3] == carry_character != 0
        and endpoint_character != 0,
        "character-three generator vanished",
    )
    intertwiner = endpoint_character * pow(carry_character, -1, PRIME) % PRIME
    inverse_p = pow(P, -1, PRIME)
    projected_carry = tuple(
        carry_character
        * inverse_p
        * pow(omega, 3 * r, PRIME)
        % PRIME
        for r in range(P)
    )
    centered_endpoint = tuple(
        pow(omega, 3 * r, PRIME) * l1 % PRIME
        for r in range(P)
    )
    require(
        tuple(intertwiner * value % PRIME for value in projected_carry)
        == centered_endpoint,
        "normalized character-three intertwiner failed",
    )

    # Node-only data are insufficient: swapping endpoint weights preserves
    # the recurrence polynomial while attaching the charged weight elsewhere.
    swapped = tuple(
        (
            alpha_r * pow(lambda_l, m, PRIME)
            - alpha_l * pow(lambda_r, m, PRIME)
        ) % PRIME
        for m in range(6)
    )
    require(
        swapped != coefficients
        and all(
            swapped[m + 2]
            == (node_sum * swapped[m + 1] - node_product * swapped[m]) % PRIME
            for m in range(4)
        ),
        "node-only orientation hostile failed",
    )

    print("ENDPOINT PRONY SPLITTER / CARRY CHARACTER-THREE INTERTWINER")
    print("multipliers=(38,64,90,116); all_gcd_91=1")
    print(
        "physical_rows=q3_origin00=q11_origin00=q11_origin12;"
        "q3_origin12=empty"
    )
    print(f"reduced_exponents={expected_exponents}; dual_field_nonzero=8/8")
    print(
        "two_nodes=(orders=182,14); recurrence_order=2;"
        f" Prony_determinant={determinant}"
    )
    print(f"first_six_specialized_coefficients={coefficients}")
    print(
        "relative_module=trivial+character3; "
        f"carry_character3={carry_character};"
        f" endpoint_character3={endpoint_character};"
        f" intertwiner={intertwiner}"
    )
    print("node_only_swapped_weight_hostile=PASS")
    print(
        "scope=endpoint factors and normalized linear interface only;"
        " common full-current bank/E3 transport OPEN"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
