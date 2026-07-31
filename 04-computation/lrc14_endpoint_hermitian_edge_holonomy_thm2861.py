#!/usr/bin/env python3
"""Exact companion for THM-2861 endpoint Hermitian edge holonomy."""

from __future__ import annotations

from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

P = 13
N = 1183
STEP = 91


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
    COMP / "lrc14_endpoint_galois_carry_torsor_thm2857.py":
        "0bae59c9b1460f37e1879a81746154593cb0699ee13b3e5e800ba0af95ea5e4c",
    ROOT / "05-knowledge/results/lrc14_endpoint_galois_carry_torsor_thm2857.out":
        "ac1194c46db2cdf43c807ece781b63971c081cc5f9070964007fdecdc20f1583",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_q3_q11_transverse_endpoint_horn_thm2847 as horn


allocation = horn.allocation


def sparse_product(left, right):
    product = {}
    for exponent_left, coefficient_left in left.items():
        for exponent_right, coefficient_right in right.items():
            exponent = (exponent_left + exponent_right) % N
            product[exponent] = (
                product.get(exponent, 0)
                + coefficient_left * coefficient_right
            )
    return {exponent: value for exponent, value in product.items() if value}


def conjugate_terms(terms):
    return {(-exponent) % N: coefficient for exponent, coefficient in terms.items()}


def shift_terms(terms, exponent):
    return {(key + exponent) % N: value for key, value in terms.items()}


def transform(values, omega, modulus):
    return tuple(
        sum(
            values[r] * pow(omega, (-frequency * r) % P, modulus)
            for r in range(P)
        )
        % modulus
        for frequency in range(P)
    )


def main():
    # Replay the actual THM-2847 horn and its q7 failure.
    (_module, full, _details, e3, clocks, _q_pairs) = allocation.build_geometry()
    period = full.T
    unit = period // P
    atom = allocation.ATOM_INTERVAL
    target = tuple(endpoint + allocation.physical.SHIFT for endpoint in atom)

    def signature(interval, s, t, clock):
        return (
            allocation.contained(interval, e3),
            allocation.contained(interval, clocks[clock]),
            horn.safe_comb_contains(
                interval, full, full.W[1], 182,
                -14 * s - 13, -14 * s + 13,
            ),
            horn.safe_comb_contains(
                interval, full, full.W[2], 182,
                -14 * t - 13, -14 * t + 13,
            ),
            horn.safe_comb_contains(
                interval, full, full.C2, 182,
                14 * s - 13, 14 * s + 13,
            ),
            horn.safe_comb_contains(
                interval, full, full.C3, 182,
                14 * t - 13, 14 * t + 13,
            ),
        )

    def full_cells(q):
        shifted = horn.circular_shift_interval(target, q * unit, period)
        return tuple(
            (s, t, clock)
            for s in allocation.COMMON_S
            for t in allocation.COMMON_T
            for clock in range(7)
            if (
                all(signature(atom, s, t, clock))
                and all(signature(target, s, t, clock))
                and all(signature(shifted, s, t, clock))
            )
        )

    common = tuple(sorted(set(full_cells(3)) & set(full_cells(11))))
    shifted_q7 = horn.circular_shift_interval(target, 7 * unit, period)
    q7_signatures = tuple(
        (cell, signature(shifted_q7, *cell)) for cell in common
    )
    q7_e3_only = tuple(
        cell
        for cell, bits in q7_signatures
        if not bits[0] and all(bits[1:])
    )
    q7_extra = tuple(
        cell
        for cell, bits in q7_signatures
        if not bits[0] and not all(bits[1:])
    )
    require(
        len(common) == 42
        and len(q7_e3_only) == 20
        and len(q7_extra) == 22,
        "q3/q11/q7 horn split changed",
    )

    endpoint_base = allocation.endpoint_base
    endpoint = endpoint_base.endpoint
    origins = (
        allocation.RIGHT_ORIGIN,
        allocation.add(allocation.RIGHT_ORIGIN, allocation.TARGET_STEP),
    )
    target_atom = ((*target, 1),)
    rows = {}
    reductions = set()
    for address in origins:
        ell = endpoint_base.REPS[address]
        present = tuple(endpoint.build_set(endpoint_base.PAT_E3, ell))
        starts = tuple(left for left, _right in present)
        for q in (3, 11, 7):
            shifted = allocation.physical.overlap.shift_weighted(
                target_atom, q * unit
            )
            restricted = allocation.indexed_weighted_intersection(
                shifted, present, starts
            )
            values = tuple(
                allocation.endpoint_sum(
                    restricted, -endpoint_base.Y0, embedding
                )
                for embedding in endpoint.MODS
            )
            rows[(address, q)] = (len(restricted), values)
            if restricted:
                exponents = tuple(
                    (
                        endpoint_base.Y0 * endpoint.RDIL * point
                    ) % endpoint.NN
                    for point in restricted[0][:2]
                )
                divisor = gcd(endpoint.NN, gcd(*exponents))
                reductions.add(
                    (
                        endpoint.NN // divisor,
                        tuple(value // divisor for value in exponents),
                    )
                )

    endpoint_value = (231164267889491750, 630230755085920022)
    require(
        rows[(origins[0], 3)] == (1, endpoint_value)
        and rows[(origins[0], 11)] == (1, endpoint_value)
        and rows[(origins[1], 3)] == (0, (0, 0))
        and rows[(origins[1], 11)] == (1, endpoint_value)
        and rows[(origins[0], 7)] == rows[(origins[1], 7)] == (0, (0, 0))
        and reductions == {(2366, (2203, 65))},
        "physical endpoint rows changed",
    )
    allocation_phase = (
        endpoint_base.Y0 * endpoint.RDIL * unit
    ) % endpoint.NN
    ancestry_phase = (
        endpoint_base.Y0 * endpoint.RDIL * period
    ) % endpoint.NN
    require(
        allocation_phase == ancestry_phase == 0,
        "canonical endpoint phase is no longer scalar-trivial",
    )

    # Exact Hermitian Galois edge in the group algebra of mu_1183.
    omega_exponent = STEP * 3

    def c_terms(r):
        return {
            624: 1,
            (510 + omega_exponent * r) % N: -1,
        }

    edge_terms = []
    reverse_terms = []
    for r in range(P):
        edge = sparse_product(c_terms((r + 1) % P), conjugate_terms(c_terms(r)))
        reverse = sparse_product(
            conjugate_terms(c_terms((r + 1) % P)), c_terms(r)
        )
        require(
            edge == shift_terms(reverse, omega_exponent),
            "Hermitian phase law failed in the exact group ring",
        )
        edge_terms.append(edge)
        reverse_terms.append(reverse)

    # Independent exact finite-field certificate: conjugation is zeta->zeta^-1.
    prime = 4733
    zeta = 625
    omega = pow(zeta, STEP, prime)
    require(
        pow(zeta, N, prime) == 1
        and pow(zeta, N // 7, prime) != 1
        and pow(zeta, N // 13, prime) != 1,
        "finite-field element lacks exact order 1183",
    )

    def evaluate(terms):
        return sum(
            coefficient * pow(zeta, exponent, prime)
            for exponent, coefficient in terms.items()
        ) % prime

    c_values = tuple(evaluate(c_terms(r)) for r in range(P))
    conjugates = tuple(evaluate(conjugate_terms(c_terms(r))) for r in range(P))
    edges = tuple(
        c_values[(r + 1) % P] * conjugates[r] % prime
        for r in range(P)
    )
    reverse_edges = tuple(
        conjugates[(r + 1) % P] * c_values[r] % prime
        for r in range(P)
    )
    antisymmetric = tuple(
        (left - right) % prime
        for left, right in zip(edges, reverse_edges)
    )
    edge_support = tuple(
        frequency
        for frequency, value in enumerate(transform(edges, omega, prime))
        if value
    )
    antisymmetric_support = tuple(
        frequency
        for frequency, value in enumerate(
            transform(antisymmetric, omega, prime)
        )
        if value
    )
    require(
        all(
            edges[r] == pow(omega, 3, prime) * reverse_edges[r] % prime
            for r in range(P)
        )
        and all(antisymmetric)
        and len(set(edges)) == P
        and edge_support == (0, 3, 10)
        and antisymmetric_support == edge_support,
        "finite-field Hermitian edge certificate failed",
    )

    # Common complex phase cancels; a constant scalar-linear carrier is sharp.
    rational_scale = 17
    scaled_edges = tuple(
        rational_scale**2 * value % prime for value in edges
    )
    require(
        scaled_edges
        == tuple(
            (rational_scale * c_values[(r + 1) % P])
            * (rational_scale * conjugates[r])
            % prime
            for r in range(P)
        ),
        "common-phase gauge cancellation failed",
    )
    constant_edges = (c_values[0] * conjugates[0] % prime,) * P
    constant_reverse = constant_edges
    require(
        all(left == right for left, right in zip(constant_edges, constant_reverse))
        and tuple(
            frequency
            for frequency, value in enumerate(
                transform(constant_edges, omega, prime)
            )
            if value
        )
        == (0,),
        "scalar-linear hostile changed",
    )

    print("ENDPOINT HERMITIAN EDGE HOLONOMY")
    print(
        f"physical_horn=(common={len(common)},"
        f"q7_E3_only={len(q7_e3_only)},q7_extra={len(q7_extra)})"
    )
    print(
        "endpoint_rows=(origin00:q3=c,q11=c,q7=0;"
        "origin12:q3=0,q11=c,q7=0)"
    )
    print(
        "endpoint_reduction=(conductor=2366,exponents=(2203,65),"
        "zeta1183=(624,510)); physical_phase_increments=(0,0)"
    )
    print(
        "canonical_carry=K0_linear_and_scalar_trivial; "
        "physical_q7_endpoint=zero_after_old_block_projection"
    )
    print(
        "Hermitian_law=c_(r+1)*bar(c_r)="
        "omega^3*bar(c_(r+1))*c_r; antisymmetric_nonzero=13/13"
    )
    print(
        f"Hermitian_edge_orbit_size={len(set(edges))}; "
        f"Fourier_support={edge_support}; "
        f"antisymmetric_support={antisymmetric_support}; bidegree=(1,1)"
    )
    print(
        "common_phase_gauge_invariant=1; orientation_reversal=conjugation; "
        "constant_scalar_T_hostile=(antisymmetric=0,support=(0,))"
    )
    print(
        f"exact_specialization=(prime={prime},zeta1183={zeta}); "
        "missing=lawful_adjacent_section_cosupport_plus_E3_transport"
    )
    print(
        "scope=coefficient_Hermitian_detector_and_physical_no-go; "
        "no_pair_twist_service_or_LRC14_conclusion=1"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
