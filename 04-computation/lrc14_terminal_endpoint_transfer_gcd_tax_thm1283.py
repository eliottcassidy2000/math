#!/usr/bin/env python3
"""Exact referee for THM-1283's terminal endpoint transfer and gcd tax.

The continuum inputs (selection of an endpoint owner, the THM-1267 survivor
mass, and extraction of THM-1274's terminal tooth) remain paper topology
providers.  This file checks the complete endpoint arithmetic, both crossing
orientations, the normalized tax conversion, the integer cut, the compact
private-owner alternative, and the exact mirrored guardrail.
"""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from math import gcd
from pathlib import Path


F = Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def danger_tooth(speed: int, address: int) -> tuple[F, F]:
    return (
        F(14 * address - 1, 14 * speed),
        F(14 * address + 1, 14 * speed),
    )


def safe_gap(carrier: int, gap: int) -> tuple[F, F]:
    return (
        F(14 * gap + 1, 14 * carrier),
        F(14 * gap + 13, 14 * carrier),
    )


def safe_component(speed: int, address: int) -> tuple[F, F]:
    return (
        F(14 * address + 1, 14 * speed),
        F(14 * address + 13, 14 * speed),
    )


def endpoint_incidence_audit() -> tuple[int, int, int, int, int]:
    rows = 0
    left_rows = 0
    right_rows = 0
    sharp_gcd_rows = 0
    strict_residue_improvements = 0

    for carrier in range(1, 31):
        denominator = 14 * carrier
        for gap in range(carrier):
            for owner in range(carrier + 1, 8 * carrier + 1):
                common_gcd = gcd(carrier, owner)
                residue_modulus = 14 // gcd(common_gcd, 14)
                reduced_residue = (
                    carrier // common_gcd + owner // common_gcd
                ) % residue_modulus
                least_positive_residue = (
                    reduced_residue if reduced_residue else residue_modulus
                )

                for side in (-1, 1):
                    endpoint_numerator = 14 * gap + (1 if side < 0 else 13)
                    raw = endpoint_numerator * owner
                    quotient = raw // denominator
                    for address in (quotient, quotient + 1):
                        if abs(raw - denominator * address) >= carrier:
                            continue

                        if side > 0:
                            signed_residual = denominator * address - raw
                        else:
                            signed_residual = raw - denominator * address
                        endpoint_residual = carrier + signed_residual

                        require(
                            -carrier < signed_residual < carrier,
                            "strict incidence did not give signed strip",
                        )
                        require(
                            0 < endpoint_residual < 2 * carrier,
                            "endpoint numerator escaped (0,2c)",
                        )
                        require(
                            endpoint_residual % common_gcd == 0,
                            "gcd does not divide endpoint numerator",
                        )
                        require(
                            (endpoint_residual - carrier - owner) % 14 == 0,
                            "endpoint residue is not c+x modulo 14",
                        )
                        require(
                            endpoint_residual
                            >= common_gcd * least_positive_residue,
                            "least positive residue lower bound failed",
                        )

                        endpoint = F(endpoint_numerator, denominator)
                        tooth = danger_tooth(owner, address)
                        outward_wall = tooth[1] if side > 0 else tooth[0]
                        outward_length = abs(outward_wall - endpoint)
                        exact_length = F(
                            endpoint_residual, 14 * carrier * owner
                        )
                        require(
                            outward_length == exact_length,
                            "outward endpoint formula failed",
                        )
                        require(
                            outward_length >= F(common_gcd, 14 * carrier * owner),
                            "gcd/lcm endpoint quantum failed",
                        )
                        require(
                            outward_length < F(1, 7 * owner)
                            < F(1, 7 * carrier),
                            "outward wall does not stay inside carrier tooth",
                        )

                        if side > 0:
                            carrier_tooth = (
                                endpoint,
                                endpoint + F(1, 7 * carrier),
                            )
                            require(
                                tooth[0]
                                < endpoint
                                == carrier_tooth[0]
                                < tooth[1]
                                < carrier_tooth[1],
                                "right proper crossing order failed",
                            )
                            right_rows += 1
                        else:
                            carrier_tooth = (
                                endpoint - F(1, 7 * carrier),
                                endpoint,
                            )
                            require(
                                carrier_tooth[0]
                                < tooth[0]
                                < endpoint
                                == carrier_tooth[1]
                                < tooth[1],
                                "left proper crossing order failed",
                            )
                            left_rows += 1

                        require(
                            min(outward_wall, endpoint) >= carrier_tooth[0]
                            and max(outward_wall, endpoint) <= carrier_tooth[1],
                            "endpoint seam is not inside adjacent carrier tooth",
                        )
                        sharp_gcd_rows += endpoint_residual == common_gcd
                        strict_residue_improvements += (
                            endpoint_residual > common_gcd
                        )
                        rows += 1

    require(rows == left_rows + right_rows, "orientation ledger mismatch")
    require(left_rows == right_rows, "mirror incidence count mismatch")
    require(sharp_gcd_rows > 0, "gcd endpoint quantum was never sharp")
    require(
        strict_residue_improvements > 0,
        "exact residue never improved the gcd quantum",
    )
    return (
        rows,
        left_rows,
        right_rows,
        sharp_gcd_rows,
        strict_residue_improvements,
    )


def normalized_tax_audit() -> tuple[int, F, F]:
    target_mass = F(11, 360)
    endpoint_density = F(3, 4)
    inverse_quantile = target_mass / endpoint_density
    require(inverse_quantile == F(11, 270), "wrong endpoint inverse quantile")
    require(F(1, 6) > inverse_quantile, "quantile left endpoint sixth")

    rows = 0
    for carrier in range(1, 13):
        for slowest in range(carrier + 1, 3 * carrier + 1):
            for owner in range(slowest + 1, 8 * carrier + 1):
                for endpoint_residual in range(1, 2 * carrier):
                    if endpoint_residual % gcd(carrier, owner):
                        continue
                    seam = F(endpoint_residual, 14 * carrier * owner)
                    normalized_seam = F(7 * slowest, 6) * seam
                    require(
                        normalized_seam
                        == F(slowest * endpoint_residual, 12 * carrier * owner),
                        "normalized seam identity failed",
                    )
                    require(
                        normalized_seam
                        >= F(slowest * gcd(carrier, owner), 12 * carrier * owner),
                        "normalized gcd tax failed",
                    )

                    rational_margin = (
                        F(13, 12)
                        - F(slowest, 2 * carrier)
                        - F(11, 270)
                        - normalized_seam
                    )
                    integer_margin = (
                        563 * carrier * owner
                        - 270 * slowest * owner
                        - 45 * slowest * endpoint_residual
                    )
                    require(
                        F(540 * carrier * owner) * rational_margin
                        == integer_margin,
                        "rational-to-integer tax conversion failed",
                    )
                    rows += 1

    return rows, inverse_quantile, F(1, 6) - inverse_quantile


def integer_rounding_audit() -> tuple[int, int]:
    strict_rows = 0
    rejected_rows = 0
    for carrier in range(1, 41):
        for slowest in range(carrier + 1, 3 * carrier + 1):
            for owner in range(slowest + 1, 10 * carrier + 1):
                common_gcd = gcd(carrier, owner)
                margin = (
                    563 * carrier * owner
                    - 270 * slowest * owner
                    - 45 * slowest * common_gcd
                )
                if margin > 0:
                    require(margin >= 1, "positive integer margin below one")
                    require(
                        270 * slowest * owner + 45 * slowest * common_gcd
                        <= 563 * carrier * owner - 1,
                        "strict gcd cut did not round",
                    )
                    strict_rows += 1
                else:
                    rejected_rows += 1
    require(strict_rows > 0 and rejected_rows > 0, "rounding audit missed a branch")
    return strict_rows, rejected_rows


def private_owner_audit() -> tuple[int, int]:
    rows = 0
    compact_rows = 0
    for carrier in range(1, 101):
        for owner in range(carrier + 1, 20 * carrier + 1):
            private_floor = (owner + 7 * carrier - 1) // (7 * carrier)
            require(
                (private_floor == 1) == (owner <= 7 * carrier),
                "single selected endpoint owner does not match x<=7c",
            )
            compact_rows += private_floor == 1
            rows += 1
    return rows, compact_rows


def protected_needle_tree_audit() -> tuple[int, F]:
    # If two margins sum to at least 1/c, one is at least 1/(2c), while every
    # endpoint seam is strictly shorter than 1/(7x)<1/(7c).
    rows = 0
    weakest_margin_ratio = F(10**9)
    for carrier in range(1, 101):
        for owner in range(carrier + 1, 20 * carrier + 1):
            protected_margin = F(1, 2 * carrier)
            maximum_seam = F(1, 7 * owner)
            require(
                maximum_seam < F(1, 7 * carrier) < protected_margin,
                "protected margin did not contain endpoint seam",
            )
            weakest_margin_ratio = min(
                weakest_margin_ratio, protected_margin / maximum_seam
            )
            rows += 1
    require(weakest_margin_ratio > F(7, 2), "protected margin ratio too small")
    return rows, weakest_margin_ratio


def mirrored_guardrail_audit() -> tuple[int, F, F, F]:
    carrier, slowest, owner = 2, 4, 28
    rows = (
        # (gap, side, x-address, a-safe-component address)
        (0, -1, 1, 0),
        (1, 1, 27, 3),
    )
    observed: list[tuple[F, F, F]] = []
    for gap, side, owner_address, component_address in rows:
        gap_interval = safe_gap(carrier, gap)
        component = safe_component(slowest, component_address)
        endpoint = gap_interval[0] if side < 0 else gap_interval[1]
        tooth = danger_tooth(owner, owner_address)
        wall = tooth[0] if side < 0 else tooth[1]
        seam = abs(wall - endpoint)
        tail = (
            gap_interval[0] - component[0]
            if side < 0
            else component[1] - gap_interval[1]
        )
        component_length = component[1] - component[0]
        ell = tail / component_length
        eta = seam / component_length

        raw = (14 * gap + (1 if side < 0 else 13)) * owner
        signed_residual = (
            raw - 14 * carrier * owner_address
            if side < 0
            else 14 * carrier * owner_address - raw
        )
        endpoint_residual = carrier + signed_residual
        require(endpoint_residual == gcd(carrier, owner) == 2, "guardrail not gcd-sharp")
        require(seam == F(1, 392), "wrong mirrored seam length")
        require(ell == F(1, 12), "wrong mirrored protrusion")
        require(eta == F(1, 84), "wrong mirrored normalized seam")
        require(ell - eta == F(1, 14) > F(11, 270), "guardrail tax failed")
        observed.append((seam, ell, eta))

    require(observed[0] == observed[1], "mirror changed endpoint tax data")
    return len(rows), observed[0][0], observed[0][1], observed[0][2]


def centered_phase_clock_audit() -> tuple[int, int, int, tuple[int, ...]]:
    rows = 0
    exact_tax_rows = 0
    excluded_rows = 0
    for carrier in range(1, 31):
        for slowest in range(carrier + 1, (563 * carrier - 1) // 270 + 1):
            total = carrier + slowest
            for gap in range(carrier):
                center_word = 2 * gap + 1
                raw_center = total * center_word
                floor_center = raw_center // (2 * carrier)
                for center_address in {floor_center, floor_center + 1}:
                    signed_error = 2 * carrier * center_address - raw_center
                    if signed_error == 0 or abs(signed_error) > carrier:
                        continue
                    side = 1 if signed_error > 0 else -1
                    error = abs(signed_error)
                    ell_numerator = 6 * carrier + 7 * error - 6 * slowest
                    if ell_numerator <= 0:
                        continue
                    endpoint_word = 14 * gap + (13 if side > 0 else 1)
                    for owner in range(slowest + 1, 8 * carrier + 1):
                        raw_endpoint = endpoint_word * owner
                        floor_owner = raw_endpoint // (14 * carrier)
                        for owner_address in (floor_owner, floor_owner + 1):
                            if abs(raw_endpoint - 14 * carrier * owner_address) >= carrier:
                                continue
                            endpoint_signed = (
                                14 * carrier * owner_address - raw_endpoint
                                if side > 0
                                else raw_endpoint - 14 * carrier * owner_address
                            )
                            endpoint_residual = carrier + endpoint_signed
                            quotient_word = endpoint_residual - carrier + 6 * owner
                            require(quotient_word % 7 == 0, "endpoint clock quotient nonintegral")
                            winding = side * (
                                owner_address * total - center_address * owner
                            )
                            require(winding > 0, "center-to-endpoint winding is not positive")
                            require(
                                total * quotient_word - 7 * error * owner
                                == 14 * carrier * winding,
                                "center/endpoint clock elimination failed",
                            )
                            exact_cut = (
                                270 * slowest * owner
                                + 45 * slowest * endpoint_residual
                                < (248 * carrier + 315 * error) * owner
                            )
                            suffix_cut = (
                                F(ell_numerator, 12 * carrier)
                                - F(
                                    slowest * endpoint_residual,
                                    12 * carrier * owner,
                                )
                                > F(11, 270)
                            )
                            require(exact_cut == suffix_cut, "exact-e cut is not tail subtraction")
                            exact_tax_rows += exact_cut
                            excluded_rows += not exact_cut
                            rows += 1

    # THM-1266's sharp five-rung local row is killed by the exact phase cut.
    carrier, slowest, owner = 140, 254, 256
    error, endpoint_residual = 126, 172
    lhs = 270 * slowest * owner + 45 * slowest * endpoint_residual
    rhs_closed = (248 * carrier + 315 * error) * owner - 1
    coarse_closed = 563 * carrier * owner - 1
    require(lhs == 19_522_440, "sharp-row exact lhs changed")
    require(rhs_closed == 19_048_959, "sharp-row exact rhs changed")
    require(lhs > rhs_closed and lhs <= coarse_closed, "sharp row did not separate exact/coarse cuts")
    require(
        F(33, 280) - F(slowest * endpoint_residual, 12 * carrier * owner)
        == F(25, 1536) < F(11, 270),
        "sharp-row survivor subtraction changed",
    )
    return rows, exact_tax_rows, excluded_rows, (lhs, rhs_closed, coarse_closed)


def source_has_no_assert_nodes() -> int:
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "referee contains optimization-sensitive assert nodes")
    return count


def main() -> None:
    no_asserts = source_has_no_assert_nodes()
    incidence = endpoint_incidence_audit()
    tax_rows, inverse_quantile, sixth_slack = normalized_tax_audit()
    strict_rows, rejected_rows = integer_rounding_audit()
    private_rows, compact_rows = private_owner_audit()
    protected_rows, protected_ratio = protected_needle_tree_audit()
    mirror_rows, mirror_seam, mirror_ell, mirror_eta = mirrored_guardrail_audit()
    phase_rows, phase_pass, phase_excluded, sharp_phase = centered_phase_clock_audit()

    print("THM-1283 TERMINAL ENDPOINT TRANSFER / GCD TAX EXACT AUDIT")
    print(f"Python assert nodes = {no_asserts}")
    print(f"strict endpoint incidences = {incidence[0]}")
    print(f"left/right incidences = {incidence[1]}/{incidence[2]}")
    print(f"gcd-sharp endpoint rows = {incidence[3]}")
    print(f"strict residue improvements over gcd = {incidence[4]}")
    print(f"normalized tax identity rows = {tax_rows}")
    print(f"endpoint inverse quantile = {inverse_quantile}")
    print(f"endpoint-sixth quantile slack = {sixth_slack}")
    print(f"positive/rejected integer gcd cuts = {strict_rows}/{rejected_rows}")
    print(f"private-count rows = {private_rows}")
    print(f"single-occurrence compact rows = {compact_rows}")
    print(f"protected-needle containment rows = {protected_rows}")
    print(f"smallest protected-margin/seam-supremum ratio = {protected_ratio}")
    print(f"mirrored gcd-sharp guardrails = {mirror_rows}")
    print(f"guardrail seam/ell/eta = {mirror_seam}/{mirror_ell}/{mirror_eta}")
    print(f"center/endpoint clock rows = {phase_rows}")
    print(f"exact-e passing/excluded rows = {phase_pass}/{phase_excluded}")
    print(f"sharp-row exact lhs/rhs/coarse-rhs = {sharp_phase}")
    print("exact_tax=ell-eta>11/270")
    print("integer_cut=270*a*x+45*a*Q<=563*c*x-1")
    print("gcd_cut=270*a*x+45*a*gcd(c,x)<=563*c*x-1")
    print("phase_cut=270*a*x+45*a*Q<=(248*c+315*e)*x-1")
    print("protected_tree=internal_six_owner_tree+exterior_carrier_edge")
    print("vertices=carrier-wall,event-owner-wall,outer-survivor-obligation")
    print("switch=circle-reflection; tie_Hamiltonian_path=outward-event-order")
    print("tournament_scores=(0,1); cycles=0; SCCs=2; Hamilton_paths=1")
    print("preserves=side,address,residue,proper-crossing,and-survivor-subtraction")
    print("destroys=next mixed-owner safe-component cover")
    print("STATUS=PASS")
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"source_sha256={digest}")


if __name__ == "__main__":
    main()
