#!/usr/bin/env python3
"""Fraction-exact replay for THM-836's two-sheet owner-shell obstruction.

The proof starts at the signed q=13 wall for the common exception pair
(13d,5d).  After multiplying by p=(5d)^(-1) mod 13, the ten quotient
speeds have the balanced residues +/-2,...,+/-6.  This script audits the
two directional first-exit coefficients of the closed 1/11-deep component,
the closed central return interval, the thin-shell owner inequality, and the
resulting modular packing condition.

All theorem-facing decisions use integers or fractions.Fraction.  No circle
sampling, floating point, optimizer, or external theorem script is used.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd


BETA = F(1, 11)
GAMMA = F(2, 143)
H_RADIUS = F(2, 169)
THRESHOLD = F(11, 13)
SIGNED_RESIDUES = (-6, -5, -4, -3, -2, 2, 3, 4, 5, 6)
EXPECTED_COEFFICIENTS = (9, 20, 31, 42, 53, 64, 75, 86, 97, 108)


def fmt(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def q_value(value: F) -> F:
    return norm(9 * value) + norm(4 * value)


def balanced_residue(value: int, modulus: int = 13) -> int:
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


def coefficient(direction: int, residue: int) -> int:
    """First 1/11-exit coefficient in the chosen time direction.

    ``direction*residue < 0`` moves the signed phase toward zero.  Otherwise
    it first moves away from zero and exits beside the next integer.
    """

    assert direction in (-1, 1)
    assert residue in SIGNED_RESIDUES
    magnitude = abs(residue)
    if direction * residue < 0:
        return 11 * magnitude - 13
    return 130 - 11 * magnitude


def coefficient_audit() -> tuple[object, ...]:
    rows = []
    endpoint_tests = 0
    eta = F(1, 143_000)
    for residue in SIGNED_RESIDUES:
        for direction in (-1, 1):
            value = coefficient(direction, residue)
            start = F(residue, 13)
            endpoint = start + direction * F(value, 143)
            before = endpoint - direction * eta
            after = endpoint + direction * eta
            assert norm(start) >= F(2, 13) > BETA
            assert norm(endpoint) == BETA
            assert norm(before) > BETA
            assert norm(after) < BETA
            rows.append((residue, direction, value))
            endpoint_tests += 4

        assert coefficient(1, residue) + coefficient(-1, residue) == 117

    for direction in (-1, 1):
        observed = tuple(sorted(coefficient(direction, residue) for residue in SIGNED_RESIDUES))
        assert observed == EXPECTED_COEFFICIENTS
    return tuple(rows), endpoint_tests


def signed_row(d: int, seed: int) -> dict[int, int]:
    """A deterministic lift of the ten normalized signed residue owners."""

    assert d > 0 and gcd(d, 13) == 1
    inverse_p = (5 * d) % 13
    answer: dict[int, int] = {}
    for index, residue in enumerate(SIGNED_RESIDUES):
        raw = (residue * inverse_p) % 13
        assert raw != 0
        lift = (seed * (index + 3) + index * index + 2 * d) % 29
        answer[residue] = raw + 13 * lift
    assert len(set(answer.values())) == 10
    return answer


def local_geometry_audit() -> tuple[object, ...]:
    rows = []
    flip_histogram: Counter[int] = Counter()
    endpoint_tests = 0
    for d in range(1, 258, 2):
        if d % 13 == 0:
            continue
        seed = (17 * d + 5) % 31
        owners = signed_row(d, seed)
        speeds = tuple(sorted(owners.values()))
        B = max(speeds)
        p = pow((5 * d) % 13, -1, 13)
        t0 = F(p, 13)

        normalized = tuple(sorted(balanced_residue(p * speed) for speed in speeds))
        assert normalized == SIGNED_RESIDUES
        assert (d * p) % 13 == 8
        assert (d * t0) % 1 == F(8, 13)
        assert all(norm(speed * t0) >= F(2, 13) > BETA for speed in speeds)

        radii: dict[int, F] = {}
        owner_sets: dict[int, tuple[int, ...]] = {}
        for direction in (-1, 1):
            candidates = {
                residue: F(coefficient(direction, residue), 143 * speed)
                for residue, speed in owners.items()
            }
            radius = min(candidates.values())
            owner_set = tuple(sorted(residue for residue, value in candidates.items() if value == radius))
            radii[direction] = radius
            owner_sets[direction] = owner_set

            endpoint = t0 + direction * radius
            eta = F(1, 1001 * B)
            assert radius > eta
            before = t0 + direction * (radius - eta)
            after = t0 + direction * (radius + eta)
            at_values = tuple(norm(speed * endpoint) for speed in speeds)
            assert min(at_values) == BETA
            assert all(norm(speed * before) > BETA for speed in speeds)
            assert any(norm(speed * after) < BETA for speed in speeds)
            endpoint_tests += 3 * len(speeds)

        central = F(2, 143 * B)
        assert max(norm(speed * central) for speed in speeds) == GAMMA
        assert max(norm(-speed * central) for speed in speeds) == GAMMA
        endpoint_tests += 2 * len(speeds)

        target_left = t0 - H_RADIUS / d
        target_right = t0 + H_RADIUS / d
        assert q_value(d * target_left) == THRESHOLD
        assert q_value(d * target_right) == THRESHOLD
        outside_step = F(1, 169_000 * d)
        assert q_value(d * (target_left - outside_step)) < THRESHOLD
        assert q_value(d * (target_right + outside_step)) < THRESHOLD
        endpoint_tests += 4

        # Optional tournament telemetry.  Right-exit time is the observable;
        # left-exit time is the switch, and (residue,speed) is the tie path.
        right_order = tuple(sorted(SIGNED_RESIDUES, key=lambda z: (F(coefficient(1, z), owners[z]), z, owners[z])))
        left_order = tuple(sorted(SIGNED_RESIDUES, key=lambda z: (F(coefficient(-1, z), owners[z]), z, owners[z])))
        right_rank = {vertex: index for index, vertex in enumerate(right_order)}
        left_rank = {vertex: index for index, vertex in enumerate(left_order)}
        flips = sum(
            (right_rank[a] < right_rank[b]) != (left_rank[a] < left_rank[b])
            for index, a in enumerate(SIGNED_RESIDUES)
            for b in SIGNED_RESIDUES[index + 1 :]
        )
        flip_histogram[flips] += 1

        rows.append(
            (
                d,
                p,
                tuple((z, owners[z]) for z in SIGNED_RESIDUES),
                radii[-1],
                radii[1],
                owner_sets[-1],
                owner_sets[1],
                central,
                right_order,
                left_order,
                flips,
            )
        )
    return len(rows), endpoint_tests, tuple(sorted(flip_histogram.items())), sha256(repr(rows).encode()).hexdigest()


def least_absolute_residue(value: int, modulus: int = 13) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def largest_with_residue(B: int, residue: int) -> int:
    residue %= 13
    assert residue != 0
    value = B - ((B - residue) % 13)
    return value if value > 0 else 0


def shell_bank_audit(limit: int = 999) -> tuple[object, ...]:
    rows = 0
    owner_feasible = 0
    packing_failures = 0
    packing_equalities = 0
    owner_equalities = 0
    s1_rows = s1_owner_feasible = 0
    s3_rows = s3_owner_feasible = 0
    s3_mod_histogram: Counter[int] = Counter()
    s5_rows = s5_owner_feasible = 0
    first_s5_witness: tuple[int, ...] | None = None
    hasher = sha256()

    for d in range(1, limit + 1, 2):
        if d % 13 == 0:
            continue
        B_min = (13 * d + 1) // 2
        for B in range(B_min, 13 * d):
            rows += 1
            s = 2 * B - 13 * d
            assert s > 0 and s % 2 == 1
            denominator = 22 * B - 26 * d
            assert denominator == 117 * d + 11 * s > 0
            threshold_numerator = 117 * d * B

            residue_plus = (3 * d) % 13
            residue_minus = (-3 * d) % 13
            top_plus = largest_with_residue(B, residue_plus)
            top_minus = largest_with_residue(B, residue_minus)
            plus_ok = top_plus * denominator >= threshold_numerator
            minus_ok = top_minus * denominator >= threshold_numerator
            feasible = plus_ok and minus_ok

            delta = least_absolute_residue(6 * d)
            g_numerator = 11 * s * (13 * d + s)
            g_denominator = 2 * (117 * d + 11 * s)
            assert g_denominator > 0

            if s == 1:
                s1_rows += 1
                assert g_numerator < g_denominator
            if s == 3:
                s3_rows += 1
                assert g_numerator > g_denominator
                assert g_numerator < 2 * g_denominator
                if delta == 1:
                    assert d % 13 in (2, 11)
                    assert B % 13 == 8
                    assert {residue_plus, residue_minus} == {6, 7}
                    # The class-6 owner is at most B-2, but G<2 makes
                    # L=B-G strictly larger than B-2.
                    assert largest_with_residue(B, 6) == B - 2
                    assert (B - 2) * denominator < threshold_numerator
                assert not feasible
            if s == 5:
                s5_rows += 1

            if feasible:
                owner_feasible += 1
                owner_equalities += int(top_plus * denominator == threshold_numerator)
                owner_equalities += int(top_minus * denominator == threshold_numerator)
                gap = abs(top_plus - top_minus)
                assert gap > 0
                assert gap >= delta
                # Both owners lie in [L,B], so their gap is at most B-L=G.
                if gap * g_denominator > g_numerator:
                    packing_failures += 1
                assert gap * g_denominator <= g_numerator
                assert delta * g_denominator <= g_numerator
                packing_equalities += int(delta * g_denominator == g_numerator)

                if s == 1:
                    s1_owner_feasible += 1
                if s == 3:
                    s3_owner_feasible += 1
                    s3_mod_histogram[d % 13] += 1
                    assert delta == 1
                    assert d % 13 in (2, 11)
                    assert gap == 1
                if s == 5:
                    s5_owner_feasible += 1
                    if first_s5_witness is None:
                        first_s5_witness = (d, B, top_minus, top_plus, g_numerator, g_denominator)

                hasher.update(
                    repr((d, B, s, top_minus, top_plus, delta, g_numerator, g_denominator)).encode()
                )

    assert s1_owner_feasible == 0
    assert s3_owner_feasible == 0
    assert s5_owner_feasible > 0
    assert set(s3_mod_histogram) <= {2, 11}
    assert packing_failures == 0
    return (
        limit,
        rows,
        owner_feasible,
        packing_failures,
        packing_equalities,
        owner_equalities,
        s1_rows,
        s1_owner_feasible,
        s3_rows,
        s3_owner_feasible,
        tuple(sorted(s3_mod_histogram.items())),
        s5_rows,
        s5_owner_feasible,
        first_s5_witness,
        hasher.hexdigest(),
    )


def main() -> None:
    coefficients, coefficient_endpoint_tests = coefficient_audit()
    geometry = local_geometry_audit()
    shells = shell_bank_audit()
    certificate = sha256(repr((coefficients, geometry, shells)).encode()).hexdigest()

    print("THM-836 TWO-SHEET OWNER-SHELL PACKING OBSTRUCTION")
    print("arithmetic=integer+fractions.Fraction floating_point=none sampled_circle=none optimizer=none")
    print("hypotheses=signed_q13_complement plus guarded_E_plus_R_subset_H_d")
    print("centre_guard=5dp=1_mod13 implies dp=8_mod13; forced deep point maps to 8/13")
    print()
    print("DIRECTIONAL_EXIT_TABLE")
    print(f"signed_residues={SIGNED_RESIDUES}")
    print(f"each_direction_coefficients={EXPECTED_COEFFICIENTS}")
    print("same_owner_complement=c_plus+c_minus=117")
    print(f"coefficient_rows={len(coefficients)} endpoint_tests={coefficient_endpoint_tests}")
    print("endpoint_semantics=E_closed boundary_at_1/11_included immediate_exterior_excluded")
    print()
    print("LOCAL_GEOMETRY_BANK")
    print(f"rows={geometry[0]} direct_endpoint_tests={geometry[1]}")
    print("central_return=[-2/(143B),2/(143B)]_closed_in_closure_R")
    print("target_component_radius=2/(169d) about forced preimage centre")
    print("necessary_each_side=min_z(c_sigma(z)/u_z)+2/B<=22/(13d)")
    print(f"telemetry_edge_flip_histogram={geometry[2]}")
    print("telemetry_observable=right_exit_time switch=left_exit_time tie=(signed_residue,speed)")
    print("telemetry_each_order=transitive score_histogram_0..9 SCCs_10_singletons HP_1")
    print(f"local_geometry_sha256={geometry[3]}")
    print()
    print("THIN_SHELL_BANK")
    print(f"odd_d_limit={shells[0]} shell_rows={shells[1]} owner_feasible_rows={shells[2]}")
    print(f"packing_failures={shells[3]} packing_equalities={shells[4]} owner_equalities={shells[5]}")
    print(f"s1_rows={shells[6]} s1_owner_feasible={shells[7]}")
    print(f"s3_rows={shells[8]} s3_owner_feasible={shells[9]} s3_mod13_histogram={shells[10]}")
    print(f"s5_rows={shells[11]} s5_owner_feasible={shells[12]} first_s5_witness={shells[13]}")
    print("packing_law=delta(d)<=11s(13d+s)/(2(117d+11s))")
    print("s1_consequence=empty")
    print("s3_consequence=empty_after_mod13_alignment")
    print("combined_consequence=13d<=2B-5")
    print("s5_guardrail=local_owner_constraints_have_survivors; no_stronger_shell_claim")
    print(f"shell_bank_sha256={shells[14]}")
    print()
    print("CARRIER_AUDIT")
    print("theorem_carrier=two_direction_obligations_x_owner_labelled_exit_times_x_mod13_gap")
    print("not_a_tournament=binary_owner_order_loses_simultaneous_left_right_minima_and_interval_packing")
    print("tournament_role=lossy_planning_telemetry_only")
    print()
    print(f"certificate_sha256={certificate}")
    print("PASS: centre guard, closed endpoints, exit coefficients, owner inequalities, and modular packing")


if __name__ == "__main__":
    main()
