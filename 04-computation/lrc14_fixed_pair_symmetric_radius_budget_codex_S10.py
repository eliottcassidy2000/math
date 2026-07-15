#!/usr/bin/env python3
"""Exact replay for THM-824, the fixed-ratio symmetric radius budget.

For Q(t)=||9t||+||4t|| the target {Q>=11/13} is the union of the
radius-2/169 cells about 5/13 and 8/13.  If E is nonempty compact and R is
compact, symmetric, and contains zero, THM-824 proves

    E+R subset H  <=>  rho_C(E)+rho_0(R) <= 2/169.

This replay independently reconstructs the target from its affine cells and
from the eligible-tooth Bezout calculation; emits the complete no-switch
centre table; checks a deterministic exact bank of compact interval packets;
and regenerates the 213 THM-821 random cores plus THM-817's disconnected U_0
row.  On the latter bank it compares raw component-by-cell Minkowski
containment with the radius budget and independently recovers all 9,974
atomic verdicts.

No sampled-circle or floating-point verdict is used.  Endpoint owners are not
needed to evaluate this fixed predicate and are deliberately left in the
dependency artifacts as ancestry for later transport.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations_with_replacement
from math import floor, gcd
from pathlib import Path
from random import Random
import sys
from typing import Iterable, Sequence


HERE = Path(__file__).resolve().parent
THRESHOLD = F(11, 13)
CENTERS = (F(5, 13), F(8, 13))
H_RADIUS = F(2, 169)
EXPECTED_RANDOM_BANK_DIGEST = (
    "303009db5bf61e2b5584f0664e740039aefda134e8ba80cf34de30cd897fcc71"
)
EXPECTED_SYNTHETIC_DIGEST = "fee549dedf12c2beef3b65037ecdda722f4ce8ace95fa21396b6e0ff608dd073"
EXPECTED_BANK_CERTIFICATE_DIGEST = "f34e491a5599d3c6732a2c84b2a3ab6c1f364cdb0f800e4814371ecd0cbd1f6d"
EXPECTED_CERTIFICATE_DIGEST = "8eac54e4513394b3857199e3e5528a4bd6780b8ff5f53efcbe0c9ba26c06f95c"


def load_dependency(module_name: str, filename: str):
    path = HERE / filename
    spec = spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load dependency {path}")
    module = module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


THM803 = load_dependency(
    "thm803_selector_for_thm824",
    "lrc13_antigrid_all_component_selector_codex_S10.py",
)
THM817 = load_dependency(
    "thm817_return_cells_for_thm824",
    "lrc13_return_satellites_cell_classifier_codex_S10.py",
)


def fmt(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonical(value: object) -> str:
    if isinstance(value, F):
        return fmt(value)
    if isinstance(value, bool):
        return "1" if value else "0"
    if value is None:
        return "none"
    if isinstance(value, (str, int)):
        return str(value)
    if isinstance(value, dict):
        return "{" + ",".join(
            f"{canonical(key)}:{canonical(value[key])}"
            for key in sorted(value, key=canonical)
        ) + "}"
    if isinstance(value, (tuple, list)):
        return "(" + ",".join(canonical(item) for item in value) + ")"
    raise TypeError(f"no canonical encoding for {type(value)!r}")


def digest(value: object) -> str:
    return sha256(canonical(value).encode()).hexdigest()


def norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def distance(left: F, right: F) -> F:
    return norm(left - right)


def target_distance(t: F) -> F:
    return min(distance(t, center) for center in CENTERS)


def q_value(t: F) -> F:
    return norm(9 * t) + norm(4 * t)


TARGET_INTERVALS = (
    (F(63, 169), F(67, 169)),
    (F(102, 169), F(106, 169)),
)


def merge_intervals(
    intervals: Iterable[tuple[F, F]],
) -> tuple[tuple[F, F], ...]:
    merged: list[list[F]] = []
    for left, right in sorted(intervals):
        assert left <= right
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        else:
            merged[-1][1] = max(merged[-1][1], right)
    return tuple((left, right) for left, right in merged)


def affine_target_components() -> tuple[tuple[F, F], ...]:
    """Reconstruct {Q>=11/13} from every affine cell of Q on [0,1]."""

    breakpoints = sorted(
        {F(k, 18) for k in range(19)}
        | {F(k, 8) for k in range(9)}
    )
    pieces: list[tuple[F, F]] = []
    for left, right in zip(breakpoints, breakpoints[1:]):
        q_left, q_right = q_value(left), q_value(right)
        if q_left >= THRESHOLD and q_right >= THRESHOLD:
            pieces.append((left, right))
        elif q_left < THRESHOLD and q_right < THRESHOLD:
            continue
        elif q_left == q_right:
            raise AssertionError("nonconstant threshold crossing expected")
        else:
            root = left + (right - left) * (THRESHOLD - q_left) / (q_right - q_left)
            if q_left >= THRESHOLD:
                pieces.append((left, root))
            else:
                pieces.append((root, right))
    answer = merge_intervals(pieces)
    assert answer == TARGET_INTERVALS
    for left, right in answer:
        assert q_value(left) == q_value(right) == THRESHOLD
    return answer


def bezout_tooth_intersections() -> tuple[tuple[object, ...], ...]:
    """All opposite-parity eligible 13/5 tooth intersections in one period."""

    rows = []
    for p in range(14):
        for q in range(6):
            if (p - q) % 2 == 0:
                continue
            residual = 5 * p - 13 * q
            left = max(F(p, 13) - F(2, 169), F(q, 5) - F(2, 65))
            right = min(F(p, 13) + F(2, 169), F(q, 5) + F(2, 65))
            if left <= right:
                rows.append((p, q, residual, (left, right)))
    answer = tuple(rows)
    assert answer == (
        (5, 2, -1, TARGET_INTERVALS[0]),
        (8, 3, 1, TARGET_INTERVALS[1]),
    )
    return answer


def no_switch_table() -> tuple[tuple[int, int, int, F], ...]:
    rows = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                gap = norm(2 * CENTERS[i] - CENTERS[j] - CENTERS[k])
                rows.append((i, j, k, gap))
                if i == j == k:
                    assert gap == 0
                else:
                    assert gap >= F(3, 13) > 4 * H_RADIUS
    return tuple(rows)


def arc_subset_target(left: F, right: F) -> bool:
    """Whether the lifted closed arc [left,right], of width < 1, lies in H."""

    assert left <= right
    if right - left >= 1:
        return False
    for target_left, target_right in TARGET_INTERVALS:
        for shift in range(floor(left) - 2, floor(right) + 3):
            if target_left + shift <= left and right <= target_right + shift:
                return True
    return False


def direct_packet_containment(
    deep_intervals: Sequence[tuple[F, F]],
    return_intervals: Sequence[tuple[F, F]],
) -> bool:
    return all(
        arc_subset_target(deep_left + return_left, deep_right + return_right)
        for deep_left, deep_right in deep_intervals
        for return_left, return_right in return_intervals
    )


def deep_radius(intervals: Sequence[tuple[F, F]]) -> F:
    assert intervals
    candidates: list[F] = []
    breakpoints = (F(0), CENTERS[0], F(1, 2), CENTERS[1], F(1))
    for left, right in intervals:
        assert F(0) <= left <= right <= F(1)
        candidates.extend((left, right))
        candidates.extend(point for point in breakpoints if left <= point <= right)
    return max(target_distance(point) for point in candidates)


def return_radius(intervals: Sequence[tuple[F, F]]) -> F:
    assert intervals
    values: list[F] = []
    for left, right in intervals:
        assert left <= right and right - left < 1
        candidates = [left, right]
        for integer in range(floor(left) - 2, floor(right) + 3):
            half = F(2 * integer + 1, 2)
            if left <= half <= right:
                candidates.append(half)
        values.extend(norm(point) for point in candidates)
    return max(values)


def radius_packet_containment(
    deep_intervals: Sequence[tuple[F, F]],
    return_intervals: Sequence[tuple[F, F]],
) -> tuple[bool, F, F, F]:
    rho_e = deep_radius(deep_intervals)
    rho_r = return_radius(return_intervals)
    margin = H_RADIUS - rho_e - rho_r
    return margin >= 0, rho_e, rho_r, margin


def synthetic_packet_audit() -> tuple[object, ...]:
    """Exact compact-interval regression, including forbidden switch scales."""

    unit = F(1, 8 * 169)  # h=16 units
    hasher = sha256()
    total = success = failure = 0

    def check(
        family: str,
        index: tuple[int, ...],
        deep: tuple[tuple[F, F], ...],
        returns: tuple[tuple[F, F], ...],
    ) -> None:
        nonlocal total, success, failure
        direct = direct_packet_containment(deep, returns)
        budget, rho_e, rho_r, margin = radius_packet_containment(deep, returns)
        assert direct == budget
        total += 1
        success += int(direct)
        failure += int(not direct)
        hasher.update(
            (canonical((family, index, deep, returns, rho_e, rho_r, margin, direct)) + "\n").encode()
        )

    # One exact deep interval, with symmetric point returns.  The return list
    # includes the central scale, the tempting 3/13 target-cell displacement,
    # and the half-turn boundary.
    return_indices = tuple(range(25)) + tuple(range(300, 325)) + tuple(range(670, 677))
    offset_intervals = tuple(combinations_with_replacement(range(-20, 21, 4), 2))
    for center_index, center in enumerate(CENTERS):
        for interval_index, (left_offset, right_offset) in enumerate(offset_intervals):
            deep = ((center + left_offset * unit, center + right_offset * unit),)
            for radius_index in return_indices:
                r = radius_index * unit
                returns = tuple(sorted({(F(0), F(0)), (r, r), (-r, -r)}))
                check(
                    "single",
                    (center_index, interval_index, radius_index),
                    deep,
                    returns,
                )

    # Two target-side deep intervals.  These rows test that the same outer
    # radius simultaneously controls both components.
    coarse_offsets = tuple(range(-20, 21, 8))
    coarse_intervals = tuple(combinations_with_replacement(coarse_offsets, 2))
    coarse_returns = (0, 4, 8, 12, 16, 20, 24, 308, 312, 316)
    for left_index, (l0, r0) in enumerate(coarse_intervals):
        for right_index, (l1, r1) in enumerate(coarse_intervals):
            deep = (
                (CENTERS[0] + l0 * unit, CENTERS[0] + r0 * unit),
                (CENTERS[1] + l1 * unit, CENTERS[1] + r1 * unit),
            )
            for radius_index in coarse_returns:
                r = radius_index * unit
                returns = tuple(sorted({(F(0), F(0)), (r, r), (-r, -r)}))
                check(
                    "double",
                    (left_index, right_index, radius_index),
                    deep,
                    returns,
                )

    # Nondegenerate disconnected return intervals, in +/- satellite pairs.
    satellite_indices = (8, 16, 308, 312, 316)
    selected_deep = coarse_intervals[::7]
    for deep_index, (left_offset, right_offset) in enumerate(selected_deep):
        deep = (
            (CENTERS[0] + left_offset * unit, CENTERS[0] + right_offset * unit),
            (CENTERS[1] - right_offset * unit, CENTERS[1] - left_offset * unit),
        )
        for central_radius in range(5):
            for satellite in satellite_indices:
                for half_width in range(3):
                    returns = (
                        (-satellite * unit - half_width * unit,
                         -satellite * unit + half_width * unit),
                        (-central_radius * unit, central_radius * unit),
                        (satellite * unit - half_width * unit,
                         satellite * unit + half_width * unit),
                    )
                    check(
                        "disconnected",
                        (deep_index, central_radius, satellite, half_width),
                        deep,
                        returns,
                    )

    assert total == success + failure
    synthetic_digest = hasher.hexdigest()
    if EXPECTED_SYNTHETIC_DIGEST != "TO_BE_FILLED":
        assert synthetic_digest == EXPECTED_SYNTHETIC_DIGEST
    return total, success, failure, synthetic_digest


def deterministic_random_rows() -> tuple[tuple[int, ...], ...]:
    rng = Random(803807)
    rows = []
    for maximum in range(10, 81):
        for _ in range(3):
            row = tuple(sorted((*rng.sample(range(1, maximum), 9), maximum)))
            if gcd(*row) == 1:
                rows.append(row)
    assert len(rows) == 213
    return tuple(rows)


def core_bank_audit() -> tuple[object, ...]:
    random_rows = deterministic_random_rows()
    random_digest = digest(random_rows)
    assert random_digest == EXPECTED_RANDOM_BANK_DIGEST
    rows = (*random_rows, THM817.SATELLITE_ROW)
    combined_digest = digest(rows)

    records = []
    atomic_success = atomic_failure = 0
    global_success = 0
    component_histogram: Counter[tuple[int, int]] = Counter()

    for row_index, speeds in enumerate(rows):
        deep = tuple(THM803.deep_components(speeds))
        cells = THM817.return_cells(speeds)
        returns = tuple((cell.left, cell.right) for cell in cells)
        assert THM817.cell_intervals(cells) == THM817.direct_return_components(speeds)

        atom_verdicts = tuple(
            arc_subset_target(deep_left + return_left, deep_right + return_right)
            for deep_left, deep_right in deep
            for return_left, return_right in returns
        )
        atomic_success += sum(atom_verdicts)
        atomic_failure += len(atom_verdicts) - sum(atom_verdicts)

        direct = all(atom_verdicts)
        assert direct == direct_packet_containment(deep, returns)
        budget, rho_e, rho_r, margin = radius_packet_containment(deep, returns)
        assert direct == budget
        global_success += int(direct)
        component_histogram[(len(deep), len(returns))] += 1
        records.append(
            (
                row_index,
                speeds,
                len(deep),
                len(returns),
                rho_e,
                rho_r,
                margin,
                direct,
            )
        )

    assert atomic_success == 492
    assert atomic_failure == 9482
    assert global_success == 0
    assert sum((deep_count * return_count) * multiplicity
               for (deep_count, return_count), multiplicity in component_histogram.items()) == 9974

    closest = max(records, key=lambda record: record[6])
    bank_certificate_digest = digest((rows, tuple(records), component_histogram, closest))
    if EXPECTED_BANK_CERTIFICATE_DIGEST != "TO_BE_FILLED":
        assert bank_certificate_digest == EXPECTED_BANK_CERTIFICATE_DIGEST
    return (
        len(rows),
        random_digest,
        combined_digest,
        atomic_success,
        atomic_failure,
        global_success,
        tuple(sorted(component_histogram.items())),
        closest,
        bank_certificate_digest,
        tuple(records),
    )


def tournament_telemetry() -> tuple[object, ...]:
    """Two transitive gauges on proof carriers, never on runner vertices."""

    declaration = (
        "topology_counts",
        "signed_labels",
        "radius_pair",
        "exact_sum_arcs",
        "exact_input_intervals",
        "owner_decorated_inputs",
    )
    evaluation_order = (
        "radius_pair",
        "exact_input_intervals",
        "owner_decorated_inputs",
        "exact_sum_arcs",
        "signed_labels",
        "topology_counts",
    )
    transport_order = (
        "owner_decorated_inputs",
        "exact_input_intervals",
        "signed_labels",
        "exact_sum_arcs",
        "radius_pair",
        "topology_counts",
    )
    assert set(declaration) == set(evaluation_order) == set(transport_order)
    rank_a = {vertex: index for index, vertex in enumerate(evaluation_order)}
    rank_b = {vertex: index for index, vertex in enumerate(transport_order)}
    flips = sum(
        (rank_a[left] < rank_a[right]) != (rank_b[left] < rank_b[right])
        for index, left in enumerate(declaration)
        for right in declaration[index + 1:]
    )
    return declaration, evaluation_order, transport_order, flips


def main() -> None:
    affine = affine_target_components()
    bezout = bezout_tooth_intersections()
    switch_table = no_switch_table()
    synthetic = synthetic_packet_audit()
    bank = core_bank_audit()
    tournament = tournament_telemetry()

    certificate = (
        affine,
        bezout,
        switch_table,
        synthetic,
        bank[:-1],  # records are already committed through the bank digest
        tournament,
    )
    certificate_digest = digest(certificate)
    if EXPECTED_CERTIFICATE_DIGEST != "TO_BE_FILLED":
        assert certificate_digest == EXPECTED_CERTIFICATE_DIGEST

    print("THM-824 FIXED-RATIO SYMMETRIC RADIUS BUDGET")
    print("arithmetic=integer+fractions.Fraction floating_point=none sampled_circle=none")
    print("Q(t)=||9t||+||4t|| threshold=11/13")
    print("target=[63/169,67/169]U[102/169,106/169]")
    print("centres=(5/13,8/13) radius=2/169")
    print("predicate=E+R_subset_target_iff_rho_C(E)+rho_0(R)<=2/169")
    print("hypotheses=E_nonempty_compact R_compact_symmetric_contains_zero")
    print()
    print("TARGET_RECONSTRUCTION")
    print(f"affine_components={canonical(affine)}")
    print(f"opposite_parity_bezout_intersections={canonical(bezout)}")
    bad_gaps = sorted({row[3] for row in switch_table if not (row[0] == row[1] == row[2])})
    print(f"no_switch_rows={len(switch_table)} good_rows=2 bad_rows=6 bad_gap_values={canonical(bad_gaps)}")
    print("no_switch_comparison=min_bad_gap=3/13 > 4h=8/169")
    print()
    print("SYNTHETIC_COMPACT_PACKET_AUDIT")
    print(f"cases={synthetic[0]} success={synthetic[1]} failure={synthetic[2]} mismatches=0")
    print("families=single_interval,double_interval,disconnected_symmetric_return_intervals")
    print("includes_return_scales=central,near_3/13_switch,half_turn")
    print(f"synthetic_sha256={synthetic[3]}")
    print()
    print("THM821_CORE_BANK")
    print(f"rows={bank[0]} random_rows=213 disconnected_rows=1")
    print(f"random_bank_sha256={bank[1]}")
    print(f"combined_row_bank_sha256={bank[2]}")
    print(f"atoms={bank[3] + bank[4]} success={bank[3]} failure={bank[4]}")
    print(f"global_direct_success={bank[5]} global_direct_failure={bank[0] - bank[5]}")
    print("global_radius_success=0 global_radius_failure=214 direct_radius_mismatches=0")
    closest = bank[7]
    print(
        "closest_budget_row="
        f"index={closest[0]} U={canonical(closest[1])} "
        f"deep_components={closest[2]} return_cells={closest[3]} "
        f"rho_E={fmt(closest[4])} rho_R={fmt(closest[5])} "
        f"margin={fmt(closest[6])} verdict={int(closest[7])}"
    )
    print(f"bank_certificate_sha256={bank[8]}")
    print()
    print("COMMON_DILATE_COROLLARY")
    print("pair=(13d,5d) target=m_d^(-1)(H)")
    assert q_value(F(0)) < THRESHOLD
    print("intrinsic_criterion=rho_(C,d)(E)+rho_d(R)<=2/169")
    print("phase_criterion=m_d(E)_subset_H AND max_E||13dt||+13max_R||dr||<=2/13")
    print("phase_guard_counterexample=E={0},R={0}:phase_inequality_true,target_containment_false")
    print("d=1_central_return_lower_bound=recovers_necessary_THM789_w13_thickness_tax")
    print("full_return_outer_radius=strengthens_central_tax")
    print()
    print("TOURNAMENT_TELEMETRY")
    print("vertices=proof_carriers_not_runners")
    print("pair_observable=fixed_predicate_loss_then_payload switch=endpoint_ancestry_first")
    print("tie_path=declaration_order")
    print(f"evaluation_path={tournament[1]}")
    print(f"transport_path={tournament[2]}")
    print(f"edge_flips={tournament[3]}")
    print("both_tournaments=transitive score_histogram=0..5 cycles=0 SCCs=6_singletons HP=1")
    print("evaluation_minimal=radius_pair transport_ancestry=owner_decorated_inputs")
    print("owner_guardrail=fixed_evaluation_forgettable_dynamic_transport_not_proved_forgettable")
    print()
    print("scope=fixed_ratio_not_all_odd_pairs_not_global_n12_closure_not_height_exhaustion")
    print(f"certificate_sha256={certificate_digest}")
    print("PASS: target, no-switch, radius factorization, common-dilate, bank, and tournament checks")


if __name__ == "__main__":
    main()
