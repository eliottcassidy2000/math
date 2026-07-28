#!/usr/bin/env python3
"""Exact referee for THM-2794.

The universe is one fixed THM-2672 slope-seven configuration,

    (rail, sector, edge, kappa, h) = (8, 0, 1, 1, 6),

at source carry 12 and missing target label 1.  For every clock, the script
builds the twelve-label common rail, common relative-present support, private
root-12 half, and all delayed carry-12 components.  It compares those open
components with the matching THM-2749 fully marked rail-eight source carrier.

The no-go is fixed-configuration and clock-matched.  It says nothing about a
different THM-2672 configuration, the full rail-eight union, or LRC(14).
"""

from fractions import Fraction
import hashlib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


DEPENDENCIES = {
    "lrc14_a12_minuscule_carry_fourier_boundary_thm2785.py":
        "e1b81f994d194bce1f1c846b85e2cec6ecf96f5fc0568c35adcd943284175478",
    "lrc14_fully_marked_root_zero_clutch_thm2749.py":
        "93b46b2701db8f72d00fa2ae131f9d9afd3200f32998959af3bb2e1fa2f56841",
}

for dependency, expected_hash in DEPENDENCIES.items():
    payload = (COMP / dependency).read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(payload).hexdigest()
    require(actual_hash == expected_hash,
            f"audited dependency changed: {dependency}")


import lrc14_a12_minuscule_carry_fourier_boundary_thm2785 as bridge
import lrc14_fully_marked_root_zero_clutch_thm2749 as marked


atlas = bridge.atlas
rebase = bridge.rebase

P = 13
R = P**6
S = P**5
T = atlas.T
RAIL = 8
SECTOR = 0
EDGE = 1
KAPPA = 1
H = 6
CARRY = 12
MISSING_LABEL = 1
ACTIVE_LABELS = tuple(
    label for label in range(P) if label != MISSING_LABEL
)

EXPECTED_FACET_CLOCK_VALUES = (
    0,
    1_509_311_690_628_117_483_066_960,
    3_259_498_009_127_699_308_489_520,
    3_259_498_009_127_699_308_489_520,
    0, 0, 0,
)
EXPECTED_MARKED_CLOCK_VALUES = (
    0,
    339_633_525_654_239_542_165_440,
    750_593_782_703_678_965_571_520,
    719_200_126_392_878_704_654_080,
    0, 0, 0,
)
EXPECTED_HALF_PIECES = (0, 239, 516, 516, 0, 0, 0)
EXPECTED_MARKED_PIECES = (0, 239, 526, 504, 0, 0, 0)
EXPECTED_COMPONENTS = (0, 165_291, 356_973, 356_973, 0, 0, 0)
EXPECTED_COMPONENT_SEAMS = (0, 94, 207, 202, 0, 0, 0)
EXPECTED_BASE_SEAMS = (0, 189, 416, 406, 0, 0, 0)


def build_common_rail(rails):
    pieces = rails[RAIL][3]
    support = [
        (left, right) for left, right, weight in pieces if weight
    ]
    common = rebase.shift_weighted(pieces, 0, T)
    for label in ACTIVE_LABELS[1:]:
        shift = 7 * label * T // R
        common = atlas.old.intersect_weighted_union(
            common, rebase.shift_union(support, shift, T)
        )
    require(common, "twelve-label common rail became empty")
    return common


def build_clock_half(module, present, common_rail, clock):
    source_present = present[clock, (-H) % P]
    common = atlas.old.intersect_weighted_union(
        common_rail, source_present
    )
    for label in ACTIVE_LABELS[1:]:
        shift = 7 * label * T // R
        common = atlas.old.intersect_weighted_union(
            common,
            rebase.shift_union(source_present, shift, T),
        )
    if not common:
        return tuple()
    return tuple(
        rebase.intersect_root_half(
            common, module.C3, EDGE, 12
        )
    )


def delayed_components(half, prefix):
    """Enumerate every strict carry-12 component in raw denominator R*T."""
    starts, lengths, _ = prefix
    period = P * T
    components = []
    for left, right, weight in half:
        if not weight:
            continue
        left_scaled = left * R
        right_scaled = right * R
        for start, length in zip(starts, lengths):
            low = CARRY * T + start
            high = low + length
            turn = (left_scaled - high) // period + 1
            while turn * period + low < right_scaled:
                component_left = max(
                    left_scaled, turn * period + low
                )
                component_right = min(
                    right_scaled, turn * period + high
                )
                if component_left < component_right:
                    components.append(
                        (component_left, component_right)
                    )
                turn += 1
    components.sort()
    require(
        all(
            components[index][1] <= components[index + 1][0]
            for index in range(len(components) - 1)
        ),
        "enumerated delayed intervals are not strict components",
    )
    return tuple(components)


def raw_marked_intervals(marked_pieces):
    return tuple(
        (left * R, right * R)
        for left, right, weight in marked_pieces
        if weight
    )


def overlap_and_one_sided_seams(left_intervals, right_intervals):
    """Return exact open overlap and oriented closure-contact counts."""
    left_index = 0
    right_index = 0
    overlap = 0
    right_to_left = 0
    left_to_right = 0
    while (
        left_index < len(left_intervals)
        and right_index < len(right_intervals)
    ):
        left_low, left_high = left_intervals[left_index]
        right_low, right_high = right_intervals[right_index]
        low = max(left_low, right_low)
        high = min(left_high, right_high)
        if low < high:
            overlap += high - low
        elif left_high == right_low:
            right_to_left += 1
        elif right_high == left_low:
            left_to_right += 1
        if left_high <= right_high:
            left_index += 1
        else:
            right_index += 1
    return overlap, right_to_left, left_to_right


def scaled_base_intervals(half):
    return tuple(
        (left * R, right * R)
        for left, right, weight in half
        if weight
    )


def main():
    require(
        (atlas.P, atlas.R, atlas.S, S)
        == (13, 4_826_809, 371_293, 371_293),
        "canonical scales changed",
    )
    require(
        ACTIVE_LABELS == (0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
        "carry-12 active-label set changed",
    )

    result = atlas.shard((RAIL, RAIL + 1))
    metadata, rows = result[5], result[6]
    require(
        result[1] == 26 and metadata == ((1, 4, 12),),
        "rail-eight primitive content or metadata changed",
    )
    root = (
        2 * CARRY + (2 * H + KAPPA) // P + (EDGE == 0)
    ) % P
    require(root == 12, "source carry stopped having private root 12")
    require(
        atlas.is_unit(
            rows[0][SECTOR][EDGE][CARRY][KAPPA][H],
            root,
            result[1],
        ),
        "selected source-carry row stopped being a primitive unit",
    )
    require(
        (2 * (6 - CARRY)) % P == MISSING_LABEL,
        "source carry stopped omitting label one",
    )

    module, _, _, _, rails, present, _ = (
        atlas.core.build_carrier_data()
    )
    require(module.C3 == 2 * S,
            "physical deep speed changed")
    prefixes = atlas.build_pair_prefixes(module)
    common_rail = build_common_rail(rails)

    source_e3 = marked.semantic.exclusive_source(module, 3)
    terminal_fork = marked.semantic.deepest_fork(module)
    semantic_prefixes = marked.build_semantic_prefixes(
        module, terminal_fork
    )
    sections = marked.semantic_sections(module, source_e3, 0, 4)
    source_vector, target_vector, marked_details = (
        marked.fully_marked_vectors(
            module,
            rails,
            present,
            semantic_prefixes,
            sections,
            RAIL,
        )
    )
    require(
        source_vector == target_vector
        and source_vector == EXPECTED_MARKED_CLOCK_VALUES,
        "fully marked rail-eight positive-control vector changed",
    )

    facet_clock_values = []
    half_piece_counts = []
    marked_piece_counts = []
    component_counts = []
    component_seams = []
    base_seams = []
    component_overlap_raw = []
    base_overlap_raw = []

    for clock in range(7):
        half = build_clock_half(
            module, present, common_rail, clock
        )
        half_piece_counts.append(len(half))
        marked_source = tuple(
            piece for piece in marked_details[clock][0] if piece[2]
        )
        marked_piece_counts.append(len(marked_source))

        if half:
            values = atlas.delayed_carry_pair(
                half, prefixes[SECTOR][clock][H], {}
            )
            facet_value = values[CARRY][KAPPA]
            components = delayed_components(
                half, prefixes[SECTOR][clock][H][KAPPA]
            )
        else:
            facet_value = 0
            components = tuple()
        facet_clock_values.append(facet_value)
        component_counts.append(len(components))

        marked_raw = raw_marked_intervals(marked_source)
        overlap, right_left, left_right = (
            overlap_and_one_sided_seams(components, marked_raw)
        )
        component_overlap_raw.append(overlap)
        component_seams.append(right_left + left_right)
        require(
            left_right == 0,
            "a marked carrier ended at a delayed-component start",
        )

        base_raw = scaled_base_intervals(half)
        base_overlap, base_right_left, base_left_right = (
            overlap_and_one_sided_seams(base_raw, marked_raw)
        )
        base_overlap_raw.append(base_overlap)
        base_seams.append(base_right_left + base_left_right)
        require(
            base_left_right == 0,
            "a marked carrier ended at a fixed-facet base start",
        )

        # Independent aggregate route: intersect the fixed-facet base with
        # the fully marked carrier before applying the delayed prefix.
        marked_support = [
            (left, right) for left, right, _weight in marked_source
        ]
        restricted = atlas.old.intersect_weighted_union(
            half, marked_support
        )
        require(
            not restricted,
            f"clock {clock} acquired a pre-prefix semantic atom",
        )

    facet_clock_values = tuple(facet_clock_values)
    half_piece_counts = tuple(half_piece_counts)
    marked_piece_counts = tuple(marked_piece_counts)
    component_counts = tuple(component_counts)
    component_seams = tuple(component_seams)
    base_seams = tuple(base_seams)
    component_overlap_raw = tuple(component_overlap_raw)
    base_overlap_raw = tuple(base_overlap_raw)

    require(
        facet_clock_values == EXPECTED_FACET_CLOCK_VALUES
        and half_piece_counts == EXPECTED_HALF_PIECES
        and marked_piece_counts == EXPECTED_MARKED_PIECES
        and component_counts == EXPECTED_COMPONENTS
        and component_seams == EXPECTED_COMPONENT_SEAMS
        and base_seams == EXPECTED_BASE_SEAMS,
        "fixed-configuration component census changed",
    )
    require(
        component_overlap_raw == (0,) * 7
        and base_overlap_raw == (0,) * 7,
        "a positive semantic overlap appeared",
    )
    require(
        sum(facet_clock_values)
        == 8_028_307_708_883_516_100_046_000
        and sum(component_counts) == 879_237
        and sum(component_seams) == 503
        and sum(base_seams) == 1_011,
        "aggregate fixed-configuration census changed",
    )
    require(
        tuple(
            clock for clock, value in enumerate(facet_clock_values)
            if value
        ) == (1, 2, 3)
        and tuple(
            clock for clock, value in enumerate(source_vector) if value
        ) == (1, 2, 3),
        "the two positive clock supports stopped matching",
    )

    print("THM2794 RAIL-EIGHT ALL-COMPONENT SEMANTIC-ATOM NO-GO EXACT REFEREE")
    print("config=(rail8,sector0,edge1,kappa1,h6) "
          "metadata=(1,4,12)")
    print("source_carry=12 private_root=12 missing_label=1 "
          f"active_labels={ACTIVE_LABELS}")
    print(f"facet_clock_values={facet_clock_values}")
    print(f"marked_clock_values={source_vector}")
    print(f"fixed_base_piece_counts={half_piece_counts}")
    print(f"marked_source_piece_counts={marked_piece_counts}")
    print(f"delayed_component_counts={component_counts}")
    print(f"component_positive_overlap_raw={component_overlap_raw}")
    print(f"component_closure_seams={component_seams}")
    print(f"base_positive_overlap_raw={base_overlap_raw}")
    print(f"base_closure_seams={base_seams}")
    print("positive_clock_supports=facet:(1,2,3) marked:(1,2,3)")
    print("total_positive_delayed_components=879237 "
          "total_component_seams=503")
    print("stronger_route=fixed-facet base intersect marked source is empty "
          "before delayed prefix on every clock")
    print("mechanism=all closure contacts are facet-right to marked-left; "
          "strict open supports do not cross")
    print("scope=one carry and one fixed rail8 configuration; "
          "no other-configuration, full-union, row, or LRC14 exclusion")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
