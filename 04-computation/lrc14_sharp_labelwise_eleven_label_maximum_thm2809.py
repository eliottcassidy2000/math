#!/usr/bin/env python3
"""Exact finite referee for THM-2809.

On each of the first fourteen THM-2749 rails, enumerate the full
2*2*2*13 THM-2640 configuration bank.  A source-label-zero packet at carry
12 can meet the marked source deep band only through two half geometries.
The exact THM-2640 unit table leaves one configuration, whose present anchor
has the opposite polarity from the marked source.  The unique wrapped
label-one configuration has future digit 12, disjoint from the marked
source's digit 6.  Together these gates close all thirteen twelve-faces.
The residual labels 2,...,12 form a positive eleven-face on all fourteen
rails, for every one of their 2^11 adjacent-edge assignments.
"""

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
    "lrc14_slope7_fixed_configuration_carry_nerve_thm2672.py":
        "83ccf3a38660a92cc990bdf304fd4ea4475339731c3e7e92ad35383ef097f361",
    "lrc14_fully_marked_root_zero_clutch_thm2749.py":
        "93b46b2701db8f72d00fa2ae131f9d9afd3200f32998959af3bb2e1fa2f56841",
}

for dependency, expected_hash in DEPENDENCIES.items():
    payload = (COMP / dependency).read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(payload).hexdigest()
    require(actual_hash == expected_hash,
            f"audited dependency changed: {dependency}")


import lrc14_predecessor_carry_private_root_atlas_thm2640 as atlas
import lrc14_slope7_fixed_configuration_carry_nerve_thm2672 as fixed
import lrc14_fully_marked_root_zero_clutch_thm2749 as marked


P = 13
SOURCE_CARRY = 12
RAIL_COUNT = 14
DENOMINATOR = 182
MARKED_BAND = (169, 181)
SOURCE_CONFIG = (0, 1, 1, 6)
LABEL_ONE_CONFIG = (1, 0, 1, 12)
EXPECTED_ELEVEN_FACE_ROWS = (
    ((1,), 399_580_256_360_672_050_023_360),
    ((6,), 74_205_644_260_590_152_069_760),
    ((2, 3), 724_908_063_903_933_297_548_160),
    ((5,), 565_104_521_676_801_927_300_480),
    ((2, 3), 1_130_171_627_188_809_393_027_840),
    ((4, 6), 682_117_240_653_421_081_629_120),
    ((2, 3), 1_267_162_127_454_119_622_485_760),
    ((5, 6), 941_819_893_732_588_135_224_960),
    ((1, 2, 3), 1_449_825_103_908_006_680_574_720),
    ((5,), 596_479_469_905_204_957_431_360),
    ((1, 3), 676_409_208_657_101_856_256_320),
    ((5,), 562_231_844_838_877_400_066_880),
    ((2,), 582_228_901_121_553_500_855_040),
    ((5,), 399_555_625_773_821_502_585_600),
)


def root(edge, kappa, h, carry):
    return (
        2 * carry
        + (2 * h + kappa) // P
        + (edge == 0)
    ) % P


def pulled_root(edge, kappa, h, label):
    carry = (SOURCE_CARRY + 7 * label) % P
    return (root(edge, kappa, h, carry) - label) % P


def private_half(edge, root_index):
    if root_index == 0:
        return (169, 182) if edge == 0 else (0, 13)
    if edge == 0:
        return 14 * root_index - 13, 14 * root_index
    return 14 * root_index, 14 * root_index + 13


def strict_intersection(left, right):
    low = max(left[0], right[0])
    high = min(left[1], right[1])
    return (low, high) if low < high else None


def prefix_intervals(prefix):
    return tuple(
        (left, left + length)
        for left, length in zip(prefix[0], prefix[1])
    )


def interval_intersection(left, right):
    out = []
    i = 0
    j = 0
    while i < len(left) and j < len(right):
        low = max(left[i][0], right[j][0])
        high = min(left[i][1], right[j][1])
        if low < high:
            out.append((low, high))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def weighted_intersection(weighted, support):
    if not weighted or not support:
        return tuple()
    return tuple(atlas.old.intersect_weighted_union(
        weighted, support
    ))


def weighted_support(weighted):
    return tuple((left, right) for left, right, *_rest in weighted)


def compatible_rows(flags, rail_index, label):
    carry = (SOURCE_CARRY + 7 * label) % P
    rows = []
    for sector in range(2):
        for edge in range(2):
            for kappa in range(2):
                for h in range(P):
                    target_root = root(edge, kappa, h, carry)
                    source_root = pulled_root(edge, kappa, h, label)
                    half = private_half(edge, source_root)
                    if (
                        flags[rail_index][sector][edge]
                             [carry][kappa][h]
                        and strict_intersection(half, MARKED_BAND)
                    ):
                        rows.append((
                            sector,
                            edge,
                            kappa,
                            h,
                            carry,
                            target_root,
                            source_root,
                            half,
                        ))
    return tuple(rows)


def main():
    require(
        (atlas.P, atlas.R, atlas.T)
        == (13, 4_826_809, 297_836_897_838_480),
        "canonical THM2672 scales changed",
    )

    shard = atlas.shard((0, RAIL_COUNT))
    content = shard[1]
    metadata = shard[5]
    rows = shard[6]
    require(
        content == 26 and len(metadata) == len(rows) == RAIL_COUNT,
        "fourteen-rail shard changed",
    )
    flags = fixed.build_flags(rows, content)

    # The slope-seven pullback removes the target-label root shift for every
    # one of the 104 configuration signatures.
    for edge in range(2):
        for kappa in range(2):
            for h in range(P):
                source_root = root(
                    edge, kappa, h, SOURCE_CARRY
                )
                for label in range(P):
                    require(
                        pulled_root(edge, kappa, h, label)
                        == source_root,
                        "slope-seven pulled-root law changed",
                    )

    # Independently of the unit flags, only these two source-coordinate half
    # types can meet the strict marked band.
    geometric_types = set()
    for edge in range(2):
        for source_root in range(P):
            half = private_half(edge, source_root)
            if strict_intersection(half, MARKED_BAND):
                geometric_types.add((edge, source_root, half))
    require(
        geometric_types
        == {
            (0, 0, (169, 182)),
            (1, 12, (168, 181)),
        },
        "marked-band-compatible half classification changed",
    )

    source_rows = []
    label_one_rows = []
    source_unit_counts = []
    for rail_index in range(RAIL_COUNT):
        source_unit_count = 0
        for sector in range(2):
            for edge in range(2):
                for kappa in range(2):
                    for h in range(P):
                        source_unit_count += int(
                            flags[rail_index][sector][edge]
                                 [SOURCE_CARRY][kappa][h]
                        )
        source_unit_counts.append(source_unit_count)

        source = compatible_rows(flags, rail_index, 0)
        label_one = compatible_rows(flags, rail_index, 1)
        require(
            source
            == ((
                *SOURCE_CONFIG,
                12,
                12,
                12,
                (168, 181),
            ),),
            f"source-label compatible configuration changed on rail "
            f"{rail_index}",
        )
        require(
            label_one
            == ((
                *LABEL_ONE_CONFIG,
                6,
                1,
                0,
                (169, 182),
            ),),
            f"label-one wrapped boundary changed on rail {rail_index}",
        )
        source_rows.append(source[0])
        label_one_rows.append(label_one[0])

    require(
        len(set(source_rows)) == len(set(label_one_rows)) == 1,
        "compatible rows stopped being rail-independent",
    )

    # The sole source-compatible row has h=6.  Its label-zero present factor
    # is F_(ell,-6)=F_(ell,7), while THM-2749's marked source construction
    # retains the strict complementary factor.
    require(
        SOURCE_CONFIG[3] == 6
        and (-SOURCE_CONFIG[3]) % P == 7,
        "source-compatible anchor changed",
    )

    # The sole label-one-compatible row has h=12,kappa=1.  Its ordinary
    # future-half prefix lies in (25/26,1).  THM-2749's complete marked source
    # prefix has h=6,kappa=1 and lies in (13/26,14/26).  The slope-seven
    # pullback fixes y={R*x}, so the comparison is in one literal coordinate.
    module, _prefixes, _whole, _masses, _rails, _present, _starts = (
        atlas.core.build_carrier_data()
    )
    ordinary_prefixes = atlas.build_pair_prefixes(module)
    terminal_fork = marked.semantic.deepest_fork(module)
    marked_prefixes = marked.build_semantic_prefixes(
        module, terminal_fork
    )
    prefix_rows = []
    for clock in range(7):
        label_one_prefix = prefix_intervals(
            ordinary_prefixes[1][clock][12][1]
        )
        marked_prefix = prefix_intervals(
            marked_prefixes[0][clock][1]
        )
        require(
            all(
                25 * atlas.T // 26 <= left < right <= atlas.T
                for left, right in label_one_prefix
            ),
            "label-one future digit escaped (25/26,1)",
        )
        require(
            all(
                13 * atlas.T // 26 <= left < right
                <= 14 * atlas.T // 26
                for left, right in marked_prefix
            ),
            "marked future digit escaped (13/26,14/26)",
        )
        meet = interval_intersection(
            label_one_prefix, marked_prefix
        )
        require(
            not meet,
            f"label-one and marked delayed prefixes met at clock {clock}",
        )
        prefix_rows.append((
            clock,
            len(label_one_prefix),
            len(marked_prefix),
            len(meet),
        ))
    prefix_rows = tuple(prefix_rows)

    # The only face avoiding both preceding killers is L={2,...,12}.
    # Each label has the same two unit choices: sector0,kappa1,h6 with
    # either adjacent edge.  On H_mark both private halves contain the
    # entire marked band, and the marked semantic prefix is a subset of
    # their common ordinary prefix.  Hence every one of the 2^11 edge
    # assignments has the same physical intersection.
    source_e3 = marked.semantic.exclusive_source(module, 3)
    sections = marked.semantic_sections(
        module, source_e3, 0, 4
    )
    marked_subset_rows = []
    for clock in range(7):
        ordinary_common = prefix_intervals(
            ordinary_prefixes[0][clock][6][1]
        )
        marked_common = prefix_intervals(
            marked_prefixes[0][clock][1]
        )
        require(
            interval_intersection(marked_common, ordinary_common)
            == marked_common,
            f"marked prefix left the common h6 prefix at clock {clock}",
        )
        marked_subset_rows.append((
            clock, len(ordinary_common), len(marked_common)
        ))
    marked_subset_rows = tuple(marked_subset_rows)

    eleven_face_rows = []
    labels = tuple(range(2, P))
    for rail_index in range(RAIL_COUNT):
        option_rows = []
        for label in labels:
            carry = (SOURCE_CARRY + 7 * label) % P
            options = []
            for edge in range(2):
                if flags[rail_index][0][edge][carry][1][6]:
                    target_root = root(edge, 1, 6, carry)
                    source_root = (
                        target_root - label
                    ) % P
                    options.append((
                        edge, target_root, source_root,
                        private_half(edge, source_root),
                    ))
            require(
                options == [
                    (0, label, 0, (169, 182)),
                    (1, (label - 1) % P, 12, (168, 181)),
                ],
                f"eleven-face edge bank changed on rail {rail_index}, "
                f"label {label}",
            )
            option_rows.append(tuple(options))

        source_vector, target_vector, details = (
            marked.fully_marked_vectors(
                module,
                _rails,
                _present,
                marked_prefixes,
                sections,
                rail_index,
            )
        )
        require(
            source_vector == target_vector,
            f"marked endpoint vectors changed on rail {rail_index}",
        )
        values = []
        for clock in range(7):
            common = tuple(
                row for row in details[clock][0] if row[2]
            )
            for label in labels:
                rail_pullback = fixed.shift_union(
                    weighted_support(_rails[rail_index][3]),
                    label * 7 * atlas.T // atlas.R,
                    atlas.T,
                )
                present_pullback = fixed.shift_union(
                    _present[clock, 7],
                    label * 7 * atlas.T // atlas.R,
                    atlas.T,
                )
                common = weighted_intersection(
                    common, rail_pullback
                )
                common = weighted_intersection(
                    common, present_pullback
                )
                if not common:
                    break
            value = (
                atlas.delayed_carry_pair(
                    common, marked_prefixes[0][clock], {}
                )[SOURCE_CARRY][1]
                if common else 0
            )
            values.append(value)
        values = tuple(values)
        positive_support = tuple(
            clock for clock, value in enumerate(values) if value
        )
        row = (positive_support, sum(values))
        require(
            row == EXPECTED_ELEVEN_FACE_ROWS[rail_index],
            f"eleven-face support or mass changed on rail {rail_index}",
        )
        require(
            len(option_rows) == 11
            and all(len(options) == 2 for options in option_rows),
            f"eleven-face edge-choice count changed on rail {rail_index}",
        )
        eleven_face_rows.append(row)
    eleven_face_rows = tuple(eleven_face_rows)

    print("THM2809 SHARP LABELWISE ELEVEN-LABEL MAXIMUM EXACT REFEREE")
    print(f"first14_rail_metadata={tuple(metadata)}")
    print("configuration_universe_per_rail=2*2*2*13=104")
    print(f"source_unit_config_counts={tuple(source_unit_counts)}")
    print(f"marked_deep_band={MARKED_BAND}/{DENOMINATOR}")
    print(f"geometric_half_types={tuple(sorted(geometric_types))}")
    print(f"unique_source_label0_row={source_rows[0]}")
    print(f"unique_wrapped_label1_row={label_one_rows[0]}")
    print("source_label0_anchor=F_(ell,7)")
    print("marked_source_anchor=F_(ell,7)^c")
    print("label1_future_half_digit=(25,26)/26")
    print("marked_future_half_digit=(13,14)/26")
    print(f"delayed_prefix_piece_rows={prefix_rows}")
    print("twelve_face_failure_census=12_label0_anchor+1_label1_digit")
    print("consequence=all13 arbitrary-configuration twelve-faces have "
          "empty full marked-source intersection")
    print(f"marked_prefix_subset_rows={marked_subset_rows}")
    print("positive_eleven_face_labels=(2,3,4,5,6,7,8,9,10,11,12)")
    print("edge_assignments_per_rail=2^11=2048")
    print(f"positive_eleven_face_rows={eleven_face_rows}")
    print("sharp_maximum_label_cardinality=11; simplex_dimension=10")
    print("scope=first14 THM2749 source rails/carry12; no target endpoint, "
          "outside-rail, row, or LRC14 exclusion")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
