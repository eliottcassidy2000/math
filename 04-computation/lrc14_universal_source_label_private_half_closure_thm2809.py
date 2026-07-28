#!/usr/bin/env python3
"""Exact finite referee for THM-2809.

On each of the first fourteen THM-2749 rails, enumerate the full
2*2*2*13 THM-2640 configuration bank.  A source-label-zero packet at carry
12 can meet the marked source deep band only through two half geometries.
The exact THM-2640 unit table leaves one configuration, whose present anchor
has the opposite polarity from the marked source.  The companion also
records the unique wrapped label-one configuration as the sharp boundary.
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


P = 13
SOURCE_CARRY = 12
RAIL_COUNT = 14
DENOMINATOR = 182
MARKED_BAND = (169, 181)
SOURCE_CONFIG = (0, 1, 1, 6)
LABEL_ONE_CONFIG = (1, 0, 1, 12)


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

    print("THM2809 UNIVERSAL SOURCE-LABEL PRIVATE-HALF CLOSURE EXACT REFEREE")
    print(f"first14_rail_metadata={tuple(metadata)}")
    print("configuration_universe_per_rail=2*2*2*13=104")
    print(f"source_unit_config_counts={tuple(source_unit_counts)}")
    print(f"marked_deep_band={MARKED_BAND}/{DENOMINATOR}")
    print(f"geometric_half_types={tuple(sorted(geometric_types))}")
    print(f"unique_source_label0_row={source_rows[0]}")
    print(f"unique_wrapped_label1_row={label_one_rows[0]}")
    print("source_label0_anchor=F_(ell,7)")
    print("marked_source_anchor=F_(ell,7)^c")
    print("consequence=every labelwise configuration mixture containing "
          "source label0 has empty marked-source intersection")
    print("boundary=faces omitting label0 are not source attachments; "
          "label1 wrapped row is not asserted positive after all factors")
    print("scope=first14 THM2749 source rails/carry12; no target endpoint, "
          "outside-rail, row, or LRC14 exclusion")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
