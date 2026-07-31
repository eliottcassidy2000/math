#!/usr/bin/env python3
"""Independent hostile audit of THM-2819.

This script reconstructs the source/target endpoint conjugacy from the
proved THM-2640, THM-2672, and THM-2749 carriers.  It deliberately does not
import the THM-2819 primary companion.  Only after all independent
mathematical checks have passed does it hash that primary artifact.

The audit pins:

* the sign and integer lift in the nonwrapping translation;
* carry and private-half relabelling;
* the target-label-zero delayed-prefix obstruction;
* the target-label-twelve odometer residual and its ordered factor death;
* strict interval endpoints;
* direct target-coordinate positive and zero rows on all fourteen rails;
* every one of the 2^11 adjacent-edge choices; and
* the exact thirteen-face census and theorem scope.
"""

from itertools import product
import hashlib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    return path.read_bytes().replace(bytes((13, 10)), bytes((10,)))


DEPENDENCIES = {
    "lrc14_predecessor_carry_private_root_atlas_thm2640.py":
        "a28b03a5903256c1c1c294ea5af389c7991fc0a5ad6908f0f25a5b0cc6e71abf",
    "lrc14_slope7_fixed_configuration_carry_nerve_thm2672.py":
        "83ccf3a38660a92cc990bdf304fd4ea4475339731c3e7e92ad35383ef097f361",
    "lrc14_fully_marked_root_zero_clutch_thm2749.py":
        "93b46b2701db8f72d00fa2ae131f9d9afd3200f32998959af3bb2e1fa2f56841",
}

for dependency, expected_hash in DEPENDENCIES.items():
    actual_hash = hashlib.sha256(lf_bytes(COMP / dependency)).hexdigest()
    require(
        actual_hash == expected_hash,
        f"audited dependency changed: {dependency}",
    )


import lrc14_predecessor_carry_private_root_atlas_thm2640 as atlas
import lrc14_slope7_fixed_configuration_carry_nerve_thm2672 as fixed
import lrc14_root_zero_overlap_clutch_20260728 as clutch
import lrc14_fully_marked_root_zero_clutch_thm2749 as marked


P = 13
R = P**6
S = P**5
T = atlas.T
SHIFT = 7 * T // R
C3 = 2 * S
SOURCE_CARRY = 12
TARGET_CARRY = 6
SOURCE_BAND = (169, 181)
TARGET_BAND = (1, 13)

EXPECTED_POSITIVE_ROWS = (
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

EXPECTED_WRAP_ORDER = (
    (784_513_810_080, 728_825_300_280, 0),
    (356_159_643_840, 330_902_579_700, 0),
    (473_353_040_160, 439_916_713_920, 0),
    (625_685_860_800, 581_278_296_060, 0),
    (620_815_595_400, 578_179_396_872, 0),
    (438_008_061_491, 407_334_158_691, 0),
    (394_663_389_120, 366_653_853_720, 0),
)


def root(edge, kappa, height, carry):
    return (
        2 * carry
        + (2 * height + kappa) // P
        + (edge == 0)
    ) % P


def private_half(edge, root_index):
    if root_index == 0:
        return (169, 182) if edge == 0 else (0, 13)
    if edge == 0:
        return 14 * root_index - 13, 14 * root_index
    return 14 * root_index, 14 * root_index + 13


def strict_meet(left, right):
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


def weighted_support(weighted):
    return tuple((left, right) for left, right, weight in weighted if weight)


def weighted_intersection(weighted, support):
    if not weighted or not support:
        return tuple()
    return tuple(atlas.old.intersect_weighted_union(weighted, support))


def weighted_mass(weighted):
    return sum(right - left for left, right, *_rest in weighted)


def normalize_weighted(pieces):
    """Canonicalize adjacent pieces solely for exact transport comparison."""
    out = []
    for left, right, weight in sorted(pieces):
        if left >= right:
            continue
        if out and out[-1][1] == left and out[-1][2] == weight:
            out[-1] = (out[-1][0], right, weight)
        else:
            out.append((left, right, weight))
    return tuple(out)


def target_rows(flags, rail_index, label):
    """Unit rows whose pulled half meets the strict marked target band."""
    carry = (TARGET_CARRY + 7 * label) % P
    rows = []
    for sector in range(2):
        for edge in range(2):
            for kappa in range(2):
                for height in range(P):
                    target_root = root(edge, kappa, height, carry)
                    # Pullback by label*tau subtracts one deep-root step per
                    # label because 182*C3*tau=196=14 mod 182.
                    pulled_root = (target_root - label) % P
                    half = private_half(edge, pulled_root)
                    if (
                        flags[rail_index][sector][edge]
                             [carry][kappa][height]
                        and strict_meet(half, TARGET_BAND)
                    ):
                        rows.append((
                            sector,
                            edge,
                            kappa,
                            height,
                            carry,
                            target_root,
                            pulled_root,
                            half,
                        ))
    return tuple(rows)


def main():
    require(
        (atlas.P, R, S, T, SHIFT, C3)
        == (
            13,
            4_826_809,
            371_293,
            297_836_897_838_480,
            431_933_040,
            742_586,
        ),
        "canonical endpoint scales changed",
    )
    require(
        R * SHIFT == 7 * T
        and 182 * C3 * SHIFT == 196 * T,
        "tau no longer has future/deep increments (0,14)",
    )

    module, _unused_prefixes, _whole, _masses, rails, present, _starts = (
        atlas.core.build_carrier_data()
    )
    require(
        module.C3 == C3 and len(rails) >= 14,
        "canonical marked carrier changed",
    )

    shard = atlas.shard((0, 14))
    content = shard[1]
    metadata = shard[5]
    raw_rows = shard[6]
    require(
        content == 26
        and len(metadata) == len(raw_rows) == 14,
        "first-fourteen THM2640 shard changed",
    )
    flags = fixed.build_flags(raw_rows, content)

    # 1. Reconstruct the nonwrapping conjugacy, including its orientation.
    # In grid coordinates preimage under x -> x+d*tau is shift by -d*SHIFT.
    nonwrap_rail_checks = 0
    nonwrap_present_checks = 0
    carry_rows = []
    for source_label in range(1, P):
        target_label = source_label - 1
        source_carry_integer = 12 + 7 * source_label
        target_carry_integer = 6 + 7 * target_label
        require(
            source_carry_integer - target_carry_integer == 13
            and source_carry_integer % P
            == target_carry_integer % P,
            "source/target carry relabelling changed",
        )
        carry_rows.append((
            source_label,
            target_label,
            source_carry_integer % P,
        ))

        for rail_index in range(14):
            base = tuple(row for row in rails[rail_index][3] if row[2])
            source_preimage = clutch.shift_weighted(
                base, -source_label * SHIFT
            )
            translated_source = clutch.shift_weighted(
                source_preimage, SHIFT
            )
            target_preimage = clutch.shift_weighted(
                base, -target_label * SHIFT
            )
            require(
                normalize_weighted(translated_source)
                == normalize_weighted(target_preimage),
                f"nonwrapping rail sign failed at rail {rail_index}, "
                f"source label {source_label}",
            )
            nonwrap_rail_checks += 1

        # Check every present bank, not merely the h=6 bank used by the
        # positive face.
        for clock in range(7):
            for present_label in range(P):
                base = tuple(present[clock, present_label])
                source_preimage = clutch.shift_union(
                    base, -source_label * SHIFT
                )
                translated_source = clutch.shift_union(
                    source_preimage, SHIFT
                )
                target_preimage = clutch.shift_union(
                    base, -target_label * SHIFT
                )
                require(
                    translated_source == target_preimage,
                    f"nonwrapping present sign failed at clock {clock}, "
                    f"present label {present_label}, source label "
                    f"{source_label}",
                )
                nonwrap_present_checks += 1

    require(
        nonwrap_rail_checks == 14 * 12
        and nonwrap_present_checks == 7 * 13 * 12,
        "nonwrapping packet-check universe changed",
    )

    # 2. Rebuild the strict target half geometry and its three distinguished
    # label families on every rail.
    geometric_types = set()
    for edge in range(2):
        for pulled_root in range(P):
            half = private_half(edge, pulled_root)
            if strict_meet(half, TARGET_BAND):
                geometric_types.add((edge, pulled_root, half))
    require(
        geometric_types
        == {
            (0, 1, (1, 14)),
            (1, 0, (0, 13)),
        },
        "strict target-band half classification changed",
    )
    require(
        strict_meet((0, 1), TARGET_BAND) is None
        and strict_meet((13, 14), TARGET_BAND) is None,
        "a boundary-only half contact was treated as positive",
    )

    target_zero_rows = []
    target_twelve_rows = []
    positive_option_rows = []
    positive_full_bank_counts = []
    for rail_index in range(14):
        zero = target_rows(flags, rail_index, 0)
        twelve = target_rows(flags, rail_index, 12)
        require(
            zero
            == ((
                1, 0, 1, 12, 6, 1, 1, (1, 14),
            ),),
            f"target-zero compatible row changed on rail {rail_index}",
        )
        require(
            twelve
            == ((
                0, 1, 1, 6, 12, 12, 0, (0, 13),
            ),),
            f"target-twelve compatible row changed on rail {rail_index}",
        )
        target_zero_rows.append(zero[0])
        target_twelve_rows.append(twelve[0])

        per_label = []
        per_label_counts = []
        for label in range(1, 12):
            carry = (TARGET_CARRY + 7 * label) % P
            expected = (
                (
                    0, 0, 1, 6, carry,
                    (label + 1) % P, 1, (1, 14),
                ),
                (
                    0, 1, 1, 6, carry,
                    label, 0, (0, 13),
                ),
            )
            actual = target_rows(flags, rail_index, label)
            selected = tuple(
                row for row in actual
                if row[0] == 0 and row[2] == 1 and row[3] == 6
            )
            require(
                selected == expected,
                f"target positive edge pair changed on rail {rail_index}, "
                f"label {label}: {selected}",
            )
            # Other unit rows can meet the bare deep band.  They do not share
            # the marked sector/height/future prefix and are not claimed as
            # edge choices of the displayed positive atom.
            require(
                set(expected).issubset(set(actual)),
                "selected target edge pair left the full compatible bank",
            )
            per_label.append(selected)
            per_label_counts.append(len(actual))
        positive_option_rows.append(tuple(per_label))
        positive_full_bank_counts.append(tuple(per_label_counts))
    require(
        len(set(target_zero_rows))
        == len(set(target_twelve_rows))
        == 1,
        "distinguished target rows stopped being rail-independent",
    )

    # Every edge half contains the full strict target band.  Enumerate the
    # entire Boolean cube, so the 2^11 claim is not inferred from a count.
    edge_assignment_count = 0
    for choices in product(range(2), repeat=11):
        common = TARGET_BAND
        for index, choice in enumerate(choices):
            half = positive_option_rows[0][index][choice][-1]
            meet = strict_meet(common, half)
            require(
                meet == TARGET_BAND,
                f"edge assignment {choices} cut the marked target band",
            )
            common = meet
        edge_assignment_count += 1
    require(
        edge_assignment_count == 2**11,
        "edge-choice cube cardinality changed",
    )

    # 3. The target-zero row dies at the same literal future coordinate.
    ordinary_prefixes = atlas.build_pair_prefixes(module)
    terminal_fork = marked.semantic.deepest_fork(module)
    marked_prefixes = marked.build_semantic_prefixes(
        module, terminal_fork
    )
    prefix_rows = []
    for clock in range(7):
        target_zero_prefix = prefix_intervals(
            ordinary_prefixes[1][clock][12][1]
        )
        target_marked_prefix = prefix_intervals(
            marked_prefixes[0][clock][1]
        )
        require(
            all(
                25 * T // 26 <= left < right <= T
                for left, right in target_zero_prefix
            ),
            "target-zero prefix escaped the strict digit (25,26)/26",
        )
        require(
            all(
                13 * T // 26 <= left < right <= 14 * T // 26
                for left, right in target_marked_prefix
            ),
            "marked target prefix escaped the strict digit (13,14)/26",
        )
        require(
            14 * T // 26 < 25 * T // 26
            and not interval_intersection(
                target_zero_prefix, target_marked_prefix
            ),
            f"target-zero delayed prefixes met at clock {clock}",
        )
        prefix_rows.append((
            clock,
            len(target_zero_prefix),
            len(target_marked_prefix),
            0,
        ))
    prefix_rows = tuple(prefix_rows)

    # 4. The target label 12 has integer lift 13 in source coordinates.
    # The wrong lift 0 is not physically equal, even though the carry agrees.
    residual = 13 * SHIFT
    require(
        residual == 91 * T // R
        and 0 < residual < T
        and residual % T != 0
        and C3 * residual == 14 * T,
        "integer-lift-13 odometer or C3 stabilizer changed",
    )

    wrap_order_rows = []
    wrong_sign_rows = []
    nonzero_lift_hostiles = []
    for clock in range(7):
        present_anchor = tuple(present[clock, 7])
        source_safe = clutch.complement(present_anchor)
        target_safe_pullback = clutch.shift_union(
            source_safe, -SHIFT
        )

        # Correct source-coordinate preimage under x -> x+13*tau.
        wrapped = clutch.shift_union(present_anchor, -residual)
        wrapped_weighted = tuple(
            (left, right, 1) for left, right in wrapped
        )
        after_source_safe = weighted_intersection(
            wrapped_weighted, source_safe
        )
        after_target_safe = weighted_intersection(
            after_source_safe, target_safe_pullback
        )
        after_deep = tuple(atlas.old.intersect_weighted_comb(
            after_target_safe, C3, 182, *SOURCE_BAND
        ))
        ordered = (
            weighted_mass(after_source_safe),
            weighted_mass(after_target_safe),
            weighted_mass(after_deep),
        )
        require(
            ordered == EXPECTED_WRAP_ORDER[clock],
            f"wrapped factor order changed at clock {clock}: {ordered}",
        )
        wrap_order_rows.append(ordered)

        # The opposite preimage sign also dies at the final C3 factor because
        # the 13-step shift is C3-trivial.  Its pre-deep target-safe mass is
        # different, so only the nonwrapping transport fixes orientation.
        wrong_wrapped = clutch.shift_union(present_anchor, residual)
        wrong_source = weighted_intersection(
            tuple((left, right, 1) for left, right in wrong_wrapped),
            source_safe,
        )
        wrong_target = weighted_intersection(
            wrong_source, target_safe_pullback
        )
        wrong_deep = tuple(atlas.old.intersect_weighted_comb(
            wrong_target, C3, 182, *SOURCE_BAND
        ))
        require(
            weighted_mass(wrong_source) == ordered[0]
            and weighted_mass(wrong_target) != ordered[1]
            and weighted_mass(wrong_deep) == 0,
            f"wrap-sign hostile stopped distinguishing clock {clock}",
        )
        wrong_sign_rows.append((
            weighted_mass(wrong_target),
            weighted_mass(wrong_deep),
        ))

        wrong_lift = tuple(present_anchor)
        symdiff_left = interval_intersection(
            wrapped, clutch.complement(wrong_lift)
        )
        symdiff_right = interval_intersection(
            wrong_lift, clutch.complement(wrapped)
        )
        symmetric_difference_mass = sum(
            right - left
            for left, right in symdiff_left + symdiff_right
        )
        require(
            symmetric_difference_mass > 0,
            f"integer lift 13 collapsed to lift 0 at clock {clock}",
        )
        nonzero_lift_hostiles.append(symmetric_difference_mass)
    wrap_order_rows = tuple(wrap_order_rows)

    # 5. Build the marked carriers and then reconstruct the target faces
    # directly in target coordinates, not by copying source coefficients.
    source_e3 = marked.semantic.exclusive_source(module, 3)
    sections = marked.semantic_sections(module, source_e3, 0, 4)
    positive_target_rows = []
    target_twelve_clock_rows = []
    marked_transport_checks = 0
    for rail_index in range(14):
        source_vector, target_vector, details = (
            marked.fully_marked_vectors(
                module,
                rails,
                present,
                marked_prefixes,
                sections,
                rail_index,
            )
        )
        require(
            source_vector == target_vector,
            f"marked source/target vector changed on rail {rail_index}",
        )
        for clock in range(7):
            require(
                normalize_weighted(clutch.shift_weighted(
                    details[clock][0], SHIFT
                ))
                == normalize_weighted(details[clock][1])
                and normalize_weighted(clutch.shift_weighted(
                    details[clock][2], SHIFT
                ))
                == normalize_weighted(details[clock][3]),
                f"marked carrier transport changed on rail {rail_index}, "
                f"clock {clock}",
            )
            marked_transport_checks += 2

        values = []
        full_values = []
        rail_support = weighted_support(rails[rail_index][3])
        for clock in range(7):
            common = tuple(
                row for row in details[clock][1] if row[2]
            )
            for label in range(1, 12):
                rail_preimage = clutch.shift_union(
                    rail_support, -label * SHIFT
                )
                present_preimage = clutch.shift_union(
                    tuple(present[clock, 7]), -label * SHIFT
                )
                common = weighted_intersection(
                    common, rail_preimage
                )
                common = weighted_intersection(
                    common, present_preimage
                )
                if not common:
                    break
            value = (
                atlas.delayed_carry_pair(
                    common, marked_prefixes[0][clock], {}
                )[TARGET_CARRY][1]
                if common else 0
            )
            values.append(value)

            # Add the sole remaining label directly in target coordinates.
            label = 12
            rail_preimage = clutch.shift_union(
                rail_support, -label * SHIFT
            )
            present_preimage = clutch.shift_union(
                tuple(present[clock, 7]), -label * SHIFT
            )
            common = weighted_intersection(common, rail_preimage)
            common = weighted_intersection(common, present_preimage)
            full_value = (
                atlas.delayed_carry_pair(
                    common, marked_prefixes[0][clock], {}
                )[TARGET_CARRY][1]
                if common else 0
            )
            full_values.append(full_value)

        values = tuple(values)
        full_values = tuple(full_values)
        row = (
            tuple(
                clock for clock, value in enumerate(values) if value
            ),
            sum(values),
        )
        require(
            row == EXPECTED_POSITIVE_ROWS[rail_index],
            f"direct target positive row changed on rail {rail_index}",
        )
        require(
            full_values == (0,) * 7,
            f"target labels 1,...,12 survived on rail {rail_index}",
        )
        positive_target_rows.append(row)
        target_twelve_clock_rows.append(full_values)
    positive_target_rows = tuple(positive_target_rows)
    target_twelve_clock_rows = tuple(target_twelve_clock_rows)
    require(
        marked_transport_checks == 14 * 7 * 2,
        "marked carrier/carry check universe changed",
    )

    # 6. Exhaust the target twelve-face combinatorics.  There are exactly
    # twelve faces containing target zero and one face omitting it.
    target_labels = set(range(P))
    faces = tuple(
        tuple(sorted(target_labels - {missing}))
        for missing in range(P)
    )
    containing_zero = tuple(face for face in faces if 0 in face)
    omitting_zero = tuple(face for face in faces if 0 not in face)
    require(
        len(faces) == 13
        and len(containing_zero) == 12
        and omitting_zero == (tuple(range(1, 13)),),
        "target twelve-face census changed",
    )

    # Only now inspect the primary artifact, and only as immutable bytes.
    primary_script = (
        COMP / "lrc14_source_target_nonwrap_odometer_sharp_eleven_thm2819.py"
    )
    primary_output = (
        ROOT / "05-knowledge/results/"
        "lrc14_source_target_nonwrap_odometer_sharp_eleven_thm2819.out"
    )
    primary_hashes = (
        hashlib.sha256(lf_bytes(primary_script)).hexdigest(),
        hashlib.sha256(lf_bytes(primary_output)).hexdigest(),
    )
    require(
        primary_hashes
        == (
            "58ca287ce8394f8008730781e5a2d851b466730ff08ffdaafc31fbbfb48d8a8c",
            "dfa30a53cfc17d00a5d6edc056952638d284614525e2c8fb1d983ba0ca528f7a",
        ),
        "THM2819 primary evidence hash changed",
    )

    print("THM2819 INDEPENDENT HOSTILE AUDIT")
    print(f"dependency_hashes={tuple(DEPENDENCIES.items())}")
    print(f"scales=(P={P},R={R},S={S},C3={C3},T={T},SHIFT={SHIFT})")
    print(f"carry_relabelling={tuple(carry_rows)}")
    print(
        f"nonwrap_checks=(rails={nonwrap_rail_checks},"
        f"present={nonwrap_present_checks},marked={marked_transport_checks})"
    )
    print(f"strict_target_half_types={tuple(sorted(geometric_types))}")
    print(f"unique_target_zero_row={target_zero_rows[0]}")
    print(f"unique_target_twelve_row={target_twelve_rows[0]}")
    print(
        f"bare_band_positive_bank_counts="
        f"{tuple(positive_full_bank_counts)}"
    )
    print(
        "selected_positive_edges=sector0,kappa1,h6; "
        "two adjacent halves per label"
    )
    print(f"target_zero_prefix_rows={prefix_rows}")
    print(
        "target_zero_gate=(25,26)/26 versus (13,14)/26; "
        "strict and disjoint"
    )
    print(
        f"integer_lift=(target12 -> source13,residual_grid={residual},"
        f"residual=91/R,C3_residual=14)"
    )
    print(f"wrap_factor_order={wrap_order_rows}")
    print(f"wrong_sign_target_safe_and_deep={tuple(wrong_sign_rows)}")
    print(
        f"lift13_vs_lift0_symdiff_masses={tuple(nonzero_lift_hostiles)}"
    )
    print(f"direct_positive_target_rows={positive_target_rows}")
    print(f"direct_target_1_to_12_rows={target_twelve_clock_rows}")
    print(f"edge_assignments_checked={edge_assignment_count}")
    print("twelve_face_census=12_target_zero_prefix+1_target12_C3_wrap")
    print("sharp_target_maximum=11_labels; simplex_dimension=10")
    print(
        "scope=first14 fully marked THM2749 target rails at carry6; "
        "no outside-rail, row, or LRC14 conclusion"
    )
    print(f"primary_hashes={primary_hashes}")
    print("VERDICT=ACCEPT FOR PROMOTION")
    print("ALL INDEPENDENT EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
