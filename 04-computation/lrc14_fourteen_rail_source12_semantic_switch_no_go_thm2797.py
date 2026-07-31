#!/usr/bin/env python3
"""Exact referee for the THM-2797 fourteen-rail semantic-switch no-go.

The universe is the first fourteen rails in the proved THM-2749 semantic
bank.  On each rail we enumerate every maximal THM-2672 configuration whose
twelve-carry unit set contains source carry 12, rebase its twelve labels to
that source carry, and compare its strict pre-prefix facet base with the
clock-matched fully marked semantic source base on the same rail.

The comparison is configurationwise.  It does not form a union over
configurations and does not claim a row or LRC(14) exclusion.
"""

from collections import Counter
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
import lrc14_slope7_rebase_facet_torsor_thm2672 as rebase
import lrc14_fully_marked_root_zero_clutch_thm2749 as marked


P = 13
R = P**6
T = atlas.T
SOURCE_CARRY = 12
RAIL_COUNT = 14

EXPECTED_COMMON_SUPPORT_HIST = (
    ((1,), 3),
    ((1, 2, 3), 3),
    ((1, 3), 3),
    ((2,), 3),
    ((2, 3), 9),
    ((4, 6), 3),
    ((5,), 11),
    ((5, 6), 4),
    ((6,), 3),
)
EXPECTED_SEAM_HIST = (
    (0, 28),
    (51, 1),
    (278, 1),
    (280, 1),
    (393, 1),
    (394, 1),
    (406, 1),
    (416, 1),
    (473, 1),
    (477, 1),
    (504, 1),
    (659, 1),
    (790, 1),
    (886, 1),
    (1011, 1),
)
EXPECTED_CONFIG_SIGNATURES = (
    (0, 0, 0, 6, 6),
    (0, 1, 1, 6, 6),
    (1, 1, 0, 0, 0),
)


def intersect_support(weighted, support):
    if not weighted or not support:
        return tuple()
    return tuple(atlas.old.intersect_weighted_union(weighted, support))


def support_of(weighted):
    return tuple((left, right) for left, right, weight in weighted if weight)


def open_overlap_and_seams(left_weighted, right_weighted):
    """Return open-overlap numerator and both oriented seam counts."""
    left = support_of(left_weighted)
    right = support_of(right_weighted)
    require(
        all(left[index][1] <= left[index + 1][0]
            for index in range(len(left) - 1)),
        "left interval bank is not sorted/disjoint",
    )
    require(
        all(right[index][1] <= right[index + 1][0]
            for index in range(len(right) - 1)),
        "right interval bank is not sorted/disjoint",
    )
    i = 0
    j = 0
    overlap = 0
    facet_right_to_marked_left = 0
    marked_right_to_facet_left = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        low = max(a, c)
        high = min(b, d)
        if low < high:
            overlap += high - low
        elif b == c:
            facet_right_to_marked_left += 1
        elif d == a:
            marked_right_to_facet_left += 1
        if b <= d:
            i += 1
        else:
            j += 1
    return (
        overlap,
        facet_right_to_marked_left,
        marked_right_to_facet_left,
    )


def common_rail(rails, rail_index, deltas, cache):
    pieces = rails[rail_index][3]
    support = support_of(pieces)
    common = tuple(piece for piece in pieces if piece[2])
    for delta in deltas:
        if delta == 0:
            continue
        key = (rail_index, delta)
        if key not in cache:
            cache[key] = fixed.shift_union(
                support, 7 * delta * T // R, T
            )
        common = intersect_support(common, cache[key])
        if not common:
            break
    return tuple(common)


def clock_layers(module, present, rail_common, deltas, edge, kappa, h,
                 clock, shifted_present):
    """Return rail, anchor-present, all-present, and private-half layers."""
    anchor_support = present[clock, (-h) % P]
    anchor = intersect_support(rail_common, anchor_support)
    all_present = anchor
    for delta in deltas:
        if delta == 0 or not all_present:
            continue
        key = (clock, h, delta)
        if key not in shifted_present:
            shifted_present[key] = fixed.shift_union(
                anchor_support, 7 * delta * T // R, T
            )
        all_present = intersect_support(
            all_present, shifted_present[key]
        )
    root = (
        2 * SOURCE_CARRY
        + (2 * h + kappa) // P
        + (edge == 0)
    ) % P
    private_half = tuple(
        rebase.intersect_root_half(
            all_present, module.C3, edge, root
        )
    )
    return rail_common, anchor, tuple(all_present), private_half, root


def first_failure(layers, marked_source):
    """Name the first facet layer disjoint from a nonempty marked base."""
    if not marked_source:
        return "marked-empty"
    names = ("rail", "anchor", "all-present", "private-half")
    marked_support = support_of(marked_source)
    for name, layer in zip(names, layers):
        if not intersect_support(layer, marked_support):
            return name
    return "survives"


def main():
    require(
        (atlas.P, atlas.R, T) == (13, 4_826_809, 297_836_897_838_480),
        "canonical THM2672 scales changed",
    )

    shard = atlas.shard((0, RAIL_COUNT))
    content = shard[1]
    metadata = shard[5]
    rows = shard[6]
    require(
        content == 26
        and len(metadata) == len(rows) == RAIL_COUNT,
        "fourteen-rail shard changed",
    )
    flags = fixed.build_flags(rows, content)

    module, _, _, _, rails, present, _starts = (
        atlas.core.build_carrier_data()
    )
    require(
        tuple((sector, clock, theta)
              for sector, clock, theta, _pieces in rails[:RAIL_COUNT])
        == metadata,
        "fourteen-rail metadata order changed",
    )
    prefixes = atlas.build_pair_prefixes(module)

    source_e3 = marked.semantic.exclusive_source(module, 3)
    terminal_fork = marked.semantic.deepest_fork(module)
    semantic_prefixes = marked.build_semantic_prefixes(
        module, terminal_fork
    )
    sections = marked.semantic_sections(module, source_e3, 0, 4)
    marked_bank = []
    for rail_index in range(RAIL_COUNT):
        source, target, details = marked.fully_marked_vectors(
            module,
            rails,
            present,
            semantic_prefixes,
            sections,
            rail_index,
        )
        require(
            source == target,
            f"marked source/target vectors differ on rail {rail_index}",
        )
        marked_bank.append((source, details))

    shifted_rail = {}
    shifted_present = {}
    delayed_cache = {}

    maximal_count = 0
    admissible_count = 0
    positive_facet_count = 0
    common_positive_clock_count = 0
    overlap_clock_count = 0
    missing_hist = Counter()
    rail_hist = Counter()
    height_hist = Counter()
    root_hist = Counter()
    failure_hist = Counter()
    common_positive_support_hist = Counter()
    seam_hist = Counter()
    total_facet_mass = 0
    total_facet_pieces = 0
    total_marked_pieces = 0
    total_facet_right_seams = 0
    total_marked_right_seams = 0
    config_rows = []

    for rail_index in range(RAIL_COUNT):
        marked_vector, marked_details = marked_bank[rail_index]
        for sector in range(2):
            for edge in range(2):
                for kappa in range(2):
                    for h in range(P):
                        carries = tuple(
                            carry for carry in range(P)
                            if flags[rail_index][sector][edge]
                                    [carry][kappa][h]
                        )
                        if len(carries) != 12:
                            continue
                        maximal_count += 1
                        if SOURCE_CARRY not in carries:
                            continue
                        admissible_count += 1
                        missing = next(
                            carry for carry in range(P)
                            if carry not in carries
                        )
                        missing_hist[missing] += 1
                        rail_hist[rail_index] += 1
                        height_hist[h] += 1
                        deltas = tuple(
                            (2 * (carry - SOURCE_CARRY)) % P
                            for carry in carries
                        )
                        require(
                            len(set(deltas)) == 12 and 0 in deltas,
                            "source-twelve label rebase changed",
                        )
                        rail_common = common_rail(
                            rails,
                            rail_index,
                            deltas,
                            shifted_rail,
                        )
                        require(
                            rail_common,
                            "a maximal source-twelve common rail vanished",
                        )

                        facet_vector = []
                        overlap_vector = []
                        config_seams = 0
                        for clock in range(7):
                            layers = clock_layers(
                                module,
                                present,
                                rail_common,
                                deltas,
                                edge,
                                kappa,
                                h,
                                clock,
                                shifted_present,
                            )
                            private_half = layers[3]
                            root_hist[layers[4]] += int(bool(private_half))
                            marked_source = tuple(
                                piece
                                for piece in marked_details[clock][0]
                                if piece[2]
                            )
                            failure_hist[
                                first_failure(layers[:4], marked_source)
                            ] += 1
                            total_facet_pieces += len(private_half)
                            total_marked_pieces += len(marked_source)

                            key = (sector, clock, h)
                            values = atlas.delayed_carry_pair(
                                private_half,
                                prefixes[sector][clock][h],
                                delayed_cache.setdefault(key, {}),
                            )
                            facet_value = values[SOURCE_CARRY][kappa]
                            facet_vector.append(facet_value)

                            overlap, forward, reverse = (
                                open_overlap_and_seams(
                                    private_half, marked_source
                                )
                            )
                            overlap_vector.append(overlap)
                            overlap_clock_count += int(overlap > 0)
                            total_facet_right_seams += forward
                            total_marked_right_seams += reverse
                            config_seams += forward + reverse

                        facet_vector = tuple(facet_vector)
                        overlap_vector = tuple(overlap_vector)
                        facet_support = tuple(
                            index for index, value
                            in enumerate(facet_vector) if value
                        )
                        marked_support = tuple(
                            index for index, value
                            in enumerate(marked_vector) if value
                        )
                        common_support = tuple(sorted(
                            set(facet_support) & set(marked_support)
                        ))
                        positive_facet_count += int(bool(facet_support))
                        common_positive_clock_count += int(
                            bool(common_support)
                        )
                        common_positive_support_hist[common_support] += 1
                        seam_hist[config_seams] += 1
                        total_facet_mass += sum(facet_vector)
                        config_rows.append((
                            rail_index,
                            metadata[rail_index],
                            sector,
                            edge,
                            kappa,
                            h,
                            missing,
                            facet_support,
                            marked_support,
                            common_support,
                            config_seams,
                        ))
                        require(
                            not any(overlap_vector),
                            "a source-twelve facet acquired a marked-base "
                            "open overlap",
                        )

    require(
        overlap_clock_count == 0,
        "a strict semantic-base overlap survived",
    )
    require(
        failure_hist["survives"] == 0,
        "a marked source survived every facet layer",
    )
    require(
        (
            maximal_count,
            admissible_count,
            positive_facet_count,
            common_positive_clock_count,
            overlap_clock_count,
        ) == (55, 42, 42, 42, 0),
        "fourteen-rail configuration census changed",
    )
    require(
        tuple(sorted(missing_hist.items())) == ((0, 14), (6, 28))
        and tuple(sorted(rail_hist.items()))
        == tuple((rail, 3) for rail in range(RAIL_COUNT))
        and tuple(sorted(height_hist.items())) == ((0, 14), (6, 28)),
        "configuration stratum census changed",
    )
    require(
        tuple(sorted(common_positive_support_hist.items()))
        == EXPECTED_COMMON_SUPPORT_HIST
        and tuple(sorted(failure_hist.items()))
        == (("anchor", 69), ("marked-empty", 225))
        and tuple(sorted(seam_hist.items())) == EXPECTED_SEAM_HIST,
        "support, first-failure, or seam census changed",
    )
    require(
        (
            total_facet_right_seams,
            total_marked_right_seams,
            total_facet_pieces,
            total_marked_pieces,
            total_facet_mass,
        ) == (
            7018,
            0,
            45449,
            30534,
            290_130_506_081_396_373_729_700_182,
        ),
        "aggregate piece, seam, or mass census changed",
    )
    by_rail = []
    for rail_index in range(RAIL_COUNT):
        rows_on_rail = tuple(
            row for row in config_rows if row[0] == rail_index
        )
        require(
            tuple(row[2:7] for row in rows_on_rail)
            == EXPECTED_CONFIG_SIGNATURES,
            f"configuration signatures changed on rail {rail_index}",
        )
        require(
            rows_on_rail[0][-1] == rows_on_rail[2][-1] == 0
            and rows_on_rail[1][-1] > 0,
            f"seam-bearing configuration changed on rail {rail_index}",
        )
        by_rail.append((
            rail_index,
            metadata[rail_index],
            tuple(
                (
                    row[7],
                    row[8],
                    row[9],
                    row[10],
                )
                for row in rows_on_rail
            ),
        ))
    by_rail = tuple(by_rail)

    print("THM2797 FOURTEEN-RAIL SOURCE-TWELVE SEMANTIC-SWITCH EXACT REFEREE")
    print(f"maximal_configs_first14={maximal_count}")
    print(f"source12_admissible_configs={admissible_count}")
    print(f"positive_source12_facets={positive_facet_count}")
    print("configs_with_common_positive_clock="
          f"{common_positive_clock_count}")
    print(f"missing_carry_hist={tuple(sorted(missing_hist.items()))}")
    print(f"rail_hist={tuple(sorted(rail_hist.items()))}")
    print(f"height_hist={tuple(sorted(height_hist.items()))}")
    print("common_positive_support_hist="
          f"{tuple(sorted(common_positive_support_hist.items()))}")
    print(f"first_failure_hist={tuple(sorted(failure_hist.items()))}")
    print(f"seam_hist={tuple(sorted(seam_hist.items()))}")
    print("total_facet_right_to_marked_left_seams="
          f"{total_facet_right_seams}")
    print("total_marked_right_to_facet_left_seams="
          f"{total_marked_right_seams}")
    print(f"total_facet_base_pieces={total_facet_pieces}")
    print(f"total_marked_base_pieces={total_marked_pieces}")
    print(f"total_source12_facet_mass={total_facet_mass}")
    print("config_signatures_each_rail="
          f"{EXPECTED_CONFIG_SIGNATURES}")
    for row in by_rail:
        print(f"rail_support_and_seam_row={row}")
    print("strict_open_overlap_clocks=0")
    print("scope=first14 matching semantic rails and each fixed maximal "
          "source12 configuration separately; no configuration union, "
          "row, or LRC14 exclusion")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
