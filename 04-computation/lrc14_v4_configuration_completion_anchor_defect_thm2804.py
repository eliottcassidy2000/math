#!/usr/bin/env python3
"""Exact referee for THM-2804's V4 configuration completion.

On each of the first fourteen semantic rails, four specific
(sector, edge, kappa) triples form the even-parity V4.  Three vertices are
the maximal source-12 configurations from THM-2797; the fourth restores
private root 12 but has only eleven unit carries.  The script checks the
finite-group statement, the exact carry-support diamond, and the fourth
vertex's positive-but-anchor-dead semantic attachment.
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

# name -> (sector, edge, kappa, h)
V4_CONFIGS = (
    ("A", 0, 0, 0, 6),
    ("B", 0, 1, 1, 6),
    ("C", 1, 1, 0, 0),
    ("D", 1, 0, 1, 0),
)
U6 = tuple(carry for carry in range(P) if carry != 6)
U0 = tuple(carry for carry in range(P) if carry != 0)
U06 = tuple(carry for carry in range(P) if carry not in (0, 6))
EXPECTED_SUPPORTS = {
    "A": U6,
    "B": U6,
    "C": U0,
    "D": U06,
}
EXPECTED_ROOTS = {"A": 12, "B": 12, "C": 11, "D": 12}
EXPECTED_D_FACET_SUPPORTS = (
    (1,),
    (6,),
    (1, 2, 3),
    (0, 5),
    (2, 3),
    (0, 4, 6),
    (0, 2, 3),
    (0, 4, 5, 6),
    (0, 1, 2, 3),
    (0, 4, 5),
    (0, 1, 3),
    (4, 5),
    (0, 2),
    (4, 5, 6),
)
EXPECTED_MARKED_SUPPORTS = (
    (1,),
    (6,),
    (2, 3),
    (5,),
    (2, 3),
    (4, 6),
    (2, 3),
    (5, 6),
    (1, 2, 3),
    (5,),
    (1, 3),
    (5,),
    (2,),
    (5, 6),
)


def xor(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


def matrix_action(matrix, vector):
    return (
        (matrix[0][0] * vector[0] + matrix[0][1] * vector[1]) % 2,
        (matrix[1][0] * vector[0] + matrix[1][1] * vector[1]) % 2,
    )


def compose(left, right):
    """Composition left after right for permutations stored as images."""
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    out = [None] * len(permutation)
    for index, image in enumerate(permutation):
        out[image] = index
    return tuple(out)


def affine_permutations():
    """Enumerate AGL(2,2) on the four (sector,edge) points."""
    points = ((0, 0), (0, 1), (1, 0), (1, 1))
    directions = ((0, 1), (1, 0), (1, 1))
    linear = []
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    if (a * d - b * c) % 2:
                        linear.append(((a, b), (c, d)))
    permutations = set()
    translations = set()
    direction_permutations = set()
    for matrix in linear:
        direction_permutations.add(tuple(
            directions.index(matrix_action(matrix, direction))
            for direction in directions
        ))
        for shift in points:
            image = tuple(
                points.index(xor(matrix_action(matrix, point), shift))
                for point in points
            )
            permutations.add(image)
            if matrix == ((1, 0), (0, 1)):
                translations.add(image)
    return (
        tuple(sorted(permutations)),
        tuple(sorted(translations)),
        tuple(sorted(direction_permutations)),
    )


def support_of(weighted):
    return tuple((left, right) for left, right, weight in weighted if weight)


def intersect_support(weighted, support):
    if not weighted or not support:
        return tuple()
    return tuple(atlas.old.intersect_weighted_union(weighted, support))


def common_rail(rails, rail_index, deltas, shifted_rail):
    pieces = rails[rail_index][3]
    support = support_of(pieces)
    common = tuple(piece for piece in pieces if piece[2])
    for delta in deltas:
        if delta == 0:
            continue
        key = (rail_index, delta)
        if key not in shifted_rail:
            shifted_rail[key] = fixed.shift_union(
                support, 7 * delta * T // R, T
            )
        common = intersect_support(common, shifted_rail[key])
        if not common:
            break
    return tuple(common)


def d_clock_layers(module, present, rail_common, deltas, clock,
                   shifted_present):
    """Build the D vertex: (sector,edge,kappa,h)=(1,0,1,0)."""
    h = 0
    anchor_support = present[clock, 0]
    anchor = intersect_support(rail_common, anchor_support)
    all_present = anchor
    for delta in deltas:
        if delta == 0 or not all_present:
            continue
        key = (clock, delta)
        if key not in shifted_present:
            shifted_present[key] = fixed.shift_union(
                anchor_support, 7 * delta * T // R, T
            )
        all_present = intersect_support(
            all_present, shifted_present[key]
        )
    private_half = tuple(
        rebase.intersect_root_half(
            all_present, module.C3, 0, 12
        )
    )
    return rail_common, anchor, tuple(all_present), private_half


def first_failure(layers, marked_source):
    if not marked_source:
        return "marked-empty"
    marked_support = support_of(marked_source)
    for name, layer in zip(
        ("rail", "anchor", "all-present", "private-half"), layers
    ):
        if not intersect_support(layer, marked_support):
            return name
    return "survives"


def open_overlap_and_seams(left_weighted, right_weighted):
    left = support_of(left_weighted)
    right = support_of(right_weighted)
    i = 0
    j = 0
    overlap = 0
    forward = 0
    reverse = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        if max(a, c) < min(b, d):
            overlap += min(b, d) - max(a, c)
        elif b == c:
            forward += 1
        elif d == a:
            reverse += 1
        if b <= d:
            i += 1
        else:
            j += 1
    return overlap, forward, reverse


def main():
    shard = atlas.shard((0, RAIL_COUNT))
    require(
        shard[1] == 26
        and len(shard[5]) == len(shard[6]) == RAIL_COUNT,
        "fourteen-rail shard changed",
    )
    metadata = shard[5]
    flags = fixed.build_flags(shard[6], shard[1])

    triples = tuple(
        (sector, edge, kappa)
        for _name, sector, edge, kappa, _h in V4_CONFIGS
    )
    require(
        set(triples) == {(0, 0, 0), (0, 1, 1),
                         (1, 1, 0), (1, 0, 1)}
        and all((sector ^ edge ^ kappa) == 0
                for sector, edge, kappa in triples),
        "configuration triples stopped being the even-parity V4",
    )
    require(
        all(h == 6 * (1 - sector)
            for _name, sector, _edge, _kappa, h in V4_CONFIGS),
        "height stopped being the sector function",
    )

    affine, translations, direction_permutations = affine_permutations()
    require(
        len(affine) == 24
        and len(translations) == 4
        and all(
            tuple(sorted(permutation)) == (0, 1, 2, 3)
            for permutation in affine
        ),
        "AGL(2,2) or its V4 translation subgroup changed",
    )
    require(
        len(direction_permutations) == 6
        and all(
            tuple(sorted(permutation)) == (0, 1, 2)
            for permutation in direction_permutations
        )
        and all(
            compose(compose(group_element, translation),
                    inverse(group_element)) in translations
            for group_element in affine
            for translation in translations
        ),
        "translation V4 stopped being normal or the quotient stopped "
        "acting as S3 on the three directions",
    )
    quotient_order = len(affine) // len(translations)
    require(
        quotient_order == 6,
        "AGL(2,2)/V4 stopped having S3 order",
    )

    support_rows = []
    for rail_index in range(RAIL_COUNT):
        row = []
        for name, sector, edge, kappa, h in V4_CONFIGS:
            carries = tuple(
                carry for carry in range(P)
                if flags[rail_index][sector][edge][carry][kappa][h]
            )
            root = (
                2 * SOURCE_CARRY
                + (2 * h + kappa) // P
                + (edge == 0)
            ) % P
            require(
                carries == EXPECTED_SUPPORTS[name]
                and root == EXPECTED_ROOTS[name],
                f"V4 support/root law changed on rail {rail_index}, "
                f"vertex {name}",
            )
            row.append((name, carries, root))
        support_rows.append(tuple(row))

    require(
        set(EXPECTED_SUPPORTS["D"])
        == set(EXPECTED_SUPPORTS["A"])
        & set(EXPECTED_SUPPORTS["C"])
        == set(EXPECTED_SUPPORTS["B"])
        & set(EXPECTED_SUPPORTS["C"])
        and set(EXPECTED_SUPPORTS["A"])
        | set(EXPECTED_SUPPORTS["C"]) == set(range(P)),
        "carry-support meet/union diamond changed",
    )
    capacities = tuple(
        len(EXPECTED_SUPPORTS[name])
        for name, *_rest in V4_CONFIGS
    )
    require(
        capacities == (12, 12, 12, 11),
        "capacity stopped singling out the D vertex",
    )
    missing_labels = tuple(
        tuple(sorted(
            (2 * (carry - SOURCE_CARRY)) % P
            for carry in range(P)
            if carry not in EXPECTED_SUPPORTS[name]
        ))
        for name, *_rest in V4_CONFIGS
    )
    require(
        missing_labels == ((1,), (1,), (2,), (1, 2)),
        "source-twelve missing-label diamond changed",
    )

    module, _, _, _, rails, present, _starts = (
        atlas.core.build_carrier_data()
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
        require(source == target,
                f"marked vector changed on rail {rail_index}")
        marked_bank.append((source, details))

    d_deltas = tuple(
        (2 * (carry - SOURCE_CARRY)) % P for carry in U06
    )
    require(
        tuple(sorted(d_deltas))
        == tuple(label for label in range(P) if label not in (1, 2)),
        "D vertex active-label set changed",
    )

    shifted_rail = {}
    shifted_present = {}
    delayed_cache = {}
    failure_hist = Counter()
    facet_supports = []
    marked_supports = []
    common_supports = []
    total_facet_mass = 0
    total_overlap = 0
    total_forward_seams = 0
    total_reverse_seams = 0
    seam_rows = []

    for rail_index in range(RAIL_COUNT):
        rail_common = common_rail(
            rails, rail_index, d_deltas, shifted_rail
        )
        require(rail_common,
                f"D common rail vanished on rail {rail_index}")
        facet_vector = []
        marked_vector, details = marked_bank[rail_index]
        rail_forward = 0
        rail_reverse = 0
        for clock in range(7):
            layers = d_clock_layers(
                module,
                present,
                rail_common,
                d_deltas,
                clock,
                shifted_present,
            )
            marked_source = tuple(
                piece for piece in details[clock][0] if piece[2]
            )
            failure_hist[
                first_failure(layers, marked_source)
            ] += 1
            values = atlas.delayed_carry_pair(
                layers[3],
                prefixes[1][clock][0],
                delayed_cache.setdefault(clock, {}),
            )
            facet_vector.append(values[SOURCE_CARRY][1])
            overlap, forward, reverse = open_overlap_and_seams(
                layers[3], marked_source
            )
            total_overlap += overlap
            total_forward_seams += forward
            total_reverse_seams += reverse
            rail_forward += forward
            rail_reverse += reverse

        facet_support = tuple(
            clock for clock, value in enumerate(facet_vector) if value
        )
        marked_support = tuple(
            clock for clock, value in enumerate(marked_vector) if value
        )
        common_support = tuple(sorted(
            set(facet_support) & set(marked_support)
        ))
        facet_supports.append(facet_support)
        marked_supports.append(marked_support)
        common_supports.append(common_support)
        total_facet_mass += sum(facet_vector)
        seam_rows.append((rail_index, rail_forward, rail_reverse))

    facet_supports = tuple(facet_supports)
    marked_supports = tuple(marked_supports)
    common_supports = tuple(common_supports)
    require(
        facet_supports == EXPECTED_D_FACET_SUPPORTS
        and marked_supports == EXPECTED_MARKED_SUPPORTS
        and all(common_supports),
        "D positive/common clock controls changed",
    )
    require(
        tuple(sorted(failure_hist.items()))
        == (("anchor", 23), ("marked-empty", 75))
        and total_overlap == 0,
        "D semantic first-failure law changed",
    )
    require(
        total_facet_mass == 103_771_702_146_779_222_198_643_000
        and total_forward_seams == total_reverse_seams == 0
        and tuple(seam_rows)
        == tuple((rail, 0, 0) for rail in range(RAIL_COUNT)),
        "D mass or strict-separation boundary changed",
    )

    print("THM2804 FOURTEEN-RAIL V4 CONFIGURATION COMPLETION EXACT REFEREE")
    print(f"v4_configs={V4_CONFIGS}")
    print("v4_bits=even_parity_subgroup_of_F2^3")
    print("height_law=h=6*(1-sector)")
    print(f"agl22_order={len(affine)} translation_v4_order="
          f"{len(translations)} abstract_quotient_order={quotient_order}")
    print(f"carry_capacities_ABCD={capacities}")
    print(f"carry_supports_ABCD={tuple(EXPECTED_SUPPORTS[name] for name,*_ in V4_CONFIGS)}")
    print(f"source12_private_roots_ABCD={tuple(EXPECTED_ROOTS[name] for name,*_ in V4_CONFIGS)}")
    print(f"missing_target_labels_ABCD={missing_labels}")
    print("carry_diamond=U_D=U_A_intersect_U_C="
          "U_B_intersect_U_C; U_A_union_U_C=F13")
    print("physical_asymmetry=D alone has capacity11; abstract "
          "Aut(V4)=S3 does not act as semantic transport")
    print(f"d_facet_supports={facet_supports}")
    print(f"marked_supports={marked_supports}")
    print(f"d_common_positive_supports={common_supports}")
    print(f"d_first_failure_hist={tuple(sorted(failure_hist.items()))}")
    print(f"d_total_source12_facet_mass={total_facet_mass}")
    print(f"d_seam_rows={tuple(seam_rows)}")
    print("d_total_overlap=0 d_forward_seams="
          f"{total_forward_seams} d_reverse_seams={total_reverse_seams}")
    print("scope=abstract V4/S3 boundary plus one eleven-label D packet "
          "on first14 rails; no physical S3 action, labelwise union, row, "
          "or LRC14 exclusion")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
