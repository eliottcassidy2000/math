#!/usr/bin/env python3
"""Independent exact hostile audit of the proposed THM-2825 collar.

This script deliberately does not import the THM-2806/2818/2825 theorem
companions.  It rebuilds the full 7 x 9 x 9 rail-eight bank from the pinned
lower physical clutch module, constructs weighted common and pulled-target
right-wing pieces, and audits the half-step metric, delayed content, and
factor/carrier/ancestry sidecars.

The key distinction is between the complete unscaled relation (all M x R
pairs in a cell) and the scale-selected map R -> M obtained from centered
circular distance.  No Python assert statement is used.
"""

from bisect import bisect_right
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
    "lrc14_full_arm_orbit_path_sheet_audit_thm2791.py":
        "1e00b6711db0d878fca70047f5f1532518084dbf6928551cd28fe51283fde543",
}
for name, expected in PINNED.items():
    payload = (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"lower physical dependency changed: {name}",
    )


import lrc14_full_arm_orbit_path_sheet_audit_thm2791 as ancestry
import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = 13
Q = 7
T = physical.T
SHIFT = physical.SHIFT
H = T // (2 * P**5)
WRAP = P * SHIFT
LENGTH = 26_444_880
CONTENT = 103_478_815_440
W1 = 27_581_135_604
W2 = 27_580_222_516
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = tuple(range(3, 12))
FACTOR_NAMES = ("E3", "clock", "q1", "q2", "c2", "c3")
EXPECTED_TYPES = {
    (1, 189, 3): 54,
    (1, 213, 2): 18,
    (1, 239, 2): 9,
    (2, 190, 3): 7,
    (2, 426, 4): 28,
    (2, 474, 2): 14,
    (2, 526, 2): 7,
    (3, 382, 5): 7,
    (3, 404, 4): 28,
    (3, 409, 3): 7,
    (3, 452, 2): 7,
    (3, 504, 2): 7,
}


def weighted_intersection(pieces, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            pieces, intervals
        )
    )


def support_of(pieces):
    return physical.overlap.merge_intervals(
        (left, right) for left, right, weight in pieces if weight
    )


def intersect_sorted(left, right):
    rows = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            rows.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(rows)


def build_geometry():
    (
        module,
        _prefixes,
        _whole,
        _masses,
        rails,
        present,
        _starts,
    ) = physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _source_vector, _target_vector, rail_pairs, details = (
        physical.overlap.overlap_vectors(
            module, pair_prefixes, rails, present, rail_index=8
        )
    )
    full_module = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full_module, 3)
    fork = physical.target.deepest_fork(full_module)
    clocks = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(Q)
    )
    q_pairs = physical.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )
    require(
        tuple(rails[8][:3]) == (1, 4, 12)
        and all(a == b for _left, _right, a, b in rail_pairs),
        "rail-eight equal-weight base changed",
    )
    return module, full_module, details, e3, clocks, q_pairs


def reconstruct_cell(details, full_module, e3, clocks, clock, s, t):
    section = physical.target.source_present_section(
        full_module, e3, clock, s, t, clocks
    )
    source_base, target_base = details[clock]
    source = weighted_intersection(source_base, section)
    target = weighted_intersection(target_base, section)
    target_pullback = physical.overlap.shift_weighted(target, -SHIFT)
    aligned = physical.overlap.intersect_weighted_profiles(
        source, target_pullback
    )
    require(
        all(a == b for _left, _right, a, b in aligned),
        f"common weights differ at {(clock, s, t)}",
    )
    common = tuple(
        (left, right, a) for left, right, a, _b in aligned
    )
    right = physical.subtract_weighted(target_pullback, common)
    require(
        not intersect_sorted(support_of(common), support_of(right)),
        f"common/right supports intersect at {(clock, s, t)}",
    )
    require(
        physical.merge_weighted(common + right)
        == physical.merge_weighted(target_pullback),
        f"common/right decomposition missed target mass at {(clock, s, t)}",
    )
    return source, target_pullback, common, right


def semantic_pair(interval, q_pair, source_cache, target_cache):
    unit = ((interval[0], interval[1], 1),)
    target_unit = physical.overlap.shift_weighted(unit, SHIFT)
    source_value = physical.relative.private.delayed_carry_pair(
        unit, q_pair, source_cache
    )[12][1]
    target_value = physical.relative.private.delayed_carry_pair(
        target_unit, q_pair, target_cache
    )[6][1]
    require(
        source_value == target_value
        and source_value in (0, CONTENT),
        f"piece has mixed or noncanonical delayed content: {interval}",
    )
    return source_value


def shifted_piece(piece, amount):
    left, right, weight = piece
    new_left = (left + amount) % T
    require(
        new_left + (right - left) <= T,
        "tested collar crossed the circle seam",
    )
    return new_left, new_left + (right - left), weight


def centered_distance(offset, circumference):
    return min(offset, circumference - offset)


def index_one(intervals):
    return intervals, tuple(right for _left, right in intervals)


def local_profile(interval, indexed):
    intervals, rights = indexed
    left, right = interval
    index = bisect_right(rights, left)
    rows = []
    while index < len(intervals) and intervals[index][0] < right:
        a = max(left, intervals[index][0])
        b = min(right, intervals[index][1])
        if a < b:
            rows.append((a - left, b - left))
        index += 1
    return tuple(rows)


def contained(interval, indexed):
    return local_profile(interval, indexed) == (
        (0, interval[1] - interval[0]),
    )


def meets(interval, indexed):
    return bool(local_profile(interval, indexed))


def section_factors(full_module, e3, clocks, clock, s, t):
    universe = ((0, T),)
    return (
        e3,
        clocks[clock],
        full_module.subtract_comb(
            universe,
            full_module.W[1],
            182,
            -14 * s - 13,
            -14 * s + 13,
        ),
        full_module.subtract_comb(
            universe,
            full_module.W[2],
            182,
            -14 * t - 13,
            -14 * t + 13,
        ),
        full_module.subtract_comb(
            universe,
            full_module.C2,
            182,
            14 * s - 13,
            14 * s + 13,
        ),
        full_module.subtract_comb(
            universe,
            full_module.C3,
            182,
            14 * t - 13,
            14 * t + 13,
        ),
    )


def indexed_factor_faces(factors):
    return (
        tuple(index_one(factor) for factor in factors),
        tuple(index_one(
            physical.overlap.shift_union(factor, -SHIFT)
        ) for factor in factors),
        tuple(index_one(factor) for factor in factors),
        tuple(index_one(
            physical.overlap.shift_union(factor, SHIFT)
        ) for factor in factors),
    )


def factor_signature(piece, faces):
    interval = piece[:2]
    target_interval = (
        interval[0] + SHIFT,
        interval[1] + SHIFT,
    )
    return (
        tuple(contained(interval, factor) for factor in faces[0]),
        tuple(contained(interval, factor) for factor in faces[1]),
        tuple(contained(target_interval, factor)
              for factor in faces[2]),
        tuple(contained(target_interval, factor)
              for factor in faces[3]),
    )


def indexed_carrier_twists(source, target_pullback):
    unit = T // P
    source_support = support_of(source)
    target_support = support_of(target_pullback)
    return (
        tuple(index_one(physical.overlap.shift_union(
            source_support, -twist * unit
        )) for twist in range(P)),
        tuple(index_one(physical.overlap.shift_union(
            target_support, -twist * unit
        )) for twist in range(P)),
    )


def carrier_signature(piece, twists):
    interval = piece[:2]
    return (
        tuple(meets(interval, carrier) for carrier in twists[0]),
        tuple(meets(interval, carrier) for carrier in twists[1]),
    )


def ancestry_setup():
    e_set = tuple(ancestry.base.build_set(
        ancestry.base.PAT_E, ancestry.base.ZELL
    ))
    q_set = tuple(ancestry.base.build_set(
        ancestry.host.PAT_QB, ancestry.base.ZELL
    ))
    u_events = tuple(sorted(set(
        ancestry.mapped_events(q_set, P**5)
        + ancestry.mapped_events(e_set, P**7)
    )))
    v_events = ancestry.mapped_events(e_set, P**5, T // P)
    return u_events, v_events


def ancestry_chamber(piece, u_events, v_events):
    left, right = piece[:2]
    u_index = bisect_right(u_events, left)
    v_index = bisect_right(v_events, left)
    require(
        0 < u_index < len(u_events)
        and 0 < v_index < len(v_events)
        and u_events[u_index] >= right
        and v_events[v_index] >= right,
        "an interval is cut by an ancestry wall",
    )
    return u_index, v_index


def main():
    require(
        (T, SHIFT, H, WRAP, LENGTH)
        == (
            297_836_897_838_480,
            431_933_040,
            401_080_680,
            5_615_129_520,
            26_444_880,
        )
        and LENGTH * 91 == 6 * H
        and WRAP == 14 * H == 7 * (2 * H),
        "scale identities changed",
    )
    module, full_module, details, e3, clocks, q_pairs = build_geometry()
    u_events, v_events = ancestry_setup()
    ancestry_chambers = set()
    semantic_cache = {}
    source_phi_cache = [{} for _clock in range(Q)]
    target_phi_cache = [{} for _clock in range(Q)]

    def content(clock, piece):
        key = (clock, piece[0], piece[1])
        if key not in semantic_cache:
            semantic_cache[key] = semantic_pair(
                piece,
                q_pairs[clock],
                source_phi_cache[clock],
                target_phi_cache[clock],
            )
        return semantic_cache[key]

    empty = nonempty = total_common = total_right = total_pairs = 0
    cell_types = Counter()
    content_pair_counts = Counter()
    right_content = Counter()
    plus_one_content = Counter()
    plus_two_content = Counter()
    wrap_landing = Counter()
    wrap_content = Counter()
    unscaled_choice_histogram = Counter()
    reverse_nearest_histogram = Counter()
    reverse_unique = reverse_tied = 0
    factor_equal = Counter()
    carrier_equal = Counter()
    ancestry_equal = Counter()
    first_factor_hostile = None
    first_carrier_hostile = None
    first_ancestry_hostile = None
    forward_targets = defaultdict(set)
    plus_two_targets = defaultdict(set)
    wrap_targets = defaultdict(set)
    sidecar_cells = []
    sidecar_clocks_done = set()
    parity_type_histogram = Counter()
    parity_unit_cells = Counter()
    signed_zero_cells = defaultdict(list)
    total_zero_cells = defaultdict(list)
    full_c2_unit_cells = Counter()

    for clock in range(Q):
        for s in COMMON_S:
            for t in COMMON_T:
                source, target_pullback, common, right = reconstruct_cell(
                    details, full_module, e3, clocks, clock, s, t
                )
                cell = (clock, s, t)
                if not common and not right:
                    empty += 1
                    continue
                require(
                    common and right,
                    f"one-sided nonempty cell appeared at {cell}",
                )
                nonempty += 1
                total_common += len(common)
                total_right += len(right)
                total_pairs += len(common) * len(right)
                cell_types[clock, len(common), len(right)] += 1
                unscaled_choice_histogram[len(common)] += len(right)

                expected_weight = W1 if clock == 1 else W2
                require(
                    {
                        piece[1] - piece[0]
                        for piece in common + right
                    } == {LENGTH}
                    and {
                        piece[2] for piece in common + right
                    } == {expected_weight}
                    and len({
                        piece[0] % H for piece in common + right
                    }) == 1,
                    f"length/weight/half-step lattice changed at {cell}",
                )
                parity_base = min(
                    piece[0] for piece in common + right
                )

                def parity_counts(pieces):
                    counts = [0, 0]
                    for piece in pieces:
                        require(
                            (piece[0] - parity_base) % H == 0,
                            f"parity base missed half-step lattice at {cell}",
                        )
                        counts[
                            ((piece[0] - parity_base) // H) % 2
                        ] += 1
                    return tuple(counts)

                m_parity = parity_counts(common)
                r_parity = parity_counts(right)
                u_parity = tuple(
                    m_parity[index] + r_parity[index]
                    for index in range(2)
                )
                parity_type_histogram[
                    clock, m_parity, r_parity
                ] += 1
                piece_mass_mod = (
                    expected_weight * LENGTH
                ) % P
                require(piece_mass_mod != 0,
                        "one physical copy ceased to be a 13-unit")
                for object_name, counts in (
                    ("M", m_parity),
                    ("R", r_parity),
                    ("M+R", u_parity),
                ):
                    colour_units = tuple(
                        count % P != 0 for count in counts
                    )
                    parity_unit_cells[
                        object_name, colour_units
                    ] += 1
                    total_character = (
                        (counts[0] + counts[1]) * piece_mass_mod
                    ) % P
                    signed_character = (
                        (counts[0] - counts[1]) * piece_mass_mod
                    ) % P
                    if total_character == 0:
                        total_zero_cells[object_name].append(
                            (cell, counts)
                        )
                    if signed_character == 0:
                        signed_zero_cells[object_name].append(
                            (cell, counts)
                        )
                    if (
                        total_character != 0
                        and signed_character != 0
                    ):
                        full_c2_unit_cells[object_name] += 1

                common_by_left = {piece[0]: piece for piece in common}
                right_by_left = {piece[0]: piece for piece in right}
                require(
                    len(common_by_left) == len(common)
                    and len(right_by_left) == len(right),
                    f"duplicate left endpoint at {cell}",
                )
                selected_sidecar_cell = clock not in sidecar_clocks_done
                if selected_sidecar_cell:
                    sidecar_clocks_done.add(clock)
                    sidecar_cells.append(cell)
                    faces = indexed_factor_faces(section_factors(
                        full_module, e3, clocks, clock, s, t
                    ))
                    twists = indexed_carrier_twists(
                        source, target_pullback
                    )
                else:
                    faces = twists = None

                common_content = {
                    piece[0]: content(clock, piece) for piece in common
                }
                right_content_cell = {
                    piece[0]: content(clock, piece) for piece in right
                }

                for r_piece in right:
                    r_value = right_content_cell[r_piece[0]]
                    right_content["live" if r_value else "dead"] += 1
                    offsets = []
                    for m_piece in common:
                        delta = (m_piece[0] - r_piece[0]) % T
                        require(
                            delta % H == 0,
                            f"nonintegral half-step relation at {cell}",
                        )
                        offset = delta // H
                        offsets.append(offset)
                        same = (
                            common_content[m_piece[0]] == r_value
                        )
                        parity_even = offset % 2 == 0
                        require(
                            same == parity_even,
                            f"content/parity equivalence failed at {cell}",
                        )
                        content_pair_counts[
                            "equal" if same else "opposite"
                        ] += 1

                    distances = tuple(
                        centered_distance(offset, full_module.C3)
                        for offset in offsets
                    )
                    minimum = min(distances)
                    nearest_indices = tuple(
                        index for index, distance in enumerate(distances)
                        if distance == minimum
                    )
                    require(
                        minimum == 1
                        and len(nearest_indices) == 1
                        and offsets[nearest_indices[0]] == 1,
                        f"centered nearest collar failed at {cell}",
                    )

                    m1 = common_by_left.get((r_piece[0] + H) % T)
                    m2 = common_by_left.get((r_piece[0] + 2 * H) % T)
                    require(
                        m1 is not None
                        and m2 is not None
                        and content(clock, m1) != r_value
                        and content(clock, m2) == r_value,
                        f"one/two-step collar law failed at {cell}",
                    )
                    forward_targets[cell].add(m1[0])
                    plus_two_targets[cell].add(m2[0])
                    plus_one_content[
                        "live" if content(clock, m1) else "dead"
                    ] += 1
                    plus_two_content[
                        "live" if content(clock, m2) else "dead"
                    ] += 1

                    wrapped = shifted_piece(r_piece, WRAP)
                    m14 = common_by_left.get(wrapped[0])
                    r14 = right_by_left.get(wrapped[0])
                    landing = (
                        "M" if m14 is not None
                        else "R" if r14 is not None
                        else "outside"
                    )
                    wrap_landing[landing] += 1
                    if m14 is not None:
                        wrap_targets[cell].add(m14[0])
                        same = content(clock, m14) == r_value
                        wrap_content[
                            "equal" if same else "opposite"
                        ] += 1

                    for step_name, target_piece in (
                        ("h", m1),
                        ("2h", m2),
                        ("14h", m14),
                    ):
                        if target_piece is None:
                            continue
                        a_equal = (
                            ancestry_chamber(r_piece, u_events, v_events)
                            == ancestry_chamber(
                                target_piece, u_events, v_events
                            )
                        )
                        ancestry_chambers.add(
                            ancestry_chamber(r_piece, u_events, v_events)
                        )
                        ancestry_chambers.add(
                            ancestry_chamber(
                                target_piece, u_events, v_events
                            )
                        )
                        ancestry_equal[step_name, a_equal] += 1
                        if selected_sidecar_cell:
                            f_equal = (
                                factor_signature(r_piece, faces)
                                == factor_signature(target_piece, faces)
                            )
                            c_equal = (
                                carrier_signature(r_piece, twists)
                                == carrier_signature(
                                    target_piece, twists
                                )
                            )
                            factor_equal[step_name, f_equal] += 1
                            carrier_equal[step_name, c_equal] += 1
                            if (
                                not f_equal
                                and first_factor_hostile is None
                            ):
                                first_factor_hostile = (
                                    cell,
                                    step_name,
                                    r_piece,
                                    target_piece,
                                    factor_signature(r_piece, faces),
                                    factor_signature(target_piece, faces),
                                )
                            if (
                                not c_equal
                                and first_carrier_hostile is None
                            ):
                                first_carrier_hostile = (
                                    cell,
                                    step_name,
                                    r_piece,
                                    target_piece,
                                    carrier_signature(r_piece, twists),
                                    carrier_signature(
                                        target_piece, twists
                                    ),
                                )
                        if not a_equal and first_ancestry_hostile is None:
                            first_ancestry_hostile = (
                                cell,
                                step_name,
                                r_piece,
                                target_piece,
                                ancestry_chamber(
                                    r_piece, u_events, v_events
                                ),
                                ancestry_chamber(
                                    target_piece, u_events, v_events
                                ),
                            )

                for m_piece in common:
                    reverse_offsets = tuple(
                        ((r_piece[0] - m_piece[0]) % T) // H
                        for r_piece in right
                    )
                    reverse_distances = tuple(
                        centered_distance(offset, full_module.C3)
                        for offset in reverse_offsets
                    )
                    minimum = min(reverse_distances)
                    tie_count = sum(
                        distance == minimum
                        for distance in reverse_distances
                    )
                    reverse_nearest_histogram[
                        minimum, tie_count
                    ] += 1
                    if tie_count == 1:
                        reverse_unique += 1
                    else:
                        reverse_tied += 1

                require(
                    len(forward_targets[cell]) == len(right)
                    and len(plus_two_targets[cell]) == len(right),
                    f"cellwise collar injection collided at {cell}",
                )

    require(
        (nonempty, empty, total_common, total_right, total_pairs)
        == (193, 374, 63_308, 587, 195_517),
        "global bank census changed",
    )
    require(
        dict(cell_types) == EXPECTED_TYPES,
        "cell-type census changed",
    )
    require(
        content_pair_counts
        == Counter({"equal": 97_661, "opposite": 97_856})
        and right_content
        == Counter({"live": 573, "dead": 14})
        and plus_one_content
        == Counter({"dead": 573, "live": 14})
        and plus_two_content
        == Counter({"live": 573, "dead": 14}),
        "delayed-content collar census changed",
    )
    require(
        sum(map(len, forward_targets.values())) == total_right
        and sum(map(len, plus_two_targets.values())) == total_right,
        "global collar injection bookkeeping changed",
    )
    require(
        (reverse_unique, reverse_tied) == (63_066, 242)
        and wrap_landing
        == Counter({"M": 527, "outside": 60})
        and wrap_content == Counter({"equal": 527}),
        "reverse-direction or direct-wrap hostile changed",
    )
    require(
        ancestry_equal
        == Counter({
            ("h", True): 587,
            ("2h", True): 587,
            ("14h", True): 527,
        })
        and factor_equal
        == Counter({
            ("h", False): 6,
            ("2h", False): 6,
            ("14h", False): 6,
        })
        and carrier_equal == factor_equal,
        "ancestry or sampled typed-sidecar boundary changed",
    )
    require(
        parity_unit_cells["M", (True, True)] == 193
        and parity_unit_cells["R", (True, True)] == 14
        and parity_unit_cells["R", (True, False)] == 179
        and full_c2_unit_cells
        == Counter({"M": 88, "R": 193, "M+R": 193})
        and len(signed_zero_cells["M"]) == 105
        and not signed_zero_cells["R"]
        and not signed_zero_cells["M+R"]
        and not total_zero_cells,
        "parity-unit or signed-augmentation census changed",
    )
    tree = ast.parse(Path(__file__).read_text())
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "audit contains a Python assert node",
    )

    print("THM-2825 INDEPENDENT COMMON/RIGHT COLLAR HOSTILE AUDIT")
    print(f"lower_dependency_sha256={tuple(PINNED.items())}")
    print(
        f"scales=T:{T},tau:{SHIFT},h:{H},13tau:{WRAP},"
        f"piece_length:{LENGTH};length=6h/91;"
        "13tau=14h=7(2h)"
    )
    print(
        f"bank_census=nonempty:{nonempty},empty:{empty},"
        f"common:{total_common},right:{total_right},"
        f"complete_bipartite_pairs:{total_pairs}"
    )
    print(f"cell_type_census={tuple(sorted(cell_types.items()))}")
    print(
        f"unscaled_common_choices_per_right={tuple(sorted(unscaled_choice_histogram.items()))};"
        "unscaled_unique=0/587;centered_nearest_unique=587/587;"
        "centered_orientation=right_to_common,+h"
    )
    print(
        f"content_pair_counts={tuple(sorted(content_pair_counts.items()))};"
        "content_equal_iff_half_step_offset_even=yes"
    )
    print(
        f"right_content={tuple(sorted(right_content.items()))};"
        f"plus_h_content={tuple(sorted(plus_one_content.items()))};"
        f"plus_2h_content={tuple(sorted(plus_two_content.items()))}"
    )
    print(
        f"reverse_common_to_right_nearest_histogram="
        f"{tuple(sorted(reverse_nearest_histogram.items()))};"
        f"reverse_unique:{reverse_unique},reverse_tied:{reverse_tied}"
    )
    print(
        f"direct_14h_landing={tuple(sorted(wrap_landing.items()))};"
        f"direct_14h_content={tuple(sorted(wrap_content.items()))}"
    )
    print(
        f"parity_type_histogram="
        f"{tuple(sorted(parity_type_histogram.items()))}"
    )
    print(
        f"parity_colour_unit_cells="
        f"{tuple(sorted(parity_unit_cells.items()))};"
        f"full_C2xC13^5_unit_cells="
        f"{tuple(sorted(full_c2_unit_cells.items()))}"
    )
    print(
        f"total_augmentation_zero_cells="
        f"{tuple((name, tuple(rows)) for name, rows in sorted(total_zero_cells.items()))};"
        f"signed_even_minus_odd_zero_cells="
        f"{tuple((name, tuple(rows)) for name, rows in sorted(signed_zero_cells.items()))}"
    )
    print(
        f"sampled_sidecar_cells={tuple(sidecar_cells)};"
        f"factor_signature_equal={tuple(sorted(factor_equal.items()))};"
        f"first_factor_hostile={first_factor_hostile}"
    )
    print(
        f"carrier_twist_signature_equal={tuple(sorted(carrier_equal.items()))};"
        f"first_carrier_hostile={first_carrier_hostile}"
    )
    print(
        f"ancestry_chamber_equal={tuple(sorted(ancestry_equal.items()))};"
        f"ancestry_chambers_evaluated={len(ancestry_chambers)};"
        f"first_ancestry_hostile={first_ancestry_hostile}"
    )
    print(
        "typing_boundary=the +h and +2h maps land in the common carrier, "
        "so their images are not in the right-domain and the collar map "
        "cannot be iterated merely because 14h=7(2h); sampled factor/"
        "carrier and full ancestry-chamber verdicts are printed above"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
