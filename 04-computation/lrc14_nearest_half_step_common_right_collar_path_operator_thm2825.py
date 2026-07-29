#!/usr/bin/env python3
"""Probe the rooted half-step paths behind the THM-2825 collar.

This is a scratch continuation of the already audited collar census.  It
imports the primary companion rather than rebuilding the physical bank from
lower modules.  The new question is whether the one-step collar extends to
directed paths

    right cofiber --(+h)--> common --(+h)--> common --> ...

and hence whether the semantic-preserving +2h map factors as one
carrier-transverse step followed by one carrier-tangent step.

No Python ``assert`` statement is used.
"""

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
    "lrc14_nearest_half_step_common_right_collar_thm2825.py":
        "bd9ffe7f6815b5c563bd483c300118fbdd683f3d9303babbab7912e031747c9a",
}
for name, expected in PINNED.items():
    payload = (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"pinned collar companion changed: {name}",
    )

PINNED_RESULTS = {
    "lrc14_nearest_half_step_common_right_collar_thm2825.out":
        "c4a31e5ee0aa5af69faa3efbe315d0968a85cba49c2d77c0ca93a229bc39fa0c",
    "lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825.out":
        "e610c5bffb720e074662f2222a50b2ce461c1b1293e946aa260faf898c4347b7",
    "lrc14_nearest_half_step_common_right_collar_independent_audit_thm2825.out":
        "b9846e35c615b6ff0183f6df47ba8e8f99796a195afb6480f922254e35cbe82d",
}
RESULTS = ROOT / "05-knowledge" / "results"
for name, expected in PINNED_RESULTS.items():
    payload = (RESULTS / name).read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"pinned collar result changed: {name}",
    )


import lrc14_nearest_half_step_common_right_collar_thm2825 as collar


P = 13
H = collar.H
T = collar.copies.T
LENGTH = collar.copies.LENGTH


def shifted_left(left, steps):
    result = left + steps * H
    require(
        0 <= result and result + LENGTH <= T,
        "tested path met the ambient circle seam",
    )
    return result


def classify(left, common_by_left, right_by_left):
    if left in common_by_left:
        return "common"
    if left in right_by_left:
        return "right"
    return "outside"


def cyclic_convolution(left, right):
    answer = [0] * P
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[(i + j) % P] += a * b
    return tuple(value % P for value in answer)


ZERO3 = ((0, 0, 0), (0, 0, 0), (0, 0, 0))


def matrix_unit(row, column):
    return tuple(
        tuple(
            1 if (i, j) == (row, column) else 0
            for j in range(3)
        )
        for i in range(3)
    )


def matrix_add(*matrices):
    return tuple(
        tuple(sum(matrix[i][j] for matrix in matrices) for j in range(3))
        for i in range(3)
    )


def matrix_scale(scalar, matrix):
    return tuple(
        tuple(scalar * matrix[i][j] for j in range(3))
        for i in range(3)
    )


def matrix_multiply(left, right):
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(3))
            for j in range(3)
        )
        for i in range(3)
    )


def diagonal_matrix(entries):
    return tuple(
        tuple(entries[i] if i == j else 0 for j in range(3))
        for i in range(3)
    )


def literal_output_rows(name):
    rows = {}
    for line in (RESULTS / name).read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        value = value.split(";", 1)[0]
        try:
            rows[key] = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            pass
    return rows


def main():
    (
        _module,
        _rails,
        _present,
        details,
        full_module,
        e3,
        clocks,
        q_pairs,
        _delayed,
        _source_weight,
        _target_weight,
        _rail_common,
    ) = collar.copies.physical_setup()

    cells = 0
    roots = 0
    common_atoms = 0
    covered_atoms = 0
    half_step_edges = 0
    transverse_edges = 0
    rooted_tangent_edges = 0
    path_length_histogram = Counter()
    terminal_histogram = Counter()
    predecessor_histogram = Counter()
    cell_profile_histogram = Counter()
    two_step_images = set()
    direct_fourteen_histogram = Counter()
    local_fourteen_histogram = Counter()
    direct_local_disagreement = Counter()
    first_uncovered = None
    first_overlap = None
    first_two_step_failure = None
    all_path_records = []
    common_component_histogram = Counter()
    common_component_predecessors = Counter()
    common_component_terminals = Counter()
    carrier_only_component_histogram = Counter()
    carrier_only_components = 0
    carrier_only_atoms = 0
    cofiber_rooted_components = 0
    max_path_length = 0
    carrier_only_start_status = Counter()
    carrier_only_unit_types = Counter()
    carrier_only_inverse_checks = 0
    semantic_caches = tuple({} for _clock in range(7))
    physical_root_counts = Counter()
    physical_m2_counts = Counter()
    m1_labels = set()
    m2_labels = set()
    semantic_triples = Counter()

    for clock in range(7):
        for sigma in collar.COMMON_S:
            for target in collar.COMMON_T:
                _source, _target, common, right = collar.cell_objects(
                    details,
                    full_module,
                    e3,
                    clocks,
                    clock,
                    sigma,
                    target,
                )
                if not right:
                    continue
                cells += 1
                roots += len(right)
                common_atoms += len(common)
                common_by_left = {piece[0]: piece for piece in common}
                right_by_left = {piece[0]: piece for piece in right}
                require(
                    len(common_by_left) == len(common)
                    and len(right_by_left) == len(right),
                    "left endpoints collided within a labelled cell",
                )

                seen = {}
                cell_lengths = []
                for root_index, rpiece in enumerate(right):
                    rleft = rpiece[0]
                    predecessor_histogram[
                        classify(
                            shifted_left(rleft, -1),
                            common_by_left,
                            right_by_left,
                        )
                    ] += 1
                    path = []
                    step = 1
                    while True:
                        left = shifted_left(rleft, step)
                        piece = common_by_left.get(left)
                        if piece is None:
                            break
                        if left in seen and first_overlap is None:
                            first_overlap = (
                                (clock, sigma, target),
                                rleft,
                                step,
                                seen[left],
                            )
                        seen[left] = (rleft, step)
                        path.append(piece)
                        step += 1

                    terminal_histogram[
                        classify(
                            shifted_left(rleft, step),
                            common_by_left,
                            right_by_left,
                        )
                    ] += 1
                    length = len(path)
                    require(length >= 2, "root lost one of its two collars")
                    max_path_length = max(max_path_length, length)
                    path_length_histogram[length] += 1
                    cell_lengths.append(length)
                    covered_atoms += length
                    half_step_edges += length
                    transverse_edges += 1
                    rooted_tangent_edges += length - 1
                    physical_root_counts[rpiece] += 1
                    cell = (clock, sigma, target)
                    m1_label = (cell, path[0][0])
                    m2_label = (cell, path[1][0])
                    m1_labels.add(m1_label)
                    m2_labels.add(m2_label)
                    semantic_values = tuple(
                        collar.semantic_value(
                            piece,
                            q_pairs[clock],
                            semantic_caches[clock],
                        )
                        != (0, 0)
                        for piece in (rpiece, path[0], path[1])
                    )
                    require(
                        semantic_values[0] == semantic_values[2]
                        and semantic_values[0] != semantic_values[1],
                        "collar triple lost its relative semantic grading",
                    )
                    semantic_triples[semantic_values] += 1

                    two_step = shifted_left(rleft, 2)
                    if two_step not in common_by_left:
                        first_two_step_failure = (
                            (clock, sigma, target),
                            rpiece,
                        )
                    else:
                        two_step_images.add(
                            ((clock, sigma, target), two_step)
                        )
                        physical_m2_counts[common_by_left[two_step]] += 1

                    direct_fourteen = classify(
                        shifted_left(rleft, 14),
                        common_by_left,
                        right_by_left,
                    )
                    local_fourteen = (
                        "common" if length >= 14 else "stopped"
                    )
                    direct_fourteen_histogram[direct_fourteen] += 1
                    local_fourteen_histogram[local_fourteen] += 1
                    direct_local_disagreement[
                        (direct_fourteen, local_fourteen)
                    ] += 1
                    all_path_records.append(
                        (
                            (clock, sigma, target),
                            root_index,
                            rleft,
                            length,
                            path[-1][0],
                            classify(
                                shifted_left(rleft, length + 1),
                                common_by_left,
                                right_by_left,
                            ),
                        )
                    )

                uncovered = sorted(set(common_by_left) - set(seen))
                if uncovered and first_uncovered is None:
                    first_uncovered = (
                        (clock, sigma, target),
                        tuple(uncovered[:10]),
                        len(uncovered),
                    )
                component_seen = set()
                rooted_starts = {
                    shifted_left(rpiece[0], 1) for rpiece in right
                }
                starts = sorted(
                    left for left in common_by_left
                    if shifted_left(left, -1) not in common_by_left
                )
                for start in starts:
                    component = []
                    left = start
                    while left in common_by_left:
                        require(
                            left not in component_seen,
                            "common component overlap",
                        )
                        component_seen.add(left)
                        component.append(left)
                        left = shifted_left(left, 1)
                    predecessor = classify(
                        shifted_left(start, -1),
                        common_by_left,
                        right_by_left,
                    )
                    terminal = classify(
                        left,
                        common_by_left,
                        right_by_left,
                    )
                    common_component_histogram[len(component)] += 1
                    common_component_predecessors[predecessor] += 1
                    common_component_terminals[terminal] += 1
                    if start in rooted_starts:
                        require(
                            predecessor == "right",
                            "cofiber-rooted component lost its root",
                        )
                        cofiber_rooted_components += 1
                    else:
                        require(
                            predecessor == "outside",
                            "carrier-only component has an unexpected predecessor",
                        )
                        carrier_only_components += 1
                        carrier_only_atoms += len(component)
                        carrier_only_component_histogram[len(component)] += 1
                        require(
                            len(component) % P == 0
                            and (len(component) // P) % 2 == 1,
                            "carrier-only component is not an odd C13 tower",
                        )
                        statuses = []
                        for component_left in component:
                            value = collar.semantic_value(
                                common_by_left[component_left],
                                q_pairs[clock],
                                semantic_caches[clock],
                            )
                            require(
                                value
                                in (
                                    (0, 0),
                                    (collar.copies.C, collar.copies.C),
                                ),
                                "carrier-only semantic value changed",
                            )
                            statuses.append(value != (0, 0))
                        require(
                            all(
                                statuses[index] != statuses[index - 1]
                                for index in range(1, len(statuses))
                            ),
                            "carrier-only component lost semantic alternation",
                        )
                        epsilon = 1 if statuses[0] else -1
                        blocks = len(component) // P
                        raw = tuple(blocks % P for _residue in range(P))
                        live = tuple(
                            sum(
                                statuses[index]
                                for index in range(residue, len(statuses), P)
                            ) % P
                            for residue in range(P)
                        )
                        alternating = tuple(
                            1 if residue % 2 == 0 else -1
                            for residue in range(P)
                        )
                        require(
                            tuple(
                                (2 * live[index] - raw[index]) % P
                                for index in range(P)
                            )
                            == tuple(
                                epsilon * value % P
                                for value in alternating
                            ),
                            "live/raw semantic contrast formula changed",
                        )
                        norm = (1,) * P
                        live_inverse = tuple(
                            (
                                epsilon * (index in (0, 1))
                                - 2 * blocks * norm[index]
                            ) % P
                            for index in range(P)
                        )
                        contrast_inverse = tuple(
                            (
                                epsilon
                                * pow(2, -1, P)
                                * (index in (0, 1))
                            ) % P
                            for index in range(P)
                        )
                        delta = (1,) + (0,) * (P - 1)
                        contrast = tuple(
                            (2 * live[index] - raw[index]) % P
                            for index in range(P)
                        )
                        require(
                            cyclic_convolution(live, live_inverse) == delta
                            and cyclic_convolution(
                                contrast, contrast_inverse
                            )
                            == delta,
                            "closed carrier-only inverse formula failed",
                        )
                        carrier_only_inverse_checks += 1
                        carrier_only_start_status[
                            "live" if statuses[0] else "dead"
                        ] += 1
                        carrier_only_unit_types[
                            (blocks, epsilon, tuple(live), live_inverse)
                        ] += 1
                require(
                    component_seen == set(common_by_left),
                    "common path decomposition missed an atom",
                )
                covered_atoms += 0
                cell_profile_histogram[
                    (
                        clock,
                        len(common),
                        len(right),
                        tuple(sorted(cell_lengths)),
                        len(uncovered),
                    )
                ] += 1

    require(cells == 193 and roots == 587, "collar input census changed")
    require(
        first_two_step_failure is None and len(two_step_images) == roots,
        "the +2h collar stopped being an injection",
    )
    require(
        transverse_edges == roots
        and rooted_tangent_edges == half_step_edges - roots,
        "carrier transverse/tangent edge split failed",
    )
    require(
        len(m1_labels) == len(m2_labels) == roots
        and not (m1_labels & m2_labels)
        and semantic_triples
        == Counter({
            (True, False, True): 573,
            (False, True, False): 14,
        }),
        "the three collar summands stopped being orthogonal and graded",
    )
    root_semantic_counts = Counter({
        "live": semantic_triples[(True, False, True)],
        "dead": semantic_triples[(False, True, False)],
    })
    require(
        root_semantic_counts == Counter({"live": 573, "dead": 14}),
        "right-root semantic imbalance changed",
    )
    require(
        cofiber_rooted_components == roots
        and carrier_only_atoms == common_atoms - covered_atoms,
        "rooted/carrier-only common decomposition changed",
    )
    rank_profile = tuple(
        (
            step,
            sum(
                count
                for length, count in path_length_histogram.items()
                if length >= step
            ),
        )
        for step in range(1, max_path_length + 1)
    )
    require(
        rank_profile[0] == (1, roots)
        and rank_profile[1] == (2, roots)
        and rank_profile[13] == (14, 527),
        "transverse power-rank controls changed",
    )
    rank_plateaus = []
    plateau_start, plateau_rank = rank_profile[0]
    previous_step = plateau_start
    for step, rank in rank_profile[1:]:
        if rank != plateau_rank:
            rank_plateaus.append(
                (plateau_start, previous_step, plateau_rank)
            )
            plateau_start = step
            plateau_rank = rank
        previous_step = step
    rank_plateaus.append(
        (plateau_start, previous_step, plateau_rank)
    )
    full_tangent_edges = common_atoms - sum(
        common_component_histogram.values()
    )
    full_half_step_edges = transverse_edges + full_tangent_edges
    require(
        full_tangent_edges == 62_623
        and full_half_step_edges == 63_210,
        "full half-step incidence ranks changed",
    )
    require(
        {
            (
                piece[0] + 2 * H,
                piece[1] + 2 * H,
                piece[2],
            ): count
            for piece, count in physical_root_counts.items()
        }
        == dict(physical_m2_counts),
        "label-forgetting multiplicities do not descend under +2h",
    )
    root_multiplicity_histogram = Counter(
        physical_root_counts.values()
    )
    matrix_units = {
        (i, j): matrix_unit(i, j)
        for i in range(3)
        for j in range(3)
    }
    matrix_unit_checks = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for ell in range(3):
                    expected = (
                        matrix_units[i, ell] if j == k else ZERO3
                    )
                    require(
                        matrix_multiply(
                            matrix_units[i, j],
                            matrix_units[k, ell],
                        )
                        == expected,
                        "M3 collar matrix-unit law failed",
                    )
                    matrix_unit_checks += 1
    ladder = matrix_add(matrix_units[1, 0], matrix_units[2, 1])
    even_collar = matrix_units[2, 0]
    ladder_square = matrix_multiply(ladder, ladder)
    ladder_cube = matrix_multiply(ladder_square, ladder)
    semantic_grading = diagonal_matrix((1, -1, 1))
    carrier_grading = diagonal_matrix((-1, 1, 1))
    transverse_step = matrix_units[1, 0]
    tangent_step = matrix_units[2, 1]
    require(
        matrix_unit_checks == 81
        and ladder_square == even_collar
        and ladder_cube == ZERO3,
        "three-state nilpotent collar ladder changed",
    )
    require(
        matrix_add(
            matrix_multiply(semantic_grading, ladder),
            matrix_multiply(ladder, semantic_grading),
        )
        == ZERO3
        and matrix_multiply(semantic_grading, even_collar)
        == matrix_multiply(even_collar, semantic_grading),
        "semantic grading stopped making B odd and B^2 even",
    )
    require(
        matrix_add(
            matrix_multiply(carrier_grading, transverse_step),
            matrix_multiply(transverse_step, carrier_grading),
        )
        == ZERO3
        and matrix_multiply(carrier_grading, tangent_step)
        == matrix_multiply(tangent_step, carrier_grading)
        and matrix_add(
            matrix_multiply(carrier_grading, even_collar),
            matrix_multiply(even_collar, carrier_grading),
        )
        == ZERO3,
        "carrier grading stopped separating transverse and tangent steps",
    )

    primary_rows = literal_output_rows(
        "lrc14_nearest_half_step_common_right_collar_thm2825.out"
    )
    physical_rows = literal_output_rows(
        "lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825.out"
    )
    independent_rows = literal_output_rows(
        "lrc14_nearest_half_step_common_right_collar_independent_audit_thm2825.out"
    )
    require(
        primary_rows["relation_parity"]
        == (((0, True), 97_661), ((1, False), 97_856)),
        "semantic half-step parity input changed",
    )
    factor_r = physical_rows["factor_census_R"]
    factor_m1 = physical_rows["factor_census_M1"]
    factor_m2 = physical_rows["factor_census_M2"]
    carrier_r = physical_rows["carrier_census_R"]
    carrier_m1 = physical_rows["carrier_census_M1"]
    carrier_m2 = physical_rows["carrier_census_M2"]
    endpoint_2h = physical_rows["endpoint_pair_census_2h"]
    factor_holes = Counter()
    for signatures, count in factor_r:
        missing = tuple(
            collar.copies.FACTOR_NAMES[index]
            for index, value in enumerate(signatures[0])
            if not value
        )
        factor_holes[missing] += count
    require(
        sum(factor_holes.values()) == roots
        and factor_m1 == factor_m2
        and factor_m2
        == (
            (
                (
                    (True,) * 6,
                    (True,) * 6,
                    (True,) * 6,
                    (True,) * 6,
                ),
                roots,
            ),
        ),
        "factor sidecar consequence changed",
    )
    require(
        carrier_r[0][1] == roots
        and carrier_m1 == carrier_m2
        and carrier_m2[0][1] == roots
        and not any(carrier_r[0][0][0])
        and carrier_r[0][0][1][0]
        and carrier_m2[0][0][0][0]
        and carrier_m2[0][0][1][0],
        "carrier transverse sidecar changed",
    )
    source_endpoint_exact = sum(
        count for row, count in endpoint_2h if row[0] == 0
    )
    target_endpoint_exact = sum(
        count for row, count in endpoint_2h if row[1] == 0
    )
    require(
        source_endpoint_exact == 74
        and target_endpoint_exact == roots
        and all(row[3] == ((0, 0),) for row, _count in endpoint_2h),
        "endpoint sidecar consequence changed",
    )
    ancestry_rows = dict(
        independent_rows["ancestry_chamber_equal"]
    )
    require(
        ancestry_rows[("2h", True)] == roots,
        "two-half-step ancestry chamber changed",
    )

    print("THM-2825 ROOTED HALF-STEP PATH PROBE")
    print(f"pinned={tuple(PINNED.items())}")
    print(
        f"cells={cells};roots={roots};common_atoms={common_atoms};"
        f"covered_atoms={covered_atoms};half_step_edges={half_step_edges}"
    )
    print(
        f"first_overlap={first_overlap};"
        f"first_uncovered={first_uncovered}"
    )
    print(
        "path_length_histogram="
        f"{tuple(sorted(path_length_histogram.items()))}"
    )
    print(
        f"terminal_histogram={tuple(sorted(terminal_histogram.items()))};"
        "predecessor_histogram="
        f"{tuple(sorted(predecessor_histogram.items()))}"
    )
    print(
        f"common_components={sum(common_component_histogram.values())};"
        f"cofiber_rooted={cofiber_rooted_components};"
        f"carrier_only={carrier_only_components};"
        f"carrier_only_atoms={carrier_only_atoms}"
    )
    print(
        "common_component_histogram="
        f"{tuple(sorted(common_component_histogram.items()))}"
    )
    print(
        "carrier_only_component_histogram="
        f"{tuple(sorted(carrier_only_component_histogram.items()))};"
        "component_predecessors="
        f"{tuple(sorted(common_component_predecessors.items()))};"
        "component_terminals="
        f"{tuple(sorted(common_component_terminals.items()))}"
    )
    print(
        "carrier_only_start_status="
        f"{tuple(sorted(carrier_only_start_status.items()))};"
        f"unit_inverse_checks={carrier_only_inverse_checks}/"
        f"{carrier_only_components}"
    )
    print(
        "carrier_only_unit_types="
        f"{tuple(sorted(carrier_only_unit_types.items()))}"
    )
    print(
        "carrier_only_group_algebra="
        "for length=13m (m odd), raw=m*N13, "
        "2live-raw=epsilon*A13, "
        "A13=sum(-1)^r X^r=2/(1+X); "
        "live^(-1)=(1+X)(epsilon-m*N13), "
        "contrast^(-1)=epsilon(1+X)/2"
    )
    print(
        f"transverse_edges={transverse_edges};"
        f"rooted_tangent_edges={rooted_tangent_edges};"
        f"full_tangent_edges={full_tangent_edges};"
        f"full_half_step_edges={full_half_step_edges};"
        f"two_step_images={len(two_step_images)}"
    )
    print(
        f"label_forgetting=distinct_physical_roots:"
        f"{len(physical_root_counts)},"
        f"distinct_physical_M2:{len(physical_m2_counts)},"
        "multiplicity_histogram:"
        f"{tuple(sorted(root_multiplicity_histogram.items()))};"
        "translation_descends_on_interval_multiplicities=yes"
    )
    print(
        f"collar_triples={roots};"
        f"M1_images={len(m1_labels)};"
        f"M2_images={len(m2_labels)};"
        "M1_intersect_M2=0;"
        f"semantic_triples={tuple(sorted(semantic_triples.items()))}"
    )
    print(
        "M3_linking=W_i^*W_j=delta_ij*I_587;"
        f"matrix_unit_checks={matrix_unit_checks}/81;"
        "E_ijE_kl=delta_jkE_il;"
        "span=M3_tensor_I587;dimension=1761"
    )
    print(
        "compressed_ladder=B=E10+E21;"
        f"B2={ladder_square};B3={ladder_cube};"
        "B2=E20=plus_2h_collar;"
        "semantic_grading=diag(1,-1,1):B_odd,B2_even;"
        "carrier_grading=diag(-1,1,1):"
        "E10_transverse,E21_tangent,E20_transverse"
    )
    print(
        "punctured_V4_grades="
        "R:(semantic,carrier)=(0,0),"
        "M1:(1,1),M2:(0,1),missing:(1,0);"
        "edge_degrees=d:(1,1),a:(1,0),ad:(0,1);"
        f"right_semantic_counts={tuple(sorted(root_semantic_counts.items()))};"
        "no semantic-flipping permutation of R"
    )
    print(
        "direct_14h="
        f"{tuple(sorted(direct_fourteen_histogram.items()))};"
        "local_N14="
        f"{tuple(sorted(local_fourteen_histogram.items()))};"
        "joint="
        f"{tuple(sorted(direct_local_disagreement.items()))}"
    )
    print(
        f"factor_holes_R={tuple(sorted(factor_holes.items()))};"
        "factor_M2=all_six_on_587"
    )
    print(
        "carrier_R=source_empty,target_delta0;"
        "carrier_M2=source_delta0,target_delta0;"
        f"endpoint_exact=source:{source_endpoint_exact},"
        f"target:{target_endpoint_exact};"
        "target_translation_delta0=587;"
        "ancestry_2h_equal=587"
    )
    print(
        "cell_profile_histogram="
        f"{tuple(sorted(cell_profile_histogram.items()))}"
    )
    print(
        "operator_law=on the free labelled interval module N=d+a with "
        "d=P*N*(1-P), a=P*N*P; "
        "d^2=d*a=0; a*d=P*N^2*(1-P) is the +2h collar; "
        "{semantic,d}={semantic,a}=0, hence [semantic,a*d]=0; "
        "[P,a*d]=a*d"
    )
    print(
        f"nilpotence=N^{max_path_length + 1}=0 sharp on rooted forest;"
        f"a^{max_path_length}=0 sharp;"
        "rank_plateaus(a^(k-1)d)="
        f"{tuple(rank_plateaus)}"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
