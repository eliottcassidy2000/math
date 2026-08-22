#!/usr/bin/env python3
"""Replace the collapsed seven-cell axis by an intrinsic Boolean square.

The actual U_full owner factor forces the desheeted base into cell zero, so
the seven-cell axis of the parent current is a rank-one delta lift
(MISTAKE-417).  The source service nevertheless has a nontrivial exact
partition inside that cell.  On the distinguished three-root spine
``(0,6,12)``, the active U roots are precisely

    {0}, {0,6}, {6,12}, {12}.

These are the two independent bits

    owner component: left/right,
    source multiplicity: singleton/doubleton.

Every state is an oriented cut S -> S^c on three roots, hence a tournament
with one missing edge.  The six nonempty proper cuts are the combinatorial
completion.  The source actually realizes the five-state Gray path
``{0},{0,6},{6},{6,12},{12}`` and globally omits ``{0,12}``; the owner factor
then excludes the realized centre ``{6}``, leaving a four-state Boolean
square.  This script retains that square before integration and tests its
Walsh x F_13 spectrum, tensor rank, marginals, and literal endpoint guards.

This is a finite-exact one-host candidate.  The Boolean square is a source
support refiner, not an exact THM-2334 address, a U_clock chronology, or a row
exclusion.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = (
    ROOT
    / "04-computation/lrc_r5_ufull_owner_node_common_base_current_probe_20260816.py"
)
PARENT_SHA256 = "6ad93cd13d47159a565f39128b9a108b5203a396ed3ed05a2892222771d45e40"
EXPECTED_STATE_SEGMENTS = (300115, 300116, 300115, 300116)
EXPECTED_STATE_MEASURES = (
    345867129050087785140,
    345867129050087785140,
    345867129050087785140,
    345867129050087785140,
)
EXPECTED_FULL_PATH_SEGMENTS = (3, 12, 3, 12, 3)
EXPECTED_FULL_PATH_TOGGLES = (6, 0, 12, 6)
EXPECTED_FULL_PATH_MEASURES = (
    345867129050087785140,
    4150405548601053421680,
    691734258100175570280,
    4150405548601053421680,
    345867129050087785140,
)
EXPECTED_RANKS = (4, 4, 3)
EXPECTED_SHAPES = (
    (52, 1, 3, 12, 36),
    (52, 1, 3, 12, 36),
    (36, 0, 0, 0, 36),
)
EXPECTED_RELATION_PROFILE = (
    79866518267205440406168,
    310652419985144092895775,
    367085997033220164935953,
    315468006625786970755333,
)
EXPECTED_RELATION_WALSH = (
    317699132065964946247468,
    576205898534886264436774,
    463338744438734120356418,
    472969917720019876075534,
)
EXPECTED_DIGESTS = (
    "8b8f2a2785b084e1578ba0512e4577ab79fd674b84b588fbdb8186f2009242c2",
    "cedec6b55700ec2854b426a4f35c73549d720af2d8a34b065a7d16cdaa57f8b5",
    "5f4b9609faaa5f148d112a7cde5cfba0ab2c1385b4c53ea9c4bcfc6e93d106fc",
    "518a002b724438f3a604732576000d61d19b03498dadcb39df32e9f1de4a5a8f",
    "84f9cbc7352413a13f517f7d2bcb90ef6fd5af5c8cbed9d3bebe8e870da25cb1",
    "1ce6928765e5433f064086a6ca945810e00878b9a8ce436dfc3c5bb369920b2b",
    "bbb69d0cc29652284cb51c292a5e9de71d8c5809dcfd5e73f4c3bf5f991544d6",
)
EXPECTED_SEMANTIC_SHA256 = "bae28345b0b1aea35b244bfbf04123414f0c8fbf9eeca98e39d2b94dd6d107ec"

P = 13
V = 4
SPINE = (0, 6, 12)
STATE_LABELS = ((0, 0), (0, 1), (1, 0), (1, 1))
STATE_SUBSETS = (
    frozenset((0,)),
    frozenset((0, 6)),
    frozenset((12,)),
    frozenset((6, 12)),
)
FULL_PATH_SUBSETS = (
    frozenset((0,)),
    frozenset((0, 6)),
    frozenset((6,)),
    frozenset((6, 12)),
    frozenset((12,)),
)
GLOBALLY_MISSING_CUT = frozenset((0, 12))
OWNER_EXCLUDED_REALIZED_CUT = frozenset((6,))
WALSH_SIGNS = tuple(
    tuple(
        -1 if ((character[0] * state[0] + character[1] * state[1]) & 1)
        else 1
        for state in STATE_LABELS
    )
    for character in STATE_LABELS
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_parent():
    require(lf_sha256(PARENT_PATH) == PARENT_SHA256, "parent source drift")
    spec = importlib.util.spec_from_file_location("owner_boolean_square_parent", PARENT_PATH)
    require(spec is not None and spec.loader is not None, "parent module loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = load_parent()
PRIME = C.JOINT_PRIME


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def source_support(point: int) -> frozenset[int]:
    return frozenset(
        root
        for root, profile in enumerate(C.context()["source_u"])
        if profile_value(profile, point)
    )


def state_of_segment(left: int, right: int) -> int:
    chamber = C.chamber_of_segment(left, right)
    support = source_support(left)
    if chamber == "left" and support == STATE_SUBSETS[0]:
        return 0
    if chamber == "left" and support == STATE_SUBSETS[1]:
        return 1
    if chamber == "right" and support == STATE_SUBSETS[2]:
        return 2
    if chamber == "right" and support == STATE_SUBSETS[3]:
        return 3
    raise RuntimeError(("untyped owner/source state", left, right, chamber, support))


def cut_certificate():
    universe = frozenset(SPINE)
    all_cuts = tuple(
        frozenset(subset)
        for size in (1, 2)
        for subset in itertools.combinations(SPINE, size)
    )
    excluded = (OWNER_EXCLUDED_REALIZED_CUT, GLOBALLY_MISSING_CUT)
    require(set(all_cuts) == set(STATE_SUBSETS) | set(excluded),
            ("six-cut family", all_cuts))
    complement = tuple(
        STATE_SUBSETS.index(universe - subset) for subset in STATE_SUBSETS
    )
    require(complement == (3, 2, 1, 0), ("complement map", complement))
    require(all(complement[state] == state ^ 3 for state in range(V)),
            "complement is not XOR 11")
    records = []
    for subset in all_cuts:
        arcs = tuple(
            (left, right)
            for left in SPINE for right in SPINE
            if left in subset and right not in subset
        )
        missing = one_way = two_way = 0
        for first_index, first in enumerate(SPINE):
            for second in SPINE[first_index + 1:]:
                forward = (first, second) in arcs
                backward = (second, first) in arcs
                if forward and backward:
                    two_way += 1
                elif forward or backward:
                    one_way += 1
                else:
                    missing += 1
        require((missing, one_way, two_way) == (1, 2, 0),
                ("cut tournament census", subset, arcs))
        records.append((tuple(sorted(subset)), arcs, (missing, one_way, two_way)))

    full_segments = {subset: 0 for subset in FULL_PATH_SUBSETS}
    full_measures = {subset: 0 for subset in FULL_PATH_SUBSETS}
    compressed_path = []
    source_boundaries = C.context()["source_boundaries"]
    for left, right in zip(source_boundaries, source_boundaries[1:]):
        support = source_support(left)
        require(support in full_segments,
                ("source support outside five-state path", left, right, support))
        if not compressed_path or compressed_path[-1] != support:
            compressed_path.append(support)
        full_segments[support] += 1
        full_measures[support] += right - left
    require(tuple(compressed_path) == FULL_PATH_SUBSETS,
            ("ordered source support path", compressed_path))
    path_toggles = tuple(
        next(iter(left ^ right))
        for left, right in zip(compressed_path, compressed_path[1:])
        if len(left ^ right) == 1
    )
    require(len(path_toggles) == len(compressed_path) - 1
            and path_toggles == EXPECTED_FULL_PATH_TOGGLES,
            ("source path is not the pinned Gray path", path_toggles))
    full_path_record = tuple(
        (tuple(sorted(subset)), full_segments[subset], full_measures[subset])
        for subset in FULL_PATH_SUBSETS
    )
    require(tuple(full_segments[subset] for subset in FULL_PATH_SUBSETS)
            == EXPECTED_FULL_PATH_SEGMENTS,
            ("full source path segments", full_path_record))
    require(tuple(full_measures[subset] for subset in FULL_PATH_SUBSETS)
            == EXPECTED_FULL_PATH_MEASURES,
            ("full source path measures", full_path_record))
    require(sum(full_measures.values()) == C.JOINT_COORDINATE,
            ("full source path measure", full_path_record))
    require(GLOBALLY_MISSING_CUT not in full_segments,
            "globally missing cut entered source path")

    boundaries = sorted(
        set(C.context()["source_boundaries"])
        | {0, C.JOINT_COORDINATE,
           C.JOINT_COORDINATE // 14,
           13 * C.JOINT_COORDINATE // 14}
    )
    measures = [0] * V
    segments = [0] * V
    for left, right in zip(boundaries, boundaries[1:]):
        if C.cell_of_segment(left, right) != 0:
            continue
        state = state_of_segment(left, right)
        measures[state] += right - left
        segments[state] += 1
    require(sum(measures) == C.JOINT_COORDINATE // 7,
            ("owner-cell measure", measures))
    return (
        tuple(records),
        (tuple(sorted(GLOBALLY_MISSING_CUT)),
         tuple(sorted(OWNER_EXCLUDED_REALIZED_CUT))),
        complement,
        WALSH_SIGNS,
        tuple(segments),
        tuple(measures),
        full_path_record,
        path_toggles,
    )


def integrate_square(alpha: int, beta: int, literal_tau: int | None = None):
    ctx = C.context()
    events, interval_count, mapped = C.endpoint_events(alpha, beta, literal_tau)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    coupled = [[0 for _state in range(V)] for _tau in tau_values]
    source_erasure = [[0 for _state in range(V)] for _tau in tau_values]
    diagonal = [[0 for _state in range(V)] for _tau in tau_values]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    state_segments = [0] * V
    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        require(C.cell_of_segment(left, right) == 0,
                ("owner escaped parent cell", alpha, beta, left, right))
        chamber = C.chamber_of_segment(left, right)
        require(chamber in ("left", "right"),
                ("active middle chamber", alpha, beta, left, right))
        state = state_of_segment(left, right)
        active_segments += 1
        state_segments[state] += 1
        jump = C.q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        u_values = tuple(
            profile_value(profile, left) for profile in ctx["source_u"]
        )
        v_values = tuple(
            profile_value(profile, left) for profile in ctx["source_v"]
        )
        if not any(u_values) or not any(v_values):
            continue
        weighted_segments += 1
        for row_index, tau in enumerate(tau_values):
            if literal_tau is None:
                selected = tuple(
                    sheet
                    for sheet in range(P)
                    if (mask >> sheet) & 1 and C.T.safe(chamber, sheet + tau)
                )
            else:
                selected = tuple(
                    sheet for sheet in range(P) if (mask >> sheet) & 1
                )
            left_value = sum(u_values[sheet] for sheet in selected)
            right_value = sum(v_values[sheet] for sheet in selected)
            diagonal_value = sum(
                u_values[sheet] * v_values[sheet] for sheet in selected
            )
            require(diagonal_value == 0,
                    ("same-root square current", alpha, beta, tau, left, right))
            coupled[row_index][state] = (
                coupled[row_index][state] + left_value * right_value * jump
            ) % PRIME
            source_erasure[row_index][state] = (
                source_erasure[row_index][state] + len(selected) ** 2 * jump
            ) % PRIME
            diagonal[row_index][state] = (
                diagonal[row_index][state] + diagonal_value * jump
            ) % PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("square endpoint mask", alpha, beta, literal_tau, mask))
    require(all(value == 0 for row in diagonal for value in row),
            ("square diagonal", alpha, beta, literal_tau))
    counts = (
        interval_count, mapped, active_segments, q_active_segments,
        weighted_segments, tuple(state_segments),
    )
    return (
        tuple(tuple(row) for row in coupled),
        tuple(tuple(row) for row in source_erasure),
        tuple(tuple(row) for row in diagonal),
        counts,
    )


def worker(alpha: int):
    zeta = C.context()["zeta"]
    coupled_rows = []
    source_erasure_rows = []
    diagonal_rows = []
    scalar_counts = [0] * 5
    state_counts = [0] * V
    for beta in range(P):
        coupled, source_erasure, diagonal, record = integrate_square(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        coupled_rows.append(tuple(
            tuple(phase * value % PRIME for value in row) for row in coupled
        ))
        source_erasure_rows.append(tuple(
            tuple(phase * value % PRIME for value in row)
            for row in source_erasure
        ))
        diagonal_rows.append(diagonal)
        scalar_counts = [
            left + right for left, right in zip(scalar_counts, record[:5])
        ]
        state_counts = [
            left + right for left, right in zip(state_counts, record[5])
        ]
    return (
        alpha,
        tuple(coupled_rows),
        tuple(source_erasure_rows),
        tuple(diagonal_rows),
        (tuple(scalar_counts), tuple(state_counts)),
    )


def inverse_state_table(gamma_states, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    table = [[0 for _relation in range(P)] for _state in range(V)]
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma_states[index]
                index += 1
                for relation in range(P):
                    phase = pow(zeta, -(alpha + tau * relation) % P, PRIME)
                    for state in range(V):
                        table[state][relation] = (
                            table[state][relation] + row[state] * phase
                        ) % PRIME
    require(index == P**3, ("square gamma size", index))
    return tuple(
        tuple(value * normalizer % PRIME for value in row) for row in table
    )


def walsh_fourier(matrix, zeta: int):
    return tuple(
        tuple(
            sum(
                matrix[state][relation]
                * WALSH_SIGNS[character][state]
                * pow(zeta, -frequency * relation % P, PRIME)
                for state in range(V) for relation in range(P)
            ) % PRIME
            for frequency in range(P)
        )
        for character in range(V)
    )


def support_shape(spectrum, rows: int):
    dc = int(spectrum[0][0] != 0)
    state_axis = sum(spectrum[character][0] != 0 for character in range(1, rows))
    residue_axis = sum(spectrum[0][frequency] != 0 for frequency in range(1, P))
    mixed = sum(
        spectrum[character][frequency] != 0
        for character in range(1, rows) for frequency in range(1, P)
    )
    return dc + state_axis + residue_axis + mixed, dc, state_axis, residue_axis, mixed


def interaction(matrix):
    rows = len(matrix)
    inv_rows = pow(rows, -1, PRIME)
    inv_p = pow(P, -1, PRIME)
    inv_total = pow(rows * P, -1, PRIME)
    row_sums = tuple(sum(row) % PRIME for row in matrix)
    column_sums = tuple(
        sum(matrix[row][column] for row in range(rows)) % PRIME
        for column in range(P)
    )
    grand = sum(row_sums) % PRIME
    answer = tuple(
        tuple(
            (
                matrix[row][column]
                - row_sums[row] * inv_p
                - column_sums[column] * inv_rows
                + grand * inv_total
            ) % PRIME
            for column in range(P)
        )
        for row in range(rows)
    )
    require(all(sum(row) % PRIME == 0 for row in answer), "square interaction rows")
    require(all(
        sum(answer[row][column] for row in range(rows)) % PRIME == 0
        for column in range(P)
    ), "square interaction columns")
    return answer


def binary_marginal(table, bit: int):
    return tuple(
        tuple(
            sum(
                table[state][relation]
                for state, label in enumerate(STATE_LABELS) if label[bit] == value
            ) % PRIME
            for relation in range(P)
        )
        for value in range(2)
    )


def binary_fourier(matrix, zeta: int):
    return tuple(
        tuple(
            sum(
                matrix[value][relation]
                * (-1 if character and value else 1)
                * pow(zeta, -frequency * relation % P, PRIME)
                for value in range(2) for relation in range(P)
            ) % PRIME
            for frequency in range(P)
        )
        for character in range(2)
    )


def main() -> None:
    ctx = C.context()
    cut_record = cut_certificate()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "square worker order")
    coupled_states = tuple(
        row
        for _alpha, coupled_rows, _erased, _diag, _counts in chunks
        for beta_rows in coupled_rows for row in beta_rows
    )
    source_erasure_states = tuple(
        row
        for _alpha, _coupled, erased_rows, _diag, _counts in chunks
        for beta_rows in erased_rows for row in beta_rows
    )
    diagonal_states = tuple(
        row
        for _alpha, _coupled, _erased, diagonal_rows, _counts in chunks
        for beta_rows in diagonal_rows for row in beta_rows
    )
    require(len(coupled_states) == source_erasure_states.__len__() == P**3,
            "square bank size")
    require(all(value == 0 for row in diagonal_states for value in row),
            "square same-root bank")

    parent_cells = tuple(
        (sum(row) % PRIME, 0, 0, 0, 0, 0, 0) for row in coupled_states
    )
    parent_erasure_cells = tuple(
        (sum(row) % PRIME, 0, 0, 0, 0, 0, 0)
        for row in source_erasure_states
    )
    require(C.digest_json(parent_cells) == C.EXPECTED_GAMMA_SHA256,
            "square marginal misses parent coupled bank")
    require(C.digest_json(parent_erasure_cells)
            == C.EXPECTED_SOURCE_ERASURE_GAMMA_SHA256,
            "square marginal misses parent source-erasure bank")

    direct_controls = ((0, 0, 0), (0, 1, 6), (1, 0, 6), (6, 6, 12), (12, 12, 0))
    direct_record = []
    for alpha, beta, tau in direct_controls:
        direct, direct_erasure, direct_diagonal, counts = integrate_square(
            alpha, beta, tau
        )
        phase = pow(ctx["zeta"], beta, PRIME)
        direct_row = tuple(phase * value % PRIME for value in direct[0])
        direct_erasure_row = tuple(
            phase * value % PRIME for value in direct_erasure[0]
        )
        index = (alpha * P + beta) * P + tau
        require(direct_row == coupled_states[index],
                ("square literal guard", alpha, beta, tau))
        require(direct_erasure_row == source_erasure_states[index],
                ("square literal erasure guard", alpha, beta, tau))
        require(all(value == 0 for value in direct_diagonal[0]),
                ("square direct diagonal", alpha, beta, tau))
        direct_record.append(
            ((alpha, beta, tau), direct_row, direct_erasure_row, counts)
        )

    table = inverse_state_table(coupled_states, ctx["zeta"])
    source_erasure_table = inverse_state_table(source_erasure_states, ctx["zeta"])
    spectrum = walsh_fourier(table, ctx["zeta"])
    source_erasure_spectrum = walsh_fourier(source_erasure_table, ctx["zeta"])
    centred = interaction(table)
    centred_spectrum = walsh_fourier(centred, ctx["zeta"])
    shapes = (
        support_shape(spectrum, V),
        support_shape(source_erasure_spectrum, V),
        support_shape(centred_spectrum, V),
    )
    coordinate_ranks = (
        C.rank_mod(table), C.rank_mod(source_erasure_table), C.rank_mod(centred)
    )

    parent_profile = tuple(sum(table[state][relation] for state in range(V)) % PRIME
                           for relation in range(P))
    require(digest_json(parent_profile) == C.EXPECTED_RESIDUE_SHA256[0],
            ("square inverse marginal", parent_profile))
    relation = 6
    relation_profile = tuple(table[state][relation] for state in range(V))
    relation_walsh = tuple(
        sum(relation_profile[state] * WALSH_SIGNS[character][state]
            for state in range(V)) % PRIME
        for character in range(V)
    )

    chamber_table = binary_marginal(table, 0)
    multiplicity_table = binary_marginal(table, 1)
    chamber_spectrum = binary_fourier(chamber_table, ctx["zeta"])
    multiplicity_spectrum = binary_fourier(multiplicity_table, ctx["zeta"])
    marginal_records = (
        (C.rank_mod(chamber_table), support_shape(chamber_spectrum, 2),
         digest_json(chamber_table), digest_json(chamber_spectrum)),
        (C.rank_mod(multiplicity_table), support_shape(multiplicity_spectrum, 2),
         digest_json(multiplicity_table), digest_json(multiplicity_spectrum)),
    )

    scalar_counts = tuple(
        sum(chunk[4][0][index] for chunk in chunks) for index in range(5)
    )
    state_counts = tuple(
        sum(chunk[4][1][state] for chunk in chunks) for state in range(V)
    )
    require(scalar_counts == C.EXPECTED_WORK_COUNTS,
            ("square work counts", scalar_counts))
    require(state_counts == EXPECTED_STATE_SEGMENTS,
            ("square state segment counts", state_counts))
    require(cut_record[5] == EXPECTED_STATE_MEASURES,
            ("square state measures", cut_record[5]))
    require(coordinate_ranks == EXPECTED_RANKS,
            ("Boolean square ranks", coordinate_ranks))
    require(shapes == EXPECTED_SHAPES, ("Boolean square shapes", shapes))
    require(relation_profile == EXPECTED_RELATION_PROFILE,
            ("fixed relation state profile", relation_profile))
    require(relation_walsh == EXPECTED_RELATION_WALSH,
            ("fixed relation Walsh profile", relation_walsh))
    require(all(record[0] == 2 and record[1] == (26, 1, 1, 12, 12)
                for record in marginal_records),
            ("Boolean bit marginals", marginal_records))

    digests = (
        digest_json(coupled_states),
        digest_json(source_erasure_states),
        digest_json(table),
        digest_json(source_erasure_table),
        digest_json(spectrum),
        digest_json(source_erasure_spectrum),
        digest_json(centred_spectrum),
    )
    require(digests == EXPECTED_DIGESTS, ("square digests", digests))
    record = (
        PARENT_SHA256,
        cut_record,
        C.JOINT_PRIME,
        C.JOINT_ROOT,
        ctx["zeta"],
        scalar_counts,
        state_counts,
        tuple(direct_record),
        digests,
        shapes,
        coordinate_ranks,
        parent_profile,
        relation,
        relation_profile,
        relation_walsh,
        marginal_records,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("square semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 U_full owner-node Boolean-square refiner probe ==")
    print(f"parent=(sha256={PARENT_SHA256},semantic={C.EXPECTED_SEMANTIC_SHA256})")
    print(f"states=(labels={STATE_LABELS},source_subsets={tuple(tuple(sorted(s)) for s in STATE_SUBSETS)},complement_XOR=3)")
    print(f"six_cut_completion=(global_missing={cut_record[1][0]},owner_excluded_realized={cut_record[1][1]},each_cut_census=(missing=1,one_way=2,two_way=0))")
    print(f"full_source_Gray_path=(subset,segments,measure)={cut_record[6]};toggles={cut_record[7]}")
    print(f"owner_cell_state_partition=(segments={cut_record[4]},measures={cut_record[5]},total={sum(cut_record[5])})")
    print(f"field=(prime={C.JOINT_PRIME},root={C.JOINT_ROOT},zeta13={ctx['zeta']})")
    print(f"work_counts={scalar_counts};state_segment_counts={state_counts}")
    print(f"literal_guard_controls={direct_controls}: PASS;parent_delta_cell_marginal=PASS;same_root=0")
    print(f"coordinate_ranks=(coupled,source_erasure,ANOVA)={coordinate_ranks}")
    print(f"spectral_shapes_(total,dc,V4axis,F13axis,mixed)=(coupled={shapes[0]},source_erasure={shapes[1]},ANOVA={shapes[2]})")
    print(f"fixed_relation=(1,0,{relation});state_profile={relation_profile};Walsh={relation_walsh}")
    print(f"bit_marginals_(chamber,multiplicity)={marginal_records}")
    print(f"digests=(gamma,erased_gamma,table,erased_table,spectrum,erased_spectrum,ANOVA_spectrum)={digests}")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT genuine rank>1 Boolean-square x F13 source/endpoint candidate on one owner-node base")
    print("scope=source support refiner,not exact C(a;X,m),not U_clock chronology,not uniform rows,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
