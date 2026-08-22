#!/usr/bin/env python3
"""Independent audit of the owner-node Boolean-square refiner.

The audited parent is the disjoint THM-2594/THM-3514 common-base checker.
This script imports that parent, not the Boolean-square candidate, and rebuilds
the four-state contraction from the parent endpoint events, source profiles,
and Q(13y) harmonic primitive.

The hierarchy is load-bearing:

* the abstract three-root cut cube has six nonempty proper cuts;
* the physical source profile realizes only the five-state Gray path
  {0}->{0,6}->{6}->{6,12}->{12};
* {0,12} is globally absent, not an OWNER-excluded physical state; and
* OWNER excludes the one realized centre {6}, clipping the two long states
  and leaving four equal-measure states forming a Boolean square.

No exact address, source/arrival identification, U_clock chronology, row
exclusion, or LRC(14) conclusion is asserted.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import importlib.util
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_refiner_"
    "independent_audit_20260816.py"
)
OUTPUT = (
    "05-knowledge/results/"
    "lrc_r5_ufull_owner_node_boolean_square_refiner_"
    "independent_audit_20260816.out"
)
PARENT_PATH = (
    ROOT
    / "04-computation/"
    "lrc_r5_ufull_owner_node_common_base_current_"
    "independent_audit_20260816.py"
)
PARENT_SHA256 = "72e819f1cae92e8969516dc79215e683c538e80da64c9e989e2c5115aebb5304"
PARENT_SEMANTIC_SHA256 = "88d4be52bcb16a52ab2656ff0c0b6bf70e33a2174652ee7e1df62376426f24e6"

P = 13
V = 4
SPINE = frozenset((0, 6, 12))
STATE_LABELS = ((0, 0), (0, 1), (1, 0), (1, 1))
STATE_SUBSETS = (
    frozenset((0,)),
    frozenset((0, 6)),
    frozenset((12,)),
    frozenset((6, 12)),
)
STATE_INDEX = {subset: index for index, subset in enumerate(STATE_SUBSETS)}
PHYSICAL_PATH = (
    frozenset((0,)),
    frozenset((0, 6)),
    frozenset((6,)),
    frozenset((6, 12)),
    frozenset((12,)),
)
TOGGLE_WORD = (6, 0, 12, 6)
PATH_MEASURES = (1, 12, 2, 12, 1)
ABSENT_CUT = frozenset((0, 12))
CONTROL_TRIPLES = ((0, 0, 0), (0, 1, 6), (1, 0, 6), (6, 6, 12), (12, 12, 0))

EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_STATE_SEGMENTS = (300115, 300116, 300115, 300116)
EXPECTED_STATE_MEASURE = 345867129050087785140
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
EXPECTED_PARENT_GAMMA_SHA256 = "b5246eb2a69f35e4dac7dabbf26b1703f21ed22bf803061399ebbf766b9a073d"
EXPECTED_PARENT_ERASED_SHA256 = "20c83a7804c9437a7ccfaae5d3bf685fc1327c4defc385c1cf213f2f9643e258"
EXPECTED_SEMANTIC_SHA256 = "af0d543232869e82ee8d0191478ba7a833954cb19b387dedb6fb6f44a6fa272c"


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
    require(lf_sha256(PARENT_PATH) == PARENT_SHA256, "parent audit hash drift")
    spec = importlib.util.spec_from_file_location("boolean_square_audited_parent", PARENT_PATH)
    require(spec is not None and spec.loader is not None, "parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.EXPECTED_SEMANTIC_SHA256 == PARENT_SEMANTIC_SHA256,
            "parent semantic drift")
    return module


R = load_parent()
PRIME = R.JOINT_PRIME


def source_hierarchy(source_u, source_v, boundaries, grid):
    merged = []
    all_roots = frozenset(range(P))
    for left, right in zip(boundaries, boundaries[1:]):
        u_support = frozenset(
            root for root, profile in enumerate(source_u)
            if R.profile_value(profile, left)
        )
        v_support = frozenset(
            root for root, profile in enumerate(source_v)
            if R.profile_value(profile, left)
        )
        require(v_support == all_roots - u_support,
                ("source complements", left, u_support, v_support))
        require(u_support <= SPINE and u_support, ("source spine", left, u_support))
        if merged and merged[-1][2] == u_support and merged[-1][1] == left:
            merged[-1] = (merged[-1][0], right, u_support)
        else:
            merged.append((left, right, u_support))

    path = tuple(row[2] for row in merged)
    measures = tuple((right - left) * 28 // grid for left, right, _support in merged)
    require(path == PHYSICAL_PATH, ("physical Gray path", path))
    require(measures == PATH_MEASURES, ("Gray path measures", measures))
    toggles = tuple(next(iter(left ^ right)) for left, right in zip(path, path[1:]))
    require(toggles == TOGGLE_WORD, ("Gray toggle word", toggles))

    abstract_cuts = tuple(
        frozenset(subset)
        for size in (1, 2)
        for subset in itertools.combinations(sorted(SPINE), size)
    )
    require(len(abstract_cuts) == 6 and set(abstract_cuts) == set(path) | {ABSENT_CUT},
            ("abstract cut completion", abstract_cuts))
    cut_census = []
    for subset in abstract_cuts:
        missing = one_way = two_way = 0
        for left, right in itertools.combinations(sorted(SPINE), 2):
            forward = left in subset and right not in subset
            backward = right in subset and left not in subset
            if forward and backward:
                two_way += 1
            elif forward or backward:
                one_way += 1
            else:
                missing += 1
        cut_census.append((missing, one_way, two_way))
    require(set(cut_census) == {(1, 2, 0)}, ("cut census", cut_census))

    owner_windows = ((0, grid // 14), (13 * grid // 14, grid))
    owner_measure = {support: 0 for support in abstract_cuts}
    for left, right, support in merged:
        for owner_left, owner_right in owner_windows:
            overlap = max(0, min(right, owner_right) - max(left, owner_left))
            owner_measure[support] += overlap
    require(owner_measure[ABSENT_CUT] == 0, "absent cut gained measure")
    require(owner_measure[frozenset((6,))] == 0, "realized centre survives OWNER")
    for subset in STATE_SUBSETS:
        require(owner_measure[subset] == grid // 28,
                ("owner square measure", subset, owner_measure[subset]))
    realized_excluded = tuple(
        subset for subset in path if owner_measure[subset] == 0
    )
    require(realized_excluded == (frozenset((6,)),),
            ("OWNER-excluded realized states", realized_excluded))

    complements = tuple(STATE_INDEX[SPINE - subset] for subset in STATE_SUBSETS)
    require(complements == tuple(index ^ 3 for index in range(V)),
            ("square complement", complements))
    for index, (label, subset) in enumerate(zip(STATE_LABELS, STATE_SUBSETS)):
        require(label[1] == int(len(subset) == 2), ("multiplicity bit", index))
        require(label[0] == int(12 in subset), ("left/right bit", index))

    path_record = tuple(
        (str(Fraction(left, grid)), str(Fraction(right, grid)), tuple(sorted(support)))
        for left, right, support in merged
    )
    owner_record = tuple(
        (tuple(sorted(subset)), owner_measure[subset]) for subset in abstract_cuts
    )
    return path_record, toggles, owner_record, tuple(cut_census), complements


def state_index(source_u_scaled, point: int) -> int:
    support = frozenset(
        root for root, profile in enumerate(source_u_scaled)
        if R.profile_value(profile, point)
    )
    require(support in STATE_INDEX, ("owner-visible source state", point, support))
    return STATE_INDEX[support]


def integrate_pair(
    alpha: int,
    beta: int,
    word,
    endpoint_grid: int,
    source_u_scaled,
    source_v_scaled,
    base_boundaries,
    harmonic,
    danger,
    literal_tau: int | None = None,
):
    events, interval_count, mapped_count = R.endpoint_events(
        word, endpoint_grid, alpha, beta, literal_tau
    )
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    coupled = [[0] * V for _tau in tau_values]
    erased = [[0] * V for _tau in tau_values]
    diagonal = [[0] * V for _tau in tau_values]
    state_counts = [0] * V
    positions = sorted(set(events) | base_boundaries)
    mask = 0
    active = q_active = weighted = 0
    primitive_left = harmonic.value(positions[0])
    for left, right in zip(positions, positions[1:]):
        mask ^= events.get(left, 0)
        primitive_right = harmonic.value(right)
        jump = (primitive_right - primitive_left) % PRIME
        primitive_left = primitive_right
        if left == right or not mask:
            continue
        active += 1
        require(R.cell_index(left, right) == 0,
                ("square escaped parent cell", alpha, beta, left, right))
        state = state_index(source_u_scaled, left)
        state_counts[state] += 1
        if not jump:
            continue
        q_active += 1
        u_values = tuple(R.profile_value(profile, left) for profile in source_u_scaled)
        v_values = tuple(R.profile_value(profile, left) for profile in source_v_scaled)
        if any(u_values) and any(v_values):
            weighted += 1
        name = R.chamber(left, right)
        require(name in ("left", "right"), ("square chamber", name))
        for row_index, tau in enumerate(tau_values):
            selected = mask if literal_tau is not None else (
                mask & R.guard_mask(name, tau, danger)
            )
            left_sum = right_sum = same_sum = count = 0
            for sheet in range(P):
                if not ((selected >> sheet) & 1):
                    continue
                left_sum += u_values[sheet]
                right_sum += v_values[sheet]
                same_sum += u_values[sheet] * v_values[sheet]
                count += 1
            require(same_sum == 0,
                    ("square same-root", alpha, beta, tau, left, right))
            coupled[row_index][state] = (
                coupled[row_index][state] + left_sum * right_sum * jump
            ) % PRIME
            erased[row_index][state] = (
                erased[row_index][state] + count * count * jump
            ) % PRIME
            diagonal[row_index][state] = (
                diagonal[row_index][state] + same_sum * jump
            ) % PRIME
    mask ^= events.get(positions[-1], 0)
    require(mask == 0, ("square endpoint mask closure", alpha, beta, literal_tau))
    require(all(value == 0 for row in diagonal for value in row),
            "square same-root row")
    return (
        tuple(tuple(row) for row in coupled),
        tuple(tuple(row) for row in erased),
        tuple(tuple(row) for row in diagonal),
        (interval_count, mapped_count, active, q_active, weighted),
        tuple(state_counts),
    )


def build_bank(
    word, endpoint_grid, source_u_scaled, source_v_scaled,
    base_boundaries, harmonic, danger,
):
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    coupled_bank = []
    erased_bank = []
    diagonal_bank = []
    totals = [0] * 5
    state_totals = [0] * V
    for alpha in range(P):
        for beta in range(P):
            coupled, erased, diagonal, counts, state_counts = integrate_pair(
                alpha, beta, word, endpoint_grid, source_u_scaled, source_v_scaled,
                base_boundaries, harmonic, danger,
            )
            phase = pow(zeta, beta, PRIME)
            for tau in range(P):
                coupled_bank.append(tuple(
                    phase * value % PRIME for value in coupled[tau]
                ))
                erased_bank.append(tuple(
                    phase * value % PRIME for value in erased[tau]
                ))
                diagonal_bank.append(tuple(
                    phase * value % PRIME for value in diagonal[tau]
                ))
            totals = [left + right for left, right in zip(totals, counts)]
            state_totals = [left + right for left, right in zip(state_totals, state_counts)]
    return (
        tuple(coupled_bank), tuple(erased_bank), tuple(diagonal_bank),
        tuple(totals), tuple(state_totals),
    )


def inverse_state_table(gamma, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    table = [[0] * P for _state in range(V)]
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                for relation_t in range(P):
                    phase = pow(zeta, -(alpha + tau * relation_t) % P, PRIME)
                    for state in range(V):
                        table[state][relation_t] = (
                            table[state][relation_t] + row[state] * phase
                        ) % PRIME
    require(index == P**3, "square inverse character count")
    return tuple(tuple(value * normalizer % PRIME for value in row) for row in table)


def walsh_fourier(matrix, zeta: int):
    return tuple(tuple(sum(
        matrix[state][relation_t]
        * (-1 if (character[0] * label[0] + character[1] * label[1]) % 2 else 1)
        * pow(zeta, -frequency * relation_t % P, PRIME)
        for state, label in enumerate(STATE_LABELS)
        for relation_t in range(P)
    ) % PRIME for frequency in range(P)) for character in STATE_LABELS)


def support_shape(spectrum, rows: int):
    dc = int(spectrum[0][0] != 0)
    state_axis = sum(spectrum[row][0] != 0 for row in range(1, rows))
    residue_axis = sum(spectrum[0][frequency] != 0 for frequency in range(1, P))
    mixed = sum(
        spectrum[row][frequency] != 0
        for row in range(1, rows) for frequency in range(1, P)
    )
    return dc + state_axis + residue_axis + mixed, dc, state_axis, residue_axis, mixed


def interaction(matrix):
    inv_v = pow(V, -1, PRIME)
    inv_p = pow(P, -1, PRIME)
    inv_vp = pow(V * P, -1, PRIME)
    row_sums = tuple(sum(row) % PRIME for row in matrix)
    column_sums = tuple(sum(matrix[state][relation_t] for state in range(V)) % PRIME
                        for relation_t in range(P))
    grand = sum(row_sums) % PRIME
    answer = tuple(tuple((
        matrix[state][relation_t]
        - row_sums[state] * inv_p
        - column_sums[relation_t] * inv_v
        + grand * inv_vp
    ) % PRIME for relation_t in range(P)) for state in range(V))
    require(all(sum(row) % PRIME == 0 for row in answer), "square ANOVA rows")
    require(all(sum(answer[state][relation_t] for state in range(V)) % PRIME == 0
                for relation_t in range(P)), "square ANOVA columns")
    return answer


def bit_marginal(table, bit: int):
    return tuple(tuple(sum(
        table[state][relation_t]
        for state, label in enumerate(STATE_LABELS) if label[bit] == value
    ) % PRIME for relation_t in range(P)) for value in (0, 1))


def binary_fourier(matrix, zeta: int):
    return tuple(tuple(sum(
        matrix[value][relation_t]
        * (-1 if character * value else 1)
        * pow(zeta, -frequency * relation_t % P, PRIME)
        for value in (0, 1) for relation_t in range(P)
    ) % PRIME for frequency in range(P)) for character in (0, 1))


def main() -> None:
    R.split_field_certificate()
    source_u, source_v, source_boundaries, profile_digest, _total, _types = R.source_profiles()
    source_grid = R.SRC.T_DEN
    hierarchy = source_hierarchy(source_u, source_v, source_boundaries, source_grid)
    word, endpoint_grid = R.endpoint_word_and_grid()
    source_u_scaled = tuple(
        R.scale_profile(profile, R.JOINT_COORDINATE, source_grid) for profile in source_u
    )
    source_v_scaled = tuple(
        R.scale_profile(profile, R.JOINT_COORDINATE, source_grid) for profile in source_v
    )
    base_boundaries = R.fixed_boundaries(source_boundaries, source_grid)
    harmonic = R.HarmonicPrimitive(word, endpoint_grid)
    danger = R.danger_arcs()

    coupled, erased, diagonal, work_counts, state_counts = build_bank(
        word, endpoint_grid, source_u_scaled, source_v_scaled,
        base_boundaries, harmonic, danger,
    )
    require(all(all(value == 0 for value in row) for row in diagonal),
            "square diagonal bank")

    parent_gamma = tuple((sum(row) % PRIME,) + (0,) * 6 for row in coupled)
    parent_erased = tuple((sum(row) % PRIME,) + (0,) * 6 for row in erased)
    require(digest_json(parent_gamma) == EXPECTED_PARENT_GAMMA_SHA256,
            "square misses audited parent gamma")
    require(digest_json(parent_erased) == EXPECTED_PARENT_ERASED_SHA256,
            "square misses audited erased parent")

    guard_records = []
    for alpha, beta, tau in CONTROL_TRIPLES:
        R.literal_guard_restoration(
            word, endpoint_grid, alpha, beta, tau, base_boundaries, danger
        )
        direct, direct_erased, direct_diagonal, _counts, _states = integrate_pair(
            alpha, beta, word, endpoint_grid, source_u_scaled, source_v_scaled,
            base_boundaries, harmonic, danger, literal_tau=tau,
        )
        phase = pow(pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME), beta, PRIME)
        index = (alpha * P + beta) * P + tau
        expected = tuple(phase * value % PRIME for value in direct[0])
        expected_erased = tuple(phase * value % PRIME for value in direct_erased[0])
        require(expected == coupled[index], ("literal square guard", alpha, beta, tau))
        require(expected_erased == erased[index],
                ("literal erased square guard", alpha, beta, tau))
        require(all(value == 0 for value in direct_diagonal[0]),
                ("literal square diagonal", alpha, beta, tau))
        guard_records.append((alpha, beta, tau))

    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    table = inverse_state_table(coupled, zeta)
    erased_table = inverse_state_table(erased, zeta)
    centred = interaction(table)
    spectrum = walsh_fourier(table, zeta)
    erased_spectrum = walsh_fourier(erased_table, zeta)
    centred_spectrum = walsh_fourier(centred, zeta)
    ranks = (R.rank_mod(table), R.rank_mod(erased_table), R.rank_mod(centred))
    shapes = (
        support_shape(spectrum, V),
        support_shape(erased_spectrum, V),
        support_shape(centred_spectrum, V),
    )
    require(ranks == EXPECTED_RANKS, ("square ranks", ranks))
    require(shapes == EXPECTED_SHAPES, ("square shapes", shapes))

    relation_t = 6
    relation_profile = tuple(table[state][relation_t] for state in range(V))
    relation_walsh = tuple(sum(
        relation_profile[state]
        * (-1 if (character[0] * label[0] + character[1] * label[1]) % 2 else 1)
        for state, label in enumerate(STATE_LABELS)
    ) % PRIME for character in STATE_LABELS)
    require(relation_profile == EXPECTED_RELATION_PROFILE,
            ("fixed square profile", relation_profile))
    require(relation_walsh == EXPECTED_RELATION_WALSH,
            ("fixed square Walsh", relation_walsh))
    require(all(value != 0 for value in relation_profile + relation_walsh),
            "fixed square nonvanishing")

    marginal_records = []
    for bit in (0, 1):
        marginal = bit_marginal(table, bit)
        marginal_spectrum = binary_fourier(marginal, zeta)
        record = (
            R.rank_mod(marginal),
            support_shape(marginal_spectrum, 2),
            digest_json(marginal),
            digest_json(marginal_spectrum),
        )
        require(record[0] == 2 and record[1] == (26, 1, 1, 12, 12),
                ("bit marginal", bit, record))
        marginal_records.append(record)

    require(work_counts == EXPECTED_WORK_COUNTS, ("square work counts", work_counts))
    require(state_counts == EXPECTED_STATE_SEGMENTS,
            ("square state segments", state_counts))
    owner_measures_source = tuple(
        next(measure for subset, measure in hierarchy[2]
             if subset == tuple(sorted(state_subset)))
        for state_subset in STATE_SUBSETS
    )
    owner_measures = tuple(
        measure * (R.JOINT_COORDINATE // source_grid)
        for measure in owner_measures_source
    )
    require(owner_measures == (EXPECTED_STATE_MEASURE,) * V,
            ("square owner measures", owner_measures))

    digests = (
        digest_json(coupled), digest_json(erased),
        digest_json(table), digest_json(erased_table),
        digest_json(spectrum), digest_json(erased_spectrum),
        digest_json(centred_spectrum),
    )
    require(digests == EXPECTED_DIGESTS, ("square candidate comparison", digests))

    parent_profile = tuple(sum(table[state][relation] for state in range(V)) % PRIME
                           for relation in range(P))
    require(parent_profile == tuple(R.inverse_table(parent_gamma, zeta)[0]),
            "inverse parent marginal")
    record = (
        profile_digest, hierarchy, tuple(guard_records), work_counts, state_counts,
        owner_measures, ranks, shapes, relation_profile, relation_walsh,
        tuple(marginal_records), parent_profile, digests,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("square audit semantic", semantic))

    print("== independent hostile audit: r=5/U_full Boolean square ==")
    print(f"parent=(script_sha256={PARENT_SHA256},semantic={PARENT_SEMANTIC_SHA256})")
    print(f"physical_gray_path=(states={hierarchy[0]},toggle_word={hierarchy[1]},measure_vector={PATH_MEASURES}/28)")
    print("hierarchy=abstract_nontrivial_cuts=6;physical_states=5;globally_absent={0,12};OWNER_excludes_only_realized_center={6}")
    print(f"abstract_cut_census_(missing,one_way,two_way)={set(hierarchy[3])};no_both_way_edges=PASS")
    print(f"owner_square=(labels={STATE_LABELS},subsets={tuple(tuple(sorted(s)) for s in STATE_SUBSETS)},measures={owner_measures},complement_XOR3={hierarchy[4]})")
    print(f"source_profile_sha256={profile_digest};same_root_pointwise_zero=PASS")
    print(f"literal_guard_controls={tuple(guard_records)}: PASS;parent_marginal=PASS")
    print(f"work_counts={work_counts};state_segment_counts={state_counts}")
    print(f"coordinate_ranks_(coupled,erased,ANOVA)={ranks}")
    print(f"spectral_shapes_(total,dc,V4axis,F13axis,mixed)=(coupled={shapes[0]},erased={shapes[1]},ANOVA={shapes[2]})")
    print(f"bit_marginals_(left_right,single_double)={tuple(marginal_records)}")
    print(f"fixed_(1,0,6)=(state_profile={relation_profile},Walsh={relation_walsh})")
    print(f"comparison_digests={digests}")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT genuine rank4 Boolean-square x F13 mixed candidate; REJECT six-cut physical-support wording")
    print("next_minimal_sidecar=inverse_owner_branch r_owner=a mod 13 on the current U leg; canonical THM-2471 Boolean sheet, cheaper than grouped exact C(a;X,m)")
    print("scope=no exact C(a;X,m),no source/arrival identification,no U_clock chronology,no uniform rows,no row exclusion,no LRC(14)")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
