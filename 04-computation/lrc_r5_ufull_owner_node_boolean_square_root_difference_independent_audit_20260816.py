#!/usr/bin/env python3
"""Independent audit of the Boolean-square x source-root-difference tensor.

The submitted root-difference candidate is not imported.  Starting from the
hash-pinned independent Boolean-square audit, this checker reopens every
selected ordered source-root pair, retains s=u-q before integration, and then
reconstructs the three-coordinate inverse tensor.

The retained difference is a source cut-arc colour.  It is not the THM-2334
grouped exact address, an inverse-ancestry digit, a deep-root label, or a clock.
"""

from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SQUARE_PATH = (
    ROOT
    / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_refiner_"
    "independent_audit_20260816.py"
)
SQUARE_SHA256 = "a75300a81efebef83683c41ac073ffa4d3268da83e96071d7b1b576b36e5bbc7"
SQUARE_SEMANTIC = "af0d543232869e82ee8d0191478ba7a833954cb19b387dedb6fb6f44a6fa272c"

P = 13
V = 4
BANK_NAMES = ("weighted", "support_only", "endpoint_only")
CONTROL_TRIPLES = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_STATE_SEGMENTS = (300115, 300116, 300115, 300116)
EXPECTED_RAW_ARCS = (
    (0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    (0, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2),
    (0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    (0, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2),
)
EXPECTED_FLATTENING_RANKS = ((4, 12, 13), (4, 12, 13), (4, 7, 13))
EXPECTED_ANOVA_RANKS = ((3, 12, 12), (3, 12, 12), (3, 6, 12))
EXPECTED_CENSUS = (676, 1, 12, 12, 144, 3, 36, 36, 432)
EXPECTED_ANOVA_CENSUS = (432, 0, 0, 0, 0, 0, 0, 0, 432)
EXPECTED_FIXED = (
    (4, 3, (10, 12, 10, 12), ((0, 1, 2), (0,), (0, 11, 12), (0,)),
     "acde8efe1065b919716bacf43fbb09afd8dd491b219b67793d8a94f5386eaf5f"),
    (4, 3, (10, 12, 10, 12), ((0, 1, 2), (0,), (0, 11, 12), (0,)),
     "fedbc2430a29f32ff39ec9e7948bde74e896c5d4d2a676fc1e9fdbe395d66c16"),
    (4, 3, (13, 13, 13, 13), ((), (), (), ()),
     "e5a9460cca848458272f5a8917f2786096dd8939ab40ce28e4fb13a5fa9576c7"),
)
EXPECTED_DIGESTS = (
    (
        "a83811099901421a16d39b1a04a9dd09bea96eba1b317208e58df1e6ab670ab6",
        "9076200eefb5d5322e8e20fdbb9fc805d17cc2deb287cc9e184e7e4a06d205d7",
        "b0439188832ee6fa92913e3873a9c3b12eb01c64179843be0a7b38dd7419ce94",
    ),
    (
        "0c6a6b831b5548d012c2954ed5fc11387657591e9de0f6ae3b55a1e593982f9f",
        "147490619b4a53c125140f3e57b0285e5153e6eeff6413410c5cfa9f769c8f7c",
        "38f44bca9e20a27dbb58962d4b750d4650e859c7418cd1bcb350476940fc0a6d",
    ),
    (
        "998d30f5c5415ebc6350b4027fdc8b4392e0e2e034e3e8345c9e8d1944ba9a6d",
        "47b2ae42cd1af50f97c288d769a10707d341ddbb12845c068abfafec6a49c1fe",
        "4f2927945813a30ca0f3cd6c24219189079249b0d836b1c042e137b41e2ddd8c",
    ),
    (
        "7acf66b6d9bd9c7363de7546c40a5f1dfba6e3e19cfe89d9f0d4bd8459b2a13d",
        "f26eb2ee2f86b1da4a1e4762c75531942a7a4f6c39acf3903c5e879f2a14aa86",
        "11014600dea82e614562bef59a203287d50f803fe1f6c766ecc13ab883d7f21d",
    ),
    (
        "6e9562dccd104a355b2fa5e35fb05a9e580f1913330ebd9aab787b32d1a2735e",
        "8618fe32231015f891ea291bd3cdd90ae7408cf39815911aae79b12b6cc86758",
        "4251ea22fea625b9f023a9714ae2a9b48fcbda2b84bbc6f7ecb8b03641ef0b86",
    ),
)
EXPECTED_SQUARE_GAMMA = "8b8f2a2785b084e1578ba0512e4577ab79fd674b84b588fbdb8186f2009242c2"
EXPECTED_ERASED_GAMMA = "cedec6b55700ec2854b426a4f35c73549d720af2d8a34b065a7d16cdaa57f8b5"
EXPECTED_SQUARE_TABLE = "5f4b9609faaa5f148d112a7cde5cfba0ab2c1385b4c53ea9c4bcfc6e93d106fc"
EXPECTED_ERASED_TABLE = "518a002b724438f3a604732576000d61d19b03498dadcb39df32e9f1de4a5a8f"
EXPECTED_SEMANTIC_SHA256 = "25bbbd871fe072915b07662c76661bf72be8d46258a6ce825ab1635ef7fe5c56"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_square():
    require(lf_sha256(SQUARE_PATH) == SQUARE_SHA256, "square audit hash drift")
    spec = importlib.util.spec_from_file_location("root_difference_square", SQUARE_PATH)
    require(spec is not None and spec.loader is not None, "square loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.EXPECTED_SEMANTIC_SHA256 == SQUARE_SEMANTIC,
            "square semantic drift")
    return module


SQ = load_square()
R = SQ.R
PRIME = R.JOINT_PRIME


def integrate_pair(
    alpha: int,
    beta: int,
    word,
    endpoint_grid: int,
    source_u,
    source_v,
    base_boundaries,
    harmonic,
    danger,
    literal_tau: int | None = None,
):
    events, interval_count, mapped_count = R.endpoint_events(
        word, endpoint_grid, alpha, beta, literal_tau
    )
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    banks = [
        [[[0] * P for _state in range(V)] for _tau in tau_values]
        for _bank in BANK_NAMES
    ]
    state_counts = [0] * V
    positions = sorted(set(events) | base_boundaries)
    mask = 0
    active = q_active = weighted_intervals = 0
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
                ("difference escaped owner cell", alpha, beta, left, right))
        state = SQ.state_index(source_u, left)
        state_counts[state] += 1
        if not jump:
            continue
        q_active += 1
        u_values = tuple(R.profile_value(profile, left) for profile in source_u)
        v_values = tuple(R.profile_value(profile, left) for profile in source_v)
        if any(u_values) and any(v_values):
            weighted_intervals += 1
        name = R.chamber(left, right)
        require(name in ("left", "right"), ("difference chamber", name))
        for row_index, tau in enumerate(tau_values):
            selected = mask if literal_tau is not None else (
                mask & R.guard_mask(name, tau, danger)
            )
            roots = tuple(root for root in range(P) if (selected >> root) & 1)
            direct_totals = [0, 0, 0]
            same_totals = [0, 0, 0]
            for u in roots:
                for q in roots:
                    difference = (u - q) % P
                    values = (
                        u_values[u] * v_values[q],
                        int(u_values[u] != 0 and v_values[q] != 0),
                        1,
                    )
                    for bank, value in enumerate(values):
                        banks[bank][row_index][state][difference] = (
                            banks[bank][row_index][state][difference]
                            + value * jump
                        ) % PRIME
                        direct_totals[bank] += value
                        if u == q:
                            same_totals[bank] += value
            expected_totals = (
                sum(u_values[root] for root in roots)
                * sum(v_values[root] for root in roots),
                sum(u_values[root] != 0 for root in roots)
                * sum(v_values[root] != 0 for root in roots),
                len(roots) ** 2,
            )
            require(tuple(direct_totals) == expected_totals,
                    ("pair marginal", alpha, beta, tau, left))
            require(same_totals[0] == 0 and same_totals[1] == 0,
                    ("same-root source cut", alpha, beta, tau, left))
    mask ^= events.get(positions[-1], 0)
    require(mask == 0, ("difference endpoint mask closure", alpha, beta, literal_tau))
    frozen = tuple(
        tuple(tuple(tuple(row) for row in tau_row) for tau_row in bank)
        for bank in banks
    )
    return frozen, (
        interval_count, mapped_count, active, q_active, weighted_intervals
    ), tuple(state_counts)


def build_banks(word, endpoint_grid, source_u, source_v, boundaries, harmonic, danger):
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    gamma = [[] for _bank in BANK_NAMES]
    totals = [0] * 5
    state_totals = [0] * V
    for alpha in range(P):
        for beta in range(P):
            local, counts, state_counts = integrate_pair(
                alpha, beta, word, endpoint_grid, source_u, source_v,
                boundaries, harmonic, danger,
            )
            phase = pow(zeta, beta, PRIME)
            for bank in range(len(BANK_NAMES)):
                for tau in range(P):
                    gamma[bank].append(tuple(
                        tuple(phase * value % PRIME for value in local[bank][tau][state])
                        for state in range(V)
                    ))
            totals = [left + right for left, right in zip(totals, counts)]
            state_totals = [left + right for left, right in zip(state_totals, state_counts)]
    return tuple(tuple(bank) for bank in gamma), tuple(totals), tuple(state_totals)


def inverse_tensor(gamma, zeta: int):
    answer = [[[0] * P for _difference in range(P)] for _state in range(V)]
    normalizer = pow(P**3, -1, PRIME)
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                for relation in range(P):
                    phase = pow(zeta, -(alpha + tau * relation) % P, PRIME)
                    for state in range(V):
                        for difference in range(P):
                            answer[state][difference][relation] = (
                                answer[state][difference][relation]
                                + row[state][difference] * phase
                            ) % PRIME
    require(index == P**3, "difference inverse character count")
    return tuple(tuple(tuple(value * normalizer % PRIME for value in line)
                       for line in sheet) for sheet in answer)


def centre_axes(tensor):
    data = [[[tensor[state][difference][relation] for relation in range(P)]
             for difference in range(P)] for state in range(V)]
    inv_v = pow(V, -1, PRIME)
    inv_p = pow(P, -1, PRIME)
    for difference in range(P):
        for relation in range(P):
            mean = sum(data[state][difference][relation] for state in range(V))
            mean = mean * inv_v % PRIME
            for state in range(V):
                data[state][difference][relation] = (
                    data[state][difference][relation] - mean
                ) % PRIME
    for state in range(V):
        for relation in range(P):
            mean = sum(data[state][difference][relation] for difference in range(P))
            mean = mean * inv_p % PRIME
            for difference in range(P):
                data[state][difference][relation] = (
                    data[state][difference][relation] - mean
                ) % PRIME
    for state in range(V):
        for difference in range(P):
            mean = sum(data[state][difference]) * inv_p % PRIME
            for relation in range(P):
                data[state][difference][relation] = (
                    data[state][difference][relation] - mean
                ) % PRIME
    answer = tuple(tuple(tuple(line) for line in sheet) for sheet in data)
    require(all(sum(answer[state][difference][relation] for state in range(V)) % PRIME == 0
                for difference in range(P) for relation in range(P)), "state centering")
    require(all(sum(answer[state][difference][relation] for difference in range(P)) % PRIME == 0
                for state in range(V) for relation in range(P)), "difference centering")
    require(all(sum(answer[state][difference]) % PRIME == 0
                for state in range(V) for difference in range(P)), "relation centering")
    return answer


def flattening_ranks(tensor):
    state = tuple(tuple(tensor[a][s][t] for s in range(P) for t in range(P))
                  for a in range(V))
    difference = tuple(tuple(tensor[a][s][t] for a in range(V) for t in range(P))
                       for s in range(P))
    relation = tuple(tuple(tensor[a][s][t] for a in range(V) for s in range(P))
                     for t in range(P))
    return R.rank_mod(state), R.rank_mod(difference), R.rank_mod(relation)


def transform(tensor, zeta: int):
    walsh = tuple(tuple(
        -1 if (character[0] * label[0] + character[1] * label[1]) % 2 else 1
        for label in SQ.STATE_LABELS
    ) for character in SQ.STATE_LABELS)
    fourier = tuple(tuple(pow(zeta, -frequency * value % P, PRIME)
                          for value in range(P)) for frequency in range(P))
    return tuple(tuple(tuple(sum(
        tensor[state][difference][relation]
        * walsh[character][state]
        * fourier[difference_frequency][difference]
        * fourier[relation_frequency][relation]
        for state in range(V)
        for difference in range(P)
        for relation in range(P)
    ) % PRIME for relation_frequency in range(P))
        for difference_frequency in range(P))
        for character in range(V))


def census(spectrum):
    counts = [0] * 9
    for character in range(V):
        for difference_frequency in range(P):
            for relation_frequency in range(P):
                if not spectrum[character][difference_frequency][relation_frequency]:
                    continue
                state = character != 0
                difference = difference_frequency != 0
                relation = relation_frequency != 0
                index = {
                    (False, False, False): 1,
                    (False, False, True): 2,
                    (False, True, False): 3,
                    (False, True, True): 4,
                    (True, False, False): 5,
                    (True, False, True): 6,
                    (True, True, False): 7,
                    (True, True, True): 8,
                }[(state, difference, relation)]
                counts[0] += 1
                counts[index] += 1
    return tuple(counts)


def marginal_table(tensor):
    return tuple(tuple(sum(tensor[state][difference][relation]
                           for difference in range(P)) % PRIME
                       for relation in range(P)) for state in range(V))


def fixed_record(tensor):
    matrix = tuple(tuple(tensor[state][difference][6] for difference in range(P))
                   for state in range(V))
    centred = SQ.interaction(matrix)
    return (
        R.rank_mod(matrix),
        R.rank_mod(centred),
        tuple(sum(value != 0 for value in row) for row in matrix),
        tuple(tuple(index for index, value in enumerate(row) if value == 0)
              for row in matrix),
        digest_json(matrix),
    )


def raw_arc_profile():
    rows = []
    for support in SQ.STATE_SUBSETS:
        row = [0] * P
        for u in support:
            for q in range(P):
                if q not in support:
                    row[(u - q) % P] += 1
        rows.append(tuple(row))
    return tuple(rows)


def main() -> None:
    R.split_field_certificate()
    source_u, source_v, source_boundaries, profile_digest, _total, _types = (
        R.source_profiles()
    )
    source_grid = R.SRC.T_DEN
    hierarchy = SQ.source_hierarchy(source_u, source_v, source_boundaries, source_grid)
    word, endpoint_grid = R.endpoint_word_and_grid()
    source_u_scaled = tuple(
        R.scale_profile(profile, R.JOINT_COORDINATE, source_grid) for profile in source_u
    )
    source_v_scaled = tuple(
        R.scale_profile(profile, R.JOINT_COORDINATE, source_grid) for profile in source_v
    )
    boundaries = R.fixed_boundaries(source_boundaries, source_grid)
    harmonic = R.HarmonicPrimitive(word, endpoint_grid)
    danger = R.danger_arcs()
    gamma, work_counts, state_counts = build_banks(
        word, endpoint_grid, source_u_scaled, source_v_scaled,
        boundaries, harmonic, danger,
    )
    require(work_counts == EXPECTED_WORK_COUNTS, ("work counts", work_counts))
    require(state_counts == EXPECTED_STATE_SEGMENTS, ("state counts", state_counts))

    gamma_marginals = tuple(tuple(
        tuple(sum(row[state]) % PRIME for state in range(V)) for row in bank
    ) for bank in gamma)
    require(digest_json(gamma_marginals[0]) == EXPECTED_SQUARE_GAMMA,
            "weighted gamma misses square")
    require(digest_json(gamma_marginals[2]) == EXPECTED_ERASED_GAMMA,
            "endpoint-only gamma misses square erasure")

    guard_records = []
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    for alpha, beta, tau in CONTROL_TRIPLES:
        R.literal_guard_restoration(
            word, endpoint_grid, alpha, beta, tau, boundaries, danger
        )
        direct, _counts, _states = integrate_pair(
            alpha, beta, word, endpoint_grid, source_u_scaled, source_v_scaled,
            boundaries, harmonic, danger, literal_tau=tau,
        )
        phase = pow(zeta, beta, PRIME)
        index = (alpha * P + beta) * P + tau
        for bank in range(len(BANK_NAMES)):
            expected = tuple(tuple(phase * value % PRIME for value in direct[bank][0][state])
                             for state in range(V))
            require(expected == gamma[bank][index],
                    ("literal root-difference guard", bank, alpha, beta, tau))
        guard_records.append((alpha, beta, tau))

    tensors = tuple(inverse_tensor(bank, zeta) for bank in gamma)
    centred = tuple(centre_axes(tensor) for tensor in tensors)
    spectra = tuple(transform(tensor, zeta) for tensor in tensors)
    centred_spectra = tuple(transform(tensor, zeta) for tensor in centred)
    ranks = tuple(flattening_ranks(tensor) for tensor in tensors)
    centred_ranks = tuple(flattening_ranks(tensor) for tensor in centred)
    censuses = tuple(census(spectrum) for spectrum in spectra)
    centred_censuses = tuple(census(spectrum) for spectrum in centred_spectra)
    require(ranks == EXPECTED_FLATTENING_RANKS, ("flattening ranks", ranks))
    require(centred_ranks == EXPECTED_ANOVA_RANKS, ("ANOVA ranks", centred_ranks))
    require(censuses == (EXPECTED_CENSUS,) * 3, ("spectral censuses", censuses))
    require(centred_censuses == (EXPECTED_ANOVA_CENSUS,) * 3,
            ("ANOVA spectral censuses", centred_censuses))

    table_marginals = tuple(marginal_table(tensor) for tensor in tensors)
    require(digest_json(table_marginals[0]) == EXPECTED_SQUARE_TABLE,
            "weighted tensor misses square")
    require(digest_json(table_marginals[2]) == EXPECTED_ERASED_TABLE,
            "endpoint-only tensor misses erased square")

    difference_support = tuple(tuple(
        difference for difference in range(P)
        if any(tensor[state][difference][relation] != 0
               for state in range(V) for relation in range(P))
    ) for tensor in tensors)
    require(difference_support == (
        tuple(range(1, P)), tuple(range(1, P)), tuple(range(P))
    ), ("difference support", difference_support))
    require(all(tensors[bank][state][0][relation] == 0
                for bank in (0, 1) for state in range(V) for relation in range(P)),
            "source-cut same-root survived")
    require(any(tensors[2][state][0][relation] != 0
                for state in range(V) for relation in range(P)),
            "endpoint-only same-root hostile failed")

    arcs = raw_arc_profile()
    require(arcs == EXPECTED_RAW_ARCS, ("raw source arc profile", arcs))
    require(R.rank_mod(arcs) == 2, "raw arc rank")
    fixed = tuple(fixed_record(tensor) for tensor in tensors)
    require(fixed == EXPECTED_FIXED, ("fixed relation records", fixed))

    digests = (
        tuple(digest_json(bank) for bank in gamma),
        tuple(digest_json(tensor) for tensor in tensors),
        tuple(digest_json(spectrum) for spectrum in spectra),
        tuple(digest_json(tensor) for tensor in centred),
        tuple(digest_json(spectrum) for spectrum in centred_spectra),
    )
    require(digests == EXPECTED_DIGESTS, ("candidate comparison digests", digests))
    record = (
        profile_digest, hierarchy, tuple(guard_records), work_counts, state_counts,
        arcs, difference_support, ranks, centred_ranks, censuses,
        centred_censuses, fixed, digests,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("root-difference semantic", semantic))

    print("== independent hostile audit: Boolean square x source root difference x relation ==")
    print(f"parent=(square_sha256={SQUARE_SHA256},square_semantic={SQUARE_SEMANTIC})")
    print("coordinates=(state in V4,source_cut_arc_difference s=u-q in F13,relation=(1,0,t) in F13)")
    print(f"source_profile_sha256={profile_digest};field=(prime={PRIME},zeta13={zeta})")
    print(f"work_counts={work_counts};state_segment_counts={state_counts};literal_controls={tuple(guard_records)}: PASS")
    print("root_pair_gate=ordered selected pairs weighted before integration;difference_orientation=s=u-q;beta_phase=zeta13^beta: PASS")
    print("marginals=(weighted->audited_square,endpoint_only->audited_source_erasure): PASS")
    print("same_root=(weighted=0,support_only=0,endpoint_only_nonzero): PASS")
    print(f"raw_cut_arc_profile=(rank={R.rank_mod(arcs)},rows={arcs})")
    print(f"difference_support_(weighted,support_only,endpoint_only)={difference_support}")
    print(f"flattening_ranks_(state,difference,relation)={ranks}")
    print(f"three_way_ANOVA_flattening_ranks={centred_ranks}")
    print("support_census_order=(total,empty,relation,difference,difference_relation,state,state_relation,state_difference,triple)")
    print(f"spectral_censuses={censuses}")
    print(f"three_way_ANOVA_spectral_censuses={centred_censuses}")
    print(f"fixed_relation=(1,0,6);records_(rank,ANOVA_rank,row_support,zero_sets,digest)={fixed}")
    print(f"digests=(gamma,tensor,spectrum,ANOVA,ANOVA_spectrum)={digests}")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT scoped finite-exact three-coordinate source-cut interaction")
    print("typing=source root difference is a lawful cut-arc colour;not THM-2334 exact address;not inverse ancestry;not chronology")
    print("scope=one owner base;no physical current,no uniform rows,no row exclusion,no LRC(14)")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
