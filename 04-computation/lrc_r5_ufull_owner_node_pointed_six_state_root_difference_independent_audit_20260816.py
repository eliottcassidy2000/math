#!/usr/bin/env python3
"""Independent audit of the pointed-six source-root refinement.

The submitted pointed-state candidate is not imported.  This checker starts
from the hash-pinned independent root-difference audit, reopens each ordered
source-root pair, and retains the actual U-tail root before integration.
"""

from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = (
    ROOT
    / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_root_difference_"
    "independent_audit_20260816.py"
)
PARENT_SHA256 = "bdb211b727550ba51c1ca83490fbda50446b6e3c0596e229cd864e741af93509"
PARENT_SEMANTIC = "25bbbd871fe072915b07662c76661bf72be8d46258a6ce825ab1635ef7fe5c56"

P = 13
POINTS = ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12))
POINT_INDEX = {point: index for index, point in enumerate(POINTS)}
STATE_FIBRES = ((0,), (1, 2), (5,), (3, 4))
CONTROL_TRIPLES = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_POINT_SEGMENTS = (41935, 51187, 51187, 51187, 51187, 41935)
EXPECTED_RANKS = ((6, 12, 13), (6, 12, 13), (4, 12, 13))
EXPECTED_ANOVA_RANKS = ((5, 12, 12), (5, 12, 12), (3, 12, 12))
EXPECTED_EDGE_RECORDS = (
    (
        (120, 120, (1, 12), "adc8dc9a773b7da8cf98c58a9324fd53e0012385449779411198a937452c49cd"),
        (144, 144, (), "a44a82621a9933055ed49f95c2f3db21e9d001ca9dad6fcd89213fe8f8290d06"),
        (120, 120, (6, 7), "8794ecd5253127ee185fdf0e80a1920463eac164c5a83bade7e4de290d4e5a86"),
        (144, 144, (), "8e8f1323e5ce4e46e40761d1139e1d988ebfd63c153a3504bd2dc11dcf9f0e75"),
        (120, 120, (1, 12), "0928acde9d91e2423440ef8facdccfe345780e4311f2d064e481c5f24e95f5ac"),
    ),
    (
        (120, 120, (1, 12), "7aa0b1f6da9ffda319accf689f9af977bf2b0acf8d588009d28bca16bfc521f2"),
        (144, 144, (), "fbe2dcafee05570f9618266d2879efd960bbd22b634bae6752826059f6c5a053"),
        (120, 120, (6, 7), "89dd0924800a2d44fa0c469d28d1e49d001d03da50fd30915836253e0ec79aba"),
        (144, 144, (), "99085cb82122a7bf3f172a1eeccb3bc482054a0e7c3c0d52d85920775a6a4afd"),
        (120, 120, (1, 12), "818f3ed81e009575090cc06a15a39e63adfc6bd8edb688cc508448c2f71bc427"),
    ),
    (
        (144, 144, (), "444390c496e60703467a17690199380849f7b25e500a3e4bcba477a754915fe0"),
        (0, 0, tuple(range(1, P)), "633ab4de4d7c45b84d18d15e75c27d5b4f00c3fa79695fc200470c897db246a8"),
        (144, 144, (), "4474d0bf355bbaf0d3ec55bd3d436f8000996130a6d9e6bb58c29f99c85717e9"),
        (0, 0, tuple(range(1, P)), "633ab4de4d7c45b84d18d15e75c27d5b4f00c3fa79695fc200470c897db246a8"),
        (144, 144, (), "910bb577eec7147c2f480de7dacdf222f409976c0920521288ea70bf94c23786"),
    ),
)
EXPECTED_FIXED = (
    (6, 5, (12, 12), "37daad2d2da203cc2ac57342765b379096aedcbf17b903bd286210a34cdfa356"),
    (6, 5, (12, 12), "43c2035d3635fa8309aa870a6f8e0e42e625da32e2573091879737a758fe1308"),
    (4, 3, (0, 0), "5cbc6f2b579968bf49df54873fc7ad55436c9137fd7bc9f5a703f5bcd8850310"),
)
EXPECTED_DIGESTS = (
    (
        "416f03a90894e526bc767a34839c2489aa4cac7051ac04687e0feacfa36d58d2",
        "a24603fe10909d5fec3f4a3668ae4a2eaffb06ba36c7d6a7121665fc8a996c38",
    ),
    (
        "9c5e227d9c142373973a562a54c6a67cac60a82da1121a028b9658920d155a19",
        "4231a2033312b13ebc3070c5084f55dff421186120831142dd8a57722f0d0c7b",
        "4958e3b37ff473644a610671763404d5aa2fe850c4e57de3c2c94eb0a89bb184",
    ),
    (
        "ae9744c987a3b2b6824feab0b64b234a239f7efa250193304e9aa5beb18e6fb5",
        "34356f952b40ee3f8ed423d732114c2bb2c1107dddf89bbe5eb3d956215d753f",
        "154bae6efd632d281fa96acafd386cde4d75ae38f460e706d02f0746371e3310",
    ),
    (
        "9d3733c61ce0dacc056d382db70ad0a10ad5aff62a8c57c3968229e21fdfad58",
        "4d9857743743203f86380b3a162c0e4314de119b95b3cbd137d0e3f180993453",
        "e59b7082fe2e25718b01c07d7d29a843d4e695075a23b85e83da960d70a8bd8d",
    ),
)
EXPECTED_PARENT_GAMMA = (
    "a83811099901421a16d39b1a04a9dd09bea96eba1b317208e58df1e6ab670ab6",
    "9076200eefb5d5322e8e20fdbb9fc805d17cc2deb287cc9e184e7e4a06d205d7",
)
EXPECTED_PARENT_TENSOR = (
    "0c6a6b831b5548d012c2954ed5fc11387657591e9de0f6ae3b55a1e593982f9f",
    "147490619b4a53c125140f3e57b0285e5153e6eeff6413410c5cfa9f769c8f7c",
)
EXPECTED_SEMANTIC_SHA256 = "66db2301f88db1ced7784868095e198e3e12f1fb79175ebd902d0f569a5decef"


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
    require(lf_sha256(PARENT_PATH) == PARENT_SHA256, "root-difference parent drift")
    spec = importlib.util.spec_from_file_location("pointed_parent", PARENT_PATH)
    require(spec is not None and spec.loader is not None, "pointed parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.EXPECTED_SEMANTIC_SHA256 == PARENT_SEMANTIC,
            "root-difference parent semantic drift")
    return module


RD = load_parent()
SQ = RD.SQ
R = RD.R
PRIME = R.JOINT_PRIME


def integrate_pair(
    alpha, beta, word, endpoint_grid, source_u, source_v,
    boundaries, harmonic, danger, literal_tau=None,
):
    events, interval_count, mapped_count = R.endpoint_events(
        word, endpoint_grid, alpha, beta, literal_tau
    )
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    banks = [[[[0] * P for _point in POINTS] for _tau in tau_values]
             for _bank in range(2)]
    point_counts = [0] * len(POINTS)
    positions = sorted(set(events) | boundaries)
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
        state = SQ.state_index(source_u, left)
        u_values = tuple(R.profile_value(profile, left) for profile in source_u)
        v_values = tuple(R.profile_value(profile, left) for profile in source_v)
        if not jump:
            continue
        q_active += 1
        if any(u_values) and any(v_values):
            weighted_intervals += 1
        for u, value in enumerate(u_values):
            if value:
                point = POINT_INDEX.get((state, u))
                require(point is not None, ("unrealized supported tail", state, u))
                point_counts[point] += 1
        name = R.chamber(left, right)
        require(name in ("left", "right"), ("pointed chamber", name))
        for row_index, tau in enumerate(tau_values):
            selected = mask if literal_tau is not None else (
                mask & R.guard_mask(name, tau, danger)
            )
            roots = tuple(root for root in range(P) if (selected >> root) & 1)
            for u in roots:
                for q in roots:
                    weighted = u_values[u] * v_values[q]
                    supported = int(u_values[u] != 0 and v_values[q] != 0)
                    if not weighted and not supported:
                        continue
                    point = POINT_INDEX.get((state, u))
                    require(point is not None, ("unrealized pointed state", state, u))
                    difference = (u - q) % P
                    banks[0][row_index][point][difference] = (
                        banks[0][row_index][point][difference] + weighted * jump
                    ) % PRIME
                    banks[1][row_index][point][difference] = (
                        banks[1][row_index][point][difference] + supported * jump
                    ) % PRIME
                    require(u != q, ("pointed same-root", alpha, beta, tau, left, u))
    mask ^= events.get(positions[-1], 0)
    require(mask == 0, ("pointed mask closure", alpha, beta, literal_tau))
    frozen = tuple(tuple(tuple(tuple(row) for row in tau_row) for tau_row in bank)
                   for bank in banks)
    return frozen, (
        interval_count, mapped_count, active, q_active, weighted_intervals
    ), tuple(point_counts)


def build_banks(word, endpoint_grid, source_u, source_v, boundaries, harmonic, danger):
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    gamma = [[], []]
    totals = [0] * 5
    point_totals = [0] * len(POINTS)
    for alpha in range(P):
        for beta in range(P):
            local, counts, point_counts = integrate_pair(
                alpha, beta, word, endpoint_grid, source_u, source_v,
                boundaries, harmonic, danger,
            )
            phase = pow(zeta, beta, PRIME)
            for bank in range(2):
                for tau in range(P):
                    gamma[bank].append(tuple(tuple(
                        phase * value % PRIME for value in local[bank][tau][point]
                    ) for point in range(len(POINTS))))
            totals = [left + right for left, right in zip(totals, counts)]
            point_totals = [left + right for left, right in zip(point_totals, point_counts)]
    return tuple(tuple(bank) for bank in gamma), tuple(totals), tuple(point_totals)


def inverse_tensor(gamma, zeta):
    size = len(POINTS)
    answer = [[[0] * P for _difference in range(P)] for _point in range(size)]
    normalizer = pow(P**3, -1, PRIME)
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                for relation in range(P):
                    phase = pow(zeta, -(alpha + tau * relation) % P, PRIME)
                    for point in range(size):
                        for difference in range(P):
                            answer[point][difference][relation] = (
                                answer[point][difference][relation]
                                + row[point][difference] * phase
                            ) % PRIME
    return tuple(tuple(tuple(value * normalizer % PRIME for value in line)
                       for line in sheet) for sheet in answer)


def parent_marginal(tensor):
    return tuple(tuple(tuple(sum(tensor[point][difference][relation]
                                 for point in STATE_FIBRES[state]) % PRIME
                             for relation in range(P)) for difference in range(P))
                 for state in range(4))


def flat_tail(parent):
    answer = []
    for state, fibre in enumerate(STATE_FIBRES):
        inv = pow(len(fibre), -1, PRIME)
        for _point in fibre:
            answer.append(tuple(tuple(parent[state][difference][relation] * inv % PRIME
                                      for relation in range(P))
                                for difference in range(P)))
    # The state-fibre order differs from geometric point order for states 2 and 3.
    by_pair = {}
    cursor = 0
    for state, fibre in enumerate(STATE_FIBRES):
        for point in fibre:
            by_pair[point] = answer[cursor]
            cursor += 1
    return tuple(by_pair[point] for point in range(len(POINTS)))


def centre_axes(tensor):
    size = len(tensor)
    data = [[[tensor[point][difference][relation] for relation in range(P)]
             for difference in range(P)] for point in range(size)]
    inv_size = pow(size, -1, PRIME)
    inv_p = pow(P, -1, PRIME)
    for difference in range(P):
        for relation in range(P):
            mean = sum(data[point][difference][relation] for point in range(size))
            mean = mean * inv_size % PRIME
            for point in range(size):
                data[point][difference][relation] = (
                    data[point][difference][relation] - mean
                ) % PRIME
    for point in range(size):
        for relation in range(P):
            mean = sum(data[point][difference][relation] for difference in range(P))
            mean = mean * inv_p % PRIME
            for difference in range(P):
                data[point][difference][relation] = (
                    data[point][difference][relation] - mean
                ) % PRIME
    for point in range(size):
        for difference in range(P):
            mean = sum(data[point][difference]) * inv_p % PRIME
            for relation in range(P):
                data[point][difference][relation] = (
                    data[point][difference][relation] - mean
                ) % PRIME
    return tuple(tuple(tuple(line) for line in sheet) for sheet in data)


def flattening_ranks(tensor):
    size = len(tensor)
    return (
        R.rank_mod(tuple(tuple(tensor[a][s][t] for s in range(P) for t in range(P))
                         for a in range(size))),
        R.rank_mod(tuple(tuple(tensor[a][s][t] for a in range(size) for t in range(P))
                         for s in range(P))),
        R.rank_mod(tuple(tuple(tensor[a][s][t] for a in range(size) for s in range(P))
                         for t in range(P))),
    )


def centre_difference_relation(tensor):
    size = len(tensor)
    inv_p = pow(P, -1, PRIME)
    data = [[[tensor[a][s][t] for t in range(P)] for s in range(P)]
            for a in range(size)]
    for a in range(size):
        for t in range(P):
            mean = sum(data[a][s][t] for s in range(P)) * inv_p % PRIME
            for s in range(P):
                data[a][s][t] = (data[a][s][t] - mean) % PRIME
        for s in range(P):
            mean = sum(data[a][s]) * inv_p % PRIME
            for t in range(P):
                data[a][s][t] = (data[a][s][t] - mean) % PRIME
    return tuple(tuple(tuple(line) for line in sheet) for sheet in data)


def edge_spectra(tensor, zeta):
    centred = centre_difference_relation(tensor)
    fourier = tuple(tuple(pow(zeta, -frequency * value % P, PRIME)
                          for value in range(P)) for frequency in range(P))
    spectra = []
    records = []
    for left, right in zip(range(len(POINTS) - 1), range(1, len(POINTS))):
        spectrum = tuple(tuple(sum(
            (centred[right][difference][relation]
             - centred[left][difference][relation])
            * fourier[difference_frequency][difference]
            * fourier[relation_frequency][relation]
            for difference in range(P) for relation in range(P)
        ) % PRIME for relation_frequency in range(P))
            for difference_frequency in range(P))
        all_support = sum(spectrum[k][ell] != 0 for k in range(P) for ell in range(P))
        mixed = sum(spectrum[k][ell] != 0
                    for k in range(1, P) for ell in range(1, P))
        zero_lines = tuple(ell for ell in range(1, P)
                           if all(spectrum[k][ell] == 0 for k in range(1, P)))
        spectra.append(spectrum)
        records.append((all_support, mixed, zero_lines, digest_json(spectrum)))
    return tuple(spectra), tuple(records)


def fixed_record(tensor):
    matrix = tuple(tuple(tensor[point][difference][6] for difference in range(P))
                   for point in range(len(POINTS)))
    inv_six = pow(len(POINTS), -1, PRIME)
    inv_p = pow(P, -1, PRIME)
    grand = sum(sum(row) for row in matrix) % PRIME
    row_sums = tuple(sum(row) % PRIME for row in matrix)
    column_sums = tuple(sum(matrix[point][difference] for point in range(len(POINTS)))
                        % PRIME for difference in range(P))
    centred = tuple(tuple((
        matrix[point][difference]
        - row_sums[point] * inv_p
        - column_sums[difference] * inv_six
        + grand * pow(len(POINTS) * P, -1, PRIME)
    ) % PRIME for difference in range(P)) for point in range(len(POINTS)))
    contrasts = ((1, 2), (3, 4))
    contrast_support = tuple(sum(
        (matrix[right][difference] - matrix[left][difference]) % PRIME != 0
        for difference in range(1, P)
    ) for left, right in contrasts)
    return R.rank_mod(matrix), R.rank_mod(centred), contrast_support, digest_json(matrix)


def main() -> None:
    R.split_field_certificate()
    source_u, source_v, source_boundaries, profile_digest, _total, _types = R.source_profiles()
    source_grid = R.SRC.T_DEN
    hierarchy = SQ.source_hierarchy(source_u, source_v, source_boundaries, source_grid)
    word, endpoint_grid = R.endpoint_word_and_grid()
    source_u_scaled = tuple(R.scale_profile(profile, R.JOINT_COORDINATE, source_grid)
                            for profile in source_u)
    source_v_scaled = tuple(R.scale_profile(profile, R.JOINT_COORDINATE, source_grid)
                            for profile in source_v)
    boundaries = R.fixed_boundaries(source_boundaries, source_grid)
    harmonic = R.HarmonicPrimitive(word, endpoint_grid)
    danger = R.danger_arcs()
    gamma, work_counts, point_counts = build_banks(
        word, endpoint_grid, source_u_scaled, source_v_scaled,
        boundaries, harmonic, danger,
    )
    require(work_counts == EXPECTED_WORK_COUNTS, ("pointed work counts", work_counts))
    require(point_counts == EXPECTED_POINT_SEGMENTS, ("point counts", point_counts))

    gamma_parents = tuple(tuple(tuple(tuple(sum(
        row[point][difference] for point in STATE_FIBRES[state]
    ) % PRIME for difference in range(P)) for state in range(4)) for row in bank)
        for bank in gamma)
    require(tuple(digest_json(bank) for bank in gamma_parents) == EXPECTED_PARENT_GAMMA,
            "pointed gamma marginal")

    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    guard_records = []
    for alpha, beta, tau in CONTROL_TRIPLES:
        R.literal_guard_restoration(word, endpoint_grid, alpha, beta, tau, boundaries, danger)
        direct, _counts, _points = integrate_pair(
            alpha, beta, word, endpoint_grid, source_u_scaled, source_v_scaled,
            boundaries, harmonic, danger, literal_tau=tau,
        )
        phase = pow(zeta, beta, PRIME)
        index = (alpha * P + beta) * P + tau
        for bank in range(2):
            expected = tuple(tuple(phase * value % PRIME for value in direct[bank][0][point])
                             for point in range(len(POINTS)))
            require(expected == gamma[bank][index],
                    ("literal pointed guard", bank, alpha, beta, tau))
        guard_records.append((alpha, beta, tau))

    tensors_real = tuple(inverse_tensor(bank, zeta) for bank in gamma)
    parents = tuple(parent_marginal(tensor) for tensor in tensors_real)
    require(tuple(digest_json(parent) for parent in parents) == EXPECTED_PARENT_TENSOR,
            "pointed tensor marginal")
    flat = flat_tail(parents[0])
    tensors = tensors_real + (flat,)
    require(parent_marginal(flat) == parents[0], "flat-tail parent marginal")
    centred = tuple(centre_axes(tensor) for tensor in tensors)
    ranks = tuple(flattening_ranks(tensor) for tensor in tensors)
    centred_ranks = tuple(flattening_ranks(tensor) for tensor in centred)
    require(ranks == EXPECTED_RANKS, ("pointed ranks", ranks))
    require(centred_ranks == EXPECTED_ANOVA_RANKS, ("pointed ANOVA ranks", centred_ranks))

    difference_support = tuple(tuple(s for s in range(P)
        if any(tensor[point][s][t] != 0 for point in range(len(POINTS)) for t in range(P)))
        for tensor in tensors)
    require(difference_support == (tuple(range(1, P)),) * 3,
            ("pointed difference support", difference_support))
    edge_data = tuple(edge_spectra(tensor, zeta) for tensor in tensors)
    edge_spectra_values = tuple(item[0] for item in edge_data)
    edge_records = tuple(item[1] for item in edge_data)
    require(edge_records == EXPECTED_EDGE_RECORDS, ("path edge records", edge_records))
    fixed = tuple(fixed_record(tensor) for tensor in tensors)
    require(fixed == EXPECTED_FIXED, ("pointed fixed records", fixed))
    digests = (
        tuple(digest_json(bank) for bank in gamma),
        tuple(digest_json(tensor) for tensor in tensors),
        tuple(digest_json(tensor) for tensor in centred),
        tuple(digest_json(spectra) for spectra in edge_spectra_values),
    )
    require(digests == EXPECTED_DIGESTS, ("pointed candidate digests", digests))
    record = (
        profile_digest, hierarchy, POINTS, STATE_FIBRES, tuple(guard_records),
        work_counts, point_counts, difference_support, ranks, centred_ranks,
        edge_records, fixed, digests,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("pointed semantic", semantic))

    print("== independent hostile audit: six pointed source cuts x root difference x relation ==")
    print(f"parent=(sha256={PARENT_SHA256},semantic={PARENT_SEMANTIC})")
    print(f"pointed_path={POINTS};state_fibres={STATE_FIBRES}")
    print(f"source_profile_sha256={profile_digest};field=(prime={PRIME},zeta13={zeta})")
    print(f"work_counts={work_counts};point_segment_counts={point_counts};literal_controls={tuple(guard_records)}: PASS")
    print("pair_typing=point=(state,u) with actual U-tail;difference=u-q;same-root zero: PASS")
    print("marginals=(weighted,support_only)->audited root-difference parent;flat_tail_same_parent: PASS")
    print(f"difference_support_(weighted,support_only,flat_tail)={difference_support}")
    print(f"flattening_ranks_(pointed,difference,relation)={ranks}")
    print(f"three_way_ANOVA_flattening_ranks={centred_ranks}")
    print(f"path_edge_support_(all,mixed,zero_relation_lines,digest)={edge_records}")
    print(f"fixed_relation=(1,0,6);records_(rank,ANOVA_rank,tail_contrast_support,digest)={fixed}")
    print(f"digests=(gamma,tensor,ANOVA,path_edge_spectra)={digests}")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT scoped six-pointed-state source-tail refinement")
    print("typing=lawful absolute source tail;not cyclic C6,not exact address,not inverse ancestry,not chronology")
    print("scope=one owner base;no physical current,no row exclusion,no LRC(14)")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
