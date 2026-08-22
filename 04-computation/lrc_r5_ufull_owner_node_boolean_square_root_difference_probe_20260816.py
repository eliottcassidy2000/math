#!/usr/bin/env python3
"""Retain source root difference inside the owner-node Boolean square.

MISTAKE-417 killed the inherited seven-cell axis, while the pinned
Boolean-square refiner recovered a genuine four-state coordinate.  The
independent common-connection audit identifies ``(u,s,ell,theta)`` as the
smallest source object surviving all failed post-marginal fits.  This probe
therefore keeps the lawful source root difference

    s = u-q in F_13

before multiplying and integrating the actual U_full endpoint factors at
the linked owner nodes ``w_u`` and ``w_q``.  It then inverts the endpoint
character bank to the thirteen refined relation residues and tests the full

    Boolean square x F_13(root difference) x F_13(relation)

tensor.  Sequential three-axis centering is the decisive hostile: formal
Fourier support alone is not treated as interaction.

The root difference is a genuine THM-2471 source-service label.  It is not a
THM-2334 exact relation address, an inverse-ancestry sheet, a chronological
U_clock lift, or a physical current.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SQUARE_PATH = (
    ROOT
    / "04-computation/lrc_r5_ufull_owner_node_boolean_square_refiner_probe_20260816.py"
)
SQUARE_SHA256 = "c7bbd2d82ed067914f39253ac52ce32f9e179a859fb55a421bd2b390f707d881"

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
EXPECTED_CENSUSES = (
    (
        (676, 1, 12, 12, 144, 3, 36, 36, 432),
        (676, 1, 12, 12, 144, 3, 36, 36, 432),
        (676, 1, 12, 12, 144, 3, 36, 36, 432),
    ),
    (
        (432, 0, 0, 0, 0, 0, 0, 0, 432),
        (432, 0, 0, 0, 0, 0, 0, 0, 432),
        (432, 0, 0, 0, 0, 0, 0, 0, 432),
    ),
    (
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    ),
)
EXPECTED_RANKS = (
    ((4, 12, 13), (4, 12, 13), (4, 7, 13)),
    ((3, 12, 12), (3, 12, 12), (3, 6, 12)),
)
EXPECTED_FIXED = (
    (4, 3, (10, 12, 10, 12), ((0, 1, 2), (0,), (0, 11, 12), (0,)),
     "acde8efe1065b919716bacf43fbb09afd8dd491b219b67793d8a94f5386eaf5f"),
    (4, 3, (10, 12, 10, 12), ((0, 1, 2), (0,), (0, 11, 12), (0,)),
     "fedbc2430a29f32ff39ec9e7948bde74e896c5d4d2a676fc1e9fdbe395d66c16"),
    (4, 3, (13, 13, 13, 13), ((), (), (), ()),
     "e5a9460cca848458272f5a8917f2786096dd8939ab40ce28e4fb13a5fa9576c7"),
)
EXPECTED_SEMANTIC_SHA256 = "38f2cbc38bbeee0a8556e0649c4cd608645cba1477e78f3d436b860b4d696278"

P = 13
V = 4


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
    observed = lf_sha256(SQUARE_PATH)
    require(observed == SQUARE_SHA256,
            ("Boolean-square source drift", str(SQUARE_PATH), observed, SQUARE_SHA256))
    spec = importlib.util.spec_from_file_location("owner_square_difference_parent", SQUARE_PATH)
    require(spec is not None and spec.loader is not None, "Boolean-square loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_square()
C = B.C
PRIME = C.JOINT_PRIME


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def zero_cube():
    return [[[0 for _difference in range(P)] for _state in range(V)] for _tau in range(P)]


def integrate_refined(alpha: int, beta: int, literal_tau: int | None = None):
    ctx = C.context()
    events, interval_count, mapped = C.endpoint_events(alpha, beta, literal_tau)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    weighted = [
        [[0 for _difference in range(P)] for _state in range(V)]
        for _tau in tau_values
    ]
    support_only = [
        [[0 for _difference in range(P)] for _state in range(V)]
        for _tau in tau_values
    ]
    endpoint_only = [
        [[0 for _difference in range(P)] for _state in range(V)]
        for _tau in tau_values
    ]
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
                ("root-difference owner escaped cell zero", alpha, beta, left, right))
        chamber = C.chamber_of_segment(left, right)
        require(chamber in ("left", "right"),
                ("root-difference active middle chamber", alpha, beta, left, right))
        state = B.state_of_segment(left, right)
        active_segments += 1
        state_segments[state] += 1
        jump = C.q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        u_values = tuple(profile_value(profile, left) for profile in ctx["source_u"])
        v_values = tuple(profile_value(profile, left) for profile in ctx["source_v"])
        u_support = tuple(root for root, value in enumerate(u_values) if value)
        v_support = tuple(root for root, value in enumerate(v_values) if value)
        if not u_support or not v_support:
            continue
        weighted_segments += 1
        require(set(u_support).isdisjoint(v_support),
                ("source supports overlap", left, right, u_support, v_support))
        for row_index, tau in enumerate(tau_values):
            if literal_tau is None:
                selected = tuple(
                    sheet
                    for sheet in range(P)
                    if (mask >> sheet) & 1 and C.T.safe(chamber, sheet + tau)
                )
            else:
                selected = tuple(sheet for sheet in range(P) if (mask >> sheet) & 1)
            selected_set = set(selected)
            selected_u = tuple(root for root in u_support if root in selected_set)
            selected_v = tuple(root for root in v_support if root in selected_set)
            for u in selected_u:
                for q in selected_v:
                    difference = (u - q) % P
                    require(difference != 0, ("same-root support pair", left, right, u, q))
                    weighted[row_index][state][difference] = (
                        weighted[row_index][state][difference]
                        + u_values[u] * v_values[q] * jump
                    ) % PRIME
                    support_only[row_index][state][difference] = (
                        support_only[row_index][state][difference] + jump
                    ) % PRIME
            for u in selected:
                for q in selected:
                    difference = (u - q) % P
                    endpoint_only[row_index][state][difference] = (
                        endpoint_only[row_index][state][difference] + jump
                    ) % PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("root-difference endpoint mask", alpha, beta, literal_tau, mask))
    counts = (
        interval_count,
        mapped,
        active_segments,
        q_active_segments,
        weighted_segments,
        tuple(state_segments),
    )
    freeze = lambda bank: tuple(
        tuple(tuple(row) for row in state_rows) for state_rows in bank
    )
    return freeze(weighted), freeze(support_only), freeze(endpoint_only), counts


def worker(alpha: int):
    zeta = C.context()["zeta"]
    weighted_rows = []
    support_rows = []
    endpoint_rows = []
    scalar_counts = [0] * 5
    state_counts = [0] * V
    for beta in range(P):
        weighted, support_only, endpoint_only, record = integrate_refined(alpha, beta)
        phase = pow(zeta, beta, PRIME)

        def twist(bank):
            return tuple(
                tuple(
                    tuple(phase * value % PRIME for value in difference_row)
                    for difference_row in state_rows
                )
                for state_rows in bank
            )

        weighted_rows.append(twist(weighted))
        support_rows.append(twist(support_only))
        endpoint_rows.append(twist(endpoint_only))
        scalar_counts = [left + right for left, right in zip(scalar_counts, record[:5])]
        state_counts = [left + right for left, right in zip(state_counts, record[5])]
    return (
        alpha,
        tuple(weighted_rows),
        tuple(support_rows),
        tuple(endpoint_rows),
        (tuple(scalar_counts), tuple(state_counts)),
    )


def inverse_tensor(gamma_rows, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    tensor = [[[0 for _relation in range(P)] for _difference in range(P)] for _state in range(V)]
    index = 0
    phases = tuple(
        tuple(pow(zeta, -(alpha + tau * relation) % P, PRIME) for relation in range(P))
        for alpha in range(P) for tau in range(P)
    )
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma_rows[index]
                phase_row = phases[alpha * P + tau]
                index += 1
                for state in range(V):
                    for difference in range(P):
                        value = row[state][difference]
                        if not value:
                            continue
                        for relation, phase in enumerate(phase_row):
                            tensor[state][difference][relation] = (
                                tensor[state][difference][relation] + value * phase
                            ) % PRIME
    require(index == P**3, ("root-difference gamma size", index))
    return tuple(
        tuple(
            tuple(value * normalizer % PRIME for value in relation_row)
            for relation_row in difference_rows
        )
        for difference_rows in tensor
    )


def fourier_tensor(tensor, zeta: int):
    roots = tuple(
        tuple(pow(zeta, -frequency * value % P, PRIME) for value in range(P))
        for frequency in range(P)
    )
    walsh_stage = [[
        [
            sum(
                tensor[state][difference][relation] * B.WALSH_SIGNS[character][state]
                for state in range(V)
            ) % PRIME
            for relation in range(P)
        ]
        for difference in range(P)
    ] for character in range(V)]
    difference_stage = [[
        [
            sum(walsh_stage[character][difference][relation] * roots[frequency][difference]
                for difference in range(P)) % PRIME
            for relation in range(P)
        ]
        for frequency in range(P)
    ] for character in range(V)]
    return tuple(
        tuple(
            tuple(
                sum(difference_stage[character][difference_frequency][relation]
                    * roots[relation_frequency][relation]
                    for relation in range(P)) % PRIME
                for relation_frequency in range(P)
            )
            for difference_frequency in range(P)
        )
        for character in range(V)
    )


def support_census(spectrum):
    bins = [0] * 8
    for character in range(V):
        for difference_frequency in range(P):
            for relation_frequency in range(P):
                if spectrum[character][difference_frequency][relation_frequency]:
                    mask = ((character != 0) << 2) | ((difference_frequency != 0) << 1) | (relation_frequency != 0)
                    bins[mask] += 1
    return (sum(bins),) + tuple(bins)


def centre_axis(tensor, axis: int):
    dimensions = (V, P, P)
    inverse = pow(dimensions[axis], -1, PRIME)
    mutable = [[[tensor[state][difference][relation] for relation in range(P)]
                for difference in range(P)] for state in range(V)]
    if axis == 0:
        for difference in range(P):
            for relation in range(P):
                mean = sum(mutable[state][difference][relation] for state in range(V)) * inverse % PRIME
                for state in range(V):
                    mutable[state][difference][relation] = (mutable[state][difference][relation] - mean) % PRIME
    elif axis == 1:
        for state in range(V):
            for relation in range(P):
                mean = sum(mutable[state][difference][relation] for difference in range(P)) * inverse % PRIME
                for difference in range(P):
                    mutable[state][difference][relation] = (mutable[state][difference][relation] - mean) % PRIME
    else:
        for state in range(V):
            for difference in range(P):
                mean = sum(mutable[state][difference][relation] for relation in range(P)) * inverse % PRIME
                for relation in range(P):
                    mutable[state][difference][relation] = (mutable[state][difference][relation] - mean) % PRIME
    return tuple(tuple(tuple(row) for row in plane) for plane in mutable)


def three_way_interaction(tensor):
    return centre_axis(centre_axis(centre_axis(tensor, 0), 1), 2)


def flattening_ranks(tensor):
    state_flat = tuple(
        tuple(tensor[state][difference][relation] for difference in range(P) for relation in range(P))
        for state in range(V)
    )
    difference_flat = tuple(
        tuple(tensor[state][difference][relation] for state in range(V) for relation in range(P))
        for difference in range(P)
    )
    relation_flat = tuple(
        tuple(tensor[state][difference][relation] for state in range(V) for difference in range(P))
        for relation in range(P)
    )
    return C.rank_mod(state_flat), C.rank_mod(difference_flat), C.rank_mod(relation_flat)


def state_relation_marginal(tensor):
    return tuple(
        tuple(sum(tensor[state][difference][relation] for difference in range(P)) % PRIME
              for relation in range(P))
        for state in range(V)
    )


def fixed_relation_matrix(tensor, relation: int):
    return tuple(
        tuple(tensor[state][difference][relation] for difference in range(P))
        for state in range(V)
    )


def cut_arc_profile():
    """Difference census of the four raw directed cuts S to F_13 minus S."""

    profile = tuple(
        tuple(
            sum(
                1
                for u in subset for q in range(P)
                if q not in subset and (u - q) % P == difference
            )
            for difference in range(P)
        )
        for subset in B.STATE_SUBSETS
    )
    expected_singleton = (0,) + (1,) * 12
    expected_doubleton = (0,) + tuple(
        1 if difference in (6, 7) else 2 for difference in range(1, P)
    )
    require(profile == (
        expected_singleton,
        expected_doubleton,
        expected_singleton,
        expected_doubleton,
    ), ("directed-cut arc profile", profile))
    require(C.rank_mod(profile) == 2, ("directed-cut arc rank", profile))
    return profile


def matrix_interaction(matrix):
    rows = len(matrix)
    columns = len(matrix[0])
    inverse_rows = pow(rows, -1, PRIME)
    inverse_columns = pow(columns, -1, PRIME)
    inverse_total = pow(rows * columns, -1, PRIME)
    row_sums = tuple(sum(row) % PRIME for row in matrix)
    column_sums = tuple(sum(matrix[row][column] for row in range(rows)) % PRIME
                        for column in range(columns))
    grand = sum(row_sums) % PRIME
    return tuple(
        tuple((matrix[row][column]
               - row_sums[row] * inverse_columns
               - column_sums[column] * inverse_rows
               + grand * inverse_total) % PRIME
              for column in range(columns))
        for row in range(rows)
    )


def main() -> None:
    ctx = C.context()
    arc_profile = cut_arc_profile()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "root-difference worker order")

    def flatten_bank(position: int):
        return tuple(
            row
            for chunk in chunks
            for beta_rows in chunk[position]
            for row in beta_rows
        )

    weighted_gamma = flatten_bank(1)
    support_gamma = flatten_bank(2)
    endpoint_gamma = flatten_bank(3)
    require(len(weighted_gamma) == len(support_gamma) == len(endpoint_gamma) == P**3,
            "root-difference bank size")

    weighted_gamma_marginal = tuple(
        tuple(sum(row[state]) % PRIME for state in range(V)) for row in weighted_gamma
    )
    endpoint_gamma_marginal = tuple(
        tuple(sum(row[state]) % PRIME for state in range(V)) for row in endpoint_gamma
    )
    require(digest_json(weighted_gamma_marginal) == B.EXPECTED_DIGESTS[0],
            "root-difference weighted marginal misses Boolean square")
    require(digest_json(endpoint_gamma_marginal) == B.EXPECTED_DIGESTS[1],
            "root-difference endpoint marginal misses square erasure")
    require(all(row[state][0] == 0 for row in weighted_gamma for state in range(V)),
            "weighted same-root slice")
    require(all(row[state][0] == 0 for row in support_gamma for state in range(V)),
            "support-only same-root slice")
    require(any(row[state][0] != 0 for row in endpoint_gamma for state in range(V)),
            "endpoint-only same-root hostile failed to fire")

    tensors = tuple(
        inverse_tensor(bank, ctx["zeta"])
        for bank in (weighted_gamma, support_gamma, endpoint_gamma)
    )
    require(digest_json(state_relation_marginal(tensors[0])) == B.EXPECTED_DIGESTS[2],
            "weighted inverse marginal misses Boolean square")
    require(digest_json(state_relation_marginal(tensors[2])) == B.EXPECTED_DIGESTS[3],
            "endpoint inverse marginal misses square erasure")

    interactions = tuple(three_way_interaction(tensor) for tensor in tensors)
    spectra = tuple(fourier_tensor(tensor, ctx["zeta"]) for tensor in tensors)
    interaction_spectra = tuple(
        fourier_tensor(tensor, ctx["zeta"]) for tensor in interactions
    )
    censuses = tuple(support_census(spectrum) for spectrum in spectra)
    interaction_censuses = tuple(
        support_census(spectrum) for spectrum in interaction_spectra
    )
    ranks = tuple(flattening_ranks(tensor) for tensor in tensors)
    interaction_ranks = tuple(flattening_ranks(tensor) for tensor in interactions)

    relation = 6
    fixed_matrices = tuple(fixed_relation_matrix(tensor, relation) for tensor in tensors)
    fixed_records = tuple(
        (
            C.rank_mod(matrix),
            C.rank_mod(matrix_interaction(matrix)),
            tuple(sum(matrix[state][difference] != 0 for difference in range(P))
                  for state in range(V)),
            tuple(tuple(difference for difference in range(P)
                        if matrix[state][difference] == 0)
                  for state in range(V)),
            digest_json(matrix),
        )
        for matrix in fixed_matrices
    )
    for record in fixed_records[:2]:
        zero_sets = record[3]
        for state in range(V):
            reflected = tuple(sorted((-difference) % P for difference in zero_sets[state]))
            require(reflected == zero_sets[state ^ 2],
                    ("fixed-relation chamber reflection", state, zero_sets))
    difference_support = tuple(
        tuple(
            difference
            for difference in range(P)
            if any(tensor[state][difference][relation_value] != 0
                   for state in range(V) for relation_value in range(P))
        )
        for tensor in tensors
    )

    scalar_counts = tuple(sum(chunk[4][0][index] for chunk in chunks) for index in range(5))
    state_counts = tuple(sum(chunk[4][1][state] for chunk in chunks) for state in range(V))
    require(scalar_counts == C.EXPECTED_WORK_COUNTS,
            ("root-difference work counts", scalar_counts))
    require(state_counts == B.EXPECTED_STATE_SEGMENTS,
            ("root-difference state counts", state_counts))

    direct_controls = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
    direct_record = []
    for alpha, beta, tau in direct_controls:
        direct = integrate_refined(alpha, beta, tau)
        index = (alpha * P + beta) * P + tau
        phase = pow(ctx["zeta"], beta, PRIME)
        for bank_index, gamma_bank in enumerate((weighted_gamma, support_gamma, endpoint_gamma)):
            direct_row = tuple(
                tuple(phase * value % PRIME for value in difference_row)
                for difference_row in direct[bank_index][0]
            )
            require(direct_row == gamma_bank[index],
                    ("root-difference literal guard", alpha, beta, tau, bank_index))
        direct_record.append(((alpha, beta, tau), direct[3]))

    digests = (
        tuple(digest_json(bank) for bank in (weighted_gamma, support_gamma, endpoint_gamma)),
        tuple(digest_json(tensor) for tensor in tensors),
        tuple(digest_json(spectrum) for spectrum in spectra),
        tuple(digest_json(tensor) for tensor in interactions),
        tuple(digest_json(spectrum) for spectrum in interaction_spectra),
    )
    if EXPECTED_DIGESTS != "TO_BE_PINNED":
        require(digests == EXPECTED_DIGESTS, ("root-difference digests", digests))
    if EXPECTED_CENSUSES != "TO_BE_PINNED":
        require((censuses, interaction_censuses, difference_support) == EXPECTED_CENSUSES,
                ("root-difference support", censuses, interaction_censuses, difference_support))
    if EXPECTED_RANKS != "TO_BE_PINNED":
        require((ranks, interaction_ranks) == EXPECTED_RANKS,
                ("root-difference ranks", ranks, interaction_ranks))
    if EXPECTED_FIXED != "TO_BE_PINNED":
        require(fixed_records == EXPECTED_FIXED,
                ("root-difference fixed relation", fixed_records))

    record = (
        SQUARE_SHA256,
        C.EXPECTED_SEMANTIC_SHA256,
        B.EXPECTED_SEMANTIC_SHA256,
        C.JOINT_PRIME,
        C.JOINT_ROOT,
        ctx["zeta"],
        scalar_counts,
        state_counts,
        tuple(direct_record),
        digests,
        censuses,
        interaction_censuses,
        ranks,
        interaction_ranks,
        difference_support,
        arc_profile,
        relation,
        fixed_records,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("root-difference semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 owner-node Boolean-square x root-difference x relation probe ==")
    print(f"parents=(square_sha={SQUARE_SHA256},square_semantic={B.EXPECTED_SEMANTIC_SHA256},owner_semantic={C.EXPECTED_SEMANTIC_SHA256})")
    print("coordinates=(state in V4,source_difference=u-q in F13,relation=(1,0,t) in F13)")
    print(f"field=(prime={C.JOINT_PRIME},root={C.JOINT_ROOT},zeta13={ctx['zeta']})")
    print(f"work_counts={scalar_counts};state_segment_counts={state_counts};literal_controls={direct_controls}: PASS")
    print("marginals=(weighted->Boolean_square,endpoint_only->source_erasure): PASS")
    print("same_root=(weighted=0,support_only=0,endpoint_only_nonzero): PASS")
    print(f"raw_cut_arc_profile=(rank=2,rows={arc_profile})")
    print(f"difference_support_(weighted,support_only,endpoint_only)={difference_support}")
    print(f"flattening_ranks_(state,difference,relation)=(weighted,support_only,endpoint_only)={ranks}")
    print(f"three_way_ANOVA_flattening_ranks={interaction_ranks}")
    print("support_census_order=(total,empty,relation,difference,difference_relation,state,state_relation,state_difference,triple)")
    print(f"spectral_censuses=(weighted,support_only,endpoint_only)={censuses}")
    print(f"three_way_ANOVA_spectral_censuses={interaction_censuses}")
    print(f"fixed_relation=(1,0,{relation});records_(rank,ANOVA_rank,row_support,zero_sets,digest)={fixed_records}")
    print(f"digests=(gamma,tensor,spectrum,ANOVA,ANOVA_spectrum)={digests}")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT three-coordinate source-root/Boolean/endpoint candidate on one owner base")
    print("scope=root difference is lawful source label;not exact address,not inverse ancestry,not chronology,not physical current,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
