#!/usr/bin/env python3
"""Independent hostile audit of the current-leg inverse-owner digit.

The submitted inverse-owner-branch probe is never imported.  Starting from
the independently audited owner-node Boolean square, this script rebuilds the
last current-leg inverse digit

    X_(u,a) = (w_u+a)/13^5,       r_owner = a mod 13,

by factoring the source fold as 13^5 = 13 * 13^4.  It folds by 13^4, pulls
the thirteen r_owner windows, and only then pulls the collision-root windows.
Their pointwise sum must recover the audited root-service profile.

The coordinate is not the source-time P_169 high digit, source root
difference, deep label, THM-2334 exact address, U_clock chronology, or a
physical current.  No row exclusion or LRC(14) conclusion is asserted.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
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
CONTROL_TRIPLES = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_STATE_COUNTS = (300115, 300116, 300115, 300116)
EXPECTED_BRANCH_COUNTS = (186244,) * P
EXPECTED_SOURCE_SHA256 = "9e2a61ea8ec00830aa93ec4daa42574a9179f5aa80a09e3388dd9af0964f4cb1"
EXPECTED_PARENT_GAMMA = "8b8f2a2785b084e1578ba0512e4577ab79fd674b84b588fbdb8186f2009242c2"
EXPECTED_PARENT_TABLE = "5f4b9609faaa5f148d112a7cde5cfba0ab2c1385b4c53ea9c4bcfc6e93d106fc"
EXPECTED_DIGESTS = (
    "344fc6f451922bfdd071cfa800444929e34e575a71344b17624cc0e73d7d9e72",
    "2a195fac7fbd60a4bbd2597bf34baf87302ead16c7c0e8c67c0667b0d320dfba",
    "8959d6734c5049dc58a6bc1ef702d20b839bed45e5d80321161a232ed367a228",
    (
        "d6d261c661a622ae315f3b713db02c23250eafbf20a38237f479f74ff99f61b0",
        "5209d2df83b48c3da26ed542054680bf39496125c2fe01aeb57617c451bcb907",
    ),
    (
        "0dca6acb11ce79478f690b509223db0ddd521aed95d8ddca7fc3cfca9e28c387",
        "197447345f85eb9bcafa5ca23a56870cb07d8ba0635314ff498abf7f9ab9d5bc",
    ),
    (
        "f2a19b58b81112fc6ec56571ebf36658731dcc026c27bfd1311ed8422742330d",
        "197447345f85eb9bcafa5ca23a56870cb07d8ba0635314ff498abf7f9ab9d5bc",
    ),
)
EXPECTED_RANKS = (((4, 4, 6), (4, 1, 4)), ((3, 4, 6), (0, 0, 0)))
EXPECTED_CENSUSES = (
    (
        (676, 1, 12, 12, 144, 3, 36, 36, 432),
        (52, 1, 12, 0, 0, 3, 36, 0, 0),
    ),
    (
        (432, 0, 0, 0, 0, 0, 0, 0, 432),
        (0, 0, 0, 0, 0, 0, 0, 0, 0),
    ),
)
EXPECTED_FIXED = (
    (4, 3, (13, 13, 13, 13),
     "7652ac0fbca27c3a41181d04e9bea962f4159b04ccc43befe3d2783af89be4a4"),
    (1, 0, (13, 13, 13, 13),
     "be782cc3ee60d8d4bc008f21a9f308731a29e180bc83d7041ff87c47eb532729"),
)
EXPECTED_SEMANTIC_SHA256 = "7063720d0e0e4847ce752102de83274ea47d7740fc435a64bae425dbd7100121"


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
    require(lf_sha256(SQUARE_PATH) == SQUARE_SHA256, "square parent hash drift")
    spec = importlib.util.spec_from_file_location("owner_digit_square_parent", SQUARE_PATH)
    require(spec is not None and spec.loader is not None, "square parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.EXPECTED_SEMANTIC_SHA256 == SQUARE_SEMANTIC,
            "square parent semantic drift")
    return module


SQ = load_square()
R = SQ.R
PRIME = R.JOINT_PRIME
WALSH_SIGNS = tuple(
    tuple(-1 if (character[0] * label[0] + character[1] * label[1]) % 2 else 1
          for label in SQ.STATE_LABELS)
    for character in SQ.STATE_LABELS
)


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def profile_mass(profile, grid: int) -> int:
    starts, values = profile
    return sum(
        values[index]
        * ((starts[index + 1] if index + 1 < len(starts) else grid) - left)
        for index, left in enumerate(starts)
    )


def radix_certificate():
    depth = P**5
    counts = [0] * P
    reflected = [0] * P
    for a in range(depth):
        quotient, digit = divmod(a, P)
        require(a == digit + P * quotient, ("radix decomposition", a))
        opposite = depth - 1 - a
        opposite_quotient, opposite_digit = divmod(opposite, P)
        require(opposite_digit == P - 1 - digit,
                ("reflected owner digit", a, digit, opposite_digit))
        require(opposite_quotient == P**4 - 1 - quotient,
                ("reflected upper word", a, quotient, opposite_quotient))
        counts[digit] += 1
        reflected[opposite_digit] += 1
    require(tuple(counts) == (P**4,) * P, ("digit census", counts))
    require(tuple(reflected) == tuple(counts), ("reflected census", reflected))
    # 13^4 X_(u,a) = (w_u + a)/13 = (w_u + digit)/13 mod 1.
    return depth, P**4, tuple(counts)


@lru_cache(maxsize=1)
def owner_digit_profiles():
    grid = R.SRC.T_DEN
    e_intervals = R.SRC.build_set(R.SRC.PAT_E, R.SRC.ZELL)
    q_intervals = R.SRC.build_set(R.SRC.PAT_QB, R.SRC.ZELL)
    packet = R.fold_weighted(
        [(left, right, 1) for left, right in e_intervals], R.SRC.RPKT, grid
    )
    f_pieces = R.intersect_profile_with_set(packet, q_intervals, grid)

    upper = R.fold_weighted(f_pieces, R.SRC.DCOLL // P, grid)
    digit_profiles = tuple(
        R.pull_profile_to_root(upper, digit, grid) for digit in range(P)
    )
    raw = tuple(
        tuple(R.pull_profile_to_root(digit_profiles[digit], root, grid)
              for digit in range(P))
        for root in range(P)
    )
    scale = R.JOINT_COORDINATE // grid
    scaled = tuple(
        tuple((tuple(point * scale for point in profile[0]), profile[1])
              for profile in root_profiles)
        for root_profiles in raw
    )

    source_u, source_v, source_boundaries, source_digest, _total, _types = (
        R.source_profiles()
    )
    source_u_scaled = tuple(
        R.scale_profile(profile, R.JOINT_COORDINATE, grid) for profile in source_u
    )
    source_v_scaled = tuple(
        R.scale_profile(profile, R.JOINT_COORDINATE, grid) for profile in source_v
    )
    boundaries = tuple(sorted(
        {0, R.JOINT_COORDINATE}
        | set(point * scale for point in source_boundaries)
        | {point for root_profiles in scaled for profile in root_profiles
           for point in profile[0]}
    ))
    require(tuple(R.JOINT_COORDINATE - point for point in reversed(boundaries))
            == boundaries, "owner-digit boundary reflection")

    reflection_intervals = 0
    square_reflection_intervals = 0
    for left, right in zip(boundaries, boundaries[1:]):
        reflected_left = R.JOINT_COORDINATE - right
        for root in range(P):
            restored = sum(profile_value(scaled[root][digit], left)
                           for digit in range(P))
            require(restored == profile_value(source_u_scaled[root], left),
                    ("rootwise digit restoration", left, root, restored))
            for digit in range(P):
                observed = profile_value(scaled[root][digit], left)
                reflected_value = profile_value(
                    scaled[P - 1 - root][P - 1 - digit], reflected_left
                )
                require(observed == reflected_value,
                        ("profile reflection", left, right, root, digit))
        reflection_intervals += 1

        support = frozenset(
            root for root, profile in enumerate(source_u_scaled)
            if profile_value(profile, left)
        )
        reflected_support = frozenset(
            root for root, profile in enumerate(source_u_scaled)
            if profile_value(profile, reflected_left)
        )
        expected_reflected_support = frozenset(P - 1 - root for root in support)
        require(reflected_support == expected_reflected_support,
                ("source support reflection", left, support, reflected_support))
        if support in SQ.STATE_INDEX:
            state = SQ.STATE_INDEX[support]
            require(reflected_support in SQ.STATE_INDEX,
                    ("reflected square state", left, reflected_support))
            require(SQ.STATE_INDEX[reflected_support] == state ^ 2,
                    ("state XOR2 reflection", left, state, reflected_support))
            square_reflection_intervals += 1

    branch_masses = tuple(
        sum(profile_mass(raw[root][digit], grid) for root in range(P))
        for digit in range(P)
    )
    require(branch_masses == tuple(reversed(branch_masses)),
            ("owner-digit mass reflection", branch_masses))
    source_record = (
        R.SRC.DCOLL,
        R.SRC.DCOLL // P,
        len(upper[0]),
        tuple(len(profile[0]) for profile in digit_profiles),
        tuple(tuple(len(profile[0]) for profile in row) for row in raw),
        len(boundaries),
        branch_masses,
        digest_json(raw),
        digest_json(scaled),
    )
    source_sha = digest_json((source_record, boundaries))
    require(source_sha == EXPECTED_SOURCE_SHA256,
            ("source decomposition", source_sha, EXPECTED_SOURCE_SHA256))
    reflection_record = (
        reflection_intervals,
        square_reflection_intervals,
        branch_masses,
    )
    return (
        scaled, boundaries, source_u_scaled, source_v_scaled,
        source_boundaries, source_digest, source_record, source_sha,
        reflection_record,
    )


@lru_cache(maxsize=1)
def endpoint_context():
    word, endpoint_grid = R.endpoint_word_and_grid()
    return word, endpoint_grid, R.HarmonicPrimitive(word, endpoint_grid), R.danger_arcs()


def integrate_pair(alpha: int, beta: int, literal_tau: int | None = None):
    (
        branches, branch_boundaries, source_u, source_v,
        _source_boundaries, _source_digest, _source_record, _source_sha,
        _reflection_record,
    ) = owner_digit_profiles()
    word, endpoint_grid, harmonic, danger = endpoint_context()
    events, interval_count, mapped_count = R.endpoint_events(
        word, endpoint_grid, alpha, beta, literal_tau
    )
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    bank = [[[0] * P for _state in range(V)] for _tau in tau_values]
    positions = sorted(
        set(events) | set(branch_boundaries)
        | {point for profile in source_v for point in profile[0]}
    )
    mask = 0
    active = q_active = weighted = 0
    state_counts = [0] * V
    branch_counts = [0] * P
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
                ("owner digit escaped cell zero", alpha, beta, left, right))
        chamber = R.chamber(left, right)
        require(chamber in ("left", "right"),
                ("owner digit entered middle chamber", alpha, beta, left, right))
        state = SQ.state_index(source_u, left)
        state_counts[state] += 1
        if not jump:
            continue
        q_active += 1
        u_values = tuple(profile_value(profile, left) for profile in source_u)
        v_values = tuple(profile_value(profile, left) for profile in source_v)
        digit_values = tuple(
            tuple(profile_value(branches[root][digit], left) for digit in range(P))
            for root in range(P)
        )
        require(all(sum(digit_values[root]) == u_values[root] for root in range(P)),
                ("local digit marginal", left))
        require(all(u_values[root] * v_values[root] == 0 for root in range(P)),
                ("pointwise same-root source", left))
        weighted += 1
        for digit in range(P):
            if any(digit_values[root][digit] for root in range(P)):
                branch_counts[digit] += 1

        for row_index, tau in enumerate(tau_values):
            selected_mask = mask if literal_tau is not None else (
                mask & R.guard_mask(chamber, tau, danger)
            )
            selected = tuple(root for root in range(P)
                             if (selected_mask >> root) & 1)
            right_value = sum(v_values[root] for root in selected)
            for digit in range(P):
                left_value = sum(digit_values[root][digit] for root in selected)
                diagonal = sum(digit_values[root][digit] * v_values[root]
                               for root in selected)
                require(diagonal == 0,
                        ("selected same-root", alpha, beta, tau, digit, left))
                bank[row_index][state][digit] = (
                    bank[row_index][state][digit]
                    + left_value * right_value * jump
                ) % PRIME
    mask ^= events.get(positions[-1], 0)
    require(mask == 0, ("endpoint mask closure", alpha, beta, literal_tau))
    rows = tuple(tuple(tuple(row) for row in state_rows) for state_rows in bank)
    counts = (
        interval_count, mapped_count, active, q_active, weighted,
        tuple(state_counts), tuple(branch_counts),
    )
    return rows, counts


def worker(alpha: int):
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    gamma = []
    scalar_counts = [0] * 5
    state_counts = [0] * V
    branch_counts = [0] * P
    for beta in range(P):
        rows, counts = integrate_pair(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        gamma.extend(
            tuple(tuple(phase * value % PRIME for value in digit_row)
                  for digit_row in state_rows)
            for state_rows in rows
        )
        scalar_counts = [left + right for left, right in
                         zip(scalar_counts, counts[:5])]
        state_counts = [left + right for left, right in
                        zip(state_counts, counts[5])]
        branch_counts = [left + right for left, right in
                         zip(branch_counts, counts[6])]
    return alpha, tuple(gamma), (
        tuple(scalar_counts), tuple(state_counts), tuple(branch_counts)
    )


def invert_gamma(gamma, zeta: int):
    tensor = [[[0] * P for _digit in range(P)] for _state in range(V)]
    normalizer = pow(P**3, -1, PRIME)
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                phases = tuple(pow(zeta, -(alpha + tau * relation) % P, PRIME)
                               for relation in range(P))
                for state in range(V):
                    for digit in range(P):
                        value = row[state][digit]
                        for relation, phase in enumerate(phases):
                            tensor[state][digit][relation] = (
                                tensor[state][digit][relation] + value * phase
                            ) % PRIME
    require(index == P**3, ("inverse gamma length", index))
    return tuple(
        tuple(tuple(value * normalizer % PRIME for value in relation_row)
              for relation_row in digit_rows)
        for digit_rows in tensor
    )


def digit_marginal(tensor):
    return tuple(tuple(
        sum(tensor[state][digit][relation] for digit in range(P)) % PRIME
        for relation in range(P)
    ) for state in range(V))


def flat_lift(parent):
    inv_p = pow(P, -1, PRIME)
    return tuple(tuple(tuple(
        parent[state][relation] * inv_p % PRIME for relation in range(P)
    ) for _digit in range(P)) for state in range(V))


def centre(tensor, axis: int):
    dimensions = (V, P, P)
    inverse = pow(dimensions[axis], -1, PRIME)
    answer = [[[tensor[state][digit][relation] for relation in range(P)]
               for digit in range(P)] for state in range(V)]
    for state in range(V):
        for digit in range(P):
            for relation in range(P):
                if axis == 0:
                    mean = sum(tensor[s][digit][relation] for s in range(V))
                elif axis == 1:
                    mean = sum(tensor[state][d][relation] for d in range(P))
                else:
                    mean = sum(tensor[state][digit][t] for t in range(P))
                answer[state][digit][relation] = (
                    tensor[state][digit][relation] - mean * inverse
                ) % PRIME
    return tuple(tuple(tuple(row) for row in plane) for plane in answer)


def three_way(tensor):
    return centre(centre(centre(tensor, 0), 1), 2)


def flattening_ranks(tensor):
    state_rows = tuple(tuple(
        tensor[state][digit][relation]
        for digit in range(P) for relation in range(P)
    ) for state in range(V))
    digit_rows = tuple(tuple(
        tensor[state][digit][relation]
        for state in range(V) for relation in range(P)
    ) for digit in range(P))
    relation_rows = tuple(tuple(
        tensor[state][digit][relation]
        for state in range(V) for digit in range(P)
    ) for relation in range(P))
    return tuple(R.rank_mod(rows) for rows in
                 (state_rows, digit_rows, relation_rows))


def fourier(tensor, zeta: int):
    roots = tuple(tuple(pow(zeta, -frequency * value % P, PRIME)
                        for value in range(P)) for frequency in range(P))
    return tuple(tuple(tuple(sum(
        tensor[state][digit][relation]
        * WALSH_SIGNS[character][state]
        * roots[digit_frequency][digit]
        * roots[relation_frequency][relation]
        for state in range(V) for digit in range(P) for relation in range(P)
    ) % PRIME for relation_frequency in range(P))
        for digit_frequency in range(P)) for character in range(V))


def support_census(spectrum):
    bins = [0] * 8
    for character in range(V):
        for digit_frequency in range(P):
            for relation_frequency in range(P):
                if spectrum[character][digit_frequency][relation_frequency]:
                    mask = (
                        ((character != 0) << 2)
                        | ((digit_frequency != 0) << 1)
                        | (relation_frequency != 0)
                    )
                    bins[mask] += 1
    return (sum(bins),) + tuple(bins)


def fixed_relation(tensor, relation: int):
    matrix = tuple(tuple(tensor[state][digit][relation] for digit in range(P))
                   for state in range(V))
    interaction = three_way(tensor)
    slice_matrix = tuple(tuple(interaction[state][digit][relation]
                               for digit in range(P)) for state in range(V))
    return (
        R.rank_mod(matrix),
        R.rank_mod(slice_matrix),
        tuple(sum(value != 0 for value in row) for row in matrix),
        digest_json(matrix),
    )


def main() -> None:
    R.split_field_certificate()
    radix = radix_certificate()
    (
        _branches, branch_boundaries, _source_u, _source_v,
        source_boundaries, source_digest, source_record, source_sha,
        reflection_record,
    ) = owner_digit_profiles()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "worker order")
    gamma = tuple(row for chunk in chunks for row in chunk[1])
    require(len(gamma) == P**3, ("gamma rows", len(gamma)))
    scalar_counts = tuple(sum(chunk[2][0][index] for chunk in chunks)
                          for index in range(5))
    state_counts = tuple(sum(chunk[2][1][state] for chunk in chunks)
                         for state in range(V))
    branch_counts = tuple(sum(chunk[2][2][digit] for chunk in chunks)
                          for digit in range(P))
    require(scalar_counts == EXPECTED_WORK_COUNTS,
            ("work counts", scalar_counts))
    require(state_counts == EXPECTED_STATE_COUNTS,
            ("state counts", state_counts))
    require(branch_counts == EXPECTED_BRANCH_COUNTS,
            ("branch counts", branch_counts))

    gamma_parent = tuple(tuple(sum(row[state]) % PRIME for state in range(V))
                         for row in gamma)
    require(digest_json(gamma_parent) == EXPECTED_PARENT_GAMMA,
            "gamma marginal misses independent square")
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    tensor = invert_gamma(gamma, zeta)
    parent = digit_marginal(tensor)
    require(digest_json(parent) == EXPECTED_PARENT_TABLE,
            "table marginal misses independent square")
    flat = flat_lift(parent)
    require(digit_marginal(flat) == parent, "flat lift marginal")

    word, endpoint_grid, _harmonic, danger = endpoint_context()
    literal_boundaries = (
        set(branch_boundaries) | R.fixed_boundaries(
            source_boundaries, R.SRC.T_DEN
        )
    )
    literal_records = []
    for alpha, beta, tau in CONTROL_TRIPLES:
        restoration = R.literal_guard_restoration(
            word, endpoint_grid, alpha, beta, tau, literal_boundaries, danger
        )
        direct, counts = integrate_pair(alpha, beta, tau)
        phase = pow(zeta, beta, PRIME)
        expected = tuple(tuple(phase * value % PRIME for value in digit_row)
                         for digit_row in direct[0])
        index = (alpha * P + beta) * P + tau
        require(expected == gamma[index],
                ("literal guard row", alpha, beta, tau))
        literal_records.append(((alpha, beta, tau), restoration, counts))

    tensors = (tensor, flat)
    interactions = tuple(three_way(value) for value in tensors)
    spectra = tuple(fourier(value, zeta) for value in tensors)
    interaction_spectra = tuple(fourier(value, zeta) for value in interactions)
    ranks = tuple(flattening_ranks(value) for value in tensors)
    interaction_ranks = tuple(flattening_ranks(value) for value in interactions)
    censuses = tuple(support_census(value) for value in spectra)
    interaction_censuses = tuple(support_census(value)
                                 for value in interaction_spectra)
    fixed = tuple(fixed_relation(value, 6) for value in tensors)
    rank_record = (ranks, interaction_ranks)
    census_record = (censuses, interaction_censuses)
    require(rank_record == EXPECTED_RANKS, ("rank record", rank_record))
    require(census_record == EXPECTED_CENSUSES,
            ("spectral census", census_record))
    require(fixed == EXPECTED_FIXED, ("fixed t=6", fixed))

    digests = (
        digest_json(gamma), digest_json(tensor), digest_json(flat),
        tuple(digest_json(value) for value in spectra),
        tuple(digest_json(value) for value in interactions),
        tuple(digest_json(value) for value in interaction_spectra),
    )
    require(digests == EXPECTED_DIGESTS, ("candidate comparison", digests))
    record = (
        SQUARE_SHA256, SQUARE_SEMANTIC, radix, source_digest, source_record,
        source_sha, reflection_record, scalar_counts, state_counts,
        branch_counts, tuple(literal_records), digests, rank_record,
        census_record, fixed,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== independent hostile audit: current-leg r_owner x Boolean square x relation ==")
    print(f"parent=(square_sha256={SQUARE_SHA256},square_semantic={SQUARE_SEMANTIC})")
    print(f"radix=(depth={radix[0]},upper_fold={radix[1]},digit_census={radix[2]}): PASS")
    print("typing=r_owner=a mod13 in T^4 X_(u,a)=(w_u+r_owner)/13 mod1;current leg only")
    print(f"source_decomposition=(upper_fold_pieces={source_record[2]},digit_window_pieces={source_record[3]},joint_boundaries={source_record[5]},branch_masses={source_record[6]},sha256={source_sha})")
    print(f"reflection=(u,r,state)->(12-u,12-r,state_XOR2);interval_record={reflection_record}: PASS")
    print(f"field=(prime={R.JOINT_PRIME},root={R.JOINT_ROOT},zeta13={zeta})")
    print(f"work_counts={scalar_counts};state_counts={state_counts};branch_counts={branch_counts}")
    print(f"literal_controls={tuple(record[0] for record in literal_records)}: PASS")
    print("marginals=(digit_sum->independently_audited_square gamma and table): PASS;same_root=0")
    print(f"flattening_ranks_(state,digit,relation)=(actual,flat)={ranks}")
    print(f"three_way_ANOVA_flattening_ranks={interaction_ranks}")
    print("support_census_order=(total,empty,relation,digit,digit_relation,state,state_relation,state_digit,triple)")
    print(f"spectral_censuses=(actual,flat)={censuses}")
    print(f"three_way_ANOVA_spectral_censuses={interaction_censuses}")
    print(f"fixed_(1,0,6)_records_(rank,threeway_slice_rank,row_support,digest)={fixed}")
    print(f"digests=(gamma,tensor,flat,spectra,ANOVA,ANOVA_spectra)={digests}")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT scoped nonredundant current-leg inverse-owner-digit sidecar")
    print("scope=not source-time P169 digit,not root difference,not deep label,not full inverse ancestry,not exact address,not U_clock chronology,not physical current,no row exclusion,no LRC(14)")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
