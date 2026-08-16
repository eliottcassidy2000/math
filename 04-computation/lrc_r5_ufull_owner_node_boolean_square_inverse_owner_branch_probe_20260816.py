#!/usr/bin/env python3
"""Retain the canonical last inverse-owner branch above the Boolean square.

For ``d=13^5`` and ``a=r+13c``, the THM-2471 current sheet satisfies

    X_(u,a) = (w_u+a)/d,
    r_owner = a mod 13.

The branch-resolved source profile is reconstructed without enumerating
``13^5`` sheets.  Split the Perron fold as ``d=13*13^4``: first fold the
current packet by ``13^4``, then take the ``r_owner`` window, then the
collision-root ``u`` window.  The resulting profiles ``U_(u,r)`` sum
pointwise to the audited ``U_u`` profiles.

The script retains ``r_owner`` together with the four owner-visible source
states before multiplying and integrating the actual U_full endpoint factors,
then inverts to the thirteen refined endpoint relation residues.  Summing the
branch coordinate must recover the independently audited Boolean-square
table.  A branch-flat lift with the same marginal is the primary hostile.

This is the canonical last inverse-owner digit only.  It is not the deep
label, a complete inverse-ancestry address, THM-2334's grouped exact address,
a U_clock chronology, or a physical current.
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
    / "04-computation/lrc_r5_ufull_owner_node_boolean_square_refiner_probe_20260816.py"
)
SQUARE_SHA256 = "c7bbd2d82ed067914f39253ac52ce32f9e179a859fb55a421bd2b390f707d881"

EXPECTED_SOURCE = "9e2a61ea8ec00830aa93ec4daa42574a9179f5aa80a09e3388dd9af0964f4cb1"
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
    (
        4,
        3,
        (13, 13, 13, 13),
        "7652ac0fbca27c3a41181d04e9bea962f4159b04ccc43befe3d2783af89be4a4",
    ),
    (
        1,
        0,
        (13, 13, 13, 13),
        "be782cc3ee60d8d4bc008f21a9f308731a29e180bc83d7041ff87c47eb532729",
    ),
)
EXPECTED_SEMANTIC_SHA256 = (
    "f43e9fb0937ae93c54c2ac0aeb1257dfd8741f9cef32ef4bf2a5f3dd2cf06388"
)

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
            ("Boolean-square source drift", observed, SQUARE_SHA256))
    spec = importlib.util.spec_from_file_location("inverse_owner_square_parent", SQUARE_PATH)
    require(spec is not None and spec.loader is not None, "inverse-owner parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_square()
C = B.C
PRIME = C.JOINT_PRIME


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


@lru_cache(maxsize=1)
def owner_branch_context():
    ctx = C.context()
    stage = C.SM
    source_grid = C.S.GRID
    require(stage.DCOLL == P**5, ("collision depth", stage.DCOLL))
    e_intervals = stage.build_set(stage.PAT_E, stage.ZELL)
    q_intervals = stage.build_set(stage.PAT_QB, stage.ZELL)
    f_pieces = C.S.B.build_f_pieces(e_intervals, q_intervals)

    upper_fold = stage.weighted_fold(f_pieces, stage.DCOLL // P, source_grid)
    branch_at_owner = tuple(
        stage.extract_window(upper_fold[0], upper_fold[1], branch, P, source_grid)
        for branch in range(P)
    )
    branch_at_root_raw = tuple(
        tuple(
            stage.extract_window(
                branch_at_owner[branch][0], branch_at_owner[branch][1], root, P, source_grid
            )
            for branch in range(P)
        )
        for root in range(P)
    )
    branch_at_root = tuple(
        tuple(
            C.scale_profile(profile, C.JOINT_COORDINATE, source_grid)
            for profile in root_profiles
        )
        for root_profiles in branch_at_root_raw
    )
    boundaries = tuple(sorted(
        {0, C.JOINT_COORDINATE}
        | {
            position
            for root_profiles in branch_at_root
            for profile in root_profiles
            for position in profile[0]
        }
        | set(ctx["source_boundaries"])
    ))

    require(
        tuple(C.JOINT_COORDINATE - point for point in reversed(boundaries))
        == boundaries,
        "inverse-owner boundary reflection",
    )

    for left, right in zip(boundaries, boundaries[1:]):
        for root in range(P):
            branch_sum = sum(
                profile_value(branch_at_root[root][branch], left)
                for branch in range(P)
            )
            total = profile_value(ctx["source_u"][root], left)
            require(branch_sum == total,
                    ("owner-branch pointwise partition", left, right, root, branch_sum, total))
            for branch in range(P):
                reflected = profile_value(
                    branch_at_root[P - 1 - root][P - 1 - branch],
                    C.JOINT_COORDINATE - right,
                )
                require(
                    profile_value(branch_at_root[root][branch], left) == reflected,
                    (
                        "inverse-owner profile reflection",
                        left,
                        right,
                        root,
                        branch,
                    ),
                )

        try:
            state = B.state_of_segment(left, right)
        except RuntimeError:
            continue
        reflected_state = B.state_of_segment(
            C.JOINT_COORDINATE - right,
            C.JOINT_COORDINATE - left,
        )
        require(reflected_state == state ^ 2,
                ("inverse-owner state reflection", left, right, state, reflected_state))

    branch_masses = tuple(
        sum(profile_mass(branch_at_root_raw[root][branch], source_grid) for root in range(P))
        for branch in range(P)
    )
    require(branch_masses == tuple(reversed(branch_masses)),
            ("inverse-owner mass reflection", branch_masses))
    source_record = (
        stage.DCOLL,
        stage.DCOLL // P,
        len(upper_fold[0]),
        tuple(len(profile[0]) for profile in branch_at_owner),
        tuple(tuple(len(profile[0]) for profile in row) for row in branch_at_root_raw),
        len(boundaries),
        branch_masses,
        digest_json(branch_at_root_raw),
        digest_json(branch_at_root),
    )
    return branch_at_root, boundaries, source_record


def integrate_branches(alpha: int, beta: int, literal_tau: int | None = None):
    ctx = C.context()
    branch_profiles, branch_boundaries, _source_record = owner_branch_context()
    events, interval_count, mapped = C.endpoint_events(alpha, beta, literal_tau)
    for boundary in branch_boundaries:
        events.setdefault(boundary, 0)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    coupled = [
        [[0 for _branch in range(P)] for _state in range(V)]
        for _tau in tau_values
    ]
    diagonal = [
        [[0 for _branch in range(P)] for _state in range(V)]
        for _tau in tau_values
    ]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    state_segments = [0] * V
    branch_segments = [0] * P
    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        require(C.cell_of_segment(left, right) == 0,
                ("inverse-owner branch escaped cell zero", alpha, beta, left, right))
        chamber = C.chamber_of_segment(left, right)
        require(chamber in ("left", "right"),
                ("inverse-owner active middle chamber", alpha, beta, left, right))
        state = B.state_of_segment(left, right)
        active_segments += 1
        state_segments[state] += 1
        jump = C.q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        u_values = tuple(profile_value(profile, left) for profile in ctx["source_u"])
        v_values = tuple(profile_value(profile, left) for profile in ctx["source_v"])
        branch_values = tuple(
            tuple(profile_value(branch_profiles[root][branch], left) for branch in range(P))
            for root in range(P)
        )
        require(all(sum(branch_values[root]) == u_values[root] for root in range(P)),
                ("inverse-owner local partition", left, right))
        if not any(u_values) or not any(v_values):
            continue
        weighted_segments += 1
        for branch in range(P):
            if any(branch_values[root][branch] for root in range(P)):
                branch_segments[branch] += 1
        for row_index, tau in enumerate(tau_values):
            if literal_tau is None:
                selected = tuple(
                    sheet
                    for sheet in range(P)
                    if (mask >> sheet) & 1 and C.T.safe(chamber, sheet + tau)
                )
            else:
                selected = tuple(sheet for sheet in range(P) if (mask >> sheet) & 1)
            right_value = sum(v_values[root] for root in selected)
            for branch in range(P):
                left_value = sum(branch_values[root][branch] for root in selected)
                diagonal_value = sum(
                    branch_values[root][branch] * v_values[root] for root in selected
                )
                require(diagonal_value == 0,
                        ("inverse-owner same-root", alpha, beta, tau, branch, left, right))
                coupled[row_index][state][branch] = (
                    coupled[row_index][state][branch] + left_value * right_value * jump
                ) % PRIME
                diagonal[row_index][state][branch] = (
                    diagonal[row_index][state][branch] + diagonal_value * jump
                ) % PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("inverse-owner endpoint mask", alpha, beta, literal_tau, mask))
    require(all(value == 0 for row in diagonal for state_row in row for value in state_row),
            ("inverse-owner diagonal bank", alpha, beta, literal_tau))
    counts = (
        interval_count,
        mapped,
        active_segments,
        q_active_segments,
        weighted_segments,
        tuple(state_segments),
        tuple(branch_segments),
    )
    freeze = lambda bank: tuple(
        tuple(tuple(row) for row in state_rows) for state_rows in bank
    )
    return freeze(coupled), freeze(diagonal), counts


def worker(alpha: int):
    zeta = C.context()["zeta"]
    coupled_rows = []
    diagonal_rows = []
    scalar_counts = [0] * 5
    state_counts = [0] * V
    branch_counts = [0] * P
    for beta in range(P):
        coupled, diagonal, record = integrate_branches(alpha, beta)
        phase = pow(zeta, beta, PRIME)

        def twist(bank):
            return tuple(
                tuple(
                    tuple(phase * value % PRIME for value in branch_row)
                    for branch_row in state_rows
                )
                for state_rows in bank
            )

        coupled_rows.append(twist(coupled))
        diagonal_rows.append(diagonal)
        scalar_counts = [left + right for left, right in zip(scalar_counts, record[:5])]
        state_counts = [left + right for left, right in zip(state_counts, record[5])]
        branch_counts = [left + right for left, right in zip(branch_counts, record[6])]
    return (
        alpha,
        tuple(coupled_rows),
        tuple(diagonal_rows),
        (tuple(scalar_counts), tuple(state_counts), tuple(branch_counts)),
    )


def inverse_tensor(gamma_rows, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    tensor = [[[0 for _relation in range(P)] for _branch in range(P)] for _state in range(V)]
    phases = tuple(
        tuple(pow(zeta, -(alpha + tau * relation) % P, PRIME) for relation in range(P))
        for alpha in range(P) for tau in range(P)
    )
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma_rows[index]
                phase_row = phases[alpha * P + tau]
                index += 1
                for state in range(V):
                    for branch in range(P):
                        value = row[state][branch]
                        if not value:
                            continue
                        for relation, phase in enumerate(phase_row):
                            tensor[state][branch][relation] = (
                                tensor[state][branch][relation] + value * phase
                            ) % PRIME
    require(index == P**3, ("inverse-owner gamma size", index))
    return tuple(
        tuple(
            tuple(value * normalizer % PRIME for value in relation_row)
            for relation_row in branch_rows
        )
        for branch_rows in tensor
    )


def branch_marginal(tensor):
    return tuple(
        tuple(sum(tensor[state][branch][relation] for branch in range(P)) % PRIME
              for relation in range(P))
        for state in range(V)
    )


def flat_branch_lift(parent):
    inverse_p = pow(P, -1, PRIME)
    return tuple(
        tuple(
            tuple(parent[state][relation] * inverse_p % PRIME for relation in range(P))
            for _branch in range(P)
        )
        for state in range(V)
    )


def centre_axis(tensor, axis: int):
    dimensions = (V, P, P)
    inverse = pow(dimensions[axis], -1, PRIME)
    mutable = [[[tensor[state][branch][relation] for relation in range(P)]
                for branch in range(P)] for state in range(V)]
    if axis == 0:
        for branch in range(P):
            for relation in range(P):
                mean = sum(mutable[state][branch][relation] for state in range(V)) * inverse % PRIME
                for state in range(V):
                    mutable[state][branch][relation] = (mutable[state][branch][relation] - mean) % PRIME
    elif axis == 1:
        for state in range(V):
            for relation in range(P):
                mean = sum(mutable[state][branch][relation] for branch in range(P)) * inverse % PRIME
                for branch in range(P):
                    mutable[state][branch][relation] = (mutable[state][branch][relation] - mean) % PRIME
    else:
        for state in range(V):
            for branch in range(P):
                mean = sum(mutable[state][branch][relation] for relation in range(P)) * inverse % PRIME
                for relation in range(P):
                    mutable[state][branch][relation] = (mutable[state][branch][relation] - mean) % PRIME
    return tuple(tuple(tuple(row) for row in plane) for plane in mutable)


def three_way_interaction(tensor):
    return centre_axis(centre_axis(centre_axis(tensor, 0), 1), 2)


def flattening_ranks(tensor):
    state_flat = tuple(
        tuple(tensor[state][branch][relation] for branch in range(P) for relation in range(P))
        for state in range(V)
    )
    branch_flat = tuple(
        tuple(tensor[state][branch][relation] for state in range(V) for relation in range(P))
        for branch in range(P)
    )
    relation_flat = tuple(
        tuple(tensor[state][branch][relation] for state in range(V) for branch in range(P))
        for relation in range(P)
    )
    return C.rank_mod(state_flat), C.rank_mod(branch_flat), C.rank_mod(relation_flat)


def fourier_tensor(tensor, zeta: int):
    roots = tuple(
        tuple(pow(zeta, -frequency * value % P, PRIME) for value in range(P))
        for frequency in range(P)
    )
    walsh_stage = [[
        [
            sum(tensor[state][branch][relation] * B.WALSH_SIGNS[character][state]
                for state in range(V)) % PRIME
            for relation in range(P)
        ]
        for branch in range(P)
    ] for character in range(V)]
    branch_stage = [[
        [
            sum(walsh_stage[character][branch][relation] * roots[frequency][branch]
                for branch in range(P)) % PRIME
            for relation in range(P)
        ]
        for frequency in range(P)
    ] for character in range(V)]
    return tuple(
        tuple(
            tuple(
                sum(branch_stage[character][branch_frequency][relation]
                    * roots[relation_frequency][relation]
                    for relation in range(P)) % PRIME
                for relation_frequency in range(P)
            )
            for branch_frequency in range(P)
        )
        for character in range(V)
    )


def support_census(spectrum):
    bins = [0] * 8
    for character in range(V):
        for branch_frequency in range(P):
            for relation_frequency in range(P):
                if spectrum[character][branch_frequency][relation_frequency]:
                    mask = ((character != 0) << 2) | ((branch_frequency != 0) << 1) | (relation_frequency != 0)
                    bins[mask] += 1
    return (sum(bins),) + tuple(bins)


def fixed_relation_record(tensor, relation: int):
    matrix = tuple(
        tuple(tensor[state][branch][relation] for branch in range(P))
        for state in range(V)
    )
    interaction = tuple(
        tuple(three_way_interaction(tensor)[state][branch][relation] for branch in range(P))
        for state in range(V)
    )
    return (
        C.rank_mod(matrix),
        C.rank_mod(interaction),
        tuple(sum(value != 0 for value in row) for row in matrix),
        digest_json(matrix),
    )


def main() -> None:
    branch_profiles, branch_boundaries, source_record = owner_branch_context()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "inverse-owner worker order")

    coupled_gamma = tuple(
        row
        for chunk in chunks
        for beta_rows in chunk[1]
        for row in beta_rows
    )
    diagonal_gamma = tuple(
        row
        for chunk in chunks
        for beta_rows in chunk[2]
        for row in beta_rows
    )
    require(len(coupled_gamma) == len(diagonal_gamma) == P**3,
            "inverse-owner bank size")
    require(all(value == 0 for row in diagonal_gamma for state_row in row for value in state_row),
            "inverse-owner same-root gamma")

    gamma_marginal = tuple(
        tuple(sum(row[state]) % PRIME for state in range(V))
        for row in coupled_gamma
    )
    require(digest_json(gamma_marginal) == B.EXPECTED_DIGESTS[0],
            "inverse-owner gamma marginal misses audited square")

    zeta = C.context()["zeta"]
    tensor = inverse_tensor(coupled_gamma, zeta)
    parent = branch_marginal(tensor)
    require(digest_json(parent) == B.EXPECTED_DIGESTS[2],
            "inverse-owner tensor marginal misses audited square")
    flat = flat_branch_lift(parent)
    require(branch_marginal(flat) == parent, "flat branch lift misses parent")
    tensors = (tensor, flat)
    interactions = tuple(three_way_interaction(value) for value in tensors)
    spectra = tuple(fourier_tensor(value, zeta) for value in tensors)
    interaction_spectra = tuple(fourier_tensor(value, zeta) for value in interactions)
    ranks = tuple(flattening_ranks(value) for value in tensors)
    interaction_ranks = tuple(flattening_ranks(value) for value in interactions)
    censuses = tuple(support_census(value) for value in spectra)
    interaction_censuses = tuple(support_census(value) for value in interaction_spectra)
    relation = 6
    fixed_records = tuple(fixed_relation_record(value, relation) for value in tensors)

    scalar_counts = tuple(sum(chunk[3][0][index] for chunk in chunks) for index in range(5))
    state_counts = tuple(sum(chunk[3][1][state] for chunk in chunks) for state in range(V))
    branch_counts = tuple(sum(chunk[3][2][branch] for chunk in chunks) for branch in range(P))
    require(all(count > 0 for count in branch_counts),
            ("empty inverse-owner branch", branch_counts))

    direct_controls = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
    direct_record = []
    for alpha, beta, tau in direct_controls:
        direct, direct_diagonal, counts = integrate_branches(alpha, beta, tau)
        phase = pow(zeta, beta, PRIME)
        direct_row = tuple(
            tuple(phase * value % PRIME for value in branch_row)
            for branch_row in direct[0]
        )
        index = (alpha * P + beta) * P + tau
        require(direct_row == coupled_gamma[index],
                ("inverse-owner literal guard", alpha, beta, tau))
        require(all(value == 0 for state_row in direct_diagonal[0] for value in state_row),
                ("inverse-owner literal diagonal", alpha, beta, tau))
        direct_record.append(((alpha, beta, tau), counts))

    digests = (
        digest_json(coupled_gamma),
        digest_json(tensor),
        digest_json(flat),
        tuple(digest_json(value) for value in spectra),
        tuple(digest_json(value) for value in interactions),
        tuple(digest_json(value) for value in interaction_spectra),
    )
    rank_record = (ranks, interaction_ranks)
    census_record = (censuses, interaction_censuses)
    source_digest = digest_json((source_record, branch_boundaries))
    if EXPECTED_SOURCE != "TO_BE_PINNED":
        require(source_digest == EXPECTED_SOURCE,
                ("inverse-owner source decomposition", source_digest, EXPECTED_SOURCE))
    if EXPECTED_DIGESTS != "TO_BE_PINNED":
        require(digests == EXPECTED_DIGESTS, ("inverse-owner digests", digests))
    if EXPECTED_RANKS != "TO_BE_PINNED":
        require(rank_record == EXPECTED_RANKS, ("inverse-owner ranks", rank_record))
    if EXPECTED_CENSUSES != "TO_BE_PINNED":
        require(census_record == EXPECTED_CENSUSES,
                ("inverse-owner censuses", census_record))
    if EXPECTED_FIXED != "TO_BE_PINNED":
        require(fixed_records == EXPECTED_FIXED,
                ("inverse-owner fixed relation", fixed_records))

    record = (
        SQUARE_SHA256,
        B.EXPECTED_SEMANTIC_SHA256,
        C.JOINT_PRIME,
        C.JOINT_ROOT,
        zeta,
        source_record,
        branch_boundaries,
        scalar_counts,
        state_counts,
        branch_counts,
        tuple(direct_record),
        digests,
        rank_record,
        census_record,
        relation,
        fixed_records,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("inverse-owner semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 Boolean-square x inverse-owner-branch x relation probe ==")
    print(f"parent=(sha256={SQUARE_SHA256},semantic={B.EXPECTED_SEMANTIC_SHA256})")
    print("coordinate=r_owner=a mod 13 on X_(u,a);construction=13^4 fold then r_owner window then collision-root window")
    print(f"source_decomposition=(upper_fold_pieces={source_record[2]},owner_window_pieces={source_record[3]},joint_boundaries={source_record[5]},branch_masses={source_record[6]},digest={source_digest})")
    print("reflection=(base y->1-y,collision_root u->12-u,owner_branch r->12-r,state->state_XOR_2): PASS")
    print(f"field=(prime={C.JOINT_PRIME},root={C.JOINT_ROOT},zeta13={zeta})")
    print(f"work_counts={scalar_counts};state_counts={state_counts};branch_counts={branch_counts};literal_controls={direct_controls}: PASS")
    print("marginals=(branch_sum->audited_Boolean_square,flat_branch_same_parent): PASS;same_root=0")
    print(f"flattening_ranks_(state,owner_branch,relation)=(actual,flat)={ranks}")
    print(f"three_way_ANOVA_flattening_ranks={interaction_ranks}")
    print("support_census_order=(total,empty,relation,branch,branch_relation,state,state_relation,state_branch,triple)")
    print(f"spectral_censuses=(actual,flat)={censuses}")
    print(f"three_way_ANOVA_spectral_censuses={interaction_censuses}")
    print(f"fixed_relation=(1,0,{relation});records_(rank,three_way_slice_rank,row_support,digest)={fixed_records}")
    print(f"digests=(gamma,tensor,flat,spectra,ANOVA,ANOVA_spectra)={digests}")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT canonical inverse-owner-branch crossing candidate on one owner base")
    print("scope=r_owner is last inverse digit only;not deep label,not full inverse ancestry,not exact C(a;X,m),not chronology,not physical current,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
