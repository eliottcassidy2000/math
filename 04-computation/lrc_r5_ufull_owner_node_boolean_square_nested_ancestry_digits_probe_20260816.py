#!/usr/bin/env python3
"""Cross the source-time b sheet with the current-leg inverse-owner digit.

On THM-2471's exact Boolean stalk, write

    a = r_owner + 13 h,       0 <= a < 13^5,
    n = c + 13 b_source,      0 <= n < 13^2.

Then

    X_(u,a) = (w_u+a)/13^5,
    Z_(u,a,n) = (X_(u,a)+n)/13^2.

Thus r_owner is an outer-current digit and b_source is an inner prescribed-
ancestry digit.  They coexist on (y,u,a,n,e') without an arrival/source
intertwiner.  The Q(X) cut remains between the inner source packet and the
outer Perron fold.

The script builds all 13x13 nonnegative source profiles and checks both exact
one-digit marginals pointwise.  The 144 double-nontrivial character profiles
have exact rank 17, so the response gate integrates a certified 17-character
basis rather than all 144 channels.  These vanish under either digit marginal
and under the additive lift with the same two marginals.  Their
state/relation-centred rank is therefore a conditional rank gate, not an
inference from the two separate rank-six/432 signatures.

The same endpoint sweep independently retains the six lawful pointed tails
``(state,u)``.  State by state, it compares the complete joint-character
response span with the pointed-tail relation span and solves explicit change-
of-coordinate equations.  This relation-output comparison can decide whether
the two nested digits add a carrier direction at that stage; it does not
identify the digits or retain the source-root difference coordinate.

This is one finite r=5 Boolean stalk.  It is not a temporal current, exact
THM-2334 address, U_clock chronology, row exclusion, or LRC(14) result.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / (
    "04-computation/lrc_r5_ufull_owner_node_boolean_square_"
    "branch_sheet_sidecar_probe_20260816.py"
)
OWNER_PATH = ROOT / (
    "04-computation/lrc_r5_ufull_owner_node_boolean_square_"
    "inverse_owner_branch_probe_20260816.py"
)
SOURCE_SHA256 = "592aa0bce31f2da5d5e2ddff7f3ffe6f1398f3a07b5ce927e0d97c9fe309ae3b"
OWNER_SHA256 = "ae1cf021ea23f325eeded42ff1dea8df903d837b1d7ab551289d62f7ab7a0348"

P = 13
V = 4
ALL_DOUBLE_PAIRS = tuple(
    (source_frequency, owner_frequency)
    for source_frequency in range(1, P)
    for owner_frequency in range(1, P)
)
CHARACTER_PAIRS = tuple(ALL_DOUBLE_PAIRS[index] for index in (
    0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 24, 25, 36, 37, 48, 49,
))
CONTROL_TRIPLES = ((0, 0, 0), (1, 0, 6), (6, 6, 12), (12, 12, 0))
POINTED_STATES = ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12))
POINT_INDEX = {pair: index for index, pair in enumerate(POINTED_STATES)}
STATE_POINTS = ((0,), (1, 2), (5,), (3, 4))
EXPECTED_SEMANTIC_SHA256 = "6e5605f58b7a94ea5ea4e8f62cfa7ee135b0d52512225f4aaee248ad6e21a9ae"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def load_module(path: Path, name: str, expected: str):
    observed = lf_sha256(path)
    require(observed == expected, (name, observed, expected))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


S = load_module(SOURCE_PATH, "nested_ancestry_source_parent", SOURCE_SHA256)
R = load_module(OWNER_PATH, "nested_ancestry_owner_parent", OWNER_SHA256)
B = S.B
PRIME = B.PRIME

require(R.P == S.P == P, "prime-cardinality mismatch")
require(R.C.JOINT_COORDINATE == B.JOINT_COORDINATE, "joint-base mismatch")
require(R.PRIME == PRIME, "split-field mismatch")


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def scale_profile(profile, scale: int):
    starts, values = profile
    return tuple(point * scale for point in starts), tuple(values)


def compress_profile(starts, values):
    out_starts = []
    out_values = []
    for start, value in zip(starts, values):
        value %= PRIME
        if out_values and out_values[-1] == value:
            continue
        out_starts.append(start)
        out_values.append(value)
    require(out_starts and out_starts[0] == 0, "compressed profile origin")
    return tuple(out_starts), tuple(out_values)


def row_basis_indices(matrix):
    basis = {}
    indices = []
    for index, source_row in enumerate(matrix):
        row = [value % PRIME for value in source_row]
        for pivot in sorted(basis):
            factor = row[pivot]
            if factor:
                pivot_row = basis[pivot]
                row = [
                    (value - factor * pivot_value) % PRIME
                    for value, pivot_value in zip(row, pivot_row)
                ]
        pivot = next((column for column, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        inverse = pow(row[pivot], -1, PRIME)
        row = [value * inverse % PRIME for value in row]
        basis[pivot] = tuple(row)
        indices.append(index)
    return tuple(indices)


@lru_cache(maxsize=1)
def joint_context():
    ctx = B.context()
    stage = ctx["M"]
    source_grid = stage.T_DEN
    require(stage.DCOLL == P**5 and stage.RPKT == P**2,
            (stage.DCOLL, stage.RPKT))
    e_intervals = stage.build_set(stage.PAT_E, stage.ZELL)
    q_intervals = stage.build_set(stage.PAT_QB, stage.ZELL)

    # Inner source ancestry: n=c+13*b_source in P_169 e, then literal Q.
    p2_by_source = S.branch_numerator_profiles(stage, e_intervals)
    f_by_source = tuple(
        B.restrict_profile(starts, values, q_intervals, source_grid)
        for starts, values in p2_by_source
    )

    # Outer current ancestry: a=r_owner+13*h in P_(13^5) f_b.
    raw = []
    upper_piece_counts = []
    for source_branch in range(P):
        upper = stage.weighted_fold(
            f_by_source[source_branch], stage.DCOLL // P, source_grid
        )
        upper_piece_counts.append(len(upper[0]))
        owner_windows = tuple(
            stage.extract_window(upper[0], upper[1], owner_branch, P, source_grid)
            for owner_branch in range(P)
        )
        # Order source_branch, collision_root, owner_branch.
        raw.append(tuple(
            tuple(
                stage.extract_window(
                    owner_windows[owner_branch][0],
                    owner_windows[owner_branch][1],
                    collision_root,
                    P,
                    source_grid,
                )
                for owner_branch in range(P)
            )
            for collision_root in range(P)
        ))
    raw = tuple(raw)
    scale = B.JOINT_COORDINATE // source_grid
    profiles = tuple(
        tuple(
            tuple(scale_profile(profile, scale) for profile in root_profiles)
            for root_profiles in source_profiles
        )
        for source_profiles in raw
    )

    source_parent = S.side_context()["windows"]  # b_source, collision_root
    owner_parent, owner_boundaries, _owner_record = R.owner_branch_context()
    boundaries = tuple(sorted(
        {0, B.JOINT_COORDINATE}
        | set(ctx["source_boundaries"])
        | set(S.side_context()["boundaries"])
        | set(owner_boundaries)
        | {
            point
            for source_profiles in profiles
            for root_profiles in source_profiles
            for profile in root_profiles
            for point in profile[0]
        }
    ))
    require(
        tuple(B.JOINT_COORDINATE - point for point in reversed(boundaries))
        == boundaries,
        "joint ancestry boundary reflection",
    )

    joint_support = [[False] * P for _source in range(P)]
    reflection_checks = 0
    for left, right in zip(boundaries, boundaries[1:]):
        reflected_left = B.JOINT_COORDINATE - right
        for collision_root in range(P):
            total = profile_value(ctx["source_u"][collision_root], left)
            for source_branch in range(P):
                source_sum = sum(
                    profile_value(profiles[source_branch][collision_root][owner], left)
                    for owner in range(P)
                )
                expected_source = profile_value(
                    source_parent[source_branch][collision_root], left
                )
                require(source_sum == expected_source,
                        ("source marginal", left, collision_root, source_branch))
            for owner_branch in range(P):
                owner_sum = sum(
                    profile_value(profiles[source][collision_root][owner_branch], left)
                    for source in range(P)
                )
                expected_owner = profile_value(
                    owner_parent[collision_root][owner_branch], left
                )
                require(owner_sum == expected_owner,
                        ("owner marginal", left, collision_root, owner_branch))
            joint_total = sum(
                profile_value(profiles[source][collision_root][owner], left)
                for source in range(P) for owner in range(P)
            )
            require(joint_total == total,
                    ("double marginal", left, collision_root, joint_total, total))

            for source_branch in range(P):
                for owner_branch in range(P):
                    value = profile_value(
                        profiles[source_branch][collision_root][owner_branch], left
                    )
                    if value:
                        joint_support[source_branch][owner_branch] = True
                    reflected = profile_value(
                        profiles[P - 1 - source_branch]
                                [P - 1 - collision_root]
                                [P - 1 - owner_branch],
                        reflected_left,
                    )
                    require(value == reflected,
                            ("joint profile reflection", left, source_branch,
                             collision_root, owner_branch, value, reflected))
                    reflection_checks += 1

    support_matrix = tuple(tuple(int(value) for value in row) for row in joint_support)
    require(support_matrix[0] == support_matrix[-1] == (0,) * P,
            ("source endpoint support", support_matrix[0], support_matrix[-1]))
    require(all(row == (1,) * P for row in support_matrix[1:-1]),
            ("interior joint support", support_matrix))

    zeta = ctx["zeta"]
    roots = tuple(
        tuple(pow(zeta, -frequency * value % P, PRIME) for value in range(P))
        for frequency in range(P)
    )
    samples = tuple(
        tuple(
            tuple(
                tuple(
                    profile_value(profiles[source][collision_root][owner], left)
                    for owner in range(P)
                )
                for source in range(P)
            )
            for left in boundaries[:-1]
        )
        for collision_root in range(P)
    )
    full_double_signatures = []
    for source_frequency, owner_frequency in ALL_DOUBLE_PAIRS:
        signature = []
        for collision_root in range(P):
            for matrix in samples[collision_root]:
                signature.append(sum(
                    roots[source_frequency][source]
                    * roots[owner_frequency][owner]
                    * matrix[source][owner]
                    for source in range(P) for owner in range(P)
                ) % PRIME)
        full_double_signatures.append(tuple(signature))
    full_double_signatures = tuple(full_double_signatures)
    full_double_profile_rank = B.rank_mod(full_double_signatures)
    full_double_basis_indices = row_basis_indices(full_double_signatures)
    require(len(full_double_basis_indices) == full_double_profile_rank,
            (full_double_basis_indices, full_double_profile_rank))
    full_double_basis_pairs = tuple(
        ALL_DOUBLE_PAIRS[index] for index in full_double_basis_indices
    )
    require(full_double_basis_pairs == CHARACTER_PAIRS,
            ("double-character basis drift", full_double_basis_pairs,
             CHARACTER_PAIRS))
    selected_signature_indices = tuple(
        ALL_DOUBLE_PAIRS.index(pair) for pair in CHARACTER_PAIRS
    )
    selected_profile_rank = B.rank_mod(tuple(
        full_double_signatures[index] for index in selected_signature_indices
    ))

    projected = []
    additive_hostile_checks = 0
    for source_frequency, owner_frequency in CHARACTER_PAIRS:
        root_profiles = []
        for collision_root in range(P):
            values = []
            for interval_index, left in enumerate(boundaries[:-1]):
                true_value = 0
                additive_value = 0
                whole = profile_value(ctx["source_u"][collision_root], left)
                for source_branch in range(P):
                    source_value = profile_value(
                        source_parent[source_branch][collision_root], left
                    )
                    source_phase = roots[source_frequency][source_branch]
                    for owner_branch in range(P):
                        owner_value = profile_value(
                            owner_parent[collision_root][owner_branch], left
                        )
                        phase = source_phase * roots[owner_frequency][owner_branch]
                        true_value += phase * samples[collision_root][interval_index][
                            source_branch
                        ][owner_branch]
                        # A=(B_b+R_r)/13-U/169 has both one-digit marginals.
                        additive = (
                            (source_value + owner_value) * pow(P, -1, PRIME)
                            - whole * pow(P * P, -1, PRIME)
                        ) % PRIME
                        additive_value += phase * additive
                values.append(true_value % PRIME)
                require(additive_value % PRIME == 0,
                        ("additive hostile double character", source_frequency,
                         owner_frequency, collision_root, left))
                additive_hostile_checks += 1
            root_profiles.append(compress_profile(boundaries[:-1], values))
        projected.append(tuple(root_profiles))
    projected = tuple(projected)
    projected_boundaries = tuple(sorted(
        {0, B.JOINT_COORDINATE}
        | {point for channel in projected for profile in channel for point in profile[0]}
    ))
    require(any(
        value
        for channel in projected for profile in channel for value in profile[1]
    ), "all true double-character profiles vanished")

    record = (
        tuple(upper_piece_counts), len(boundaries), len(projected_boundaries),
        support_matrix, reflection_checks, additive_hostile_checks,
        full_double_profile_rank, selected_profile_rank,
        full_double_basis_indices, full_double_basis_pairs,
        digest(raw), digest(profiles), digest(full_double_signatures), digest(projected),
    )
    return {
        "profiles": profiles,
        "boundaries": boundaries,
        "projected": projected,
        "projected_boundaries": projected_boundaries,
        "support": support_matrix,
        "record": record,
    }


def integrate_projected(alpha: int, beta: int, literal_tau: int | None = None):
    ctx = B.context()
    joint = joint_context()
    events, interval_count, mapped = B.endpoint_events(alpha, beta, literal_tau)
    for boundary in joint["projected_boundaries"]:
        events.setdefault(boundary, 0)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    rows = [
        [[0 for _state in range(V)] for _channel in CHARACTER_PAIRS]
        for _tau in tau_values
    ]
    diagonal = [
        [[0 for _state in range(V)] for _channel in CHARACTER_PAIRS]
        for _tau in tau_values
    ]
    pointed = [[0 for _point in POINTED_STATES] for _tau in tau_values]
    positions = sorted(events)
    mask = 0
    active = q_active = weighted = 0
    state_counts = [0] * V
    for position, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position + 1]
        if mask == 0 or left == right:
            continue
        chamber = B.chamber_of_segment(left, right)
        state = B.state_of_segment(left, right)
        active += 1
        state_counts[state] += 1
        jump = B.q_phase_jump(left, right)
        if not jump:
            continue
        q_active += 1
        whole_u, v_values = B.source_values(left)
        support = tuple(root for root, value in enumerate(whole_u) if value)
        if not support or not any(v_values):
            continue
        weighted += 1
        projected_by_root = tuple(
            tuple(profile_value(joint["projected"][channel][root], left)
                  for root in range(P))
            for channel in range(len(CHARACTER_PAIRS))
        )
        for row_index, tau in enumerate(tau_values):
            if literal_tau is None:
                selected = tuple(
                    sheet for sheet in range(P)
                    if (mask >> sheet) & 1 and B.guard_safe(sheet, chamber, tau)
                )
            else:
                selected = tuple(sheet for sheet in range(P) if (mask >> sheet) & 1)
            selected_support = tuple(root for root in support if root in selected)
            right_value = sum(v_values[root] for root in selected)
            if not selected_support or not right_value:
                continue
            factor = right_value * jump % PRIME
            for root in selected_support:
                point = POINT_INDEX[(state, root)]
                pointed[row_index][point] = (
                    pointed[row_index][point] + whole_u[root] * factor
                ) % PRIME
            for channel in range(len(CHARACTER_PAIRS)):
                left_value = sum(
                    projected_by_root[channel][root] for root in selected_support
                ) % PRIME
                same_root = sum(
                    projected_by_root[channel][root] * v_values[root]
                    for root in selected_support
                ) % PRIME
                require(same_root == 0,
                        ("projected same root", alpha, beta, tau, channel, left))
                rows[row_index][channel][state] = (
                    rows[row_index][channel][state] + left_value * factor
                ) % PRIME
                diagonal[row_index][channel][state] = (
                    diagonal[row_index][channel][state] + same_root * jump
                ) % PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("endpoint mask", alpha, beta, literal_tau, mask))
    require(all(value == 0 for row in diagonal for channel in row for value in channel),
            ("projected diagonal bank", alpha, beta, literal_tau))
    freeze = lambda bank: tuple(
        tuple(tuple(state_row) for state_row in row) for row in bank
    )
    return (
        freeze(rows), freeze(diagonal), tuple(tuple(row) for row in pointed),
        (interval_count, mapped, active, q_active, weighted, tuple(state_counts)),
    )


def worker(alpha: int):
    zeta = B.context()["zeta"]
    rows = []
    diagonal_rows = []
    all_pointed_rows = []
    scalar_counts = [0] * 5
    state_counts = [0] * V
    for beta in range(P):
        coupled, diagonal, pointed, counts = integrate_projected(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        rows.append(tuple(
            tuple(tuple(phase * value % PRIME for value in state_row)
                  for state_row in row)
            for row in coupled
        ))
        diagonal_rows.append(diagonal)
        pointed_rows = tuple(
            tuple(phase * value % PRIME for value in row) for row in pointed
        )
        all_pointed_rows.append(pointed_rows)
        scalar_counts = [left + right for left, right in zip(scalar_counts, counts[:5])]
        state_counts = [left + right for left, right in zip(state_counts, counts[5])]
    return (
        alpha, tuple(rows), tuple(diagonal_rows), tuple(all_pointed_rows),
        (tuple(scalar_counts), tuple(state_counts)),
    )


def inverse_tensor(gamma, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    tensor = [
        [[0 for _relation in range(P)] for _state in range(V)]
        for _channel in CHARACTER_PAIRS
    ]
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                phases = tuple(
                    pow(zeta, -(alpha + tau * relation) % P, PRIME)
                    for relation in range(P)
                )
                for channel in range(len(CHARACTER_PAIRS)):
                    for state in range(V):
                        value = row[channel][state]
                        if not value:
                            continue
                        for relation, phase in enumerate(phases):
                            tensor[channel][state][relation] = (
                                tensor[channel][state][relation] + value * phase
                            ) % PRIME
    require(index == P**3, ("gamma size", index))
    return tuple(
        tuple(
            tuple(value * normalizer % PRIME for value in row)
            for row in channel
        )
        for channel in tensor
    )


def inverse_pointed(gamma, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    tensor = [[0 for _relation in range(P)] for _point in POINTED_STATES]
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                for relation in range(P):
                    phase = pow(zeta, -(alpha + tau * relation) % P, PRIME)
                    for point, value in enumerate(row):
                        tensor[point][relation] = (
                            tensor[point][relation] + value * phase
                        ) % PRIME
    require(index == P**3, ("pointed gamma size", index))
    return tuple(
        tuple(value * normalizer % PRIME for value in row) for row in tensor
    )


def walsh_relation_centre(tensor):
    inverse_p = pow(P, -1, PRIME)
    answer = []
    for channel in range(len(CHARACTER_PAIRS)):
        channel_rows = []
        for character in range(1, V):
            row = []
            for relation in range(P):
                total = 0
                for state in range(V):
                    parity = (
                        (B.STATE_LABELS[state][0] & B.STATE_LABELS[character][0])
                        ^ (B.STATE_LABELS[state][1] & B.STATE_LABELS[character][1])
                    )
                    total += (-1 if parity else 1) * tensor[channel][state][relation]
                row.append(total % PRIME)
            mean = sum(row) % PRIME * inverse_p % PRIME
            channel_rows.append(tuple((value - mean) % PRIME for value in row))
        answer.append(tuple(channel_rows))
    return tuple(answer)


def conditional_spectrum(four_way, zeta: int):
    spectrum = []
    for channel in range(len(CHARACTER_PAIRS)):
        channel_rows = []
        for character in range(V - 1):
            row = tuple(
                sum(
                    four_way[channel][character][relation]
                    * pow(zeta, -frequency * relation % P, PRIME)
                    for relation in range(P)
                ) % PRIME
                for frequency in range(1, P)
            )
            channel_rows.append(row)
        spectrum.append(tuple(channel_rows))
    return tuple(spectrum)


def channel_subset_record(four_way):
    per_channel = tuple(
        B.rank_mod(tuple(four_way[channel][character] for character in range(V - 1)))
        for channel in range(len(CHARACTER_PAIRS))
    )
    target = B.rank_mod(tuple(
        four_way[channel][character]
        for channel in range(len(CHARACTER_PAIRS))
        for character in range(V - 1)
    ))
    witnesses = []
    for size in range(1, len(CHARACTER_PAIRS) + 1):
        for subset in combinations(range(len(CHARACTER_PAIRS)), size):
            matrix = tuple(
                four_way[channel][character]
                for channel in subset for character in range(V - 1)
            )
            if B.rank_mod(matrix) == target:
                witnesses.append(subset)
        if witnesses:
            break
    return per_channel, target, tuple(witnesses)


def solve_square(matrix, vector):
    size = len(matrix)
    rows = [
        [value % PRIME for value in matrix[index]] + [vector[index] % PRIME]
        for index in range(size)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if rows[row][column]), None)
        require(pivot is not None, ("singular coordinate matrix", matrix))
        rows[column], rows[pivot] = rows[pivot], rows[column]
        inverse = pow(rows[column][column], -1, PRIME)
        rows[column] = [value * inverse % PRIME for value in rows[column]]
        for row in range(size):
            if row == column:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    (value - factor * pivot_value) % PRIME
                    for value, pivot_value in zip(rows[row], rows[column])
                ]
    return tuple(rows[index][-1] for index in range(size))


def pointed_rowspace_record(tensor, pointed_tensor):
    records = []
    coordinates = []
    for state in range(V):
        joint_rows = tuple(
            tensor[channel][state] for channel in range(len(CHARACTER_PAIRS))
        )
        point_rows = tuple(pointed_tensor[point] for point in STATE_POINTS[state])
        point_count = len(point_rows)
        ranks = (
            point_count,
            B.rank_mod(joint_rows),
            B.rank_mod(point_rows),
            B.rank_mod(joint_rows + point_rows),
        )
        records.append(ranks)
        require(ranks == (point_count,) * 4,
                ("nested digits leave pointed relation row space", state, ranks))

        transposed = tuple(
            tuple(point_rows[row][column] for row in range(point_count))
            for column in range(P)
        )
        pivot_columns = row_basis_indices(transposed)
        require(len(pivot_columns) == point_count,
                ("pointed coordinate pivots", state, pivot_columns))
        square = tuple(
            tuple(point_rows[row][column] for row in range(point_count))
            for column in pivot_columns
        )
        state_coordinates = []
        for target in joint_rows:
            coefficients = solve_square(
                square, tuple(target[column] for column in pivot_columns)
            )
            reconstruction = tuple(
                sum(coefficients[row] * point_rows[row][column]
                    for row in range(point_count)) % PRIME
                for column in range(P)
            )
            require(reconstruction == target,
                    ("pointed coordinate reconstruction", state, coefficients))
            state_coordinates.append(coefficients)
        coordinates.append((pivot_columns, tuple(state_coordinates)))
    joint_rows = tuple(
        tensor[channel][state]
        for channel in range(len(CHARACTER_PAIRS)) for state in range(V)
    )
    point_rows = tuple(pointed_tensor[point] for point in range(len(POINTED_STATES)))
    global_record = (
        len(POINTED_STATES), B.rank_mod(joint_rows), B.rank_mod(point_rows),
        B.rank_mod(joint_rows + point_rows),
    )
    require(global_record == (len(POINTED_STATES),) * 4,
            ("nested digits leave global pointed relation row space",
             global_record))
    return tuple(records), global_record, tuple(coordinates)


def main() -> None:
    joint = joint_context()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    gamma = tuple(
        row
        for chunk in chunks
        for beta_rows in chunk[1]
        for row in beta_rows
    )
    diagonal = tuple(
        row
        for chunk in chunks
        for beta_rows in chunk[2]
        for row in beta_rows
    )
    pointed_gamma = tuple(
        row
        for chunk in chunks
        for beta_rows in chunk[3]
        for row in beta_rows
    )
    require(len(gamma) == len(diagonal) == len(pointed_gamma) == P**3,
            (len(gamma), len(diagonal), len(pointed_gamma)))
    require(all(value == 0 for row in diagonal for channel in row for value in channel),
            "global projected diagonal")

    direct_records = []
    zeta = B.context()["zeta"]
    for alpha, beta, tau in CONTROL_TRIPLES:
        direct, direct_diagonal, direct_pointed, counts = integrate_projected(
            alpha, beta, tau
        )
        phase = pow(zeta, beta, PRIME)
        direct_row = tuple(
            tuple(phase * value % PRIME for value in state_row)
            for state_row in direct[0]
        )
        index = (alpha * P + beta) * P + tau
        require(direct_row == gamma[index], ("literal guard", alpha, beta, tau))
        require(all(value == 0 for channel in direct_diagonal[0] for value in channel),
                ("literal diagonal", alpha, beta, tau))
        direct_pointed_row = tuple(
            phase * value % PRIME for value in direct_pointed[0]
        )
        require(direct_pointed_row == pointed_gamma[index],
                ("literal pointed", alpha, beta, tau))
        direct_records.append(((alpha, beta, tau), counts))

    tensor = inverse_tensor(gamma, zeta)
    pointed_tensor = inverse_pointed(pointed_gamma, zeta)
    pointed_parent = tuple(
        tuple(
            sum(pointed_tensor[point][relation] for point in STATE_POINTS[state])
            % PRIME
            for relation in range(P)
        )
        for state in range(V)
    )
    require(digest(pointed_parent) == B.EXPECTED_DIGESTS[2],
            ("pointed relation marginal", digest(pointed_parent),
             B.EXPECTED_DIGESTS[2]))
    pointed_record, pointed_global_record, pointed_coordinates = pointed_rowspace_record(
        tensor, pointed_tensor
    )
    four_way = walsh_relation_centre(tensor)
    spectrum = conditional_spectrum(four_way, zeta)
    raw_matrix = tuple(
        tensor[channel][state]
        for channel in range(len(CHARACTER_PAIRS)) for state in range(V)
    )
    four_way_matrix = tuple(
        four_way[channel][character]
        for channel in range(len(CHARACTER_PAIRS))
        for character in range(V - 1)
    )
    raw_rank = B.rank_mod(raw_matrix)
    conditional_rank = B.rank_mod(four_way_matrix)
    raw_channel_matrix = tuple(
        tuple(
            tensor[channel][state][relation]
            for state in range(V) for relation in range(P)
        )
        for channel in range(len(CHARACTER_PAIRS))
    )
    four_way_channel_matrix = tuple(
        tuple(
            four_way[channel][character][relation]
            for character in range(V - 1) for relation in range(P)
        )
        for channel in range(len(CHARACTER_PAIRS))
    )
    raw_channel_rank = B.rank_mod(raw_channel_matrix)
    four_way_channel_rank = B.rank_mod(four_way_channel_matrix)
    support = sum(
        value != 0
        for channel in spectrum for character in channel for value in character
    )
    channel_ranks, subset_target, minimal_witnesses = channel_subset_record(four_way)
    fixed_relation = 6
    fixed_matrix = tuple(
        tuple(tensor[channel][state][fixed_relation] for state in range(V))
        for channel in range(len(CHARACTER_PAIRS))
    )
    fixed_rank = B.rank_mod(fixed_matrix)
    require(raw_rank >= conditional_rank > 0, (raw_rank, conditional_rank))
    full_profile_rank = joint["record"][6]
    selected_profile_rank = joint["record"][7]
    require(selected_profile_rank == full_profile_rank,
            ("selected characters miss double-profile span",
             selected_profile_rank, full_profile_rank))
    require(raw_channel_rank <= full_profile_rank and
            four_way_channel_rank <= full_profile_rank,
            ("response cannot exceed source-profile basis",
             full_profile_rank, raw_channel_rank, four_way_channel_rank))
    require(
        (full_profile_rank, selected_profile_rank, raw_rank, conditional_rank,
         raw_channel_rank, four_way_channel_rank, support, fixed_rank)
        == (17, 17, 6, 6, 4, 4, 612, 4),
        ("rank/support drift", full_profile_rank, selected_profile_rank,
         raw_rank, conditional_rank, raw_channel_rank,
         four_way_channel_rank, support, fixed_rank),
    )
    require(subset_target == conditional_rank and minimal_witnesses,
            ("conditional subset witness", subset_target, conditional_rank,
             channel_ranks, minimal_witnesses))

    scalar_counts = tuple(sum(chunk[4][0][index] for chunk in chunks) for index in range(5))
    state_counts = tuple(sum(chunk[4][1][state] for chunk in chunks) for state in range(V))
    require(scalar_counts == S.EXPECTED_WORK_COUNTS,
            ("work count drift", scalar_counts, S.EXPECTED_WORK_COUNTS))

    digests = (
        digest(gamma), digest(tensor), digest(pointed_gamma), digest(pointed_tensor),
        digest(pointed_coordinates), digest(four_way), digest(spectrum),
        digest(fixed_matrix),
    )
    semantic_surface = (
        SOURCE_SHA256, OWNER_SHA256, CHARACTER_PAIRS, joint["record"],
        tuple(direct_records), scalar_counts, state_counts, raw_rank,
        conditional_rank, raw_channel_rank, four_way_channel_rank, support,
        channel_ranks, minimal_witnesses, pointed_record, pointed_global_record,
        pointed_coordinates, fixed_relation, fixed_rank, digests,
    )
    semantic = digest(semantic_surface)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 nested source-sheet x inverse-owner digit conditional probe ==")
    print(f"parents=(source_sha={SOURCE_SHA256},owner_sha={OWNER_SHA256})")
    print("stage_word=(n=c+13*b_source before Q(X);a=r_owner+13*h in outer P_(13^5);X=(w_u+a)/13^5;Z=(X+n)/169)")
    print("typing=both labels coexist on THM-2471 Boolean stalk (y,u,a,n,e');no source/arrival atom or temporal intertwiner used")
    print(f"joint_source=(upper_piece_counts,boundaries,projected_boundaries,support,reflection_checks,additive_hostile_checks,all_double_profile_rank,selected_profile_rank,basis_indices,basis_pairs,raw_sha,profile_sha,all_double_sha,projected_sha)={joint['record']}")
    print("marginals=(sum_r_owner->audited b_source profiles,sum_b_source->current-leg r_owner candidate profiles,double_sum->U_u): POINTWISE PASS")
    print("reflection=(y,u,b_source,r_owner)->(1-y,12-u,12-b_source,12-r_owner): PASS")
    print(f"conditional_character_bank={CHARACTER_PAIRS};both_digit_frequencies_nonzero=True")
    print(f"work_counts={scalar_counts};state_counts={state_counts};literal_controls={CONTROL_TRIPLES}: PASS;same_root=0")
    print(f"conditional_ranks=(raw_character_state_x_relation={raw_rank},pure_four_way_character_state_x_relation={conditional_rank})")
    print(f"channel_to_output_ranks=(source_profile_basis={full_profile_rank},raw_state_x_relation={raw_channel_rank},pure_four_way_state_x_relation={four_way_channel_rank},pure_kernel_dimension={full_profile_rank-four_way_channel_rank})")
    print("inherited_rank_ceiling=relation rows equal the pointed-six carrier exactly;the joint interaction supplies four amplitude-coordinate directions inside it")
    print(f"per_character_pair_four_way_ranks={channel_ranks}")
    print(f"minimal_conditional_rank_witnesses=(count={len(minimal_witnesses)},first_indices={minimal_witnesses[0]},first_pairs={tuple(CHARACTER_PAIRS[i] for i in minimal_witnesses[0])},sha256={digest(minimal_witnesses)})")
    print(f"pure_four_way_relation_fourier_support={support}/{len(CHARACTER_PAIRS)*(V-1)*(P-1)}")
    print(f"pointed_relation_rowspace_(point_count,joint_rank,point_rank,union_rank)={pointed_record}: PASS")
    print(f"pointed_global_relation_rowspace_(point_count,joint_rank,point_rank,union_rank)={pointed_global_record}: PASS")
    print(f"pointed_coordinate_map=(sha256={digest(pointed_coordinates)},pivot_columns={tuple(record[0] for record in pointed_coordinates)}): PASS")
    print(f"fixed_relation_(1,0,{fixed_relation})_double_character_x_state_rank={fixed_rank}")
    print("hostile=additive lift A_(b,r)=B_b/13+R_r/13-U/169 has both exact one-digit marginals and zero every double-nontrivial character: PASS")
    print(f"digests_(gamma,tensor,pointed_gamma,pointed_tensor,pointed_coordinates,four_way,spectrum,fixed)={digests}")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=FINITE-EXACT lawful nested ancestry construction with conditional-rank and pointed-carrier gates on one r=5 owner base")
    print("scope=relation-output carrier comparison only;root difference marginalized;not equality of digits,not exact C(a;X,m),not chronology,not physical current,no row exclusion,no LRC(14)")
    print("commands=python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_nested_ancestry_digits_probe_20260816.py;python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_nested_ancestry_digits_probe_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
