#!/usr/bin/env python3
"""Retain one inverse-branch ancestry sidecar over the audited Boolean square.

The sidecar is the high base-13 digit

    b = floor(n/13),  0 <= n < 169,

of the inverse branch in the left P_169 source packet.  THM-2471's temporal
transition identifies this digit with the source guard sheet.  The program
splits the source transfer before P_(13^5), couples the thirteen pieces to the
actual endpoint integrand, and retains (b, Boolean-square state, residue).

This is source-time finite ancestry information, not an exact address
C(a;X,m), arrival atom, temporal intertwiner, U_clock chronology, physical
current, row exclusion, or LRC(14) theorem.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / (
    "04-computation/lrc_r5_ufull_owner_node_boolean_square_"
    "independent_audit_20260816.py"
)
BASE_SHA256 = "e585c1a18d846c4ab87fc159ad411bad90859be831831cea1c0e8f326d9b440f"
BASE_SEMANTIC = "b0996af3f1760b2118187490c93e0e01b322cb57fd6b25d3bf3688778b7e664c"

P = 13
V = 4
ROWS = P * V
EXPECTED_PARENT_GAMMA = "8b8f2a2785b084e1578ba0512e4577ab79fd674b84b588fbdb8186f2009242c2"
CONTROL_TRIPLES = ((0, 0, 0), (1, 0, 6), (12, 12, 0))
EXPECTED_NONEMPTY = (0,) + (1,) * 11 + (0,)
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_FULL_SHAPE = (
    ("dc", 1), ("branch", 12), ("state", 3), ("residue", 12),
    ("branch_state", 36), ("branch_residue", 144),
    ("state_residue", 36), ("triple", 432),
)
EXPECTED_THREE_WAY_SHAPE = (
    ("dc", 0), ("branch", 0), ("state", 0), ("residue", 0),
    ("branch_state", 0), ("branch_residue", 0),
    ("state_residue", 0), ("triple", 432),
)
EXPECTED_DIGESTS = (
    "ba0143302119c93b7b7b4f778478c14c265efb8c0bcb9c2945e9eaae9009b493",
    "63d9304d28689ffa04b1568d538c67f3da59ecdc7e0d1c6535ad323c800ac0e5",
    "504e91431cf4d55f6e3cb5aa7c6e6f2c34db538cb692c67857a5c7627b252261",
    "61cf3999a1ef902db958851e99b8ea923c20099d81464c159613c0cf0054efdf",
    "8874b55a71e975e20c7b2b29efd9eec0ea93dcaecc0209d43c61712bbc9c0acf",
    "5d02353fd1cfa642be290264353e5e02441a2b69fa4144be72a7ae5ed9776396",
    "d8d388655fe34b27122c92da84c3cfb29d96303451334fbb4fe9d8c8fc13f2b8",
    "a70c2937fa755de450cc913d031b7e5f1eaf4b756900e79832db247582f76cba",
)
EXPECTED_SEMANTIC_SHA256 = "c82d9aba055f2f5ccbea75e092dd7a45b9a9f2a4a1b91762f32ad5852cf9c552"

_SIDE_CONTEXT = None


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def load_base():
    actual = lf_sha256(BASE_PATH)
    require(actual == BASE_SHA256, (actual, BASE_SHA256))
    spec = importlib.util.spec_from_file_location("branch_sheet_boolean_square_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "base loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules["branch_sheet_boolean_square_base"] = module
    spec.loader.exec_module(module)
    require(module.EXPECTED_SEMANTIC_SHA256 == BASE_SEMANTIC, "base semantic")
    return module


B = load_base()


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def sum_profile_values(profiles, point: int) -> int:
    return sum(profile_value(profile, point) for profile in profiles)


def inverse_branch_source_sheet_certificate():
    cells = 91 * P**2
    branch_counts = [0] * P
    for cell in range(cells):
        midpoint_numerator = 2 * cell + 1
        # s=91t=(2*cell+1)/(2*169).  Its guard sheet is floor(s/7).
        source_sheet = midpoint_numerator // (2 * P**2 * 7)
        branch = midpoint_numerator // (2 * 91)
        high_digit = branch // P
        require(0 <= branch < P**2, (cell, branch))
        require(source_sheet == high_digit, (cell, source_sheet, branch, high_digit))
        branch_counts[high_digit] += 1
    require(tuple(branch_counts) == (91 * P,) * P, branch_counts)
    return cells, tuple(branch_counts)


def branch_numerator_profiles(M, e_intervals):
    """Split sum_n 1_E((y+n)/169) by b=floor(n/13)."""
    grid = M.T_DEN
    unit = grid // M.RPKT
    groups = [[] for _branch in range(P)]
    input_measure = 0
    lifted_measure = 0
    for interval_left, interval_right in e_intervals:
        input_measure += interval_right - interval_left
        first = interval_left // unit
        last = (interval_right - 1) // unit
        for inverse_branch in range(first, last + 1):
            left = max(interval_left, inverse_branch * unit)
            right = min(interval_right, (inverse_branch + 1) * unit)
            if left >= right:
                continue
            mapped_left = M.RPKT * left - inverse_branch * grid
            mapped_right = M.RPKT * right - inverse_branch * grid
            require(0 <= mapped_left < mapped_right <= grid,
                    (inverse_branch, mapped_left, mapped_right))
            groups[inverse_branch // P].append((mapped_left, mapped_right, 1))
            lifted_measure += mapped_right - mapped_left
    require(lifted_measure == M.RPKT * input_measure,
            (lifted_measure, M.RPKT * input_measure))
    return tuple(M.weighted_fold(group, 1, grid) for group in groups)


def build_branch_windows():
    ctx = B.context()
    M = ctx["M"]
    e_intervals = M.build_set(M.PAT_E, M.ZELL)
    q_intervals = M.build_set(M.PAT_QB, M.ZELL)
    p2_branches = branch_numerator_profiles(M, e_intervals)
    f_branches = tuple(
        B.restrict_profile(starts, values, q_intervals, M.T_DEN)
        for starts, values in p2_branches
    )
    u_branches = tuple(
        M.weighted_fold(pieces, M.DCOLL, M.T_DEN) for pieces in f_branches
    )
    raw_windows = tuple(
        tuple(M.extract_window(starts, values, root, P, M.T_DEN) for root in range(P))
        for starts, values in u_branches
    )
    scale = B.JOINT_COORDINATE // M.T_DEN
    windows = tuple(
        tuple(
            (tuple(point * scale for point in starts), tuple(values))
            for starts, values in branch_roots
        )
        for branch_roots in raw_windows
    )
    boundaries = tuple(sorted(
        {0, B.JOINT_COORDINATE}
        | {point for branch in windows for profile in branch for point in profile[0]}
        | set(ctx["source_boundaries"])
    ))

    # Exact linear restoration of every one of the thirteen root profiles.
    for root in range(P):
        profiles = tuple(windows[branch][root] for branch in range(P))
        check_points = sorted(
            {0, B.JOINT_COORDINATE}
            | {point for profile in profiles for point in profile[0]}
            | set(ctx["source_u"][root][0])
        )
        require(all(
            sum_profile_values(profiles, point)
            == profile_value(ctx["source_u"][root], point)
            for point in check_points[:-1]
        ), ("branch restoration", root))
    nonempty = tuple(
        int(any(value for profile in windows[branch] for value in profile[1]))
        for branch in range(P)
    )
    require(nonempty == EXPECTED_NONEMPTY, nonempty)
    return windows, boundaries, digest(windows), nonempty


def side_context():
    global _SIDE_CONTEXT
    if _SIDE_CONTEXT is None:
        windows, boundaries, profile_digest, nonempty = build_branch_windows()
        _SIDE_CONTEXT = {
            "windows": windows,
            "boundaries": boundaries,
            "profile_digest": profile_digest,
            "nonempty": nonempty,
        }
    return _SIDE_CONTEXT


def branch_values_on_support(point: int, support):
    windows = side_context()["windows"]
    return {
        root: tuple(profile_value(windows[branch][root], point) for branch in range(P))
        for root in support
    }


def geometric_branch_state_support():
    side = side_context()
    support = [[False] * V for _branch in range(P)]
    components = [[0] * V for _branch in range(P)]
    active_previous = [[False] * V for _branch in range(P)]
    for left, right in zip(side["boundaries"], side["boundaries"][1:]):
        midpoint_twice = left + right
        in_owner = (
            B.Q * midpoint_twice < B.JOINT_COORDINATE
            or B.Q * midpoint_twice > 13 * B.JOINT_COORDINATE
        )
        current = [[False] * V for _branch in range(P)]
        if in_owner:
            state = B.state_of_segment(left, right)
            for branch in range(P):
                branch_nonzero = any(
                    profile_value(side["windows"][branch][root], left)
                    for root in range(P)
                )
                if branch_nonzero:
                    support[branch][state] = True
                    current[branch][state] = True
                    if not active_previous[branch][state]:
                        components[branch][state] += 1
        active_previous = current
    support_counts_by_state = tuple(
        sum(support[branch][state] for branch in range(P)) for state in range(V)
    )
    support_counts_by_branch = tuple(sum(row) for row in support)
    require(all(count > 1 for count in support_counts_by_state), support_counts_by_state)
    require(len(set(tuple(row) for row in support)) > 1, "branch support rows collapsed")
    return (
        tuple(tuple(int(value) for value in row) for row in support),
        tuple(tuple(row) for row in components),
        support_counts_by_state,
        support_counts_by_branch,
    )


def integrate_sidecar(alpha: int, beta: int, literal_tau: int | None = None):
    events, interval_count, mapped = B.endpoint_events(alpha, beta, literal_tau)
    for boundary in side_context()["boundaries"]:
        events.setdefault(boundary, 0)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    result = [[0] * ROWS for _tau in tau_values]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    for position, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position + 1]
        if mask == 0 or left == right:
            continue
        chamber = B.chamber_of_segment(left, right)
        state = B.state_of_segment(left, right)
        active_segments += 1
        jump = B.q_phase_jump(left, right)
        if jump == 0:
            continue
        q_active_segments += 1
        whole_u, v_values = B.source_values(left)
        root_support = tuple(root for root, value in enumerate(whole_u) if value)
        branch_by_root = branch_values_on_support(left, root_support)
        require(all(
            sum(branch_by_root[root]) == whole_u[root] for root in root_support
        ), ("local branch restoration", left, root_support))
        weighted_segments += 1
        for row_index, tau in enumerate(tau_values):
            if literal_tau is None:
                selected = tuple(
                    sheet for sheet in range(P)
                    if (mask >> sheet) & 1 and B.guard_safe(sheet, chamber, tau)
                )
            else:
                selected = tuple(sheet for sheet in range(P) if (mask >> sheet) & 1)
            selected_source_roots = tuple(root for root in root_support if root in selected)
            right_value = sum(v_values[sheet] for sheet in selected)
            left_by_branch = tuple(
                sum(branch_by_root[root][branch] for root in selected_source_roots)
                for branch in range(P)
            )
            require(sum(left_by_branch) == sum(whole_u[root] for root in selected_source_roots),
                    ("selected branch restoration", alpha, beta, tau, left))
            for root in selected_source_roots:
                require(v_values[root] == 0, ("branch diagonal", root, left))
            for branch, left_value in enumerate(left_by_branch):
                result[row_index][branch * V + state] = (
                    result[row_index][branch * V + state]
                    + left_value * right_value * jump
                ) % B.PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("mask closure", alpha, beta, literal_tau, mask))
    return (
        tuple(tuple(row) for row in result),
        (interval_count, mapped, active_segments, q_active_segments, weighted_segments),
    )


def worker(alpha: int):
    zeta = B.context()["zeta"]
    rows = []
    counts = [0] * 5
    for beta in range(P):
        values, record = integrate_sidecar(alpha, beta)
        phase = pow(zeta, beta, B.PRIME)
        rows.append(tuple(
            tuple(phase * value % B.PRIME for value in row) for row in values
        ))
        counts = [left + right for left, right in zip(counts, record)]
    return alpha, tuple(rows), tuple(counts)


def inverse_table(gamma, zeta: int):
    normalizer = pow(P**3, -1, B.PRIME)
    table = [[0] * P for _row in range(ROWS)]
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                for residue in range(P):
                    phase = pow(zeta, -(alpha + tau * residue) % P, B.PRIME)
                    for coordinate in range(ROWS):
                        table[coordinate][residue] = (
                            table[coordinate][residue] + row[coordinate] * phase
                        ) % B.PRIME
    require(index == P**3, index)
    return tuple(tuple(value * normalizer % B.PRIME for value in row) for row in table)


def triple_fourier(table, zeta: int):
    answer = []
    for branch_frequency in range(P):
        branch_block = []
        for channel in range(V):
            channel_row = []
            for residue_frequency in range(P):
                total = 0
                for branch in range(P):
                    branch_phase = pow(zeta, (-branch_frequency * branch) % P, B.PRIME)
                    for state in range(V):
                        parity = (
                            (B.STATE_LABELS[state][0] & B.STATE_LABELS[channel][0])
                            ^ (B.STATE_LABELS[state][1] & B.STATE_LABELS[channel][1])
                        )
                        sign = -1 if parity else 1
                        row = table[branch * V + state]
                        total += branch_phase * sign * sum(
                            row[residue]
                            * pow(zeta, (-residue_frequency * residue) % P, B.PRIME)
                            for residue in range(P)
                        )
                channel_row.append(total % B.PRIME)
            branch_block.append(tuple(channel_row))
        answer.append(tuple(branch_block))
    return tuple(answer)


def spectrum_shape(spectrum):
    counts = {
        "dc": 0,
        "branch": 0,
        "state": 0,
        "residue": 0,
        "branch_state": 0,
        "branch_residue": 0,
        "state_residue": 0,
        "triple": 0,
    }
    for branch_frequency in range(P):
        for channel in range(V):
            for residue_frequency in range(P):
                if spectrum[branch_frequency][channel][residue_frequency] == 0:
                    continue
                flags = (branch_frequency != 0, channel != 0, residue_frequency != 0)
                key = {
                    (False, False, False): "dc",
                    (True, False, False): "branch",
                    (False, True, False): "state",
                    (False, False, True): "residue",
                    (True, True, False): "branch_state",
                    (True, False, True): "branch_residue",
                    (False, True, True): "state_residue",
                    (True, True, True): "triple",
                }[flags]
                counts[key] += 1
    return tuple((key, counts[key]) for key in counts)


def three_way_interaction(table):
    tensor = [
        [[table[branch * V + state][residue] for residue in range(P)]
         for state in range(V)]
        for branch in range(P)
    ]
    inv_b = pow(P, -1, B.PRIME)
    inv_s = pow(V, -1, B.PRIME)
    inv_r = inv_b
    answer = [[[0] * P for _state in range(V)] for _branch in range(P)]
    grand = sum(value for branch in tensor for state in branch for value in state) % B.PRIME
    grand = grand * inv_b % B.PRIME * inv_s % B.PRIME * inv_r % B.PRIME
    for branch in range(P):
        for state in range(V):
            for residue in range(P):
                p_b = sum(tensor[other][state][residue] for other in range(P)) % B.PRIME * inv_b % B.PRIME
                p_s = sum(tensor[branch][other][residue] for other in range(V)) % B.PRIME * inv_s % B.PRIME
                p_r = sum(tensor[branch][state]) % B.PRIME * inv_r % B.PRIME
                p_bs = sum(
                    tensor[other_b][other_s][residue]
                    for other_b in range(P) for other_s in range(V)
                ) % B.PRIME * inv_b % B.PRIME * inv_s % B.PRIME
                p_br = sum(
                    tensor[other_b][state][other_r]
                    for other_b in range(P) for other_r in range(P)
                ) % B.PRIME * inv_b % B.PRIME * inv_r % B.PRIME
                p_sr = sum(
                    tensor[branch][other_s][other_r]
                    for other_s in range(V) for other_r in range(P)
                ) % B.PRIME * inv_s % B.PRIME * inv_r % B.PRIME
                answer[branch][state][residue] = (
                    tensor[branch][state][residue] - p_b - p_s - p_r
                    + p_bs + p_br + p_sr - grand
                ) % B.PRIME
    require(all(
        sum(answer[branch][state][residue] for branch in range(P)) % B.PRIME == 0
        for state in range(V) for residue in range(P)
    ), "three-way branch marginals")
    require(all(
        sum(answer[branch][state][residue] for state in range(V)) % B.PRIME == 0
        for branch in range(P) for residue in range(P)
    ), "three-way state marginals")
    require(all(
        sum(answer[branch][state]) % B.PRIME == 0
        for branch in range(P) for state in range(V)
    ), "three-way residue marginals")
    return tuple(
        tuple(tuple(row) for row in branch) for branch in answer
    )


def flatten_tensor(tensor):
    return tuple(tensor[branch][state] for branch in range(P) for state in range(V))


def main() -> None:
    temporal_certificate = inverse_branch_source_sheet_certificate()
    geometry = geometric_branch_state_support()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    gamma = tuple(row for chunk in chunks for beta_rows in chunk[1] for row in beta_rows)
    require(len(gamma) == P**3, len(gamma))
    parent_gamma = tuple(
        tuple(
            sum(row[branch * V + state] for branch in range(P)) % B.PRIME
            for state in range(V)
        )
        for row in gamma
    )
    parent_hash = digest(parent_gamma)
    require(parent_hash == EXPECTED_PARENT_GAMMA, parent_hash)

    literal_records = []
    for alpha, beta, tau in CONTROL_TRIPLES:
        direct, counts = integrate_sidecar(alpha, beta, tau)
        phase = pow(B.context()["zeta"], beta, B.PRIME)
        direct_row = tuple(phase * value % B.PRIME for value in direct[0])
        index = (alpha * P + beta) * P + tau
        require(direct_row == gamma[index], ("literal guard", alpha, beta, tau))
        literal_records.append(((alpha, beta, tau), counts))

    table = inverse_table(gamma, B.context()["zeta"])
    table_rank = B.rank_mod(table)
    parent_table = tuple(
        tuple(
            sum(table[branch * V + state][residue] for branch in range(P)) % B.PRIME
            for residue in range(P)
        )
        for state in range(V)
    )
    require(digest(parent_table) == B.EXPECTED_DIGESTS[2], "parent table marginal")
    require(B.rank_mod(parent_table) == 4, "parent rank")
    require(table_rank > 4, ("sidecar rank did not increase", table_rank))

    spectrum = triple_fourier(table, B.context()["zeta"])
    shape = spectrum_shape(spectrum)
    three_way = three_way_interaction(table)
    three_way_flat = flatten_tensor(three_way)
    three_way_rank = B.rank_mod(three_way_flat)
    three_way_spectrum = triple_fourier(three_way_flat, B.context()["zeta"])
    three_way_shape = spectrum_shape(three_way_spectrum)
    triple_modes = dict(three_way_shape)["triple"]
    require(three_way_rank > 0 and triple_modes > 0,
            ("three-way interaction vanished", three_way_rank, three_way_shape))

    fixed_residue = tuple(
        tuple(table[branch * V + state][6] for state in range(V))
        for branch in range(P)
    )
    fixed_rank = B.rank_mod(fixed_residue)
    fixed_support = sum(value != 0 for row in fixed_residue for value in row)
    branch_marginal = tuple(
        tuple(
            sum(table[branch * V + state][residue] for state in range(V)) % B.PRIME
            for residue in range(P)
        )
        for branch in range(P)
    )
    branch_marginal_rank = B.rank_mod(branch_marginal)
    require(branch_marginal_rank > 1, branch_marginal_rank)

    work_counts = tuple(sum(chunk[2][index] for chunk in chunks) for index in range(5))
    digests = (
        side_context()["profile_digest"], digest(gamma), digest(table), digest(spectrum),
        digest(three_way), digest(three_way_spectrum), digest(fixed_residue),
        digest(branch_marginal),
    )
    require(work_counts == EXPECTED_WORK_COUNTS, work_counts)
    require((table_rank, branch_marginal_rank) == (6, 3),
            (table_rank, branch_marginal_rank))
    require(shape == EXPECTED_FULL_SHAPE, shape)
    require((three_way_rank, three_way_shape) == (6, EXPECTED_THREE_WAY_SHAPE),
            (three_way_rank, three_way_shape))
    require((fixed_rank, fixed_support) == (3, 44), (fixed_rank, fixed_support))
    require(digests == EXPECTED_DIGESTS, digests)
    proof = (
        "inverse_branch_high_digit_equals_source_guard_sheet",
        "thirteen_nonnegative_branch_profiles_restore_every_source_root_profile",
        "branch_state_geometry_is_not_a_function_of_the_Boolean_square",
        "branch_marginal_recovers_the_audited_Boolean_square_exactly",
        "adjoining_branch_sheet_strictly_raises_output_rank",
        "nonzero_three_way_branch_state_residue_interaction",
    )
    boundary = (
        "source_time_inverse_branch_sheet_not_arrival_atom",
        "finite_ancestry_sidecar_not_exact_C_a_X_m",
        "no_temporal_intertwiner_U_clock_uniform_rows_row_exclusion_or_LRC14",
    )
    semantic_surface = (
        BASE_SHA256, BASE_SEMANTIC, temporal_certificate, side_context()["nonempty"], geometry,
        tuple(literal_records), parent_hash, work_counts, table_rank,
        shape, three_way_rank, three_way_shape, fixed_rank, fixed_support,
        branch_marginal_rank, digests, proof, boundary,
    )
    semantic = hashlib.sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("R5 BOOLEAN SQUARE + INVERSE-BRANCH SHEET -- SIDECAR PROBE")
    print("status=FINITE_EXACT_NONREDUNDANT_SOURCE_ANCESTRY_SIDECAR;LRC14_OPEN")
    print(f"base=(sha256={BASE_SHA256},semantic={BASE_SEMANTIC})")
    print(f"sidecar=high_base13_digit_b=floor(n/13);temporal_source_sheet_certificate={temporal_certificate}:PASS")
    print(f"branch_profile=(sha256={side_context()['profile_digest']},nonempty={side_context()['nonempty']},rootwise_restoration=True)")
    print(f"geometric_branch_x_state_support=(matrix,components,counts_by_state,counts_by_branch)={geometry}")
    print(f"work_counts_(endpoint_intervals,mapped,active,q_active,weighted)={work_counts}")
    print(f"literal_guard_controls={CONTROL_TRIPLES}:PASS;Boolean_square_marginal_sha256={parent_hash}:PASS")
    print(f"rank_jump=(Boolean_square=4,branch_x_square={table_rank},branch_marginal={branch_marginal_rank})")
    print(f"full_F13branch_x_V4_x_F13residue_support={shape}")
    print(f"three_way_interaction=(rank={three_way_rank},support={three_way_shape})")
    print(f"fixed_relation_1_0_6=(branch_x_state_rank={fixed_rank},nonzero_entries={fixed_support}/52,sha256={digests[6]})")
    print(f"digests_(branch_profiles,gamma,table,spectrum,three_way,three_way_spectrum,fixed,branch_marginal)={digests}")
    print(f"proof={proof}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("reproducibility=no_randomness;no_elapsed_fields;normal_and_O_transcripts_must_match")
    print("commands=python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_branch_sheet_sidecar_probe_20260816.py;python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_branch_sheet_sidecar_probe_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
