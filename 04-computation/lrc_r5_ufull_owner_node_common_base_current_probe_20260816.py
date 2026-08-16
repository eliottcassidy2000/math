#!/usr/bin/env python3
"""Exact owner-node common-base source/endpoint coupling at the r=5 window.

The source side is the literal THM-2471 root service

    U_u(y) = (P_(13^5) f)((y+u)/13),
    V_q(y) = (P_(13^5) e)((y+q)/13),

with f=1_Q P^2 e.  The endpoint side is the actual U_full Boolean set on
the same owner nodes w_u=(y+u)/13.  For every U_full character this script
forms, before integration,

    Q(13y) exp(2*pi*i*57122*y)
      [sum_u U_u(y) E_u(y)] [sum_q V_q(y) E_q(y)],

restores the shifted endpoint guard literally, and records the seven
THM-2594 cells.  A load-bearing hostile then checks the geometric collapse:
the endpoint owner condition becomes ``||y||<1/14``, so only cell zero can
occur.  Thus neither endpoint leg is a preintegrated AX/BY scalar, but the
resulting ``7 x 13`` array is a separable delta-cell lift rather than a
genuine two-coordinate cell interaction.

The two inherited grids do not embed in one another.  We therefore work on

    C = lcm(T_source, 13*T_endpoint)

and in a Lucas-certified split field for

    L = lcm(T_source, 13^2*T_endpoint).

This is a finite-exact one-base owner-node current candidate.  It is not an
exact THM-2334 address C(a;X,m), does not identify source time with arrival
time, supplies no U_clock transplant or chronology theorem, and proves no
row exclusion or LRC(14).
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from hashlib import sha256
import importlib.util
import json
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = (
    ROOT
    / "04-computation/lrc_r5_source_aligned_guard_atom_branch_sidecar_probe_20260816.py"
)
SOURCE_SHA256 = "22c5c748392817ccc36889a007c65bd5f44b26c10638df6f6aac48e917547f41"
TARGET_PATH = (
    ROOT
    / "04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py"
)
TARGET_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"

P = 13
Q = 7
SOURCE_PROFILE_SHA256 = "2de29f9be7fd16ceb4be5d15f7a71aa3d09f2907ec39f7af0b2017eadf3c18d2"
SOURCE_TOTAL_NUMERATOR = 168908463464745122312762880
JOINT_COORDINATE = 9684279613402457983920
JOINT_ORDER = 125895634974231953790960
JOINT_PRIME = 755373809845391722745761
JOINT_GENERATOR = 23
JOINT_ROOT = 148035889
JOINT_FACTORS = (
    (2, 4), (3, 3), (5, 1), (7, 2), (11, 1), (13, 8),
    (53, 1), (61, 1), (131, 1), (313, 1),
)
PRIME_FACTORS = (
    (2, 5), (3, 4), (5, 1), (7, 2), (11, 1), (13, 8),
    (53, 1), (61, 1), (131, 1), (313, 1),
)
LUCAS_WITNESSES = tuple((prime, JOINT_GENERATOR) for prime, _power in PRIME_FACTORS)
X_FREQUENCY = 4394
Y_FREQUENCY = P * X_FREQUENCY
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_GAMMA_SHA256 = "b5246eb2a69f35e4dac7dabbf26b1703f21ed22bf803061399ebbf766b9a073d"
EXPECTED_SOURCE_ERASURE_GAMMA_SHA256 = (
    "20c83a7804c9437a7ccfaae5d3bf685fc1327c4defc385c1cf213f2f9643e258"
)
EXPECTED_TABLE_SHA256 = (
    "d30230ee4e48146651d50af82514fa4d6f5acd08360303430e47e0a0fccf4411",
    "a59f045bf3e5f76d38b5e93b2ec4c98ba059dd5bc8d2a0807604771332814049",
)
EXPECTED_SPECTRUM_SHA256 = (
    "98f2fcde147575ecfb07d55246f4d4197d8e310c11b7daebf8213c7b7a54e60c",
    "4f83193667214bd64587c554e60ae26bd68e2978c85d81305b371dc7bc261f19",
    "070fbba00107a58f84cbd705d26fab2afc83da40a90f4deead13abd57bee46d8",
)
EXPECTED_RESIDUE_SHA256 = (
    "c30a063cc3d44f808b1da5b19ea00a6529ed6709f8c055fc8a59b6e7c07caf91",
    "9ddd65dc7b046e2fb6bd304929494bb15a2936779cabec34c932deb82929b6d6",
)
EXPECTED_ROLE_VALUES = (
    125385278409587426725290,
    657486478079327229022863,
    223272610175651920448188,
)
EXPECTED_SOURCE_ERASURE_VALUES = (
    471060960989539924924555,
    594285905723663416170205,
    632148865111268231500111,
)
EXPECTED_RELATION_PROFILE = (317699132065964946247468, 0, 0, 0, 0, 0, 0)
EXPECTED_SOURCE_ERASURE_RELATION_PROFILE = (
    454454155013190282848607, 0, 0, 0, 0, 0, 0,
)
EXPECTED_SHAPES = (
    (91, 1, 6, 12, 72),
    (91, 1, 6, 12, 72),
    (72, 0, 0, 0, 72),
)
EXPECTED_SEMANTIC_SHA256 = "98a27d4540648377c544d8e1b86c3dd3df7bb16d3431f6a9471f4844ba2e6b9f"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(body).hexdigest()


def digest_integers(values) -> str:
    return sha256(
        ",".join(str(value) for value in values).encode("ascii")
    ).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    require(lf_sha256(path) == expected_hash, (name, "source hash drift"))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "module loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_module(SOURCE_PATH, "owner_node_source_service", SOURCE_SHA256)
T = load_module(TARGET_PATH, "owner_node_endpoint_engine", TARGET_SHA256)
SM = S.M
EM = T.M
CTX = None


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def source_service_profiles():
    """Rebuild the two unsplit THM-2471 root-service families."""
    grid = S.GRID
    e_intervals = SM.build_set(SM.PAT_E, SM.ZELL)
    q_intervals = SM.build_set(SM.PAT_QB, SM.ZELL)
    f_pieces = S.B.build_f_pieces(e_intervals, q_intervals)
    u_whole = SM.weighted_fold(f_pieces, SM.DCOLL, grid)
    v_whole = SM.weighted_fold(
        [(left, right, 1) for left, right in e_intervals], SM.DCOLL, grid
    )
    source_u = tuple(
        SM.extract_window(u_whole[0], u_whole[1], root, P, grid)
        for root in range(P)
    )
    source_v = tuple(
        SM.extract_window(v_whole[0], v_whole[1], root, P, grid)
        for root in range(P)
    )
    profile_digest = digest_json((source_u, source_v))
    require(profile_digest == SOURCE_PROFILE_SHA256,
            ("source profile drift", profile_digest))

    boundaries = sorted(
        {0, grid}
        | {position for profile in source_u + source_v for position in profile[0]}
    )
    total = 0
    diagonal = 0
    type_measure = Counter()
    type_segments = Counter()
    type_examples = {}
    for left, right in zip(boundaries, boundaries[1:]):
        u_values = tuple(profile_value(profile, left) for profile in source_u)
        v_values = tuple(profile_value(profile, left) for profile in source_v)
        total += sum(u_values) * sum(v_values) * (right - left)
        diagonal_piece = sum(u * v for u, v in zip(u_values, v_values))
        diagonal += diagonal_piece * (right - left)
        require(diagonal_piece == 0, ("same-root source overlap", left, right))

        u_support = tuple(index for index, value in enumerate(u_values) if value)
        v_support = tuple(index for index, value in enumerate(v_values) if value)
        missing = one_way = two_way = 0
        for first in range(P):
            for second in range(first + 1, P):
                forward = first in u_support and second in v_support
                backward = second in u_support and first in v_support
                if forward and backward:
                    two_way += 1
                elif forward or backward:
                    one_way += 1
                else:
                    missing += 1
        graph_type = (len(u_support), len(v_support), missing, one_way, two_way)
        type_measure[graph_type] += right - left
        type_segments[graph_type] += 1
        type_examples.setdefault(graph_type, (u_support, v_support))

    require(total == SOURCE_TOTAL_NUMERATOR, ("source anchor numerator", total))
    require(diagonal == 0, ("source diagonal", diagonal))
    expected_graph_types = {
        (1, 12, 66, 12, 0): 42548128262640,
        (2, 11, 56, 22, 0): 255288769575840,
    }
    require(dict(type_measure) == expected_graph_types,
            ("source graph-type measure", type_measure))
    graph_record = tuple(
        (
            graph_type,
            type_segments[graph_type],
            type_measure[graph_type],
            type_examples[graph_type],
        )
        for graph_type in sorted(type_measure)
    )
    return source_u, source_v, boundaries, profile_digest, graph_record


def scale_profile(profile, numerator: int, denominator: int):
    require(numerator % denominator == 0, ("profile scale", numerator, denominator))
    scale = numerator // denominator
    starts, values = profile
    return tuple(position * scale for position in starts), tuple(values)


def context():
    global CTX
    if CTX is None:
        (
            word,
            endpoint_grid,
            endpoint_order,
            _old_prime,
            _old_root,
            _old_zeta,
            endpoint_q,
            _endpoint_q_starts,
            _embeddings,
            _tabs,
            atom_intervals,
        ) = T.context()
        source_grid = S.GRID
        require(word == (1, 183, 27, 131, 53, 313, 13, 2197, 742586),
                ("U_full word", word))
        require(endpoint_order == P * P * endpoint_grid,
                ("endpoint order", endpoint_order, endpoint_grid))
        require(JOINT_COORDINATE == lcm(source_grid, P * endpoint_grid),
                "joint coordinate drift")
        require(JOINT_ORDER == lcm(source_grid, endpoint_order),
                "joint order drift")
        require(JOINT_ORDER % JOINT_COORDINATE == 0,
                "coordinate does not divide field order")
        EM.verify_lucas_certificate(
            JOINT_PRIME, PRIME_FACTORS, LUCAS_WITNESSES, "joint owner-node prime"
        )
        EM.verify_embedding(
            JOINT_PRIME, JOINT_ROOT, JOINT_ORDER, JOINT_FACTORS,
            "joint owner-node root",
        )
        require(pow(JOINT_GENERATOR, 6, JOINT_PRIME) == JOINT_ROOT,
                "joint generator/root mismatch")

        coordinate_scale = JOINT_ORDER // JOINT_COORDINATE
        endpoint_scale = JOINT_COORDINATE // endpoint_grid
        source_scale = JOINT_COORDINATE // source_grid
        require((coordinate_scale, endpoint_scale, source_scale) == (13, 20020, 32515379),
                ("joint scales", coordinate_scale, endpoint_scale, source_scale))
        q_intervals = tuple(
            (left * endpoint_scale, right * endpoint_scale)
            for left, right in endpoint_q
        )
        q_starts = tuple(left for left, _right in q_intervals)
        q_left_powers = tuple(
            pow(
                JOINT_ROOT,
                X_FREQUENCY * left * coordinate_scale % JOINT_ORDER,
                JOINT_PRIME,
            )
            for left, _right in q_intervals
        )
        q_right_powers = tuple(
            pow(
                JOINT_ROOT,
                X_FREQUENCY * right * coordinate_scale % JOINT_ORDER,
                JOINT_PRIME,
            )
            for _left, right in q_intervals
        )

        source_u_raw, source_v_raw, source_boundaries_raw, profile_digest, graph_record = (
            source_service_profiles()
        )
        source_u = tuple(
            scale_profile(profile, JOINT_COORDINATE, source_grid)
            for profile in source_u_raw
        )
        source_v = tuple(
            scale_profile(profile, JOINT_COORDINATE, source_grid)
            for profile in source_v_raw
        )
        source_boundaries = tuple(
            position * source_scale for position in source_boundaries_raw
        )
        zeta = pow(JOINT_ROOT, JOINT_ORDER // P, JOINT_PRIME)
        eta = pow(JOINT_ROOT, JOINT_ORDER // Q, JOINT_PRIME)
        require(pow(zeta, P, JOINT_PRIME) == 1 and zeta != 1,
                "joint order-13 root")
        require(pow(eta, Q, JOINT_PRIME) == 1 and eta != 1,
                "joint order-7 root")
        CTX = {
            "word": word,
            "endpoint_grid": endpoint_grid,
            "atom_intervals": atom_intervals,
            "endpoint_scale": endpoint_scale,
            "coordinate_scale": coordinate_scale,
            "q_intervals": q_intervals,
            "q_starts": q_starts,
            "q_left_powers": q_left_powers,
            "q_right_powers": q_right_powers,
            "source_u": source_u,
            "source_v": source_v,
            "source_boundaries": source_boundaries,
            "profile_digest": profile_digest,
            "graph_record": graph_record,
            "zeta": zeta,
            "eta": eta,
        }
    return CTX


def cell_of_segment(left: int, right: int) -> int:
    band = (Q * (left + right)) // JOINT_COORDINATE
    require(0 <= band < 2 * Q, ("cell band", left, right, band))
    if band in (0, 2 * Q - 1):
        return 0
    return (band + 1) // 2


def chamber_of_segment(left: int, right: int) -> str:
    midpoint_twice = left + right
    if Q * midpoint_twice < 2 * JOINT_COORDINATE:
        return "left"
    if Q * midpoint_twice < 12 * JOINT_COORDINATE:
        return "middle"
    return "right"


@lru_cache(maxsize=400000)
def q_endpoint_jump(left: int, right: int) -> int:
    """Endpoint numerator for Q(13y) exp(2*pi*i*57122*y)."""
    ctx = context()
    scaled_left = P * left
    branch, start = divmod(scaled_left, JOINT_COORDINATE)
    span = P * (right - left)
    stop = start + span
    require(0 <= branch < P and 0 < span <= JOINT_COORDINATE,
            ("inverse branch", left, right, branch, span))
    require(stop <= JOINT_COORDINATE,
            ("segment crosses inverse branch", left, right, start, stop))

    q_intervals = ctx["q_intervals"]
    q_starts = ctx["q_starts"]
    index = bisect_right(q_starts, start) - 1
    if index < 0:
        index = 0
    elif q_intervals[index][1] <= start:
        index += 1
    total = 0
    scale = ctx["coordinate_scale"]
    while index < len(q_intervals):
        q_left, q_right = q_intervals[index]
        if q_left >= stop:
            break
        overlap_left = max(start, q_left)
        overlap_right = min(stop, q_right)
        if overlap_left < overlap_right:
            if overlap_left == q_left:
                value_left = ctx["q_left_powers"][index]
            else:
                value_left = pow(
                    JOINT_ROOT,
                    X_FREQUENCY * overlap_left * scale % JOINT_ORDER,
                    JOINT_PRIME,
                )
            if overlap_right == q_right:
                value_right = ctx["q_right_powers"][index]
            else:
                value_right = pow(
                    JOINT_ROOT,
                    X_FREQUENCY * overlap_right * scale % JOINT_ORDER,
                    JOINT_PRIME,
                )
            total = (total + value_left - value_right) % JOINT_PRIME
        index += 1

    # x=13y-branch and 4394*branch is integral, so the branch phase is one.
    require(Y_FREQUENCY == P * X_FREQUENCY, "frequency descent")
    return total


def endpoint_events(alpha: int, beta: int, literal_tau: int | None):
    ctx = context()
    word = ctx["word"]
    endpoint_grid = ctx["endpoint_grid"]
    if literal_tau is None:
        pattern = dict(EM.PATTERN_E)
        require(pattern.pop(EM.GUARD) == "guard_safe",
                "removed endpoint factor is not the guard")
        tau = 0
    else:
        pattern = EM.PATTERN_E
        tau = literal_tau
    shift = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    intervals = EM.fast_build_set(word, endpoint_grid, pattern, shift)
    groups = T.partition_two_pointer(intervals, ctx["atom_intervals"])
    events: dict[int, int] = {}
    mapped = 0
    endpoint_scale = ctx["endpoint_scale"]
    for atom_index, pieces in enumerate(groups):
        sheet, _chamber = T.ATOMS[atom_index]
        bit = 1 << sheet
        for interval_left, interval_right in pieces:
            y_left = (P * interval_left - sheet * endpoint_grid) * endpoint_scale
            y_right = (P * interval_right - sheet * endpoint_grid) * endpoint_scale
            require(0 <= y_left < y_right <= JOINT_COORDINATE,
                    ("desheeted endpoint interval", atom_index, y_left, y_right))
            events[y_left] = events.get(y_left, 0) ^ bit
            events[y_right] = events.get(y_right, 0) ^ bit
            mapped += 1

    boundaries = {0, JOINT_COORDINATE}
    boundaries.update(index * (JOINT_COORDINATE // P) for index in range(P + 1))
    boundaries.update(
        index * (JOINT_COORDINATE // (2 * Q)) for index in range(2 * Q + 1)
    )
    boundaries.update(ctx["source_boundaries"])
    for boundary in boundaries:
        events.setdefault(boundary, 0)
    return events, len(intervals), mapped


def integrate_profile(alpha: int, beta: int, literal_tau: int | None = None):
    """Integrate all guard shifts, or one independently guarded control."""
    ctx = context()
    events, interval_count, mapped = endpoint_events(alpha, beta, literal_tau)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    coupled = [[0 for _ell in range(Q)] for _tau in tau_values]
    source_erasure = [[0 for _ell in range(Q)] for _tau in tau_values]
    diagonal = [[0 for _ell in range(Q)] for _tau in tau_values]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        chamber = chamber_of_segment(left, right)
        require(chamber != "middle",
                ("owner support leaked into middle chamber", alpha, beta, left, right))
        active_segments += 1
        jump = q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        u_values = tuple(profile_value(profile, left) for profile in ctx["source_u"])
        v_values = tuple(profile_value(profile, left) for profile in ctx["source_v"])
        if not any(u_values) or not any(v_values):
            continue
        weighted_segments += 1
        ell = cell_of_segment(left, right)
        # OWNER has speed 13 and is an ``in`` factor.  On the desheeted
        # branch 13t=y+sheet, hence OWNER is exactly ||y||<1/14=cell_0.
        # This is a geometric zero over Q, not merely a zero in JOINT_PRIME.
        require(ell == 0,
                ("endpoint owner escaped cell zero", alpha, beta, left, right, ell))
        for row_index, tau in enumerate(tau_values):
            if literal_tau is None:
                selected = tuple(
                    sheet
                    for sheet in range(P)
                    if (mask >> sheet) & 1 and T.safe(chamber, sheet + tau)
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
                    ("same-root current fired", alpha, beta, tau, left, right))
            coupled[row_index][ell] = (
                coupled[row_index][ell] + left_value * right_value * jump
            ) % JOINT_PRIME
            source_erasure[row_index][ell] = (
                source_erasure[row_index][ell] + len(selected) ** 2 * jump
            ) % JOINT_PRIME
            diagonal[row_index][ell] = (
                diagonal[row_index][ell] + diagonal_value * jump
            ) % JOINT_PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("endpoint mask did not close", alpha, beta, literal_tau, mask))
    require(all(value == 0 for row in diagonal for value in row),
            ("same-root row", alpha, beta, literal_tau))
    counts = (
        interval_count, mapped, active_segments, q_active_segments, weighted_segments
    )
    return (
        tuple(tuple(row) for row in coupled),
        tuple(tuple(row) for row in source_erasure),
        tuple(tuple(row) for row in diagonal),
        counts,
    )


def worker(alpha: int):
    zeta = context()["zeta"]
    coupled_rows = []
    source_erasure_rows = []
    diagonal_rows = []
    counts = [0] * 5
    for beta in range(P):
        coupled, source_erasure, diagonal, record = integrate_profile(alpha, beta)
        phase = pow(zeta, beta, JOINT_PRIME)
        coupled_rows.append(tuple(
            tuple(phase * value % JOINT_PRIME for value in row)
            for row in coupled
        ))
        source_erasure_rows.append(tuple(
            tuple(phase * value % JOINT_PRIME for value in row)
            for row in source_erasure
        ))
        diagonal_rows.append(diagonal)
        counts = [left + right for left, right in zip(counts, record)]
    return (
        alpha,
        tuple(coupled_rows),
        tuple(source_erasure_rows),
        tuple(diagonal_rows),
        tuple(counts),
    )


def inverse_cell_table(gamma_cells, zeta: int):
    normalizer = pow(P**3, -1, JOINT_PRIME)
    table = [[0 for _t in range(P)] for _ell in range(Q)]
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma_cells[index]
                index += 1
                for relation_t in range(P):
                    phase = pow(zeta, -(alpha + tau * relation_t) % P, JOINT_PRIME)
                    for ell in range(Q):
                        table[ell][relation_t] = (
                            table[ell][relation_t] + row[ell] * phase
                        ) % JOINT_PRIME
    require(index == P**3, ("character bank size", index))
    return tuple(
        tuple(value * normalizer % JOINT_PRIME for value in row)
        for row in table
    )


def fourier_2d(matrix, eta: int, zeta: int):
    return tuple(
        tuple(
            sum(
                matrix[ell][relation_t]
                * pow(eta, -h * ell % Q, JOINT_PRIME)
                * pow(zeta, -k * relation_t % P, JOINT_PRIME)
                for ell in range(Q) for relation_t in range(P)
            ) % JOINT_PRIME
            for k in range(P)
        )
        for h in range(Q)
    )


def support_shape(spectrum) -> tuple[int, int, int, int, int]:
    dc = int(spectrum[0][0] != 0)
    cell_axis = sum(spectrum[h][0] != 0 for h in range(1, Q))
    residue_axis = sum(spectrum[0][k] != 0 for k in range(1, P))
    mixed = sum(
        spectrum[h][k] != 0 for h in range(1, Q) for k in range(1, P)
    )
    return dc + cell_axis + residue_axis + mixed, dc, cell_axis, residue_axis, mixed


def rank_mod(matrix) -> int:
    rows = [[value % JOINT_PRIME for value in row] for row in matrix]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column]), None
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, JOINT_PRIME)
        rows[rank] = [value * inverse % JOINT_PRIME for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (value - factor * pivot_value) % JOINT_PRIME
                for value, pivot_value in zip(rows[row], rows[rank])
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def matrix_interaction(matrix):
    inv_q = pow(Q, -1, JOINT_PRIME)
    inv_p = pow(P, -1, JOINT_PRIME)
    inv_qp = pow(Q * P, -1, JOINT_PRIME)
    row_sums = tuple(sum(row) % JOINT_PRIME for row in matrix)
    column_sums = tuple(
        sum(matrix[ell][relation_t] for ell in range(Q)) % JOINT_PRIME
        for relation_t in range(P)
    )
    grand = sum(row_sums) % JOINT_PRIME
    answer = tuple(
        tuple(
            (
                matrix[ell][relation_t]
                - row_sums[ell] * inv_p
                - column_sums[relation_t] * inv_q
                + grand * inv_qp
            ) % JOINT_PRIME
            for relation_t in range(P)
        )
        for ell in range(Q)
    )
    require(all(sum(row) % JOINT_PRIME == 0 for row in answer),
            "interaction row sums")
    require(all(
        sum(answer[ell][relation_t] for ell in range(Q)) % JOINT_PRIME == 0
        for relation_t in range(P)
    ), "interaction column sums")
    return answer


def inverse_value(table, relation_t: int) -> int:
    return sum(table[ell][relation_t] for ell in range(Q)) % JOINT_PRIME


def main() -> None:
    ctx = context()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "worker order")
    coupled_cells = tuple(
        row
        for _alpha, coupled_rows, _source_erasure_rows, _diagonal_rows, _counts in chunks
        for beta_rows in coupled_rows
        for row in beta_rows
    )
    source_erasure_cells = tuple(
        row
        for _alpha, _coupled_rows, source_erasure_rows, _diagonal_rows, _counts in chunks
        for beta_rows in source_erasure_rows
        for row in beta_rows
    )
    diagonal_cells = tuple(
        row
        for _alpha, _coupled_rows, _source_erasure_rows, diagonal_rows, _counts in chunks
        for beta_rows in diagonal_rows
        for row in beta_rows
    )
    require(len(coupled_cells) == len(source_erasure_cells) == len(diagonal_cells) == P**3,
            "gamma cell bank size")
    require(all(value == 0 for row in diagonal_cells for value in row),
            "same-root gamma bank")

    direct_controls = ((0, 0, 0), (0, 1, 6), (1, 0, 6), (6, 6, 12), (12, 12, 0))
    direct_record = []
    for alpha, beta, tau in direct_controls:
        direct_coupled, direct_source_erasure, direct_diagonal, counts = integrate_profile(
            alpha, beta, tau
        )
        phase = pow(ctx["zeta"], beta, JOINT_PRIME)
        direct_coupled_row = tuple(
            phase * value % JOINT_PRIME for value in direct_coupled[0]
        )
        direct_source_erasure_row = tuple(
            phase * value % JOINT_PRIME for value in direct_source_erasure[0]
        )
        index = (alpha * P + beta) * P + tau
        require(direct_coupled_row == coupled_cells[index],
                ("literal coupled guard restoration", alpha, beta, tau))
        require(direct_source_erasure_row == source_erasure_cells[index],
                ("literal source-erasure guard restoration", alpha, beta, tau))
        require(all(value == 0 for value in direct_diagonal[0]),
                ("literal same-root hostile", alpha, beta, tau))
        direct_record.append(
            ((alpha, beta, tau), direct_coupled_row, direct_source_erasure_row, counts)
        )

    zeta = ctx["zeta"]
    eta = ctx["eta"]
    coupled_table = inverse_cell_table(coupled_cells, zeta)
    source_erasure_table = inverse_cell_table(source_erasure_cells, zeta)
    coupled_spectrum = fourier_2d(coupled_table, eta, zeta)
    source_erasure_spectrum = fourier_2d(source_erasure_table, eta, zeta)
    interaction = matrix_interaction(coupled_table)
    interaction_spectrum = fourier_2d(interaction, eta, zeta)
    require(all(row[ell] == 0 for row in coupled_cells for ell in range(1, Q)),
            "character bank escaped the owner cell")
    require(all(row[ell] == 0
                for row in source_erasure_cells for ell in range(1, Q)),
            "source-erasure bank escaped the owner cell")
    require(all(value == 0 for row in coupled_table[1:] for value in row),
            "inverse table escaped the owner cell")
    require(all(value == 0
                for row in source_erasure_table[1:] for value in row),
            "source-erasure inverse escaped the owner cell")
    coordinate_ranks = (
        rank_mod(coupled_table),
        rank_mod(source_erasure_table),
        rank_mod(interaction),
    )
    require(coordinate_ranks == (1, 1, 1),
            ("separable coordinate ranks", coordinate_ranks))
    residue_profile = coupled_table[0]
    residue_spectrum = tuple(
        sum(
            residue_profile[relation_t]
            * pow(zeta, -k * relation_t % P, JOINT_PRIME)
            for relation_t in range(P)
        ) % JOINT_PRIME
        for k in range(P)
    )
    require(all(value != 0 for value in residue_spectrum),
            ("one-dimensional residue spectrum", residue_spectrum))
    require(all(
        coupled_spectrum[h][k] == residue_spectrum[k]
        for h in range(Q) for k in range(P)
    ), "delta-cell Fourier factorization")
    residue_mean = sum(residue_profile) * pow(P, -1, JOINT_PRIME) % JOINT_PRIME
    cell_mean = pow(Q, -1, JOINT_PRIME)
    require(all(
        interaction[ell][relation_t]
        == (
            ((int(ell == 0) - cell_mean) % JOINT_PRIME)
            * ((residue_profile[relation_t] - residue_mean) % JOINT_PRIME)
        ) % JOINT_PRIME
        for ell in range(Q) for relation_t in range(P)
    ), "ANOVA is not the expected rank-one outer product")
    shapes = (
        support_shape(coupled_spectrum),
        support_shape(source_erasure_spectrum),
        support_shape(interaction_spectrum),
    )

    relation_t = 6
    relation_profile = tuple(
        coupled_table[ell][relation_t] for ell in range(Q)
    )
    relation_spectrum = tuple(
        sum(
            relation_profile[ell] * pow(eta, -h * ell % Q, JOINT_PRIME)
            for ell in range(Q)
        ) % JOINT_PRIME
        for h in range(Q)
    )
    source_erasure_relation_profile = tuple(
        source_erasure_table[ell][relation_t] for ell in range(Q)
    )
    source_erasure_relation_spectrum = tuple(
        sum(
            source_erasure_relation_profile[ell]
            * pow(eta, -h * ell % Q, JOINT_PRIME)
            for ell in range(Q)
        ) % JOINT_PRIME
        for h in range(Q)
    )

    role_values = (
        inverse_value(coupled_table, 1),
        inverse_value(coupled_table, 0),
    )
    bridge = (role_values[0] - role_values[1]) % JOINT_PRIME
    source_erasure_values = (
        inverse_value(source_erasure_table, 1),
        inverse_value(source_erasure_table, 0),
    )
    source_erasure_bridge = (
        source_erasure_values[0] - source_erasure_values[1]
    ) % JOINT_PRIME
    require(bridge != 0, ("coupled bridge vanished", role_values))
    require(all(value != 0 for value in relation_spectrum),
            ("fixed-relation F7 mode vanished", relation_spectrum))

    counts = tuple(
        sum(chunk[4][index] for chunk in chunks) for index in range(5)
    )
    gamma_digest = digest_json(coupled_cells)
    source_erasure_gamma_digest = digest_json(source_erasure_cells)
    table_digests = (
        digest_json(coupled_table), digest_json(source_erasure_table)
    )
    spectrum_digests = (
        digest_json(coupled_spectrum), digest_json(source_erasure_spectrum),
        digest_json(interaction_spectrum),
    )
    residue_digests = (
        digest_json(residue_profile), digest_json(residue_spectrum)
    )
    require(counts == EXPECTED_WORK_COUNTS, ("work counts", counts))
    require(gamma_digest == EXPECTED_GAMMA_SHA256,
            ("coupled gamma digest", gamma_digest))
    require(source_erasure_gamma_digest == EXPECTED_SOURCE_ERASURE_GAMMA_SHA256,
            ("source-erasure gamma digest", source_erasure_gamma_digest))
    require(table_digests == EXPECTED_TABLE_SHA256,
            ("table digests", table_digests))
    require(spectrum_digests == EXPECTED_SPECTRUM_SHA256,
            ("spectrum digests", spectrum_digests))
    require(residue_digests == EXPECTED_RESIDUE_SHA256,
            ("residue digests", residue_digests))
    require((*role_values, bridge) == EXPECTED_ROLE_VALUES,
            ("coupled inverse values", role_values, bridge))
    require((*source_erasure_values, source_erasure_bridge)
            == EXPECTED_SOURCE_ERASURE_VALUES,
            ("source-erasure inverse values", source_erasure_values,
             source_erasure_bridge))
    require(relation_profile == EXPECTED_RELATION_PROFILE,
            ("fixed relation profile", relation_profile))
    require(source_erasure_relation_profile
            == EXPECTED_SOURCE_ERASURE_RELATION_PROFILE,
            ("source-erasure fixed relation profile",
             source_erasure_relation_profile))
    require(relation_spectrum == (EXPECTED_RELATION_PROFILE[0],) * Q,
            ("fixed relation F7 spectrum", relation_spectrum))
    require(source_erasure_relation_spectrum
            == (EXPECTED_SOURCE_ERASURE_RELATION_PROFILE[0],) * Q,
            ("source-erasure fixed relation F7 spectrum",
             source_erasure_relation_spectrum))
    require(shapes == EXPECTED_SHAPES, ("spectral shapes", shapes))
    record = (
        SOURCE_SHA256,
        TARGET_SHA256,
        ctx["word"],
        S.GRID,
        ctx["endpoint_grid"],
        JOINT_COORDINATE,
        JOINT_ORDER,
        JOINT_PRIME,
        JOINT_GENERATOR,
        JOINT_ROOT,
        eta,
        zeta,
        X_FREQUENCY,
        Y_FREQUENCY,
        ctx["profile_digest"],
        ctx["graph_record"],
        counts,
        tuple(direct_record),
        gamma_digest,
        source_erasure_gamma_digest,
        table_digests,
        spectrum_digests,
        shapes,
        coordinate_ranks,
        residue_digests,
        role_values,
        bridge,
        source_erasure_values,
        source_erasure_bridge,
        relation_t,
        relation_profile,
        relation_spectrum,
        source_erasure_relation_profile,
        source_erasure_relation_spectrum,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 U_full owner-node common-base source/endpoint current probe ==")
    print(f"dependencies=((source,{SOURCE_SHA256}),(endpoint,{TARGET_SHA256}))")
    print(f"grids=(source={S.GRID},endpoint={ctx['endpoint_grid']},joint_coordinate={JOINT_COORDINATE},joint_order={JOINT_ORDER})")
    print(f"field=(prime={JOINT_PRIME},generator={JOINT_GENERATOR},root={JOINT_ROOT},eta7={eta},zeta13={zeta})")
    print(f"frequency_descent=(x_frequency={X_FREQUENCY},y_frequency={Y_FREQUENCY})")
    print(f"source_profile=(sha256={ctx['profile_digest']},total_numerator={SOURCE_TOTAL_NUMERATOR},same_root=0)")
    print(f"source_directed_pair_types_(|U|,|V|,missing,one_way,two_way)={ctx['graph_record']}")
    print(f"work_counts=(endpoint_intervals,mapped_intervals,active_segments,q_active_segments,source_weighted_segments)={counts}")
    print(f"literal_guard_restoration_controls={direct_controls}: PASS")
    print("same_root_source_current=pointwise_zero_before_integration: PASS")
    print(f"gamma_sha256=(coupled={gamma_digest},source_erasure={source_erasure_gamma_digest})")
    print(f"inverse_roles=(q_H={role_values[0]},q_q5={role_values[1]},bridge={bridge})")
    print(f"source_erasure_inverse=(q_H={source_erasure_values[0]},q_q5={source_erasure_values[1]},bridge={source_erasure_bridge})")
    print(f"spectral_shapes_(total,dc,F7axis,F13axis,mixed)=(coupled={shapes[0]},source_erasure={shapes[1]},coupled_ANOVA={shapes[2]})")
    print(f"cell_geometry=(owner_implies_cell_0,character_support_by_cell={(P**3, 0, 0, 0, 0, 0, 0)},coordinate_ranks={coordinate_ranks})")
    print(f"separable_factorization=table=delta_cell0*residue_profile;ANOVA=(delta_cell0-1/7)*(residue_profile-mean);residue_sha256={residue_digests}")
    print(f"fixed_relation_class=(1,0,{relation_t});seven_cell_profile={relation_profile}")
    print(f"fixed_relation_F7_spectrum={relation_spectrum}")
    print(f"source_erasure_fixed_relation_profile={source_erasure_relation_profile}")
    print(f"source_erasure_fixed_relation_F7_spectrum={source_erasure_relation_spectrum}")
    print(f"table_sha256=(coupled={table_digests[0]},source_erasure={table_digests[1]})")
    print(f"spectrum_sha256=(coupled={spectrum_digests[0]},source_erasure={spectrum_digests[1]},coupled_ANOVA={spectrum_digests[2]})")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT one-common-base owner-node current candidate with nonzero one-dimensional residue signal; the 7x13 support is a rank-one delta-cell lift, not genuine cell/residue mixing")
    print("scope=no exact C(a;X,m),no arrival/source-time identification,no U_clock chronology,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
