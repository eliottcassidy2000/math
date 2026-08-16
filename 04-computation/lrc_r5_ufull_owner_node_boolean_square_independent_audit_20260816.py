#!/usr/bin/env python3
"""Clean-room audit of the r=5 owner-node Boolean-square refiner.

The submitted Boolean-square and common-owner scripts are not imported.  The
source service is rebuilt from the hash-pinned THM-2594 interval engine and
the actual endpoint factors are rebuilt from the independent THM-3514 atom
engine.  They are then coupled on a newly derived common coordinate by an
event sweep that lifts Q(13y) directly, rather than using the candidate's
inverse-branch jump routine.

The audit first applies MISTAKE-417's geometric gates: the new coordinate must
occupy more than one state over Q and its output matrix must have rank greater
than one.  Fourier support is interpreted only after those gates pass.  This
is one finite exact host computation, not an exact address, U_clock chronology,
uniform-row current, row exclusion, or LRC(14) theorem.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from math import gcd, lcm
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
ENDPOINT_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
SOURCE_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
ENDPOINT_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"

P = 13
Q = 7
V = 4
STATE_LABELS = ((0, 0), (0, 1), (1, 0), (1, 1))
STATE_SUBSETS = (
    frozenset((0,)),
    frozenset((0, 6)),
    frozenset((12,)),
    frozenset((6, 12)),
)
ROOT_SPINE = frozenset((0, 6, 12))
GRAY_PATH = (
    frozenset((0,)),
    frozenset((0, 6)),
    frozenset((6,)),
    frozenset((6, 12)),
    frozenset((12,)),
)
GRAY_MEASURE_NUMERATORS = (1, 12, 2, 12, 1)
CONTROL_TRIPLES = ((0, 0, 0), (0, 1, 6), (1, 0, 6), (6, 6, 12), (12, 12, 0))

JOINT_COORDINATE = 9684279613402457983920
JOINT_ORDER = 125895634974231953790960
PRIME = 755373809845391722745761
GENERATOR = 23
ROOT_OF_UNITY = 148035889
PRIME_FACTORS = (
    (2, 5), (3, 4), (5, 1), (7, 2), (11, 1), (13, 8),
    (53, 1), (61, 1), (131, 1), (313, 1),
)
SOURCE_PROFILE_SHA256 = "2de29f9be7fd16ceb4be5d15f7a71aa3d09f2907ec39f7af0b2017eadf3c18d2"
EXPECTED_STATE_WEIGHT_SEGMENTS = (3, 2, 3, 2)
EXPECTED_STATE_MEASURES = (JOINT_COORDINATE // 28,) * V
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_STATE_SEGMENTS = (300115, 300116, 300115, 300116)
EXPECTED_PARENT_GAMMA = "b5246eb2a69f35e4dac7dabbf26b1703f21ed22bf803061399ebbf766b9a073d"
EXPECTED_PARENT_ERASURE_GAMMA = (
    "20c83a7804c9437a7ccfaae5d3bf685fc1327c4defc385c1cf213f2f9643e258"
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
EXPECTED_BIT_MARGINALS = (
    (2, (26, 1, 1, 12, 12),
     "43450c53365624033ec2c098fb52812934d94de5daf1c3396d4e59b2e03bc680",
     "d48f20b5fe98a19bcfd20a55781dce268a3127a7fc76e8d6752f1c5f656ee896"),
    (2, (26, 1, 1, 12, 12),
     "daebb6ead49585eba1e118db3333af16179388a3d987bd626c2f151e3980db6d",
     "051da2cd94ed11bf9fb4fb708a1d15b8206a3a26729f35b5dc4f68609ef5efa6"),
)
EXPECTED_SEMANTIC_SHA256 = "b0996af3f1760b2118187490c93e0e01b322cb57fd6b25d3bf3688778b7e664c"

_CONTEXT = None


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    actual = lf_sha256(path)
    require(actual == expected_hash, (name, actual, expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def is_small_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return value == divisor
        divisor += 1
    return True


def split_field_certificate() -> tuple[object, ...]:
    factor_product = 1
    for prime, power in PRIME_FACTORS:
        require(is_small_prime(prime), ("nonprime factor", prime))
        factor_product *= prime**power
    require(factor_product == PRIME - 1 == 6 * JOINT_ORDER, "factorization")
    require(pow(GENERATOR, PRIME - 1, PRIME) == 1, "Fermat witness")
    lucas = []
    for prime, _power in PRIME_FACTORS:
        residue = pow(GENERATOR, (PRIME - 1) // prime, PRIME)
        require(gcd(residue - 1, PRIME) == 1, ("Lucas witness", prime, residue))
        lucas.append((prime, residue))
    require(pow(GENERATOR, 6, PRIME) == ROOT_OF_UNITY, "root derivation")
    require(pow(ROOT_OF_UNITY, JOINT_ORDER, PRIME) == 1, "root order upper bound")
    for prime, _power in PRIME_FACTORS:
        if JOINT_ORDER % prime == 0:
            require(pow(ROOT_OF_UNITY, JOINT_ORDER // prime, PRIME) != 1,
                    ("root exact order", prime))
    zeta = pow(ROOT_OF_UNITY, JOINT_ORDER // P, PRIME)
    require(pow(zeta, P, PRIME) == 1 and zeta != 1, "order 13 root")
    return zeta, tuple(lucas)


def profile_value(profile, point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def restrict_profile(starts, values, intervals, grid: int):
    pieces = []
    for interval_left, interval_right in intervals:
        index = bisect_right(starts, interval_left) - 1
        while True:
            profile_right = starts[index + 1] if index + 1 < len(starts) else grid
            left = max(interval_left, starts[index])
            right = min(interval_right, profile_right)
            if left < right and values[index]:
                pieces.append((left, right, values[index]))
            if profile_right >= interval_right:
                break
            index += 1
    return pieces


def build_source_profiles(M):
    e_intervals = M.build_set(M.PAT_E, M.ZELL)
    q_intervals = M.build_set(M.PAT_QB, M.ZELL)
    p2_starts, p2_values = M.weighted_fold(
        [(left, right, 1) for left, right in e_intervals], M.RPKT, M.T_DEN
    )
    f_pieces = restrict_profile(p2_starts, p2_values, q_intervals, M.T_DEN)
    u_whole = M.weighted_fold(f_pieces, M.DCOLL, M.T_DEN)
    v_whole = M.weighted_fold(
        [(left, right, 1) for left, right in e_intervals], M.DCOLL, M.T_DEN
    )
    source_u_raw = tuple(
        M.extract_window(u_whole[0], u_whole[1], root, P, M.T_DEN)
        for root in range(P)
    )
    source_v_raw = tuple(
        M.extract_window(v_whole[0], v_whole[1], root, P, M.T_DEN)
        for root in range(P)
    )
    source_digest = digest((source_u_raw, source_v_raw))
    require(source_digest == SOURCE_PROFILE_SHA256, source_digest)
    scale = JOINT_COORDINATE // M.T_DEN
    source_u = tuple(
        (tuple(point * scale for point in starts), tuple(values))
        for starts, values in source_u_raw
    )
    source_v = tuple(
        (tuple(point * scale for point in starts), tuple(values))
        for starts, values in source_v_raw
    )
    boundaries = tuple(sorted(
        {0, JOINT_COORDINATE}
        | {point for profile in source_u + source_v for point in profile[0]}
        | {JOINT_COORDINATE // 14, 13 * JOINT_COORDINATE // 14}
    ))
    return source_u, source_v, boundaries, source_digest


def partition_intervals(intervals, atoms):
    groups = [[] for _atom in atoms]
    atom_index = 0
    input_measure = 0
    output_measure = 0
    for interval_left, interval_right in intervals:
        input_measure += interval_right - interval_left
        while atoms[atom_index][3] <= interval_left:
            atom_index += 1
        scan = atom_index
        while scan < len(atoms) and atoms[scan][2] < interval_right:
            _sheet, _chamber, atom_left, atom_right = atoms[scan]
            left = max(interval_left, atom_left)
            right = min(interval_right, atom_right)
            if left < right:
                groups[scan].append((left, right))
                output_measure += right - left
            if atom_right >= interval_right:
                break
            scan += 1
    require(input_measure == output_measure, (input_measure, output_measure))
    return tuple(tuple(group) for group in groups)


def context():
    global _CONTEXT
    if _CONTEXT is None:
        M = load_module(SOURCE_PATH, "boolean_square_clean_source", SOURCE_SHA256)
        A = load_module(ENDPOINT_PATH, "boolean_square_clean_endpoint", ENDPOINT_SHA256)
        (
            word, endpoint_grid, endpoint_order, _old_prime, _old_root, _old_zeta,
            endpoint_q, _q_starts, _embeddings, _tabs, atom_intervals,
        ) = A.context()
        require(word == (1, 183, 27, 131, 53, 313, 13, 2197, 742586), word)
        require(word[A.M.OWNER] == P, "owner speed")
        require(endpoint_order == P**2 * endpoint_grid, "endpoint order")
        require(JOINT_COORDINATE == lcm(M.T_DEN, P * endpoint_grid), "joint coordinate")
        require(JOINT_ORDER == lcm(M.T_DEN, endpoint_order), "joint order")
        zeta, lucas = split_field_certificate()
        source_u, source_v, source_boundaries, source_digest = build_source_profiles(M)

        q_scale = JOINT_COORDINATE // (P * endpoint_grid)
        q_lift = tuple(
            ((branch * endpoint_grid + left) * q_scale,
             (branch * endpoint_grid + right) * q_scale)
            for branch in range(P)
            for left, right in endpoint_q
        )
        require(all(left < right for left, right in q_lift), "Q lift orientation")
        require(all(old[1] <= new[0] for old, new in zip(q_lift, q_lift[1:])),
                "Q lift overlap")
        q_lift_starts = tuple(left for left, _right in q_lift)
        endpoint_scale = JOINT_COORDINATE // endpoint_grid
        phase_scale = JOINT_ORDER // JOINT_COORDINATE
        require((endpoint_scale, phase_scale) == (20020, 13),
                (endpoint_scale, phase_scale))
        _CONTEXT = {
            "M": M, "A": A, "word": word, "endpoint_grid": endpoint_grid,
            "atom_intervals": atom_intervals, "endpoint_scale": endpoint_scale,
            "source_u": source_u, "source_v": source_v,
            "source_boundaries": source_boundaries, "source_digest": source_digest,
            "q_lift": q_lift, "q_lift_starts": q_lift_starts,
            "phase_scale": phase_scale, "zeta": zeta, "lucas": lucas,
        }
    return _CONTEXT


def chamber_of_segment(left: int, right: int) -> str:
    midpoint_twice = left + right
    if Q * midpoint_twice < JOINT_COORDINATE:
        return "left"
    if Q * midpoint_twice > 13 * JOINT_COORDINATE:
        return "right"
    raise RuntimeError(("OWNER segment outside cell zero", left, right))


def guard_safe(sheet: int, chamber: str, tau: int) -> int:
    midpoint_numerator = {"left": 1, "middle": 7, "right": 13}[chamber]
    doubled = (14 * (sheet + tau) + midpoint_numerator) % 182
    return int(min(doubled, 182 - doubled) > 26)


def source_values(point: int):
    ctx = context()
    u_values = tuple(profile_value(profile, point) for profile in ctx["source_u"])
    v_values = tuple(profile_value(profile, point) for profile in ctx["source_v"])
    require(all(not (u and v) for u, v in zip(u_values, v_values)),
            ("same-root source overlap", point))
    return u_values, v_values


def state_of_segment(left: int, right: int) -> int:
    chamber = chamber_of_segment(left, right)
    u_values, _v_values = source_values(left)
    support = frozenset(index for index, value in enumerate(u_values) if value)
    require(support <= ROOT_SPINE, ("root spine", left, support))
    mapping = {
        ("left", STATE_SUBSETS[0]): 0,
        ("left", STATE_SUBSETS[1]): 1,
        ("right", STATE_SUBSETS[2]): 2,
        ("right", STATE_SUBSETS[3]): 3,
    }
    require((chamber, support) in mapping, ("state map", left, right, chamber, support))
    return mapping[(chamber, support)]


def gray_path_and_owner_certificate():
    roots = tuple(sorted(ROOT_SPINE))
    ctx = context()
    boundaries = ctx["source_boundaries"]
    full_runs = []
    for left, right in zip(boundaries, boundaries[1:]):
        u_values, _v_values = source_values(left)
        support = frozenset(index for index, value in enumerate(u_values) if value)
        require(support and support < ROOT_SPINE, ("proper realized cut", left, support))
        if full_runs and full_runs[-1][0] == support and full_runs[-1][2] == left:
            full_runs[-1] = (support, full_runs[-1][1], right)
        else:
            full_runs.append((support, left, right))
    realized_path = tuple(run[0] for run in full_runs)
    realized_measures = tuple(right - left for _support, left, right in full_runs)
    expected_gray_measures = tuple(
        JOINT_COORDINATE * numerator // 28 for numerator in GRAY_MEASURE_NUMERATORS
    )
    require(realized_path == GRAY_PATH, realized_path)
    require(realized_measures == expected_gray_measures, realized_measures)
    require(frozenset((0, 12)) not in realized_path, realized_path)
    toggles = tuple(
        tuple(sorted(left ^ right)) for left, right in zip(realized_path, realized_path[1:])
    )
    require(toggles == ((6,), (0,), (12,), (6,)), toggles)

    cut_census = []
    for subset in realized_path:
        missing = one_way = two_way = 0
        for first_index, first in enumerate(roots):
            for second in roots[first_index + 1:]:
                forward = first in subset and second not in subset
                backward = second in subset and first not in subset
                if forward and backward:
                    two_way += 1
                elif forward or backward:
                    one_way += 1
                else:
                    missing += 1
        require((missing, one_way, two_way) == (1, 2, 0),
                (subset, missing, one_way, two_way))
        cut_census.append((tuple(sorted(subset)), missing, one_way, two_way))

    intervals = [[] for _state in range(V)]
    weight_segments = [0] * V
    for left, right in zip(boundaries, boundaries[1:]):
        midpoint_twice = left + right
        if not (Q * midpoint_twice < JOINT_COORDINATE
                or Q * midpoint_twice > 13 * JOINT_COORDINATE):
            continue
        state = state_of_segment(left, right)
        weight_segments[state] += 1
        if intervals[state] and intervals[state][-1][1] == left:
            intervals[state][-1] = (intervals[state][-1][0], right)
        else:
            intervals[state].append((left, right))
    components = tuple(len(group) for group in intervals)
    weight_segments_tuple = tuple(weight_segments)
    measures = tuple(sum(right - left for left, right in group) for group in intervals)
    require(weight_segments_tuple == EXPECTED_STATE_WEIGHT_SEGMENTS,
            weight_segments_tuple)
    require(components == (1, 1, 1, 1), components)
    require(measures == EXPECTED_STATE_MEASURES, measures)
    require(sum(measures) == JOINT_COORDINATE // Q, measures)
    complements = tuple(
        STATE_SUBSETS.index(ROOT_SPINE - subset) for subset in STATE_SUBSETS
    )
    require(complements == (3, 2, 1, 0), complements)
    require(all(complements[state] == (state ^ 3) for state in range(V)), complements)
    return (
        tuple(tuple(sorted(support)) for support in realized_path),
        realized_measures,
        toggles,
        tuple(cut_census),
        weight_segments_tuple,
        components,
        measures,
        complements,
    )


def endpoint_events(alpha: int, beta: int, literal_tau: int | None):
    ctx = context()
    A = ctx["A"]
    pattern = dict(A.M.PATTERN_E)
    if literal_tau is None:
        require(pattern.pop(A.M.GUARD) == "guard_safe", "guard deletion")
        tau = 0
    else:
        tau = literal_tau
    shift = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    intervals = A.M.fast_build_set(ctx["word"], ctx["endpoint_grid"], pattern, shift)
    groups = partition_intervals(intervals, ctx["atom_intervals"])
    events: dict[int, int] = {0: 0, JOINT_COORDINATE: 0}
    mapped = 0
    for atom_index, pieces in enumerate(groups):
        sheet, chamber, _atom_left, _atom_right = ctx["atom_intervals"][atom_index]
        require(chamber != "middle" or not pieces, ("middle atom", atom_index))
        bit = 1 << sheet
        for interval_left, interval_right in pieces:
            y_left = (P * interval_left - sheet * ctx["endpoint_grid"]) * ctx["endpoint_scale"]
            y_right = (P * interval_right - sheet * ctx["endpoint_grid"]) * ctx["endpoint_scale"]
            require(0 <= y_left < y_right <= JOINT_COORDINATE,
                    ("desheeted endpoint", atom_index, y_left, y_right))
            events[y_left] = events.get(y_left, 0) ^ bit
            events[y_right] = events.get(y_right, 0) ^ bit
            mapped += 1
    for boundary in ctx["source_boundaries"]:
        events.setdefault(boundary, 0)
    return events, len(intervals), mapped


def q_phase_jump(left: int, right: int) -> int:
    ctx = context()
    q_lift = ctx["q_lift"]
    q_starts = ctx["q_lift_starts"]
    index = bisect_right(q_starts, left) - 1
    if index < 0:
        index = 0
    elif q_lift[index][1] <= left:
        index += 1
    total = 0
    frequency = 57122 * ctx["phase_scale"]
    require(frequency == 742586, frequency)
    while index < len(q_lift) and q_lift[index][0] < right:
        q_left, q_right = q_lift[index]
        overlap_left = max(left, q_left)
        overlap_right = min(right, q_right)
        if overlap_left < overlap_right:
            total = (
                total
                + pow(ROOT_OF_UNITY, frequency * overlap_left % JOINT_ORDER, PRIME)
                - pow(ROOT_OF_UNITY, frequency * overlap_right % JOINT_ORDER, PRIME)
            ) % PRIME
        index += 1
    return total


def integrate_square(alpha: int, beta: int, literal_tau: int | None = None):
    events, interval_count, mapped = endpoint_events(alpha, beta, literal_tau)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    coupled = [[0] * V for _tau in tau_values]
    erased = [[0] * V for _tau in tau_values]
    diagonal = [[0] * V for _tau in tau_values]
    state_counts = [0] * V
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    for index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[index + 1]
        if mask == 0 or left == right:
            continue
        chamber = chamber_of_segment(left, right)
        state = state_of_segment(left, right)
        active_segments += 1
        state_counts[state] += 1
        jump = q_phase_jump(left, right)
        if jump == 0:
            continue
        q_active_segments += 1
        u_values, v_values = source_values(left)
        require(any(u_values) and any(v_values), ("source void", left, right))
        weighted_segments += 1
        for row_index, tau in enumerate(tau_values):
            if literal_tau is None:
                selected = tuple(
                    sheet for sheet in range(P)
                    if (mask >> sheet) & 1 and guard_safe(sheet, chamber, tau)
                )
            else:
                selected = tuple(sheet for sheet in range(P) if (mask >> sheet) & 1)
            left_value = sum(u_values[sheet] for sheet in selected)
            right_value = sum(v_values[sheet] for sheet in selected)
            diagonal_value = sum(u_values[sheet] * v_values[sheet] for sheet in selected)
            require(diagonal_value == 0,
                    ("same-root square", alpha, beta, tau, left, right))
            coupled[row_index][state] = (
                coupled[row_index][state] + left_value * right_value * jump
            ) % PRIME
            erased[row_index][state] = (
                erased[row_index][state] + len(selected) ** 2 * jump
            ) % PRIME
            diagonal[row_index][state] = (
                diagonal[row_index][state] + diagonal_value * jump
            ) % PRIME
    mask ^= events[positions[-1]]
    require(mask == 0, ("endpoint mask closure", alpha, beta, literal_tau, mask))
    require(all(value == 0 for row in diagonal for value in row), "diagonal bank")
    counts = (interval_count, mapped, active_segments, q_active_segments, weighted_segments)
    return (
        tuple(tuple(row) for row in coupled),
        tuple(tuple(row) for row in erased),
        tuple(tuple(row) for row in diagonal),
        counts,
        tuple(state_counts),
    )


def worker(alpha: int):
    zeta = context()["zeta"]
    coupled_rows = []
    erased_rows = []
    diagonal_rows = []
    counts = [0] * 5
    state_counts = [0] * V
    for beta in range(P):
        coupled, erased, diagonal, record, states = integrate_square(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        coupled_rows.append(tuple(
            tuple(phase * value % PRIME for value in row) for row in coupled
        ))
        erased_rows.append(tuple(
            tuple(phase * value % PRIME for value in row) for row in erased
        ))
        diagonal_rows.append(diagonal)
        counts = [left + right for left, right in zip(counts, record)]
        state_counts = [left + right for left, right in zip(state_counts, states)]
    return alpha, tuple(coupled_rows), tuple(erased_rows), tuple(diagonal_rows), tuple(counts), tuple(state_counts)


def inverse_state_table(gamma_states, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    table = [[0] * P for _state in range(V)]
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma_states[index]
                index += 1
                for residue in range(P):
                    phase = pow(zeta, -(alpha + tau * residue) % P, PRIME)
                    for state in range(V):
                        table[state][residue] = (
                            table[state][residue] + row[state] * phase
                        ) % PRIME
    require(index == P**3, index)
    return tuple(tuple(value * normalizer % PRIME for value in row) for row in table)


def walsh_fourier(table, zeta: int):
    return tuple(
        tuple(
            sum(
                table[state][residue]
                * (-1 if ((STATE_LABELS[state][0] & STATE_LABELS[channel][0])
                           ^ (STATE_LABELS[state][1] & STATE_LABELS[channel][1])) else 1)
                * pow(zeta, (-frequency * residue) % P, PRIME)
                for state in range(V) for residue in range(P)
            ) % PRIME
            for frequency in range(P)
        )
        for channel in range(V)
    )


def support_shape(spectrum, rows: int):
    dc = int(spectrum[0][0] != 0)
    row_axis = sum(spectrum[row][0] != 0 for row in range(1, rows))
    residue_axis = sum(spectrum[0][frequency] != 0 for frequency in range(1, P))
    mixed = sum(
        spectrum[row][frequency] != 0
        for row in range(1, rows) for frequency in range(1, P)
    )
    return dc + row_axis + residue_axis + mixed, dc, row_axis, residue_axis, mixed


def rank_mod(matrix) -> int:
    rows = [list(value % PRIME for value in row) for row in matrix]
    columns = len(rows[0]) if rows else 0
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, PRIME)
        rows[rank] = [value * inverse % PRIME for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [
                (value - factor * pivot_value) % PRIME
                for value, pivot_value in zip(rows[row], rows[rank])
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def interaction(table):
    inverse_rows = pow(len(table), -1, PRIME)
    inverse_columns = pow(P, -1, PRIME)
    row_means = tuple(sum(row) % PRIME * inverse_columns % PRIME for row in table)
    column_means = tuple(
        sum(table[row][column] for row in range(len(table))) % PRIME * inverse_rows % PRIME
        for column in range(P)
    )
    grand = sum(sum(row) for row in table) % PRIME
    grand = grand * inverse_rows % PRIME * inverse_columns % PRIME
    answer = tuple(
        tuple(
            (table[row][column] - row_means[row] - column_means[column] + grand) % PRIME
            for column in range(P)
        )
        for row in range(len(table))
    )
    require(all(sum(row) % PRIME == 0 for row in answer), "interaction rows")
    require(all(sum(answer[row][column] for row in range(len(answer))) % PRIME == 0
                for column in range(P)), "interaction columns")
    return answer


def binary_marginal(table, bit: int):
    return tuple(
        tuple(
            sum(table[state][residue] for state in range(V)
                if STATE_LABELS[state][bit] == value) % PRIME
            for residue in range(P)
        )
        for value in range(2)
    )


def binary_fourier(table, zeta: int):
    return tuple(
        tuple(
            sum(
                table[value][residue]
                * (-1 if channel and value else 1)
                * pow(zeta, (-frequency * residue) % P, PRIME)
                for value in range(2) for residue in range(P)
            ) % PRIME
            for frequency in range(P)
        )
        for channel in range(2)
    )


def main() -> None:
    ctx = context()
    (
        gray_path, gray_measures, gray_toggles, cut_record,
        state_weight_segments, state_components, state_measures, complements,
    ) = gray_path_and_owner_certificate()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    gamma = tuple(
        row for chunk in chunks for beta_rows in chunk[1] for row in beta_rows
    )
    erased_gamma = tuple(
        row for chunk in chunks for beta_rows in chunk[2] for row in beta_rows
    )
    diagonal_gamma = tuple(
        row for chunk in chunks for beta_rows in chunk[3] for row in beta_rows
    )
    require(len(gamma) == len(erased_gamma) == len(diagonal_gamma) == P**3,
            "gamma size")
    require(all(value == 0 for row in diagonal_gamma for value in row), "diagonal")

    parent_gamma = tuple((sum(row) % PRIME,) + (0,) * 6 for row in gamma)
    parent_erased = tuple((sum(row) % PRIME,) + (0,) * 6 for row in erased_gamma)
    parent_hashes = (digest(parent_gamma), digest(parent_erased))
    require(parent_hashes == (EXPECTED_PARENT_GAMMA, EXPECTED_PARENT_ERASURE_GAMMA),
            parent_hashes)

    literal_records = []
    for alpha, beta, tau in CONTROL_TRIPLES:
        direct, erased, diagonal, counts, states = integrate_square(alpha, beta, tau)
        phase = pow(ctx["zeta"], beta, PRIME)
        direct_row = tuple(phase * value % PRIME for value in direct[0])
        erased_row = tuple(phase * value % PRIME for value in erased[0])
        index = (alpha * P + beta) * P + tau
        require(direct_row == gamma[index], ("literal coupled", alpha, beta, tau))
        require(erased_row == erased_gamma[index], ("literal erased", alpha, beta, tau))
        require(all(value == 0 for value in diagonal[0]), "literal diagonal")
        literal_records.append(((alpha, beta, tau), counts, states))

    table = inverse_state_table(gamma, ctx["zeta"])
    erased_table = inverse_state_table(erased_gamma, ctx["zeta"])
    centred = interaction(table)
    spectrum = walsh_fourier(table, ctx["zeta"])
    erased_spectrum = walsh_fourier(erased_table, ctx["zeta"])
    centred_spectrum = walsh_fourier(centred, ctx["zeta"])
    ranks = (rank_mod(table), rank_mod(erased_table), rank_mod(centred))
    shapes = (
        support_shape(spectrum, V),
        support_shape(erased_spectrum, V),
        support_shape(centred_spectrum, V),
    )
    require(ranks == (4, 4, 3), ranks)
    require(shapes == ((52, 1, 3, 12, 36), (52, 1, 3, 12, 36),
                       (36, 0, 0, 0, 36)), shapes)

    relation_profile = tuple(table[state][6] for state in range(V))
    relation_walsh = tuple(
        sum(
            relation_profile[state]
            * (-1 if ((STATE_LABELS[state][0] & STATE_LABELS[channel][0])
                       ^ (STATE_LABELS[state][1] & STATE_LABELS[channel][1])) else 1)
            for state in range(V)
        ) % PRIME
        for channel in range(V)
    )
    require(relation_profile == EXPECTED_RELATION_PROFILE, relation_profile)
    require(relation_walsh == EXPECTED_RELATION_WALSH, relation_walsh)

    bit_records = []
    for bit in range(2):
        marginal = binary_marginal(table, bit)
        marginal_spectrum = binary_fourier(marginal, ctx["zeta"])
        record = (
            rank_mod(marginal), support_shape(marginal_spectrum, 2),
            digest(marginal), digest(marginal_spectrum),
        )
        bit_records.append(record)
    bit_records_tuple = tuple(bit_records)
    require(bit_records_tuple == EXPECTED_BIT_MARGINALS, bit_records_tuple)

    work_counts = tuple(sum(chunk[4][index] for chunk in chunks) for index in range(5))
    state_segments = tuple(sum(chunk[5][state] for chunk in chunks) for state in range(V))
    require(work_counts == EXPECTED_WORK_COUNTS, work_counts)
    require(state_segments == EXPECTED_STATE_SEGMENTS, state_segments)
    digests = (
        digest(gamma), digest(erased_gamma), digest(table), digest(erased_table),
        digest(spectrum), digest(erased_spectrum), digest(centred_spectrum),
    )
    require(digests == EXPECTED_DIGESTS, digests)

    proof = (
        "five_state_source_support_Gray_path_with_toggles_6_0_12_6",
        "source_support_0_12_never_occurs_and_OWNER_excludes_only_type_6",
        "each_realized_support_is_a_directed_cut_with_one_missing_edge",
        "four_states_have_equal_exact_measure_and_complement_XOR_3",
        "direct_Q13y_lift_and_independent_event_sweep",
        "five_literal_guard_controls_and_pointwise_zero_same_root",
        "state_marginal_recovers_corrected_delta_cell_parent",
        "rank4_and_centred_rank3_precede_Fourier_interpretation",
        "both_binary_marginals_rank2_and_all_36_mixed_modes_nonzero",
    )
    boundary = (
        "typed_source_support_refiner_on_one_r5_owner_node_host",
        "not_exact_address_or_U_clock_chronology_or_uniform_rows",
        "not_row_exclusion_or_LRC14",
    )
    semantic_surface = (
        SOURCE_SHA256, ENDPOINT_SHA256, ctx["source_digest"], JOINT_COORDINATE,
        JOINT_ORDER, PRIME, GENERATOR, ROOT_OF_UNITY, ctx["zeta"], ctx["lucas"],
        gray_path, gray_measures, gray_toggles, cut_record,
        state_weight_segments, state_components, state_measures, complements, work_counts,
        state_segments, tuple(literal_records), parent_hashes, digests, ranks,
        shapes, relation_profile, relation_walsh, bit_records_tuple, proof, boundary,
    )
    semantic = hashlib.sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("R5 U_FULL OWNER-NODE BOOLEAN SQUARE -- INDEPENDENT HOSTILE AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT_GENUINE_RANK_GT_1_REFINER;LRC14_OPEN")
    print(f"dependencies=(THM2594_engine={SOURCE_SHA256},THM3514_engine={ENDPOINT_SHA256})")
    print(f"joint_field=(coordinate={JOINT_COORDINATE},order={JOINT_ORDER},prime={PRIME},generator={GENERATOR},root={ROOT_OF_UNITY},zeta13={ctx['zeta']})")
    print(f"source_profile=(hash={ctx['source_digest']},root_spine={tuple(sorted(ROOT_SPINE))},same_root_pointwise_zero=True)")
    print(f"source_support_Gray_path=(states={gray_path},toggles={gray_toggles},measures={gray_measures},absent=(0,12))")
    print(f"realized_cut_certificate={cut_record}")
    print(f"owner_state_partition_over_Q=(visible_subsets={tuple(tuple(sorted(s)) for s in STATE_SUBSETS)},excluded_type=(6,),weight_segments={state_weight_segments},components={state_components},measures={state_measures},total={sum(state_measures)},complement_XOR=3)")
    print(f"work_counts={work_counts};state_segment_counts={state_segments}")
    print(f"literal_guard_controls={CONTROL_TRIPLES}:PASS;parent_delta_cell_marginal={parent_hashes}:PASS")
    print(f"coordinate_ranks=(coupled,source_erasure,ANOVA)={ranks}")
    print(f"spectral_shapes_(total,dc,V4axis,F13axis,mixed)={shapes}")
    print(f"fixed_relation_1_0_6=(profile={relation_profile},Walsh={relation_walsh})")
    print(f"bit_marginals_(component,multiplicity)={bit_records_tuple}")
    print(f"digests_(gamma,erased_gamma,table,erased_table,spectrum,erased_spectrum,ANOVA_spectrum)={digests}")
    print(f"proof={proof}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("reproducibility=no_candidate_imports;no_randomness;no_elapsed_fields;normal_and_O_transcripts_must_match")
    print("commands=python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_independent_audit_20260816.py;python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
