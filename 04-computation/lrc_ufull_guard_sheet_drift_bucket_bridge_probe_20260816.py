#!/usr/bin/env python3
"""Exact 39-atom/117-drift decomposition of the frozen U_full bridge.

The refined THM-3479 endpoint bank varies over

    ell=(tau,-alpha,-beta,0,0,0,0,alpha,beta) in F_13^9.

Only the H ``guard_safe`` factor depends on tau.  This companion removes
that one factor, cuts every remaining E interval by the fixed decomposition

    s=7*a+r,  a in F_13,  r in [0,1), [1,6), [6,7),

and computes AX and BY before restoring the guard multiplier g_C(a+tau).
Thus the full 13^3 character bank is reconstructed from 39 actual guard
atoms.  Expanding AX*BY then groups the q_H and q_q5 inverse transforms into
the 117 kernel types (C,D,d), d=a_R-a_L, with the common left sheet retained
through its Fourier phase.

This is an exact endpoint-factor decomposition and an implementation of the
previous guard-sheet address law.  It is not a common-ancestry support
relation: the Cartesian product of the two endpoint atom tables remains in
use.  No root/owner/word/source sheet, horizon, chronology, physical current,
grouped C(a;X,m), all-unit projector, row exclusion, or LRC(14) conclusion is
constructed.
"""

from __future__ import annotations

import ast
from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_ufull_guard_sheet_drift_bucket_bridge_probe_20260816.py"
OUTPUT = "05-knowledge/results/lrc_ufull_guard_sheet_drift_bucket_bridge_probe_20260816.out"

BRIDGE_PATH = ROOT / "04-computation/lrc_half_twist_relation_current_bridge_thm3479.py"
BRIDGE_SHA256 = "ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b"
GUARD_PATH = ROOT / "04-computation/lrc_ufull_guard_sheet_joint_drift_address_probe_20260816.py"
GUARD_SHA256 = "79b78637b2cc0ff54051fde02a6651ef10c8694a8d7a865ae403696370125179"
REFINED_PATH = ROOT / "04-computation/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.py"
REFINED_SHA256 = "ee2105742abee578a9c41ff7ec954a07ada324fccc2c643429e7ac6e6e6f8fc2"
REFINED_OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.out"
REFINED_OUTPUT_SHA256 = "10a98351cc59615a5b6d2b8f555e0936d1a39566d9906127edc2b0fbc3918e73"

P = 13
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)
EXPECTED_GAMMA_SHA256 = "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
EXPECTED_VALUES = (320618948602619577408, 503604956476841920373)
EXPECTED_BRIDGE = 389266878372286537904
EXPECTED_SEMANTIC_SHA256 = "0b31a992ba23cecd05f28ae353133531f41cc6d84a4a935c34a12d77fd3db590"

CHAMBERS = (
    ("left", 0, 1, frozenset((12, 0, 1))),
    ("middle", 1, 6, frozenset((11, 12, 0, 1))),
    ("right", 6, 7, frozenset((11, 12, 0))),
)
CHAMBER_NAMES = tuple(row[0] for row in CHAMBERS)
CHAMBER_DANGER = {name: danger for name, _low, _high, danger in CHAMBERS}
ATOMS = tuple((sheet, name) for sheet in range(P) for name in CHAMBER_NAMES)
BUCKETS = tuple(
    (left, right, drift)
    for left in CHAMBER_NAMES
    for right in CHAMBER_NAMES
    for drift in range(P)
)
BUCKET_INDEX = {bucket: index for index, bucket in enumerate(BUCKETS)}
CONTROL_PAIRS = frozenset(((0, 0), (1, 0), (0, 1), (12, 12)))


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_integers(values: tuple[int, ...]) -> str:
    return hashlib.sha256(
        ",".join(str(value) for value in values).encode("ascii")
    ).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def load_bridge_module():
    require(lf_sha256(BRIDGE_PATH) == BRIDGE_SHA256, "bridge source hash drift")
    spec = importlib.util.spec_from_file_location("thm3479_guard_buckets", BRIDGE_PATH)
    require(spec is not None and spec.loader is not None, "bridge module loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_bridge_module()
CTX = None


def safe(chamber: str, sheet: int) -> int:
    return int(sheet % P not in CHAMBER_DANGER[chamber])


def context():
    global CTX
    if CTX is None:
        word = M.to_current(M.U_FULL_REL)
        t_den = 182 * M.lcm_tuple(word)
        nn = M.R_DILATION * t_den
        prime, root, prime_factors, lucas_witnesses = M.FULL_EMBEDDINGS[0]
        M.verify_lucas_certificate(
            prime, prime_factors, lucas_witnesses, "U_full guard-bucket prime"
        )
        M.verify_embedding(
            prime, root, nn, M.FULL_NN_FACTORS, "U_full guard-bucket embedding"
        )
        require(word[M.GUARD] == 1, ("guard speed", word[M.GUARD]))
        require(word[M.OWNER] == 13, ("owner speed", word[M.OWNER]))
        require(t_den % 91 == 0, ("guard unit", t_den))
        guard_unit = t_den // 91
        atom_intervals = tuple(
            (
                sheet,
                chamber,
                (7 * sheet + low) * guard_unit,
                (7 * sheet + high) * guard_unit,
            )
            for sheet in range(P)
            for chamber, low, high, _danger in CHAMBERS
        )
        require(atom_intervals[0][2] == 0, "atom partition left endpoint")
        require(atom_intervals[-1][3] == t_den, "atom partition right endpoint")
        require(
            all(left[3] == right[2] for left, right in zip(atom_intervals, atom_intervals[1:])),
            "atom partition gap",
        )
        atom_starts = tuple(row[2] for row in atom_intervals)
        zero = (0,) * 9
        q_intervals = M.fast_build_set(word, t_den, M.PATTERN_QA, zero)
        q_starts = [left for left, _right in q_intervals]
        embeddings = ((prime, root),)
        tabs = M.fast_make_tabs(q_intervals, M.X_FREQUENCY, nn, embeddings)
        zeta = pow(root, nn // P, prime)
        require(pow(zeta, P, prime) == 1 and zeta != 1, "bad order-13 root")
        tau_h = []
        tau_q5 = []
        for left_sheet, left_chamber in ATOMS:
            h_row = []
            q5_row = []
            for right_sheet, right_chamber in ATOMS:
                products = tuple(
                    safe(left_chamber, left_sheet + tau)
                    * safe(right_chamber, right_sheet + tau)
                    for tau in range(P)
                )
                q5_row.append(sum(products) % prime)
                h_row.append(
                    sum(
                        value * pow(zeta, -tau % P, prime)
                        for tau, value in enumerate(products)
                    )
                    % prime
                )
            tau_h.append(tuple(h_row))
            tau_q5.append(tuple(q5_row))
        CTX = (
            word,
            t_den,
            nn,
            prime,
            root,
            zeta,
            q_intervals,
            q_starts,
            embeddings,
            tabs,
            atom_intervals,
            atom_starts,
            tuple(tau_h),
            tuple(tau_q5),
        )
    return CTX


def split_by_guard_atom(
    intervals: list[tuple[int, int]],
    atom_intervals: tuple[tuple[int, str, int, int], ...],
    atom_starts: tuple[int, ...],
) -> tuple[tuple[tuple[int, int], ...], ...]:
    groups: list[list[tuple[int, int]]] = [[] for _atom in ATOMS]
    input_length = 0
    output_length = 0
    for interval_left, interval_right in intervals:
        require(interval_left < interval_right, ("bad E interval", interval_left, interval_right))
        input_length += interval_right - interval_left
        atom_index = bisect_right(atom_starts, interval_left) - 1
        require(atom_index >= 0, ("atom lookup", interval_left))
        while atom_index < len(atom_intervals):
            _sheet, _chamber, atom_left, atom_right = atom_intervals[atom_index]
            if atom_left >= interval_right:
                break
            left = max(interval_left, atom_left)
            right = min(interval_right, atom_right)
            if left < right:
                groups[atom_index].append((left, right))
                output_length += right - left
            if atom_right >= interval_right:
                break
            atom_index += 1
    require(output_length == input_length, ("atom partition changed measure", input_length, output_length))
    return tuple(tuple(group) for group in groups)


def atom_endpoint_values(
    groups: tuple[tuple[tuple[int, int], ...], ...],
) -> tuple[tuple[int, ...], tuple[int, ...], int]:
    (
        word,
        t_den,
        nn,
        prime,
        _root,
        _zeta,
        q_intervals,
        q_starts,
        embeddings,
        tabs,
        _atom_intervals,
        _atom_starts,
        _tau_h,
        _tau_q5,
    ) = context()
    ax_values = []
    by_values = []
    overlap = 0
    y_frequency = M.X_FREQUENCY + word[M.TARGET_B]
    for intervals in groups:
        if intervals:
            ax, atom_overlap = M.fast_x_sweep(
                list(intervals),
                q_intervals,
                q_starts,
                M.X_FREQUENCY,
                t_den,
                nn,
                embeddings,
                tabs,
            )
            by = M.fast_endpoint_sum(
                list(intervals), -y_frequency, nn, embeddings
            )
            ax_values.append(ax[0])
            by_values.append(by[0])
            overlap += atom_overlap
        else:
            ax_values.append(0)
            by_values.append(0)
    require(len(ax_values) == len(by_values) == len(ATOMS), "atom value count")
    return tuple(ax_values), tuple(by_values), overlap


def direct_guarded_values(alpha: int, beta: int, tau: int) -> tuple[int, int]:
    (
        word,
        t_den,
        nn,
        _prime,
        _root,
        _zeta,
        q_intervals,
        q_starts,
        embeddings,
        tabs,
        _atom_intervals,
        _atom_starts,
        _tau_h,
        _tau_q5,
    ) = context()
    ell = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    e_intervals = M.fast_build_set(word, t_den, M.PATTERN_E, ell)
    ax, _overlap = M.fast_x_sweep(
        e_intervals,
        q_intervals,
        q_starts,
        M.X_FREQUENCY,
        t_den,
        nn,
        embeddings,
        tabs,
    )
    by = M.fast_endpoint_sum(
        e_intervals,
        -(M.X_FREQUENCY + word[M.TARGET_B]),
        nn,
        embeddings,
    )
    return ax[0], by[0]


def worker(alpha: int) -> tuple[object, ...]:
    (
        word,
        t_den,
        _nn,
        prime,
        root,
        zeta,
        _q_intervals,
        _q_starts,
        _embeddings,
        _tabs,
        atom_intervals,
        atom_starts,
        tau_h,
        tau_q5,
    ) = context()
    no_guard_pattern = dict(M.PATTERN_E)
    require(no_guard_pattern.pop(M.GUARD) == "guard_safe", "guard factor absent")
    gamma_rows = []
    h_buckets = [0] * len(BUCKETS)
    q5_buckets = [0] * len(BUCKETS)
    raw_count = 0
    split_count = 0
    nonempty_atom_count = 0
    overlap_length = 0
    control_rows = []
    alpha_character = pow(zeta, -alpha % P, prime)
    for beta in range(P):
        ell = (0, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
        e_unguarded = M.fast_build_set(word, t_den, no_guard_pattern, ell)
        groups = split_by_guard_atom(e_unguarded, atom_intervals, atom_starts)
        middle_indices = tuple(
            index
            for index, (_sheet, chamber) in enumerate(ATOMS)
            if chamber == "middle"
        )
        require(
            all(not groups[index] for index in middle_indices),
            ("owner speed-13 gate leaked into middle chamber", alpha, beta),
        )
        ax_values, by_values, overlap = atom_endpoint_values(groups)
        raw_count += len(e_unguarded)
        split_count += sum(len(group) for group in groups)
        nonempty_atom_count += sum(bool(group) for group in groups)
        overlap_length += overlap
        phase = pow(root, beta * (context()[2] // P) % context()[2], prime)

        reconstructed = []
        for tau in range(P):
            ax = sum(
                value * safe(chamber, sheet + tau)
                for value, (sheet, chamber) in zip(ax_values, ATOMS)
            ) % prime
            by = sum(
                value * safe(chamber, sheet + tau)
                for value, (sheet, chamber) in zip(by_values, ATOMS)
            ) % prime
            reconstructed.append((ax, by))
            gamma_rows.append(phase * ax % prime * by % prime)
            if (alpha, beta) in CONTROL_PAIRS:
                direct = direct_guarded_values(alpha, beta, tau)
                require((ax, by) == direct, ("parent replay", alpha, beta, tau, (ax, by), direct))
                control_rows.append((beta, tau, ax, by))

        for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
            left_value = ax_values[left_index]
            if left_value == 0:
                continue
            for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
                right_value = by_values[right_index]
                if right_value == 0:
                    continue
                drift = (right_sheet - left_sheet) % P
                bucket = BUCKET_INDEX[(left_chamber, right_chamber, drift)]
                coefficient = (
                    phase * alpha_character % prime * left_value % prime * right_value
                ) % prime
                h_buckets[bucket] = (
                    h_buckets[bucket] + coefficient * tau_h[left_index][right_index]
                ) % prime
                q5_buckets[bucket] = (
                    q5_buckets[bucket] + coefficient * tau_q5[left_index][right_index]
                ) % prime

        # The atom product expansion is checked before any inverse transform.
        for tau, (ax, by) in enumerate(reconstructed):
            direct_product = phase * ax % prime * by % prime
            require(
                direct_product == gamma_rows[-P + tau],
                ("atom product", alpha, beta, tau),
            )

    return (
        alpha,
        raw_count,
        split_count,
        nonempty_atom_count,
        overlap_length,
        tuple(gamma_rows),
        tuple(h_buckets),
        tuple(q5_buckets),
        tuple(control_rows),
    )


def inverse_value(
    gamma: tuple[int, ...], q: tuple[int, int, int], prime: int, zeta: int
) -> int:
    total = 0
    index = 0
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                exponent = -(alpha * q[0] + beta * q[1] + tau * q[2]) % P
                total = (total + gamma[index] * pow(zeta, exponent, prime)) % prime
                index += 1
    return total * pow(P**3, -1, prime) % prime


def rank_mod(rows: tuple[tuple[int, ...], ...], prime: int) -> int:
    matrix = [list(value % prime for value in row) for row in rows]
    if not matrix:
        return 0
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [value * inverse % prime for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank:
                continue
            factor = matrix[index][column]
            if factor:
                matrix[index] = [
                    (value - factor * pivot_value) % prime
                    for value, pivot_value in zip(matrix[index], matrix[rank])
                ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def root_alignment_certificate(zeta: int, prime: int) -> tuple[object, ...]:
    """Classify the common C13-equivariant guard-sheet/root label gauges."""

    maps = tuple(
        tuple((address + gauge) % P for address in range(P))
        for gauge in range(P)
    )
    require(len(set(maps)) == P, "root-alignment gauge count")
    for gauge, mapping in enumerate(maps):
        require(len(set(mapping)) == P, ("root gauge not bijective", gauge))
        for address in range(P):
            for tau in range(P):
                require(
                    mapping[(address + tau) % P]
                    == (mapping[address] + tau) % P,
                    ("root gauge not equivariant", gauge, address, tau),
                )
        # Any equivariant map f obeys f(a)=f(0)+a by taking tau=a.
        require(
            all(mapping[address] == (mapping[0] + address) % P for address in range(P)),
            ("equivariant classification", gauge),
        )
        for frequency in range(1, P):
            gauge_phase = pow(zeta, frequency * gauge % P, prime)
            for address in range(P):
                require(
                    pow(zeta, frequency * mapping[address] % P, prime)
                    == gauge_phase
                    * pow(zeta, frequency * address % P, prime)
                    % prime,
                    ("root gauge phase", gauge, frequency, address),
                )
        for left in range(P):
            for right in range(P):
                require(
                    (mapping[right] - mapping[left]) % P == (right - left) % P,
                    ("common-gauge drift", gauge, left, right),
                )
    independent_gauge_hostile = (
        ((1 + 1) - (0 + 0)) % P,
        (1 - 0) % P,
    )
    require(
        independent_gauge_hostile[0] != independent_gauge_hostile[1],
        "independent endpoint gauges preserved drift",
    )
    return (
        len(maps),
        "u=a+c",
        "guard pullback g_C(a+tau) and root translation phi(y,u)+tau/13=phi(y,u+tau) both use +tau",
        "primitive orbit coefficient changes by zeta^(k*c)",
        "u_R-u_L=a_R-a_L under one common c",
        "independent c_L,c_R shift the drift by c_R-c_L",
        "label torsor only; chamber and physical-root realization remain absent",
        independent_gauge_hostile,
    )


def security_certificate(path: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            bad.append("Assert")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in {"eval", "exec", "compile", "__import__"}
        ):
            bad.append(node.func.id)
    require(not bad, ("security", bad))
    return len(tuple(ast.walk(tree))), tuple(bad)


def main() -> None:
    require(lf_sha256(GUARD_PATH) == GUARD_SHA256, "guard-sheet source hash drift")
    require(lf_sha256(REFINED_PATH) == REFINED_SHA256, "refined source hash drift")
    require(
        lf_sha256(REFINED_OUTPUT_PATH) == REFINED_OUTPUT_SHA256,
        "refined output hash drift",
    )
    refined_output = REFINED_OUTPUT_PATH.read_text(encoding="utf-8")
    require(str(EXPECTED_BRIDGE) in refined_output, "frozen bridge absent")
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(row[0] for row in chunks) == tuple(range(P)), "worker order")
    gamma = tuple(value for row in chunks for value in row[5])
    require(len(gamma) == P**3, ("gamma size", len(gamma)))
    gamma_hash = digest_integers(gamma)
    require(gamma_hash == EXPECTED_GAMMA_SHA256, ("gamma digest", gamma_hash))

    (
        word,
        t_den,
        nn,
        prime,
        root,
        zeta,
        q_intervals,
        _q_starts,
        _embeddings,
        _tabs,
        atom_intervals,
        _atom_starts,
        tau_h,
        tau_q5,
    ) = context()
    h_value = inverse_value(gamma, Q_H, prime, zeta)
    q5_value = inverse_value(gamma, Q_Q5, prime, zeta)
    bridge = (h_value - q5_value) % prime
    require((h_value, q5_value) == EXPECTED_VALUES, ("inverse values", h_value, q5_value))
    require(bridge == EXPECTED_BRIDGE, ("bridge", bridge))
    root_alignment = root_alignment_certificate(zeta, prime)

    normalizer = pow(P**3, -1, prime)
    h_buckets = tuple(
        sum(row[6][index] for row in chunks) % prime * normalizer % prime
        for index in range(len(BUCKETS))
    )
    q5_buckets = tuple(
        sum(row[7][index] for row in chunks) % prime * normalizer % prime
        for index in range(len(BUCKETS))
    )
    bridge_buckets = tuple(
        (left - right) % prime for left, right in zip(h_buckets, q5_buckets)
    )
    require(sum(h_buckets) % prime == h_value, "H bucket reconstruction")
    require(sum(q5_buckets) % prime == q5_value, "q5 bucket reconstruction")
    require(sum(bridge_buckets) % prime == bridge, "bridge bucket reconstruction")

    restriction_predicates = (
        ("same_sheet", lambda bucket: bucket[2] == 0),
        ("same_chamber", lambda bucket: bucket[0] == bucket[1]),
        (
            "same_guard_atom",
            lambda bucket: bucket[0] == bucket[1] and bucket[2] == 0,
        ),
    )
    restriction_rows = tuple(
        (
            name,
            sum(
                value
                for bucket, value in zip(BUCKETS, h_buckets)
                if predicate(bucket)
            )
            % prime,
            sum(
                value
                for bucket, value in zip(BUCKETS, q5_buckets)
                if predicate(bucket)
            )
            % prime,
            sum(
                value
                for bucket, value in zip(BUCKETS, bridge_buckets)
                if predicate(bucket)
            )
            % prime,
        )
        for name, predicate in restriction_predicates
    )
    require(
        all(row[3] != bridge for row in restriction_rows),
        ("guard equality unexpectedly recovered frozen bridge", restriction_rows),
    )

    # The owner gate leaves two chamber states per endpoint.  Therefore each
    # drift layer is canonically F_2^2: its four states are a K4 carrier, and
    # the Walsh basis is (constant,left,right,mixed).  This is a representation
    # statement; no orientation, hence no tournament, is manufactured.
    corner_states = (
        ("left", "left"),
        ("left", "right"),
        ("right", "left"),
        ("right", "right"),
    )
    corner_rows = tuple(
        tuple(
            bridge_buckets[BUCKET_INDEX[(left, right, drift)]]
            for drift in range(P)
        )
        for left, right in corner_states
    )
    walsh_rows = (
        tuple(sum(corner_rows[index][drift] for index in range(4)) % prime for drift in range(P)),
        tuple((corner_rows[0][drift] + corner_rows[1][drift] - corner_rows[2][drift] - corner_rows[3][drift]) % prime for drift in range(P)),
        tuple((corner_rows[0][drift] - corner_rows[1][drift] + corner_rows[2][drift] - corner_rows[3][drift]) % prime for drift in range(P)),
        tuple((corner_rows[0][drift] - corner_rows[1][drift] - corner_rows[2][drift] + corner_rows[3][drift]) % prime for drift in range(P)),
    )
    walsh_names = ("constant", "left", "right", "mixed")
    walsh_nonzero = tuple(sum(value != 0 for value in row) for row in walsh_rows)
    walsh_zero_drifts = tuple(
        (name, tuple(index for index, value in enumerate(row) if value == 0))
        for name, row in zip(walsh_names, walsh_rows)
    )
    walsh_spectra = tuple(
        tuple(
            sum(
                row[drift] * pow(zeta, -frequency * drift % P, prime)
                for drift in range(P)
            ) % prime
            for frequency in range(P)
        )
        for row in walsh_rows
    )
    walsh_spectral_nonzero = tuple(
        sum(value != 0 for value in row) for row in walsh_spectra
    )
    walsh_zero_frequencies = tuple(
        (name, tuple(index for index, value in enumerate(row) if value == 0))
        for name, row in zip(walsh_names, walsh_spectra)
    )
    corner_rank = rank_mod(corner_rows, prime)
    require(corner_rank == rank_mod(walsh_rows, prime), "Walsh rank changed")

    kernel_h_nonzero = sum(value != 0 for row in tau_h for value in row)
    kernel_q5_nonzero = sum(value != 0 for row in tau_q5 for value in row)
    require(kernel_h_nonzero == kernel_q5_nonzero == len(ATOMS) ** 2, "pair kernel zero")
    controls = tuple(item for row in chunks for item in row[8])
    require(len(controls) == len(CONTROL_PAIRS) * P, ("control count", len(controls)))
    control_hash = digest_json(controls)
    h_bucket_hash = digest_json(tuple(zip(BUCKETS, h_buckets)))
    q5_bucket_hash = digest_json(tuple(zip(BUCKETS, q5_buckets)))
    bridge_bucket_hash = digest_json(tuple(zip(BUCKETS, bridge_buckets)))
    zero_h = tuple(bucket for bucket, value in zip(BUCKETS, h_buckets) if value == 0)
    zero_q5 = tuple(bucket for bucket, value in zip(BUCKETS, q5_buckets) if value == 0)
    zero_bridge = tuple(bucket for bucket, value in zip(BUCKETS, bridge_buckets) if value == 0)
    expected_zero_buckets = tuple(
        bucket
        for bucket in BUCKETS
        if "middle" in bucket[:2]
    )
    require(
        zero_h == zero_q5 == zero_bridge == expected_zero_buckets,
        ("owner-boundary support pattern", zero_h, zero_q5, zero_bridge),
    )
    require(len(expected_zero_buckets) == 5 * P, "zero bucket count")
    security = security_certificate(ROOT / SCRIPT)
    semantic = (
        BRIDGE_SHA256,
        GUARD_SHA256,
        REFINED_SHA256,
        REFINED_OUTPUT_SHA256,
        (prime, root, t_den, nn),
        tuple(word),
        tuple(atom_intervals),
        tuple(BUCKETS),
        tuple(row[1] for row in chunks),
        tuple(row[2] for row in chunks),
        tuple(row[3] for row in chunks),
        tuple(row[4] for row in chunks),
        gamma_hash,
        (h_value, q5_value, bridge),
        h_buckets,
        q5_buckets,
        bridge_buckets,
        restriction_rows,
        root_alignment,
        expected_zero_buckets,
        corner_rows,
        walsh_rows,
        walsh_spectra,
        control_hash,
        security,
        "Cartesian endpoint atoms only; no common ancestry or current",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic hash", semantic_hash))

    print("LRC U_FULL GUARD-SHEET DRIFT-BUCKET BRIDGE DECOMPOSITION")
    print("status=FINITE-EXACT endpoint factorization; actual guard address, Cartesian endpoint pairing; LRC(14) OPEN")
    print(f"dependency_hashes={(BRIDGE_SHA256, GUARD_SHA256, REFINED_SHA256, REFINED_OUTPUT_SHA256)}")
    print(f"embedding=(prime={prime},root={root},order={nn})")
    print(f"common_guard_partition=(atoms={len(ATOMS)},sheets={P},chambers={CHAMBER_NAMES},half_open=((0,1),(1,6),(6,7)))")
    print("owner_boundary_gate=(speed=13,condition=distance(13*t,Z)<1/14,residual_support_half_open=([0,1/2),[13/2,7)),active_guard_atoms=26,middle_atoms_empty=True)")
    print(f"pair_kernel_types=(full_atom_pairs={len(ATOMS)**2},drift_types={len(BUCKETS)},H_nonzero={kernel_h_nonzero},q5_nonzero={kernel_q5_nonzero})")
    print(f"unguarded_interval_counts_by_alpha={tuple(row[1] for row in chunks)} total={sum(row[1] for row in chunks)}")
    print(f"atom_split_interval_counts_by_alpha={tuple(row[2] for row in chunks)} total={sum(row[2] for row in chunks)}")
    print(f"nonempty_atom_tables_by_alpha={tuple(row[3] for row in chunks)} total={sum(row[3] for row in chunks)}")
    print(f"overlap_lengths_by_alpha={tuple(row[4] for row in chunks)}")
    print(f"direct_parent_controls=(pairs={tuple(sorted(CONTROL_PAIRS))},tau_each={P},count={len(controls)},sha256={control_hash})")
    print(f"reconstructed_gamma_sha256={gamma_hash}")
    print(f"inverse_values=((q_H,{h_value}),(q_q5,{q5_value}),bridge={bridge})")
    print(f"bucket_nonzero_counts=(H={len(BUCKETS)-len(zero_h)},q5={len(BUCKETS)-len(zero_q5)},bridge={len(BUCKETS)-len(zero_bridge)})/{len(BUCKETS)}")
    print("bucket_support=exactly (left/right)x(left/right)xF13; every one of 52 types is nonzero for H, q5, and their bridge")
    print(f"bucket_zero_types_sha256={digest_json(expected_zero_buckets)}; schema=all 13 drifts for the five chamber pairs containing middle")
    print(f"bucket_sha256=(H={h_bucket_hash},q5={q5_bucket_hash},bridge={bridge_bucket_hash})")
    print(f"guard_equality_hostiles=(name,H,q5,bridge)={restriction_rows}; none recovers the frozen bridge")
    print(f"root_alignment_torsor={root_alignment}")
    print(f"K4xF13_corner_table=(rank={corner_rank},vertex_rows_sha256={digest_json(corner_rows)},walsh_rows_sha256={digest_json(walsh_rows)})")
    print(f"K4_walsh_support=(names={walsh_names},drift_nonzero={walsh_nonzero},zero_drifts={walsh_zero_drifts})")
    print(f"K4_walsh_drift_spectrum=(nonzero={walsh_spectral_nonzero},zero_frequencies={walsh_zero_frequencies},sha256={digest_json(walsh_spectra)})")
    print("K4_scope=four left/right chamber-pair states are the F2^2/K4 carrier at each drift; no pairwise orientation or tournament is inherited")
    print("factorization=remove only H guard; owner speed-13 in-gate restricts to 26 boundary atoms; restore g_C(a+tau); expand AX*BY; sum into 52 supported (C,D,d) buckets with common-sheet Fourier phase")
    print("positive_control=39-atom reconstruction matches the frozen full 13^3 gamma bank and both inverse role values exactly")
    print("hostile_boundary=117 drift buckets classify guard covariance only; their Cartesian endpoint pairs have no THM-2471 common-base/root/owner/word/source/horizon support relation")
    print("next_test=retain these guard buckets while adding one independently lawful ancestry-support predicate before multiplying endpoint factors")
    print("nonconsequence=no physical current, grouped exact-address coefficient, all-unit projector, scalar-row exclusion, or LRC(14)")
    print(f"security_ast={security}")
    print(f"semantic_sha256={semantic_hash}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
