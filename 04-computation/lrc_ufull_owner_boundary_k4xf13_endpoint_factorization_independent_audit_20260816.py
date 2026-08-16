#!/usr/bin/env python3
"""Independent audit of the U_full owner-boundary K4 x F_13 factorization.

The trusted input is the promoted THM-3479 endpoint engine.  This audit does
not import the candidate guard/bucket companion.  It independently

* derives the speed-13 owner support, including the strict-set versus
  half-open interval-representative boundary convention;
* removes only the H guard, partitions the resulting E-set by 39 common
  guard atoms, and restores the guard to reconstruct the frozen 13^3 bank;
* contracts the left/right atom tables centrally into 117 chamber/drift
  buckets, rather than accumulating buckets inside the endpoint workers;
* audits the supported 4 x 13 table, its four Walsh channels, and all drift
  Fourier modes; and
* keeps Cartesian endpoint pairing distinct from THM-2471 ancestry.

All endpoint values live in one certified split-field image.  Nonvanishing
there proves nonvanishing of the corresponding cyclotomic integer.  This is
not a physical current, an ancestry realization, a row exclusion, or an
LRC(14) result.
"""

from __future__ import annotations

import ast
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
OUTPUT = (
    "05-knowledge/results/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.out"
)
PARENT_PATH = ROOT / "04-computation/lrc_half_twist_relation_current_bridge_thm3479.py"
PARENT_SHA256 = "ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b"
CANDIDATE_COMMIT = "e92508714abec0f054e033ef933e04b61ae3f1e1"

P = 13
EXPECTED_GAMMA_SHA256 = "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
EXPECTED_VALUES = (320618948602619577408, 503604956476841920373)
EXPECTED_BRIDGE = 389266878372286537904
EXPECTED_SEMANTIC_SHA256 = "d52c9f0a56c14a83e1e6b175c7b725314c99f09d44509bc8582847a5857f7da6"

# The chamber endpoints are integer coordinates on s=91t.  They are encoded
# half-open because the endpoint engine integrates indicator functions and
# records interval endpoint differences.  Literal strict membership differs
# from this representative only at the two owner-boundary singleton points.
CHAMBERS = (
    ("left", 0, 1),
    ("middle", 1, 6),
    ("right", 6, 7),
)
CHAMBER_NAMES = tuple(name for name, _left, _right in CHAMBERS)
DANGER = {
    "left": frozenset((12, 0, 1)),
    "middle": frozenset((11, 12, 0, 1)),
    "right": frozenset((11, 12, 0)),
}
ATOMS = tuple((sheet, chamber) for sheet in range(P) for chamber in CHAMBER_NAMES)
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}
BUCKETS = tuple(
    (left, right, drift)
    for left in CHAMBER_NAMES
    for right in CHAMBER_NAMES
    for drift in range(P)
)
BUCKET_INDEX = {bucket: index for index, bucket in enumerate(BUCKETS)}
ACTIVE_CORNERS = (
    ("left", "left"),
    ("left", "right"),
    ("right", "left"),
    ("right", "right"),
)
CONTROL_PAIRS = frozenset(((0, 0), (0, 1), (1, 0), (6, 6), (12, 0), (12, 12)))


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
    return hashlib.sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_parent():
    require(lf_sha256(PARENT_PATH) == PARENT_SHA256, "THM-3479 source hash drift")
    spec = importlib.util.spec_from_file_location("thm3479_owner_boundary_audit", PARENT_PATH)
    require(spec is not None and spec.loader is not None, "parent module loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_parent()
CTX = None


def safe(chamber: str, sheet: int) -> int:
    return int(sheet % P not in DANGER[chamber])


def direct_guard_safe(sheet: int, residual_numerator: int, tau: int) -> int:
    """Guard value at r=residual_numerator/2, away from chamber boundaries."""
    doubled_s = (14 * (sheet + tau) + residual_numerator) % 182
    distance = min(doubled_s, 182 - doubled_s)
    return int(distance > 26)


def boundary_certificate() -> tuple[object, ...]:
    # Midpoint route derives each danger arc without consulting the declared
    # table.  Residual numerators 1, 7, 13 represent r=1/2,7/2,13/2.
    midpoint_numerators = {"left": 1, "middle": 7, "right": 13}
    derived = {}
    for chamber in CHAMBER_NAMES:
        danger = frozenset(
            sheet
            for sheet in range(P)
            if not direct_guard_safe(sheet, midpoint_numerators[chamber], 0)
        )
        require(danger == DANGER[chamber], ("guard danger arc", chamber, danger))
        derived[chamber] = tuple(sorted(danger))

    # At the chamber cuts, use the endpoint engine's half-open unsafe window
    # [-13-7*tau,13-7*tau) on the 91-circle.  Hence r=1 belongs to middle and
    # r=6 belongs to right.  This is an a.e. representative, not a claim that
    # all strict inequalities include those boundary points.
    def half_open_guard_danger(sheet: int, residual: int) -> bool:
        value = (7 * sheet + residual) % 91
        return value < 13 or value >= 78

    for residual, chamber in ((0, "left"), (1, "middle"), (6, "right")):
        danger = frozenset(
            sheet for sheet in range(P) if half_open_guard_danger(sheet, residual)
        )
        require(danger == DANGER[chamber], ("guard chamber cut", residual, danger))

    # The owner has speed 13 and shift zero.  Since s=91t=7a+r,
    # 13t=a+r/7 mod 1.  Literal strict support is
    # [0,1/2) union (13/2,7).  The fast interval engine chooses the
    # measure-equivalent half-open representative [0,1/2) union [13/2,7).
    strict_probes = (
        (0, True),
        (2, False),   # r=1/2, encoded in quarter-r units
        (24, False),
        (26, False),  # r=13/2 is excluded literally
        (27, True),   # r=27/4 lies in the right owner sliver
    )
    for four_r, expected in strict_probes:
        distance = min(four_r, 28 - four_r)
        literal = distance < 2
        require(literal == expected, ("owner strict probe", four_r, literal))

    # Compare PATTERN_E's owner primitive with the independently written
    # expected half-open windows.  In units T_DEN/182, each centre 14*a has
    # [14*a-1,14*a+1), cyclically.  This is the exact source-level convention
    # from fast_in_comb, not an inference from the candidate output.
    word = M.to_current(M.U_FULL_REL)
    t_den = 182 * M.lcm_tuple(word)
    unit = t_den // 182
    owner_intervals = M.fast_build_set(
        word, t_den, {M.OWNER: "in"}, (0,) * 9
    )
    expected_owner_intervals = []
    for sheet in range(P):
        left = (14 * sheet - 1) * unit
        right = (14 * sheet + 1) * unit
        if left < 0:
            expected_owner_intervals.extend(((0, right), (t_den + left, t_den)))
        else:
            expected_owner_intervals.append((left, right))
    expected_owner_intervals.sort()
    require(
        owner_intervals == expected_owner_intervals,
        ("owner half-open intervals", owner_intervals, expected_owner_intervals),
    )
    half_open_singleton_difference = (
        "r=13/2 in each sheet",
        P,
        "engine_in/literal_out",
    )
    return (
        tuple((name, derived[name]) for name in CHAMBER_NAMES),
        "literal_owner=[0,1/2)U(13/2,7)",
        "interval_representative=[0,1/2)U[13/2,7)",
        half_open_singleton_difference,
    )


def context():
    global CTX
    if CTX is None:
        word = M.to_current(M.U_FULL_REL)
        t_den = 182 * M.lcm_tuple(word)
        nn = M.R_DILATION * t_den
        prime, root, prime_factors, lucas_witnesses = M.FULL_EMBEDDINGS[0]
        M.verify_lucas_certificate(
            prime, prime_factors, lucas_witnesses, "owner-boundary audit prime"
        )
        M.verify_embedding(
            prime, root, nn, M.FULL_NN_FACTORS, "owner-boundary audit embedding"
        )
        require(word == (1, 183, 27, 131, 53, 313, 13, 2197, 742586), word)
        require(word[M.OWNER] == 13 and word[M.GUARD] == 1, "owner/guard speeds")
        require(t_den % 91 == 0, ("guard grid", t_den))
        zero = (0,) * 9
        q_intervals = M.fast_build_set(word, t_den, M.PATTERN_QA, zero)
        q_starts = [left for left, _right in q_intervals]
        embeddings = ((prime, root),)
        tabs = M.fast_make_tabs(q_intervals, M.X_FREQUENCY, nn, embeddings)
        zeta = pow(root, nn // P, prime)
        require(pow(zeta, P, prime) == 1 and zeta != 1, "order-13 root")
        guard_unit = t_den // 91
        atom_intervals = tuple(
            (
                sheet,
                chamber,
                (7 * sheet + low) * guard_unit,
                (7 * sheet + high) * guard_unit,
            )
            for sheet in range(P)
            for chamber, low, high in CHAMBERS
        )
        require(atom_intervals[0][2] == 0, "atom partition start")
        require(atom_intervals[-1][3] == t_den, "atom partition stop")
        require(
            all(left[3] == right[2] for left, right in zip(atom_intervals, atom_intervals[1:])),
            "atom partition gap",
        )
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
        )
    return CTX


def partition_two_pointer(
    intervals: list[tuple[int, int]],
    atom_intervals: tuple[tuple[int, str, int, int], ...],
) -> tuple[tuple[tuple[int, int], ...], ...]:
    """Intersect sorted E intervals with the contiguous atom partition."""
    groups: list[list[tuple[int, int]]] = [[] for _atom in ATOMS]
    atom_index = 0
    input_measure = 0
    output_measure = 0
    for interval_left, interval_right in intervals:
        require(interval_left < interval_right, ("bad E interval", interval_left, interval_right))
        input_measure += interval_right - interval_left
        while atom_intervals[atom_index][3] <= interval_left:
            atom_index += 1
        scan = atom_index
        while scan < len(atom_intervals) and atom_intervals[scan][2] < interval_right:
            _sheet, _chamber, atom_left, atom_right = atom_intervals[scan]
            left = max(interval_left, atom_left)
            right = min(interval_right, atom_right)
            if left < right:
                groups[scan].append((left, right))
                output_measure += right - left
            if atom_right >= interval_right:
                break
            scan += 1
    require(input_measure == output_measure, ("atom measure", input_measure, output_measure))
    return tuple(tuple(group) for group in groups)


def atom_endpoint_tables(
    groups: tuple[tuple[tuple[int, int], ...], ...],
) -> tuple[tuple[int, ...], tuple[int, ...], int]:
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
    ) = context()
    ax_values = []
    by_values = []
    overlap = 0
    y_frequency = M.X_FREQUENCY + word[M.TARGET_B]
    for group in groups:
        if not group:
            ax_values.append(0)
            by_values.append(0)
            continue
        intervals = list(group)
        ax, atom_overlap = M.fast_x_sweep(
            intervals,
            q_intervals,
            q_starts,
            M.X_FREQUENCY,
            t_den,
            nn,
            embeddings,
            tabs,
        )
        by = M.fast_endpoint_sum(intervals, -y_frequency, nn, embeddings)
        ax_values.append(ax[0])
        by_values.append(by[0])
        overlap += atom_overlap
    require(len(ax_values) == len(by_values) == len(ATOMS), "atom table size")
    return tuple(ax_values), tuple(by_values), overlap


def direct_guarded(alpha: int, beta: int, tau: int) -> tuple[int, int]:
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
    ) = context()
    ell = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    intervals = M.fast_build_set(word, t_den, M.PATTERN_E, ell)
    ax, _overlap = M.fast_x_sweep(
        intervals,
        q_intervals,
        q_starts,
        M.X_FREQUENCY,
        t_den,
        nn,
        embeddings,
        tabs,
    )
    by = M.fast_endpoint_sum(
        intervals,
        -(M.X_FREQUENCY + word[M.TARGET_B]),
        nn,
        embeddings,
    )
    return ax[0], by[0]


def worker(alpha: int) -> tuple[object, ...]:
    (
        word,
        _t_den,
        _nn,
        _prime,
        _root,
        _zeta,
        _q_intervals,
        _q_starts,
        _embeddings,
        _tabs,
        atom_intervals,
    ) = context()
    no_guard = dict(M.PATTERN_E)
    require(no_guard.pop(M.GUARD) == "guard_safe", "removed factor is not H guard")
    rows = []
    direct_rows = []
    raw_count = 0
    split_count = 0
    nonempty_count = 0
    overlap_count = 0
    for beta in range(P):
        ell = (0, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
        unguarded = M.fast_build_set(word, context()[1], no_guard, ell)
        groups = partition_two_pointer(unguarded, atom_intervals)
        middle = tuple(
            index for index, (_sheet, chamber) in enumerate(ATOMS) if chamber == "middle"
        )
        require(all(not groups[index] for index in middle), ("middle owner leak", alpha, beta))
        ax_values, by_values, overlap = atom_endpoint_tables(groups)
        rows.append((beta, ax_values, by_values))
        raw_count += len(unguarded)
        split_count += sum(len(group) for group in groups)
        nonempty_count += sum(bool(group) for group in groups)
        overlap_count += overlap
        if (alpha, beta) in CONTROL_PAIRS:
            for tau in range(P):
                direct_rows.append((beta, tau, direct_guarded(alpha, beta, tau)))
    return (
        alpha,
        tuple(rows),
        raw_count,
        split_count,
        nonempty_count,
        overlap_count,
        tuple(direct_rows),
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
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [value * inverse % prime for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank:
                continue
            factor = matrix[row][column]
            if factor:
                matrix[row] = [
                    (value - factor * pivot_value) % prime
                    for value, pivot_value in zip(matrix[row], matrix[rank])
                ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def security_certificate(path: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    forbidden = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            forbidden.append("Assert")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in {"eval", "exec", "compile", "__import__"}
        ):
            forbidden.append(node.func.id)
    require(not forbidden, ("security", forbidden))
    return len(tuple(ast.walk(tree))), tuple(forbidden)


def main() -> None:
    boundary = boundary_certificate()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")

    (
        _word,
        t_den,
        nn,
        prime,
        root,
        zeta,
        _q_intervals,
        _q_starts,
        _embeddings,
        _tabs,
        atom_intervals,
    ) = context()
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(len(tables) == P**2, ("table count", len(tables)))

    # Fresh reconstruction route: atom tables are materialized first, then the
    # guard is restored centrally.  The candidate instead accumulated its
    # drift buckets inside each endpoint worker.
    gamma_rows = []
    reconstructed_controls = []
    table_index = 0
    for alpha in range(P):
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, ("beta order", alpha, beta, stored_beta))
            phase = pow(zeta, beta, prime)
            for tau in range(P):
                ax = sum(
                    value * safe(chamber, sheet + tau)
                    for value, (sheet, chamber) in zip(ax_values, ATOMS)
                ) % prime
                by = sum(
                    value * safe(chamber, sheet + tau)
                    for value, (sheet, chamber) in zip(by_values, ATOMS)
                ) % prime
                gamma_rows.append(phase * ax % prime * by % prime)
                if (alpha, beta) in CONTROL_PAIRS:
                    reconstructed_controls.append((alpha, beta, tau, ax, by))
    gamma = tuple(gamma_rows)
    require(len(gamma) == P**3, ("gamma size", len(gamma)))
    gamma_hash = digest_integers(gamma)
    require(gamma_hash == EXPECTED_GAMMA_SHA256, ("frozen gamma digest", gamma_hash))

    direct_controls = tuple(
        (chunk[0], beta, tau, values[0], values[1])
        for chunk in chunks
        for beta, tau, values in chunk[6]
    )
    require(
        tuple(reconstructed_controls) == direct_controls,
        ("direct guarded reconstruction", reconstructed_controls[:2], direct_controls[:2]),
    )

    h_value = inverse_value(gamma, (1, 0, 1), prime, zeta)
    q5_value = inverse_value(gamma, (1, 0, 0), prime, zeta)
    bridge = (h_value - q5_value) % prime
    require((h_value, q5_value) == EXPECTED_VALUES, ("role values", h_value, q5_value))
    require(bridge == EXPECTED_BRIDGE, ("bridge", bridge))

    # Direct pair kernels; no guard-sidecar source is imported.
    h_kernels = []
    q5_kernels = []
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
        h_kernels.append(tuple(h_row))
        q5_kernels.append(tuple(q5_row))
    h_kernels_t = tuple(h_kernels)
    q5_kernels_t = tuple(q5_kernels)
    require(
        all(value != 0 for row in h_kernels_t for value in row)
        and all(value != 0 for row in q5_kernels_t for value in row),
        "abstract pair-kernel zero",
    )

    h_buckets = [0] * len(BUCKETS)
    q5_buckets = [0] * len(BUCKETS)
    table_index = 0
    for alpha in range(P):
        alpha_weight = pow(zeta, -alpha % P, prime)
        for beta in range(P):
            _stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            phase_weight = pow(zeta, beta, prime) * alpha_weight % prime
            for left_index, (left_sheet, left_chamber) in enumerate(ATOMS):
                left_value = ax_values[left_index]
                if left_value == 0:
                    continue
                for right_index, (right_sheet, right_chamber) in enumerate(ATOMS):
                    right_value = by_values[right_index]
                    if right_value == 0:
                        continue
                    drift = (right_sheet - left_sheet) % P
                    bucket_index = BUCKET_INDEX[(left_chamber, right_chamber, drift)]
                    coefficient = phase_weight * left_value % prime * right_value % prime
                    h_buckets[bucket_index] = (
                        h_buckets[bucket_index]
                        + coefficient * h_kernels_t[left_index][right_index]
                    ) % prime
                    q5_buckets[bucket_index] = (
                        q5_buckets[bucket_index]
                        + coefficient * q5_kernels_t[left_index][right_index]
                    ) % prime
    normalizer = pow(P**3, -1, prime)
    h_bucket_tuple = tuple(value * normalizer % prime for value in h_buckets)
    q5_bucket_tuple = tuple(value * normalizer % prime for value in q5_buckets)
    bridge_buckets = tuple(
        (left - right) % prime for left, right in zip(h_bucket_tuple, q5_bucket_tuple)
    )
    require(sum(h_bucket_tuple) % prime == h_value, "H bucket reconstruction")
    require(sum(q5_bucket_tuple) % prime == q5_value, "q5 bucket reconstruction")
    require(sum(bridge_buckets) % prime == bridge, "bridge bucket reconstruction")

    expected_zero = tuple(bucket for bucket in BUCKETS if "middle" in bucket[:2])
    zero_h = tuple(bucket for bucket, value in zip(BUCKETS, h_bucket_tuple) if value == 0)
    zero_q5 = tuple(bucket for bucket, value in zip(BUCKETS, q5_bucket_tuple) if value == 0)
    zero_bridge = tuple(bucket for bucket, value in zip(BUCKETS, bridge_buckets) if value == 0)
    require(zero_h == zero_q5 == zero_bridge == expected_zero, "52-bucket support law")
    require(len(BUCKETS) - len(expected_zero) == 52, ("supported count", len(expected_zero)))

    corner_rows = tuple(
        tuple(
            bridge_buckets[BUCKET_INDEX[(left, right, drift)]]
            for drift in range(P)
        )
        for left, right in ACTIVE_CORNERS
    )
    walsh_rows = (
        tuple(sum(corner_rows[row][drift] for row in range(4)) % prime for drift in range(P)),
        tuple((corner_rows[0][drift] + corner_rows[1][drift] - corner_rows[2][drift] - corner_rows[3][drift]) % prime for drift in range(P)),
        tuple((corner_rows[0][drift] - corner_rows[1][drift] + corner_rows[2][drift] - corner_rows[3][drift]) % prime for drift in range(P)),
        tuple((corner_rows[0][drift] - corner_rows[1][drift] - corner_rows[2][drift] + corner_rows[3][drift]) % prime for drift in range(P)),
    )
    walsh_names = ("constant", "left", "right", "mixed")
    corner_rank = rank_mod(corner_rows, prime)
    walsh_rank = rank_mod(walsh_rows, prime)
    require(corner_rank == walsh_rank == 4, ("K4 rank", corner_rank, walsh_rank))
    require(all(value != 0 for row in corner_rows for value in row), "corner drift zero")
    require(all(value != 0 for row in walsh_rows for value in row), "Walsh drift zero")
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
    require(
        all(value != 0 for row in walsh_spectra for value in row),
        "Walsh drift-Fourier zero",
    )

    # Equality hostiles use only guard geometry and therefore cannot be called
    # ancestry.  None reconstructs the full bridge.
    restrictions = (
        ("same_sheet", lambda bucket: bucket[2] == 0),
        ("same_chamber", lambda bucket: bucket[0] == bucket[1]),
        ("same_guard_atom", lambda bucket: bucket[2] == 0 and bucket[0] == bucket[1]),
    )
    restriction_rows = tuple(
        (
            name,
            sum(value for bucket, value in zip(BUCKETS, bridge_buckets) if predicate(bucket)) % prime,
        )
        for name, predicate in restrictions
    )
    require(all(value != bridge for _name, value in restriction_rows), restriction_rows)

    control_hash = digest_json(direct_controls)
    security = security_certificate(ROOT / SCRIPT)
    semantic = (
        CANDIDATE_COMMIT,
        PARENT_SHA256,
        boundary,
        (prime, root, t_den, nn),
        tuple(atom_intervals),
        tuple(chunk[2] for chunk in chunks),
        tuple(chunk[3] for chunk in chunks),
        tuple(chunk[4] for chunk in chunks),
        tuple(chunk[5] for chunk in chunks),
        gamma_hash,
        (h_value, q5_value, bridge),
        h_bucket_tuple,
        q5_bucket_tuple,
        bridge_buckets,
        corner_rows,
        walsh_rows,
        walsh_spectra,
        restriction_rows,
        control_hash,
        security,
        "K4 carrier without orientation",
        "Cartesian endpoint pairs without THM-2471 ancestry support",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic hash", semantic_hash))

    print("LRC U_FULL OWNER-BOUNDARY K4xF13 ENDPOINT FACTORIZATION INDEPENDENT AUDIT")
    print("status=FINITE-EXACT INDEPENDENT AUDIT; endpoint factorization only; LRC(14) OPEN")
    print(f"candidate_commit={CANDIDATE_COMMIT}; imported_candidate_code=False")
    print(f"parent_dependency={PARENT_PATH.name}:{PARENT_SHA256}")
    print(f"embedding=(prime={prime},root={root},order={nn})")
    print(f"owner_support=(literal=[0,1/2)U(13/2,7),interval_representative=[0,1/2)U[13/2,7),singleton_difference=13 left endpoints r=13/2)")
    print("guard_partition=(39=13*3 atoms; chambers=[0,1),[1,6),[6,7); r=1 middle, r=6 right in half-open engine convention)")
    print(f"guard_danger_arcs={tuple((name,tuple(sorted(DANGER[name]))) for name in CHAMBER_NAMES)}")
    print("owner_gate=(26 active left/right atoms; all 13 middle atoms empty before Fourier contraction)")
    print(f"unguarded_interval_counts_by_alpha={tuple(chunk[2] for chunk in chunks)} total={sum(chunk[2] for chunk in chunks)}")
    print(f"atom_split_interval_counts_by_alpha={tuple(chunk[3] for chunk in chunks)} total={sum(chunk[3] for chunk in chunks)}")
    print(f"nonempty_atom_tables_by_alpha={tuple(chunk[4] for chunk in chunks)} total={sum(chunk[4] for chunk in chunks)}")
    print(f"direct_guarded_controls=(pairs={tuple(sorted(CONTROL_PAIRS))},tau_each={P},count={len(direct_controls)},sha256={control_hash})")
    print(f"frozen_bank_reconstruction=(size={len(gamma)},sha256={gamma_hash})")
    print(f"inverse_values=((q_H,{h_value}),(q_q5,{q5_value}),bridge={bridge})")
    print(f"abstract_guard_pair_kernels=(H_nonzero={sum(value != 0 for row in h_kernels_t for value in row)},q5_nonzero={sum(value != 0 for row in q5_kernels_t for value in row)})/{len(ATOMS)**2}")
    print("drift_bucket_support=(52/117 exactly; (left/right)^2 x F13; all H,q5,bridge supported buckets nonzero)")
    print(f"bucket_hashes=(H={digest_json(tuple(zip(BUCKETS,h_bucket_tuple)))},q5={digest_json(tuple(zip(BUCKETS,q5_bucket_tuple)))},bridge={digest_json(tuple(zip(BUCKETS,bridge_buckets)))})")
    print(f"K4xF13_table=(rank={corner_rank},corner_sha256={digest_json(corner_rows)},walsh_sha256={digest_json(walsh_rows)})")
    print(f"Walsh_channels=(names={walsh_names},drift_support={(13,13,13,13)},Fourier_support={(13,13,13,13)},spectrum_sha256={digest_json(walsh_spectra)})")
    print("K4_scope=four F2^2 chamber-pair states are vertices of an undirected complete carrier; no intrinsic pairwise orientation, hence no tournament")
    print(f"guard_equality_hostiles={restriction_rows}; none reconstructs the frozen bridge")
    print("pairing_scope=AX and BY atoms are multiplied by Cartesian product; no common (y,u,v,a,b,e') stalk, root, owner, word, source, horizon, or support predicate from THM-2471")
    print("nonconsequence=no ancestry/bispectrum realization, physical current, grouped C(a;X,m), all-unit B(q), scalar-row exclusion, or LRC(14)")
    print(f"security_ast={security}")
    print(f"semantic_sha256={semantic_hash}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
