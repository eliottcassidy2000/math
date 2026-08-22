#!/usr/bin/env python3
"""Test the canonical desheeted common-base coupling for U_full endpoints.

Write the endpoint guard coordinate as

    91 t = 7 a + r,       a in F_13, 0 <= r < 7,
    y = r/7 = 13 t - a,   t = (y+a)/13.

This is exactly the THM-2471 owner-node form w_a=(y+a)/13.  Two endpoint
points share the proposed outer base when their desheeted residual y agrees;
they need not be the same circle point.  For every U_full character, this
probe pushes all thirteen endpoint sheets to y, inserts Q(13y), multiplies
the two endpoint indicators before summing, and retains THM-2594's seven
cells on y.

The same-sheet sector must reproduce the previously audited common-point
diagonal exactly.  The full square of the sheet count adds the genuinely new
cross-sheet/common-residual sector.  This is an actual one-base endpoint
integrand, but it supplies no collision horizon, source packet, chronology,
grouped exact-address coefficient, row exclusion, or LRC(14) conclusion.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TARGET_PATH = (
    ROOT
    / "04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py"
)
TARGET_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
POINT_PATH = ROOT / "04-computation/lrc_endpoint_ufull_atom_bridge_kernel_probe_20260816.py"
POINT_SHA256 = "d1182de4d777bab20a8d423cf942151ac3149014b67d9c34883cbce37a7b0a9f"
EXPECTED_POINT_GAMMA_SHA256 = "771545a5cb1f0f03459b8d351de668ad950ece5fcb985fa61d599d643de3303f"
EXPECTED_POINT_VALUES = (
    633668780131603861,
    405160484437854840264,
    167726070588785644466,
)
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1199656, 186244)
EXPECTED_FULL_GAMMA_SHA256 = "c16e6c4e82b73aced92f5d9ab61779d535263deed17a9caa90fce0ca6c947dcf"
EXPECTED_CROSS_GAMMA_SHA256 = "f02e10902bd8cc5ec9b588f059c7ebb72eac14b5f9ee46f6055f2f2ee0a530c6"
EXPECTED_TABLE_SHA256 = (
    "64adfc3850b8a8924546a75e051116aa10e3a62661ca409175aad174a9ebb926",
    "b80423dad8b155bd02e79d98c981ec3ada15934d91cb72ebacc385ded6b7ed6d",
    "50d44c5db5cea5842ceada123d06c54f0221561d08744c83eaf3aebd55f2420c",
)
EXPECTED_SPECTRUM_SHA256 = (
    "398a2c2ccb57714aa9d735d3bca4f24c2212396a8be9e2bd07e8e4ec57b6c05a",
    "17e08d5c0a5d107bad1d9828e754ed197d2737ee4dc7c1fcbfdeccb410992f72",
    "14c5490ee52e4283848290feee3b0995a44791fee30a9aaef6ad2a945b9442e8",
)
EXPECTED_FULL_VALUES = (
    95783417057771114126,
    124341028951542154618,
    543695274352737840377,
)
EXPECTED_CROSS_VALUES = (
    95149748277639510265,
    291433430760196195223,
    375969203763952195911,
)
EXPECTED_RELATION_CELL_PROFILE = (
    289814661037836286866, 0, 0, 0, 0, 0, 0,
)
EXPECTED_SEMANTIC_SHA256 = "97694f50b071dbf875802223dd894d7e1df86ef65a194c22d588505c59eea507"

P = 13
Q = 7
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)


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


def load_target():
    require(lf_sha256(TARGET_PATH) == TARGET_SHA256, "target source drift")
    require(lf_sha256(POINT_PATH) == POINT_SHA256, "point-hostile source drift")
    spec = importlib.util.spec_from_file_location("ufull_desheet_target", TARGET_PATH)
    require(spec is not None and spec.loader is not None, "target module loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_target()
M = T.M
WORKER_CONTEXT = None
SAFE_COUNTS = None


def worker_context():
    global WORKER_CONTEXT
    if WORKER_CONTEXT is None:
        (
            word,
            t_den,
            nn,
            prime,
            root,
            zeta,
            q_intervals,
            q_starts,
            embeddings,
            _tabs,
            atom_intervals,
        ) = T.context()
        require(word[M.TARGET_B] % (P * P) == 0,
                ("target frequency is not 13^2-divisible", word[M.TARGET_B]))
        diagonal_frequency = -word[M.TARGET_B]
        diagonal_tabs = M.fast_make_tabs(
            q_intervals, diagonal_frequency, nn, embeddings
        )[0]
        require(t_den % (2 * Q * P) == 0, "base grid misses 7/13 boundaries")
        WORKER_CONTEXT = (
            word,
            t_den,
            nn,
            prime,
            root,
            zeta,
            q_intervals,
            q_starts,
            diagonal_frequency,
            diagonal_tabs,
            atom_intervals,
        )
    return WORKER_CONTEXT


def safe_count_table():
    global SAFE_COUNTS
    if SAFE_COUNTS is None:
        tables = {}
        for chamber in ("left", "middle", "right"):
            tables[chamber] = tuple(
                tuple(
                    sum(
                        ((mask >> sheet) & 1) * T.safe(chamber, sheet + tau)
                        for sheet in range(P)
                    )
                    for tau in range(P)
                )
                for mask in range(1 << P)
            )
        SAFE_COUNTS = tables
    return SAFE_COUNTS


def cell_of_segment(left: int, right: int, t_den: int) -> int:
    midpoint_numerator = 14 * (left + right)
    band = midpoint_numerator // (2 * t_den)
    require(0 <= band < 14, ("cell band", left, right, band))
    if band in (0, 13):
        return 0
    return (band + 1) // 2


def chamber_of_segment(left: int, right: int, t_den: int) -> str:
    midpoint = left + right
    if Q * midpoint < 2 * t_den:
        return "left"
    if Q * midpoint < 12 * t_den:
        return "middle"
    return "right"


def q_endpoint_jump(left: int, right: int) -> int:
    """Endpoint numerator on Q(13y) for one y segment."""
    (
        word,
        t_den,
        nn,
        prime,
        root,
        _zeta,
        q_intervals,
        q_starts,
        diagonal_frequency,
        diagonal_tabs,
        _atom_intervals,
    ) = worker_context()
    scaled_left = P * left
    start = scaled_left % t_den
    span = P * (right - left)
    require(0 < span <= t_den, ("desheet span", left, right, span))
    stop = start + span
    require(stop <= t_den, ("desheet branch crossing", left, right, start, stop))

    index = bisect_right(q_starts, start) - 1
    if index < 0:
        index = 0
    elif q_intervals[index][1] <= start:
        index += 1
    total = 0
    while index < len(q_intervals):
        q_left, q_right = q_intervals[index]
        if q_left >= stop:
            break
        overlap_left = max(start, q_left)
        overlap_right = min(stop, q_right)
        if overlap_left < overlap_right:
            if overlap_left == q_left:
                value_left = diagonal_tabs[0][index]
            else:
                value_left = pow(
                    root, (-diagonal_frequency * overlap_left) % nn, prime
                )
            if overlap_right == q_right:
                value_right = diagonal_tabs[1][index]
            else:
                value_right = pow(
                    root, (-diagonal_frequency * overlap_right) % nn, prime
                )
            total = (total + value_left - value_right) % prime
        index += 1

    # If y=(x+j)/13, the discarded branch phase is
    # exp(2*pi*i*(word[TARGET_B]/13^2)*j)=1 because that quotient is integral.
    require(word[M.TARGET_B] % (P * P) == 0, "branch phase hostile")
    return total


def desheet_profile(alpha: int, beta: int):
    """Return full-pair and same-sheet cell rows for all guard shifts tau."""
    (
        word,
        t_den,
        _nn,
        prime,
        _root,
        _zeta,
        _q_intervals,
        _q_starts,
        _diagonal_frequency,
        _diagonal_tabs,
        atom_intervals,
    ) = worker_context()
    no_guard = dict(M.PATTERN_E)
    require(no_guard.pop(M.GUARD) == "guard_safe", "removed factor is not guard")
    shift = (0, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    unguarded = M.fast_build_set(word, t_den, no_guard, shift)
    groups = T.partition_two_pointer(unguarded, atom_intervals)
    require(all(not groups[index]
                for index, (_sheet, chamber) in enumerate(T.ATOMS)
                if chamber == "middle"),
            ("middle owner leak", alpha, beta))

    events: dict[int, int] = {}
    mapped_intervals = 0
    for atom_index, intervals in enumerate(groups):
        sheet, _chamber = T.ATOMS[atom_index]
        bit = 1 << sheet
        for interval_left, interval_right in intervals:
            y_left = P * interval_left - sheet * t_den
            y_right = P * interval_right - sheet * t_den
            require(0 <= y_left < y_right <= t_den,
                    ("desheet interval", atom_index, y_left, y_right))
            events[y_left] = events.get(y_left, 0) ^ bit
            events[y_right] = events.get(y_right, 0) ^ bit
            mapped_intervals += 1

    # Split at every x=13y branch, every 1/14 cell band, and both chamber walls.
    boundaries = {0, t_den}
    boundaries.update(index * (t_den // P) for index in range(P + 1))
    boundaries.update(index * (t_den // (2 * Q)) for index in range(2 * Q + 1))
    for boundary in boundaries:
        events.setdefault(boundary, 0)

    full = [[0 for _ell in range(Q)] for _tau in range(P)]
    same = [[0 for _ell in range(Q)] for _tau in range(P)]
    mask = 0
    positions = sorted(events)
    tables = safe_count_table()
    active_segments = 0
    q_active_segments = 0
    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        chamber = chamber_of_segment(left, right, t_den)
        require(chamber != "middle", ("middle active segment", left, right, mask))
        active_segments += 1
        jump = q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        ell = cell_of_segment(left, right, t_den)
        counts = tables[chamber][mask]
        for tau, count in enumerate(counts):
            full[tau][ell] = (full[tau][ell] + count * count * jump) % prime
            same[tau][ell] = (same[tau][ell] + count * jump) % prime
    mask ^= events[positions[-1]]
    require(mask == 0, ("desheet profile did not close", alpha, beta, mask))
    return (
        tuple(tuple(row) for row in full),
        tuple(tuple(row) for row in same),
        len(unguarded),
        mapped_intervals,
        active_segments,
        q_active_segments,
    )


def direct_guarded_profile(alpha: int, beta: int, tau: int):
    """Independent guard-restored control using the literal shifted E set."""
    (
        word,
        t_den,
        _nn,
        prime,
        _root,
        _zeta,
        _q_intervals,
        _q_starts,
        _diagonal_frequency,
        _diagonal_tabs,
        atom_intervals,
    ) = worker_context()
    shift = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    guarded = M.fast_build_set(word, t_den, M.PATTERN_E, shift)
    groups = T.partition_two_pointer(guarded, atom_intervals)
    events: dict[int, int] = {}
    for atom_index, intervals in enumerate(groups):
        sheet, _chamber = T.ATOMS[atom_index]
        bit = 1 << sheet
        for interval_left, interval_right in intervals:
            y_left = P * interval_left - sheet * t_den
            y_right = P * interval_right - sheet * t_den
            require(0 <= y_left < y_right <= t_den,
                    ("direct desheet interval", atom_index, y_left, y_right))
            events[y_left] = events.get(y_left, 0) ^ bit
            events[y_right] = events.get(y_right, 0) ^ bit
    boundaries = {0, t_den}
    boundaries.update(index * (t_den // P) for index in range(P + 1))
    boundaries.update(index * (t_den // (2 * Q)) for index in range(2 * Q + 1))
    for boundary in boundaries:
        events.setdefault(boundary, 0)

    full = [0 for _ell in range(Q)]
    same = [0 for _ell in range(Q)]
    mask = 0
    positions = sorted(events)
    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        count = mask.bit_count()
        jump = q_endpoint_jump(left, right)
        if not jump:
            continue
        ell = cell_of_segment(left, right, t_den)
        full[ell] = (full[ell] + count * count * jump) % prime
        same[ell] = (same[ell] + count * jump) % prime
    mask ^= events[positions[-1]]
    require(mask == 0, ("direct guarded profile did not close", alpha, beta, tau))
    phase = pow(worker_context()[5], beta, prime)
    return (
        tuple(phase * value % prime for value in full),
        tuple(phase * value % prime for value in same),
        len(guarded),
    )


def worker(alpha: int):
    (
        _word,
        _t_den,
        _nn,
        prime,
        _root,
        zeta,
        _q,
        _qs,
        _df,
        _dt,
        _atoms,
    ) = worker_context()
    full_rows = []
    same_rows = []
    counts = [0, 0, 0, 0]
    for beta in range(P):
        full, same, *record = desheet_profile(alpha, beta)
        phase = pow(zeta, beta, prime)
        full_rows.append(tuple(
            tuple(phase * value % prime for value in row) for row in full
        ))
        same_rows.append(tuple(
            tuple(phase * value % prime for value in row) for row in same
        ))
        counts = [left + right for left, right in zip(counts, record)]
    return alpha, tuple(full_rows), tuple(same_rows), tuple(counts)


def inverse_cell_table(gamma_cells, zeta: int, prime: int):
    """Invert only the thirteen q=(1,0,t) refined residue classes."""
    normalizer = pow(P**3, -1, prime)
    table = [[0 for _t in range(P)] for _ell in range(Q)]
    index = 0
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                row = gamma_cells[index]
                index += 1
                for t in range(P):
                    phase = pow(zeta, -(alpha + tau * t) % P, prime)
                    for ell in range(Q):
                        table[ell][t] = (
                            table[ell][t] + row[ell] * phase
                        ) % prime
    require(index == P**3, ("gamma count", index))
    return tuple(
        tuple(value * normalizer % prime for value in row) for row in table
    )


def fourier_2d(matrix, eta: int, zeta: int, prime: int):
    return tuple(
        tuple(
            sum(
                matrix[ell][t]
                * pow(eta, -h * ell % Q, prime)
                * pow(zeta, -k * t % P, prime)
                for ell in range(Q) for t in range(P)
            ) % prime
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


def matrix_interaction(matrix, prime: int):
    inv_q = pow(Q, -1, prime)
    inv_p = pow(P, -1, prime)
    inv_qp = pow(Q * P, -1, prime)
    row_sums = tuple(sum(row) % prime for row in matrix)
    col_sums = tuple(
        sum(matrix[ell][t] for ell in range(Q)) % prime for t in range(P)
    )
    grand = sum(row_sums) % prime
    answer = tuple(
        tuple((matrix[ell][t] - row_sums[ell] * inv_p
               - col_sums[t] * inv_q + grand * inv_qp) % prime
              for t in range(P))
        for ell in range(Q)
    )
    require(all(sum(row) % prime == 0 for row in answer), "interaction rows")
    require(all(sum(answer[ell][t] for ell in range(Q)) % prime == 0
                for t in range(P)), "interaction columns")
    return answer


def inverse_value_from_table(table, t: int) -> int:
    return sum(table[ell][t] for ell in range(Q)) % worker_context()[3]


def main() -> None:
    (
        word,
        t_den,
        nn,
        prime,
        root,
        zeta,
        _q_intervals,
        _q_starts,
        _diagonal_frequency,
        _diagonal_tabs,
        _atom_intervals,
    ) = worker_context()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")

    full_cells = tuple(
        row
        for _alpha, full_rows, _same_rows, _counts in chunks
        for beta_rows in full_rows
        for row in beta_rows
    )
    same_cells = tuple(
        row
        for _alpha, _full_rows, same_rows, _counts in chunks
        for beta_rows in same_rows
        for row in beta_rows
    )
    require(len(full_cells) == len(same_cells) == P**3, "character cell bank size")

    direct_controls = ((0, 0, 0), (0, 1, 6), (1, 0, 6), (6, 6, 12), (12, 12, 0))
    direct_record = []
    for alpha, beta, tau in direct_controls:
        direct_full, direct_same, interval_count = direct_guarded_profile(
            alpha, beta, tau
        )
        index = (alpha * P + beta) * P + tau
        require(direct_full == full_cells[index],
                ("direct full guard restoration", alpha, beta, tau))
        require(direct_same == same_cells[index],
                ("direct same guard restoration", alpha, beta, tau))
        direct_record.append(
            ((alpha, beta, tau), direct_full, direct_same, interval_count)
        )
    full_gamma = tuple(sum(row) % prime for row in full_cells)
    same_gamma = tuple(sum(row) % prime for row in same_cells)
    cross_cells = tuple(
        tuple((full - same) % prime for full, same in zip(full_row, same_row))
        for full_row, same_row in zip(full_cells, same_cells)
    )
    cross_gamma = tuple(sum(row) % prime for row in cross_cells)

    same_hash = digest_integers(same_gamma)
    require(same_hash == EXPECTED_POINT_GAMMA_SHA256,
            ("same-sheet point recovery", same_hash))
    require(all((same + cross) % prime == full
                for same, cross, full in zip(same_gamma, cross_gamma, full_gamma)),
            "same plus cross reconstruction")

    full_table = inverse_cell_table(full_cells, zeta, prime)
    same_table = inverse_cell_table(same_cells, zeta, prime)
    cross_table = inverse_cell_table(cross_cells, zeta, prime)
    require(all((same_table[ell][t] + cross_table[ell][t]) % prime
                == full_table[ell][t]
                for ell in range(Q) for t in range(P)),
            "inverse cell reconstruction")

    same_h = inverse_value_from_table(same_table, 1)
    same_q5 = inverse_value_from_table(same_table, 0)
    same_bridge = (same_h - same_q5) % prime
    require((same_h, same_q5, same_bridge) == EXPECTED_POINT_VALUES,
            ("point inverse recovery", same_h, same_q5, same_bridge))

    full_h = inverse_value_from_table(full_table, 1)
    full_q5 = inverse_value_from_table(full_table, 0)
    full_bridge = (full_h - full_q5) % prime
    cross_h = inverse_value_from_table(cross_table, 1)
    cross_q5 = inverse_value_from_table(cross_table, 0)
    cross_bridge = (cross_h - cross_q5) % prime
    require((same_bridge + cross_bridge) % prime == full_bridge,
            "bridge sector reconstruction")
    require((full_h, full_q5, full_bridge) == EXPECTED_FULL_VALUES,
            ("full common-residual values", full_h, full_q5, full_bridge))
    require((cross_h, cross_q5, cross_bridge) == EXPECTED_CROSS_VALUES,
            ("cross-sheet values", cross_h, cross_q5, cross_bridge))

    eta = pow(root, nn // Q, prime)
    require(pow(eta, Q, prime) == 1 and eta != 1, "order-seven root")
    full_spectrum = fourier_2d(full_table, eta, zeta, prime)
    same_spectrum = fourier_2d(same_table, eta, zeta, prime)
    cross_spectrum = fourier_2d(cross_table, eta, zeta, prime)
    full_interaction = matrix_interaction(full_table, prime)
    full_interaction_spectrum = fourier_2d(full_interaction, eta, zeta, prime)

    relation_t = 6
    relation_cell_profile = tuple(full_table[ell][relation_t] for ell in range(Q))
    relation_cell_spectrum = tuple(
        sum(relation_cell_profile[ell] * pow(eta, -h * ell % Q, prime)
            for ell in range(Q)) % prime
        for h in range(Q)
    )
    require(relation_cell_profile == EXPECTED_RELATION_CELL_PROFILE,
            ("fixed relation cell profile", relation_cell_profile))
    require(all(value == EXPECTED_RELATION_CELL_PROFILE[0]
                for value in relation_cell_spectrum),
            ("fixed relation F7 spectrum", relation_cell_spectrum))

    counts = tuple(
        sum(chunk[3][index] for chunk in chunks) for index in range(4)
    )
    require(counts == EXPECTED_WORK_COUNTS, ("work counts", counts))
    shapes = (
        support_shape(full_spectrum),
        support_shape(same_spectrum),
        support_shape(cross_spectrum),
    )
    require(shapes == ((91, 1, 6, 12, 72),) * 3,
            ("common-residual spectral closure", shapes))
    interaction_shape = support_shape(full_interaction_spectrum)
    require(interaction_shape == (72, 0, 0, 0, 72),
            ("common-residual mixed spectrum", interaction_shape))
    full_gamma_digest = digest_integers(full_gamma)
    cross_gamma_digest = digest_integers(cross_gamma)
    table_digests = (
        digest_json(full_table), digest_json(same_table), digest_json(cross_table)
    )
    spectrum_digests = (
        digest_json(full_spectrum),
        digest_json(same_spectrum),
        digest_json(cross_spectrum),
    )
    require(full_gamma_digest == EXPECTED_FULL_GAMMA_SHA256,
            ("full gamma", full_gamma_digest))
    require(cross_gamma_digest == EXPECTED_CROSS_GAMMA_SHA256,
            ("cross gamma", cross_gamma_digest))
    require(table_digests == EXPECTED_TABLE_SHA256,
            ("table digests", table_digests))
    require(spectrum_digests == EXPECTED_SPECTRUM_SHA256,
            ("spectrum digests", spectrum_digests))
    record = (
        TARGET_SHA256,
        POINT_SHA256,
        word,
        t_den,
        nn,
        prime,
        root,
        eta,
        zeta,
        counts,
        tuple(direct_record),
        digest_json(full_cells),
        same_hash,
        cross_gamma_digest,
        full_gamma_digest,
        table_digests[0],
        table_digests[1],
        table_digests[2],
        full_h,
        full_q5,
        full_bridge,
        same_h,
        same_q5,
        same_bridge,
        cross_h,
        cross_q5,
        cross_bridge,
        shapes[0],
        shapes[1],
        shapes[2],
        interaction_shape,
        spectrum_digests[0],
        spectrum_digests[1],
        spectrum_digests[2],
        relation_t,
        relation_cell_profile,
        relation_cell_spectrum,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== U_full desheeted common-residual base probe ==")
    print(f"dependencies=((endpoint,{TARGET_SHA256}),(point_hostile,{POINT_SHA256}))")
    print(f"map=91t=7a+r -> (sheet=a,base=y=r/7=13t-a); t=(y+a)/13")
    print(f"frequency_descent=(target_B={word[M.TARGET_B]}=13^2*{word[M.TARGET_B] // (P*P)},base_frequency={word[M.TARGET_B] // P})")
    print(f"split_embedding=(prime={prime},zeta7={eta},zeta13={zeta})")
    print(f"work_counts=(unguarded_intervals,mapped_intervals,active_segments,q_active_segments)={counts}")
    print(f"literal_guarded_profile_controls={direct_controls}: PASS")
    print(f"same_sheet_point_gamma_sha256={same_hash}; inherited_point_recovery=PASS")
    print(f"gamma_sha256=(full={full_gamma_digest},cross_sheet={cross_gamma_digest})")
    print(f"same_sheet_inverse=((q_H,{same_h}),(q_q5,{same_q5}),bridge={same_bridge})")
    print(f"cross_sheet_inverse=((q_H,{cross_h}),(q_q5,{cross_q5}),bridge={cross_bridge})")
    print(f"full_common_residual_inverse=((q_H,{full_h}),(q_q5,{full_q5}),bridge={full_bridge})")
    print(f"sector_reconstruction=(same+cross=full): PASS")
    print(f"cell_x_refined_residue_spectrum_shapes_(total,dc,F7axis,F13axis,mixed)=(full={shapes[0]},same={shapes[1]},cross={shapes[2]})")
    print(f"full_cell_residue_ANOVA_spectrum_shape={interaction_shape}")
    print(f"fixed_relation_class=(1,0,{relation_t}); seven_cell_profile={relation_cell_profile}")
    print(f"fixed_relation_F7_spectrum={relation_cell_spectrum}")
    print(f"table_sha256=(full={table_digests[0]},same={table_digests[1]},cross={table_digests[2]})")
    print(f"spectrum_sha256=(full={spectrum_digests[0]},same={spectrum_digests[1]},cross={spectrum_digests[2]})")
    print(f"semantic_sha256={semantic}")
    print("scope=actual U_full endpoint integrands on one desheeted base with seven cells; no THM-2471 collision/horizons/source packet,exact C(a;X,m),chronology,row exclusion,LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
