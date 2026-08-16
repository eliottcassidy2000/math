#!/usr/bin/env python3
"""Retain a third inverse-cylinder digit on the pointed-six response bundle.

On the fixed THM-2471 ``r=5`` sheet this script writes

    a = r0 + 13*r1 + 169*r2 + 2197*c

and retains ``(point=(Boolean state,U-root),r0,r1,r2,s=u-q,t)``.  The
expensive endpoint transform is factored once through exact refined source
cells: for every ``(cell,u,s)`` we compute its complete alpha/beta/tau
character response, then stream all 2,197 address fibres through that sparse
operator.  This is algebraically the same endpoint integration as the
two-digit parent, but avoids a dense 6*13^5 character bank.

The mandatory gate is exact marginalization over r2 to the hash-pinned
two-digit pointed tensor.  For every prior address (r0,r1), the script then
solves the unique restriction on its live pointed lines; a full Mat_6 map is
unique precisely when all six parent lines survive.  It tests diagonal
preservation of the six directed P4 arc lines, projective normalization,
stationary/lower-memory ansatzes, the adjacent-pair cocycle, and the
corresponding 78-state (six arcs times thirteen digits) de Bruijn path law.
Support-normalized and flat parent-preserving splits are hostile controls.

Everything here is a static finite response identity.  No chronology,
complete address, arrival/current typing, ancestry support map, row exclusion,
physical current, or LRC(14) consequence is asserted.
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
PARENT_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_two_current_digits_pointed_root_difference_carrier_transition_probe_20260816.py"
)
P4_SIDECAR_PATH = (
    ROOT / "04-computation/"
    "lrc_r5_pointed_bidirected_p4_cocycle_sidecar_20260816.py"
)
P4_SIDECAR_RAW_SHA256 = (
    "8b03b8fcf146c33d9203da68572fcc477f32f0fae2042f7c51c0f570a6bc3583"
)
P4_SIDECAR_SEMANTIC_SHA256 = (
    "b2ba313f88fbab0d36e95a63ade832492743ed6007fa0480434083d7dab0ecd3"
)
# Pin the literal tracked bytes.  The earlier reflection printed the stale
# value below as an "LF SHA"; the b1baa781a blob is already LF-only, so raw and
# LF-normalized hashes are both the actual 9d1671... value.
PARENT_RAW_SHA256 = "9d1671e0f823fdbaa9ab79915ba05dbb4dda4c6eabb97fe4484baf4c2e3205f2"
PARENT_LF_SHA256 = PARENT_RAW_SHA256
PARENT_STALE_DOCUMENTED_SHA256 = (
    "bc8727733804da38b9e7c691e2e9ff02de9d70d398916af3661816b9ae36c279"
)
PARENT_SEMANTIC_SHA256 = (
    "38725dc1d7129b326634c99bd70e1eb414590dc24fb83bd9522e2095e41f204c"
)
PARENT_TENSOR_SHA256 = (
    "f08ced17b1f727c8032a692e59df174de14ed06a9bc821e72788c8f347b28986"
)
PARENT_K2_MATRIX_DIGEST_SHA256 = (
    "185b2fb843d37a6e6f73be48375e40e93029ff507392df16f0058892184a1db2"
)

P = 13
V = 4
DIRECT_CONTROLS = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
ARC_PAIRS = ((0, 1), (2, 3), (4, 5))

# Pinned after the first exact normal replay.
EXPECTED_SOURCE_SHA256 = "c7e08cb20c445801bfab7ed87cd87a02084bb577323e0afe732f890f440ad77c"
EXPECTED_OPERATOR_SHA256 = "651e4a7ce39041fadfc89e56613ba976daea4b3ef2d67b5886e91fed729c5455"
EXPECTED_TRANSITION_SHA256 = "c5e3a65c92ee52aeb10f99cfccd1960abba71f7b71b1fafa19fd27c42111416d"
EXPECTED_COCYCLE_SHA256 = "22f0774a722e7e7f16e577e70921f9cfc053bdd8bd6c455749a5db8a8a5ef5f0"
EXPECTED_SEMANTIC_SHA256 = "3d1527fb4ce4931680e50d7135b9d1129c1816e3a9158645523e2728ddc71ec2"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def raw_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_parent():
    raw = raw_sha256(PARENT_PATH)
    normalized = lf_sha256(PARENT_PATH)
    require((raw, normalized) == (PARENT_RAW_SHA256, PARENT_LF_SHA256),
            ("two-digit pointed parent drift", raw, normalized))
    spec = importlib.util.spec_from_file_location("third_digit_pointed_parent", PARENT_PATH)
    require(spec is not None and spec.loader is not None, "third-digit parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B2 = load_parent()
T2 = B2.T2
C = B2.C
B = B2.B
PRIME = B2.PRIME
POINTS = B2.POINTS
POINT_INDEX = B2.POINT_INDEX
STATE_POINTS = B2.STATE_POINTS

require(B2.EXPECTED_SEMANTIC_SHA256 == PARENT_SEMANTIC_SHA256,
        "two-digit pointed semantic drift")
require(raw_sha256(P4_SIDECAR_PATH) == P4_SIDECAR_RAW_SHA256,
        ("bidirected-P4 typing sidecar drift", raw_sha256(P4_SIDECAR_PATH),
         P4_SIDECAR_RAW_SHA256))
require(POINTS == ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12)),
        ("pointed arc order", POINTS))
require(tuple(tuple(POINTS[index][1] for index in pair) for pair in ARC_PAIRS)
        == ((0, 0), (6, 6), (12, 12)), "bidirected P4 arc-pair labels")


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
def three_digit_context():
    """Build exact half-open mod-2197 profiles and their pointed cells."""

    ctx = C.context()
    stage = C.SM
    source_grid = C.S.GRID
    joint_grid = C.JOINT_COORDINATE
    require(stage.DCOLL == P**5, ("collision depth", stage.DCOLL))

    e_intervals = stage.build_set(stage.PAT_E, stage.ZELL)
    q_intervals = stage.build_set(stage.PAT_QB, stage.ZELL)
    f_pieces = C.S.B.build_f_pieces(e_intervals, q_intervals)

    # Contravariant extraction: highest retained digit first.  Every window
    # is the half-open chart [j/p,(j+1)/p), as enforced by extract_window.
    base_fold = stage.weighted_fold(f_pieces, stage.DCOLL // P**3, source_grid)
    r2_windows = tuple(
        stage.extract_window(base_fold[0], base_fold[1], r2, P, source_grid)
        for r2 in range(P)
    )
    r1_windows = tuple(
        tuple(
            stage.extract_window(
                r2_windows[r2][0], r2_windows[r2][1], r1, P, source_grid
            )
            for r2 in range(P)
        )
        for r1 in range(P)
    )
    address_windows = tuple(
        tuple(
            tuple(
                stage.extract_window(
                    r1_windows[r1][r2][0], r1_windows[r1][r2][1],
                    r0, P, source_grid,
                )
                for r2 in range(P)
            )
            for r1 in range(P)
        )
        for r0 in range(P)
    )

    # Full half-open endpoint audit: nested extraction equals the one-shot
    # mod-2197 window, including every boundary convention.
    for r0 in range(P):
        for r1 in range(P):
            for r2 in range(P):
                combined = stage.extract_window(
                    base_fold[0], base_fold[1], r0 + P * r1 + P**2 * r2,
                    P**3, source_grid,
                )
                require(address_windows[r0][r1][r2] == combined,
                        ("nested/combined mod-2197 window", r0, r1, r2))

    # Root extraction is composed in one shot.  The equality with an explicit
    # final half-open root window is checked on all 13^4 profiles.
    root_raw = tuple(
        tuple(
            tuple(
                tuple(
                    stage.extract_window(
                        base_fold[0], base_fold[1],
                        root + P * r0 + P**2 * r1 + P**3 * r2,
                        P**4, source_grid,
                    )
                    for r2 in range(P)
                )
                for r1 in range(P)
            )
            for r0 in range(P)
        )
        for root in range(P)
    )
    for root in range(P):
        for r0 in range(P):
            for r1 in range(P):
                for r2 in range(P):
                    nested_root = stage.extract_window(
                        address_windows[r0][r1][r2][0],
                        address_windows[r0][r1][r2][1],
                        root, P, source_grid,
                    )
                    require(root_raw[root][r0][r1][r2] == nested_root,
                            ("root/nested mod-13^4 window", root, r0, r1, r2))

    root_profiles = tuple(
        tuple(
            tuple(
                tuple(C.scale_profile(profile, joint_grid, source_grid)
                      for profile in r2_rows)
                for r2_rows in r1_rows
            )
            for r1_rows in r0_rows
        )
        for r0_rows in root_raw
    )
    parent_profiles, parent_boundaries, _parent_cells, _parent_record = (
        T2.two_digit_context()
    )
    boundaries = tuple(sorted(
        {0, joint_grid}
        | set(ctx["source_boundaries"])
        | set(parent_boundaries)
        | {
            position
            for root_rows in root_profiles
            for r0_rows in root_rows
            for r1_rows in r0_rows
            for profile in r1_rows
            for position in profile[0]
        }
    ))
    require(tuple(joint_grid - point for point in reversed(boundaries)) == boundaries,
            "three-digit source-boundary reflection")

    cells = []
    active_r2_counts = set()
    full_address_support = [set() for _state in range(V)]
    for left, right in zip(boundaries, boundaries[1:]):
        address_values = tuple(
            tuple(
                tuple(
                    tuple(profile_value(root_profiles[root][r0][r1][r2], left)
                          for r2 in range(P))
                    for r1 in range(P)
                )
                for r0 in range(P)
            )
            for root in range(P)
        )
        u_values = tuple(profile_value(profile, left) for profile in ctx["source_u"])
        v_values = tuple(profile_value(profile, left) for profile in ctx["source_v"])

        for root in range(P):
            for r0 in range(P):
                for r1 in range(P):
                    child_sum = sum(address_values[root][r0][r1])
                    expected = profile_value(parent_profiles[root][r0][r1], left)
                    require(child_sum == expected,
                            ("r2 source-cell marginal", left, right, root,
                             r0, r1, child_sum, expected))
                    if child_sum:
                        active_r2_counts.add(sum(value != 0 for value in
                                                 address_values[root][r0][r1]))
            require(
                sum(address_values[root][r0][r1][r2]
                    for r0 in range(P) for r1 in range(P) for r2 in range(P))
                == u_values[root],
                ("mod-2197 total marginal", left, right, root),
            )
            require(
                all(address_values[root][r0][r1][r2] * v_values[root] == 0
                    for r0 in range(P) for r1 in range(P) for r2 in range(P)),
                ("three-digit same-root source overlap", left, right, root),
            )
            for r0 in range(P):
                for r1 in range(P):
                    for r2 in range(P):
                        reflected = profile_value(
                            root_profiles[P - 1 - root][P - 1 - r0]
                                         [P - 1 - r1][P - 1 - r2],
                            joint_grid - right,
                        )
                        require(address_values[root][r0][r1][r2] == reflected,
                                ("three-digit profile reflection", left, right,
                                 root, r0, r1, r2))

        state = None
        try:
            state = B.state_of_segment(left, right)
        except RuntimeError:
            pass
        if state is not None:
            reflected_state = B.state_of_segment(joint_grid - right,
                                                 joint_grid - left)
            require(reflected_state == state ^ 2,
                    ("three-digit state reflection", left, right, state,
                     reflected_state))
            for root, value in enumerate(u_values):
                if not value:
                    continue
                require((state, root) in POINT_INDEX,
                        ("unsupported pointed source tail", left, state, root))
                full_address_support[state].add(sum(
                    address_values[root][r0][r1][r2] != 0
                    for r0 in range(P) for r1 in range(P) for r2 in range(P)
                ))

        u_support_mask = sum(1 << root for root, value in enumerate(u_values) if value)
        cells.append((state, u_values, v_values, u_support_mask, address_values))

    address_masses = tuple(
        sum(profile_mass(root_raw[root][r0][r1][r2], source_grid)
            for root in range(P))
        for r0 in range(P) for r1 in range(P) for r2 in range(P)
    )
    reflected_masses = tuple(
        address_masses[
            (P - 1 - r0) * P**2 + (P - 1 - r1) * P + (P - 1 - r2)
        ]
        for r0 in range(P) for r1 in range(P) for r2 in range(P)
    )
    require(address_masses == reflected_masses,
            "three-digit address-mass reflection")

    source_record = (
        stage.DCOLL,
        stage.DCOLL // P**3,
        len(f_pieces),
        len(base_fold[0]),
        sum(len(profile[0]) for profile in r2_windows),
        sum(len(profile[0]) for row in r1_windows for profile in row),
        sum(len(profile[0]) for slab in address_windows
            for row in slab for profile in row),
        sum(len(profile[0]) for root_rows in root_raw
            for r0_rows in root_rows for r1_rows in r0_rows for profile in r1_rows),
        len(boundaries),
        tuple(sorted(active_r2_counts)),
        tuple(tuple(sorted(values)) for values in full_address_support),
        address_masses,
        digest_json(root_raw),
        digest_json(root_profiles),
    )
    return root_profiles, boundaries, tuple(cells), source_record


@lru_cache(maxsize=64)
def operator_cells(boundaries):
    ctx = C.context()
    cells = []
    for left, right in zip(boundaries, boundaries[1:]):
        state = None
        try:
            state = B.state_of_segment(left, right)
        except RuntimeError:
            pass
        u_values = tuple(profile_value(profile, left) for profile in ctx["source_u"])
        v_values = tuple(profile_value(profile, left) for profile in ctx["source_v"])
        cells.append((state, u_values, v_values))
    return tuple(cells)


def endpoint_coefficients(alpha: int, beta: int, boundaries,
                          literal_tau: int | None = None):
    """Endpoint sweep before inserting any inverse-cylinder amplitudes."""

    cells = operator_cells(boundaries)
    events, interval_count, mapped = C.endpoint_events(alpha, beta, literal_tau)
    for boundary in boundaries:
        events.setdefault(boundary, 0)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    coefficients = [dict() for _tau in tau_values]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    state_segments = [0] * V
    updates = 0

    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        require(C.cell_of_segment(left, right) == 0,
                ("three-digit endpoint escaped cell zero", alpha, beta, left, right))
        chamber = C.chamber_of_segment(left, right)
        require(chamber in ("left", "right"),
                ("three-digit endpoint middle chamber", alpha, beta, left, right))
        cell_index = bisect_right(boundaries, left) - 1
        require(0 <= cell_index < len(cells)
                and right <= boundaries[cell_index + 1],
                ("endpoint crosses refined source cell", left, right, cell_index))
        state, u_values, v_values = cells[cell_index]
        require(state is not None, ("active untyped source cell", left, right))
        active_segments += 1
        state_segments[state] += 1

        jump = C.q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        u_support = tuple(root for root, value in enumerate(u_values) if value)
        v_support = tuple(root for root, value in enumerate(v_values) if value)
        if not u_support or not v_support:
            continue
        weighted_segments += 1
        require(set(u_support).isdisjoint(v_support),
                ("refined source supports overlap", left, right, u_support, v_support))

        for tau_index, tau in enumerate(tau_values):
            selected_mask = mask if literal_tau is not None else (
                mask & T2.SAFE_MASKS[chamber][tau]
            )
            selected_u = tuple(root for root in u_support
                               if (selected_mask >> root) & 1)
            selected_v = tuple(root for root in v_support
                               if (selected_mask >> root) & 1)
            for root in selected_u:
                require((state, root) in POINT_INDEX,
                        ("unrealized refined pointed source", left, state, root))
                for q_root in selected_v:
                    difference = (root - q_root) % P
                    require(difference != 0,
                            ("refined pointed same-root pair", left, root, q_root))
                    key = (cell_index, root, difference)
                    coefficients[tau_index][key] = (
                        coefficients[tau_index].get(key, 0)
                        + v_values[q_root] * jump
                    ) % PRIME
                    updates += 1

    mask ^= events[positions[-1]]
    require(mask == 0, ("three-digit endpoint mask", alpha, beta,
                        literal_tau, mask))
    counts = (
        interval_count, mapped, active_segments, q_active_segments,
        weighted_segments, tuple(state_segments), updates,
        tuple(len(row) for row in coefficients),
    )
    return tuple(coefficients), counts


def operator_worker(payload):
    alpha, boundaries = payload
    zeta = C.context()["zeta"]
    accum = [dict() for _tau in range(P)]
    scalar_counts = [0] * 5
    state_counts = [0] * V
    updates = 0
    key_counts = [0] * P
    row_digests = []
    for beta in range(P):
        coefficients, counts = endpoint_coefficients(alpha, beta, boundaries)
        phase = pow(zeta, (beta - alpha) % P, PRIME)
        row_digests.append(digest_json(tuple(
            tuple(sorted(row.items())) for row in coefficients
        )))
        for tau, row in enumerate(coefficients):
            target = accum[tau]
            for key, value in row.items():
                target[key] = (target.get(key, 0) + phase * value) % PRIME
        scalar_counts = [left + right for left, right in
                         zip(scalar_counts, counts[:5])]
        state_counts = [left + right for left, right in
                        zip(state_counts, counts[5])]
        updates += counts[6]
        key_counts = [left + right for left, right in
                      zip(key_counts, counts[7])]
    return (
        alpha,
        tuple(tuple(sorted(row.items())) for row in accum),
        digest_json(tuple(row_digests)),
        (tuple(scalar_counts), tuple(state_counts), updates, tuple(key_counts)),
    )


def build_sparse_operator(boundaries):
    zeta = C.context()["zeta"]
    tau_accum = [dict() for _tau in range(P)]
    alpha_digests = []
    work_records = []
    payloads = tuple((alpha, boundaries) for alpha in range(P))
    with ProcessPoolExecutor(max_workers=4) as pool:
        for expected_alpha, chunk in enumerate(pool.map(operator_worker, payloads)):
            alpha, rows, row_digest, work = chunk
            require(alpha == expected_alpha, ("operator worker order", alpha,
                                               expected_alpha))
            alpha_digests.append(row_digest)
            work_records.append(work)
            for tau, items in enumerate(rows):
                target = tau_accum[tau]
                for key, value in items:
                    target[key] = (target.get(key, 0) + value) % PRIME

    keys = tuple(sorted(set().union(*(set(row) for row in tau_accum))))
    roots = tuple(
        tuple(pow(zeta, -tau * relation % P, PRIME) for tau in range(P))
        for relation in range(P)
    )
    normalizer = pow(P**3, -1, PRIME)
    operator = {}
    for key in keys:
        tau_values = tuple(tau_accum[tau].get(key, 0) for tau in range(P))
        relation_row = tuple(
            sum(tau_values[tau] * roots[relation][tau] for tau in range(P))
            * normalizer % PRIME
            for relation in range(P)
        )
        if any(relation_row):
            operator[key] = relation_row
    record = (
        tuple(alpha_digests), tuple(work_records),
        tuple(len(row) for row in tau_accum), len(keys), len(operator),
        digest_json(tuple(tuple(sorted(row.items())) for row in tau_accum)),
        digest_json(tuple(sorted(operator.items()))),
    )
    return operator, record


def freeze_parent(bank):
    return tuple(
        tuple(
            tuple(tuple(tuple(line) for line in sheet) for sheet in block)
            for block in slab
        )
        for slab in bank
    )


def operator_index(operator):
    index = {}
    for (cell, root, difference), relation_row in operator.items():
        index.setdefault((cell, root), []).append((difference, relation_row))
    return {
        key: tuple(sorted(rows))
        for key, rows in index.items()
    }


def point_sources(cells):
    sources = [[] for _point in POINTS]
    for cell_index, (state, u_values, _v_values, _mask, addresses) in enumerate(cells):
        if state is None:
            continue
        for root, value in enumerate(u_values):
            if not value:
                continue
            point = POINT_INDEX.get((state, root))
            require(point is not None,
                    ("unindexed three-digit pointed source", cell_index, state, root))
            sources[point].append((cell_index, root, addresses[root]))
    return tuple(tuple(rows) for rows in sources)


def expand_parent_coefficients(coefficients, cells):
    """Insert sum_r2 source amplitudes into raw pre-integration rows."""

    bank = [[[[[0 for _s in range(P)] for _r1 in range(P)]
              for _r0 in range(P)] for _point in POINTS]
            for _tau in coefficients]
    for tau_index, coefficient_map in enumerate(coefficients):
        for (cell_index, root, difference), scalar in coefficient_map.items():
            if not scalar:
                continue
            state, _u_values, _v_values, _mask, addresses = cells[cell_index]
            require(state is not None, ("untyped parent coefficient", cell_index))
            point = POINT_INDEX[(state, root)]
            for r0 in range(P):
                for r1 in range(P):
                    amplitude = sum(addresses[root][r0][r1]) % PRIME
                    if amplitude:
                        bank[tau_index][point][r0][r1][difference] = (
                            bank[tau_index][point][r0][r1][difference]
                            + amplitude * scalar
                        ) % PRIME
    return tuple(B2.freeze_joint_row(row) for row in bank)


def direct_parent_controls(boundaries, cells):
    records = []
    for alpha, beta, tau in DIRECT_CONTROLS:
        coefficients, counts = endpoint_coefficients(
            alpha, beta, boundaries, literal_tau=tau
        )
        rebuilt = expand_parent_coefficients(coefficients, cells)
        parent_banks, parent_counts = B2.integrate_joint(
            alpha, beta, literal_tau=tau
        )
        require(rebuilt == parent_banks[0],
                ("literal r2 parent marginal", alpha, beta, tau,
                 digest_json(rebuilt), digest_json(parent_banks[0])))
        records.append((
            (alpha, beta, tau), digest_json(rebuilt),
            counts, parent_counts,
        ))
    return tuple(records)


def assemble_parent_tensor(operator, cells):
    """Rebuild the entire b1baa781a tensor through the refined operator."""

    indexed = operator_index(operator)
    sources = point_sources(cells)
    tensor = [[[[[0 for _relation in range(P)] for _difference in range(P)]
                for _r1 in range(P)] for _r0 in range(P)]
              for _point in POINTS]
    for point, entries in enumerate(sources):
        for cell_index, root, addresses in entries:
            rows = indexed.get((cell_index, root), ())
            for r0 in range(P):
                for r1 in range(P):
                    amplitude = sum(addresses[r0][r1]) % PRIME
                    if not amplitude:
                        continue
                    for difference, relation_row in rows:
                        target = tensor[point][r0][r1][difference]
                        for relation, value in enumerate(relation_row):
                            target[relation] = (
                                target[relation] + amplitude * value
                            ) % PRIME
    frozen = freeze_parent(tensor)
    require(all(frozen[point][r0][r1][0][relation] == 0
                for point in range(len(POINTS))
                for r0 in range(P) for r1 in range(P)
                for relation in range(P)), "parent same-root line")
    return frozen, indexed, sources


def parent_rows_at(parent, r0: int, r1: int):
    return tuple(
        tuple(parent[point][r0][r1][difference][relation]
              for difference in range(P) for relation in range(P))
        for point in range(len(POINTS))
    )


def child_rows_at(indexed, sources, parent_rows, r0: int, r1: int,
                  r2: int, mode: str):
    if mode == "flat":
        inverse = pow(P, -1, PRIME)
        return tuple(tuple(value * inverse % PRIME for value in row)
                     for row in parent_rows)

    rows = [[0 for _column in range(P * P)] for _point in POINTS]
    for point, entries in enumerate(sources):
        target = rows[point]
        for cell_index, root, addresses in entries:
            address_row = addresses[r0][r1]
            if mode == "actual":
                amplitude = address_row[r2] % PRIME
            elif mode == "support":
                total = sum(address_row) % PRIME
                active = tuple(index for index, value in enumerate(address_row) if value)
                require(bool(active) == bool(total),
                        ("support/total mismatch", cell_index, root, r0, r1))
                amplitude = (
                    total * pow(len(active), -1, PRIME) % PRIME
                    if r2 in active else 0
                )
            else:
                raise RuntimeError(("unknown third-digit mode", mode))
            if not amplitude:
                continue
            for difference, relation_row in indexed.get((cell_index, root), ()):
                offset = difference * P
                for relation, value in enumerate(relation_row):
                    target[offset + relation] = (
                        target[offset + relation] + amplitude * value
                    ) % PRIME
    return tuple(tuple(row) for row in rows)


def prepare_parent_solver(parent_rows):
    require(B2.rank_rows(parent_rows) == len(POINTS),
            ("prior-address carrier rank", B2.rank_rows(parent_rows)))
    pivots = B2.pivot_columns(parent_rows)
    require(len(pivots) == len(POINTS), ("prior-address pivots", pivots))
    minor = tuple(tuple(parent_rows[row][column] for column in pivots)
                  for row in range(len(POINTS)))
    return pivots, B2.matrix_inverse(minor)


def solve_prepared(parent_rows, child_rows, prepared):
    pivots, inverse = prepared
    child_minor = tuple(tuple(child_rows[row][column] for column in pivots)
                        for row in range(len(POINTS)))
    matrix = B2.matrix_multiply(child_minor, inverse)
    for row in range(len(POINTS)):
        for column in range(len(parent_rows[0])):
            expected = sum(matrix[row][middle] * parent_rows[middle][column]
                           for middle in range(len(POINTS))) % PRIME
            if expected != child_rows[row][column]:
                return matrix, (row, column, expected, child_rows[row][column])
    return matrix, None


def solve_live_lines(parent_rows, child_rows):
    """Solve the diagonal restriction on every nonzero pointed line.

    At a rank-deficient prior address, columns belonging to zero parent rows
    are unobservable, so a full Mat_6 map is not unique.  The returned matrix
    is the canonical minimal extension: the unique scalar on each live line
    and zero on each dead line.
    """

    size = len(POINTS)
    matrix = [[0 for _column in range(size)] for _row in range(size)]
    active = []
    for point in range(size):
        parent_row = parent_rows[point]
        child_row = child_rows[point]
        pivot = next((column for column, value in enumerate(parent_row) if value), None)
        if pivot is None:
            if any(child_row):
                return tuple(tuple(row) for row in matrix), tuple(active), (
                    "dead-parent-live-child", point,
                    next(column for column, value in enumerate(child_row) if value),
                )
            continue
        active.append(point)
        scalar = child_row[pivot] * pow(parent_row[pivot], -1, PRIME) % PRIME
        for column, value in enumerate(parent_row):
            if scalar * value % PRIME != child_row[column]:
                return tuple(tuple(row) for row in matrix), tuple(active), (
                    "line-not-proportional", point, column,
                    scalar * value % PRIME, child_row[column],
                )
        matrix[point][point] = scalar
    require(B2.rank_rows(parent_rows) == len(active),
            ("live-line rank mismatch", B2.rank_rows(parent_rows), active))
    return tuple(tuple(row) for row in matrix), tuple(active), None


def is_diagonal(matrix) -> bool:
    return all(matrix[row][column] == 0
               for row in range(len(POINTS)) for column in range(len(POINTS))
               if row != column)


def is_arc_pair_block(matrix) -> bool:
    block_of = {
        point: block
        for block, pair in enumerate(ARC_PAIRS)
        for point in pair
    }
    return all(matrix[row][column] == 0
               for row in range(len(POINTS)) for column in range(len(POINTS))
               if block_of[row] != block_of[column])


def identity_matrix():
    return tuple(tuple(int(row == column) for column in range(len(POINTS)))
                 for row in range(len(POINTS)))


def support_projector(active):
    active_set = set(active)
    return tuple(
        tuple(int(row == column and row in active_set)
              for column in range(len(POINTS)))
        for row in range(len(POINTS))
    )


def matrix_sum(matrices):
    return tuple(
        tuple(sum(matrix[row][column] for matrix in matrices) % PRIME
              for column in range(len(POINTS)))
        for row in range(len(POINTS))
    )


def histogram(values):
    counts = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return tuple(sorted(counts.items()))


def transition_bank(parent, indexed, sources, mode: str):
    matrices = [[[None for _r2 in range(P)] for _r1 in range(P)]
                for _r0 in range(P)]
    parent_ranks = []
    contained = full_unique = arc_block = diagonal = scalar = 0
    projective_identity = projective_support = 0
    first_nonunique = first_failure = first_off_diagonal = None
    rank_values = []
    support_values = []
    conditional_ranks = []
    global_r2_diagonals = [[] for _r2 in range(P)]

    for r0 in range(P):
        for r1 in range(P):
            parent_rows = parent_rows_at(parent, r0, r1)
            parent_rank = B2.rank_rows(parent_rows)
            parent_ranks.append(parent_rank)
            prepared = prepare_parent_solver(parent_rows) if parent_rank == len(POINTS) else None
            active_parent = tuple(point for point, row in enumerate(parent_rows) if any(row))
            require(len(active_parent) == parent_rank,
                    ("pointed parent lines/rank", r0, r1, active_parent, parent_rank))
            if parent_rank < len(POINTS) and first_nonunique is None:
                first_nonunique = (r0, r1, parent_rank, active_parent)
            address_matrices = []
            diagonal_rows = []
            for r2 in range(P):
                child_rows = child_rows_at(
                    indexed, sources, parent_rows, r0, r1, r2, mode
                )
                matrix, solved_active, witness = solve_live_lines(parent_rows, child_rows)
                require(solved_active == active_parent,
                        ("r2 live-line support", r0, r1, r2,
                         solved_active, active_parent))
                if prepared is not None:
                    full_matrix, full_witness = solve_prepared(
                        parent_rows, child_rows, prepared
                    )
                    require((full_witness, full_matrix) == (None, matrix),
                            ("full/live-line solver mismatch", r0, r1, r2,
                             full_witness, digest_json(full_matrix),
                             digest_json(matrix)))
                matrices[r0][r1][r2] = matrix
                if witness is None:
                    contained += 1
                    full_unique += int(parent_rank == len(POINTS))
                elif first_failure is None:
                    first_failure = (r0, r1, r2, witness,
                                     B2.rank_rows(parent_rows + child_rows))
                matrix_rank, matrix_support = B2.matrix_rank_and_support(matrix)
                rank_values.append(matrix_rank)
                support_values.append(matrix_support)
                if is_arc_pair_block(matrix):
                    arc_block += 1
                if is_diagonal(matrix):
                    diagonal += 1
                elif first_off_diagonal is None:
                    first_off_diagonal = (r0, r1, r2, digest_json(matrix))
                diagonal_row = tuple(matrix[index][index]
                                     for index in range(len(POINTS)))
                diagonal_rows.append(diagonal_row)
                global_r2_diagonals[r2].extend(diagonal_row)
                if len(set(diagonal_row)) == 1 and is_diagonal(matrix):
                    scalar += 1
                address_matrices.append(matrix)
            summed = matrix_sum(address_matrices)
            if summed == identity_matrix():
                projective_identity += 1
            if summed == support_projector(active_parent):
                projective_support += 1
            inverse = pow(P, -1, PRIME)
            mean = tuple(sum(row[column] for row in diagonal_rows) * inverse % PRIME
                         for column in range(len(POINTS)))
            contrasts = tuple(
                tuple((value - mean[column]) % PRIME
                      for column, value in enumerate(row))
                for row in diagonal_rows
            )
            conditional_ranks.append((B2.rank_rows(diagonal_rows),
                                      B2.rank_rows(contrasts)))

    frozen = tuple(tuple(tuple(row) for row in slab) for slab in matrices)
    matrix_digests = tuple(
        digest_json(frozen[r0][r1][r2])
        for r0 in range(P) for r1 in range(P) for r2 in range(P)
    )
    global_mean = tuple(
        sum(global_r2_diagonals[r2][column] for r2 in range(P))
        * pow(P, -1, PRIME) % PRIME
        for column in range(len(global_r2_diagonals[0]))
    )
    global_contrasts = tuple(
        tuple((value - global_mean[column]) % PRIME
              for column, value in enumerate(global_r2_diagonals[r2]))
        for r2 in range(P)
    )
    record = (
        mode,
        histogram(parent_ranks),
        contained, full_unique, arc_block, diagonal, scalar,
        projective_identity, projective_support,
        histogram(rank_values), histogram(support_values),
        histogram(conditional_ranks),
        (B2.rank_rows(tuple(tuple(row) for row in global_r2_diagonals)),
         B2.rank_rows(global_contrasts)),
        len(set(matrix_digests)),
        first_nonunique, first_failure, first_off_diagonal,
        digest_json(matrix_digests),
        digest_json(tuple(
            tuple(
                tuple(tuple(frozen[r0][r1][r2][point][point]
                            for point in range(len(POINTS)))
                      for r2 in range(P))
                for r1 in range(P)
            )
            for r0 in range(P)
        )),
    )
    return frozen, record


def first_unequal(items, labels):
    first = items[0]
    for index, item in enumerate(items[1:], 1):
        if item != first:
            return labels[0], labels[index], digest_json(first), digest_json(item)
    return None


def memory_record(bank, parent=None, full_rank_only: bool = False):
    def eligible(r0, r1):
        return (not full_rank_only
                or B2.rank_rows(parent_rows_at(parent, r0, r1)) == len(POINTS))

    fixed_r2_counts = []
    first_stationary = None
    for r2 in range(P):
        labels = [(r0, r1, r2) for r0 in range(P) for r1 in range(P)
                  if eligible(r0, r1)]
        rows = [bank[r0][r1][r2] for r0, r1, _r2 in labels]
        fixed_r2_counts.append(len({digest_json(row) for row in rows}))
        if first_stationary is None and rows:
            first_stationary = first_unequal(rows, labels)

    last_digit_counts = []
    first_last_digit = None
    for r1 in range(P):
        for r2 in range(P):
            labels = [(r0, r1, r2) for r0 in range(P) if eligible(r0, r1)]
            rows = [bank[r0][r1][r2] for r0, _r1, _r2 in labels]
            last_digit_counts.append(len({digest_json(row) for row in rows}))
            if first_last_digit is None and rows:
                first_last_digit = first_unequal(rows, labels)

    first_digit_counts = []
    first_first_digit = None
    for r0 in range(P):
        for r2 in range(P):
            labels = [(r0, r1, r2) for r1 in range(P) if eligible(r0, r1)]
            rows = [bank[r0][r1][r2] for _r0, r1, _r2 in labels]
            first_digit_counts.append(len({digest_json(row) for row in rows}))
            if first_first_digit is None and rows:
                first_first_digit = first_unequal(rows, labels)

    return (
        "full_rank_130" if full_rank_only else "canonical_minimal_169",
        tuple(fixed_r2_counts), all(count == 1 for count in fixed_r2_counts),
        histogram(last_digit_counts), all(count <= 1 for count in last_digit_counts),
        histogram(first_digit_counts), all(count <= 1 for count in first_digit_counts),
        first_stationary, first_last_digit, first_first_digit,
    )


def source_ratio_record(cells, mode: str, bank):
    sources = point_sources(cells)
    size_values = []
    dead = comparable = matches = exceptions = 0
    first_exception = first_mismatch = None
    for point, entries in enumerate(sources):
        for r0 in range(P):
            for r1 in range(P):
                for r2 in range(P):
                    ratios = set()
                    for _cell_index, _root, addresses in entries:
                        row = addresses[r0][r1]
                        total = sum(row) % PRIME
                        if not total:
                            continue
                        if mode == "actual":
                            child = row[r2] % PRIME
                        elif mode == "support":
                            active = tuple(index for index, value in enumerate(row)
                                           if value)
                            child = (total * pow(len(active), -1, PRIME) % PRIME
                                     if r2 in active else 0)
                        else:
                            raise RuntimeError(("source ratio mode", mode))
                        ratios.add(child * pow(total, -1, PRIME) % PRIME)
                    size_values.append(len(ratios))
                    if not ratios:
                        dead += 1
                    elif len(ratios) == 1:
                        comparable += 1
                        expected = next(iter(ratios))
                        observed = bank[r0][r1][r2][point][point]
                        if expected == observed:
                            matches += 1
                        elif first_mismatch is None:
                            first_mismatch = (point, r0, r1, r2, expected, observed)
                    else:
                        exceptions += 1
                        if first_exception is None:
                            first_exception = (
                                point, POINTS[point], r0, r1, r2, len(ratios),
                                digest_json(tuple(sorted(ratios))),
                            )
    return (
        mode, histogram(size_values), dead, comparable, matches, exceptions,
        first_exception, first_mismatch,
    )


def rebuild_k2(parent):
    """Recover the b1baa781a conditional maps without importing its tensor."""

    one_digit = [[[[0 for _relation in range(P)] for _difference in range(P)]
                  for _r0 in range(P)] for _point in POINTS]
    for point in range(len(POINTS)):
        for r0 in range(P):
            for difference in range(P):
                for relation in range(P):
                    one_digit[point][r0][difference][relation] = sum(
                        parent[point][r0][r1][difference][relation]
                        for r1 in range(P)
                    ) % PRIME

    bank = [[None for _r1 in range(P)] for _r0 in range(P)]
    projective = diagonal = arc_block = 0
    first_failure = None
    for r0 in range(P):
        base_rows = tuple(
            tuple(one_digit[point][r0][difference][relation]
                  for difference in range(P) for relation in range(P))
            for point in range(len(POINTS))
        )
        prepared = prepare_parent_solver(base_rows)
        matrices = []
        for r1 in range(P):
            child_rows = parent_rows_at(parent, r0, r1)
            matrix, witness = solve_prepared(base_rows, child_rows, prepared)
            if witness is not None and first_failure is None:
                first_failure = (r0, r1, witness)
            bank[r0][r1] = matrix
            matrices.append(matrix)
            diagonal += int(is_diagonal(matrix))
            arc_block += int(is_arc_pair_block(matrix))
        projective += int(matrix_sum(matrices) == identity_matrix())

    frozen = tuple(tuple(row) for row in bank)
    matrix_digest = digest_json(tuple(
        digest_json(frozen[r0][r1])
        for r0 in range(P) for r1 in range(P)
    ))
    require(matrix_digest == PARENT_K2_MATRIX_DIGEST_SHA256,
            ("reconstructed K2 matrix bank", matrix_digest,
             PARENT_K2_MATRIX_DIGEST_SHA256))
    record = (
        diagonal, arc_block, projective, first_failure,
        len({digest_json(frozen[r0][r1])
             for r0 in range(P) for r1 in range(P)}),
        matrix_digest,
    )
    return frozen, record


def matrix_difference(left, right):
    for row in range(len(left)):
        for column in range(len(left[row])):
            if left[row][column] != right[row][column]:
                return row, column, left[row][column], right[row][column]
    return None


def cocycle_record(k2, k3):
    """Test literal conditional and type-correct cumulative path laws."""

    literal_matches = cumulative_matches = 0
    literal_by_parent_rank = {}
    cumulative_by_parent_rank = {}
    first_literal = first_cumulative = None
    path_matches_by_arc = [0] * len(POINTS)
    for r0 in range(P):
        for r1 in range(P):
            parent_rank = B2.rank_rows(k2[r0][r1])
            for r2 in range(P):
                first_edge = k2[r0][r1]
                second_edge = k2[r1][r2]
                predicted = B2.matrix_multiply(second_edge, first_edge)
                conditional = k3[r0][r1][r2]
                if conditional == predicted:
                    literal_matches += 1
                    literal_by_parent_rank[parent_rank] = (
                        literal_by_parent_rank.get(parent_rank, 0) + 1
                    )
                elif first_literal is None:
                    first_literal = (
                        r0, r1, r2, matrix_difference(conditional, predicted),
                        digest_json(conditional), digest_json(predicted),
                    )

                cumulative = B2.matrix_multiply(conditional, first_edge)
                if cumulative == predicted:
                    cumulative_matches += 1
                    cumulative_by_parent_rank[parent_rank] = (
                        cumulative_by_parent_rank.get(parent_rank, 0) + 1
                    )
                elif first_cumulative is None:
                    first_cumulative = (
                        r0, r1, r2, matrix_difference(cumulative, predicted),
                        digest_json(cumulative), digest_json(predicted),
                    )
                for point in range(len(POINTS)):
                    if cumulative[point][point] == predicted[point][point]:
                        path_matches_by_arc[point] += 1

    # The literal product has the wrong conditional normalization in general:
    # summing over r2 returns K2(r0,r1), not I.  Verify that type invoice.
    product_sum_to_k2 = product_sum_to_identity = 0
    for r0 in range(P):
        for r1 in range(P):
            products = tuple(
                B2.matrix_multiply(k2[r1][r2], k2[r0][r1])
                for r2 in range(P)
            )
            summed = matrix_sum(products)
            product_sum_to_k2 += int(summed == k2[r0][r1])
            product_sum_to_identity += int(summed == identity_matrix())

    edge_table = tuple(
        tuple(tuple(k2[r0][r1][point][point]
                    for r1 in range(P))
              for r0 in range(P))
        for point in range(len(POINTS))
    )
    return (
        literal_matches, P**3, first_literal,
        cumulative_matches, P**3, first_cumulative,
        tuple(sorted(literal_by_parent_rank.items())),
        tuple(sorted(cumulative_by_parent_rank.items())),
        tuple(path_matches_by_arc),
        product_sum_to_k2, product_sum_to_identity,
        digest_json(edge_table),
    )


def arc_digit_rows(k2, k3=None, fixed_r2=None):
    """Compress exact trajectories to six arc blocks times thirteen digits.

    Columns are the same 78 arc-by-r0 coordinates.  The compression is exact
    because the six parent rows at every r0 are independent and all K maps
    used here are diagonal.
    """

    rows = []
    for point in range(len(POINTS)):
        for r1 in range(P):
            row = [0] * (len(POINTS) * P)
            for r0 in range(P):
                value = k2[r0][r1][point][point]
                if k3 is not None:
                    value = value * k3[r0][r1][fixed_r2][point][point] % PRIME
                row[point * P + r0] = value
            rows.append(tuple(row))
    return tuple(rows)


def block_preserving_78(matrix):
    return all(
        matrix[left][right] == 0
        for left in range(len(POINTS) * P)
        for right in range(len(POINTS) * P)
        if left // P != right // P
    )


def trajectory_78_record(k2, k3):
    """Test the unrestricted static 78-state lift separately from path law."""

    parent_rows = arc_digit_rows(k2)
    parent_rank = B2.rank_rows(parent_rows)
    child_ranks = []
    union_ranks = []
    matrices = []
    first_failure = None
    if parent_rank == len(parent_rows):
        parent_inverse = B2.matrix_inverse(parent_rows)
    else:
        parent_inverse = None

    for r2 in range(P):
        child_rows = arc_digit_rows(k2, k3, r2)
        child_rank = B2.rank_rows(child_rows)
        union_rank = B2.rank_rows(parent_rows + child_rows)
        child_ranks.append(child_rank)
        union_ranks.append(union_rank)
        if parent_inverse is not None:
            matrix = B2.matrix_multiply(child_rows, parent_inverse)
            require(B2.matrix_multiply(matrix, parent_rows) == child_rows,
                    ("78-state coordinate solve", r2))
            matrices.append(matrix)
        elif union_rank > parent_rank and first_failure is None:
            first_failure = (r2, parent_rank, child_rank, union_rank)

    all_unique = parent_inverse is not None and len(matrices) == P
    projective = False
    block_preserving = 0
    matrix_digest = "NONE"
    rank_support = ()
    if all_unique:
        projective = matrix_sum_78(matrices) == identity_matrix_78()
        block_preserving = sum(block_preserving_78(matrix) for matrix in matrices)
        matrix_digest = digest_json(tuple(digest_json(matrix) for matrix in matrices))
        rank_support = histogram(tuple(
            (B2.rank_rows(matrix), sum(value != 0 for row in matrix for value in row))
            for matrix in matrices
        ))
    return (
        parent_rank, tuple(child_ranks), tuple(union_ranks),
        all_unique, block_preserving, projective,
        len({digest_json(matrix) for matrix in matrices}),
        rank_support, first_failure, matrix_digest,
    )


def full_trajectory_rows(parent, k3=None, fixed_r2=None):
    rows = []
    for point in range(len(POINTS)):
        for r1 in range(P):
            row = []
            for r0 in range(P):
                scalar = (1 if k3 is None else
                          k3[r0][r1][fixed_r2][point][point])
                row.extend(
                    scalar * parent[point][r0][r1][difference][relation] % PRIME
                    for difference in range(P) for relation in range(P)
                )
            rows.append(tuple(row))
    return tuple(rows)


def trajectory_full_controls(parent, k3, compressed_record):
    parent_rows = full_trajectory_rows(parent)
    parent_rank = B2.rank_rows(parent_rows)
    require(parent_rank == compressed_record[0],
            ("full/compressed 78 parent rank", parent_rank,
             compressed_record[0]))
    records = []
    for r2 in (0, 6, 12):
        child_rows = full_trajectory_rows(parent, k3, r2)
        child_rank = B2.rank_rows(child_rows)
        union_rank = B2.rank_rows(parent_rows + child_rows)
        require((child_rank, union_rank)
                == (compressed_record[1][r2], compressed_record[2][r2]),
                ("full/compressed 78 child rank", r2, child_rank,
                 union_rank, compressed_record[1][r2],
                 compressed_record[2][r2]))
        records.append((r2, child_rank, union_rank))
    return parent_rank, tuple(records)


def identity_matrix_78():
    size = len(POINTS) * P
    return tuple(tuple(int(row == column) for column in range(size))
                 for row in range(size))


def matrix_sum_78(matrices):
    size = len(POINTS) * P
    return tuple(
        tuple(sum(matrix[row][column] for matrix in matrices) % PRIME
              for column in range(size))
        for row in range(size)
    )


def main() -> None:
    _profiles, boundaries, cells, source_record = three_digit_context()
    source_digest = digest_json((source_record, boundaries))
    if EXPECTED_SOURCE_SHA256 != "TO_BE_PINNED":
        require(source_digest == EXPECTED_SOURCE_SHA256,
                ("three-digit source digest", source_digest,
                 EXPECTED_SOURCE_SHA256))

    direct_records = direct_parent_controls(boundaries, cells)
    operator, operator_record = build_sparse_operator(boundaries)
    operator_digest = digest_json(operator_record)
    if EXPECTED_OPERATOR_SHA256 != "TO_BE_PINNED":
        require(operator_digest == EXPECTED_OPERATOR_SHA256,
                ("sparse operator digest", operator_digest,
                 EXPECTED_OPERATOR_SHA256))

    parent, indexed, sources = assemble_parent_tensor(operator, cells)
    parent_digest = digest_json(parent)
    require(parent_digest == PARENT_TENSOR_SHA256,
            ("r2 tensor marginal to b1baa781a", parent_digest,
             PARENT_TENSOR_SHA256))
    k2, k2_record = rebuild_k2(parent)

    mode_names = ("actual", "support", "flat")
    banks = []
    transition_records = []
    memory_records = []
    for mode in mode_names:
        bank, record = transition_bank(parent, indexed, sources, mode)
        banks.append(bank)
        transition_records.append(record)
        memory_records.append((
            memory_record(bank),
            memory_record(bank, parent, full_rank_only=True),
        ))

    actual, support, flat = banks
    expected_parent_ranks = ((0, 35), (2, 2), (4, 2), (6, 130))
    require(all(record[1] == expected_parent_ranks
                and record[2] == P**3
                and record[3] == 130 * P
                and record[4] == P**3
                and record[5] == P**3
                and record[7] == 130
                and record[8] == P**2
                for record in transition_records),
            ("three-digit live-line/projective gate", transition_records))
    inverse = pow(P, -1, PRIME)
    for r0 in range(P):
        for r1 in range(P):
            active = tuple(point for point, row in
                           enumerate(parent_rows_at(parent, r0, r1)) if any(row))
            expected_flat = tuple(
                tuple(inverse if row == column and row in active else 0
                      for column in range(len(POINTS)))
                for row in range(len(POINTS))
            )
            require(all(flat[r0][r1][r2] == expected_flat for r2 in range(P)),
                    ("flat K3 live-support projector", r0, r1, active))

    source_ratios = (
        source_ratio_record(cells, "actual", actual),
        source_ratio_record(cells, "support", support),
    )
    cocycle = cocycle_record(k2, actual)
    require(cocycle[0] == 35 * P and cocycle[3] == 35 * P
            and cocycle[6] == ((0, 35 * P),)
            and cocycle[7] == ((0, 35 * P),),
            ("adjacent-pair cocycle must be vacuous-only", cocycle))
    trajectory_78 = trajectory_78_record(k2, actual)
    full_trajectory_controls = trajectory_full_controls(
        parent, actual, trajectory_78
    )

    transition_digest = digest_json((
        tuple(transition_records), tuple(memory_records), source_ratios,
    ))
    cocycle_digest = digest_json((cocycle, trajectory_78,
                                  full_trajectory_controls, k2_record))
    if EXPECTED_TRANSITION_SHA256 != "TO_BE_PINNED":
        require(transition_digest == EXPECTED_TRANSITION_SHA256,
                ("three-digit transitions", transition_digest,
                 EXPECTED_TRANSITION_SHA256))
    if EXPECTED_COCYCLE_SHA256 != "TO_BE_PINNED":
        require(cocycle_digest == EXPECTED_COCYCLE_SHA256,
                ("three-digit cocycle", cocycle_digest,
                 EXPECTED_COCYCLE_SHA256))

    record = (
        PARENT_RAW_SHA256, PARENT_LF_SHA256, PARENT_STALE_DOCUMENTED_SHA256,
        PARENT_SEMANTIC_SHA256,
        P4_SIDECAR_RAW_SHA256, P4_SIDECAR_SEMANTIC_SHA256,
        PARENT_TENSOR_SHA256, PARENT_K2_MATRIX_DIGEST_SHA256,
        PRIME, C.JOINT_ROOT, C.context()["zeta"], POINTS, ARC_PAIRS,
        source_digest, direct_records, operator_record, parent_digest,
        k2_record, tuple(transition_records), tuple(memory_records),
        source_ratios, cocycle, trajectory_78, full_trajectory_controls,
        transition_digest, cocycle_digest,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("three-digit semantic", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 third inverse-cylinder digit x pointed-six diagonal bundle ==")
    print(f"parent=(commit=b1baa781a,raw_sha={PARENT_RAW_SHA256},lf_normalized_sha={PARENT_LF_SHA256},stale_reflection_sha={PARENT_STALE_DOCUMENTED_SHA256},semantic={PARENT_SEMANTIC_SHA256},tensor={PARENT_TENSOR_SHA256})")
    print(f"typing_parent=(bidirected_P4_sidecar_sha={P4_SIDECAR_RAW_SHA256},semantic={P4_SIDECAR_SEMANTIC_SHA256},P7_H1=0)")
    print(f"coordinates=(point=(state,u) in {POINTS},arc_pairs={ARC_PAIRS},r0,r1,r2,s=u-q!=0,relation=(1,0,t))")
    print("construction=a=r0+13*r1+169*r2+2197*c;contravariant half-open r2/r1/r0/root windows;refined-cell endpoint operator before address expansion")
    print(f"field=(prime={PRIME},root={C.JOINT_ROOT},zeta13={C.context()['zeta']})")
    print(f"source=(fold={source_record[1]},pieces=(f={source_record[2]},base={source_record[3]},r2={source_record[4]},r1={source_record[5]},address={source_record[6]},root={source_record[7]}),boundaries={source_record[8]},active_r2_counts={source_record[9]},full_support_by_state={source_record[10]},digest={source_digest})")
    print(f"literal_preintegration_parent_controls={direct_records}: PASS")
    print(f"sparse_operator=(record={operator_record},digest={operator_digest})")
    print(f"mandatory_r2_marginal=(source_cell_exact,tensor_sha={parent_digest},expected_b1baa781a={PARENT_TENSOR_SHA256}): PASS")
    print(f"reconstructed_two_digit_maps={k2_record}")
    print("transition_record=(mode,parent_rank_hist,linewise_contained,full_Mat6_unique,P4_arc_pair_block,canonical_diagonal,canonical_scalar,sum_r2_I6,sum_r2_live_projector,rank_hist,nnz_hist,conditional_amplitude_rank_hist,global_r2_amplitude_rank,distinct,first_nonunique,first_failure,first_offdiag,matrix_digest,diagonal_digest)")
    print(f"third_digit_transitions={tuple(transition_records)}")
    print("memory_record_per_mode=(canonical_minimal_all_169,unique_full_rank_130);each=(scope,unique_at_fixed_r2,stationary_K_r2,unique_across_r0_at_fixed_r1_r2,K_(r1,r2),unique_across_r1_at_fixed_r0_r2,K_(r0,r2),first_witnesses)")
    print(f"stationary_and_lower_memory={tuple(zip(mode_names, memory_records))}")
    print(f"source_cell_child_over_parent={source_ratios}")
    print("cocycle_record=(literal_conditional_matches,total,first;type_correct_cumulative_matches,total,first;literal_matches_by_parent_rank;cumulative_matches_by_parent_rank;per_arc_matches;product_sums_to_K2;product_sums_to_I;edge_table_digest)")
    print(f"adjacent_pair_cocycle={cocycle}")
    print("trajectory_78=(parent_rank,child_ranks,union_ranks,all_unique,arc_block_count,sum_M_r2_I,distinct,(rank,nnz)_hist,first_failure,matrix_digest)")
    print(f"static_78_state_lift={trajectory_78}")
    print(f"full_2197_column_trajectory_rank_controls={full_trajectory_controls}: PASS")
    print(f"verification_digests=(transition={transition_digest},cocycle={cocycle_digest})")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT static third-cylinder pointed response/transition probe on one owner base")
    print("typing=no chronology,no complete address,no arrival/current typing,no ancestry support map,no physical current,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
