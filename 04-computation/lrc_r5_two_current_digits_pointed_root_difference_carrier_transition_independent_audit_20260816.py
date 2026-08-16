#!/usr/bin/env python3
"""Independent hostile audit of the two-current pointed diagonal bundle.

The submitted five-coordinate candidate is never imported.  This audit uses
two independently audited parents: the exact mod-169 current-cylinder source
profiles and the six-pointed ordered-root endpoint engine.  On each endpoint
segment it aggregates literal ordered pairs (u,q) by exact source cell before
expanding the 169 address weights.  This is an independent pre-integration
route to

    T(point, r0, r1, s=u-q, relation).

The transition matrices below are static coordinate maps between finite
response tables.  They are not clocks, complete addresses, arrival atoms,
physical currents, row exclusions, or an LRC(14) argument.
"""

from __future__ import annotations

from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TWO_DIGIT_AUDIT_PATH = ROOT / (
    "04-computation/lrc_r5_ufull_owner_node_boolean_square_"
    "two_digit_current_ancestry_independent_audit_20260816.py"
)
POINTED_AUDIT_PATH = ROOT / (
    "04-computation/lrc_r5_ufull_owner_node_pointed_six_state_"
    "root_difference_independent_audit_20260816.py"
)
CANDIDATE_PATH = ROOT / (
    "04-computation/lrc_r5_two_current_digits_pointed_root_difference_"
    "carrier_transition_probe_20260816.py"
)

TWO_DIGIT_AUDIT_SHA256 = (
    "126be106a34a990f22e66b26d79ea3568ee4f394419936ed47f1bfbe0656788f"
)
POINTED_AUDIT_SHA256 = (
    "8c7cb5f98b15a768d4f4d6060074e0815a8f089f857ec4f3c55a0e7d877e1fec"
)
CANDIDATE_SHA256_OBSERVED = (
    "9d1671e0f823fdbaa9ab79915ba05dbb4dda4c6eabb97fe4484baf4c2e3205f2"
)
CANDIDATE_SHA256_PUBLISHED = (
    "bc8727733804da38b9e7c691e2e9ff02de9d70d398916af3661816b9ae36c279"
)

P = 13
V = 4
CONTROLS = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
INVERSION_PROBES = (
    (0, 0, 1, 1),
    (1, 6, 7, 6),
    (3, 4, 9, 7),
    (5, 12, 11, 12),
)
PROBE_RELATIONS = (0, 6, 12)

ONE_DIGIT_ROOT_GAMMA_SHA256 = (
    "843fe780896a23471f9e844383cc03cbdc70836a138a6addc218e87e319649de"
)
ONE_DIGIT_ROOT_TENSOR_SHA256 = (
    "ccfdb2373578190156dca9e9012c0a98ac50cdbacea84f1ff21d5ec0ac94db5b"
)
POINTED_GAMMA_SHA256 = (
    "416f03a90894e526bc767a34839c2489aa4cac7051ac04687e0feacfa36d58d2"
)
POINTED_TENSOR_SHA256 = (
    "9c5e227d9c142373973a562a54c6a67cac60a82da1121a028b9658920d155a19"
)
EXPECTED_CANDIDATE_TENSOR_SHA256 = (
    "f08ced17b1f727c8032a692e59df174de14ed06a9bc821e72788c8f347b28986",
    "6aa28175baa458ae8646e03b4b23234febb3a5ef60714ec7e35d623e4d2267a0",
    "e4f96fe95854d1f78649c9bbe71c91a5276ab29813934c90e11c6cd616bc7964",
    "e36cedbce618771be7ac5e94d74fa6be948ad8236c3145222877ee274a6b5412",
)

# Pinned after the first exact normal replay.
EXPECTED_SEMANTIC_SHA256 = (
    "b58d98c283a2dc42111365a3d61af0948757feda3b1ff11ca65cd7d562b15a56"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_module(path: Path, expected: str, name: str):
    observed = lf_sha256(path)
    require(observed == expected, (name, "source drift", observed, expected))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


TD = load_module(
    TWO_DIGIT_AUDIT_PATH, TWO_DIGIT_AUDIT_SHA256, "diagonal_audit_two_digit_parent"
)
PA = load_module(
    POINTED_AUDIT_PATH, POINTED_AUDIT_SHA256, "diagonal_audit_pointed_parent"
)

R = TD.R
SQ = TD.SQ
PRIME = TD.PRIME
POINTS = PA.POINTS
POINT_INDEX = PA.POINT_INDEX
STATE_FIBRES = PA.STATE_FIBRES

require(POINTS == ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12)),
        ("point order", POINTS))
require(PRIME == PA.PRIME, "split-field mismatch")
require(TD.EXPECTED_SEMANTIC_SHA256
        == "f2af726e9d5abd1487e841623ce2f62ca647c86e5a1a68e41eda4d9dda6c81ac",
        "two-digit parent semantic drift")
require(PA.EXPECTED_SEMANTIC_SHA256
        == "66db2301f88db1ced7784868095e198e3e12f1fb79175ebd902d0f569a5decef",
        "pointed parent semantic drift")
require(lf_sha256(CANDIDATE_PATH) == CANDIDATE_SHA256_OBSERVED,
        "submitted candidate source drift")


def freeze_row(bank):
    return tuple(
        tuple(tuple(tuple(line) for line in block) for block in slab)
        for slab in bank
    )


def freeze_tensor(bank):
    return tuple(
        tuple(
            tuple(tuple(tuple(line) for line in block) for block in slab)
            for slab in sheet
        )
        for sheet in bank
    )


def zero_row():
    return [[[[0 for _difference in range(P)] for _r1 in range(P)]
             for _r0 in range(P)] for _point in POINTS]


@lru_cache(maxsize=1)
def source_data():
    (
        profiles, boundaries, cells, one_digit, source_u, source_v,
        source_boundaries, source_profile_digest, source_record, source_sha,
        reflection_record,
    ) = TD.source_context()
    require(source_sha == TD.EXPECTED_SOURCE_SHA256,
            ("two-digit source parent", source_sha, TD.EXPECTED_SOURCE_SHA256))
    owner_cells = []
    owner_state_map = {
        ("left", frozenset((0,))): 0,
        ("left", frozenset((0, 6))): 1,
        ("right", frozenset((12,))): 2,
        ("right", frozenset((6, 12))): 3,
    }
    joint = R.JOINT_COORDINATE
    for left, right, cell in zip(boundaries, boundaries[1:], cells):
        support_state, u_values, v_values, support_mask, addresses = cell
        doubled_midpoint = left + right
        chamber = (
            "left" if 7 * doubled_midpoint < 2 * joint else
            "right" if 7 * doubled_midpoint >= 12 * joint else
            None
        )
        support = frozenset(root for root, value in enumerate(u_values) if value)
        owner_state = owner_state_map.get((chamber, support))
        if owner_state is not None:
            require(owner_state == support_state,
                    ("owner/support state agreement", left, right,
                     chamber, support, owner_state, support_state))
        owner_cells.append((owner_state, u_values, v_values, support_mask, addresses))
    return {
        "profiles": profiles,
        "boundaries": boundaries,
        "cells": tuple(owner_cells),
        "one_digit": one_digit,
        "source_u": source_u,
        "source_v": source_v,
        "source_boundaries": source_boundaries,
        "source_profile_digest": source_profile_digest,
        "source_record": source_record,
        "source_sha": source_sha,
        "reflection": reflection_record,
    }


def source_support_record(cells):
    sizes = set()
    patterns = set()
    for state, u_values, _v_values, _mask, addresses in cells:
        if state is None:
            continue
        for root, value in enumerate(u_values):
            if not value:
                continue
            pattern = tuple(
                sum(addresses[root][r0][r1] != 0 for r1 in range(P))
                for r0 in range(P)
            )
            sizes.update(pattern)
            patterns.add(pattern)
    return tuple(sorted(sizes)), tuple(sorted(patterns))


def split_value(addresses, root: int, r0: int, r1: int, normalized: bool) -> int:
    row = addresses[root][r0]
    if not normalized:
        return row[r1] % PRIME
    active = tuple(index for index, value in enumerate(row) if value)
    if r1 not in active:
        return 0
    parent = sum(row) % PRIME
    require(parent != 0 and active, ("support split", root, r0, r1))
    return parent * pow(len(active), -1, PRIME) % PRIME


def source_ratio_record(cells):
    all_kernels = []
    records = []
    for normalized in (False, True):
        kernels = [[[None for _r1 in range(P)] for _r0 in range(P)]
                   for _point in POINTS]
        size_histogram = {}
        exceptions = []
        for point, (state, root) in enumerate(POINTS):
            for r0 in range(P):
                for r1 in range(P):
                    ratios = set()
                    for cell_state, u_values, _v_values, _mask, addresses in cells:
                        if cell_state != state or not u_values[root]:
                            continue
                        parent = sum(addresses[root][r0]) % PRIME
                        child = split_value(addresses, root, r0, r1, normalized)
                        if not parent:
                            require(child == 0,
                                    ("child without source parent", point, r0, r1))
                            continue
                        ratios.add(child * pow(parent, -1, PRIME) % PRIME)
                    require(ratios, ("empty source ratio fibre", point, r0, r1))
                    size_histogram[len(ratios)] = (
                        size_histogram.get(len(ratios), 0) + 1
                    )
                    if len(ratios) == 1:
                        kernels[point][r0][r1] = next(iter(ratios))
                    else:
                        exceptions.append(
                            (point, POINTS[point], r0, r1, tuple(sorted(ratios)))
                        )
        compact_exceptions = tuple(
            (point, pair, r0, tuple(
                r1 for q_point, q_pair, q_r0, r1, _ratios in exceptions
                if (q_point, q_pair, q_r0) == (point, pair, r0)
            ))
            for point, pair, r0 in sorted(
                {(point, pair, r0) for point, pair, r0, _r1, _ratios in exceptions}
            )
        )
        record = (
            tuple(sorted(size_histogram.items())), len(exceptions),
            compact_exceptions, digest_json(tuple(exceptions)),
        )
        all_kernels.append(tuple(
            tuple(tuple(row) for row in point_rows) for point_rows in kernels
        ))
        records.append(record)
    return tuple(all_kernels), tuple(records)


def pair_coefficient_maps(alpha: int, beta: int, literal_tau: int | None = None):
    data = source_data()
    boundaries = data["boundaries"]
    cells = data["cells"]
    word, endpoint_grid, harmonic, danger = TD.endpoint_context()
    events, interval_count, mapped_count = R.endpoint_events(
        word, endpoint_grid, alpha, beta, literal_tau
    )
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    maps = [dict() for _tau in tau_values]
    positions = sorted(set(events) | set(boundaries))
    mask = 0
    active = q_active = weighted = ordered_updates = 0
    state_counts = [0] * V
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
                ("joint audit escaped cell zero", alpha, beta, left, right))
        chamber = R.chamber(left, right)
        require(chamber in ("left", "right"),
                ("joint audit entered middle chamber", chamber))
        cell_index = bisect_right(boundaries, left) - 1
        require(0 <= cell_index < len(cells)
                and right <= boundaries[cell_index + 1],
                ("joint audit source crossing", left, right, cell_index))
        state, u_values, v_values, _support_mask, addresses = cells[cell_index]
        require(state is not None, ("untyped active source cell", left, right))
        state_counts[state] += 1
        if not jump:
            continue
        q_active += 1
        u_support = tuple(root for root, value in enumerate(u_values) if value)
        v_support = tuple(root for root, value in enumerate(v_values) if value)
        require(u_support and v_support, ("empty weighted root side", left))
        require(set(u_support).isdisjoint(v_support),
                ("same-root source overlap", left, u_support, v_support))
        weighted += 1
        for row_index, tau in enumerate(tau_values):
            selected = mask if literal_tau is not None else (
                mask & R.guard_mask(chamber, tau, danger)
            )
            for u in u_support:
                if not ((selected >> u) & 1):
                    continue
                require(sum(addresses[u][r0][r1]
                            for r0 in range(P) for r1 in range(P))
                        == u_values[u],
                        ("address restoration", cell_index, u))
                for q in v_support:
                    if not ((selected >> q) & 1):
                        continue
                    require(u != q, ("same-root selected pair", left, u, q))
                    scalar = v_values[q] * jump % PRIME
                    key = (cell_index, u, q)
                    maps[row_index][key] = (
                        maps[row_index].get(key, 0) + scalar
                    ) % PRIME
                    ordered_updates += 1
    mask ^= events.get(positions[-1], 0)
    require(mask == 0, ("joint endpoint mask closure", alpha, beta, literal_tau))
    return tuple(maps), (
        interval_count, mapped_count, active, q_active, weighted,
        tuple(state_counts), ordered_updates, tuple(len(row) for row in maps),
    )


@lru_cache(maxsize=100000)
def expanded_pair_vector(cell_index: int, u: int, q: int):
    state, _u_values, _v_values, _mask, addresses = source_data()["cells"][cell_index]
    point = POINT_INDEX.get((state, u))
    require(point is not None, ("unrealized pointed root", cell_index, state, u))
    difference = (u - q) % P
    require(difference != 0, ("zero root difference", cell_index, u, q))
    entries = []
    for r0 in range(P):
        for r1 in range(P):
            actual = split_value(addresses, u, r0, r1, False)
            support = split_value(addresses, u, r0, r1, True)
            if actual or support:
                entries.append((point, r0, r1, difference, actual, support))
    return tuple(entries)


def expand_pair_map(coefficient_map):
    actual = zero_row()
    support = zero_row()
    for (cell_index, u, q), scalar in coefficient_map.items():
        if not scalar:
            continue
        for point, r0, r1, difference, left_actual, left_support in (
            expanded_pair_vector(cell_index, u, q)
        ):
            actual[point][r0][r1][difference] = (
                actual[point][r0][r1][difference] + left_actual * scalar
            ) % PRIME
            support[point][r0][r1][difference] = (
                support[point][r0][r1][difference] + left_support * scalar
            ) % PRIME
    return freeze_row(actual), freeze_row(support)


def scale_row(row, scalar: int):
    return tuple(
        tuple(
            tuple(tuple(scalar * value % PRIME for value in line) for line in block)
            for block in slab
        )
        for slab in row
    )


def row_two_digit_marginal(row):
    return tuple(
        tuple(tuple(
            sum(row[point][r0][r1][difference]
                for point in STATE_FIBRES[state] for difference in range(P)) % PRIME
            for r1 in range(P)
        ) for r0 in range(P))
        for state in range(V)
    )


def row_one_digit_pointed(row):
    return tuple(
        tuple(tuple(
            sum(row[point][r0][r1][difference] for r1 in range(P)) % PRIME
            for difference in range(P)
        ) for r0 in range(P))
        for point in range(len(POINTS))
    )


def row_one_digit_state(row):
    pointed = row_one_digit_pointed(row)
    return tuple(
        tuple(tuple(
            sum(pointed[point][r0][difference] for point in STATE_FIBRES[state])
            % PRIME
            for difference in range(P)
        ) for r0 in range(P))
        for state in range(V)
    )


def row_pointed_parent(row):
    return tuple(
        tuple(
            sum(row[point][r0][r1][difference]
                for r0 in range(P) for r1 in range(P)) % PRIME
            for difference in range(P)
        )
        for point in range(len(POINTS))
    )


def add_row_to_core(core, tau: int, row, scalar: int = 1):
    for point in range(len(POINTS)):
        for r0 in range(P):
            for r1 in range(P):
                for difference in range(P):
                    core[tau][point][r0][r1][difference] = (
                        core[tau][point][r0][r1][difference]
                        + scalar * row[point][r0][r1][difference]
                    ) % PRIME


def compact_json_rows(rows) -> str:
    return json.dumps(tuple(rows), separators=(",", ":"), sort_keys=True)


def worker(alpha: int):
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    tau_actual = [zero_row() for _tau in range(P)]
    tau_support = [zero_row() for _tau in range(P)]
    two_digit_blocks = []
    one_digit_rows = []
    pointed_rows = []
    row_digest_pairs = []
    controls = {}
    probe_values = [[] for _probe in INVERSION_PROBES]
    scalar_counts = [0] * 5
    state_counts = [0] * V
    ordered_updates = 0
    key_counts = [0] * P
    for beta in range(P):
        maps, counts = pair_coefficient_maps(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        two_digit_rows = []
        for tau in range(P):
            raw_actual, raw_support = expand_pair_map(maps[tau])
            actual = scale_row(raw_actual, phase)
            support = scale_row(raw_support, phase)
            require(row_one_digit_pointed(actual) == row_one_digit_pointed(support),
                    ("support-normalized gamma parent", alpha, beta, tau))
            two_digit_rows.append(row_two_digit_marginal(actual))
            one_digit_rows.append(row_one_digit_state(actual))
            pointed_rows.append(row_pointed_parent(actual))
            row_digest_pairs.append((digest_json(actual), digest_json(support)))
            add_row_to_core(tau_actual, tau, actual)
            add_row_to_core(tau_support, tau, support)
            if (alpha, beta, tau) in CONTROLS:
                controls[(alpha, beta, tau)] = (actual, support)
            for probe_index, (point, r0, r1, difference) in enumerate(INVERSION_PROBES):
                probe_values[probe_index].append(
                    actual[point][r0][r1][difference]
                )
        two_digit_blocks.append(tuple(two_digit_rows))
        scalar_counts = [left + right for left, right in zip(scalar_counts, counts[:5])]
        state_counts = [left + right for left, right in zip(state_counts, counts[5])]
        ordered_updates += counts[6]
        key_counts = [left + right for left, right in zip(key_counts, counts[7])]
    require(digest_json(tuple(two_digit_blocks)) == TD.EXPECTED_ALPHA_GAMMA_SHA256[alpha],
            ("two-digit alpha gamma marginal", alpha,
             digest_json(tuple(two_digit_blocks)), TD.EXPECTED_ALPHA_GAMMA_SHA256[alpha]))
    return (
        alpha,
        tuple(freeze_row(row) for row in tau_actual),
        tuple(freeze_row(row) for row in tau_support),
        compact_json_rows(one_digit_rows),
        compact_json_rows(pointed_rows),
        tuple(row_digest_pairs),
        tuple(sorted(controls.items())),
        tuple(tuple(values) for values in probe_values),
        (tuple(scalar_counts), tuple(state_counts), ordered_updates,
         tuple(key_counts)),
    )


def combine_json_blocks(blocks) -> str:
    pieces = []
    for block in blocks:
        require(block.startswith("[") and block.endswith("]"), "JSON block shape")
        if len(block) > 2:
            pieces.append(block[1:-1])
    return "[" + ",".join(pieces) + "]"


def inverse_core(tau_core, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    tensor = [[[[[0 for _relation in range(P)] for _difference in range(P)]
                for _r1 in range(P)] for _r0 in range(P)]
              for _point in POINTS]
    for point in range(len(POINTS)):
        for r0 in range(P):
            for r1 in range(P):
                for difference in range(P):
                    for relation in range(P):
                        tensor[point][r0][r1][difference][relation] = (
                            sum(
                                tau_core[tau][point][r0][r1][difference]
                                * pow(zeta, -tau * relation % P, PRIME)
                                for tau in range(P)
                            ) * normalizer
                        ) % PRIME
    return freeze_tensor(tensor)


def tensor_two_digit_marginal(tensor):
    return tuple(
        tuple(tuple(tuple(
            sum(tensor[point][r0][r1][difference][relation]
                for point in STATE_FIBRES[state] for difference in range(P)) % PRIME
            for relation in range(P)
        ) for r1 in range(P)) for r0 in range(P))
        for state in range(V)
    )


def tensor_one_digit_pointed(tensor):
    return tuple(
        tuple(tuple(tuple(
            sum(tensor[point][r0][r1][difference][relation]
                for r1 in range(P)) % PRIME
            for relation in range(P)
        ) for difference in range(P)) for r0 in range(P))
        for point in range(len(POINTS))
    )


def tensor_one_digit_state(pointed):
    return tuple(
        tuple(tuple(tuple(
            sum(pointed[point][r0][difference][relation]
                for point in STATE_FIBRES[state]) % PRIME
            for relation in range(P)
        ) for difference in range(P)) for r0 in range(P))
        for state in range(V)
    )


def tensor_pointed_parent(pointed):
    return tuple(
        tuple(tuple(
            sum(pointed[point][r0][difference][relation] for r0 in range(P))
            % PRIME
            for relation in range(P)
        ) for difference in range(P))
        for point in range(len(POINTS))
    )


def r1_flat_lift(parent):
    inverse = pow(P, -1, PRIME)
    return tuple(
        tuple(
            tuple(
                tuple(tuple(parent[point][r0][difference][relation] * inverse % PRIME
                            for relation in range(P))
                      for difference in range(P))
                for _r1 in range(P)
            )
            for r0 in range(P)
        )
        for point in range(len(POINTS))
    )


def point_difference_flat_lift(two_digit):
    answer = [[[[[0 for _relation in range(P)] for _difference in range(P)]
                for _r1 in range(P)] for _r0 in range(P)]
              for _point in POINTS]
    for state, fibre in enumerate(STATE_FIBRES):
        inverse = pow(len(fibre) * P, -1, PRIME)
        for point in fibre:
            for r0 in range(P):
                for r1 in range(P):
                    for difference in range(P):
                        for relation in range(P):
                            answer[point][r0][r1][difference][relation] = (
                                two_digit[state][r0][r1][relation] * inverse
                            ) % PRIME
    return freeze_tensor(answer)


def rank_mod(matrix) -> int:
    basis = {}
    for source in matrix:
        row = [value % PRIME for value in source]
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
        basis[pivot] = tuple(value * inverse % PRIME for value in row)
    return len(basis)


def flatten_st(block):
    return tuple(value for line in block for value in line)


def base_rows(base, points=None):
    chosen = tuple(range(len(POINTS))) if points is None else tuple(points)
    return tuple(flatten_st(base[point]) for point in chosen)


def parent_rows(parent, points=None, fixed_r0=None):
    chosen = tuple(range(len(POINTS))) if points is None else tuple(points)
    r0_values = range(P) if fixed_r0 is None else (fixed_r0,)
    return tuple(
        flatten_st(parent[point][r0])
        for point in chosen for r0 in r0_values
    )


def child_rows(tensor, points=None, fixed_r0=None, fixed_r1=None):
    chosen = tuple(range(len(POINTS))) if points is None else tuple(points)
    r0_values = range(P) if fixed_r0 is None else (fixed_r0,)
    r1_values = range(P) if fixed_r1 is None else (fixed_r1,)
    return tuple(
        flatten_st(tensor[point][r0][r1])
        for point in chosen for r0 in r0_values for r1 in r1_values
    )


def r1_amplitude_ranks(tensor):
    raw = tuple(
        tuple(
            tensor[point][r0][r1][difference][relation]
            for point in range(len(POINTS)) for r0 in range(P)
            for difference in range(P) for relation in range(P)
        )
        for r1 in range(P)
    )
    inverse = pow(P, -1, PRIME)
    means = tuple(sum(raw[r1][column] for r1 in range(P)) * inverse % PRIME
                  for column in range(len(raw[0])))
    contrast = tuple(
        tuple((value - means[column]) % PRIME
              for column, value in enumerate(raw[r1]))
        for r1 in range(P)
    )
    return rank_mod(raw), rank_mod(contrast)


def response_rank_record(tensor):
    one = tensor_one_digit_pointed(tensor)
    base = tensor_pointed_parent(one)
    base_matrix = base_rows(base)
    one_matrix = parent_rows(one)
    two_matrix = child_rows(tensor)
    global_record = (
        rank_mod(base_matrix), rank_mod(one_matrix),
        rank_mod(base_matrix + one_matrix), rank_mod(two_matrix),
        rank_mod(one_matrix + two_matrix), r1_amplitude_ranks(tensor),
    )
    state_records = tuple(
        (
            len(fibre),
            rank_mod(base_rows(base, fibre)),
            rank_mod(parent_rows(one, fibre)),
            rank_mod(child_rows(tensor, fibre)),
            rank_mod(parent_rows(one, fibre) + child_rows(tensor, fibre)),
        )
        for fibre in STATE_FIBRES
    )
    fixed_r0 = tuple(
        (
            rank_mod(parent_rows(one, fixed_r0=r0)),
            rank_mod(child_rows(tensor, fixed_r0=r0)),
            rank_mod(parent_rows(one, fixed_r0=r0)
                     + child_rows(tensor, fixed_r0=r0)),
        )
        for r0 in range(P)
    )
    return global_record, state_records, fixed_r0


def matrix_inverse(matrix):
    size = len(matrix)
    rows = [
        [value % PRIME for value in matrix[index]]
        + [int(index == column) for column in range(size)]
        for index in range(size)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if rows[row][column]), None)
        require(pivot is not None, ("singular matrix", matrix))
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
    return tuple(tuple(row[size:]) for row in rows)


def matrix_multiply(left, right):
    return tuple(
        tuple(sum(left[row][middle] * right[middle][column]
                  for middle in range(len(right))) % PRIME
              for column in range(len(right[0])))
        for row in range(len(left))
    )


def pivot_columns(rows):
    width = len(rows[0])
    chosen = []
    rank = 0
    for column in range(width):
        trial = tuple(tuple(row[index] for index in chosen + [column]) for row in rows)
        trial_rank = rank_mod(trial)
        if trial_rank > rank:
            chosen.append(column)
            rank = trial_rank
            if rank == len(rows):
                break
    return tuple(chosen)


def coordinate_matrix(parent, child):
    require(len(parent) == len(child) == len(POINTS), "coordinate row count")
    require(rank_mod(parent) == len(POINTS), "parent coordinate rank")
    pivots = pivot_columns(parent)
    require(len(pivots) == len(POINTS), ("coordinate pivots", pivots))
    parent_minor = tuple(tuple(parent[row][column] for column in pivots)
                         for row in range(len(POINTS)))
    child_minor = tuple(tuple(child[row][column] for column in pivots)
                        for row in range(len(POINTS)))
    matrix = matrix_multiply(child_minor, matrix_inverse(parent_minor))
    predicted = matrix_multiply(matrix, parent)
    witness = next((
        (row, column, predicted[row][column], child[row][column])
        for row in range(len(POINTS)) for column in range(len(parent[0]))
        if predicted[row][column] != child[row][column]
    ), None)
    return matrix, witness


def stationary_record(tensor):
    one = tensor_one_digit_pointed(tensor)
    parent = tuple(
        tuple(one[point][r0][difference][relation]
              for r0 in range(P) for difference in range(P) for relation in range(P))
        for point in range(len(POINTS))
    )
    parent_rank = rank_mod(parent)
    records = []
    for r1 in range(P):
        child = tuple(
            tuple(tensor[point][r0][r1][difference][relation]
                  for r0 in range(P) for difference in range(P)
                  for relation in range(P))
            for point in range(len(POINTS))
        )
        child_rank = rank_mod(child)
        union_rank = rank_mod(parent + child)
        if parent_rank == len(POINTS):
            matrix, witness = coordinate_matrix(parent, child)
        else:
            matrix, witness = None, ("parent-rank", parent_rank)
        exists = union_rank == parent_rank
        records.append((child_rank, union_rank, exists, witness))
    return parent_rank, tuple(records)


def matrix_rank_support(matrix):
    return rank_mod(matrix), sum(value != 0 for row in matrix for value in row)


def address_bundle_record(tensor, source_kernel=None):
    one = tensor_one_digit_pointed(tensor)
    matrices = []
    contained = unique = diagonal = state_block = scalar = 0
    rank_hist = {}
    support_hist = {}
    first_failure = None
    source_comparable = source_matches = 0
    first_source_mismatch = None
    fixed_r0_sums = []
    for r0 in range(P):
        parent = tuple(flatten_st(one[point][r0]) for point in range(len(POINTS)))
        parent_rank = rank_mod(parent)
        r0_matrices = []
        for r1 in range(P):
            child = tuple(flatten_st(tensor[point][r0][r1])
                          for point in range(len(POINTS)))
            union_rank = rank_mod(parent + child)
            if union_rank == parent_rank:
                contained += 1
            if parent_rank != len(POINTS):
                if first_failure is None:
                    first_failure = (r0, r1, parent_rank, union_rank)
                matrices.append(None)
                r0_matrices.append(None)
                continue
            matrix, witness = coordinate_matrix(parent, child)
            require(witness is None, ("address-conditioned residual", r0, r1, witness))
            unique += 1
            matrices.append(matrix)
            r0_matrices.append(matrix)
            is_state_block = all(
                matrix[left][right] == 0
                for left, (left_state, _left_root) in enumerate(POINTS)
                for right, (right_state, _right_root) in enumerate(POINTS)
                if left_state != right_state
            )
            is_diagonal = all(
                matrix[left][right] == 0
                for left in range(len(POINTS)) for right in range(len(POINTS))
                if left != right
            )
            state_block += int(is_state_block)
            diagonal += int(is_diagonal)
            scalar += int(is_diagonal and len({matrix[i][i] for i in range(len(POINTS))}) == 1)
            rank, support = matrix_rank_support(matrix)
            rank_hist[rank] = rank_hist.get(rank, 0) + 1
            support_hist[support] = support_hist.get(support, 0) + 1
            if source_kernel is not None and is_diagonal:
                for point in range(len(POINTS)):
                    source_value = source_kernel[point][r0][r1]
                    if source_value is None:
                        continue
                    source_comparable += 1
                    if source_value == matrix[point][point]:
                        source_matches += 1
                    elif first_source_mismatch is None:
                        first_source_mismatch = (
                            point, r0, r1, source_value, matrix[point][point]
                        )
        if all(matrix is not None for matrix in r0_matrices):
            summed = tuple(tuple(
                sum(matrix[row][column] for matrix in r0_matrices) % PRIME
                for column in range(len(POINTS))
            ) for row in range(len(POINTS)))
            identity = tuple(tuple(int(row == column) for column in range(len(POINTS)))
                             for row in range(len(POINTS)))
            require(summed == identity, ("fixed-r0 projective identity", r0))
            fixed_r0_sums.append(summed)
    actual_matrices = tuple(matrix for matrix in matrices if matrix is not None)
    distinct = len(set(actual_matrices))
    unique_by_r1 = tuple(len({matrices[r0 * P + r1] for r0 in range(P)
                              if matrices[r0 * P + r1] is not None})
                         for r1 in range(P))
    return (
        contained, unique, state_block, diagonal, scalar,
        len(fixed_r0_sums), tuple(sorted(rank_hist.items())),
        tuple(sorted(support_hist.items())), distinct, unique_by_r1,
        first_failure, source_comparable, source_matches,
        first_source_mismatch, digest_json(actual_matrices),
    ), tuple(matrices)


def direct_probe_record(chunks, tensor, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    records = []
    by_alpha = {chunk[0]: chunk[7] for chunk in chunks}
    for probe_index, probe in enumerate(INVERSION_PROBES):
        values = []
        for relation in PROBE_RELATIONS:
            total = 0
            for alpha in range(P):
                sequence = by_alpha[alpha][probe_index]
                require(len(sequence) == P * P, ("probe sequence", alpha, probe))
                index = 0
                for _beta in range(P):
                    for tau in range(P):
                        total += (
                            sequence[index]
                            * pow(zeta, -(alpha + tau * relation) % P, PRIME)
                        )
                        index += 1
            value = total % PRIME * normalizer % PRIME
            point, r0, r1, difference = probe
            require(value == tensor[point][r0][r1][difference][relation],
                    ("direct inversion probe", probe, relation, value,
                     tensor[point][r0][r1][difference][relation]))
            values.append(value)
        records.append((probe, tuple(values)))
    return tuple(records)


def main() -> None:
    data = source_data()
    cells = data["cells"]
    source_support = source_support_record(cells)
    source_kernels, source_ratios = source_ratio_record(cells)

    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    tau_actual = [zero_row() for _tau in range(P)]
    tau_support = [zero_row() for _tau in range(P)]
    chunks_by_alpha = {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(worker, alpha): alpha for alpha in range(P)}
        for future in as_completed(futures):
            chunk = future.result()
            alpha = chunk[0]
            chunks_by_alpha[alpha] = chunk
            alpha_phase = pow(zeta, -alpha % P, PRIME)
            for tau in range(P):
                add_row_to_core(tau_actual, tau, chunk[1][tau], alpha_phase)
                add_row_to_core(tau_support, tau, chunk[2][tau], alpha_phase)
    chunks = tuple(chunks_by_alpha[alpha] for alpha in range(P))

    one_json = combine_json_blocks(tuple(chunk[3] for chunk in chunks))
    pointed_json = combine_json_blocks(tuple(chunk[4] for chunk in chunks))
    one_gamma_digest = sha256(one_json.encode("ascii")).hexdigest()
    pointed_gamma_digest = sha256(pointed_json.encode("ascii")).hexdigest()
    require(one_gamma_digest == ONE_DIGIT_ROOT_GAMMA_SHA256,
            ("one-digit/root gamma marginal", one_gamma_digest))
    require(pointed_gamma_digest == POINTED_GAMMA_SHA256,
            ("pointed gamma marginal", pointed_gamma_digest))

    controls = dict(item for chunk in chunks for item in chunk[6])
    literal_records = []
    for alpha, beta, tau in CONTROLS:
        maps, counts = pair_coefficient_maps(alpha, beta, tau)
        raw_actual, raw_support = expand_pair_map(maps[0])
        phase = pow(zeta, beta, PRIME)
        direct = (scale_row(raw_actual, phase), scale_row(raw_support, phase))
        require(direct == controls[(alpha, beta, tau)],
                ("literal guard control", alpha, beta, tau))
        literal_records.append(((alpha, beta, tau), counts))

    actual = inverse_core(tau_actual, zeta)
    support = inverse_core(tau_support, zeta)
    two_digit = tensor_two_digit_marginal(actual)
    require(digest_json(two_digit) == TD.EXPECTED_DIGESTS[0],
            ("two-digit tensor marginal", digest_json(two_digit)))
    one_pointed = tensor_one_digit_pointed(actual)
    one_state = tensor_one_digit_state(one_pointed)
    pointed = tensor_pointed_parent(one_pointed)
    require(digest_json(one_state) == ONE_DIGIT_ROOT_TENSOR_SHA256,
            ("one-digit/root tensor marginal", digest_json(one_state)))
    require(digest_json(pointed) == POINTED_TENSOR_SHA256,
            ("pointed tensor marginal", digest_json(pointed)))
    support_one = tensor_one_digit_pointed(support)
    require(support_one == one_pointed,
            "support normalization misses complete one-digit pointed parent")
    require(all(actual[point][r0][r1][0][relation] == 0
                for point in range(len(POINTS)) for r0 in range(P)
                for r1 in range(P) for relation in range(P)),
            "same-root tensor sector")

    flat = r1_flat_lift(one_pointed)
    point_flat = point_difference_flat_lift(two_digit)
    require(tensor_one_digit_pointed(flat) == one_pointed, "r1-flat parent")
    require(tensor_two_digit_marginal(point_flat) == two_digit,
            "point/difference-flat two-digit parent")
    tensors = (actual, support, flat, point_flat)
    names = ("actual", "r1_support_normalized", "r1_flat", "point_difference_flat")
    tensor_digests = tuple(digest_json(tensor) for tensor in tensors)
    require(tensor_digests == EXPECTED_CANDIDATE_TENSOR_SHA256,
            ("submitted tensor comparison", tensor_digests,
             EXPECTED_CANDIDATE_TENSOR_SHA256))

    rank_records = tuple(response_rank_record(tensor) for tensor in tensors)
    require(rank_records[0][0] == (6, 6, 6, 6, 6, (13, 12)),
            ("actual amplitude/carrier rank", rank_records[0][0]))
    require(rank_records[1][0] == (6, 6, 6, 6, 6, (4, 3)),
            ("support amplitude/carrier rank", rank_records[1][0]))
    require(rank_records[2][0] == (6, 6, 6, 6, 6, (1, 0)),
            ("flat amplitude/carrier rank", rank_records[2][0]))
    require(rank_records[0][1] == (
        (1, 1, 1, 1, 1), (2, 2, 2, 2, 2),
        (1, 1, 1, 1, 1), (2, 2, 2, 2, 2),
    ), ("actual statewise carrier", rank_records[0][1]))
    require(all(record == (6, 6, 6) for record in rank_records[0][2]),
            ("fixed-r0 carrier", rank_records[0][2]))

    stationary = tuple(stationary_record(tensor) for tensor in tensors)
    actual_stationary_shapes = tuple(
        (child_rank, union_rank, exists)
        for child_rank, union_rank, exists, _witness in stationary[0][1]
    )
    require(actual_stationary_shapes == (
        (4, 10, False),
        *((6, 12, False),) * 5,
        (4, 10, False),
        *((6, 12, False),) * 5,
        (4, 10, False),
    ), ("stationary K_r1 obstruction", actual_stationary_shapes))
    require(stationary[0][1][0][3]
            == (2, 2041, 0, 148404725567356333796897),
            ("first stationary witness", stationary[0][1][0][3]))
    require(all(record[2] for record in stationary[2][1]),
            "flat stationary positive control")

    address_records = []
    matrix_banks = []
    for index, tensor in enumerate(tensors):
        kernel = source_kernels[index] if index < 2 else None
        record, matrices = address_bundle_record(tensor, kernel)
        address_records.append(record)
        matrix_banks.append(matrices)
    address_records = tuple(address_records)
    require(address_records[0][:10] == (
        169, 169, 169, 169, 35, 13,
        ((0, 35), (2, 2), (4, 2), (6, 130)),
        ((0, 35), (2, 2), (4, 2), (6, 130)),
        121,
        (2, 13, 13, 13, 13, 13, 3, 13, 13, 13, 13, 13, 2),
    ), ("actual diagonal bundle", address_records[0]))
    require(address_records[0][11:14] == (954, 954, None),
            ("actual source/output ratios", address_records[0][11:14]))
    require(address_records[1][11:14] == (1014, 1014, None),
            ("support source/output ratios", address_records[1][11:14]))
    inverse_13 = pow(P, -1, PRIME)
    flat_expected = tuple(tuple(inverse_13 if row == column else 0
                                for column in range(len(POINTS)))
                          for row in range(len(POINTS)))
    require(all(matrix == flat_expected for matrix in matrix_banks[2]),
            "r1-flat diagonal identity splitter")

    require(source_ratios[0][0] == ((1, 954), (2, 40), (3, 20)),
            ("actual source ratio histogram", source_ratios[0]))
    require(source_ratios[0][1] == 60,
            ("exceptional source fibre count", source_ratios[0]))
    expected_exception_fibres = (
        (0, (0, 0), 12, (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)),
        (1, (1, 0), 12, (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)),
        (2, (1, 6), 9, (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)),
        (3, (3, 6), 3, (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)),
        (4, (3, 12), 0, (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)),
        (5, (2, 12), 0, (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)),
    )
    require(source_ratios[0][2] == expected_exception_fibres,
            ("exceptional source fibres", source_ratios[0][2]))
    require(source_ratios[1][0] == ((1, 1014),) and source_ratios[1][1] == 0,
            ("support source proportionality", source_ratios[1]))

    probe_record = direct_probe_record(chunks, actual, zeta)
    work_records = tuple(chunk[8] for chunk in chunks)
    totals = tuple(sum(record[0][index] for record in work_records)
                   for index in range(5))
    state_totals = tuple(sum(record[1][state] for record in work_records)
                         for state in range(V))
    require(totals == TD.EXPECTED_WORK_COUNTS,
            ("endpoint work totals", totals, TD.EXPECTED_WORK_COUNTS))
    require(state_totals == TD.EXPECTED_STATE_COUNTS,
            ("endpoint state totals", state_totals, TD.EXPECTED_STATE_COUNTS))

    joint_merkle = digest_json(tuple(
        pair for chunk in chunks for pair in chunk[5]
    ))
    semantic_surface = (
        TWO_DIGIT_AUDIT_SHA256, TD.EXPECTED_SEMANTIC_SHA256,
        POINTED_AUDIT_SHA256, PA.EXPECTED_SEMANTIC_SHA256,
        CANDIDATE_SHA256_OBSERVED, data["source_sha"], source_support,
        source_ratios, tuple(literal_records), one_gamma_digest,
        pointed_gamma_digest, joint_merkle, tensor_digests, rank_records,
        stationary, address_records, probe_record, totals, state_totals,
        digest_json(work_records),
    )
    semantic = digest_json(semantic_surface)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("audit semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== independent hostile audit: two current digits x pointed root-difference diagonal bundle ==")
    print(f"parents=(two_digit_independent_sha={TWO_DIGIT_AUDIT_SHA256},pointed_independent_sha={POINTED_AUDIT_SHA256})")
    print(f"candidate_source=(observed_sha={CANDIDATE_SHA256_OBSERVED},published_sha={CANDIDATE_SHA256_PUBLISHED},published_metadata_matches=False)")
    print("construction=ordered_(u,q)_pair coefficients aggregated by exact source cell before independent 169-address expansion;candidate module not imported")
    print(f"coordinates=(point={POINTS},r0,r1,s=u-q,relation);field=(prime={PRIME},zeta13={zeta})")
    print(f"source=(sha256={data['source_sha']},r1_support={source_support})")
    print(f"source_ratio_records_(actual,support)={source_ratios}")
    print(f"work=(totals={totals},state_totals={state_totals},ordered_pair_updates={sum(record[2] for record in work_records)},keys_by_tau={tuple(sum(record[3][tau] for record in work_records) for tau in range(P))})")
    print(f"literal_controls={CONTROLS}: PASS;direct_inversion_probes={probe_record}: PASS;same_root_s0=0")
    print(f"preintegration_marginals=(two_digit_alpha_blocks={TD.EXPECTED_ALPHA_GAMMA_SHA256},one_digit_root={one_gamma_digest},pointed_root={pointed_gamma_digest}): PASS")
    print("inverse_marginals=(two_digit,one_digit_root,pointed_root): PASS")
    print(f"tensor_digests_(actual,support,r1_flat,point_difference_flat)={tensor_digests}: MATCH submitted candidate")
    print(f"response_rank_records={tuple(zip(names, rank_records))}")
    print(f"stationary_K_r1_records={tuple(zip(names, stationary))}")
    print(f"address_diagonal_bundle_records={tuple(zip(names, address_records))}")
    print("fixed_r0_projective_sums=sum_r1 K_(r0,r1)=I6 for all 13 r0: PASS")
    print("hostiles=(support amplitude 4/3 with same one-digit-pointed parent;r1-flat 1/0 with K=13^-1 I6;point/difference-flat same two-digit parent but carrier base rank4): PASS")
    print("exceptional_source_profiles=60 in six listed (point,r0) fibres;all integrated residuals annihilated by the common endpoint operator: PASS")
    print(f"joint_gamma_merkle={joint_merkle};semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("verdict=ACCEPT finite-exact five-coordinate diagonal-bundle candidate;repair stale published script hash")
    print("typing=static address-conditioned partition of identity only;not clock,not complete address,not arrival ancestry,not physical current,no row exclusion,no LRC(14)")
    print("commands=python -B 04-computation/lrc_r5_two_current_digits_pointed_root_difference_carrier_transition_independent_audit_20260816.py;python -B -O 04-computation/lrc_r5_two_current_digits_pointed_root_difference_carrier_transition_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
