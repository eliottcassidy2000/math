#!/usr/bin/env python3
"""Independent hostile audit of two current inverse digits over the r=5 square.

The submitted two-digit probe is never imported.  This audit starts from the
independently audited one-digit current-leg tensor, reconstructs the source
packet from its canonical THM-2471 dependencies, and retains

    a = r0 + 13*r1 + 169*c,       0 <= r0,r1 < 13.

The source profile is folded by 13^3, pulled first through the r1 window,
then through r0, and finally through collision root u.  A separately written
generic modulus-window routine proves literal equality with the direct
mod-169 half-open window.  Endpoint coefficients are integrated before the
169-address expansion; this is an exact Fubini regrouping on the 33 source
cells, not an import of the candidate implementation.

The conditional-rank hostile has a narrow conclusion.  It refutes scalar
state/relation-independent factorizations K(r0,r1)T1(state,r0,t).  It does
not refute hidden-state, matrix-valued, state-dependent, nonlinear, or
temporal recurrences.  No complete address, chronology, physical current,
row exclusion, or LRC(14) conclusion is asserted.
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
    ROOT
    / "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_inverse_owner_branch_"
    "independent_audit_20260816.py"
)
PARENT_SHA256 = "4b0bd05ffa6195ff484433329e334d771bc27e7cd380136b50b45e7248bb98ba"
PARENT_SEMANTIC = "7063720d0e0e4847ce752102de83274ea47d7740fc435a64bae425dbd7100121"

P = 13
V = 4
CONTROL_TRIPLES = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
INVERSION_PROBES = ((0, 0, 1), (1, 6, 7), (2, 4, 9), (3, 12, 11))
PROBE_RELATIONS = (0, 6, 12)
EXPECTED_SOURCE_SHA256 = "dd2f48375fc38b419babdd0ed13365fb56918829a663270c4b23d56635d6e097"
EXPECTED_WORK_COUNTS = (7107008, 7108460, 1200462, 186244, 186244)
EXPECTED_STATE_COUNTS = (300115, 300116, 300115, 300116)
EXPECTED_COEFFICIENT_UPDATES = 1248644
EXPECTED_KEY_COUNTS = (310, 310, 1186, 1186, 1186, 902, 618, 618, 902, 1186, 1186, 1186, 310)
EXPECTED_ALPHA_GAMMA_SHA256 = (
    "9eb5ee6775d877e319272190d1d39df9a1bea1cb60f7947cdc4f34c9f7f362e1",
    "a4e0cc5521ec8dae63f513ec0e34ac96d550d2597676689ea0f2bc3ea036d9d5",
    "6d796c3fc5d0a1f58ab701f0fccc183d3c7954a6b61f27120a7735ddb3e74dc1",
    "6c523a119e6fef952b29b5c4cb7f77bede76da74e63fc759d81a06c41434c740",
    "93f13337db351986d01d772b02868794135021472f1c7cd03d061403be642c12",
    "216c33d878bfad06dddd75e06a4589773d7a961bde423f1504940a3b1e4108f7",
    "bc85bce549d7ab8e0abaec1ee488dd77efdeaa1126320b9d00fc54f6ff697592",
    "ffefd0093190c730018653856e43214f9869fc0d34e4dd5edb71fc1885cbf3da",
    "b865880bbbb274ffc43b82ef6fd3ef7cbc79826a8f4ed7e43969df227a6b74a5",
    "762b0a14df8ed52dd942b3e1a927864c7ea85c0cf27f14770c856baa02b449ca",
    "5b0f2f8c3957dd8d3a5fce6d45de993376603d7d713e9c5ef9303a0b2cd83f08",
    "fa5098b6c025587ce3741b53bcf9f2196931e1ce80274ad8b47a73bf0f91b694",
    "32c3d17ca29cadcb6f0c3fc4dde9d5690cbf1d903c514df1e1d24dab0ce3ca67",
)
EXPECTED_DIGESTS = (
    "53eb6e618d0669bdb27841a1800c46e32456bb6c6d3698c590ae0d5e68822033",
    "fd1b837d9de3e4f9e586d29b69ed6726364ce97535d2e48f441c6ccd694250de",
    "2a195fac7fbd60a4bbd2597bf34baf87302ead16c7c0e8c67c0667b0d320dfba",
    "5f4b9609faaa5f148d112a7cde5cfba0ab2c1385b4c53ea9c4bcfc6e93d106fc",
    "0ba2b6072c22ed3ef8058f54c428923b633892849b4cc3ddad0d137479cceca8",
    (
        "50ac8fbd386ee78aead3ad370e6945c521a44cd426934f611531edf4c54a4dba",
        "9f8cd9688152ba01c408eccbf5979960c60a693c52e226096bef896e8ef7338d",
    ),
    (
        "4351f36a0727317990ea67b969f9b68dd668a9f5e904f24b0cdf44a95dcfcca6",
        "eb4ad01f374e55b0bca39a47290aa5a9bf96953f1c5aff4693995d0a3b1421e5",
    ),
    (
        "e7bd15c1765c8cfe38d6ac8f53d9e55d7fc30a788c812e1fb09fbddacd4ddffa",
        "9f8cd9688152ba01c408eccbf5979960c60a693c52e226096bef896e8ef7338d",
    ),
)
EXPECTED_RANKS = (
    ((4, 13, 13, 6, 4), (4, 4, 1, 6, 4)),
    ((3, 12, 12, 6, 4), (0, 0, 0, 0, 0)),
)
EXPECTED_CONDITIONAL_RANKS = (
    (
        (4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4),
        (4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4),
    ),
    (
        (1,) * P,
        (0,) * P,
    ),
)
EXPECTED_CENSUSES = (
    (
        (8788, 1, 12, 12, 144, 12, 144, 144, 1728,
         3, 36, 36, 432, 36, 432, 432, 5184),
        (676, 1, 12, 0, 0, 12, 144, 0, 0,
         3, 36, 0, 0, 36, 432, 0, 0),
    ),
    (
        (5184, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 5184),
        (0,) * 17,
    ),
)
EXPECTED_FIXED = (
    (
        4,
        (132, 133, 132, 133),
        (4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4),
        "d4c0e037f48c7f7c5759d497798c208e3e70fd94600d827d80d22c568043522c",
    ),
    (
        4,
        (169, 169, 169, 169),
        (1,) * P,
        "9b9e8372ae008139744bac343d9a89b7003691d46ddb4c9d5ced071e594b8a9f",
    ),
)
EXPECTED_PROBES = (
    ((0, 0, 1), (322383479003395571142375, 652209513571172027017012, 380488686199075270671641)),
    ((1, 6, 7), (697962230161895224492129, 321725971461886166285600, 242345101589067781680478)),
    ((2, 4, 9), (9243580176343212298593, 85755045745196838326775, 624201457125133858214241)),
    ((3, 12, 11), (546928898629614439108454, 48211157869431777318692, 258819341437113154438961)),
)
EXPECTED_SEMANTIC_SHA256 = "f2af726e9d5abd1487e841623ce2f62ca647c86e5a1a68e41eda4d9dda6c81ac"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_parent():
    require(lf_sha256(PARENT_PATH) == PARENT_SHA256, "one-digit parent hash drift")
    spec = importlib.util.spec_from_file_location("two_digit_independent_parent", PARENT_PATH)
    require(spec is not None and spec.loader is not None, "one-digit parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.EXPECTED_SEMANTIC_SHA256 == PARENT_SEMANTIC,
            "one-digit parent semantic drift")
    return module


ONE = load_parent()
SQ = ONE.SQ
R = ONE.R
PRIME = R.JOINT_PRIME
WALSH_SIGNS = ONE.WALSH_SIGNS


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


def pull_window(profile, index: int, modulus: int, grid: int):
    """Pull one exact half-open inverse window and map it back to [0,grid)."""
    require(grid % modulus == 0, ("window modulus", modulus, grid))
    require(0 <= index < modulus, ("window index", index, modulus))
    starts, values = profile
    width = grid // modulus
    lower = index * width
    upper = lower + width
    cursor = lower
    source_index = bisect_right(starts, lower) - 1
    out_starts = []
    out_values = []
    while cursor < upper:
        stop = starts[source_index + 1] if source_index + 1 < len(starts) else grid
        stop = min(stop, upper)
        mapped = modulus * cursor - index * grid
        value = values[source_index]
        if not out_starts or out_values[-1] != value:
            out_starts.append(mapped)
            out_values.append(value)
        cursor = stop
        source_index += 1
    require(out_starts and out_starts[0] == 0,
            ("window starts", index, modulus, out_starts[:1]))
    return tuple(out_starts), tuple(out_values)


def freeze3(bank):
    return tuple(tuple(tuple(row) for row in plane) for plane in bank)


def freeze4(bank):
    return tuple(
        tuple(tuple(tuple(row) for row in plane) for plane in slab)
        for slab in bank
    )


@lru_cache(maxsize=1)
def source_context():
    grid = R.SRC.T_DEN
    joint = R.JOINT_COORDINATE
    e_intervals = R.SRC.build_set(R.SRC.PAT_E, R.SRC.ZELL)
    q_intervals = R.SRC.build_set(R.SRC.PAT_QB, R.SRC.ZELL)
    packet = R.fold_weighted(
        [(left, right, 1) for left, right in e_intervals], R.SRC.RPKT, grid
    )
    f_pieces = R.intersect_profile_with_set(packet, q_intervals, grid)

    base_fold = R.fold_weighted(f_pieces, R.SRC.DCOLL // P**2, grid)
    r1_windows = tuple(pull_window(base_fold, r1, P, grid) for r1 in range(P))
    address_windows = tuple(
        tuple(pull_window(r1_windows[r1], r0, P, grid) for r1 in range(P))
        for r0 in range(P)
    )
    for r0 in range(P):
        for r1 in range(P):
            direct = pull_window(base_fold, r0 + P * r1, P**2, grid)
            require(address_windows[r0][r1] == direct,
                    ("nested/direct mod169 window", r0, r1))

    raw = tuple(
        tuple(
            tuple(pull_window(address_windows[r0][r1], root, P, grid)
                  for r1 in range(P))
            for r0 in range(P)
        )
        for root in range(P)
    )
    scale = joint // grid
    profiles = tuple(
        tuple(
            tuple((tuple(point * scale for point in profile[0]), profile[1])
                  for profile in r1_rows)
            for r1_rows in r0_rows
        )
        for r0_rows in raw
    )

    (
        one_digit, _one_boundaries, source_u, source_v,
        source_boundaries, source_profile_digest, _one_record, _one_sha,
        _one_reflection,
    ) = ONE.owner_digit_profiles()
    boundaries = tuple(sorted(
        {0, joint}
        | {point for profile in source_u + source_v for point in profile[0]}
        | {point for root_rows in profiles for r0_rows in root_rows
           for profile in r0_rows for point in profile[0]}
    ))
    require(tuple(joint - point for point in reversed(boundaries)) == boundaries,
            "two-digit boundary reflection")

    cells = []
    state_address_support = [set() for _state in range(V)]
    reflection_intervals = 0
    square_intervals = 0
    for left, right in zip(boundaries, boundaries[1:]):
        address_values = tuple(
            tuple(tuple(profile_value(profiles[root][r0][r1], left)
                        for r1 in range(P)) for r0 in range(P))
            for root in range(P)
        )
        u_values = tuple(profile_value(profile, left) for profile in source_u)
        v_values = tuple(profile_value(profile, left) for profile in source_v)
        for root in range(P):
            for r0 in range(P):
                observed = sum(address_values[root][r0])
                expected = profile_value(one_digit[root][r0], left)
                require(observed == expected,
                        ("r1 profile marginal", left, root, r0, observed, expected))
            require(sum(address_values[root][r0][r1]
                        for r0 in range(P) for r1 in range(P)) == u_values[root],
                    ("two-digit profile total", left, root))
            require(all(address_values[root][r0][r1] * v_values[root] == 0
                        for r0 in range(P) for r1 in range(P)),
                    ("two-digit pointwise same-root", left, root))
            for r0 in range(P):
                for r1 in range(P):
                    reflected = profile_value(
                        profiles[P - 1 - root][P - 1 - r0][P - 1 - r1],
                        joint - right,
                    )
                    require(address_values[root][r0][r1] == reflected,
                            ("two-digit reflection", left, root, r0, r1))
        reflection_intervals += 1

        support = frozenset(root for root, value in enumerate(u_values) if value)
        state = SQ.STATE_INDEX.get(support)
        reflected_left = joint - right
        reflected_support = frozenset(
            root for root, profile in enumerate(source_u)
            if profile_value(profile, reflected_left)
        )
        require(reflected_support == frozenset(P - 1 - root for root in support),
                ("root-set reflection", left, support, reflected_support))
        if state is not None:
            require(SQ.STATE_INDEX[reflected_support] == state ^ 2,
                    ("state XOR2", left, state, reflected_support))
            square_intervals += 1
            for root in range(P):
                count = sum(address_values[root][r0][r1] != 0
                            for r0 in range(P) for r1 in range(P))
                if count:
                    state_address_support[state].add(count)
        support_mask = sum(1 << root for root, value in enumerate(u_values) if value)
        cells.append((state, u_values, v_values, support_mask, address_values))

    require(all(values == {132} for values in state_address_support),
            ("active-root address support", state_address_support))
    address_masses = tuple(
        sum(profile_mass(raw[root][r0][r1], grid) for root in range(P))
        for r0 in range(P) for r1 in range(P)
    )
    require(address_masses == tuple(
        address_masses[(P - 1 - r0) * P + (P - 1 - r1)]
        for r0 in range(P) for r1 in range(P)
    ), "address-mass reflection")
    source_record = (
        R.SRC.DCOLL,
        R.SRC.DCOLL // P**2,
        len(f_pieces),
        len(base_fold[0]),
        tuple(len(profile[0]) for profile in r1_windows),
        sum(len(profile[0]) for row in address_windows for profile in row),
        sum(len(profile[0]) for root_rows in raw
            for row in root_rows for profile in row),
        len(boundaries),
        tuple(tuple(sorted(values)) for values in state_address_support),
        address_masses,
        digest_json(raw),
        digest_json(profiles),
    )
    source_sha = digest_json((source_record, boundaries))
    require(source_sha == EXPECTED_SOURCE_SHA256,
            ("source digest", source_sha, EXPECTED_SOURCE_SHA256))
    reflection_record = (reflection_intervals, square_intervals, address_masses)
    return (
        profiles, boundaries, tuple(cells), one_digit, source_u, source_v,
        source_boundaries, source_profile_digest, source_record, source_sha,
        reflection_record,
    )


@lru_cache(maxsize=1)
def endpoint_context():
    return ONE.endpoint_context()


@lru_cache(maxsize=200000)
def sparse_left_vector(cell_index: int, selected_u_mask: int):
    _profiles, _boundaries, cells, *_rest = source_context()
    state, _u_values, _v_values, _support_mask, address_values = cells[cell_index]
    require(state is not None, ("untyped left vector", cell_index))
    roots = tuple(root for root in range(P) if (selected_u_mask >> root) & 1)
    entries = []
    for r0 in range(P):
        for r1 in range(P):
            value = sum(address_values[root][r0][r1] for root in roots)
            if value:
                entries.append((state, r0, r1, value))
    return tuple(entries)


def coefficient_maps(alpha: int, beta: int, literal_tau: int | None = None):
    _profiles, boundaries, cells, *_rest = source_context()
    word, endpoint_grid, harmonic, danger = endpoint_context()
    events, interval_count, mapped_count = R.endpoint_events(
        word, endpoint_grid, alpha, beta, literal_tau
    )
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    maps = [dict() for _tau in tau_values]
    positions = sorted(set(events) | set(boundaries))
    mask = 0
    active = q_active = weighted = 0
    state_counts = [0] * V
    update_count = 0
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
                ("two-digit escaped cell zero", alpha, beta, left, right))
        chamber = R.chamber(left, right)
        require(chamber in ("left", "right"),
                ("two-digit entered middle chamber", alpha, beta, left, right))
        cell_index = bisect_right(boundaries, left) - 1
        require(0 <= cell_index < len(cells)
                and right <= boundaries[cell_index + 1],
                ("source-cell crossing", left, right, cell_index))
        state, u_values, v_values, support_mask, _address_values = cells[cell_index]
        require(state is not None, ("active untyped source state", left, right))
        state_counts[state] += 1
        if not jump:
            continue
        q_active += 1
        require(any(u_values) and any(v_values), ("empty weighted source", left))
        weighted += 1
        for row_index, tau in enumerate(tau_values):
            selected_mask = mask if literal_tau is not None else (
                mask & R.guard_mask(chamber, tau, danger)
            )
            selected_u_mask = selected_mask & support_mask
            right_value = sum(v_values[root] for root in range(P)
                              if (selected_mask >> root) & 1)
            if not selected_u_mask or not right_value:
                continue
            require(all(u_values[root] * v_values[root] == 0
                        for root in range(P) if (selected_mask >> root) & 1),
                    ("selected same-root", alpha, beta, tau, left))
            key = (cell_index, selected_u_mask)
            maps[row_index][key] = (
                maps[row_index].get(key, 0) + right_value * jump
            ) % PRIME
            update_count += 1
    mask ^= events.get(positions[-1], 0)
    require(mask == 0, ("endpoint mask closure", alpha, beta, literal_tau))
    counts = (
        interval_count, mapped_count, active, q_active, weighted,
        tuple(state_counts), update_count, tuple(len(row) for row in maps),
    )
    return tuple(maps), counts


def expand_map(coefficient_map):
    bank = [[[0] * P for _r0 in range(P)] for _state in range(V)]
    for (cell_index, selected_u_mask), scalar in coefficient_map.items():
        if not scalar:
            continue
        for state, r0, r1, value in sparse_left_vector(cell_index, selected_u_mask):
            bank[state][r0][r1] = (
                bank[state][r0][r1] + value * scalar
            ) % PRIME
    return freeze3(bank)


def phase_row(row, phase: int):
    return tuple(tuple(tuple(phase * value % PRIME for value in r1_row)
                       for r1_row in r0_rows) for r0_rows in row)


def r1_row_marginal(row):
    return tuple(tuple(sum(row[state][r0]) % PRIME for r0 in range(P))
                 for state in range(V))


def worker(alpha: int):
    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    beta_blocks = []
    parent_rows = []
    alpha_tau_sum = [[[[0 for _r1 in range(P)] for _r0 in range(P)]
                      for _state in range(V)] for _tau in range(P)]
    scalar_counts = [0] * 5
    state_counts = [0] * V
    update_count = 0
    key_counts = [0] * P
    controls = {}
    probe_values = [[[] for _relation in (0,)] for _probe in INVERSION_PROBES]

    for beta in range(P):
        maps, counts = coefficient_maps(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        tau_rows = []
        for tau in range(P):
            phased = phase_row(expand_map(maps[tau]), phase)
            tau_rows.append(phased)
            parent_rows.append(r1_row_marginal(phased))
            if (alpha, beta, tau) in CONTROL_TRIPLES:
                controls[(alpha, beta, tau)] = phased
            for state in range(V):
                for r0 in range(P):
                    for r1 in range(P):
                        alpha_tau_sum[tau][state][r0][r1] = (
                            alpha_tau_sum[tau][state][r0][r1]
                            + phased[state][r0][r1]
                        ) % PRIME
            for probe_index, (state, r0, r1) in enumerate(INVERSION_PROBES):
                probe_values[probe_index][0].append(phased[state][r0][r1])
        beta_blocks.append(tuple(tau_rows))
        scalar_counts = [left + right for left, right in
                         zip(scalar_counts, counts[:5])]
        state_counts = [left + right for left, right in
                        zip(state_counts, counts[5])]
        update_count += counts[6]
        key_counts = [left + right for left, right in zip(key_counts, counts[7])]

    frozen_tau_sum = freeze4(alpha_tau_sum)
    probes = tuple((probe, tuple(probe_values[index][0]))
                   for index, probe in enumerate(INVERSION_PROBES))
    return (
        alpha,
        frozen_tau_sum,
        tuple(parent_rows),
        digest_json(beta_blocks),
        (tuple(scalar_counts), tuple(state_counts), update_count, tuple(key_counts)),
        tuple(sorted(controls.items())),
        probes,
    )


def invert_worker_sums(chunks, zeta: int):
    tau_core = [[[[0 for _r1 in range(P)] for _r0 in range(P)]
                 for _state in range(V)] for _tau in range(P)]
    for chunk in chunks:
        alpha = chunk[0]
        alpha_phase = pow(zeta, -alpha % P, PRIME)
        for tau in range(P):
            for state in range(V):
                for r0 in range(P):
                    for r1 in range(P):
                        tau_core[tau][state][r0][r1] = (
                            tau_core[tau][state][r0][r1]
                            + alpha_phase * chunk[1][tau][state][r0][r1]
                        ) % PRIME
    frozen_core = freeze4(tau_core)
    normalizer = pow(P**3, -1, PRIME)
    tensor = [[[[0] * P for _r1 in range(P)] for _r0 in range(P)]
              for _state in range(V)]
    for state in range(V):
        for r0 in range(P):
            for r1 in range(P):
                for relation in range(P):
                    tensor[state][r0][r1][relation] = (
                        sum(tau_core[tau][state][r0][r1]
                            * pow(zeta, -tau * relation % P, PRIME)
                            for tau in range(P))
                        * normalizer
                    ) % PRIME
    return freeze4(tensor), frozen_core


def marginal_r1(tensor):
    return tuple(tuple(tuple(
        sum(tensor[state][r0][r1][relation] for r1 in range(P)) % PRIME
        for relation in range(P)
    ) for r0 in range(P)) for state in range(V))


def marginal_r0(first_digit):
    return tuple(tuple(
        sum(first_digit[state][r0][relation] for r0 in range(P)) % PRIME
        for relation in range(P)
    ) for state in range(V))


def flat_r1_lift(first_digit):
    inverse = pow(P, -1, PRIME)
    return tuple(tuple(tuple(tuple(
        first_digit[state][r0][relation] * inverse % PRIME
        for relation in range(P)
    ) for _r1 in range(P)) for r0 in range(P)) for state in range(V))


def centre_axis(tensor, axis: int):
    dimensions = (V, P, P, P)
    inverse = pow(dimensions[axis], -1, PRIME)
    answer = [[[[tensor[state][r0][r1][relation] for relation in range(P)]
                for r1 in range(P)] for r0 in range(P)] for state in range(V)]
    for state in range(V):
        for r0 in range(P):
            for r1 in range(P):
                for relation in range(P):
                    if axis == 0:
                        total = sum(tensor[s][r0][r1][relation] for s in range(V))
                    elif axis == 1:
                        total = sum(tensor[state][d0][r1][relation] for d0 in range(P))
                    elif axis == 2:
                        total = sum(tensor[state][r0][d1][relation] for d1 in range(P))
                    else:
                        total = sum(tensor[state][r0][r1][t] for t in range(P))
                    answer[state][r0][r1][relation] = (
                        tensor[state][r0][r1][relation] - total * inverse
                    ) % PRIME
    return freeze4(answer)


def four_way_interaction(tensor):
    answer = tensor
    for axis in range(4):
        answer = centre_axis(answer, axis)
    return answer


def axis_ranks(tensor):
    matrices = (
        tuple(tuple(tensor[state][r0][r1][relation]
                    for r0 in range(P) for r1 in range(P) for relation in range(P))
              for state in range(V)),
        tuple(tuple(tensor[state][r0][r1][relation]
                    for state in range(V) for r1 in range(P) for relation in range(P))
              for r0 in range(P)),
        tuple(tuple(tensor[state][r0][r1][relation]
                    for state in range(V) for r0 in range(P) for relation in range(P))
              for r1 in range(P)),
        tuple(tuple(tensor[state][r0][r1][relation]
                    for state in range(V) for r0 in range(P) for r1 in range(P))
              for relation in range(P)),
        tuple(tuple(tensor[state][r0][r1][relation]
                    for state in range(V) for relation in range(P))
              for r0 in range(P) for r1 in range(P)),
    )
    return tuple(R.rank_mod(matrix) for matrix in matrices)


def conditional_ranks(tensor):
    raw = []
    contrast = []
    inverse = pow(P, -1, PRIME)
    for r0 in range(P):
        matrix = tuple(tuple(tensor[state][r0][r1][relation]
                             for state in range(V) for relation in range(P))
                       for r1 in range(P))
        means = tuple(sum(matrix[r1][column] for r1 in range(P)) * inverse % PRIME
                      for column in range(V * P))
        centred = tuple(tuple((matrix[r1][column] - means[column]) % PRIME
                              for column in range(V * P)) for r1 in range(P))
        raw.append(R.rank_mod(matrix))
        contrast.append(R.rank_mod(centred))
    return tuple(raw), tuple(contrast)


def fourier_tensor(tensor, zeta: int):
    roots = tuple(tuple(pow(zeta, -frequency * value % P, PRIME)
                        for value in range(P)) for frequency in range(P))
    state_stage = [[[[sum(
        tensor[state][r0][r1][relation] * WALSH_SIGNS[character][state]
        for state in range(V)
    ) % PRIME for relation in range(P)] for r1 in range(P)] for r0 in range(P)]
        for character in range(V)]
    r0_stage = [[[[sum(
        state_stage[character][r0][r1][relation] * roots[f0][r0]
        for r0 in range(P)
    ) % PRIME for relation in range(P)] for r1 in range(P)] for f0 in range(P)]
        for character in range(V)]
    r1_stage = [[[[sum(
        r0_stage[character][f0][r1][relation] * roots[f1][r1]
        for r1 in range(P)
    ) % PRIME for relation in range(P)] for f1 in range(P)] for f0 in range(P)]
        for character in range(V)]
    return tuple(tuple(tuple(tuple(sum(
        r1_stage[character][f0][f1][relation] * roots[ft][relation]
        for relation in range(P)
    ) % PRIME for ft in range(P)) for f1 in range(P)) for f0 in range(P))
        for character in range(V))


def support_census(spectrum):
    bins = [0] * 16
    for character in range(V):
        for f0 in range(P):
            for f1 in range(P):
                for ft in range(P):
                    if spectrum[character][f0][f1][ft]:
                        mask = (
                            ((character != 0) << 3)
                            | ((f0 != 0) << 2)
                            | ((f1 != 0) << 1)
                            | (ft != 0)
                        )
                        bins[mask] += 1
    return (sum(bins),) + tuple(bins)


def fixed_relation_record(tensor, relation: int):
    matrix = tuple(tuple(tensor[state][r0][r1][relation]
                         for r0 in range(P) for r1 in range(P))
                   for state in range(V))
    conditionals = tuple(R.rank_mod(tuple(
        tuple(tensor[state][r0][r1][relation] for state in range(V))
        for r1 in range(P)
    )) for r0 in range(P))
    return (
        R.rank_mod(matrix),
        tuple(sum(value != 0 for value in row) for row in matrix),
        conditionals,
        digest_json(matrix),
    )


def direct_inversion_probes(chunks, tensor, zeta: int):
    normalizer = pow(P**3, -1, PRIME)
    records = []
    for probe_index, probe in enumerate(INVERSION_PROBES):
        state, r0, r1 = probe
        values = []
        for relation in PROBE_RELATIONS:
            total = 0
            for chunk in chunks:
                alpha = chunk[0]
                observed_probe, rows = chunk[6][probe_index]
                require(observed_probe == probe,
                        ("probe order", observed_probe, probe))
                require(len(rows) == P**2, ("probe row count", probe, len(rows)))
                index = 0
                for _beta in range(P):
                    for tau in range(P):
                        total = (
                            total + rows[index]
                            * pow(zeta, -(alpha + tau * relation) % P, PRIME)
                        ) % PRIME
                        index += 1
            value = total * normalizer % PRIME
            require(value == tensor[state][r0][r1][relation],
                    ("direct inversion", probe, relation, value,
                     tensor[state][r0][r1][relation]))
            values.append(value)
        records.append((probe, tuple(values)))
    return tuple(records)


def main() -> None:
    R.split_field_certificate()
    (
        _profiles, boundaries, _cells, _one_digit, _source_u, _source_v,
        source_boundaries, source_profile_digest, source_record, source_sha,
        reflection_record,
    ) = source_context()
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "worker order")

    alpha_digests = tuple(chunk[3] for chunk in chunks)
    require(alpha_digests == EXPECTED_ALPHA_GAMMA_SHA256,
            ("alpha gamma digests", alpha_digests))
    gamma_parent = tuple(row for chunk in chunks for row in chunk[2])
    require(digest_json(gamma_parent) == ONE.EXPECTED_DIGESTS[0],
            "gamma r1 marginal misses audited r0 parent")
    scalar_counts = tuple(sum(chunk[4][0][index] for chunk in chunks)
                          for index in range(5))
    state_counts = tuple(sum(chunk[4][1][state] for chunk in chunks)
                         for state in range(V))
    update_count = sum(chunk[4][2] for chunk in chunks)
    key_counts = tuple(sum(chunk[4][3][tau] for chunk in chunks)
                       for tau in range(P))
    require(scalar_counts == EXPECTED_WORK_COUNTS,
            ("work counts", scalar_counts))
    require(state_counts == EXPECTED_STATE_COUNTS,
            ("state counts", state_counts))
    require(update_count == EXPECTED_COEFFICIENT_UPDATES,
            ("coefficient updates", update_count))
    require(key_counts == EXPECTED_KEY_COUNTS, ("key counts", key_counts))

    zeta = pow(R.JOINT_ROOT, R.JOINT_ORDER // P, PRIME)
    tensor, tau_core = invert_worker_sums(chunks, zeta)
    first_digit = marginal_r1(tensor)
    square = marginal_r0(first_digit)
    require(digest_json(first_digit) == ONE.EXPECTED_DIGESTS[1],
            "tensor r1 marginal misses audited r0 parent")
    require(digest_json(square) == ONE.EXPECTED_PARENT_TABLE,
            "tensor r0 marginal misses audited square")
    flat = flat_r1_lift(first_digit)
    require(marginal_r1(flat) == first_digit, "flat r1 marginal")

    control_records = []
    word, endpoint_grid, _harmonic, danger = endpoint_context()
    literal_boundaries = (
        set(boundaries) | R.fixed_boundaries(source_boundaries, R.SRC.T_DEN)
    )
    stored_controls = dict(item for chunk in chunks for item in chunk[5])
    for alpha, beta, tau in CONTROL_TRIPLES:
        restoration = R.literal_guard_restoration(
            word, endpoint_grid, alpha, beta, tau, literal_boundaries, danger
        )
        maps, counts = coefficient_maps(alpha, beta, tau)
        phase = pow(zeta, beta, PRIME)
        direct = phase_row(expand_map(maps[0]), phase)
        require(direct == stored_controls[(alpha, beta, tau)],
                ("literal guard", alpha, beta, tau))
        control_records.append(((alpha, beta, tau), restoration, counts))

    probes = direct_inversion_probes(chunks, tensor, zeta)
    require(probes == EXPECTED_PROBES, ("inversion probes", probes))
    tensors = (tensor, flat)
    interactions = tuple(four_way_interaction(value) for value in tensors)
    ranks = tuple(axis_ranks(value) for value in tensors)
    interaction_ranks = tuple(axis_ranks(value) for value in interactions)
    conditional = tuple(conditional_ranks(value) for value in tensors)
    require((ranks, interaction_ranks) == EXPECTED_RANKS,
            ("axis ranks", ranks, interaction_ranks))
    require(conditional == EXPECTED_CONDITIONAL_RANKS,
            ("conditional ranks", conditional))
    require(all(rank > 1 for rank in conditional[0][0]),
            ("scalar-kernel hostile", conditional[0][0]))

    spectra = tuple(fourier_tensor(value, zeta) for value in tensors)
    interaction_spectra = tuple(fourier_tensor(value, zeta)
                                for value in interactions)
    censuses = tuple(support_census(value) for value in spectra)
    interaction_censuses = tuple(support_census(value)
                                 for value in interaction_spectra)
    require((censuses, interaction_censuses) == EXPECTED_CENSUSES,
            ("spectral censuses", censuses, interaction_censuses))
    fixed = tuple(fixed_relation_record(value, 6) for value in tensors)
    require(fixed == EXPECTED_FIXED, ("fixed t6", fixed))

    digests = (
        digest_json(tensor),
        digest_json(tau_core),
        digest_json(first_digit),
        digest_json(square),
        digest_json(flat),
        tuple(digest_json(value) for value in interactions),
        tuple(digest_json(value) for value in spectra),
        tuple(digest_json(value) for value in interaction_spectra),
    )
    require(digests == EXPECTED_DIGESTS, ("candidate digests", digests))
    rank_record = (ranks, interaction_ranks)
    census_record = (censuses, interaction_censuses)
    record = (
        PARENT_SHA256, PARENT_SEMANTIC, source_profile_digest, source_record,
        source_sha, reflection_record, scalar_counts, state_counts,
        update_count, key_counts, tuple(control_records), probes,
        alpha_digests, digests, rank_record, conditional, census_record, fixed,
        "scalar state/relation-independent K(r0,r1) factorization only",
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== independent hostile audit: two current inverse digits x Boolean square x relation ==")
    print(f"parent=(one_digit_sha256={PARENT_SHA256},semantic={PARENT_SEMANTIC})")
    print("digits=a=r0+13*r1+169*c;construction=independent 13^3 fold -> r1 -> r0 -> root")
    print("half_open_nested_window_equals_direct_mod169_window: PASS")
    print(f"source=(f_pieces={source_record[2]},fold_pieces={source_record[3]},r1_pieces={source_record[4]},address_pieces={source_record[5]},root_pieces={source_record[6]},boundaries={source_record[7]},active_root_support={source_record[8]},sha256={source_sha})")
    print(f"reflection=(u,r0,r1,state)->(12-u,12-r0,12-r1,state_XOR2);intervals={reflection_record[:2]};mass_sha256={digest_json(reflection_record[2])}: PASS")
    print(f"field=(prime={R.JOINT_PRIME},root={R.JOINT_ROOT},zeta13={zeta})")
    print(f"sparse_work=(counts={scalar_counts},state_counts={state_counts},updates={update_count},keys_by_tau={key_counts})")
    print(f"literal_controls={tuple(row[0] for row in control_records)}: PASS;same_root=0")
    print(f"direct_2197_term_inversion_probes={probes}: PASS")
    print("marginal_chain=(sum_r1->independently_audited_r0 gamma/tensor,sum_r0->audited_square): PASS")
    print(f"axis_ranks_(state,r0,r1,relation,address)=(actual,flat)={ranks}")
    print(f"four_way_ANOVA_axis_ranks={interaction_ranks}")
    print(f"conditional_r1_by_(state,relation)_ranks_(raw,contrast)=(actual,flat)={conditional}")
    print("amplitude_carrier_typing=conditional ranks are r1-to-(state,relation) amplitude ranks;pointed-six coordinate marginalized;no carrier or recurrence conclusion")
    print("hostile_scope=actual rank>1 refutes only scalar state/relation-independent K(r0,r1)*T1 factorization")
    print("nonclaim=does_not_refute_hidden_state,matrix_valued,state_dependent,nonlinear,or_temporal_recurrence")
    print("spectral_census_order=(total,16 bins by state/r0/r1/relation mask 0000..1111)")
    print(f"spectral_censuses_(actual,flat)={censuses}")
    print(f"four_way_ANOVA_spectral_censuses={interaction_censuses}")
    print(f"fixed_(1,0,6)_records={fixed}")
    print(f"alpha_gamma_digests={alpha_digests}")
    print(f"digests=(tensor,tau_core,r0_parent,square,flat,ANOVA,spectra,ANOVA_spectra)={digests}")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT scoped two-static-current-digit projective cylinder and scalar-factor obstruction")
    print("scope=not complete address,not general autonomous-recurrence no-go,not U_clock chronology,not arrival ancestry,not physical current,no row exclusion,no LRC(14)")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
