#!/usr/bin/env python3
"""Retain two consecutive current-leg inverse digits over the Boolean square.

At the THM-2471 r=5 stalk, write

    X_(u,a) = (w_u+a)/13^5,
    a = r0 + 13*r1 + 13^2*c.

Here r0 is the already audited current-leg owner digit and r1 is its
immediate predecessor.  The exact construction folds by 13^3, takes the r1
window, then the r0 window, and finally the collision-root window.  Half-open
windows make this identical to taking the single mod-169 window r0+13*r1.

The implementation is sparse in source geometry rather than dense in endpoint
segments.  Endpoint contributions are first accumulated by the 34 exact
source cells and the selected U-root mask; only then are the 169 address
weights expanded.  Summing r1 must recover the hash-pinned one-digit candidate
row by row, and summing r0 next must recover the audited Boolean square.

The decisive hostile is stronger than a flat split.  Any scalar,
state/relation-independent one-step digit factorization has the form

    T(state,r0,r1,t) = K(r0,r1) T1(state,r0,t),

so for each fixed r0 its r1-by-(state,t) conditional matrix has rank at most
one.  The script compares those ranks with the exact depth-two tensor.

This is a finite exact two-cylinder probe.  It is not a complete inverse
address, a U_clock chronology, an arrival atom, a physical current, a row
exclusion, or a proof of LRC(14).
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
CURRENT_PATH = (
    ROOT
    / "04-computation/lrc_r5_ufull_owner_node_boolean_square_inverse_owner_branch_probe_20260816.py"
)
CURRENT_SHA256 = "ae1cf021ea23f325eeded42ff1dea8df903d837b1d7ab551289d62f7ab7a0348"

P = 13
V = 4
DIRECT_CONTROLS = ((0, 0, 0), (1, 0, 6), (6, 6, 12))
INVERSION_PROBES = ((0, 0, 1), (1, 6, 7), (2, 4, 9), (3, 12, 11))
PROBE_RELATIONS = (0, 6, 12)

# Pinned after the first exact normal replay.
EXPECTED_SOURCE_SHA256 = "dd2f48375fc38b419babdd0ed13365fb56918829a663270c4b23d56635d6e097"
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
    ((1,) * P, (0,) * P),
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
EXPECTED_SEMANTIC_SHA256 = "61743457afc0cff984c87affa7f2e67bf3a21e08a401ea69b679319f2f51e826"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_current():
    observed = lf_sha256(CURRENT_PATH)
    require(observed == CURRENT_SHA256,
            ("one-digit current source drift", observed, CURRENT_SHA256))
    spec = importlib.util.spec_from_file_location("two_digit_current_parent", CURRENT_PATH)
    require(spec is not None and spec.loader is not None, "two-digit parent loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


D = load_current()
B = D.B
C = D.C
PRIME = C.JOINT_PRIME


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


def freeze4(bank):
    return tuple(
        tuple(tuple(tuple(row) for row in plane) for plane in slab)
        for slab in bank
    )


@lru_cache(maxsize=1)
def two_digit_context():
    ctx = C.context()
    stage = C.SM
    source_grid = C.S.GRID
    joint_grid = C.JOINT_COORDINATE
    require(stage.DCOLL == P**5, ("collision depth", stage.DCOLL))

    e_intervals = stage.build_set(stage.PAT_E, stage.ZELL)
    q_intervals = stage.build_set(stage.PAT_QB, stage.ZELL)
    f_pieces = C.S.B.build_f_pieces(e_intervals, q_intervals)

    # a=r0+13*r1+13^2*c.  Window composition is contravariant: the
    # higher retained digit r1 is extracted before the lower digit r0.
    base_fold = stage.weighted_fold(f_pieces, stage.DCOLL // P**2, source_grid)
    r1_windows = tuple(
        stage.extract_window(base_fold[0], base_fold[1], r1, P, source_grid)
        for r1 in range(P)
    )
    address_windows = tuple(
        tuple(
            stage.extract_window(
                r1_windows[r1][0], r1_windows[r1][1], r0, P, source_grid
            )
            for r1 in range(P)
        )
        for r0 in range(P)
    )

    # The nested half-open windows must literally equal one mod-169 window.
    for r0 in range(P):
        for r1 in range(P):
            combined = stage.extract_window(
                base_fold[0], base_fold[1], r0 + P * r1, P**2, source_grid
            )
            require(address_windows[r0][r1] == combined,
                    ("nested/combined mod-169 window", r0, r1))

    root_raw = tuple(
        tuple(
            tuple(
                stage.extract_window(
                    address_windows[r0][r1][0],
                    address_windows[r0][r1][1],
                    root,
                    P,
                    source_grid,
                )
                for r1 in range(P)
            )
            for r0 in range(P)
        )
        for root in range(P)
    )
    root_profiles = tuple(
        tuple(
            tuple(
                C.scale_profile(profile, joint_grid, source_grid)
                for profile in r1_row
            )
            for r1_row in r0_rows
        )
        for r0_rows in root_raw
    )

    current_profiles, _current_boundaries, _current_record = D.owner_branch_context()
    boundaries = tuple(sorted(
        {0, joint_grid}
        | set(ctx["source_boundaries"])
        | {
            position
            for root_rows in root_profiles
            for r0_rows in root_rows
            for profile in r0_rows
            for position in profile[0]
        }
    ))
    require(
        tuple(joint_grid - point for point in reversed(boundaries)) == boundaries,
        "two-digit source-boundary reflection",
    )

    cells = []
    state_address_support = [set() for _state in range(V)]
    for left, right in zip(boundaries, boundaries[1:]):
        address_values = tuple(
            tuple(
                tuple(profile_value(root_profiles[root][r0][r1], left)
                      for r1 in range(P))
                for r0 in range(P)
            )
            for root in range(P)
        )
        u_values = tuple(profile_value(profile, left) for profile in ctx["source_u"])
        v_values = tuple(profile_value(profile, left) for profile in ctx["source_v"])

        for root in range(P):
            for r0 in range(P):
                first_digit = sum(address_values[root][r0])
                expected_first = profile_value(current_profiles[root][r0], left)
                require(first_digit == expected_first,
                        ("r1 cylinder marginal", left, right, root, r0,
                         first_digit, expected_first))
            require(
                sum(address_values[root][r0][r1]
                    for r0 in range(P) for r1 in range(P)) == u_values[root],
                ("mod-169 total marginal", left, right, root),
            )
            require(
                all(address_values[root][r0][r1] * v_values[root] == 0
                    for r0 in range(P) for r1 in range(P)),
                ("two-digit same-root source overlap", left, right, root),
            )
            for r0 in range(P):
                for r1 in range(P):
                    reflected = profile_value(
                        root_profiles[P - 1 - root][P - 1 - r0][P - 1 - r1],
                        joint_grid - right,
                    )
                    require(address_values[root][r0][r1] == reflected,
                            ("two-digit profile reflection", left, right,
                             root, r0, r1))

        state = None
        try:
            state = B.state_of_segment(left, right)
        except RuntimeError:
            pass
        if state is not None:
            reflected_state = B.state_of_segment(joint_grid - right, joint_grid - left)
            require(reflected_state == state ^ 2,
                    ("two-digit state reflection", left, right, state, reflected_state))
            for root in range(P):
                count = sum(
                    address_values[root][r0][r1] != 0
                    for r0 in range(P) for r1 in range(P)
                )
                if count:
                    state_address_support[state].add(count)

        u_support_mask = sum(1 << root for root, value in enumerate(u_values) if value)
        cells.append((state, u_values, v_values, u_support_mask, address_values))

    require(all(values == {132} for values in state_address_support),
            ("two-digit geometric support", state_address_support))

    address_masses = tuple(
        sum(profile_mass(root_raw[root][r0][r1], source_grid) for root in range(P))
        for r0 in range(P) for r1 in range(P)
    )
    reflected_masses = tuple(
        address_masses[(P - 1 - r0) * P + (P - 1 - r1)]
        for r0 in range(P) for r1 in range(P)
    )
    require(address_masses == reflected_masses,
            "two-digit address-mass reflection")

    source_record = (
        stage.DCOLL,
        stage.DCOLL // P**2,
        len(f_pieces),
        len(base_fold[0]),
        tuple(len(profile[0]) for profile in r1_windows),
        sum(len(profile[0]) for row in address_windows for profile in row),
        sum(len(profile[0]) for root_rows in root_raw
            for row in root_rows for profile in row),
        len(boundaries),
        tuple(tuple(sorted(values)) for values in state_address_support),
        address_masses,
        digest_json(root_raw),
        digest_json(root_profiles),
    )
    return root_profiles, boundaries, tuple(cells), source_record


SAFE_MASKS = {
    chamber: tuple(
        sum(1 << sheet for sheet in range(P) if C.T.safe(chamber, sheet + tau))
        for tau in range(P)
    )
    for chamber in ("left", "right")
}


@lru_cache(maxsize=200000)
def selection_data(cell_index: int, selected_mask: int):
    _profiles, _boundaries, cells, _record = two_digit_context()
    _state, _u_values, v_values, u_support_mask, _address_values = cells[cell_index]
    right_value = sum(
        v_values[root] for root in range(P) if (selected_mask >> root) & 1
    )
    return selected_mask & u_support_mask, right_value


def expand_coefficients(coefficients):
    _profiles, _boundaries, cells, _record = two_digit_context()
    bank = [[[[0 for _r1 in range(P)] for _r0 in range(P)]
             for _state in range(V)] for _tau in coefficients]
    for tau_index, coefficient_map in enumerate(coefficients):
        for (cell_index, selected_u_mask), scalar in coefficient_map.items():
            if scalar == 0 or selected_u_mask == 0:
                continue
            state, _u_values, _v_values, _u_support, address_values = cells[cell_index]
            require(state is not None, ("untyped expanded source cell", cell_index))
            roots = tuple(
                root for root in range(P) if (selected_u_mask >> root) & 1
            )
            for r0 in range(P):
                for r1 in range(P):
                    left_value = sum(address_values[root][r0][r1] for root in roots)
                    if left_value:
                        bank[tau_index][state][r0][r1] = (
                            bank[tau_index][state][r0][r1] + left_value * scalar
                        ) % PRIME
    return freeze4(bank)


def integrate_two_digits(alpha: int, beta: int, literal_tau: int | None = None):
    ctx = C.context()
    _profiles, boundaries, cells, _source_record = two_digit_context()
    events, interval_count, mapped = C.endpoint_events(alpha, beta, literal_tau)
    for boundary in boundaries:
        events.setdefault(boundary, 0)
    tau_values = tuple(range(P)) if literal_tau is None else (literal_tau,)
    coefficients = [dict() for _tau in tau_values]
    mask = 0
    positions = sorted(events)
    active_segments = q_active_segments = weighted_segments = 0
    state_segments = [0] * V
    coefficient_updates = 0

    for position_index, left in enumerate(positions[:-1]):
        mask ^= events[left]
        right = positions[position_index + 1]
        if left == right or mask == 0:
            continue
        require(C.cell_of_segment(left, right) == 0,
                ("two-digit branch escaped cell zero", alpha, beta, left, right))
        chamber = C.chamber_of_segment(left, right)
        require(chamber in ("left", "right"),
                ("two-digit active middle chamber", alpha, beta, left, right))
        cell_index = bisect_right(boundaries, left) - 1
        require(0 <= cell_index < len(cells)
                and right <= boundaries[cell_index + 1],
                ("endpoint segment crosses source cell", left, right, cell_index))
        state, u_values, v_values, _u_support_mask, _address_values = cells[cell_index]
        require(state is not None, ("untyped active source cell", left, right))
        active_segments += 1
        state_segments[state] += 1

        jump = C.q_endpoint_jump(left, right)
        if not jump:
            continue
        q_active_segments += 1
        if not any(u_values) or not any(v_values):
            continue
        weighted_segments += 1

        for tau_index, tau in enumerate(tau_values):
            selected_mask = mask if literal_tau is not None else mask & SAFE_MASKS[chamber][tau]
            selected_u_mask, right_value = selection_data(cell_index, selected_mask)
            if selected_u_mask == 0 or right_value == 0:
                continue
            key = (cell_index, selected_u_mask)
            coefficients[tau_index][key] = (
                coefficients[tau_index].get(key, 0) + right_value * jump
            ) % PRIME
            coefficient_updates += 1

    mask ^= events[positions[-1]]
    require(mask == 0, ("two-digit endpoint mask", alpha, beta, literal_tau, mask))
    bank = expand_coefficients(coefficients)
    counts = (
        interval_count,
        mapped,
        active_segments,
        q_active_segments,
        weighted_segments,
        tuple(state_segments),
        coefficient_updates,
        tuple(len(row) for row in coefficients),
    )
    return bank, counts


def marginal_r1(tensor):
    return tuple(
        tuple(
            tuple(sum(tensor[state][r0][r1][relation] for r1 in range(P)) % PRIME
                  for relation in range(P))
            for r0 in range(P)
        )
        for state in range(V)
    )


def marginal_r0(first_digit_tensor):
    return tuple(
        tuple(sum(first_digit_tensor[state][r0][relation] for r0 in range(P)) % PRIME
              for relation in range(P))
        for state in range(V)
    )


def marginal_gamma_r1(row):
    return tuple(
        tuple(sum(row[state][r0]) % PRIME for r0 in range(P))
        for state in range(V)
    )


def phase_bank(bank, phase: int):
    return tuple(
        tuple(
            tuple(phase * value % PRIME for value in r1_row)
            for r1_row in state_rows
        )
        for state_rows in bank
    )


def worker(alpha: int):
    zeta = C.context()["zeta"]
    alpha_rows = []
    marginal_rows = []
    tau_sum = [[[[0 for _r1 in range(P)] for _r0 in range(P)]
                for _state in range(V)] for _tau in range(P)]
    scalar_counts = [0] * 5
    state_counts = [0] * V
    update_count = 0
    key_counts = [0] * P
    control_rows = {}

    for beta in range(P):
        bank, counts = integrate_two_digits(alpha, beta)
        phase = pow(zeta, beta, PRIME)
        phased = tuple(phase_bank(row, phase) for row in bank)
        alpha_rows.append(phased)
        for tau in range(P):
            marginal_rows.append(marginal_gamma_r1(phased[tau]))
            if (alpha, beta, tau) in DIRECT_CONTROLS:
                control_rows[(alpha, beta, tau)] = phased[tau]
            for state in range(V):
                for r0 in range(P):
                    for r1 in range(P):
                        tau_sum[tau][state][r0][r1] = (
                            tau_sum[tau][state][r0][r1]
                            + phased[tau][state][r0][r1]
                        ) % PRIME
        scalar_counts = [left + right for left, right in zip(scalar_counts, counts[:5])]
        state_counts = [left + right for left, right in zip(state_counts, counts[5])]
        update_count += counts[6]
        key_counts = [left + right for left, right in zip(key_counts, counts[7])]

    probe_rows = tuple(
        (
            probe,
            tuple(
                alpha_rows[beta][tau][probe[0]][probe[1]][probe[2]]
                for beta in range(P) for tau in range(P)
            ),
        )
        for probe in INVERSION_PROBES
    )
    return (
        alpha,
        freeze4(tau_sum),
        tuple(marginal_rows),
        digest_json(alpha_rows),
        (tuple(scalar_counts), tuple(state_counts), update_count, tuple(key_counts)),
        tuple(sorted(control_rows.items())),
        probe_rows,
    )


def inverse_compressed(chunks, zeta: int):
    # beta is already twisted and summed inside each alpha worker.  Sum alpha
    # with zeta^-alpha, then perform only the remaining tau DFT.
    tau_core = [[[[0 for _r1 in range(P)] for _r0 in range(P)]
                 for _state in range(V)] for _tau in range(P)]
    for chunk in chunks:
        alpha = chunk[0]
        phase = pow(zeta, -alpha % P, PRIME)
        alpha_bank = chunk[1]
        for tau in range(P):
            for state in range(V):
                for r0 in range(P):
                    for r1 in range(P):
                        tau_core[tau][state][r0][r1] = (
                            tau_core[tau][state][r0][r1]
                            + phase * alpha_bank[tau][state][r0][r1]
                        ) % PRIME

    normalizer = pow(P**3, -1, PRIME)
    roots = tuple(
        tuple(pow(zeta, -tau * relation % P, PRIME) for tau in range(P))
        for relation in range(P)
    )
    tensor = [[[[0 for _relation in range(P)] for _r1 in range(P)]
               for _r0 in range(P)] for _state in range(V)]
    for state in range(V):
        for r0 in range(P):
            for r1 in range(P):
                values = tuple(tau_core[tau][state][r0][r1] for tau in range(P))
                for relation in range(P):
                    tensor[state][r0][r1][relation] = (
                        sum(values[tau] * roots[relation][tau] for tau in range(P))
                        * normalizer
                    ) % PRIME
    return freeze4(tensor), freeze4(tau_core)


def flat_second_digit_lift(first_digit):
    inverse = pow(P, -1, PRIME)
    return tuple(
        tuple(
            tuple(
                tuple(first_digit[state][r0][relation] * inverse % PRIME
                      for relation in range(P))
                for _r1 in range(P)
            )
            for r0 in range(P)
        )
        for state in range(V)
    )


def centre_axis(tensor, axis: int):
    dimensions = (V, P, P, P)
    inverse = pow(dimensions[axis], -1, PRIME)
    mutable = [[[[tensor[state][r0][r1][relation] for relation in range(P)]
                 for r1 in range(P)] for r0 in range(P)] for state in range(V)]
    if axis == 0:
        for r0 in range(P):
            for r1 in range(P):
                for relation in range(P):
                    mean = sum(mutable[state][r0][r1][relation]
                               for state in range(V)) * inverse % PRIME
                    for state in range(V):
                        mutable[state][r0][r1][relation] = (
                            mutable[state][r0][r1][relation] - mean
                        ) % PRIME
    elif axis == 1:
        for state in range(V):
            for r1 in range(P):
                for relation in range(P):
                    mean = sum(mutable[state][r0][r1][relation]
                               for r0 in range(P)) * inverse % PRIME
                    for r0 in range(P):
                        mutable[state][r0][r1][relation] = (
                            mutable[state][r0][r1][relation] - mean
                        ) % PRIME
    elif axis == 2:
        for state in range(V):
            for r0 in range(P):
                for relation in range(P):
                    mean = sum(mutable[state][r0][r1][relation]
                               for r1 in range(P)) * inverse % PRIME
                    for r1 in range(P):
                        mutable[state][r0][r1][relation] = (
                            mutable[state][r0][r1][relation] - mean
                        ) % PRIME
    else:
        for state in range(V):
            for r0 in range(P):
                for r1 in range(P):
                    mean = sum(mutable[state][r0][r1][relation]
                               for relation in range(P)) * inverse % PRIME
                    for relation in range(P):
                        mutable[state][r0][r1][relation] = (
                            mutable[state][r0][r1][relation] - mean
                        ) % PRIME
    return freeze4(mutable)


def four_way_interaction(tensor):
    answer = tensor
    for axis in range(4):
        answer = centre_axis(answer, axis)
    return answer


def axis_ranks(tensor):
    state_flat = tuple(
        tuple(tensor[state][r0][r1][relation]
              for r0 in range(P) for r1 in range(P) for relation in range(P))
        for state in range(V)
    )
    r0_flat = tuple(
        tuple(tensor[state][r0][r1][relation]
              for state in range(V) for r1 in range(P) for relation in range(P))
        for r0 in range(P)
    )
    r1_flat = tuple(
        tuple(tensor[state][r0][r1][relation]
              for state in range(V) for r0 in range(P) for relation in range(P))
        for r1 in range(P)
    )
    relation_flat = tuple(
        tuple(tensor[state][r0][r1][relation]
              for state in range(V) for r0 in range(P) for r1 in range(P))
        for relation in range(P)
    )
    address_flat = tuple(
        tuple(tensor[state][r0][r1][relation]
              for state in range(V) for relation in range(P))
        for r0 in range(P) for r1 in range(P)
    )
    return (
        C.rank_mod(state_flat),
        C.rank_mod(r0_flat),
        C.rank_mod(r1_flat),
        C.rank_mod(relation_flat),
        C.rank_mod(address_flat),
    )


def conditional_ranks(tensor):
    ranks = []
    contrast_ranks = []
    inverse = pow(P, -1, PRIME)
    for r0 in range(P):
        matrix = tuple(
            tuple(tensor[state][r0][r1][relation]
                  for state in range(V) for relation in range(P))
            for r1 in range(P)
        )
        contrasts = []
        means = tuple(
            sum(matrix[r1][column] for r1 in range(P)) * inverse % PRIME
            for column in range(V * P)
        )
        for r1 in range(P):
            contrasts.append(tuple(
                (matrix[r1][column] - means[column]) % PRIME
                for column in range(V * P)
            ))
        ranks.append(C.rank_mod(matrix))
        contrast_ranks.append(C.rank_mod(tuple(contrasts)))
    return tuple(ranks), tuple(contrast_ranks)


def fourier_tensor(tensor, zeta: int):
    roots = tuple(
        tuple(pow(zeta, -frequency * value % P, PRIME) for value in range(P))
        for frequency in range(P)
    )
    state_stage = [[[[
        sum(tensor[state][r0][r1][relation] * B.WALSH_SIGNS[character][state]
            for state in range(V)) % PRIME
        for relation in range(P)] for r1 in range(P)] for r0 in range(P)]
        for character in range(V)]
    r0_stage = [[[[
        sum(state_stage[character][r0][r1][relation] * roots[frequency][r0]
            for r0 in range(P)) % PRIME
        for relation in range(P)] for r1 in range(P)] for frequency in range(P)]
        for character in range(V)]
    r1_stage = [[[[
        sum(r0_stage[character][f0][r1][relation] * roots[frequency][r1]
            for r1 in range(P)) % PRIME
        for relation in range(P)] for frequency in range(P)] for f0 in range(P)]
        for character in range(V)]
    return tuple(
        tuple(
            tuple(
                tuple(
                    sum(r1_stage[character][f0][f1][relation]
                        * roots[ft][relation] for relation in range(P)) % PRIME
                    for ft in range(P)
                )
                for f1 in range(P)
            )
            for f0 in range(P)
        )
        for character in range(V)
    )


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
    state_matrix = tuple(
        tuple(tensor[state][r0][r1][relation]
              for r0 in range(P) for r1 in range(P))
        for state in range(V)
    )
    conditional = tuple(
        C.rank_mod(tuple(
            tuple(tensor[state][r0][r1][relation] for state in range(V))
            for r1 in range(P)
        ))
        for r0 in range(P)
    )
    return (
        C.rank_mod(state_matrix),
        tuple(sum(value != 0 for value in row) for row in state_matrix),
        conditional,
        digest_json(state_matrix),
    )


def brute_inversion_probe(chunks, tensor, zeta: int):
    """Direct 2,197-term inversion at held-out entries.

    This deliberately does not use the compressed alpha/beta/tau route in
    ``inverse_compressed``.  Worker probe rows already contain the beta twist.
    """
    normalizer = pow(P**3, -1, PRIME)
    records = []
    for probe_index, probe in enumerate(INVERSION_PROBES):
        state, r0, r1 = probe
        values = []
        for relation in PROBE_RELATIONS:
            total = 0
            for chunk in chunks:
                alpha = chunk[0]
                recorded_probe, rows = chunk[6][probe_index]
                require(recorded_probe == probe,
                        ("brute inversion probe order", recorded_probe, probe))
                index = 0
                for _beta in range(P):
                    for tau in range(P):
                        phase = pow(zeta, -(alpha + tau * relation) % P, PRIME)
                        total = (total + rows[index] * phase) % PRIME
                        index += 1
                require(index == P**2, ("brute inversion row size", probe, index))
            value = total * normalizer % PRIME
            require(value == tensor[state][r0][r1][relation],
                    ("compressed/direct inversion mismatch", probe, relation,
                     value, tensor[state][r0][r1][relation]))
            values.append(value)
        records.append((probe, tuple(values)))
    return tuple(records)


def main() -> None:
    _profiles, boundaries, _cells, source_record = two_digit_context()
    source_digest = digest_json((source_record, boundaries))
    if EXPECTED_SOURCE_SHA256 != "TO_BE_PINNED":
        require(source_digest == EXPECTED_SOURCE_SHA256,
                ("two-digit source drift", source_digest, EXPECTED_SOURCE_SHA256))

    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
            "two-digit worker order")

    marginal_gamma = tuple(
        row for chunk in chunks for row in chunk[2]
    )
    marginal_gamma_digest = digest_json(marginal_gamma)
    require(marginal_gamma_digest == D.EXPECTED_DIGESTS[0],
            ("two-digit gamma misses one-digit parent", marginal_gamma_digest))
    alpha_gamma_digests = tuple(chunk[3] for chunk in chunks)
    if EXPECTED_ALPHA_GAMMA_SHA256 != "TO_BE_PINNED":
        require(alpha_gamma_digests == EXPECTED_ALPHA_GAMMA_SHA256,
                ("two-digit alpha gamma drift", alpha_gamma_digests))

    zeta = C.context()["zeta"]
    tensor, tau_core = inverse_compressed(chunks, zeta)
    inversion_probe_record = brute_inversion_probe(chunks, tensor, zeta)
    first_digit = marginal_r1(tensor)
    square = marginal_r0(first_digit)
    require(digest_json(first_digit) == D.EXPECTED_DIGESTS[1],
            "two-digit tensor misses one-digit parent")
    require(digest_json(square) == B.EXPECTED_DIGESTS[2],
            "two-digit marginal chain misses Boolean square")

    flat = flat_second_digit_lift(first_digit)
    require(marginal_r1(flat) == first_digit,
            "flat second-digit lift misses one-digit parent")
    tensors = (tensor, flat)
    interactions = tuple(four_way_interaction(value) for value in tensors)
    spectra = tuple(fourier_tensor(value, zeta) for value in tensors)
    interaction_spectra = tuple(fourier_tensor(value, zeta) for value in interactions)
    ranks = tuple(axis_ranks(value) for value in tensors)
    interaction_ranks = tuple(axis_ranks(value) for value in interactions)
    conditionals = tuple(conditional_ranks(value) for value in tensors)
    require(all(rank > 1 for rank in conditionals[0][0]),
            ("scalar state/relation-blind digit-factor hostile did not separate",
             conditionals[0]))
    require(conditionals[1][0] == (1,) * P and conditionals[1][1] == (0,) * P,
            ("flat recurrence hostile", conditionals[1]))
    censuses = tuple(support_census(value) for value in spectra)
    interaction_censuses = tuple(support_census(value) for value in interaction_spectra)

    relation = 6
    fixed_records = tuple(fixed_relation_record(value, relation) for value in tensors)
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
    rank_record = (ranks, interaction_ranks)
    census_record = (censuses, interaction_censuses)
    if EXPECTED_DIGESTS != "TO_BE_PINNED":
        require(digests == EXPECTED_DIGESTS, ("two-digit digests", digests))
    if EXPECTED_RANKS != "TO_BE_PINNED":
        require(rank_record == EXPECTED_RANKS, ("two-digit ranks", rank_record))
    if EXPECTED_CONDITIONAL_RANKS != "TO_BE_PINNED":
        require(conditionals == EXPECTED_CONDITIONAL_RANKS,
                ("two-digit conditional ranks", conditionals))
    if EXPECTED_CENSUSES != "TO_BE_PINNED":
        require(census_record == EXPECTED_CENSUSES,
                ("two-digit spectral censuses", census_record))
    if EXPECTED_FIXED != "TO_BE_PINNED":
        require(fixed_records == EXPECTED_FIXED,
                ("two-digit fixed relation", fixed_records))

    regular_controls = dict(
        item for chunk in chunks for item in chunk[5]
    )
    direct_record = []
    for alpha, beta, tau in DIRECT_CONTROLS:
        direct_bank, direct_counts = integrate_two_digits(alpha, beta, tau)
        phase = pow(zeta, beta, PRIME)
        direct_row = phase_bank(direct_bank[0], phase)
        require(direct_row == regular_controls[(alpha, beta, tau)],
                ("two-digit literal endpoint guard", alpha, beta, tau))
        parent_direct, parent_diagonal, _parent_counts = D.integrate_branches(
            alpha, beta, tau
        )
        parent_row = tuple(
            tuple(phase * value % PRIME for value in r0_row)
            for r0_row in parent_direct[0]
        )
        require(marginal_gamma_r1(direct_row) == parent_row,
                ("literal control misses one-digit parent", alpha, beta, tau))
        require(all(value == 0 for state_row in parent_diagonal[0]
                    for value in state_row),
                ("literal parent diagonal", alpha, beta, tau))
        direct_record.append(((alpha, beta, tau), direct_counts))

    scalar_counts = tuple(
        sum(chunk[4][0][index] for chunk in chunks) for index in range(5)
    )
    state_counts = tuple(
        sum(chunk[4][1][state] for chunk in chunks) for state in range(V)
    )
    update_count = sum(chunk[4][2] for chunk in chunks)
    key_counts = tuple(
        sum(chunk[4][3][tau] for chunk in chunks) for tau in range(P)
    )
    work_record = (scalar_counts, state_counts, update_count, key_counts,
                   tuple(direct_record))

    record = (
        CURRENT_SHA256,
        D.EXPECTED_SEMANTIC_SHA256,
        B.EXPECTED_SEMANTIC_SHA256,
        C.JOINT_PRIME,
        C.JOINT_ROOT,
        zeta,
        source_record,
        boundaries,
        source_digest,
        alpha_gamma_digests,
        marginal_gamma_digest,
        work_record,
        inversion_probe_record,
        digests,
        rank_record,
        conditionals,
        census_record,
        relation,
        fixed_records,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("two-digit semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 Boolean-square x two current inverse digits x relation probe ==")
    print(f"parent=(sha256={CURRENT_SHA256},semantic={D.EXPECTED_SEMANTIC_SHA256})")
    print("digits=(a=r0+13*r1+169*c);factorization=13^3 fold -> r1 window -> r0 window -> root window")
    print("half_open_nested_window_equals_single_mod169_window: PASS")
    print(f"source=(fold_pieces={source_record[3]},r1_window_pieces={source_record[4]},address_pieces={source_record[5]},root_pieces={source_record[6]},joint_boundaries={source_record[7]},support_per_active_root={source_record[8]},digest={source_digest})")
    print("reflection=(y->1-y,u->12-u,r0->12-r0,r1->12-r1,state->state_XOR_2): PASS")
    print(f"field=(prime={C.JOINT_PRIME},root={C.JOINT_ROOT},zeta13={zeta})")
    print(f"sparse_work=(parent_counts={scalar_counts},state_counts={state_counts},coefficient_updates={update_count},postaggregate_keys_by_tau={key_counts})")
    print(f"literal_controls={DIRECT_CONTROLS}: PASS;same_root=0")
    print(f"direct_2197_term_inversion_probes={inversion_probe_record}: PASS")
    print("marginal_chain=(sum_r1->audited_r_owner,sum_r0->audited_Boolean_square): PASS")
    print(f"axis_ranks_(state,r0,r1,relation,address)=(actual,flat)={ranks}")
    print(f"four_way_ANOVA_axis_ranks={interaction_ranks}")
    print(f"conditional_r1_by_(state,relation)_ranks_(raw,contrast)=(actual,flat)={conditionals}")
    print("scalar_state_relation_blind_factor_hostile=rank<=1 at each fixed r0;actual rank>1 for all r0: PASS")
    print("spectral_census_order=(total,16 bins by mask state/r0/r1/relation from 0000 to 1111)")
    print(f"spectral_censuses_(actual,flat)={censuses}")
    print(f"four_way_ANOVA_spectral_censuses={interaction_censuses}")
    print(f"fixed_relation=(1,0,{relation});records_(state_rank,state_support,conditional_ranks,digest)={fixed_records}")
    print(f"alpha_gamma_digests={alpha_gamma_digests}")
    print(f"digests=(tensor,tau_core,r0_parent,square,flat,ANOVA,spectra,ANOVA_spectra)={digests}")
    print(f"semantic_sha256={semantic}")
    print("status=FINITE-EXACT two-digit current-ancestry projective/conditional-rank candidate on one owner base")
    print("scope=two current digits only;not complete address,not U_clock chronology,not arrival ancestry,not physical current,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
