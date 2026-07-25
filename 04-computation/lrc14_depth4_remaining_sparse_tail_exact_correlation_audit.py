#!/usr/bin/env python3
"""Exact sparse-tail audit for scalar profiles 114, 124, and 134.

This is a neutral computation artifact.  It imports the exact modular
convolution engine from

    lrc14_depth4_hinge_gram_exact_correlation_audit.py

and contains no floating-point arithmetic.

On the signed unit group G_5=(Z/13^5 Z)^*/{+/-1}, write every exponent as

    e = r + j |G_4|,       0 <= j < 13.

For a terminal exponent q, its root-fibre row is

    h_q(r)=sum_j C(r+j|G_4|) D_5(q+r+j|G_4|) in {0,1,2}.

The full singleton capacity F_q is sum_r h_q(r).  A blocker of depth
lambda has a periodic active row D_(5-lambda)(r+a), so all singleton
owner capacities at a fixed q are one short exact cyclic correlation of
h_q with that row.

The hinge thresholds leave sparse tails:

* 114: theta=17495, 1,835 (terminal,owner,value) triples.  Only 177,699
  of 86,889,153 unordered owner pairs share any high terminal label; a
  nine-sheet overlap floor settles all but two of those pairs, and exact
  meet/Hellinger rows settle the two exceptions.
* 124: theta=17544.  No terminal label is high for both a depth-one and
  a depth-two owner.
* 134: theta=17524.  No depth-three singleton row has a high entry.

All convolutions are exact finite-field NTTs.  Every endpoint inequality
is strict on the power-of-thirteen carrier.  No Python ``assert`` is load
bearing, so checks remain active under ``python -O``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import struct
from pathlib import Path

import numpy as np
from numba import njit, prange


P = 13
N = P**5
FULL_ORDER = 6 * P**4
BASE_ORDER = 6 * P**3
GENERATOR = 2

THETA_114 = 17495
THETA_124 = 17544
THETA_134 = 17524
TERMINAL_COUNT = 5


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ENGINE_PATH = (
    Path(__file__).resolve().parent
    / "lrc14_depth4_hinge_gram_exact_correlation_audit.py"
)
ENGINE_SPEC = importlib.util.spec_from_file_location(
    "depth4_exact_engine", ENGINE_PATH
)
require(ENGINE_SPEC is not None and ENGINE_SPEC.loader is not None,
        "could not load the exact correlation engine")
ENGINE = importlib.util.module_from_spec(ENGINE_SPEC)
ENGINE_SPEC.loader.exec_module(ENGINE)

MODULUS = ENGINE.NTT_MODULUS
PRIMITIVE_ROOT = ENGINE.NTT_ROOT


def norm_numerators(values, modulus):
    residues = values % modulus
    return np.minimum(residues, modulus - residues)


def canonical_label_from_exponent(exponent, modulus):
    value = pow(GENERATOR, int(exponent), modulus)
    return min(value, modulus - value)


def next_power_of_two(value):
    answer = 1
    while answer < value:
        answer *= 2
    return answer


def bit_reversal(size):
    bits = size.bit_length() - 1
    answer = np.empty(size, dtype=np.int64)
    for value in range(size):
        source = value
        target = 0
        for _ in range(bits):
            target = (target << 1) | (source & 1)
            source >>= 1
        answer[value] = target
    return answer


@njit
def small_ntt_from_natural(
    source,
    target,
    reversal,
    stage_roots,
    inverse_scale,
):
    """Natural-order exact NTT for any supported power-of-two length."""
    size = reversal.size
    for index in range(size):
        target[index] = source[reversal[index]]

    length = 2
    stage = 0
    while length <= size:
        half = length // 2
        root_step = stage_roots[stage]
        for block in range(0, size, length):
            root_power = 1
            for offset in range(half):
                left = target[block + offset]
                right = (
                    target[block + offset + half] * root_power
                ) % MODULUS
                total = left + right
                if total >= MODULUS:
                    total -= MODULUS
                difference = left - right
                if difference < 0:
                    difference += MODULUS
                target[block + offset] = total
                target[block + offset + half] = difference
                root_power = (root_power * root_step) % MODULUS
        stage += 1
        length *= 2

    if inverse_scale != 0:
        for index in range(size):
            target[index] = (
                target[index] * inverse_scale
            ) % MODULUS


@njit(parallel=True)
def small_exact_cyclic_correlations(
    rows,
    period,
    padded_length,
    danger_transform,
    reversal,
    forward_roots,
    inverse_roots,
    inverse_scale,
):
    """Exact correlations sum_r rows[r] danger[r+a]."""
    count = rows.shape[0]
    answers = np.empty((count, period), dtype=np.int32)
    for row_index in prange(count):
        padded = np.zeros(padded_length, dtype=np.int64)
        transformed = np.empty(padded_length, dtype=np.int64)
        padded[0] = rows[row_index, 0]
        for exponent in range(1, period):
            padded[exponent] = rows[row_index, period - exponent]

        small_ntt_from_natural(
            padded, transformed, reversal, forward_roots, 0
        )
        for index in range(padded_length):
            transformed[index] = (
                transformed[index] * danger_transform[index]
            ) % MODULUS
        small_ntt_from_natural(
            transformed,
            padded,
            reversal,
            inverse_roots,
            inverse_scale,
        )
        for shift in range(period):
            coefficient = padded[shift]
            if shift + period < padded_length:
                coefficient += padded[shift + period]
            answers[row_index, shift] = coefficient
    return answers


def reduced_danger(depth):
    reduced_power = 5 - depth
    modulus = P**reduced_power
    period = 6 * P ** (reduced_power - 1)
    residues = np.empty(period, dtype=np.int64)
    residues[0] = 1
    for exponent in range(1, period):
        residues[exponent] = (
            residues[exponent - 1] * GENERATOR
        ) % modulus
    danger = 14 * norm_numerators(residues, modulus) < modulus
    return modulus, period, danger


def small_ntt_data(period, danger):
    padded_length = next_power_of_two(2 * period - 1)
    require(
        (MODULUS - 1) % padded_length == 0,
        "small NTT length does not divide the field order",
    )
    reversal = bit_reversal(padded_length)
    stages = padded_length.bit_length() - 1
    forward_roots = np.asarray(
        [
            pow(
                PRIMITIVE_ROOT,
                (MODULUS - 1) // (1 << stage),
                MODULUS,
            )
            for stage in range(1, stages + 1)
        ],
        dtype=np.int64,
    )
    inverse_roots = np.asarray(
        [
            pow(int(root), MODULUS - 2, MODULUS)
            for root in forward_roots
        ],
        dtype=np.int64,
    )
    inverse_scale = pow(padded_length, MODULUS - 2, MODULUS)

    padded = np.zeros(padded_length, dtype=np.int64)
    padded[:period] = danger.astype(np.int64)
    danger_transform = np.empty(padded_length, dtype=np.int64)
    small_ntt_from_natural(
        padded, danger_transform, reversal, forward_roots, 0
    )
    return (
        padded_length,
        danger_transform,
        reversal,
        forward_roots,
        inverse_roots,
        inverse_scale,
    )


def correlate_small(rows, period, data):
    (
        padded_length,
        danger_transform,
        reversal,
        forward_roots,
        inverse_roots,
        inverse_scale,
    ) = data
    return small_exact_cyclic_correlations(
        rows,
        period,
        padded_length,
        danger_transform,
        reversal,
        forward_roots,
        inverse_roots,
        inverse_scale,
    )


def full_carrier():
    representatives = np.empty(FULL_ORDER, dtype=np.int64)
    representatives[0] = 1
    for exponent in range(1, FULL_ORDER):
        representatives[exponent] = (
            representatives[exponent - 1] * GENERATOR
        ) % N
    norms = norm_numerators(representatives, N)
    guard = 7 * norms > N
    danger = 14 * norms < N
    require(int(guard.sum()) == 122405, "bad signed guard size")
    require(N % 7 != 0 and N % 14 != 0, "torsion endpoint risk")
    return representatives, guard, danger


def exact_full_capacity(guard, full_ntt_data):
    (
        reversal,
        forward_roots,
        inverse_roots,
        danger_transform,
    ) = full_ntt_data
    return ENGINE.exact_cyclic_correlations(
        guard[None, :],
        danger_transform,
        reversal,
        forward_roots,
        inverse_roots,
    )[0]


def root_rows_for_candidates(candidate_exponents, guard, danger):
    """Build h_q(r) exactly for a batch of terminal exponents."""
    doubled_danger = np.concatenate((danger, danger))
    guard_fibres = guard.reshape(13, BASE_ORDER)
    shifted = np.asarray(
        [
            doubled_danger[
                int(exponent):int(exponent) + FULL_ORDER
            ]
            for exponent in candidate_exponents
        ],
        dtype=np.bool_,
    )
    counts = (
        shifted.reshape(-1, 13, BASE_ORDER)
        & guard_fibres[None, :, :]
    ).sum(axis=1, dtype=np.int16)
    require(int(counts.max()) <= 2, "root row exceeded two sheets")
    return counts


def digest_little_endian(array, dtype):
    return hashlib.sha256(
        np.asarray(array).astype(dtype, copy=False).tobytes(order="C")
    ).hexdigest()


def owner_loss_rows(guard_fibre_counts, danger_by_depth, ntt_by_depth):
    """Weighted one-owner losses L_a at each depth."""
    losses = {}
    for depth in (1, 2, 3):
        period = danger_by_depth[depth][1]
        aggregated = guard_fibre_counts.reshape(
            -1, period
        ).sum(axis=0, dtype=np.int16)
        losses[depth] = correlate_small(
            aggregated[None, :],
            period,
            ntt_by_depth[depth],
        )[0]
    return losses


def exact_depth_one_overlap(danger_one, ntt_one):
    """J(d)=|A_0 intersection A_d|."""
    return correlate_small(
        danger_one.astype(np.int16)[None, :],
        danger_one.size,
        ntt_one,
    )[0]


def collect_sparse_tails(
    full_capacity,
    candidates,
    guard,
    danger,
    danger_by_depth,
    ntt_by_depth,
):
    """All exact singleton tails needed by the three profiles."""
    requested = {
        (1, THETA_114): [],
        (1, THETA_124): [],
        (1, THETA_134): [],
        (2, THETA_124): [],
        (3, THETA_134): [],
    }
    maxima = {key: {} for key in requested}
    batch_size = 26
    for start in range(0, len(candidates), batch_size):
        batch = candidates[start:start + batch_size]
        root_rows = root_rows_for_candidates(batch, guard, danger)
        for depth in (1, 2, 3):
            period = danger_by_depth[depth][1]
            aggregated = root_rows.reshape(
                len(batch), -1, period
            ).sum(axis=1, dtype=np.int16)
            removals = correlate_small(
                aggregated,
                period,
                ntt_by_depth[depth],
            )

            thresholds = {
                1: (THETA_114, THETA_124, THETA_134),
                2: (THETA_124,),
                3: (THETA_134,),
            }[depth]
            for local_index, terminal_exponent in enumerate(batch):
                singleton = (
                    int(full_capacity[terminal_exponent])
                    - removals[local_index].astype(np.int64)
                )
                for threshold in thresholds:
                    maxima[(depth, threshold)][
                        int(terminal_exponent)
                    ] = int(singleton.max())
                    active = np.flatnonzero(singleton > threshold)
                    bucket = requested[(depth, threshold)]
                    for owner_exponent in active:
                        bucket.append(
                            (
                                int(terminal_exponent),
                                int(owner_exponent),
                                int(singleton[owner_exponent]),
                            )
                        )
    for key in requested:
        requested[key].sort()
    return requested, maxima


def tail_by_terminal(triples):
    answer = {}
    for terminal, owner, value in triples:
        answer.setdefault(terminal, []).append((owner, value))
    return answer


def floor_114_global_minimum(weight, losses_one, overlap):
    """Minimize W0-L_a-L_b+9J(b-a) over all unordered pairs."""
    period = len(losses_one)
    minimum = None
    record = None
    minimizers = 0
    indices = np.arange(period, dtype=np.int64)
    for difference in range(period):
        right = (indices + difference) % period
        values = (
            weight
            - losses_one.astype(np.int64)
            - losses_one[right].astype(np.int64)
            + 9 * int(overlap[difference])
        )
        local_minimum = int(values.min())
        active = np.flatnonzero(values == local_minimum)
        for left in active:
            right_exponent = int(right[left])
            first = min(int(left), right_exponent)
            second = max(int(left), right_exponent)
            candidate = (first, second, local_minimum, difference)
            if minimum is None or local_minimum < minimum:
                minimum = local_minimum
                record = candidate
                minimizers = 1
            elif local_minimum == minimum:
                minimizers += 1
                if candidate < record:
                    record = candidate
    # Every unordered pair appears twice except the diagonal.  We only need
    # the value and first canonical witness; count unique minimizers directly.
    unique_minimizers = set()
    for difference in range(period):
        right = (indices + difference) % period
        values = (
            weight
            - losses_one.astype(np.int64)
            - losses_one[right].astype(np.int64)
            + 9 * int(overlap[difference])
        )
        for left in np.flatnonzero(values == minimum):
            pair = tuple(sorted((int(left), int(right[left]))))
            unique_minimizers.add(pair)
    return minimum, record, tuple(sorted(unique_minimizers))


def build_114_pair_records(
    triples,
    weight,
    losses_one,
    overlap,
    guard_fibre_counts,
    danger_one,
):
    """Sparse shared-tail records and the complete deterministic digest."""
    by_terminal = tail_by_terminal(triples)
    packed_danger = np.packbits(
        danger_one[None, :], axis=1, bitorder="little"
    )[0]
    base_bits = int.from_bytes(
        packed_danger.tobytes(), "little"
    )
    packed_ten = np.packbits(
        (guard_fibre_counts == 10)[None, :],
        axis=1,
        bitorder="little",
    )[0]
    ten_sheet_bits = int.from_bytes(
        packed_ten.tobytes(), "little"
    )
    full_bits = (1 << BASE_ORDER) - 1
    active_bits = []
    for shift in range(BASE_ORDER):
        if shift == 0:
            rotated = base_bits
        else:
            rotated = (
                (base_bits >> shift)
                | ((base_bits & ((1 << shift) - 1))
                   << (BASE_ORDER - shift))
            )
        active_bits.append(rotated & full_bits)

    high_values = {}
    occurrence_count = 0
    for terminal in sorted(by_terminal):
        entries = by_terminal[terminal]
        for left_index, (left_owner, left_value) in enumerate(entries):
            for right_index in range(left_index, len(entries)):
                right_owner, right_value = entries[right_index]
                pair = (
                    min(left_owner, right_owner),
                    max(left_owner, right_owner),
                )
                high_values.setdefault(pair, []).append(
                    min(left_value, right_value)
                )
                occurrence_count += 1

    digest = hashlib.sha256()
    failures = []
    minimum_margin = None
    minimum_record = None
    for pair in sorted(high_values):
        left, right = pair
        difference = (right - left) % BASE_ORDER
        floor_weight = (
            weight
            - int(losses_one[left])
            - int(losses_one[right])
            + 9 * int(overlap[difference])
        )
        correction = (
            active_bits[left]
            & active_bits[right]
            & ten_sheet_bits
        ).bit_count()
        exact_weight = floor_weight + correction
        high = sorted(high_values[pair], reverse=True)[:TERMINAL_COUNT]
        high_count = len(high)
        bound = (
            sum(high)
            + (TERMINAL_COUNT - high_count) * THETA_114
        )
        margin = exact_weight - bound
        padded = tuple(high + [0] * (TERMINAL_COUNT - high_count))
        digest.update(
            struct.pack(
                "<11i",
                left,
                right,
                exact_weight,
                bound,
                margin,
                high_count,
                *padded,
            )
        )
        record = (
            left,
            right,
            exact_weight,
            bound,
            margin,
            high_count,
            padded,
        )
        if margin <= 0:
            failures.append(record)
        if minimum_margin is None or margin < minimum_margin:
            minimum_margin = margin
            minimum_record = record
    return (
        len(high_values),
        occurrence_count,
        minimum_margin,
        minimum_record,
        tuple(failures),
        digest.hexdigest(),
    )


def direct_owner_singleton_row(
    guard,
    danger_one,
    owner_exponent,
    full_ntt_data,
):
    base_classes = np.arange(FULL_ORDER, dtype=np.int64) % BASE_ORDER
    residual = guard & ~danger_one[
        (base_classes + owner_exponent) % BASE_ORDER
    ]
    (
        reversal,
        forward_roots,
        inverse_roots,
        danger_transform,
    ) = full_ntt_data
    capacities = ENGINE.exact_cyclic_correlations(
        residual[None, :],
        danger_transform,
        reversal,
        forward_roots,
        inverse_roots,
    )[0]
    return residual, capacities


def exact_pair_row(
    guard,
    left_depth,
    left_exponent,
    right_depth,
    right_exponent,
    danger_by_depth,
    full_ntt_data,
):
    exponents = np.arange(FULL_ORDER, dtype=np.int64)
    left_period = danger_by_depth[left_depth][1]
    right_period = danger_by_depth[right_depth][1]
    left_danger = danger_by_depth[left_depth][2]
    right_danger = danger_by_depth[right_depth][2]
    residual = (
        guard
        & ~left_danger[
            (exponents % left_period + left_exponent) % left_period
        ]
        & ~right_danger[
            (exponents % right_period + right_exponent) % right_period
        ]
    )
    (
        reversal,
        forward_roots,
        inverse_roots,
        danger_transform,
    ) = full_ntt_data
    capacities = ENGINE.exact_cyclic_correlations(
        residual[None, :],
        danger_transform,
        reversal,
        forward_roots,
        inverse_roots,
    )[0]
    return residual, capacities


def top_values(array, count=5):
    indices = np.argpartition(array, -count)[-count:]
    return tuple(
        sorted((int(array[index]) for index in indices), reverse=True)
    )


def verify_direct_top_capacities(
    name,
    representatives,
    residual,
    capacities,
    count=5,
):
    """Replay the selected NTT winners from the residue definition."""
    indices = np.argpartition(capacities, -count)[-count:]
    for index in indices:
        products = (
            int(representatives[index]) * representatives
        ) % N
        dangerous = 14 * norm_numerators(products, N) < N
        direct = int(np.sum(residual & dangerous))
        require(
            direct == int(capacities[index]),
            f"{name} direct terminal capacity changed",
        )


def exception_114_control(
    representatives,
    guard,
    danger_one,
    danger_by_depth,
    full_ntt_data,
    left,
    right,
    theta,
):
    _, left_row = direct_owner_singleton_row(
        guard, danger_one, left, full_ntt_data
    )
    _, right_row = direct_owner_singleton_row(
        guard, danger_one, right, full_ntt_data
    )
    meet = np.minimum(left_row, right_row)
    meet_top = top_values(meet)

    left_tail = np.maximum(left_row.astype(np.int64) - theta, 0)
    right_tail = np.maximum(right_row.astype(np.int64) - theta, 0)
    common = np.flatnonzero((left_tail > 0) & (right_tail > 0))
    hellinger = sum(
        int(ENGINE.ceiling_square_root(
            int(left_tail[index] * right_tail[index])
        ))
        for index in common
    )
    bound = TERMINAL_COUNT * theta + hellinger

    residual, pair_capacities = exact_pair_row(
        guard,
        1,
        left,
        1,
        right,
        danger_by_depth,
        full_ntt_data,
    )
    pair_top = top_values(pair_capacities)
    verify_direct_top_capacities(
        "114 exception",
        representatives,
        residual,
        pair_capacities,
    )
    return (
        int(residual.sum()),
        meet_top,
        sum(meet_top),
        theta,
        len(common),
        hellinger,
        bound,
        int(residual.sum()) - bound,
        pair_top,
        sum(pair_top),
        int(residual.sum()) - sum(pair_top),
    )


def mixed_summary_digest(
    full_capacity,
    threshold,
    left_triples,
    right_triples,
    left_maxima,
    right_maxima,
):
    left = tail_by_terminal(left_triples)
    right = tail_by_terminal(right_triples)
    digest = hashlib.sha256()
    candidates = np.flatnonzero(full_capacity > threshold)
    cross_products = 0
    maximum_left = 0
    maximum_right = 0
    for terminal in candidates:
        left_entries = left.get(int(terminal), ())
        right_entries = right.get(int(terminal), ())
        max_left = int(left_maxima[int(terminal)])
        max_right = int(right_maxima[int(terminal)])
        count_left = len(left_entries)
        count_right = len(right_entries)
        cross_products += count_left * count_right
        maximum_left = max(maximum_left, count_left)
        maximum_right = max(maximum_right, count_right)
        digest.update(
            struct.pack(
                "<6i",
                int(terminal),
                int(full_capacity[terminal]),
                max_left,
                max_right,
                count_left,
                count_right,
            )
        )
    return (
        len(candidates),
        sum(len(entries) for entries in left.values()),
        len(left),
        maximum_left,
        sum(len(entries) for entries in right.values()),
        len(right),
        maximum_right,
        cross_products,
        digest.hexdigest(),
    )


def periodic_active_bitsets(danger, period):
    """All shifted periodic blocker masks on the G4 phase carrier."""
    phase_indices = np.arange(BASE_ORDER, dtype=np.int64) % period
    answer = []
    for shift in range(period):
        row = danger[(phase_indices + shift) % period]
        packed = np.packbits(
            row[None, :], axis=1, bitorder="little"
        )[0]
        answer.append(int.from_bytes(packed.tobytes(), "little"))
    return tuple(answer)


def exact_mixed_weight_minimum(
    weight,
    left_losses,
    right_losses,
    left_bits,
    right_bits,
    guard_fibre_counts,
):
    """Complete exact typed-pair W table via nine/ten-sheet bitsets."""
    ten_packed = np.packbits(
        (guard_fibre_counts == 10)[None, :],
        axis=1,
        bitorder="little",
    )[0]
    ten_bits = int.from_bytes(ten_packed.tobytes(), "little")
    digest = hashlib.sha256()
    minimum = None
    record = None
    minimizers = 0
    pair_count = 0
    for left, left_mask in enumerate(left_bits):
        for right, right_mask in enumerate(right_bits):
            intersection = left_mask & right_mask
            overlap = (
                9 * intersection.bit_count()
                + (intersection & ten_bits).bit_count()
            )
            residual = (
                weight
                - int(left_losses[left])
                - int(right_losses[right])
                + overlap
            )
            digest.update(struct.pack("<3i", left, right, residual))
            pair_count += 1
            candidate = (left, right, residual, overlap)
            if minimum is None or residual < minimum:
                minimum = residual
                record = candidate
                minimizers = 1
            elif residual == minimum:
                minimizers += 1
                if candidate < record:
                    record = candidate
    return (
        pair_count,
        minimum,
        record,
        minimizers,
        digest.hexdigest(),
    )


def mixed_hostile_control(
    representatives,
    guard,
    danger_by_depth,
    full_ntt_data,
    left_depth,
    left_exponent,
    right_depth,
    right_exponent,
):
    left_danger = danger_by_depth[left_depth][2]
    right_danger = danger_by_depth[right_depth][2]
    # Full singleton rows are exact correlations on their residual signals.
    exponents = np.arange(FULL_ORDER, dtype=np.int64)
    left_period = danger_by_depth[left_depth][1]
    right_period = danger_by_depth[right_depth][1]
    (
        reversal,
        forward_roots,
        inverse_roots,
        danger_transform,
    ) = full_ntt_data
    left_residual = guard & ~left_danger[
        (exponents % left_period + left_exponent) % left_period
    ]
    right_residual = guard & ~right_danger[
        (exponents % right_period + right_exponent) % right_period
    ]
    singleton_rows = ENGINE.exact_cyclic_correlations(
        np.asarray((left_residual, right_residual), dtype=np.bool_),
        danger_transform,
        reversal,
        forward_roots,
        inverse_roots,
    )
    meet = np.minimum(singleton_rows[0], singleton_rows[1])
    meet_top = top_values(meet)

    pair_residual, pair_capacities = exact_pair_row(
        guard,
        left_depth,
        left_exponent,
        right_depth,
        right_exponent,
        danger_by_depth,
        full_ntt_data,
    )
    actual_top = top_values(pair_capacities)
    verify_direct_top_capacities(
        "mixed hostile",
        representatives,
        pair_residual,
        pair_capacities,
    )
    return (
        int(pair_residual.sum()),
        int(meet.max()),
        meet_top,
        sum(meet_top),
        int(pair_residual.sum()) - sum(meet_top),
        actual_top,
        sum(actual_top),
        int(pair_residual.sum()) - sum(actual_top),
    )


def run():
    representatives, guard, danger = full_carrier()
    full_ntt_data = ENGINE.build_ntt_data(danger)
    full_capacity = exact_full_capacity(guard, full_ntt_data)
    require(
        digest_little_endian(full_capacity, "<i4")
        == "2ef7b86380667c3da2ef22f4031ea38737f458facb6111260bf1d5b67956efb7",
        "full capacity digest changed",
    )

    guard_fibre_counts = guard.reshape(
        13, BASE_ORDER
    ).sum(axis=0, dtype=np.uint8)
    require(
        set(map(int, np.unique(guard_fibre_counts))) == {9, 10},
        "guard fibres were not nine/ten sheets",
    )
    require(
        digest_little_endian(guard_fibre_counts, "<u1")
        == "702d0205453ac1948544e942420f9ea38640b568c4dd010304f71b4ba8ab5dba",
        "guard-fibre digest changed",
    )

    danger_by_depth = {
        depth: reduced_danger(depth) for depth in (1, 2, 3)
    }
    ntt_by_depth = {
        depth: small_ntt_data(
            danger_by_depth[depth][1],
            danger_by_depth[depth][2],
        )
        for depth in (1, 2, 3)
    }
    losses = owner_loss_rows(
        guard_fibre_counts, danger_by_depth, ntt_by_depth
    )
    overlap_one = exact_depth_one_overlap(
        danger_by_depth[1][2], ntt_by_depth[1]
    )
    require(
        digest_little_endian(losses[1], "<i4")
        == "c5db84957db4c9667e3c8b1fbd04e39cb02bcf59e777658428d8c387c191c111",
        "depth-one owner-loss digest changed",
    )
    require(
        digest_little_endian(overlap_one, "<i4")
        == "382288fd86c362d6ada87eb147b64bca677b7a1b05469fe83874fb2cc57f4be5",
        "depth-one overlap digest changed",
    )

    candidates = np.flatnonzero(full_capacity > THETA_114)
    require(len(candidates) == 6294, "bad broad candidate count")
    tails, maxima = collect_sparse_tails(
        full_capacity,
        candidates,
        guard,
        danger,
        danger_by_depth,
        ntt_by_depth,
    )

    triples_114 = tails[(1, THETA_114)]
    triples_digest = hashlib.sha256()
    for triple in triples_114:
        triples_digest.update(struct.pack("<III", *triple))
    require(
        len(triples_114) == 1835
        and triples_digest.hexdigest()
        == "0db7ecf9fc595119e9d3feb9d84da32625ca541f0fe58dfbe9ec132bb5ee1d65",
        "114 tail triples changed",
    )

    floor_minimum = floor_114_global_minimum(
        int(guard.sum()), losses[1], overlap_one
    )
    require(
        floor_minimum[0] == 87478
        and floor_minimum[1][:3] == (4055, 13181, 87478)
        and (4055, 13181) in floor_minimum[2],
        "114 global floor minimum changed",
    )
    pair_records = build_114_pair_records(
        triples_114,
        int(guard.sum()),
        losses[1],
        overlap_one,
        guard_fibre_counts,
        danger_by_depth[1][2],
    )
    require(
        pair_records[0] == 177699
        and pair_records[1] == 178781
        and pair_records[5]
        == "d3e70461c1ba01d3407d45bbcb1f84da6ee285a058cbc997c5cdaa0a56be4476",
        "114 sparse pair table changed",
    )
    failure_pairs = tuple(
        record[:2] for record in pair_records[4]
    )
    require(
        failure_pairs == ((0, 1998), (1997, 13181)),
        "114 padded exceptions changed",
    )

    exception_one = exception_114_control(
        representatives,
        guard,
        danger_by_depth[1][2],
        danger_by_depth,
        full_ntt_data,
        0,
        1998,
        17000,
    )
    exception_two = exception_114_control(
        representatives,
        guard,
        danger_by_depth[1][2],
        danger_by_depth,
        full_ntt_data,
        1997,
        13181,
        17000,
    )
    require(
        exception_one[:8]
        == (
            87708,
            (18045, 17260, 16947, 16946, 16852),
            86050,
            17000,
            2,
            1742,
            86742,
            966,
        ),
        "first 114 exception meet/Hellinger row changed",
    )
    require(
        exception_one[9:] == (78180, 9528),
        "first 114 actual pair row changed",
    )
    require(
        exception_two[:8]
        == (
            87865,
            (18673, 17260, 17072, 16946, 16852),
            86803,
            17000,
            3,
            2371,
            87371,
            494,
        ),
        "second 114 exception meet/Hellinger row changed",
    )
    require(
        exception_two[9:] == (78572, 9293),
        "second 114 actual pair row changed",
    )

    summary_124 = mixed_summary_digest(
        full_capacity,
        THETA_124,
        tails[(1, THETA_124)],
        tails[(2, THETA_124)],
        maxima[(1, THETA_124)],
        maxima[(2, THETA_124)],
    )
    summary_134 = mixed_summary_digest(
        full_capacity,
        THETA_134,
        tails[(1, THETA_134)],
        tails[(3, THETA_134)],
        maxima[(1, THETA_134)],
        maxima[(3, THETA_134)],
    )
    require(
        summary_124
        == (
            819,
            250,
            129,
            43,
            3,
            3,
            1,
            0,
            "5065137672e85a192361941e1c2cdf8880d2bea8d4b4a467a90de1a605929bc2",
        ),
        "124 sparse summary changed",
    )
    require(
        summary_134
        == (
            1273,
            357,
            183,
            69,
            0,
            0,
            0,
            0,
            "066bf47aba202c38ed94821f1d97e3e0e5b92e767ba4ade0fc331042ebde846c",
        ),
        "134 sparse summary changed",
    )

    active_one = periodic_active_bitsets(
        danger_by_depth[1][2], danger_by_depth[1][1]
    )
    active_two = periodic_active_bitsets(
        danger_by_depth[2][2], danger_by_depth[2][1]
    )
    active_three = periodic_active_bitsets(
        danger_by_depth[3][2], danger_by_depth[3][1]
    )
    exact_weights_124 = exact_mixed_weight_minimum(
        int(guard.sum()),
        losses[1],
        losses[2],
        active_one,
        active_two,
        guard_fibre_counts,
    )
    exact_weights_134 = exact_mixed_weight_minimum(
        int(guard.sum()),
        losses[1],
        losses[3],
        active_one,
        active_three,
        guard_fibre_counts,
    )
    require(
        exact_weights_124
        == (
            13366548,
            87725,
            (0, 0, 87725, 1440),
            2,
            "dac30421571468ee979ecab38dea5b5fb1fcacefdbfd4f049c0e2dc33a3b343e",
        ),
        "124 exact weight minimum changed",
    )
    require(
        exact_weights_134
        == (
            1028196,
            87623,
            (0, 3, 87623, 2890),
            13,
            "45934a3a144c6f3f29fd58a1a3992aed2b98f9106b930ba3fd84b2e74a158817",
        ),
        "134 exact weight minimum changed",
    )
    require(
        exact_weights_124[1] - 5 * THETA_124 == 5,
        "124 uniform margin changed",
    )
    require(
        exact_weights_134[1] - 5 * THETA_134 == 3,
        "134 uniform margin changed",
    )

    hostile_124 = mixed_hostile_control(
        representatives,
        guard,
        danger_by_depth,
        full_ntt_data,
        1,
        0,
        2,
        0,
    )
    hostile_134 = mixed_hostile_control(
        representatives,
        guard,
        danger_by_depth,
        full_ntt_data,
        1,
        0,
        3,
        3,
    )
    require(
        (
            hostile_124[0],
            hostile_124[1],
            hostile_124[3],
            hostile_124[4],
            hostile_124[6],
            hostile_124[7],
        )
        == (87725, 17520, 85364, 2361, 76389, 11336),
        "124 hostile control changed",
    )
    require(
        (
            hostile_134[0],
            hostile_134[1],
            hostile_134[3],
            hostile_134[4],
            hostile_134[6],
            hostile_134[7],
        )
        == (87623, 17263, 84611, 3012, 76262, 11361),
        "134 hostile control changed",
    )

    pair_count_114 = BASE_ORDER * (BASE_ORDER + 1) // 2
    pair_count_124 = BASE_ORDER * danger_by_depth[2][1]
    pair_count_134 = BASE_ORDER * danger_by_depth[3][1]
    require(pair_count_114 == 86889153, "bad 114 pair count")
    require(pair_count_124 == 13366548, "bad 124 pair count")
    require(pair_count_134 == 1028196, "bad 134 pair count")

    labels = {
        exponent: canonical_label_from_exponent(exponent, P**4)
        for exponent in (0, 1997, 1998, 4055, 13181)
    }
    print("EXACT REMAINING DEPTH-FOUR SPARSE-TAIL AUDIT")
    print("carrier N / signed G5 / signed G4 =",
          (N, FULL_ORDER, BASE_ORDER))
    print("guard signed mass / root histogram =",
          (int(guard.sum()),
           tuple((value, int(np.sum(guard_fibre_counts == value)))
                 for value in (9, 10))))
    print("arithmetic = exact finite-field/integer; floating point = none")
    print("F / h / L1 / J digests =",
          (
              digest_little_endian(full_capacity, "<i4"),
              digest_little_endian(guard_fibre_counts, "<u1"),
              digest_little_endian(losses[1], "<i4"),
              digest_little_endian(overlap_one, "<i4"),
          ))
    print("114 threshold / candidates / triples / triple digest =",
          (THETA_114, len(candidates), len(triples_114),
           triples_digest.hexdigest()))
    print("114 global floor minimum =", floor_minimum)
    print("114 sparse pair audit =", pair_records)
    print("114 exception exponent->label map =", labels)
    print("114 exception (1,12) =", exception_one)
    print("114 exception (6,14280) =", exception_two)
    print("124 threshold / summary / exact weight table =",
          (THETA_124, summary_124, exact_weights_124))
    print("124 hostile control =", hostile_124)
    print("134 threshold / summary / exact weight table =",
          (THETA_134, summary_134, exact_weights_134))
    print("134 hostile control =", hostile_134)
    print("pair universes =", (pair_count_114, pair_count_124,
                                pair_count_134))
    print("VERDICT = profiles 114, 124, and 134 have positive deficits")


if __name__ == "__main__":
    run()
