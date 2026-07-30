#!/usr/bin/env python3
"""Exact all-label closure of the non-tail projected k=2 descent rows.

The global scalar atlases on 1800 <= z1 <= 1835 leave ten exact-suffix
rows: one at 1824, one at 1812, six at 1810, one at 1805, and one at 1800.
Every body here has ruler L=11760.  This verifier applies one data-driven
pipeline to all ten rows:

  residue rays -> denominator multisets -> crude all-divisor capacities
  -> common K5 Hunter status -> exact scalar-slack label packets
  -> lossless projected residual.

The first, second, and last rows die already in the quotient/status stages.
The remaining status states have a positive scalar cutoff on every residue
ray, so their literal universe of later distinct nonaligned labels is finite
without imposing a horizon.  Every scalar-eligible literal packet has
projected residual at least 25/91, contradicting containment in the union of
the two aligned danger combs.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[1]
UNIFORM = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
)
PROJECTED = (
    ROOT
    / "04-computation"
    / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
)
BAND_1810 = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1810_1835_thm2941.out"
)
BAND_1800 = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1800_1809_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_exact_descent_1800_1824_closure_thm2941.out"
)
EXPECTED_UNIFORM_SHA256 = (
    "dfa4788297b8c31fc9b5dce1afadf29d20b267cb4159fa95dadb9346b1980b36"
)
EXPECTED_PROJECTED_SHA256 = (
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
)
EXPECTED_BAND_1810_SHA256 = (
    "d197eb6179a3f7c7da08d4389fde988c0bd1fbc5db8cfaf8e30435ace3c7d87f"
)
EXPECTED_BAND_1800_SHA256 = (
    "a652db146760a151572ca2ff8f093cf297cf3a6322df441e530d9da3fb24ba0a"
)
EXPECTED_PROFILE_SHA256 = (
    "17c26ab5bb9d2c4bf76ea2b82c29a6387c7ab78add6839100a9b68da337ca542"
)
EXPECTED_SEMANTIC_SHA256 = (
    "fb22f354d62253bfe8d46ac50bb7d8541bcf50ae9277bd9bba542b8d911dcbc9"
)

QUANTIFIER = "distinct later nonaligned labels"

HOSTILE = (1, 4, 8, 10, 12, 14)
SLICES = (
    (1824, HOSTILE, "STATUS"),
    (1812, HOSTILE, "STATUS"),
    (1810, (1, 2, 8, 10, 12, 14), "PROJECT"),
    (1810, HOSTILE, "PROJECT"),
    (1810, (1, 6, 8, 10, 12, 14), "PROJECT"),
    (1810, (2, 4, 8, 10, 12, 14), "PROJECT"),
    (1810, (2, 6, 8, 10, 12, 14), "PROJECT"),
    (1810, (4, 6, 8, 10, 12, 14), "PROJECT"),
    (1805, HOSTILE, "PROJECT"),
    (1800, HOSTILE, "STATUS"),
)
HIGH_CASE_KEYS = (
    (1810, (2, 6, 8, 9, 10, 14)),
    (1810, (2, 8, 9, 10, 12, 14)),
    (1807, (1, 2, 10, 12, 13, 14)),
    (1807, (1, 4, 10, 12, 13, 14)),
    (1807, (1, 8, 10, 12, 13, 14)),
    (1800, (1, 2, 10, 12, 13, 14)),
)

# (scalar states, crude kills, common-status kills, status survivors,
#  scalar-eligible literal packets).  The z=1805 row is discovered below
#  and pinned after the first independent run.
EXPECTED_COUNTS = {
    (1824, HOSTILE): (38, 15, 23, 0, 0),
    (1812, HOSTILE): (11, 4, 7, 0, 0),
    (1810, (1, 2, 8, 10, 12, 14)): (15, 3, 6, 6, 8),
    (1810, HOSTILE): (304, 86, 216, 2, 2),
    (1810, (1, 6, 8, 10, 12, 14)): (47, 13, 30, 4, 6),
    (1810, (2, 4, 8, 10, 12, 14)): (94, 11, 23, 60, 105),
    (1810, (2, 6, 8, 10, 12, 14)): (27, 0, 12, 15, 20),
    (1810, (4, 6, 8, 10, 12, 14)): (6, 0, 2, 4, 4),
    (1805, HOSTILE): (2, 0, 0, 2, 2),
    (1800, HOSTILE): (14, 5, 9, 0, 0),
}

ALIGNED_TWO_UNION_CAP = F(25, 91)
SUFFIX_SLOTS = 4


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(UNIFORM) == EXPECTED_UNIFORM_SHA256, "uniform engine changed")
require(file_sha256(PROJECTED) == EXPECTED_PROJECTED_SHA256, "projected engine changed")
require(file_sha256(BAND_1810) == EXPECTED_BAND_1810_SHA256, "1810 atlas changed")
require(file_sha256(BAND_1800) == EXPECTED_BAND_1800_SHA256, "1800 atlas changed")
USPEC = spec_from_file_location("k2_exact_descent_uniform", UNIFORM)
require(USPEC is not None and USPEC.loader is not None, "cannot load uniform engine")
U = module_from_spec(USPEC)
USPEC.loader.exec_module(U)
PSPEC = spec_from_file_location("k2_exact_descent_projected", PROJECTED)
require(PSPEC is not None and PSPEC.loader is not None, "cannot load projected engine")
P = module_from_spec(PSPEC)
PSPEC.loader.exec_module(P)
require(P.A.RULER == U.suffix.A.RULER, "uniform/projected master rulers disagree")


def spanning_trees(vertex_count):
    edges = tuple(combinations(range(vertex_count), 2))
    result = []
    for chosen in combinations(edges, vertex_count - 1):
        seen = {0}
        changed = True
        while changed:
            changed = False
            for left, right in chosen:
                if left in seen and right not in seen:
                    seen.add(right)
                    changed = True
                elif right in seen and left not in seen:
                    seen.add(left)
                    changed = True
        if len(seen) == vertex_count:
            result.append(tuple(sorted(chosen)))
    require(
        len(result) == vertex_count ** (vertex_count - 2),
        ("Cayley tree count changed", vertex_count, len(result)),
    )
    return tuple(sorted(result))


def prufer_trees(vertex_count):
    """Independent K5 tree atlas from Prüfer words."""
    trees = set()
    for word in product(range(vertex_count), repeat=vertex_count - 2):
        degree = [1] * vertex_count
        for vertex in word:
            degree[vertex] += 1
        edges = []
        for vertex in word:
            leaf = next(index for index, value in enumerate(degree) if value == 1)
            edges.append(tuple(sorted((leaf, vertex))))
            degree[leaf] -= 1
            degree[vertex] -= 1
        leaves = tuple(index for index, value in enumerate(degree) if value == 1)
        require(len(leaves) == 2, ("Pruefer terminal leaves", word, leaves))
        edges.append(tuple(sorted(leaves)))
        trees.add(tuple(sorted(edges)))
    return tuple(sorted(trees))


TREES5 = spanning_trees(5)
require(TREES5 == prufer_trees(5), "edge and Pruefer K5 tree atlases disagree")


@lru_cache(maxsize=None)
def hunter_cap5(M, ds, es):
    sizes = tuple((M // d) * e for d, e in zip(ds, es))
    overlaps = {
        (i, j): U.pair_lower(M, ds[i], es[i], ds[j], es[j])
        for i, j in combinations(range(5), 2)
    }
    invoice = max(sum(overlaps[edge] for edge in tree) for tree in TREES5)
    return min(M, sum(sizes) - invoice)


def hunter_status_data5(D, ds, q):
    M = D // q
    inner_ds = []
    lows = []
    marginals = []
    for d in ds:
        common = gcd(d, q)
        low, remainder = divmod((d + 6) // 7, common)
        inner_ds.append(d // common)
        lows.append(low)
        marginals.append((q // common) * remainder)
    capacities = tuple(
        hunter_cap5(
            M,
            tuple(inner_ds),
            tuple(lows[index] + ((pattern >> index) & 1) for index in range(5)),
        )
        for pattern in range(32)
    )
    return tuple(marginals), capacities


@lru_cache(maxsize=None)
def common_status_feasible5(q, marginals, capacities, histogram):
    """Float-discover a status verdict; certify every kill over Q."""
    tail_rows = []
    tail_rhs = []
    thresholds = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(count for load, count in histogram if load >= threshold)
        good = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(good):
            continue
        thresholds.append(threshold)
        if not any(good):
            return False, ((threshold,), (F(1),), (F(0),) * 6)
        tail_rows.append(good)
        tail_rhs.append(demand)
    equality_rows = [
        (1,) * 32,
        *[
            tuple((pattern >> index) & 1 for pattern in range(32))
            for index in range(5)
        ],
    ]
    equality_rhs = (q, *marginals)
    if not tail_rows:
        return True, None
    primal = linprog(
        np.zeros(32),
        A_ub=-np.asarray(tail_rows, dtype=float),
        b_ub=-np.asarray(tail_rhs, dtype=float),
        A_eq=np.asarray(equality_rows, dtype=float),
        b_eq=np.asarray(equality_rhs, dtype=float),
        bounds=(0, None),
        method="highs",
    )
    if primal.success:
        return True, None
    require(primal.status == 2, ("unexpected primal status", primal.status))
    tail_count = len(tail_rows)
    dual_rows = []
    dual_rhs = []
    for pattern in range(32):
        dual_rows.append(
            tuple(tail_rows[row][pattern] for row in range(tail_count))
            + tuple(-equality_rows[row][pattern] for row in range(6))
        )
        dual_rhs.append(0)
    dual_rows.append(tuple(-value for value in tail_rhs) + tuple(equality_rhs))
    dual_rhs.append(-1)
    dual = linprog(
        np.zeros(tail_count + 6),
        A_ub=np.asarray(dual_rows, dtype=float),
        b_ub=np.asarray(dual_rhs, dtype=float),
        bounds=[(0, None)] * tail_count + [(None, None)] * 6,
        method="highs",
    )
    if not dual.success:
        return True, None
    alpha = tuple(F(float(value)).limit_denominator(1_000_000) for value in dual.x[:tail_count])
    z = tuple(F(float(value)).limit_denominator(1_000_000) for value in dual.x[tail_count:])
    slacks = tuple(
        sum(z[row] * equality_rows[row][pattern] for row in range(6))
        - sum(alpha[row] * tail_rows[row][pattern] for row in range(tail_count))
        for pattern in range(32)
    )
    contradiction = (
        sum(z[row] * equality_rhs[row] for row in range(6))
        - sum(alpha[row] * tail_rhs[row] for row in range(tail_count))
    )
    if (
        all(value >= 0 for value in alpha)
        and all(value >= 0 for value in slacks)
        and contradiction < 0
    ):
        return False, (tuple(thresholds), alpha, z)
    return True, None


def projected_safe_lower_bound(cells, L, labels, stop_at_cap=True):
    common_danger = ((F(0), F(1)),)
    cells_used = 0
    for cell in cells:
        local_union = P.merge_fraction(
            [
                interval
                for label in labels
                for interval in P.phase_danger(cell, label, L)
            ]
        )
        common_danger = P.intersect_fraction(common_danger, local_union)
        cells_used += 1
        if stop_at_cap and P.interval_mass(common_danger) <= 1 - ALIGNED_TWO_UNION_CAP:
            break
    return 1 - P.interval_mass(common_danger), cells_used, common_danger


def direct_projected_mass(body, L, labels):
    carrier = tuple(
        (F(left, P.A.RULER), F(right, P.A.RULER))
        for left, right in P.A.carrier_for(body)
    )
    removed = P.merge_fraction(
        [interval for label in labels for interval in P.danger_fraction(label)]
    )
    residual = P.subtract_fraction(carrier, removed)
    projected = []
    for left, right in residual:
        scaled_left = L * left
        scaled_right = L * right
        for integer in range(P.floor_fraction(scaled_left), P.ceil_fraction(scaled_right)):
            piece_left = max(scaled_left, F(integer)) - integer
            piece_right = min(scaled_right, F(integer + 1)) - integer
            if piece_left < piece_right:
                projected.append((piece_left, piece_right))
    return P.interval_mass(P.merge_fraction(projected))


def atlas_partition():
    keys = []
    high_keys = []
    for path in (BAND_1810, BAND_1800):
        for line in path.read_text().splitlines():
            if not line.startswith("SURVIVOR;"):
                continue
            fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
            key = (int(fields["z1"]), tuple(map(int, fields["E"].split(","))))
            keys.append(key)
            if "HIGH-TAIL" in line:
                high_keys.append(key)
    return tuple(sorted(keys)), tuple(sorted(high_keys))


ATLAS_KEYS, ATLAS_HIGH_KEYS = atlas_partition()
require(
    ATLAS_KEYS
    == tuple(
        sorted(tuple((first, body) for first, body, _mode in SLICES) + HIGH_CASE_KEYS)
    ),
    "exact/high cases do not partition the scalar-atlas survivors",
)
require(
    ATLAS_HIGH_KEYS == tuple(sorted(HIGH_CASE_KEYS)),
    "HIGH-TAIL atlas classification disagrees with the high-ray cases",
)


def ray_and_status(first, body, carrier, h, lower, L):
    first_delta = U.delta(carrier, h, first)
    first_d = L // gcd(L, first)
    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require(
            (residue + L) * U.delta(carrier, h, residue + L) == amplitude,
            (first, body, "ray recurrence", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(amplitudes[L - residue] == -amplitudes[residue] for residue in range(1, L)),
        (first, body, "ray antipode"),
    )
    require(L * U.delta(carrier, h, L) == 0, (first, body, "aligned ray is nonzero"))
    require(
        signs[1] == signs[-1] and sum(signs.values()) == L - 1,
        (first, body, "ray sign census", signs),
    )
    ray_digest = sha256(repr(tuple(amplitudes)).encode()).hexdigest()
    divisors = tuple(d for d in U.support.divisors(L) if d > 1)

    @lru_cache(maxsize=None)
    def top_values(d, multiplicity):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + L - 1) // L) * L
            candidates.extend(
                (amplitude / (label + offset * L), label + offset * L, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        require(len(candidates) >= multiplicity, (first, body, "missing ray", d))
        return tuple(candidates[:multiplicity])

    trials = comb(len(divisors) + SUFFIX_SLOTS - 1, SUFFIX_SLOTS)
    scalar = []
    for tail in combinations_with_replacement(divisors, SUFFIX_SLOTS):
        upper = first_delta
        labels = []
        for d, multiplicity in Counter(tail).items():
            chosen = top_values(d, multiplicity)
            upper += sum((value for value, _label, _residue in chosen), F())
            labels.extend(label for _value, label, _residue in chosen)
        if upper >= lower:
            scalar.append(
                (
                    tuple(sorted((first_d, *tail))),
                    upper,
                    (first, *sorted(labels)),
                )
            )
    require(
        len({ds for ds, _upper, _labels in scalar}) == len(scalar),
        (first, body, "duplicate denominator state"),
    )

    actual_L, ranges = U.support.safe_cell_ranges(body)
    require(actual_L == L, (first, body, "safe-cell ruler"))
    arcs_cache = {}
    crude_kills = []
    crude_survivors = []
    for ds, upper, labels in scalar:
        D = lcm(*ds)
        if D not in arcs_cache:
            arcs_cache[D] = U.fibre.projected_support_arcs(D, ranges)
        arcs = arcs_cache[D]
        witness = None
        for q in U.support.divisors(D):
            histogram = U.fibre.residue_load_histogram(arcs, q)
            target = max(load for load, count in histogram if count)
            capacity = sum(U.fibre_cap(D, d, q) for d in ds)
            if target > capacity:
                witness = (q, D // q, target, capacity)
                break
        row = (ds, upper, labels, witness)
        (crude_survivors if witness is None else crude_kills).append(row)

    status_kills = []
    status_survivors = []
    for ds, upper, labels, _crude_witness in crude_survivors:
        D = lcm(*ds)
        arcs = arcs_cache[D]
        witness = None
        for M in U.support.divisors(D):
            q = D // M
            marginals, capacities = hunter_status_data5(D, ds, q)
            histogram = U.fibre.residue_load_histogram(arcs, q)
            feasible, certificate = common_status_feasible5(
                q, marginals, capacities, histogram
            )
            if not feasible:
                require(certificate is not None, (first, body, ds, "missing witness"))
                witness = (
                    q,
                    M,
                    marginals,
                    tuple(sorted(set(capacities))),
                    histogram,
                    certificate,
                )
                break
        row = (ds, upper, labels, witness)
        (status_survivors if witness is None else status_kills).append(row)

    stages = (scalar, crude_kills, status_kills, status_survivors)
    stage_digests = tuple(sha256(repr(tuple(rows)).encode()).hexdigest() for rows in stages)
    return (
        tuple(amplitudes),
        ray_digest,
        len(divisors),
        trials,
        first_delta,
        first_d,
        tuple(scalar),
        tuple(crude_kills),
        tuple(status_kills),
        tuple(status_survivors),
        stage_digests,
    )


def projected_packets(first, body, carrier, h, lower, L, amplitudes, first_delta, states):
    cells = P.body_cells(P.A.carrier_for(body), L)

    @lru_cache(maxsize=None)
    def ray_heads(d, multiplicity):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + L - 1) // L) * L
            candidates.extend(
                (amplitude / (label + offset * L), label + offset * L, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        return tuple(candidates[:multiplicity])

    @lru_cache(maxsize=None)
    def eligible_group(d, multiplicity, slack):
        heads = ray_heads(d, multiplicity)
        top_sum = sum((row[0] for row in heads), F())
        threshold = heads[-1][0] - slack
        require(
            threshold > 0,
            (first, body, "scalar slack did not make ray finite", d, multiplicity, slack),
        )
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude <= 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + L - 1) // L) * L
            while amplitude / label >= threshold:
                candidates.append((amplitude / label, label, residue))
                label += L
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        choices = []
        for chosen in combinations(candidates, multiplicity):
            labels = tuple(sorted(row[1] for row in chosen))
            if len(set(labels)) != multiplicity:
                continue
            value = sum((row[0] for row in chosen), F())
            deficit = top_sum - value
            if deficit <= slack:
                choices.append((deficit, labels, value))
        choices.sort()
        require(
            choices and choices[0][0] == 0,
            (first, body, "maximizing group choice missing", d, multiplicity),
        )
        return top_sum, threshold, tuple(candidates), tuple(choices)

    packet_count = 0
    killed_count = 0
    survivor_packets = []
    minimum_margin = None
    maximum_cells_used = 0
    minimum_row = None
    state_audit = []
    first_d = L // gcd(L, first)
    for state_index, (ds, upper, maximizing_labels, witness) in enumerate(states, 1):
        require(witness is None, (first, body, "status survivor has witness", ds))
        slack = upper - lower
        tail = list(ds)
        tail.remove(first_d)
        groups = []
        for d, multiplicity in sorted(Counter(tail).items()):
            top_sum, threshold, candidates, choices = eligible_group(d, multiplicity, slack)
            groups.append((d, multiplicity, top_sum, threshold, len(candidates), choices))
        expected_upper = first_delta + sum((group[2] for group in groups), F())
        require(
            expected_upper == upper,
            (first, body, "state upper changed", state_index, expected_upper, upper),
        )
        state_packets = 0
        state_kills = 0
        state_minimum = None
        for selection in product(*(group[5] for group in groups)):
            total_deficit = sum((choice[0] for choice in selection), F())
            if total_deficit > slack:
                continue
            labels = (first, *sorted(label for choice in selection for label in choice[1]))
            require(
                len(labels) == 5 and len(set(labels)) == 5,
                (first, body, "literal labels not distinct", labels),
            )
            actual_upper = first_delta + sum((choice[2] for choice in selection), F())
            require(
                actual_upper == upper - total_deficit and actual_upper >= lower,
                (first, body, "literal scalar arithmetic", labels),
            )
            projected_lower, cells_used, _common = projected_safe_lower_bound(
                cells, L, labels
            )
            margin = projected_lower - ALIGNED_TWO_UNION_CAP
            packet_count += 1
            state_packets += 1
            maximum_cells_used = max(maximum_cells_used, cells_used)
            minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
            state_minimum = margin if state_minimum is None else min(state_minimum, margin)
            if minimum_row is None or (margin, labels) < (minimum_row[0], minimum_row[1]):
                minimum_row = (margin, labels, projected_lower, cells_used)
            if margin > 0:
                killed_count += 1
                state_kills += 1
            else:
                survivor_packets.append((state_index, ds, labels, actual_upper, projected_lower))
        require(state_packets > 0, (first, body, "status state has no literal packet", ds))
        state_audit.append(
            (
                state_index,
                ds,
                upper,
                maximizing_labels,
                slack,
                tuple(
                    (d, multiplicity, top_sum, threshold, candidate_count, len(choices))
                    for d, multiplicity, top_sum, threshold, candidate_count, choices in groups
                ),
                state_packets,
                state_kills,
                state_minimum,
            )
        )
    require(
        packet_count == killed_count and not survivor_packets,
        (first, body, "projected packet survived", survivor_packets),
    )
    direct_mass = None
    if minimum_row is not None:
        _margin, control_labels, prefix_mass, _prefix_cells = minimum_row
        full_cell_mass, full_cells, _common = projected_safe_lower_bound(
            cells, L, control_labels, stop_at_cap=False
        )
        direct_mass = direct_projected_mass(body, L, control_labels)
        require(
            full_cells == len(cells)
            and direct_mass == full_cell_mass
            and direct_mass >= prefix_mass,
            (first, body, "independent projected control"),
        )
    state_digest = sha256(repr(tuple(state_audit)).encode()).hexdigest()
    return (
        len(cells),
        packet_count,
        killed_count,
        minimum_margin,
        maximum_cells_used,
        minimum_row,
        direct_mass,
        state_digest,
        tuple(state_audit),
    )


def profile(item):
    first, body, mode = item
    carrier = U.suffix.A.carrier_for(body)
    require(P.A.carrier_for(body) == carrier, (first, body, "carrier engines disagree"))
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    require(L == 11760, (first, body, "unexpected exact-body ruler", L))
    require(P.A.RULER % L == 0, (first, body, "body ruler left master ruler"))
    high_floor = max(15, (U.suffix.PROJECTED_RATIOS[2] * L).numerator // (U.suffix.PROJECTED_RATIOS[2] * L).denominator + 1)
    require(first >= high_floor, (first, body, "unexpected high-wall obligation"))
    (
        amplitudes,
        ray_digest,
        divisor_count,
        trials,
        first_delta,
        first_d,
        scalar,
        crude_kills,
        status_kills,
        states,
        stage_digests,
    ) = ray_and_status(first, body, carrier, h, lower, L)
    counts_prefix = (len(scalar), len(crude_kills), len(status_kills), len(states))
    if mode == "STATUS":
        require(not states, (first, body, "status-only row survived", len(states)))
        projected = (0, 0, 0, None, 0, None, None, sha256(b"()").hexdigest(), ())
    else:
        projected = projected_packets(
            first, body, carrier, h, lower, L, amplitudes, first_delta, states
        )
        require(projected[4] <= 2, (first, body, "projected prefix exceeded two cells"))
    counts = (*counts_prefix, projected[1])
    expected = EXPECTED_COUNTS[(first, body)]
    if expected is not None:
        require(counts == expected, (first, body, "stage counts changed", counts))
    return (
        first,
        body,
        mode,
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        ray_digest,
        divisor_count,
        trials,
        counts,
        stage_digests,
        *projected[:-1],
    )


def render(profiles):
    require(tuple((row[0], row[1], row[2]) for row in profiles) == SLICES, "slice universe changed")
    require(all(row[12][3] == 0 or row[15] == row[16] for row in profiles), "closure count mismatch")
    total_scalar = sum(row[12][0] for row in profiles)
    total_crude = sum(row[12][1] for row in profiles)
    total_status = sum(row[12][2] for row in profiles)
    total_residual = sum(row[12][3] for row in profiles)
    total_packets = sum(row[15] for row in profiles)
    total_kills = sum(row[16] for row in profiles)
    positive_margins = tuple(row[17] for row in profiles if row[17] is not None)
    require(total_packets == total_kills, "a projected packet survived globally")
    require(all(margin > 0 for margin in positive_margins), "nonpositive projected margin")
    require(
        (
            total_scalar,
            total_crude,
            total_status,
            total_residual,
            total_packets,
            total_kills,
        )
        == (558, 137, 328, 93, 147, 147),
        "global descent ledger changed",
    )
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        SLICES,
        HIGH_CASE_KEYS,
        QUANTIFIER,
        ALIGNED_TWO_UNION_CAP,
        total_scalar,
        total_crude,
        total_status,
        total_residual,
        total_packets,
        total_kills,
        positive_margins,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exact descent closure 1800..1824",
        f"uniform_engine_sha256={file_sha256(UNIFORM)}",
        f"projected_engine_sha256={file_sha256(PROJECTED)}",
        f"scalar_band_1810_sha256={file_sha256(BAND_1810)}",
        f"scalar_band_1800_sha256={file_sha256(BAND_1800)}",
        f"scope=ten exact-suffix rows;all {QUANTIFIER};no finite label horizon",
        (
            "pipeline=residue rays;denominator multisets;all-divisor crude;"
            "common-K5 status;scalar-slack literal packets;projected residual"
        ),
        "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);two-aligned cap=25/91",
        (
            f"global_counts=scalar:{total_scalar};crude:{total_crude};status:{total_status};"
            f"status_survivors:{total_residual};literal_packets:{total_packets};"
            f"projected_kills:{total_kills};survivors:0"
        ),
    ]
    for row in profiles:
        (
            first,
            body,
            mode,
            h,
            components,
            L,
            lower,
            first_delta,
            first_d,
            ray_digest,
            divisor_count,
            trials,
            counts,
            stage_digests,
            cell_count,
            packet_count,
            killed_count,
            minimum_margin,
            maximum_cells,
            minimum_row,
            direct_mass,
            state_digest,
        ) = row
        lines.append(
            f"SLICE;z1={first};E={','.join(map(str, body))};mode={mode};"
            f"h={ftext(h)};r={components};L={L};lower={ftext(lower)};"
            f"delta1={ftext(first_delta)};first_d={first_d};"
            f"ray_sha256={ray_digest};denominators={divisor_count};trials={trials};"
            f"counts={counts};stage_sha256={stage_digests};body_cells={cell_count};"
            f"packets={packet_count};kills={killed_count};"
            f"min_margin={ftext(minimum_margin)};max_prefix_cells={maximum_cells};"
            f"direct_control_mass={ftext(direct_mass)};state_sha256={state_digest};"
            f"minimum_row={minimum_row}"
        )
    lines.extend(
        (
            "conclusion=all ten exact-suffix scalar exceptions on 1800..1835 are empty uniformly",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(len(SLICES), mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    if args.workers == 1:
        profiles = [profile(item) for item in SLICES]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, SLICES))
    order = {(first, body): index for index, (first, body, _mode) in enumerate(SLICES)}
    profiles.sort(key=lambda row: order[(row[0], row[1])])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
