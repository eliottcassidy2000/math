#!/usr/bin/env python3
"""Close both projected k=2 atlas rows at ``z1=1732`` exactly.

The proof is reconstructed directly from the current canonical THM-2941
uniform-ray/status and projected-cell engines.  It deliberately does not
import a prior k=2 descent wrapper.  The only earlier closure used for the
splice is the hash-pinned ``z1=1736`` theorem transcript.

For each of the two atlas bodies this verifier:

* checks the complete period ray, including recurrence, antipodes, zero rays,
  and balanced positive/negative directions;
* maximizes every denominator multiset over all later labels, so there is no
  label horizon;
* applies exact crude capacities and rationally replayed common-K5 Farkas
  witnesses;
* expands every residual scalar state through its positive exact slack; and
* compares the projected-cell calculation with direct full-carrier
  subtraction on every literal packet.

The first body is exhausted by ``33+155`` capacity/status witnesses.  The
second leaves twenty scalar states, each with one literal packet, and all
twenty packets have projected residual strictly above ``25/91``.
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
UNIFORM = ROOT / "04-computation/lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
PROJECTED = ROOT / "04-computation/lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1680_1742_thm2941.py"
ATLAS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
PRIOR_SOURCE = ROOT / "04-computation/lrc14_j7_k2_z1736_hybrid_closure_thm2941.py"
PRIOR_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_z1736_hybrid_closure_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1732_ray_status_projected_closure_thm2941.out"

DEPENDENCIES = (UNIFORM, PROJECTED, ATLAS_SOURCE, ATLAS_OUTPUT, PRIOR_SOURCE, PRIOR_OUTPUT)
EXPECTED_DEPENDENCY_HASHES = (
    "34ab29162ed33d90093e6d2bf781def36c420a1cd6596158b5d6579a3a8f3f46",
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662",
    "89016f939c961fa979ec5b30812981456df5bfb2af3066f1f1b38e5a83f1a412",
    "4a36611b26585964e185bbaa3d583be3f1c67a7b608cca785920266bc217a779",
    "5965eddea10a7e2c2d2b70d94052f6d69b3593c865fd72d3a8f1c8052cf1f96f",
    "548fd0d318e09ae4fa2da1844a2df5a50d8917c3fbfc48ab1b72ddf46f3d9678",
)

FIRST = 1732
SUFFIX_SLOTS = 4
SCALAR_ETA = F(1, 91)
ALIGNED_TWO_UNION_CAP = F(25, 91)
QUANTIFIER = "all distinct later nonaligned labels"
CASES = (
    (FIRST, (1, 4, 8, 10, 12, 14)),
    (FIRST, (2, 4, 8, 10, 12, 14)),
)
EXPECTED_HEIGHTS = (
    (1683, 1),
    (1694, 10),
    (1702, 3),
    (1708, 14),
    (1722, 11),
    (1724, 2),
    (1732, 2),
    (1736, 15),
)
# scalar states, crude kills, common-status kills, status survivors, packets
EXPECTED_COUNTS = {
    CASES[0]: (188, 33, 155, 0, 0),
    CASES[1]: (24, 1, 3, 20, 20),
}
EMPTY_SHA256 = "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"
EXPECTED_STAGE_DIGESTS = {
    CASES[0]: (
        "c4ce1fde263280dda337184ded196f4f31fae14ed1d0adb197478847d54d7467",
        "b26ca570b73cfa3c4f867e4895cc4078179e232dd5167d34f592e519be79cd85",
        "e9411b800f9dd248b14669be866015a9251c1c0f40ff3e0720ada710a326abae",
        EMPTY_SHA256,
    ),
    CASES[1]: (
        "14b8157df7fb8300a5a4d41bafa2bd3786d4995e8e9b2a90f4daa05e49d7d70d",
        "678cd0cd9a69539e98a0618175e0354641c600e6d0ec63135f589974c8f5cf94",
        "78293ffc40054b2b6a37384a9f1c53205b4009b812fc47285671b70664a6d571",
        "be0363f9e947ad7ab2bf768614433e933b2712e02bfb1276b2d42a2dd46ab03e",
    ),
}
EXPECTED_PACKET_AUDIT = (
    F(903353, 7131943),
    2,
    (F(903353, 7131943), (1732, 1736, 1750, 1810, 1836), F(31458, 78373), 2),
    F(15, 13423),
    F(1),
    "70462ba71eb90368300dc9c17e75a82eb70fca2ee04b1fa405f163bbc8018d02",
)
# Filled only after a full development replay, then frozen for fail-closed use.
EXPECTED_PROFILE_SHA256 = None
EXPECTED_SEMANTIC_SHA256 = None


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(
    tuple(file_sha256(path) for path in DEPENDENCIES) == EXPECTED_DEPENDENCY_HASHES,
    "dependency changed",
)
U = load_module("k2_z1732_uniform", UNIFORM)
P = load_module("k2_z1732_projected", PROJECTED)
require(P.A.RULER == U.suffix.A.RULER, "canonical rulers disagree")


def atlas_partition():
    keys = []
    high = []
    heights = Counter()
    for line in ATLAS_OUTPUT.read_text(encoding="utf-8").splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
        key = (int(fields["z1"]), tuple(map(int, fields["E"].split(","))))
        keys.append(key)
        heights[key[0]] += 1
        if "HIGH-TAIL" in line:
            high.append(key)
    return tuple(sorted(keys)), tuple(sorted(high)), tuple(sorted(heights.items()))


ATLAS_KEYS, ATLAS_HIGH, ATLAS_HEIGHTS = atlas_partition()
require(ATLAS_HEIGHTS == EXPECTED_HEIGHTS, "atlas height profile changed")
require(tuple(key for key in ATLAS_KEYS if key[0] == FIRST) == CASES, "z1732 rows changed")
require(not any(key[0] == FIRST for key in ATLAS_HIGH), "z1732 unexpectedly forced-high")
require(not any(1725 <= key[0] <= 1731 for key in ATLAS_KEYS), "splice interval occupied")
require(sum(key[0] == 1724 for key in ATLAS_KEYS) == 2, "next occupied height changed")
require(
    "consequence=projected k=2 first drift label z1<=1732" in PRIOR_OUTPUT.read_text(),
    "prior cap transcript changed",
)


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
    require(len(result) == vertex_count ** (vertex_count - 2), "Cayley count changed")
    return tuple(sorted(result))


def prufer_trees(vertex_count):
    result = set()
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
        require(len(leaves) == 2, "Pruefer terminal changed")
        edges.append(tuple(sorted(leaves)))
        result.add(tuple(sorted(edges)))
    return tuple(sorted(result))


TREES5 = spanning_trees(5)
require(TREES5 == prufer_trees(5), "independent K5 tree atlases disagree")


@lru_cache(maxsize=None)
def hunter_cap5(M, ds, es):
    sizes = tuple((M // d) * e for d, e in zip(ds, es))
    overlaps = {
        (left, right): U.pair_lower(M, ds[left], es[left], ds[right], es[right])
        for left, right in combinations(range(5), 2)
    }
    invoice = max(sum(overlaps[edge] for edge in tree) for tree in TREES5)
    return min(M, sum(sizes) - invoice)


def hunter_status_data5(D, ds, q):
    M = D // q
    inner = []
    lows = []
    marginals = []
    for d in ds:
        common = gcd(d, q)
        low, remainder = divmod((d + 6) // 7, common)
        inner.append(d // common)
        lows.append(low)
        marginals.append((q // common) * remainder)
    capacities = tuple(
        hunter_cap5(
            M,
            tuple(inner),
            tuple(lows[index] + ((pattern >> index) & 1) for index in range(5)),
        )
        for pattern in range(32)
    )
    return tuple(marginals), capacities


@lru_cache(maxsize=None)
def common_status_feasible5(q, marginals, capacities, histogram):
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
        *[tuple((pattern >> index) & 1 for pattern in range(32)) for index in range(5)],
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
    for pattern in range(32):
        dual_rows.append(
            tuple(tail_rows[row][pattern] for row in range(tail_count))
            + tuple(-equality_rows[row][pattern] for row in range(6))
        )
    dual_rows.append(tuple(-value for value in tail_rhs) + tuple(equality_rhs))
    dual = linprog(
        np.zeros(tail_count + 6),
        A_ub=np.asarray(dual_rows, dtype=float),
        b_ub=np.asarray((0,) * 32 + (-1,), dtype=float),
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
    contradiction = sum(z[row] * equality_rhs[row] for row in range(6)) - sum(
        alpha[row] * tail_rhs[row] for row in range(tail_count)
    )
    if all(value >= 0 for value in alpha) and all(value >= 0 for value in slacks) and contradiction < 0:
        return False, (tuple(thresholds), alpha, z)
    return True, None


def projected_safe_lower_bound(cells, L, labels, stop_at_cap=True):
    common_danger = ((F(0), F(1)),)
    cells_used = 0
    for cell in cells:
        local_union = P.merge_fraction(
            [interval for label in labels for interval in P.phase_danger(cell, label, L)]
        )
        common_danger = P.intersect_fraction(common_danger, local_union)
        cells_used += 1
        if stop_at_cap and P.interval_mass(common_danger) <= 1 - ALIGNED_TWO_UNION_CAP:
            break
    return 1 - P.interval_mass(common_danger), cells_used, common_danger


def projected_prefix(cells, L, labels, count):
    require(1 <= count <= len(cells), "bad projected prefix")
    common = ((F(0), F(1)),)
    for index, cell in enumerate(cells[:count], 1):
        local = P.merge_fraction(
            [interval for label in labels for interval in P.phase_danger(cell, label, L)]
        )
        common = P.intersect_fraction(common, local)
    return 1 - P.interval_mass(common), index, common


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


def ray_and_status(first, body, carrier, h, lower, L):
    first_delta = U.delta(carrier, h, first)
    first_d = L // gcd(L, first)
    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require(
            (residue + L) * U.delta(carrier, h, residue + L) == amplitude,
            (body, "ray recurrence"),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(amplitudes[L - residue] == -amplitudes[residue] for residue in range(1, L)),
        (body, "ray antipode"),
    )
    require(L * U.delta(carrier, h, L) == 0, (body, "zero ray"))
    require(signs[1] == signs[-1] and sum(signs.values()) == L - 1, (body, "ray signs"))
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
            # An antipodal negative ray only lowers the scalar excess, so it
            # cannot occur in a maximizing necessary-state representative.
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
        require(len(candidates) >= multiplicity, (body, "missing ray", d))
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
            scalar.append((tuple(sorted((first_d, *tail))), upper, (first, *sorted(labels))))
    require(len({row[0] for row in scalar}) == len(scalar), (body, "duplicate state"))

    actual_L, ranges = U.support.safe_cell_ranges(body)
    require(actual_L == L, (body, "safe-cell ruler"))
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
    for ds, upper, labels, _witness in crude_survivors:
        D = lcm(*ds)
        arcs = arcs_cache[D]
        witness = None
        for M in U.support.divisors(D):
            q = D // M
            marginals, capacities = hunter_status_data5(D, ds, q)
            histogram = U.fibre.residue_load_histogram(arcs, q)
            feasible, certificate = common_status_feasible5(q, marginals, capacities, histogram)
            if not feasible:
                require(certificate is not None, (body, ds, "missing exact Farkas witness"))
                witness = (q, M, marginals, tuple(sorted(set(capacities))), histogram, certificate)
                break
        row = (ds, upper, labels, witness)
        (status_survivors if witness is None else status_kills).append(row)
    stages = (scalar, crude_kills, status_kills, status_survivors)
    stage_digests = tuple(sha256(repr(tuple(rows)).encode()).hexdigest() for rows in stages)
    return (
        tuple(amplitudes),
        tuple(signs[index] for index in (1, -1, 0)),
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
        # This is the exact finiteness gate.  Any later member of any ray is
        # below threshold and cannot fit inside the state's total slack.
        require(threshold > 0, (body, "nonfinite ray cutoff", d, multiplicity))
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
        require(choices and choices[0][0] == 0, (body, "maximizer missing", d))
        return top_sum, threshold, tuple(candidates), tuple(choices)

    packet_count = 0
    killed_count = 0
    direct_controls = 0
    minimum_margin = None
    maximum_cells = 0
    minimum_row = None
    state_audit = []
    first_d = L // gcd(L, first)
    for state_index, (ds, upper, maximizing_labels, witness) in enumerate(states, 1):
        require(witness is None, (body, "status witness on survivor"))
        slack = upper - lower
        tail = list(ds)
        tail.remove(first_d)
        groups = []
        for d, multiplicity in sorted(Counter(tail).items()):
            top_sum, threshold, candidates, choices = eligible_group(d, multiplicity, slack)
            groups.append((d, multiplicity, top_sum, threshold, len(candidates), choices))
        require(first_delta + sum((group[2] for group in groups), F()) == upper, (body, "upper"))
        state_packets = 0
        state_kills = 0
        state_minimum = None
        for selection in product(*(group[5] for group in groups)):
            deficit = sum((choice[0] for choice in selection), F())
            if deficit > slack:
                continue
            labels = (first, *sorted(label for choice in selection for label in choice[1]))
            require(len(labels) == 5 and len(set(labels)) == 5, (body, "literal labels"))
            actual_upper = first_delta + sum((choice[2] for choice in selection), F())
            require(actual_upper == upper - deficit and actual_upper >= lower, (body, "scalar"))
            projected, cells_used, _common = projected_safe_lower_bound(cells, L, labels)
            full_mass, full_cells, _full_common = projected_safe_lower_bound(cells, L, labels, False)
            direct_mass = direct_projected_mass(body, L, labels)
            require(
                full_cells == len(cells) and direct_mass == full_mass >= projected,
                (body, labels, "direct projection disagreement"),
            )
            margin = projected - ALIGNED_TWO_UNION_CAP
            require(margin > 0 and direct_mass > ALIGNED_TWO_UNION_CAP, (body, labels, "packet survived"))
            packet_count += 1
            killed_count += 1
            direct_controls += 1
            state_packets += 1
            state_kills += 1
            maximum_cells = max(maximum_cells, cells_used)
            minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
            state_minimum = margin if state_minimum is None else min(state_minimum, margin)
            if minimum_row is None or (margin, labels) < (minimum_row[0], minimum_row[1]):
                minimum_row = (margin, labels, projected, cells_used)
        require(state_packets > 0 and state_packets == state_kills, (body, "empty or live state"))
        state_audit.append(
            (
                state_index,
                ds,
                upper,
                maximizing_labels,
                slack,
                tuple(
                    (d, multiplicity, top, threshold, candidate_count, len(choices))
                    for d, multiplicity, top, threshold, candidate_count, choices in groups
                ),
                state_packets,
                state_kills,
                state_minimum,
            )
        )
    require(packet_count == killed_count == direct_controls, (body, "packet ledger"))
    hostile_one_cell = None
    direct_minimum = None
    if minimum_row is not None:
        _margin, labels, _prefix_mass, _cells_used = minimum_row
        hostile_one_cell, one_cell_count, _one_cell_common = projected_prefix(cells, L, labels, 1)
        require(one_cell_count == 1 and hostile_one_cell <= ALIGNED_TWO_UNION_CAP, (body, "hostile control"))
        direct_minimum = direct_projected_mass(body, L, labels)
    return (
        len(cells),
        packet_count,
        killed_count,
        direct_controls,
        minimum_margin,
        maximum_cells,
        minimum_row,
        hostile_one_cell,
        direct_minimum,
        sha256(repr(tuple(state_audit)).encode()).hexdigest(),
        tuple(state_audit),
    )


def case_profile(case):
    first, body = case
    carrier = U.suffix.A.carrier_for(body)
    require(P.A.carrier_for(body) == carrier, (body, "carrier disagreement"))
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * SCALAR_ETA
    L = 14 * lcm(*body)
    require(U.suffix.A.RULER % L == 0, (body, "ruler"))
    data = ray_and_status(first, body, carrier, h, lower, L)
    (
        amplitudes,
        ray_signs,
        ray_digest,
        divisor_count,
        trials,
        first_delta,
        first_d,
        scalar,
        crude,
        status,
        states,
        stage_digests,
    ) = data
    if states:
        projected = projected_packets(
            first, body, carrier, h, lower, L, amplitudes, first_delta, states
        )
    else:
        projected = (
            len(P.body_cells(carrier, L)),
            0,
            0,
            0,
            None,
            0,
            None,
            None,
            None,
            EMPTY_SHA256,
            (),
        )
    counts = (len(scalar), len(crude), len(status), len(states), projected[1])
    require(counts == EXPECTED_COUNTS[case], (case, "counts", counts))
    require(stage_digests == EXPECTED_STAGE_DIGESTS[case], (case, "stage digest"))
    if case == CASES[0]:
        require(
            projected[1:] == (0, 0, 0, None, 0, None, None, None, EMPTY_SHA256, ()),
            "empty route",
        )
    else:
        packet_audit = (
            projected[4],
            projected[5],
            projected[6],
            projected[7],
            projected[8],
            projected[9],
        )
        require(packet_audit == EXPECTED_PACKET_AUDIT, (case, "packet audit", packet_audit))
        require(projected[1] == projected[2] == projected[3] == 20, (case, "packet controls"))
    return (
        first,
        body,
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        ray_signs,
        ray_digest,
        divisor_count,
        trials,
        counts,
        stage_digests,
        *projected[:-1],
    )


def render(profiles):
    require(len(profiles) == len(CASES), "profile count changed")
    totals = tuple(sum(row[12][index] for row in profiles) for index in range(5))
    require(totals == (212, 34, 158, 20, 20), ("global counts", totals))
    require(sum(row[16] for row in profiles) == 20, "projected kill total changed")
    require(sum(row[17] for row in profiles) == 20, "direct-control total changed")
    prior_ledger = 2_239_774
    final_ledger = prior_ledger - len(CASES)
    profile_digest = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_digest == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        EXPECTED_DEPENDENCY_HASHES,
        QUANTIFIER,
        CASES,
        EXPECTED_HEIGHTS,
        EXPECTED_COUNTS,
        EXPECTED_STAGE_DIGESTS,
        EXPECTED_PACKET_AUDIT,
        totals,
        prior_ledger,
        final_ledger,
        profile_digest,
    )
    semantic_digest = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 ray/status/projected closure at z1=1732",
        *(f"dependency_sha256={path.name}:{file_sha256(path)}" for path in DEPENDENCIES),
        f"scope=both atlas rows;{QUANTIFIER};no finite label horizon",
        f"atlas=height_counts:{EXPECTED_HEIGHTS};z1732_rows:2;empty:1725..1731;next:1724",
        "universe=first label 1732;four ordered-distinct later nonaligned labels;all denominator multisets and primitive ray directions",
        "negative_ray_control=exact recurrence plus antipodal sign balance;negative amplitudes cannot improve scalar excess",
        "routes=exact ray maxima -> crude fibre capacity -> rational common-K5 Farkas -> finite positive-slack packets -> lossless projected residual",
        f"global_counts=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};status_survivors:{totals[3]};literal_packets:{totals[4]};projected_survivors:0",
    ]
    for row in profiles:
        (
            first,
            body,
            h,
            components,
            L,
            lower,
            first_delta,
            first_d,
            ray_signs,
            ray_digest,
            divisor_count,
            trials,
            counts,
            stage_digests,
            body_cells,
            packet_count,
            killed_count,
            direct_controls,
            minimum_margin,
            maximum_cells,
            minimum_row,
            hostile_one_cell,
            direct_minimum,
            state_digest,
        ) = row
        lines.append(
            f"CASE;z1={first};E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};"
            f"lower={ftext(lower)};delta1={ftext(first_delta)};first_d={first_d};"
            f"ray_signs={ray_signs};ray_sha256={ray_digest};denominators={divisor_count};trials={trials};"
            f"counts={counts};stage_sha256={stage_digests};body_cells={body_cells};packets={packet_count};"
            f"kills={killed_count};direct_controls={direct_controls};min_margin={ftext(minimum_margin)};"
            f"max_prefix_cells={maximum_cells};minimum_row={minimum_row};hostile_one_cell={ftext(hostile_one_cell)};"
            f"direct_minimum_mass={ftext(direct_minimum)};state_sha256={state_digest};conclusion=EMPTY"
        )
    lines.extend(
        (
            "global_z1732_rows=2;empty=2;survivors=0;partition_complete=yes",
            f"necessary_ledger={prior_ledger}-2={final_ledger}",
            "splice=no occupied height at 1725..1731;no interpolation",
            "consequence=projected k=2 first drift label z1<=1724",
            f"profile_sha256={profile_digest}",
            f"semantic_sha256={semantic_digest}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(2, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1 and args.hash_seed >= 0, "bad replay options")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    if args.workers == 1:
        profiles = [case_profile(case) for case in CASES]
    else:
        with mp.get_context("spawn").Pool(min(args.workers, len(CASES))) as pool:
            profiles = list(pool.imap(case_profile, CASES, chunksize=1))
    order = {case: index for index, case in enumerate(CASES)}
    profiles.sort(key=lambda row: order[(row[0], row[1])])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
