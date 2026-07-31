#!/usr/bin/env python3
"""Independent hostile referee for the projected k=2, z1=1784 closure.

This program intentionally imports neither the primary z=1784 verifier nor
its ray/status/projected parent.  It reconstructs the exact carrier and
singleton integral from the defining danger combs, builds denominator rays,
enumerates all denominator multisets and literal scalar-eligible packets,
and implements the quotient and projected interval tests locally.

In particular, the former HIGH-TAIL body receives exactly the same all-label
enumeration as the other five rows: no projected-wall/high-label predicate is
read or applied.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
ATLAS = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1750_1799_thm2941.out"
)
Z1788 = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1788_ray_status_closure_thm2941.out"
)
Z1790 = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1790_exact_descent_closure_thm2941.out"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1784_independent_referee_thm2941.out"
)
EXPECTED_ATLAS_SHA256 = "00973f3679ebe5cf9e401a0311a9a19181ffcf4edb7c30f4c81b71f0d4a73d36"
EXPECTED_Z1788_SHA256 = "542970df3d4c901ac1e6120e69b768b9a7999cd4c778663d2bce99aafe8fa62c"
EXPECTED_Z1790_SHA256 = "56e1a1613ddd019cc97a1bdc5cda2ee2430ca8f87cc93d625ffdad1d304b9836"
EXPECTED_PROFILE_SHA256 = "c98cbb2b4e18457c521edb7f4d4f06668681866130190fa6d22b39bd566780b2"
EXPECTED_SEMANTIC_SHA256 = "d3e0d76ece1fe58845fa6043f67dc960b4843fbad42028cf71c63cdd2d66bf4e"

FIRST = 1784
MASTER = 14 * lcm(*range(1, 15))
LOWER_RATIO = Q(1, 91)
TWO_ALIGNED_CAP = Q(25, 91)
TAIL_SLOTS = 4
BODIES = (
    (1, 4, 8, 10, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (2, 8, 9, 10, 12, 14),
    (4, 6, 8, 10, 12, 14),
)
HIGH_TAIL_BODY = (2, 8, 9, 10, 12, 14)
EXPECTED = {
    BODIES[0]: (96, 28, 68, 0, 0, None, 0),
    BODIES[1]: (5, 4, 0, 1, 3, Q(1026, 16471), 2),
    BODIES[2]: (10, 1, 0, 9, 9, Q(66, 91), 1),
    BODIES[3]: (20, 3, 0, 17, 26, Q(723886, 3632447), 2),
    BODIES[4]: (7, 0, 0, 7, 7, Q(66, 91), 4),
    BODIES[5]: (2, 1, 0, 1, 1, Q(3980, 20293), 1),
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


require(digest(ATLAS) == EXPECTED_ATLAS_SHA256, "scalar atlas changed")
require(digest(Z1788) == EXPECTED_Z1788_SHA256, "z1788 closure output changed")
require(digest(Z1790) == EXPECTED_Z1790_SHA256, "z1790 closure output changed")


def ftext(value: Q | None) -> str:
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


def floor_q(value: Q) -> int:
    return value.numerator // value.denominator


def ceil_q(value: Q) -> int:
    return -((-value.numerator) // value.denominator)


def divisors(number: int) -> tuple[int, ...]:
    low: list[int] = []
    high: list[int] = []
    trial = 1
    while trial * trial <= number:
        if number % trial == 0:
            low.append(trial)
            if trial * trial != number:
                high.append(number // trial)
        trial += 1
    return tuple(low + high[::-1])


def merge_integer(rows: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    ordered = sorted((left, right) for left, right in rows if left < right)
    merged: list[list[int]] = []
    for left, right in ordered:
        if merged and left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def carrier_for(body: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    """Complement of the six open danger combs on the common integer ruler."""
    danger: list[tuple[int, int]] = []
    for speed in body:
        period = MASTER // speed
        radius = MASTER // (14 * speed)
        danger.append((0, radius))
        danger.extend(
            (index * period - radius, index * period + radius)
            for index in range(1, speed)
        )
        danger.append((MASTER - radius, MASTER))
    carrier: list[tuple[int, int]] = []
    cursor = 0
    for left, right in merge_integer(danger):
        if cursor < left:
            carrier.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < MASTER:
        carrier.append((cursor, MASTER))
    require(carrier and carrier[0][0] > 0 and carrier[-1][1] < MASTER, body)
    return tuple(carrier)


def primitive_one_comb(numerator: int) -> int:
    """MASTER times the integral of 1_{||t||<1/14} up to numerator/MASTER."""
    cycles, residue = divmod(numerator, MASTER)
    radius = MASTER // 14
    return (
        cycles * (MASTER // 7)
        + min(residue, radius)
        + max(0, residue - 13 * radius)
    )


def singleton_mass(carrier: tuple[tuple[int, int], ...], label: int) -> Q:
    numerator = sum(
        primitive_one_comb(label * right) - primitive_one_comb(label * left)
        for left, right in carrier
    )
    return Q(numerator, MASTER * label)


def excess(carrier: tuple[tuple[int, int], ...], h: Q, label: int) -> Q:
    return singleton_mass(carrier, label) - h / 7


def body_ranges(
    carrier: tuple[tuple[int, int], ...], L: int
) -> tuple[tuple[int, int], ...]:
    require(MASTER % L == 0, (L, "body ruler does not divide master"))
    scale = MASTER // L
    require(
        all(left % scale == right % scale == 0 for left, right in carrier),
        (L, "carrier endpoint off grid"),
    )
    return tuple((left // scale, right // scale) for left, right in carrier)


def cyclic_projection(
    modulus: int, ranges: tuple[tuple[int, int], ...]
) -> tuple[tuple[int, int], ...]:
    pieces: list[tuple[int, int]] = []
    for left, right in ranges:
        length = right - left
        if length >= modulus:
            return ((0, modulus),)
        start = left % modulus
        if start + length <= modulus:
            pieces.append((start, start + length))
        else:
            pieces.extend(((start, modulus), (0, start + length - modulus)))
    return merge_integer(pieces)


def load_histogram(
    arcs: tuple[tuple[int, int], ...], modulus: int
) -> tuple[tuple[int, int], ...]:
    base = 0
    events: defaultdict[int, int] = defaultdict(int)
    events[0] = 0
    events[modulus] = 0
    for left, right in arcs:
        quotient, remainder = divmod(right - left, modulus)
        base += quotient
        start = left % modulus
        finish = start + remainder
        if finish <= modulus:
            if start < finish:
                events[start] += 1
                events[finish] -= 1
        else:
            events[start] += 1
            events[modulus] -= 1
            events[0] += 1
            events[finish - modulus] -= 1
    histogram: defaultdict[int, int] = defaultdict(int)
    running = base
    points = sorted(events)
    for left, right in zip(points, points[1:]):
        running += events[left]
        histogram[running] += right - left
    result = tuple(sorted(histogram.items()))
    require(sum(count for _load, count in result) == modulus, result)
    require(
        sum(load * count for load, count in result)
        == sum(right - left for left, right in arcs),
        result,
    )
    return result


def fibre_capacity(D: int, d: int, q: int) -> int:
    teeth = (d + 6) // 7
    common = gcd(d, q)
    height = D // lcm(d, q)
    return height * ((teeth + common - 1) // common)


def prufer_tree_atlas(vertices: int) -> tuple[tuple[tuple[int, int], ...], ...]:
    trees: set[tuple[tuple[int, int], ...]] = set()
    for word in product(range(vertices), repeat=vertices - 2):
        degree = [1] * vertices
        for vertex in word:
            degree[vertex] += 1
        edges: list[tuple[int, int]] = []
        for vertex in word:
            leaf = next(index for index, value in enumerate(degree) if value == 1)
            edges.append(tuple(sorted((leaf, vertex))))
            degree[leaf] -= 1
            degree[vertex] -= 1
        leaves = [index for index, value in enumerate(degree) if value == 1]
        require(len(leaves) == 2, (word, leaves))
        edges.append(tuple(sorted(leaves)))
        trees.add(tuple(sorted(edges)))
    require(len(trees) == vertices ** (vertices - 2), len(trees))
    return tuple(sorted(trees))


TREES5 = prufer_tree_atlas(5)


def forced_pair_overlap(M: int, d: int, e: int, f: int, z: int) -> int:
    common = gcd(d, f)
    ae, re = divmod(e, common)
    az, rz = divmod(z, common)
    one_period = common * ae * az + ae * rz + az * re + max(0, re + rz - common)
    return (M // lcm(d, f)) * one_period


@lru_cache(maxsize=None)
def hunter_union_cap(M: int, ds: tuple[int, ...], es: tuple[int, ...]) -> int:
    sizes = tuple((M // d) * e for d, e in zip(ds, es))
    overlaps = {
        edge: forced_pair_overlap(M, ds[edge[0]], es[edge[0]], ds[edge[1]], es[edge[1]])
        for edge in combinations(range(5), 2)
    }
    invoice = max(sum(overlaps[edge] for edge in tree) for tree in TREES5)
    return min(M, sum(sizes) - invoice)


def status_signature(
    D: int, ds: tuple[int, ...], q: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    M = D // q
    inner: list[int] = []
    lows: list[int] = []
    marginals: list[int] = []
    for d in ds:
        common = gcd(d, q)
        low, remainder = divmod((d + 6) // 7, common)
        inner.append(d // common)
        lows.append(low)
        marginals.append((q // common) * remainder)
    capacities = tuple(
        hunter_union_cap(
            M,
            tuple(inner),
            tuple(lows[index] + ((pattern >> index) & 1) for index in range(5)),
        )
        for pattern in range(32)
    )
    return tuple(marginals), capacities


@lru_cache(maxsize=None)
def simultaneous_status_feasible(
    q: int,
    marginals: tuple[int, ...],
    capacities: tuple[int, ...],
    histogram: tuple[tuple[int, int], ...],
) -> tuple[bool, object]:
    """Every negative verdict is checked by an exact rational Farkas witness."""
    tail_rows: list[tuple[int, ...]] = []
    tail_rhs: list[int] = []
    thresholds: list[int] = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(count for load, count in histogram if load >= threshold)
        admissible = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(admissible):
            continue
        thresholds.append(threshold)
        if not any(admissible):
            return False, ((threshold,), (Q(1),), (Q(0),) * 6)
        tail_rows.append(admissible)
        tail_rhs.append(demand)

    equal_rows = [
        (1,) * 32,
        *[
            tuple((pattern >> index) & 1 for pattern in range(32))
            for index in range(5)
        ],
    ]
    equal_rhs = (q, *marginals)
    if not tail_rows:
        return True, None

    primal = linprog(
        np.zeros(32),
        A_ub=-np.asarray(tail_rows, dtype=float),
        b_ub=-np.asarray(tail_rhs, dtype=float),
        A_eq=np.asarray(equal_rows, dtype=float),
        b_eq=np.asarray(equal_rhs, dtype=float),
        bounds=(0, None),
        method="highs",
    )
    if primal.success:
        return True, None
    require(primal.status == 2, ("unexpected primal status", primal.status))

    tail_count = len(tail_rows)
    dual_rows: list[tuple[int, ...]] = []
    for pattern in range(32):
        dual_rows.append(
            tuple(tail_rows[row][pattern] for row in range(tail_count))
            + tuple(-equal_rows[row][pattern] for row in range(6))
        )
    dual_rows.append(tuple(-value for value in tail_rhs) + tuple(equal_rhs))
    dual = linprog(
        np.zeros(tail_count + 6),
        A_ub=np.asarray(dual_rows, dtype=float),
        b_ub=np.asarray((0,) * 32 + (-1,), dtype=float),
        bounds=[(0, None)] * tail_count + [(None, None)] * 6,
        method="highs",
    )
    if not dual.success:
        return True, None
    alpha = tuple(Q(float(value)).limit_denominator(1_000_000) for value in dual.x[:tail_count])
    z = tuple(Q(float(value)).limit_denominator(1_000_000) for value in dual.x[tail_count:])
    slacks = tuple(
        sum(z[row] * equal_rows[row][pattern] for row in range(6))
        - sum(alpha[row] * tail_rows[row][pattern] for row in range(tail_count))
        for pattern in range(32)
    )
    contradiction = (
        sum(z[row] * equal_rhs[row] for row in range(6))
        - sum(alpha[row] * tail_rhs[row] for row in range(tail_count))
    )
    if all(value >= 0 for value in alpha) and all(value >= 0 for value in slacks) and contradiction < 0:
        return False, (tuple(thresholds), alpha, z, contradiction)
    return True, None


def merge_q(rows: list[tuple[Q, Q]]) -> tuple[tuple[Q, Q], ...]:
    ordered = sorted((left, right) for left, right in rows if left < right)
    merged: list[list[Q]] = []
    for left, right in ordered:
        if merged and left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def intersect_q(
    first: tuple[tuple[Q, Q], ...], second: tuple[tuple[Q, Q], ...]
) -> tuple[tuple[Q, Q], ...]:
    rows: list[tuple[Q, Q]] = []
    i = 0
    j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            rows.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(rows)


def mass(rows: tuple[tuple[Q, Q], ...]) -> Q:
    return sum((right - left for left, right in rows), Q())


def cell_danger(cell: int, label: int, L: int) -> tuple[tuple[Q, Q], ...]:
    start = Q(label * cell, L)
    finish = Q(label * (cell + 1), L)
    rows: list[tuple[Q, Q]] = []
    for tooth in range(floor_q(start) - 1, ceil_q(finish) + 2):
        left = max(Q(), Q(L, label) * (Q(tooth) - Q(1, 14)) - cell)
        right = min(Q(1), Q(L, label) * (Q(tooth) + Q(1, 14)) - cell)
        if left < right:
            rows.append((left, right))
    return merge_q(rows)


def projected_prefix(
    cells: tuple[int, ...], L: int, labels: tuple[int, ...], stop: bool
) -> tuple[Q, int]:
    common: tuple[tuple[Q, Q], ...] = ((Q(), Q(1)),)
    used = 0
    for cell in cells:
        local = merge_q(
            [interval for label in labels for interval in cell_danger(cell, label, L)]
        )
        common = intersect_q(common, local)
        used += 1
        if stop and mass(common) <= 1 - TWO_ALIGNED_CAP:
            break
    return 1 - mass(common), used


def danger_on_circle(label: int) -> tuple[tuple[Q, Q], ...]:
    radius = Q(1, 14 * label)
    rows = [(Q(), radius), (Q(1) - radius, Q(1))]
    rows.extend(
        (Q(index, label) - radius, Q(index, label) + radius)
        for index in range(1, label)
    )
    return merge_q(rows)


def subtract_q(
    source: tuple[tuple[Q, Q], ...], removed: tuple[tuple[Q, Q], ...]
) -> tuple[tuple[Q, Q], ...]:
    result: list[tuple[Q, Q]] = []
    j = 0
    for left, right in source:
        cursor = left
        while j < len(removed) and removed[j][1] <= left:
            j += 1
        k = j
        while k < len(removed) and removed[k][0] < right:
            if cursor < removed[k][0]:
                result.append((cursor, min(right, removed[k][0])))
            cursor = max(cursor, removed[k][1])
            if cursor >= right:
                break
            k += 1
        if cursor < right:
            result.append((cursor, right))
    return tuple(result)


def physical_projected_mass(
    carrier: tuple[tuple[int, int], ...], L: int, labels: tuple[int, ...]
) -> Q:
    carrier_q = tuple((Q(left, MASTER), Q(right, MASTER)) for left, right in carrier)
    removed = merge_q(
        [interval for label in labels for interval in danger_on_circle(label)]
    )
    residual = subtract_q(carrier_q, removed)
    projected: list[tuple[Q, Q]] = []
    for left, right in residual:
        scaled_left = L * left
        scaled_right = L * right
        for integer in range(floor_q(scaled_left), ceil_q(scaled_right)):
            piece_left = max(scaled_left, Q(integer)) - integer
            piece_right = min(scaled_right, Q(integer + 1)) - integer
            if piece_left < piece_right:
                projected.append((piece_left, piece_right))
    return mass(merge_q(projected))


def profile(body: tuple[int, ...]) -> tuple[object, ...]:
    carrier = carrier_for(body)
    h = Q(sum(right - left for left, right in carrier), MASTER)
    lower = h * LOWER_RATIO
    L = 14 * lcm(*body)
    ranges = body_ranges(carrier, L)
    cells = tuple(cell for left, right in ranges for cell in range(left, right))
    first_delta = excess(carrier, h, FIRST)
    first_d = L // gcd(L, FIRST)

    amplitudes = [Q()]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * excess(carrier, h, residue)
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(
            (residue + L) * excess(carrier, h, residue + L) == amplitudes[residue]
            for residue in range(1, L)
        ),
        (body, "ray recurrence"),
    )
    require(
        all(amplitudes[L - residue] == -amplitudes[residue] for residue in range(1, L)),
        (body, "ray antipodes"),
    )
    require(L * excess(carrier, h, L) == 0, (body, "aligned ray"))
    require(signs[1] == signs[-1], (body, signs))

    ds_universe = tuple(d for d in divisors(L) if d > 1)

    @lru_cache(maxsize=None)
    def heads(d: int, multiplicity: int) -> tuple[tuple[Q, int, int], ...]:
        candidates: list[tuple[Q, int, int]] = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= FIRST:
                label += ((FIRST + 1 - label + L - 1) // L) * L
            candidates.extend(
                (amplitude / (label + offset * L), label + offset * L, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        require(len(candidates) >= multiplicity, (body, d, multiplicity))
        return tuple(candidates[:multiplicity])

    scalar: list[tuple[tuple[int, ...], Q, tuple[int, ...]]] = []
    trials = 0
    for tail in combinations_with_replacement(ds_universe, TAIL_SLOTS):
        trials += 1
        upper = first_delta
        labels: list[int] = []
        for d, multiplicity in Counter(tail).items():
            selected = heads(d, multiplicity)
            upper += sum((value for value, _label, _residue in selected), Q())
            labels.extend(label for _value, label, _residue in selected)
        if upper >= lower:
            scalar.append((tuple(sorted((first_d, *tail))), upper, (FIRST, *sorted(labels))))
    require(trials == comb(len(ds_universe) + 3, 4), (body, trials))
    require(len({row[0] for row in scalar}) == len(scalar), (body, "duplicate state"))

    arc_cache: dict[int, tuple[tuple[int, int], ...]] = {}
    crude_kills: list[tuple[object, ...]] = []
    crude_survivors: list[tuple[object, ...]] = []
    for ds, upper, labels in scalar:
        D = lcm(*ds)
        arcs = arc_cache.setdefault(D, cyclic_projection(D, ranges))
        witness = None
        for q in divisors(D):
            histogram = load_histogram(arcs, q)
            target = max(load for load, count in histogram if count)
            capacity = sum(fibre_capacity(D, d, q) for d in ds)
            if target > capacity:
                witness = (q, target, capacity)
                break
        row = (ds, upper, labels, witness)
        (crude_kills if witness else crude_survivors).append(row)

    status_kills: list[tuple[object, ...]] = []
    survivors: list[tuple[object, ...]] = []
    for ds, upper, labels, _ in crude_survivors:
        D = lcm(*ds)
        arcs = arc_cache[D]
        witness = None
        for M in divisors(D):
            q = D // M
            marginals, capacities = status_signature(D, ds, q)
            histogram = load_histogram(arcs, q)
            feasible, certificate = simultaneous_status_feasible(
                q, marginals, capacities, histogram
            )
            if not feasible:
                require(certificate is not None, (body, ds, q))
                witness = (q, M, certificate)
                break
        row = (ds, upper, labels, witness)
        (status_kills if witness else survivors).append(row)

    packet_count = 0
    projected_kills = 0
    minimum_margin: Q | None = None
    maximum_prefix = 0
    minimum_labels: tuple[int, ...] | None = None
    packet_rows: list[tuple[object, ...]] = []

    @lru_cache(maxsize=None)
    def literal_choices(d: int, multiplicity: int, slack: Q) -> tuple[tuple[Q, tuple[int, ...], Q], ...]:
        top = heads(d, multiplicity)
        top_sum = sum((row[0] for row in top), Q())
        threshold = top[-1][0] - slack
        require(threshold > 0, (body, d, multiplicity, slack, threshold))
        candidates: list[tuple[Q, int, int]] = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude <= 0:
                continue
            label = residue
            if label <= FIRST:
                label += ((FIRST + 1 - label + L - 1) // L) * L
            while amplitude / label >= threshold:
                candidates.append((amplitude / label, label, residue))
                label += L
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        choices: list[tuple[Q, tuple[int, ...], Q]] = []
        for chosen in combinations(candidates, multiplicity):
            chosen_labels = tuple(sorted(row[1] for row in chosen))
            if len(set(chosen_labels)) != multiplicity:
                continue
            value = sum((row[0] for row in chosen), Q())
            deficit = top_sum - value
            if deficit <= slack:
                choices.append((deficit, chosen_labels, value))
        choices.sort()
        require(choices and choices[0][0] == 0, (body, d, multiplicity))
        return tuple(choices)

    for state_index, (ds, upper, maximizing_labels, witness) in enumerate(survivors, 1):
        require(witness is None, (body, ds))
        slack = upper - lower
        tail = list(ds)
        tail.remove(first_d)
        groups = [
            (d, multiplicity, literal_choices(d, multiplicity, slack))
            for d, multiplicity in sorted(Counter(tail).items())
        ]
        state_packets = 0
        for selection in product(*(group[2] for group in groups)):
            deficit = sum((choice[0] for choice in selection), Q())
            if deficit > slack:
                continue
            labels = (FIRST, *sorted(label for choice in selection for label in choice[1]))
            require(
                len(labels) == 5 and len(set(labels)) == 5 and all(label > FIRST for label in labels[1:]),
                (body, labels),
            )
            actual = first_delta + sum((choice[2] for choice in selection), Q())
            require(actual == upper - deficit and actual >= lower, (body, labels, actual))
            projected, prefix = projected_prefix(cells, L, labels, True)
            margin = projected - TWO_ALIGNED_CAP
            packet_count += 1
            state_packets += 1
            maximum_prefix = max(maximum_prefix, prefix)
            if minimum_margin is None or (margin, labels) < (minimum_margin, minimum_labels):
                minimum_margin = margin
                minimum_labels = labels
            if margin > 0:
                projected_kills += 1
            packet_rows.append((state_index, ds, labels, actual, margin, prefix))
        require(state_packets > 0, (body, ds, maximizing_labels))

    direct_control = None
    full_cell_control = None
    if minimum_labels is not None:
        full_cell_control, used = projected_prefix(cells, L, minimum_labels, False)
        require(used == len(cells), (body, used, len(cells)))
        direct_control = physical_projected_mass(carrier, L, minimum_labels)
        require(direct_control == full_cell_control, (body, direct_control, full_cell_control))

    counts = (len(scalar), len(crude_kills), len(status_kills), len(survivors), packet_count)
    expected = EXPECTED[body]
    require(counts == expected[:5], (body, "count mismatch", counts, expected[:5]))
    require(projected_kills == packet_count, (body, "projected survivor"))
    require(minimum_margin == expected[5], (body, "margin mismatch", minimum_margin, expected[5]))
    require(maximum_prefix == expected[6], (body, "prefix mismatch", maximum_prefix, expected[6]))

    # The exact Farkas representative is solver-selected and therefore not a
    # canonical part of the infeasible status instance (MISTAKE-331).  Keep
    # verifying it above, but hash only q and M at this evidence boundary.
    canonical_status_kills = tuple(
        (ds, upper, labels, witness[:-1])
        for ds, upper, labels, witness in status_kills
    )
    require(
        not status_kills
        or sha256(repr(tuple(status_kills)).encode()).hexdigest()
        != sha256(repr(canonical_status_kills).encode()).hexdigest(),
        (body, "raw LP certificate entered status digest"),
    )
    stage_payloads = (scalar, crude_kills, canonical_status_kills, survivors, packet_rows)
    stage_hashes = tuple(sha256(repr(tuple(rows)).encode()).hexdigest() for rows in stage_payloads)
    return (
        body,
        h,
        len(carrier),
        L,
        len(cells),
        lower,
        first_delta,
        first_d,
        len(ds_universe),
        trials,
        counts,
        projected_kills,
        minimum_margin,
        maximum_prefix,
        minimum_labels,
        full_cell_control,
        direct_control,
        body == HIGH_TAIL_BODY,
        False,
        stage_hashes,
    )


def atlas_audit() -> tuple[object, ...]:
    rows: list[tuple[int, tuple[int, ...], bool]] = []
    for line in ATLAS.read_text().splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
        body = tuple(map(int, fields["E"].split(",")))
        first = int(fields["z1"])
        rows.append((first, body, "HIGH-TAIL" in fields["suffix"]))
    z1784 = tuple(body for first, body, _high in rows if first == FIRST)
    high = tuple(body for first, body, is_high in rows if first == FIRST and is_high)
    require(z1784 == BODIES, ("atlas z1784 rows", z1784))
    require(high == (HIGH_TAIL_BODY,), ("atlas high row", high))
    require("global_z1788_rows=1;empty=1;survivors=0" in Z1788.read_text(), "z1788 closure missing")
    require("global_z1790_rows=8;empty=8;survivors=0" in Z1790.read_text(), "z1790 closure missing")
    require("all_exact_controls=PASS" in Z1788.read_text(), "z1788 controls missing")
    require("all_exact_controls=PASS" in Z1790.read_text(), "z1790 controls missing")
    next_occupied = max(first for first, _body, _high in rows if first < FIRST)
    require(next_occupied == 1780, next_occupied)
    return len(rows), tuple(sorted(Counter(first for first, _body, _high in rows).items())), next_occupied


def render(profiles: tuple[tuple[object, ...], ...]) -> str:
    require(tuple(row[0] for row in profiles) == BODIES, "body order")
    atlas_rows, atlas_heights, next_occupied = atlas_audit()
    totals = tuple(sum(row[10][index] for row in profiles) for index in range(5))
    total_kills = sum(row[11] for row in profiles)
    require(totals == (140, 37, 68, 35, 46), totals)
    require(total_kills == totals[-1] == 46, total_kills)
    minimum = min(row[12] for row in profiles if row[12] is not None)
    require(minimum == Q(1026, 16471), minimum)
    require(sum(int(row[17]) for row in profiles) == 1, "HIGH-TAIL row count")
    require(not any(row[18] for row in profiles), "a high side filter was applied")
    profile_hash = sha256(repr(profiles).encode()).hexdigest()
    require(
        profile_hash == EXPECTED_PROFILE_SHA256,
        f"profile digest changed: expected {EXPECTED_PROFILE_SHA256}; got {profile_hash}",
    )
    semantic_payload = (
        FIRST,
        BODIES,
        totals,
        total_kills,
        minimum,
        next_occupied,
        tuple((row[0], row[10], row[12], row[13], row[17], row[18]) for row in profiles),
        digest(ATLAS),
        digest(Z1788),
        digest(Z1790),
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    require(
        semantic_hash == EXPECTED_SEMANTIC_SHA256,
        f"semantic digest changed: expected {EXPECTED_SEMANTIC_SHA256}; got {semantic_hash}",
    )
    lines = [
        "LRC14 projected k=2 z1784 independent hostile referee",
        f"referee_source_sha256={digest(SOURCE)}",
        "independence=no import of primary z1784 or parent ray/status/projected engines",
        "construction=defining integer danger combs;local singleton primitive;local quotient/status/projected geometry",
        "scope=all distinct later nonaligned labels;no finite label horizon;no HIGH-TAIL side condition",
        f"atlas_sha256={digest(ATLAS)}",
        f"z1788_output_sha256={digest(Z1788)}",
        f"z1790_output_sha256={digest(Z1790)}",
        f"atlas_rows={atlas_rows};height_counts={atlas_heights};next_occupied_after_z1784={next_occupied}",
        f"global_counts=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};packets:{totals[4]};projected_kills:{total_kills}",
        f"global_minimum_positive_margin={ftext(minimum)}",
    ]
    for row in profiles:
        (
            body, h, components, L, cells, lower, delta1, first_d,
            divisor_count, trials, counts, kills, min_margin, max_prefix,
            min_labels, full_control, direct_control, atlas_high, high_filter,
            stage_hashes,
        ) = row
        lines.append(
            "CASE;"
            f"E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};cells={cells};"
            f"lower={ftext(lower)};delta1={ftext(delta1)};first_d={first_d};"
            f"denominators={divisor_count};trials={trials};counts={counts};kills={kills};"
            f"min_margin={ftext(min_margin)};max_prefix_cells={max_prefix};"
            f"minimum_labels={min_labels};full_cell_control={ftext(full_control)};"
            f"physical_control={ftext(direct_control)};atlas_HIGH_TAIL={atlas_high};"
            f"high_side_filter_applied={high_filter};stage_sha256={stage_hashes}"
        )
    lines.extend(
        (
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(6, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    if args.workers == 1:
        rows = [profile(body) for body in BODIES]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.map(profile, BODIES))
    rows.sort(key=lambda row: row[0])
    output = render(tuple(rows))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
