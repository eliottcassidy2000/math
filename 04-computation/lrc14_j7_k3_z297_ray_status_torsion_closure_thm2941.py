#!/usr/bin/env python3
"""Close the projected k=3, z1=297 slice by located p-torsion.

The lossless scalar atlas has seven body rows at first drift 297.  Exact
periodic rays and the common 16-status Hunter transport leave 71 denominator
states on three bodies.  The pinned projected-wall branch requires at least
one later label at or above floor(13L/132)+1.  A duplicate-permitting exact
ray upper bound puts every packet with at least two high labels strictly below
the scalar wall, so every possible packet has exactly one high label.

All labels below the projected wall are then finite for a structural reason,
not an imposed horizon.  Negative-amplitude low labels are retained.  Exact
enumeration leaves 73 (state, high denominator, low pair) cases involving
only seven body-distinct low pairs.  For every case the set S of distinct
fixed-safe cell residues modulo the high denominator d satisfies |S|>d/r
for some divisor 2<=r<=7 of d.  Two residues therefore share a coset modulo
d/r.  Their nonzero difference has effective torsion order at most r, so
every high label z=(L/d)u+hL, gcd(u,d)=1, separates their phases by at least
1/r>=1/7.  The two strict-open radius-1/14 danger sets are pointwise
disjoint, including the endpoint-touching r=7 boundary.  The projected
residual has full mass one and cannot be covered by the three aligned combs,
whose union cap is 36/91.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import Counter, defaultdict
from fractions import Fraction as Q
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_PATH = ROOT / "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
ATLAS = ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.out"
UPSTREAM = ROOT / "04-computation/lrc14_j7_k3_z306_z302_z298_ray_status_descent_closure_thm2941.py"
UPSTREAM_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z306_z302_z298_ray_status_descent_closure_thm2941.out"
FIRST = 297
PROJECTED_WALL = Q(13, 132)
THREE_ALIGNED_CAP = Q(36, 91)
EXPECTED_RAY_SHA256 = "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2"
EXPECTED_ATLAS_SOURCE_SHA256 = "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded"
EXPECTED_ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
EXPECTED_UPSTREAM_SHA256 = "1e06537c7a9f194b93c6913ed87559027e18104f46bfad4347aa964cfe6b29f2"
EXPECTED_UPSTREAM_OUTPUT_SHA256 = "2a5845c36712c5b91991bed27ddf578150e11b73fe61972494618765e3c476ef"
EXPECTED_SEMANTIC_SHA256 = "52d02b6f0cec464f24b0e3076f1dee5bd954a41d07d8bc42b1b1bc51d9cac768"
INHERITED_LEDGER = 375_765
FINAL_LEDGER = 375_758
EXPECTED_BODIES = (
    (1, 4, 10, 11, 12, 14),
    (1, 8, 10, 11, 12, 13),
    (1, 8, 10, 11, 12, 14),
    (1, 10, 11, 12, 13, 14),
    (2, 6, 8, 10, 11, 14),
    (2, 8, 10, 11, 12, 13),
    (2, 8, 10, 11, 12, 14),
)
EXPECTED_COUNTS = {
    EXPECTED_BODIES[0]: (73, 0, 72, 1),
    EXPECTED_BODIES[1]: (317, 195, 122, 0),
    EXPECTED_BODIES[2]: (545, 76, 423, 46),
    EXPECTED_BODIES[3]: (0, 0, 0, 0),
    EXPECTED_BODIES[4]: (0, 0, 0, 0),
    EXPECTED_BODIES[5]: (0, 0, 0, 0),
    EXPECTED_BODIES[6]: (237, 0, 213, 24),
}
EXPECTED_GAPS = {
    EXPECTED_BODIES[0]: Q(3426764729, 1136083728780),
    EXPECTED_BODIES[2]: Q(30957017759, 11803465764090),
    EXPECTED_BODIES[6]: Q(1791956807, 582618246210),
}
EXPECTED_TERMINAL = {
    EXPECTED_BODIES[0]: (1, 1, 1),
    EXPECTED_BODIES[2]: (45, 48, 5),
    EXPECTED_BODIES[6]: (24, 24, 1),
}
EXPECTED_QUALIFYING_HISTOGRAM = ((2, 49), (3, 14), (4, 6), (6, 4))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(RAY_PATH) == EXPECTED_RAY_SHA256, "ray engine changed")
require(file_sha256(ATLAS_SOURCE) == EXPECTED_ATLAS_SOURCE_SHA256, "atlas source changed")
require(file_sha256(ATLAS) == EXPECTED_ATLAS_SHA256, "atlas output changed")
require(file_sha256(UPSTREAM) == EXPECTED_UPSTREAM_SHA256, "z298 closure changed")
require(
    file_sha256(UPSTREAM_OUTPUT) == EXPECTED_UPSTREAM_OUTPUT_SHA256,
    "z298 closure output changed",
)
require(
    "projected k3 cap<=297;next exact frontier=297"
    in UPSTREAM_OUTPUT.read_text(),
    "upstream cap statement changed",
)
ray = load("z297_terminal_ray", RAY_PATH)
ray.FIRST = FIRST


def atlas_data():
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);"
        r"z1=([0-9]+);.*;suffix=(.*);lower="
    )
    counts = Counter()
    frontier_rows = []
    for line in ATLAS.read_text().splitlines():
        match = pattern.match(line)
        if not match:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        L = int(match.group(2))
        high = int(match.group(3))
        first = int(match.group(4))
        suffix = match.group(5)
        counts[first] += 1
        if first == FIRST:
            derived_high = (
                PROJECTED_WALL.numerator * L // PROJECTED_WALL.denominator + 1
            )
            require(
                high == derived_high and FIRST < high,
                (body, L, high, derived_high),
            )
            require(
                "HIGH-TAIL:" in suffix,
                (body, "projected high-wall branch absent"),
            )
            frontier_rows.append((body, L, high))
    bodies = tuple(row[0] for row in frontier_rows)
    require(bodies == EXPECTED_BODIES, bodies)
    require(counts[FIRST] == 7, counts[FIRST])
    require(
        counts[296] == counts[295] == 0 and counts[294] == 1,
        tuple(sorted(counts.items())),
    )
    return (
        bodies,
        tuple(frontier_rows),
        tuple(sorted(counts.items())),
        294,
        counts[294],
    )


def independent_farkas_check(q, marginals, capacities, histogram, certificate):
    """Rebuild the common-table inequalities and audit the witness over Q."""
    thresholds, alpha, z = certificate
    tail_rows = []
    tail_rhs = []
    rebuilt_thresholds = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(count for load, count in histogram if load >= threshold)
        good = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(good):
            continue
        rebuilt_thresholds.append(threshold)
        tail_rows.append(good)
        tail_rhs.append(demand)
    equality_rows = [
        (1,) * 16,
        *[
            tuple((pattern >> index) & 1 for pattern in range(16))
            for index in range(4)
        ],
    ]
    equality_rhs = (q, *marginals)
    require(tuple(rebuilt_thresholds) == tuple(thresholds), "threshold mismatch")
    require(len(alpha) == len(tail_rows) and len(z) == 5, "Farkas shape")
    require(all(value >= 0 for value in alpha), "negative Farkas alpha")
    slacks = tuple(
        sum(z[row] * equality_rows[row][pattern] for row in range(5))
        - sum(alpha[row] * tail_rows[row][pattern] for row in range(len(alpha)))
        for pattern in range(16)
    )
    contradiction = sum(z[row] * equality_rhs[row] for row in range(5))
    contradiction -= sum(
        alpha[row] * tail_rhs[row] for row in range(len(alpha))
    )
    require(all(value >= 0 for value in slacks), "negative Farkas column")
    require(contradiction < 0, "nonnegative Farkas contradiction")
    return contradiction


def exact_common_status_screen(stream, states):
    """Replay the pinned screen while genuinely caching geometric fibres.

    The upstream routine is mathematically exact, but its ``setdefault``
    argument eagerly rebuilds an already cached support arc.  This local
    spelling makes the finite universe and the cache key explicit; it changes
    no witness or rejection order.
    """
    crude_kills = {}
    status_kills = {}
    survivors = []
    arcs_cache = {}
    histogram_cache = {}
    for ds in sorted(states):
        killed, witness = stream.pattern_screen(ds)
        if killed:
            crude_kills[ds] = witness
            continue
        D = lcm(*ds)
        if D not in arcs_cache:
            arcs_cache[D] = ray.fibre.projected_support_arcs(D, stream.ranges)
        arcs = arcs_cache[D]
        status_witness = None
        for M in ray.support.divisors(D):
            q = D // M
            marginals, capacities = ray.local.hunter_status_data(D, ds, q)
            key = (D, q)
            if key not in histogram_cache:
                histogram_cache[key] = ray.fibre.residue_load_histogram(arcs, q)
            histogram = histogram_cache[key]
            feasible, certificate = ray.local.common_status_feasible(
                q,
                marginals,
                capacities,
                histogram,
            )
            if not feasible:
                status_witness = (
                    q,
                    M,
                    marginals,
                    tuple(sorted(set(capacities))),
                    histogram,
                    certificate,
                )
                break
        if status_witness is None:
            survivors.append(ds)
        else:
            status_kills[ds] = status_witness
    return crude_kills, status_kills, tuple(survivors)


def evaluate(body):
    stream = ray.Stream(body)
    require(
        stream.high_floor
        == PROJECTED_WALL.numerator * stream.L // PROJECTED_WALL.denominator + 1,
        (body, "high floor changed"),
    )
    trials, states, checks, signs = ray.ray_quotient_states(stream)
    crude, status, survivors = exact_common_status_screen(stream, states)
    require(set(states) == set(crude) | set(status) | set(survivors), body)
    require(not (set(crude) & set(status)), (body, "kill stages overlap"))
    require(
        (len(states), len(crude), len(status), len(survivors))
        == EXPECTED_COUNTS[body],
        (body, len(states), len(crude), len(status), len(survivors)),
    )
    for ds, witness in crude.items():
        gap, q, M, target, capacity = witness
        D = lcm(*ds)
        require(M == D // q, (body, ds, "crude cofactor"))
        require(target == dict(stream.target_data(D))[q], (body, ds, "target"))
        require(
            capacity == sum(ray.local.fibre_cap(D, d, q) for d in ds),
            (body, ds, "capacity"),
        )
        require(gap == target - capacity and gap > 0, (body, ds, witness))
    contradictions = []
    verified_instance_digest = hashlib.sha256()
    representative_instance = None
    arcs_cache = {}
    histogram_cache = {}
    for ds, witness in sorted(status.items()):
        q, M, marginals, cap_set, histogram, certificate = witness
        D = lcm(*ds)
        require(M == D // q, (body, ds, "status cofactor"))
        rebuilt_marginals, capacities = ray.local.hunter_status_data(D, ds, q)
        require(rebuilt_marginals == marginals, (body, ds, "marginals"))
        require(tuple(sorted(set(capacities))) == cap_set, (body, ds, "caps"))
        if D not in arcs_cache:
            arcs_cache[D] = ray.fibre.projected_support_arcs(D, stream.ranges)
        arcs = arcs_cache[D]
        key = (D, q)
        if key not in histogram_cache:
            histogram_cache[key] = ray.fibre.residue_load_histogram(arcs, q)
        require(
            histogram_cache[key] == histogram,
            (body, ds, "histogram"),
        )
        contradictions.append(
            independent_farkas_check(q, marginals, capacities, histogram, certificate)
        )
        instance = witness[:-1]
        if representative_instance is None:
            representative_instance = (body, ds, instance)
        verified_instance_digest.update(f"{body}|{ds}|{instance}\n".encode())
    packets = []
    for ds in survivors:
        state = states[ds]
        labels = state["labels"]
        require(tuple(sorted(stream.L // gcd(stream.L, label) for label in labels)) == ds, (body, ds, labels))
        require(labels[0] == FIRST and len(labels) == len(set(labels)) == 4, labels)
        require(
            sum(label >= stream.high_floor for label in labels[1:]) == 1,
            (body, ds, labels, "maximizer high grammar"),
        )
        packets.append((ds, labels, state["total"], state["excess"]))
    sign_totals = {
        sign: sum(count for (_d, candidate_sign), count in signs.items() if candidate_sign == sign)
        for sign in (-1, 0, 1)
    }
    require(sign_totals[-1] == sign_totals[1], (body, sign_totals))
    return (
        body,
        stream.L,
        stream.high_floor,
        stream.first_d,
        trials,
        checks,
        tuple(sorted(sign_totals.items())),
        tuple(sorted(states)),
        tuple(sorted(crude)),
        tuple(sorted(status)),
        tuple(packets),
        tuple(sorted(Counter(witness[1] for witness in status.values()).items())),
        verified_instance_digest.hexdigest(),
        len(contradictions),
        representative_instance,
    )


def delta(stream, label):
    return ray.suffix.A.singleton_coverage(stream.carrier, label) - stream.h / 7


def first_on_ray(residue, modulus, threshold):
    if residue >= threshold:
        return residue
    return residue + ((threshold - residue + modulus - 1) // modulus) * modulus


def suffix_slots(ds, first_d):
    slots = list(ds)
    slots.remove(first_d)
    require(len(slots) == 3, (ds, first_d, slots))
    return tuple(slots)


def build_literal_tables(stream, needed):
    """Keep every finite low label, including zero/negative amplitudes."""
    low = {d: [] for d in needed}
    for label in range(stream.first + 1, stream.high_floor):
        d = stream.L // gcd(stream.L, label)
        if d in low:
            low[d].append((delta(stream, label), label))
    sign_census = Counter()
    for d in low:
        sign_census.update(
            (value > 0) - (value < 0) for value, _label in low[d]
        )
        low[d].sort(key=lambda item: (-item[0], item[1]))
        low[d] = tuple(low[d])

    # Above the wall, each unit residue ray is K/z.  Its exact maximum is at
    # the first high point when K>0; K<=0 has rigorous upper bound zero.  This
    # finite unit-ray quotient seals every arbitrarily large label, so no
    # search horizon is present.
    high = {}
    recurrence_checks = 0
    for d in needed:
        step = stream.L // d
        best = (Q(0), None, None, None)
        for unit in range(1, d):
            if gcd(unit, d) != 1:
                continue
            residue = step * unit
            amplitude = residue * delta(stream, residue)
            require(
                (residue + stream.L) * delta(stream, residue + stream.L)
                == amplitude,
                (stream.body, d, unit, "ray recurrence"),
            )
            recurrence_checks += 1
            label = first_on_ray(residue, stream.L, stream.high_floor)
            value = amplitude / label if amplitude > 0 else Q(0)
            candidate = (value, label, unit, amplitude)
            if candidate[0] > best[0] or (
                candidate[0] == best[0]
                and candidate[1] is not None
                and (best[1] is None or candidate[1] < best[1])
            ):
                best = candidate
        require(best[1] is not None, (stream.body, d, "empty high class"))
        high[d] = best
    return low, high, tuple(sorted(sign_census.items())), recurrence_checks


def duplicate_two_high_gap(stream, residuals, low, high):
    """Upper-bound >=2-high packets even while allowing repeated labels."""
    required = stream.lower - stream.first_delta
    best = None
    witness = None
    for ds in residuals:
        slots = suffix_slots(ds, stream.first_d)
        for mask in range(8):
            if mask.bit_count() < 2:
                continue
            value = Q()
            packet = []
            possible = True
            for index, d in enumerate(slots):
                if (mask >> index) & 1:
                    row = high[d]
                elif low[d]:
                    row = (*low[d][0], None, None)
                else:
                    possible = False
                    break
                value += row[0]
                packet.append((d, row[1], bool((mask >> index) & 1), row[0]))
            if possible and (best is None or value > best):
                best = value
                witness = (ds, mask, tuple(packet))
    require(best is not None, stream.body)
    return required - best, witness


def zero_high_scalar_passes(stream, residuals, low):
    """Hostile control excluded only by the inherited projected high gate."""
    required = stream.lower - stream.first_delta
    rows = []
    for ds in residuals:
        slots = suffix_slots(ds, stream.first_d)
        if not all(low[d] for d in slots):
            continue
        # Duplicates are deliberately allowed: this can only enlarge the
        # hostile zero-high upper and is not used to reject a state.
        value = sum((low[d][0][0] for d in slots), Q())
        if value >= required:
            rows.append((ds, tuple(low[d][0][1] for d in slots), value - required))
    return tuple(rows)


def finite_low_pairs(first, second, same, threshold):
    rows = []
    if same:
        for index, left in enumerate(first):
            if (
                index + 1 >= len(first)
                or left[0] + first[index + 1][0] < threshold
            ):
                break
            for right in first[index + 1 :]:
                if left[0] + right[0] < threshold:
                    break
                if left[1] != right[1]:
                    rows.append((left, right))
    else:
        for left in first:
            if not second or left[0] + second[0][0] < threshold:
                break
            for right in second:
                if left[0] + right[0] < threshold:
                    break
                if left[1] != right[1]:
                    rows.append((left, right))
    return tuple(rows)


def one_high_cases(stream, residuals, low, high):
    required = stream.lower - stream.first_delta
    cases = set()
    for ds in residuals:
        slots = suffix_slots(ds, stream.first_d)
        for high_index, high_d in enumerate(slots):
            low_ds = tuple(
                d for index, d in enumerate(slots) if index != high_index
            )
            if not low[low_ds[0]] or not low[low_ds[1]]:
                continue
            threshold = required - high[high_d][0]
            pairs = finite_low_pairs(
                low[low_ds[0]],
                low[low_ds[1]],
                low_ds[0] == low_ds[1],
                threshold,
            )
            for left, right in pairs:
                low_rows = tuple(
                    sorted(((low_ds[0], left[1]), (low_ds[1], right[1])))
                )
                excess = left[0] + right[0] + high[high_d][0] - required
                cases.add((ds, high_d, low_rows, excess))
    return tuple(sorted(cases))


def cell_clean(cell, label, L):
    residue = (label * cell) % L
    # The closed body cell misses the strict danger comb exactly when both
    # weak endpoint inequalities hold.  Every fixed label is <13L/132<L.
    return 14 * residue >= L and 14 * (residue + label) <= 13 * L


def fixed_safe_cells(stream, low_labels):
    labels = (stream.first, *low_labels)
    require(all(label < stream.high_floor < stream.L for label in labels), labels)
    return tuple(
        cell
        for left, right in stream.ranges
        for cell in range(left, right)
        if all(cell_clean(cell, label, stream.L) for label in labels)
    )


def torsion_pigeonhole(clean_cells, d):
    """Return the least-r density certificate and a concrete torsion pair."""
    cell_for_residue = {}
    for cell in clean_cells:
        cell_for_residue.setdefault(cell % d, cell)
    residues = tuple(sorted(cell_for_residue))
    for qualifying_order in range(2, 8):
        if d % qualifying_order or len(residues) <= d // qualifying_order:
            continue
        quotient = d // qualifying_order
        cosets = defaultdict(list)
        for residue in residues:
            cosets[residue % quotient].append(residue)
        crowded = next(rows for rows in cosets.values() if len(rows) >= 2)
        first_residue, second_residue = crowded[:2]
        shift = (second_residue - first_residue) % d
        effective_order = d // gcd(d, shift)
        require(
            2 <= effective_order <= qualifying_order <= 7,
            (d, qualifying_order, shift, effective_order),
        )
        phase_numerators = tuple(
            min(unit, effective_order - unit)
            for unit in range(1, effective_order)
            if gcd(unit, effective_order) == 1
        )
        require(
            phase_numerators
            and 7 * min(phase_numerators) >= effective_order,
            (d, qualifying_order, effective_order, phase_numerators),
        )
        first_cell = cell_for_residue[first_residue]
        second_cell = cell_for_residue[second_residue]
        require((second_cell - first_cell) % d == shift, "cell shift changed")
        return (
            qualifying_order,
            effective_order,
            len(residues),
            quotient,
            first_cell,
            second_cell,
            first_residue,
            second_residue,
            shift,
            min(phase_numerators),
            len(clean_cells),
        )
    return (None, None, len(residues), None, None, None, None, None, None, None, len(clean_cells))


def analyze_residual_body(body, residuals):
    stream = ray.Stream(body)
    require(stream.first == FIRST < stream.high_floor, (body, "high gate inactive"))
    needed = {d for ds in residuals for d in suffix_slots(ds, stream.first_d)}
    low, high, low_signs, high_checks = build_literal_tables(stream, needed)
    require(dict(low_signs).get(-1, 0) > 0, (body, "negative lows were lost"))
    gap, gap_witness = duplicate_two_high_gap(stream, residuals, low, high)
    require(gap == EXPECTED_GAPS[body] > 0, (body, gap))
    zero_rows = zero_high_scalar_passes(stream, residuals, low)
    cases = one_high_cases(stream, residuals, low, high)
    distinct_low_pairs = {
        tuple(sorted(label for _d, label in low_rows))
        for _ds, _high_d, low_rows, _excess in cases
    }
    expected_zero, expected_cases, expected_pairs = EXPECTED_TERMINAL[body]
    require(
        (len(zero_rows), len(cases), len(distinct_low_pairs))
        == (expected_zero, expected_cases, expected_pairs),
        (body, len(zero_rows), len(cases), len(distinct_low_pairs)),
    )
    clean_cache = {}
    witness_cache = {}
    rows = []
    qualifying_histogram = Counter()
    effective_histogram = Counter()
    for ds, high_d, low_rows, scalar_excess in cases:
        labels = tuple(sorted(label for _d, label in low_rows))
        if labels not in clean_cache:
            clean_cache[labels] = fixed_safe_cells(stream, labels)
        key = (labels, high_d)
        if key not in witness_cache:
            witness_cache[key] = torsion_pigeonhole(clean_cache[labels], high_d)
        witness = witness_cache[key]
        require(witness[0] is not None, (body, ds, high_d, labels, witness))
        qualifying_histogram[witness[0]] += 1
        effective_histogram[witness[1]] += 1
        require(
            all(
                cell_clean(cell, label, stream.L)
                for cell in (witness[4], witness[5])
                for label in (stream.first, *labels)
            ),
            (body, ds, "fixed comb entered witness cell"),
        )
        rows.append((ds, high_d, low_rows, scalar_excess, witness))
    return (
        body,
        stream.L,
        stream.high_floor,
        stream.lower - stream.first_delta,
        len(residuals),
        gap,
        gap_witness,
        len(zero_rows),
        len(cases),
        len(distinct_low_pairs),
        low_signs,
        high_checks,
        tuple(sorted(qualifying_histogram.items())),
        tuple(sorted(effective_histogram.items())),
        tuple(rows),
    )


def ftext(value):
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


def render(records, terminal_records, frontier_rows, atlas_counts, next_height, next_count):
    totals = tuple(sum(len(row[index]) for row in records) for index in (7, 8, 9, 10))
    require(totals == (1172, 271, 830, 71), totals)
    residual_bodies = tuple(row[0] for row in records if row[10])
    require(
        residual_bodies
        == (
            (1, 4, 10, 11, 12, 14),
            (1, 8, 10, 11, 12, 14),
            (2, 8, 10, 11, 12, 14),
        ),
        residual_bodies,
    )
    require(tuple(row[0] for row in terminal_records) == residual_bodies, "terminal order")
    total_zero = sum(row[7] for row in terminal_records)
    total_cases = sum(row[8] for row in terminal_records)
    total_pairs = sum(row[9] for row in terminal_records)
    require((total_zero, total_cases, total_pairs) == (70, 73, 7), (total_zero, total_cases, total_pairs))
    qualifying_histogram = Counter()
    effective_histogram = Counter()
    for row in terminal_records:
        qualifying_histogram.update(dict(row[12]))
        effective_histogram.update(dict(row[13]))
    require(
        tuple(sorted(qualifying_histogram.items())) == EXPECTED_QUALIFYING_HISTOGRAM,
        qualifying_histogram,
    )
    require(sum(qualifying_histogram.values()) == total_cases, qualifying_histogram)
    require(THREE_ALIGNED_CAP < 1, THREE_ALIGNED_CAP)
    require(len(frontier_rows) == 7 and next_height == 294 and next_count == 1, (frontier_rows, next_height, next_count))
    require(INHERITED_LEDGER - len(frontier_rows) == FINAL_LEDGER, "ledger arithmetic")
    representative_instance = next(
        row[14] for row in records if row[14] is not None
    )
    semantic_payload = (
        FIRST,
        PROJECTED_WALL,
        THREE_ALIGNED_CAP,
        records,
        terminal_records,
        frontier_rows,
        atlas_counts,
        next_height,
        next_count,
        INHERITED_LEDGER,
        FINAL_LEDGER,
        EXPECTED_UPSTREAM_OUTPUT_SHA256,
        representative_instance,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=3 z1=297 ray/status and located-torsion closure",
        f"ray_source_sha256={file_sha256(RAY_PATH)}",
        f"atlas_source_sha256={file_sha256(ATLAS_SOURCE)}",
        f"atlas_output_sha256={file_sha256(ATLAS)}",
        f"upstream_z298_source_sha256={file_sha256(UPSTREAM)}",
        f"upstream_z298_output_sha256={file_sha256(UPSTREAM_OUTPUT)}",
        "scope=seven exact z297 scalar-atlas body rows;three distinct later nonaligned labels;no finite label horizon",
        "projected_high_gate=each pinned atlas row is HIGH-TAIL;first=297<floor(13L/132)+1;therefore at least one later label is high",
        "ray_law=(z+L)delta(z+L)=zdelta(z);all arbitrarily large high labels reduced to finitely many unit directions",
        "low_policy=every literal label 298<=z<floor(13L/132)+1 is retained, including negative-amplitude labels",
        f"frontier_totals=states:{totals[0]};crude_kills:{totals[1]};status_kills:{totals[2]};residual_states:{totals[3]}",
        "independent_exact_farkas_checks=830/830:PASS",
        f"residual_bodies={residual_bodies}",
        "logical_split=zero-high excluded only by inherited projected gate;duplicate-permitting two-or-more-high upper has positive exact gap on all three bodies;therefore exactly one high",
        f"zero_high_hostile_control={total_zero}/71 denominator states pass scalar arithmetic without a high label",
        f"finite_reduction=states:71;one_high_cases:{total_cases};body-distinct_low_pairs:{total_pairs}",
        "torsion_pigeonhole=S is a set of distinct mod-d residues represented by complete fixed-safe body cells;least divisor r|d with 2<=r<=7 and |S|>d/r forces two residues in one mod-(d/r) coset;effective phase order<=r",
        f"least_qualifying_r_histogram={dict(sorted(qualifying_histogram.items()))};effective_order_histogram={dict(sorted(effective_histogram.items()))}",
        "phase_separation=high z=(L/d)u+hL with gcd(u,d)=1;height cancels and multiplication by u preserves effective order;the nonzero circular gap is >=1/r>=1/7,so two strict-open radius-1/14 danger conditions cannot coexist (at r=7 only excluded endpoints meet)",
        f"projected_consequence=all {total_cases} cases have full projected drift-safe mass 1>three-aligned union cap {ftext(THREE_ALIGNED_CAP)}",
        f"representative_infeasible_instance={representative_instance};exact_farkas=VERIFIED;solver_basis_not_frozen",
    ]
    for row in records:
        (
            body, L, high, first_d, trials, checks, signs, states, crude,
            status, packets, mhist, instance_digest, verified_count,
            _representative_instance,
        ) = row
        require(verified_count == len(status), (body, verified_count, len(status)))
        lines.append(
            f"BODY;E={body};L={L};high={high};d1={first_d};trials={trials};checks={checks};"
            f"signs={dict(signs)};states={len(states)};crude={len(crude)};status={len(status)};"
            f"residual={len(packets)};M={dict(mhist)};"
            f"verified_farkas_instance_sha256={instance_digest};"
            f"verified_farkas_checks={verified_count};all_negative=1;"
            "solver_basis_not_frozen;contradiction_magnitudes_not_frozen"
        )
        for packet in packets:
            lines.append(f"RESIDUAL;E={body};row={packet}")
    for row in terminal_records:
        (
            body, L, high, required, residual_count, gap, gap_witness,
            zero_count, case_count, pair_count, low_signs, high_checks,
            qualifying, effective, cases,
        ) = row
        lines.append(
            f"TERMINAL;E={body};L={L};high={high};required_suffix={ftext(required)};"
            f"residual_states={residual_count};two_high_gap={ftext(gap)};"
            f"gap_witness={gap_witness};zero_high_passes={zero_count};"
            f"one_high_cases={case_count};low_pairs={pair_count};"
            f"low_sign_census={dict(low_signs)};high_ray_checks={high_checks};"
            f"least_r={dict(qualifying)};effective_orders={dict(effective)}"
        )
        for case in cases:
            lines.append(f"CASE;E={body};row={case}")
    lines.extend(
        (
            "atlas_exact_check=z297 rows:7;z295..296 empty;next occupied z294 with1 body",
            f"ledger_decrement={INHERITED_LEDGER}-7={FINAL_LEDGER};decrement counts body rows,not 71 quotient states",
            "conclusion=all seven projected k3 z1=297 body rows are empty;with pinned z298 closure,projected k3 cap<=294",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=7)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    bodies, frontier_rows, atlas_counts, next_height, next_count = atlas_data()
    if args.processes == 1:
        records = tuple(evaluate(body) for body in bodies)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(bodies))) as pool:
            records = tuple(sorted(pool.map(evaluate, bodies)))
    residual_ledger = {
        row[0]: tuple(packet[0] for packet in row[10])
        for row in records
        if row[10]
    }
    terminal_records = tuple(
        analyze_residual_body(body, residuals)
        for body, residuals in residual_ledger.items()
    )
    payload = render(
        records,
        terminal_records,
        frontier_rows,
        atlas_counts,
        next_height,
        next_count,
    )
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
