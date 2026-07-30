#!/usr/bin/env python3
"""Close every projected ``k=3,z1=312`` literal packet exactly.

The exact ray/status frontier leaves 70 denominator states on four bodies.
This verifier partitions every later label by its literal low/high status and
enumerates *all signs* of the finite low rays.  Every scalar-viable pattern
has exactly one high label.  The bounded branch is checked by the lossless
projected-cell identity.  On an unbounded branch the two low labels already
clear the scalar wall; two body-safe cells missed by the fixed packet move
the arbitrary H phase by a half-turn, except for denominators 1911 and 9555,
where a third-turn does the same job.  The argument covers every height and
both amplitude signs on the H ray.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FRONTIER_PATH = ROOT / "04-computation" / "lrc14_j7_k3_z312_ray_status_frontier_thm2941.py"
FRONTIER_OUTPUT_PATH = ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z312_ray_status_frontier_thm2941.out"
PROJECTION_PATH = ROOT / "04-computation" / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
DEFAULT_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z312_torsion_h_drift_closure_thm2941.out"
EXPECTED_FRONTIER_SHA256 = "202588c1232fc507fa376ae57dbf630d340c54e5b6be727e6e5be5b83cdaa20d"
EXPECTED_FRONTIER_OUTPUT_SHA256 = "75c1f42d86654efbb719c96d882e934ecf64b486c3325ad885e989767745223f"
EXPECTED_PROJECTION_SHA256 = "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
EXPECTED_SEMANTIC_SHA256 = "29b61bc296b7dfd63ef7db685cb446b684a54ae9af10ab873b55b7539395f211"
FIRST = 312
CAP = F(36, 91)
EXPECTED_COUNTS = {
    "states": 70,
    "mask_tests": 490,
    "viable_patterns": 76,
    "finite_packets": 27491,
    "infinite_bases": 73,
}
EXPECTED_TORSION_ORDERS = {2: 71, 3: 2}
EXPECTED_THIRD_TURN_SHIFTS = {1911: 637, 9555: 3185}
EXPECTED_NEGATIVE_LOW_CANDIDATES = 15_528
EXPECTED_FINITE_MINIMUM = (
    F(2255, 2_786_147),
    (1, 8, 10, 11, 12, 14),
    (1617, 1848, 5390, 129360),
    (312, 350, 400, 30617),
    F(12137, 30617),
    1,
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if EXPECTED_FRONTIER_SHA256 is not None:
    require(file_sha256(FRONTIER_PATH) == EXPECTED_FRONTIER_SHA256, "frontier changed")
if EXPECTED_FRONTIER_OUTPUT_SHA256 is not None:
    require(file_sha256(FRONTIER_OUTPUT_PATH) == EXPECTED_FRONTIER_OUTPUT_SHA256, "frontier transcript changed")
require(file_sha256(PROJECTION_PATH) == EXPECTED_PROJECTION_SHA256, "projection engine changed")
frontier = load("z312_exact_frontier", FRONTIER_PATH)
P = load("z312_lossless_projection", PROJECTION_PATH)
ray = frontier.ray


def residuals_from_frontier():
    current = None
    residuals = {}
    for line in FRONTIER_OUTPUT_PATH.read_text().splitlines():
        if line.startswith("E="):
            body_text = line.split(";", 1)[0][2:]
            current = tuple(ast.literal_eval(body_text))
        elif line.startswith("  residual_denominators="):
            rows = tuple(ast.literal_eval(line.split("=", 1)[1]))
            require(current is not None and rows, "orphan residual ledger")
            residuals[current] = rows
    require(tuple(map(len, residuals.values())) == (25, 41, 2, 2), ("residual ledger changed", residuals))
    return residuals


def first_on_ray(residue, L, threshold):
    if residue >= threshold:
        return residue
    return residue + ((threshold - residue + L - 1) // L) * L


def denominator(L, label):
    return L // gcd(L, label)


def ray_tables(stream, denominators):
    require(
        stream.first < stream.high_floor < stream.L,
        ("low/high partition changed", stream.body, stream.high_floor, stream.L),
    )
    low = defaultdict(list)
    high = defaultdict(list)
    sign_counts = Counter()
    direction_counts = {}
    for d in sorted(denominators):
        step = stream.L // d
        amplitudes = {}
        for unit in range(1, d):
            if gcd(unit, d) != 1:
                continue
            residue = step * unit
            z = first_on_ray(residue, stream.L, FIRST + 1)
            amplitude = z * ray.local.delta(stream.carrier, stream.h, z)
            next_z = z + stream.L
            require(
                next_z * ray.local.delta(stream.carrier, stream.h, next_z) == amplitude,
                ("ray recurrence failed", stream.body, d, unit),
            )
            amplitudes[unit] = amplitude
            sign_counts[(d, (amplitude > 0) - (amplitude < 0))] += 1
            if z < stream.high_floor:
                low[d].append((amplitude / z, z, amplitude, unit))
            high_z = first_on_ray(residue, stream.L, stream.high_floor)
            high[d].append((amplitude / high_z, high_z, amplitude, unit))
        require(amplitudes, ("empty denominator class", d))
        for unit, amplitude in amplitudes.items():
            require(
                amplitude + amplitudes[(-unit) % d] == 0,
                ("antipodal amplitude failed", stream.body, d, unit),
            )
        direction_counts[d] = len(amplitudes)
        low[d].sort(key=lambda row: (row[0], -row[1]), reverse=True)
        high[d].sort(key=lambda row: (row[0], -row[1]), reverse=True)
        require(
            high[d][0][0] >= 0,
            ("high-ray maximum became a remote negative supremum", stream.body, d),
        )
    return dict(low), dict(high), tuple(sorted(sign_counts.items())), tuple(sorted(direction_counts.items()))


def max_pattern(ds, high_mask, low, high):
    groups = Counter((d, bool(high_mask[index])) for index, d in enumerate(ds))
    total = F(0)
    labels = []
    for (d, is_high), multiplicity in groups.items():
        rows = (high if is_high else low).get(d, ())
        if len(rows) < multiplicity:
            return None
        total += sum((row[0] for row in rows[:multiplicity]), F(0))
        labels.extend(row[1] for row in rows[:multiplicity])
    require(len(labels) == len(set(labels)), "pattern maximum repeats a label")
    return total, tuple(labels)


def grouped_choices(ds, table):
    grouped = []
    for d, multiplicity in sorted(Counter(ds).items()):
        grouped.append(tuple(combinations(table.get(d, ()), multiplicity)))
    for selected_groups in product(*grouped):
        selected = tuple(row for group in selected_groups for row in group)
        labels = tuple(row[1] for row in selected)
        require(len(labels) == len(set(labels)), "low choice repeats a label")
        yield selected


def high_labels_meeting(rows, deficit, L):
    require(deficit > 0, "finite H deficit is not positive")
    for _value, head, amplitude, unit in rows:
        if amplitude <= 0:
            continue
        last_floor = (amplitude / deficit).numerator // (amplitude / deficit).denominator
        if last_floor < head:
            continue
        for level in range((last_floor - head) // L + 1):
            label = head + level * L
            value = amplitude / label
            if value >= deficit:
                yield value, label, amplitude, unit


@lru_cache(maxsize=None)
def phase_danger(cell, label, L):
    return P.phase_danger(cell, label, L)


def spread_cells(cells):
    selected = []
    seen = set()
    count = len(cells)
    for power in range(9):
        denominator_grid = 1 << power
        for odd in range(1, 2 * denominator_grid, 2):
            index = min(count - 1, odd * count // (2 * denominator_grid))
            cell = cells[index]
            if cell not in seen:
                seen.add(cell)
                selected.append(cell)
    selected.extend(cell for cell in cells if cell not in seen)
    require(len(selected) == len(cells), "spread-cell permutation changed size")
    return tuple(selected)


def projected_lower(cells, L, labels, *, stop=True):
    common = ((F(0), F(1)),)
    for used, cell in enumerate(cells, 1):
        local = P.merge_fraction(
            [
                interval
                for label in labels
                for interval in phase_danger(cell, label, L)
            ]
        )
        common = P.intersect_fraction(common, local)
        if stop and 1 - P.interval_mass(common) > CAP:
            return 1 - P.interval_mass(common), used, common
    return 1 - P.interval_mass(common), len(cells), common


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


@lru_cache(maxsize=None)
def torsion_data(d):
    units = tuple(unit for unit in range(1, d) if gcd(unit, d) == 1)
    for order in range(2, 7):
        if d % order:
            continue
        shift = d // order
        phases = tuple((unit * shift) % d for unit in units)
        if all(7 * min(phase, d - phase) > d for phase in phases):
            return order, shift, tuple(sorted(set(phases))), len(units)
    return None


def torsion_pair(cells, L, fixed, d, safe_cache):
    data = torsion_data(d)
    require(data is not None, ("no small torsion separation", d))
    order, shift, phases, unit_count = data
    for label in fixed:
        if label not in safe_cache:
            safe_cache[label] = frozenset(
                cell for cell in cells if not phase_danger(cell, label, L)
            )
    safe = set.intersection(*(set(safe_cache[label]) for label in fixed))
    by_residue = {}
    for cell in sorted(safe):
        by_residue.setdefault(cell % d, cell)
    for residue, first_cell in sorted(by_residue.items()):
        second_cell = by_residue.get((residue + shift) % d)
        if second_cell is not None:
            require((second_cell - first_cell) % d == shift, "torsion cell shift changed")
            return first_cell, second_cell, order, shift, phases, unit_count
    return None


def digest(value):
    return hashlib.sha256(repr(value).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    residuals = residuals_from_frontier()
    frontier_keys = {(body, ds) for body, rows in residuals.items() for ds in rows}
    closed_keys = set()
    counts = Counter()
    order_histogram = Counter()
    order_denominators = defaultdict(set)
    torsion_formulae = set()
    negative_low_candidates = negative_low_admissible = 0
    finite_digest = hashlib.sha256()
    torsion_digest = hashlib.sha256()
    global_minimum = None
    global_max_cells = 0
    torsion_controls = {}
    body_rows = []
    for body, states in residuals.items():
        stream = ray.Stream(body)
        tails = []
        for ds in states:
            tail = list(ds)
            tail.remove(stream.first_d)
            tail = tuple(sorted(tail))
            require(len(tail) == len(set(tail)) == 3, ("repeated tail denominator", body, ds))
            tails.append(tail)
        denominators = {d for tail in tails for d in tail}
        low, high, signs, directions = ray_tables(stream, denominators)
        negative_low_candidates += sum(
            1 for rows in low.values() for row in rows if row[2] < 0
        )
        raw_cells = P.body_cells(P.A.carrier_for(body), stream.L)
        cells = spread_cells(raw_cells)
        safe_cache = {}
        body_finite = body_tails = 0
        for ds, tail in zip(states, tails):
            required = stream.lower - stream.first_delta
            viable = []
            counts["mask_tests"] += 7
            for mask_integer in range(1, 8):
                mask = tuple(bool(mask_integer & (1 << index)) for index in range(3))
                maximum = max_pattern(tail, mask, low, high)
                if maximum is not None and maximum[0] >= required:
                    viable.append((mask, maximum))
            require(viable and all(sum(mask) == 1 for mask, _maximum in viable), ("non-single-H pattern survived", body, ds, viable))
            patterns = tuple(
                sorted(
                    {
                        (tail[index], tuple(tail[j] for j in range(3) if j != index))
                        for mask, _maximum in viable
                        for index in range(3)
                        if mask[index]
                    }
                )
            )
            counts["viable_patterns"] += len(patterns)
            state_packets = set()
            state_tails = set()
            for high_d, low_ds in patterns:
                for selected in grouped_choices(low_ds, low):
                    base = sum((row[0] for row in selected), F(0))
                    low_labels = tuple(row[1] for row in selected)
                    fixed = (FIRST, *sorted(low_labels))
                    if base >= required:
                        key = (high_d, fixed)
                        if key in state_tails:
                            continue
                        state_tails.add(key)
                        witness = torsion_pair(raw_cells, stream.L, fixed, high_d, safe_cache)
                        require(witness is not None, ("unclosed H tail", body, ds, high_d, fixed))
                        j, k, order, shift, phases, unit_count = witness
                        representative = high[high_d][0][1]
                        require(
                            representative >= stream.high_floor
                            and denominator(stream.L, representative) == high_d
                            and len(set((*fixed, representative))) == 4,
                            ("bad H representative", body, ds, representative),
                        )
                        packet = (*fixed, representative)
                        local_j = P.merge_fraction(
                            [interval for label in packet for interval in phase_danger(j, label, stream.L)]
                        )
                        local_k = P.merge_fraction(
                            [interval for label in packet for interval in phase_danger(k, label, stream.L)]
                        )
                        require(P.intersect_fraction(local_j, local_k) == (), "torsion sections overlap")
                        tail_row = (
                            body,
                            ds,
                            high_d,
                            fixed,
                            base - required,
                            j,
                            k,
                            order,
                            shift,
                            phases,
                            unit_count,
                            representative,
                        )
                        torsion_digest.update(f"{tail_row}\n".encode())
                        order_histogram[order] += 1
                        order_denominators[order].add(high_d)
                        torsion_formulae.add((high_d, order, shift, phases, unit_count))
                        counts["infinite_bases"] += 1
                        body_tails += 1
                        torsion_controls.setdefault(order, tail_row)
                        if any(row[2] < 0 for row in selected):
                            negative_low_admissible += 1
                        continue
                    deficit = required - base
                    if deficit > high[high_d][0][0]:
                        continue
                    for hrow in high_labels_meeting(high[high_d], deficit, stream.L):
                        labels = (FIRST, *sorted((*low_labels, hrow[1])))
                        if len(labels) != len(set(labels)) or labels in state_packets:
                            continue
                        require(
                            tuple(sorted(denominator(stream.L, label) for label in labels)) == ds,
                            ("literal denominator mismatch", body, ds, labels),
                        )
                        total = stream.first_delta + base + hrow[0]
                        require(total >= stream.lower and max(labels[1:]) >= stream.high_floor, "literal scalar gate failed")
                        mass, used, _common = projected_lower(cells, stream.L, labels)
                        require(mass > CAP, ("finite packet survives projection", body, ds, labels, mass))
                        state_packets.add(labels)
                        packet_row = (body, ds, labels, total - stream.lower, mass, used)
                        finite_digest.update(f"{packet_row}\n".encode())
                        candidate = (mass - CAP, body, ds, labels, mass, used)
                        if global_minimum is None or candidate < global_minimum:
                            global_minimum = candidate
                        global_max_cells = max(global_max_cells, used)
                        counts["finite_packets"] += 1
                        body_finite += 1
                        if any(row[2] < 0 for row in selected):
                            negative_low_admissible += 1
            require(state_packets or state_tails, ("empty literal state", body, ds))
            closed_keys.add((body, ds))
        body_rows.append((body, len(states), body_finite, body_tails, signs, directions, len(raw_cells)))
    counts["states"] = len(closed_keys)
    require(dict(counts) == EXPECTED_COUNTS, ("literal closure counts changed", counts))
    require(dict(order_histogram) == EXPECTED_TORSION_ORDERS, ("torsion orders changed", order_histogram))
    require(order_denominators[3] == {1911, 9555}, ("third-turn denominators changed", order_denominators))
    third_turn_shifts = {
        d: shift
        for d, order, shift, phases, _unit_count in torsion_formulae
        if order == 3
        and phases == (shift, 2 * shift)
    }
    require(
        third_turn_shifts == EXPECTED_THIRD_TURN_SHIFTS,
        ("third-turn phase formulas changed", third_turn_shifts),
    )
    require(
        (negative_low_candidates, negative_low_admissible)
        == (EXPECTED_NEGATIVE_LOW_CANDIDATES, 0),
        ("negative-low audit changed", negative_low_candidates, negative_low_admissible),
    )
    require(closed_keys == frontier_keys and len(frontier_keys) == 70, "closure union is not the exact frontier")
    require(
        global_minimum == EXPECTED_FINITE_MINIMUM
        and global_minimum[0] > 0
        and global_max_cells == 3,
        ("finite strict-margin ledger changed", global_minimum, global_max_cells),
    )
    _margin, control_body, _control_ds, control_labels, _prefix_mass, _used = global_minimum
    control_stream = ray.Stream(control_body)
    control_cells = P.body_cells(P.A.carrier_for(control_body), control_stream.L)
    full_mass, full_used, _common = projected_lower(control_cells, control_stream.L, control_labels, stop=False)
    direct_mass = direct_projected_mass(control_body, control_stream.L, control_labels)
    require(full_used == len(control_cells) and full_mass == direct_mass, "finite direct control failed")
    torsion_direct_controls = {}
    for order, row in sorted(torsion_controls.items()):
        body, _ds, _high_d, fixed, _surplus, _j, _k, _order, _shift, _phases, _units, representative = row
        stream = ray.Stream(body)
        labels = (*fixed, representative)
        cells = P.body_cells(P.A.carrier_for(body), stream.L)
        mass, _used, common = projected_lower(cells, stream.L, labels, stop=False)
        direct = direct_projected_mass(body, stream.L, labels)
        require(common == () and mass == direct == 1, ("torsion direct control failed", order, mass, direct))
        torsion_direct_controls[order] = (body, labels, mass, direct)
    semantic_payload = (
        tuple(body_rows),
        tuple(sorted(counts.items())),
        tuple(sorted(order_histogram.items())),
        tuple((order, tuple(sorted(ds))) for order, ds in sorted(order_denominators.items())),
        tuple(sorted(torsion_formulae)),
        negative_low_candidates,
        negative_low_admissible,
        finite_digest.hexdigest(),
        torsion_digest.hexdigest(),
        global_minimum,
        global_max_cells,
        full_mass,
        direct_mass,
        tuple(sorted(torsion_direct_controls.items())),
        tuple(sorted(closed_keys)),
    )
    semantic_hash = digest(semantic_payload)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=3 z1=312 finite-packet and torsion H-drift closure",
        f"frontier_source_sha256={file_sha256(FRONTIER_PATH)}",
        f"frontier_output_sha256={file_sha256(FRONTIER_OUTPUT_PATH)}",
        f"projection_source_sha256={file_sha256(PROJECTION_PATH)}",
        "scope=all 70 exact frontier residual states;all low-ray amplitude signs",
        "literal_partition=490 nonempty low/high masks;76 viable patterns;every viable pattern has exactly one high label",
        f"finite_packets={counts['finite_packets']};finite_kills={counts['finite_packets']};max_cells_used={global_max_cells}",
        f"finite_strict_minimum={global_minimum}",
        f"finite_packet_sha256={finite_digest.hexdigest()}",
        f"negative_low_candidates={negative_low_candidates};negative_low_admissible={negative_low_admissible}",
        f"unbounded_bases={counts['infinite_bases']};torsion_orders={dict(order_histogram)}",
        f"torsion_denominators={dict((order, tuple(sorted(ds))) for order, ds in sorted(order_denominators.items()))}",
        f"torsion_formulae={tuple(sorted(torsion_formulae))}",
        "H_phase_law=z=(L/d)u+mL => phase_k-phase_j=u(k-j)/d mod 1",
        "half_turn_scope=71 bases;third_turn_scope=d=1911,9555;all H heights and signs",
        f"torsion_witness_sha256={torsion_digest.hexdigest()}",
        f"finite_direct_control=full_cell_mass:{full_mass};direct_subtraction_mass:{direct_mass}",
        f"torsion_direct_controls={torsion_direct_controls}",
        "closure_union=70/70 exact frontier residual states",
        "next_lossless_scalar_atlas_height=306",
        "consequence=z1=312 slice empty;with the pinned lossless atlas, projected k=3 cap <=306",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
