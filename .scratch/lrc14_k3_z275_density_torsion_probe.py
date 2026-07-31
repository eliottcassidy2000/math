#!/usr/bin/env python3
"""Generic exact residual reduction and torsion-density probe."""

from __future__ import annotations

import importlib.util
import multiprocessing as mp
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_PATH = ROOT / "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
FIRST = 275
SURVIVOR_BODIES = (
    (1, 5, 7, 9, 11, 13),
    (1, 5, 9, 11, 13, 14),
    (1, 8, 10, 11, 12, 14),
    (2, 5, 7, 9, 11, 13),
    (2, 5, 9, 11, 13, 14),
)


def load_ray(first=FIRST):
    spec = importlib.util.spec_from_file_location(f"ray{first}_terminal", RAY_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(RAY_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.FIRST = first
    return module


def delta(ray, stream, label):
    return ray.suffix.A.singleton_coverage(stream.carrier, label) - stream.h / 7


def first_on_ray(residue, modulus, threshold):
    if residue >= threshold:
        return residue
    return residue + ((threshold - residue + modulus - 1) // modulus) * modulus


def suffix_slots(ds, first_d):
    slots = list(ds)
    slots.remove(first_d)
    if len(slots) != 3:
        raise RuntimeError((ds, first_d, slots))
    return tuple(slots)


def build_tables(ray, stream, needed):
    low = {d: [] for d in needed}
    low_signs = Counter()
    for label in range(stream.first + 1, stream.high_floor):
        d = stream.L // gcd(stream.L, label)
        if d not in low:
            continue
        value = delta(ray, stream, label)
        low_signs[(d, (value > 0) - (value < 0))] += 1
        low[d].append((value, label))
    for d in low:
        low[d].sort(key=lambda item: (-item[0], item[1]))
        low[d] = tuple(low[d])

    high = {}
    high_units = {}
    recurrence_checks = 0
    for d in needed:
        step = stream.L // d
        best = (Q(0), None, None, None)
        unit_rows = []
        for unit in range(1, d):
            if gcd(unit, d) != 1:
                continue
            residue = step * unit
            amplitude = residue * delta(ray, stream, residue)
            if (residue + stream.L) * delta(ray, stream, residue + stream.L) != amplitude:
                raise RuntimeError((stream.body, d, unit, "ray recurrence"))
            recurrence_checks += 1
            label = first_on_ray(residue, stream.L, stream.high_floor)
            unit_rows.append((unit, residue, amplitude, label))
            value = amplitude / label if amplitude > 0 else Q(0)
            candidate = (value, label, unit, amplitude)
            if candidate[0] > best[0] or (
                candidate[0] == best[0]
                and candidate[1] is not None
                and (best[1] is None or candidate[1] < best[1])
            ):
                best = candidate
        high[d] = best
        high_units[d] = tuple(unit_rows)
    return low, high, high_units, tuple(sorted(low_signs.items())), recurrence_checks


def duplicate_two_high_gap(stream, residuals, low, high):
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
    if best is None:
        raise RuntimeError((stream.body, "missing two-high upper"))
    return required - best, witness


def zero_high_passes(stream, residuals, low):
    required = stream.lower - stream.first_delta
    passes = []
    for ds in residuals:
        slots = suffix_slots(ds, stream.first_d)
        if not all(low[d] for d in slots):
            continue
        value = sum((low[d][0][0] for d in slots), Q())
        if value >= required:
            passes.append((ds, tuple(low[d][0][1] for d in slots), value - required))
    return tuple(passes)


def finite_pairs(first, second, same, threshold):
    rows = []
    if same:
        for index, left in enumerate(first):
            if index + 1 >= len(first) or left[0] + first[index + 1][0] < threshold:
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
            low_ds = tuple(d for index, d in enumerate(slots) if index != high_index)
            if not low[low_ds[0]] or not low[low_ds[1]]:
                continue
            threshold = required - high[high_d][0]
            pairs = finite_pairs(
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
    return 14 * residue >= L and 14 * (residue + label) <= 13 * L


def clean_cells(stream, low_labels):
    labels = (stream.first, *low_labels)
    return tuple(
        cell
        for left, right in stream.ranges
        for cell in range(left, right)
        if all(cell_clean(cell, label, stream.L) for label in labels)
    )


def density_certificate(clean, d):
    residues = sorted({cell % d for cell in clean})
    for r in range(2, 8):
        if d % r:
            continue
        quotient = d // r
        buckets = {}
        for residue in residues:
            key = residue % quotient
            if key in buckets:
                first = buckets[key]
                difference = (residue - first) % d
                effective = d // gcd(d, difference)
                if not (2 <= effective <= r <= 7):
                    raise RuntimeError((d, r, first, residue, effective))
                return (
                    r,
                    effective,
                    len(residues),
                    quotient,
                    first,
                    residue,
                    difference,
                    len(clean),
                )
            buckets[key] = residue
        if len(residues) > quotient:
            raise RuntimeError((d, r, len(residues), quotient, "pigeonhole"))
    return (None, None, len(residues), None, None, None, None, len(clean))


def forced_density_certificate(clean, d):
    """Return the least torsion order forced by residue cardinality alone."""
    residues = sorted({cell % d for cell in clean})
    for r in range(2, 8):
        if d % r or len(residues) <= d // r:
            continue
        quotient = d // r
        buckets = {}
        for residue in residues:
            key = residue % quotient
            if key in buckets:
                first = buckets[key]
                difference = (residue - first) % d
                effective = d // gcd(d, difference)
                if not (2 <= effective <= r <= 7):
                    raise RuntimeError((d, r, first, residue, effective))
                return (
                    r,
                    effective,
                    len(residues),
                    quotient,
                    first,
                    residue,
                    difference,
                    len(clean),
                )
            buckets[key] = residue
        raise RuntimeError((d, r, len(residues), quotient, "pigeonhole"))
    return (None, None, len(residues), None, None, None, None, len(clean))


def ray_resolved_phase_certificate(clean, d, unit):
    """Locate a fixed-safe cell whose high-ray phase misses the open comb."""
    for cell in clean:
        phase = (unit * (cell % d)) % d
        if 14 * min(phase, d - phase) >= d:
            return (cell, cell % d, phase)
    return None


def ray_resolved_gate(ray, stream, cases, clean_cache, high_units):
    """Retain the primitive ray before taking its scalar envelope.

    On a primitive unit ray the high contribution is ``A/z`` and its phase
    on a fixed cell is independent of the height.  Positive ``A`` is maximal
    at the first high point, zero ``A`` is attained, and negative ``A`` tends
    monotonically to zero from below.  These three signs give an exact test
    for whether *some* height on the ray can meet the scalar wall.
    """
    required = stream.lower - stream.first_delta
    active = 0
    safe = 0
    failures = []
    sign_histogram = Counter()
    result_cache = {}
    for ds, high_d, low_rows, _excess in cases:
        labels = tuple(sorted(label for _d, label in low_rows))
        clean = clean_cache[labels]
        key = (labels, high_d)
        if key in result_cache:
            key_active, key_safe, key_signs, key_failures = result_cache[key]
            active += key_active
            safe += key_safe
            sign_histogram.update(dict(key_signs))
            failures.extend(
                (ds, *failure[1:]) for failure in key_failures
            )
            continue
        low_sum = sum((delta(ray, stream, label) for label in labels), Q())
        slack = low_sum - required
        key_active = 0
        key_safe = 0
        key_signs = Counter()
        key_failures = []
        for unit, _residue, amplitude, first_high in high_units[high_d]:
            sign = (amplitude > 0) - (amplitude < 0)
            if amplitude > 0:
                can_reach = slack + amplitude / first_high >= 0
            elif amplitude == 0:
                can_reach = slack >= 0
            else:
                can_reach = slack > 0
            if not can_reach:
                continue
            key_active += 1
            key_signs[sign] += 1
            certificate = ray_resolved_phase_certificate(clean, high_d, unit)
            if certificate is None:
                key_failures.append(
                    (ds, high_d, labels, unit, amplitude, first_high, slack)
                )
            else:
                key_safe += 1
        result_cache[key] = (
            key_active,
            key_safe,
            tuple(sorted(key_signs.items())),
            tuple(key_failures),
        )
        active += key_active
        safe += key_safe
        sign_histogram.update(key_signs)
        failures.extend(key_failures)
    return active, safe, tuple(sorted(sign_histogram.items())), tuple(failures)


def evaluate(task):
    if (
        isinstance(task, tuple)
        and len(task) in (2, 3)
        and isinstance(task[0], int)
        and isinstance(task[1], tuple)
    ):
        first, body = task[:2]
        run_ray_resolved = task[2] if len(task) == 3 else True
    else:
        first, body = FIRST, task
        run_ray_resolved = True
    ray = load_ray(first)
    stream = ray.Stream(body)
    trials, states, checks, _signs = ray.ray_quotient_states(stream)
    crude, status, survivors = ray.common_status_screen(stream, states)
    residuals = tuple(sorted(survivors))
    if not residuals:
        return (
            body,
            stream.L,
            stream.high_floor,
            trials,
            checks,
            len(states),
            len(crude),
            len(status),
            0,
            None,
            None,
            0,
            0,
            0,
            (),
            (),
            (),
            0,
            (),
            (),
            (),
            0,
            0,
            (),
            (),
        )
    needed = {
        d
        for ds in residuals
        for d in suffix_slots(ds, stream.first_d)
    }
    low, high, high_units, low_signs, recurrence_checks = build_tables(
        ray, stream, needed
    )
    gap, gap_witness = duplicate_two_high_gap(stream, residuals, low, high)
    # If the first label is already at or above the scalar high floor, then
    # every later suffix label is high as well: there are no finite low rows.
    # The one-high gate is indeed unavailable, but a strict all-high upper gap
    # is already a complete closure.  Reject only when that stronger direct
    # closure fails; never fall through to the one-high logic in this regime.
    all_high_terminal = bool(residuals and stream.first >= stream.high_floor)
    if all_high_terminal:
        if any(low[d] for d in needed):
            raise RuntimeError(
                (stream.body, stream.first, stream.high_floor, "unexpected low label")
            )
        if not gap > 0:
            raise RuntimeError(
                (
                    stream.body,
                    stream.first,
                    stream.high_floor,
                    gap,
                    "all-high terminal unresolved",
                )
            )
    zero = zero_high_passes(stream, residuals, low)
    cases = one_high_cases(stream, residuals, low, high)
    if all_high_terminal and (zero or cases):
        raise RuntimeError((stream.body, zero, cases, "one-high logic leaked"))
    clean_cache = {}
    certificate_cache = {}
    forced_cache = {}
    failures = []
    forced_failures = []
    histogram = Counter()
    forced_histogram = Counter()
    forced_effective_histogram = Counter()
    pair_set = set()
    for ds, high_d, low_rows, excess in cases:
        labels = tuple(sorted(label for _d, label in low_rows))
        pair_set.add(labels)
        if labels not in clean_cache:
            clean_cache[labels] = clean_cells(stream, labels)
        key = (labels, high_d)
        if key not in certificate_cache:
            certificate_cache[key] = density_certificate(clean_cache[labels], high_d)
            forced_cache[key] = forced_density_certificate(
                clean_cache[labels], high_d
            )
        certificate = certificate_cache[key]
        forced = forced_cache[key]
        if certificate[0] is None:
            failures.append((ds, high_d, low_rows, excess, certificate))
        else:
            histogram[(certificate[0], certificate[1])] += 1
        if forced[0] is None:
            forced_failures.append((ds, high_d, low_rows, excess, forced))
        else:
            forced_histogram[forced[0]] += 1
            forced_effective_histogram[forced[1]] += 1
    if run_ray_resolved:
        ray_active, ray_safe, ray_signs, ray_failures = ray_resolved_gate(
            ray, stream, cases, clean_cache, high_units
        )
    else:
        ray_active, ray_safe, ray_signs, ray_failures = 0, 0, (), ()
    return (
        body,
        stream.L,
        stream.high_floor,
        trials,
        checks,
        len(states),
        len(crude),
        len(status),
        len(residuals),
        gap,
        gap_witness,
        len(zero),
        len(cases),
        len(pair_set),
        tuple(sorted(histogram.items())),
        tuple(failures),
        low_signs,
        recurrence_checks,
        tuple(sorted(forced_histogram.items())),
        tuple(sorted(forced_effective_histogram.items())),
        tuple(forced_failures),
        ray_active,
        ray_safe,
        ray_signs,
        ray_failures,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=int, default=FIRST)
    parser.add_argument(
        "--body",
        action="append",
        default=[],
        help="comma-separated six-entry body; repeat for every residual body",
    )
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--compact", action="store_true")
    parser.add_argument("--skip-ray-resolved", action="store_true")
    args = parser.parse_args()
    bodies = (
        tuple(tuple(map(int, body.split(","))) for body in args.body)
        if args.body
        else SURVIVOR_BODIES
    )
    tasks = tuple(
        (args.first, body, not args.skip_ray_resolved) for body in bodies
    )
    records = []
    with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
        iterator = pool.imap_unordered(evaluate, tasks, chunksize=1)
        for record in iterator:
            records.append(record)
            if args.compact:
                print(
                    "BODY_SUMMARY",
                    {
                        "body": record[0],
                        "L": record[1],
                        "states": record[5],
                        "crude": record[6],
                        "status": record[7],
                        "residuals": record[8],
                        "two_high_gap": record[9],
                        "zero_high_passes": record[11],
                        "one_high_cases": record[12],
                        "low_pairs": record[13],
                        "histogram": dict(record[14]),
                        "failure_count": len(record[15]),
                        "failure_sample": record[15][:3],
                        "recurrence_checks": record[17],
                        "forced_histogram": dict(record[18]),
                        "forced_effective_histogram": dict(record[19]),
                        "forced_failure_count": len(record[20]),
                        "forced_failure_sample": record[20][:3],
                        "ray_resolved_active": record[21],
                        "ray_resolved_safe": record[22],
                        "ray_resolved_signs": dict(record[23]),
                        "ray_resolved_failure_count": len(record[24]),
                        "ray_resolved_failure_sample": record[24][:3],
                    },
                    flush=True,
                )
            else:
                print("BODY_RECORD", record, flush=True)
    records = tuple(records)
    print(
        "TOTAL",
        {
            "residuals": sum(record[8] for record in records),
            "zero_high_passes": sum(record[11] for record in records),
            "one_high_cases": sum(record[12] for record in records),
            "low_pairs": sum(record[13] for record in records),
            "failures": sum(len(record[15]) for record in records),
            "histogram": dict(
                sum((Counter(dict(record[14])) for record in records), Counter())
            ),
            "recurrence_checks": sum(record[17] for record in records),
            "forced_histogram": dict(
                sum((Counter(dict(record[18])) for record in records), Counter())
            ),
            "forced_effective_histogram": dict(
                sum((Counter(dict(record[19])) for record in records), Counter())
            ),
            "forced_failures": sum(len(record[20]) for record in records),
            "ray_resolved_active": sum(record[21] for record in records),
            "ray_resolved_safe": sum(record[22] for record in records),
            "ray_resolved_failures": sum(len(record[24]) for record in records),
        },
    )


if __name__ == "__main__":
    main()
