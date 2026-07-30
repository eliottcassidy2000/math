#!/usr/bin/env python3
"""Independent exact referee for the z312 finite-low-pair torsion closure.

The inherited projected high-wall gate requires at least one of the three
later labels to lie at or above the body's high floor.  This condition is
essential: scalar arithmetic alone admits a zero-high packet on 69 of the 70
residual denominator states.  The new argument only excludes two-or-more
high labels; together with the inherited gate this leaves exactly one high.

This referee reconstructs the scalar split directly from unit rays, enumerates
all finite low pairs, and uses an integer whole-cell test independent of the
primary projected-cell implementation.  It then checks all order-two and
order-three torsion pairs pointwise.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction as F
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FRONTIER = ROOT / "04-computation/lrc14_j7_k3_z312_ray_status_frontier_thm2941.py"
FRONTIER_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z312_ray_status_frontier_thm2941.out"
DEFAULT_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z312_finite_low_pair_torsion_independent_thm2941.out"
EXPECTED_FRONTIER_SHA256 = "d6a80c15c5c4d8ef2ea8be9fc886c40e70e3189123b5d0b3fce48765fa301977"
EXPECTED_FRONTIER_OUTPUT_SHA256 = "d03fc39ed1f5f64cd2be4e7e28f5cf23e8d7adc0c6a737abc6944bdb7672515f"
EXPECTED_SEMANTIC_SHA256 = "4667f8f0ef18add7466f849c6f902048d436692aa04c1f9037c2a3623709c91e"
EXPECTED_GAPS = {
    (1, 8, 10, 11, 12, 14): F(271403663, 168333225060),
    (1, 8, 10, 12, 13, 14): F(22539649297, 15003760917600),
    (1, 8, 11, 12, 13, 14): F(295936150144, 96567353117235),
    (2, 8, 10, 11, 12, 14): F(41681149, 16799037255),
}
EXPECTED_CASE_COUNTS = {
    (1, 8, 10, 11, 12, 14): 26,
    (1, 8, 10, 12, 13, 14): 50,
    (1, 8, 11, 12, 13, 14): 2,
    (2, 8, 10, 11, 12, 14): 2,
}
EXPECTED_LOW_PAIR_COUNTS = {
    (1, 8, 10, 11, 12, 14): 6,
    (1, 8, 10, 12, 13, 14): 12,
    (1, 8, 11, 12, 13, 14): 1,
    (2, 8, 10, 11, 12, 14): 1,
}
EXPECTED_ZERO_HIGH_SCALAR_PASSES = {
    (1, 8, 10, 11, 12, 14): 24,
    (1, 8, 10, 12, 13, 14): 41,
    (1, 8, 11, 12, 13, 14): 2,
    (2, 8, 10, 11, 12, 14): 2,
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(FRONTIER) == EXPECTED_FRONTIER_SHA256, "frontier source changed")
require(
    file_sha256(FRONTIER_OUTPUT) == EXPECTED_FRONTIER_OUTPUT_SHA256,
    "frontier output changed",
)
frontier = load("z312_independent_frontier", FRONTIER)
ray = frontier.ray
ray.FIRST = frontier.FIRST
A = ray.suffix.A


def first_on_ray(residue, modulus, threshold):
    if residue >= threshold:
        return residue
    return residue + ((threshold - residue + modulus - 1) // modulus) * modulus


def delta(stream, label):
    return A.singleton_coverage(stream.carrier, label) - stream.h / 7


def suffix_denominators(ds, first_d):
    result = list(ds)
    result.remove(first_d)
    require(len(result) == len(set(result)) == 3, ("suffix multiplicity changed", ds))
    return tuple(sorted(result))


def low_tables(stream, needed):
    rows = {d: [] for d in needed}
    for label in range(stream.first + 1, stream.high_floor):
        d = stream.L // gcd(stream.L, label)
        if d in rows:
            rows[d].append((delta(stream, label), label))
    for d in rows:
        rows[d].sort(key=lambda item: (item[0], -item[1]), reverse=True)
        rows[d] = tuple(rows[d])
    return rows


def high_uppers(stream, needed):
    """Direct unit-ray construction, independent of the frontier top tables."""
    rows = {d: (F(0), None, None) for d in needed}
    for residue in range(1, stream.L):
        d = stream.L // gcd(stream.L, residue)
        if d not in rows:
            continue
        amplitude = residue * delta(stream, residue)
        if amplitude <= 0:
            continue
        label = first_on_ray(residue, stream.L, stream.high_floor)
        value = amplitude / label
        candidate = (value, label, residue)
        incumbent = rows[d]
        if candidate[0] > incumbent[0] or (
            candidate[0] == incumbent[0] and candidate[1] < incumbent[1]
        ):
            rows[d] = candidate
    require(all(value >= 0 for value, _label, _residue in rows.values()), "negative upper")
    return rows


def two_high_gap(stream, residuals, low, high):
    required = stream.lower - stream.first_delta
    best = None
    witness = None
    for ds in residuals:
        suffix = suffix_denominators(ds, stream.first_d)
        for mask in range(1 << 3):
            if mask.bit_count() < 2:
                continue
            value = F(0)
            packet = []
            possible = True
            for index, d in enumerate(suffix):
                if (mask >> index) & 1:
                    item = high[d]
                elif low[d]:
                    item = low[d][0]
                else:
                    possible = False
                    break
                value += item[0]
                packet.append((d, item[1], bool((mask >> index) & 1)))
            if not possible:
                continue
            if best is None or value > best:
                best = value
                witness = (ds, tuple(packet))
    require(best is not None, "missing two-high upper")
    return required - best, witness


def one_high_cases(stream, residuals, low, high):
    required = stream.lower - stream.first_delta
    cases = []
    for ds in residuals:
        suffix = suffix_denominators(ds, stream.first_d)
        for high_d in suffix:
            low_ds = tuple(d for d in suffix if d != high_d)
            if not low[low_ds[0]] or not low[low_ds[1]]:
                continue
            threshold = required - high[high_d][0]
            first_rows = low[low_ds[0]]
            second_rows = low[low_ds[1]]
            for first_value, first_label in first_rows:
                if first_value + second_rows[0][0] < threshold:
                    break
                for second_value, second_label in second_rows:
                    total = first_value + second_value
                    if total < threshold:
                        break
                    require(first_label != second_label, "low labels collided")
                    cases.append(
                        (
                            ds,
                            high_d,
                            (low_ds[0], first_label),
                            (low_ds[1], second_label),
                            total + high[high_d][0] - required,
                        )
                    )
    require(len(cases) == len(set(cases)), "duplicate literal case")
    return tuple(cases)


def zero_high_scalar_passes(stream, residuals, low):
    """Hostile control: zero-high scalar maxima, excluded only by high wall."""
    required = stream.lower - stream.first_delta
    rows = []
    for ds in residuals:
        suffix = suffix_denominators(ds, stream.first_d)
        if not all(low[d] for d in suffix):
            continue
        labels = tuple(low[d][0][1] for d in suffix)
        value = sum((low[d][0][0] for d in suffix), F(0))
        if value >= required:
            rows.append((ds, labels, value - required))
    return tuple(rows)


def cell_is_clean(cell, label, modulus):
    residue = (label * cell) % modulus
    return 14 * residue >= modulus and 14 * (residue + label) <= 13 * modulus


def torsion_pair(stream, labels, high_d):
    torsion = 2 if high_d % 2 == 0 else 3
    require(high_d % torsion == 0 and torsion in (2, 3), ("bad torsion", high_d))
    shift = high_d // torsion
    clean = []
    for left, right in stream.ranges:
        for cell in range(left, right):
            if all(cell_is_clean(cell, label, stream.L) for label in labels):
                clean.append(cell)
    by_residue = {}
    for cell in clean:
        residue = cell % high_d
        by_residue.setdefault(residue, cell)
    for residue, first in by_residue.items():
        target = (residue + shift) % high_d
        if target in by_residue:
            second = by_residue[target]
            require((second - first) % high_d == shift, "torsion shift changed")
            phase_gaps = {
                (unit * shift) % high_d
                for unit in range(1, high_d)
                if gcd(unit, high_d) == 1
            }
            require(
                min(min(gap, high_d - gap) for gap in phase_gaps) * 7 > high_d,
                ("danger arcs may meet", high_d, phase_gaps),
            )
            require(
                all(cell_is_clean(c, label, stream.L) for c in (first, second) for label in labels),
                "fixed comb entered witness cell",
            )
            return torsion, first, second, len(clean), tuple(sorted(phase_gaps))
    return torsion, None, None, len(clean), ()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    semantic = hashlib.sha256()
    total_cases = 0
    total_failures = 0
    total_zero_high = 0
    torsion_histogram = Counter()
    lines = [
        "LRC14 projected k=3 z1=312 finite-low-pair torsion independent referee",
        f"frontier_source_sha256={file_sha256(FRONTIER)}",
        f"frontier_output_sha256={file_sha256(FRONTIER_OUTPUT)}",
        "inherited_gate=first label 312 lies below every body high floor;projected high-wall theorem requires at least one later high label",
        "hostile_control=zero-high scalar arithmetic is tested separately and is not used as a proof step",
    ]
    for body, residuals in frontier.EXPECTED_RESIDUALS.items():
        stream = ray.Stream(body)
        require(
            stream.first == frontier.FIRST and stream.first < stream.high_floor,
            ("inherited high-wall gate inactive", body, stream.first, stream.high_floor),
        )
        needed = {
            d
            for ds in residuals
            for d in suffix_denominators(ds, stream.first_d)
        }
        low = low_tables(stream, needed)
        high = high_uppers(stream, needed)
        zero_high = zero_high_scalar_passes(stream, residuals, low)
        require(
            len(zero_high) == EXPECTED_ZERO_HIGH_SCALAR_PASSES[body],
            ("zero-high hostile ledger changed", body, len(zero_high)),
        )
        gap, gap_witness = two_high_gap(stream, residuals, low, high)
        require(gap == EXPECTED_GAPS[body] > 0, ("two-high gap changed", body, gap))
        cases = one_high_cases(stream, residuals, low, high)
        require(len(cases) == EXPECTED_CASE_COUNTS[body], ("case count changed", body, len(cases)))
        low_pairs = {
            tuple(sorted((case[2][1], case[3][1])))
            for case in cases
        }
        require(
            len(low_pairs) == EXPECTED_LOW_PAIR_COUNTS[body],
            ("low-pair count changed", body, len(low_pairs)),
        )
        failures = 0
        witness_cache = {}
        for case in cases:
            ds, high_d, first_low, second_low, excess = case
            labels = (stream.first, first_low[1], second_low[1])
            cache_key = (labels, high_d)
            if cache_key not in witness_cache:
                witness_cache[cache_key] = torsion_pair(stream, labels, high_d)
            witness = witness_cache[cache_key]
            torsion, first_cell, second_cell, clean_count, phase_gaps = witness
            torsion_histogram[torsion] += 1
            if first_cell is None:
                failures += 1
            row = (body, ds, high_d, first_low, second_low, excess, witness)
            semantic.update(f"{row}\n".encode())
            lines.append(
                f"case=E:{body};ds:{ds};high_d:{high_d};lows:{first_low},{second_low};"
                f"scalar_excess:{ftext(excess)};torsion:{torsion};cells:{first_cell},{second_cell};"
                f"clean_cells:{clean_count};phase_gaps:{phase_gaps}"
            )
        require(failures == 0, ("torsion pair failure", body, failures))
        total_cases += len(cases)
        total_failures += failures
        total_zero_high += len(zero_high)
        lines.append(
            f"body={body};two_high_gap={ftext(gap)};gap_witness={gap_witness};"
            f"zero_high_scalar_passes={len(zero_high)}/{len(residuals)};"
            f"cases={len(cases)};distinct_low_pairs={len(low_pairs)};pair_failures={failures}"
        )
    require(total_cases == 80 and total_failures == 0, "global case ledger changed")
    require(total_zero_high == 69, ("global zero-high hostile ledger changed", total_zero_high))
    require(torsion_histogram == Counter({2: 78, 3: 2}), ("torsion histogram changed", torsion_histogram))
    semantic_hash = semantic.hexdigest()
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines.extend(
        (
            "logical_split=zero high excluded by inherited projected high-wall gate;two-or-more high excluded by positive exact gaps;therefore exactly one high",
            "zero_high_hostile_total=69/70;body_split=24/41/2/2;this demonstrates the inherited gate is nonredundant",
            f"totals=cases:{total_cases};pair_failures:{total_failures};torsion:{dict(sorted(torsion_histogram.items()))}",
            f"semantic_sha256={semantic_hash}",
            "all_independent_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
