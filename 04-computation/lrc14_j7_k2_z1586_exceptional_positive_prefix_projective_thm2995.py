#!/usr/bin/env python3
"""Compressed zero-high and exact one-high audit for the exceptional row."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "04-computation/lrc14_j7_k2_1600_1679_zero_high_and_high_order_phase_closure_thm2980.py"
SCALAR = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1580_1599_thm2995.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_z1586_exceptional_positive_prefix_projective_thm2995.out"
EXPECTED_ENGINE_SHA256 = "7cc93dfe2e9365aa1e313d97ab19e76bce0943fd9f6a243e8c1b33c80d439ef4"
EXPECTED_SCALAR_SHA256 = "01577b566a51d37c96c2fed783f133e27403ad3566c3df002aa5f129883f2ba4"
KEY = (1586, (1, 10, 11, 12, 13, 14))


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def load_engine():
    require(lf_sha256(ENGINE) == EXPECTED_ENGINE_SHA256, ("engine changed", lf_sha256(ENGINE)))
    spec = importlib.util.spec_from_file_location("thm2995_exceptional_prefix_projective_engine", ENGINE)
    require(spec is not None and spec.loader is not None, ENGINE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_row():
    require(lf_sha256(SCALAR) == EXPECTED_SCALAR_SHA256, ("scalar output changed", lf_sha256(SCALAR)))
    for line in SCALAR.read_text(encoding="utf-8").splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
        key = (int(fields["z1"]), tuple(map(int, fields["E"].split(","))))
        if key == KEY:
            return {
                "body": key[1],
                "L": int(fields["L"]),
                "floor": int(fields["largest_floor"]),
                "first": key[0],
                "first_delta": Fraction(fields["delta1"]),
                "lower": Fraction(fields["lower"]),
            }
    raise RuntimeError(KEY)


def kappa(denominator: int) -> int:
    return (denominator + 6) // 7


def semantic_digest(lines) -> str:
    return hashlib.sha256(("\n".join(lines) + "\n").encode("ascii")).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    m = load_engine()
    row = parse_row()
    body, modulus, first, floor = row["body"], row["L"], row["first"], row["floor"]
    required = row["lower"] - row["first_delta"]
    values = m.amplitude(body, modulus)
    lows = []
    for label in range(first + 1, floor):
        raw = int(values[label - 1])
        if raw > 0:
            lows.append((m.ray_value(raw, label), label, raw))
    lows.sort(key=lambda item: (-item[0], item[1]))

    maxima = {}
    rays = {}
    high_start = max(first + 1, floor)
    for residue, raw_value in enumerate(values, 1):
        raw = int(raw_value)
        if raw <= 0:
            continue
        denominator = modulus // gcd(modulus, residue)
        divisor = modulus // denominator
        unit = residue // divisor
        require(gcd(unit, denominator) == 1, (denominator, unit))
        label = m.first_at_or_after(residue, high_start, modulus)
        value = m.ray_value(raw, label)
        ray = (value, label, unit, residue, raw)
        rays.setdefault(denominator, []).append(ray)
        old = maxima.get(denominator)
        if old is None or value > old[0] or (value == old[0] and label < old[1]):
            maxima[denominator] = (value, label, residue, raw)

    best_d, best_high = max(maxima.items(), key=lambda item: (item[1][0], -item[0]))
    two_positive_gap = required - lows[0][0] - lows[1][0]
    three_gap = required - sum((item[0] for item in lows[:3]), m.F())
    mixed_gap = required - best_high[0] - lows[0][0] - lows[1][0]
    zero_gap = required - sum((item[0] for item in lows[:4]), m.F())
    two_gap = required - 2 * best_high[0] - lows[0][0] - lows[1][0]
    require(
        best_high[0] <= lows[1][0]
        and zero_gap < 0
        and two_positive_gap > 0
        and mixed_gap > 0
        and two_gap > 0,
        (best_high[0], lows[1][0], zero_gap, two_positive_gap, mixed_gap, two_gap),
    )

    # Exact zero-high quotient: order four lows, puncture the fourth, and
    # collapse its heights/units to the leading triple and denominator.
    low_candidates = tuple((value, label, label % modulus, raw) for value, label, raw in lows)
    low_values = tuple(item[0] for item in low_candidates)
    low_labels = tuple(item[1] for item in low_candidates)
    prefix_denominators = defaultdict(set)
    zero_packets = 0
    n = len(low_candidates)
    for i in range(n - 3):
        if sum(low_values[i : i + 4], m.F()) < required:
            break
        for j in range(i + 1, n - 2):
            if low_values[i] + sum(low_values[j : j + 3], m.F()) < required:
                break
            for k in range(j + 1, n - 1):
                partial = low_values[i] + low_values[j] + low_values[k]
                if partial + low_values[k + 1] < required:
                    break
                prefix = tuple(sorted((low_labels[i], low_labels[j], low_labels[k])))
                for ell in range(k + 1, n):
                    if partial + low_values[ell] < required:
                        break
                    denominator = modulus // gcd(modulus, low_labels[ell])
                    require(denominator > 1, (low_labels[ell], modulus))
                    prefix_denominators[prefix].add(denominator)
                    zero_packets += 1

    cells = m.base.complete_body_cells(m.base.carrier_for(body), modulus)
    cells_array = np.asarray(cells, dtype=np.int64)
    safe_cache = {first: m.vector_safe_mask(cells_array, modulus, first)}
    zero_tests = 0
    zero_coarse_passes = 0
    zero_exact_folds = 0
    zero_failures = []
    zero_records = []
    zero_min_slack = None
    zero_min_witness = None
    for prefix in sorted(prefix_denominators):
        mask = safe_cache[first]
        for label in prefix:
            if label not in safe_cache:
                safe_cache[label] = m.vector_safe_mask(cells_array, modulus, label)
            mask &= safe_cache[label]
        fixed_cells = mask.bit_count()
        for denominator in sorted(prefix_denominators[prefix]):
            zero_tests += 1
            if fixed_cells * denominator > kappa(denominator) * modulus:
                zero_coarse_passes += 1
                zero_records.append(
                    "|".join((",".join(map(str, prefix)), str(denominator), "coarse", str(fixed_cells)))
                )
                continue
            presence = m.vector_residue_presence(mask, modulus, denominator)
            cardinality = presence.bit_count()
            slack = cardinality - kappa(denominator)
            zero_exact_folds += 1
            witness = (prefix, denominator, cardinality, kappa(denominator), fixed_cells)
            zero_records.append("|".join((",".join(map(str, prefix)), str(denominator), str(cardinality), str(fixed_cells))))
            if zero_min_slack is None or slack < zero_min_slack:
                zero_min_slack = slack
                zero_min_witness = witness
            if slack <= 0:
                zero_failures.append(witness)

    # The exactly-one-high projective classes are still audited directly.
    classes = small = full = units = 0
    denominator_failures = Counter()
    small_records = []
    phase_records = []
    minimum_full_margin = None
    minimum_pair_margin = None
    for denominator, high in sorted(maxima.items()):
        threshold = required - high[0]
        for labels, low_total in m.finite_triples(lows, threshold):
            classes += 1
            mask = safe_cache[first]
            for label in labels:
                if label not in safe_cache:
                    safe_cache[label] = m.vector_safe_mask(cells_array, modulus, label)
                mask &= safe_cache[label]
            orders = m.base.fast_torsion_orders(mask, modulus, denominator)
            if orders:
                small += 1
                small_records.append(
                    f"d={denominator};lows={','.join(map(str,labels))};orders={','.join(str(item[0]) for item in orders)}"
                )
                continue
            full += 1
            denominator_failures[denominator] += 1
            residues_mask = m.vector_residue_presence(mask, modulus, denominator)
            require(residues_mask == (1 << denominator) - 1, (denominator, labels, mask.bit_count()))
            remaining = required - low_total
            relevant = tuple(ray for ray in rays[denominator] if ray[0] >= remaining)
            require(relevant, (denominator, labels, remaining))
            class_full_margin = None
            class_pair_margin = None
            for _value, _label, unit, _residue, _raw in relevant:
                units += 1
                span, largest_gap = m.cyclic_span(range(denominator), denominator, unit)
                full_margin = 7 * span - denominator
                require(largest_gap == 1 and full_margin == 6 * denominator - 7, (denominator, unit))
                t = (denominator + 6) // 7
                second = (pow(unit, -1, denominator) * t) % denominator
                pair_span, pair_gap = m.cyclic_span((0, second), denominator, unit)
                pair_margin = 7 * t - denominator
                require(
                    (unit * second) % denominator == t
                    and pair_span == min(t, denominator - t) == t
                    and pair_gap == denominator - t
                    and pair_margin >= 0,
                    (denominator, unit, second, t, pair_span, pair_gap),
                )
                class_full_margin = full_margin if class_full_margin is None else min(class_full_margin, full_margin)
                class_pair_margin = pair_margin if class_pair_margin is None else min(class_pair_margin, pair_margin)
                minimum_full_margin = full_margin if minimum_full_margin is None else min(minimum_full_margin, full_margin)
                minimum_pair_margin = pair_margin if minimum_pair_margin is None else min(minimum_pair_margin, pair_margin)
            phase_records.append(
                f"d={denominator};lows={','.join(map(str,labels))};safe={mask.bit_count()};"
                f"remaining={remaining};units={len(relevant)};full_margin={class_full_margin};pair_margin={class_pair_margin}"
            )

    require(zero_packets > 0 and zero_tests > 0 and not zero_failures, (zero_packets, zero_tests, zero_failures[:1]))
    require(
        classes == small + full and (minimum_pair_margin is None or minimum_pair_margin >= 0),
        (classes, small, full, minimum_pair_margin),
    )
    payload = "\n".join(
        (
            "LRC14 projected k2 z1586 exceptional positive prefix-projective audit",
            f"dependency={ENGINE.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_ENGINE_SHA256}",
            f"dependency={SCALAR.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_SCALAR_SHA256}",
            f"E={','.join(map(str,body))};L={modulus};floor={floor};lows={len(lows)};denominators={len(maxima)}",
            f"two_positive_gap={two_positive_gap};three_gap={three_gap};mixed_gap={mixed_gap};zero_gap={zero_gap};two_gap={two_gap}",
            f"best_high_d={best_d};best_high_label={best_high[1]};best_high_le_second_low=1",
            f"zero_packets={zero_packets};prefixes={len(prefix_denominators)};prefix_denominator_tests={zero_tests};"
            f"coarse_passes={zero_coarse_passes};exact_folds={zero_exact_folds};failures={len(zero_failures)};"
            f"exact_fold_min_slack={zero_min_slack};min_witness={zero_min_witness}",
            f"zero_record_sha256={semantic_digest(zero_records)}",
            f"one_high_classes={classes};small_torsion={small};full_projection={full};unit_instances={units}",
            f"full_projection_denominators={dict(sorted(denominator_failures.items()))};minimum_full_margin={minimum_full_margin};minimum_pair_margin={minimum_pair_margin}",
            f"small_class_sha256={semantic_digest(small_records)};phase_record_sha256={semantic_digest(phase_records)}",
            "all_exact_checks=PASS",
            "",
        )
    )
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
