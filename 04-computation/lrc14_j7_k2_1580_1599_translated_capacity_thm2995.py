#!/usr/bin/env python3
"""Exact translated-capacity closure for the two non-cheap 1580..1599 rows."""

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
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_1580_1599_translated_capacity_thm2995.out"
EXPECTED_ENGINE_SHA256 = "7cc93dfe2e9365aa1e313d97ab19e76bce0943fd9f6a243e8c1b33c80d439ef4"
EXPECTED_SCALAR_SHA256 = "01577b566a51d37c96c2fed783f133e27403ad3566c3df002aa5f129883f2ba4"
LARGE_KEY = (1586, (1, 8, 10, 12, 13, 14))
EXCEPTIONAL_KEY = (1586, (1, 10, 11, 12, 13, 14))


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def load_engine():
    require(lf_sha256(ENGINE) == EXPECTED_ENGINE_SHA256, ("engine changed", lf_sha256(ENGINE)))
    spec = importlib.util.spec_from_file_location("thm2995_translated_capacity_engine", ENGINE)
    require(spec is not None and spec.loader is not None, ENGINE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_rows():
    require(lf_sha256(SCALAR) == EXPECTED_SCALAR_SHA256, ("scalar output changed", lf_sha256(SCALAR)))
    rows = []
    for line in SCALAR.read_text(encoding="utf-8").splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
        rows.append(
            {
                "body": tuple(map(int, fields["E"].split(","))),
                "L": int(fields["L"]),
                "floor": int(fields["largest_floor"]),
                "first": int(fields["z1"]),
                "first_delta": Fraction(fields["delta1"]),
                "lower": Fraction(fields["lower"]),
                "high_tail": "HIGH-TAIL" in fields["suffix"],
            }
        )
    require(len(rows) == 26, len(rows))
    return tuple(sorted(rows, key=lambda row: (row["first"], row["body"])))


def kappa(denominator: int) -> int:
    return (denominator + 6) // 7


def semantic_digest(lines) -> str:
    return hashlib.sha256(("\n".join(lines) + "\n").encode("ascii")).hexdigest()


def large_row_audit(m, row):
    required, cutoff, candidates, _heads = m.finite_candidates(row)
    require(cutoff > 0 and len(candidates) == 12_588, (cutoff, len(candidates)))
    values = tuple(item[0] for item in candidates)
    labels = tuple(item[1] for item in candidates)
    modulus = row["L"]
    prefix_denominators = defaultdict(set)
    prefix_tail_counts = Counter()
    total_packets = 0
    n = len(candidates)
    for i in range(n - 3):
        if sum(values[i : i + 4], m.F()) < required:
            break
        for j in range(i + 1, n - 2):
            if values[i] + sum(values[j : j + 3], m.F()) < required:
                break
            for k in range(j + 1, n - 1):
                partial = values[i] + values[j] + values[k]
                if partial + values[k + 1] < required:
                    break
                prefix = tuple(sorted((labels[i], labels[j], labels[k])))
                for ell in range(k + 1, n):
                    if partial + values[ell] < required:
                        break
                    denominator = modulus // gcd(modulus, labels[ell])
                    require(denominator > 1, (labels[ell], modulus))
                    prefix_denominators[prefix].add(denominator)
                    prefix_tail_counts[prefix] += 1
                    total_packets += 1

    cells = m.base.complete_body_cells(m.base.carrier_for(row["body"]), modulus)
    cells_array = np.asarray(cells, dtype=np.int64)
    label_cache = {row["first"]: m.vector_safe_mask(cells_array, modulus, row["first"])}
    failures = []
    records = []
    tests = 0
    min_slack = None
    min_witness = None
    fixed_cell_min = None
    for prefix in sorted(prefix_denominators):
        mask = label_cache[row["first"]]
        for label in prefix:
            if label not in label_cache:
                label_cache[label] = m.vector_safe_mask(cells_array, modulus, label)
            mask &= label_cache[label]
        fixed_cells = mask.bit_count()
        fixed_cell_min = fixed_cells if fixed_cell_min is None else min(fixed_cell_min, fixed_cells)
        for denominator in sorted(prefix_denominators[prefix]):
            presence = m.vector_residue_presence(mask, modulus, denominator)
            cardinality = presence.bit_count()
            slack = cardinality - kappa(denominator)
            tests += 1
            witness = (prefix, denominator, cardinality, kappa(denominator), fixed_cells)
            records.append("|".join((",".join(map(str, prefix)), str(denominator), str(cardinality), str(fixed_cells))))
            if min_slack is None or slack < min_slack:
                min_slack = slack
                min_witness = witness
            if slack <= 0:
                failures.append(witness)
    require(
        total_packets == 29_372
        and len(prefix_denominators) == 249
        and tests == 5_507
        and not failures
        and fixed_cell_min == 23_688
        and min_slack == 6,
        (total_packets, len(prefix_denominators), tests, len(failures), fixed_cell_min, min_slack),
    )
    return {
        "cutoff": cutoff,
        "packets": total_packets,
        "prefixes": len(prefix_denominators),
        "tests": tests,
        "tail_min": min(prefix_tail_counts.values()),
        "tail_max": max(prefix_tail_counts.values()),
        "fixed_min": fixed_cell_min,
        "min_slack": min_slack,
        "min_witness": min_witness,
        "digest": semantic_digest(records),
    }


def exceptional_row_audit(m, row):
    modulus, first, body = row["L"], row["first"], row["body"]
    required = row["lower"] - row["first_delta"]
    values = m.amplitude(body, modulus)
    lows = sorted(
        (
            (m.ray_value(int(values[label - 1]), label), label)
            for label in range(first + 1, row["floor"])
            if int(values[label - 1]) > 0
        ),
        key=lambda item: (-item[0], item[1]),
    )
    triples = tuple(m.finite_triples(tuple(lows), required))
    denominators = sorted(
        {
            modulus // gcd(modulus, residue)
            for residue, raw_value in enumerate(values, 1)
            if int(raw_value) <= 0
        }
    )
    cells = m.base.complete_body_cells(m.base.carrier_for(body), modulus)
    cells_array = np.asarray(cells, dtype=np.int64)
    safe_cache = {first: m.vector_safe_mask(cells_array, modulus, first)}
    failures = []
    records = []
    min_slack = None
    min_witness = None
    tests = cell_gate_passes = 0
    for labels, total in triples:
        require(total >= required, (labels, total, required))
        mask = safe_cache[first]
        for label in labels:
            if label not in safe_cache:
                safe_cache[label] = m.vector_safe_mask(cells_array, modulus, label)
            mask &= safe_cache[label]
        fixed_cells = mask.bit_count()
        for denominator in denominators:
            presence = m.vector_residue_presence(mask, modulus, denominator)
            cardinality = presence.bit_count()
            slack = cardinality - kappa(denominator)
            tests += 1
            if fixed_cells * denominator > kappa(denominator) * modulus:
                cell_gate_passes += 1
            witness = (labels, denominator, cardinality, kappa(denominator), fixed_cells)
            records.append("|".join((",".join(map(str, labels)), str(denominator), str(cardinality), str(fixed_cells))))
            if min_slack is None or slack < min_slack:
                min_slack = slack
                min_witness = witness
            if slack <= 0:
                failures.append(witness)
    require(
        len(lows) == 35_697
        and len(triples) == 5
        and len(denominators) == 191
        and tests == 955
        and not failures
        and cell_gate_passes == 889
        and min_slack == 1,
        (len(lows), len(triples), len(denominators), tests, len(failures), cell_gate_passes, min_slack),
    )
    return {
        "required": required,
        "lows": len(lows),
        "triples": len(triples),
        "denominators": len(denominators),
        "tests": tests,
        "cell_gate_passes": cell_gate_passes,
        "min_slack": min_slack,
        "min_witness": min_witness,
        "digest": semantic_digest(records),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    m = load_engine()
    rows = parse_rows()
    large = next(row for row in rows if (row["first"], row["body"]) == LARGE_KEY)
    exceptional = next(row for row in rows if (row["first"], row["body"]) == EXCEPTIONAL_KEY)
    large_result = large_row_audit(m, large)
    exceptional_result = exceptional_row_audit(m, exceptional)
    lines = [
        "LRC14 projected k2 1580..1599 translated-capacity closure",
        f"dependency={ENGINE.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_ENGINE_SHA256}",
        f"dependency={SCALAR.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_SCALAR_SHA256}",
        "capacity=kappa(d)=ceil(d/7);strict_open=1;complete_cell_projection=THM-2984",
        "LARGE;"
        f"z=1586;E=1,8,10,12,13,14;eta={large_result['cutoff']};"
        f"packets={large_result['packets']};prefixes={large_result['prefixes']};tests={large_result['tests']};"
        f"tail_min={large_result['tail_min']};tail_max={large_result['tail_max']};"
        f"fixed_cell_min={large_result['fixed_min']};min_slack={large_result['min_slack']};"
        f"min_witness={large_result['min_witness']};record_sha256={large_result['digest']}",
        "EXCEPTIONAL;"
        f"z=1586;E=1,10,11,12,13,14;required={exceptional_result['required']};"
        f"positive_lows={exceptional_result['lows']};triples={exceptional_result['triples']};"
        f"denominators={exceptional_result['denominators']};tests={exceptional_result['tests']};"
        f"cell_gate_passes={exceptional_result['cell_gate_passes']};"
        f"min_slack={exceptional_result['min_slack']};min_witness={exceptional_result['min_witness']};"
        f"record_sha256={exceptional_result['digest']}",
        "all_exact_checks=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
