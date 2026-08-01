#!/usr/bin/env python3
"""Exact finite-packet torsion closure for the 24 cheap 1580..1599 rows."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "04-computation/lrc14_j7_k2_1600_1679_zero_high_and_high_order_phase_closure_thm2980.py"
SCALAR = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1580_1599_thm2995.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_1580_1599_cheap_finite_torsion_thm2995.out"
EXPECTED_ENGINE_SHA256 = "7cc93dfe2e9365aa1e313d97ab19e76bce0943fd9f6a243e8c1b33c80d439ef4"
EXPECTED_SCALAR_SHA256 = "01577b566a51d37c96c2fed783f133e27403ad3566c3df002aa5f129883f2ba4"
SKIPPED_KEYS = {
    (1586, (1, 8, 10, 12, 13, 14)),
    (1586, (1, 10, 11, 12, 13, 14)),
}
EXPECTED_HEIGHT_PACKETS = {1581: 498, 1586: 6606, 1588: 1, 1590: 95, 1594: 43, 1595: 2181, 1599: 2554}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def load_engine():
    require(lf_sha256(ENGINE) == EXPECTED_ENGINE_SHA256, ("engine changed", lf_sha256(ENGINE)))
    spec = importlib.util.spec_from_file_location("thm2995_cheap_torsion_engine", ENGINE)
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


def semantic_digest(lines) -> str:
    return hashlib.sha256(("\n".join(lines) + "\n").encode("ascii")).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    m = load_engine()
    records = []
    witness_records = []
    packet_count = failure_count = q2_count = 0
    witness_hist = Counter()
    height_packets = Counter()
    skipped = set()
    for row in parse_rows():
        required, cutoff, candidates, _heads = m.finite_candidates(row)
        key = (row["first"], row["body"])
        if cutoff <= 0 or len(candidates) > 1000:
            skipped.add(key)
            continue
        require(cutoff > 0, (key, cutoff))
        modulus = row["L"]
        cells = m.base.complete_body_cells(m.base.carrier_for(row["body"]), modulus)
        labels = {row["first"]}
        labels.update(item[1] for item in candidates)
        safe_masks = {label: m.base.safe_mask(cells, modulus, label) for label in sorted(labels)}
        local_packets = local_failures = local_q2 = 0
        for suffix, total in m.admissible_suffixes(candidates, required):
            packet = (row["first"], *suffix)
            require(len(set(packet)) == 5 and total >= required, (key, packet, total, required))
            local_packets += 1
            packet_count += 1
            height_packets[row["first"]] += 1
            witness = None
            for index, label in enumerate(packet):
                fixed = packet[:index] + packet[index + 1 :]
                mask = safe_masks[fixed[0]]
                for other in fixed[1:]:
                    mask &= safe_masks[other]
                denominator, unit = m.denominator_and_unit(label, modulus)
                orders = m.base.fast_torsion_orders(mask, modulus, denominator)
                if orders:
                    witness = (index, label, denominator, unit, mask.bit_count(), orders[0])
                    break
            if witness is None:
                local_failures += 1
                failure_count += 1
                continue
            require(witness[-1][0] == 2, (key, packet, witness))
            local_q2 += 1
            q2_count += 1
            witness_hist[(witness[0], witness[2])] += 1
            witness_records.append(
                f"{row['first']}|{','.join(map(str,row['body']))}|{','.join(map(str,packet))}|{repr(witness)}"
            )
        require(local_packets > 0 and local_failures == 0 and local_q2 == local_packets, (key, local_packets, local_failures, local_q2))
        records.append(
            "ROW;"
            f"z={row['first']};E={','.join(map(str,row['body']))};"
            f"candidates={len(candidates)};packets={local_packets};failures={local_failures};q2={local_q2}"
        )

    require(skipped == SKIPPED_KEYS, skipped)
    require(
        len(records) == 24
        and packet_count == 11_978
        and failure_count == 0
        and q2_count == packet_count
        and dict(height_packets) == EXPECTED_HEIGHT_PACKETS,
        (len(records), packet_count, failure_count, q2_count, dict(height_packets)),
    )
    lines = [
        "LRC14 projected k2 1580..1599 cheap finite-packet torsion closure",
        f"dependency={ENGINE.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_ENGINE_SHA256}",
        f"dependency={SCALAR.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_SCALAR_SHA256}",
        f"rows={len(records)};packets={packet_count};failures={failure_count};q2={q2_count};skipped=2",
        f"height_packets={dict(sorted(height_packets.items()))}",
        f"q2_witness_index_denominator={dict(sorted(witness_hist.items()))}",
        f"row_record_sha256={semantic_digest(records)}",
        f"packet_witness_sha256={semantic_digest(witness_records)}",
        *records,
        "all_exact_checks=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
