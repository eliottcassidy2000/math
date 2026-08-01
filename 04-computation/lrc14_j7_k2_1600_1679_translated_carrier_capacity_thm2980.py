#!/usr/bin/env python3
"""Exact translated-capacity audit for the 27 THM-2980 repair carriers."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from hashlib import sha256
from math import gcd
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc14_j7_k2_1600_1679_zero_high_and_high_order_phase_closure_thm2980.py"
EXPECTED_SOURCE_SHA256 = "7cc93dfe2e9365aa1e313d97ab19e76bce0943fd9f6a243e8c1b33c80d439ef4"
TARGET_BODIES = {
    (1, 8, 10, 12, 13, 14): 18,
    (1, 10, 11, 12, 13, 14): 9,
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


def load_source():
    spec = importlib.util.spec_from_file_location("thm2980_cardinality_shadow", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(SOURCE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def kappa(denominator: int) -> int:
    """Sharp capacity of an arbitrarily translated open 1/7-band."""
    return (denominator + 6) // 7


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    source_hash = lf_sha256(SOURCE)
    require(source_hash == EXPECTED_SOURCE_SHA256, "THM-2980 closure source hash changed")

    # MISTAKE-334 hostile: the centered beta(28)=3 bound is false after an
    # arbitrary translation.  At theta=-3/56 the four residues 0,1,2,3 all
    # lie strictly inside the translated danger band.  The correct capacity
    # is kappa(28)=4.
    translated_hostile = tuple(range(4))
    require(kappa(28) == 4, "translated capacity changed at d=28")
    require(
        all(abs(2 * residue - 3) < 4 for residue in translated_hostile),
        "translated d=28 hostile changed",
    )

    m = load_source()
    lines = [
        "LRC14 PROJECTED K2 NONPOSITIVE TRANSLATED-CARRIER CAPACITY",
        f"dependency={SOURCE.relative_to(ROOT).as_posix()};lf_sha256={source_hash}",
        "capacity=kappa(d)=ceil(d/7);phase_scope=arbitrary_translation;strict_open=1",
        "hostile=d28,u1,theta=-3/56,S=0,1,2,3;centered_beta=3;translated_kappa=4;all_danger=1",
    ]
    total_tests = total_failures = 0
    global_min_slack = None
    global_witness = None
    for row in m.parse_rows():
        body = row["body"]
        if row["first"] != 1612 or body not in TARGET_BODIES:
            continue
        modulus, first = row["L"], row["first"]
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
        require(len(triples) == TARGET_BODIES[body], (body, len(triples)))
        denominators = set()
        for residue, raw_value in enumerate(values, 1):
            if int(raw_value) <= 0:
                denominator = modulus // gcd(modulus, residue)
                require(denominator > 1, (body, residue, denominator))
                denominators.add(denominator)

        cells = m.base.complete_body_cells(m.base.carrier_for(body), modulus)
        cells_array = np.asarray(cells, dtype=np.int64)
        safe_cache = {first: m.vector_safe_mask(cells_array, modulus, first)}
        local_min_slack = None
        local_witness = None
        local_failures = 0
        local_cell_gate_passes = 0
        local_tests = 0
        for labels, low_total in triples:
            require(low_total >= required, (body, labels, low_total, required))
            mask = safe_cache[first]
            for label in labels:
                if label not in safe_cache:
                    safe_cache[label] = m.vector_safe_mask(cells_array, modulus, label)
                mask &= safe_cache[label]
            fixed_cells = mask.bit_count()
            for denominator in sorted(denominators):
                presence = m.vector_residue_presence(mask, modulus, denominator)
                cardinality = presence.bit_count()
                slack = cardinality - kappa(denominator)
                local_tests += 1
                if fixed_cells * denominator > kappa(denominator) * modulus:
                    local_cell_gate_passes += 1
                if local_min_slack is None or slack < local_min_slack:
                    local_min_slack = slack
                    local_witness = (labels, denominator, cardinality, kappa(denominator), fixed_cells)
                if global_min_slack is None or slack < global_min_slack:
                    global_min_slack = slack
                    global_witness = (
                        body,
                        labels,
                        denominator,
                        cardinality,
                        kappa(denominator),
                        fixed_cells,
                    )
                if slack <= 0:
                    local_failures += 1
        total_tests += local_tests
        total_failures += local_failures
        lines.append(
            "ROW;"
            f"E={','.join(map(str, body))};triples={len(triples)};denominators={len(denominators)};"
            f"carrier_denominator_tests={local_tests};failures={local_failures};"
            f"cell_count_gate_passes={local_cell_gate_passes};"
            f"min_slack={local_min_slack};"
            f"min_witness={local_witness}",
        )
    require(total_tests == 3861, total_tests)
    require(total_failures == 0, total_failures)
    require(global_min_slack == 1, global_min_slack)
    lines.append(
        f"TOTAL;carrier_denominator_tests={total_tests};failures={total_failures};"
        f"global_min_slack={global_min_slack};global_witness={global_witness}",
    )
    lines.extend(
        [
            "centered_beta_not_used=1;translated_capacity_is_load_bearing=1",
            "exhaustive_unit_replay_role=HOSTILE_REGRESSION_NOT_CONCEPTUAL_PROOF",
            "all_exact_checks=PASS",
        ]
    )
    transcript = "\n".join(lines) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
