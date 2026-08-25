#!/usr/bin/env python3
"""Reclassify THM-4013 fibres by center-prediction error mass.

Scratch only.  This imports the canonical producer but changes no canonical
artifact.  Its universe and feature caps are exactly those of THM-4013.
"""

from __future__ import annotations

from collections import defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
PRODUCER = (
    ROOT
    / "04-computation"
    / "rule30_depth15_history_adaptive_route_thm4013.py"
)
EXPECTED = {
    "none": (
        (5637, 4435, 26149, 13237, 43),
        (15235, 9050, 20917, 23701, 15),
        (45916, 8105, 8661, 48213, 4),
        (61838, 1750, 1760, 62015, 2),
        (64912, 302, 307, 64921, 2),
        (65371, 71, 74, 65387, 2),
        (65450, 41, 43, 65449, 2),
    ),
    "target": (
        (10074, 23, 23, 65489, 1),
        (24284, 18, 18, 65499, 1),
        (54020, 5, 5, 65525, 1),
        (63587, 2, 2, 65531, 1),
        (65214, 0, 0, 65535, 0),
        (65442, 0, 0, 65535, 0),
        (65491, 0, 0, 65535, 0),
    ),
    "carry": (
        (5664, 4433, 26136, 13263, 43),
        (15256, 9046, 20905, 23725, 15),
        (45923, 8102, 8658, 48219, 4),
        (61840, 1749, 1759, 62017, 2),
        (64912, 302, 307, 64921, 2),
        (65371, 71, 74, 65387, 2),
        (65450, 41, 43, 65449, 2),
    ),
    "full": (
        (10097, 0, 0, 65535, 0),
        (24302, 0, 0, 65535, 0),
        (54025, 0, 0, 65535, 0),
        (63589, 0, 0, 65535, 0),
        (65214, 0, 0, 65535, 0),
        (65442, 0, 0, 65535, 0),
        (65491, 0, 0, 65535, 0),
    ),
}


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def load_producer():
    spec = importlib.util.spec_from_file_location("rule30_thm4013_mass_source", PRODUCER)
    require(spec is not None and spec.loader is not None, "producer loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def carry_shadow(record, mode: str):
    state = record.carry_state
    if mode == "none":
        return ()
    if state[0] == "pre":
        return ("pre",)
    if mode == "target":
        return ("visible", state[1])
    if mode == "carry":
        return ("visible", state[2])
    require(mode == "full", ("carry mode", mode))
    return state


def key(record, mode: str, history: int):
    return (
        record.gap,
        record.owner_bank,
        record.owner_portrait,
        record.odd_block,
        carry_shadow(record, mode),
        record.projective[-history:] if history else (),
    )


def census(records, mode: str, history: int):
    fibres: dict[object, list[int]] = defaultdict(lambda: [0, 0])
    for record in records:
        fibres[key(record, mode, history)][record.center] += 1
    error = sum(min(pair) for pair in fibres.values())
    mismatches = sum(pair[0] > 0 and pair[1] > 0 for pair in fibres.values())
    variation = sum(abs(pair[1] - pair[0]) for pair in fibres.values())
    maximum_minority = max((min(pair) for pair in fibres.values()), default=0)
    require(variation == len(records) - 2 * error, ("mass identity", mode, history))
    return len(fibres), mismatches, error, variation, maximum_minority


def main() -> None:
    producer = load_producer()
    dependency = producer.load_dependency()
    records, _rows, _low, _v, _gaps, controls = producer.build_records(dependency)
    require(len(records) == 65535, "canonical record universe")
    for record in records:
        state = record.carry_state
        predicted = 0 if state[0] == "pre" else state[1] ^ state[2]
        require(predicted == record.center, ("full carry determines center", record.n))

    tables = {
        mode: tuple(census(records, mode, history) for history in range(7))
        for mode in ("none", "target", "carry", "full")
    }
    require(tables == EXPECTED, "frozen reclassification table")
    semantic_sha256 = hashlib.sha256(
        json.dumps((controls, tables), separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("RULE30_PRIZE23_OBSERVER_MASS_RECLASSIFICATION_SCRATCH")
    print("universe=THM4013_records_n=1..65535;projective_history=0..6")
    print(f"producer={PRODUCER.relative_to(ROOT).as_posix()}")
    print(f"dependency_sha256={producer.DEPENDENCY_SHA256}")
    print(
        "columns=(fibres,center_mismatch_fibres,optimal_center_errors,"
        "conditional_L1,max_fibre_minority)"
    )
    for mode in ("none", "target", "carry", "full"):
        print(f"carry_mode_{mode}={tables[mode]}")
    print(f"semantic_sha256={semantic_sha256}")
    print(
        "interpretation=THM4013 target mismatches are not center mismatches; full "
        "carry determines center; target mode deletes only the lower carry bit"
    )
    print("scope=FINITE-EXACT canonical capped observer family; no all-scale or prize claim")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
