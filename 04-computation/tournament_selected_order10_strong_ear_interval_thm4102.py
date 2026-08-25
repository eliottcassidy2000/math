#!/usr/bin/env python3
"""Selected order-ten strong-ear interval compiler for THM-4102.

Start with the deterministic first labelled strong order-nine witness retained
for each of the 1,482 values in THM-4097, and evaluate every nonconstant ear
over those parents.  This is a proved construction subset, not a complete
order-ten census.  Every executable gate uses ``require`` so ``python -O``
removes no verification.
"""

from __future__ import annotations

import importlib.util
import json
from hashlib import sha256
from math import isqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PARENT_COMPILER = HERE / "tournament_order9_strong_ear_spectrum_thm4097.py"
PARENT_COMPILER_SHA256 = (
    "610ca5850b272e0e75c574f2c1a710a0b96c75cc7191b1e1f1a03dfbdd1378d6"
)
KEY_VALUES = (125, 249, 2887, 2933, 14649, 14653, 14655, 15055, 15551, 15621)


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_parent_compiler():
    require(file_sha256(PARENT_COMPILER) == PARENT_COMPILER_SHA256,
            "THM-4097 compiler hash")
    spec = importlib.util.spec_from_file_location("thm4097_for_thm4102", PARENT_COMPILER)
    require(spec is not None and spec.loader is not None, "parent compiler import")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def first_missing_prime(values: set[int], start: int) -> int:
    candidate = start if start % 2 else start + 1
    while True:
        if is_prime(candidate) and candidate not in values:
            return candidate
        candidate += 2


def first_missing_seven_prime(values: set[int], start_prime: int) -> int:
    candidate = start_prime if start_prime % 2 else start_prime + 1
    while True:
        if is_prime(candidate) and 7 * candidate not in values:
            return candidate
        candidate += 2


def main() -> None:
    parent = load_parent_compiler()
    engine = parent.load_engine()
    representatives, counts = engine.generate(8)
    require(counts == parent.A000568, "A000568 class counts through order eight")

    # Reproduce THM-4097's deterministic first witness for every strong value.
    selected: dict[int, list[int]] = {}
    selected_codes: dict[int, int] = {}
    for old in representatives[8]:
        if not engine.is_strong(8, old):
            continue
        state = parent.boundary_state(old)
        for signature in range(1, (1 << 8) - 1):
            value = parent.insertion_h(state, signature)
            if value not in selected:
                child = engine.extend(old, 8, signature)
                selected[value] = child
                selected_codes[value] = parent.raw_code(child)

    historical, all_classes, strong_classes = parent.read_historical_strong_values()
    require(set(selected) == historical, "selected parents cover exact THM-4097 spectrum")
    require(len(selected) == 1_482, "one selected parent per order-nine value")
    require((all_classes, strong_classes) == (191_536, 178_133),
            "historical order-nine controls")

    values: set[int] = set()
    witnesses: dict[int, tuple[int, int, int]] = {}
    parent_failures: list[tuple[int, int, bool]] = []
    ears_checked = 0
    for parent_h, old in sorted(selected.items()):
        direct_h = engine.Hcount(9, old)
        strong = engine.is_strong(9, old)
        if direct_h != parent_h or not strong:
            parent_failures.append((parent_h, direct_h, strong))
        state = parent.boundary_state(old)
        for signature in range(1, (1 << 9) - 1):
            ears_checked += 1
            value = parent.insertion_h(state, signature)
            values.add(value)
            if value not in witnesses:
                child = engine.extend(old, 9, signature)
                witnesses[value] = (parent_h, signature, parent.raw_code(child))

    require(not parent_failures, "all selected order-nine parents direct/strong")
    require(ears_checked == 1_482 * 510 == 755_820, "complete selected-ear universe")
    require(len(values) == 7_566, "selected order-ten value count")
    require((min(values), max(values)) == (125, 15_621), "selected image extrema")

    intervals = parent.odd_intervals(values)
    require((249, 14_649, 7_201) in intervals, "primary solid interval")
    require((15_055, 15_551, 249) in intervals, "secondary solid interval")
    require(14_651 not in values and 14_653 in values and 14_655 in values,
            "primary right boundary and two sporadic successors")

    # Directly re-evaluate one retained labelled witness for every value.
    witness_failures: list[tuple[int, int, bool]] = []
    key_rows: dict[str, dict[str, int]] = {}
    for value, (parent_h, signature, code) in sorted(witnesses.items()):
        child = parent.decode_raw(code, 10)
        direct_h = engine.Hcount(10, child)
        strong = engine.is_strong(10, child)
        if direct_h != value or not strong:
            witness_failures.append((value, direct_h, strong))
        if value in KEY_VALUES:
            key_rows[str(value)] = {
                "parent_h": parent_h,
                "signature": signature,
                "cut_weight": signature.bit_count(),
                "code": code,
            }
    require(not witness_failures, "direct DP/strongness check for retained witnesses")
    require(set(key_rows) == {str(value) for value in KEY_VALUES}, "key boundary witnesses")
    require(key_rows["2887"]["code"] == 25_133_469_073_343, "stable H=2887 code")
    require(key_rows["2933"]["code"] == 34_960_494_755_519, "stable H=2933 code")

    ordinary_missing = first_missing_prime(values, 2_887)
    seven_missing_prime = first_missing_seven_prime(values, 419)
    require(ordinary_missing == 14_657, "next ordinary-prime target")
    require(seven_missing_prime == 2_111 and 7 * seven_missing_prime == 14_777,
            "next seven-prime target")
    require(14_651 == 49 * 13 * 23, "multiplicative bridge over selected hole")

    longest = sorted(intervals, key=lambda row: row[2], reverse=True)[:12]
    ledger = {
        "parent_compiler_sha256": PARENT_COMPILER_SHA256,
        "selected_parent_codes": sorted(selected_codes.items()),
        "ears_checked": ears_checked,
        "values": sorted(values),
        "intervals": intervals,
        "key_rows": key_rows,
        "lane_next": [ordinary_missing, seven_missing_prime],
        "direct_checks": [len(selected), len(witnesses)],
    }
    semantic = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("SELECTED ORDER-TEN STRONG-EAR INTERVAL")
    print("selected_order9_parents=", len(selected))
    print("nonconstant_ears_checked=", ears_checked)
    print("selected_value_count=", len(values), "min=", min(values), "max=", max(values))
    print("primary_solid_interval= [249,14649] count=7201")
    print("secondary_solid_interval= [15055,15551] count=249")
    print("first_intervals=", intervals[:12])
    print("longest_intervals=", longest)
    print("direct_parent_checks=", len(selected), "failures=", len(parent_failures))
    print("direct_witness_checks=", len(witnesses), "failures=", len(witness_failures))
    print("key_witnesses=", key_rows)
    print("global_allowed_prefix_through=14655 multiplicative_bridge_14651=49*13*23")
    print("ordinary_prime_lane_through=14653 next_unforced=", ordinary_missing)
    print("seven_prime_lane_p_through=2099 next_unforced_value=", 7 * seven_missing_prime)
    print("semantic_sha256=", semantic)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
