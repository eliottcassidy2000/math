#!/usr/bin/env python3
"""Exact order-nine strong-ear H-spectrum compiler for THM-4097.

The universe is the complete A000568 tower through order eight from the frozen
S5 engine.  For every strong order-eight representative and every nonconstant
new-vertex cut, the script evaluates the exact start/end/exposed-slot boundary
polynomial.  Moon's nonseparating-vertex theorem makes this an exhaustive
order-nine strong-spectrum computation, without canonicalizing order nine.

Every executable gate uses ``require`` so ``python -O`` removes no check.
"""

from __future__ import annotations

import importlib.util
import json
from hashlib import sha256
from math import isqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ENGINE = HERE / "strong_H_spectrum_m8_isoclass_monad_s5.py"
HISTOGRAM = ROOT / "05-knowledge/results/h_spectrum_n9_histogram_monad_s6.tsv"
LEGACY_VALUES = ROOT / "05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out"

ENGINE_SHA256 = "6ab922de4a8b6f6c15ee0ca7e0b036c3821b3e800dbdf961de72194e73346419"
HISTOGRAM_SHA256 = "e7d5594879d4c3af739cb94ca8cfd944879c4d586747d993dd6687e60126552f"
LEGACY_VALUES_SHA256 = "27fbef5b06fcf0369eeb602e513c3802ea171492e1292a3f6afa3efeadef9f55"

A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
KEY_TARGETS = (75, 81, 85, 613, 623, 2881, 2885, 3357)


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_engine():
    require(file_sha256(ENGINE) == ENGINE_SHA256, "frozen order-eight engine hash")
    spec = importlib.util.spec_from_file_location("strong_s5_thm4097", ENGINE)
    require(spec is not None and spec.loader is not None, "engine import specification")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def all_end_counts(adj: list[int]) -> list[list[int]]:
    """dp[mask][v] = Hamiltonian paths of T[mask] ending at v."""
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for mask in range(1, 1 << n):
        for vertex, count in enumerate(dp[mask]):
            if count == 0:
                continue
            available = adj[vertex] & ~mask
            while available:
                bit = available & -available
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += count
                available -= bit
    return dp


def converse(adj: list[int]) -> list[int]:
    n = len(adj)
    out = [0] * n
    for source in range(n):
        available = adj[source]
        while available:
            bit = available & -available
            target = bit.bit_length() - 1
            out[target] |= 1 << source
            available -= bit
    return out


def boundary_state(adj: list[int]) -> tuple[list[int], list[int], list[list[int]]]:
    """Return start, end, and exposed-slot Q counts by subset DP."""
    n = len(adj)
    full = (1 << n) - 1
    end_dp = all_end_counts(adj)
    reverse_end_dp = all_end_counts(converse(adj))
    starts = list(reverse_end_dp[full])
    ends = list(end_dp[full])
    exposed = [[0] * n for _ in range(n)]

    for left_end in range(n):
        for right_start in range(n):
            if left_end == right_start:
                continue
            free = full ^ (1 << left_end) ^ (1 << right_start)
            subset = free
            total = 0
            while True:
                left = subset | (1 << left_end)
                right = full ^ left
                total += end_dp[left][left_end] * reverse_end_dp[right][right_start]
                if subset == 0:
                    break
                subset = (subset - 1) & free
            exposed[left_end][right_start] = total
    return starts, ends, exposed


def insertion_h(
    state: tuple[list[int], list[int], list[list[int]]], signature: int
) -> int:
    starts, ends, exposed = state
    n = len(starts)
    total = 0
    for vertex in range(n):
        total += starts[vertex] if (signature >> vertex) & 1 else ends[vertex]
    for left_end in range(n):
        if (signature >> left_end) & 1:
            continue
        for right_start in range(n):
            if (signature >> right_start) & 1:
                total += exposed[left_end][right_start]
    return total


def raw_code(adj: list[int]) -> int:
    code = 0
    bit_index = 0
    for source in range(len(adj)):
        for target in range(source + 1, len(adj)):
            code |= ((adj[source] >> target) & 1) << bit_index
            bit_index += 1
    return code


def decode_raw(code: int, n: int) -> list[int]:
    adj = [0] * n
    bit_index = 0
    for source in range(n):
        for target in range(source + 1, n):
            if (code >> bit_index) & 1:
                adj[source] |= 1 << target
            else:
                adj[target] |= 1 << source
            bit_index += 1
    return adj


def odd_intervals(values: set[int]) -> list[tuple[int, int, int]]:
    ordered = sorted(values)
    require(bool(ordered), "nonempty spectrum")
    intervals: list[tuple[int, int, int]] = []
    first = last = ordered[0]
    for value in ordered[1:]:
        if value == last + 2:
            last = value
        else:
            intervals.append((first, last, (last - first) // 2 + 1))
            first = last = value
    intervals.append((first, last, (last - first) // 2 + 1))
    return intervals


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def next_odd_prime(after: int) -> int:
    candidate = after + 1 + (after % 2)
    while not is_prime(candidate):
        candidate += 2
    return candidate


def read_historical_strong_values() -> tuple[set[int], int, int]:
    require(file_sha256(HISTOGRAM) == HISTOGRAM_SHA256, "historical histogram hash")
    values: set[int] = set()
    all_classes = 0
    strong_classes = 0
    lines = HISTOGRAM.read_text(encoding="ascii").splitlines()
    require(lines and lines[0] == "H\tcount_all\tcount_strong", "histogram header")
    for line in lines[1:]:
        value_text, all_text, strong_text = line.split("\t")
        value = int(value_text)
        all_count = int(all_text)
        strong_count = int(strong_text)
        all_classes += all_count
        strong_classes += strong_count
        if strong_count:
            values.add(value)
    return values, all_classes, strong_classes


def main() -> None:
    engine = load_engine()
    require(file_sha256(LEGACY_VALUES) == LEGACY_VALUES_SHA256, "legacy strong-values hash")

    representatives, counts = engine.generate(8)
    require(counts == A000568, "A000568 class counts through order eight")

    ear_values: set[int] = set()
    witnesses: dict[int, tuple[int, int, int]] = {}
    strong_parents = 0
    cuts_checked = 0
    for parent_index, parent in enumerate(representatives[8]):
        if not engine.is_strong(8, parent):
            continue
        strong_parents += 1
        state = boundary_state(parent)
        for signature in range(1, (1 << 8) - 1):
            cuts_checked += 1
            value = insertion_h(state, signature)
            ear_values.add(value)
            if value not in witnesses:
                child = engine.extend(parent, 8, signature)
                witnesses[value] = (parent_index, signature, raw_code(child))

    require(strong_parents == 6008, "strong order-eight parent count")
    require(cuts_checked == 1_526_032, "complete nonconstant ear count")
    require(len(ear_values) == 1482, "order-nine strong spectrum size")
    require(min(ear_values) == 75 and max(ear_values) == 3357, "spectrum extrema")

    historical, all_classes, historical_strong_classes = read_historical_strong_values()
    legacy = {int(value) for value in LEGACY_VALUES.read_text(encoding="ascii").split()}
    require(all_classes == 191_536, "historical all-class total")
    require(historical_strong_classes == 178_133, "historical strong-class total")
    require(ear_values == historical, "ear spectrum equals historical n=9 histogram")
    require(ear_values == legacy, "ear spectrum equals independent legacy value dump")

    intervals = odd_intervals(ear_values)
    require(intervals[:3] == [(75, 75, 1), (81, 81, 1), (85, 2881, 1399)],
            "solid order-nine interval")
    require(2883 not in ear_values, "first post-interval hostile")
    tail = sorted(value for value in ear_values if value > 2881)
    require(len(tail) == 81, "post-interval tail size")

    direct_failures = []
    witness_rows = {}
    for value, (parent_index, signature, code) in sorted(witnesses.items()):
        child = decode_raw(code, 9)
        direct_value = engine.Hcount(9, child)
        strong = engine.is_strong(9, child)
        if direct_value != value or not strong:
            direct_failures.append((value, direct_value, strong))
        if value in KEY_TARGETS:
            witness_rows[str(value)] = {
                "parent_index": parent_index,
                "signature": signature,
                "cut_weight": signature.bit_count(),
                "code": code,
            }
    require(not direct_failures, "direct DP/strongness check for every retained witness")
    require(set(witness_rows) == {str(value) for value in KEY_TARGETS}, "key witnesses")

    # The smaller unlock for 613 is useful because it is already an order-eight atom.
    witness_613_order8 = None
    for parent_index, parent in enumerate(representatives[7]):
        if not engine.is_strong(7, parent):
            continue
        state = boundary_state(parent)
        for signature in range(1, (1 << 7) - 1):
            if insertion_h(state, signature) == 613:
                child = engine.extend(parent, 7, signature)
                witness_613_order8 = {
                    "parent_index": parent_index,
                    "signature": signature,
                    "cut_weight": signature.bit_count(),
                    "code": raw_code(child),
                }
                require(engine.Hcount(8, child) == 613, "order-eight H=613 witness")
                require(engine.is_strong(8, child), "order-eight H=613 strongness")
                break
        if witness_613_order8 is not None:
            break
    require(witness_613_order8 is not None, "order-eight H=613 unlock")
    require(witness_613_order8["code"] == 251_585_423, "stable order-eight H=613 code")

    ordinary_missing = [
        value for value in range(85, 2882, 2) if is_prime(value) and value not in ear_values
    ]
    seven_missing = [
        prime for prime in range(13, 410, 2)
        if is_prime(prime) and 7 * prime not in ear_values
    ]
    require(not ordinary_missing, "ordinary-prime interval through 2879")
    require(not seven_missing, "seven-prime interval through p=409")
    next_ordinary = next_odd_prime(2881)
    next_seven_prime = next_odd_prime(409)
    require(next_ordinary == 2887, "next ordinary-prime target")
    require(next_seven_prime == 419 and 7 * next_seven_prime == 2933,
            "next seven-prime target")

    ledger = {
        "A000568": A000568,
        "strong_parents": strong_parents,
        "cuts_checked": cuts_checked,
        "spectrum": sorted(ear_values),
        "intervals": intervals,
        "tail": tail,
        "witness_rows": witness_rows,
        "witness_613_order8": witness_613_order8,
        "historical_class_totals": [all_classes, historical_strong_classes],
        "next_targets": [next_ordinary, 7 * next_seven_prime],
        "input_hashes": [ENGINE_SHA256, HISTOGRAM_SHA256, LEGACY_VALUES_SHA256],
    }
    semantic = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("ORDER-NINE STRONG-EAR SPECTRUM")
    print("A000568_counts=", [counts[n] for n in range(1, 9)])
    print("strong_order8_parents=", strong_parents)
    print("nonconstant_ears_checked=", cuts_checked)
    print("strong_spectrum_size=", len(ear_values), "min=", min(ear_values), "max=", max(ear_values))
    print("first_intervals=", intervals[:3])
    print("solid_odd_interval= [85,2881] count=1399")
    print("post_interval_tail_count=", len(tail))
    print("post_interval_tail=", tail)
    print("historical_all_classes=", all_classes, "historical_strong_classes=", historical_strong_classes)
    print("ear_histogram_legacy_sets_equal=", ear_values == historical == legacy)
    print("direct_witness_checks=", len(witnesses), "failures=", len(direct_failures))
    print("key_witnesses=", witness_rows)
    print("order8_H613_witness=", witness_613_order8)
    print("ordinary_prime_lane_through=2879 next_unforced=", next_ordinary)
    print("seven_prime_lane_p_through=409 next_unforced_value=", 7 * next_seven_prime)
    print("semantic_sha256=", semantic)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
