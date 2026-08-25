#!/usr/bin/env python3
"""Independent audit for THM-4097.

This path does not import the strong-ear compiler.  It compares two historical
order-nine censuses, characterizes their common strong spectrum, and verifies
selected labelled witnesses by both Held--Karp DP and literal permutation
enumeration.
"""

from __future__ import annotations

import json
from hashlib import sha256
from itertools import permutations
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
HISTOGRAM = ROOT / "05-knowledge/results/h_spectrum_n9_histogram_monad_s6.tsv"
LEGACY_VALUES = ROOT / "05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out"

HISTOGRAM_SHA256 = "e7d5594879d4c3af739cb94ca8cfd944879c4d586747d993dd6687e60126552f"
LEGACY_VALUES_SHA256 = "27fbef5b06fcf0369eeb602e513c3802ea171492e1292a3f6afa3efeadef9f55"

WITNESSES = {
    75: (9, 9_897_508_671),
    81: (9, 10_165_944_127),
    85: (9, 49_392_123_839),
    613: (8, 251_585_423),
    623: (9, 63_853_559_615),
    2881: (9, 34_308_081_455),
    2885: (9, 63_870_568_207),
    3357: (9, 68_164_491_031),
}


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


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
    require(code >> bit_index == 0, f"code fits order {n}")
    return adj


def h_dp(adj: list[int]) -> int:
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
    return sum(dp[-1])


def h_literal(adj: list[int]) -> int:
    n = len(adj)
    total = 0
    for order in permutations(range(n)):
        if all((adj[order[index]] >> order[index + 1]) & 1 for index in range(n - 1)):
            total += 1
    return total


def reachable(adj: list[int], reverse: bool) -> int:
    n = len(adj)
    seen = 1
    frontier = 1
    while frontier:
        nxt = 0
        for source in range(n):
            if not ((frontier >> source) & 1):
                continue
            if reverse:
                for target in range(n):
                    if (adj[target] >> source) & 1:
                        nxt |= 1 << target
            else:
                nxt |= adj[source]
        nxt &= ~seen
        seen |= nxt
        frontier = nxt
    return seen


def is_strong(adj: list[int]) -> bool:
    full = (1 << len(adj)) - 1
    return reachable(adj, False) == full and reachable(adj, True) == full


def odd_intervals(values: set[int]) -> list[tuple[int, int, int]]:
    ordered = sorted(values)
    first = last = ordered[0]
    out: list[tuple[int, int, int]] = []
    for value in ordered[1:]:
        if value == last + 2:
            last = value
        else:
            out.append((first, last, (last - first) // 2 + 1))
            first = last = value
    out.append((first, last, (last - first) // 2 + 1))
    return out


def main() -> None:
    require(file_sha256(HISTOGRAM) == HISTOGRAM_SHA256, "histogram hash")
    require(file_sha256(LEGACY_VALUES) == LEGACY_VALUES_SHA256, "legacy-values hash")

    histogram_values: set[int] = set()
    all_classes = 0
    strong_classes = 0
    lines = HISTOGRAM.read_text(encoding="ascii").splitlines()
    require(lines[0] == "H\tcount_all\tcount_strong", "histogram header")
    for line in lines[1:]:
        value_text, all_text, strong_text = line.split("\t")
        value = int(value_text)
        all_count = int(all_text)
        strong_count = int(strong_text)
        all_classes += all_count
        strong_classes += strong_count
        if strong_count:
            histogram_values.add(value)
    legacy_values = {int(value) for value in LEGACY_VALUES.read_text(encoding="ascii").split()}

    require(histogram_values == legacy_values, "independent historical strong spectra agree")
    require(all_classes == 191_536 and strong_classes == 178_133, "class totals")
    require(len(histogram_values) == 1482, "strong value count")
    intervals = odd_intervals(histogram_values)
    require(intervals[:3] == [(75, 75, 1), (81, 81, 1), (85, 2881, 1399)],
            "solid interval characterization")
    require(2883 not in histogram_values and 2885 in histogram_values,
            "first post-interval boundary")

    witness_rows = []
    for expected_h, (n, code) in WITNESSES.items():
        adj = decode_raw(code, n)
        dp_value = h_dp(adj)
        literal_value = h_literal(adj)
        strong = is_strong(adj)
        require(dp_value == expected_h, f"DP witness H={expected_h}")
        require(literal_value == expected_h, f"literal witness H={expected_h}")
        require(strong, f"strong witness H={expected_h}")
        witness_rows.append((expected_h, n, code, dp_value, literal_value, strong))

    ledger = {
        "input_hashes": [HISTOGRAM_SHA256, LEGACY_VALUES_SHA256],
        "class_totals": [all_classes, strong_classes],
        "strong_values": sorted(histogram_values),
        "intervals": intervals,
        "witness_rows": witness_rows,
    }
    semantic = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("THM-4097 INDEPENDENT AUDIT")
    print("historical_sets_equal=", histogram_values == legacy_values)
    print("class_totals=", (all_classes, strong_classes))
    print("strong_spectrum_size=", len(histogram_values))
    print("first_intervals=", intervals[:3])
    print("first_post_interval_boundary=", (2883 in histogram_values, 2885 in histogram_values))
    for row in witness_rows:
        print("witness=", row)
    print("semantic_sha256=", semantic)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
