#!/usr/bin/env python3
"""Exact incremental profile subtraction for one physical first-label spike.

This imports the frozen bridge mechanics, reconstructs the exact THM-2941
physical rows, and classifies only the ``(body,d1)`` groups occurring at the
requested ``z1``.  It is intended for fast reconciliation when a new scalar
spike closure lands: the reported stage counts and profile occurrences can be
subtracted from a previously frozen full-ledger summary without recomputing
all 5,338 conditional-GF pairs.

The script does not prove that the requested spike is empty.  That remains an
independently pinned scalar/ray/status input.
"""

from __future__ import annotations

import argparse
from collections import defaultdict

import bridge


EXPECTED_SPIKE_ROWS = {
    378: 9,
    364: 25,
    350: 53,
    336: 8,
}


def main(z1: int, workers: int) -> None:
    bridge.require(z1 >= 15, "the external first label must be at least 15")
    _profiles, rows = bridge.reconstruct_physical_rows(workers)
    spike_rows = tuple(row for row in rows if row[2] == z1)
    expected = EXPECTED_SPIKE_ROWS.get(z1)
    if expected is not None:
        bridge.require(len(spike_rows) == expected, ("spike row count changed", z1))
    bridge.require(spike_rows, ("empty physical spike", z1))

    support_rows, body_divisor_rows, support_body_count = bridge.build_support_rows()
    grouped_rows = defaultdict(list)
    for row in spike_rows:
        grouped_rows[(row[0], row[3])].append(row)
    classifications, pair_count, coefficient_semantic = bridge.classify_all_groups(
        grouped_rows, support_rows
    )
    summary = bridge.summarize(spike_rows, classifications, set())

    bridge.require_inputs_unchanged("spike_delta_end")
    lines = [
        "LRC14 k=3 physical denominator spike delta",
        *[
            f"{label}_sha256={bridge.INPUT_SHA256[label]}"
            for label in sorted(bridge.INPUT_SHA256)
        ],
        f"z1={z1}",
        f"physical_spike_rows={len(spike_rows)}",
        f"physical_body_d1_groups={len(grouped_rows)}",
        f"body_divisor_rows={body_divisor_rows}",
        f"support_body_count={support_body_count}",
        f"conditional_GF_pairs={pair_count}",
        f"conditional_coefficient_semantic_sha256={coefficient_semantic}",
        *bridge.render_summary(f"spike_z{z1}", summary),
        "scope=exact denominator-screen delta only; spike emptiness is not proved here",
        "input_snapshot_guard=start/after_module_load/spike_delta_end:PASS",
        "all_exact_controls=PASS",
    ]
    print("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--z1", type=int, required=True)
    parser.add_argument("--workers", type=int, default=min(8, bridge.mp.cpu_count() or 1))
    arguments = parser.parse_args()
    main(arguments.z1, arguments.workers)
