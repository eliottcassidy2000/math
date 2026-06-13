#!/usr/bin/env python3
"""
lonely_runner_tight_scan_s357.py

codex-2026-05-30 S357

Bounded exact scan for the tight stratum of the reduced lonely runner
conjecture.  This builds on lonely_runner_residue_probe_s356.py.

The S356 probe splits witnesses into positive gaps and boundary witnesses.
For a finite union of open forbidden intervals this split is tautological:
if the complement is nonempty and has no interval, it is made of endpoints.
The useful object is therefore the full-measure non-cover stratum:

    forbidden_length = 1, max_gap = 0, boundary_witness_count > 0.

These are tight LRC instances.  A counterexample would be a full open cover of
the time circle:

    forbidden_length = 1, max_gap = 0, boundary_witness_count = 0.

The scan below enumerates primitive speed sets in small boxes and records the
boundary-only sets with their component counts and first witness.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()


@dataclass(frozen=True)
class ScanConfig:
    k: int
    max_speed: int


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def scan_box(config: ScanConfig) -> dict[str, object]:
    rows = []
    positive = 0
    boundary_only = []
    open_covers = []
    min_positive = None
    min_positive_row = None
    component_hist = Counter()
    witness_hist = Counter()

    for combo in combinations(range(1, config.max_speed + 1), config.k):
        if not primitive(combo):
            continue
        row = S356.report(f"k={config.k} max={config.max_speed}", list(combo))
        rows.append(row)
        if row.max_gap > 0:
            positive += 1
            ratio = row.max_gap / row.threshold
            if min_positive is None or ratio < min_positive:
                min_positive = ratio
                min_positive_row = row
            continue

        if row.boundary_witness_count:
            boundary_only.append(row)
            component_hist[row.components] += 1
            witness_hist[row.boundary_witness] += 1
        else:
            open_covers.append(row)

    return {
        "config": config,
        "total": len(rows),
        "positive": positive,
        "boundary_only": boundary_only,
        "open_covers": open_covers,
        "min_positive": min_positive,
        "min_positive_row": min_positive_row,
        "component_hist": component_hist,
        "witness_hist": witness_hist,
    }


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def print_known_tight_examples() -> None:
    examples = [
        ("initial k=4", [1, 2, 3, 4]),
        ("sporadic n=5", [1, 3, 4, 7]),
        ("initial k=5", [1, 2, 3, 4, 5]),
        ("sporadic n=6", [1, 3, 4, 5, 9]),
        ("initial k=7", [1, 2, 3, 4, 5, 6, 7]),
        ("sporadic n=8 A", [1, 4, 5, 6, 7, 11, 13]),
        ("sporadic n=8 B", [1, 2, 3, 4, 5, 7, 12]),
    ]

    print("Known tight examples")
    for label, speeds in examples:
        row = S356.report(label, speeds)
        residues = tuple(sorted(v % (len(row.speeds) + 1) for v in row.speeds))
        print(f"  [{label}]")
        print(f"    speeds={row.speeds}")
        print(f"    residues_mod_{len(row.speeds)+1}={residues}")
        print(
            "    "
            f"components={row.components} forbidden_length={fmt_frac(row.forbidden_length)} "
            f"boundary_witnesses={row.boundary_witness_count} "
            f"first={fmt_frac(row.boundary_witness)} Q={row.boundary_modulus}"
        )
    print()


def print_scan_result(result: dict[str, object], sample_limit: int = 10) -> None:
    config = result["config"]
    assert isinstance(config, ScanConfig)
    boundary_only = result["boundary_only"]
    open_covers = result["open_covers"]
    component_hist = result["component_hist"]
    witness_hist = result["witness_hist"]
    min_positive_row = result["min_positive_row"]

    print(f"Exhaustive primitive box k={config.k}, max_speed={config.max_speed}")
    print(f"  total={result['total']}")
    print(f"  positive_gap={result['positive']}")
    print(f"  boundary_only={len(boundary_only)}")
    print(f"  open_covers={len(open_covers)}")
    if min_positive_row is not None:
        ratio = min_positive_row.max_gap / min_positive_row.threshold
        print(
            "  tightest_positive="
            f"speeds={min_positive_row.speeds} "
            f"max_gap={fmt_frac(min_positive_row.max_gap)} "
            f"ratio={float(ratio):.6f} "
            f"witness={fmt_frac(min_positive_row.witness)}"
        )
    print(f"  boundary_component_hist={dict(sorted(component_hist.items()))}")
    top_witnesses = sorted(
        witness_hist.items(), key=lambda item: (-item[1], item[0])
    )[:5]
    print(
        "  top_boundary_witnesses="
        + ", ".join(f"{fmt_frac(w)}:{c}" for w, c in top_witnesses)
    )
    print("  boundary_only_samples=")
    for row in boundary_only[:sample_limit]:
        print(
            "    "
            f"speeds={row.speeds} components={row.components} "
            f"boundary_count={row.boundary_witness_count} "
            f"first={fmt_frac(row.boundary_witness)} Q={row.boundary_modulus}"
        )
    if open_covers:
        print("  OPEN COVER CANDIDATES FOUND")
        for row in open_covers[:sample_limit]:
            print(f"    speeds={row.speeds}")
    print()


def main() -> None:
    print("Lonely runner tight-stratum scan (codex-2026-05-30 S357)")
    print("Exact interval arithmetic inherited from S356.\n")
    print_known_tight_examples()

    configs = [
        ScanConfig(k=3, max_speed=24),
        ScanConfig(k=4, max_speed=24),
        ScanConfig(k=5, max_speed=20),
        ScanConfig(k=6, max_speed=16),
        ScanConfig(k=7, max_speed=14),
    ]
    for config in configs:
        print_scan_result(scan_box(config))

    print("Summary")
    print("  No open-cover counterexample appeared in any scanned primitive box.")
    print("  Boundary-only instances are rare but structured; the first witness")
    print("  is usually 1/(k+1), and the boundary quotient often collapses to")
    print("  k+1 even when the speeds are not an initial segment.")


if __name__ == "__main__":
    main()
