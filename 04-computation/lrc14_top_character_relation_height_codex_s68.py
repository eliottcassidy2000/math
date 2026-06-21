#!/usr/bin/env python3
"""HYP-2717 relation-height scout for the LRC14 multi-block gap.

This is a small exact addendum to the miss-zeta/top-character route.  It scans
two coherent 4-block carriers and compares:

* raw carrier separation,
* the primitive exact relation height for (M1,M2),
* the exact product-vs-actual cover error,
* the share of cap-product slack consumed by that error.

No proof is claimed.  The goal is to check the HYP-2717 warning that raw gap is
not the right variable: exact integer relation height is the carrier datum that
survives the line integral.
"""

from __future__ import annotations

import importlib.util
import itertools
from collections import Counter
from fractions import Fraction as F
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent
MISS_ZETA_PATH = HERE / "lrc14_multiblock_miss_zeta_layers_codex_20260621.py"
spec = importlib.util.spec_from_file_location("miss_zeta_layers", MISS_ZETA_PATH)
miss_zeta = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(miss_zeta)


CAP_9 = F(1979, 4004)
BLOCK = (0, 1, 2, 3)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def relation_vector(m1: int, m2: int) -> tuple[int, int]:
    g = gcd(m1, m2)
    return (m2 // g, -m1 // g)


def relation_height(m1: int, m2: int) -> int:
    a, b = relation_vector(m1, m2)
    return abs(a) + abs(b)


def row_report(m1: int, m2: int) -> dict[str, object]:
    offsets = ((m1, BLOCK), (m2, BLOCK))
    row = miss_zeta.row_from_blocks(offsets)
    actual_z = miss_zeta.zeta_from_hit_profile(miss_zeta.actual_hit_profile(row))
    product_z = miss_zeta.product_miss_zeta_shared_x((BLOCK, BLOCK))
    actual = miss_zeta.cover_from_zeta(actual_z)
    product = miss_zeta.cover_from_zeta(product_z)
    error = product - actual
    cap_product = CAP_9 - product
    ratio = abs(error) / cap_product
    h = relation_height(m1, m2)
    return {
        "offsets": (m1, m2),
        "row": row,
        "gap": m2 - (m1 + max(BLOCK)),
        "relation": relation_vector(m1, m2),
        "height": h,
        "actual": actual,
        "product": product,
        "error": error,
        "cap_product": cap_product,
        "ratio": ratio,
        "height_error": abs(error) * h,
        "primitive": miss_zeta.primitive(row),
    }


def print_report(report: dict[str, object]) -> None:
    print(f"offsets={report['offsets']} row={report['row']}")
    print(
        f"  raw_gap={report['gap']} relation={report['relation']} "
        f"height={report['height']} primitive={report['primitive']}"
    )
    print(f"  actual={fmt(report['actual'])}")
    print(f"  product={fmt(report['product'])}")
    print(f"  product-actual={fmt(report['error'])}")
    print(f"  cap-product={fmt(report['cap_product'])}")
    print(f"  |error|/(cap-product)={fmt(report['ratio'])}")
    print(f"  height*|error|={fmt(report['height_error'])}")


def tournament(reports: list[dict[str, object]]) -> None:
    print("\nTOURNAMENT ANALYSIS")
    print("  vertices: tested two-carrier rows")
    print("  pairwise observable: larger slack-error ratio is more dangerous;")
    print("    ties broken by lower relation height, then lower raw gap")
    print("  switch/gauge: relation-height carrier basis after top-character quotient")
    scores = Counter()
    edges = set()
    for i, j in itertools.combinations(range(len(reports)), 2):
        ai = (reports[i]["ratio"], -reports[i]["height"], -reports[i]["gap"])
        aj = (reports[j]["ratio"], -reports[j]["height"], -reports[j]["gap"])
        if ai >= aj:
            edges.add((i, j))
            scores[i] += 1
        else:
            edges.add((j, i))
            scores[j] += 1
    cycles = 0
    for a, b, c in itertools.combinations(range(len(reports)), 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    score_hist = Counter(scores[i] for i in range(len(reports)))
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  proof-pressure path:")
    order = sorted(
        range(len(reports)),
        key=lambda i: (reports[i]["ratio"], -reports[i]["height"], -reports[i]["gap"]),
        reverse=True,
    )
    for i in order:
        r = reports[i]
        print(
            f"    offsets={r['offsets']} H={r['height']} raw_gap={r['gap']} "
            f"ratio={fmt(r['ratio'])} error={fmt(r['error'])}"
        )


def main() -> None:
    print("HYP-2717 relation-height scout for two 4-block carriers")
    print("Exact arithmetic via the HYP-2715/HYP-2716 miss-zeta engine.\n")
    cases = [
        (14, 28),  # large raw gap but tiny relation height
        (14, 30),
        (30, 60),  # much larger raw gap, same tiny relation height
        (30, 80),
        (15, 31),
        (20, 43),
        (30, 47),
        (40, 61),
        (60, 89),
        (80, 121),
    ]
    reports = [row_report(*case) for case in cases]
    for report in reports:
        print_report(report)
    tournament(reports)
    print("\nSYNTHESIS")
    print("  Raw separation is not the right carrier variable: rows with large gaps")
    print("  can still have very low exact relation height, e.g. offsets (30,60).")
    print("  The cap-budget ratios remain small in this bank, but they organize by")
    print("  relation-height and arithmetic phase more naturally than by raw gap.")
    print("  This supports the HYP-2717 proof split: finite-ledger low-height exact")
    print("  relations first, then a top-character Fourier tail for high-height rows.")


if __name__ == "__main__":
    main()

