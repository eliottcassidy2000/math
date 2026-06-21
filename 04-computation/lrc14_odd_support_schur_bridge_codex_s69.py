#!/usr/bin/env python3
"""Bridge HYP-2720 odd-support envelope with HYP-2719 Schur support-size lever.

HYP-2719 says the Weyl/carrier error is organized by Fourier relation support
size, with support-3 additive triangles as the leading cross-block lever.
HYP-2720 says that, after THM-561 factorial inversion, the origin atom has an
odd-support L1 envelope but not signed odd dominance.

This exact scout puts both diagnostics on the same split-row bank.  It is not a
proof; it checks whether factorial odd support is merely a shadow of Schur
triangle count, or whether the two ledgers should be kept as separate axes.
"""

from __future__ import annotations

import importlib.util
import itertools
import math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ODD_PATH = HERE / "lrc14_odd_support_weyl_error_codex_s69.py"
spec = importlib.util.spec_from_file_location("odd_support", ODD_PATH)
odd_support = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(odd_support)


BLOCK4 = (0, 1, 2, 3)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def row_ids(offset_blocks: tuple[tuple[int, tuple[int, ...]], ...]) -> dict[int, int]:
    ids = {0: -1}
    for idx, (offset, block) in enumerate(offset_blocks):
        for d in block:
            ids[offset + d] = idx
    return ids


def schur_counts(row: tuple[int, ...], ids: dict[int, int]) -> dict[str, int]:
    row_set = set(row)
    counts = Counter()
    for a in row:
        for b in row:
            c = a + b
            if c not in row_set:
                continue
            counts["all_ordered"] += 1
            if a > 0 and b > 0:
                counts["positive_ordered"] += 1
            if len({a, b, c}) == 3:
                counts["distinct_ordered"] += 1
            blocks = {ids.get(x, -2) for x in (a, b, c) if ids.get(x, -2) >= 0}
            if len(blocks) >= 2:
                counts["cross_ordered"] += 1
                if a > 0 and b > 0:
                    counts["cross_positive_ordered"] += 1
                if len({a, b, c}) == 3:
                    counts["cross_distinct_ordered"] += 1
            if a == 0 or b == 0:
                counts["apex_copy_ordered"] += 1
    return dict(counts)


def corr(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def report(name: str, offset_blocks: tuple[tuple[int, tuple[int, ...]], ...]) -> dict[str, object]:
    r = odd_support.row_report(name, offset_blocks)
    ids = row_ids(offset_blocks)
    counts = schur_counts(r["row"], ids)  # type: ignore[arg-type]
    r.update(counts)
    return r


def main() -> None:
    print("HYP-2719/HYP-2720 bridge: factorial odd support vs Schur support-size")
    print("Exact origin-atom diagnostics plus integer Schur/additive-triangle counts.\n")
    cases = [
        ("two 4-blocks, ratio 2:1", ((14, BLOCK4), (28, BLOCK4))),
        ("two 4-blocks, moderate gap", ((14, BLOCK4), (30, BLOCK4))),
        ("two 4-blocks, wider gap", ((30, BLOCK4), (80, BLOCK4))),
        ("two 4-blocks, high relation phase", ((15, BLOCK4), (31, BLOCK4))),
        ("two 4-blocks, positive Q0", ((40, BLOCK4), (61, BLOCK4))),
        ("two 4-blocks, positive Q0 high", ((80, BLOCK4), (121, BLOCK4))),
        ("5+3 split", ((20, (0, 1, 2, 3, 4)), (55, (0, 1, 2)))),
        (
            "3+3+2 split",
            ((18, (0, 1, 2)), (45, (0, 1, 2)), (90, (0, 1))),
        ),
        (
            "five 2-blocks",
            (
                (15, (0, 1)),
                (30, (0, 1)),
                (46, (0, 1)),
                (63, (0, 1)),
                (81, (0, 1)),
            ),
        ),
        (
            "seven singleton carriers",
            (
                (19, (0,)),
                (31, (0,)),
                (44, (0,)),
                (58, (0,)),
                (73, (0,)),
                (89, (0,)),
                (106, (0,)),
            ),
        ),
    ]
    rows = [report(name, blocks) for name, blocks in cases]
    print(
        "  row                                  odd_share  Q0          "
        "allSchur posSchur cross+ crossDistinct signedOdd"
    )
    for r in rows:
        print(
            f"  {r['name'][:34]:34s} "
            f"{float(r['odd_l1_share']):9.6f} "
            f"{float(r['q0']):+11.7f} "
            f"{r.get('all_ordered', 0):8d} "
            f"{r.get('positive_ordered', 0):8d} "
            f"{r.get('cross_positive_ordered', 0):6d} "
            f"{r.get('cross_distinct_ordered', 0):13d} "
            f"{str(r['signed_odd_dominates']):>9s}"
        )

    odd = [float(r["odd_l1_share"]) for r in rows]
    q0 = [abs(float(r["q0"])) for r in rows]
    all_s = [float(r.get("all_ordered", 0)) for r in rows]
    pos_s = [float(r.get("positive_ordered", 0)) for r in rows]
    cross_pos = [float(r.get("cross_positive_ordered", 0)) for r in rows]
    cross_dist = [float(r.get("cross_distinct_ordered", 0)) for r in rows]
    print("\nCORRELATIONS (Pearson on this small exact bank)")
    print(f"  corr(odd_share, all_ordered_schur)      = {corr(odd, all_s):+.6f}")
    print(f"  corr(odd_share, positive_schur)         = {corr(odd, pos_s):+.6f}")
    print(f"  corr(odd_share, cross_positive_schur)   = {corr(odd, cross_pos):+.6f}")
    print(f"  corr(odd_share, cross_distinct_schur)   = {corr(odd, cross_dist):+.6f}")
    print(f"  corr(|Q0|, cross_positive_schur)        = {corr(q0, cross_pos):+.6f}")

    print("\nTOURNAMENT ANALYSIS")
    print("  vertices: rows; observable: odd_share then cross-positive Schur count")
    print("  switch/gauge: compare factorial support with Fourier support-size proxy")
    scores = Counter()
    edges = set()
    for i, j in itertools.combinations(range(len(rows)), 2):
        ai = (rows[i]["odd_l1_share"], rows[i].get("cross_positive_ordered", 0), rows[i]["q0"])
        aj = (rows[j]["odd_l1_share"], rows[j].get("cross_positive_ordered", 0), rows[j]["q0"])
        if ai >= aj:
            edges.add((i, j))
            scores[i] += 1
        else:
            edges.add((j, i))
            scores[j] += 1
    cycles = 0
    for a, b, c in itertools.combinations(range(len(rows)), 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    print(f"  score_hist={dict(sorted(Counter(scores[i] for i in range(len(rows))).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  path:")
    for i in sorted(
        range(len(rows)),
        key=lambda idx: (
            rows[idx]["odd_l1_share"],
            rows[idx].get("cross_positive_ordered", 0),
            rows[idx]["q0"],
        ),
        reverse=True,
    ):
        r = rows[i]
        print(
            f"    {r['name']}: odd_share={fmt(r['odd_l1_share'])} "
            f"cross+={r.get('cross_positive_ordered', 0)} Q0={fmt(r['q0'])}"
        )

    print("\nSYNTHESIS")
    print("  The two support notions are not the same.  Schur/additive triangles are")
    print("  a Fourier relation-support proxy; odd/even W_j is the factorial-origin")
    print("  support after miss-zeta inversion.  In this bank, odd_share is not a")
    print("  monotone function of cross-Schur counts, so HYP-2720 should not replace")
    print("  HYP-2719.  The useful split is sequential: HYP-2719 selects low-support")
    print("  carrier packets; HYP-2720 taxes their origin atom in factorial parity.")


if __name__ == "__main__":
    main()

