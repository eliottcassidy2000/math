#!/usr/bin/env python3
"""HYP-2721/S70: q0 origin atom versus the Vitali atom cone.

The tournament Vitali atom is a low-observable-preserving move whose target
functional still changes through a hidden packet channel.  This scout tests the
same pattern for the LRC14 carrier-product gap:

    Q_0 = ProductCover - ActualCover = sum_j (-1)^j W_j.

Here W_j are the factorial miss-zeta discrepancies and Q_t are the missed-count
atom discrepancies.  The universal finite-difference basis

    B_j(t) = (-1)^(j-t) binom(j,t), 0 <= t <= j,

is the atom move created by a unit W_j.  It preserves all lower factorial
moments W_i for i<j but moves the origin atom by (-1)^j.  That is the abstract
"Vitali atom cone" for the q0 target.
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


def fmt(q: F | int) -> str:
    if isinstance(q, int):
        return str(q)
    return f"{q} ({float(q):.9f})"


def cases() -> list[tuple[str, tuple[tuple[int, tuple[int, ...]], ...]]]:
    return [
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


def finite_difference_basis() -> list[tuple[int, ...]]:
    basis = []
    for j in range(7):
        move = []
        for t in range(7):
            if t <= j:
                move.append(((-1) ** (j - t)) * math.comb(j, t))
            else:
                move.append(0)
        basis.append(tuple(move))
    return basis


def low_eval(W: tuple[F, ...], r: int) -> F:
    return sum(((-1) ** j) * W[j] for j in range(r + 1))


def tail_eval(W: tuple[F, ...], r: int) -> F:
    return sum(((-1) ** j) * W[j] for j in range(r + 1, 7))


def low_l1(W: tuple[F, ...], r: int) -> F:
    return sum(abs(W[j]) for j in range(r + 1))


def tail_l1(W: tuple[F, ...], r: int) -> F:
    return sum(abs(W[j]) for j in range(r + 1, 7))


def enrich(name: str, offset_blocks: tuple[tuple[int, tuple[int, ...]], ...]) -> dict[str, object]:
    r = odd_support.row_report(name, offset_blocks)
    ids = row_ids(offset_blocks)
    r.update(schur_counts(r["row"], ids))  # type: ignore[arg-type]
    W = r["W"]
    atom = r["atom"]
    q0 = atom[0]
    atom_l1 = sum(abs(x) for x in atom)
    nonorigin_l1 = sum(abs(x) for x in atom[1:])
    r.update(
        {
            "atom_l1": atom_l1,
            "nonorigin_l1": nonorigin_l1,
            "origin_share": abs(q0) / atom_l1 if atom_l1 else F(0),
            "tax_ratio": nonorigin_l1 / abs(q0) if q0 else None,
            "low2_eval": low_eval(W, 2),
            "tail2_eval": tail_eval(W, 2),
            "low2_l1": low_l1(W, 2),
            "tail2_l1": tail_l1(W, 2),
            "tail2_pressure": tail_l1(W, 2) / abs(q0) if q0 else None,
            "low3_eval": low_eval(W, 3),
            "tail3_eval": tail_eval(W, 3),
            "low3_l1": low_l1(W, 3),
            "tail3_l1": tail_l1(W, 3),
            "tail3_pressure": tail_l1(W, 3) / abs(q0) if q0 else None,
            "cone_margin": nonorigin_l1 - abs(q0),
        }
    )
    assert r["low2_eval"] + r["tail2_eval"] == q0
    assert r["low3_eval"] + r["tail3_eval"] == q0
    return r


def print_basis() -> None:
    print("UNIVERSAL FACTORIAL-ATOM VITALI MOVES")
    print("  B_j is the atom profile created by unit W_j; it preserves W_i for i<j.")
    print("  columns are atom coordinates Q_0..Q_6")
    for j, move in enumerate(finite_difference_basis()):
        print(f"  B_{j}: {move}  q0={move[0]}")


def print_row_table(rows: list[dict[str, object]]) -> None:
    print("\nROW-LEVEL Q0 / VITALI-CONE DIAGNOSTICS")
    print(
        "  row                                  Q0          origin%  nonorigin/|Q0| "
        "tail2L1/|Q0| tail3L1/|Q0| cross+"
    )
    for r in rows:
        tax = r["tax_ratio"]
        tail2 = r["tail2_pressure"]
        tail3 = r["tail3_pressure"]
        print(
            f"  {str(r['name'])[:34]:34s} "
            f"{float(r['q0']):+11.7f} "
            f"{float(r['origin_share']):8.4f} "
            f"{float(tax) if tax is not None else 0:15.4f} "
            f"{float(tail2) if tail2 is not None else 0:13.4f} "
            f"{float(tail3) if tail3 is not None else 0:13.4f} "
            f"{r.get('cross_positive_ordered', 0):6d}"
        )
        print(
            f"      low2_eval={fmt(r['low2_eval'])}; tail2_eval={fmt(r['tail2_eval'])}; "
            f"cone_margin={fmt(r['cone_margin'])}"
        )


def aggregate(rows: list[dict[str, object]]) -> None:
    print("\nAGGREGATE")
    q0_l1 = sum(abs(r["q0"]) for r in rows)  # type: ignore[arg-type]
    atom_l1 = sum(r["atom_l1"] for r in rows)  # type: ignore[arg-type]
    nonorigin_l1 = sum(r["nonorigin_l1"] for r in rows)  # type: ignore[arg-type]
    tail2_l1 = sum(r["tail2_l1"] for r in rows)  # type: ignore[arg-type]
    tail3_l1 = sum(r["tail3_l1"] for r in rows)  # type: ignore[arg-type]
    low2_abs_eval = sum(abs(r["low2_eval"]) for r in rows)  # type: ignore[arg-type]
    tail2_abs_eval = sum(abs(r["tail2_eval"]) for r in rows)  # type: ignore[arg-type]
    print(f"  sum |Q0|={fmt(q0_l1)}")
    print(f"  sum atom_L1={fmt(atom_l1)}")
    print(f"  sum nonorigin_L1={fmt(nonorigin_l1)}")
    print(f"  aggregate origin share={fmt(q0_l1 / atom_l1)}")
    print(f"  aggregate nonorigin tax ratio={fmt(nonorigin_l1 / q0_l1)}")
    print(f"  sum tail2_L1={fmt(tail2_l1)}; tail2_L1/sum|Q0|={fmt(tail2_l1 / q0_l1)}")
    print(f"  sum tail3_L1={fmt(tail3_l1)}; tail3_L1/sum|Q0|={fmt(tail3_l1 / q0_l1)}")
    print(f"  sum |low2_eval|={fmt(low2_abs_eval)}")
    print(f"  sum |tail2_eval|={fmt(tail2_abs_eval)}")
    print(
        "  cone_margin positive rows="
        f"{sum(1 for r in rows if r['cone_margin'] > 0)}/{len(rows)}"
    )

    q0 = [abs(float(r["q0"])) for r in rows]
    origin_share = [float(r["origin_share"]) for r in rows]
    tail2_pressure = [float(r["tail2_pressure"]) for r in rows]
    cross_pos = [float(r.get("cross_positive_ordered", 0)) for r in rows]
    print("\nCORRELATIONS (small exact bank)")
    print(f"  corr(|Q0|, cross_positive_schur)       = {corr(q0, cross_pos):+.6f}")
    print(f"  corr(origin_share, cross_positive)     = {corr(origin_share, cross_pos):+.6f}")
    print(f"  corr(tail2_pressure, cross_positive)   = {corr(tail2_pressure, cross_pos):+.6f}")


def tournament(rows: list[dict[str, object]]) -> None:
    print("\nTOURNAMENT ANALYSIS")
    print("  vertices: split rows")
    print("  pairwise observable: tail2_L1/|Q0|, then nonorigin atom tax, then support-3 Schur proxy")
    print("  switch/gauge: low-factorial observer -> full atom cone -> relation-support proxy")
    scores = Counter()
    edges = set()
    for i, j in itertools.combinations(range(len(rows)), 2):
        ai = (
            rows[i]["tail2_pressure"] or F(0),
            rows[i]["tax_ratio"] or F(0),
            rows[i].get("cross_positive_ordered", 0),
        )
        aj = (
            rows[j]["tail2_pressure"] or F(0),
            rows[j]["tax_ratio"] or F(0),
            rows[j].get("cross_positive_ordered", 0),
        )
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
    print("  Hamiltonian pressure path:")
    for i in sorted(
        range(len(rows)),
        key=lambda idx: (
            rows[idx]["tail2_pressure"] or F(0),
            rows[idx]["tax_ratio"] or F(0),
            rows[idx].get("cross_positive_ordered", 0),
        ),
        reverse=True,
    ):
        r = rows[i]
        print(
            f"    {r['name']}: tail2_pressure={fmt(r['tail2_pressure'])} "
            f"tax={fmt(r['tax_ratio'])} cross+={r.get('cross_positive_ordered', 0)} "
            f"Q0={fmt(r['q0'])}"
        )


def synthesis() -> None:
    print("\nSYNTHESIS")
    print("  The q0 origin atom is Vitali-like only after the factorial transform.")
    print("  Low factorial moments are the analogue of the lambda/c3/c5-visible")
    print("  data: they can be held fixed by higher finite-difference moves B_j,")
    print("  while Q0 still changes.  The observed split rows all have positive")
    print("  nonorigin atom-cone margin, so bounding the whole missed-count law is")
    print("  stronger than necessary.  HYP-2719 should choose the support-size")
    print("  carrier packets; HYP-2720/HYP-2721 should tax those packets in the")
    print("  factorial-origin basis before scalarizing to Q0.")


def main() -> None:
    print("HYP-2721/S70 q0 origin atom versus Vitali atom cone")
    print("Exact arithmetic inherited from the S69 miss-zeta engine; no proof claimed.\n")
    print_basis()
    rows = [enrich(name, blocks) for name, blocks in cases()]
    print_row_table(rows)
    aggregate(rows)
    tournament(rows)
    synthesis()


if __name__ == "__main__":
    main()
