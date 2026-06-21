#!/usr/bin/env python3
"""Forced-zero/additive-energy tie-break audit for LRC14 cover rows.

This is a compact follow-up to HYP-2735.  Additive energy is a strong selector
for structured cover-risk rows, but it is not a monotone scalar proof.  The
useful next lemma is narrower: once additive energy places a row in the
arithmetic-progression stratum, the LRC convention (0 forced, or positive block
anchored at 1) supplies the missing offset/span tie-break.

Tournament Analysis:
  vertices: scalar proof lenses on a curated exact row bank
  observable: pairwise agreement with exact p0=measS7 risk ordering
  switch/gauge: fewer inversions means the scalar lens better preserves cover risk
  tie Hamiltonian path: lenses sorted by (inversions, tied pairs, label)
"""

from __future__ import annotations

import itertools
from collections import Counter
from fractions import Fraction
from math import gcd

P = 7


def breakpoints(row: tuple[int, ...]) -> list[Fraction]:
    pts = {Fraction(0), Fraction(1)}
    for e in row:
        if e == 0:
            continue
        for t in range(P * e + 1):
            pts.add(Fraction(t, P * e))
    return sorted(pts)


def meas_s7(row: tuple[int, ...]) -> Fraction:
    pts = breakpoints(row)
    total = Fraction(0)
    nonzero = [e for e in row if e != 0]
    for a, b in zip(pts, pts[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        sectors = {0} if 0 in row else set()
        for e in nonzero:
            sectors.add(int(P * ((e * mid) % 1)))
            if len(sectors) == P:
                total += b - a
                break
    return total


def additive_energy(row: tuple[int, ...]) -> int:
    sums: Counter[int] = Counter(a + b for a in row for b in row)
    return sum(v * v for v in sums.values())


def sumset_size(row: tuple[int, ...]) -> int:
    return len({a + b for a in row for b in row})


def ordered_schur(row: tuple[int, ...]) -> int:
    rowset = set(row)
    return sum(1 for a in row for b in row if a + b in rowset)


def ap3_count(row: tuple[int, ...]) -> int:
    rowset = set(row)
    return sum(
        1
        for a, c in itertools.combinations(row, 2)
        if (a + c) % 2 == 0 and (a + c) // 2 in rowset
    )


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for e in row:
        if e:
            g = gcd(g, e)
    return g == 1


def fmt(fr: Fraction) -> str:
    return f"{fr} ({float(fr):.6f})"


def row_stats(label: str, row: tuple[int, ...]) -> dict[str, object]:
    return {
        "label": label,
        "row": row,
        "p0": meas_s7(row),
        "energy": additive_energy(row),
        "neg_sumset": -sumset_size(row),
        "neg_span": -(row[-1] - row[0]),
        "schur": ordered_schur(row),
        "ap3": ap3_count(row),
    }


def ap_stratum_audit() -> None:
    print("=" * 88)
    print("AP stratum: forced-zero dilation equality and positive-offset tie-break")
    print("=" * 88)
    for k in (8, 9, 10):
        forced = tuple(range(k))
        base = row_stats(f"forced C{k}", forced)
        print(f"k={k} forced-zero consecutive {forced}: p0={fmt(base['p0'])}, AE={base['energy']}")
        for d in (2, 3, 5):
            dil = tuple(d * i for i in range(k))
            cur = row_stats(f"{d}*forced C{k}", dil)
            print(
                f"  dilation d={d:<2} p0={fmt(cur['p0'])} "
                f"AE={cur['energy']} equal_p0={cur['p0'] == base['p0']}"
            )

    print("\nPositive APs with k=8 and equal additive energy:")
    positive_rows = [
        ("anchored [1..8]", tuple(range(1, 9))),
        ("dilation [2,4,..,16]", tuple(range(2, 17, 2))),
        ("shifted odd [1,3,..,15]", tuple(range(1, 16, 2))),
        ("shifted [3,5,..,17]", tuple(range(3, 18, 2))),
    ]
    ref = row_stats(*positive_rows[0])
    for label, row in positive_rows:
        cur = row_stats(label, row)
        print(
            f"  {label:22s} p0={fmt(cur['p0'])} AE={cur['energy']} "
            f"same_AE_as_anchor={cur['energy'] == ref['energy']}"
        )
    print("  conclusion: additive energy finds the AP stratum; offset/dilation data is still needed.")

    print("\nStrict additive-energy inversion from the current HYP-2735 landscape:")
    high = (14, 15, 23, 27, 28, 32, 33, 34)
    low = (12, 14, 15, 20, 24, 28, 30, 32)
    for label, row in [("higher AE / lower p0", high), ("lower AE / higher p0", low)]:
        cur = row_stats(label, row)
        print(f"  {label:22s} row={row} p0={fmt(cur['p0'])} AE={cur['energy']}")


def inversion_count(rows: list[dict[str, object]], lens: str) -> tuple[int, int]:
    inversions = 0
    tied_pairs = 0
    for a, b in itertools.combinations(rows, 2):
        score_a = a[lens]
        score_b = b[lens]
        if score_a == score_b:
            tied_pairs += 1
            continue
        predicts_a = score_a > score_b
        actual_a = a["p0"] > b["p0"]
        if predicts_a != actual_a:
            inversions += 1
    return inversions, tied_pairs


def scalar_lens_tournament() -> None:
    bank = [
        ("forced C8", tuple(range(8))),
        ("positive C8", tuple(range(1, 9))),
        ("positive AP d2", tuple(range(2, 17, 2))),
        ("shifted odd AP", tuple(range(1, 16, 2))),
        ("dyadic", (0, 1, 2, 4, 8, 16, 32, 64)),
        ("two-block", (0, 1, 2, 3, 40, 41, 42, 43)),
        ("Mian-Chowla", (0, 1, 3, 7, 12, 20, 30, 44)),
        ("HYP-2735 inv high", (14, 15, 23, 27, 28, 32, 33, 34)),
        ("HYP-2735 inv low", (12, 14, 15, 20, 24, 28, 30, 32)),
    ]
    rows = [row_stats(label, row) for label, row in bank]
    print("\n" + "=" * 88)
    print("Curated exact row bank")
    print("=" * 88)
    for r in sorted(rows, key=lambda x: x["p0"], reverse=True):
        print(
            f"  {r['label']:18s} p0={float(r['p0']):.6f} "
            f"AE={r['energy']:4d} sumset={-r['neg_sumset']:3d} "
            f"span={-r['neg_span']:3d} schur={r['schur']:3d} ap3={r['ap3']:3d}"
        )

    lenses = ["energy", "neg_sumset", "neg_span", "schur", "ap3"]
    stats = [(lens, *inversion_count(rows, lens)) for lens in lenses]
    stats.sort(key=lambda t: (t[1], t[2], t[0]))
    beats = {
        (a, b): (dict((x, (i, ties)) for x, i, ties in stats)[a], a)
        < (dict((x, (i, ties)) for x, i, ties in stats)[b], b)
        for a, b in itertools.permutations(lenses, 2)
    }
    cycles = 0
    for a, b, c in itertools.combinations(lenses, 3):
        if beats[(a, b)] and beats[(b, c)] and beats[(c, a)]:
            cycles += 1
        if beats[(a, c)] and beats[(c, b)] and beats[(b, a)]:
            cycles += 1
    hp_count = sum(
        1
        for perm in itertools.permutations(lenses)
        if all(beats[(perm[i], perm[i + 1])] for i in range(len(perm) - 1))
    )

    print("\nTournament Analysis over scalar lenses:")
    for lens, inv, tied in stats:
        print(f"  {lens:10s} inversions={inv:2d} tied_pairs={tied:2d}")
    print(f"  directed_3_cycles={cycles}  hamiltonian_path_count={hp_count}")
    print("  tie Hamiltonian path: " + " > ".join(lens for lens, _, _ in stats))


def main() -> None:
    ap_stratum_audit()
    scalar_lens_tournament()


if __name__ == "__main__":
    main()
