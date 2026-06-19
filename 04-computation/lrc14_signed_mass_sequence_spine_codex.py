#!/usr/bin/env python3
"""
LRC(14) signed-mass sequence spine.

The current support-6 residual has a striking numerical signature:
large absolute mass on boundary/cusp faces, but tiny signed mass after the
seven-sector kernel is applied.  This script collects the recurring integer
and fractional sequences that make that signature proof-useful instead of
anecdotal.

Inputs:
  * THM-538/S12 kernel code:
      04-computation/lrc14_support6_residue_cusp_codex_s12.py
  * S12 exact support-6 shell ledger:
      05-knowledge/results/lrc14_support6_residue_cusp_codex_s12.out
  * HYP-2598/HYP-2608 sequence data, repeated here as exact fractions.

Tournament Analysis declaration.
  Vertex set: proof-obligation sequence families, not runners:
    support_floor, residue_constant_decay, cusp_signed_ratio,
    empty_window_degree_drop, universal_center_survivors,
    low_height_wall_ledger, coimage_reciprocal_tail.
  Pairwise observable: (exactness, proof leverage, portability in k/d,
    remaining gap pressure).  The switch is lexicographic, with the listed
    Hamiltonian path breaking ties.

Assumption challenge.
  I considered runners, offsets, arcs, residues, boundary faces, Fourier modes,
  relation hyperplanes, empty-window regions, and proof obligations as possible
  tournament vertices.  The chosen quotient preserves the LRC predicate
  "support-6 correction stays below the cap margin" and the recurrence of the
  obstruction as runner count / ambient dimension changes.  It destroys
  witness-time geometry, so it cannot by itself construct a lonely time.
"""

from __future__ import annotations

import itertools
import math
import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402

S12_OUT = ROOT / "05-knowledge/results/lrc14_support6_residue_cusp_codex_s12.out"


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def fstr(fr: Fraction) -> str:
    return f"{fr.numerator}/{fr.denominator}"


def support_floor_sequence(max_d: int = 13) -> None:
    section("SUPPORT-FLOOR SEQUENCE")
    print(
        "THM-538 is an exact-support statement for the signed Fourier expansion "
        "over the relation lattice.  It says all genuine support <=5 relation "
        "terms vanish, and the first live layer is support 6."
    )
    print(
        "Guardrail: do not test this by plugging zero coordinates into the raw "
        "ambient K function; that includes lower-dimensional boundary/main pieces. "
        "The proof object is the exact-support coefficient after the signed "
        "inclusion-exclusion projection."
    )
    print(f"{'d':>3} {'support 0':>10} {'1':>4} {'2':>4} {'3':>4} {'4':>4} {'5':>4} {'6':>4}")
    for d in range(6, max_d + 1):
        row = [0, 0, 0, 0, 0, 0, 1]
        print(f"{d:>3} {row[0]:>10} " + " ".join(f"{x:>4}" for x in row[1:]))
    print(
        "Sequence readout: the persistent small detail is the 0,0,0,0,0,0,1 "
        "floor.  As d/k changes, lower faces remain killed and the live object "
        "is always the six-fold residue coimage."
    )


def residue_constant_decay(max_d: int = 13) -> None:
    section("RESIDUE-CONSTANT DECAY SEQUENCE")
    blunt = 64.0 * (0.697303**6)
    print(
        "C_d is defined by K(n_1,...,n_6,0,...)=C_d(n mod 7)/(n_1...n_6)."
    )
    print(f"Blunt free-product constant 64*c1^6 = {blunt:.9g}.")
    print(f"{'d':>3} {'max |C_d|':>14} {'ratio':>11} {'argmax residue':>20} {'worst 1-marginal':>20}")
    max_sequence: list[float] = []
    marginal_sequence: list[float] = []
    for d in range(6, max_d + 1):
        best_abs = -1.0
        best_r: tuple[int, ...] | None = None
        marginals: dict[tuple[int, tuple[int, ...]], complex] = defaultdict(complex)
        for residues in itertools.product(range(1, 7), repeat=6):
            c = s12.residue_coeff(residues, d)
            ac = abs(c)
            if ac > best_abs:
                best_abs = ac
                best_r = residues
            for i in range(6):
                fixed = residues[:i] + residues[i + 1 :]
                marginals[(i, fixed)] += c
        worst_marginal = max(abs(v) for v in marginals.values())
        max_sequence.append(best_abs)
        marginal_sequence.append(worst_marginal)
        print(
            f"{d:>3} {best_abs:>14.9g} {best_abs / blunt:>11.5g} "
            f"{str(best_r):>20} {worst_marginal:>20.9g}"
        )
    print("max |C_d| sequence:", ", ".join(f"{x:.9g}" for x in max_sequence))
    print("worst one-coordinate marginal sequence:", ", ".join(f"{x:.9g}" for x in marginal_sequence))
    print(
        "The max constants fall almost geometrically, while the marginals remain "
        "nonzero.  This is the exact shape of the warning: cancellation is a "
        "relation-hyperplane/coimage phenomenon, not a free coordinate marginal."
    )


def parse_float(token: str) -> float:
    if token == "inf":
        return math.inf
    match = re.match(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?", token)
    if not match:
        raise ValueError(f"could not parse numeric token {token!r}")
    return float(match.group(0))


def parse_s12_cases() -> list[dict[str, object]]:
    text = S12_OUT.read_text()
    chunks = re.split(r"\n=+\nCASE: ", text)
    cases = []
    for chunk in chunks[1:]:
        name, rest = chunk.split("\n", 1)
        name = name.strip()
        cumulative: list[dict[str, float | int]] = []
        boundary: list[dict[str, float | int]] = []
        mode = None
        for line in rest.splitlines():
            if line.startswith("Cumulative exact signed"):
                mode = "cumulative"
                continue
            if line.startswith("Boundary-face split"):
                mode = "boundary"
                continue
            parts = line.split()
            if mode == "cumulative" and len(parts) == 5 and parts[0].isdigit():
                cumulative.append(
                    {
                        "H": int(parts[0]),
                        "relations": int(parts[1]),
                        "signed": parse_float(parts[2]),
                        "absK": parse_float(parts[3]),
                        "ratio": parse_float(parts[4]),
                    }
                )
            elif mode == "boundary" and len(parts) == 5 and parts[0].isdigit():
                boundary.append(
                    {
                        "touch": int(parts[0]),
                        "relations": int(parts[1]),
                        "signed": parse_float(parts[2]),
                        "absK": parse_float(parts[3]),
                        "ratio": parse_float(parts[4]),
                    }
                )
        if cumulative:
            cases.append({"name": name, "cumulative": cumulative, "boundary": boundary})
    return cases


def cusp_ratio_sequences() -> None:
    section("BOUNDARY/CUSP ABSOLUTE-VERSUS-SIGNED SEQUENCES")
    cases = parse_s12_cases()
    print(
        "Parsed the exact S12 support-6 relation ledger.  The displayed numbers "
        "are absK/|signed| ratios; larger means the absolute envelope is mostly "
        "cusp mass cancelled by the signed kernel."
    )
    print(f"{'case':<49} {'final H':>7} {'final ratio':>13} {'one-face ratio':>16} {'one-face count':>15}")
    final_ratios = []
    one_face_ratios = []
    for case in cases:
        final = case["cumulative"][-1]  # type: ignore[index]
        boundary = case["boundary"]  # type: ignore[assignment]
        one_face = next((row for row in boundary if row["touch"] == 1), None)
        final_ratios.append(float(final["ratio"]))
        if one_face is not None:
            one_face_ratios.append(float(one_face["ratio"]))
        print(
            f"{str(case['name'])[:49]:<49} {int(final['H']):>7} "
            f"{float(final['ratio']):>13.6g} "
            f"{float(one_face['ratio']) if one_face else math.nan:>16.6g} "
            f"{int(one_face['relations']) if one_face else 0:>15}"
        )
    print("final ratio sequence:", ", ".join(f"{x:.6g}" for x in final_ratios))
    print("one-face boundary ratio sequence:", ", ".join(f"{x:.6g}" for x in one_face_ratios))
    print(
        "Integer side: the one-face boundary relation counts are the large mass.  "
        "Fractional side: the signed/absolute ratios are the small coimage error."
    )


def center_survivor_sequences() -> None:
    section("UNIVERSAL-CENTER SURVIVOR SEQUENCE")
    survivors = [math.comb(7, s) + math.comb(9, s) - math.comb(5, s) for s in range(14)]
    totals = [math.comb(13, s) for s in range(14)]
    mixed = [t - u for t, u in zip(totals, survivors)]
    print("survivor(s)=C(7,s)+C(9,s)-C(5,s), s=0..13")
    print("survivor sequence:", survivors)
    print("mixed residual sequence:", mixed)
    print("k=13-s readout (cluster size k):")
    print(f"{'k':>3} {'survivor P':>12} {'mixed P':>10} {'total P':>9}")
    for k in range(3, 14):
        s = 13 - k
        print(f"{k:>3} {survivors[s]:>12} {mixed[s]:>10} {totals[s]:>9}")
    print(
        "The fixed denominator centers solve all-odd or 3-free small parts; the "
        "mixed parity/triadic sequence is exactly where colored resonance remains."
    )


def empty_window_sequences() -> None:
    section("EMPTY-WINDOW CERTIFICATE SEQUENCES")
    rows = [
        (8, 4, Fraction(19633, 30870), Fraction(3637, 5880), Fraction(431, 24696)),
        (9, 3, Fraction(127, 245), Fraction(2025, 4004), Fraction(1769, 140140)),
        (10, 3, Fraction(10765, 24696), Fraction(36, 91), Fraction(12937, 321048)),
        (11, 2, Fraction(21673, 70560), Fraction(25, 91), Fraction(29749, 917280)),
        (12, 1, Fraction(89533, 543312), Fraction(1, 7), Fraction(11917, 543312)),
    ]
    print("AP threshold-clearing certificate sequence from HYP-2608:")
    print(f"{'k':>3} {'degree':>6} {'Phi(AP)':>16} {'threshold':>16} {'margin':>16} {'margin dec':>12}")
    for k, degree, phi, threshold, margin in rows:
        print(
            f"{k:>3} {degree:>6} {fstr(phi):>16} {fstr(threshold):>16} "
            f"{fstr(margin):>16} {float(margin):>12.8f}"
        )
    print("degree sequence:", [r[1] for r in rows])
    print("margin numerator sequence:", [r[4].numerator for r in rows])
    print(
        "The degree drop 4,3,3,2,1 is the region-side analogue of the support "
        "floor: as k grows, the required proof obligation moves to fewer moments."
    )


def hamiltonian_paths_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [Counter() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v, count in list(dp[mask].items()):
            for w in range(n):
                if not (mask & (1 << w)) and adj[v][w]:
                    dp[mask | (1 << w)][w] += count
    return sum(dp[-1].values())


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in itertools.combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w, ok in enumerate(adj[v]):
            if ok and w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)
    seen.clear()
    sizes = []
    for start in reversed(order):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for w, ok in enumerate(radj[v]):
                if ok and w not in seen:
                    seen.add(w)
                    stack.append(w)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS: SEQUENCE-FAMILY QUOTIENTS")
    names = [
        "support_floor",
        "residue_constant_decay",
        "cusp_signed_ratio",
        "empty_window_degree_drop",
        "universal_center_survivors",
        "low_height_wall_ledger",
        "coimage_reciprocal_tail",
    ]
    metrics = {
        # exactness, leverage, portability, negative remaining-gap pressure
        "support_floor": (4, 5, 5, -1),
        "residue_constant_decay": (3, 4, 5, -2),
        "cusp_signed_ratio": (3, 5, 4, -2),
        "empty_window_degree_drop": (4, 3, 4, -1),
        "universal_center_survivors": (4, 3, 5, -1),
        "low_height_wall_ledger": (2, 5, 3, -4),
        "coimage_reciprocal_tail": (1, 5, 5, -5),
    }
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    flips = 0
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            adj[i][j] = metrics[a] > metrics[b] or (metrics[a] == metrics[b] and i < j)
            if i < j and not adj[i][j]:
                flips += 1
    scores = {names[i]: sum(adj[i]) for i in range(n)}
    print("Hamiltonian proof path:", sorted(names, key=lambda x: scores[x], reverse=True))
    print("score histogram:", dict(sorted(Counter(scores.values()).items())))
    print("directed 3-cycles:", count_directed_3cycles(adj))
    print("SCC sizes:", scc_sizes(adj))
    print("Hamiltonian path count:", hamiltonian_paths_count(adj))
    print("edge flips against listed order:", flips)
    print(
        "Assumption challenged: sequence-family vertices preserve the analytic "
        "tail predicate and the recursion across d/k, while runner vertices hide "
        "the support floor and the coimage cancellation."
    )


def final_reading() -> None:
    section("READING FOR THE LRC(14) PROOF ROUTE")
    print(
        "1. The support-floor sequence is the hard cancellation: all lower support "
        "faces vanish, so boundary/cusp mass must be read through six-fold residues."
    )
    print(
        "2. The residue constants decay with ambient dimension, but coordinate "
        "marginals stay nonzero; the proof object is the signed coimage over "
        "sum e_i n_i=0, not a product of one-coordinate cancellations."
    )
    print(
        "3. The S12 cusp ratio sequence shows why absolute Minkowski volume is the "
        "wrong measure.  The integer relation counts are huge exactly where the "
        "signed reciprocal mass is tiny."
    )
    print(
        "4. The center-survivor and empty-window sequences are the companion finite "
        "spines.  They say which runner-count regimes need a high-degree local "
        "certificate, and which small-speed residues pass to colored recurrence."
    )
    print(
        "5. LRC(14) remains open.  The sharpened sub-problem is a tail theorem for "
        "residue-addressed reciprocal sums after low-height wall deletion."
    )


def main() -> None:
    section("LRC(14) SIGNED-MASS SEQUENCE SPINE")
    print(
        "Goal: collect the integer/fractional sequences behind the observation "
        "'large absolute mass but tiny signed mass', especially on boundary/cusp faces."
    )
    support_floor_sequence()
    residue_constant_decay()
    cusp_ratio_sequences()
    center_survivor_sequences()
    empty_window_sequences()
    tournament_analysis()
    final_reading()


if __name__ == "__main__":
    main()
