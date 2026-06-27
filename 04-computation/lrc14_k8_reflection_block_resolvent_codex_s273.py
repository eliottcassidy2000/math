#!/usr/bin/env python3
"""Exact k=8 reflection-block scout for the LRC14 bounded-core node.

This is an executable research note, not a proof.  It connects:

* HYP-3085: pairwise sector co-emptiness matrix M and the Z/2 reflection block;
* HYP-3122/HYP-3132: the phi4 / biquadratic resolvent view of the k=8 dip;
* HYP-3133: A000568 as a finite middle quotient for edge witnesses.

Tournament Analysis declaration:
  vertices: proof carriers / quotient coordinates, not runners;
  pairwise observable: majority vote over retained LRC payload axes;
  switch: orient A -> B when A wins more axes, with a Hamiltonian tie path;
  tie path: exact core block, biquadratic fold, SPEC certificate, A000568
            middle shadow, Worpitzky coefficient vocabulary, raw spectrum.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Dict, Iterable, List, Sequence, Tuple

import sympy as sp


def sector_of(x: Fraction) -> int:
    return int((x % 1) * 7)


def pairwise_matrix(speeds: Iterable[int]) -> List[List[Fraction]]:
    """M[i][j] = measure of times sectors i and j are both empty."""
    speeds = sorted(set(speeds))
    cuts = {Fraction(0), Fraction(1)}
    for e in speeds:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            cuts.add(Fraction(m, 7 * e))
    cuts = sorted(cuts)
    matrix = [[Fraction(0) for _ in range(7)] for _ in range(7)]
    for left, right in zip(cuts, cuts[1:]):
        if right <= left:
            continue
        mid = (left + right) / 2
        weight = right - left
        covered = {sector_of(e * mid) for e in speeds}
        empty = [s for s in range(7) if s not in covered]
        for i in empty:
            for j in empty:
                matrix[i][j] += weight
    return matrix


def fmt(frac: Fraction) -> str:
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def sp_matrix(rows: Sequence[Sequence[Fraction]]) -> sp.Matrix:
    return sp.Matrix(
        [[sp.Rational(x.numerator, x.denominator) for x in row] for row in rows]
    )


def charpoly_expr(rows: Sequence[Sequence[Fraction]]) -> str:
    x = sp.Symbol("x")
    expr = sp.factor(sp_matrix(rows).charpoly(x).as_expr())
    return str(expr)


def eigenvals_float(rows: Sequence[Sequence[Fraction]]) -> List[float]:
    vals = sp_matrix(rows).evalf().eigenvals()
    out: List[float] = []
    for val, mult in vals.items():
        out.extend([float(sp.re(val))] * mult)
    return sorted(out)


def reflection_core_blocks(matrix: List[List[Fraction]]) -> Dict[str, object]:
    """Fold sectors 1..5 by s -> 6-s, leaving sector 6 as boundary leakage.

    The invariant core has orbits A={1,5}, B={2,4}, C={3}.  For a reflection-
    symmetric vector with values (a,b,c), M acts through a rational 3x3 block.
    The antisymmetric block is on A- and B-differences.
    """
    # Symmetric quotient: row representatives 1,2,3; column orbits {1,5},{2,4},{3}.
    sym_rows = []
    for i in (1, 2, 3):
        sym_rows.append(
            [
                matrix[i][1] + matrix[i][5],
                matrix[i][2] + matrix[i][4],
                matrix[i][3],
            ]
        )

    # Antisymmetric quotient: representatives 1,2; columns e1-e5, e2-e4.
    anti_rows = []
    for i in (1, 2):
        anti_rows.append(
            [
                matrix[i][1] - matrix[i][5],
                matrix[i][2] - matrix[i][4],
            ]
        )

    core = [[matrix[i][j] for j in range(1, 6)] for i in range(1, 6)]
    shell_2x2 = [row[:2] for row in sym_rows[:2]]
    center_coupling = [sym_rows[0][2], sym_rows[1][2], sym_rows[2][0], sym_rows[2][1]]
    boundary_vector = [matrix[i][6] for i in range(1, 6)]
    core_s2 = sum(matrix[i][j] for i in range(1, 6) for j in range(i + 1, 6))
    boundary_s2 = sum(matrix[i][6] for i in range(1, 6))
    total_s2 = sum(matrix[i][j] for i in range(7) for j in range(i + 1, 7))
    return {
        "core": core,
        "sym": sym_rows,
        "shell_2x2": shell_2x2,
        "center_coupling": center_coupling,
        "anti": anti_rows,
        "boundary_vector": boundary_vector,
        "core_s2": core_s2,
        "boundary_s2": boundary_s2,
        "total_s2": total_s2,
    }


def reflection_ok(matrix: List[List[Fraction]]) -> bool:
    rho = {1: 5, 5: 1, 2: 4, 4: 2, 3: 3}
    for i in range(1, 6):
        for j in range(1, 6):
            if matrix[i][j] != matrix[rho[i]][rho[j]]:
                return False
    return True


def fold_polynomials() -> List[Dict[str, object]]:
    """All reflected 4-root shell choices around center 3 in sectors 0..6."""
    rows = []
    t = sp.Symbol("t")
    y = sp.Symbol("y")
    for shells in combinations((1, 2, 3), 2):
        y_roots = [d * d for d in shells]
        y_poly = sp.expand((y - y_roots[0]) * (y - y_roots[1]))
        t_roots = sorted([3 - shells[0], 3 + shells[0], 3 - shells[1], 3 + shells[1]])
        t_poly = sp.expand(sp.prod(t - r for r in t_roots))
        rows.append(
            {
                "shells": shells,
                "roots_t": tuple(t_roots),
                "roots_y": tuple(y_roots),
                "y_poly": str(y_poly),
                "t_poly": str(t_poly),
                "y_sum": sum(y_roots),
                "y_product": y_roots[0] * y_roots[1],
                "discriminant": (y_roots[1] - y_roots[0]) ** 2,
                "uses_boundary_shell": 3 in shells,
                "is_hyp3132_inner_fold": shells == (1, 2),
            }
        )
    return rows


AXES = (
    "preserves_status",
    "preserves_route",
    "exact_arithmetic",
    "folds_k8_obligation",
    "connects_spec_floor",
    "connects_a000568",
    "formalization_ready",
    "scalar_guardrail",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: Dict[str, int]


def carrier(name: str, values: Sequence[int]) -> Carrier:
    return Carrier(name, dict(zip(AXES, values)))


CARRIERS = [
    carrier("exact_reflection_3x3_core_block", (9, 7, 10, 10, 7, 6, 9, 9)),
    carrier("hyp3132_inner_biquadratic_fold", (8, 7, 10, 10, 6, 7, 8, 9)),
    carrier("hyp3129_spec_resonance_certificate", (10, 8, 9, 6, 10, 5, 8, 9)),
    carrier("hyp3133_a000568_middle_shadow", (6, 6, 7, 5, 6, 10, 7, 9)),
    carrier("worpitzky_path_moment_vocabulary", (4, 4, 8, 6, 4, 7, 6, 8)),
    carrier("antisymmetric_2x2_nonmax_block", (5, 3, 9, 7, 3, 5, 8, 8)),
    carrier("boundary_sector6_leak_vector", (7, 5, 10, 7, 6, 4, 8, 8)),
    carrier("raw_numeric_spectrum", (2, 1, 3, 3, 2, 2, 3, 1)),
]

TIE_PATH = [c.name for c in CARRIERS]


def compare(a: Carrier, b: Carrier) -> str:
    av = bv = 0
    for axis in AXES:
        if a.scores[axis] > b.scores[axis]:
            av += 1
        elif b.scores[axis] > a.scores[axis]:
            bv += 1
    if av > bv:
        return a.name
    if bv > av:
        return b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_summary(vertices: Sequence[Carrier]) -> Dict[str, object]:
    adj = {v.name: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        winner = compare(a, b)
        loser = b.name if winner == a.name else a.name
        adj[winner].add(loser)
    score_hist: Dict[int, int] = {}
    for outs in adj.values():
        score_hist[len(outs)] = score_hist.get(len(outs), 0) + 1
    cycles = 0
    for a, b, c in combinations(sorted(adj), 3):
        if (b in adj[a] and c in adj[b] and a in adj[c]) or (
            a in adj[b] and b in adj[c] and c in adj[a]
        ):
            cycles += 1
    path = sorted(vertices, key=lambda v: (-len(adj[v.name]), TIE_PATH.index(v.name)))
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": cycles,
        "selected_path": [v.name for v in path],
    }


def print_matrix(name: str, rows: Sequence[Sequence[Fraction]]) -> None:
    print(f"{name}:")
    for row in rows:
        print("  [" + ", ".join(f"{fmt(x):>12}" for x in row) + "]")


def main() -> None:
    print("HYP-3139 / codex-2026-06-27-S273")
    print("Exact k=8 reflection-block resolvent scout; executable evidence, not a proof.")
    print()

    print("1. REFLECTED SHELL POLYNOMIALS AROUND CENTER 3")
    for row in fold_polynomials():
        marker = " <-- HYP-3132 inner fold" if row["is_hyp3132_inner_fold"] else ""
        print(
            f"shells={row['shells']} roots_t={row['roots_t']} "
            f"y_poly={row['y_poly']} discr={row['discriminant']} "
            f"boundary={row['uses_boundary_shell']}{marker}"
        )
    print(
        "Readout: the HYP-3132 quartic is the unique reflected 4-root fold using "
        "only the two inner shells.  Boundary-shell folds remain solvable, but "
        "they are not the k=8 bounded-core node."
    )
    print()

    print("2. EXACT PAIRWISE CO-EMPTINESS REFLECTION BLOCKS")
    for k in (8, 9, 10):
        matrix = pairwise_matrix(range(k))
        blocks = reflection_core_blocks(matrix)
        print(f"\n--- consec_{k} ---")
        print(f"core_reflection_exact={reflection_ok(matrix)}")
        print(f"S2_core_1..5={fmt(blocks['core_s2'])}")
        print(f"S2_boundary_with_sector6={fmt(blocks['boundary_s2'])}")
        print(f"S2_total={fmt(blocks['total_s2'])}")
        print(
            "core_fraction_of_S2="
            f"{float(blocks['core_s2'] / blocks['total_s2']):.6f}"
        )
        print_matrix("reflection_symmetric_3x3_block_on_(1+5,2+4,3)", blocks["sym"])
        print(f"sym_charpoly={charpoly_expr(blocks['sym'])}")
        print(f"sym_eigenvalues_float={[round(v, 9) for v in eigenvals_float(blocks['sym'])]}")
        print_matrix("inner_shell_2x2_page_on_(1+5,2+4)", blocks["shell_2x2"])
        print(f"shell_charpoly={charpoly_expr(blocks['shell_2x2'])}")
        print(f"shell_eigenvalues_float={[round(v, 9) for v in eigenvals_float(blocks['shell_2x2'])]}")
        print(
            "center_shell_coupling=["
            + ", ".join(fmt(x) for x in blocks["center_coupling"])
            + "]"
        )
        print_matrix("antisymmetric_2x2_block_on_(1-5,2-4)", blocks["anti"])
        print(f"anti_charpoly={charpoly_expr(blocks['anti'])}")
        print(f"anti_eigenvalues_float={[round(v, 9) for v in eigenvals_float(blocks['anti'])]}")
        print(
            "boundary_sector6_vector=["
            + ", ".join(fmt(x) for x in blocks["boundary_vector"])
            + "]"
        )

    print()
    print("3. PROOF-ROUTE INTERPRETATION")
    print(
        "The k=8 matrix has an exact reflection core on sectors 1..5.  The "
        "two reflected shell pairs (1,5) and (2,4) give a genuine 2x2 page, "
        "matching the degree-2 biquadratic reduction in HYP-3132; the fixed "
        "center sector 3 is an extra scalar coupling.  The all-positive "
        "S2/Perron witness lives in the full reflection-symmetric 3x3 block, "
        "while the antisymmetric 2x2 block is a nonmaximizing oscillation.  "
        "Sector 6 is a boundary leakage vector, not part of the HYP-3132 "
        "inner-shell biquadratic.  Therefore the next proof target can be made "
        "finite: bound the inner 2x2 shell page, control the center coupling "
        "and boundary vector, then attach the phi4 sign for the S3/S4 k=8 "
        "correction."
    )
    print(
        "Worpitzky/Eulerian path moments remain useful as coefficient-expansion "
        "vocabulary, but this scout says the load-bearing certificate is smaller: "
        "a reflected shell 2x2 page, a center-coupling scalar page, and a "
        "one-sector boundary page."
    )
    print()

    print("4. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    summary = tournament_summary(CARRIERS)
    print(f"vertices={','.join(c.name for c in CARRIERS)}")
    print(f"axes={','.join(AXES)}")
    print(f"score_hist={summary['score_hist']}")
    print(f"directed_3cycles={summary['directed_3cycles']}")
    print("selected_path=" + " -> ".join(summary["selected_path"]))
    print()

    print("5. CANDIDATE LEMMA")
    print(
        "In the bounded-core k=8 row, after removing the pinned sector 0, the "
        "sectors 1..5 split under s->6-s into an inner reflected 2x2 shell "
        "page, a fixed-center coupling, an antisymmetric 2x2 oscillation, and "
        "a boundary-sector-6 leakage vector.  The HYP-3132 biquadratic roots "
        "are exactly the two inner reflected shells.  A proof of the 2x2 shell "
        "bound, plus center/boundary ceilings and the known phi4 sign, would "
        "close the k=8 dip without invoking any infinite tournament quotient "
        "or general quartic machinery."
    )


if __name__ == "__main__":
    main()
