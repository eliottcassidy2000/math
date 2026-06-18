#!/usr/bin/env python3
"""
lrc14_colored_resonance_deficit_codex.py

Codex 2026-06-18: proof-oriented scout for the LRC(14) colored
arithmetic discrepancy bound.

HYP-2594 proved the crude interval lower bound

    actual_count(q=14V) >= V*Sigma(P,E) - K(P,E).

This script pushes toward the sharper arithmetic theorem.  It records two
objects:

  1. A Fourier/color-resonance identity explaining why most micro-components
     should cancel after summing the 14 residue colors.
  2. Exact deficit decompositions by the coarse safe components of G_P, testing
     additive bounds such as C*(k+c_GP)+1 rather than raw K.

Tournament Analysis declaration.
  Vertex set: candidate discrepancy bounds.
  Pairwise observable: (violation count, worst normalized pressure, median
    bound size) on the reproducible stress bank.
  Gauge/tie path: fewer violations wins; then lower worst pressure; then
    smaller median bound.  Ties follow the printed candidate order.
  Fingerprints: score histogram, directed 3-cycles, SCCs, Hamiltonian paths.

Assumption challenge.
  I considered phase colors, micro-interval endpoints, coarse G_P components,
  Fourier resonances, and proof-obligation bounds as vertices.  The micro-
  endpoint quotient preserves the trivial K-bound but destroys cancellation.
  The coarse-component/resonance quotient preserves the actual q=14V counting
  predicate and the color sum, while forgetting individual endpoint origins.
  The challenged assumption is that every interval endpoint costs one unit of
  discrepancy; the color-resonance identity says only compatible V-resonances
  should survive.
"""

from __future__ import annotations

import itertools
import random
import statistics
import sys
from collections import Counter
from fractions import Fraction as F

import lrc14_colored_discrepancy_bound_codex as disc
import lrc14_global_threshold_ladder_codex as ladder
import lrc14_phase_color_reservoir_codex as pc

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

RNG = random.Random(20260618 + 2595)


def ffloat(x: F) -> str:
    return f"{float(x):.6f}"


def first_covering_lift(P: tuple[int, ...], E: tuple[int, ...], extra: int = 180) -> int | None:
    lo = max(E) + 14
    for V in range(lo, lo + extra + 1):
        S = pc.build_S(P, E, V)
        if len(S) == 13 and min(V - e for e in E) > 13 and pc.gcd_all(S) == 1 and pc.is_covering(S):
            return V
    return None


def closed_actual(P: tuple[int, ...], E: tuple[int, ...], V: int) -> int:
    return pc.actual_crt_count(pc.build_S(P, E, V), V)


def row_metrics(label: str, P: tuple[int, ...], E: tuple[int, ...], V: int) -> dict:
    sig, K, _ = disc.sigma_K(P, E)
    open_count = disc.exact_open_count(P, E, V)
    actual = closed_actual(P, E, V)
    gp = ladder.safe_set(P)
    cgp = max(1, len(gp))
    k = len(E)
    color_errors = []
    interval_errors = []
    max_cell = (F(-10**9), None)
    positive_cell_sum = F(0)

    for b in range(14):
        err = F(0)
        comps = 0
        for I in gp:
            sub = ladder.intersect([I], pc.phase_safe_set(E, b))
            err += V * ladder.measure(sub) - disc.open_grid_count(sub, b, V)
            comps += len(sub)
        color_errors.append((err, b, comps))

    for idx, I in enumerate(gp):
        err = F(0)
        for b in range(14):
            sub = ladder.intersect([I], pc.phase_safe_set(E, b))
            meas = ladder.measure(sub)
            count = disc.open_grid_count(sub, b, V)
            cell_err = V * meas - count
            if cell_err > max_cell[0]:
                max_cell = (cell_err, (idx, b, len(sub), meas, count))
            if cell_err > 0:
                positive_cell_sum += cell_err
            err += cell_err
        interval_errors.append((err, idx))

    return {
        "label": label,
        "P": P,
        "E": E,
        "V": V,
        "k": k,
        "cgp": cgp,
        "Sigma": sig,
        "K": K,
        "open": open_count,
        "actual": actual,
        "open_deficit": sig * V - open_count,
        "actual_deficit": sig * V - actual,
        "boundary_bonus": actual - open_count,
        "max_color_error": max(color_errors),
        "min_color_error": min(color_errors),
        "max_interval_error": max(interval_errors),
        "max_cell_error": max_cell,
        "positive_cell_sum": positive_cell_sum,
    }


def deterministic_rows() -> list[tuple[str, tuple[int, ...], tuple[int, ...], int]]:
    rows = [
        ("HYP2594_stored_max_deficit", (4, 6, 10), (0, 49, 85, 97, 196, 205, 207, 269, 284, 326), 349),
        ("P_empty_large_open_deficit", (), (0, 40, 83, 196, 230, 390, 455, 505, 553, 587, 590, 635, 660), 686),
        ("P_empty_boundary_bonus", (), (0, 93, 163, 185, 218, 272, 283, 322, 329, 473, 486, 542, 619), 658),
        ("constant_100_failure", (1, 2, 11), (0, 84, 293, 301, 355, 416, 485, 665, 843, 886), 1203),
        ("structured_worst_K_over_Sigma", (1, 2, 3, 7, 8, 9, 10, 11, 12, 13), (0, 1, 6), 244),
        ("quarter_min", (1, 2, 3), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10), 140),
        ("via_zero_k7", (1, 2, 3, 6, 12, 13), (0, 2, 3, 4, 5, 6, 8), 98),
        ("hard_core_84m_m1", tuple(list(range(1, 12)) + [13]), (0,), 84),
        ("hard_core_84m_m5", tuple(list(range(1, 12)) + [13]), (0,), 420),
    ]
    return rows


def stress_rows(limit: int = 140) -> list[tuple[str, tuple[int, ...], tuple[int, ...], int]]:
    rows = []
    seen = set()
    attempts = 0
    while len(rows) < limit and attempts < 20 * limit:
        attempts += 1
        k = RNG.choice([3, 4, 5, 6, 8, 10, 11, 12, 13])
        psz = 13 - k
        P = tuple(sorted(RNG.sample(range(1, 14), psz)))
        spread = RNG.choice([3 * k, 5 * k, 8 * k, 13 * k, 21 * k, 34 * k])
        E = (0, *tuple(sorted(RNG.sample(range(1, spread + 1), k - 1))))
        if pc.gcd_all(P + E) != 1:
            continue
        key = (P, E)
        if key in seen:
            continue
        seen.add(key)
        V = first_covering_lift(P, E, 220)
        if V is None:
            continue
        rows.append((f"stress_{len(rows):03d}", P, E, V))
    return rows


def candidate_bounds(row: dict) -> dict[str, F]:
    k = row["k"]
    c = row["cgp"]
    return {
        "100": F(100, 1),
        "7(k+c)+1": F(7 * (k + c) + 1, 1),
        "8(k+c)+1": F(8 * (k + c) + 1, 1),
        "8kc+1": F(8 * k * c + 1, 1),
        "K": F(row["K"], 1),
    }


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in itertools.combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = set()
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


def hamiltonian_paths_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [Counter() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v, cnt in list(dp[mask].items()):
            if not cnt:
                continue
            for w in range(n):
                if not (mask & (1 << w)) and adj[v][w]:
                    dp[mask | (1 << w)][w] += cnt
    return sum(dp[(1 << n) - 1].values())


def bound_tournament(rows: list[dict]) -> dict:
    names = list(candidate_bounds(rows[0]).keys())
    stats = {}
    for name in names:
        violations = 0
        pressures = []
        sizes = []
        for row in rows:
            deficit = max(row["actual_deficit"], F(0))
            bound = candidate_bounds(row)[name]
            sizes.append(float(bound))
            pressures.append(float(deficit / bound) if bound else 0.0)
            if deficit > bound:
                violations += 1
        stats[name] = {
            "violations": violations,
            "max_pressure": max(pressures),
            "median_bound": statistics.median(sizes),
        }
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            key_a = (stats[a]["violations"], stats[a]["max_pressure"], stats[a]["median_bound"], i)
            key_b = (stats[b]["violations"], stats[b]["max_pressure"], stats[b]["median_bound"], j)
            if key_a < key_b:
                adj[i][j] = True
    scores = [sum(row) for row in adj]
    return {
        "names": names,
        "stats": stats,
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles3": count_directed_3cycles(adj),
        "scc": scc_sizes(adj),
        "hamiltonian_paths": hamiltonian_paths_count(adj),
        "leaders": sorted(zip(scores, names), reverse=True),
    }


def print_fourier_identity() -> None:
    print("\n" + "=" * 88)
    print("STEP 1. Rigorous Fourier/color-resonance identity")
    print("=" * 88)
    print("Let h(x)=1_{||x||>=1/14}, with boundary values taken as 1/2")
    print("for the Fourier identity.  For a fixed color b define")
    print("  F_b(x)=prod_{p in P} h(px) * prod_{e in E} h(ex-b/14).")
    print("The half-boundary colored count is")
    print("  N(P,E,V)=sum_b sum_{t=0}^{V-1} F_b(t/V+b/(14V)).")
    print("Writing Fourier coefficients in x, the grid quadrature gives")
    print("  N - V*Sigma = V * sum_{ell != 0} sum_b Fhat_b(ell V) exp(2pi i ell b/14).")
    print("Expanding each h-factor, a term can survive the color sum only if")
    print("  sum_p p*n_p + sum_e e*n_e = ell*V")
    print("  and  sum_e n_e == ell (mod 14).")
    print("This is the colored resonance condition.  Raw interval endpoints count")
    print("all micro-events; the actual discrepancy only sees these compatible")
    print("V-resonances.  The closed CRT count differs from the half-boundary")
    print("count by the explicit binding-boundary bonus, which the diagnostics print.")


def print_row(row: dict) -> None:
    bounds = candidate_bounds(row)
    deficit = row["actual_deficit"]
    print(f"  {row['label']}: P={row['P']} Emax={max(row['E'])} k={row['k']} V={row['V']} cGP={row['cgp']}")
    print(
        f"    Sigma={row['Sigma']} ({ffloat(row['Sigma'])}) K={row['K']} "
        f"actual={row['actual']} open={row['open']} boundary_bonus={row['boundary_bonus']}"
    )
    print(
        f"    actual_deficit={deficit} ({float(deficit):.3f}) "
        f"open_deficit={float(row['open_deficit']):.3f} "
        f"positive_cell_sum={float(row['positive_cell_sum']):.3f}"
    )
    print(
        f"    max_color_error={float(row['max_color_error'][0]):.3f} at b={row['max_color_error'][1]} "
        f"max_interval_error={float(row['max_interval_error'][0]):.3f} at comp={row['max_interval_error'][1]}"
    )
    print(
        "    bounds "
        + " ".join(
            f"{name}:{float(bound):.1f}{'!' if deficit > bound else ''}"
            for name, bound in bounds.items()
            if name != "K"
        )
    )


def main() -> None:
    print("=" * 88)
    print("LRC(14) colored resonance deficit scout")
    print("=" * 88)
    print("Goal: replace raw K by a color-resonance/coarse-component discrepancy bound.")
    print_fourier_identity()

    print("\n" + "=" * 88)
    print("STEP 2. Named and adversarial rows: coarse decomposition")
    print("=" * 88)
    named = [row_metrics(*args) for args in deterministic_rows()]
    for row in named:
        print_row(row)

    print("\n" + "=" * 88)
    print("STEP 3. Reproducible stress bank for candidate bounds")
    print("=" * 88)
    stress = [row_metrics(*args) for args in stress_rows()]
    rows = named + stress
    print(f"  stress rows={len(stress)} total rows={len(rows)}")
    print(f"  zero actual witnesses={sum(1 for r in rows if r['actual'] == 0)}")
    for name in candidate_bounds(rows[0]):
        violations = []
        worst = None
        for row in rows:
            deficit = max(row["actual_deficit"], F(0))
            bound = candidate_bounds(row)[name]
            pressure = deficit / bound if bound else F(0)
            if worst is None or pressure > worst[0]:
                worst = (pressure, row, bound)
            if deficit > bound:
                violations.append(row)
        print(
            f"  bound {name:10s}: violations={len(violations):3d} "
            f"max_pressure={float(worst[0]):.6f} "
            f"worst={worst[1]['label']} deficit={float(worst[1]['actual_deficit']):.3f} "
            f"bound={float(worst[2]):.1f} k={worst[1]['k']} cGP={worst[1]['cgp']} K={worst[1]['K']}"
        )

    worst_def = max(rows, key=lambda r: r["actual_deficit"])
    worst_pressure = max(rows, key=lambda r: max(r["actual_deficit"], F(0)) / candidate_bounds(r)["8(k+c)+1"])
    print("\n  Worst actual deficit row:")
    print_row(worst_def)
    print("\n  Worst pressure against 8(k+c)+1:")
    print_row(worst_pressure)

    tour = bound_tournament(rows)
    print("\n  Candidate-bound tournament:")
    for name in tour["names"]:
        st = tour["stats"][name]
        print(
            f"    {name:10s}: violations={st['violations']} "
            f"max_pressure={st['max_pressure']:.6f} median_bound={st['median_bound']:.1f}"
        )
    print(
        f"    score_hist={tour['score_hist']} cycles3={tour['cycles3']} "
        f"scc={tour['scc']} hamiltonian_paths={tour['hamiltonian_paths']}"
    )
    print(f"    leaders={tour['leaders']}")

    print("\n" + "=" * 88)
    print("TAKEAWAY")
    print("=" * 88)
    print("  1. The exact discrepancy is governed by a colored resonance condition,")
    print("     not by raw micro-component count K.")
    print("  2. Direct rows refute too-small constant bounds, but additive bounds")
    print("     C*(k+cGP)+1 are the right-looking next theorem target.")
    print("  3. A proof should bound the surviving Fourier/color resonances or,")
    print("     equivalently, prove endpoint cancellation inside each coarse G_P")
    print("     component after summing the 14 phase colors.")
    print("  4. This is still not a proof of LRC(14), but it is a precise lemma:")
    print("     colored discrepancy <= C*(k+cGP)+1 plus the Sigma floor would leave")
    print("     only a finite low-V exact check.")


if __name__ == "__main__":
    main()
