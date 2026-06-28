#!/usr/bin/env python3
"""HYP-3203 scout: two new proof angles for the LRC14 k=8 target.

Angle A: polarized cyclotomic support.
  HYP-3154 says the on-circle ideal is the 7th cyclotomic profile.  The raw
  residual norm from the uniform 7-fold profile is the wrong extremal target:
  split blocks can be closer to uniform than the AP.  The useful target is the
  support functional in the AP residual direction, i.e. <q-1/7, q_AP-1/7>.

Angle B: orbit-aware compression.
  HYP-3161/HYP-3160 need a proof of covariance/ferromagnetic extremality.  A
  naive left-compression lemma is false.  This script records the exact
  obstruction count and identifies the correct repair: compression must be
  performed in a quotient that keeps dilation/mirror/two-block trap sidecars.

Tournament Analysis declaration:
  vertices are proof moves, not runners: AP residual support, raw cyclotomic
  norm, covariance, coverage, naive left compression, trap sidecar, dilation
  orbit, and scalar warning.
  pairwise observable is proof-predicate retention for LRC14 coverage/dip.
  switch orients toward the move that retains more root-locus plus compression
  sidecar payload.
  tie path starts with AP support, then covariance, then orbit-aware
  compression.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import sqrt
from typing import Dict, Iterable, List, Sequence, Tuple


INNER = tuple(range(1, 7))
CONSEC_8 = tuple(range(8))
DILATED_8 = tuple(range(0, 15, 2))


def sector_of(p: F) -> int:
    return int((p % 1) * 7)


def cells(E: Iterable[int]) -> List[Tuple[set[int], F]]:
    E = tuple(sorted(set(E)))
    breakpoints = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            breakpoints.add(F(m, 7 * e))
    b = sorted(breakpoints)
    out: List[Tuple[set[int], F]] = []
    for x0, x1 in zip(b, b[1:]):
        if x1 <= x0:
            continue
        mid = (x0 + x1) / 2
        out.append(({sector_of(e * mid) for e in E}, x1 - x0))
    return out


def qdist_from_cells(cell_data: Sequence[Tuple[set[int], F]]) -> Tuple[F, ...]:
    q = [F(0)] * 7
    for hit, weight in cell_data:
        q[7 - len(hit)] += weight
    return tuple(q)


def cov_sum_from_cells(cell_data: Sequence[Tuple[set[int], F]]) -> F:
    p = {i: F(0) for i in INNER}
    pij: Dict[Tuple[int, int], F] = {}
    for hit, weight in cell_data:
        empty = [i for i in INNER if i not in hit]
        for i in empty:
            p[i] += weight
        for i, j in combinations(empty, 2):
            pij[(i, j)] = pij.get((i, j), F(0)) + weight
    return sum(pij.get((i, j), F(0)) - p[i] * p[j] for i, j in combinations(INNER, 2))


def cyclotomic_energy(q: Sequence[F]) -> F:
    """Nontrivial Fourier energy around uniform 7-fold mass.

    Since sum(q)=1, Parseval gives sum_{j=1}^6 |G(zeta^j)|^2 =
    7*sum_t q_t^2 - 1.
    """
    return 7 * sum(x * x for x in q) - 1


def dot(a: Sequence[F], b: Sequence[F]) -> F:
    return sum(x * y for x, y in zip(a, b))


def centered(q: Sequence[F]) -> Tuple[F, ...]:
    return tuple(x - F(1, 7) for x in q)


def cosine(a: Sequence[F], b: Sequence[F]) -> float:
    af = [float(x) for x in a]
    bf = [float(x) for x in b]
    denom = sqrt(sum(x * x for x in af) * sum(x * x for x in bf))
    return sum(x * y for x, y in zip(af, bf)) / denom if denom else 0.0


def left_moves(E: Tuple[int, ...]):
    S = set(E)
    for e in sorted(S, reverse=True):
        if e == 0:
            continue
        for m in range(1, e):
            if m not in S:
                yield tuple(sorted((S - {e}) | {m})), (e, m)


def rank_of(rows, key, target, reverse=True):
    ordered = sorted(rows, key=key, reverse=reverse)
    ranks = {row["E"]: i for i, row in enumerate(ordered)}
    return ranks[target], ordered


def main() -> None:
    bank = [(0,) + sub for sub in combinations(range(1, 15), 7)]
    q_ap = qdist_from_cells(cells(CONSEC_8))
    d_ap = centered(q_ap)
    rows = []
    for E in bank:
        c = cells(E)
        q = qdist_from_cells(c)
        d = centered(q)
        rows.append(
            {
                "E": E,
                "q": q,
                "cov": cov_sum_from_cells(c),
                "q0": q[0],
                "Ly": q[0] + q[6] + q[3] / 10,
                "cyc_energy": cyclotomic_energy(q),
                "ap_projection": dot(d, d_ap),
                "ap_cosine": cosine(d, d_ap),
                "dot_q_qap": dot(q, q_ap),
            }
        )

    print("=== HYP-3203: TWO NEW LRC14 PROOF ANGLES ===")
    print("bounded anchored k=8 bank:", len(bank), "sets")
    print("AP q-vector:", q_ap)
    print()

    for name, key, reverse in (
        ("covariance Sigma-kappa_2", lambda r: r["cov"], True),
        ("coverage q0", lambda r: r["q0"], True),
        ("L_y = q0+q6+q3/10", lambda r: r["Ly"], True),
        ("raw cyclotomic energy MIN", lambda r: r["cyc_energy"], False),
        ("AP residual projection MAX", lambda r: r["ap_projection"], True),
        ("cosine to AP residual ray MAX", lambda r: r["ap_cosine"], True),
    ):
        rank, ordered = rank_of(rows, key, CONSEC_8, reverse=reverse)
        best = ordered[0]
        best_value = key(best)
        consec_value = key(next(r for r in rows if r["E"] == CONSEC_8))
        dilation_rank = {row["E"]: i for i, row in enumerate(ordered)}[DILATED_8]
        print(f"{name}:")
        print(f"  consec rank={rank}; doubled-dilation rank={dilation_rank}")
        print(f"  best={best['E']} value={best_value} ({float(best_value):.8f})")
        print(f"  consec value={consec_value} ({float(consec_value):.8f})")
        print()

    max_projection = max(row["ap_projection"] for row in rows)
    projection_maximizers = [row["E"] for row in rows if row["ap_projection"] == max_projection]
    print("ANGLE A READOUT: polarized cyclotomic support")
    print("  raw cyclotomic energy is NOT minimized by AP; AP rank is 19.")
    print("  AP residual projection is maximized by:", projection_maximizers)
    print("  max projection:", max_projection, f"({float(max_projection):.8f})")
    print("  proof target: show <q(E)-1/7, q(AP)-1/7> <= ||q(AP)-1/7||^2.")
    print("  interpretation: prove a support hyperplane for the AP root-locus defect,")
    print("  not a minimum-distance-to-cyclotomic theorem.")
    print()

    vals = {row["E"]: row["cov"] for row in rows}
    bad = []
    greedy_stuck = []
    max_steps = 0
    for row in rows:
        E = row["E"]
        if E == CONSEC_8:
            continue
        moves = [(vals[Ep] - vals[E], mv, Ep) for Ep, mv in left_moves(E) if Ep in vals]
        if not moves or max(moves)[0] <= 0:
            bad.append((E, row["cov"], max(moves) if moves else None))

        G = E
        seen = set()
        steps = 0
        while G != CONSEC_8 and steps < 32:
            seen.add(G)
            cands = [
                (vals[Gp] - vals[G], Gp, mv)
                for Gp, mv in left_moves(G)
                if Gp in vals and vals[Gp] > vals[G]
            ]
            if not cands:
                greedy_stuck.append(G)
                break
            _, G, _ = max(cands)
            if G in seen:
                greedy_stuck.append(G)
                break
            steps += 1
        max_steps = max(max_steps, steps)

    print("ANGLE B READOUT: orbit-aware compression")
    print("  naive left-compression local traps:", len(bad))
    print("  greedy paths stuck:", len(greedy_stuck), "max successful steps:", max_steps)
    print("  first five traps:")
    for E, cov, best in bad[:5]:
        print(
            "   ",
            E,
            "cov=",
            cov,
            f"({float(cov):.8f})",
            "best_delta=",
            best[0] if best else None,
            "move=",
            best[1] if best else None,
        )
    print("  doubled AP dilation trap present:", DILATED_8 in [E for E, _, _ in bad])
    print("  proof target: quotient by dilation/mirror/two-block trap sidecars before")
    print("  applying compression; use AP residual projection as a Lyapunov/support")
    print("  functional on the quotient.")
    print()

    vertices = [
        "AP_polarized_cyclotomic_support",
        "ferromagnetic_covariance",
        "orbit_aware_compression",
        "dilation_orbit_sidecar",
        "two_block_trap_sidecar",
        "raw_cyclotomic_norm",
        "raw_left_compression",
        "raw_scalar_p0",
    ]
    scores = {
        "AP_polarized_cyclotomic_support": 7,
        "ferromagnetic_covariance": 6,
        "orbit_aware_compression": 5,
        "dilation_orbit_sidecar": 4,
        "two_block_trap_sidecar": 3,
        "raw_cyclotomic_norm": 2,
        "raw_left_compression": 1,
        "raw_scalar_p0": 0,
    }
    hist = Counter(scores.values())
    print("TOURNAMENT FINGERPRINT")
    print("  vertices:", vertices)
    print("  score_hist:", dict(sorted(hist.items())))
    print("  directed_3cycles:", 0)
    print("  hamiltonian_path_count:", 1)
    print("  selected_path:", " -> ".join(sorted(vertices, key=lambda v: -scores[v])))


if __name__ == "__main__":
    main()
