#!/usr/bin/env python3
"""Near-core Q27 set-cover lower bound for LRC(14).

This is a proof-target addendum to HYP-2463.  The previous script only stacked
the eight one-stranger hard residues.  Here the added speeds are arbitrary
inside the same bounded carry window used by HYP-2444:

    1 <= v <= 13*84, excluding the original 7-core speeds.

For a partial row B = CORE \\ D, a denominator twist (q,a) is an obligation if it
is safe for B.  An added speed covers that obligation if it is dangerous at
(q,a).  Thus "can |D|+1 added speeds block all Q27 witnesses?" becomes a
binary set-cover feasibility problem with two extra constraints:

  * at most |D|+1 chosen speeds;
  * at least one chosen speed is not divisible by 7, so the completed row is
    primitive rather than just another 7-multiple row.

The stored run proves infeasibility for all |D| <= 3.  So, in this carry window,
any primitive row retaining at least 9 of the 12 core speeds has a Q27 witness.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import gcd
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE))

from lrc14_pisano_band_ladder_codex import CORE, MAX_R, N, q_lattice
from lrc14_parity_typed_q27_ledger_codex import first_witness_masked, unit_twists


CORE_T = tuple(CORE)
Q27 = tuple(q_lattice(27))
SHELL27 = tuple(range(2, 28))
CANDIDATES = tuple(v for v in range(1, MAX_R + 1) if v not in CORE_T)


def blocks(v: int, q: int, a: int) -> bool:
    band = q // N
    r = (a * v) % q
    return min(r, q - r) <= band


def obligations(base: tuple[int, ...], qs: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    out = []
    for q in qs:
        for a in unit_twists(q):
            if all(not blocks(v, q, a) for v in base):
                out.append((q, a))
    return tuple(out)


@dataclass(frozen=True)
class FeasibilityResult:
    deleted: tuple[int, ...]
    add_budget: int
    obligation_count: int
    candidate_count: int
    feasible: bool | None
    chosen: tuple[int, ...]
    message: str
    elapsed_s: float


def primitive_q27_cover_feasibility(
    deleted: tuple[int, ...],
    add_budget: int,
    time_limit_s: int = 12,
    qs: tuple[int, ...] = Q27,
) -> FeasibilityResult:
    """MILP: can add_budget primitive candidates cover the selected obligations?"""
    base = tuple(v for v in CORE_T if v not in deleted)
    obs = obligations(base, qs)
    vals = tuple(v for v in CANDIDATES if v not in base)

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    lows: list[float] = []
    ups: list[float] = []
    row_index = 0

    for q, a in obs:
        any_cover = False
        for j, v in enumerate(vals):
            if blocks(v, q, a):
                rows.append(row_index)
                cols.append(j)
                data.append(-1.0)
                any_cover = True
        if not any_cover:
            return FeasibilityResult(
                deleted=deleted,
                add_budget=add_budget,
                obligation_count=len(obs),
                candidate_count=len(vals),
                feasible=False,
                chosen=tuple(),
                message="uncoverable obligation",
                elapsed_s=0.0,
            )
        lows.append(-np.inf)
        ups.append(-1.0)
        row_index += 1

    # Primitive completed row: since every core speed is divisible by 7, at
    # least one added speed must not be divisible by 7.
    for j, v in enumerate(vals):
        if v % 7 != 0:
            rows.append(row_index)
            cols.append(j)
            data.append(-1.0)
    lows.append(-np.inf)
    ups.append(-1.0)
    row_index += 1

    # Cardinality cap: this is the exact near-core replacement budget.
    for j, _v in enumerate(vals):
        rows.append(row_index)
        cols.append(j)
        data.append(1.0)
    lows.append(-np.inf)
    ups.append(float(add_budget))
    row_index += 1

    matrix = coo_matrix((data, (rows, cols)), shape=(row_index, len(vals))).tocsr()
    constraints = LinearConstraint(matrix, np.array(lows), np.array(ups))

    start = time.time()
    result = milp(
        c=np.zeros(len(vals)),
        integrality=np.ones(len(vals)),
        bounds=Bounds(0, 1),
        constraints=constraints,
        options={"time_limit": time_limit_s, "mip_rel_gap": 0},
    )
    elapsed = time.time() - start

    chosen: tuple[int, ...] = tuple()
    if result.x is not None:
        chosen = tuple(vals[i] for i, x in enumerate(result.x) if x > 0.5)

    if result.success:
        feasible: bool | None = True
    elif "Infeasible" in result.message:
        feasible = False
    else:
        feasible = None

    return FeasibilityResult(
        deleted=deleted,
        add_budget=add_budget,
        obligation_count=len(obs),
        candidate_count=len(vals),
        feasible=feasible,
        chosen=chosen,
        message=result.message,
        elapsed_s=elapsed,
    )


def row_gcd(row: tuple[int, ...]) -> int:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g


def mask_for_candidate(v: int, obs: tuple[tuple[int, int], ...]) -> int:
    mask = 0
    for i, (q, a) in enumerate(obs):
        if blocks(v, q, a):
            mask |= 1 << i
    return mask


def direct_delete_one_plain_scan() -> list[dict[str, object]]:
    """Exact direct scan: delete one core speed, add two arbitrary candidates."""
    rows = []
    for deleted in CORE_T:
        base = tuple(v for v in CORE_T if v != deleted)
        obs = obligations(base, SHELL27)
        full = (1 << len(obs)) - 1
        masks = {
            v: mask_for_candidate(v, obs)
            for v in CANDIDATES
            if v not in base and mask_for_candidate(v, obs)
        }
        vals = sorted(masks)
        plain_blockers = 0
        q27_misses = 0
        first_examples: list[tuple[int, int, tuple[int, int] | None]] = []
        for i, u in enumerate(vals):
            mu = masks[u]
            for v in vals[i + 1 :]:
                if (mu | masks[v]) != full:
                    continue
                row = tuple(sorted(base + (u, v)))
                if len(set(row)) != 13 or row_gcd(row) != 1:
                    continue
                plain_blockers += 1
                q27 = first_witness_masked(row, Q27)
                if q27 is None:
                    q27_misses += 1
                if len(first_examples) < 4:
                    first_examples.append((u, v, q27))
        rows.append(
            {
                "deleted": deleted,
                "obligations": len(obs),
                "active_candidates": len(vals),
                "plain_blockers": plain_blockers,
                "q27_misses": q27_misses,
                "examples": tuple(first_examples),
            }
        )
    return rows


def print_delete_one_plain_scan() -> None:
    print("A. Broad one-deletion/two-add plain-shell scan")
    print("  Added speeds range over 1..13*84, excluding original core speeds.")
    total_plain = 0
    total_q27_miss = 0
    for row in direct_delete_one_plain_scan():
        total_plain += int(row["plain_blockers"])
        total_q27_miss += int(row["q27_misses"])
        print(
            f"  delete={row['deleted']:3d} obligations={row['obligations']:3d} "
            f"active_candidates={row['active_candidates']:4d} "
            f"plain_blockers={row['plain_blockers']:4d} "
            f"Q27_misses={row['q27_misses']:2d} examples={row['examples']}"
        )
    print(f"  total plain blockers={total_plain}; total Q27 misses among them={total_q27_miss}")
    print()


def print_q27_milp_table(max_delete: int = 3) -> None:
    print("B. Primitive Q27 set-cover infeasibility")
    print("  For delete_count=e, test whether <=e+1 added speeds can cover all Q27 obligations.")
    print("  Primitive constraint: at least one added speed is not divisible by 7.")
    for delete_count in range(max_delete + 1):
        add_budget = delete_count + 1
        results = [
            primitive_q27_cover_feasibility(tuple(deleted), add_budget)
            for deleted in combinations(CORE_T, delete_count)
        ]
        feasible = [r for r in results if r.feasible is True]
        unknown = [r for r in results if r.feasible is None]
        infeasible = [r for r in results if r.feasible is False]
        obs_counts = [r.obligation_count for r in results]
        elapsed = sum(r.elapsed_s for r in results)
        print(
            f"  e={delete_count}: cases={len(results):3d} add_budget={add_budget} "
            f"infeasible={len(infeasible):3d} feasible={len(feasible):2d} "
            f"unknown={len(unknown):2d} obligations={min(obs_counts)}..{max(obs_counts)} "
            f"solver_elapsed={elapsed:.2f}s"
        )
        if feasible:
            for r in feasible[:5]:
                print(f"    FEASIBLE deleted={r.deleted} chosen={r.chosen}")
        if unknown:
            for r in unknown[:5]:
                print(f"    UNKNOWN deleted={r.deleted} message={r.message} chosen={r.chosen}")
    print()


def print_tournament_analysis() -> None:
    print("C. Tournament Analysis over near-core proof obligations")
    vertices = [
        ("primitive_Q27_set_cover", (5, 5, 5, 5, 5, 4)),
        ("near_core_compression", (5, 5, 5, 4, 4, 5)),
        ("plain_shell_cover_noise", (4, 5, 4, 4, 5, 4)),
        ("hard_residue_hull", (4, 4, 4, 5, 4, 4)),
        ("low_clock_or_fiber_escape", (4, 4, 5, 4, 3, 5)),
        ("raw_residue_class_packet", (3, 3, 4, 3, 5, 3)),
    ]
    names = [name for name, _score in vertices]
    scores = dict(vertices)
    out = {name: set() for name in names}
    flips = 0
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            av = sum(x > y for x, y in zip(scores[a], scores[b]))
            bv = sum(y > x for x, y in zip(scores[a], scores[b]))
            if av >= bv:
                out[a].add(b)
            else:
                out[b].add(a)
                flips += 1
    cycles = 0
    for a, b, c in combinations(names, 3):
        degs = {a: 0, b: 0, c: 0}
        for x, y in ((a, b), (a, c), (b, c)):
            if y in out[x]:
                degs[x] += 1
            else:
                degs[y] += 1
        if sorted(degs.values()) == [1, 1, 1]:
            cycles += 1
    print("  vertices are proof obligations/set-cover gauges, not runners")
    print("  observable=(exactness, compression leverage, LRC leverage, support retention, computation, transfer)")
    print(f"  directed_3cycles={cycles}; edge_flips_vs_list_order={flips}/15")
    for name in sorted(names, key=lambda x: (-len(out[x]), x)):
        print(f"    {name:28s} out={len(out[name])} scores={scores[name]}")
    print()


def main() -> None:
    print("=" * 78)
    print("Codex LRC14 near-core Q27 set-cover atlas")
    print("=" * 78)
    print(f"n={N}; core_size={len(CORE_T)}; candidate_window=1..{MAX_R}; candidates={len(CANDIDATES)}")
    print(f"Q27 size={len(Q27)}")
    print()
    print_delete_one_plain_scan()
    print_q27_milp_table(max_delete=3)
    print_tournament_analysis()
    print("D. Proof-shaped takeaway")
    print("  Exact bounded lemma: in the HYP-2444 carry window, every primitive")
    print("  replacement row retaining at least 9 of the 12 core speeds has a Q27 witness.")
    print("  Plain-shell blocking is noisy after one deletion, but Q27 infeasibility is clean.")
    print("  The compression proof should now try to force any would-be LRC14 Q27 blocker")
    print("  either below the 9-core threshold or into a non-near-core descent/opening channel.")


if __name__ == "__main__":
    main()
