#!/usr/bin/env python3
"""Eight-core Q27 set-cover boundary scout for LRC(14).

This is the first below-nine-core extension of HYP-2465/HYP-2469.

The model is deliberately the same one used in
``lrc14_near_core_q27_setcover_codex.py``:

* start with the 7-core ``CORE=7*{1,...,12}``;
* delete a set ``D`` of core speeds;
* add at most ``|D|+1`` non-core speeds in the carry window ``1..13*84``;
* require primitivity by forcing at least one added non-7-multiple;
* ask whether the added speeds can cover every Q27 twist that is safe for the
  retained core.

For ``|D|=4`` this asks whether a primitive row retaining exactly eight core
speeds can block Q27.  If all 495 deletion sets are infeasible, HYP-2465
upgrades from "retain at least 9 core speeds" to "retain at least 8 core
speeds", leaving only rows that delete at least five core speeds, leave the
carry window, or evade the replacement-normalization model.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from statistics import mean
import sys
import time

from lrc14_near_core_q27_setcover_codex import (
    CANDIDATES,
    CORE_T,
    MAX_R,
    N,
    Q27,
    blocks,
    obligations,
    primitive_q27_cover_feasibility,
)
from lrc14_pisano_band_ladder_codex import (
    bprime_any,
    first_witness_in,
    min_witness_modulus,
    q_lattice,
)
from lrc14_lebesgue_wall_s676 import depth_sweep


def residue_signature(deleted: tuple[int, ...]) -> str:
    """Compact address of deleted core speeds."""
    k_values = tuple(v // 7 for v in deleted)
    mod3 = "".join(str(k % 3) for k in k_values)
    mod13 = ",".join(str(v % 13) for v in deleted)
    return f"k={k_values}; k_mod3={mod3}; v_mod13=({mod13})"


def coverage_snapshot(deleted: tuple[int, ...], top_n: int = 6) -> dict[str, object]:
    """Local coverage diagnostics for one deletion set.

    These diagnostics are not the proof certificate.  They help identify which
    candidate channels are trying hardest before the exact MILP rules them out.
    """
    base = tuple(v for v in CORE_T if v not in deleted)
    obs = obligations(base, Q27)
    rows: list[tuple[int, int]] = []
    best_13 = (0, 0)
    best_non13 = (0, 0)
    best_non7 = (0, 0)
    for v in CANDIDATES:
        count = sum(1 for q, a in obs if blocks(v, q, a))
        if count == 0:
            continue
        rows.append((count, v))
        if v % 13 == 0 and count > best_13[0]:
            best_13 = (count, v)
        if v % 13 != 0 and count > best_non13[0]:
            best_non13 = (count, v)
        if v % 7 != 0 and count > best_non7[0]:
            best_non7 = (count, v)
    rows.sort(reverse=True)
    max_single = rows[0][0] if rows else 0
    lower_bound = (len(obs) + max_single - 1) // max_single if max_single else None
    return {
        "deleted": deleted,
        "obligations": len(obs),
        "active_candidates": len(rows),
        "max_single_lower_bound": lower_bound,
        "best": tuple((v, c) for c, v in rows[:top_n]),
        "best_13": (best_13[1], best_13[0]),
        "best_non13": (best_non13[1], best_non13[0]),
        "best_non7": (best_non7[1], best_non7[0]),
    }


def completed_row(deleted: tuple[int, ...], chosen: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(tuple(v for v in CORE_T if v not in deleted) + tuple(chosen)))


def opening_diagnostics(deleted: tuple[int, ...], chosen: tuple[int, ...]) -> dict[str, object]:
    row = completed_row(deleted, chosen)
    sweep = depth_sweep(row)
    return {
        "row": row,
        "plain_min_q": min_witness_modulus(list(row), 1200),
        "q41": first_witness_in(list(row), q_lattice(41)),
        "bprime_any": bprime_any(list(row)),
        "p0": sweep.p0,
        "components": sweep.positive_safe_components,
        "safe_points": len(sweep.safe_points),
    }


def run_boundary(
    delete_count: int,
    add_budget: int,
    time_limit: int,
    limit: int | None,
    start_index: int,
    stop_on_feasible: bool,
) -> list[object]:
    combos = list(combinations(CORE_T, delete_count))
    combos = combos[start_index:]
    if limit is not None:
        combos = combos[:limit]
    results = []
    start = time.time()
    for index, deleted in enumerate(combos, start=1):
        result = primitive_q27_cover_feasibility(tuple(deleted), add_budget, time_limit)
        results.append(result)
        print(
            f"progress {index:3d}/{len(combos):3d}: deleted={deleted} "
            f"obs={result.obligation_count} feasible={result.feasible} "
            f"elapsed={result.elapsed_s:.2f}s",
            file=sys.stderr,
            flush=True,
        )
        if result.feasible is True:
            print(f"first feasible set found: deleted={deleted} chosen={result.chosen}", file=sys.stderr)
            if stop_on_feasible:
                break
    print(f"progress total wall time {time.time() - start:.2f}s", file=sys.stderr)
    return results


def print_summary(results: list[object], delete_count: int, add_budget: int, time_limit: int) -> None:
    feasible = [r for r in results if r.feasible is True]
    unknown = [r for r in results if r.feasible is None]
    infeasible = [r for r in results if r.feasible is False]
    obs = [r.obligation_count for r in results]
    elapsed = [r.elapsed_s for r in results]
    print("A. Exact Q27 set-cover boundary")
    print(f"  delete_count={delete_count}; add_budget={add_budget}; time_limit={time_limit}s")
    print(f"  candidate_window=1..{MAX_R}; n={N}; Q27_size={len(Q27)}")
    print(
        f"  cases={len(results)} infeasible={len(infeasible)} "
        f"feasible={len(feasible)} unknown={len(unknown)}"
    )
    print(
        f"  obligations={min(obs)}..{max(obs)} mean={mean(obs):.2f}; "
        f"solver_elapsed={sum(elapsed):.2f}s; max_case={max(elapsed):.2f}s"
    )
    if feasible:
        print("  feasible examples:")
        for r in feasible[:8]:
            diag = opening_diagnostics(r.deleted, r.chosen)
            print(f"    deleted={r.deleted} chosen={r.chosen} message={r.message}")
            print(
                f"      row={diag['row']} plain_min_q={diag['plain_min_q']} "
                f"Q41={diag['q41']} Bprime={diag['bprime_any']} "
                f"p0={diag['p0']} components={diag['components']} safe_points={diag['safe_points']}"
            )
    if unknown:
        print("  unknown examples:")
        for r in unknown[:8]:
            print(f"    deleted={r.deleted} chosen={r.chosen} message={r.message}")
    print()

    print("B. Hardest deletion addresses")
    ranked_by_obs = sorted(results, key=lambda r: (-r.obligation_count, r.deleted))[:8]
    ranked_by_time = sorted(results, key=lambda r: (-r.elapsed_s, r.deleted))[:8]
    print("  highest obligation loads:")
    for r in ranked_by_obs:
        print(f"    deleted={r.deleted} obligations={r.obligation_count} {residue_signature(r.deleted)}")
    print("  slowest solver cases:")
    for r in ranked_by_time:
        print(f"    deleted={r.deleted} elapsed={r.elapsed_s:.2f}s obligations={r.obligation_count}")
    print()

    print("C. Coverage snapshots for high-load cases")
    for r in ranked_by_obs[:3]:
        snap = coverage_snapshot(r.deleted)
        print(
            f"  deleted={snap['deleted']} obligations={snap['obligations']} "
            f"active_candidates={snap['active_candidates']} "
            f"single_cover_lb={snap['max_single_lower_bound']}"
        )
        print(f"    best={snap['best']}")
        print(
            f"    best_13={snap['best_13']} best_non13={snap['best_non13']} "
            f"best_non7={snap['best_non7']}"
        )
    print()


def print_tournament_analysis() -> None:
    print("D. Tournament Analysis and assumption challenge")
    vertices = [
        ("Q27_obligation_hypergraph", (5, 5, 5, 5, 4, 5)),
        ("deleted_core_addresses", (5, 4, 5, 4, 5, 4)),
        ("candidate_speed_cover_masks", (5, 4, 4, 5, 5, 3)),
        ("13_clock_resource", (4, 5, 4, 4, 4, 4)),
        ("divisor_fiber_residue_classes", (4, 4, 4, 4, 5, 4)),
        ("support_load_owner_channels", (4, 5, 5, 3, 3, 5)),
        ("raw_runner_geometry", (3, 3, 3, 5, 5, 3)),
        ("outside_window_normalizer", (3, 5, 5, 3, 3, 5)),
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
    print("  observable=(exact finite control, descent leverage, retained side-channel,")
    print("              geometric fidelity, computation, proof-transfer leverage)")
    print(f"  directed_3cycles={cycles}; edge_flips_vs_list_order={flips}/28")
    for name in sorted(names, key=lambda x: (-len(out[x]), x)):
        print(f"    {name:31s} out={len(out[name])} scores={scores[name]}")
    print()
    print("  Assumption challenge:")
    print("    considered vertices: runners, gaps, fixed circle sections, section boundaries,")
    print("    wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,")
    print("    deletion addresses, candidate masks, and proof obligations.")
    print("    selected quotient: Q27 obligations with deleted-core/candidate-cover side channels.")
    print("    preserves: exact bounded Q27 feasibility and the Church-style retained channel.")
    print("    destroys: outside-window geometry, time-continuous owner motion, and arbitrary")
    print("    multi-stranger interactions not reducible to the replacement-normalization model.")
    print()


def print_takeaway(results: list[object], delete_count: int) -> None:
    feasible = [r for r in results if r.feasible is True]
    unknown = [r for r in results if r.feasible is None]
    print("E. Proof-shaped takeaway")
    if not feasible and not unknown and len(results) == 495 and delete_count == 4:
        print("  Exact bounded upgrade: in the carry window, every primitive replacement row")
        print("  retaining at least 8 of the 12 core speeds has a Q27 witness.")
        print("  Therefore a normalized no-Q27 row must delete at least 5 core speeds.")
        print("  The Church-style descent target narrows from below-nine-core to below-eight-core,")
        print("  unless the row leaves the carry window or the replacement normalization fails.")
    elif feasible:
        print("  Feasible Q27 blockers appear at four core deletions, so the universal")
        print("  eight-core Q27 set-cover lemma is false as stated.  The useful theorem")
        print("  should instead classify every feasible packet by its opening channel:")
        print("  plain-shell witness, Bprime/owner opening, positive measure, or descent.")
    elif unknown:
        print("  Some cases were unknown at this time limit.  Rerun with a larger limit before")
        print("  treating this boundary as a certificate.")
    else:
        print("  This was a partial run.  Use the summary as a scout, not a certificate.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--delete-count", type=int, default=4)
    parser.add_argument("--add-budget", type=int, default=5)
    parser.add_argument("--time-limit", type=int, default=12)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--stop-on-feasible", action="store_true")
    args = parser.parse_args()

    print("=" * 78)
    print("Codex LRC14 eight-core Q27 set-cover boundary")
    print("=" * 78)
    print()
    results = run_boundary(
        args.delete_count,
        args.add_budget,
        args.time_limit,
        args.limit,
        args.start_index,
        args.stop_on_feasible,
    )
    print_summary(results, args.delete_count, args.add_budget, args.time_limit)
    print_tournament_analysis()
    print_takeaway(results, args.delete_count)


if __name__ == "__main__":
    main()
