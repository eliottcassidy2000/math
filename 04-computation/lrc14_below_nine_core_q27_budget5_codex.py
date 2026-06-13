#!/usr/bin/env python3
"""Below-nine-core Q27 budget-5 scout for LRC(14).

HYP-2465 proves that, inside the HYP-2444 carry window 1..13*84,
every primitive replacement row retaining at least 9 of the 12 speeds

    CORE = 7*{1,...,12}

has a Q27 witness.  HYP-2469 names the next live portal as the first
below-nine-core layer: delete four core speeds and add five arbitrary
window speeds.

This script runs that exact first portal as a set-cover problem:

    B = CORE \\ D, |D| = 4
    obligations = Q27 twists safe for B
    added speeds cover an obligation by being dangerous there

Can five primitive added speeds cover all Q27 obligations?  If not, every
primitive 13-speed row in the carry window retaining at least 8 core speeds
has a Q27 witness.

Tournament Analysis / assumption challenge:
  Vertices are proof portals, not runners.  Candidate vertex sets considered
  include runners, deleted core speeds, added speeds, Q27 twists, support-load
  obligations, Bprime owner exits, outside-window carry fibers, and proof
  obligations.  The selected quotient preserves the finite descent decision
  tree and destroys raw time geometry and arbitrary speed magnitudes.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import sys
import time


HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE))

from lrc14_near_core_q27_setcover_codex import (  # noqa: E402
    CANDIDATES,
    CORE_T,
    Q27,
    blocks,
    obligations,
    primitive_q27_cover_feasibility,
)
from lrc14_pisano_band_ladder_codex import antipodal_class  # noqa: E402


HARD = (260, 351, 442, 611, 702, 793, 962, 1053)
ANCHOR_CORE = (7, 21, 49)


@dataclass(frozen=True)
class CaseResult:
    deleted: tuple[int, ...]
    feasible: bool | None
    obligation_count: int
    candidate_count: int
    elapsed_s: float
    message: str
    chosen: tuple[int, ...]

    @property
    def deleted_quotients(self) -> tuple[int, ...]:
        return tuple(v // 7 for v in self.deleted)

    @property
    def anchor_signature(self) -> str:
        kept = tuple(v for v in ANCHOR_CORE if v not in self.deleted)
        return "keeps_" + "_".join(str(v) for v in kept) if kept else "deletes_7_21_49"


def speed_tag(v: int) -> str:
    parts: list[str] = []
    if v in HARD:
        parts.append("hard")
    if v % 13 == 0:
        parts.append("clock13")
    if v % 7 == 0:
        parts.append("mult7")
    c27 = antipodal_class(v, 27)
    if v % 27 == 0:
        parts.append("zero27")
    elif c27 == 10:
        parts.append("pm10")
    elif c27 is None:
        parts.append("nonunit27")
    else:
        parts.append(f"c27={c27}")
    if v % 91 == 0:
        parts.append("q91")
    return "+".join(parts)


def top_covering_speeds(deleted: tuple[int, ...], limit: int = 8) -> tuple[tuple[int, int, str], ...]:
    base = tuple(v for v in CORE_T if v not in deleted)
    obs = obligations(base, Q27)
    vals = tuple(v for v in CANDIDATES if v not in base)
    scored: list[tuple[int, int, str]] = []
    for v in vals:
        count = sum(1 for q, a in obs if blocks(v, q, a))
        if count:
            scored.append((count, v, speed_tag(v)))
    scored.sort(key=lambda x: (-x[0], x[1]))
    return tuple((v, count, tag) for count, v, tag in scored[:limit])


def run_delete4(time_limit_s: int = 20, progress_every: int = 25) -> list[CaseResult]:
    results: list[CaseResult] = []
    cases = list(combinations(CORE_T, 4))
    start = time.time()
    for idx, deleted in enumerate(cases, 1):
        res = primitive_q27_cover_feasibility(tuple(deleted), 5, time_limit_s=time_limit_s)
        results.append(
            CaseResult(
                deleted=tuple(deleted),
                feasible=res.feasible,
                obligation_count=res.obligation_count,
                candidate_count=res.candidate_count,
                elapsed_s=res.elapsed_s,
                message=res.message,
                chosen=res.chosen,
            )
        )
        if progress_every and idx % progress_every == 0:
            infeasible = sum(1 for r in results if r.feasible is False)
            feasible = sum(1 for r in results if r.feasible is True)
            unknown = sum(1 for r in results if r.feasible is None)
            print(
                f"  progress {idx:3d}/{len(cases)} "
                f"infeasible={infeasible:3d} feasible={feasible:2d} unknown={unknown:2d} "
                f"elapsed={time.time() - start:.1f}s",
                flush=True,
            )
    return results


def tournament_fingerprint() -> tuple[dict[int, int], int, int, list[tuple[str, int, tuple[int, ...]]]]:
    vertices: list[tuple[str, tuple[int, ...]]] = [
        ("Q27_budget5_setcover", (5, 5, 5, 5, 4, 5)),
        ("eight_core_window_barrier", (5, 5, 5, 4, 4, 5)),
        ("support_load_descent", (4, 5, 5, 4, 3, 5)),
        ("outside_window_normalizer", (3, 4, 5, 4, 3, 5)),
        ("owner_Bprime_exception", (4, 4, 4, 5, 3, 4)),
        ("plain_shell_scalar", (2, 3, 3, 2, 5, 2)),
        ("raw_denominator_search", (1, 2, 3, 1, 5, 1)),
    ]
    n = len(vertices)
    out = [0] * n
    flips = 0
    for i in range(n):
        for j in range(i + 1, n):
            if vertices[i][1] >= vertices[j][1]:
                out[i] += 1
            else:
                out[j] += 1
                flips += 1
    cycles = 0
    for i, j, k in combinations(range(n), 3):
        edges = {
            (i, j): vertices[i][1] >= vertices[j][1],
            (i, k): vertices[i][1] >= vertices[k][1],
            (j, k): vertices[j][1] >= vertices[k][1],
        }
        # With vertices ordered by declaration, a 3-cycle occurs exactly when
        # neither all three edges follow declaration order nor all can be made
        # transitive by one middle vertex.
        wins = defaultdict(set)
        for a, b in ((i, j), (i, k), (j, k)):
            if edges[(a, b)]:
                wins[a].add(b)
            else:
                wins[b].add(a)
        if sorted(len(wins[x]) for x in (i, j, k)) == [1, 1, 1]:
            cycles += 1
    hist = dict(sorted(Counter(out).items()))
    ranked = sorted(
        [(name, out[idx], score) for idx, (name, score) in enumerate(vertices)],
        key=lambda row: (-row[1], row[0]),
    )
    return hist, cycles, flips, ranked


def print_summary(results: list[CaseResult]) -> None:
    feasible = [r for r in results if r.feasible is True]
    unknown = [r for r in results if r.feasible is None]
    infeasible = [r for r in results if r.feasible is False]
    obs = [r.obligation_count for r in results]
    elapsed = [r.elapsed_s for r in results]

    print()
    print("A. Exact delete-count 4 Q27 budget-5 result")
    print("  row form: (CORE \\ D) union A, |D|=4, |A|<=5, primitive")
    print(f"  cases={len(results)} infeasible={len(infeasible)} feasible={len(feasible)} unknown={len(unknown)}")
    print(f"  obligations={min(obs)}..{max(obs)} avg={sum(obs) / len(obs):.1f}")
    print(f"  solver_elapsed_total={sum(elapsed):.2f}s max_case={max(elapsed):.2f}s")
    if feasible:
        print(f"  feasible examples={[r.deleted for r in feasible[:5]]}")
    if unknown:
        print(f"  unknown examples={[(r.deleted, r.message) for r in unknown[:5]]}")
    print()

    print("B. Deletion-shape diagnostics")
    by_anchor = Counter(r.anchor_signature for r in results)
    print(f"  anchor signatures={dict(sorted(by_anchor.items()))}")
    by_first = Counter(r.deleted_quotients[0] for r in results)
    print(f"  first deleted quotient histogram={dict(sorted(by_first.items()))}")
    by_sum = Counter(sum(r.deleted_quotients) for r in results)
    print(f"  deleted quotient sum range={min(by_sum)}..{max(by_sum)}")
    print("  quotient-sum extremes:")
    for r in sorted(results, key=lambda x: (sum(x.deleted_quotients), x.deleted))[:3]:
        print(f"    low sum {sum(r.deleted_quotients):2d}: D={r.deleted} obligations={r.obligation_count}")
    for r in sorted(results, key=lambda x: (-sum(x.deleted_quotients), x.deleted))[:3]:
        print(f"    high sum {sum(r.deleted_quotients):2d}: D={r.deleted} obligations={r.obligation_count}")
    print()

    print("C. Hardest and extremal cases")
    notable: list[CaseResult] = []
    notable.extend(sorted(results, key=lambda r: -r.elapsed_s)[:3])
    notable.extend(sorted(results, key=lambda r: -r.obligation_count)[:2])
    notable.extend(sorted(results, key=lambda r: r.obligation_count)[:2])
    seen: set[tuple[int, ...]] = set()
    for r in notable:
        if r.deleted in seen:
            continue
        seen.add(r.deleted)
        print(
            f"  D={r.deleted} quotients={r.deleted_quotients} "
            f"obligations={r.obligation_count} elapsed={r.elapsed_s:.2f}s"
        )
        print(f"    top single covers={top_covering_speeds(r.deleted, 6)}")
    print()

    print("D. Proof-shaped consequence")
    if len(infeasible) == len(results):
        print("  In the carry window 1..1092, no primitive row retaining >=8 core speeds")
        print("  can block Q27.  HYP-2465's nine-core barrier extends one full layer.")
        print("  The below-nine-core portal has moved to a below-eight-core portal:")
        print("  any bounded Q27 blocker must delete at least five core speeds, leave")
        print("  the carry window, or trigger a named side-channel exception/descent.")
    else:
        print("  The first below-nine-core layer is not closed by raw Q27 budget.")
        print("  Inspect feasible/unknown cases as candidate exceptions or typed-budget gates.")
    print()

    hist, cycles, flips, ranked = tournament_fingerprint()
    print("E. Tournament Analysis over proof portals")
    print("  observable=(exactness, support retention, LRC leverage, side-channel, computability, descent)")
    print("  switch/gauge=rank proof portals by retained-channel strength, not scalar q size")
    print("  tie Hamiltonian path=declaration order")
    print(f"  score histogram={hist}; directed_3cycles={cycles}; edge_flips_vs_declaration={flips}/21")
    for name, out, score in ranked:
        print(f"    {out}: {name:28s} score={score}")


def main() -> None:
    print("=" * 80)
    print("Codex LRC14 below-nine-core Q27 budget-5 scout")
    print("=" * 80)
    print(f"core={CORE_T}")
    print(f"Q27 size={len(Q27)}; candidate_window=1..{13 * 84}; candidates={len(CANDIDATES)}")
    print("Running exact MILP feasibility for all |D|=4 cases...")
    results = run_delete4()
    print_summary(results)


if __name__ == "__main__":
    main()
