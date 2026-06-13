#!/usr/bin/env python3
"""Codex 2026-06-12: marked ladder support for the LRC(14) proof route.

This script is a support-gate reframing of HYP-2438.  The scalar question
"does a clock q work?" is replaced by a marked support question:

    if q is blocked, which runners paid for blocking every unit twist a/q?

For total runner count n=14, delta=1/14.  At denominator q and unit twist a,
a speed v is dangerous when

    ||a v / q|| <= 1/14.

Thus q is a witness iff some unit a avoids every runner's danger band.  If q is
blocked, the set of dangerous runners over all unit twists is a tiny hitting-set
problem.  The minimum hitting-set size is a local resource cost for blocking q.

The point is not to brute-force LRC(14).  It is to find the missing marked
support object, analogous to the recent Type-II-code and Erdos-Moser support
gates: scalar shadows are cheap; marked support is load-bearing.

Tournament Analysis / assumption challenge:
  Vertices are proof routes/marked blocker types, not runners.  Candidate
  vertex sets considered: runners, denominators q, unit twists a, danger bands,
  blocker hitting sets, apex speeds, stranger speeds, safe intervals, and proof
  obligations.  The selected quotient preserves "how expensive is it to block
  this ladder rung?" and destroys raw time geometry except for the band index.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path
import random


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "05-knowledge" / "results" / "lrc14_marked_ladder_support_codex.out"

N = 14
K = N - 1
DELTA = Fraction(1, N)
DIVISORS = (1, 2, 7, 14)

AP = tuple(range(1, N))
VSTAR = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)
TWOP = tuple(2 * v for v in AP)
EVADER_RS = (611, 702, 793, 962, 1053)


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def row_gcd(row: tuple[int, ...]) -> int:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g


def norm_residue(r: int, q: int) -> int:
    r %= q
    return min(r, q - r)


def band(q: int) -> int:
    return q // N


def denominators(m_max: int) -> tuple[int, ...]:
    return tuple(sorted({d * m for d in DIVISORS for m in range(1, m_max + 1) if d * m > 1}))


def shell_denominators(q_max: int) -> tuple[int, ...]:
    return tuple(range(2, q_max + 1))


def unit_twists(q: int) -> tuple[int, ...]:
    return tuple(a for a in range(1, q) if gcd(a, q) == 1)


def blocks_twist(v: int, a: int, q: int) -> bool:
    return norm_residue(a * v, q) <= band(q)


def safe_twists(row: tuple[int, ...], q: int) -> tuple[int, ...]:
    return tuple(a for a in unit_twists(q) if all(not blocks_twist(v, a, q) for v in row))


def first_witness(row: tuple[int, ...], qs: tuple[int, ...]) -> tuple[int, int] | None:
    for q in qs:
        safe = safe_twists(row, q)
        if safe:
            return q, safe[0]
    return None


@dataclass(frozen=True)
class BlockLedger:
    q: int
    twist_count: int
    band: int
    safe_count: int
    min_blockers: int
    canonical_cover: tuple[int, ...]
    universal_blockers: tuple[int, ...]
    cover_twist_counts: tuple[tuple[int, int], ...]

    @property
    def blocked(self) -> bool:
        return self.safe_count == 0

    @property
    def has_universal(self) -> bool:
        return bool(self.universal_blockers)

    @property
    def cover_speeds(self) -> tuple[int, ...]:
        return self.canonical_cover


def block_ledger(row: tuple[int, ...], q: int) -> BlockLedger:
    twists = unit_twists(q)
    safe_count = 0
    masks: list[int] = []
    for v in row:
        mask = 0
        for j, a in enumerate(twists):
            if blocks_twist(v, a, q):
                mask |= 1 << j
        masks.append(mask)
    full = (1 << len(twists)) - 1
    for j, _a in enumerate(twists):
        bit = 1 << j
        if not any(mask & bit for mask in masks):
            safe_count += 1

    universal = tuple(row[i] for i, mask in enumerate(masks) if mask == full)
    if safe_count:
        min_size = 0
        cover: tuple[int, ...] = tuple()
    else:
        min_size = K + 1
        cover = tuple()
        indices = range(len(row))
        for size in range(1, len(row) + 1):
            found = False
            for combo in combinations(indices, size):
                merged = 0
                for i in combo:
                    merged |= masks[i]
                if merged == full:
                    min_size = size
                    cover = tuple(row[i] for i in combo)
                    found = True
                    break
            if found:
                break

    cover_counts = tuple(
        sorted(
            ((row[i], masks[i].bit_count()) for i in range(len(row)) if masks[i]),
            key=lambda item: (-item[1], item[0]),
        )
    )
    return BlockLedger(
        q=q,
        twist_count=len(twists),
        band=band(q),
        safe_count=safe_count,
        min_blockers=min_size if safe_count == 0 else 0,
        canonical_cover=cover,
        universal_blockers=universal,
        cover_twist_counts=cover_counts,
    )


@dataclass(frozen=True)
class RowSummary:
    name: str
    row: tuple[int, ...]
    m_max: int
    witness: tuple[int, int] | None
    blocked_count: int
    blocked_by_band: tuple[tuple[int, int], ...]
    min_blocker_hist: tuple[tuple[int, int], ...]
    universal_count: int
    no_universal_count: int
    max_min_blockers: int
    max_min_qs: tuple[int, ...]
    cover_load: tuple[tuple[int, int], ...]


def summarize_row(name: str, row: tuple[int, ...], m_max: int) -> RowSummary:
    qs = denominators(m_max)
    ledgers = [block_ledger(row, q) for q in qs]
    blocked = [ledger for ledger in ledgers if ledger.blocked]
    blocked_by_band = Counter(ledger.band for ledger in blocked)
    hist = Counter(ledger.min_blockers for ledger in blocked)
    universal_count = sum(1 for ledger in blocked if ledger.has_universal)
    no_universal = len(blocked) - universal_count
    max_min = max((ledger.min_blockers for ledger in blocked), default=0)
    max_qs = tuple(ledger.q for ledger in blocked if ledger.min_blockers == max_min)
    load: Counter[int] = Counter()
    for ledger in blocked:
        for v in ledger.cover_speeds:
            load[v] += 1
    return RowSummary(
        name=name,
        row=row,
        m_max=m_max,
        witness=first_witness(row, qs),
        blocked_count=len(blocked),
        blocked_by_band=tuple(sorted(blocked_by_band.items())),
        min_blocker_hist=tuple(sorted(hist.items())),
        universal_count=universal_count,
        no_universal_count=no_universal,
        max_min_blockers=max_min,
        max_min_qs=max_qs,
        cover_load=tuple(load.most_common(8)),
    )


def print_row_summary(summary: RowSummary) -> None:
    w = "none" if summary.witness is None else f"{summary.witness[1]}/{summary.witness[0]}"
    print(f"  {summary.name:22s} m<={summary.m_max:<2d} witness={w:>7s} blocked={summary.blocked_count:3d}")
    print(
        f"    bands={dict(summary.blocked_by_band)} "
        f"min_cover_hist={dict(summary.min_blocker_hist)}"
    )
    print(
        f"    universal_q={summary.universal_count} "
        f"cover_q={summary.no_universal_count} "
        f"max_min_blockers={summary.max_min_blockers} at q={list(summary.max_min_qs[:10])}"
    )
    print(f"    canonical cover load top={summary.cover_load}")


def family_row(r: int) -> tuple[int, ...]:
    return tuple(sorted([7 * k for k in range(1, 13)] + [r]))


def evader_scan(limit: int = 1200) -> list[tuple[int, tuple[int, int] | None, tuple[int, int] | None, tuple[int, int] | None]]:
    shell27 = shell_denominators(27)
    shell41 = shell_denominators(41)
    fiber27 = denominators(27)
    out = []
    for r in range(1, limit):
        row = family_row(r)
        if len(set(row)) != K or row_gcd(row) != 1 or r % 7 == 0:
            continue
        w_shell27 = first_witness(row, shell27)
        w_shell41 = first_witness(row, shell41)
        w_fiber27 = first_witness(row, fiber27)
        if w_shell27 is None and w_shell41 is not None:
            out.append((r, w_shell27, w_shell41, w_fiber27))
    return out


def random_structured_rows(seed: int = 2443, trials: int = 600) -> list[tuple[str, tuple[int, ...]]]:
    rng = random.Random(seed)
    rows: list[tuple[str, tuple[int, ...]]] = []
    for trial in range(trials):
        core_size = rng.randint(7, 12)
        core = [7 * k for k in rng.sample(range(1, 20), core_size)]
        if not any(v % N == 0 for v in core):
            core[0] = 14 * rng.randint(1, 10)
        rest = rng.sample([v for v in range(1, 220) if v % 7], K - core_size)
        row = tuple(sorted(core + rest))
        if len(set(row)) == K and row_gcd(row) == 1:
            rows.append((f"rand{trial:03d}", row))
    return rows


def search_high_support() -> None:
    print("C. Random structured rows ranked by marked ladder pressure")
    q27 = denominators(27)
    q41 = denominators(41)
    scored = []
    for name, row in random_structured_rows():
        s27 = summarize_row(name, row, 27)
        s41 = summarize_row(name, row, 41)
        no_witness27 = s27.witness is None
        score = (
            int(no_witness27),
            s27.no_universal_count,
            s27.max_min_blockers,
            s27.blocked_count,
            -min((q for q, _a in [s41.witness] if q is not None), default=10**9),
        )
        if no_witness27 or s27.no_universal_count >= 25:
            scored.append((score, name, row, s27, s41))
    scored.sort(reverse=True)
    print(f"  candidates retained={len(scored)} out of {len(random_structured_rows())}")
    for _score, name, row, s27, s41 in scored[:8]:
        print(f"  {name}: row={row}")
        print_row_summary(s27)
        print_row_summary(s41)
    print()
    print(
        "  Reading: the hard rows are not rows with many clocks blocked by one "
        "apex multiple.  The interesting rows have many cover-blocked q's, where "
        "two or more marked runners share the cost.  Those are the rows where an "
        "owner-private deletion or Bprime(any runner) lemma has actual content."
    )
    print()


@dataclass(frozen=True)
class RouteVertex:
    name: str
    scores: tuple[int, ...]


def tournament_edges(vertices: tuple[RouteVertex, ...]) -> dict[str, set[str]]:
    edges = {v.name: set() for v in vertices}
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            av = sum(1 for x, y in zip(a.scores, b.scores) if x > y)
            bv = sum(1 for x, y in zip(a.scores, b.scores) if y > x)
            if av >= bv:
                edges[a.name].add(b.name)
            else:
                edges[b.name].add(a.name)
    return edges


def count_directed_triangles(vertices: tuple[RouteVertex, ...], edges: dict[str, set[str]]) -> int:
    names = [v.name for v in vertices]
    total = 0
    for a, b, c in combinations(names, 3):
        if (
            (b in edges[a] and c in edges[b] and a in edges[c])
            or (c in edges[a] and b in edges[c] and a in edges[b])
        ):
            total += 1
    return total


def scc_sizes(vertices: tuple[RouteVertex, ...], edges: dict[str, set[str]]) -> list[int]:
    names = [v.name for v in vertices]
    reverse = {name: set() for name in names}
    for a, outs in edges.items():
        for b in outs:
            reverse[b].add(a)

    def reach(start: str, graph: dict[str, set[str]]) -> set[str]:
        seen = {start}
        todo = deque([start])
        while todo:
            x = todo.popleft()
            for y in graph[x]:
                if y not in seen:
                    seen.add(y)
                    todo.append(y)
        return seen

    remaining = set(names)
    sizes = []
    while remaining:
        start = min(remaining)
        comp = reach(start, edges) & reach(start, reverse)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(vertices: tuple[RouteVertex, ...], edges: dict[str, set[str]]) -> int:
    names = tuple(v.name for v in vertices)
    total = 0

    def rec(path: tuple[str, ...], rest: tuple[str, ...]) -> None:
        nonlocal total
        if not rest:
            total += 1
            return
        last = path[-1]
        for i, nxt in enumerate(rest):
            if nxt in edges[last]:
                rec(path + (nxt,), rest[:i] + rest[i + 1 :])

    for i, start in enumerate(names):
        rec((start,), names[:i] + names[i + 1 :])
    return total


def print_tournament() -> None:
    print("D. Tournament Analysis over proof routes")
    vertices = (
        RouteVertex("marked_ladder_setcover", (5, 5, 5, 4, 5, 5)),
        RouteVertex("owner_private_deletion", (4, 5, 4, 5, 4, 5)),
        RouteVertex("Bprime_any_runner", (4, 4, 3, 5, 3, 5)),
        RouteVertex("apex_opening_modes", (5, 5, 3, 4, 4, 4)),
        RouteVertex("raw_Q_search", (3, 4, 5, 1, 3, 3)),
        RouteVertex("raw_scalar_floor", (2, 1, 2, 1, 5, 1)),
    )
    edges = tournament_edges(vertices)
    hist = Counter(len(edges[v.name]) for v in vertices)
    ordered = sorted(vertices, key=lambda v: (-len(edges[v.name]), v.name))
    print("  vertices=proof routes/support quotients, not runners")
    print("  observable=(marking, exactness, coverage, owner-content, quotient-stability, action)")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={count_directed_triangles(vertices, edges)}")
    print(f"  scc_sizes={scc_sizes(vertices, edges)}")
    print(f"  hamiltonian_paths={hamiltonian_paths(vertices, edges)}")
    print("  outdegree order:")
    for v in ordered:
        print(f"    {v.name:24s} out={len(edges[v.name])} scores={v.scores}")
    print()


def main() -> None:
    print("=" * 78)
    print("Codex LRC14 marked ladder support gate")
    print("=" * 78)
    print("HYP-2443 / T787")
    print(f"n={N}; delta={fmt_frac(DELTA)}; divisors={DIVISORS}")
    print()

    print("A. Baseline and known hard-shape rows")
    rows = [
        ("AP", AP),
        ("Vstar", VSTAR),
        ("2AP", TWOP),
        ("S(611)", family_row(611)),
        ("S(702)", family_row(702)),
        ("S(1053)", family_row(1053)),
    ]
    for name, row in rows:
        print(f"  row {name}: gcd={row_gcd(row)} mult14={any(v % N == 0 for v in row)}")
        print_row_summary(summarize_row(name, row, 27))
        print_row_summary(summarize_row(name, row, 41))
    print()

    print("B. Single-stranger evader scan 7*{1..12}+r")
    evaders = evader_scan()
    print(f"  pure-shell q<=27 blocked but q<=41 caught rows for r<1200: {len(evaders)}")
    print("  first 16:")
    for r, _w27, w41, wfiber in evaders[:16]:
        assert w41 is not None and wfiber is not None
        print(
            f"    r={r:4d}  shell41={w41[1]}/{w41[0]}  "
            f"fiber27={wfiber[1]}/{wfiber[0]}  r_mod27={r%27:2d} r_mod13={r%13:2d}"
        )
    print("  distinguished HYP-2438 evaders:")
    for r in EVADER_RS:
        row = family_row(r)
        w_shell41 = first_witness(row, shell_denominators(41))
        w_fiber27 = first_witness(row, denominators(27))
        print(
            f"    r={r:4d}  shell41={w_shell41[1]}/{w_shell41[0] if w_shell41 else 'none'}  "
            f"fiber27={w_fiber27[1]}/{w_fiber27[0] if w_fiber27 else 'none'}  r_mod27={r%27:2d}"
        )
    print(
        "  Reading: the S(r) rows are failures of the pure band-1 shell, not of the "
        "fibered m<=27 ladder.  The fibered address q=91 already sees the 7-core "
        "as the proven n=13 problem; shell q=40/41 is the non-fibered rung-up rescue."
    )
    print()

    search_high_support()
    print_tournament()

    print("E. Proof target sharpened")
    print(
        "  HYP-2438 can be split into two marked-support lemmas:\n"
        "    L1 finite ladder: if Q41 has a witness, the band-ladder closes directly.\n"
        "    L2 concentrated support: if Q41 is blocked, the marked set-cover ledger has\n"
        "       either a universal blocker or repeated high-load cover runners; that marked\n"
        "       runner should feed Bprime(any runner) or an owner-private deletion opening.\n"
        "  This is a proof obligation, not yet a proof.  The new invariant is the minimum\n"
        "  blocker support per q plus cross-q cover load; it is the address coordinate\n"
        "  missing from raw Res_27/Q search."
    )


if __name__ == "__main__":
    main()
