#!/usr/bin/env python3
"""
lrc_n14_dihedral_spine_pivots_s576.py

codex-2026-06-03 S576

Use the recent polygon/dihedral and obligation-hypergraph ideas on the n=14
proof route.

For n=14, C=2n-1=27.  The unit-shell gate U_a has exactly nine obligations:
  a in {1,2,4,5,7,8,10,11,13}.

Thus any strict sub-edge candidate must spend at least nine of its thirteen
runners hitting unit antipodal shells.  The remaining four runners are the
"mesh slack" layer where the composite nonunit holes, endpoint blockers, and
small-denominator shields must live.

This script audits the canonical polygon spine
  {1,2,4,5,7,8,10,11,13}
plus four nonunit/zero-residue slack runners (multiples of 3) up to a bound.
It records:
  * the dihedral quotient D/U/N obligation cover,
  * private quotient obligations/pivots,
  * exact maximin M(S),
  * whether any row falls below the LRC floor 1/14.

Tournament Analysis:
  Vertices are slack rows over the fixed unit spine.
  Pair observable is (maximin, full-cover flag, number of private quotient
  obligations, slack shell multiplicity, maximum speed).
  Switch/gauge: harder = smaller maximin, then full cover, then fewer private
  pivots, then smaller max speed.  Tie Hamiltonian path is lexicographic slack.

Assumption challenge:
  Vertices could be runners, time words, dihedral orbits, D/U/N obligations,
  unit shells, slack residues, or endpoint blockers.  This script uses slack
  rows over the mandatory unit-shell polygon spine.  It preserves the n=14
  unit-shell proof obligation and exact maximin, but destroys arbitrary choices
  of large unit-shell representatives.  So it is a spine normal-form probe,
  not a full n=14 proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


N = 14
K = N - 1
C = 2 * N - 1
FLOOR = Fraction(1, N)
EDGE = Fraction(2, C)
UNIT_SPINE = (1, 2, 4, 5, 7, 8, 10, 11, 13)
SLACK_BOUND = 42

Obligation = tuple[str, int]


@dataclass(frozen=True)
class RowReport:
    speeds: tuple[int, ...]
    slack: tuple[int, ...]
    full_cover: bool
    uncovered: tuple[Obligation, ...]
    private_count: int
    private_layer_hist: tuple[tuple[str, int], ...]
    maximin: Fraction | None
    witness: Fraction | None
    active: tuple[int, ...]
    slack_shells: tuple[int, ...]


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def score_time(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm(Fraction(v) * t) for v in speeds)


def exact_maximin(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    candidates: set[Fraction] = set()
    for i, a in enumerate(speeds):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                candidates.add(t)
        for b in speeds[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    candidates.add(Fraction(m, den))

    best = Fraction(0)
    best_t = Fraction(0)
    for t in candidates:
        score = score_time(speeds, t)
        if (score, -t) > (best, -best_t):
            best = score
            best_t = t
    return best, best_t


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "-"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def active_runners(speeds: tuple[int, ...], t: Fraction | None) -> tuple[int, ...]:
    if t is None:
        return tuple()
    score = score_time(speeds, t)
    return tuple(v for v in speeds if norm(Fraction(v) * t) == score)


def unit_shell(a: int) -> int:
    r = a % C
    if r == 0:
        return 0
    return min(r, C - r)


def obligation_universe_quotient() -> list[Obligation]:
    obligations: list[Obligation] = []
    obligations.extend(("D", q) for q in range(2, N))
    obligations.extend(("U", a) for a in UNIT_SPINE)
    obligations.extend(("N", j) for j in range(1, N // 2 + 1))
    return obligations


def covers(v: int, obligation: Obligation) -> bool:
    layer, label = obligation
    if layer == "D":
        return v % label == 0
    if layer == "U":
        return unit_shell(v) == label and gcd(label, C) == 1
    if layer == "N":
        # N_j and N_{14-j} are reflection-paired and have identical blockers.
        r = (v * label) % N
        return min(r, N - r) <= 1
    raise ValueError(obligation)


def incidence(speeds: tuple[int, ...]) -> dict[Obligation, tuple[int, ...]]:
    out: dict[Obligation, tuple[int, ...]] = {}
    for obligation in obligation_universe_quotient():
        out[obligation] = tuple(v for v in speeds if covers(v, obligation))
    return out


def layer_hist(obligations: tuple[Obligation, ...] | list[Obligation]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted(Counter(layer for layer, _ in obligations).items()))


def slack_shells(slack: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(unit_shell(v) for v in slack))


def analyze(slack: tuple[int, ...]) -> RowReport:
    speeds = tuple(sorted(UNIT_SPINE + slack))
    inc = incidence(speeds)
    uncovered = tuple(o for o, owners in inc.items() if not owners)
    private = tuple(o for o, owners in inc.items() if len(owners) == 1)
    full = not uncovered
    maximin = witness = None
    active: tuple[int, ...] = tuple()
    if full:
        maximin, witness = exact_maximin(speeds)
        active = active_runners(speeds, witness)
    return RowReport(
        speeds=speeds,
        slack=slack,
        full_cover=full,
        uncovered=uncovered,
        private_count=len(private),
        private_layer_hist=layer_hist(private),
        maximin=maximin,
        witness=witness,
        active=active,
        slack_shells=slack_shells(slack),
    )


def tournament_fingerprint(reports: list[RowReport]) -> dict[str, object]:
    # Summarize only full covers; non-covers have immediate clock witnesses.
    rows = [r for r in reports if r.full_cover]
    sample = sorted(rows, key=lambda r: (r.maximin or Fraction(9), r.private_count, max(r.speeds), r.slack))[:24]

    def key(row: RowReport) -> tuple[Fraction, int, int, tuple[int, ...]]:
        assert row.maximin is not None
        return (-row.maximin, -row.private_count, -max(row.speeds), tuple(-v for v in row.slack))

    n = len(sample)
    adj = [[False] * n for _ in range(n)]
    for i, left in enumerate(sample):
        for j, right in enumerate(sample):
            if i == j:
                continue
            adj[i][j] = key(left) > key(right) or (key(left) == key(right) and i < j)

    scores = [sum(row) for row in adj]
    cycles = 0
    for i, j, k in combinations(range(n), 3):
        cycles += int(
            (adj[i][j] and adj[j][k] and adj[k][i])
            or (adj[i][k] and adj[k][j] and adj[j][i])
        )

    def reach(start: int) -> set[int]:
        seen = {start}
        todo = deque([start])
        while todo:
            u = todo.popleft()
            for v, edge in enumerate(adj[u]):
                if edge and v not in seen:
                    seen.add(v)
                    todo.append(v)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        u = next(iter(remaining))
        ru = reach(u)
        comp = {v for v in remaining if v in ru and u in reach(v)}
        sccs.append(len(comp))
        remaining -= comp

    return {
        "sampled_vertices": [r.slack for r in sample],
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": cycles,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": 1 if cycles == 0 and len(set(scores)) == n else "not_counted",
    }


def print_report(label: str, row: RowReport) -> None:
    print(f"[{label}]")
    print(f"  speeds={row.speeds}")
    print(f"  slack={row.slack} slack_shells={row.slack_shells}")
    print(f"  full_DUN_quotient_cover={row.full_cover} uncovered={row.uncovered or '-'}")
    print(f"  private_quotient_obligations={row.private_count} {dict(row.private_layer_hist)}")
    print(f"  M={fmt_frac(row.maximin)} witness={fmt_frac(row.witness)} active={row.active or '-'}")
    print()


def main() -> None:
    print("S576 n=14 dihedral unit-spine private-pivot audit")
    print("=" * 78)
    print(f"n={N}; C={C}; floor={fmt_frac(FLOOR)}; edge={fmt_frac(EDGE)}")
    print(f"mandatory unit-shell spine={UNIT_SPINE}")
    print(f"dihedral quotient obligations={len(obligation_universe_quotient())}")
    print("  D: 12 fixed denominator gates; U: 9 unit shells; N: 7 reflected n-clock gates")
    print()

    named = {
        "AP floor": (3, 6, 9, 12),
        "V* floor": (3, 6, 9, 24),
        "zero-slack variant": (3, 6, 9, 27),
        "double-3 variant": (3, 6, 24, 30),
    }
    for label, slack in named.items():
        print_report(label, analyze(tuple(sorted(slack))))

    slack_candidates = tuple(v for v in range(3, SLACK_BOUND + 1) if v % 3 == 0 and v not in UNIT_SPINE)
    reports: list[RowReport] = []
    uncovered_hist: Counter[tuple[Obligation, ...]] = Counter()
    for slack in combinations(slack_candidates, 4):
        row = analyze(slack)
        reports.append(row)
        if not row.full_cover:
            uncovered_hist[row.uncovered] += 1

    full = [r for r in reports if r.full_cover]
    below_floor = [r for r in full if r.maximin is not None and r.maximin < FLOOR]
    floor_rows = [r for r in full if r.maximin == FLOOR]
    open_gap = [r for r in full if r.maximin is not None and FLOOR < r.maximin < EDGE]
    edge_or_above = [r for r in full if r.maximin is not None and r.maximin >= EDGE]
    min_m = min((r.maximin for r in full if r.maximin is not None), default=None)
    min_rows = [r for r in full if r.maximin == min_m]

    print("Canonical unit-spine slack scan")
    print("-" * 78)
    print(f"  slack_bound={SLACK_BOUND}; slack_candidates={len(slack_candidates)}; rows={len(reports)}")
    print(f"  full_DUN_quotient_covers={len(full)}")
    print(f"  below_floor={len(below_floor)}")
    print(f"  floor_rows={len(floor_rows)} open_gap_rows={len(open_gap)} edge_or_above={len(edge_or_above)}")
    print(f"  min_M={fmt_frac(min_m)} min_count={len(min_rows)}")
    print(f"  private_count_hist={dict(sorted(Counter(r.private_count for r in full).items()))}")
    print(f"  floor_slack_shell_hist={dict(Counter(r.slack_shells for r in floor_rows).most_common(8))}")
    print("  common_uncovered_for_non_covers:")
    for obs, count in uncovered_hist.most_common(6):
        preview = ",".join(f"{a}{b}" for a, b in obs[:8])
        more = "" if len(obs) <= 8 else f",...(+{len(obs)-8})"
        print(f"    {count:5d}: {preview}{more}")
    print("  min_rows_sample:")
    for row in min_rows[:8]:
        print(
            f"    slack={row.slack} shells={row.slack_shells} "
            f"M={fmt_frac(row.maximin)} t={fmt_frac(row.witness)} active={row.active}"
        )
    print("  open_gap_sample:")
    for row in open_gap[:8]:
        print(
            f"    slack={row.slack} shells={row.slack_shells} "
            f"M={fmt_frac(row.maximin)} t={fmt_frac(row.witness)} active={row.active}"
        )
    print()

    print("Tournament Analysis")
    print("-" * 78)
    print("  vertices: full-cover slack rows over the fixed unit spine")
    print("  observable: M, private quotient pivots, max speed, slack tuple")
    print("  switch: harder = smaller M, fewer pivots, smaller max speed")
    print(f"  fingerprints: {tournament_fingerprint(reports)}")
    print()

    print("Synthesis")
    print("-" * 78)
    print("The n=14 unit-shell gate forces a nine-runner polygon spine before any")
    print("counterexample can enter the strict sub-edge residual.  In the canonical")
    print(f"spine scan through slack <={SLACK_BOUND}, every full D/U/N quotient cover has")
    print("M>=1/14; the minimum rows are exactly floor rows.")
    print("This does not prove n=14: large representatives of the nine unit shells")
    print("are not scanned.  It does isolate a proof subproblem:")
    print("  normalize the unit-shell spine or prove an exchange that lowers it while")
    print("  preserving private U pivots; then the four nonunit slack runners carry")
    print("  all remaining D/N, endpoint-blocker, and composite-hole debt.")


if __name__ == "__main__":
    main()
