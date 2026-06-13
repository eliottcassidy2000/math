#!/usr/bin/env python3
"""
lrc_n14_unit_spine_exchange_s578.py

codex-2026-06-03 S578

Push on the HYP-2096 unit-spine exchange lemma.

Question:
  If n=14 and a row has one runner in each of the nine unit shells mod 27,
  can large unit-shell representatives be lowered to the canonical spine
  (1,2,4,5,7,8,10,11,13), unless a cheap HYP-2095 unblocked small-pair witness
  already proves loneliness?

This script audits two bounded exchange regimes:

  1. one-unit-lift exhaustive:
       all canonical slack quadruples through 42, and one unit-shell
       representative replaced by any same-shell representative through 81.

  2. named local/extreme unit stress:
       all one- and two-unit lifts through 81, plus several all-shell extreme
       patterns, over selected slack layers AP, V*, first open-gap, and two
       non-cover controls.

For every full D/U/N quotient cover it asks:
  * does HYP-2095 already give an unblocked small pair?
  * if not, does lowering all unit representatives to the canonical spine
    preserve the D/U/N quotient cover?
  * if not, what quotient obligations are lost?
  * among the no-cheap residual, does exact M ever fall below 1/14?

Tournament Analysis:
  Vertices are exchange patterns, not runners: one-lift shell labels and named
  local/extreme stress families.  Pair observable is
    (bad lowering count, below-floor count, no-cheap residual count,
     full-cover count, max representative).
  Switch/gauge: harder = more bad lowerings, then more below-floor rows, then
  more no-cheap residual.  Tie Hamiltonian path is lexicographic by label.

Assumption challenge:
  Vertices could be runners, small pairs, unit shells, D/U/N obligations,
  fixed round classes, or exchange moves.  This script uses exchange patterns.
  It preserves the predicate relevant to the unit-spine exchange lemma:
  whether large unit representatives are needed after cheap pair witnesses are
  removed.  It destroys the full 64 fixed-round class identity, so it is still
  a local exchange audit, not the S578 fixed-class certificate table.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd


N = 14
K = N - 1
C = 2 * N - 1
FLOOR = Fraction(1, N)
EDGE = Fraction(2, C)
UNIT_SPINE = (1, 2, 4, 5, 7, 8, 10, 11, 13)
ONE_LIFT_UNIT_BOUND = 81
ONE_LIFT_SLACK_BOUND = 42
LOCAL_UNIT_BOUND = 81

Obligation = tuple[str, int]


@dataclass(frozen=True)
class Witness:
    pair: tuple[int, int]
    denominator: int
    numerator: int

    @property
    def time(self) -> Fraction:
        return Fraction(self.numerator, self.denominator)


@dataclass(frozen=True)
class AuditSummary:
    label: str
    rows: int
    full_cover: int
    cheap_witness: int
    no_cheap: int
    lowers_to_full: int
    bad_lowering: int
    below_floor: int
    floor_or_above: int
    min_m: Fraction | None
    max_rep: int
    witness_den_hist: tuple[tuple[int, int], ...]
    witness_examples: tuple[str, ...]
    lost_hist: tuple[tuple[tuple[Obligation, ...], int], ...]
    examples: tuple[str, ...]

    @property
    def burden(self) -> tuple[int, int, int, int, int]:
        return (self.bad_lowering, self.below_floor, self.no_cheap, self.full_cover, self.max_rep)


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "-"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


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


def unit_shell(v: int) -> int:
    r = v % C
    if r == 0:
        return 0
    return min(r, C - r)


def unit_reps(shell: int, bound: int) -> tuple[int, ...]:
    return tuple(v for v in range(1, bound + 1) if unit_shell(v) == shell and gcd(v, C) == 1)


def obligation_universe_quotient() -> tuple[Obligation, ...]:
    return (
        tuple(("D", q) for q in range(2, N))
        + tuple(("U", a) for a in UNIT_SPINE)
        + tuple(("N", j) for j in range(1, N // 2 + 1))
    )


OBLIGATIONS = obligation_universe_quotient()


def covers(v: int, obligation: Obligation) -> bool:
    layer, label = obligation
    if layer == "D":
        return v % label == 0
    if layer == "U":
        return unit_shell(v) == label and gcd(label, C) == 1
    if layer == "N":
        r = (v * label) % N
        return min(r, N - r) <= 1
    raise ValueError(obligation)


def uncovered_obligations(speeds: tuple[int, ...]) -> tuple[Obligation, ...]:
    return tuple(ob for ob in OBLIGATIONS if not any(covers(v, ob) for v in speeds))


def is_full_cover(speeds: tuple[int, ...]) -> bool:
    return not uncovered_obligations(speeds)


def canonical_row(slack: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(UNIT_SPINE + slack))


def reduced_sum(a: int, b: int) -> int:
    return (a + b) // gcd(a, b)


@lru_cache(maxsize=None)
def pair_safe_numerators(a: int, b: int) -> tuple[int, ...]:
    den = a + b
    return tuple(
        m
        for m in range(1, den)
        if norm(Fraction(a * m, den)) >= FLOOR and norm(Fraction(b * m, den)) >= FLOOR
    )


def unblocked_small_pair(speeds: tuple[int, ...]) -> Witness | None:
    row = tuple(sorted(speeds))
    for a, b in combinations(row, 2):
        if reduced_sum(a, b) > N:
            continue
        den = a + b
        for m in pair_safe_numerators(a, b):
            t = Fraction(m, den)
            if all(norm(Fraction(c) * t) >= FLOOR for c in row):
                return Witness((a, b), den, m)
    return None


def active_runners(speeds: tuple[int, ...], t: Fraction) -> tuple[int, ...]:
    score = score_time(speeds, t)
    return tuple(v for v in speeds if norm(Fraction(v) * t) == score)


def format_obligations(obs: tuple[Obligation, ...]) -> str:
    if not obs:
        return "-"
    return ",".join(f"{layer}{label}" for layer, label in obs)


def exact_for_residual(
    speeds: tuple[int, ...],
    min_seen: Fraction | None,
) -> tuple[Fraction, Fraction, tuple[int, ...], Fraction]:
    m, t = exact_maximin(speeds)
    active = active_runners(speeds, t)
    new_min = m if min_seen is None else min(min_seen, m)
    return m, t, active, new_min


def summarize(
    label: str,
    rows: int,
    full_cover: int,
    cheap: int,
    no_cheap: int,
    lowers_to_full: int,
    bad_lowering: int,
    below_floor: int,
    floor_or_above: int,
    min_m: Fraction | None,
    max_rep: int,
    witness_den_hist: Counter[int],
    witness_examples: list[str],
    lost_hist: Counter[tuple[Obligation, ...]],
    examples: list[str],
) -> AuditSummary:
    return AuditSummary(
        label=label,
        rows=rows,
        full_cover=full_cover,
        cheap_witness=cheap,
        no_cheap=no_cheap,
        lowers_to_full=lowers_to_full,
        bad_lowering=bad_lowering,
        below_floor=below_floor,
        floor_or_above=floor_or_above,
        min_m=min_m,
        max_rep=max_rep,
        witness_den_hist=tuple(witness_den_hist.most_common(8)),
        witness_examples=tuple(witness_examples[:6]),
        lost_hist=tuple(lost_hist.most_common(6)),
        examples=tuple(examples[:6]),
    )


def audit_one_lift() -> AuditSummary:
    slack_candidates = tuple(
        v for v in range(3, ONE_LIFT_SLACK_BOUND + 1) if v % 3 == 0 and v not in UNIT_SPINE
    )
    reps_by_shell = {
        shell: tuple(v for v in unit_reps(shell, ONE_LIFT_UNIT_BOUND) if v != shell)
        for shell in UNIT_SPINE
    }

    canonical_full_cache: dict[tuple[int, ...], bool] = {}
    rows = full_cover = cheap = no_cheap = 0
    lowers_to_full = bad_lowering = below_floor = floor_or_above = 0
    min_m: Fraction | None = None
    max_rep = 0
    witness_den_hist: Counter[int] = Counter()
    witness_examples: list[str] = []
    lost_hist: Counter[tuple[Obligation, ...]] = Counter()
    examples: list[str] = []

    for slack in combinations(slack_candidates, 4):
        slack = tuple(sorted(slack))
        lowered = canonical_row(slack)
        lowered_uncovered = uncovered_obligations(lowered)
        canonical_full_cache[slack] = not lowered_uncovered
        for i, shell in enumerate(UNIT_SPINE):
            for rep in reps_by_shell[shell]:
                units = list(UNIT_SPINE)
                units[i] = rep
                speeds = tuple(sorted(tuple(units) + slack))
                rows += 1
                max_rep = max(max_rep, rep)
                if not is_full_cover(speeds):
                    continue
                full_cover += 1
                witness = unblocked_small_pair(speeds)
                if witness is not None:
                    cheap += 1
                    witness_den_hist[witness.denominator] += 1
                    if len(witness_examples) < 6:
                        witness_examples.append(
                            f"slack={slack} lift {shell}->{rep} "
                            f"pair={witness.pair} t={fmt_frac(witness.time)}"
                        )
                    continue
                no_cheap += 1
                m, t, active, min_m = exact_for_residual(speeds, min_m)
                if m < FLOOR:
                    below_floor += 1
                else:
                    floor_or_above += 1
                if canonical_full_cache[slack]:
                    lowers_to_full += 1
                else:
                    bad_lowering += 1
                    lost_hist[lowered_uncovered] += 1
                    if len(examples) < 6:
                        examples.append(
                            "slack="
                            f"{slack} lift {shell}->{rep} M={fmt_frac(m)} t={fmt_frac(t)} "
                            f"active={active} lowered_loses={format_obligations(lowered_uncovered)}"
                        )

    return summarize(
        "one_unit_lift_all_slack",
        rows,
        full_cover,
        cheap,
        no_cheap,
        lowers_to_full,
        bad_lowering,
        below_floor,
        floor_or_above,
        min_m,
        max_rep,
        witness_den_hist,
        witness_examples,
        lost_hist,
        examples,
    )


def local_unit_patterns() -> tuple[tuple[int, ...], ...]:
    reps_by_shell = [
        tuple(v for v in unit_reps(shell, LOCAL_UNIT_BOUND) if v != shell) for shell in UNIT_SPINE
    ]
    patterns: set[tuple[int, ...]] = {UNIT_SPINE}

    for i, reps in enumerate(reps_by_shell):
        for rep in reps:
            units = list(UNIT_SPINE)
            units[i] = rep
            patterns.add(tuple(units))

    for i, j in combinations(range(len(UNIT_SPINE)), 2):
        for rep_i in reps_by_shell[i]:
            for rep_j in reps_by_shell[j]:
                units = list(UNIT_SPINE)
                units[i] = rep_i
                units[j] = rep_j
                patterns.add(tuple(units))

    patterns.add(tuple(C - shell for shell in UNIT_SPINE))
    patterns.add(tuple(shell + C for shell in UNIT_SPINE))
    patterns.add(tuple((C - shell) + C for shell in UNIT_SPINE))
    patterns.add(tuple(max(unit_reps(shell, LOCAL_UNIT_BOUND)) for shell in UNIT_SPINE))
    return tuple(sorted(patterns))


def audit_named_local_unit_patterns() -> list[AuditSummary]:
    named_slacks = {
        "AP_slack": (3, 6, 9, 12),
        "Vstar_slack": (3, 6, 9, 24),
        "first_open_gap_slack": (3, 6, 9, 36),
        "zero_slack_control": (3, 6, 9, 27),
        "double3_control": (3, 6, 24, 30),
    }
    patterns = local_unit_patterns()
    out: list[AuditSummary] = []

    for label, slack_raw in named_slacks.items():
        slack = tuple(sorted(slack_raw))
        lowered = canonical_row(slack)
        lowered_uncovered = uncovered_obligations(lowered)
        lowered_full = not lowered_uncovered
        rows = full_cover = cheap = no_cheap = 0
        lowers_to_full = bad_lowering = below_floor = floor_or_above = 0
        min_m: Fraction | None = None
        max_rep = 0
        witness_den_hist: Counter[int] = Counter()
        witness_examples: list[str] = []
        lost_hist: Counter[tuple[Obligation, ...]] = Counter()
        examples: list[str] = []

        for units in patterns:
            speeds = tuple(sorted(units + slack))
            rows += 1
            max_rep = max(max_rep, max(units))
            if not is_full_cover(speeds):
                continue
            full_cover += 1
            witness = unblocked_small_pair(speeds)
            if witness is not None:
                cheap += 1
                witness_den_hist[witness.denominator] += 1
                if len(witness_examples) < 6:
                    witness_examples.append(
                        f"units={units} pair={witness.pair} t={fmt_frac(witness.time)}"
                    )
                continue
            no_cheap += 1
            m, t, active, min_m = exact_for_residual(speeds, min_m)
            if m < FLOOR:
                below_floor += 1
            else:
                floor_or_above += 1
            if lowered_full:
                lowers_to_full += 1
            else:
                bad_lowering += 1
                lost_hist[lowered_uncovered] += 1
                if len(examples) < 6:
                    examples.append(
                        f"units={units} M={fmt_frac(m)} t={fmt_frac(t)} "
                        f"active={active} lowered_loses={format_obligations(lowered_uncovered)}"
                    )

        out.append(
            summarize(
                label,
                rows,
                full_cover,
                cheap,
                no_cheap,
                lowers_to_full,
                bad_lowering,
                below_floor,
                floor_or_above,
                min_m,
                max_rep,
                witness_den_hist,
                witness_examples,
                lost_hist,
                examples,
            )
        )
    return out


def tournament_fingerprint(rows: list[AuditSummary]) -> dict[str, object]:
    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    for i, left in enumerate(rows):
        for j, right in enumerate(rows):
            if i == j:
                continue
            adj[i][j] = (left.burden, left.label) > (right.burden, right.label)

    scores = [sum(row) for row in adj]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        ):
            c3 += 1

    def reach(start: int) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for w, edge in enumerate(adj[v]):
                if edge and w not in seen:
                    seen.add(w)
                    q.append(w)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        start = next(iter(remaining))
        forward = reach(start)
        comp = {v for v in remaining if v in forward and start in reach(v)}
        sccs.append(len(comp))
        remaining -= comp

    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hardness_order": [r.label for r in sorted(rows, key=lambda x: (x.burden, x.label), reverse=True)],
    }


def print_summary(summary: AuditSummary) -> None:
    print(f"[{summary.label}]")
    print(f"  rows={summary.rows} full_cover={summary.full_cover}")
    print(
        f"  cheap_unblocked_pair={summary.cheap_witness} no_cheap_residual={summary.no_cheap}"
    )
    print(f"  cheap_witness_den_hist={dict(summary.witness_den_hist)}")
    print("  cheap_witness_examples:")
    if summary.witness_examples:
        for example in summary.witness_examples:
            print(f"    {example}")
    else:
        print("    -")
    print(
        f"  lowers_to_full={summary.lowers_to_full} bad_lowering={summary.bad_lowering}"
    )
    print(
        f"  below_floor={summary.below_floor} floor_or_above={summary.floor_or_above} "
        f"min_M={fmt_frac(summary.min_m)} max_rep={summary.max_rep}"
    )
    print("  lost_obligation_hist:")
    if summary.lost_hist:
        for obs, count in summary.lost_hist:
            print(f"    {count:6d}: {format_obligations(obs)}")
    else:
        print("    -")
    print("  examples:")
    if summary.examples:
        for example in summary.examples:
            print(f"    {example}")
    else:
        print("    -")
    print()


def main() -> None:
    print("S578 n=14 unit-spine exchange audit")
    print("=" * 78)
    print(f"n={N}; C={C}; floor={fmt_frac(FLOOR)}; edge={fmt_frac(EDGE)}")
    print(f"canonical unit spine={UNIT_SPINE}")
    print(
        "unit reps per shell through 81="
        f"{[len(unit_reps(shell, LOCAL_UNIT_BOUND)) for shell in UNIT_SPINE]}"
    )
    print(
        "one-lift reps per shell through 81 excluding canonical="
        f"{[len(tuple(v for v in unit_reps(shell, ONE_LIFT_UNIT_BOUND) if v != shell)) for shell in UNIT_SPINE]}"
    )
    print(f"named local/extreme unit patterns={len(local_unit_patterns())}")
    print()

    summaries = [audit_one_lift()] + audit_named_local_unit_patterns()
    for summary in summaries:
        print_summary(summary)

    print("Tournament Analysis")
    print("-" * 78)
    print("  vertices: exchange audit patterns")
    print("  observable: bad lowerings, below-floor rows, no-cheap residual, full covers")
    print("  switch: harder = more bad lowerings, then below-floor, then no-cheap")
    print(f"  fingerprints: {tournament_fingerprint(summaries)}")
    print()

    print("Synthesis")
    print("-" * 78)
    print("Bounded evidence supports the exchange lemma in a sharper form:")
    print("  after HYP-2095 removes rows with an unblocked small pair, the audited")
    print("  full D/U/N unit-lift residual is empty.")
    print("In the one-lift exhaustive scan, every full cover already has a cheap")
    print("unblocked-pair witness.  In the named local/extreme unit stress,")
    print("AP/V*/open-gap families again have no bad lowerings or below-floor")
    print("residuals.")
    print("The two non-cover controls show why D9/D12-style composite debt cannot be")
    print("repaired by unit representatives: no full-cover residual survives there.")
    print("Next proof move: formalize this as a local exchange rule: a lifted unit")
    print("representative that is not part of an unblocked small-pair witness can be")
    print("lowered shellwise, because unit representatives do not repair the gcd-3")
    print("D-obligations that make the canonical slack layer fail.")


if __name__ == "__main__":
    main()
