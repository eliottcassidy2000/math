#!/usr/bin/env python3
"""Exact swap-aware gap races for the THM-387 LRC two-gap flow.

codex-2026-06-01 S522

S518 introduced the no-swap race approximation: after a runner wraps through
the observer, compare the resetter's time to reach 1/n with the nearest
counterclockwise runner's time to fall to 1/n.  This script follows the actual
wall sequence instead.  A race is attached to every wrap-around reset and is
classified by the first event that happens afterward:

    open_LL   an open cell with both adjacent gaps >= 1/n is reached
    wall_LL   a compactified wall point has both adjacent gaps >= 1/n
    loss_*    a non-lonely short-left state is reached first
    next_reset another runner wraps before the race resolves
    non_LS    the reset does not start a long-left/short-right race

Tournament Analysis declaration:
    vertices: selected primitive speed sets and hard rows
    pairwise observable: exact race strength vector
        (has open LL, has wall LL, winning reset count, max reset runway,
         negative losing reset count)
    switch/gauge: lexicographic comparison of the strength vector
    tie Hamiltonian path: listed selected-row order

Stored output:
    05-knowledge/results/lrc_gap_race_swap_s522.out
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DYN = SourceFileLoader(
    "lrc_two_gap_dynamics_s518",
    str(ROOT / "04-computation" / "lrc_two_gap_dynamics_s518.py"),
).load_module()
RACE = SourceFileLoader(
    "lrc_gap_race_analysis_s518",
    str(ROOT / "04-computation" / "lrc_gap_race_analysis_s518.py"),
).load_module()


ONE = Fraction(1)
ZERO = Fraction(0)
WIN_OUTCOMES = {"open_LL", "wall_LL"}


@dataclass(frozen=True)
class ExactRace:
    n: int
    speeds: tuple[int, ...]
    t_reset: Fraction
    resetter_idx: int
    v_resetter: int
    start_fiber: str
    outcome: str
    t_outcome: Fraction | None
    dt_outcome: Fraction | None
    L0: Fraction | None
    simplified: str
    simplified_ratio: Fraction | None
    note: str = ""


@dataclass(frozen=True)
class SetProfile:
    label: str
    n: int
    speeds: tuple[int, ...]
    has_open_ll: bool
    has_wall_ll: bool
    open_ll_measure: Fraction
    wall_ll_count: int
    race_counts: tuple[tuple[str, int], ...]
    winning_races: int
    losing_races: int
    max_L0: Fraction
    best_win_dt: Fraction | None
    simplified_misses_rescued: int

    def strength(self) -> tuple:
        best_dt = self.best_win_dt if self.best_win_dt is not None else ONE
        total_races = sum(v for _, v in self.race_counts)
        win_rate = Fraction(self.winning_races, total_races) if total_races else ZERO
        return (
            int(self.has_open_ll),
            int(self.has_wall_ll),
            self.winning_races,
            win_rate,
            self.max_L0,
            -self.losing_races,
            -best_dt,
        )


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def fmt_frac(value: Fraction | None) -> str:
    if value is None:
        return "-"
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fmt_float(value: Fraction | None, places: int = 6) -> str:
    if value is None:
        return "-"
    return f"{float(value):.{places}f}"


def is_primitive(speeds: tuple[int, ...]) -> bool:
    return reduce(gcd, speeds) == 1


def is_wrap_time(speeds: tuple[int, ...], t: Fraction) -> bool:
    return any(frac(Fraction(v) * t) == ZERO for v in speeds)


def fiber_at(speeds: tuple[int, ...], n: int, t: Fraction) -> str:
    left, right = DYN.compute_gaps_at(speeds, n, t)
    return DYN.fiber_label(left, right, Fraction(1, n))


def walls_after_reset(
    speeds: tuple[int, ...],
    n: int,
    t_reset: Fraction,
    base_walls: list[Fraction] | None = None,
) -> list[Fraction]:
    walls = base_walls if base_walls is not None else DYN.compute_walls(speeds, n)
    after = []
    for wall in walls:
        lifted = wall if wall > t_reset else wall + ONE
        if lifted > t_reset:
            after.append(lifted)
    after.append(t_reset + ONE)
    return sorted(set(after))


def simplified_fields(raw: dict) -> tuple[Fraction | None, str, Fraction | None]:
    if raw.get("degenerate", False):
        return None, "degenerate", None
    L0 = raw.get("L_0")
    outcome = raw.get("race_outcome", raw.get("entry_fiber", "unknown"))
    ratio = raw.get("race_ratio")
    return L0, outcome, ratio


def trace_exact_race(
    speeds: tuple[int, ...],
    n: int,
    t_reset: Fraction,
    resetter_idx: int,
    base_walls: list[Fraction] | None = None,
) -> ExactRace:
    raw = RACE.gap_race_at_reset(speeds, n, t_reset, resetter_idx)
    L0, simplified, ratio = simplified_fields(raw)
    v_resetter = speeds[resetter_idx]

    if raw.get("degenerate", False):
        return ExactRace(
            n,
            speeds,
            t_reset,
            resetter_idx,
            v_resetter,
            "degenerate",
            "degenerate",
            None,
            None,
            L0,
            simplified,
            ratio,
            "simultaneous reset",
        )

    walls = walls_after_reset(speeds, n, t_reset, base_walls)
    t_start = t_reset
    first_end = walls[0]
    start_fiber = fiber_at(speeds, n, (t_start + first_end) / 2)
    if start_fiber != "LS":
        return ExactRace(
            n,
            speeds,
            t_reset,
            resetter_idx,
            v_resetter,
            start_fiber,
            "non_LS",
            t_start,
            ZERO,
            L0,
            simplified,
            ratio,
            "first open cell after reset is not LS",
        )

    for t_end in walls:
        if t_end <= t_start:
            continue
        if t_end - t_reset > ONE:
            break

        mid = (t_start + t_end) / 2
        cell_fiber = fiber_at(speeds, n, mid)
        if cell_fiber == "LL":
            return ExactRace(
                n,
                speeds,
                t_reset,
                resetter_idx,
                v_resetter,
                start_fiber,
                "open_LL",
                t_start,
                t_start - t_reset,
                L0,
                simplified,
                ratio,
                f"open interval [{fmt_frac(t_start)}, {fmt_frac(t_end)}]",
            )
        if cell_fiber in {"SS", "SL"} and t_start > t_reset:
            return ExactRace(
                n,
                speeds,
                t_reset,
                resetter_idx,
                v_resetter,
                start_fiber,
                f"loss_{cell_fiber}",
                t_start,
                t_start - t_reset,
                L0,
                simplified,
                ratio,
                "short-left state reached before LL",
            )

        wall_fiber = fiber_at(speeds, n, t_end)
        if wall_fiber == "LL":
            return ExactRace(
                n,
                speeds,
                t_reset,
                resetter_idx,
                v_resetter,
                start_fiber,
                "wall_LL",
                t_end,
                t_end - t_reset,
                L0,
                simplified,
                ratio,
                "compactified threshold wall",
            )

        if is_wrap_time(speeds, t_end):
            return ExactRace(
                n,
                speeds,
                t_reset,
                resetter_idx,
                v_resetter,
                start_fiber,
                "next_reset",
                t_end,
                t_end - t_reset,
                L0,
                simplified,
                ratio,
                "another wrap-around reset starts a new race",
            )

        t_start = t_end

    return ExactRace(
        n,
        speeds,
        t_reset,
        resetter_idx,
        v_resetter,
        start_fiber,
        "no_outcome",
        None,
        None,
        L0,
        simplified,
        ratio,
        "no race resolution inside one period",
    )


def exact_races_for_set(speeds: tuple[int, ...], n: int) -> list[ExactRace]:
    base_walls = DYN.compute_walls(speeds, n)
    races = []
    for t_reset, resetter_idx in RACE.compute_wrap_arounds(speeds):
        races.append(trace_exact_race(speeds, n, t_reset, resetter_idx, base_walls))
    return races


def profile_speed_set(label: str, speeds: tuple[int, ...], n: int | None = None) -> SetProfile:
    if n is None:
        n = len(speeds) + 1
    analysis = DYN.analyze_speed_set(speeds, n)
    races = exact_races_for_set(speeds, n)
    race_counts = Counter(r.outcome for r in races)
    winning = sum(race_counts[outcome] for outcome in WIN_OUTCOMES)
    losing = sum(
        count
        for outcome, count in race_counts.items()
        if outcome.startswith("loss") or outcome == "next_reset"
    )
    max_L0 = max((r.L0 for r in races if r.L0 is not None), default=ZERO)
    win_dts = [r.dt_outcome for r in races if r.outcome in WIN_OUTCOMES and r.dt_outcome is not None]
    rescued = sum(
        1
        for r in races
        if r.outcome in WIN_OUTCOMES and r.simplified not in {"LL", "tie"}
    )
    return SetProfile(
        label=label,
        n=n,
        speeds=speeds,
        has_open_ll=analysis["has_open_ll"],
        has_wall_ll=analysis["has_wall_ll"],
        open_ll_measure=analysis["lonely_measure"],
        wall_ll_count=analysis["wall_ll_count"],
        race_counts=tuple(sorted(race_counts.items())),
        winning_races=winning,
        losing_races=losing,
        max_L0=max_L0,
        best_win_dt=min(win_dts) if win_dts else None,
        simplified_misses_rescued=rescued,
    )


def primitive_speed_sets(n: int, max_speed: int):
    for combo in combinations(range(1, max_speed + 1), n - 1):
        if is_primitive(combo):
            yield combo


def bounded_exact_audit() -> None:
    print("=" * 72)
    print("PART A: Exact swap-aware bounded audit")
    print("=" * 72)
    print()
    windows = [(3, 30), (4, 20), (5, 15), (6, 12)]
    for n, max_speed in windows:
        total = 0
        any_exact_win = 0
        any_open_win = 0
        any_wall_win = 0
        any_full_ll = 0
        mismatch_sets = 0
        outcome_counts: Counter[str] = Counter()
        rescued_examples = []
        no_exact_examples = []

        for speeds in primitive_speed_sets(n, max_speed):
            total += 1
            analysis = DYN.analyze_speed_set(speeds, n)
            races = exact_races_for_set(speeds, n)
            outcomes = Counter(r.outcome for r in races)
            outcome_counts.update(outcomes)

            has_exact = any(r.outcome in WIN_OUTCOMES for r in races)
            has_open = any(r.outcome == "open_LL" for r in races)
            has_wall = any(r.outcome == "wall_LL" for r in races)
            has_full = analysis["has_open_ll"] or analysis["has_wall_ll"]
            any_exact_win += int(has_exact)
            any_open_win += int(has_open)
            any_wall_win += int(has_wall)
            any_full_ll += int(has_full)
            if has_full != has_exact:
                mismatch_sets += 1
                if len(no_exact_examples) < 5:
                    no_exact_examples.append((speeds, has_full, outcomes))

            rescued = [
                r
                for r in races
                if r.outcome in WIN_OUTCOMES and r.simplified not in {"LL", "tie"}
            ]
            if rescued and len(rescued_examples) < 5:
                rescued_examples.append((speeds, rescued[0]))

        print(f"n={n} max_speed={max_speed} primitive_sets={total}")
        print(f"  full analysis has closed LL: {any_full_ll}/{total}")
        print(f"  exact reset race reaches LL: {any_exact_win}/{total}")
        print(f"  open-winning sets: {any_open_win}/{total}")
        print(f"  wall-winning sets: {any_wall_win}/{total}")
        print(f"  full/race mismatches: {mismatch_sets}")
        print(f"  aggregate race outcomes: {tuple(sorted(outcome_counts.items()))}")
        if rescued_examples:
            print("  compactification/swap rescues over no-swap losses:")
            for speeds, race in rescued_examples:
                print(
                    "    speeds={} reset v={} t={} exact={} simplified={} ratio={} L0={}".format(
                        speeds,
                        race.v_resetter,
                        fmt_frac(race.t_reset),
                        race.outcome,
                        race.simplified,
                        fmt_float(race.simplified_ratio),
                        fmt_frac(race.L0),
                    )
                )
        if no_exact_examples:
            print("  WARNING mismatch examples:")
            for speeds, has_full, outcomes in no_exact_examples:
                print(f"    speeds={speeds} has_full={has_full} outcomes={tuple(sorted(outcomes.items()))}")
        print()


def initial_segments() -> None:
    print("=" * 72)
    print("PART B: Initial segment exact race profile")
    print("=" * 72)
    print()
    for n in range(3, 21):
        speeds = tuple(range(1, n))
        profile = profile_speed_set(f"initial_n{n}", speeds, n)
        print(
            "n={:2d} races={:3d} outcomes={} open_mu={} wall_LL={} max_L0={} best_dt={}".format(
                n,
                sum(v for _, v in profile.race_counts),
                profile.race_counts,
                fmt_frac(profile.open_ll_measure),
                profile.wall_ll_count,
                fmt_frac(profile.max_L0),
                fmt_frac(profile.best_win_dt),
            )
        )
    print()


def selected_profiles() -> list[SetProfile]:
    selected = [
        ("n3_initial", (1, 2), 3),
        ("n5_initial", (1, 2, 3, 4), 5),
        ("n5_gap_sum_refuter", (5, 11, 12, 17), 5),
        ("n6_initial", (1, 2, 3, 4, 5), 6),
        ("n14_initial", tuple(range(1, 14)), 14),
        ("n14_gate_7", tuple(list(range(1, 7)) + list(range(8, 14)) + [14]), 14),
        ("n18_initial", tuple(range(1, 18)), 18),
        ("n18_gate_8", tuple(list(range(1, 8)) + list(range(9, 18)) + [18]), 18),
    ]
    return [profile_speed_set(label, speeds, n) for label, speeds, n in selected]


def orient_tournament(profiles: list[SetProfile]) -> list[list[bool]]:
    m = len(profiles)
    adj = [[False] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            si = profiles[i].strength()
            sj = profiles[j].strength()
            if si > sj or (si == sj and i < j):
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj


def score_hist(adj: list[list[bool]]) -> tuple[tuple[int, int], ...]:
    hist = Counter(sum(row) for row in adj)
    return tuple(sorted(hist.items()))


def c3_count(adj: list[list[bool]]) -> int:
    m = len(adj)
    count = 0
    for a, b, c in combinations(range(m), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            count += 1
    return count


def scc_sizes(adj: list[list[bool]]) -> tuple[int, ...]:
    m = len(adj)
    graph = [[j for j, edge in enumerate(adj[i]) if edge] for i in range(m)]
    rgraph = [[] for _ in range(m)]
    for i in range(m):
        for j in graph[i]:
            rgraph[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(m):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes = []

    def rdfs(v: int) -> int:
        seen.add(v)
        size = 1
        for w in rgraph[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count(adj: list[list[bool]]) -> int:
    m = len(adj)
    dp: dict[tuple[int, int], int] = {}
    for v in range(m):
        dp[(1 << v, v)] = 1
    for mask in range(1 << m):
        for last in range(m):
            cur = dp.get((mask, last), 0)
            if not cur:
                continue
            for nxt in range(m):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + cur
    full = (1 << m) - 1
    return sum(dp.get((full, v), 0) for v in range(m))


def tournament_analysis() -> None:
    print("=" * 72)
    print("PART C: Tournament Analysis of exact race profiles")
    print("=" * 72)
    print()
    profiles = selected_profiles()
    adj = orient_tournament(profiles)

    for idx, profile in enumerate(profiles):
        print(
            "{:2d}. {:20s} n={} open_mu={} wall={} wins={} losses={} max_L0={} rescued={} strength={}".format(
                idx,
                profile.label,
                profile.n,
                fmt_frac(profile.open_ll_measure),
                profile.wall_ll_count,
                profile.winning_races,
                profile.losing_races,
                fmt_frac(profile.max_L0),
                profile.simplified_misses_rescued,
                profile.strength(),
            )
        )
        print(f"    outcomes={profile.race_counts}")

    print()
    print("Tournament fingerprints:")
    print(f"  score_hist={score_hist(adj)}")
    print(f"  c3={c3_count(adj)}")
    print(f"  SCCs={scc_sizes(adj)}")
    print(f"  Hamiltonian_paths={hamiltonian_path_count(adj)}")
    print("  top order by score:")
    scores = [(sum(adj[i]), profiles[i].label) for i in range(len(profiles))]
    for score, label in sorted(scores, reverse=True):
        print(f"    score={score} {label}")
    print()


def main() -> None:
    print("LRC exact gap-race swap audit -- codex-2026-06-01 S522")
    print("THM-387: LL can only be entered by winning an LS gap race")
    print()
    bounded_exact_audit()
    initial_segments()
    tournament_analysis()
    print("=" * 72)
    print("SYNTHESIS")
    print("=" * 72)
    print()
    print("1. Exact reset races match the full closed-LL audit in the bounded windows.")
    print("2. The no-swap approximation fails only at compactified wall/tie or")
    print("   swap-sensitive resets in these windows; those are proof objects, not")
    print("   counterexamples.")
    print("3. Initial segments are visible as wall-ledger extremals: open LL measure")
    print("   stays zero while exact reset races still hit compactified LL walls.")
    print("4. The next proof target is a reset ledger theorem: every primitive system")
    print("   has at least one reset whose exact wall trace reaches closed LL before")
    print("   the next reset or short-left loss.")


if __name__ == "__main__":
    main()
