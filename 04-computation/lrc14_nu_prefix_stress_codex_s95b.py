#!/usr/bin/env python3
"""
lrc14_nu_prefix_stress_codex_s95b.py

Follow-up stress test for HYP-2866.  The S95 audit found that individual
least-denominator buckets are not monotone, but every bounded-core prefix

    B_Q(E) = sum_{least denominator q <= Q} width_q(E)

was dominated by the consecutive block C_k in the exact `[0,14]` bank.  This
script looks for the first violation outside that bank.

It deliberately reuses the exact S95 dense-interval engine.  This is not a
proof; it is a falsification-oriented scout over a larger exact window plus
structured wide families and reproducible random wide rows.

Tournament Analysis declaration.
  Vertices: stress-bank generators, not runners.
  Pairwise observable:
    (zero_prefix_violations, zero_total_violations, smallest_prefix_slack,
     smallest_total_slack, max_density_ratio_seen)
  Gauge: lexicographic comparison, with smaller density ratio preferred only
    after violation/slack coordinates.  Ties follow the declared bank order.
  Fingerprints: score histogram, directed 3-cycles, SCC sizes, Hamiltonian path
    count.

Assumption challenge.
  I considered runners, individual q-buckets, prefix ledgers, high-q tails,
  random rows, one-tail AP collars, two-tail collars, two-block rows, Fourier
  channels, and proof obligations as vertices.  The chosen vertices are stress
  banks because this script's predicate is falsification coverage: did a whole
  family find a prefix violation?  This preserves the HYP-2866 obstruction
  search and destroys row-level geometry except through the exact witness rows
  reported when a bank is tight.
"""

from __future__ import annotations

import importlib.util
import random
import sys
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

HERE = Path(__file__).resolve().parent
S95_PATH = HERE / "lrc14_nu_denom_center_budget_codex_s95.py"
SPEC = importlib.util.spec_from_file_location("s95_denom_budget", S95_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {S95_PATH}")
s95 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(s95)

RNG = random.Random(20260622)
EXACT_W = 15
RANDOM_ROWS_PER_K = 60
RANDOM_W = 45
TAIL_MAX = 50


def measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((b - a for a, b in intervals), F(0))


def widths_for(E: tuple[int, ...]) -> dict[int | None, F]:
    intervals = s95.dense_intervals(E)
    if not intervals:
        return {}
    return s95.denom_budget_from_intervals(intervals)[1]


def prefix(widths: dict[int | None, F], qcut: int) -> F:
    return sum((w for q, w in widths.items() if q is not None and q <= qcut), F(0))


CONSEC_DATA: dict[int, tuple[F, dict[int | None, F], dict[int, F]]] = {}
for k in range(8, 14):
    C = tuple(range(k))
    widths = widths_for(C)
    D = sum(widths.values(), F(0))
    prefixes = {qcut: prefix(widths, qcut) for qcut in range(7, k + 1)}
    CONSEC_DATA[k] = (D, widths, prefixes)


class BankStats:
    def __init__(self, name: str) -> None:
        self.name = name
        self.checked = 0
        self.prefix_violations: list[tuple[int, int, F, tuple[int, ...], dict[int | None, F]]] = []
        self.total_violations: list[tuple[int, F, tuple[int, ...], dict[int | None, F]]] = []
        self.min_prefix_slack: tuple[F | None, int | None, int | None, tuple[int, ...] | None] = (None, None, None, None)
        self.min_total_slack: tuple[F | None, int | None, tuple[int, ...] | None] = (None, None, None)
        self.max_ratio: tuple[F, int | None, tuple[int, ...] | None] = (F(0), None, None)

    def record(self, k: int, E: tuple[int, ...]) -> None:
        if len(E) != k or E[0] != 0 or not s95.primitive(E):
            return
        self.checked += 1
        widths = widths_for(E)
        D = sum(widths.values(), F(0))
        consec_D, _, consec_prefixes = CONSEC_DATA[k]
        total_slack = consec_D - D
        if total_slack < 0:
            self.total_violations.append((k, total_slack, E, widths))
        old_total, _, _ = self.min_total_slack
        if old_total is None or total_slack < old_total:
            self.min_total_slack = (total_slack, k, E)
        ratio = D / consec_D if consec_D else F(0)
        if ratio > self.max_ratio[0] and E != tuple(range(k)):
            self.max_ratio = (ratio, k, E)
        for qcut in range(7, k + 1):
            slack = consec_prefixes[qcut] - prefix(widths, qcut)
            if slack < 0:
                self.prefix_violations.append((k, qcut, slack, E, widths))
            old_prefix, _, _, _ = self.min_prefix_slack
            if old_prefix is None or slack < old_prefix:
                self.min_prefix_slack = (slack, k, qcut, E)

    def score_tuple(self) -> tuple[int, int, F, F, F]:
        p_ok = 1 if not self.prefix_violations else 0
        t_ok = 1 if not self.total_violations else 0
        p_slack = self.min_prefix_slack[0] if self.min_prefix_slack[0] is not None else F(0)
        t_slack = self.min_total_slack[0] if self.min_total_slack[0] is not None else F(0)
        # Smaller max ratio is safer; negate so lexicographic larger is better.
        return (p_ok, t_ok, p_slack, t_slack, -self.max_ratio[0])


def exact_window_bank() -> BankStats:
    stats = BankStats(f"exact_window_W{EXACT_W}")
    for k in range(8, 14):
        print(f"[bank exact] k={k}, W={EXACT_W}")
        for tail in combinations(range(1, EXACT_W + 1), k - 1):
            stats.record(k, (0,) + tail)
    return stats


def one_tail_bank() -> BankStats:
    stats = BankStats("one_tail_consecutive_collar")
    for k in range(8, 14):
        base = tuple(range(k - 1))
        for far in range(k, TAIL_MAX + 1):
            stats.record(k, tuple(sorted(base + (far,))))
            stats.record(k, tuple(sorted((0,) + tuple(range(2, k + 1)) + (far,))))
    return stats


def two_tail_bank() -> BankStats:
    stats = BankStats("two_tail_consecutive_collar")
    for k in range(8, 14):
        base = tuple(range(k - 2))
        for far1 in range(k - 1, min(35, TAIL_MAX)):
            for far2 in range(far1 + 1, min(far1 + 9, TAIL_MAX + 1)):
                stats.record(k, tuple(sorted(base + (far1, far2))))
    return stats


def two_block_bank() -> BankStats:
    stats = BankStats("two_block_rows")
    for k in range(8, 14):
        for left_len in range(2, k - 1):
            right_len = k - left_len
            left = tuple(range(left_len))
            for gap in range(1, 36):
                start = left_len + gap
                right = tuple(range(start, start + right_len))
                stats.record(k, left + right)
    return stats


def random_wide_bank() -> BankStats:
    stats = BankStats("random_wide_rows")
    for k in range(8, 14):
        seen: set[tuple[int, ...]] = set()
        while len(seen) < RANDOM_ROWS_PER_K:
            E = tuple(sorted((0,) + tuple(RNG.sample(range(1, RANDOM_W + 1), k - 1))))
            if E in seen or not s95.primitive(E):
                continue
            seen.add(E)
            stats.record(k, E)
    return stats


def print_bank(stats: BankStats) -> None:
    print(f"\n--- {stats.name} ---")
    print(f"checked={stats.checked}")
    print(f"prefix violations={len(stats.prefix_violations)}")
    print(f"total D violations={len(stats.total_violations)}")
    p_slack, p_k, p_q, p_E = stats.min_prefix_slack
    t_slack, t_k, t_E = stats.min_total_slack
    ratio, r_k, r_E = stats.max_ratio
    print(f"min prefix slack={p_slack} ({s95.frac_float(p_slack) if p_slack is not None else 'n/a'}), k={p_k}, Q={p_q}, E={p_E}")
    print(f"min total slack={t_slack} ({s95.frac_float(t_slack) if t_slack is not None else 'n/a'}), k={t_k}, E={t_E}")
    print(f"max nonconsec D/consec ratio={ratio} ({s95.frac_float(ratio)}), k={r_k}, E={r_E}")
    if stats.prefix_violations[:3]:
        print("first prefix violations:")
        for k, qcut, slack, E, widths in stats.prefix_violations[:3]:
            print(f"  k={k}, Q={qcut}, slack={slack}, E={E}, widths={s95.fmt_widths(widths, 8)}")
    if stats.total_violations[:3]:
        print("first total violations:")
        for k, slack, E, widths in stats.total_violations[:3]:
            print(f"  k={k}, slack={slack}, E={E}, widths={s95.fmt_widths(widths, 8)}")


def strongly_connected_components(vertices: list[str], winner: dict[tuple[str, str], str]) -> list[list[str]]:
    adj = {v: [] for v in vertices}
    radj = {v: [] for v in vertices}
    for (a, b), w in winner.items():
        loser = b if w == a else a
        adj[w].append(loser)
        radj[loser].append(w)
    order: list[str] = []
    seen: set[str] = set()

    def dfs(v: str) -> None:
        seen.add(v)
        for u in adj[v]:
            if u not in seen:
                dfs(u)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)
    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for u in radj[v]:
            if u not in seen:
                rdfs(u, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(sorted(comp))
    return comps


def count_directed_triangles(vertices: list[str], winner: dict[tuple[str, str], str]) -> int:
    total = 0
    for a, b, c in combinations(vertices, 3):
        out = Counter(
            [
                winner[(min(a, b), max(a, b))],
                winner[(min(a, c), max(a, c))],
                winner[(min(b, c), max(b, c))],
            ]
        )
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def edge(winner: dict[tuple[str, str], str], a: str, b: str) -> bool:
    return winner[(min(a, b), max(a, b))] == a


def count_hamiltonian_paths(vertices: list[str], winner: dict[tuple[str, str], str]) -> int:
    return sum(
        1
        for path in permutations(vertices)
        if all(edge(winner, path[i], path[i + 1]) for i in range(len(path) - 1))
    )


def tournament_analysis(stats_list: list[BankStats]) -> None:
    vertices = [s.name for s in stats_list]
    by_name = {s.name: s for s in stats_list}
    order = {name: i for i, name in enumerate(vertices)}
    wins: dict[tuple[str, str], str] = {}
    score = Counter()
    for a, b in combinations(sorted(vertices), 2):
        sa, sb = by_name[a].score_tuple(), by_name[b].score_tuple()
        if sa > sb:
            w = a
        elif sb > sa:
            w = b
        else:
            w = a if order[a] < order[b] else b
        wins[(a, b)] = w
        score[w] += 1
    hist = Counter(score[v] for v in vertices)
    print("\nTournament Analysis: stress-bank coverage")
    print("  score histogram:", dict(sorted(hist.items())))
    print("  directed 3-cycles:", count_directed_triangles(vertices, wins))
    print("  SCC sizes:", [len(c) for c in strongly_connected_components(vertices, wins)])
    print("  Hamiltonian path count:", count_hamiltonian_paths(vertices, wins))
    ranked = sorted(vertices, key=lambda v: (score[v], by_name[v].score_tuple()), reverse=True)
    print("  leading Hamiltonian path:", " -> ".join(ranked))


def main() -> None:
    print("=" * 78)
    print("HYP-2866 denominator-prefix stress test")
    print("=" * 78)
    print(
        f"exact_W={EXACT_W}, one/two-tail max={TAIL_MAX}, "
        f"random rows per k={RANDOM_ROWS_PER_K}, random_W={RANDOM_W}"
    )
    bank_fns = [
        exact_window_bank,
        one_tail_bank,
        two_tail_bank,
        two_block_bank,
        random_wide_bank,
    ]
    banks = []
    for fn in bank_fns:
        print(f"\n[run bank] {fn.__name__}")
        banks.append(fn())
    for stats in banks:
        print_bank(stats)
    total_prefix = sum(len(s.prefix_violations) for s in banks)
    total_total = sum(len(s.total_violations) for s in banks)
    print("\nSummary:")
    print(f"  aggregate checked={sum(s.checked for s in banks)}")
    print(f"  aggregate prefix violations={total_prefix}")
    print(f"  aggregate total D violations={total_total}")
    if total_prefix == 0 and total_total == 0:
        print("  No violations found; HYP-2866 survives this wider stress bank.")
    else:
        print("  Violation found; inspect rows above before using HYP-2866.")
    tournament_analysis(banks)
    print("\nDONE.")


if __name__ == "__main__":
    main()
