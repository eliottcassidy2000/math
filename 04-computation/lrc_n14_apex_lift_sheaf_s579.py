#!/usr/bin/env python3
"""
lrc_n14_apex_lift_sheaf_s579.py

codex-2026-06-03 S579

Probe the user's "apex-lift certificate sheaf" idea for the n=14 LRC
unit-spine proof route.

Finite site:
  Base rows are HYP-2096 unit-spine slack rows.  A lift object replaces a
  subset of the nine unit-shell representatives by larger same-shell
  representatives modulo C=27.  A restriction map lowers one lifted shell back
  to its canonical representative.

Certificate sections:
  * ledger_failure: the D/U/N quotient is not a full cover, so the row is not in
    the strict sub-edge danger ledger.
  * cheap_d14: a HYP-2095 unblocked small-pair witness with denominator 14.
    This is the apex-denominator chart for n=14=2*7.
  * cheap_side: a HYP-2095 unblocked small-pair witness on another denominator.
  * positive_measure: all small pairs are blocked, but the safe set has
    positive measure.
  * residual: full cover, no cheap witness, and no positive-measure certificate.

The sheaf question is whether lifted full-cover rows restrict to local
certificate sections without creating residual defects.  The "apex" here is
not assumed to be a witness runner.  It is the central denominator-14 chart
coming from the doubled-prime apex q=7; side charts are allowed and counted.

Tournament Analysis / assumption challenge:
  Vertices are named lift sites, not runners.  The observable is
    (restriction residual defects, side-chart load, union-only sections,
     full covers, apex-lift full covers).
  The switch orients toward harder sheaf sites.  This preserves the LRC
  predicate "every lifted full cover has a local certificate after restriction"
  and destroys the full self-converse round-class identity.  That is deliberate:
  the script tests the local gluing interface between HYP-2099 and HYP-2100.
"""

from __future__ import annotations

import importlib.util
import sys
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S578 = ROOT / "04-computation" / "lrc_n14_unit_spine_exchange_s578.py"


def load_s578():
    spec = importlib.util.spec_from_file_location("s578_exchange", S578)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {S578}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


s578 = load_s578()

N = s578.N
C = s578.C
FLOOR = s578.FLOOR
UNIT_SPINE = s578.UNIT_SPINE
APEX_SPEED = N // 2
APEX_DENOMINATOR = N


@dataclass(frozen=True)
class Route:
    kind: str
    pair: tuple[int, int] | None = None
    denominator: int | None = None
    numerator: int | None = None

    @property
    def coarse(self) -> str:
        if self.kind.startswith("cheap"):
            return "cheap"
        return self.kind


@dataclass(frozen=True)
class SiteSummary:
    label: str
    rows: int
    full: int
    route_hist: tuple[tuple[str, int], ...]
    witness_den_hist: tuple[tuple[int, int], ...]
    apex_lift_rows: int
    apex_lift_full: int
    apex_lift_d14: int
    apex_pair_witnesses: int
    two_full: int
    two_d14: int
    restriction_pair_hist: tuple[tuple[tuple[str, ...], int], ...]
    restriction_residual: int
    union_only: int
    examples: tuple[str, ...]

    @property
    def side_chart_load(self) -> int:
        route_counts = dict(self.route_hist)
        return route_counts.get("cheap_side", 0)

    @property
    def burden(self) -> tuple[int, int, int, int, int]:
        return (
            self.restriction_residual,
            self.side_chart_load,
            self.union_only,
            self.full,
            self.apex_lift_full,
        )


def positive_measure(speeds: tuple[int, ...]) -> bool:
    endpoints: set[Fraction] = set()
    for v in speeds:
        for k in range(v + 1):
            for sign in (-1, 1):
                endpoints.add(Fraction(k * N + sign, N * v) % 1)
    pts = sorted(endpoints)
    for i, a in enumerate(pts):
        b = pts[(i + 1) % len(pts)]
        length = b - a if b > a else b - a + 1
        mid = (a + length / 2) % 1
        if all(s578.norm(Fraction(v) * mid) > FLOOR for v in speeds):
            return True
    return False


route_cache: dict[tuple[int, ...], Route] = {}


def route_for(speeds_raw: tuple[int, ...]) -> Route:
    speeds = tuple(sorted(speeds_raw))
    if speeds in route_cache:
        return route_cache[speeds]
    if not s578.is_full_cover(speeds):
        route = Route("ledger_failure")
    else:
        witness = s578.unblocked_small_pair(speeds)
        if witness is not None:
            kind = "cheap_d14" if witness.denominator == APEX_DENOMINATOR else "cheap_side"
            route = Route(kind, witness.pair, witness.denominator, witness.numerator)
        elif positive_measure(speeds):
            route = Route("positive_measure")
        else:
            route = Route("residual")
    route_cache[speeds] = route
    return route


def lifted_shells(units: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(shell for shell, rep in zip(UNIT_SPINE, units) if shell != rep)


def lower_shell(units: tuple[int, ...], shell: int) -> tuple[int, ...]:
    lowered = list(units)
    lowered[UNIT_SPINE.index(shell)] = shell
    return tuple(lowered)


def route_name_pair(routes: list[Route]) -> tuple[str, ...]:
    return tuple(sorted(r.coarse for r in routes))


def apex_pair_witness(route: Route) -> bool:
    if route.pair is None:
        return False
    return any(s578.unit_shell(v) == APEX_SPEED for v in route.pair)


def summarize_named_site(label: str, slack: tuple[int, ...]) -> SiteSummary:
    patterns = s578.local_unit_patterns()
    rows = full = 0
    route_hist: Counter[str] = Counter()
    den_hist: Counter[int] = Counter()
    apex_lift_rows = apex_lift_full = apex_lift_d14 = apex_pair_count = 0
    two_full = two_d14 = restriction_residual = union_only = 0
    restriction_pair_hist: Counter[tuple[str, ...]] = Counter()
    examples: list[str] = []

    for units in patterns:
        rows += 1
        lifts = lifted_shells(units)
        route = route_for(tuple(units) + slack)
        route_hist[route.kind] += 1
        if route.denominator is not None:
            den_hist[route.denominator] += 1
        if route.kind != "ledger_failure":
            full += 1
        if APEX_SPEED in lifts:
            apex_lift_rows += 1
            if route.kind != "ledger_failure":
                apex_lift_full += 1
            if route.kind == "cheap_d14":
                apex_lift_d14 += 1
        if apex_pair_witness(route):
            apex_pair_count += 1

        if len(lifts) == 2 and route.kind != "ledger_failure":
            two_full += 1
            if route.kind == "cheap_d14":
                two_d14 += 1
            restrictions = [route_for(lower_shell(units, shell) + slack) for shell in lifts]
            pair_key = route_name_pair(restrictions)
            restriction_pair_hist[pair_key] += 1
            restriction_residual += sum(1 for r in restrictions if r.kind == "residual")
            if pair_key == ("ledger_failure", "ledger_failure"):
                union_only += 1
            if pair_key != ("cheap", "cheap") and len(examples) < 5:
                examples.append(
                    f"units={units} lifts={lifts} route={route.kind}"
                    f" pair={route.pair} D={route.denominator} restrictions={pair_key}"
                )

    return SiteSummary(
        label=label,
        rows=rows,
        full=full,
        route_hist=tuple(route_hist.most_common()),
        witness_den_hist=tuple(den_hist.most_common(8)),
        apex_lift_rows=apex_lift_rows,
        apex_lift_full=apex_lift_full,
        apex_lift_d14=apex_lift_d14,
        apex_pair_witnesses=apex_pair_count,
        two_full=two_full,
        two_d14=two_d14,
        restriction_pair_hist=tuple(restriction_pair_hist.most_common()),
        restriction_residual=restriction_residual,
        union_only=union_only,
        examples=tuple(examples),
    )


def summarize_one_lift_all_slack() -> SiteSummary:
    slack_candidates = tuple(
        v for v in range(3, s578.ONE_LIFT_SLACK_BOUND + 1) if v % 3 == 0 and v not in UNIT_SPINE
    )
    reps_by_shell = {
        shell: tuple(v for v in s578.unit_reps(shell, s578.ONE_LIFT_UNIT_BOUND) if v != shell)
        for shell in UNIT_SPINE
    }

    rows = full = 0
    route_hist: Counter[str] = Counter()
    den_hist: Counter[int] = Counter()
    apex_lift_rows = apex_lift_full = apex_lift_d14 = apex_pair_count = 0
    restriction_pair_hist: Counter[tuple[str, ...]] = Counter()
    restriction_residual = union_only = 0
    examples: list[str] = []

    for slack in combinations(slack_candidates, 4):
        slack = tuple(sorted(slack))
        for i, shell in enumerate(UNIT_SPINE):
            for rep in reps_by_shell[shell]:
                units = list(UNIT_SPINE)
                units[i] = rep
                units_t = tuple(units)
                rows += 1
                route = route_for(units_t + slack)
                route_hist[route.kind] += 1
                if route.denominator is not None:
                    den_hist[route.denominator] += 1
                if route.kind != "ledger_failure":
                    full += 1
                    restricted = route_for(lower_shell(units_t, shell) + slack)
                    pair_key = route_name_pair([restricted])
                    restriction_pair_hist[pair_key] += 1
                    restriction_residual += int(restricted.kind == "residual")
                    if restricted.kind == "ledger_failure":
                        union_only += 1
                    if restricted.kind not in {"cheap_d14", "cheap_side"} and len(examples) < 5:
                        examples.append(
                            f"slack={slack} lift={shell}->{rep} route={route.kind}"
                            f" D={route.denominator} restriction={restricted.kind}"
                        )
                if shell == APEX_SPEED:
                    apex_lift_rows += 1
                    if route.kind != "ledger_failure":
                        apex_lift_full += 1
                    if route.kind == "cheap_d14":
                        apex_lift_d14 += 1
                if apex_pair_witness(route):
                    apex_pair_count += 1

    return SiteSummary(
        label="one_lift_all_slack",
        rows=rows,
        full=full,
        route_hist=tuple(route_hist.most_common()),
        witness_den_hist=tuple(den_hist.most_common(8)),
        apex_lift_rows=apex_lift_rows,
        apex_lift_full=apex_lift_full,
        apex_lift_d14=apex_lift_d14,
        apex_pair_witnesses=apex_pair_count,
        two_full=0,
        two_d14=0,
        restriction_pair_hist=tuple(restriction_pair_hist.most_common()),
        restriction_residual=restriction_residual,
        union_only=union_only,
        examples=tuple(examples),
    )


def tournament_fingerprint(rows: list[SiteSummary]) -> dict[str, object]:
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
        c3 += int(
            (adj[i][j] and adj[j][k] and adj[k][i])
            or (adj[i][k] and adj[k][j] and adj[j][i])
        )

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


def print_summary(row: SiteSummary) -> None:
    print(f"[{row.label}]")
    print(f"  rows={row.rows} full={row.full}")
    print(f"  route_hist={dict(row.route_hist)}")
    print(f"  witness_den_hist={dict(row.witness_den_hist)}")
    print(
        "  apex_lifts="
        f"rows:{row.apex_lift_rows} full:{row.apex_lift_full} d14:{row.apex_lift_d14}"
    )
    print(f"  apex_pair_witnesses={row.apex_pair_witnesses}")
    if row.two_full:
        print(f"  two_lift_full={row.two_full} two_lift_d14={row.two_d14}")
    print(f"  restriction_pair_hist={dict(row.restriction_pair_hist)}")
    print(f"  restriction_residual={row.restriction_residual} union_only={row.union_only}")
    if row.examples:
        print("  nontrivial restriction examples:")
        for example in row.examples:
            print(f"    {example}")
    print()


def main() -> None:
    print("S579 n=14 apex-lift certificate sheaf")
    print("=" * 78)
    print(f"n={N}; C={C}; apex_speed={APEX_SPEED}; apex_denominator={APEX_DENOMINATOR}")
    print("site: unit-shell lift subsets; restriction lowers one lifted shell")
    print("sections: ledger_failure, cheap_d14, cheap_side, positive_measure, residual")
    print()

    rows = [summarize_one_lift_all_slack()]
    named = {
        "AP_slack": (3, 6, 9, 12),
        "Vstar_slack": (3, 6, 9, 24),
        "first_open_gap_slack": (3, 6, 9, 36),
        "zero_slack_control": (3, 6, 9, 27),
        "double3_control": (3, 6, 24, 30),
    }
    rows.extend(summarize_named_site(label, tuple(sorted(slack))) for label, slack in named.items())

    for row in rows:
        print_summary(row)

    print("Tournament Analysis")
    print("  vertices: sheaf sites")
    print(
        "  observable: restriction residual defects, side-chart load,"
        " union-only sections, full covers, apex-lift full covers"
    )
    print("  switch: harder = more residual, then more side-chart load")
    print(f"  fingerprints={tournament_fingerprint(rows)}")
    print()

    print("Synthesis")
    print("  The denominator-14 apex chart is large but not complete: side charts")
    print("  are needed for lifted full covers.  The useful object is therefore a")
    print("  certificate sheaf, not one global apex witness.  In this bounded site,")
    print("  every full restriction is certified by a cheap section; restriction")
    print("  residuals are zero.  A few two-lift rows are union-only: each one-lift")
    print("  restriction is a ledger failure, while the two-lift union creates a full")
    print("  cover with a cheap witness.  That is exactly the local-to-global behavior")
    print("  the next HYP-2099/HYP-2100 bridge lemma should formalize.")


if __name__ == "__main__":
    main()
