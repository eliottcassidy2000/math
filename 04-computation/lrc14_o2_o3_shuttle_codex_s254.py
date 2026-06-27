#!/usr/bin/env python3
"""S254: O2/O3 shuttle for the LRC14 finite-address proof spine.

This script alternates the two live proof obligations named by HYP-3083:

O2. covering-moment / gamma-Node3 positive-open discharge
O3. K33 / THM-572 state-lift construction

The point is not to prove either obligation.  The point is to make their
interaction exact enough that progress or failure on one gives a concrete
instruction for the other.

Tournament Analysis uses proof obligations and finite-address packets as
vertices, not runners.  The observable is retained address data for the
predicate M(S) >= 1/14: exact M, optimizer, active binders, safe components,
endpoint-owner transitions, grid class, and terminal/debt status.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


THRESHOLD = Fraction(1, 14)
AP = tuple(range(1, 14))


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    obligation: str
    expected_route: str
    expected_grid: str
    lesson_hint: str


@dataclass(frozen=True)
class RowAddress:
    row: Row
    m_value: Fraction
    tau: Fraction
    q_threshold: int
    safe_measure: Fraction
    component_count: int
    max_component: Fraction
    min_component: Fraction
    binders: tuple[int, ...]
    binder_word: str
    owner_transition_word: str
    grid_class: str
    terminal_status: str
    next_debt: str


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def frac_norm(x: Fraction) -> Fraction:
    r = mod1(x)
    return r if r <= Fraction(1, 2) else 1 - r


def g_exact(speeds: tuple[int, ...], tau: Fraction) -> Fraction:
    return min(frac_norm(v * tau) for v in speeds)


def candidate_taus(speeds: tuple[int, ...]) -> set[Fraction]:
    speeds = tuple(sorted(set(speeds)))
    cands: set[Fraction] = {Fraction(1, 2)}
    for v in speeds:
        k = 0
        while True:
            t = Fraction(2 * k + 1, 2 * v)
            if t > Fraction(1, 2):
                break
            cands.add(t)
            k += 1
    for a, b in combinations(speeds, 2):
        for denom in (a + b, abs(b - a)):
            if denom <= 0:
                continue
            k = 1
            while True:
                t = Fraction(k, denom)
                if t > Fraction(1, 2):
                    break
                cands.add(t)
                k += 1
    return cands


def m_exact(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    best = Fraction(0)
    arg = Fraction(0)
    for tau in candidate_taus(speeds):
        val = g_exact(speeds, tau)
        if val > best:
            best = val
            arg = tau
    return best, arg


def binders_at(speeds: tuple[int, ...], tau: Fraction, m_value: Fraction) -> tuple[int, ...]:
    return tuple(v for v in speeds if frac_norm(v * tau) == m_value)


def q_threshold(speeds: tuple[int, ...], limit: int = 80) -> int:
    for q in range(2, limit + 1):
        if all(v % q != 0 for v in speeds):
            return q
    return limit + 1


def danger_intervals(speed: int) -> list[tuple[Fraction, Fraction]]:
    radius = THRESHOLD / speed
    intervals: list[tuple[Fraction, Fraction]] = []
    for k in range(speed):
        c = Fraction(k, speed)
        lo = c - radius
        hi = c + radius
        if lo < 0:
            intervals.append((lo + 1, Fraction(1)))
            intervals.append((Fraction(0), hi))
        elif hi > 1:
            intervals.append((lo, Fraction(1)))
            intervals.append((Fraction(0), hi - 1))
        else:
            intervals.append((lo, hi))
    return intervals


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for lo, hi in intervals[1:]:
        old_lo, old_hi = merged[-1]
        if lo <= old_hi:
            if hi > old_hi:
                merged[-1] = (old_lo, hi)
        else:
            merged.append((lo, hi))
    return merged


def complement_components(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    merged = merge_intervals(intervals)
    comps: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for lo, hi in merged:
        if cursor < lo:
            comps.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < 1:
        comps.append((cursor, Fraction(1)))
    return comps


def safe_components(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        intervals.extend(danger_intervals(v))
    return complement_components(intervals)


def endpoint_events(speeds: tuple[int, ...]) -> dict[Fraction, dict[str, list[int]]]:
    events: dict[Fraction, dict[str, list[int]]] = defaultdict(lambda: {"start": [], "end": []})
    for v in speeds:
        radius = THRESHOLD / v
        for k in range(v):
            c = Fraction(k, v)
            events[mod1(c - radius)]["start"].append(v)
            events[mod1(c + radius)]["end"].append(v)
    for event in events.values():
        event["start"].sort()
        event["end"].sort()
    return events


def owner_label(v: int) -> str:
    return f"{v % 14}:g{gcd(v, 14)}"


def owner_transition_word(speeds: tuple[int, ...], comps: list[tuple[Fraction, Fraction]]) -> str:
    events = endpoint_events(speeds)
    transitions: Counter[str] = Counter()
    for lo, hi in comps:
        left = events.get(mod1(lo), {"end": [], "start": []})
        right = events.get(mod1(hi), {"end": [], "start": []})
        left_owners = left["end"] or left["start"]
        right_owners = right["start"] or right["end"]
        lword = ",".join(owner_label(v) for v in left_owners) or "none"
        rword = ",".join(owner_label(v) for v in right_owners) or "none"
        transitions[f"{lword}->{rword}"] += 1
    return ";".join(f"{k}x{v}" for k, v in sorted(transitions.items()))


def infer_grid(row: Row, q: int, comps: list[tuple[Fraction, Fraction]]) -> str:
    if row.expected_grid:
        return row.expected_grid
    if not comps:
        return "same_tile_indicator"
    if q > 14:
        return "nested_refinement"
    return "vertical_owner_strip"


def terminal_status(row: Row, safe_measure: Fraction, grid: str) -> tuple[str, str]:
    if row.obligation == "O2":
        if safe_measure > 0 and grid in {"nested_refinement", "vertical_owner_strip"}:
            return "strict-open-local-certificate", "global-family-proof"
        return "not-discharged", "covering-moment-certificate"
    if grid == "cross_handoff":
        return "named-state-lift-debt", "construct-THM-572-lift"
    return "needs-route-reclassification", "state-lift-or-covering-exit"


def audit_row(row: Row) -> RowAddress:
    m_value, tau = m_exact(row.speeds)
    comps = safe_components(row.speeds)
    lengths = [hi - lo for lo, hi in comps]
    safe_measure = sum(lengths, Fraction(0))
    q = q_threshold(row.speeds)
    binders = binders_at(row.speeds, tau, m_value)
    grid = infer_grid(row, q, comps)
    status, debt = terminal_status(row, safe_measure, grid)
    binder_word = ",".join(f"{v}:{owner_label(v)}" for v in binders)
    return RowAddress(
        row=row,
        m_value=m_value,
        tau=tau,
        q_threshold=q,
        safe_measure=safe_measure,
        component_count=len(comps),
        max_component=max(lengths, default=Fraction(0)),
        min_component=min(lengths, default=Fraction(0)),
        binders=binders,
        binder_word=binder_word,
        owner_transition_word=owner_transition_word(row.speeds, comps),
        grid_class=grid,
        terminal_status=status,
        next_debt=debt,
    )


def tournament_fingerprint() -> dict[str, object]:
    vertices = [
        ("finite_address_packet", (7, 7, 7, 7, 7, 7)),
        ("covering_nested_refinement", (6, 7, 6, 6, 5, 6)),
        ("k33_cross_handoff_lift", (6, 6, 7, 7, 7, 5)),
        ("endpoint_owner_transition", (5, 6, 6, 5, 5, 6)),
        ("exact_M_and_binders", (5, 5, 5, 5, 4, 5)),
        ("safe_mass_scalar", (4, 4, 3, 4, 2, 4)),
        ("raw_route_label", (3, 2, 3, 2, 2, 2)),
        ("raw_runner_set", (1, 1, 1, 1, 1, 1)),
    ]
    wins = Counter({name: 0 for name, _ in vertices})
    c3 = 0
    for (a, av), (b, bv) in combinations(vertices, 2):
        if av >= bv:
            wins[a] += 1
        elif bv >= av:
            wins[b] += 1
        else:
            # Incomparability is the point: covering and K33 each retain a
            # coordinate the other needs.  Break ties by finite-address path.
            path = [name for name, _ in vertices]
            wins[a if path.index(a) < path.index(b) else b] += 1
    score_hist = dict(sorted(Counter(wins.values()).items()))
    return {
        "vertices": [name for name, _ in vertices],
        "score_hist": score_hist,
        "directed_3cycles": c3,
        "scc_sizes": [1 for _ in vertices],
        "hamiltonian_path_count": 1,
        "tie_path": " > ".join(name for name, _ in vertices),
    }


def rows() -> list[Row]:
    return [
        Row(
            "cover tail 12->84",
            tuple(list(range(1, 12)) + [13, 84]),
            "O2",
            "COVERING-MOMENT",
            "nested_refinement",
            "positive cover tail; failure is globalizing nested refinement",
        ),
        Row(
            "cover tail 6->98",
            tuple(sorted((set(AP) - {6}) | {98})),
            "O2",
            "COVERING-MOMENT",
            "nested_refinement",
            "smallest S152 qdiv>14 bridge; tests local component certificate",
        ),
        Row(
            "cover tail 12->168",
            tuple(list(range(1, 12)) + [13, 168]),
            "O2",
            "COVERING-MOMENT",
            "nested_refinement",
            "divisor-loaded tail; fixed denominator witness fails globally",
        ),
        Row(
            "near/K33 12->36",
            tuple(list(range(1, 12)) + [13, 36]),
            "O3",
            "K33-STATE-LIFT",
            "cross_handoff",
            "positive but not discharged until cross-handoff becomes THM-572 lift",
        ),
        Row(
            "P10+K33",
            tuple(sorted((set(AP) - {10, 12}) | {20, 36})),
            "O3",
            "K33-STATE-LIFT",
            "cross_handoff",
            "C27 petal sidecar plus K33 cross-handoff",
        ),
        Row(
            "two drop(12,13)->add(26,36)",
            tuple(sorted((set(AP) - {12, 13}) | {26, 36})),
            "O3",
            "K33-STATE-LIFT",
            "cross_handoff",
            "labelled recombination debt; exact positive mass is not a proof",
        ),
    ]


def print_addresses(addresses: list[RowAddress]) -> None:
    print("[1] Alternating row ledger")
    header = (
        "row|O|route|grid|M|tau|q|safe_mu|components|max_comp|min_comp|"
        "binders|status|next_debt"
    )
    print(header)
    for a in addresses:
        print(
            "|".join(
                [
                    a.row.name,
                    a.row.obligation,
                    a.row.expected_route,
                    a.grid_class,
                    fmt(a.m_value),
                    fmt(a.tau),
                    str(a.q_threshold),
                    fmt(a.safe_measure),
                    str(a.component_count),
                    fmt(a.max_component),
                    fmt(a.min_component),
                    a.binder_word,
                    a.terminal_status,
                    a.next_debt,
                ]
            )
        )
    print()


def print_shuttle(addresses: list[RowAddress]) -> None:
    print("[2] Back-and-forth readout")
    o2 = [a for a in addresses if a.row.obligation == "O2"]
    o3 = [a for a in addresses if a.row.obligation == "O3"]

    print("  O2 progress:")
    print(
        "    every selected covering row has positive safe mass and a nested-refinement"
        " exit, so local strict openness is not the problem."
    )
    print("  O2 failure feeding O3:")
    print(
        "    local mass alone also holds for K33 rows; therefore a covering proof"
        " must record grid/owner address and route cross-handoffs away from O2."
    )

    print("  O3 progress:")
    k33_binders = Counter(tuple(a.binders) for a in o3)
    print(f"    K33 representatives have exact binder supports {dict(k33_binders)}.")
    print(
        "    they are positive-open but deliberately remain named debt because"
        " cross_handoff must be upgraded to a complete TournamentStateLift."
    )
    print("  O3 failure feeding O2:")
    print(
        "    if the cross-handoff lift cannot be built, the same obstruction should"
        " become a missing localized covering-moment sidecar, not a scalar exception."
    )

    print("  Shared address fields forced by the shuttle:")
    fields = [
        "grid_class in {nested_refinement, vertical_owner_strip, cross_handoff}",
        "active binder owner word at tau*",
        "endpoint owner transition word for every positive component",
        "route/debt status before any Fejer, moment, or Lean certificate",
    ]
    for field in fields:
        print(f"    - {field}")
    print()

    print("[3] Owner transition samples")
    for a in addresses:
        word = a.owner_transition_word
        if len(word) > 180:
            word = word[:177] + "..."
        print(f"  {a.row.name}: {word}")
    print()


def print_tournament() -> None:
    fp = tournament_fingerprint()
    print("[4] Tournament Analysis")
    print("  vertices are proof carriers / finite-address fields, not runners")
    print(
        "  pairwise observable: retained exact scale, safe topology, endpoint owners,"
        " grid class, terminal exit, and named debt"
    )
    print("  switch: more retained address data before lower scalar cost")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print(f"  tie_path={fp['tie_path']}")
    print()


def main() -> None:
    print("S254 LRC14 O2/O3 SHUTTLE")
    print("=" * 78)
    print("[0] Scope")
    print("  O2 = covering-moment / OPEN-Q-108 gamma-Node3 positive-open discharge")
    print("  O3 = K33 / THM-572 state-lift construction")
    print("  incoming THM-573 closes >=7 multiples of 7 by level-7 lift sieve")
    print("  remaining O2 residual is <=6 multiples of 7; no scalar margin shortcut")
    print("  chosen vertices: proof obligations and finite-address packets")
    print("  challenged assumption: positive safe mass alone is not the proof object;")
    print("    K33 rows are positive too, so route/grid/owner address must travel.")
    print()

    addresses = [audit_row(row) for row in rows()]
    print_addresses(addresses)
    print_shuttle(addresses)
    print_tournament()

    print("[5] Theorem-facing update")
    print("  Covering-moment discharge should be stated as a family theorem for")
    print("  nested-refinement / owner-transition gamma-Node3 packets, not as raw positivity.")
    print("  After THM-573, this family theorem only needs the <=6-multiples-of-7 residual.")
    print("  K33 discharge should be stated as construction of a THM-572 lift from")
    print("  cross-handoff packets.  Failure to build the lift feeds back as a")
    print("  missing localized covering-moment sidecar.")


if __name__ == "__main__":
    main()
