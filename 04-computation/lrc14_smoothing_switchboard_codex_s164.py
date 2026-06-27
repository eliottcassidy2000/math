#!/usr/bin/env python3
"""HYP-2984: LRC14 smoothing switchboard.

This is a proof-interface audit, not a proof of LRC14.  It joins the newest
analytic-sieve/Kaczynski lane to the existing Fejer and Ramanujan packet
certificates:

* HYP-2981: selected packet-anchored Fejer interval certificates.
* HYP-2982: Phi/G analytic packet weights and squarefree-blindness warning.
* HYP-2983: exponential-sum/smoothing/Kaczynski proof-template synthesis.
* HYP-2979: Ramanujan exact-period primitive witnesses.

The output is a finite route matrix for named packet families.  Its purpose is
to make the admissible-smoothing lemma less vague: each packet is assigned a
kernel route, a retained side channel, and a forbidden scalarization.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s162 = load_module(
    "s164_fejer_scaffold",
    REPO / "04-computation" / "lrc14_packet_fejer_interval_scaffold_codex_s162.py",
)
s2979 = load_module(
    "s164_ramanujan_projector",
    REPO / "04-computation" / "lrc14_ramanujan_exact_period_projector_codex_20260624.py",
)


def qstr(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


@dataclass(frozen=True)
class AuditRow:
    name: str
    source: str
    route: str
    speeds: tuple[int, ...]
    fejer_degree: int | None = None


@dataclass(frozen=True)
class RouteRecord:
    row: AuditRow
    q_threshold: int
    safe_mu: Fraction
    component_count: int
    largest_width: Fraction
    first_weak_q: int | None
    first_strict_q: int | None
    q14_weak: int
    q14_strict: int
    q14_boundary_only: int
    fejer_certified: bool
    fejer_degree: int | None
    switch_route: str
    smoothing_kernel: str
    retained_side_channel: str
    forbidden_forgetting: str


def one(drop: int, add: int) -> tuple[int, ...]:
    return tuple(sorted((set(s2979.AP) - {drop}) | {add}))


def many(drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(s2979.AP) - set(drops)) | set(adds)))


def audit_rows() -> list[AuditRow]:
    rows: dict[tuple[int, ...], AuditRow] = {}

    for row in s2979.named_rows():
        rows[row.speeds] = AuditRow(row.name, "HYP-2979 named", row.route, row.speeds)

    # Add the selected interval-certificate rows, preserving the Fejer degree.
    for row in s162.default_rows():
        route = {
            "K33 state-lift": "K33-STATE-LIFT",
            "two-block splice": "BOUNDARY-PETAL",
            "covering comb": "COVERING-MOMENT",
            "two-swap": "Q-WITNESS",
            "one-swap AP bank": "ONE-SWAP",
        }.get(row.source_family, row.source_family)
        rows[row.speeds] = AuditRow(
            row.name,
            "HYP-2981 Fejer selected",
            route,
            row.speeds,
            row.degree,
        )

    # A few route-probing rows that exercise the switchboard decisions.
    extras = [
        AuditRow("petal splice P13+GW", "S164 probe", "BOUNDARY-PETAL", many((12, 13), (24, 26))),
        AuditRow("q14 liar 12->96", "S164 probe", "Q14-FRONT", one(12, 96)),
        AuditRow("committed wall 84", "S164 probe", "LATE-Q-COVERING", one(12, 84)),
        AuditRow("late repeated-prime 12->280", "S164 probe", "Q-WITNESS-LATE-FACTOR", one(12, 280)),
    ]
    for row in extras:
        rows.setdefault(row.speeds, row)

    return list(rows.values())


def fejer_certified(row: AuditRow, safe_center: Fraction, terms: int = 36) -> bool:
    if row.fejer_degree is None:
        return False
    if safe_center == 0:
        return False
    iv = s162.fejer_interval(row.speeds, row.fejer_degree, safe_center, terms)
    return iv.hi < 0


def route_record(row: AuditRow) -> RouteRecord:
    safe_mu, largest_width, center, component_count = s162.safe_center(row.speeds)
    q_threshold = s2979.q_threshold(row.speeds)
    first_weak = s2979.first_witness_q(row.speeds, qmax=42, strict=False)
    first_strict = s2979.first_witness_q(row.speeds, qmax=42, strict=True)
    q14_weak, q14_strict, q14_boundary, _ = s2979.phase_summary(row.speeds, 14)
    fejer_ok = fejer_certified(row, center)

    if safe_mu == 0:
        switch_route = "AP_GW_BOUNDARY_ATOM"
        smoothing_kernel = "closed-boundary endpoint kernel"
        retained = "endpoint zero-credit pairs plus Kaczynski approach class"
        forbidden = "open-measure averaging or raw Fejer negativity"
    elif fejer_ok:
        switch_route = "INTERVAL_FEJER_CERTIFICATE"
        smoothing_kernel = "packet-centered Fejer kernel"
        retained = "packet key, rational center, degree, interval upper bound"
        forbidden = "floating Fejer value without interval/packet label"
    elif "K33" in row.route:
        switch_route = "K33_STATE_LIFT"
        smoothing_kernel = "state-lift boundary/resonance kernel"
        retained = "K33 owner packet plus HYP-2908/THM-572 debt"
        forbidden = "large-sieve scalar tail without state-lift side channel"
    elif "COVERING" in row.route or "LATE-Q" in row.route:
        switch_route = "COVERING_LIFT_OR_BOUNDARY_MOMENT"
        smoothing_kernel = "multi-chart lift / boundary-moment kernel"
        retained = "lift chart, endpoint owners, exact safe mass, late q label"
        forbidden = "squarefree mu^2/phi collapse of repeated-prime q"
    elif first_strict is not None:
        switch_route = "RAMANUJAN_PRE_SPLIT_THEN_FEJER"
        smoothing_kernel = "exact-period Ramanujan pre-split plus Fejer handoff"
        retained = "first strict q, primitive phase packet, endpoint labels"
        forbidden = "qdiv-only or Ramanujan scalar profile"
    else:
        switch_route = "UNCLASSIFIED_RESONANCE_DEBT"
        smoothing_kernel = "adaptive Kaczynski defect kernel"
        retained = "all labels until new side channel is found"
        forbidden = "declaring a density estimate theorem-safe"

    return RouteRecord(
        row=row,
        q_threshold=q_threshold,
        safe_mu=safe_mu,
        component_count=component_count,
        largest_width=largest_width,
        first_weak_q=first_weak,
        first_strict_q=first_strict,
        q14_weak=q14_weak,
        q14_strict=q14_strict,
        q14_boundary_only=q14_boundary,
        fejer_certified=fejer_ok,
        fejer_degree=row.fejer_degree,
        switch_route=switch_route,
        smoothing_kernel=smoothing_kernel,
        retained_side_channel=retained,
        forbidden_forgetting=forbidden,
    )


def beats(edges: dict[tuple[int, int], int], i: int, j: int) -> bool:
    if i < j:
        return edges[(i, j)] == 1
    return edges[(j, i)] == 0


def score_hist(edges: dict[tuple[int, int], int], n: int) -> dict[int, int]:
    scores = Counter({i: 0 for i in range(n)})
    for i, j in combinations(range(n), 2):
        scores[i if beats(edges, i, j) else j] += 1
    return dict(sorted(Counter(scores.values()).items()))


def directed_3cycles(edges: dict[tuple[int, int], int], n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        out = Counter()
        for i, j in ((a, b), (a, c), (b, c)):
            out[i if beats(edges, i, j) else j] += 1
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def hamiltonian_count(edges: dict[tuple[int, int], int], n: int) -> int:
    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for size in range(2, n + 1):
        nxt: dict[tuple[int, int], int] = {}
        for mask in range(1 << n):
            if mask.bit_count() != size:
                continue
            for last in range(n):
                if not ((mask >> last) & 1):
                    continue
                prev_mask = mask ^ (1 << last)
                total = 0
                for prev in range(n):
                    if ((prev_mask >> prev) & 1) and beats(edges, prev, last):
                        total += dp.get((prev_mask, prev), 0)
                if total:
                    nxt[(mask, last)] = total
        dp.update(nxt)
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def carrier_tournament():
    vertices = [
        ("raw_scalar_smoothing_choice", (1, 0, 0, 0, 0, 0, 0), "chooses a kernel but forgets route labels"),
        ("large_sieve_middle_bound", (1, 1, 0, 0, 0, 0, 0), "bounds a family but needs residual labels"),
        ("ramanujan_exact_period_presplit", (1, 1, 1, 0, 0, 0, 0), "keeps primitive q packets before smoothing"),
        ("kaczynski_boundary_defect", (1, 1, 1, 1, 0, 0, 0), "keeps approach class for boundary ambiguity"),
        ("boundary_moment_lift_chart", (1, 1, 1, 1, 1, 0, 0), "keeps lift chart and endpoint owners"),
        ("interval_fejer_certificate", (1, 1, 1, 1, 1, 1, 0), "keeps center, degree, and interval sign"),
        ("labelled_smoothing_switchboard", (1, 1, 1, 1, 1, 1, 1), "chooses kernel after packet route is known"),
    ]
    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(len(vertices)), 2):
        edges[(i, j)] = 1 if vertices[i][1] > vertices[j][1] else 0
    return vertices, edges


def main() -> None:
    records = [route_record(row) for row in audit_rows()]
    records.sort(key=lambda r: (r.switch_route, r.row.name))

    print("HYP-2984: LRC14 smoothing switchboard")
    print("=" * 78)
    print("Rows are packet/proof obligations, not runners.")
    print("Each route declares a smoothing kernel and the side channel it must retain.")
    print()

    route_counts = Counter(rec.switch_route for rec in records)
    print("Route counts")
    for route, count in sorted(route_counts.items()):
        print(f"  {route:<38s} {count}")
    print()

    print("Packet route matrix")
    header = (
        "row | route | qthr | safe_mu | comps | q14 weak/strict/bdy | "
        "first weak/strict q | Fejer | switch"
    )
    print(header)
    print("-" * len(header))
    for rec in records:
        fejer = "-"
        if rec.fejer_degree is not None:
            fejer = f"d={rec.fejer_degree}:{'cert' if rec.fejer_certified else 'fail'}"
        print(
            f"{rec.row.name} | {rec.row.route} | {rec.q_threshold} | "
            f"{qstr(rec.safe_mu)} | {rec.component_count} | "
            f"{rec.q14_weak}/{rec.q14_strict}/{rec.q14_boundary_only} | "
            f"{rec.first_weak_q or '-'}/{rec.first_strict_q or '-'} | "
            f"{fejer} | {rec.switch_route}"
        )

    print()
    print("Switch contracts")
    for rec in records:
        print(f"- {rec.row.name}:")
        print(f"    kernel={rec.smoothing_kernel}")
        print(f"    retain={rec.retained_side_channel}")
        print(f"    forbidden_forgetting={rec.forbidden_forgetting}")

    vertices, edges = carrier_tournament()
    n = len(vertices)
    path = " > ".join(name for name, _score, _why in reversed(vertices))
    print()
    print("Tournament Analysis")
    print("vertices: smoothing/proof carriers, not runners")
    print(
        "pairwise observable: retained predicate payload = kernel, q packet,"
        " boundary approach, lift chart, interval sign, and route handoff"
    )
    print("tie Hamiltonian path: listed carrier order")
    print(f"score_hist={score_hist(edges, n)}")
    print(f"directed_3cycles={directed_3cycles(edges, n)}")
    print(f"hamiltonian_paths={hamiltonian_count(edges, n)}")
    print(f"Hamiltonian path: {path}")
    for name, score, why in vertices:
        print(f"  {name:<34s} score={score} :: {why}")

    print()
    print("Proof-pass readout")
    print(
        "The creative compression is a switchboard lemma: choose Fejer,"
        " Ramanujan, boundary-moment, large-sieve, or Kaczynski smoothing only"
        " after the packet route is known.  AP/GW are boundary atoms; selected"
        " positive hard rows already have Fejer interval certificates; covering"
        " and K33 rows keep lift/state debt; late q packets keep exact-period"
        " labels before any squarefree large-sieve normalization."
    )


if __name__ == "__main__":
    main()
