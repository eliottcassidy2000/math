#!/usr/bin/env python3
"""S164: kernel homotopy and boundary-defect ledger for LRC14.

This is a small exact proof-interface scout between two recent lanes:

* HYP-2981: packet-anchored Fejer interval certificates;
* HYP-2982/HYP-2983: analytic smoothing, large-sieve, and Kaczynski boundary
  guardrails.

The point is not to introduce another scalar kernel.  The point is to say when
changing a kernel is theorem-safe.  A kernel deformation is allowed only if it
preserves a labelled packet certificate or records the boundary defect it
creates.  Exact regular-open safe components give a kernel-stability radius;
AP/Goddyn-Wong have no open component and instead produce taut endpoint atoms.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s147 = load_module(
    "s164_baire_haar",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)
s155 = load_module(
    "s164_taut_bridge",
    REPO / "04-computation" / "lrc14_taut_bridge_graph_codex_s155.py",
)


@dataclass(frozen=True)
class RowSpec:
    name: str
    family: str
    route: str
    speeds: tuple[int, ...]


@dataclass(frozen=True)
class KernelRecord:
    row: RowSpec
    M: Fraction
    qdiv: int
    safe_mu: Fraction
    component_count: int
    max_width: Fraction
    stability_radius: Fraction
    taut_count: int
    zero_sum_pairs: int
    bridge_count: int
    state: str
    route: str


def fmt(x: Fraction | None) -> str:
    if x is None:
        return "-"
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g == 1


def row_replace(holes: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(holes)) | set(adds)))


def named_rows() -> list[RowSpec]:
    return [
        RowSpec("AP", "AP/GW equality", "boundary-defect atom", AP),
        RowSpec(
            "GW 12->24",
            "AP/GW equality",
            "boundary-defect atom",
            tuple(list(range(1, 12)) + [13, 24]),
        ),
        RowSpec(
            "near/K33 12->36",
            "K33 state-lift",
            "open-stable, then Fejer/state-lift",
            row_replace((12,), (36,)),
        ),
        RowSpec(
            "petal 10->20",
            "unit petal",
            "open-stable petal discharge",
            row_replace((10,), (20,)),
        ),
        RowSpec(
            "petal 13->26",
            "unit petal",
            "open-stable petal discharge",
            row_replace((13,), (26,)),
        ),
        RowSpec(
            "P10+GW",
            "two-block splice",
            "open-stable unit/GW splice",
            row_replace((10, 12), (20, 24)),
        ),
        RowSpec(
            "P10+K33",
            "two-block K33",
            "open-stable before K33 label",
            row_replace((10, 12), (20, 36)),
        ),
        RowSpec(
            "residue liar 12->26",
            "q-witness",
            "coarser q-threshold witness",
            row_replace((12,), (26,)),
        ),
        RowSpec(
            "magnitude liar 12->96",
            "q14 loose",
            "same apex residues, off-apex open-stable",
            row_replace((12,), (96,)),
        ),
        RowSpec(
            "covering 12->84",
            "covering comb",
            "open-stable covering packet",
            row_replace((12,), (84,)),
        ),
        RowSpec(
            "covering 12->168",
            "covering comb",
            "open-stable small-margin covering packet",
            row_replace((12,), (168,)),
        ),
        RowSpec(
            "few-apex 6->28",
            "few-apex",
            "open-stable few-apex packet",
            row_replace((6,), (28,)),
        ),
    ]


def kernel_record(row: RowSpec) -> KernelRecord:
    exact = s147.exact_row_measure(row.speeds)
    comps = tuple(exact["safe_components"])
    safe_mu = exact["safe_measure"]
    max_width = max((hi - lo for lo, hi in comps), default=Fraction(0))
    stability = max_width / 2 if max_width > 0 else Fraction(0)
    audit = s155.audit_row(row.name, row.family, row.speeds)
    taut_count = len(audit.taut_vertices)
    zero_sum_pairs = sum(len(t.zero_sum_pairs) for t in audit.taut_vertices)
    bridge_count = len(audit.positive_bridges)
    if safe_mu > 0:
        state = "open-stable"
    elif taut_count:
        state = "boundary-defect"
    else:
        state = "covered-residual"
    return KernelRecord(
        row=row,
        M=audit.M,
        qdiv=audit.qdiv,
        safe_mu=safe_mu,
        component_count=len(comps),
        max_width=max_width,
        stability_radius=stability,
        taut_count=taut_count,
        zero_sum_pairs=zero_sum_pairs,
        bridge_count=bridge_count,
        state=state,
        route=audit.route,
    )


def scan_single_swaps(limit: int) -> dict[str, object]:
    seen: set[tuple[int, ...]] = set()
    counts: Counter[str] = Counter()
    min_positive: tuple[Fraction, str, tuple[int, ...]] | None = None
    zero_rows: list[tuple[str, tuple[int, ...]]] = []
    total = 0
    for remove in AP:
        for add in range(14, limit + 1):
            row = row_replace((remove,), (add,))
            if len(row) != 13 or row in seen or not primitive(row):
                continue
            seen.add(row)
            total += 1
            data = s147.exact_row_measure(row)
            safe_mu = data["safe_measure"]
            if safe_mu > 0:
                counts["open-stable"] += 1
                label = f"{remove}->{add}"
                if min_positive is None or safe_mu < min_positive[0]:
                    min_positive = (safe_mu, label, row)
            else:
                counts["zero-open"] += 1
                zero_rows.append((f"{remove}->{add}", row))
    return {
        "limit": limit,
        "total": total,
        "counts": counts,
        "min_positive": min_positive,
        "zero_rows": zero_rows,
    }


RETENTION_SCORES = {
    # strict predicate, boundary predicate, packet labels, kernel audit,
    # formal interval path, scalar-decoy resistance
    "open_component_certificate": (8, 4, 7, 8, 8, 8),
    "boundary_defect_atom": (4, 9, 8, 8, 6, 8),
    "packet_fejer_certificate": (7, 5, 9, 7, 9, 8),
    "kaczynski_approach_class": (5, 8, 7, 8, 5, 7),
    "analytic_smoothing_kernel": (5, 4, 5, 9, 6, 5),
    "ramanujan_exact_period_packet": (5, 6, 8, 5, 6, 7),
    "raw_kernel_scalar": (2, 1, 1, 2, 2, 1),
}

TIE_PATH = (
    "open_component_certificate",
    "boundary_defect_atom",
    "packet_fejer_certificate",
    "kaczynski_approach_class",
    "analytic_smoothing_kernel",
    "ramanujan_exact_period_packet",
    "raw_kernel_scalar",
)


def winner(a: str, b: str) -> str:
    va = RETENTION_SCORES[a]
    vb = RETENTION_SCORES[b]
    score = sum(x > y for x, y in zip(va, vb)) - sum(x < y for x, y in zip(va, vb))
    if score > 0:
        return a
    if score < 0:
        return b
    return a if TIE_PATH.index(a) < TIE_PATH.index(b) else b


def strongly_connected_components(vertices: tuple[str, ...], edges: dict[str, set[str]]) -> list[set[str]]:
    def reach(start: str, graph: dict[str, set[str]]) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in graph[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    rev = {v: set() for v in vertices}
    for v, outs in edges.items():
        for w in outs:
            rev[w].add(v)
    unseen = set(vertices)
    out: list[set[str]] = []
    while unseen:
        v = next(iter(unseen))
        comp = reach(v, edges) & reach(v, rev)
        out.append(comp)
        unseen -= comp
    return out


def tournament_fingerprint() -> dict[str, object]:
    vertices = tuple(RETENTION_SCORES)
    edges = {v: set() for v in vertices}
    wins = Counter()
    for a, b in combinations(vertices, 2):
        w = winner(a, b)
        loser = b if w == a else a
        edges[w].add(loser)
        wins[w] += 1
    c3 = 0
    for a, b, c in combinations(vertices, 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            c3 += 1
        if c in edges[a] and b in edges[c] and a in edges[b]:
            c3 += 1
    sccs = strongly_connected_components(vertices, edges)
    order = sorted(vertices, key=lambda v: (-wins[v], TIE_PATH.index(v)))
    return {
        "score_hist": dict(sorted(Counter(wins[v] for v in vertices).items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted((len(c) for c in sccs), reverse=True),
        "hamiltonian_path": order,
    }


def print_records(records: list[KernelRecord]) -> None:
    print("[1] Named packet kernel-homotopy ledger")
    print(
        "  "
        + "row".ljust(24)
        + "state".ljust(17)
        + "M".rjust(8)
        + "qdiv".rjust(6)
        + "safe_mu".rjust(12)
        + "comp".rjust(6)
        + "max_w".rjust(12)
        + "eps<".rjust(12)
        + "taut".rjust(6)
        + " "
        + "route".rjust(18)
    )
    for rec in records:
        print(
            "  "
            + rec.row.name.ljust(24)
            + rec.state.ljust(17)
            + fmt(rec.M).rjust(8)
            + str(rec.qdiv).rjust(6)
            + fmt(rec.safe_mu).rjust(12)
            + str(rec.component_count).rjust(6)
            + fmt(rec.max_width).rjust(12)
            + fmt(rec.stability_radius).rjust(12)
            + str(rec.taut_count).rjust(6)
            + " "
            + rec.route.rjust(18)
        )
    print()
    print("  Interpretation:")
    print("    eps< is a certificate-safe support radius: a symmetric kernel with")
    print("    support below that value can be centered inside the largest strict")
    print("    safe component without crossing a danger boundary.")
    print("    Boundary-defect rows have eps<=0 and must retain taut endpoint atoms.")
    print()


def print_boundary_atoms(records: list[KernelRecord]) -> None:
    print("[2] Boundary-defect atoms")
    for rec in records:
        if rec.state != "boundary-defect":
            continue
        audit = s155.audit_row(rec.row.name, rec.row.family, rec.row.speeds)
        print(
            f"  {rec.row.name}: taut_vertices={rec.taut_count} "
            f"zero_sum_pairs={rec.zero_sum_pairs}"
        )
        for taut in audit.taut_vertices:
            left = s155.owners_text(taut.left_end)
            right = s155.owners_text(taut.right_start)
            print(
                f"    t={fmt(taut.t):>5s} left={left:20s} "
                f"right={right:20s} pair_mod14={list(taut.pair_sum_mod14)}"
            )
    print()


def print_single_scan(scan: dict[str, object]) -> None:
    print("[3] Single-swap kernel stability sanity scan")
    print(f"  add_limit={scan['limit']} rows={scan['total']}")
    print(f"  counts={dict(sorted(scan['counts'].items()))}")
    min_positive = scan["min_positive"]
    if min_positive:
        mu, label, row = min_positive
        print(f"  min_positive_safe_mu={fmt(mu)} at {label} row={row}")
    zero_rows = scan["zero_rows"]
    if zero_rows:
        print("  zero-open rows:")
        for label, row in zero_rows:
            print(f"    {label}: row={row}")
    print()


def print_tournament() -> None:
    fp = tournament_fingerprint()
    print("[4] Tournament Analysis")
    print("  vertex sets challenged:")
    print("    runners, arcs, kernels, safe components, boundary events,")
    print("    endpoint-owner pairs, packet fibers, smoothing defects, proof obligations.")
    print("  chosen vertices:")
    print("    proof carriers in the kernel-deformation step.")
    print("  pairwise observable:")
    print("    retention of strict witness, boundary equality, packet labels,")
    print("    kernel auditability, formal interval path, and scalar-decoy resistance.")
    print("  switch/gauge:")
    print("    componentwise majority of the six retention scores; ties follow the")
    print("    declared Hamiltonian path.")
    print(
        f"  fingerprint: score_hist={fp['score_hist']} "
        f"directed_3cycles={fp['directed_3cycles']} scc_sizes={fp['scc_sizes']}"
    )
    print("  Hamiltonian path:")
    for idx, vertex in enumerate(fp["hamiltonian_path"], start=1):
        print(f"    {idx}. {vertex}")
    print()


def print_readout() -> None:
    print("[5] Readout")
    print("  The kernel-homotopy rule is a proof-interface guardrail, not a proof of")
    print("  LRC14.  Open components are stable under small kernel support changes and")
    print("  can feed HYP-2981-style Fejer interval certificates.  AP/GW have no such")
    print("  support radius; they are boundary-defect atoms and must remain labelled.")
    print("  A future smoothing proof must therefore state exactly which labels survive")
    print("  the homotopy and which boundary defects are emitted.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=160)
    args = parser.parse_args()

    print("S164 LRC14 KERNEL HOMOTOPY BOUNDARY-DEFECT LEDGER")
    print("=" * 78)
    print("[0] Scope")
    print("  Exact regular-open safe components provide kernel-stability radii.")
    print("  Isolated threshold-safe points provide boundary-defect atoms.")
    print("  Quotients may forget raw runners only after this predicate is retained.")
    print()

    records = [kernel_record(row) for row in named_rows()]
    print_records(records)
    print_boundary_atoms(records)
    print_single_scan(scan_single_swaps(args.single_limit))
    print_tournament()
    print_readout()


if __name__ == "__main__":
    main()
