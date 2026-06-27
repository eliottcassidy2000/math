#!/usr/bin/env python3
"""S168: residual-section characterization and packet-grid verifier.

This is an executable proof-interface audit, not an LRC14 proof.

The previous cocycle passes say a quotient is theorem-safe only when its lost
coordinate is exact, dual-annihilated, descended, boundary-equality, or emitted
as a named residual.  This script makes that statement concrete for the
HYP-2963 packet bank by assigning each packet to a residual section and to a
HYP-2989 Haar-product grid class.

The main verification is conservative:

* AP/GW boundary packets must be same-tile boundary atoms.
* Non-AP/GW hard packets must land in owner-strip, cross-handoff, or nested
  refinement grid exits.
* Any hard, zero-open, non-AP/GW packet not covered by those exits is reported
  as a candidate F7 harmonic residual section.

Tournament Analysis uses residual sections and grid exits as vertices, not
runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import argparse
import os
import sys


ROOT = Path(__file__).resolve().parents[1]
THRESHOLD = Fraction(1, 14)
LOW_227 = Fraction(2, 27)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module {name} from {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


packets_mod = load_module(
    "s168_hyp2963_packets",
    ROOT / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)
tile_mod = load_module(
    "s168_hyp2989_tiles",
    ROOT / "04-computation" / "lrc14_haar_product_tile_discrepancy_codex_s165.py",
)


@dataclass(frozen=True)
class SectionRecord:
    section: str
    route: str
    grid_class: str
    cohomology_role: str
    exit_status: str
    residual_if_failed: str
    zeta_readout: str
    vector: tuple[int, ...]


@dataclass(frozen=True)
class PacketGridRecord:
    name: str
    source_family: str
    m_value: Fraction
    q_threshold: int
    strict_safe_mu: Fraction
    packet_route: str
    packet_family: str
    branch: str
    transfer: str
    section: str
    grid_class: str
    cohomology_role: str
    exit_status: str
    zeta_readout: str


SECTION_ORDER = [
    "direct_q_witness_section",
    "ap_gw_boundary_section",
    "unit_petal_c27_strip_section",
    "open_haar_witness_section",
    "covering_boundary_moment_section",
    "k33_state_lift_section",
    "candidate_f7_harmonic_section",
    "raw_scalar_shadow",
]


SECTION_RULES: dict[str, SectionRecord] = {
    "direct_q_witness_section": SectionRecord(
        section="direct_q_witness_section",
        route="Q-WITNESS",
        grid_class="orthogonal_zero",
        cohomology_role="exact 0-cochain witness",
        exit_status="discharged by q<=13 strict witness",
        residual_if_failed="none; failure would contradict the exact q-threshold computation",
        zeta_readout="mixed Haar mode irrelevant after denominator witness",
        vector=(3, 3, 0, 1, 0, 0, 0, 3, 3),
    ),
    "ap_gw_boundary_section": SectionRecord(
        section="ap_gw_boundary_section",
        route="BOUNDARY-AP-GW",
        grid_class="same_tile_indicator",
        cohomology_role="declared boundary cohomology",
        exit_status="allowed equality atom",
        residual_if_failed="boundary-owner/C27 labels must distinguish the atom",
        zeta_readout="zeta vanishes on same-tile boundary atom",
        vector=(3, 3, 3, 3, 0, 0, 2, 3, 3),
    ),
    "unit_petal_c27_strip_section": SectionRecord(
        section="unit_petal_c27_strip_section",
        route="BOUNDARY-PETAL-SPORADIC",
        grid_class="horizontal_owner_strip",
        cohomology_role="owner-strip coboundary with C27 label",
        exit_status="discharged by unit-petal/C27 strip",
        residual_if_failed="C27 transfer must be retained before scalarizing",
        zeta_readout="one-axis owner strip exposes the mixed mode",
        vector=(3, 3, 2, 3, 0, 1, 1, 3, 3),
    ),
    "open_haar_witness_section": SectionRecord(
        section="open_haar_witness_section",
        route="COVERING-MOMENT",
        grid_class="vertical_owner_strip",
        cohomology_role="positive regular-open section",
        exit_status="discharged by positive Haar-open witness",
        residual_if_failed="exact interval certificate must name the owner strip",
        zeta_readout="one-axis owner strip carries a nonzero witness address",
        vector=(3, 2, 2, 3, 0, 1, 0, 3, 3),
    ),
    "covering_boundary_moment_section": SectionRecord(
        section="covering_boundary_moment_section",
        route="COVERING-MOMENT",
        grid_class="nested_refinement",
        cohomology_role="descended boundary-moment section",
        exit_status="routed to smaller exact-period / moment chart",
        residual_if_failed="covering rows must emit boundary-moment or F7 debt",
        zeta_readout="same-direction nesting keeps the residual in a smaller packet",
        vector=(3, 2, 1, 3, 1, 3, 1, 2, 3),
    ),
    "k33_state_lift_section": SectionRecord(
        section="k33_state_lift_section",
        route="K33-STATE-LIFT",
        grid_class="cross_handoff",
        cohomology_role="named THM-572 state-lift residual",
        exit_status="not discharged here; routed to state-lift theorem",
        residual_if_failed="this is the named K33/THM-572 residual debt",
        zeta_readout="opposite-coordinate nesting creates a cross handoff",
        vector=(3, 3, 1, 3, 3, 1, 3, 1, 3),
    ),
    "candidate_f7_harmonic_section": SectionRecord(
        section="candidate_f7_harmonic_section",
        route="SOURCE-SPECTRUM-UNKNOWN",
        grid_class="cross_handoff",
        cohomology_role="candidate non-exact harmonic residual",
        exit_status="unresolved; promote immediately if present",
        residual_if_failed="genuine F7 harmonic/source-spectrum section",
        zeta_readout="unpaired mixed mode survives all known exits",
        vector=(3, 3, 1, 2, 2, 2, 3, 0, 3),
    ),
    "raw_scalar_shadow": SectionRecord(
        section="raw_scalar_shadow",
        route="RAW",
        grid_class="orthogonal_zero",
        cohomology_role="negative control",
        exit_status="not a proof exit",
        residual_if_failed="scalar quotient forgot the packet grid",
        zeta_readout="zeta coordinate erased",
        vector=(1, 1, 0, 0, 0, 0, 0, 0, 0),
    ),
}


def fmt(frac: Fraction) -> str:
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def section_for_packet(packet) -> SectionRecord:
    if packet.route == "Q-WITNESS":
        return SECTION_RULES["direct_q_witness_section"]
    if packet.route == "BOUNDARY-AP-GW":
        return SECTION_RULES["ap_gw_boundary_section"]
    if packet.route == "BOUNDARY-PETAL-SPORADIC":
        return SECTION_RULES["unit_petal_c27_strip_section"]
    if packet.route == "K33-STATE-LIFT":
        return SECTION_RULES["k33_state_lift_section"]
    if packet.route == "SOURCE-SPECTRUM-UNKNOWN":
        return SECTION_RULES["candidate_f7_harmonic_section"]
    if packet.route == "COUNTEREXAMPLE":
        return SECTION_RULES["candidate_f7_harmonic_section"]
    if packet.route == "COVERING-MOMENT":
        if packet.family == "covering boundary-moment family":
            return SECTION_RULES["covering_boundary_moment_section"]
        return SECTION_RULES["open_haar_witness_section"]
    # Shell/magnitude-liar routes do not appear in the default HYP-2963 bank
    # after q-witness/open-Haar routing, but they are scalar negative controls.
    return SECTION_RULES["raw_scalar_shadow"]


def build_records(args: argparse.Namespace) -> list[PacketGridRecord]:
    rows = packets_mod.build_bank(
        args.single_limit,
        args.two_swap_limit,
        args.alias_depth,
        args.lcm_tail_max,
    )
    packets = packets_mod.compute_packets(rows, workers=args.workers)
    records: list[PacketGridRecord] = []
    for packet in packets:
        section = section_for_packet(packet)
        records.append(
            PacketGridRecord(
                name=packet.name,
                source_family=packet.source_family,
                m_value=packet.M,
                q_threshold=packet.q_threshold,
                strict_safe_mu=packet.strict_safe_mu,
                packet_route=packet.route,
                packet_family=packet.family,
                branch=packet.branch,
                transfer=packet.transfer,
                section=section.section,
                grid_class=section.grid_class,
                cohomology_role=section.cohomology_role,
                exit_status=section.exit_status,
                zeta_readout=section.zeta_readout,
            )
        )
    return records


def validate_records(records: list[PacketGridRecord]) -> list[str]:
    errors: list[str] = []
    valid_grid = set(tile_mod.CLASS_ORDER)
    for rec in records:
        if rec.grid_class not in valid_grid:
            errors.append(f"{rec.name}: invalid grid class {rec.grid_class}")
        if rec.packet_route == "BOUNDARY-AP-GW" and rec.grid_class != "same_tile_indicator":
            errors.append(f"{rec.name}: AP/GW boundary not mapped to same_tile_indicator")
        if rec.packet_route != "BOUNDARY-AP-GW" and rec.q_threshold >= 14:
            if rec.grid_class == "same_tile_indicator":
                errors.append(f"{rec.name}: non-AP/GW hard packet collapsed to same tile")
        if (
            rec.q_threshold >= 14
            and rec.packet_route != "BOUNDARY-AP-GW"
            and rec.strict_safe_mu == 0
            and rec.section != "candidate_f7_harmonic_section"
        ):
            errors.append(f"{rec.name}: zero-open hard packet not promoted to F7 candidate")
    return errors


def tile_census() -> tuple[Counter[str], dict[str, Counter[int]]]:
    return tile_mod.ordered_product_census(3)


def orient_section(a: str, b: str) -> tuple[str, str]:
    weights = (5, 4, 4, 5, 4, 3, 4, 5, 6)
    va = SECTION_RULES[a].vector
    vb = SECTION_RULES[b].vector
    sa = sum(w * x for w, x in zip(weights, va))
    sb = sum(w * x for w, x in zip(weights, vb))
    if sa > sb:
        return a, b
    if sb > sa:
        return b, a
    return (a, b) if SECTION_ORDER.index(a) < SECTION_ORDER.index(b) else (b, a)


def tournament_edges(vertices: list[str]) -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    for a, b in combinations(vertices, 2):
        edges.add(orient_section(a, b))
    return edges


def has_edge(edges: set[tuple[str, str]], a: str, b: str) -> bool:
    return (a, b) in edges


def directed_3cycles(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    total = 0
    for a, b, c in combinations(vertices, 3):
        if (
            has_edge(edges, a, b)
            and has_edge(edges, b, c)
            and has_edge(edges, c, a)
        ) or (
            has_edge(edges, a, c)
            and has_edge(edges, c, b)
            and has_edge(edges, b, a)
        ):
            total += 1
    return total


def strongly_connected_components(vertices: list[str], edges: set[tuple[str, str]]) -> list[list[str]]:
    graph = {v: [] for v in vertices}
    rev = {v: [] for v in vertices}
    for a, b in edges:
        graph[a].append(b)
        rev[b].append(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    comps: list[list[str]] = []
    seen.clear()

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(comp)
    return comps


def hamiltonian_path_count(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    index = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp: dict[tuple[int, int], int] = {}
    for v in vertices:
        dp[(1 << index[v], index[v])] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if not ways:
                continue
            a = vertices[last]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                b = vertices[nxt]
                if has_edge(edges, a, b):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + ways
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def print_header(args: argparse.Namespace) -> None:
    print("LRC14 residual-section packet-grid verifier")
    print("=" * 78)
    print("namespace: HYP-2996 / T1080 / codex-2026-06-24-S168")
    print(
        f"packet bank: single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}, workers={args.workers}"
    )
    print()


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, fixed circle sections, section boundaries, wall-crossing")
    print("    events, residues, cover arcs, Fourier modes, Haar rectangles, matroid")
    print("    circuits, proof obligations, residual sections, and packet-grid exits.")
    print("  chosen vertices:")
    print("    residual sections and Haar-product grid exits over HYP-2963 packets.")
    print("  preserved LRC predicate:")
    print("    exact M(S)>=1/14, q-threshold, open-vs-boundary status, C27/K33 owner")
    print("    labels, covering route, and the named F7/THM-572 residual predicate.")
    print("  destroyed information:")
    print("    raw runner identity and scalar M buckets are deliberately discarded only")
    print("    after their lost coordinates have named section records.")
    print()


def print_section_schema() -> None:
    print("[1] Residual section schema")
    print(f"  {'section':34s} {'grid':24s} {'role':34s} exit")
    for name in SECTION_ORDER:
        rule = SECTION_RULES[name]
        print(f"  {name:34s} {rule.grid_class:24s} {rule.cohomology_role:34s} {rule.exit_status}")
    print()


def print_grid_baseline() -> None:
    counts, signs = tile_census()
    total = sum(counts.values())
    print("[2] HYP-2989 Haar-product grid baseline")
    print(f"  depth<=3 rectangles={(2 ** 4 - 1) ** 2}; ordered products={total}")
    print(f"  {'grid class':24s} {'tile count':>10s} {'signs':>22s}")
    for name in tile_mod.CLASS_ORDER:
        sign_text = ", ".join(f"{k}:{v}" for k, v in sorted(signs[name].items()))
        print(f"  {name:24s} {counts[name]:10d} {sign_text:>22s}")
    print()


def print_packet_verification(records: list[PacketGridRecord]) -> None:
    errors = validate_records(records)
    by_section = Counter(r.section for r in records)
    by_grid = Counter(r.grid_class for r in records)
    by_route = Counter(r.packet_route for r in records)
    hard = [r for r in records if r.q_threshold >= 14]
    hard_non_boundary = [r for r in hard if r.packet_route != "BOUNDARY-AP-GW"]
    exposed_hard = [
        r
        for r in hard_non_boundary
        if r.grid_class in {"vertical_owner_strip", "horizontal_owner_strip", "cross_handoff", "nested_refinement"}
    ]
    f7 = [r for r in records if r.section == "candidate_f7_harmonic_section"]
    zero_open_non_boundary = [
        r
        for r in records
        if r.q_threshold >= 14 and r.packet_route != "BOUNDARY-AP-GW" and r.strict_safe_mu == 0
    ]

    print("[3] HYP-2963 packet-grid verification")
    print(f"  audited packets={len(records)}")
    print(f"  hard q>=14 packets={len(hard)}")
    print(f"  hard non-AP/GW packets={len(hard_non_boundary)}")
    print(f"  hard non-AP/GW packets with owner/cross/nested grid exit={len(exposed_hard)}")
    print(f"  zero-open hard non-AP/GW packets={len(zero_open_non_boundary)}")
    print(f"  candidate F7 harmonic sections={len(f7)}")
    print(f"  validation errors={len(errors)}")
    if errors:
        for err in errors[:12]:
            print(f"    ERROR {err}")
        if len(errors) > 12:
            print(f"    ... {len(errors) - 12} more errors")
    print("  route counts:")
    for route, count in sorted(by_route.items()):
        print(f"    {route:24s} {count}")
    print("  section counts:")
    for section in SECTION_ORDER:
        if by_section[section]:
            print(f"    {section:34s} {by_section[section]}")
    print("  packet grid counts:")
    for grid in tile_mod.CLASS_ORDER:
        if by_grid[grid]:
            print(f"    {grid:24s} {by_grid[grid]}")
    print()


def representative_records(records: list[PacketGridRecord]) -> list[PacketGridRecord]:
    selected: list[PacketGridRecord] = []
    wanted_names = {
        "AP",
        "GW 12->24",
        "near/K33 12->36",
        "petal 10->20",
        "petal 13->26",
        "P10+GW",
        "P10+K33",
        "covering comb 12->84",
        "covering repair drop13 add182",
        "scale-separated AP tail 200",
    }
    for rec in records:
        if rec.name in wanted_names:
            selected.append(rec)
    for section in SECTION_ORDER:
        section_rows = [r for r in records if r.section == section and r not in selected]
        selected.extend(section_rows[:2])
    selected = sorted(
        selected,
        key=lambda r: (SECTION_ORDER.index(r.section), r.m_value, r.name),
    )
    return selected


def print_representatives(records: list[PacketGridRecord]) -> None:
    print("[4] Representative packet-grid records")
    print(f"  {'name':34s} {'M':>8s} {'q':>3s} {'mu':>9s} {'section':34s} grid")
    for rec in representative_records(records)[:32]:
        print(
            f"  {rec.name[:34]:34s} {fmt(rec.m_value):>8s} {rec.q_threshold:3d} "
            f"{fmt(rec.strict_safe_mu):>9s} {rec.section:34s} {rec.grid_class}"
        )
        print(f"      branch={rec.branch}; transfer={rec.transfer}")
        print(f"      exit={rec.exit_status}; zeta={rec.zeta_readout}")
    print()


def print_tournament_analysis(records: list[PacketGridRecord]) -> None:
    vertices = [v for v in SECTION_ORDER if v == "raw_scalar_shadow" or any(r.section == v for r in records)]
    edges = tournament_edges(vertices)
    out = Counter()
    for a, _ in edges:
        out[a] += 1
    score_hist = Counter(out[v] for v in vertices)
    print("[5] Tournament Analysis")
    print("  vertices: residual sections / packet-grid exits, not runners")
    print("  pairwise observable:")
    print("    retained LRC predicate, exact scale, boundary data, Haar grid class,")
    print("    state-lift visibility, moment-route visibility, named residual status,")
    print("    proof exit, and anti-scalarization guard.")
    print("  switch:")
    print("    orient toward the stronger weighted retention vector.")
    print("  tie Hamiltonian path:")
    print("    " + " -> ".join(vertices))
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={directed_3cycles(vertices, edges)}")
    print(f"  SCC_sizes={[len(c) for c in strongly_connected_components(vertices, edges)]}")
    print(f"  Hamiltonian_path_count={hamiltonian_path_count(vertices, edges)}")
    print(f"  {'section':34s} {'out':>3s} vector")
    for v in vertices:
        print(f"  {v:34s} {out[v]:3d} {SECTION_RULES[v].vector}")
    print()


def print_theorem_readout(records: list[PacketGridRecord]) -> None:
    f7 = [r for r in records if r.section == "candidate_f7_harmonic_section"]
    hard_non_boundary = [
        r
        for r in records
        if r.q_threshold >= 14 and r.packet_route != "BOUNDARY-AP-GW"
    ]
    print("[6] Theorem-facing readout")
    print("  Residual section characterization:")
    print("    F0/q-witness packets are exact 0-cochain exits.")
    print("    AP/GW are boundary cohomology on same-tile atoms.")
    print("    Unit-petal/C27 rows are owner-strip coboundaries with C27 labels retained.")
    print("    K33 rows are cross-handoff residuals routed to THM-572 state lift.")
    print("    Covering rows either expose a vertical owner strip or descend to a")
    print("    boundary-moment chart by nested refinement.")
    print("    F7 is now a testable section: hard, zero-open, non-AP/GW, and not")
    print("    discharged by owner-strip/cross/nested grid exits.")
    print()
    print("  Verification status in this HYP-2963 bank:")
    print(f"    hard non-AP/GW residual packets={len(hard_non_boundary)}")
    print(f"    candidate F7 sections={len(f7)}")
    if not f7:
        print("    no F7 section appears in the audited bank.")
    print("  This remains a bounded proof-interface audit, not a global proof.")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    args = parser.parse_args()

    print_header(args)
    print_assumption_challenge()
    print_section_schema()
    print_grid_baseline()
    records = build_records(args)
    print_packet_verification(records)
    print_representatives(records)
    print_tournament_analysis(records)
    print_theorem_readout(records)


if __name__ == "__main__":
    main()
