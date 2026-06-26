#!/usr/bin/env python3
"""Status-topology bridge for LRC14.

HYP-3024 shows that a coarse Erdos-Turan plus Henselian-unit gate has no
mixed boundary/open fibers on the full HYP-2963 packet bank, while HYP-3025
shows that AP/Goddyn-Wong boundary rows are closed arc-Cech cover cycles.

This script audits the overlap.  It focuses on the only residue-terminal
fibers that still mix boundary and open packets, then on the route-mixed
fibers that survive the coarse ET+unit gate.

Tournament Analysis declaration:
  vertices: proof gates, not runners;
  pairwise observable: predicate preservation, topology exactness, arithmetic
    compression, route scheduling, quotient-defect visibility, and proof cost;
  switch/gauge: majority comparison of observable vectors;
  tie Hamiltonian path: arc boundary gate > coarse ET+unit status gate >
    magnitude route splitter > barcode packet scheduler > raw residue word.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


fz = load_module(
    "lrc14_fiber_zipper_convergence_s194",
    REPO / "04-computation" / "lrc14_fiber_zipper_convergence_codex_s188.py",
)
ac = load_module(
    "lrc14_arc_cech_nerve_s194",
    REPO / "04-computation" / "lrc14_arc_cech_nerve_carrier_codex_s174.py",
)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


@dataclass(frozen=True)
class TopologySig:
    closed_b0: int
    closed_b1: int
    open_b0: int
    open_b1: int
    safe_topes: int
    safe_measure: Fraction
    quotient_defect: int
    owner_sums_mod14: tuple[int, ...]

    @property
    def compact(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.closed_b0,
            self.closed_b1,
            self.open_b0,
            self.open_b1,
            self.safe_topes,
            self.quotient_defect,
        )

    @property
    def apgw_boundary_cycle(self) -> bool:
        return (
            self.closed_b0 == 1
            and self.closed_b1 == 1
            and self.open_b0 == 6
            and self.open_b1 == 0
            and self.safe_topes == 0
            and self.owner_sums_mod14 == (0, 0, 0, 0, 0, 0)
        )


def topology_sig(packet) -> TopologySig:
    rec = ac.row_record(packet.name, tuple(packet.speeds))
    return TopologySig(
        rec.closed_arc_betti.beta0,
        rec.closed_arc_betti.beta1,
        rec.arc_betti.beta0,
        rec.arc_betti.beta1,
        rec.audit.safe_topes,
        rec.audit.safe_measure,
        rec.runner_quotient_defect,
        rec.boundary_pair_sums,
    )


def group_by_split(packets: list, split_name: str) -> dict[object, list]:
    split = next(s for s in fz.SPLITS if s.name == split_name)
    groups: dict[object, list] = defaultdict(list)
    for packet in packets:
        groups[split.key_func(packet)].append(packet)
    return groups


def mixed_status_groups(groups: dict[object, list]) -> list[list]:
    return [rows for rows in groups.values() if len({fz.status(p) for p in rows}) > 1]


def mixed_route_groups(groups: dict[object, list]) -> list[list]:
    return [rows for rows in groups.values() if len({p.route for p in rows}) > 1]


def short_row(packet, sig: TopologySig) -> str:
    return (
        f"{packet.name:<34} status={fz.status(packet):<8} route={packet.route:<16} "
        f"M={fmt(packet.M):>7} mu={fmt(packet.strict_safe_mu):>12} "
        f"closed=({sig.closed_b0},{sig.closed_b1}) open=({sig.open_b0},{sig.open_b1}) "
        f"topes={sig.safe_topes:<2d} qdef={sig.quotient_defect:<2d} "
        f"apgw_cycle={sig.apgw_boundary_cycle}"
    )


def audit_residue_status(packets: list) -> dict[str, object]:
    groups = group_by_split(packets, "residue_terminal_fiber")
    mixed = sorted(mixed_status_groups(groups), key=lambda rows: (-len(rows), rows[0].name))

    boundary_good = 0
    boundary_bad = 0
    open_beta1 = 0
    open_min_topes = None
    topology_classes = Counter()
    lines: list[str] = []

    for index, rows in enumerate(mixed):
        rows = sorted(rows, key=lambda p: (fz.status(p), p.route, p.name))
        lines.append(
            f"  residue_mixed_status_fiber[{index}] size={len(rows)} "
            f"status={dict(Counter(fz.status(p) for p in rows))} "
            f"routes={dict(Counter(p.route for p in rows))}"
        )
        for packet in rows:
            sig = topology_sig(packet)
            topology_classes[sig.compact] += 1
            if fz.status(packet) == "boundary":
                if sig.apgw_boundary_cycle:
                    boundary_good += 1
                else:
                    boundary_bad += 1
            else:
                if sig.closed_b1:
                    open_beta1 += 1
                open_min_topes = sig.safe_topes if open_min_topes is None else min(open_min_topes, sig.safe_topes)
            lines.append("    " + short_row(packet, sig))

    return {
        "groups": len(groups),
        "mixed_status": len(mixed),
        "mixed_route": len(mixed_route_groups(groups)),
        "max_route": max((len(rows) for rows in mixed_route_groups(groups)), default=0),
        "boundary_good": boundary_good,
        "boundary_bad": boundary_bad,
        "open_beta1": open_beta1,
        "open_min_topes": open_min_topes,
        "topology_classes": topology_classes,
        "lines": lines,
    }


def audit_coarse_route(packets: list) -> dict[str, object]:
    groups = group_by_split(packets, "coarse_et_unit_gate")
    mixed = sorted(mixed_route_groups(groups), key=lambda rows: (-len(rows), rows[0].name))

    status_counter = Counter()
    route_counter = Counter()
    topology_classes = Counter()
    open_beta1 = 0
    open_min_topes = None
    fiber_lines: list[str] = []

    for index, rows in enumerate(mixed):
        rows = sorted(rows, key=lambda p: (p.route, p.M, p.name))
        classes = Counter()
        for packet in rows:
            sig = topology_sig(packet)
            status_counter[fz.status(packet)] += 1
            route_counter[packet.route] += 1
            topology_classes[sig.compact] += 1
            classes[sig.compact] += 1
            if fz.status(packet) == "open":
                if sig.closed_b1:
                    open_beta1 += 1
                open_min_topes = sig.safe_topes if open_min_topes is None else min(open_min_topes, sig.safe_topes)
        fiber_lines.append(
            f"  coarse_route_mixed_fiber[{index}] size={len(rows)} "
            f"status={dict(Counter(fz.status(p) for p in rows))} "
            f"routes={dict(Counter(p.route for p in rows))} "
            f"topology_classes={dict(classes)}"
        )
        for packet in rows:
            fiber_lines.append("    " + short_row(packet, topology_sig(packet)))

    return {
        "groups": len(groups),
        "mixed_status": len(mixed_status_groups(groups)),
        "mixed_route": len(mixed),
        "max_route": max((len(rows) for rows in mixed), default=0),
        "status_counter": status_counter,
        "route_counter": route_counter,
        "topology_classes": topology_classes,
        "open_beta1": open_beta1,
        "open_min_topes": open_min_topes,
        "lines": fiber_lines,
    }


@dataclass(frozen=True)
class Gate:
    name: str
    predicate: int
    topology: int
    arithmetic: int
    route: int
    quotient_defect: int
    proof_cost: int

    @property
    def values(self) -> tuple[int, ...]:
        return (
            self.predicate,
            self.topology,
            self.arithmetic,
            self.route,
            self.quotient_defect,
            8 - self.proof_cost,
        )


GATES = (
    Gate("arc_boundary_cycle_gate", 5, 5, 1, 2, 5, 2),
    Gate("coarse_et_unit_status_gate", 5, 2, 5, 3, 2, 3),
    Gate("magnitude_route_splitter", 4, 1, 4, 5, 2, 4),
    Gate("barcode_packet_scheduler", 4, 3, 2, 5, 3, 5),
    Gate("raw_residue_terminal_word", 1, 0, 1, 1, 0, 1),
)


def gate_tournament() -> dict[str, object]:
    edges: set[tuple[str, str]] = set()
    names = [gate.name for gate in GATES]
    order = {name: i for i, name in enumerate(names)}
    score = Counter()
    for a, b in combinations(GATES, 2):
        aw = sum(x > y for x, y in zip(a.values, b.values))
        bw = sum(x < y for x, y in zip(a.values, b.values))
        if aw > bw or (aw == bw and order[a.name] < order[b.name]):
            edges.add((a.name, b.name))
            score[a.name] += 1
        else:
            edges.add((b.name, a.name))
            score[b.name] += 1

    c3 = 0
    for a, b, c in combinations(names, 3):
        out = {name: 0 for name in (a, b, c)}
        for x, y in combinations((a, b, c), 2):
            if (x, y) in edges:
                out[x] += 1
            else:
                out[y] += 1
        if sorted(out.values()) == [1, 1, 1]:
            c3 += 1

    n = len(names)
    edge_idx = {(names.index(a), names.index(b)) for a, b in edges}
    hp = 0
    for perm in permutations(range(n)):
        if all((perm[i], perm[i + 1]) in edge_idx for i in range(n - 1)):
            hp += 1

    ranking = sorted(names, key=lambda name: (-score[name], order[name]))
    return {
        "score_hist": dict(sorted(Counter(score[name] for name in names).items())),
        "directed_3cycles": c3,
        "hamiltonian_path_count": hp,
        "score_order": ranking,
    }


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, individual danger arcs, endpoint cells, boundary")
    print("    cocircuits, residues, Fourier/ET clocks, p-adic unit roots,")
    print("    magnitude fibers, barcode packets, and proof obligations.")
    print("  chosen vertices:")
    print("    proof gates: arc boundary cycle, coarse ET+unit status,")
    print("    magnitude route split, barcode packet scheduler, residue word.")
    print("  preserved LRC predicate:")
    print("    boundary/open status for the question whether a strict safe")
    print("    interval exists at threshold 1/14.")
    print("  destroyed information:")
    print("    raw runner identity, exact route labels, and full endpoint geometry")
    print("    unless a later gate reattaches them.")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    rows = fz.lp.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    packets = fz.lp.compute_packets(rows, args.workers)

    residue = audit_residue_status(packets)
    coarse = audit_coarse_route(packets)

    print("=== LRC14 status-topology gate bridge S194 ===")
    print(
        "bank=HYP-2963 default rows "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"packets={len(packets)}")
    print()
    print_assumption_challenge()

    print("[1] Residue-terminal boundary/open collisions")
    print(
        f"  fibers={residue['groups']} mixed_status={residue['mixed_status']} "
        f"mixed_route={residue['mixed_route']} max_route_mixed={residue['max_route']}"
    )
    print(
        f"  boundary_apgw_cycle={residue['boundary_good']} "
        f"boundary_noncycle={residue['boundary_bad']} "
        f"open_with_closed_beta1={residue['open_beta1']} "
        f"open_min_safe_topes={residue['open_min_topes']}"
    )
    print("  topology_classes=" + str(dict(residue["topology_classes"])))
    for line in residue["lines"]:
        print(line)
    print()

    print("[2] Coarse ET+unit route-mixed fibers")
    print(
        f"  fibers={coarse['groups']} mixed_status={coarse['mixed_status']} "
        f"mixed_route={coarse['mixed_route']} max_route_mixed={coarse['max_route']}"
    )
    print(f"  aggregate_status={dict(coarse['status_counter'])}")
    print(f"  aggregate_routes={dict(coarse['route_counter'])}")
    print(
        f"  open_with_closed_beta1={coarse['open_beta1']} "
        f"open_min_safe_topes={coarse['open_min_topes']}"
    )
    print("  topology_classes=" + str(dict(coarse["topology_classes"])))
    for line in coarse["lines"]:
        print(line)
    print()

    print("[3] Tournament Analysis over proof gates")
    fp = gate_tournament()
    print("  vertices_are=proof gates, not runners")
    print("  score_hist=" + str(fp["score_hist"]))
    print("  directed_3cycles=" + str(fp["directed_3cycles"]))
    print("  hamiltonian_path_count=" + str(fp["hamiltonian_path_count"]))
    print("  score_order=" + " > ".join(fp["score_order"]))
    print()

    print("[4] Proof readout")
    print("  1. The residue-terminal quotient has exactly two boundary/open")
    print("     collisions on the tested full bank.  They are AP and GW against")
    print("     open impostors.")
    print("  2. In those collisions, the boundary rows are exactly the AP/GW")
    print("     arc-Cech cycle: closed arc beta=(1,1), open arc beta=(6,0),")
    print("     zero safe topes, and six owner-current sums 0 mod 14.")
    print("  3. Every open cohabitant in those two fibers has closed arc beta1=0")
    print("     and at least four safe topes.  This gives a topology-first")
    print("     separator for equality atoms before arithmetic route labels.")
    print("  4. The coarse ET+unit gate has no boundary/open mixing.  Its 15")
    print("     remaining route-mixed fibers contain 38 packets, all open, all")
    print("     closed beta1=0.  Route mixing here is harmless for the direct")
    print("     LRC predicate and can be scheduled later by magnitude/covering")
    print("     certificates.")
    print("  5. Candidate theorem: primitive zero-open packets must carry the")
    print("     AP/GW arc-boundary cycle, or else the arc topology emits positive")
    print("     safe topes and the coarse ET+unit gate is allowed to forget the")
    print("     exact route while preserving status.")


if __name__ == "__main__":
    main()
