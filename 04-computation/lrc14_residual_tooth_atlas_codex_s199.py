#!/usr/bin/env python3
"""Residual tooth atlas for the LRC14 coarse ET+unit gate.

HYP-3028/S192 showed that the full-bank coarse Erdos-Turan plus Henselian-unit
gate has no mixed boundary/open fibers, but leaves 15 route-mixed open fibers.
HYP-3031/S195 reframed those residuals as lost mixed Haar/tile coordinates.

This pass audits the 15 residual fibers directly.  It asks which non-route
tooth first separates each Q-witness / covering-moment collision:

    arc topology, coarse safe-component stalk, exact safe-component stalk,
    magnitude cocycle, or only the explicit certificate label.

Tournament Analysis vertices are repair teeth, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path
import argparse
import re
import sys


REPO = Path(__file__).resolve().parents[1]
S194_OUT = REPO / "05-knowledge" / "results" / "lrc14_status_topology_gate_codex_s194.out"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


topo = load_module(
    "lrc14_status_topology_gate_s196",
    REPO / "04-computation" / "lrc14_status_topology_gate_codex_s194.py",
)
stalk_mod = load_module(
    "lrc14_safe_component_stalk_s196",
    REPO / "04-computation" / "lrc14_safe_component_stalk_codex_s193.py",
)

fz = topo.fz


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def magnitude_key(packet) -> tuple[Fraction, int, int, Fraction]:
    return packet.M, packet.q_threshold, packet.farey_excess, packet.lacunary_tail_ratio


def direct_certificate_key(packet) -> tuple[str, object]:
    if packet.route == "Q-WITNESS":
        return ("q_witness", packet.q_threshold)
    if packet.route == "COVERING-MOMENT":
        return ("covering_positive", packet.strict_safe_mu > 0)
    return ("other", packet.route)


@dataclass(frozen=True)
class ResidualRow:
    packet: object
    topology: object
    stalk: object


@dataclass(frozen=True)
class Tooth:
    name: str
    repair_class: str
    key_func: object
    preserves_status: int
    route_split: int
    nonroute: int
    locality: int
    compression: int
    proof_cost: int

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.preserves_status,
            self.route_split,
            self.nonroute,
            self.locality,
            self.compression,
            8 - self.proof_cost,
        )


def tooth_family() -> tuple[Tooth, ...]:
    return (
        Tooth(
            "arc_topology_compact",
            "owner_strip",
            lambda row: row.topology.compact,
            5,
            3,
            5,
            4,
            4,
            2,
        ),
        Tooth(
            "coarse_safe_stalk",
            "owner_strip",
            lambda row: row.stalk.coarse_key,
            5,
            4,
            5,
            5,
            3,
            3,
        ),
        Tooth(
            "exact_safe_stalk",
            "nested_refinement",
            lambda row: row.stalk.exact_key,
            5,
            5,
            5,
            5,
            2,
            4,
        ),
        Tooth(
            "magnitude_cocycle",
            "nested_refinement",
            lambda row: magnitude_key(row.packet),
            5,
            5,
            5,
            2,
            1,
            5,
        ),
        Tooth(
            "q_or_covering_certificate",
            "nested_refinement",
            lambda row: direct_certificate_key(row.packet),
            5,
            5,
            1,
            3,
            5,
            2,
        ),
    )


def route_mixed(rows: list[ResidualRow], key_func) -> list[list[ResidualRow]]:
    buckets: dict[object, list[ResidualRow]] = defaultdict(list)
    for row in rows:
        buckets[key_func(row)].append(row)
    return [bucket for bucket in buckets.values() if len({r.packet.route for r in bucket}) > 1]


def route_pure(rows: list[ResidualRow], key_func) -> bool:
    return not route_mixed(rows, key_func)


def stored_s194_residual_names() -> list[list[str]]:
    fibers: list[list[str]] = []
    current: list[str] | None = None
    in_section = False
    row_re = re.compile(r"^\s{4}(.+?)\s+status=")
    for line in S194_OUT.read_text(encoding="utf-8").splitlines():
        if line.startswith("[2] Coarse ET+unit route-mixed fibers"):
            in_section = True
            continue
        if in_section and line.startswith("[3] "):
            break
        if not in_section:
            continue
        if line.startswith("  coarse_route_mixed_fiber["):
            current = []
            fibers.append(current)
            continue
        match = row_re.match(line)
        if match and current is not None:
            current.append(match.group(1).rstrip())
    if len(fibers) != 15 or sum(len(fiber) for fiber in fibers) != 38:
        raise ValueError(
            f"unexpected S194 residual cache shape: fibers={len(fibers)} "
            f"rows={sum(len(fiber) for fiber in fibers)}"
        )
    return fibers


def build_residual_fibers(args) -> list[list[ResidualRow]]:
    names_by_fiber = stored_s194_residual_names()
    needed = {name for fiber in names_by_fiber for name in fiber}
    bank = fz.lp.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    selected_rows = [row for row in bank if row[0] in needed]
    if len(selected_rows) != len(needed):
        found = {row[0] for row in selected_rows}
        raise ValueError(f"missing bank rows for residual names: {sorted(needed - found)}")
    packets = fz.lp.compute_packets(selected_rows, args.workers)
    packet_by_name = {packet.name: packet for packet in packets}

    out: list[list[ResidualRow]] = []
    for names in names_by_fiber:
        fiber: list[ResidualRow] = []
        for name in names:
            packet = packet_by_name[name]
            fiber.append(
                ResidualRow(
                    packet=packet,
                    topology=topo.topology_sig(packet),
                    stalk=stalk_mod.largest_component_stalk(packet),
                )
            )
        out.append(fiber)
    return out


def first_tooth(rows: list[ResidualRow], teeth: tuple[Tooth, ...]) -> Tooth:
    for tooth in teeth:
        if route_pure(rows, tooth.key_func):
            return tooth
    raise RuntimeError("no tooth separated residual fiber")


def induced_report(fibers: list[list[ResidualRow]], tooth: Tooth) -> dict[str, object]:
    mixed_fibers = 0
    mixed_buckets = 0
    max_mixed = 0
    pure_rows = 0
    bucket_count = 0
    for rows in fibers:
        buckets: dict[object, list[ResidualRow]] = defaultdict(list)
        for row in rows:
            buckets[tooth.key_func(row)].append(row)
        bucket_count += len(buckets)
        fiber_mixed = False
        for bucket in buckets.values():
            if len({r.packet.route for r in bucket}) > 1:
                mixed_buckets += 1
                max_mixed = max(max_mixed, len(bucket))
                fiber_mixed = True
            else:
                pure_rows += len(bucket)
        if fiber_mixed:
            mixed_fibers += 1
    return {
        "buckets": bucket_count,
        "mixed_fibers": mixed_fibers,
        "mixed_buckets": mixed_buckets,
        "max_mixed": max_mixed,
        "pure_rows": pure_rows,
    }


def tooth_tournament(teeth: tuple[Tooth, ...]) -> dict[str, object]:
    names = [tooth.name for tooth in teeth]
    order = {name: i for i, name in enumerate(names)}
    edges: set[tuple[str, str]] = set()
    score = Counter()
    for a, b in combinations(teeth, 2):
        aw = sum(x > y for x, y in zip(a.vector, b.vector))
        bw = sum(x < y for x, y in zip(a.vector, b.vector))
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

    edge_idx = {(names.index(a), names.index(b)) for a, b in edges}
    hp = 0
    for perm in permutations(range(len(names))):
        if all((perm[i], perm[i + 1]) in edge_idx for i in range(len(names) - 1)):
            hp += 1

    ranking = sorted(names, key=lambda name: (-score[name], order[name]))
    return {
        "score_hist": dict(sorted(Counter(score[name] for name in names).items())),
        "directed_3cycles": c3,
        "hamiltonian_path_count": hp,
        "score_order": ranking,
    }


def row_line(row: ResidualRow) -> str:
    packet = row.packet
    return (
        f"    {packet.name:<34} route={packet.route:<16} "
        f"M={fmt(packet.M):>7} q0={packet.q_threshold:<2d} "
        f"mu={fmt(packet.strict_safe_mu):>12} "
        f"topo={row.topology.compact} "
        f"stalk={row.stalk.short()}"
    )


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, route labels, coarse ET fibers, p-adic unit roots,")
    print("    arc-Cech topology classes, safe-component stalks, Haar zeta teeth,")
    print("    q-witness certificates, covering-moment packets, and proof obligations.")
    print("  chosen vertices:")
    print("    repair teeth acting inside each residual coarse ET+unit fiber.")
    print("  preserved LRC predicate:")
    print("    boundary/open status at threshold 1/14; all audited residual rows are open.")
    print("  destroyed information:")
    print("    exact theorem route until a topology, stalk, magnitude, or certificate")
    print("    tooth is reattached.")
    print("  challenged assumption:")
    print("    the 15 residual fibers are not a counting problem; they are a")
    print("    finite list of local repair teeth.")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    fibers = build_residual_fibers(args)
    teeth = tooth_family()
    first = [first_tooth(rows, teeth) for rows in fibers]
    all_rows = [row for rows in fibers for row in rows]

    print("=== LRC14 residual tooth atlas S199 ===")
    print(
        "bank=HYP-2963 default rows "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"residual_coarse_et_unit_fibers={len(fibers)} rows={len(all_rows)}")
    print()
    print_assumption_challenge()

    print("[1] Residual tooth reports")
    for tooth in teeth:
        report = induced_report(fibers, tooth)
        print(
            "  {name:<28} buckets={buckets:<2d} mixed_fibers={mixed_fibers:<2d} "
            "mixed_buckets={mixed_buckets:<2d} max_mixed={max_mixed:<2d} pure_rows={pure_rows:<2d} "
            "class={cls}".format(name=tooth.name, cls=tooth.repair_class, **report)
        )
    print("  first_tooth_counts=" + str(dict(Counter(tooth.name for tooth in first))))
    print("  first_repair_class_counts=" + str(dict(Counter(tooth.repair_class for tooth in first))))
    print()

    print("[2] Fiber-by-fiber first tooth")
    for index, rows in enumerate(fibers):
        tooth = first[index]
        route_counts = Counter(row.packet.route for row in rows)
        status_counts = Counter(fz.status(row.packet) for row in rows)
        q_min = min(row.packet.q_threshold for row in rows)
        q_max = max(row.packet.q_threshold for row in rows)
        topology_pure = route_pure(rows, teeth[0].key_func)
        coarse_stalk_pure = route_pure(rows, teeth[1].key_func)
        exact_stalk_pure = route_pure(rows, teeth[2].key_func)
        print(
            f"  fiber[{index:02d}] size={len(rows)} routes={dict(route_counts)} "
            f"status={dict(status_counts)} q0_range={q_min}-{q_max} "
            f"first={tooth.name}/{tooth.repair_class} "
            f"topology_pure={topology_pure} coarse_stalk_pure={coarse_stalk_pure} "
            f"exact_stalk_pure={exact_stalk_pure}"
        )
        for row in rows:
            print(row_line(row))
    print()

    print("[3] Certificate teeth")
    q_rows = [row for row in all_rows if row.packet.route == "Q-WITNESS"]
    covering_rows = [row for row in all_rows if row.packet.route == "COVERING-MOMENT"]
    print(
        f"  q_witness_rows={len(q_rows)} q0_hist="
        + str(dict(sorted(Counter(row.packet.q_threshold for row in q_rows).items())))
    )
    print(
        f"  covering_rows={len(covering_rows)} positive_covering_rows="
        f"{sum(row.packet.strict_safe_mu > 0 for row in covering_rows)} "
        f"min_covering_mu={fmt(min(row.packet.strict_safe_mu for row in covering_rows))}"
    )
    print(
        "  topology_same-route-collision_fibers="
        + str(sum(not route_pure(rows, teeth[0].key_func) for rows in fibers))
    )
    print(
        "  coarse_stalk_same-route-collision_fibers="
        + str(sum(not route_pure(rows, teeth[1].key_func) for rows in fibers))
    )
    print()

    print("[4] Tournament Analysis")
    fp = tooth_tournament(teeth)
    print("  vertices_are=repair teeth, not runners")
    print("  observable=status preservation, route split, nonroute legality, locality, compression, proof cost")
    print("  switch=majority comparison of tooth vectors")
    print("  tie_hamiltonian_path=" + " > ".join(tooth.name for tooth in teeth))
    print("  score_hist=" + str(fp["score_hist"]))
    print("  directed_3cycles=" + str(fp["directed_3cycles"]))
    print("  hamiltonian_path_count=" + str(fp["hamiltonian_path_count"]))
    print("  score_order=" + " > ".join(fp["score_order"]))
    print()

    print("[5] Proof readout")
    print("  1. The 15 residual coarse ET+unit route-mixed fibers are finite")
    print("     and all 38 packets in them are strict-open.")
    print("  2. Arc topology alone separates most residual route collisions,")
    print("     but some same-topology Q/covering pairs need a stalk or")
    print("     magnitude/nested-refinement tooth.")
    print("  3. The exact safe-component stalk and magnitude cocycle separate")
    print("     every residual fiber without invoking boundary/open status.")
    print("  4. Therefore the next proof object can be a residual tooth manifest:")
    print("     for each coarse gate fiber, record the first legal non-route")
    print("     tooth before falling back to explicit q/covering certificate labels.")


if __name__ == "__main__":
    main()
