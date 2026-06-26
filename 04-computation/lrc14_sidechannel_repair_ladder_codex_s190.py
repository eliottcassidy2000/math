#!/usr/bin/env python3
"""LRC14 side-channel repair ladder for automatic quotient failures.

HYP-3016 and HYP-3017 show that automatic words are useful sidecars but unsafe
quotients: one word fiber can contain AP/Goddyn-Wong boundary atoms, q-witness
rows, petal rows, and covering rows.  This script asks the next question:

    which extra side-channel repairs the mixed fibers, and which still leaks?

The audit is intentionally theorem-facing.  It compares quotient joins over the
HYP-2963 labelled-packet bank and reports whether each join is pure with
respect to exact open/boundary status and proof route.  Tournament Analysis is
over quotient-repair carriers, not runners.
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
from typing import Callable, Iterable


REPO = Path(__file__).resolve().parents[1]
CLASSIFIER_PATH = (
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py"
)


def load_classifier():
    spec = spec_from_file_location("lrc14_packet_classifier_s186", CLASSIFIER_PATH)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules["lrc14_packet_classifier_s186"] = module
    spec.loader.exec_module(module)
    return module


lp = load_classifier()
Packet = lp.Packet


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def open_status(packet: Packet) -> str:
    return "boundary" if packet.strict_safe_mu == 0 else "open"


def tail_bucket(r: Fraction) -> str:
    if r <= 2:
        return "<=2"
    if r <= 3:
        return "(2,3]"
    if r <= 5:
        return "(3,5]"
    return ">5"


def mu_bucket(mu: Fraction) -> str:
    if mu == 0:
        return "0"
    if mu <= Fraction(1, 1000):
        return "(0,1e-3]"
    if mu <= Fraction(1, 100):
        return "(1e-3,1e-2]"
    return ">1e-2"


KeyFunc = Callable[[Packet], object]


@dataclass(frozen=True)
class Quotient:
    name: str
    kind: str
    key: KeyFunc
    keeps_magnitude: int
    keeps_boundary_topology: int
    keeps_packet_labels: int
    noncircular: int
    proof_cost: int
    note: str


@dataclass(frozen=True)
class QuotientStats:
    name: str
    kind: str
    fibers: int
    mixed_route_fibers: int
    mixed_status_fibers: int
    mixed_family_fibers: int
    max_fiber_size: int
    max_mixed_size: int
    max_route_count: int
    largest_mixed_key: object
    largest_mixed_routes: tuple[tuple[str, int], ...]
    largest_mixed_statuses: tuple[tuple[str, int], ...]
    sample: tuple[str, ...]


def quotient_bank() -> tuple[Quotient, ...]:
    return (
        Quotient(
            "automatic_word",
            "negative control",
            lambda p: p.automatic_word,
            0,
            0,
            0,
            6,
            1,
            "Moser/fibbinary word only",
        ),
        Quotient(
            "automatic_counts",
            "negative control",
            lambda p: p.automatic_counts,
            0,
            0,
            0,
            6,
            1,
            "forgets order and magnitude",
        ),
        Quotient(
            "word_plus_q_class",
            "scalar sector",
            lambda p: (p.automatic_word, p.q_class),
            1,
            0,
            0,
            6,
            2,
            "word with coarse q sector",
        ),
        Quotient(
            "word_plus_q_threshold",
            "scalar sector",
            lambda p: (p.automatic_word, p.q_threshold),
            2,
            0,
            0,
            6,
            3,
            "word with exact q-threshold",
        ),
        Quotient(
            "word_plus_q_factor_power",
            "arithmetic sidecar",
            lambda p: (p.automatic_word, p.q_factorization, p.power_lift_guard),
            1,
            0,
            0,
            6,
            3,
            "Fermat-Catalan / factor sidecar",
        ),
        Quotient(
            "word_plus_tail_bucket",
            "magnitude proxy",
            lambda p: (p.automatic_word, tail_bucket(p.lacunary_tail_ratio)),
            1,
            0,
            0,
            6,
            2,
            "coarse lacunary tail magnitude",
        ),
        Quotient(
            "word_plus_tail_exact",
            "magnitude proxy",
            lambda p: (p.automatic_word, p.lacunary_tail_ratio),
            2,
            0,
            0,
            6,
            4,
            "exact tail magnitude cocycle",
        ),
        Quotient(
            "word_plus_M",
            "magnitude cocycle",
            lambda p: (p.automatic_word, p.M),
            4,
            0,
            0,
            5,
            4,
            "exact LRC scale retained",
        ),
        Quotient(
            "word_plus_M_q",
            "magnitude cocycle",
            lambda p: (p.automatic_word, p.M, p.q_threshold),
            5,
            0,
            0,
            5,
            5,
            "exact scale plus q-threshold",
        ),
        Quotient(
            "word_plus_boundary_topology",
            "geometric sidecar",
            lambda p: (
                p.automatic_word,
                p.q_threshold,
                p.strict_components,
                p.boundary_count,
                mu_bucket(p.strict_safe_mu),
            ),
            2,
            5,
            0,
            4,
            5,
            "safe-component topology and coarse mass",
        ),
        Quotient(
            "word_plus_packet_labels",
            "labelled packet",
            lambda p: (
                p.automatic_word,
                p.packet_route,
                p.packet_rank,
                p.state_lift,
                p.transfer,
                tuple(p.unknown_pairs),
            ),
            2,
            2,
            6,
            5,
            5,
            "C27/K33/transfer labels reattached",
        ),
        Quotient(
            "word_plus_guarded_packet_no_route",
            "guarded nonroute join",
            lambda p: (
                p.automatic_word,
                p.M,
                p.q_threshold,
                p.strict_components,
                p.boundary_count,
                p.packet_route,
                p.packet_rank,
                p.state_lift,
                p.transfer,
                p.power_lift_guard,
                p.lacunary_tail_ratio,
            ),
            6,
            6,
            6,
            5,
            7,
            "maximal non-route packet sidecar",
        ),
        Quotient(
            "word_plus_filter_exit",
            "circular label check",
            lambda p: (p.automatic_word, p.automatic_filter_exit),
            1,
            1,
            3,
            1,
            2,
            "uses route-derived exit; diagnostic only",
        ),
        Quotient(
            "word_plus_route",
            "circular label check",
            lambda p: (p.automatic_word, p.route),
            1,
            1,
            6,
            0,
            1,
            "route is the target label; sanity check only",
        ),
    )


def groups_for(packets: Iterable[Packet], key: KeyFunc) -> dict[object, list[Packet]]:
    groups: dict[object, list[Packet]] = defaultdict(list)
    for packet in packets:
        groups[key(packet)].append(packet)
    return groups


def stats_for(packets: list[Packet], quotient: Quotient) -> QuotientStats:
    groups = groups_for(packets, quotient.key)
    mixed_route: list[tuple[object, list[Packet]]] = []
    mixed_status: list[tuple[object, list[Packet]]] = []
    mixed_family: list[tuple[object, list[Packet]]] = []
    max_route_count = 0
    for key, rows in groups.items():
        route_count = len({p.route for p in rows})
        max_route_count = max(max_route_count, route_count)
        if route_count > 1:
            mixed_route.append((key, rows))
        if len({open_status(p) for p in rows}) > 1:
            mixed_status.append((key, rows))
        if len({p.family for p in rows}) > 1:
            mixed_family.append((key, rows))

    if mixed_route:
        largest_key, largest_rows = max(
            mixed_route,
            key=lambda item: (len(item[1]), len({p.route for p in item[1]}), str(item[0])),
        )
    else:
        largest_key, largest_rows = None, []

    sample_rows = sorted(
        largest_rows,
        key=lambda p: (lp.ROUTE_ORDER.index(p.route), p.M, p.name),
    )[:8]
    return QuotientStats(
        name=quotient.name,
        kind=quotient.kind,
        fibers=len(groups),
        mixed_route_fibers=len(mixed_route),
        mixed_status_fibers=len(mixed_status),
        mixed_family_fibers=len(mixed_family),
        max_fiber_size=max((len(rows) for rows in groups.values()), default=0),
        max_mixed_size=len(largest_rows),
        max_route_count=max_route_count,
        largest_mixed_key=largest_key,
        largest_mixed_routes=tuple(sorted(Counter(p.route for p in largest_rows).items())),
        largest_mixed_statuses=tuple(sorted(Counter(open_status(p) for p in largest_rows).items())),
        sample=tuple(p.name for p in sample_rows),
    )


def key_preview(key: object, limit: int = 100) -> str:
    text = repr(key)
    if len(text) > limit:
        return text[: limit - 3] + "..."
    return text


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, fixed circle sections, section boundaries, wall crossings,")
    print("    residues, cover arcs, Fourier modes, endpoint owners, automatic states,")
    print("    safe-component bars, packet labels, and proof obligations.")
    print("  chosen vertices:")
    print("    quotient-repair carriers joining automatic words to one side-channel.")
    print("  preserved predicate tested:")
    print("    open/boundary status and theorem route are checked for fiber constancy.")
    print("  challenged assumption:")
    print("    a finite automaton quotient is not proof data until every mixed fiber has")
    print("    a declared magnitude, topology, endpoint, packet-label, or residual exit.")
    print()


def print_repair_ladder(stats: list[QuotientStats], quotients: tuple[Quotient, ...]) -> None:
    q_by_name = {q.name: q for q in quotients}
    print("[1] Quotient repair ladder")
    print(
        "  name                              fibers mixed_route mixed_status "
        "mixed_family max_fiber max_mixed routes note"
    )
    for st in stats:
        q = q_by_name[st.name]
        print(
            f"  {st.name:32s} {st.fibers:6d} {st.mixed_route_fibers:11d} "
            f"{st.mixed_status_fibers:12d} {st.mixed_family_fibers:12d} "
            f"{st.max_fiber_size:9d} {st.max_mixed_size:9d} "
            f"{st.max_route_count:6d} {q.note}"
        )
    print()


def print_largest_mixed(stats: list[QuotientStats]) -> None:
    print("[2] Largest mixed-route fibers by quotient")
    for st in stats:
        if st.mixed_route_fibers == 0:
            print(f"  {st.name}: route-pure")
            continue
        print(
            f"  {st.name}: key={key_preview(st.largest_mixed_key)} "
            f"rows={st.max_mixed_size} routes={dict(st.largest_mixed_routes)} "
            f"status={dict(st.largest_mixed_statuses)}"
        )
        print(f"    sample={list(st.sample)}")
    print()


def print_ap_gw_fibers(packets: list[Packet], quotients: tuple[Quotient, ...]) -> None:
    ap = next(p for p in packets if p.name == "AP")
    gw = next(p for p in packets if p.name == "GW 12->24")
    word = ap.automatic_word
    print("[3] AP/Goddyn-Wong word fiber repair")
    print(f"  boundary word={word}")
    raw_rows = [p for p in packets if p.automatic_word == word]
    print(f"  raw word rows={len(raw_rows)} routes={dict(Counter(p.route for p in raw_rows))}")
    for quotient in quotients:
        if quotient.name not in {
            "automatic_word",
            "word_plus_q_threshold",
            "word_plus_q_factor_power",
            "word_plus_tail_bucket",
            "word_plus_tail_exact",
            "word_plus_M",
            "word_plus_M_q",
            "word_plus_boundary_topology",
            "word_plus_packet_labels",
            "word_plus_guarded_packet_no_route",
        }:
            continue
        groups = groups_for(packets, quotient.key)
        ap_rows = groups[quotient.key(ap)]
        gw_rows = groups[quotient.key(gw)]
        ap_routes = dict(Counter(p.route for p in ap_rows))
        gw_routes = dict(Counter(p.route for p in gw_rows))
        same = quotient.key(ap) == quotient.key(gw)
        print(
            f"  {quotient.name:32s} same_AP_GW={same!s:5s} "
            f"AP_rows={len(ap_rows):4d} AP_routes={ap_routes} "
            f"GW_rows={len(gw_rows):4d} GW_routes={gw_routes}"
        )
    print()


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: tuple[int, ...]


def carrier_vector(stat: QuotientStats, quotient: Quotient) -> tuple[int, ...]:
    route_purity = 10000 - stat.mixed_route_fibers
    status_purity = 10000 - stat.mixed_status_fibers
    family_purity = 10000 - stat.mixed_family_fibers
    # Smaller fiber count is more discriminating, but after purity; cap so it
    # does not swamp the declared theorem-facing side-channel scores.
    discriminates = min(9999, stat.fibers)
    cost = 10 - quotient.proof_cost
    return (
        route_purity,
        status_purity,
        family_purity,
        quotient.keeps_magnitude,
        quotient.keeps_boundary_topology,
        quotient.keeps_packet_labels,
        quotient.noncircular,
        cost,
        discriminates,
    )


def tournament_mask(carriers: list[Carrier]) -> int:
    mask = 0
    bit = 0
    for i, j in combinations(range(len(carriers)), 2):
        if carriers[i].vector >= carriers[j].vector:
            mask |= 1 << bit
        bit += 1
    return mask


def print_tournament(stats: list[QuotientStats], quotients: tuple[Quotient, ...]) -> None:
    q_by_name = {q.name: q for q in quotients}
    carriers = [
        Carrier(st.name, carrier_vector(st, q_by_name[st.name]))
        for st in stats
    ]
    mask = tournament_mask(carriers)
    fp = lp.s138.tournament_fingerprint(mask, len(carriers))
    ordered = sorted(carriers, key=lambda c: c.vector, reverse=True)
    print("[4] Tournament Analysis")
    print("  vertices are quotient-repair carriers, not runners.")
    print("  pairwise observable:")
    print("    route purity, open/boundary purity, family purity, retained magnitude,")
    print("    retained topology, retained packet labels, noncircularity, proof cost,")
    print("    and discriminating fiber count.")
    print("  switch/gauge:")
    print("    A -> B iff A's observable vector is lexicographically >= B's.")
    print("  tie Hamiltonian path:")
    print("    " + " -> ".join(q.name for q in quotients))
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print("  high-retention path:")
    print("    " + " > ".join(c.name for c in ordered))
    print()


def print_theorem_readout(stats: list[QuotientStats]) -> None:
    by_name = {st.name: st for st in stats}
    print("[5] Proof-angle readout")
    print("  automatic words alone are far from route-pure.")
    print(
        "  q-threshold and q-factor/power fields reduce but do not eliminate mixed "
        "fibers."
    )
    if by_name["word_plus_M"].mixed_status_fibers == 0:
        print("  exact M attached to the automatic word already separates open/boundary status.")
    if by_name["word_plus_M"].mixed_route_fibers:
        print(
            "  exact M is not enough for theorem-route purity; it still needs packet "
            "labels or topology."
        )
    if by_name["word_plus_guarded_packet_no_route"].mixed_route_fibers == 0:
        print(
            "  the maximal non-route guarded packet signature is route-pure in this bank."
        )
    else:
        print(
            "  even the maximal non-route guarded packet signature still has mixed route "
            "fibers; inspect [2]."
        )
    print("  Candidate zipper theorem:")
    print(
        "    In any fixed automatic word fiber, the first nonzero repair cochain among "
        "exact M/q, boundary topology, C27/K33 transfer labels, tail magnitude, "
        "convergence/arc-Cech sidecars, and harmonic dual data either opens a strict "
        "component, descends to a known family, is killed by a dual certificate, or "
        "emits F7/THM-572 residual debt."
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    args = parser.parse_args()

    rows = lp.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    packets = lp.compute_packets(rows, workers=args.workers)
    quotients = quotient_bank()
    stats = [stats_for(packets, q) for q in quotients]

    print("LRC14 side-channel repair ladder")
    print("=" * 78)
    print(
        f"parameters: single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max} workers={args.workers}"
    )
    print(f"packets={len(packets)} routes={dict(Counter(p.route for p in packets))}")
    print()
    print_assumption_challenge()
    print_repair_ladder(stats, quotients)
    print_largest_mixed(stats)
    print_ap_gw_fibers(packets, quotients)
    print_tournament(stats, quotients)
    print_theorem_readout(stats)


if __name__ == "__main__":
    main()
