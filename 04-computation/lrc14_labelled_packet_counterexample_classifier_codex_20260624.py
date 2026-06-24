#!/usr/bin/env python3
"""Labelled-packet counterexample classifier for LRC14.

This script turns the S148 adversarial gauntlet into a theorem-facing packet
ledger.  The target is not a new scalar test.  The target is a finite list of
families and sporadics that any primitive LRC14 counterexample must enter after
the known q-threshold and AP-neighborhood reductions.

The packet design follows the boundary-moment bridge:

    scalar count sector        exact M, q-threshold, Farey excess
    boundary sector            strict Haar/Baire mass and boundary debt
    labelled packet sector     C27 transfer, S145 component route, K33 lift
    covering sector            q-threshold > 14 / multiple-of-14 pressure

Tournament Analysis is included over proof routes rather than runners.  The
pairwise observable is retention of exact scale, boundary status, packet owner,
state-lift address, covering information, and anti-scalarization guard.  Ties
are oriented by the declared Hamiltonian path in ROUTE_ORDER.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import os
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
THRESHOLD = Fraction(1, 14)
LOW_227 = Fraction(2, 27)
LOW_341 = Fraction(3, 41)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s138 = load_module(
    "lp_s138_frontier",
    REPO / "04-computation" / "lrc14_c27_two_swap_frontier_codex_s138.py",
)
s145 = load_module(
    "lp_s145_packet",
    REPO / "04-computation" / "lrc14_measurable_rank_recombination_codex_s145.py",
)
s147 = load_module(
    "lp_s147_haar",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)
s148 = load_module(
    "lp_s148_gauntlet",
    REPO / "04-computation" / "lrc14_adversarial_counterexample_gauntlet_codex_s148.py",
)


@dataclass(frozen=True)
class Packet:
    name: str
    source_family: str
    speeds: tuple[int, ...]
    M: Fraction
    denoms: tuple[int, ...]
    q_threshold: int
    branch: str
    farey_excess: int
    strict_safe_mu: Fraction
    strict_components: int
    boundary_count: int
    transfer: str
    packet_route: str
    packet_rank: int
    state_lift: bool
    unknown_pairs: tuple[tuple[int, int], ...]
    q_class: str
    route: str
    family: str
    theorem_role: str


ROUTE_ORDER = (
    "COUNTEREXAMPLE",
    "Q-WITNESS",
    "BOUNDARY-AP-GW",
    "BOUNDARY-PETAL-SPORADIC",
    "K33-STATE-LIFT",
    "COVERING-MOMENT",
    "SHELL-ALIAS-LOOSE",
    "MAGNITUDE-LIAR-LOOSE",
    "SOURCE-SPECTRUM-UNKNOWN",
)


ROUTE_VECTORS: dict[str, tuple[int, ...]] = {
    "COUNTEREXAMPLE": (6, 6, 6, 6, 6, 6),
    "Q-WITNESS": (5, 6, 2, 2, 2, 4),
    "BOUNDARY-AP-GW": (6, 6, 6, 5, 2, 6),
    "BOUNDARY-PETAL-SPORADIC": (5, 5, 6, 4, 2, 5),
    "K33-STATE-LIFT": (5, 5, 6, 6, 4, 6),
    "COVERING-MOMENT": (5, 4, 5, 5, 6, 6),
    "SHELL-ALIAS-LOOSE": (4, 4, 4, 3, 3, 6),
    "MAGNITUDE-LIAR-LOOSE": (4, 5, 3, 3, 4, 6),
    "SOURCE-SPECTRUM-UNKNOWN": (3, 3, 3, 3, 3, 3),
}


def fmt(frac: Fraction) -> str:
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def unique_rows(rows: list[tuple[str, str, tuple[int, ...]]]) -> list[tuple[str, str, tuple[int, ...]]]:
    seen: set[tuple[int, ...]] = set()
    out = []
    for name, family, speeds in rows:
        speeds = tuple(sorted(set(speeds)))
        if len(speeds) == 13 and primitive(speeds) and speeds not in seen:
            seen.add(speeds)
            out.append((name, family, speeds))
    return out


def handbuilt_rows() -> list[tuple[str, str, tuple[int, ...]]]:
    return [
        ("residue liar 12->26", "magnitude-liar", tuple(list(range(1, 12)) + [13, 26])),
        ("magnitude liar 12->96", "magnitude-liar", tuple(list(range(1, 12)) + [13, 96])),
        ("covering comb 12->84", "covering", tuple(list(range(1, 12)) + [13, 84])),
        ("covering repair drop13 add182", "covering", tuple(list(range(1, 13)) + [182])),
        ("scale-separated AP tail 200", "covering", tuple(list(range(1, 13)) + [200])),
        ("drop6 fattening core add180", "drop-core", tuple([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 180])),
        ("drop6 fattening core add210", "drop-core", tuple([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 210])),
    ]


def build_bank(single_limit: int, two_swap_limit: int, alias_depth: int, lcm_tail_max: int) -> list[tuple[str, str, tuple[int, ...]]]:
    rows: list[tuple[str, str, tuple[int, ...]]] = []
    rows.extend(s148.named_adversaries(alias_depth, lcm_tail_max))
    rows.extend(handbuilt_rows())
    rows.extend((f"single {s138.row_label(row)}", "single-swap", row) for row in s148.ap_single_swap_rows(single_limit))
    rows.extend((f"two {s138.row_label(row)}", "two-swap", row) for row in s138.s124.bank_two_swaps(two_swap_limit))
    return unique_rows(rows)


def q_class(q_threshold: int) -> str:
    if q_threshold <= 13:
        return "direct q<=13"
    if q_threshold == 14:
        return "apex q=14"
    return "covering q>14"


def audit_row(name: str, source_family: str, speeds: tuple[int, ...], M: Fraction, denoms: tuple[int, ...], ap_mask: int) -> Packet:
    row = s138.audit_row(speeds, M, denoms, ap_mask)
    packet = s145.classify_packet(row)
    safe = s147.exact_row_measure(speeds)
    strict_safe_mu = safe["safe_measure"]
    strict_components = len(safe["safe_components"])
    boundary_count = safe["boundary_count"]
    route, family, role = classify_route(
        source_family,
        row.label,
        M,
        row.q_threshold,
        strict_safe_mu,
        packet.route,
        packet.ph_rank,
        packet.state_lift,
        packet.unknown_pairs,
    )
    return Packet(
        name=name,
        source_family=source_family,
        speeds=speeds,
        M=M,
        denoms=denoms,
        q_threshold=row.q_threshold,
        branch=row.branch,
        farey_excess=row.farey_excess,
        strict_safe_mu=strict_safe_mu,
        strict_components=strict_components,
        boundary_count=boundary_count,
        transfer=row.transfer,
        packet_route=packet.route,
        packet_rank=packet.ph_rank,
        state_lift=packet.state_lift,
        unknown_pairs=packet.unknown_pairs,
        q_class=q_class(row.q_threshold),
        route=route,
        family=family,
        theorem_role=role,
    )


def classify_route(
    source_family: str,
    row_label: str,
    M: Fraction,
    q_threshold: int,
    strict_safe_mu: Fraction,
    packet_route: str,
    packet_rank: int,
    state_lift: bool,
    unknown_pairs: tuple[tuple[int, int], ...],
) -> tuple[str, str, str]:
    if M < THRESHOLD:
        return (
            "COUNTEREXAMPLE",
            "sporadic counterexample",
            "actual disproof row; none appeared in the audited bank",
        )
    if q_threshold <= 13:
        return (
            "Q-WITNESS",
            "infinite q-witness family",
            "direct denominator witness gives M(S) >= 1/q > 1/14",
        )
    if row_label in {"AP", "GW 12->24"} and M == THRESHOLD and strict_safe_mu == 0:
        return (
            "BOUNDARY-AP-GW",
            "AP/Goddyn-Wong boundary family",
            "endpoint-only labelled boundary atom",
        )
    if state_lift:
        return (
            "K33-STATE-LIFT",
            "K33/state-lift sporadic family",
            "nonunit packet must feed HYP-2908 / THM-572 state-lift",
        )
    if "petal" in packet_route or "GW discharge" in packet_route:
        return (
            "BOUNDARY-PETAL-SPORADIC",
            "unit-petal sporadic family",
            "unit C27 p=2 packet has positive open mass or rank-0 discharge",
        )
    if q_threshold > 14 or source_family in {"covering", "drop-core"}:
        return (
            "COVERING-MOMENT",
            "covering boundary-moment family",
            "covering branch; route to adaptive exact-period boundary-moment map",
        )
    if strict_safe_mu > 0:
        return (
            "COVERING-MOMENT",
            "open-Haar loose family",
            "positive regular-open witness interval; not a counterexample",
        )
    if source_family == "shell-alias":
        return (
            "SHELL-ALIAS-LOOSE",
            "shell-alias loose family",
            "same coarse C27 shell can be loose; exact labels are mandatory",
        )
    if source_family in {"magnitude-liar", "impostor", "scale-separated", "lcm-tail"}:
        return (
            "MAGNITUDE-LIAR-LOOSE",
            "magnitude-liar loose family",
            "raw residues/tournament shadows are magnitude-blind",
        )
    if unknown_pairs or packet_rank >= 99:
        return (
            "SOURCE-SPECTRUM-UNKNOWN",
            "unclassified sporadic",
            "zero-open, unlabeled packet; theorem must label or eliminate it",
        )
    return (
        "SOURCE-SPECTRUM-UNKNOWN",
        "unclassified sporadic",
        "zero-open, q>=14 packet not routed by current labels",
    )


def compute_packets(rows: list[tuple[str, str, tuple[int, ...]]], workers: int) -> list[Packet]:
    ap_shells, _ = s138.shell_profile(AP)
    ap_mask = s138.tournament_mask([s138.shell_priority(shell) for shell in ap_shells])
    exact = s138.compute_exact([speeds for _, _, speeds in rows], workers=workers, progress_every=0)
    packets = [
        audit_row(name, source_family, speeds, *exact[speeds], ap_mask)
        for name, source_family, speeds in rows
    ]
    return sorted(packets, key=lambda p: (ROUTE_ORDER.index(p.route), p.M, p.name))


def tournament_mask(names: list[str]) -> int:
    vectors = [ROUTE_VECTORS[name] for name in names]
    mask = 0
    bit = 0
    for i, j in combinations(range(len(names)), 2):
        vi = vectors[i]
        vj = vectors[j]
        if vi >= vj:
            mask |= 1 << bit
        bit += 1
    return mask


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, safe components, denominator walls, C27 shells,")
    print("    source-spectrum packets, K33 incidence, covering combs, Fourier/moment")
    print("    modes, and proof obligations.")
    print("  chosen vertices:")
    print("    labelled packet routes for possible counterexamples.")
    print("  preserved LRC predicate:")
    print("    exact M(S) >= 1/14, q-threshold, strict/open Haar mass, boundary debt,")
    print("    C27 owner transfer, and K33/state-lift labels.")
    print("  challenged assumption:")
    print("    neither a raw tournament isomorphism class nor a scalar M bucket is allowed")
    print("    to forget magnitude, boundary owners, or packet labels before discharge.")
    print()


def print_packet_summary(packets: list[Packet]) -> None:
    by_route = Counter(p.route for p in packets)
    by_family = Counter(p.family for p in packets)
    below = [p for p in packets if p.M < THRESHOLD]
    tight = [p for p in packets if p.M == THRESHOLD]
    low = [p for p in packets if p.M <= LOW_227]
    unknown = [p for p in packets if p.route == "SOURCE-SPECTRUM-UNKNOWN"]
    print("[1] Labelled packet census")
    print(f"  audited rows={len(packets)}")
    print(f"  below threshold={len(below)}")
    print(f"  tight rows={len(tight)}: {[p.name for p in tight[:10]]}")
    print(f"  M<=2/27 low packets={len(low)}")
    print(f"  unknown packets={len(unknown)}")
    print("  route counts:")
    for route in ROUTE_ORDER:
        if by_route[route]:
            print(f"    {route:24s} {by_route[route]}")
    print("  family counts:")
    for family, count in sorted(by_family.items()):
        print(f"    {family:34s} {count}")
    print()


def print_low_frontier(packets: list[Packet]) -> None:
    print("[2] Low-frontier families and sporadics")
    rows = [p for p in packets if p.M <= LOW_227 or p.strict_safe_mu == 0 or p.state_lift]
    rows = sorted(rows, key=lambda p: (p.M, p.route, p.name))
    print(f"  selected rows={len(rows)}")
    print(
        f"  {'name':36s} {'M':>8s} {'q0':>3s} {'mu':>9s} "
        f"{'route':24s} {'family'}"
    )
    for p in rows[:36]:
        print(
            f"  {p.name[:36]:36s} {fmt(p.M):>8s} {p.q_threshold:3d} "
            f"{fmt(p.strict_safe_mu):>9s} {p.route:24s} {p.family}"
        )
        print(f"      branch={p.branch}; transfer={p.transfer}; packet={p.packet_route}")
    if len(rows) > 36:
        print(f"  ... {len(rows) - 36} additional selected rows omitted")
    print()


def print_theorem_shape(packets: list[Packet]) -> None:
    unknown = [p for p in packets if p.route == "SOURCE-SPECTRUM-UNKNOWN"]
    covering = [p for p in packets if p.route == "COVERING-MOMENT"]
    state_lift = [p for p in packets if p.route == "K33-STATE-LIFT"]
    petals = [p for p in packets if p.route == "BOUNDARY-PETAL-SPORADIC"]
    print("[3] Candidate labelled packet theorem")
    print("  For every primitive 13-speed residual S, a source-spectrum packet P(S)")
    print("  should land in exactly one theorem bucket:")
    print("    (A) q-witness family: q_threshold(S)<=13.")
    print("    (B) AP/GW boundary atoms: q=14, strict Haar interior 0, labelled packet AP/GW.")
    print("    (C) unit-petal sporadics: C27 p=2 packet, positive open mass or rank-0 strip.")
    print("    (D) K33 sporadics: nonunit packet producing TournamentStateLift.")
    print("    (E) covering-moment family: q>14 or multiple-of-14 pressure; exact-period")
    print("        boundary maps to positive gK8/L_y moment slack.")
    print("    (F) new sporadic: any SOURCE-SPECTRUM-UNKNOWN packet.  In this audited bank,")
    print(f"        its count is {len(unknown)}.")
    print()
    print("  Computed pressure in the audited bank:")
    print(f"    covering-moment rows={len(covering)}")
    print(f"    K33/state-lift rows={len(state_lift)}")
    print(f"    unit-petal rows={len(petals)}")
    print("  This is a classification scaffold, not a proof of the global reduction.")
    print()


def print_swap_chain_analogy() -> None:
    print("[4] Binary fixed-margin swap-chain analogy")
    print("  Fu-Qin-Wang arXiv:2606.22636 proves an inverse-polynomial spectral")
    print("  gap for binary fixed-margin swap chains by comparing to two-row")
    print("  heat-bath moves, reducing to three rows, then splitting a scalar count")
    print("  sector from Johnson-harmonic non-scalar sectors.")
    print("  LRC14 import:")
    print("    scalar count sector       -> q_threshold, exact M, Farey excess")
    print("    three-row reduction       -> AP/GW/petal/K33 local packet atlas")
    print("    Johnson harmonic sectors  -> owner-labelled C27/K33/boundary modes")
    print("    heat-bath comparison      -> packet-preserving local replacement moves")
    print("  The proof lesson is conservative: prove comparisons only after conditioning")
    print("  on the labels that make scalar count sectors honest.")
    print()


def print_tournament_analysis() -> None:
    print("[5] Tournament Analysis")
    names = list(ROUTE_ORDER)
    mask = tournament_mask(names)
    fp = s138.tournament_fingerprint(mask, len(names))
    print("  vertices are counterexample-classification routes, not runners.")
    print("  pair observable:")
    print("    exact scale retention, boundary retention, packet-owner retention,")
    print("    state-lift visibility, covering visibility, anti-scalar guard.")
    print("  switch:")
    print("    route A -> route B iff A's retention vector is lexicographically >= B's.")
    print("  tie Hamiltonian path:")
    print("    " + " -> ".join(names))
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print()


def print_machine_readable(packets: list[Packet], emit_all: bool) -> None:
    print("[6] Machine-readable packet rows")
    selected = packets
    if not emit_all:
        chosen: list[Packet] = []
        for route in ROUTE_ORDER:
            route_rows = [p for p in packets if p.route == route]
            chosen.extend(route_rows[:8])
        for p in packets:
            if p.M <= LOW_227 and p not in chosen:
                chosen.append(p)
        selected = chosen
        print("  representative rows only; pass --emit-all-rows for the full packet ledger")
    print("  name|source|M|q|mu|route|family|role")
    for p in selected:
        print(
            f"  {p.name}|{p.source_family}|{fmt(p.M)}|{p.q_threshold}|"
            f"{fmt(p.strict_safe_mu)}|{p.route}|{p.family}|{p.theorem_role}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    parser.add_argument("--emit-all-rows", action="store_true")
    args = parser.parse_args()

    print("LRC14 labelled-packet counterexample classifier")
    print("=" * 78)
    print(
        f"single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}, workers={args.workers}"
    )
    print_assumption_challenge()

    rows = build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    packets = compute_packets(rows, args.workers)
    print_packet_summary(packets)
    print_low_frontier(packets)
    print_theorem_shape(packets)
    print_swap_chain_analogy()
    print_tournament_analysis()
    print_machine_readable(packets, args.emit_all_rows)


if __name__ == "__main__":
    main()
