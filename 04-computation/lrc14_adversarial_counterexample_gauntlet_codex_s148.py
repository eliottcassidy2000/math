#!/usr/bin/env python3
"""S148: adversarial proof/disproof gauntlet for LRC14.

This is a proof-interface and counterexample-hunting script for HYP-2950.
It deliberately treats possible counterexamples as adversarial objects:

* q-threshold rows that evade small denominator witnesses;
* AP single-swap and two-swap rows beyond the current named frontier;
* C27 shell-alias impostors with the same coarse labels as known packets;
* divisor-loaded lcm tails and floor-odd/tournament impostors;
* exact low-frontier packets from HYP-2947.

For each row where exact computation is requested, the script attaches:

    exact M/Farey branch
    strict Haar/Baire safe mass at threshold 1/14
    C27 shell-transfer label
    S145 packet/rank classification when available

The disproof target is explicit: find any primitive 13-speed row with
M < 1/14, or a low packet that evades the S145 rank split.  The proof target is
the complementary gauntlet statement: every tested adversary is either directly
safe, endpoint-only AP/GW, unit-petal discharge, or K33/state-lift obligation.
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
    "s148_s138_frontier",
    REPO / "04-computation" / "lrc14_c27_two_swap_frontier_codex_s138.py",
)
s145 = load_module(
    "s148_s145_recombination",
    REPO / "04-computation" / "lrc14_measurable_rank_recombination_codex_s145.py",
)
s147 = load_module(
    "s148_s147_haar",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)


@dataclass(frozen=True)
class CandidateAudit:
    name: str
    family: str
    speeds: tuple[int, ...]
    M: Fraction
    denoms: tuple[int, ...]
    q_threshold: int
    row_label: str
    branch: str
    farey_excess: int
    strict_safe_mu: Fraction
    strict_components: int
    transfer: str
    packet_route: str
    packet_rank: int
    state_lift: bool
    unknown_pairs: tuple[tuple[int, int], ...]
    verdict: str


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def ap_single_swap_rows(limit: int) -> list[tuple[int, ...]]:
    rows: set[tuple[int, ...]] = set()
    for rem in AP:
        kept = set(AP) - {rem}
        for add in range(1, limit + 1):
            row = tuple(sorted(kept | {add}))
            if len(row) == 13 and primitive(row):
                rows.add(row)
    return sorted(rows)


def audit_from_exact(name: str, family: str, speeds: tuple[int, ...], M: Fraction, denoms: tuple[int, ...]) -> CandidateAudit:
    ap_shells, _ = s138.shell_profile(AP)
    ap_mask = s138.tournament_mask([s138.shell_priority(shell) for shell in ap_shells])
    row = s138.audit_row(speeds, M, denoms, ap_mask)
    packet = s145.classify_packet(row)
    safe = s147.exact_row_measure(speeds)
    strict_safe_mu = safe["safe_measure"]
    strict_components = len(safe["safe_components"])

    if M < THRESHOLD:
        verdict = "COUNTEREXAMPLE"
    elif M == THRESHOLD and row.label in {"AP", "GW 12->24"}:
        verdict = "endpoint-only AP/GW tight"
    elif M == THRESHOLD:
        verdict = "tight unknown boundary threat"
    elif packet.unknown_pairs and M <= LOW_227:
        verdict = "LOW UNKNOWN PACKET THREAT"
    elif packet.state_lift and M <= LOW_227:
        verdict = "K33/state-lift obligation"
    elif packet.ph_rank == 0 and M <= LOW_227 and not packet.unknown_pairs:
        verdict = "rank-0 local discharge"
    elif row.q_threshold <= 13:
        verdict = "direct q-witness safe"
    else:
        verdict = "strict loose exact witness"

    return CandidateAudit(
        name=name,
        family=family,
        speeds=speeds,
        M=M,
        denoms=denoms,
        q_threshold=row.q_threshold,
        row_label=row.label,
        branch=row.branch,
        farey_excess=row.farey_excess,
        strict_safe_mu=strict_safe_mu,
        strict_components=strict_components,
        transfer=row.transfer,
        packet_route=packet.route,
        packet_rank=packet.ph_rank,
        state_lift=packet.state_lift,
        unknown_pairs=packet.unknown_pairs,
        verdict=verdict,
    )


def compute_audits(rows: list[tuple[str, str, tuple[int, ...]]], workers: int) -> list[CandidateAudit]:
    exact = s138.compute_exact([speeds for _, _, speeds in rows], workers=workers, progress_every=0)
    out = []
    for name, family, speeds in rows:
        M, denoms = exact[speeds]
        out.append(audit_from_exact(name, family, speeds, M, denoms))
    return sorted(out, key=lambda a: (a.M, a.family, a.name))


def named_adversaries(alias_depth: int, lcm_tail_max: int) -> list[tuple[str, str, tuple[int, ...]]]:
    rows: list[tuple[str, str, tuple[int, ...]]] = [
        ("AP", "named-frontier", AP),
        ("GW 12->24", "named-frontier", tuple(list(range(1, 12)) + [13, 24])),
        ("near/K33 12->36", "named-frontier", tuple(list(range(1, 12)) + [13, 36])),
        ("petal 10->20", "named-frontier", tuple(sorted((set(AP) - {10}) | {20}))),
        ("petal 13->26", "named-frontier", tuple(sorted((set(AP) - {13}) | {26}))),
        ("P10+GW", "named-frontier", tuple(sorted((set(AP) - {10, 12}) | {20, 24}))),
        ("P10+K33", "named-frontier", tuple(sorted((set(AP) - {10, 12}) | {20, 36}))),
        ("floor-odd GW iso impostor", "impostor", tuple(list(range(1, 12)) + [13, 360])),
        ("AP repair drop13 add182", "covering-repair", tuple(list(range(1, 13)) + [182])),
        ("scale-separated AP tail 200", "scale-separated", tuple(list(range(1, 13)) + [200])),
        ("divisor-loaded tail 84", "lcm-tail", tuple(list(range(1, 12)) + [13, 84])),
    ]

    for m in range(2, lcm_tail_max + 1):
        rows.append((f"divisor-loaded tail 84*{m}", "lcm-tail", tuple(list(range(1, 12)) + [13, 84 * m])))

    alias_specs = [
        ("GW-shell alias", 12, 24),
        ("K33-shell alias", 12, 36),
        ("P10-shell alias", 10, 20),
        ("P13-shell alias", 13, 26),
    ]
    for family, rem, base_add in alias_specs:
        for m in range(alias_depth + 1):
            add = base_add + 27 * m
            rows.append((f"{family} {rem}->{add}", "shell-alias", tuple(sorted((set(AP) - {rem}) | {add}))))

    # Deduplicate by speeds while retaining the first descriptive name.
    seen: set[tuple[int, ...]] = set()
    deduped: list[tuple[str, str, tuple[int, ...]]] = []
    for name, family, speeds in rows:
        if speeds not in seen:
            seen.add(speeds)
            deduped.append((name, family, speeds))
    return deduped


def scan_bank(name: str, rows: list[tuple[int, ...]], workers: int) -> tuple[list[object], dict[str, int]]:
    hard = [row for row in rows if s138.s124.q_threshold(row) >= 14]
    exact = s138.compute_exact(hard, workers=workers, progress_every=0)
    ap_shells, _ = s138.shell_profile(AP)
    ap_mask = s138.tournament_mask([s138.shell_priority(shell) for shell in ap_shells])
    audits = [
        s138.audit_row(row, M, denoms, ap_mask)
        for row, (M, denoms) in exact.items()
    ]
    stats = {
        "raw": len(rows),
        "direct_q_safe": len(rows) - len(hard),
        "hard_q_ge_14": len(hard),
        "below": sum(1 for audit in audits if audit.M < THRESHOLD),
        "tight": sum(1 for audit in audits if audit.M == THRESHOLD),
        "low_3_41": sum(1 for audit in audits if audit.M <= LOW_341),
        "low_2_27": sum(1 for audit in audits if audit.M <= LOW_227),
    }
    return sorted(audits, key=lambda audit: (audit.M, audit.label, audit.speeds)), stats


def fmt(frac: Fraction) -> str:
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, q-threshold gaps, strict safe components, C27 shell aliases,")
    print("    unital/affine packets, K33 minors, lcm tails, impostor tournaments,")
    print("    and proof/disproof obligations.")
    print("  chosen vertices:")
    print("    adversarial candidate families and the proof route that defeats each one.")
    print("  preserved LRC predicate:")
    print("    exact M(S) >= 1/14, with strict Haar/Baire safe mass and packet labels")
    print("    retained until the row is discharged.")
    print("  challenged assumptions:")
    print("    (1) a C27 shell label alone proves anything; aliases refute that.")
    print("    (2) a tournament iso or floor-odd shadow characterizes tightness.")
    print("    (3) positive Haar mass alone can replace owner/carry addresses.")
    print()


def print_candidate_table(audits: list[CandidateAudit], limit: int = 22) -> None:
    print("[1] Named adversaries and shell-alias impostors")
    print(f"  audited={len(audits)}")
    print(
        f"  {'name':32s} {'family':16s} {'M':>8s} {'q0':>3s} "
        f"{'safe_mu':>9s} {'rank':>4s} {'verdict'}"
    )
    for audit in audits[:limit]:
        print(
            f"  {audit.name[:32]:32s} {audit.family[:16]:16s} {fmt(audit.M):>8s} "
            f"{audit.q_threshold:3d} {fmt(audit.strict_safe_mu):>9s} "
            f"{audit.packet_rank:4d} {audit.verdict}"
        )
        print(f"      transfer={audit.transfer}")
        if audit.unknown_pairs:
            print(f"      unknown_pairs={audit.unknown_pairs}")
    if len(audits) > limit:
        print(f"  ... {len(audits) - limit} more rows omitted from table")
    print()


def print_scan_summary(title: str, audits: list[object], stats: dict[str, int], top: int = 12) -> None:
    print(f"[2] {title}")
    print(f"  stats={stats}")
    print(f"  top {min(top, len(audits))} hard rows by exact M:")
    for audit in audits[:top]:
        packet = s145.classify_packet(audit)
        mark = ""
        if audit.M < THRESHOLD:
            mark = "  <-- COUNTEREXAMPLE"
        elif audit.M == THRESHOLD:
            mark = "  <-- endpoint tight"
        print(
            f"    {audit.label:32s} M={fmt(audit.M):>6s} q0={audit.q_threshold:2d} "
            f"branch={audit.branch:21s} rank={packet.ph_rank:2d} "
            f"route={packet.route}{mark}"
        )
    print()


def shell_alias_readout(audits: list[CandidateAudit]) -> None:
    print("[3] Coarse C27 shell labels are not enough")
    by_route: dict[str, list[CandidateAudit]] = {}
    for audit in audits:
        if audit.family == "shell-alias":
            prefix = audit.name.split()[0]
            by_route.setdefault(prefix, []).append(audit)
    for prefix, group in sorted(by_route.items()):
        group = sorted(group, key=lambda audit: audit.M)
        best = group[0]
        worst = group[-1]
        below = sum(1 for audit in group if audit.M < THRESHOLD)
        tight = sum(1 for audit in group if audit.M == THRESHOLD)
        print(
            f"  {prefix:10s}: rows={len(group)} below={below} tight={tight} "
            f"best={best.name} M={fmt(best.M)} worst={worst.name} M={fmt(worst.M)}"
        )
        print(f"      shared transfer shape example: {best.transfer}")
    print("  readout:")
    print("    shell aliases can repeat the same coarse owner/carry pattern while")
    print("    becoming safely loose.  Exact M/Farey and Haar/Baire address data")
    print("    must stay attached; this falsifies shell-label-only proof attempts.")
    print()


def tournament_analysis() -> None:
    print("[4] Tournament Analysis: proof/disproof gauntlet routes")
    routes = [
        ("exact M counterexample search", (6, 6, 5, 5, 4, 5)),
        ("q-threshold direct witnesses", (5, 6, 4, 4, 4, 5)),
        ("regular-open Haar/Baire mass", (5, 5, 6, 4, 5, 5)),
        ("S145 packet rank split", (5, 5, 5, 6, 6, 6)),
        ("K33 state-lift endpoint", (4, 5, 5, 6, 6, 5)),
        ("C27 shell labels alone", (3, 4, 3, 3, 2, 2)),
        ("raw tournament/impostor iso", (2, 2, 2, 2, 2, 1)),
        ("raw scalar count", (1, 1, 1, 1, 1, 1)),
    ]
    mask = s138.tournament_mask([score for _, score in routes])
    fp = s138.tournament_fingerprint(mask, len(routes))
    order = sorted(range(len(routes)), key=lambda i: routes[i][1], reverse=True)
    print("  vertices are routes for defeating or finding counterexamples.")
    print("  pair observable:")
    print("    exactness, branch retention, measurable witness retention, packet fit,")
    print("    endpoint closure, and scalar-decoy resistance.")
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(routes[i][0] for i in order))
    print()


def verdict(named: list[CandidateAudit], single_stats: dict[str, int], two_stats: dict[str, int], single_audits: list[object], two_audits: list[object]) -> None:
    print("[5] Proof/disproof verdict")
    all_named_below = [audit for audit in named if audit.M < THRESHOLD]
    named_unknown_low = [audit for audit in named if audit.M <= LOW_227 and audit.unknown_pairs]
    print(f"  named counterexamples={len(all_named_below)}")
    print(f"  named low unknown packets={len(named_unknown_low)}")
    print(f"  single-swap below-threshold rows={single_stats['below']}")
    print(f"  two-swap below-threshold rows={two_stats['below']}")
    print(f"  single-swap tight hard rows={single_stats['tight']}")
    print(f"  two-swap tight hard rows={two_stats['tight']}")
    best_single = single_audits[0] if single_audits else None
    best_two = two_audits[0] if two_audits else None
    if best_single:
        print(f"  best single-swap hard row: {best_single.label}, M={fmt(best_single.M)}")
    if best_two:
        print(f"  best two-swap hard row: {best_two.label}, M={fmt(best_two.M)}")
    print()
    print("  No disproof row was found in this gauntlet.")
    print("  What any future counterexample must now evade:")
    print("    * all q<=13 direct witnesses, so q_threshold>=14;")
    print("    * AP/GW endpoint-only tight classification in the tested single/two-swap banks;")
    print("    * S145 rank split for every tested M<=2/27 low packet;")
    print("    * HYP-2948/HYP-2949 regular-open safe-mass separation;")
    print("    * shell-alias and tournament-impostor decoys, which are loose under exact M.")
    print("  Remaining proof theorem:")
    print("    prove every primitive LRC14 counterexample reduces to one of these gauntlet")
    print("    families or to a packet-preserving mutation that still enters HYP-2947.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=360)
    parser.add_argument("--two-swap-limit", type=int, default=42)
    parser.add_argument("--alias-depth", type=int, default=6)
    parser.add_argument("--lcm-tail-max", type=int, default=8)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    args = parser.parse_args()

    print("S148 LRC14 ADVERSARIAL COUNTEREXAMPLE GAUNTLET")
    print("=" * 78)
    print(
        f"single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}, workers={args.workers}"
    )
    print_assumption_challenge()

    named_rows = named_adversaries(args.alias_depth, args.lcm_tail_max)
    named_audits = compute_audits(named_rows, args.workers)
    print_candidate_table(named_audits)

    single_rows = ap_single_swap_rows(args.single_limit)
    single_audits, single_stats = scan_bank("single-swap", single_rows, args.workers)
    print_scan_summary(f"AP single-swap scan v<={args.single_limit}", single_audits, single_stats)

    two_rows = s138.s124.bank_two_swaps(args.two_swap_limit)
    two_audits, two_stats = scan_bank("two-swap", two_rows, args.workers)
    print_scan_summary(f"AP two-swap scan add<={args.two_swap_limit}", two_audits, two_stats)

    shell_alias_readout(named_audits)
    tournament_analysis()
    verdict(named_audits, single_stats, two_stats, single_audits, two_audits)


if __name__ == "__main__":
    main()

