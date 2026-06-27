#!/usr/bin/env python3
"""S154: boundary-moment packet ledger for the LRC14 labelled theorem route.

This script extends the HYP-2963 labelled-packet counterexample audit with an
explicit exact-period sector ledger for the HYP-2954 boundary-moment bridge.

The ledger is a conservative proxy, not the final gK8 theorem.  For a row S and
denominator D it scans unit packets a in (Z/DZ)^* and records:

  * whether the LRC14 threshold test at t=a/D is covered, boundary, or strict;
  * which of the six 1/7-sector slots are hit by the phases s*a/D;
  * the missed-depth histogram q_0,...,q_6;
  * the moment readout L_y = 10*q_0 + q_3 + 10*q_6.

The purpose is to make the live COVERING-MOMENT bucket theorem-facing.  Any
strict counterexample would have to survive qdiv, Haar-open migration, C27/K33
labelling, and then appear here as a zero-open, all-covered, non-K33 packet
with no positive moment exit.  The audited bounded bank finds none.

Tournament Analysis declaration:
  vertices: proof packet routes / moment sectors, not runners.
  pairwise observable: which layer preserves strict-counterexample status
    before scalarization.
  switch: componentwise retention vector, ties by the declared Hamiltonian path.
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
THRESHOLD = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


lp = load_module(
    "s154_labelled_packet",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)
s138 = load_module(
    "s154_s138_frontier",
    REPO / "04-computation" / "lrc14_c27_two_swap_frontier_codex_s138.py",
)


@dataclass(frozen=True)
class MomentLedger:
    name: str
    route: str
    family: str
    q_threshold: int
    M: Fraction
    strict_safe_mu: Fraction
    denominator: int
    units: int
    covered_units: int
    boundary_units: int
    strict_units: int
    min_slack_num: int
    missed_hist: tuple[int, ...]
    ly_count: int
    ly_norm: Fraction
    extreme_norm: Fraction
    sector_support_hist: tuple[int, ...]
    kernel_flag: str


ROUTE_ORDER = (
    "Q-WITNESS",
    "BOUNDARY-AP-GW",
    "BOUNDARY-PETAL-SPORADIC",
    "K33-STATE-LIFT",
    "COVERING-MOMENT",
    "SHELL-ALIAS-LOOSE",
    "MAGNITUDE-LIAR-LOOSE",
    "SOURCE-SPECTRUM-UNKNOWN",
    "COUNTEREXAMPLE",
)


LAYER_VECTORS = {
    "qdiv_gate": (8, 5, 2, 2, 1, 5),
    "exact_M_Farey": (7, 7, 2, 2, 2, 6),
    "Haar_Baire_front": (6, 6, 6, 3, 2, 6),
    "C27_K33_labels": (4, 4, 4, 8, 4, 7),
    "fixed_margin_packet": (5, 5, 5, 6, 6, 7),
    "boundary_moment_ledger": (4, 4, 5, 6, 8, 8),
    "TournamentStateLift": (3, 3, 3, 7, 7, 8),
    "raw_scalar_or_residue": (1, 1, 1, 1, 1, 1),
}


def choose_denominator(packet: object, cap: int) -> int:
    """Choose a small exact-period chart that still sees the packet branch."""
    denoms = tuple(getattr(packet, "denoms", ())) or ()
    candidates: list[int] = []
    if getattr(packet, "q_threshold", 0) >= 2:
        candidates.append(int(packet.q_threshold))
    candidates.extend(int(d) for d in denoms if int(d) >= 2)
    candidates.append(14)
    # Add 27 for C27/K33 packets and 7*13 for covering tails when affordable.
    transfer = getattr(packet, "transfer", "")
    if "g3" in transfer or "g9" in transfer or "K33" in getattr(packet, "packet_route", ""):
        candidates.append(27)
    if getattr(packet, "q_threshold", 0) > 14:
        candidates.append(91)

    for d in candidates:
        if 2 <= d <= cap:
            return d
    return min(max(candidates[0], 2), cap)


def sector_mask_for(speeds: tuple[int, ...], a: int, D: int) -> int:
    """Mark six nonzero seventh-sectors visited by phases s*a/D."""
    mask = 0
    for s in speeds:
        r = (s * a) % D
        sector = (7 * r) // D
        if 1 <= sector <= 6:
            mask |= 1 << (sector - 1)
    return mask


def unit_moment_ledger(packet: object, D: int) -> MomentLedger:
    speeds = tuple(packet.speeds)
    covered = boundary = strict = 0
    min_slack_num: int | None = None
    missed_hist = [0] * 7
    support_hist = [0] * 7

    for a in range(1, D):
        if gcd(a, D) != 1:
            continue
        best_dist = min(min((s * a) % D, D - ((s * a) % D)) for s in speeds)
        slack_num = 14 * best_dist - D
        min_slack_num = slack_num if min_slack_num is None else min(min_slack_num, slack_num)
        if slack_num < 0:
            covered += 1
        elif slack_num == 0:
            boundary += 1
        else:
            strict += 1

        mask = sector_mask_for(speeds, a, D)
        support = mask.bit_count()
        missed = 6 - support
        support_hist[support] += 1
        missed_hist[missed] += 1

    units = covered + boundary + strict
    ly_count = 10 * missed_hist[0] + missed_hist[3] + 10 * missed_hist[6]
    ly_norm = Fraction(ly_count, units) if units else Fraction(0)
    extreme_norm = Fraction(missed_hist[0] + missed_hist[6], units) if units else Fraction(0)

    kernel_flag = classify_kernel(packet, covered, boundary, strict, tuple(missed_hist))
    return MomentLedger(
        packet.name,
        packet.route,
        packet.family,
        packet.q_threshold,
        packet.M,
        packet.strict_safe_mu,
        D,
        units,
        covered,
        boundary,
        strict,
        min_slack_num or 0,
        tuple(missed_hist),
        ly_count,
        ly_norm,
        extreme_norm,
        tuple(support_hist),
        kernel_flag,
    )


def classify_kernel(packet: object, covered: int, boundary: int, strict: int, missed_hist: tuple[int, ...]) -> str:
    if packet.M < THRESHOLD:
        return "ACTUAL-COUNTEREXAMPLE"
    if packet.strict_safe_mu == 0 and packet.route == "BOUNDARY-AP-GW":
        return "AP/GW-zero-open-equality"
    if packet.state_lift:
        return "named-K33-state-lift"
    if packet.route == "BOUNDARY-PETAL-SPORADIC":
        return "unit-petal-named"
    if packet.route == "Q-WITNESS":
        return "q-witness-discharge"
    if packet.strict_safe_mu > 0:
        return "positive-Haar-open"
    if strict == 0 and boundary == 0 and packet.route not in {"BOUNDARY-AP-GW", "K33-STATE-LIFT"}:
        if missed_hist[0] == 0 and missed_hist[6] == 0:
            return "ZERO-OPEN-ALL-COVERED-MOMENT-BLIND"
        return "ZERO-OPEN-ALL-COVERED-MOMENT-EXTREME"
    if packet.route == "SOURCE-SPECTRUM-UNKNOWN":
        return "SOURCE-SPECTRUM-UNKNOWN"
    return "labelled-nonbad"


def build_packets(args: argparse.Namespace) -> list[object]:
    if args.full_bank:
        rows = lp.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    else:
        rows = curated_bank(args.alias_depth, args.lcm_tail_max)
    packets = lp.compute_packets(rows, args.workers)
    return packets


def curated_bank(alias_depth: int, lcm_tail_max: int) -> list[tuple[str, str, tuple[int, ...]]]:
    """Small default bank covering each named theorem bucket."""
    rows: list[tuple[str, str, tuple[int, ...]]] = []
    rows.extend(lp.s148.named_adversaries(alias_depth, lcm_tail_max))
    rows.extend(lp.handbuilt_rows())
    rows.extend(
        [
            ("direct-q decoy 12->28", "fixed-margin-decoy", tuple(list(range(1, 12)) + [13, 28])),
            ("direct-q decoy 12->56", "fixed-margin-decoy", tuple(list(range(1, 12)) + [13, 56])),
            ("covering lcm family 12->168", "covering", tuple(list(range(1, 12)) + [13, 168])),
            (
                "K33 larger marker drop(12,13)->add(26,36)",
                "K33-marker",
                tuple(sorted((set(lp.AP) - {12, 13}) | {26, 36})),
            ),
            (
                "unit petal splice drop(10,13)->add(20,26)",
                "unit-petal",
                tuple(sorted((set(lp.AP) - {10, 13}) | {20, 26})),
            ),
            (
                "S138 splice drop(10,12)->add(20,24)",
                "S138-splice",
                tuple(sorted((set(lp.AP) - {10, 12}) | {20, 24})),
            ),
            (
                "S138 splice drop(10,12)->add(20,36)",
                "S138-splice",
                tuple(sorted((set(lp.AP) - {10, 12}) | {20, 36})),
            ),
        ]
    )
    return lp.unique_rows(rows)


def selected_packets(packets: list[object], sample_per_route: int, include_all_covering: bool) -> list[object]:
    chosen: list[object] = []
    seen: set[tuple[int, ...]] = set()

    def add(packet: object) -> None:
        key = tuple(packet.speeds)
        if key not in seen:
            seen.add(key)
            chosen.append(packet)

    for route in ROUTE_ORDER:
        route_rows = [p for p in packets if p.route == route]
        route_rows = sorted(route_rows, key=lambda p: (p.M, p.name))
        for p in route_rows[:sample_per_route]:
            add(p)
    for p in packets:
        if p.M <= Fraction(2, 27) or p.strict_safe_mu == 0 or p.state_lift:
            add(p)
    if include_all_covering:
        for p in packets:
            if p.route == "COVERING-MOMENT":
                add(p)
    return chosen


def tournament_fingerprint(names: list[str]) -> dict[str, object]:
    mask = 0
    bit = 0
    vectors = [LAYER_VECTORS[name] for name in names]
    for i, j in combinations(range(len(names)), 2):
        vi, vj = vectors[i], vectors[j]
        si = sum(1 for a, b in zip(vi, vj) if a >= b)
        sj = sum(1 for a, b in zip(vi, vj) if b >= a)
        if si > sj or (si == sj and i < j):
            mask |= 1 << bit
        bit += 1
    return s138.tournament_fingerprint(mask, len(names))


def fmt(frac: Fraction) -> str:
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, residues, exact-period unit packets, missed-sector states,")
    print("    qdiv gates, C27/K33 owner labels, boundary-moment vectors, and proof")
    print("    obligations.")
    print("  chosen vertices:")
    print("    labelled packet routes plus exact-period missed-depth moment ledgers.")
    print("  preserved LRC predicate:")
    print("    strict-counterexample status through qdiv, exact M, Haar-open mass,")
    print("    boundary debt, C27/K33 labels, and all-covered unit-packet tests.")
    print("  destroyed information:")
    print("    continuous torus geometry between probed denominators and the true gK8")
    print("    feasible moment region.  This is a scout ledger, not the final bridge.")
    print("  challenged assumption:")
    print("    the COVERING-MOMENT bucket can remain a word.  It must emit a moment")
    print("    signature, and any all-covered non-K33 signature is the next target.")
    print()


def print_ledger_summary(packets: list[object], ledgers: list[MomentLedger]) -> None:
    print("[1] Boundary-moment ledger census")
    print(f"  source packets audited={len(packets)}")
    print(f"  moment ledgers emitted={len(ledgers)}")
    print(f"  below-threshold packets={sum(1 for p in packets if p.M < THRESHOLD)}")
    print(f"  zero-open packets={sum(1 for p in packets if p.strict_safe_mu == 0)}")
    print("  route counts:")
    for route, count in sorted(Counter(p.route for p in packets).items()):
        print(f"    {route:24s} {count}")
    print("  kernel flags among emitted ledgers:")
    for flag, count in sorted(Counter(l.kernel_flag for l in ledgers).items()):
        print(f"    {flag:38s} {count}")
    print()


def print_route_moment_stats(ledgers: list[MomentLedger]) -> None:
    print("[2] Route-level moment statistics")
    by_route: dict[str, list[MomentLedger]] = defaultdict(list)
    for ledger in ledgers:
        by_route[ledger.route].append(ledger)
    for route in ROUTE_ORDER:
        rows = by_route.get(route, [])
        if not rows:
            continue
        all_covered = sum(1 for r in rows if r.strict_units == 0 and r.boundary_units == 0)
        min_ly = min(r.ly_norm for r in rows)
        max_ly = max(r.ly_norm for r in rows)
        min_extreme = min(r.extreme_norm for r in rows)
        max_extreme = max(r.extreme_norm for r in rows)
        print(
            f"  {route:24s} n={len(rows):4d} all-covered@D={all_covered:3d} "
            f"Ly_norm=[{fmt(min_ly)}, {fmt(max_ly)}] "
            f"extreme=[{fmt(min_extreme)}, {fmt(max_extreme)}]"
        )
    print()


def print_representatives(ledgers: list[MomentLedger], limit: int) -> None:
    print("[3] Representative ledgers")
    ordered = sorted(
        ledgers,
        key=lambda r: (
            r.kernel_flag != "ZERO-OPEN-ALL-COVERED-MOMENT-BLIND",
            r.kernel_flag != "ZERO-OPEN-ALL-COVERED-MOMENT-EXTREME",
            r.route,
            r.ly_norm,
            r.name,
        ),
    )
    print(
        f"  {'name':34s} {'route':24s} {'D':>4s} {'state(c/b/s)':>14s} "
        f"{'missed q0..q6':>23s} {'Ly':>8s} {'flag'}"
    )
    for r in ordered[:limit]:
        state = f"{r.covered_units}/{r.boundary_units}/{r.strict_units}"
        print(
            f"  {r.name[:34]:34s} {r.route:24s} {r.denominator:4d} "
            f"{state:>14s} {str(r.missed_hist):>23s} {fmt(r.ly_norm):>8s} "
            f"{r.kernel_flag}"
        )
    if len(ordered) > limit:
        print(f"  ... {len(ordered) - limit} additional emitted ledgers omitted")
    print()


def print_theorem_readout(ledgers: list[MomentLedger]) -> None:
    bad = [
        r for r in ledgers
        if r.kernel_flag.startswith("ZERO-OPEN-ALL-COVERED")
        or r.kernel_flag == "SOURCE-SPECTRUM-UNKNOWN"
        or r.kernel_flag == "ACTUAL-COUNTEREXAMPLE"
    ]
    print("[4] Labelled packet theorem readout")
    print("  Desired theorem shape:")
    print("    Every primitive LRC14 row emits a fixed-margin labelled packet.")
    print("    If it is strict-bad, qdiv>14 and its packet must sit in a covering")
    print("    boundary-moment fiber.  That fiber either has positive Haar/moment")
    print("    exit, carries a named K33 state-lift label, or is an explicit new")
    print("    Johnson-harmonic sector.")
    print(f"  audited dangerous moment-kernel rows={len(bad)}")
    if bad:
        for r in bad[:12]:
            print(f"    {r.name}: route={r.route}, D={r.denominator}, flag={r.kernel_flag}")
    else:
        print("    none in the emitted bounded ledger.")
    print("  Next rigorous target:")
    print("    replace this finite sector proxy by the true gK8/L_y feasible-region")
    print("    map B_D and prove every COVERING-MOMENT fixed-margin fiber has positive")
    print("    image unless it constructs HYP-2908/THM-572.")
    print()


def print_arxiv_bridge() -> None:
    print("[5] arXiv:2606.22636 proof-shape import")
    print("  Fu-Qin-Wang prove a spectral gap for binary fixed-margin swap chains by")
    print("  keeping fixed margins, comparing local swap dynamics to two-row heat-bath")
    print("  moves, reducing to three rows, and splitting scalar count sectors from")
    print("  Johnson harmonic sectors.")
    print("  LRC translation used here:")
    print("    fixed margins        -> labelled packet signatures")
    print("    swap fiber           -> packet-preserving mutations")
    print("    scalar sector        -> qdiv / exact M / Haar-open status")
    print("    Johnson sectors      -> C27, K33, source-spectrum, boundary-moment")
    print("    three-row reduction  -> triples of proof obligations, not runner triples")
    print()


def print_tournament_analysis() -> None:
    print("[6] Tournament Analysis")
    names = list(LAYER_VECTORS)
    fp = tournament_fingerprint(names)
    print("  vertices are proof layers / packet sectors, not runners.")
    print("  pairwise observable:")
    print("    which layer preserves strict-counterexample status and owner labels")
    print("    before scalarization.")
    print("  tie Hamiltonian path:")
    print("    " + " -> ".join(names))
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} hp={fp['hp']}")
    print(f"  scc={fp['scc']}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=120)
    parser.add_argument("--two-swap-limit", type=int, default=30)
    parser.add_argument("--alias-depth", type=int, default=3)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    parser.add_argument("--sample-per-route", type=int, default=10)
    parser.add_argument("--denom-cap", type=int, default=140)
    parser.add_argument("--include-all-covering", action="store_true")
    parser.add_argument("--full-bank", action="store_true")
    parser.add_argument("--representative-limit", type=int, default=64)
    args = parser.parse_args()

    print("LRC14 boundary-moment packet ledger")
    print("=" * 78)
    print(
        f"single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}, "
        f"full_bank={args.full_bank}, "
        f"workers={args.workers}, denom_cap={args.denom_cap}"
    )
    print_assumption_challenge()

    packets = build_packets(args)
    chosen = selected_packets(packets, args.sample_per_route, args.include_all_covering)
    ledgers = [unit_moment_ledger(p, choose_denominator(p, args.denom_cap)) for p in chosen]

    print_ledger_summary(packets, ledgers)
    print_route_moment_stats(ledgers)
    print_representatives(ledgers, args.representative_limit)
    print_theorem_readout(ledgers)
    print_arxiv_bridge()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
