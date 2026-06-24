#!/usr/bin/env python3
"""S152: boundary-gap packet bridge for the LRC14 covering branch.

The current labelled-packet route says a strict LRC14 counterexample must
survive the q-witness gate, hence it must cover every divisor d<=14.  The
remaining danger is the F6 bucket from HYP-2956/HYP-2962/HYP-2963:

    qdiv(S)>14, strict safe open set empty, no petal/K33/state-lift label.

This script makes that bucket falsifiable by replacing the phrase "positive
Haar front" with exact boundary-gap packets.  A positive strict component is
an interval (lo, hi) between two adjacent danger-arc boundaries.  Its length is
an exact rational certificate for M(S)>1/14.  A zero-open covering kernel must
therefore pinch every such boundary bridge.

Tournament Analysis is included with vertices equal to proof coordinates, not
runners.  The selected quotient preserves the LRC witness predicate by keeping
q-covering, exact boundary gaps, endpoint owners, C27/K33 labels, and the
zero-open kernel obligation together.
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


s138 = load_module(
    "s152_s138_frontier",
    REPO / "04-computation" / "lrc14_c27_two_swap_frontier_codex_s138.py",
)
s147 = load_module(
    "s152_s147_haar",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)
s2963 = load_module(
    "s152_hyp2963_packets",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)


@dataclass(frozen=True)
class BoundaryBridge:
    lo: Fraction
    hi: Fraction
    length: Fraction
    left_end_owners: tuple[int, ...]
    right_start_owners: tuple[int, ...]

    @property
    def left_label(self) -> str:
        return owner_label(self.left_end_owners)

    @property
    def right_label(self) -> str:
        return owner_label(self.right_start_owners)

    @property
    def transition(self) -> str:
        return f"{self.left_label}->{self.right_label}"


@dataclass(frozen=True)
class RowBridgeAudit:
    name: str
    source_family: str
    speeds: tuple[int, ...]
    q_threshold: int
    safe_mu: Fraction
    component_count: int
    max_gap: Fraction
    min_gap: Fraction
    bridge_transitions: tuple[str, ...]
    current_l1: int
    owner_mod14_support: tuple[int, ...]
    covering_signature: str


PROOF_VERTICES = (
    "q_cover_gate",
    "exact_boundary_gap",
    "endpoint_owner_current",
    "source_spectrum_pullback",
    "C27_unital_packet",
    "K33_state_lift_packet",
    "gK8_Ly_moment_image",
    "fixed_margin_swap_fiber",
    "raw_runner_row",
)

PROOF_SCORES = {
    "q_cover_gate": (7, 4, 2, 2, 2, 4),
    "exact_boundary_gap": (7, 7, 5, 4, 4, 6),
    "endpoint_owner_current": (6, 6, 7, 5, 5, 6),
    "source_spectrum_pullback": (7, 6, 7, 7, 6, 7),
    "C27_unital_packet": (5, 5, 7, 6, 6, 6),
    "K33_state_lift_packet": (5, 5, 7, 7, 7, 7),
    "gK8_Ly_moment_image": (6, 6, 6, 6, 7, 7),
    "fixed_margin_swap_fiber": (6, 5, 6, 6, 5, 7),
    "raw_runner_row": (2, 2, 2, 1, 1, 1),
}

TIE_PATH = (
    "source_spectrum_pullback",
    "exact_boundary_gap",
    "endpoint_owner_current",
    "gK8_Ly_moment_image",
    "K33_state_lift_packet",
    "C27_unital_packet",
    "fixed_margin_swap_fiber",
    "q_cover_gate",
    "raw_runner_row",
)


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def owner_label(owners: tuple[int, ...]) -> str:
    if not owners:
        return "circle"
    parts = []
    for v in owners:
        parts.append(f"{v % 14}:g{gcd(v, 14)}")
    return ",".join(parts)


def endpoint_events(speeds: tuple[int, ...]) -> dict[Fraction, dict[str, list[int]]]:
    events: dict[Fraction, dict[str, list[int]]] = defaultdict(lambda: {"start": [], "end": []})
    for v in speeds:
        for k in range(v):
            center = Fraction(k, v)
            start = mod1(center - THRESHOLD / v)
            end = mod1(center + THRESHOLD / v)
            events[start]["start"].append(v)
            events[end]["end"].append(v)
    for data in events.values():
        data["start"].sort()
        data["end"].sort()
    return events


def boundary_bridges(speeds: tuple[int, ...]) -> tuple[Fraction, tuple[BoundaryBridge, ...]]:
    data = s147.exact_row_measure(speeds)
    events = endpoint_events(speeds)
    bridges: list[BoundaryBridge] = []
    for lo, hi in data["safe_components"]:
        key_lo = mod1(lo)
        key_hi = mod1(hi)
        left = tuple(events.get(key_lo, {}).get("end", []))
        right = tuple(events.get(key_hi, {}).get("start", []))
        # Split intervals can put a component at 0 or 1.  If the directional
        # owner is absent, keep the whole endpoint owner packet instead of
        # pretending the owner is unknown.
        if not left:
            both = events.get(key_lo, {})
            left = tuple(sorted(set(both.get("start", []) + both.get("end", []))))
        if not right:
            both = events.get(key_hi, {})
            right = tuple(sorted(set(both.get("start", []) + both.get("end", []))))
        bridges.append(
            BoundaryBridge(
                lo=lo,
                hi=hi,
                length=hi - lo,
                left_end_owners=left,
                right_start_owners=right,
            )
        )
    return data["safe_measure"], tuple(bridges)


def current_l1(bridges: tuple[BoundaryBridge, ...]) -> int:
    current: Counter[tuple[int, int]] = Counter()
    for bridge in bridges:
        for v in bridge.left_end_owners:
            current[(v % 14, gcd(v, 14))] += 1
        for v in bridge.right_start_owners:
            current[(v % 14, gcd(v, 14))] -= 1
    return sum(abs(value) for value in current.values())


def owner_mod14_support(bridges: tuple[BoundaryBridge, ...]) -> tuple[int, ...]:
    support = set()
    for bridge in bridges:
        support.update(v % 14 for v in bridge.left_end_owners)
        support.update(v % 14 for v in bridge.right_start_owners)
    return tuple(sorted(support))


def covering_signature(speeds: tuple[int, ...]) -> str:
    hits = []
    for d in range(2, 15):
        owners = tuple(v for v in speeds if v % d == 0)
        hits.append(f"{d}:{len(owners)}")
    return " ".join(hits)


def audit_row(name: str, source_family: str, speeds: tuple[int, ...]) -> RowBridgeAudit:
    q = s138.s124.q_threshold(speeds)
    safe_mu, bridges = boundary_bridges(speeds)
    lengths = tuple(bridge.length for bridge in bridges)
    return RowBridgeAudit(
        name=name,
        source_family=source_family,
        speeds=speeds,
        q_threshold=q,
        safe_mu=safe_mu,
        component_count=len(bridges),
        max_gap=max(lengths, default=Fraction(0)),
        min_gap=min(lengths, default=Fraction(0)),
        bridge_transitions=tuple(sorted(Counter(bridge.transition for bridge in bridges))),
        current_l1=current_l1(bridges),
        owner_mod14_support=owner_mod14_support(bridges),
        covering_signature=covering_signature(speeds),
    )


def tournament_mask(names: tuple[str, ...]) -> int:
    tie_rank = {name: len(TIE_PATH) - i for i, name in enumerate(TIE_PATH)}
    mask = 0
    bit = 0
    for i, j in combinations(range(len(names)), 2):
        a = names[i]
        b = names[j]
        if PROOF_SCORES[a] > PROOF_SCORES[b]:
            orient = True
        elif PROOF_SCORES[a] < PROOF_SCORES[b]:
            orient = False
        else:
            orient = tie_rank[a] >= tie_rank[b]
        if orient:
            mask |= 1 << bit
        bit += 1
    return mask


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, danger arcs, safe components, component endpoints, endpoint")
    print("    owner events, divisor clocks, C27 shells, K33 incidence, exact-period")
    print("    residues, Fourier/moment modes, fixed-margin fibers, and proof obligations.")
    print("  chosen vertices:")
    print("    proof coordinates of a covering boundary-gap packet.")
    print("  preserved LRC predicate:")
    print("    qdiv>14 candidate status plus exact strict safe intervals; any positive")
    print("    bridge length is a rigorous witness that M(S)>1/14.")
    print("  challenged assumption:")
    print("    the covering branch is not a raw row family; it is an endpoint-collision")
    print("    problem.  A counterexample must pinch every labelled boundary bridge.")
    print()


def print_summary(audits: list[RowBridgeAudit]) -> None:
    zero_open = [a for a in audits if a.safe_mu == 0]
    positive = [a for a in audits if a.safe_mu > 0]
    net_zero = [a for a in audits if a.current_l1 == 0]
    transition_hist = Counter()
    mod14_hist = Counter()
    for audit in audits:
        transition_hist.update(audit.bridge_transitions)
        mod14_hist.update(audit.owner_mod14_support)
    print("[1] Covering boundary-gap audit")
    print(f"  qdiv>14 covering rows audited={len(audits)}")
    print(f"  positive strict-open rows={len(positive)}")
    print(f"  zero-open covering rows={len(zero_open)}")
    print(f"  rows with zero net endpoint current={len(net_zero)}")
    if positive:
        min_mu = min(positive, key=lambda a: (a.safe_mu, a.name))
        min_max_gap = min(positive, key=lambda a: (a.max_gap, a.name))
        print(f"  smallest safe_mu={fmt(min_mu.safe_mu)} at {min_mu.name}")
        print(f"  smallest max bridge={fmt(min_max_gap.max_gap)} at {min_max_gap.name}")
    print("  endpoint owner residues mod 14 used by bridge certificates:")
    print("    " + ", ".join(f"{k}:{v}" for k, v in sorted(mod14_hist.items())))
    print("  top transition labels:")
    for transition, count in transition_hist.most_common(12):
        print(f"    {transition:38s} {count}")
    print()


def print_tight_rows(audits: list[RowBridgeAudit], limit: int) -> None:
    print("[2] Tightest covering bridge certificates")
    rows = sorted(audits, key=lambda a: (a.safe_mu, a.max_gap, a.name))
    print(
        f"  {'name':38s} {'q':>3s} {'mu':>10s} {'comps':>5s} "
        f"{'min_gap':>10s} {'max_gap':>10s} {'I1':>4s} owner_mod14"
    )
    for audit in rows[:limit]:
        print(
            f"  {audit.name[:38]:38s} {audit.q_threshold:3d} "
            f"{fmt(audit.safe_mu):>10s} {audit.component_count:5d} "
            f"{fmt(audit.min_gap):>10s} {fmt(audit.max_gap):>10s} "
            f"{audit.current_l1:4d} {audit.owner_mod14_support}"
        )
        print(f"      q-cover={audit.covering_signature}")
        print(f"      transitions={'; '.join(audit.bridge_transitions[:6])}")
    print()


def print_named_components(audits: list[RowBridgeAudit], rows: list[tuple[str, str, tuple[int, ...]]]) -> None:
    names = {
        "divisor-loaded tail 84",
        "divisor-loaded tail 84*2",
        "AP repair drop13 add182",
        "covering repair drop13 add182",
        "covering comb 12->84",
        "single swap 12->84",
        "single swap 13->182",
    }
    lookup = {audit.name: audit for audit in audits}
    row_lookup = {name: speeds for name, _, speeds in rows}
    print("[3] Named bridge decompositions")
    shown = 0
    for name in sorted(names):
        if name not in lookup or name not in row_lookup:
            continue
        safe_mu, bridges = boundary_bridges(row_lookup[name])
        print(f"  {name}: safe_mu={fmt(safe_mu)} components={len(bridges)}")
        for bridge in sorted(bridges, key=lambda b: (-b.length, b.lo))[:8]:
            print(
                f"    [{fmt(bridge.lo)}, {fmt(bridge.hi)}] len={fmt(bridge.length)} "
                f"{bridge.transition}"
            )
        shown += 1
    if not shown:
        print("  no named rows from the display list were present in this bank")
    print()


def print_theorem_readout(audits: list[RowBridgeAudit]) -> None:
    zero_open = [a for a in audits if a.safe_mu == 0]
    net_zero = [a for a in audits if a.current_l1 == 0]
    print("[4] Boundary-gap packet theorem readout")
    print("  Exact lemma now available for any finite row S:")
    print("    S is LRC14-safe with strict slack iff at least one boundary bridge")
    print("    between adjacent danger arcs has positive rational length.")
    print("  Therefore a strict counterexample in the covering branch must satisfy:")
    print("    qdiv(S)>14 and every labelled boundary bridge is pinched to length 0")
    print("    after merging danger arcs.")
    print("  Audited conclusion:")
    print(f"    zero-open qdiv>14 packets in this bank = {len(zero_open)}.")
    print(f"    zero net first-current packets in this bank = {len(net_zero)}.")
    print("    Thus the bridge cannot be a raw divergence proof; the useful token is")
    print("    the localized transition packet / second moment of the boundary gaps.")
    print("  Proof target promoted from a search condition:")
    print("    show the localized endpoint-owner transition packet of any qdiv>14")
    print("    zero-open row either has positive gK8/L_y moment image or carries a")
    print("    K33/H=7 state-lift label.")
    print("  This is the boundary-moment bridge in packet form.")
    print()


def print_tournament_analysis() -> None:
    print("[5] Tournament Analysis")
    names = PROOF_VERTICES
    mask = tournament_mask(names)
    fp = s138.tournament_fingerprint(mask, len(names))
    print("  vertices are proof coordinates, not runners.")
    print("  pair observable:")
    print("    which coordinate preserves the covering LRC predicate and prevents a")
    print("    positive boundary gap from being scalarized away.")
    print("  switch/gauge:")
    print("    lexicographic retention of q-covering, exact gap, endpoint owner,")
    print("    source identity, state-lift visibility, moment image, and anti-scalar guard.")
    print("  tie Hamiltonian path:")
    print("    " + " > ".join(TIE_PATH))
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--top", type=int, default=14)
    args = parser.parse_args()

    print("S152 LRC14 boundary-gap packet bridge")
    print("=" * 78)
    print(
        f"single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}"
    )
    print_assumption_challenge()

    rows = s2963.build_bank(
        single_limit=args.single_limit,
        two_swap_limit=args.two_swap_limit,
        alias_depth=args.alias_depth,
        lcm_tail_max=args.lcm_tail_max,
    )
    covering_rows = [
        (name, family, speeds)
        for name, family, speeds in rows
        if s138.s124.q_threshold(speeds) > 14
    ]
    audits = [audit_row(name, family, speeds) for name, family, speeds in covering_rows]
    audits.sort(key=lambda a: (a.safe_mu, a.max_gap, a.name))

    print_summary(audits)
    print_tight_rows(audits, args.top)
    print_named_components(audits, covering_rows)
    print_theorem_readout(audits)
    print_tournament_analysis()


if __name__ == "__main__":
    main()
