#!/usr/bin/env python3
"""S145: measurable rank recombination packet for the LRC14 low frontier.

This script is a proof-interface experiment for HYP-2947.  It reuses the exact
S138 AP two-swap frontier and then attaches the newer POKE labels:

  * exact M/Farey branch;
  * C27 shell-transfer component;
  * q=3 unital local-chart role;
  * affine-depth/order signature;
  * K33/Kuratowski state-lift flag;
  * Borel/Haar-safe witness-address status;
  * PH-style bad-child rank placeholder.

The aim is not to prove LRC14 outright.  The aim is to test whether the known
M<=2/27 local frontier already recombines into the desired three-way split:

  tight AP/GW | unit-petal discharge | K33/state-lift obligation.

Tournament Analysis vertices are proof interfaces, not runners or arcs.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path
import argparse
import os
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
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
    "s145_s138_frontier",
    REPO / "04-computation" / "lrc14_c27_two_swap_frontier_codex_s138.py",
)
s142 = load_module(
    "s145_s142_affine",
    REPO / "04-computation" / "lrc14_affine_depth_unital_chain_codex_s142.py",
)


@dataclass(frozen=True)
class PacketComponent:
    key: str
    label: str
    removed: int
    added: int
    hole_shell: int
    double_shell: int
    depth: int
    unital_block: tuple[str, ...]
    role: str
    measurable_address: str


COMPONENTS: dict[tuple[int, int], PacketComponent] = {
    (12, 24): PacketComponent(
        "GW",
        "GW 12->24",
        12,
        24,
        12,
        3,
        s142.COMPONENTS["GW"].depth,
        s142.COMPONENTS["GW"].block,
        "tight AP/GW transfer",
        "Borel-coded after AP/GW owner labels",
    ),
    (12, 36): PacketComponent(
        "K33",
        "near-miss 12->36",
        12,
        36,
        12,
        9,
        s142.COMPONENTS["K33"].depth,
        s142.COMPONENTS["K33"].block,
        "nonunit K33/Kuratowski obligation",
        "requires nonunit owner/K33 address",
    ),
    (10, 20): PacketComponent(
        "P10",
        "petal 10->20",
        10,
        20,
        10,
        7,
        s142.COMPONENTS["P10"].depth,
        s142.COMPONENTS["P10"].block,
        "unit C27 petal discharge",
        "unit-petal Borel address",
    ),
    (13, 26): PacketComponent(
        "P13",
        "petal 13->26",
        13,
        26,
        13,
        1,
        s142.COMPONENTS["P13"].depth,
        s142.COMPONENTS["P13"].block,
        "unit C27 petal discharge",
        "unit-petal Borel address",
    ),
}


@dataclass(frozen=True)
class PacketAudit:
    row: object
    removed: tuple[int, ...]
    added: tuple[int, ...]
    components: tuple[PacketComponent, ...]
    unknown_pairs: tuple[tuple[int, int], ...]
    route: str
    measurable_code: str
    ph_rank: int
    state_lift: bool
    depth_sum: int
    suffixes: tuple[int, ...]


def removed_added(speeds: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    s = set(speeds)
    return tuple(sorted(set(AP) - s)), tuple(sorted(s - set(AP)))


def match_components(removed: tuple[int, ...], added: tuple[int, ...]) -> tuple[tuple[PacketComponent, ...], tuple[tuple[int, int], ...]]:
    """Match known low-frontier replacement atoms.

    The low-frontier replacements are rigid in the S136/S138 bank.  For
    unfamiliar rows we greedily record unknown removed/added pairs so the
    script fails visibly rather than hiding a new atom.
    """

    remaining = list(added)
    comps: list[PacketComponent] = []
    unknown: list[tuple[int, int]] = []
    for r in removed:
        hits = [a for a in remaining if (r, a) in COMPONENTS]
        if hits:
            a = hits[0]
            remaining.remove(a)
            comps.append(COMPONENTS[(r, a)])
        elif remaining:
            a = remaining.pop(0)
            unknown.append((r, a))
        else:
            unknown.append((r, -1))
    for a in remaining:
        unknown.append((-1, a))
    return tuple(comps), tuple(unknown)


def affine_profile(keys: tuple[str, ...]) -> tuple[int, tuple[int, ...]]:
    if not keys:
        return 0, ()
    depths = [s142.COMPONENTS[k].depth for k in keys]
    _, suffixes, _, depth_sum = s142.affine_depth_profile(depths)
    return depth_sum, tuple(suffixes)


def classify_packet(row: object) -> PacketAudit:
    removed, added = removed_added(row.speeds)
    comps, unknown = match_components(removed, added)
    keys = tuple(c.key for c in comps)
    depth_sum, suffixes = affine_profile(keys)

    if row.label == "AP":
        route = "tight AP floor"
        code = "Haar-safe base event"
        rank = 0
        state_lift = False
    elif unknown:
        route = "unknown low-frontier atom"
        code = "not Borel-safe until labelled"
        rank = 99
        state_lift = False
    elif keys == ("GW",):
        route = "tight GW floor"
        code = "AP/GW-calibrated Borel address"
        rank = 0
        state_lift = False
    elif set(keys) <= {"P10", "P13"}:
        route = "unit-petal discharge"
        code = "unit C27 petal address"
        rank = 0
        state_lift = False
    elif keys == ("K33",):
        route = "K33 state-lift obligation"
        code = "nonunit K33/Kuratowski address"
        rank = 1
        state_lift = True
    elif set(keys) == {"P10", "GW"}:
        route = "petal-plus-tight-GW discharge strip"
        code = "product Borel code: unit petal x AP/GW"
        rank = 0
        state_lift = False
    elif set(keys) == {"P10", "K33"}:
        route = "petal-plus-K33 state-lift strip"
        code = "product Borel code: unit petal x K33"
        rank = 1
        state_lift = True
    else:
        route = "labelled but unproved recombination"
        code = "Borel-coded, proof route open"
        rank = sum(1 for k in keys if k == "K33")
        state_lift = rank > 0

    return PacketAudit(row, removed, added, comps, unknown, route, code, rank, state_lift, depth_sum, suffixes)


def compute_frontier(limit: int, workers: int, progress_every: int) -> list[object]:
    rows = s138.s124.bank_two_swaps(limit)
    q14 = [S for S in rows if s138.s124.q_threshold(S) >= 14]
    ap_shells, _ = s138.shell_profile(AP)
    ap_mask = s138.tournament_mask([s138.shell_priority(s) for s in ap_shells])
    exact = s138.compute_exact(q14, workers=workers, progress_every=progress_every)
    audits = [
        s138.audit_row(S, M, denoms, ap_mask)
        for S, (M, denoms) in exact.items()
    ]
    return sorted(audits, key=lambda a: (a.M, a.label, a.speeds))


def component_table(packets: list[PacketAudit]) -> None:
    print("[2] Measurable recombination table for M<=2/27")
    print(f"  rows={len(packets)}")
    for p in packets:
        comp = ",".join(c.key for c in p.components) or "-"
        print(
            f"  {p.row.label:32s} M={str(p.row.M):>5s} "
            f"branch={p.row.branch:21s} comps={comp:8s} "
            f"depth_sum={p.depth_sum:2d} suffixes={list(p.suffixes)!s:10s} "
            f"rank={p.ph_rank:2d} lift={str(p.state_lift):5s} route={p.route}"
        )
        print(f"      code={p.measurable_code}")
        print(f"      transfer={p.row.transfer}")


def virtual_depth14_chain() -> None:
    print()
    print("[3] Virtual linked-chain depth-14 audit")
    keys = ("GW", "K33", "P10")
    print("  components: GW -> K33 -> P10")
    report = s142.sequence_report(keys)
    print(f"  depths={report['depths']} suffixes={report['suffixes']} depth_sum={report['depth_sum']} beta={report['beta']}")
    print("  all permutations of the same component multiset:")
    for perm in sorted(set(permutations(keys)), key=lambda x: s142.sequence_report(x)["depth_sum"]):
        r = s142.sequence_report(perm)
        mark = " <-- unique LRC14 depth" if r["depth_sum"] == 14 else ""
        print(f"    {perm}: depth_sum={r['depth_sum']:2d} suffixes={r['suffixes']}{mark}")
    print("  readout:")
    print("    depth 14 is not a scalar equality.  It appears only when the")
    print("    AP/GW block is followed by the nonunit K33 block and then by the")
    print("    unit petal discharge address.  This is the order-sensitive state-lift")
    print("    packet that HYP-2947 wants to construct from any surviving bad atom.")


def quotient_audit() -> None:
    print()
    print("[4] Haar/Borel quotient audit")
    rows = [
        ("global phase", "compact group translation on T=R/Z", "safe"),
        ("C27 shell labels", "finite labelled residue packet", "keep label"),
        ("q=3 unital chart", "branch-local pair-completion chart", "keep chart"),
        ("Kpq/K33 minor flag", "graph-minor predicate, not a measure quotient", "keep flag"),
        ("owner-private deletion", "HYP-2248 address tax", "mandatory"),
        ("raw row scalar M only", "forgets witness address", "unsafe alone"),
    ]
    for name, meaning, verdict in rows:
        print(f"  {name:24s} | {verdict:12s} | {meaning}")
    print("  preserved predicate:")
    print("    positive Haar measure of GOOD ∩ G_P after exact packet labels survive.")
    print("  destroyed only after proof:")
    print("    runner identity and raw residue labels, once a measurable invariant")
    print("    action proves they do not affect the witness event.")


def proof_channel_tournament() -> None:
    print()
    print("[5] Tournament Analysis: measurable proof interfaces")
    channels = [
        ("exact M/Farey branch", (6, 6, 6, 5, 5, 5)),
        ("Haar/Borel witness carrier", (6, 5, 6, 6, 5, 5)),
        ("C27 owner/carry shell code", (5, 6, 5, 5, 5, 5)),
        ("q=3 unital local chart", (5, 5, 5, 5, 5, 5)),
        ("affine order/depth code", (5, 5, 5, 4, 5, 5)),
        ("K33/Kuratowski state-lift flag", (5, 4, 5, 6, 4, 5)),
        ("PH bad-child rank", (4, 5, 5, 5, 6, 5)),
        ("raw scalar count", (1, 1, 1, 1, 1, 1)),
    ]
    mask = s138.tournament_mask([score for _, score in channels])
    fp = s138.tournament_fingerprint(mask, len(channels))
    order = sorted(range(len(channels)), key=lambda i: channels[i][1], reverse=True)
    print("  vertices are proof interfaces, not runners.")
    print("  pair observable:")
    print("    branch retention, address retention, measurability, state-lift fit,")
    print("    extension-rank control, anti-scalar guard.")
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(channels[i][0] for i in order))


def proof_attempt_summary(packets: list[PacketAudit], all_frontier: list[object]) -> None:
    print()
    print("[6] Proof-attempt verdict")
    route_counts = Counter(p.route for p in packets)
    rank_counts = Counter(p.ph_rank for p in packets)
    unknown = [p for p in packets if p.unknown_pairs or p.ph_rank >= 99]
    lifts = [p for p in packets if p.state_lift]
    print(f"  M<=3/41 rows in scan: {sum(1 for a in all_frontier if a.M <= LOW_341)}")
    print(f"  M<=2/27 rows in scan: {len(packets)}")
    print(f"  route_counts={dict(sorted(route_counts.items()))}")
    print(f"  rank_counts={dict(sorted(rank_counts.items()))}")
    print(f"  state_lift_obligations={[p.row.label for p in lifts]}")
    print(f"  unknown_low_frontier_atoms={len(unknown)}")
    if unknown:
        for p in unknown:
            print(f"    UNKNOWN {p.row.label}: removed={p.removed} added={p.added} unknown={p.unknown_pairs}")
    else:
        print("  No unknown atom appears in the exact M<=2/27 two-swap frontier.")
    print()
    print("  What this would prove if globalized:")
    print("    after AP-tail and q-threshold reductions, any low-frontier bad atom")
    print("    recombines into either AP/GW tightness, a unit-petal discharge, or")
    print("    a K33/state-lift obligation.  THM-572 already closes the endpoint")
    print("    once the K33 packet is functorially lifted into tournament H=7.")
    print("  Remaining hard theorem:")
    print("    replace the finite two-swap ceiling by a structural theorem that every")
    print("    primitive LRC14 counterexample reaches this measurable packet language.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=40)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    parser.add_argument("--progress-every", type=int, default=0)
    args = parser.parse_args()

    print("S145 LRC14 MEASURABLE RANK RECOMBINATION")
    print("=" * 78)
    print(f"limit={args.limit}, workers={args.workers}")
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, residues, cover arcs, C27 shells, unital blocks, Kpq minors,")
    print("    affine words, owner-private addresses, Borel event codes, and proof obligations.")
    print("  chosen vertices:")
    print("    measurable proof-interface packets on the exact low frontier.")
    print("  preserved LRC predicate:")
    print("    positive Haar measure of GOOD ∩ G_P after exact M/Farey and packet labels.")
    print("  challenged assumption:")
    print("    a scalar low value M or a raw shell transfer is enough.  It is not;")
    print("    the witness address and state-lift route must be retained.")

    print()
    print("[1] Exact frontier recomputation")
    audits = compute_frontier(args.limit, args.workers, args.progress_every)
    packets = [classify_packet(a) for a in audits if a.M <= LOW_227]
    print(f"  q>=14 exact rows audited={len(audits)}")
    print(f"  M<=3/41 rows={sum(1 for a in audits if a.M <= LOW_341)}")
    print(f"  M<=2/27 rows={len(packets)}")

    component_table(packets)
    virtual_depth14_chain()
    quotient_audit()
    proof_channel_tournament()
    proof_attempt_summary(packets, audits)


if __name__ == "__main__":
    main()
