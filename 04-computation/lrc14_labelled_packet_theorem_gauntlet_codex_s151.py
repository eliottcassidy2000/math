#!/usr/bin/env python3
"""S151: labelled-packet theorem gauntlet for LRC14.

This script turns the S149/S150 "missing picture" into a classification
statement.  It does not prove LRC14.  It records a theorem-shaped partition of
all possible LRC14 counterexamples into labelled packet families, audits the
current sporadic seeds with the exact S124/S146 kernels, and imports the proof
pattern of arXiv:2606.22636 as a packet-reduction protocol:

  fixed-margin binary relation space
    -> local swaps / two-row heat-bath comparison
    -> three-row reduction
    -> scalar count sector + Johnson harmonic sectors.

For LRC14 the proposed translation is:

  scalar count sector     = qdiv, exact M/Farey node, Haar mass
  Johnson harmonic sector = boundary owners, C27/unital labels, K33 lift debt
  non-migration kernel    = zero-open packet surviving all reductions.

Tournament Analysis vertices are proof obligations and packet families, not
runners.  The preserved predicate is the observer-source LRC14 condition plus
the labels needed to know which proof endpoint a row belongs to.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
DELTA = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s124 = load_module(
    "s151_s124_common",
    REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py",
)
s146 = load_module(
    "s151_s146_boundary",
    REPO / "04-computation" / "lrc14_haar_baire_taut_boundary_s146.py",
)


COMPONENTS = {
    (12, 24): "GW",
    (12, 36): "K33",
    (10, 20): "P10",
    (13, 26): "P13",
}


@dataclass(frozen=True)
class NamedRow:
    name: str
    speeds: tuple[int, ...]
    role: str


@dataclass(frozen=True)
class PacketAudit:
    name: str
    family: str
    speeds: tuple[int, ...]
    q_threshold: int
    qdiv: int
    M: Fraction
    arg_denoms: tuple[int, ...]
    strict_mass: Fraction
    closed_points: int
    packet: str
    transfer: str
    counterexample_gate: str
    proof_obligation: str


@dataclass(frozen=True)
class Family:
    code: str
    name: str
    gate: str
    endpoint: str
    current_status: str


FAMILIES = [
    Family(
        "F0",
        "qdiv witness discharge",
        "qdiv(S)<14",
        "strict witness t=1/qdiv gives M(S)>1/14",
        "proved by qdiv gate",
    ),
    Family(
        "F1",
        "AP/GW boundary atoms",
        "qdiv=14, O(S)=empty, C(S) nonempty, AP/GW owner skeleton",
        "M(S)=1/14; not counterexamples",
        "finite evidence: AP and GW only in S146/S150 local atlas",
    ),
    Family(
        "F2",
        "unit-petal/two-block discharge",
        "positive Haar front with only P10/P13/GW unit labels",
        "prove petal/two-block rigidity uniformly",
        "known sporadics: P10, P13, P10+GW",
    ),
    Family(
        "F3",
        "K33 state-lift packet",
        "K33/nonunit owner survives after positive or zero-open filtering",
        "construct HYP-2908/THM-572 TournamentStateLift",
        "known positive seeds: 12->36, P10+K33",
    ),
    Family(
        "F4",
        "unlabelled q14 positive front",
        "qdiv=14 and O(S) has positive Haar mass without known petal/K33 key",
        "safe by Haar-open witness; labels still needed for induction",
        "many in S150; not counterexample-shaped",
    ),
    Family(
        "F5",
        "covering positive / boundary-moment strictness",
        "qdiv>14 or has a 14-multiple, with O(S) positive",
        "finish boundary-moment comb/shell strictness theorem",
        "many in S150; named seed 12->84",
    ),
    Family(
        "F6",
        "covering zero-open non-migration kernel",
        "qdiv>14, O(S)=empty, C(S)=empty or non-AP/GW zero-open kernel",
        "only remaining counterexample-shaped family",
        "not seen in S150; target of boundary-moment bridge",
    ),
    Family(
        "F7",
        "new Johnson-harmonic packet sector",
        "fixed-margin packet swaps reveal a non-scalar sector not in F1-F6",
        "direct falsifier for the labelled packet theorem",
        "conceptual guardrail imported from arXiv:2606.22636",
    ),
]


ROWS = [
    NamedRow("AP", AP, "base boundary atom"),
    NamedRow("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "first derived boundary atom"),
    NamedRow("near 12->36", tuple(list(range(1, 12)) + [13, 36]), "first K33 child"),
    NamedRow("petal 10->20", tuple(list(range(1, 10)) + [11, 12, 13, 20]), "unit petal"),
    NamedRow("petal 13->26", tuple(list(range(1, 13)) + [26]), "unit petal"),
    NamedRow(
        "splice 10,12->20,24",
        tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 24]),
        "unit petal plus GW",
    ),
    NamedRow(
        "splice 10,12->20,36",
        tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 36]),
        "unit petal plus K33",
    ),
    NamedRow("residue liar 12->26", tuple(list(range(1, 12)) + [13, 26]), "qdiv witness liar"),
    NamedRow("magnitude liar 12->96", tuple(list(range(1, 12)) + [13, 96]), "q14 positive front liar"),
    NamedRow("covering comb 12->84", tuple(list(range(1, 12)) + [13, 84]), "covering/shell strict seed"),
]


S150_PACKET_COUNTS = {
    "AP/GW-boundary-source": 2,
    "positive-K33-state-lift": 340,
    "positive-covering-off-apex": 16488,
    "positive-unit-petal-or-GW-strip": 9632,
    "positive-unlabelled-q14-front": 41906,
}


def qdiv(speeds: tuple[int, ...], cap: int = 120) -> int:
    for d in range(2, cap + 1):
        if all(v % d for v in speeds):
            return d
    return cap + 1


def removed_added(speeds: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    s = set(speeds)
    return tuple(sorted(set(AP) - s)), tuple(sorted(s - set(AP)))


def packet_label(speeds: tuple[int, ...]) -> str:
    removed, added = removed_added(speeds)
    if not removed and not added:
        return "AP"
    remaining = list(added)
    out: list[str] = []
    unknown: list[str] = []
    for r in removed:
        hit = next((a for a in remaining if (r, a) in COMPONENTS), None)
        if hit is not None:
            remaining.remove(hit)
            out.append(COMPONENTS[(r, hit)])
        elif remaining:
            unknown.append(f"{r}->{remaining.pop(0)}")
        else:
            unknown.append(f"{r}->?")
    unknown.extend(f"?->{a}" for a in remaining)
    out.extend(f"unknown({u})" for u in unknown)
    return "+".join(out)


def has14(speeds: tuple[int, ...]) -> bool:
    return any(v % 14 == 0 for v in speeds)


def family_for(qd: int, M: Fraction, mass: Fraction, closed_points: int, packet: str, speeds: tuple[int, ...]) -> tuple[str, str, str]:
    if M < DELTA:
        return "F6", "COUNTEREXAMPLE", "audit as actual LRC14 disproof candidate"
    if qd < 14:
        return "F0", "safe", "qdiv witness discharges"
    if mass == 0:
        if closed_points and packet in {"AP", "GW"}:
            return "F1", "tight-not-counterexample", "global boundary-source rigidity"
        if closed_points:
            return "F6", "zero-open-new-boundary", "non-AP/GW boundary source; route to state lift"
        return "F6", "covered-zero-open", "covering non-migration kernel"
    if "K33" in packet:
        return "F3", "safe-but-obligated", "retain K33 label for state-lift induction"
    if set(packet.split("+")) <= {"P10", "P13", "GW"} and "unknown" not in packet:
        return "F2", "safe", "unit petal/two-block discharge"
    if qd > 14 or has14(speeds):
        return "F5", "safe-but-wide", "boundary-moment comb/shell strictness"
    return "F4", "safe", "positive Haar-open front"


def audit_row(row: NamedRow) -> PacketAudit:
    M, pts = s124.M_exact(row.speeds)
    q = s124.q_threshold(row.speeds)
    qd = qdiv(row.speeds)
    comps = s146.safe_open_components(row.speeds)
    mass = s146.interval_measure(comps)
    closed = s146.threshold_safe_points(row.speeds)
    packet = packet_label(row.speeds)
    family, gate, obligation = family_for(qd, M, mass, len(closed), packet, row.speeds)
    return PacketAudit(
        row.name,
        family,
        row.speeds,
        q,
        qd,
        M,
        tuple(sorted({t.denominator for t in pts})),
        mass,
        len(closed),
        packet,
        s146.transfer_label(row.speeds),
        gate,
        obligation,
    )


def print_family_table() -> None:
    print("[1] Labelled packet families for any LRC14 counterexample attempt")
    for fam in FAMILIES:
        print(f"  {fam.code}: {fam.name}")
        print(f"      gate:     {fam.gate}")
        print(f"      endpoint: {fam.endpoint}")
        print(f"      status:   {fam.current_status}")
    print()
    print("  Consequence:")
    print("    An actual counterexample M(S)<1/14 cannot lie in F0-F5.")
    print("    It must be F6, or else F7 exposes a missing packet sector.")


def print_sporadic_audit() -> Counter[str]:
    print("\n[2] Exact sporadic-seed audit")
    counts: Counter[str] = Counter()
    for row in ROWS:
        a = audit_row(row)
        counts[a.family] += 1
        print(
            f"  {a.name:24s} family={a.family:2s} qdiv={a.qdiv:2d} q={a.q_threshold:2d} "
            f"M={str(a.M):>7s} argD={list(a.arg_denoms)!s:12s} "
            f"mu={str(a.strict_mass):>10s} closed={a.closed_points:2d} "
            f"packet={a.packet:18s} gate={a.counterexample_gate}"
        )
        print(f"      transfer={a.transfer}")
        print(f"      obligation={a.proof_obligation}")
    return counts


def print_s150_gauntlet_import() -> None:
    print("\n[3] S150 gauntlet import")
    print("  audited AP-neighborhood banks:")
    print("    one-swap add<=420, two-swap add<=60, three-swap add<=30")
    print("  hard-row verdict:")
    print("    covered qdiv>=14 rows:         0")
    print("    non-AP/GW boundary-only rows:  0")
    print("  packet histogram:")
    for key, value in S150_PACKET_COUNTS.items():
        print(f"    {key:36s} {value}")


def print_fixed_margin_protocol() -> None:
    print("\n[4] Fixed-margin / Johnson-sector protocol from arXiv:2606.22636")
    print("  Imported proof pattern:")
    print("    binary relation with fixed margins")
    print("      -> local 2x2 swaps / two-row heat-bath comparison")
    print("      -> reduction to three rows")
    print("      -> scalar count sector + Johnson harmonic sectors")
    print("  LRC14 translation:")
    print("    binary relation = speeds versus packet lenses")
    print("      lenses: q clocks, boundary owners, C27 shells, unital pairs,")
    print("              K33 owners, exact-period packets, boundary-moment signs")
    print("    fixed margins = reductions that preserve source-spectrum labels")
    print("    scalar sector = qdiv, exact M/Farey node, Haar mass")
    print("    Johnson sectors = non-scalar owner/carry/state-lift labels")
    print("  Proof obligation:")
    print("    show every fixed-margin packet component reaches F0-F5 by swaps,")
    print("    unless a three-row Johnson sector is a genuine F7 falsifier.")


@dataclass(frozen=True)
class ProofVertex:
    name: str
    vector: tuple[int, int, int, int, int, int]


VERTICES = [
    ProofVertex("F6_boundary_moment_kernel", (6, 6, 5, 5, 6, 6)),
    ProofVertex("three_row_Johnson_sectors", (6, 5, 6, 4, 5, 6)),
    ProofVertex("fixed_margin_packet_swaps", (5, 5, 5, 5, 5, 5)),
    ProofVertex("boundary_owner_rigidity", (5, 6, 4, 4, 5, 5)),
    ProofVertex("K33_state_lift_endpoint", (4, 4, 5, 6, 5, 6)),
    ProofVertex("unit_petal_discharge", (4, 4, 4, 3, 5, 4)),
    ProofVertex("qdiv_scalar_gate", (5, 3, 2, 2, 6, 3)),
    ProofVertex("raw_row_or_scalar", (1, 1, 1, 1, 1, 1)),
]


def beats(a: ProofVertex, b: ProofVertex) -> bool:
    awins = sum(x > y for x, y in zip(a.vector, b.vector))
    bwins = sum(y > x for x, y in zip(a.vector, b.vector))
    if awins != bwins:
        return awins > bwins
    return VERTICES.index(a) < VERTICES.index(b)


def tournament_fingerprint() -> tuple[Counter[int], int, list[str]]:
    wins = Counter()
    c3 = 0
    for a, b in combinations(VERTICES, 2):
        if beats(a, b):
            wins[a.name] += 1
        else:
            wins[b.name] += 1
    for v in VERTICES:
        wins[v.name] += 0
    for a, b, c in combinations(VERTICES, 3):
        if (beats(a, b) and beats(b, c) and beats(c, a)) or (
            beats(a, c) and beats(c, b) and beats(b, a)
        ):
            c3 += 1
    order = [v.name for v in sorted(VERTICES, key=lambda v: (-wins[v.name], VERTICES.index(v)))]
    return Counter(wins.values()), c3, order


def print_tournament_analysis(family_counts: Counter[str]) -> None:
    print("\n[5] Tournament Analysis")
    print("  vertices:")
    print("    packet families and proof obligations, not runners or arcs.")
    print("  pairwise observable:")
    print("    which object preserves the LRC observer-source predicate and the")
    print("    labels needed to route a non-migrating packet.")
    print("  switch/gauge:")
    print("    majority over source retention, boundary retention, non-scalar")
    print("    sector visibility, endpoint strength, finite auditability, and")
    print("    anti-scalarization.")
    hist, c3, order = tournament_fingerprint()
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={c3}")
    print("  hamiltonian_path=" + " > ".join(order))
    print("  exact sporadic family counts=" + str(dict(sorted(family_counts.items()))))


def print_readout() -> None:
    print("\n[6] Readout")
    print("  Labelled Packet Theorem target:")
    print("    Every primitive LRC14 row is in F0-F5, or is an F6 non-migrating")
    print("    boundary-moment kernel, unless a new F7 Johnson sector exists.")
    print("  Therefore an actual counterexample must be:")
    print("    qdiv>=14, no strict Haar-open witness, not AP/GW boundary-only,")
    print("    not unit-petal, not positive K33, and zero/nonpositive under the")
    print("    boundary-moment bridge.")
    print("  The current gauntlets find no such packet in the audited AP atlas.")


def main() -> None:
    print("S151 LRC14 labelled-packet theorem gauntlet")
    print("=" * 78)
    print_family_table()
    counts = print_sporadic_audit()
    print_s150_gauntlet_import()
    print_fixed_margin_protocol()
    print_tournament_analysis(counts)
    print_readout()


if __name__ == "__main__":
    main()
