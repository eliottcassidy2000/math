#!/usr/bin/env python3
"""S149 missing-picture bridge atlas for the LRC14 proof route.

This script is intentionally small and connective.  It reuses exact local
kernels from the existing repo history:

* S124 exact M/q-threshold/AP-Goddyn-Wong filters;
* S146 Haar-Baire strict-open and boundary-owner fronts.

The goal is not to reprove every predecessor.  The goal is to expose the
missing theorem in a computable ledger:

    reduced LRC14 residual
      -> exact M/Farey branch
      -> Haar-Baire open or boundary front
      -> C27/unital/K33 address
      -> discharge or TournamentStateLift.

Tournament Analysis declaration
-------------------------------
Vertices:
    proof interfaces/quotients, not runners.
Pairwise observable:
    branch retention, boundary retention, C27/unital retention,
    K33/state-lift fit, finite checkability, theorem endpoint strength,
    anti-scalarization guard.
Switch:
    componentwise score sum; ties follow the declared Hamiltonian path.
Tie Hamiltonian path:
    exact_M_Farey -> Haar_Baire_front -> boundary_owner_skeleton
    -> C27_transfer -> unital_pair_chart -> affine_depth
    -> K33_incidence -> TournamentStateLift -> raw_apex_tournament
    -> raw_scalar_M.

Challenged assumption:
    The right final vertices are not runners, arcs, or a single tournament
    iso class.  They are typed proof obligations that preserve the LRC witness
    predicate through the exact scale, measurable front, and owner address.
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
    "s149_s124_common",
    REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py",
)
s146 = load_module(
    "s149_s146_boundary",
    REPO / "04-computation" / "lrc14_haar_baire_taut_boundary_s146.py",
)


@dataclass(frozen=True)
class NamedRow:
    name: str
    speeds: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class BridgeAudit:
    name: str
    speeds: tuple[int, ...]
    q_threshold: int
    covering14: bool
    has14: bool
    M: Fraction
    arg_denoms: tuple[int, ...]
    strict_mass: Fraction
    open_components: int
    boundary_points: int
    transfer: str
    packet: str
    route: str
    missing_arrow: str


ROWS = [
    NamedRow("AP", AP, "base endpoint atom"),
    NamedRow("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "first derived endpoint atom"),
    NamedRow("near 12->36", tuple(list(range(1, 12)) + [13, 36]), "first K33/Farey child"),
    NamedRow("petal 10->20", tuple(list(range(1, 10)) + [11, 12, 13, 20]), "unit C27 petal"),
    NamedRow("petal 13->26", tuple(list(range(1, 13)) + [26]), "unit C27 petal"),
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
    NamedRow("residue liar 12->26", tuple(list(range(1, 12)) + [13, 26]), "AP residues but q=12"),
    NamedRow("magnitude liar 12->96", tuple(list(range(1, 12)) + [13, 96]), "AP residues and q=14 but loose"),
    NamedRow("covering comb 12->84", tuple(list(range(1, 12)) + [13, 84]), "14-covering strict row"),
]


COMPONENTS = {
    (12, 24): "GW",
    (12, 36): "K33",
    (10, 20): "P10",
    (13, 26): "P13",
}


def removed_added(speeds: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    s = set(speeds)
    return tuple(sorted(set(AP) - s)), tuple(sorted(s - set(AP)))


def packet_label(speeds: tuple[int, ...]) -> str:
    removed, added = removed_added(speeds)
    if not removed and not added:
        return "AP"
    remaining = list(added)
    parts: list[str] = []
    unknown: list[str] = []
    for r in removed:
        hits = [a for a in remaining if (r, a) in COMPONENTS]
        if hits:
            a = hits[0]
            remaining.remove(a)
            parts.append(COMPONENTS[(r, a)])
        elif remaining:
            a = remaining.pop(0)
            unknown.append(f"{r}->{a}")
        else:
            unknown.append(f"{r}->?")
    for a in remaining:
        unknown.append(f"?->{a}")
    if unknown:
        parts.extend(f"unknown({u})" for u in unknown)
    return "+".join(parts)


def covering14(speeds: tuple[int, ...]) -> bool:
    return all(any(v % d == 0 for v in speeds) for d in range(2, 15))


def has14(speeds: tuple[int, ...]) -> bool:
    return any(v % 14 == 0 for v in speeds)


def route_for(row: NamedRow, q: int, cov14: bool, h14: bool, M: Fraction, strict_mass: Fraction, packet: str) -> tuple[str, str]:
    if M < DELTA:
        return "COUNTEREXAMPLE-SHAPED", "no route; would refute LRC14"
    if M == DELTA and strict_mass == 0:
        if packet in {"AP", "GW"}:
            return "boundary atom", "prove global boundary-only residuals reduce to AP/GW"
        return "new boundary packet", "feed this exact packet to HYP-2908/THM-572"
    if h14 or cov14:
        return "covering strictness", "finish shell-collapse/comb-margin theorem for 14-multiple rows"
    if q < 14:
        return "q-witness loose", "q-threshold already explains looseness"
    if "K33" in packet:
        return "K33 state-lift obligation", "construct TournamentStateLift from this nonunit owner packet"
    if packet in {"P10", "P13", "P10+GW", "GW+P10"}:
        return "unit-petal discharge", "prove petal/two-block rigidity uniformly"
    if "unknown" in packet:
        return "magnitude-aware open front", "retain spectrum/Farey scale; raw apex tournament forgot magnitude"
    return "positive Haar open front", "propagate Haar-Baire front with retained owner labels"


def audit_row(row: NamedRow) -> BridgeAudit:
    M, pts = s124.M_exact(row.speeds)
    q = s124.q_threshold(row.speeds)
    comps = s146.safe_open_components(row.speeds)
    strict_mass = s146.interval_measure(comps)
    boundary = s146.threshold_safe_points(row.speeds)
    packet = packet_label(row.speeds)
    route, missing = route_for(row, q, covering14(row.speeds), has14(row.speeds), M, strict_mass, packet)
    return BridgeAudit(
        row.name,
        row.speeds,
        q,
        covering14(row.speeds),
        has14(row.speeds),
        M,
        tuple(sorted({t.denominator for t in pts})),
        strict_mass,
        len(comps),
        len(boundary),
        s146.transfer_label(row.speeds),
        packet,
        route,
        missing,
    )


def print_bridge_table() -> None:
    print("[1] Named-row bridge audit")
    audits = [audit_row(row) for row in ROWS]
    for a in audits:
        print(
            f"  {a.name:24s} q={a.q_threshold:2d} cov14={str(a.covering14):5s} "
            f"has14={str(a.has14):5s} M={str(a.M):>7s} argD={list(a.arg_denoms)!s:14s} "
            f"strict_mu={str(a.strict_mass):>8s} comps={a.open_components:2d} "
            f"closed={a.boundary_points:2d} packet={a.packet:18s} route={a.route}"
        )
        print(f"      transfer={a.transfer}")
        print(f"      missing_arrow={a.missing_arrow}")


def print_boundary_database() -> None:
    print("\n[2] Reused S146 local boundary database")
    one_boundary, one_positive, one_covered = s146.classify_rows(s146.one_swap_rows(160))
    two_boundary, two_positive, two_covered = s146.classify_rows(s146.two_swap_rows(40))
    print(
        "  one-swap add<=160: "
        f"boundary_only={len(one_boundary)} positive_open={len(one_positive)} covered={len(one_covered)}"
    )
    print("    boundary rows=" + ", ".join(row.name for row, _ in one_boundary))
    print(
        "  two-swap add<=40:  "
        f"boundary_only={len(two_boundary)} positive_open={len(two_positive)} covered={len(two_covered)}"
    )
    print("    smallest positive two-swap rows:")
    for mass, row in two_positive[:5]:
        print(f"      {row.name:22s} mass={mass} transfer={row.note}")


@dataclass(frozen=True)
class Interface:
    name: str
    vector: tuple[int, int, int, int, int, int, int]
    note: str


INTERFACES = [
    Interface("exact_M_Farey", (5, 3, 2, 2, 5, 4, 5), "binding scale and Farey branch"),
    Interface("Haar_Baire_front", (4, 5, 2, 2, 5, 3, 4), "open-vs-boundary measurable front"),
    Interface("boundary_owner_skeleton", (3, 5, 3, 2, 4, 3, 4), "endpoint owner pairs"),
    Interface("C27_transfer", (3, 3, 5, 3, 5, 3, 5), "hole/double owner-carry shell"),
    Interface("unital_pair_chart", (3, 2, 5, 4, 4, 3, 5), "branch-local pair completion"),
    Interface("affine_depth", (2, 2, 4, 4, 4, 3, 4), "order-sensitive depth-14 marker"),
    Interface("K33_incidence", (3, 2, 3, 5, 4, 4, 5), "minor-transitive nonunit packet"),
    Interface("TournamentStateLift", (2, 2, 3, 5, 3, 5, 5), "machine-checked contradiction endpoint once constructed"),
    Interface("raw_apex_tournament", (1, 1, 1, 1, 4, 1, 2), "residue order only; magnitude blind"),
    Interface("raw_scalar_M", (4, 1, 1, 1, 5, 2, 1), "useful value, unsafe quotient"),
]


def beats(a: Interface, b: Interface) -> bool:
    sa = sum(a.vector)
    sb = sum(b.vector)
    if sa != sb:
        return sa > sb
    return INTERFACES.index(a) < INTERFACES.index(b)


def directed_3cycles(vertices: list[Interface]) -> int:
    count = 0
    for a, b, c in combinations(vertices, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        ac = beats(a, c)
        cb = beats(c, b)
        ba = beats(b, a)
        if (ab and bc and ca) or (ac and cb and ba):
            count += 1
    return count


def print_interface_tournament() -> None:
    print("\n[3] Proof-interface tournament")
    criteria = (
        "branch_retention",
        "boundary_retention",
        "C27_unital_retention",
        "K33_state_lift_fit",
        "finite_checkability",
        "theorem_endpoint",
        "anti_scalar_guard",
    )
    print("  criteria=" + ", ".join(criteria))
    ranking = sorted(INTERFACES, key=lambda x: (sum(x.vector), -INTERFACES.index(x)), reverse=True)
    wins = {v.name: 0 for v in INTERFACES}
    for a, b in combinations(INTERFACES, 2):
        if beats(a, b):
            wins[a.name] += 1
        else:
            wins[b.name] += 1
    for v in ranking:
        print(f"  {v.name:24s} vector={v.vector} score={wins[v.name]} note={v.note}")
    print(f"  directed_3cycles={directed_3cycles(INTERFACES)}")
    print("  ranking=" + " > ".join(v.name for v in ranking))


def print_readout() -> None:
    print("\n[4] Missing-picture readout")
    print("  The repo history is no longer missing a scalar invariant; it is missing a functor.")
    print("  Needed functor:")
    print("    reduced residual -> (exact scale, Haar/Baire front, C27/unital/K33 address)")
    print("                       -> discharge or TournamentStateLift.")
    print("  Local evidence:")
    print("    boundary-only local rows are AP and GW; their owner skeletons agree.")
    print("    positive-open rows are safe, but K33-labelled rows remain state-lift obligations.")
    print("    raw apex tournament and raw scalar M each forget a load-bearing label.")
    print("  Sharp theorem target:")
    print("    every primitive residual that is not discharged by q-threshold or Haar-open")
    print("    strictness must either be AP/GW boundary-only or construct the HYP-2908")
    print("    TournamentStateLift packet.")


def main() -> None:
    print("LRC14 missing-picture bridge atlas (S149)")
    print("=" * 72)
    print_bridge_table()
    print_boundary_database()
    print_interface_tournament()
    print_readout()


if __name__ == "__main__":
    main()
