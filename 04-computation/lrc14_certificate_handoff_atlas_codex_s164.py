#!/usr/bin/env python3
"""S164: LRC14 certificate handoff atlas.

This is a creative proof-pass artifact, not a proof.  The recent LRC14 work has
many strong local certificates: q-witnesses, Haar/Baire open fronts, twist
ladders, danger-count duals, Ramanujan exact-period packets, Fejer interval
certificates, analytic-sieve/Kaczynski packets, and the HYP-2908 state-lift
endpoint.  The missing theorem is a gluing theorem: a packet may move from one
certificate carrier to another only if the LRC predicate and all load-bearing
labels survive the quotient.

The script records that handoff atlas explicitly.  Tournament Analysis vertices
are proof obligations / certificate carriers, not runners.  The pairwise
observable is the retention vector:

    LRC predicate, exact scale, phase/period, topology, endpoint owners,
    packet family, dual certificate, formal checkability, residual routing.

Ties follow the listed Hamiltonian path.  The output is intended to make the
next proof pass more surgical: prove the open arrows rather than searching for
another scalar invariant.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: tuple[int, ...]
    role: str
    may_forget: tuple[str, ...]
    must_retain: tuple[str, ...]


@dataclass(frozen=True)
class PacketHandoff:
    packet: str
    current_carrier: str
    exit_carrier: str
    status: str
    retained_predicate: str
    destroyed_without_guard: tuple[str, ...]
    open_arrow: str


CARRIERS = [
    Carrier(
        "labelled_source_packet",
        (5, 5, 5, 5, 5, 5, 4, 4, 5),
        "root sheaf: qdiv/Farey/Haar/endpoint/C27/K33/state labels",
        ("raw runner order after packet key is fixed",),
        ("qdiv", "exact_M_Farey", "strict_vs_boundary", "endpoint_owners", "route"),
    ),
    Carrier(
        "exact_interval_fejer",
        (5, 5, 4, 4, 3, 5, 5, 5, 4),
        "formal dual certificate target: Q_d(center)<0 on packet fibers",
        ("floating decimal shadow", "raw endpoint list after center/component is fixed"),
        ("packet_key", "center", "degree", "atom_formula", "interval_upper"),
    ),
    Carrier(
        "ramanujan_exact_period",
        (4, 4, 5, 3, 2, 4, 4, 4, 4),
        "primitive phase projector and late-q pre-splitter",
        ("speed names inside a primitive residue profile",),
        ("q", "unit_residue", "weak_vs_strict_witness", "route_label"),
    ),
    Carrier(
        "endpoint_bridge_graph",
        (4, 4, 3, 5, 5, 4, 4, 4, 4),
        "taut/positive bridge layer over danger endpoints",
        ("global scalar safe mass after bridges are labelled",),
        ("left_owner", "right_owner", "bridge_sign", "winding_or_potential"),
    ),
    Carrier(
        "twist_ladder_blocker",
        (4, 4, 4, 3, 2, 3, 5, 5, 4),
        "primal rational witness or finite speed-vs-twist blocker",
        ("continuous endpoint geometry after a twist witness is found",),
        ("denominator_ladder", "failed_twist_hypergraph", "first_q"),
    ),
    Carrier(
        "danger_count_moment_dual",
        (4, 3, 2, 4, 1, 3, 5, 5, 3),
        "count-distribution Farkas/Delsarte dual",
        ("endpoint ownership", "phase position"),
        ("count_histogram", "dual_polynomial", "degree", "sign"),
    ),
    Carrier(
        "analytic_sieve_kaczynski",
        (3, 4, 4, 3, 2, 4, 3, 3, 5),
        "large-sieve/circle-method/Kaczynski boundary middle certificate",
        ("raw prime or Mobius scalar",),
        ("smoothing_kernel", "exceptional_set", "approach_class", "exact_period_labels"),
    ),
    Carrier(
        "tournament_state_lift",
        (5, 3, 2, 2, 4, 5, 4, 5, 5),
        "terminal forbidden-H=7 state-lift endpoint",
        ("metric values after the conflict packet is built",),
        ("packet_value", "H(T)=7", "conflict_tournament", "realizability_map"),
    ),
    Carrier(
        "raw_scalar_shadow",
        (1, 1, 1, 1, 1, 1, 1, 1, 1),
        "unsafe summary: tau/sigma/raw mass/raw tournament class",
        tuple(),
        ("lost_label_certificate",),
    ),
]

TIE_PATH = [carrier.name for carrier in CARRIERS]


HANDOFFS = [
    PacketHandoff(
        "qdiv<14 / direct q-witness",
        "labelled_source_packet",
        "twist_ladder_blocker",
        "closed in the current grammar",
        "existence of a strict rational witness",
        ("endpoint-owner route",),
        "global reduction must put every direct-witness packet here before scalarizing",
    ),
    PacketHandoff(
        "AP/Goddyn-Wong boundary atom",
        "endpoint_bridge_graph",
        "labelled_source_packet",
        "equality atom, not a strict counterexample",
        "zero strict-open mass with q=14 boundary owner skeleton",
        ("C27/GW hidden acceleration", "closed boundary support"),
        "tight-locus census still must show no other global equality atom",
    ),
    PacketHandoff(
        "K33 / 12->36 and K33 splices",
        "labelled_source_packet",
        "exact_interval_fejer",
        "selected packet certified negative; state-lift endpoint remains theorem target",
        "positive safe interval or forbidden K33/H=7 conflict packet",
        ("state-lift owner labels", "K33 incidence address"),
        "prove K33 family compression or construct the HYP-2908 state lift",
    ),
    PacketHandoff(
        "unit petal and P10+GW strip",
        "labelled_source_packet",
        "exact_interval_fejer",
        "selected P10+GW certified negative; C27 petal route available",
        "positive safe interval in the unit-petal packet",
        ("C27 unit-hole label", "two-block splice order"),
        "prove the petal/two-block family template, not just named rows",
    ),
    PacketHandoff(
        "covering boundary-moment packets",
        "endpoint_bridge_graph",
        "exact_interval_fejer",
        "selected covering rows certified; full family needs multi-chart positivity",
        "strict safe component or positive boundary-moment image",
        ("apex lift section", "moment chart", "endpoint transition packet"),
        "build multi-chart gK8/L_y or Fejer family templates for all covering packets",
    ),
    PacketHandoff(
        "lcm-tail / finite-denominator wall",
        "ramanujan_exact_period",
        "analytic_sieve_kaczynski",
        "q<=42 twist evidence closes audited lcm tails; global fixed basis impossible",
        "adaptive exact-period witness or labelled Kaczynski resonance",
        ("prime-power denominator labels", "approach class"),
        "prove dynamic ladder or admissible smoothed exponential-sum lemma",
    ),
    PacketHandoff(
        "SOURCE-SPECTRUM-UNKNOWN / F7",
        "labelled_source_packet",
        "tournament_state_lift",
        "only honest falsifier bucket",
        "simultaneous invisibility to all known primal/dual certificates",
        ("which certificate failed", "residual sector label"),
        "define F7 and prove it is empty or constructs TournamentStateLift",
    ),
]

FIXED_MARGIN_SWAP_IMPORT = [
    (
        "fixed margins are labels",
        "binary row/column margins behave like packet fibers: a comparison "
        "argument may move between chains only while preserving the margin "
        "fiber or recording its kernel.",
    ),
    (
        "reduce to a low-row core",
        "the swap-chain proof reduces arbitrary margins to three rows; the "
        "LRC14 analogue is to reduce arbitrary primitive rows to a finite "
        "packet core before applying Fejer, twist, or state-lift certificates.",
    ),
    (
        "split count and harmonic sectors",
        "the three-row analysis decomposes by column count and Johnson "
        "harmonic sectors; the LRC14 F7 bucket should be defined as a named "
        "harmonic residual sector, not as an anonymous failure of existing "
        "certificates.",
    ),
]


def carrier_index() -> dict[str, int]:
    return {carrier.name: i for i, carrier in enumerate(CARRIERS)}


def beats(a: Carrier, b: Carrier) -> bool:
    aw = sum(1 for x, y in zip(a.vector, b.vector) if x > y)
    bw = sum(1 for x, y in zip(a.vector, b.vector) if x < y)
    if aw != bw:
        return aw > bw
    return carrier_index()[a.name] < carrier_index()[b.name]


def tournament_edges() -> dict[str, set[str]]:
    edges = {c.name: set() for c in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        if beats(a, b):
            edges[a.name].add(b.name)
        else:
            edges[b.name].add(a.name)
    return edges


def score_hist(edges: dict[str, set[str]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for outs in edges.values():
        hist[len(outs)] = hist.get(len(outs), 0) + 1
    return dict(sorted(hist.items()))


def directed_3cycles(edges: dict[str, set[str]]) -> int:
    count = 0
    names = [c.name for c in CARRIERS]
    for a, b, c in combinations(names, 3):
        ab = b in edges[a]
        bc = c in edges[b]
        ca = a in edges[c]
        ba = a in edges[b]
        cb = b in edges[c]
        ac = c in edges[a]
        if (ab and bc and ca) or (ba and cb and ac):
            count += 1
    return count


def sccs(edges: dict[str, set[str]]) -> list[list[str]]:
    names = [c.name for c in CARRIERS]
    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    low: dict[str, int] = {}
    out: list[list[str]] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in edges[v]:
            if w not in indices:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp: list[str] = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            out.append(comp)

    for name in names:
        if name not in indices:
            strongconnect(name)
    return sorted((sorted(comp, key=carrier_index().get) for comp in out), key=lambda c: carrier_index()[c[0]])


def hamiltonian_path_count(edges: dict[str, set[str]]) -> tuple[int, list[str]]:
    n = len(CARRIERS)
    idx = carrier_index()
    names = [c.name for c in CARRIERS]

    @lru_cache(maxsize=None)
    def dp(last: int, mask: int) -> int:
        if mask == (1 << n) - 1:
            return 1
        total = 0
        last_name = names[last]
        for nxt, nxt_name in enumerate(names):
            if mask & (1 << nxt):
                continue
            if nxt_name in edges[last_name]:
                total += dp(nxt, mask | (1 << nxt))
        return total

    total = sum(dp(i, 1 << i) for i in range(n))
    path = sorted(names, key=idx.get)
    return total, path


def open_arrow_summary() -> list[str]:
    return [
        "O1 source-kernel exclusion: every primitive row enters this atlas.",
        "O2 formal interval backend: replace prototype pi/trig intervals by checkable balls.",
        "O3 family compression: selected Fejer/twist certificates lift to packet templates.",
        "O4 admissible smoothing: analytic-sieve/Kaczynski quotients retain approach labels.",
        "O5 state-lift construction: any zero-open non-AP/GW residual builds HYP-2908/THM-572.",
        "O6 F7 definition: if a Johnson-harmonic residual sector exists, name its preserved predicate.",
    ]


def main() -> None:
    edges = tournament_edges()
    hp_count, hp_path = hamiltonian_path_count(edges)
    comps = sccs(edges)

    print("S164 LRC14 CERTIFICATE HANDOFF ATLAS")
    print("=" * 72)
    print("Assumption challenge: vertices are proof carriers / obligations, not runners.")
    print("Pairwise observable: retained LRC predicate payload across handoffs.")
    print()
    print("Carrier tournament")
    print("------------------")
    print(f"score_hist: {score_hist(edges)}")
    print(f"directed_3cycles: {directed_3cycles(edges)}")
    print(f"SCC_sizes: {[len(c) for c in comps]}")
    print(f"Hamiltonian_path_count: {hp_count}")
    print("tie/retention Hamiltonian path:")
    print("  " + " > ".join(hp_path))
    print()

    print("Carrier roles")
    print("-------------")
    for carrier in CARRIERS:
        print(f"- {carrier.name}: vector={carrier.vector}")
        print(f"  role={carrier.role}")
        print(f"  must_retain={', '.join(carrier.must_retain)}")
        if carrier.may_forget:
            print(f"  may_forget={', '.join(carrier.may_forget)}")
    print()

    print("Packet handoff table")
    print("--------------------")
    for handoff in HANDOFFS:
        print(f"- packet={handoff.packet}")
        print(f"  {handoff.current_carrier} -> {handoff.exit_carrier}")
        print(f"  status={handoff.status}")
        print(f"  retained_predicate={handoff.retained_predicate}")
        print(f"  destroyed_without_guard={', '.join(handoff.destroyed_without_guard)}")
        print(f"  open_arrow={handoff.open_arrow}")
    print()

    print("Open arrows for a proof")
    print("-----------------------")
    for arrow in open_arrow_summary():
        print(f"- {arrow}")
    print()

    print("Fixed-margin swap-chain import")
    print("------------------------------")
    print("Source pattern: arXiv:2606.22636 reduces fixed-margin swap chains to a")
    print("low-row heat-bath core and then splits count sectors from Johnson harmonics.")
    for name, transfer in FIXED_MARGIN_SWAP_IMPORT:
        print(f"- {name}: {transfer}")
    print()

    print("Zipper theorem candidate")
    print("------------------------")
    print(
        "If O1-O6 hold, then every primitive LRC14 row either has a strict "
        "witness/dual certificate in the atlas, is the AP/GW equality atom, "
        "or constructs the forbidden HYP-2908/THM-572 state lift.  This would "
        "close LRC14.  The atlas does not prove O1-O6; it localizes them."
    )


if __name__ == "__main__":
    main()
