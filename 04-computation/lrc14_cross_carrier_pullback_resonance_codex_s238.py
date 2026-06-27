#!/usr/bin/env python3
"""S238: cross-carrier pullback resonance scout for remaining LRC14 angles.

This is a proof-interface computation, not a proof of LRC14.  It takes the
CPI carrier-pullback index seriously as an object of computation: proof
carriers, hidden coordinates, sidecars, analytic clocks, automata, topology,
and formal exits are the vertices.  Runners and routes remain possible
secondary coordinates, but they are not assumed to be the universal vertex set.

The script asks two concrete questions:

1. Which cross-disciplinary portfolios cover the remaining proof obligations?
2. Which single carriers are dangerous because they preserve a useful signal
   while destroying the coordinate that the obligation actually needs?

Tournament Analysis vertices are proof carriers / CPI pullback rows, not
runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations


Axis = str


CORE_AXES: tuple[Axis, ...] = (
    "status",
    "route",
    "exact_scale",
    "topology",
    "owner",
    "period_deck",
    "analytic_certificate",
    "automaton_partial_cube",
    "crt_padic",
    "observer_cut",
    "hodge_cycle",
    "formal_exit",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    source: str
    discipline: str
    axes: frozenset[Axis]
    destroys: frozenset[Axis]
    cost: int
    guardrail: str


@dataclass(frozen=True)
class Obligation:
    name: str
    required: dict[Axis, int]
    reading: str


def fs(*items: Axis) -> frozenset[Axis]:
    return frozenset(items)


CARRIERS: tuple[Carrier, ...] = (
    Carrier(
        "labelled_packet_sheaf",
        "HYP-2963/CPI-020",
        "packet_sheaf",
        fs("status", "route", "packet_key", "family_label", "qdiv", "exact_scale", "destroyed_coordinate"),
        fs("topology", "owner", "analytic_certificate"),
        2,
        "unlabelled packets invite false analogies, but packet labels alone do not close a certificate",
    ),
    Carrier(
        "closed_arc_cech_topology",
        "HYP-3025/CPI-081",
        "topology",
        fs("status", "topology", "arc_cech_beta", "open_tope", "quotient_defect", "boundary_cycle"),
        fs("route", "period_deck", "exact_scale"),
        2,
        "Cech topology separates boundary/open shape but forgets arithmetic period unless sidecars return",
    ),
    Carrier(
        "safe_stalk_barcode_normal_fan",
        "HYP-3029/HYP-3015/HYP-3018/CPI-076,CPI-078",
        "stalk_geometry",
        fs(
            "status",
            "topology",
            "owner",
            "barcode_stalk",
            "normal_fan",
            "safe_component",
            "largest_stalk",
        ),
        fs("route", "period_deck"),
        2,
        "barcode shape is not a theorem route unless active support and owners remain attached",
    ),
    Carrier(
        "endpoint_owner_strip_current",
        "HYP-3042/HYP-3045/CPI-034",
        "owner_current",
        fs("status", "owner", "endpoint_owner_strip", "boundary_current", "owner_transfer", "route"),
        fs("period_deck", "exact_scale"),
        2,
        "owner strips can forget primitive clocks such as the AP-tail mod-13 split",
    ),
    Carrier(
        "ramanujan_primitive_period_deck",
        "HYP-3036/CPI-043",
        "periodic_arithmetic",
        fs("route", "period_deck", "prime_power_label", "primitive_safe_deck", "exact_period", "owner"),
        fs("topology"),
        3,
        "period projectors must be endpoint-aware; squarefree divisor shadows merge live prime powers",
    ),
    Carrier(
        "mobius_squarefree_blindness_report",
        "HYP-3032/CPI-053",
        "analytic_sieve",
        fs("analytic_capacity", "squarefree_blindness", "squarefree_status", "phi_density"),
        fs("prime_power_label", "owner", "route", "period_deck"),
        1,
        "mu^2/phi is a capacity meter, not a final equality atom; it zeros repeated-prime routes",
    ),
    Carrier(
        "haar_zeta_square",
        "HYP-3038/HYP-2991/CPI-036,CPI-060",
        "haar_harmonic",
        fs("status", "owner", "haar_zeta", "fixed_margin", "product_tiling", "analytic_certificate"),
        fs("route", "period_deck"),
        2,
        "row/column margins collide unless the mixed zeta coordinate is retained",
    ),
    Carrier(
        "fejer_interval_backend",
        "HYP-2981/CPI-041",
        "analytic_certificate",
        fs("status", "analytic_certificate", "fejer_interval", "packet_key", "sign_bound", "exact_scale"),
        fs("route", "owner"),
        3,
        "floating Fejer values are invalid without packet key, center, degree, and interval backend",
    ),
    Carrier(
        "crt_padic_residual_tree",
        "HYP-3033/CPI-056",
        "arithmetic_crt",
        fs("crt_padic", "unit_root", "zero_root_debt", "crt_recombination", "period_deck", "exact_scale"),
        fs("topology", "owner"),
        3,
        "denominator-only p-adic trees lose endpoint owners and route family",
    ),
    Carrier(
        "automatic_partial_cube_sidecar",
        "HYP-3063/HYP-3023/CPI-054,CPI-067,CPI-077",
        "automata_geometry",
        fs(
            "automaton_partial_cube",
            "theta_class",
            "gated_subcube",
            "simplex_face",
            "native_transition",
            "bit_phase",
        ),
        fs("status", "exact_scale", "owner", "route"),
        2,
        "Moser/fibbinary membership is telemetry; it mixes status and route without exact packet sidecars",
    ),
    Carrier(
        "observer_extension_cut_payload",
        "HYP-3054/HYP-3056/CPI-071",
        "observer_cut",
        fs("observer_cut", "source_mark", "deletion_fiber", "cut_signature", "root_object", "destroyed_coordinate"),
        fs("status", "route", "exact_scale"),
        2,
        "observer cuts name the lost payload but do not certify boundary/open status by themselves",
    ),
    Carrier(
        "rectangle_hourglass_diagonal_flow",
        "HYP-3065/HYP-3053/S217",
        "diagonal_flow",
        fs("observer_cut", "rectangle_residue", "hourglass_residue", "line_potential", "edge_sector", "fixed_path_flow"),
        fs("owner", "route", "exact_scale"),
        2,
        "rectangle residues are controlled-forgetting sidecars, not standalone LRC packets",
    ),
    Carrier(
        "hodge_cycle_lift",
        "HYP-3066/HYP-3071/CPI-033,CPI-061,CPI-065",
        "cohomology",
        fs("hodge_cycle", "formal_exit", "cycle_generators", "image_status", "cochain", "state_lift"),
        fs("exact_scale", "packet_key"),
        3,
        "rank/image language is legal only after actual packet cochains are emitted",
    ),
    Carrier(
        "median_route_center_control",
        "HYP-3069/HYP-3070",
        "median_sidecar",
        fs("route", "median_center", "sidecar_page", "root_object", "owner", "formal_exit", "destroyed_coordinate"),
        fs("exact_scale", "analytic_certificate"),
        3,
        "a raw route vocabulary center is not proof-compatible until legal sidecar pages are attached",
    ),
    Carrier(
        "toeplitz_square_scale_gate",
        "HYP-3064/CPI-042",
        "noncollapse_geometry",
        fs("exact_scale", "toeplitz_gate", "d4_orbit", "midpoint_balance", "psd_dual", "analytic_certificate"),
        fs("route", "owner"),
        3,
        "approximate four-witness shadows must keep positive scale or route to named debt",
    ),
    Carrier(
        "roth_minkowski_low_height_wall",
        "HYP-3062",
        "diophantine_geometry",
        fs("exact_scale", "low_height_wall", "relation_lattice", "covolume", "successive_minima", "crt_padic"),
        fs("route", "status"),
        3,
        "volume/height pressure is only a sidecar unless the deleted anti-cosets are named",
    ),
    Carrier(
        "k33_state_lift_incidence",
        "THM-572/CPI-025",
        "state_lift",
        fs("state_lift", "route", "k33_incidence", "forbidden_lift", "formal_exit", "hodge_cycle"),
        fs("topology", "exact_scale"),
        3,
        "positive K33 incidence without state address is unfinished residual debt",
    ),
    Carrier(
        "carrier_fusion_switchboard",
        "HYP-3026/CPI-082",
        "fusion",
        fs(
            "status",
            "route",
            "topology",
            "owner",
            "barcode_stalk",
            "exact_scale",
            "fusion_signature",
            "crt_padic",
            "haar_zeta",
        ),
        fs("formal_exit"),
        4,
        "a fused signature can still hide which subcarrier is load-bearing",
    ),
    Carrier(
        "sidechannel_repair_ladder",
        "HYP-3027/HYP-3028/HYP-3039",
        "repair_ladder",
        fs("status", "route", "destroyed_coordinate", "repair_stage", "first_repair_cochain", "formal_exit"),
        fs("exact_scale", "topology"),
        3,
        "a repair stage is only admissible when it names the first lost coordinate",
    ),
    Carrier(
        "exact_farey_kpq_scale",
        "HYP-2934/CPI-049,CPI-050,CPI-068",
        "farey_incidence",
        fs("exact_scale", "qdiv", "binding_scale", "farey_excess", "kpq_incidence", "route"),
        fs("topology", "owner"),
        2,
        "product/farey shadows are unsafe without endpoint and packet-family sidecars",
    ),
    Carrier(
        "pair_good_blocker_grammar",
        "HYP-3021/HYP-3022/CPI-080",
        "local_switch",
        fs("status", "pair_blocker", "blocker_tooth", "generator_class", "barcode_stalk", "normal_fan", "owner"),
        fs("route", "exact_scale"),
        2,
        "raw decoy count is misleading; generator and barcode/normal-fan sidecars must stay attached",
    ),
    Carrier(
        "hyperbolic_power_guard",
        "HYP-3058/CPI-072",
        "curvature_power",
        fs("automaton_partial_cube", "valuation_budget", "power_guard", "safe_component", "exact_scale", "geometry_regime"),
        fs("route", "owner", "topology"),
        2,
        "the reciprocal-sum <1 condition is hyperbolic curvature debt, not an LRC quotient",
    ),
)


OBLIGATIONS: tuple[Obligation, ...] = (
    Obligation(
        "residual_route_mixed_fibers",
        {
            "status": 3,
            "route": 3,
            "topology": 2,
            "barcode_stalk": 2,
            "owner": 2,
            "period_deck": 1,
            "destroyed_coordinate": 1,
        },
        "After coarse ET+unit, route-mixed fibers need topology/stalk first teeth plus owner and clock sidecars.",
    ),
    Obligation(
        "q23_squarefree_petal_covering",
        {
            "exact_scale": 3,
            "haar_zeta": 2,
            "owner": 2,
            "period_deck": 2,
            "squarefree_blindness": 2,
            "route": 1,
        },
        "The q=23 petal/covering pair needs zeta, owner strips, exact M, and a squarefree blindness warning.",
    ),
    Obligation(
        "automatic_partial_cube_route_purity",
        {
            "automaton_partial_cube": 3,
            "exact_scale": 3,
            "topology": 2,
            "owner": 2,
            "route": 2,
            "status": 1,
        },
        "Moser/fibbinary/partial-cube telemetry is useful only after exact scale, topology, owner, route, and status return.",
    ),
    Obligation(
        "two_tail_ap_clock_search",
        {
            "period_deck": 3,
            "exact_scale": 3,
            "owner": 2,
            "route": 2,
            "crt_padic": 1,
        },
        "Search two-tail AP-core residuals by primitive clock, fixed point, owner strip, exact M, and CRT sidecar.",
    ),
    Obligation(
        "k33_f7_state_lift_exit",
        {
            "state_lift": 3,
            "hodge_cycle": 3,
            "formal_exit": 2,
            "route": 2,
            "topology": 1,
            "exact_scale": 1,
        },
        "K33/F7 debt must become a named state-lift cochain with formal exit and retained packet geometry.",
    ),
    Obligation(
        "observer_extension_rootless_payload",
        {
            "observer_cut": 3,
            "destroyed_coordinate": 2,
            "rectangle_residue": 2,
            "root_object": 2,
            "formal_exit": 1,
            "exact_scale": 1,
        },
        "A000568/rootless extension defects need observer cuts, rectangle/hourglass residues, root objects, and exact scale.",
    ),
    Obligation(
        "carrier_pullback_formalization",
        {
            "formal_exit": 3,
            "destroyed_coordinate": 3,
            "status": 2,
            "route": 2,
            "median_center": 2,
        },
        "A quotient is admissible only after the preserved predicate, destroyed coordinate, and legal exit are explicit.",
    ),
    Obligation(
        "analytic_certificate_backend",
        {
            "analytic_certificate": 3,
            "fejer_interval": 2,
            "period_deck": 2,
            "haar_zeta": 2,
            "exact_scale": 2,
            "owner": 1,
        },
        "Fejer/Ramanujan/Haar backends need packet-keyed intervals, primitive periods, zeta, exact M, and owners.",
    ),
    Obligation(
        "multi_scale_binding_spectrum",
        {
            "exact_scale": 3,
            "binding_scale": 3,
            "route": 2,
            "qdiv": 1,
            "observer_cut": 1,
        },
        "The spectrum/binding-scale route needs exact Farey scale with route and observer-cut sidecars.",
    ),
)


TIE_PATH: tuple[str, ...] = (
    "carrier_fusion_switchboard",
    "labelled_packet_sheaf",
    "sidechannel_repair_ladder",
    "hodge_cycle_lift",
    "median_route_center_control",
    "safe_stalk_barcode_normal_fan",
    "closed_arc_cech_topology",
    "endpoint_owner_strip_current",
    "exact_farey_kpq_scale",
    "ramanujan_primitive_period_deck",
    "haar_zeta_square",
    "fejer_interval_backend",
    "crt_padic_residual_tree",
    "observer_extension_cut_payload",
    "rectangle_hourglass_diagonal_flow",
    "automatic_partial_cube_sidecar",
    "hyperbolic_power_guard",
    "pair_good_blocker_grammar",
    "toeplitz_square_scale_gate",
    "roth_minkowski_low_height_wall",
    "k33_state_lift_incidence",
    "mobius_squarefree_blindness_report",
)


def core_word(carrier: Carrier) -> str:
    return "".join("1" if axis in carrier.axes else "0" for axis in CORE_AXES)


def obligation_score(axes: frozenset[Axis], obligation: Obligation) -> tuple[int, int, tuple[Axis, ...]]:
    got = sum(weight for axis, weight in obligation.required.items() if axis in axes)
    total = sum(obligation.required.values())
    missing = tuple(axis for axis in obligation.required if axis not in axes)
    return got, total, missing


def portfolio_axes(carriers: tuple[Carrier, ...]) -> frozenset[Axis]:
    axes: set[Axis] = set()
    for carrier in carriers:
        axes.update(carrier.axes)
    return frozenset(axes)


def discipline_count(carriers: tuple[Carrier, ...]) -> int:
    return len({carrier.discipline for carrier in carriers})


def portfolio_cost(carriers: tuple[Carrier, ...]) -> int:
    return sum(carrier.cost for carrier in carriers)


def portfolio_covers(carriers: tuple[Carrier, ...], obligation: Obligation) -> bool:
    axes = portfolio_axes(carriers)
    return all(axis in axes for axis in obligation.required)


def best_portfolios_for_obligation(obligation: Obligation, max_size: int = 4) -> list[tuple[Carrier, ...]]:
    found: list[tuple[Carrier, ...]] = []
    for size in range(1, max_size + 1):
        for combo in combinations(CARRIERS, size):
            if portfolio_covers(combo, obligation):
                found.append(combo)
        if found:
            break
    return sorted(found, key=lambda combo: (-discipline_count(combo), portfolio_cost(combo), [c.name for c in combo]))[:5]


def global_requirement_axes() -> frozenset[Axis]:
    axes: set[Axis] = set()
    for obligation in OBLIGATIONS:
        axes.update(obligation.required)
    return frozenset(axes)


def global_portfolios(max_size: int = 9) -> list[tuple[Carrier, ...]]:
    target = global_requirement_axes()
    found: list[tuple[Carrier, ...]] = []
    for size in range(1, max_size + 1):
        for combo in combinations(CARRIERS, size):
            if target <= portfolio_axes(combo):
                found.append(combo)
        if found:
            break
    return sorted(found, key=lambda combo: (-discipline_count(combo), portfolio_cost(combo), [c.name for c in combo]))[:10]


def carrier_payload(carrier: Carrier) -> tuple[int, int, int, int, int, int]:
    full_obligations = 0
    weighted = 0
    critical = 0
    for obligation in OBLIGATIONS:
        got, total, _ = obligation_score(carrier.axes, obligation)
        weighted += got
        if got == total:
            full_obligations += 1
        critical += sum(1 for axis in obligation.required if axis in carrier.axes)
    return (
        full_obligations,
        weighted,
        critical,
        len(carrier.axes),
        -len(carrier.destroys),
        -carrier.cost,
    )


def orient(a: Carrier, b: Carrier) -> tuple[str, str]:
    pa = carrier_payload(a)
    pb = carrier_payload(b)
    if pa > pb:
        return a.name, b.name
    if pb > pa:
        return b.name, a.name
    rank = {name: i for i, name in enumerate(TIE_PATH)}
    return (a.name, b.name) if rank[a.name] < rank[b.name] else (b.name, a.name)


def tournament_edges() -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    for a, b in combinations(CARRIERS, 2):
        edges.add(orient(a, b))
    return edges


def score_hist(edges: set[tuple[str, str]]) -> dict[int, int]:
    scores = Counter()
    for winner, _ in edges:
        scores[winner] += 1
    return dict(sorted(Counter(scores[carrier.name] for carrier in CARRIERS).items()))


def directed_3cycles(edges: set[tuple[str, str]]) -> int:
    edge_set = set(edges)
    count = 0
    for a, b, c in combinations([carrier.name for carrier in CARRIERS], 3):
        if ((a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set) or (
            (a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set
        ):
            count += 1
    return count


def scc_sizes(edges: set[tuple[str, str]]) -> list[int]:
    names = [carrier.name for carrier in CARRIERS]
    graph = {name: set() for name in names}
    rev = {name: set() for name in names}
    for a, b in edges:
        graph[a].add(b)
        rev[b].add(a)
    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)
    seen.clear()
    sizes: list[int] = []
    for name in reversed(order):
        if name in seen:
            continue
        stack = [name]
        seen.add(name)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for w in rev[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edges: set[tuple[str, str]]) -> int:
    names = [carrier.name for carrier in CARRIERS]
    n = len(names)
    index = {name: i for i, name in enumerate(names)}
    outmask = [0] * n
    for a, b in edges:
        outmask[index[a]] |= 1 << index[b]
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp.get((mask, last), 0)
            if not value:
                continue
            candidates = outmask[last] & ~mask
            while candidates:
                bit = candidates & -candidates
                nxt = bit.bit_length() - 1
                dp[(mask | bit, nxt)] = dp.get((mask | bit, nxt), 0) + value
                candidates ^= bit
    full = (1 << n) - 1
    return sum(dp.get((full, i), 0) for i in range(n))


def topo_path(edges: set[tuple[str, str]]) -> list[str]:
    names = [carrier.name for carrier in CARRIERS]
    indeg = {name: 0 for name in names}
    graph = {name: set() for name in names}
    for a, b in edges:
        graph[a].add(b)
        indeg[b] += 1
    queue = deque(sorted([name for name in names if indeg[name] == 0], key=TIE_PATH.index))
    order: list[str] = []
    while queue:
        v = queue.popleft()
        order.append(v)
        for w in sorted(graph[v], key=TIE_PATH.index):
            indeg[w] -= 1
            if indeg[w] == 0:
                queue.append(w)
    return order


def resonance_score(combo: tuple[Carrier, ...]) -> tuple[int, int, int, int, tuple[str, ...]]:
    axes = portfolio_axes(combo)
    individual_best = max(len(carrier.axes) for carrier in combo)
    shared = sum(1 for axis in axes if sum(axis in carrier.axes for carrier in combo) > 1)
    return (
        len(axes) - individual_best,
        discipline_count(combo),
        -shared,
        -portfolio_cost(combo),
        tuple(carrier.name for carrier in combo),
    )


def top_resonances(size: int, limit: int = 8) -> list[tuple[Carrier, ...]]:
    combos = list(combinations(CARRIERS, size))
    return sorted(combos, key=resonance_score, reverse=True)[:limit]


def fmt_combo(combo: tuple[Carrier, ...]) -> str:
    return " + ".join(carrier.name for carrier in combo)


def fmt_missing(missing: tuple[Axis, ...]) -> str:
    return "none" if not missing else ",".join(missing)


def main() -> None:
    print("=== LRC14 cross-carrier pullback resonance scout S238 ===")
    print("source=CPI carrier-pullback index + HYP-2963 residual-stack summaries")
    print(f"carriers={len(CARRIERS)} obligations={len(OBLIGATIONS)} core_axes={len(CORE_AXES)}")
    print()

    print("[0] Assumption challenge")
    print("  considered vertex sets:")
    print("    runners, gaps, route labels, proof obligations, CPI carrier rows, hidden coordinates,")
    print("    endpoint-owner strips, primitive decks, analytic clocks, automaton states,")
    print("    partial-cube Theta classes, CRT roots, observer cuts, rectangle/hourglass residues,")
    print("    Hodge cochains, median-center pages, state-lift sectors, and formal exits.")
    print("  chosen tournament vertices:")
    print("    proof-carrier pullbacks / CPI rows, not runners.")
    print("  preserved LRC predicates:")
    print("    boundary/open status, theorem route, exact scale, endpoint/topology, primitive period,")
    print("    analytic certificate, automaton telemetry, observer-cut payload, cycle image status,")
    print("    and formal discharge/debt labels, depending on the carrier.")
    print("  destroyed information:")
    print("    recorded carrier-by-carrier; a quotient is not legal unless the destroyed coordinate is")
    print("    reconstructed, dual-annihilated, descended, AP/GW boundary, or named THM-572/F7 debt.")
    print()

    print("[1] Duodecimal incident words over the core payload alphabet")
    print("  alphabet=" + ",".join(CORE_AXES))
    for carrier in CARRIERS:
        print(
            f"  {carrier.name:34s} word={core_word(carrier)} "
            f"discipline={carrier.discipline:22s} source={carrier.source}"
        )
    print("  reading: the useful object is the incident word plus the lost-coordinate ledger,")
    print("  not the raw carrier name.  Several attractive carriers are sparse in the core word.")
    print()

    print("[2] Single-carrier blindness report")
    for carrier in CARRIERS:
        dangerous = []
        for obligation in OBLIGATIONS:
            _, total, missing = obligation_score(carrier.axes, obligation)
            destroyed_needed = sorted(set(missing) & set(carrier.destroys))
            if destroyed_needed:
                dangerous.append(f"{obligation.name}:{','.join(destroyed_needed)}")
        print(
            f"  {carrier.name:34s} payload={carrier_payload(carrier)} "
            f"destroys={sorted(carrier.destroys)}"
        )
        print(f"    guardrail={carrier.guardrail}")
        if dangerous:
            print(f"    destroys_needed_by={'; '.join(dangerous)}")
        else:
            print("    destroys_needed_by=none in this finite obligation list")
    print()

    print("[3] Obligation coverage and minimal cross-carrier portfolios")
    for obligation in OBLIGATIONS:
        print(f"  obligation={obligation.name}")
        print(f"    reading={obligation.reading}")
        top_singles = sorted(
            (
                (
                    obligation_score(carrier.axes, obligation)[0],
                    -carrier.cost,
                    carrier.name,
                    obligation_score(carrier.axes, obligation)[2],
                )
                for carrier in CARRIERS
            ),
            reverse=True,
        )[:5]
        for got, neg_cost, name, missing in top_singles:
            total = sum(obligation.required.values())
            print(f"    single {name:34s} score={got}/{total} missing={fmt_missing(missing)}")
        portfolios = best_portfolios_for_obligation(obligation)
        if portfolios:
            for combo in portfolios[:3]:
                print(
                    f"    cover size={len(combo)} disciplines={discipline_count(combo)} "
                    f"cost={portfolio_cost(combo)} :: {fmt_combo(combo)}"
                )
        else:
            print("    cover none with size<=4")
    print()

    print("[4] Global set-cover portfolios across all obligation axes")
    print(f"  target_axes={sorted(global_requirement_axes())}")
    global_covers = global_portfolios()
    if not global_covers:
        print("  no cover found with size<=9")
    for combo in global_covers:
        print(
            f"  cover size={len(combo)} disciplines={discipline_count(combo)} "
            f"cost={portfolio_cost(combo)} :: {fmt_combo(combo)}"
        )
    print("  reading: the first full covers appear only at size 9 in this finite atlas,")
    print("  so no small universal scalar-like bundle is visible.  The recurring skeleton is")
    print("  packet/fusion + owner/topology/stalk + period/analytic + observer/formal/state-lift.")
    print()

    print("[5] Cross-disciplinary resonance pairs")
    for combo in top_resonances(2):
        axes = sorted(portfolio_axes(combo))
        print(
            f"  score={resonance_score(combo)[:4]} disciplines={discipline_count(combo)} "
            f":: {fmt_combo(combo)}"
        )
        print(f"    union_axes={axes}")
    print()

    print("[6] Cross-disciplinary resonance triples")
    for combo in top_resonances(3):
        axes = sorted(portfolio_axes(combo))
        print(
            f"  score={resonance_score(combo)[:4]} disciplines={discipline_count(combo)} "
            f":: {fmt_combo(combo)}"
        )
        print(f"    union_axes={axes}")
    print()

    edges = tournament_edges()
    print("[7] Tournament Analysis")
    print("  vertices_are=proof-carrier pullbacks / CPI rows, not runners")
    print("  pairwise_observable=(full_obligation_count, weighted_axis_coverage, critical_axis_hits, payload_width, -destroyed_count, -cost)")
    print("  gauge=orient toward more retained proof payload; ties use the CPI pullback tie path")
    print(f"  score_hist={score_hist(edges)}")
    print(f"  directed_3cycles={directed_3cycles(edges)}")
    print(f"  scc_sizes={scc_sizes(edges)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(edges)}")
    print("  tie_hamiltonian_path=" + " > ".join(topo_path(edges)))
    print()

    print("[8] Two creative proof angles extracted")
    print("  angle A: blindness-pair proof")
    print("    For every scalar or low-payload carrier, pair it with the carrier that restores")
    print("    exactly the coordinate it destroys.  Examples: mu^2/phi pairs with Ramanujan/CRT")
    print("    for prime powers; automata pair with exact Farey scale plus arc-Cech/owner;")
    print("    observer cuts pair with formal exits plus exact scale.  The proof target is a")
    print("    no-unpaid-forgetting lemma.")
    print("  angle B: resonance-portfolio proof")
    print("    Work with minimal portfolios per obligation instead of one master quotient.")
    print("    The covering family, q=23 square, automatic fibers, and F7/K33 debt each need")
    print("    a different three- or four-carrier packet.  Glue them by the formal exit ledger,")
    print("    not by forcing all rows through one scalar invariant.")
    print()

    print("[9] Next exact job")
    print("  Emit actual HYP-2963 packet rows with these fields:")
    print("    carrier_pullback_row_id, core_incident_word, preserved_lrc_predicate,")
    print("    destroyed_coordinate, required_sidecar, blindness_pair_id, resonance_portfolio_id,")
    print("    status_mixing_result, route_mixing_result, and legal_exit_status.")
    print("  Then test whether the listed portfolios make the residual coarse fibers")
    print("  status-pure and route-pure before any new theorem debt is named.")


if __name__ == "__main__":
    main()
