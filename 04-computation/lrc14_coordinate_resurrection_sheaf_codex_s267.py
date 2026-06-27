#!/usr/bin/env python3
"""HYP-3118 coordinate-resurrection scout for LRC14 proof carriers.

This is an abstract executable synthesis, not an LRC proof.

The point is to test the post-HYP-3116 question:

    once a quotient has a nonempty missing-input vector, which minimal
    sidecar, adjoint section, or sheaf-gluing datum resurrects the coordinate
    needed by the next proof obligation?

Tournament Analysis uses reframes as vertices, not runners, roots, gates,
constants, or scalar signals.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations


CORE = frozenset({"finite_address", "observer_gluing", "endpoint_owner", "uniformity"})
LIVE = frozenset(
    {
        "root_ear_payload",
        "relation_lattice_shape",
        "component_bound",
        "cocycle_exactness",
        "state_lift",
        "pde_weak_form",
    }
)


@dataclass(frozen=True)
class Sidecar:
    name: str
    source: str
    coords: frozenset[str]
    note: str
    cost: int = 1


@dataclass(frozen=True)
class Shortcut:
    name: str
    source: str
    kept: frozenset[str]
    destroyed: frozenset[str]
    challenged_assumption: str


@dataclass(frozen=True)
class Reframe:
    name: str
    axes: dict[str, int]
    preserves: str
    destroys: str
    next_signal: str


@dataclass(frozen=True)
class Inspiration:
    name: str
    external_pattern: str
    lrc_reframe: str
    resurrects: tuple[str, ...]
    new_signal: str


SIDECARS = [
    Sidecar(
        "labelled_packet_decision_tree",
        "HYP-3117/HYP-2963",
        frozenset({"finite_address", "endpoint_owner", "packet_label", "route_purity"}),
        "old labelled packets resurrect finite address and endpoint owner together",
    ),
    Sidecar(
        "observer_gluing_certificate",
        "HYP-3098/HYP-3107",
        frozenset({"observer_gluing", "chart_overlap", "terminal_exit"}),
        "observer chart compatibility and terminal discharge interface",
    ),
    Sidecar(
        "ear_decomposition_gluing_grammar",
        "ear-decomposition theorems/HYP-3108",
        frozenset(
            {
                "observer_gluing",
                "ear_decomposition",
                "ear_parity",
                "nested_ear_order",
                "strong_connectivity_witness",
            }
        ),
        "treats directed/odd/nested ears as named gluing grammars for proof-state charts",
    ),
    Sidecar(
        "proof_circuit_uniformity_bus",
        "HYP-3116/HYP-3117",
        frozenset({"uniformity", "proof_circuit_basis", "missing_input_vector"}),
        "turns low-depth shortcuts into uniform proof-state circuits",
    ),
    Sidecar(
        "savitch_midpoint_certificate_ladder",
        "Savitch/HYP-3108/HYP-3117",
        frozenset({"uniformity", "proof_circuit_basis", "savitch_midpoint", "reachability_halving"}),
        "uses midpoint reachability recursion as a uniform certificate ladder, not a stored path",
    ),
    Sidecar(
        "lee_yang_ear_payload_gate",
        "HYP-3112/HYP-3109",
        frozenset({"root_ear_payload", "root_zero_locus", "ear_payload_A_vector", "root_motion"}),
        "restores root motion and one-runner ear payload, not only root count",
    ),
    Sidecar(
        "pgf_root_curve_section",
        "HYP-3108/HYP-3109/HYP-3112",
        frozenset(
            {
                "root_ear_payload",
                "root_zero_locus",
                "pgf_curve",
                "root_argument_path",
                "lee_yang_extremality",
                "discriminant_wall",
            }
        ),
        "promotes the whole PGF curve and zero trajectory over any single p0 scalar",
    ),
    Sidecar(
        "root_lattice_two_map_portfolio",
        "HYP-3113",
        frozenset({"root_ear_payload", "relation_lattice_shape", "savitch_midpoint", "ear_certificate"}),
        "couples analytic root portfolio to certificate ladder",
    ),
    Sidecar(
        "bravais_relation_shape_atlas",
        "Bravais/HYP-3115",
        frozenset({"relation_lattice_shape", "bravais_shape", "basis_reduction_word", "low_height_wall"}),
        "records lattice shape class and low-height relation walls before using volume pressure",
    ),
    Sidecar(
        "minkowski_relation_wall_gate",
        "HYP-3115/HYP-3062",
        frozenset({"relation_lattice_shape", "low_height_wall", "bravais_shape"}),
        "relation-rich rows need wall labels and shape, not raw covolume",
    ),
    Sidecar(
        "normal_fan_component_packet",
        "HYP-3101/HYP-3096",
        frozenset({"component_bound", "normal_fan_chamber", "cech_barcode"}),
        "component-bound proof route for normalized lonely-set packets",
    ),
    Sidecar(
        "first_obstruction_cocycle_gate",
        "HYP-3102",
        frozenset({"cocycle_exactness", "obstruction_syndrome", "dual_annihilator"}),
        "hidden payload difference is exact, generated, killed, descended, or debt",
    ),
    Sidecar(
        "state_lift_K33_gate",
        "HYP-3094/THM-572",
        frozenset({"state_lift", "active_binder", "K33_cross_handoff"}),
        "routes cross-handoff packets through named state-lift debt or construction",
    ),
    Sidecar(
        "pde_weak_form_operator",
        "HYP-3111/S264b",
        frozenset({"pde_weak_form", "boundary_data", "stiffness_matrix", "zero_mode"}),
        "weak-form carrier preserves boundary/operator data before energy analogies",
    ),
    Sidecar(
        "phi4_quartic_density_stress_test",
        "phi4_2/HYP-3111",
        frozenset({"pde_weak_form", "boundary_data", "zero_mode", "quartic_density_moment_curve"}),
        "quartic density is a stress test for moment/root analogies only with boundary data retained",
    ),
    Sidecar(
        "finite_interval_margin_packet",
        "HYP-3114",
        frozenset({"approximation_margin", "finite_grid_hit", "endpoint_margin"}),
        "resurrects rational-grid witness transfer only after an interval exists",
    ),
    Sidecar(
        "route_center_median_ledger",
        "HYP-3069/HYP-3070",
        frozenset({"route_center", "median_center", "sidecar_closure"}),
        "legal route centers require sidecar closure before median arguments",
    ),
]


SHORTCUTS = [
    Shortcut(
        "raw_p0_scalar",
        "HYP-3103/HYP-3116",
        frozenset({"scalar_value"}),
        frozenset(CORE | LIVE),
        "scalar value preserves neither mechanism nor proof exit",
    ),
    Shortcut(
        "one_literal_apex7_threshold",
        "HYP-3115/HYP-3116",
        frozenset({"finite_bank_threshold", "apex7_angle"}),
        frozenset({"finite_address", "observer_gluing", "endpoint_owner", "uniformity"}),
        "finite classifier may be a uniformity warning, not a proof",
    ),
    Shortcut(
        "root_count_only",
        "HYP-3109/HYP-3112",
        frozenset({"root_count"}),
        frozenset({"root_zero_locus", "ear_payload_A_vector", "finite_address", "observer_gluing"}),
        "real-root count forgets root motion and the extension payload",
    ),
    Shortcut(
        "ear_payload_without_owner",
        "HYP-3112",
        frozenset({"root_ear_payload", "ear_payload_A_vector"}),
        frozenset({"endpoint_owner", "finite_address", "observer_gluing", "uniformity"}),
        "ear payload is not terminal until owner and proof route survive",
    ),
    Shortcut(
        "raw_minkowski_volume",
        "HYP-3111/HYP-3115",
        frozenset({"relation_lattice_shape", "covolume"}),
        frozenset({"finite_address", "observer_gluing", "endpoint_owner", "uniformity", "low_height_wall"}),
        "volume pressure needs a declared body and route predicate",
    ),
    Shortcut(
        "raw_ising_energy",
        "HYP-3111/HYP-3115",
        frozenset({"root_ear_payload", "ising_energy"}),
        frozenset({"finite_address", "observer_gluing", "endpoint_owner", "uniformity", "boundary_data"}),
        "spin energy forgets which domain wall is a legal root-collision exit",
    ),
    Shortcut(
        "direct_time_interval_only",
        "HYP-3114/HYP-3088",
        frozenset({"approximation_margin", "finite_grid_hit"}),
        frozenset({"finite_address", "observer_gluing", "endpoint_owner", "uniformity", "normalized_clock"}),
        "raw time intervals do not preserve the normalized slow/ruler clock",
    ),
    Shortcut(
        "first_obstruction_only",
        "HYP-3102",
        frozenset({"cocycle_exactness", "obstruction_syndrome"}),
        frozenset({"finite_address", "observer_gluing", "endpoint_owner", "uniformity"}),
        "knowing the first obstruction is not the same as discharging it",
    ),
    Shortcut(
        "component_bound_only",
        "HYP-3101",
        frozenset({"component_bound", "cech_barcode"}),
        frozenset({"finite_address", "observer_gluing", "endpoint_owner", "uniformity"}),
        "component control must feed a terminal observer certificate",
    ),
    Shortcut(
        "pde_energy_only",
        "HYP-3111/S264b",
        frozenset({"pde_weak_form", "stiffness_matrix"}),
        frozenset({"finite_address", "observer_gluing", "endpoint_owner", "uniformity", "boundary_data"}),
        "energy analogies collapse without boundary and zero-mode sidecars",
    ),
    Shortcut(
        "proof_circuit_compiler_without_payload",
        "HYP-3117",
        frozenset({"proof_circuit_basis", "missing_input_vector", "uniformity"}),
        frozenset({"finite_address", "observer_gluing", "endpoint_owner"}),
        "compiler detects missing coordinates before it supplies them",
    ),
]


INSPIRATIONS = [
    Inspiration(
        "Savitch midpoint recursion",
        "reachability can be certified by recursively naming midpoints while reusing space",
        "proof reachability should store a recomputable midpoint certificate, not the whole route",
        ("uniformity", "savitch_midpoint", "proof_circuit_basis"),
        "midpoint_certificate_depth_profile",
    ),
    Inspiration(
        "Bravais lattice shape",
        "lattice behavior is controlled by shape class and basis choices, not only covolume",
        "Minkowski pressure needs relation-shape labels and low-height wall ownership",
        ("relation_lattice_shape", "bravais_shape", "low_height_wall"),
        "bravais_shape_wall_signature",
    ),
    Inspiration(
        "Lee-Yang extremality",
        "zeros move as a structured locus under parameter changes",
        "replace p0-only maximization by the PGF root curve, nearest-zero radius, and discriminant walls",
        ("root_ear_payload", "root_zero_locus", "pgf_curve", "discriminant_wall"),
        "pgf_zero_trajectory_signature",
    ),
    Inspiration(
        "phi4 quartic density",
        "exp(-lambda*S^4-b*S^2) is controlled by moments, wells, and boundary regimes",
        "quartic energy analogies are legal only after boundary data and zero modes are sidecars",
        ("pde_weak_form", "boundary_data", "zero_mode", "quartic_density_moment_curve"),
        "quartic_moment_wall_profile",
    ),
    Inspiration(
        "strong connectivity iff ear decomposition",
        "global directed connectivity can be certified by an ear-building grammar",
        "observer gluing should be tested by explicit ear certificates rather than asserted globally",
        ("observer_gluing", "ear_decomposition", "strong_connectivity_witness"),
        "observer_ear_certificate_type",
    ),
    Inspiration(
        "factor-critical iff odd ear decomposition",
        "matching-critical structure is exposed by parity-constrained ears",
        "payload ownership may require parity debt, not just a connected chart",
        ("endpoint_owner", "ear_parity", "parity_debt"),
        "odd_ear_payload_parity_debt",
    ),
    Inspiration(
        "series-parallel iff nested ear decomposition",
        "forbidden-minor simplicity appears as nested ear order",
        "proof routes may become tractable only after a nested-ear or branchwidth sidecar is declared",
        ("observer_gluing", "nested_ear_order", "series_parallel_shadow"),
        "nested_ear_branchwidth_shadow",
    ),
]


REFRAMES = [
    Reframe(
        "coordinate_resurrection_sheaf",
        {
            "predicate": 3,
            "repair": 3,
            "minimality": 3,
            "composition": 3,
            "computable": 3,
            "falsifiable": 3,
            "unifies": 3,
            "proof_exit": 3,
        },
        "destroyed coordinates plus legal sections",
        "none essential if sections are proved",
        "instantiate sheaf sections on HYP-2963 packet fibers",
    ),
    Reframe(
        "missing_input_matroid",
        {
            "predicate": 3,
            "repair": 3,
            "minimality": 3,
            "composition": 2,
            "computable": 3,
            "falsifiable": 3,
            "unifies": 2,
            "proof_exit": 2,
        },
        "circuits of mutually missing coordinates",
        "geometric meaning unless labelled by sidecars",
        "compute circuits/cocircuits of repair hypergraph",
    ),
    Reframe(
        "adjoint_quotient_section_pair",
        {
            "predicate": 3,
            "repair": 3,
            "minimality": 2,
            "composition": 3,
            "computable": 2,
            "falsifiable": 2,
            "unifies": 3,
            "proof_exit": 3,
        },
        "which quotient has which legal right adjoint section",
        "finite enumeration unless section is concrete",
        "write quotient/section pairs for root, lattice, time, and circuit maps",
    ),
    Reframe(
        "formal_concept_lattice",
        {
            "predicate": 2,
            "repair": 2,
            "minimality": 3,
            "composition": 2,
            "computable": 3,
            "falsifiable": 3,
            "unifies": 2,
            "proof_exit": 2,
        },
        "sidecar extents and common coordinate intents",
        "dynamics of next operation unless transitions are added",
        "mine stable coordinate intents across current carrier library",
    ),
    Reframe(
        "repair_hypergraph_cover",
        {
            "predicate": 2,
            "repair": 3,
            "minimality": 3,
            "composition": 2,
            "computable": 3,
            "falsifiable": 3,
            "unifies": 2,
            "proof_exit": 2,
        },
        "minimal sidecar covers of shortcut missing vectors",
        "semantic legality if cover is untyped",
        "attach cover certificate to every shortcut audit",
    ),
    Reframe(
        "energy_entropy_dual",
        {
            "predicate": 2,
            "repair": 2,
            "minimality": 1,
            "composition": 2,
            "computable": 2,
            "falsifiable": 2,
            "unifies": 3,
            "proof_exit": 1,
        },
        "PDE/Ising/Lee-Yang entropy-energy analogies",
        "endpoint owner and finite address unless resurrected",
        "test boundary-data sidecars against root and component packets",
    ),
    Reframe(
        "interval_margin_functor",
        {
            "predicate": 2,
            "repair": 2,
            "minimality": 2,
            "composition": 2,
            "computable": 3,
            "falsifiable": 3,
            "unifies": 1,
            "proof_exit": 2,
        },
        "positive intervals and rational-grid resurrection",
        "normalized proof clock unless included",
        "move HYP-3114 from raw time to normalized slow/ruler intervals",
    ),
    Reframe(
        "pde_weak_form_compiler",
        {
            "predicate": 2,
            "repair": 2,
            "minimality": 1,
            "composition": 2,
            "computable": 2,
            "falsifiable": 2,
            "unifies": 2,
            "proof_exit": 2,
        },
        "operator/boundary data as uniform proof sidecar",
        "discrete packet address without finite-address bridge",
        "compile Lee-Yang and component data into weak-form rows",
    ),
    Reframe(
        "raw_scalar_tuning",
        {
            "predicate": 0,
            "repair": 0,
            "minimality": 1,
            "composition": 0,
            "computable": 3,
            "falsifiable": 1,
            "unifies": 0,
            "proof_exit": 0,
        },
        "cheap ranking values",
        "almost every coordinate that HYP-3116 says is essential",
        "demote to scout feature unless repair cover is attached",
    ),
]


TIE_PATH = [
    "coordinate_resurrection_sheaf",
    "missing_input_matroid",
    "adjoint_quotient_section_pair",
    "repair_hypergraph_cover",
    "formal_concept_lattice",
    "interval_margin_functor",
    "pde_weak_form_compiler",
    "energy_entropy_dual",
    "raw_scalar_tuning",
]


def proof_satisfied(coords: set[str]) -> bool:
    return CORE.issubset(coords) and bool(LIVE & coords)


def minimal_repair(shortcut: Shortcut, max_size: int = 5) -> list[tuple[int, tuple[Sidecar, ...], frozenset[str]]]:
    best: list[tuple[int, tuple[Sidecar, ...], frozenset[str]]] = []
    start = set(shortcut.kept)
    if proof_satisfied(start):
        return [(0, tuple(), frozenset(start))]
    for r in range(1, max_size + 1):
        for cover in combinations(SIDECARS, r):
            coords = set(start)
            cost = 0
            for sidecar in cover:
                coords |= sidecar.coords
                cost += sidecar.cost
            if proof_satisfied(coords):
                best.append((cost, cover, frozenset(coords)))
        if best:
            best.sort(key=lambda item: (item[0], len(item[1]), [s.name for s in item[1]]))
            return best
    return []


def all_minimal_global_covers() -> list[tuple[Sidecar, ...]]:
    covers: list[tuple[Sidecar, ...]] = []
    for r in range(1, 6):
        for cover in combinations(SIDECARS, r):
            coords: set[str] = set()
            for sidecar in cover:
                coords |= sidecar.coords
            if proof_satisfied(coords):
                covers.append(cover)
        if covers:
            return sorted(covers, key=lambda cov: [s.name for s in cov])
    return covers


def concept_lattice() -> list[tuple[tuple[str, ...], tuple[str, ...]]]:
    concepts: set[tuple[tuple[str, ...], tuple[str, ...]]] = set()
    sidecar_names = [s.name for s in SIDECARS]
    sidecar_by_name = {s.name: s for s in SIDECARS}
    all_coords = set().union(*(s.coords for s in SIDECARS))
    for r in range(len(SIDECARS) + 1):
        for subset_names in combinations(sidecar_names, r):
            if subset_names:
                intent = set(all_coords)
                for name in subset_names:
                    intent &= sidecar_by_name[name].coords
            else:
                intent = set(all_coords)
            extent = tuple(name for name in sidecar_names if intent.issubset(sidecar_by_name[name].coords))
            closed_intent = set(all_coords)
            if extent:
                for name in extent:
                    closed_intent &= sidecar_by_name[name].coords
            else:
                closed_intent = set()
            concepts.add((tuple(sorted(extent)), tuple(sorted(closed_intent))))
    return sorted(concepts, key=lambda ei: (-len(ei[1]), -len(ei[0]), ei[0]))


def edge_winner(a: Reframe, b: Reframe) -> str:
    axes = sorted(set(a.axes) | set(b.axes))
    av = bv = 0
    for axis in axes:
        da = a.axes.get(axis, 0)
        db = b.axes.get(axis, 0)
        if da > db:
            av += 1
        elif db > da:
            bv += 1
    if av > bv:
        return a.name
    if bv > av:
        return b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_edges() -> dict[str, set[str]]:
    by_name = {r.name: r for r in REFRAMES}
    edges = {r.name: set() for r in REFRAMES}
    for a, b in combinations([r.name for r in REFRAMES], 2):
        winner = edge_winner(by_name[a], by_name[b])
        loser = b if winner == a else a
        edges[winner].add(loser)
    return edges


def directed_3cycles(edges: dict[str, set[str]]) -> int:
    count = 0
    names = sorted(edges)
    for a, b, c in combinations(names, 3):
        triples = [(a, b, c), (a, c, b)]
        for x, y, z in triples:
            if y in edges[x] and z in edges[y] and x in edges[z]:
                count += 1
                break
    return count


def sccs(edges: dict[str, set[str]]) -> list[list[str]]:
    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    low: dict[str, int] = {}
    comps: list[list[str]] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in sorted(edges[v]):
            if w not in indices:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            comps.append(sorted(comp))

    for v in sorted(edges):
        if v not in indices:
            strongconnect(v)
    return comps


def hamiltonian_path_count(edges: dict[str, set[str]]) -> tuple[int, tuple[str, ...] | None]:
    names = [r.name for r in REFRAMES]
    count = 0
    first = None
    for path in permutations(names):
        if all(path[i + 1] in edges[path[i]] for i in range(len(path) - 1)):
            count += 1
            if first is None:
                first = path
    return count, first


def main() -> None:
    print("HYP-3118 coordinate-resurrection sheaf scout -- codex S267")
    print("This is evidence and route selection, not a proof.")
    print()

    print("=" * 100)
    print("REPAIR HYPERGRAPH: SIDECARS AS COORDINATE-RESURRECTION EDGES")
    print("=" * 100)
    for sidecar in SIDECARS:
        print(f"- {sidecar.name} [{sidecar.source}]")
        print(f"  coords={','.join(sorted(sidecar.coords))}")
        print(f"  note={sidecar.note}")
    print()

    print("proof-route predicate used by scout:")
    print("  CORE = finite_address + observer_gluing + endpoint_owner + uniformity")
    print("  plus at least one LIVE sidecar:")
    print("  " + ", ".join(sorted(LIVE)))
    print()

    print("=" * 100)
    print("TWO COMPREHENSIVE MAPS")
    print("=" * 100)
    print("MAP A: quotient -> destroyed coordinate vector -> first proof demand")
    for shortcut in SHORTCUTS:
        core_gap = sorted(CORE - shortcut.kept)
        live_present = sorted(LIVE & shortcut.kept)
        live_gap = sorted(LIVE - shortcut.kept)
        demand = "core repair"
        if live_present:
            demand = "terminalize existing live section"
        elif live_gap:
            demand = "choose one live section after core repair"
        print(f"- {shortcut.name}")
        print(f"  core_gap={','.join(core_gap) if core_gap else 'none'}")
        print(f"  live_present={','.join(live_present) if live_present else 'none'}")
        print(f"  live_options={','.join(live_gap)}")
        print(f"  first_demand={demand}")
    print()
    print("MAP B: external pattern -> LRC resurrection signal")
    for inspiration in INSPIRATIONS:
        print(f"- {inspiration.name}")
        print(f"  external_pattern={inspiration.external_pattern}")
        print(f"  lrc_reframe={inspiration.lrc_reframe}")
        print(f"  resurrects={','.join(inspiration.resurrects)}")
        print(f"  new_signal={inspiration.new_signal}")
    print()

    print("=" * 100)
    print("SHORTCUT REPAIR COVERS")
    print("=" * 100)
    missing_counter: Counter[str] = Counter()
    first_cover_counter: Counter[str] = Counter()
    all_cover_counter: Counter[str] = Counter()
    cover_size_hist: Counter[int] = Counter()
    for shortcut in SHORTCUTS:
        missing = (CORE | LIVE) - shortcut.kept
        for coord in sorted(missing):
            missing_counter[coord] += 1
        repairs = minimal_repair(shortcut)
        print(f"- {shortcut.name} [{shortcut.source}]")
        print(f"  kept={','.join(sorted(shortcut.kept))}")
        print(f"  challenged_assumption={shortcut.challenged_assumption}")
        if not repairs:
            print("  no repair cover found within size bound")
            continue
        cost, cover, coords = repairs[0]
        cover_names = [s.name for s in cover]
        for name in cover_names:
            first_cover_counter[name] += 1
        for _, alt_cover, _ in repairs:
            for name in sorted({s.name for s in alt_cover}):
                all_cover_counter[name] += 1
        cover_size_hist[len(cover)] += 1
        print(f"  minimal_cover_size={len(cover)} cost={cost}")
        print(f"  minimal_cover_count_at_size={len(repairs)}")
        print(f"  cover={ ' + '.join(cover_names) }")
        if len(repairs) > 1:
            for alt_index, (_, alt_cover, _) in enumerate(repairs[1:4], start=2):
                print(f"  alternative_cover_{alt_index}={ ' + '.join(s.name for s in alt_cover) }")
        resurrected = sorted((set(coords) - set(shortcut.kept)) & (CORE | LIVE))
        print(f"  resurrected={','.join(resurrected)}")
    print()
    print("missing_coordinate_frequency_across_shortcuts:")
    for coord, count in sorted(missing_counter.items(), key=lambda item: (-item[1], item[0])):
        print(f"  {coord}: {count}")
    print("repair_sidecar_frequency_in_first_minimal_covers:")
    for name, count in sorted(first_cover_counter.items(), key=lambda item: (-item[1], item[0])):
        print(f"  {name}: {count}")
    print("repair_sidecar_frequency_across_all_minimal_covers:")
    for name, count in sorted(all_cover_counter.items(), key=lambda item: (-item[1], item[0])):
        print(f"  {name}: {count}")
    print(f"cover_size_hist={dict(sorted(cover_size_hist.items()))}")
    print()

    print("=" * 100)
    print("GLOBAL PROOF-ROUTE COVERS")
    print("=" * 100)
    covers = all_minimal_global_covers()
    print(f"minimal_global_cover_count={len(covers)} size={len(covers[0]) if covers else None}")
    for cover in covers[:12]:
        coords = set().union(*(s.coords for s in cover))
        print("  " + " + ".join(s.name for s in cover))
        print(f"    live={','.join(sorted(coords & LIVE))}")
    print()

    print("=" * 100)
    print("FORMAL CONCEPT LATTICE: SIDECAR EXTENTS AND COORDINATE INTENTS")
    print("=" * 100)
    concepts = concept_lattice()
    informative = [c for c in concepts if c[0] and c[1]]
    print(f"concept_count={len(concepts)} informative={len(informative)}")
    for extent, intent in informative[:16]:
        print(f"  extent_size={len(extent):2d} intent_size={len(intent):2d}")
        print(f"    extent={', '.join(extent)}")
        print(f"    intent={', '.join(intent)}")
    print()

    print("=" * 100)
    print("TOURNAMENT ANALYSIS ON ABSTRACT REFRAMES")
    print("=" * 100)
    edges = tournament_edges()
    scores = {name: len(out) for name, out in edges.items()}
    score_hist = Counter(scores.values())
    print("axes=predicate, repair, minimality, composition, computable, falsifiable, unifies, proof_exit")
    print("tie_path=" + " > ".join(TIE_PATH))
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"scores={dict(sorted(scores.items(), key=lambda kv: (-kv[1], kv[0])))}")
    print(f"directed_3cycles={directed_3cycles(edges)}")
    print(f"SCCs={sccs(edges)}")
    hp_count, first_hp = hamiltonian_path_count(edges)
    print(f"Hamiltonian_path_count={hp_count}")
    if first_hp:
        print("selected_retention_path=" + " > ".join(first_hp))
    print()
    print("vertex ledger:")
    by_name = {r.name: r for r in REFRAMES}
    for name in sorted(scores, key=lambda n: (-scores[n], TIE_PATH.index(n))):
        r = by_name[name]
        print(f"- {name}")
        print(f"  preserves: {r.preserves}")
        print(f"  destroys: {r.destroys}")
        print(f"  next_signal: {r.next_signal}")
    print()

    print("=" * 100)
    print("GENERATED HYPOTHESES AND REFRAMES")
    print("=" * 100)
    hypotheses = [
        (
            "H1 coordinate resurrection",
            "Every failed scalar shortcut can be made productive by naming the first destroyed coordinate and a minimal repair cover; the cover is the theorem-facing object.",
        ),
        (
            "H2 proof sheaf",
            "Finite-address, observer-gluing, endpoint-owner, and uniformity are the base sheaf stalk; root/ear, lattice, component, cocycle, state-lift, and PDE sidecars are local sections.",
        ),
        (
            "H3 adjoint quotient",
            "A quotient is legal exactly when it has a computable right-adjoint section on the coordinate demanded by the next proof operation.",
        ),
        (
            "H4 matroid debt",
            "The recurring missing-input vectors form circuits of a matroid-like repair object; finite_address/observer_gluing/endpoint_owner/uniformity behave like coloops in the current model.",
        ),
        (
            "H5 energy-root duality",
            "PDE weak forms and Ising energy are dual to Lee-Yang root confinement only after boundary and zero-mode coordinates are resurrected.",
        ),
        (
            "H6 approximation as section",
            "Irrational/transcendental approximation is not a separate proof route; it is a section from positive interval margin to finite grid address.",
        ),
        (
            "H7 concept-lattice atlas",
            "The stable intents in the sidecar concept lattice should become packet schema fields before another broad analogy is added.",
        ),
        (
            "H8 PGF curve over scalar",
            "The single p0 value is a lossy quotient of a whole PGF root curve; extremality should track zero trajectories, discriminant walls, and root-ear payload transfer.",
        ),
        (
            "H9 ear grammar as gluing",
            "Directed, odd, and nested ear decompositions are not analogies in bulk: they are three possible observer-gluing certificate grammars with different parity and nesting debts.",
        ),
        (
            "H10 Bravais wall ownership",
            "Minkowski/Bravais signals are strongest when they name the relation-lattice shape and low-height wall owner; raw volume is an unowned pressure signal.",
        ),
        (
            "H11 Savitch proof compression",
            "A finite LRC route may need a Savitch-style midpoint ladder: recomputable certificate segments replace storing the entire proof path.",
        ),
        (
            "H12 phi4 stress discipline",
            "The quartic density analogy is useful as a stress test for moment walls and zero modes, but illegal as a proof route until boundary/operator data are retained.",
        ),
    ]
    for label, text in hypotheses:
        print(f"- {label}: {text}")

    print()
    print("new signals:")
    for signal in [
        "destroyed_coordinate_vector",
        "coordinate_resurrection_cover",
        "adjoint_section_status",
        "repair_cover_rank",
        "concept_lattice_intent_id",
        "core_stalk_presence",
        "live_section_type",
        "bottleneck_coordinate_frequency",
        "matroid_circuit_debt",
        "section_gluing_failure_mode",
        "pgf_zero_trajectory_signature",
        "discriminant_wall_crossing_count",
        "observer_ear_certificate_type",
        "odd_ear_payload_parity_debt",
        "nested_ear_branchwidth_shadow",
        "bravais_shape_wall_signature",
        "midpoint_certificate_depth_profile",
        "quartic_moment_wall_profile",
    ]:
        print(f"  - {signal}")


if __name__ == "__main__":
    main()
