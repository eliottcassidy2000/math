#!/usr/bin/env python3
"""Hidden-connection accelerator index for the LRC14 proof stack.

This is a source-audit and proof-carrier index, not a proof of LRC14.

The script collects connector lemmas that were easy to miss because different
sessions used different names for the same proof carrier.  It verifies local
source markers, ranks the connectors, and emits a Tournament Analysis whose
vertices are proof obligations rather than runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from textwrap import fill


REPO = Path(__file__).resolve().parents[1]
NAMESPACE = "HYP-3046 / T1127 / LTI-194 / LTT-092"


@dataclass(frozen=True)
class SourceRef:
    path: str
    markers: tuple[str, ...]


@dataclass(frozen=True)
class Accelerator:
    name: str
    title: str
    sources: tuple[SourceRef, ...]
    hidden_detail: str
    target: str
    action: str
    sidecars: tuple[str, ...]
    guardrail: str
    destroys: str
    scores: tuple[int, int, int, int, int, int, int]


ACCELERATORS: tuple[Accelerator, ...] = (
    Accelerator(
        name="section_grid_exits_are_capacitor_exit_codes",
        title="Residual section exits are the capacitor exit alphabet",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-2996-lrc14-residual-section-packet-grid-verification.md",
                ("horizontal_owner_strip", "vertical_owner_strip", "nested_refinement", "cross_handoff"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-2992-lrc14-haar-product-tile-discrepancy.md",
                ("same_tile_indicator", "cross_handoff", "nested_refinement"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3037-lrc14-residual-capacitor-flow-cuts.md",
                ("exit class: nested_refinement", "exit class: cross_handoff", "same_tile_boundary"),
            ),
        ),
        hidden_detail=(
            "The HYP-3037 capacitor labels are not new.  They are the same exit "
            "alphabet already counted in HYP-2996 and in the exact Haar-product "
            "tile census: owner strips, same-tile boundary atoms, cross handoffs, "
            "and nested refinement."
        ),
        target="Promote a single `zeta_exit_class` / `residual_section_exit` field.",
        action=(
            "Prove capacitor exits by the old residual-section/Haar dictionary "
            "before inventing any new residual labels."
        ),
        sidecars=("zeta_exit_class", "residual_section_exit", "haar_product_exit_class"),
        guardrail="Do not treat capacitor min-cuts as a separate taxonomy until the old exit alphabet fails.",
        destroys="Raw section row identity and raw Haar rectangle positions; keeps the proof exit code.",
        scores=(5, 5, 5, 5, 5, 5, 5),
    ),
    Accelerator(
        name="endpoint_owner_transfer_is_B18Z6_address_lift",
        title="Endpoint-owner transfer is the address lift inside B18Z6",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-3045-lrc14-endpoint-owner-transfer.md",
                (
                    "B18Z6",
                    "external endpoint-owner strip",
                    "owner-transfer carrier",
                    "endpoint_owner_transfer_delta",
                ),
            ),
            SourceRef(
                "05-knowledge/results/lrc14_endpoint_owner_transfer_codex_s208.out",
                (
                    "coarse_endpoint=B18Z6",
                    "external_owner_strip=split",
                    "owner_transfer_carrier",
                    "residue_delta",
                ),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3042-lrc14-owner-strip-filtration.md",
                ("endpoint-owner strip current", "first_surviving_filtration_page"),
            ),
        ),
        hidden_detail=(
            "The endpoint-owner transfer carrier from HYP-3045 is the concrete "
            "address lift hidden inside the coarse B18Z6 residual surface.  It "
            "splits the q=23 petal/covering diagonal and both residual capacitor "
            "pairs by external owner tokens before route labels are used."
        ),
        target="Prove B18Z6 coarse endpoint equality -> external owner-current lift -> residual capacitor split.",
        action=(
            "Use endpoint-owner transfer deltas as the next bridge between "
            "HYP-3038 nested-refinement squares, HYP-3042 owner-strip filtration, "
            "and HYP-3037 capacitor pairs."
        ),
        sidecars=(
            "coarse_endpoint_word",
            "external_endpoint_owner_strip",
            "endpoint_owner_transfer_delta",
            "endpoint_owner_residue_delta",
            "safe_component_owner_stalk",
        ),
        guardrail="B18Z6 is only a shadow; a proof must keep, reconstruct, or annihilate the owner current.",
        destroys="Full packet labels and nonlargest safe interiors; keeps the non-route endpoint owner address.",
        scores=(5, 5, 5, 5, 5, 5, 4),
    ),
    Accelerator(
        name="q23_square_is_nested_refinement_subroutine",
        title="The q=23 square is a concrete nested-refinement subroutine",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-3038-lrc14-q23-drop-add-haar-square.md",
                ("diagonal_doubling_match", "endpoint-owner strip", "nested_refinement_to_q23_diagonal_then_owner_strip"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3037-lrc14-residual-capacitor-flow-cuts.md",
                ("M_q_petal_covering_capacitor", "exit class: nested_refinement"),
            ),
            SourceRef(
                "00-navigation/LRC-TECHNIQUE-INDEX.md",
                ("drop_add_square_id", "exact_M_zeta", "endpoint_owner_strip"),
            ),
        ),
        hidden_detail=(
            "The first displayed nested-refinement capacitor already has a "
            "two-coordinate fixed-margin model.  HYP-3038 shows that the q=23 "
            "petal/covering pair is a drop/add square: off-diagonal rows open, "
            "while the diagonal needs endpoint-owner strips."
        ),
        target="Use HYP-3038 as the local normal form for `nested_refinement` exits.",
        action=(
            "Turn every future nested-refinement report into a request for a "
            "drop/add square, diagonal q-layer, and endpoint-owner split."
        ),
        sidecars=("drop_add_square_id", "diagonal_doubling_match", "exact_M_zeta", "endpoint_owner_strip"),
        guardrail="Exact M alone still mixes the q=23 diagonal; retain owners before declaring a route.",
        destroys="Row/column margins and raw analytic q; keeps mixed-coordinate zeta and owner strip.",
        scores=(5, 4, 4, 5, 5, 4, 5),
    ),
    Accelerator(
        name="owner_strip_is_normal_fan_owner_deletion",
        title="Owner-strip repairs are normal-fan and path-lift owner support",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-3034-lrc14-arc-boundary-path-lift.md",
                ("owner-essential", "owner_deletion_persistence", "closed_arc_h1_owner_support"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3035-lrc14-residual-tooth-atlas.md",
                ("owner strip", "Arc-topology owner strips", "Coarse safe-component stalk owner strips"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3042-lrc14-owner-strip-filtration.md",
                ("endpoint-owner strip current", "owner_strip_page", "first_surviving_filtration_page"),
            ),
            SourceRef(
                "00-navigation/LRC-TECHNIQUE-INDEX.md",
                ("Active-bottleneck normal fan", "closed_arc_h1_owner_support", "owner_deletion_beta1_word"),
            ),
        ),
        hidden_detail=(
            "`owner_strip` is not just a repair label.  The same coordinate is "
            "visible as active bottleneck normal-fan support, safe-component "
            "stalk ownership, HYP-3042's endpoint-owner filtration page, and "
            "owner-deletion persistence of the AP/GW closed H1 representative."
        ),
        target="Prove owner-essential arc cycle -> normal-fan support -> endpoint-owner filtration -> owner-strip descent.",
        action=(
            "Compare `closed_arc_h1_owner_support`, `largest_safe_stalk_owners`, "
            "`normal_fan_active_owner_support`, and `first_surviving_filtration_page` "
            "on the residual ledger."
        ),
        sidecars=(
            "closed_arc_h1_owner_support",
            "largest_safe_stalk_owners",
            "normal_fan_active_owner_support",
            "first_surviving_filtration_page",
        ),
        guardrail="A coarse endpoint word such as `B18Z6` is unsafe without the owner identities.",
        destroys="Raw arc labels and full barcodes; keeps owner support and deletion persistence.",
        scores=(5, 4, 5, 5, 5, 4, 4),
    ),
    Accelerator(
        name="topology_exceptions_are_owner_stalk_primitive_teeth",
        title="Topology exceptions are owner-stalk and primitive-deck teeth",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-3044-lrc14-residual-topology-exception-teeth.md",
                (
                    "topology_exception_fibers=2",
                    "coarse_stalk_splits_all_exceptions=True",
                    "primitive_deck_2_13_splits_all_exceptions=True",
                ),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3035-lrc14-residual-tooth-atlas.md",
                (
                    "arc_topology_compact  13 fibers",
                    "coarse_safe_stalk      2 fibers",
                    "first_tooth_counts={'arc_topology_compact': 13, 'coarse_safe_stalk': 2}",
                ),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3036-lrc14-ramanujan-route-scheduler.md",
                ("primitive_safe_deck_2_13", "first_primitive_safe_q_2_13", "q=14"),
            ),
        ),
        hidden_detail=(
            "HYP-3044 shows that the two compact-topology failures from HYP-3035 "
            "are not new residual families.  They are single-swap owner-stalk "
            "collars at drops 9 and 11, and the primitive q<=13 deck independently "
            "splits the same Q-witness versus covering routes."
        ),
        target="Prove compact-topology failure -> owner-stalk collar -> primitive-deck split.",
        action=(
            "Add topology-exception sidecars before route labels: exception id, "
            "drop speed, stalk key, and first primitive safe q."
        ),
        sidecars=(
            "residual_topology_exception_id",
            "topology_exception_drop",
            "topology_exception_stalk_key",
            "topology_exception_first_primitive_q",
            "topology_then_owner_stalk_rule",
        ),
        guardrail="Do not promote compact topology alone when the two known failures are owner-stalk collars.",
        destroys="Raw topology bucket identity; keeps owner-stalk and primitive-period repair teeth.",
        scores=(5, 5, 5, 5, 5, 4, 4),
    ),
    Accelerator(
        name="primitive_deck_is_exact_period_packet_atlas",
        title="Primitive safe decks are exact-period packet atlases",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-2886-lrc14-exact-period-packet-atlas.md",
                ("exact-period units", "multiplicativity defect", "finite bases are not global closures"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3036-lrc14-ramanujan-route-scheduler.md",
                ("primitive_safe_deck_2_13", "first_primitive_safe_q_2_13", "q=14"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/INDEX.md",
                ("HYP-2979", "Ramanujan exact-period", "no strict q<=42 = 2"),
            ),
        ),
        hidden_detail=(
            "The HYP-3036 primitive safe deck is HYP-2886's exact-period "
            "unit atlas specialized to the post-status residual fibers.  Raw "
            "Ramanujan traces are diagnostic only after the safe-phase inequality "
            "has been evaluated."
        ),
        target="Replace scalar Ramanujan traces by exact-period safe-unit decks.",
        action=(
            "Store exact-period safe counts, first safe q, and CRT multiplicativity "
            "defect before using any Ramanujan or divisor scalar."
        ),
        sidecars=("primitive_safe_deck_2_13", "first_primitive_safe_q_2_13", "crt_multiplicativity_defect"),
        guardrail="Do not collapse exact-period packets to phi(q), c_q(v), or a squarefree capacity.",
        destroys="Individual unit residues only after safe counts are retained; keeps exact-period packet status.",
        scores=(5, 5, 5, 4, 5, 5, 4),
    ),
    Accelerator(
        name="q14_guardrail_is_covering_grid_theorem",
        title="The q=14 guardrail is the old covering-grid threshold",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/INDEX.md",
                ("THM-523", "q(S)=14", "q>14"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3036-lrc14-ramanujan-route-scheduler.md",
                ("The `q=14` layer", "boundary/covering-moment layer", "primitive_safe_deck_2_13"),
            ),
        ),
        hidden_detail=(
            "HYP-3036's q=14 warning is exactly the old THM-523/HYP-2917 "
            "divisibility-threshold split: q<14 gives a direct witness, q=14 "
            "is tight boundary, and q>14 is the covering hard core."
        ),
        target="Route zero primitive deck rows through q(S) and covering-grid debt.",
        action=(
            "Add the divisibility-threshold field next to the primitive safe deck "
            "so q=14 never masquerades as a direct Q-witness split."
        ),
        sidecars=("divisibility_threshold_qS", "covering_grid_debt", "q14_boundary_layer"),
        guardrail="Positive mass at q=14 belongs to boundary/covering scheduling, not q<=13 witness scheduling.",
        destroys="Raw denominator optimism; keeps the old non-covering/covering theorem split.",
        scores=(4, 5, 4, 4, 5, 5, 5),
    ),
    Accelerator(
        name="capacitor_min_cut_is_sidechannel_repair_cochain",
        title="Capacitor min-cuts are first nonzero repair cochains",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-3027-lrc14-sidechannel-repair-ladder.md",
                ("repair ladder", "first nonzero repair cochain", "guarded non-route signature"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3031-lrc14-haar-tile-repair-ladder-synthesis.md",
                ("repair ladder", "zeta", "first tested non-route splitter"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3037-lrc14-residual-capacitor-flow-cuts.md",
                ("first cut", "residual_capacitor_id", "first_cut_stage"),
            ),
        ),
        hidden_detail=(
            "The capacitor language is a flow picture for the same object HYP-3027 "
            "called the first nonzero repair cochain, with HYP-3031 naming the "
            "local Haar zeta coordinate."
        ),
        target="Prove `first_cut_stage` by the side-channel repair theorem.",
        action=(
            "Keep `residual_capacitor_id` only as an index; discharge it by the "
            "first nonzero legal repair cochain."
        ),
        sidecars=("residual_capacitor_id", "first_cut_stage", "first_nonzero_repair_cochain"),
        guardrail="A min-cut is not a scalar rank; exact scale and topology are incomparable cheap cuts.",
        destroys="Flow-graph plate identity; keeps repair stage and side-channel class.",
        scores=(5, 4, 5, 5, 4, 4, 5),
    ),
    Accelerator(
        name="pair_good_decoys_are_blocker_decks",
        title="Pair-good decoys are generated blocker decks",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-3022-lrc14-pair-good-decoy-barcode-normal-fan.md",
                ("pair_good_decoy_blocker_deck", "pair_good_decoy_generation_key", "barcode_relation"),
            ),
            SourceRef(
                "00-navigation/LRC-TECHNIQUE-INDEX.md",
                ("Pair-good decoy barcode/normal-fan refinement", "blocker-deck grammar", "Replaces raw decoy counts"),
            ),
        ),
        hidden_detail=(
            "Large pair-good decoy counts are generated switch objects, not proof "
            "objects.  The useful data is the lane/shell/blocker deck and its "
            "barcode/normal-fan relation."
        ),
        target="Discharge false pair switches by blocker-deck grammar.",
        action=(
            "Attach decoy generation fields to residual ledgers before using "
            "pair-good predicates in a proof carrier."
        ),
        sidecars=("pair_good_decoy_generation_key", "pair_good_decoy_blocker_deck", "pair_good_decoy_normal_fan_relation"),
        guardrail="Raw decoy counts reward repeated shadows of the same blocker deck.",
        destroys="Individual decoy event multiplicity; keeps generator grammar and blocker owners.",
        scores=(4, 4, 3, 5, 5, 4, 3),
    ),
    Accelerator(
        name="perfect_divisor_lanes_explain_prime_power_blindness",
        title="Perfect-number lanes explain prime-power blindness",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-3013-lrc14-perfect-number-packet-merge.md",
                ("abundancy_defect", "prime_q_flag", "q=27=3^3"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3032-lrc14-analytic-sieve-clock-bridge.md",
                ("squarefree blindness", "fibbinary_first13", "q=27=3^3"),
            ),
            SourceRef(
                "00-navigation/LRC-TECHNIQUE-INDEX.md",
                ("Perfect-number divisor packet merge", "abundancy_defect", "squarefree blindness"),
            ),
        ),
        hidden_detail=(
            "The q=27/C27 and q=25/fibbinary failures of squarefree analytic "
            "capacity are not analytic accidents.  They are prime-power and "
            "divisor-lattice packet lanes that HYP-3013 already said must be "
            "kept before product analogies are used."
        ),
        target="Repair analytic squarefree blindness with divisor-lattice packet fields.",
        action=(
            "Preserve prime/composite q, factorization, abundancy defect, and "
            "exact period before invoking sieve clocks."
        ),
        sidecars=("prime_q_flag", "divisor_lattice_factorization", "abundancy_defect", "squarefree_blindness_report"),
        guardrail="Squarefree capacity is a meter, not an LRC quotient.",
        destroys="Raw product scalar and mu2/phi-only view; keeps prime-power exact-period lanes.",
        scores=(3, 4, 3, 4, 5, 3, 4),
    ),
    Accelerator(
        name="cocycle_exactness_is_global_contract",
        title="All recent sidecars are cochains in one exactness contract",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-2995-lrc14-cocycle-carrier-atlas.md",
                ("omega_Q", "exact as a coboundary", "named F7/THM-572 residual"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3006-lrc14-cocycle-sheaf-exactness.md",
                ("exactness at `C1`", "emitted local cocycles", "named residual summand"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/HYP-3031-lrc14-haar-tile-repair-ladder-synthesis.md",
                ("mixed cocycle", "annihilation, descent, boundary", "named residual debt"),
            ),
        ),
        hidden_detail=(
            "The recent sidecars are cochains: zeta switches, endpoint currents, "
            "Ramanujan phases, Fejer debts, pair tensions, owner strips, and "
            "capacitor cuts are all instances of `omega_Q` for a quotient Q."
        ),
        target="Use cocycle exactness as the global proof contract.",
        action=(
            "For each protected status fiber, require `omega_Q` to be retained, "
            "reconstructed, exact, dual-annihilated, descended, boundary-equal, "
            "or emitted as named F7/THM-572 debt."
        ),
        sidecars=("omega_Q_class", "coboundary_certificate", "dual_annihilator", "named_residual_summand"),
        guardrail="A quotient that only says it separates packets but not what happened to the emitted cochain is incomplete.",
        destroys="Local vocabulary differences; keeps the exactness obligation.",
        scores=(4, 5, 5, 3, 5, 5, 2),
    ),
    Accelerator(
        name="exact_period_atlas_prevents_finite_denominator_shortcuts",
        title="Exact-period packets turn finite denominator failures into charts",
        sources=(
            SourceRef(
                "05-knowledge/hypotheses/HYP-2886-lrc14-exact-period-packet-atlas.md",
                ("finite bases are not global closures", "divload_B90", "Multiplicativity defect"),
            ),
            SourceRef(
                "05-knowledge/hypotheses/INDEX.md",
                ("THM-566 refutes", "finite bases are atlas charts", "q-witness lemma"),
            ),
            SourceRef(
                "00-navigation/LRC-TECHNIQUE-INDEX.md",
                ("Ramanujan exact-period projectors", "Exact-period denominators", "Restore prime-power"),
            ),
        ),
        hidden_detail=(
            "The old bounded-denominator no-go and the newer primitive-deck "
            "success are compatible: finite bases are charts.  They accelerate "
            "a proof only when denominator deaths and exact-period survivors "
            "are stored as packet fields."
        ),
        target="Use exact-period atlas charts as schedulers, not shortcuts.",
        action=(
            "When a finite witness basis fails, record the killed denominators, "
            "first surviving exact period, and route handoff."
        ),
        sidecars=("killed_denominator_set", "first_surviving_exact_period", "finite_exception_budget", "exact_period_route_handoff"),
        guardrail="No fixed finite denominator list proves the global theorem.",
        destroys="Universal bounded-D hope; keeps chart-local exact-period data.",
        scores=(3, 5, 4, 4, 5, 5, 3),
    ),
)


TIE_PATH = tuple(acc.name for acc in ACCELERATORS)


@lru_cache(maxsize=None)
def read_source(path: str) -> str | None:
    full = REPO / path
    if not full.exists():
        return None
    return full.read_text(encoding="utf-8", errors="replace")


def audit_source(source: SourceRef) -> tuple[list[str], list[str], bool]:
    text = read_source(source.path)
    if text is None:
        return [], list(source.markers), False
    hits = [marker for marker in source.markers if marker in text]
    misses = [marker for marker in source.markers if marker not in text]
    return hits, misses, True


def audit_accelerator(acc: Accelerator) -> tuple[int, int, list[str]]:
    hits = 0
    misses = 0
    lines: list[str] = []
    for source in acc.sources:
        source_hits, source_misses, exists = audit_source(source)
        hits += len(source_hits)
        misses += len(source_misses)
        status = "ok" if exists else "missing-file"
        lines.append(
            f"{source.path} [{status}] hit={len(source_hits)}/{len(source.markers)}"
            + (f" miss={','.join(source_misses)}" if source_misses else "")
        )
    return hits, misses, lines


def rank_score(acc: Accelerator) -> tuple[int, int, int]:
    hits, misses, _lines = audit_accelerator(acc)
    return (sum(acc.scores) + hits - 2 * misses, hits, -TIE_PATH.index(acc.name))


def pair_winner(a: Accelerator, b: Accelerator) -> Accelerator:
    a_votes = 0
    b_votes = 0
    for av, bv in zip(a.scores, b.scores):
        if av > bv:
            a_votes += 1
        elif bv > av:
            b_votes += 1
    if a_votes > b_votes:
        return a
    if b_votes > a_votes:
        return b
    return a if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b


def build_tournament(vertices: tuple[Accelerator, ...]) -> dict[str, set[str]]:
    edges: dict[str, set[str]] = {acc.name: set() for acc in vertices}
    for a, b in combinations(vertices, 2):
        winner = pair_winner(a, b)
        loser = b if winner is a else a
        edges[winner.name].add(loser.name)
    return edges


def directed_3cycles(vertices: tuple[Accelerator, ...], edges: dict[str, set[str]]) -> list[tuple[str, str, str]]:
    cycles: list[tuple[str, str, str]] = []
    for a, b, c in combinations([v.name for v in vertices], 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            cycles.append((a, b, c))
        elif c in edges[a] and b in edges[c] and a in edges[b]:
            cycles.append((a, c, b))
    return cycles


def scc_sizes(vertices: tuple[Accelerator, ...], edges: dict[str, set[str]]) -> list[int]:
    names = [v.name for v in vertices]
    reverse: dict[str, set[str]] = {name: set() for name in names}
    for a, outs in edges.items():
        for b in outs:
            reverse[b].add(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in edges[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str, bag: list[str]) -> None:
        seen.add(v)
        bag.append(v)
        for w in reverse[v]:
            if w not in seen:
                rdfs(w, bag)

    for name in reversed(order):
        if name not in seen:
            bag: list[str] = []
            rdfs(name, bag)
            sizes.append(len(bag))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(vertices: tuple[Accelerator, ...], edges: dict[str, set[str]]) -> int:
    names = [v.name for v in vertices]
    n = len(names)
    idx = {name: i for i, name in enumerate(names)}
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if names[nxt] in edges[names[last]]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, i), 0) for i in range(n))


def greedy_path(vertices: tuple[Accelerator, ...], edges: dict[str, set[str]]) -> list[str]:
    remaining = set(v.name for v in vertices)
    path: list[str] = []
    while remaining:
        best = max(
            remaining,
            key=lambda name: (len(edges[name] & remaining), -TIE_PATH.index(name)),
        )
        path.append(best)
        remaining.remove(best)
    return path


def wrap(text: str, indent: str = "  ", width: int = 104) -> str:
    return fill(text, width=width, initial_indent=indent, subsequent_indent=indent)


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print(
        wrap(
            "Alternate vertex sets considered: runners, gaps, residual fibers, old hypothesis files, "
            "sidecars, proof obligations, capacitor pairs, Haar exits, denominator layers, cover arcs, "
            "Fourier modes, matroid cocircuits, and wall-crossing events."
        )
    )
    print(
        wrap(
            "Chosen vertices: hidden connector lemmas/proof obligations.  This preserves the LRC predicate "
            "of boundary/open status plus theorem-route schedulability and tells which older proof object "
            "can discharge which new residual obligation."
        )
    )
    print(
        wrap(
            "Destroyed by the quotient: row identity, raw runner names, exact times, and sometimes route "
            "labels.  Retained: status, exit class, sidecar field, and the proof obligation."
        )
    )
    print()


def print_source_audit(ranked: list[Accelerator]) -> None:
    print("[1] Source marker audit")
    total_hits = 0
    total_misses = 0
    for acc in ranked:
        hits, misses, lines = audit_accelerator(acc)
        total_hits += hits
        total_misses += misses
        print(f"{acc.name}: markers_hit={hits} missing={misses}")
        for line in lines:
            print(f"  - {line}")
    print(f"marker_totals: hit={total_hits} missing={total_misses}")
    print()


def print_ranked(ranked: list[Accelerator]) -> None:
    print("[2] Ranked accelerators")
    for i, acc in enumerate(ranked, 1):
        print(f"{i}. {acc.name} :: {acc.title}")
        print(f"  score_vector={acc.scores} audit_rank={rank_score(acc)[0]}")
        print(wrap(f"hidden_detail: {acc.hidden_detail}"))
        print(wrap(f"target: {acc.target}"))
        print(wrap(f"action: {acc.action}"))
        print(f"  sidecars={', '.join(acc.sidecars)}")
        print(wrap(f"guardrail: {acc.guardrail}"))
        print(wrap(f"destroys: {acc.destroys}"))
    print()


def print_tournament(vertices: tuple[Accelerator, ...]) -> None:
    edges = build_tournament(vertices)
    score_hist = Counter(len(outs) for outs in edges.values())
    cycles = directed_3cycles(vertices, edges)
    path = greedy_path(vertices, edges)
    print("[3] Tournament Analysis")
    print("vertices=hidden connector lemmas / proof obligations, not runners or packets")
    print(
        "observable=(recent_stack_reach, legacy_evidence, route_power, sidecar_readiness, "
        "scalar_guardrail, family_transfer, low_proof_cost)"
    )
    print("switch=orient A->B by majority of observable coordinates; tie path is source order")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={len(cycles)}")
    if cycles:
        print("cycle_examples=" + " ; ".join(" > ".join(cyc) + " > " + cyc[0] for cyc in cycles[:5]))
    print(f"scc_sizes={scc_sizes(vertices, edges)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(vertices, edges)}")
    print("high_retention_path=" + " > ".join(path))
    print()


def print_readout(ranked: list[Accelerator]) -> None:
    print("[4] Proof acceleration readout")
    print(
        wrap(
            "Highest leverage: collapse HYP-3037 capacitor exit classes to the HYP-2996/HYP-2992 "
            "residual-section/Haar exit alphabet, read HYP-3045 as the B18Z6 endpoint-owner "
            "address lift, then use HYP-3038 as the concrete q=23 nested-refinement normal form."
        )
    )
    print(
        wrap(
            "Owner-strip work should route through HYP-3042's endpoint-owner filtration page plus "
            "HYP-3045 endpoint-owner transfer and HYP-3018/HYP-3034 owner supports, not through "
            "coarse endpoint counts.  HYP-3044's two topology exceptions should route through "
            "owner-stalk and primitive-deck teeth.  The q=14 primitive-deck guardrail should route "
            "through THM-523/HYP-2917.  Pair-good decoys should be blocker decks, and analytic "
            "squarefree blindness should be repaired by exact-period/divisor-lattice fields."
        )
    )
    sidecars: list[str] = []
    for acc in ranked:
        for field in acc.sidecars:
            if field not in sidecars:
                sidecars.append(field)
    print("packet_sidecar_todo=" + ", ".join(sidecars))
    print()


def main() -> None:
    print("=== LRC14 hidden connection accelerators S209 ===")
    print(f"namespace={NAMESPACE}")
    print("script=04-computation/lrc14_hidden_connection_accelerators_codex_s209.py")
    print()
    ranked = sorted(ACCELERATORS, key=rank_score, reverse=True)
    print_assumption_challenge()
    print_source_audit(ranked)
    print_ranked(ranked)
    print_tournament(tuple(ranked))
    print_readout(ranked)


if __name__ == "__main__":
    main()
