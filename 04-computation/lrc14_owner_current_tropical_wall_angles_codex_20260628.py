#!/usr/bin/env python3
"""HYP-3402 scout: endpoint-owner currents plus tropical height walls.

This is a proof-angle scout, not a proof of LRC14.  It deliberately avoids
rerunning the HYP-3300 observability/Morse program, the HYP-3301 sheaf/cusp
program, or the HYP-3310/HYP-3311 residue-word instantiation.  The question is:
what should be tried after the first actual-packet ambiguity is repaired by a
nonunit residue word, but before residue-only data is mistaken for a theorem?

Two different carriers are tested:

1. endpoint-owner boundary current: treat endpoint owners and theorem exits as
   sources/sinks in a finite current ledger;
2. tropical height/discriminant wall: treat same-residue or same-v2 covering
   flexes as valuation bends in a Newton/secondary fan.

Tournament Analysis is over proof carriers and hidden sidecars, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


OBLIGATION_WEIGHT = {
    "endpoint_owner_memory": 9,
    "off_grid_floor": 8,
    "height_flex_legality": 8,
    "residue_magnitude_split": 7,
    "open_boundary_status": 7,
    "finite_checkability": 7,
    "quotient_descent_legality": 7,
    "green_current_certificate": 6,
    "normal_fan_wall": 6,
    "state_lift_exit": 6,
    "analytic_zero_control": 5,
    "bulk_floor_margin": 5,
    "odd_sign_sidecar": 5,
    "formal_packet_schema": 4,
}


KEY_FILES = [
    "00-navigation/LRC-LENS-MAP.md",
    "00-navigation/LRC-TECHNIQUE-INDEX.md",
    "00-navigation/LRC-TOURNAMENT-TECHNIQUE-INDEX.md",
    "00-navigation/OPEN-QUESTIONS.md",
    "05-knowledge/hypotheses/INDEX.md",
    "05-knowledge/hypotheses/HYP-3311-lrc14-actual-packet-sheaf-instantiation.md",
    "05-knowledge/hypotheses/HYP-3311-lrc14-crt-galois-sidecar-audit.md",
    "05-knowledge/hypotheses/HYP-3400-lrc14-shadow-charge-conservation-atlas.md",
    "07-reflections/lrc14-actual-packet-sheaf-instantiation-codex-20260628.md",
    "07-reflections/the-census-factors-via-crt-7-adic-residue-c3-skeleton-times-2-adic-magnitude-doubling-hinge.md",
    "07-reflections/lrc14-shadow-charge-conservation-atlas-codex-20260628.md",
    "comms/POKE-COORDINATION.md",
]


@dataclass(frozen=True)
class Carrier:
    name: str
    family: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    cross_terms: tuple[str, ...]
    theorem_test: str
    risk: str
    priority: int


CARRIERS = [
    Carrier(
        name="endpoint_owner_boundary_current",
        family="Hodge/Farkas/flow proof carrier",
        preserves=frozenset(
            {
                "endpoint_owner_memory",
                "open_boundary_status",
                "quotient_descent_legality",
                "finite_checkability",
                "green_current_certificate",
                "state_lift_exit",
                "formal_packet_schema",
            }
        ),
        destroys=frozenset({"height_flex_legality", "analytic_zero_control"}),
        cross_terms=(
            "endpoint owner",
            "owner current",
            "Farkas",
            "Green current",
            "Thomson",
            "Hodge",
            "boundary cocircuit",
            "safe component",
        ),
        theorem_test=(
            "Turn each residual fiber into a signed owner-current boundary. "
            "A legal quotient must conserve the current, route it to a dual "
            "Farkas/Green certificate, stop at AP/GW boundary H1, or name the "
            "first owner-current leak."
        ),
        risk="owner currents can classify exits without proving the positive floor",
        priority=16,
    ),
    Carrier(
        name="tropical_height_discriminant_wall",
        family="Newton polygon / secondary fan proof carrier",
        preserves=frozenset(
            {
                "height_flex_legality",
                "off_grid_floor",
                "residue_magnitude_split",
                "normal_fan_wall",
                "finite_checkability",
                "analytic_zero_control",
                "quotient_descent_legality",
            }
        ),
        destroys=frozenset({"endpoint_owner_memory", "odd_sign_sidecar"}),
        cross_terms=(
            "height flex",
            "same-residue",
            "v2",
            "Newton",
            "tropical",
            "secondary fan",
            "Plucker",
            "discriminant",
            "off-grid floor",
        ),
        theorem_test=(
            "Treat covering flex as a min-plus bend locus.  Same-residue or "
            "same-v2 rows must either cross a tropical wall with positive "
            "off-grid floor, land on the AP/GW 12->24 hinge, or expose a named "
            "height/discriminant debt."
        ),
        risk="a tropical wall without endpoint owners can prove the wrong exit",
        priority=15,
    ),
    Carrier(
        name="owner_valuation_bicurrent",
        family="two-layer current plus valuation carrier",
        preserves=frozenset(
            {
                "endpoint_owner_memory",
                "height_flex_legality",
                "off_grid_floor",
                "residue_magnitude_split",
                "quotient_descent_legality",
                "finite_checkability",
            }
        ),
        destroys=frozenset({"analytic_zero_control"}),
        cross_terms=(
            "endpoint owner",
            "v2",
            "height",
            "cover signature",
            "residue word",
            "owner current",
        ),
        theorem_test=(
            "Pair the endpoint-owner boundary current with the first valuation "
            "owner of the nonunit cover.  Use it as the minimal fallback when "
            "residue words or v2 words alone stop separating enlarged packets."
        ),
        risk="may be too close to a sidecar product unless it yields a theorem",
        priority=4,
    ),
    Carrier(
        name="cluster_mutation_wall_crossing",
        family="cluster/scattering diagram carrier",
        preserves=frozenset(
            {
                "normal_fan_wall",
                "height_flex_legality",
                "residue_magnitude_split",
                "quotient_descent_legality",
                "finite_checkability",
            }
        ),
        destroys=frozenset({"endpoint_owner_memory", "bulk_floor_margin"}),
        cross_terms=("cluster", "mutation", "scattering", "wall crossing", "Plucker", "fan"),
        theorem_test=(
            "Model local packet changes as mutations in a finite scattering "
            "diagram.  Illegal loops should have nonzero wall charge, while "
            "legal AP/GW and positive-open routes have consistent wall products."
        ),
        risk="beautiful wall language can hide the original LRC predicate",
        priority=9,
    ),
    Carrier(
        name="ramanujan_projector_backend",
        family="exact-period harmonic certificate carrier",
        preserves=frozenset(
            {
                "analytic_zero_control",
                "bulk_floor_margin",
                "finite_checkability",
                "quotient_descent_legality",
                "formal_packet_schema",
            }
        ),
        destroys=frozenset({"endpoint_owner_memory", "height_flex_legality"}),
        cross_terms=("Ramanujan", "exact-period", "projector", "period", "Haar", "Fejer"),
        theorem_test=(
            "Rebuild positive-open exits by exact-period projectors, but require "
            "owner and height sidecars before using a scalar harmonic certificate."
        ),
        risk="harmonic positivity is blind to the first owner/height leak",
        priority=6,
    ),
    Carrier(
        name="boolean_minterm_certificate_dag",
        family="circuit/Farkas minimal-certificate carrier",
        preserves=frozenset(
            {
                "finite_checkability",
                "formal_packet_schema",
                "quotient_descent_legality",
                "open_boundary_status",
            }
        ),
        destroys=frozenset({"endpoint_owner_memory", "height_flex_legality", "off_grid_floor"}),
        cross_terms=("circuit", "minterm", "Farkas", "missing input", "Boolean", "proof DAG"),
        theorem_test=(
            "Translate every exit into a minimal certificate minterm; any "
            "missing endpoint/height input becomes an explicit lower-bound gate."
        ),
        risk="certificate size is not proof content unless the minterms carry packet fields",
        priority=3,
    ),
    Carrier(
        name="raw_residue_word_reuse",
        family="do-not-repeat residue-only carrier",
        preserves=frozenset({"residue_magnitude_split", "finite_checkability"}),
        destroys=frozenset(
            {
                "endpoint_owner_memory",
                "height_flex_legality",
                "off_grid_floor",
                "analytic_zero_control",
                "odd_sign_sidecar",
            }
        ),
        cross_terms=("residue word", "mod14", "C6", "C3"),
        theorem_test="Useful bank-local sidecar, but no longer a new proof angle after HYP-3311.",
        risk="known to be bank-local; same-residue height moves remain live",
        priority=-20,
    ),
]


@dataclass(frozen=True)
class MixedRow:
    name: str
    kernel: str
    residue_word: tuple[int, ...]
    v2_word: tuple[int, ...]
    owner_current: tuple[str, ...]
    tropical_wall: tuple[str, ...]
    bicurrent: tuple[str, ...]

    @property
    def coarse_base(self) -> tuple[str, ...]:
        return ("eq14", "six_unit_boundary", "strict_safe_positive", "no_state_lift")


MIXED_ROWS = [
    MixedRow(
        "drop6 fattening core add180",
        "positive-Haar-open",
        (2, 4, 7, 8, 10, 12, 12),
        (0, 1, 1, 2, 2, 2, 3),
        ("covering_floor_source", "drop6_owner", "add180_sink"),
        ("height_bend", "drop6_to_180", "apex7_active"),
        ("drop6_owner", "v2_owner_12", "covering_floor"),
    ),
    MixedRow(
        "floor-odd GW iso impostor",
        "positive-Haar-open",
        (2, 4, 6, 7, 8, 10, 10),
        (0, 1, 1, 1, 2, 3, 3),
        ("covering_floor_source", "odd_floor_owner", "gw_shadow_sink"),
        ("odd_floor_bend", "gw_shadow", "apex7_active"),
        ("odd_floor_owner", "v2_owner_10", "covering_floor"),
    ),
    MixedRow(
        "magnitude liar 12->96",
        "positive-Haar-open",
        (2, 4, 6, 7, 8, 10, 12),
        (0, 1, 1, 1, 2, 3, 5),
        ("covering_floor_source", "12_owner", "96_height_sink"),
        ("height_bend", "12_to_96", "v2_depth5"),
        ("12_owner", "v2_owner_96", "covering_floor"),
    ),
    MixedRow(
        "P10+GW",
        "unit-petal-named",
        (2, 4, 6, 6, 7, 8, 10),
        (0, 1, 1, 2, 2, 3, 3),
        ("unit_petal_source", "P10_owner", "GW_boundary_sink"),
        ("unit_petal_bend", "P10_plus_GW", "unit_boundary"),
        ("P10_owner", "v2_owner_10", "unit_petal"),
    ),
    MixedRow(
        "petal 10->20",
        "unit-petal-named",
        (2, 4, 6, 6, 7, 8, 12),
        (0, 1, 1, 2, 2, 2, 3),
        ("unit_petal_source", "10_owner", "20_boundary_sink"),
        ("unit_petal_bend", "10_to_20", "unit_boundary"),
        ("10_owner", "v2_owner_20", "unit_petal"),
    ),
    MixedRow(
        "petal 13->26",
        "unit-petal-named",
        (2, 4, 6, 7, 8, 10, 12, 12),
        (0, 1, 1, 1, 1, 2, 2, 3),
        ("unit_petal_source", "13_owner", "26_boundary_sink"),
        ("unit_petal_bend", "13_to_26", "unit_boundary"),
        ("13_owner", "v2_owner_26", "unit_petal"),
    ),
    MixedRow(
        "unit petal splice drop(10,13)->add(20,26)",
        "unit-petal-named",
        (2, 4, 6, 6, 7, 8, 12, 12),
        (0, 1, 1, 1, 2, 2, 2, 3),
        ("unit_petal_source", "10_13_owner_pair", "20_26_boundary_sink"),
        ("unit_petal_bend", "splice_10_13_to_20_26", "unit_boundary"),
        ("10_13_owner_pair", "v2_owner_20_26", "unit_petal"),
    ),
]


def corpus_text() -> str:
    chunks = []
    for rel in KEY_FILES:
        path = ROOT / rel
        if path.exists():
            chunks.append(path.read_text(encoding="utf-8", errors="ignore").lower())
    return "\n".join(chunks)


def corpus_hits(carrier: Carrier, text: str) -> int:
    return sum(text.count(term.lower()) for term in carrier.cross_terms)


def score(carrier: Carrier, text: str) -> int:
    base = sum(OBLIGATION_WEIGHT[key] for key in carrier.preserves)
    penalty = sum(OBLIGATION_WEIGHT[key] for key in carrier.destroys)
    hit_bonus = min(corpus_hits(carrier, text), 80) // 2
    return base - penalty + carrier.priority + hit_bonus


def fibers(rows: list[MixedRow], key_name: str) -> dict[tuple[object, ...], list[MixedRow]]:
    out: dict[tuple[object, ...], list[MixedRow]] = defaultdict(list)
    for row in rows:
        out[getattr(row, key_name)].append(row)
    return out


def mixed_fiber_count(rows: list[MixedRow], key_name: str) -> int:
    total = 0
    for fiber in fibers(rows, key_name).values():
        if len({row.kernel for row in fiber}) > 1:
            total += 1
    return total


def max_target_width(rows: list[MixedRow], key_name: str) -> int:
    width = 0
    for fiber in fibers(rows, key_name).values():
        width = max(width, len({row.kernel for row in fiber}))
    return width


def hamiltonian_paths(order: list[str], edges: set[tuple[str, str]]) -> int:
    # The score plus tie gauge induces one total priority order, so only that
    # path is tested for this report.
    return int(all((a, b) in edges for a, b in zip(order, order[1:])))


def directed_triangles(names: list[str], edges: set[tuple[str, str]]) -> int:
    count = 0
    for a, b, c in combinations(names, 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            count += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            count += 1
    return count


def minimal_sidecar_sets(rows: list[MixedRow], sidecars: list[str]) -> list[tuple[str, ...]]:
    good: list[tuple[str, ...]] = []
    for size in range(1, len(sidecars) + 1):
        for combo in combinations(sidecars, size):
            grouped: dict[tuple[object, ...], list[MixedRow]] = defaultdict(list)
            for row in rows:
                key = row.coarse_base
                for sidecar in combo:
                    key += (getattr(row, sidecar),)
                grouped[key].append(row)
            if all(len({row.kernel for row in fiber}) == 1 for fiber in grouped.values()):
                good.append(combo)
        if good:
            return good
    return []


def main() -> None:
    text = corpus_text()
    ranked = sorted(CARRIERS, key=lambda carrier: (-score(carrier, text), carrier.name))
    scores = {carrier.name: score(carrier, text) for carrier in CARRIERS}
    names = [carrier.name for carrier in ranked]
    edges = {(a.name, b.name) for a, b in combinations(ranked, 2)}
    sidecars = ["residue_word", "v2_word", "owner_current", "tropical_wall", "bicurrent"]

    print("HYP-3402 LRC14 OWNER-CURRENT + TROPICAL-WALL PROOF-ANGLE SCOUT")
    print("=" * 78)
    print("status=synthesis / theorem-target router; not an LRC14 proof")
    print("new_handles=HYP-3402,T1363,LTI-363,LTT-263")
    print()
    print("SELECTED REMAINING ANGLES")
    selected = [
        carrier
        for carrier in ranked
        if carrier.name
        in {"endpoint_owner_boundary_current", "tropical_height_discriminant_wall"}
    ]
    for carrier in selected:
        print(f"  {carrier.name} score={scores[carrier.name]}")
        print(f"    family={carrier.family}")
        print(f"    preserves={','.join(sorted(carrier.preserves))}")
        print(f"    destroys={','.join(sorted(carrier.destroys))}")
        print(f"    test={carrier.theorem_test}")
        print(f"    risk={carrier.risk}")
    print()
    print("FULL CARRIER RANKING")
    for carrier in ranked:
        print(
            f"  {carrier.name:34s} score={scores[carrier.name]:4d} "
            f"corpus_hits={corpus_hits(carrier, text):4d}"
        )
    print()
    print("COMPOSITE BRIDGE")
    bridge = next(carrier for carrier in ranked if carrier.name == "owner_valuation_bicurrent")
    print(f"  {bridge.name} score={scores[bridge.name]}")
    print("    use only after one of the two separate angles finds the first leak")
    print("    or when enlarged packets need endpoint-owner and valuation owner together.")
    print()

    print("HYP-3311 MIXED-FIBER SIDECAR AUDIT")
    print("  source=the unique actual-packet coarse fiber from HYP-3311")
    print(f"  rows={len(MIXED_ROWS)} kernels={dict(Counter(row.kernel for row in MIXED_ROWS))}")
    for sidecar in ["coarse_base", *sidecars]:
        print(
            f"  {sidecar:15s} fibers={len(fibers(MIXED_ROWS, sidecar)):2d} "
            f"mixed_kernel_fibers={mixed_fiber_count(MIXED_ROWS, sidecar):1d} "
            f"max_target_width={max_target_width(MIXED_ROWS, sidecar)}"
        )
    print()
    print("  minimal single/few sidecar repairs over coarse_base:")
    for combo in minimal_sidecar_sets(MIXED_ROWS, sidecars):
        print(f"    {combo}")
    print()

    print("OWNER-CURRENT READOUT")
    print("  owner_current treats theorem exits as signed boundary currents:")
    print("    unit-petal-named rows are unit-boundary inflows;")
    print("    positive-Haar-open rows are covering-floor outflows.")
    print("  theorem target: the first endpoint-owner leak in an enlarged HYP-2963")
    print("  bank must become a Farkas/Green dual certificate, AP/GW boundary H1,")
    print("  forbidden H7 lift, or named residual owner-current debt.")
    print()

    print("TROPICAL-WALL READOUT")
    print("  v2_word alone repeats the HYP-3311 warning:")
    print("    drop6 fattening core add180 and petal 10->20 share the same v2 word")
    print("    but have different theorem exits.")
    print("  tropical_wall adds the first bending owner / wall type and separates")
    print("  that collision without using the full residue word.")
    print("  theorem target: same-residue or same-v2 height flex must cross a")
    print("  Newton/secondary-fan wall with positive off-grid floor, land on the")
    print("  AP/GW 12->24 hinge, or emit height-discriminant debt.")
    print()

    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof carriers and hidden sidecars, not runners/arcs")
    print(
        "  pairwise_observable=weighted retained proof obligations minus destroyed "
        "endpoint/height/off-grid coordinates"
    )
    print("  switch/gauge=higher proof-carrier score; ties use fewer destroyed scarce fields")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={directed_triangles(names, edges)}")
    print(f"  hamiltonian_path_count={hamiltonian_paths(names, edges)}")
    print("  priority_path=" + " -> ".join(names))
    print()

    print("CONCLUSION")
    print(
        "  The next two non-repeating proof angles are endpoint-owner boundary "
        "currents and tropical height/discriminant walls.  They attack the two "
        "places HYP-3311 leaves open: endpoint-owner loss and same-residue or "
        "same-v2 height flex outside the curated bank.  Residue words remain a "
        "good sidecar, but HYP-3402 demotes residue-only reuse because it is "
        "already bank-local evidence, not the global theorem."
    )


if __name__ == "__main__":
    main()
