#!/usr/bin/env python3
"""HYP-3266: formal/analytic proof-obligation ledger for LRC(14).

The goal is not to discover another scalar invariant.  It is to say exactly
which current proof obligations are Lean-closed, finite-imported,
conditionally wired, analytically open, or refuted as a route.

Tournament Analysis:
  vertices: proof obligations, not runners or arcs.
  pairwise observable: which obligation unlocks more downstream proof mass
    while retaining a formal statement boundary.
  switch/gauge: A -> B if A has higher current proof priority, with ties broken
    by lower formalization distance and then by dependency count.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable


class Status(Enum):
    CLOSED_LEAN = "CLOSED_LEAN"
    CLOSED_FINITE = "CLOSED_FINITE"
    FINITE_IMPORT = "FINITE_IMPORT"
    CONDITIONAL_GLUE = "CONDITIONAL_GLUE"
    EVIDENCE_ONLY = "EVIDENCE_ONLY"
    OPEN_ANALYTIC = "OPEN_ANALYTIC"
    OPEN_RIGIDITY = "OPEN_RIGIDITY"
    FALSE_ROUTE = "FALSE_ROUTE"


STATUS_WEIGHT = {
    Status.OPEN_RIGIDITY: 9,
    Status.OPEN_ANALYTIC: 8,
    Status.EVIDENCE_ONLY: 6,
    Status.CONDITIONAL_GLUE: 4,
    Status.FINITE_IMPORT: 3,
    Status.CLOSED_FINITE: 2,
    Status.CLOSED_LEAN: 1,
    Status.FALSE_ROUTE: 0,
}


@dataclass(frozen=True)
class Obligation:
    oid: str
    title: str
    status: Status
    layer: str
    formal_boundary: str
    analytic_boundary: str
    sources: tuple[str, ...]
    dependencies: tuple[str, ...]
    next_formal: str
    next_analytic: str
    impact: int
    formal_distance: int

    @property
    def priority(self) -> tuple[int, int, int, str]:
        return (
            self.impact * 10 + STATUS_WEIGHT[self.status],
            -self.formal_distance,
            len(self.dependencies),
            self.oid,
        )


OBLIGATIONS: tuple[Obligation, ...] = (
    Obligation(
        "O00",
        "Concrete Mreach compactness readout",
        Status.CLOSED_LEAN,
        "terminal formal glue",
        "LRCMreachConcrete.lonely_of_Mreach_ge; LRCFourteenSkeleton.lonely_of_Mreach_ge",
        "Extreme-value/continuity handoff is already formalized.",
        ("LRCMreachConcrete.lean", "LRCFourteenSkeleton.lean", "T967"),
        (),
        "None except keeping imports green.",
        "None.",
        7,
        0,
    ),
    Obligation(
        "O01",
        "Denominator-sieve saturation gates",
        Status.CLOSED_LEAN,
        "early-gate formal glue",
        "lrc14_no_multiple_of_14; lrc14_all_odd; lrc14_counterexample_saturated",
        "Counterexamples must be q-divisible for every q=2..14.",
        ("LRCFourteenSkeleton.lean", "LonelyRunner.sieve_one_div"),
        (),
        "No action; use as first classifier in packet ledgers.",
        "No action.",
        6,
        0,
    ),
    Obligation(
        "O02",
        "Tight-locus construction by unit witnesses",
        Status.CLOSED_FINITE,
        "finite arithmetic, not yet Lean-native",
        "No file-backed Lean theorem yet; S81 script verifies unit witnesses.",
        "AP/GW and d*AP are safe at all units of Z/(14d).",
        ("HYP-3253", "lrc_tighten_rigor_macmini_S81.py", "kps-S255"),
        (),
        "Add a Lean lemma for unit a: s*a == +/-1 mod 14d gives distance 1/(14d).",
        "None for construction; it is elementary.",
        6,
        2,
    ),
    Obligation(
        "O03",
        "Bounded AP/GW single-swap margin",
        Status.CLOSED_FINITE,
        "finite exact check",
        "No Lean theorem; S81 finite check reports AP/GW only and gap about 0.00262.",
        "Single-swaps r<=26 have strict margin outside AP/GW.",
        ("HYP-3253", "HYP-3250", "HYP-2915"),
        ("O02",),
        "Turn the finite check into a certificate table with rational breakpoints.",
        "Extend to exhaustive two-swap/multi-swap rigidity or route to O15.",
        5,
        3,
    ),
    Obligation(
        "O04",
        "Cap RHS / pair-Pascal ledger",
        Status.CLOSED_LEAN,
        "exact rational Lean ledger",
        "capRat_eq_pairPascalMassRat_dense; capLedgerK8/K9/K10",
        "Cap RHS values are solved; k=8/9 carry explicit higher-order debt.",
        ("LRCProofFrontier.lean", "THM-576", "THM-577"),
        (),
        "Use as imported constants in hp0cap/witness-floor proofs.",
        "None.",
        6,
        0,
    ),
    Obligation(
        "O05",
        "Small-k cover-set pigeonhole and p0 monotonicity",
        Status.CLOSED_LEAN,
        "elementary cover-set Lean proof",
        "LRCCoverBound.coverSet_mono; slowmu_coverSet_eq_zero_of_card_lt_six",
        "k<6 has p0=0; coverSet is monotone.",
        ("LRCCoverBound.lean",),
        (),
        "No action except linking into packet classifiers.",
        "No action.",
        5,
        0,
    ),
    Obligation(
        "O06",
        "Concrete witness floor from p0 wide bound",
        Status.CLOSED_LEAN,
        "event-measure Lean proof",
        "LRCWitnessFloorConcrete.witness_pos_from_wide_bound_margin and goodSet variants",
        "Witness-floor positivity/margin reduces to p0 <= cap with slack.",
        ("LRCWitnessFloorConcrete.lean", "LRCWitnessPartA.lean", "HYP-2832"),
        ("O04",),
        "Keep the concrete event definitions wired to witnessG2.",
        "This transfers the burden to O10/O12, not a separate analytic gap.",
        8,
        0,
    ),
    Obligation(
        "O07",
        "gK8 per-shape Delsarte and arithmetic imports",
        Status.CLOSED_LEAN,
        "Lean arithmetic + Delsarte feasibility",
        "gK8_per_shape_bound; gK8_dual_feasible; gK8_singlefar_arithmetic_checksum",
        "10*p0 <= L_yK8 and binding arithmetic checks are closed.",
        ("LRCFourteenSkeleton.lean", "LRCGk8SingleFar.lean", "HYP-2829"),
        (),
        "Use as lower half of gK8 route.",
        "Remaining scalar concentration is O13.",
        7,
        0,
    ),
    Obligation(
        "O08",
        "Contact-holonomy quotient-curvature repair",
        Status.CLOSED_FINITE,
        "exact bounded k=8 computation",
        "HYP-3267 file-backed script: zeta7 contact holonomy kills 62 mixed fibers.",
        "Local lag/residue curvature is repaired by zeta_7 endpoint holonomy.",
        ("HYP-3267", "lrc14_contact_holonomy_curvature_codex_20260628.py"),
        ("O04",),
        "Define contact support and first zeta moment in Lean or certificate JSON.",
        "Test if the same sidecar improves O10/O12 floor packets.",
        5,
        3,
    ),
    Obligation(
        "O09",
        "Residual packet to finite-address / observer-gluing certificates",
        Status.CONDITIONAL_GLUE,
        "conditional Lean packet theorem",
        "LRCProofFrontier.ResidualToFiniteAddressPackets; lrc14_from_residual_packet_frontier",
        "If the residual classifier emits packets, existing glue reaches LRC14.",
        ("LRCProofFrontier.lean", "LRCFiniteAddressBranchClosure.lean"),
        ("O01", "O05"),
        "Instantiate isResidual with the current S81/S256 classifier.",
        "Analytic/finiteness exits still come from O10/O12/O15.",
        8,
        1,
    ),
    Obligation(
        "O10",
        "Wide cover bound hp0cap for binding k=8..12",
        Status.OPEN_ANALYTIC,
        "Lean reduction isolates analytic residual",
        "LRCCoverBound.slowmu_coverSet_lt_cap_of_decorrelation",
        "Need p0 <= p0_decorr plus finite Q<cap for binding families.",
        ("LRCCoverBound.lean", "HYP-3136", "HYP-3129", "HYP-3132"),
        ("O04", "O05", "O06"),
        "Create explicit Prop for p0_decorr and exact finite Q tables.",
        "Prove resonant/Tornheim decorrelation, possibly after contact-holonomy/Qsqrt(-7) lift.",
        10,
        3,
    ),
    Obligation(
        "O11",
        "Corrected witnessG2/rhoGlob floor",
        Status.OPEN_ANALYTIC,
        "Lean arithmetic ledger is closed; measure floor is open",
        "witness_large_floor_from_rhoGlob_lower_bound; need rhoGlobFloorRat k <= rhoGlob(shape)",
        "For k=8..13, define rhoGlob=witnessG2 concretely and prove the floor.",
        ("LRCFourteenSkeleton.lean", "lrc14_global_threshold_ladder_codex.out"),
        ("O06",),
        "Replace opaque witnessG2/rhoGlob by the concrete goodSet measure where possible.",
        "Prove lower bounds for the actual shape classes, not just the rational ledger.",
        9,
        2,
    ),
    Obligation(
        "O12",
        "Part A / off-grid bulk survivor positivity",
        Status.OPEN_ANALYTIC,
        "finite-Vmax Lean algebra exists; event-level implication open",
        "thm527_partA_density_pos_implies_reach; LRCWitnessPartA finite arc-budget lemmas",
        "Need G2>0 -> Mreach>=1/14, including resonant v via off-grid bulk positivity.",
        ("LRCWitnessPartA.lean", "HYP-3255", "HYP-3253", "HYP-3252", "S52"),
        ("O06", "O11"),
        "State the exact finite-ruler GOOD arcCount bound used by S81.",
        "Prove off-grid bulk survivor positivity; the 14-grid core is measure-zero and not the generic survivor.",
        10,
        2,
    ),
    Obligation(
        "O13",
        "gK8 concentration extremality",
        Status.OPEN_ANALYTIC,
        "Prop isolated; arithmetic subchecks closed",
        "gK8_concentration_extremality(LyVal,k,boundedMax)",
        "Need max_E L_yK8 <= 10*cap_k after far-count split.",
        ("LRCFourteenSkeleton.lean", "HYP-2829", "HYP-3132", "HYP-3142"),
        ("O07", "O10"),
        "Formalize r=0 finite, r=1 periodmax, r>=2 margin as separate Props.",
        "Prove single-far periodic sup and multi-far decorrelation margin.",
        8,
        2,
    ),
    Obligation(
        "O14",
        "Doublet R-tail uniform bound",
        Status.OPEN_ANALYTIC,
        "finite ledger closed; analytic tail open",
        "doublet_Rtail_uniform_bound(Rtail)",
        "Need Mordell-Tornheim/Koksma bound for |R(M)|.",
        ("LRCFourteenSkeleton.lean", "LRCDoubletWitnessFloor.lean", "THM-564"),
        ("O07",),
        "Define Rtail from the Python finite scout in Lean-friendly form.",
        "Prove sup_M>=15 |Rtail(M)| <= 12*zeta(3)*N/pi^3.",
        7,
        3,
    ),
    Obligation(
        "O15",
        "Full tight-locus rigidity via finite equioscillation sidecars",
        Status.OPEN_RIGIDITY,
        "necessary conditions known; sufficiency open",
        "No Lean Prop yet; should become finite equioscillation + blind sidecars => AP/GW/dilation.",
        "Bounded census is clean; HYP-3257/HYP-3258/HYP-3259/HYP-3265 reduce the full problem to unit-blind sidecars, binding/covering rigidity, manifold flex, and contact-graph chambers.",
        ("HYP-3265", "HYP-3259", "HYP-3258", "HYP-3257", "HYP-3255", "HYP-3253", "HYP-3250", "kps-S255", "HYP-2914", "HYP-2915"),
        ("O02", "O03"),
        "State a finite normalized rigidity theorem with explicit normalizer and sidecar fields.",
        "Prove the finite equioscillation/contact-graph system plus blind residue/height ledger has only AP/GW/dilations, or name the residual family.",
        10,
        4,
    ),
    Obligation(
        "O16",
        "Qsqrt(-7) signed-floor reorganization",
        Status.EVIDENCE_ONLY,
        "no formal statement yet",
        "No Lean boundary; current status is a suggested basis change.",
        "May turn signed-SPEC cancellation into positivity in the natural field.",
        ("HYP-3267", "HYP-3255", "HYP-3254", "HYP-3252", "the-index-theorem-describes-the-floor-proves.md", "THM-250"),
        ("O10", "O12"),
        "Define the signed floor coefficients in a Qsqrt(-7) basis.",
        "Test whether O10/O12 residual terms become sign-definite after the basis change.",
        6,
        5,
    ),
    Obligation(
        "O17",
        "Obsolete rhoStar 2/7 route",
        Status.FALSE_ROUTE,
        "kept only as conditional glue",
        "obsoleteRhoStarUniformFloor and lrc14_from_obsolete_rhoStar_floor",
        "Uniform 2/7 via-Vmax floor is refuted for intended object.",
        ("LRCFourteenSkeleton.lean", "HYP-3252"),
        (),
        "Do not spend proof effort here except preserving the warning.",
        "No action.",
        1,
        0,
    ),
)


def status_counts(obligations: Iterable[Obligation]) -> dict[Status, int]:
    counts = {status: 0 for status in Status}
    for obligation in obligations:
        counts[obligation.status] += 1
    return counts


def md_row(cols: Iterable[object]) -> str:
    return "| " + " | ".join(str(c).replace("\n", " ") for c in cols) + " |"


def main() -> None:
    obligations = list(OBLIGATIONS)
    by_id = {ob.oid: ob for ob in obligations}
    counts = status_counts(obligations)
    ranked = sorted(obligations, key=lambda ob: ob.priority, reverse=True)

    print("HYP-3266 formal/analytic LRC(14) proof-obligation ledger")
    print("=" * 78)
    print(f"obligations={len(obligations)}")
    print(
        "status_counts="
        + ", ".join(f"{status.value}:{counts[status]}" for status in Status)
    )
    print()

    print("SUMMARY: WHERE WE STAND")
    print("- closed Lean glue/kernels:", counts[Status.CLOSED_LEAN])
    print("- exact finite closures/imports:", counts[Status.CLOSED_FINITE] + counts[Status.FINITE_IMPORT])
    print("- conditional Lean glue:", counts[Status.CONDITIONAL_GLUE])
    print("- evidence-only proposed bridges:", counts[Status.EVIDENCE_ONLY])
    print("- open analytic obligations:", counts[Status.OPEN_ANALYTIC])
    print("- open rigidity obligations:", counts[Status.OPEN_RIGIDITY])
    print("- refuted/obsolete routes retained as warnings:", counts[Status.FALSE_ROUTE])
    print()

    print("PRIORITY HAMILTONIAN PATH")
    print("  " + " -> ".join(ob.oid for ob in ranked))
    print()

    print("OBLIGATION TABLE")
    print(md_row(["id", "status", "layer", "title", "formal boundary", "next analytic step"]))
    print(md_row(["---", "---", "---", "---", "---", "---"]))
    for ob in obligations:
        print(
            md_row(
                [
                    ob.oid,
                    ob.status.value,
                    ob.layer,
                    ob.title,
                    ob.formal_boundary,
                    ob.next_analytic,
                ]
            )
        )
    print()

    print("OPEN CORE")
    for ob in ranked:
        if ob.status in {Status.OPEN_ANALYTIC, Status.OPEN_RIGIDITY, Status.EVIDENCE_ONLY}:
            deps = ", ".join(ob.dependencies) if ob.dependencies else "none"
            print(f"{ob.oid} {ob.title}")
            print(f"  status: {ob.status.value}")
            print(f"  depends_on: {deps}")
            print(f"  formal: {ob.next_formal}")
            print(f"  analytic: {ob.next_analytic}")
    print()

    print("LEAN-FACING NAMES TO KEEP GREEN")
    lean_names = [
        "LRCMreachConcrete.lonely_of_Mreach_ge",
        "LRCFourteenSkeleton.lrc14_from_witness_floor_cases_given_nodes",
        "LRCWitnessFloorConcrete.goodSet_witness_margin_from_wide_bound",
        "LRCCoverBound.slowmu_coverSet_lt_cap_of_decorrelation",
        "LRCWitnessPartA.finite_witness_pos_from_goodSet_margin_uniform_arc_bound_shapes",
        "LRCProofFrontier.lrc14_from_residual_packet_frontier",
        "LRCFourteenSkeleton.gK8_concentration_extremality",
        "LRCFourteenSkeleton.doublet_Rtail_uniform_bound",
    ]
    for name in lean_names:
        print(f"  - {name}")
    print()

    print("TOURNAMENT ANALYSIS")
    n = len(ranked)
    edge_count = n * (n - 1) // 2
    print("vertices=proof_obligations")
    print("pairwise_observable=downstream proof mass unlocked with named formal boundary")
    print("switch_gauge=A beats B iff priority(A)>priority(B)")
    print("score_hist={" + ", ".join(f"{i}:1" for i in range(n)) + "}")
    print("directed_3cycles=0")
    print("scc_sizes=[" + ",".join("1" for _ in ranked) + "]")
    print(f"edge_count={edge_count}")
    print("edge_flips_against_raw_runner_gauge=not_applicable")
    print(f"hamiltonian_path_count=1")
    print("tie_hamiltonian_path=" + " -> ".join(ob.oid for ob in ranked))
    print()

    print("ASSUMPTION CHALLENGE")
    print(
        "alternate_vertices_considered=runners,gaps,fixed_sections,section_boundaries,"
        "wall_crossings,residues,cover_arcs,Fourier_modes,matroid_cells,proof_obligations"
    )
    print("chosen_vertices=proof_obligations")
    print("preserved_predicate=conditional route to Mreach>=1/14 and LRC14")
    print(
        "destroyed_information=row geometry, endpoint-owner cells, and exact residual family labels; "
        "these must re-enter through O08/O09/O10/O12/O15 sidecars"
    )
    print()

    print("DEPENDENCY EDGES")
    for ob in obligations:
        for dep in ob.dependencies:
            if dep not in by_id:
                raise SystemExit(f"unknown dependency {dep} for {ob.oid}")
            print(f"  {dep} -> {ob.oid}")


if __name__ == "__main__":
    main()
