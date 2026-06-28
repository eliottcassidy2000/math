#!/usr/bin/env python3
"""Creative proof-reframe lead atlas for LRC14.

This script is a research router, not a proof.  It pulls one hard fact from
the current actual-packet sheaf instantiation:

    the HYP-3301 coarse packet has one mixed theorem-exit fiber on the curated
    bank, and the HYP-3310 nonunit residue word kills it while v2 alone does not.

It then generates and ranks proof leads that try to turn that local repair into
global rigor.  Tournament Analysis is included over proof leads rather than
runners or arcs.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Callable
import itertools
import sys


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relpath: str):
    path = ROOT / relpath
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


actual_packet = load_module(
    "creative_actual_packet",
    "04-computation/lrc14_actual_packet_sheaf_instantiation_codex_20260628.py",
)


@dataclass(frozen=True)
class Lead:
    lead_id: str
    title: str
    core_move: str
    preserved: tuple[str, ...]
    destroyed_if_naive: tuple[str, ...]
    first_test: str
    theorem_gate: str
    risk: str
    tags: tuple[str, ...]
    scores: dict[str, int]

    @property
    def total(self) -> int:
        weights = {
            "actual_packet_pressure": 4,
            "finite_check": 3,
            "globalization_path": 3,
            "sidecar_discipline": 2,
            "covering_layer": 2,
            "novelty": 1,
            "formalizable": 2,
            "risk_penalty": -2,
        }
        return sum(weights[key] * value for key, value in self.scores.items())


def actual_packet_facts() -> dict[str, object]:
    rows = actual_packet.build_rows()
    coarse_mixed = actual_packet.mixed_fibers(rows, lambda row: row.coarse_sheaf_base)
    residue_mixed = actual_packet.mixed_fibers(
        rows, lambda row: row.coarse_sheaf_base + (row.nonunit_residue_word,)
    )
    v2_mixed = actual_packet.mixed_fibers(
        rows, lambda row: row.coarse_sheaf_base + (row.nonunit_v2_word,)
    )
    qgt14 = [row for row in rows if row.q_threshold > 14]
    return {
        "row_count": len(rows),
        "coarse_mixed_count": len(coarse_mixed),
        "coarse_mixed_size": len(coarse_mixed[0]) if coarse_mixed else 0,
        "residue_mixed_count": len(residue_mixed),
        "v2_mixed_count": len(v2_mixed),
        "qgt14_count": len(qgt14),
        "qgt14_kernel_flags": sorted({row.kernel_flag for row in qgt14}),
        "coarse_mixed_names": [row.name for row in coarse_mixed[0]] if coarse_mixed else [],
    }


def lead_catalog() -> list[Lead]:
    return [
        Lead(
            "R01",
            "Residue-word breakpoint theorem",
            "Enlarge the actual-packet bank until nonunit residue-word exactness first fails; prove no failure, or let the first failure name the missing sidecar.",
            ("coarse sheaf base", "nonunit residue word", "kernel exit"),
            ("height", "endpoint owner", "off-grid floor"),
            "Run the HYP-3301 actual-packet test over progressively larger HYP-2963 residual samples and log the first mixed residue fiber.",
            "A finite theorem: residue-word fibers are route-pure unless a named height/owner/floor obstruction appears.",
            "Could only move the problem to a larger finite bank unless the first-failure class has a symbolic parametrization.",
            ("actual-packet", "first-failure", "covering", "sheaf"),
            {
                "actual_packet_pressure": 3,
                "finite_check": 3,
                "globalization_path": 3,
                "sidecar_discipline": 3,
                "covering_layer": 3,
                "novelty": 2,
                "formalizable": 3,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R02",
            "Denominator-curvature transport",
            "Treat blocker ledgers under q -> 2q, q -> 7q, and q -> lcm(q,14) as a discrete connection; nonzero transport defect is the witness-producing curvature.",
            ("blocker support", "qdiv", "fiber map", "runner owner"),
            ("scalar blocked/unblocked denominator",),
            "For each qgt14 positive-open packet, compare minimal blocker supports under q -> 2q and q -> 7q against the image of the q support.",
            "Zero curvature across enough ladder maps forces AP/GW/GW-hinge monodromy; positive curvature gives a finite witness.",
            "Needs support-hypergraph data beyond the current curated sheaf rows.",
            ("p-curvature", "denominator", "monodromy", "qdiv"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 2,
                "globalization_path": 3,
                "sidecar_discipline": 3,
                "covering_layer": 2,
                "novelty": 3,
                "formalizable": 2,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R03",
            "Haar-Baire owner-strip zipper",
            "Replace scalar safe mass by a typed Haar product table over endpoint walls and proof clocks; vanishing owner-strip and cross-handoff coefficients should force AP/GW boundary atoms.",
            ("regular-open safe events", "endpoint owner", "Haar coefficient type"),
            ("exact period", "C27/K33 transfer", "height"),
            "Attach typed Haar rectangles to the seven-row mixed fiber and test which coefficient separates petal-named from positive-Haar-open exits.",
            "If no positive open witness exists, a nonzero labelled Haar coefficient or AP/GW boundary H1 must survive.",
            "Haar coefficients can become another analogy unless tied to exact endpoint owners.",
            ("Haar", "Baire", "zipper", "endpoint"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 3,
                "globalization_path": 2,
                "sidecar_discipline": 3,
                "covering_layer": 2,
                "novelty": 3,
                "formalizable": 2,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R04",
            "Rank-one covering-flex Hessian",
            "Freeze the unit C3 skeleton and nonunit residue word, then study the height variables as a constrained flex space; prove the nullspace is one-dimensional and generated by 12->24.",
            ("unit contact graph", "nonunit residue word", "height vector"),
            ("endpoint topology", "off-grid safe mass"),
            "On residue-pure families, compute the linearized boundary-moment Jacobian in nonunit height coordinates and identify first zero modes.",
            "Only the AP/Goddyn-Wong hinge survives the height-flex equations; every other mode opens strict mass or a Phi14d witness.",
            "Most technically demanding; needs symbolic inequalities after finite evidence.",
            ("height-flex", "Jacobian", "covering", "12-24"),
            {
                "actual_packet_pressure": 3,
                "finite_check": 2,
                "globalization_path": 3,
                "sidecar_discipline": 2,
                "covering_layer": 3,
                "novelty": 2,
                "formalizable": 2,
                "risk_penalty": 2,
            },
        ),
        Lead(
            "R05",
            "Shadow-charge Farkas ledger",
            "Turn HYP-3400 into a finite dual certificate: every scalar shadow must preserve, transfer, or debt C3 charge, quadratic charge, and height/flex charge.",
            ("C3 charge", "Qsqrt7 charge", "height/flex charge"),
            ("raw scalar score", "raw tournament rank"),
            "Build a small incidence matrix with rows = quotient attempts and columns = charges; compute minimal hitting sidecars on HYP-3311 plus actual-packet rows.",
            "A no-naked-quotient lemma becomes a Farkas/cover certificate over named charges.",
            "Can become bookkeeping unless one row is tied to a real theorem exit.",
            ("shadow-charge", "Farkas", "controlled-forgetting"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 3,
                "globalization_path": 2,
                "sidecar_discipline": 3,
                "covering_layer": 2,
                "novelty": 2,
                "formalizable": 3,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R06",
            "Taut-wave interval routing",
            "Use the any-angle path-planning analogy literally: nodes are safe intervals with boundary debt; edges are taut owner-preserving shortcuts.",
            ("safe interval", "boundary owner", "route exit"),
            ("raw grid denominator",),
            "For qgt14 rows, compute whether each exact-period boundary has a taut shortcut to positive-Haar-open mass before it can enter zero-open.",
            "Every exact-period cover routes to a positive boundary-moment image, AP/GW boundary, or named K33/H7 debt.",
            "Needs interval endpoint data imported from packet scripts.",
            ("ANYA", "CWave", "interval", "qdiv"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 2,
                "globalization_path": 2,
                "sidecar_discipline": 2,
                "covering_layer": 3,
                "novelty": 3,
                "formalizable": 2,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R07",
            "Petal-positive separator polynomial",
            "The seven-row actual mixed fiber should have a small algebraic separator in nonunit residue counts plus one owner term, distinguishing petal-named from positive-open.",
            ("nonunit residue counts", "petal owner", "kernel exit"),
            ("v2 height", "global route"),
            "Fit all integer linear separators on the seven-row fiber, then test on the larger HYP-2963 labelled bank.",
            "A tiny separator could become the local lemma that discharges the only current coarse mixed fiber.",
            "A fitted separator may fail outside the curated bank.",
            ("mixed-fiber", "separator", "petal", "positive-open"),
            {
                "actual_packet_pressure": 3,
                "finite_check": 3,
                "globalization_path": 1,
                "sidecar_discipline": 2,
                "covering_layer": 2,
                "novelty": 2,
                "formalizable": 2,
                "risk_penalty": 2,
            },
        ),
        Lead(
            "R08",
            "Ramanujan-period projector breakpoint",
            "Project exact-period rows by Ramanujan sums before residue-word comparison; the first failure should distinguish exact-period boundary from residue coincidence.",
            ("exact period", "Ramanujan projector", "nonunit residue word"),
            ("endpoint owner", "height"),
            "For the qgt14 positive-open rows, compare residue exactness before and after exact-period Ramanujan projection.",
            "Residue-word exactness lifts to exact-period exactness, or the projector emits the missing owner/height debt.",
            "Requires careful handling of zero-residue and nonprimitive periods.",
            ("Ramanujan", "period", "projector", "residue"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 2,
                "globalization_path": 2,
                "sidecar_discipline": 3,
                "covering_layer": 2,
                "novelty": 2,
                "formalizable": 2,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R09",
            "Unlabelled-tournament sidecar metagraph",
            "Use A000568-style quotient classes as a stress test: every sidecar deletion or addition is an edge in a metagraph whose SCCs should be proof-obligation classes.",
            ("sidecar set", "mixed-fiber status", "deletion edge"),
            ("raw runner identity",),
            "Build deletion decks for sidecar bundles on actual-packet rows; identify minimal deletions that recreate the unique mixed fiber.",
            "Minimal deletion obstructions become exact first-obstruction cocycles.",
            "Metagraphs can get large unless restricted to current packet fields.",
            ("tournament", "A000568", "deletion", "metagraph"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 3,
                "globalization_path": 1,
                "sidecar_discipline": 3,
                "covering_layer": 1,
                "novelty": 2,
                "formalizable": 2,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R10",
            "Bulk-core charge conservation theorem",
            "Fuse Vitali bulk/core with HYP-3400: positive safe mass cannot disappear at the core unless Phi14/Phi14d witness charge turns on.",
            ("strict safe mass", "Phi witness", "bulk-core side"),
            ("residue word", "height"),
            "Check the actual-packet positive-Haar-open rows for the first Phi14d or finite-trap charge that accounts for their positive mass.",
            "Every positive-open residual descends to a charge-conserving witness or a named finite trap.",
            "Less directly tied to the one mixed fiber than residue/flex approaches.",
            ("Vitali", "bulk-core", "Phi14d", "charge"),
            {
                "actual_packet_pressure": 1,
                "finite_check": 2,
                "globalization_path": 3,
                "sidecar_discipline": 2,
                "covering_layer": 2,
                "novelty": 2,
                "formalizable": 2,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R11",
            "Colored-resonance half-boundary sieve",
            "Fuse the residue-word sheaf with the HYP-2593/HYP-2595 colored CRT discrepancy machinery; a mixed fiber should expose a surviving half-boundary resonance rather than raw component count.",
            ("color class", "residue word", "component boundary", "CRT witness count"),
            ("endpoint owner", "exact denominator period"),
            "On each enlarged residue-word fiber, compute the colored resonance deficit and ask whether the pessimistic component term collapses to an owner-labelled half-boundary term.",
            "A uniform colored deficit bound plus the residue-word floor certifies qdiv>14 positive-open rows without a new zero-open kernel.",
            "Needs a formal bridge from finite color-resonance evidence to a familywise discrepancy theorem.",
            ("colored", "resonance", "CRT", "half-boundary"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 3,
                "globalization_path": 3,
                "sidecar_discipline": 3,
                "covering_layer": 3,
                "novelty": 2,
                "formalizable": 2,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R12",
            "Mayer-gas parity reality reduction",
            "Treat the boundary-moment obstruction as the signed real shadow of an OCF/Mayer gas: support-parity signs must be retained until F7* orbit sums become sign-definite or emit debt.",
            ("support parity", "real Gauss-sum orbit", "relation lattice"),
            ("unsigned count", "QR/NQR magnitude split"),
            "For first-failure candidate rows, group the signed lattice correction by F7* dilation orbit and test sign coherence against the dominant residue weights.",
            "The complex correction reduces to a real signed sum whose nonpositive orbit packets discharge wide residuals, while incoherent orbits become named analytic sidecar debt.",
            "The reality part is exact, but magnitude cancellation is the hard analytic tail.",
            ("Mayer", "OCF", "Gauss", "parity"),
            {
                "actual_packet_pressure": 1,
                "finite_check": 2,
                "globalization_path": 3,
                "sidecar_discipline": 3,
                "covering_layer": 2,
                "novelty": 3,
                "formalizable": 2,
                "risk_penalty": 2,
            },
        ),
        Lead(
            "R13",
            "Analytic-lifting stability ledger",
            "Use the S293 analytic-lifting frame as a guarded sidecar: bulk discrepancy estimates may act only after Krasner-stable 2-adic/7-adic local lifts preserve the packet coordinates.",
            ("bulk discrepancy", "Krasner-stable local lift", "2-adic/7-adic valuation"),
            ("packet owner", "residue-height coupling"),
            "Attach p-adic lift-stability flags to same-residue height moves such as 12->24 and 12->96, then test whether unstable moves are exactly the positive-open ones.",
            "Bulk estimates discharge large-height tails while local stability prevents height/flex roots from sliding into a hidden AP/GW-like zero-open kernel.",
            "The analytic inputs must be stated as precise inequalities before this can be more than a ledger.",
            ("BDH", "Krasner", "p-adic", "lifting"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 1,
                "globalization_path": 3,
                "sidecar_discipline": 3,
                "covering_layer": 3,
                "novelty": 2,
                "formalizable": 1,
                "risk_penalty": 2,
            },
        ),
        Lead(
            "R14",
            "Endpoint deletion-cut Menger theorem",
            "Turn sidecar deletion decks into a min-cut problem: a quotient may forget endpoint owners only if every owner-deletion cut still has a route to the same theorem exit.",
            ("endpoint owner", "sidecar deletion deck", "route cut"),
            ("raw unlabelled tournament class", "fiber multiplicity"),
            "Build the deletion deck of actual-packet fields and compute minimal cuts that recreate the seven-row mixed fiber.",
            "Minimal mixed-fiber cuts are exactly first-obstruction cocycles; if every cut is owner-labelled, endpoint resurrection closes the quotient.",
            "May duplicate the metagraph route unless it proves a real cut/flow duality.",
            ("Menger", "deletion", "endpoint", "cocycle"),
            {
                "actual_packet_pressure": 3,
                "finite_check": 3,
                "globalization_path": 2,
                "sidecar_discipline": 3,
                "covering_layer": 2,
                "novelty": 2,
                "formalizable": 2,
                "risk_penalty": 1,
            },
        ),
        Lead(
            "R15",
            "Worpitzky-Faulhaber finite-difference ladder",
            "Read height/flex motion through odd Faulhaber moments and Worpitzky finite differences; the legal flex should have vanishing low odd moments until the 12->24 hinge appears.",
            ("odd moment", "finite difference", "height-flex variable"),
            ("endpoint owner", "color resonance"),
            "For residue-pure covering families, compute low-order odd boundary-moment differences along height moves and search for a sign-regular ladder.",
            "Moment sign regularity makes the covering-flex Hessian rank-one theorem symbolic instead of purely computational.",
            "Needs careful normalization so moment positivity is not another scalar shadow.",
            ("Worpitzky", "Faulhaber", "moment", "height-flex"),
            {
                "actual_packet_pressure": 2,
                "finite_check": 2,
                "globalization_path": 3,
                "sidecar_discipline": 2,
                "covering_layer": 3,
                "novelty": 3,
                "formalizable": 2,
                "risk_penalty": 2,
            },
        ),
    ]


def procedural_recombinations(leads: list[Lead]) -> list[tuple[str, str, str, int]]:
    objects = ["residue_word", "height_flex", "endpoint_owner", "boundary_period"]
    operators = ["first_failure", "curvature_transport", "Haar_zipper", "charge_dual"]
    certificates = ["finite_mixed_fiber", "positive_open_mass", "AP_GW_boundary", "named_debt"]
    known_tags = {tag for lead in leads for tag in lead.tags}
    out = []
    for obj, op, cert in itertools.product(objects, operators, certificates):
        score = 0
        if obj in ("residue_word", "height_flex"):
            score += 2
        if op in ("first_failure", "curvature_transport"):
            score += 2
        if cert in ("finite_mixed_fiber", "positive_open_mass"):
            score += 2
        if any(part in known_tags for part in (obj, op, cert)):
            score += 1
        out.append((obj, op, cert, score))
    return sorted(out, key=lambda item: (-item[3], item[:3]))[:12]


def tournament(leads: list[Lead]) -> dict[str, object]:
    ordered = sorted(leads, key=lambda lead: (-lead.total, lead.lead_id))
    scores = {lead.lead_id: lead.total for lead in leads}
    score_hist: dict[int, int] = {}
    for score in scores.values():
        score_hist[score] = score_hist.get(score, 0) + 1
    return {
        "vertices": len(leads),
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": 0,
        "hamiltonian_path_count": 1,
        "priority_path": [lead.lead_id for lead in ordered],
    }


def print_report() -> None:
    facts = actual_packet_facts()
    leads = lead_catalog()
    tour = tournament(leads)
    ordered = sorted(leads, key=lambda lead: (-lead.total, lead.lead_id))

    print("HYP-3404 creative LRC14 proof-reframe lead atlas")
    print("=" * 78)
    print("actual_packet_anchor:")
    print(f"  rows={facts['row_count']}")
    print(
        "  coarse_mixed_fibers="
        f"{facts['coarse_mixed_count']} size={facts['coarse_mixed_size']}"
    )
    print(f"  residue_word_mixed_fibers={facts['residue_mixed_count']}")
    print(f"  v2_word_mixed_fibers={facts['v2_mixed_count']}")
    print(
        "  qgt14="
        f"{facts['qgt14_count']} kernel_flags={facts['qgt14_kernel_flags']}"
    )
    print("  unique_mixed_fiber_names=")
    for name in facts["coarse_mixed_names"]:
        print(f"    - {name}")
    print()

    print("RANKED LEADS")
    for lead in ordered:
        print(f"{lead.lead_id} score={lead.total:2d} :: {lead.title}")
        print(f"  move: {lead.core_move}")
        print(f"  first_test: {lead.first_test}")
        print(f"  theorem_gate: {lead.theorem_gate}")
        print(f"  risk: {lead.risk}")
        print()

    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof-reframe leads, not runners/arcs")
    print(
        "  pairwise_observable=weighted proof leverage: actual-packet pressure, "
        "finite checkability, globalization path, sidecar discipline, covering-layer relevance"
    )
    print("  switch/gauge=A -> B iff A has larger weighted score; ties use lead_id Hamiltonian path")
    print(f"  vertices_count={tour['vertices']}")
    print(f"  score_hist={tour['score_hist']}")
    print(f"  directed_3cycles={tour['directed_3cycles']}")
    print(f"  hamiltonian_path_count={tour['hamiltonian_path_count']}")
    print("  priority_path=" + " -> ".join(tour["priority_path"]))
    print()

    print("PROCEDURAL RECOMBINATION BUCKET")
    for obj, op, cert, score in procedural_recombinations(leads):
        print(f"  {obj:16s} x {op:20s} -> {cert:18s} score={score}")
    print()

    print("ASSUMPTION CHALLENGE")
    print("  Alternate vertices considered: packet fibers, sidecar bundles, quotient maps,")
    print("  denominator transports, endpoint events, Haar rectangles, charge reservoirs,")
    print("  deletion decks, and theorem exits.")
    print("  Preserved predicate: whether a primitive packet exits by q-witness, AP/GW")
    print("  boundary, petal/K33 named debt, positive-open covering mass, or a new")
    print("  zero-open kernel.")
    print("  Destroyed if naively quotiented: height/flex, endpoint owner, exact period,")
    print("  off-grid floor, and the transverse Qsqrt(-7) sidecar.")
    print()

    print("BOTTOM LINE")
    print("  The strongest new route is not another scalar.  It is a first-failure")
    print("  theorem: enlarge the actual-packet sheaf until residue-word exactness fails;")
    print("  prove no such failure, or show the first failure is exactly height/flex,")
    print("  owner-current, tropical off-grid floor, exact-period, or named K33/H7")
    print("  debt.  HYP-3402 and HYP-3403 make owner-current/tropical-wall plus")
    print("  shadow-charge packet gluing the first handoffs before harder Hessian")
    print("  or curvature routes.  This would convert")
    print("  the current creative stack into a rigorous labelled-packet theorem target.")


if __name__ == "__main__":
    print_report()
