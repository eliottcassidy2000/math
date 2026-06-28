#!/usr/bin/env python3
"""HYP-3440 scout: Bring/SC/Menger motifs and recursive Charal signatures.

This is a proof-angle router, not an LRC14 proof.  It starts from the current
post-HYP-3406 state:

* nonunit residue is only a curated-bank separator;
* residue plus height still leaves one enlarged-bank mixed theorem-exit fiber;
* residue plus endpoint-owner support separates the scanned enlarged banks.

The requested outside motifs are treated as typed sidecars.  They are useful
only when they name the LRC coordinate they preserve and the coordinate they
would otherwise destroy.

This HYP-3440 scout executes the requested Bring/SC/Menger slice after
HYP-3412 measures the special-function sidecars and HYP-3414 turns the owner
cut into a finite clause/transversal calculus.  After rebasing over the
HYP-3415/HYP-3418 critical-path correction, HYP-3424's covering-floor duality
transfer, the paired HYP-3425 Helly/energy-sheet scouts, HYP-3426's
one-branch mirror endpoint-support audit, HYP-3427's wall-signature atlas, and
HYP-3428's two-adic descent loss ledger, plus HYP-3429's component-spine
endpoint certificate, it treats the endpoint-cut geometry as a
sidecar/certificate router, not as the completion theorem.  After rebasing over
HYP-3430-HYP-3437, this router is even narrower: its Menger language should feed
the branch-cover and overlap-tax certificate stack rather than compete with it.
The LRC14 completion route is the covering decorrelation floor, now sharpened
to a 2-adic even-speed descent problem.  After the HYP-3423 guardrail, topology
and boundary language are allowed to certify q-uniform residue/order data only;
q-specific magnitude or floor claims must restore arithmetic, owner-current, or
two-adic sidecars.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


KEY_FILES = [
    "00-navigation/LRC-LENS-MAP.md",
    "00-navigation/LRC-TECHNIQUE-INDEX.md",
    "00-navigation/LRC-TOURNAMENT-TECHNIQUE-INDEX.md",
    "00-navigation/OPEN-QUESTIONS.md",
    "00-navigation/SESSION-LOG.md",
    "05-knowledge/hypotheses/INDEX.md",
    "05-knowledge/hypotheses/HYP-3404-lrc14-creative-reframe-lead-atlas.md",
    "05-knowledge/hypotheses/HYP-3405-lrc14-ap-collar-finite-lemma-certificate.md",
    "05-knowledge/hypotheses/HYP-3406-lrc14-expanded-residue-owner-repair.md",
    "05-knowledge/hypotheses/HYP-3407-lrc14-boundary-uniformization-cut-stability.md",
    "05-knowledge/hypotheses/HYP-3408-lrc14-exotic-guardrail-reframe-atlas.md",
    "05-knowledge/hypotheses/HYP-3409-lrc14-recursive-sidecar-pattern-atlas.md",
    "05-knowledge/hypotheses/HYP-3410-lrc14-bring-schwarz-bdh-menger-charal-recursion.md",
    "05-knowledge/hypotheses/HYP-3412-lrc14-special-function-cut-signature-recursion.md",
    "05-knowledge/hypotheses/HYP-3414-lrc14-owner-cut-resurrection-calculus.md",
    "05-knowledge/hypotheses/HYP-3416-lrc14-recursive-quotient-ladder.md",
    "05-knowledge/hypotheses/HYP-3417-lrc14-owner-cut-dual-certificate-synthesis.md",
    "05-knowledge/hypotheses/HYP-3419-lrc14-charal-owner-cut-recursion-prototype.md",
    "05-knowledge/hypotheses/HYP-3420-lrc14-owner-cut-chiral-transcendence-synthesis.md",
    "05-knowledge/hypotheses/HYP-3421-lrc14-offgrid-resonance-transparency-rprime-closure.md",
    "05-knowledge/hypotheses/HYP-3422-lrc14-two-adic-offgrid-relocation.md",
    "05-knowledge/hypotheses/HYP-3423-lrc14-quniform-topology-arithmetic-break-guardrail.md",
    "05-knowledge/hypotheses/HYP-3424-lrc14-covering-floor-duality-transfer.md",
    "05-knowledge/hypotheses/HYP-3425-lrc14-two-branch-obstruction-helly.md",
    "05-knowledge/hypotheses/HYP-3425-lrc14-additive-energy-sheet-sidecar.md",
    "05-knowledge/hypotheses/HYP-3426-lrc14-one-branch-mirror-endpoint-support.md",
    "05-knowledge/hypotheses/HYP-3427-lrc14-two-branch-wall-signature-atlas.md",
    "05-knowledge/hypotheses/HYP-3428-lrc14-two-adic-descent-loss-ledger.md",
    "05-knowledge/hypotheses/HYP-3429-lrc14-component-spine-certificate.md",
    "05-knowledge/hypotheses/HYP-3430-lrc14-euler-mascheroni-harmonic-intercept.md",
    "05-knowledge/hypotheses/HYP-3431-lrc14-canonical-corridor-fence.md",
    "05-knowledge/hypotheses/HYP-3432-lrc14-euler-mascheroni-wall-budget.md",
    "05-knowledge/hypotheses/HYP-3433-lrc14-euler-mascheroni-endpoint-spine-finite-part.md",
    "05-knowledge/hypotheses/HYP-3434-lrc14-gamma-harmonic-sieve-remainder.md",
    "05-knowledge/hypotheses/HYP-3435-lrc14-two-adic-branch-cover-certificate.md",
    "05-knowledge/hypotheses/HYP-3436-lrc14-minimal-bad-core-cover-extractor.md",
    "05-knowledge/hypotheses/HYP-3437-lrc14-overlap-menger-cut-certificate.md",
    "05-knowledge/results/lrc14_expanded_residue_owner_repair_codex_20260628.out",
    "05-knowledge/results/lrc14_special_function_cut_signature_recursion_codex_20260628.out",
    "05-knowledge/results/lrc14_owner_cut_resurrection_calculus_codex_20260628.out",
    "05-knowledge/results/lrc14_boundary_uniformization_cut_stability_codex_20260628.out",
    "05-knowledge/results/lrc14_exotic_guardrail_reframe_atlas_codex_20260628.out",
    "05-knowledge/results/lrc14_recursive_sidecar_pattern_atlas_codex_20260628.out",
    "05-knowledge/results/lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628.out",
    "05-knowledge/results/lrc14_owner_cut_dual_certificate_synthesis_codex_20260628.out",
    "05-knowledge/results/lrc14_charal_owner_cut_recursion_prototype_codex_20260628.out",
    "05-knowledge/results/lrc14_transcendence_cut_chiral_synthesis_codex_20260628.out",
    "05-knowledge/results/lrc14_offgrid_resonance_transparency_rprime_closure_codex_20260628.out",
    "05-knowledge/results/lrc14_two_adic_offgrid_relocation_codex_20260628.out",
    "05-knowledge/results/lrc14_quniform_topology_arithmetic_break_codex_20260628.out",
    "05-knowledge/results/lrc14_covering_floor_duality_transfer_codex_20260628.out",
    "05-knowledge/results/lrc14_two_branch_obstruction_helly_codex_20260628.out",
    "05-knowledge/results/lrc14_additive_energy_sheet_sidecar_codex_20260628.out",
    "05-knowledge/results/lrc14_one_branch_mirror_endpoint_support_codex_20260628.out",
    "05-knowledge/results/lrc14_two_branch_wall_signature_atlas_codex_20260628.out",
    "05-knowledge/results/lrc14_two_adic_descent_loss_ledger_codex_20260628.out",
    "05-knowledge/results/lrc14_component_spine_certificate_codex_20260628.out",
    "05-knowledge/results/lrc14_euler_mascheroni_harmonic_intercept_codex_20260628.out",
    "05-knowledge/results/lrc14_canonical_corridor_fence_codex_20260628.out",
    "05-knowledge/results/lrc14_euler_mascheroni_wall_budget_codex_20260628.out",
    "05-knowledge/results/lrc14_euler_mascheroni_endpoint_spine_ledger_codex_20260628.out",
    "05-knowledge/results/lrc14_overlap_tax_harmonic_sieve_remainder_codex_20260628.out",
    "05-knowledge/results/lrc14_two_adic_branch_cover_certificate_codex_20260628.out",
    "05-knowledge/results/lrc14_overlap_menger_cut_certificate_codex_20260628.out",
    "04-computation/lrc14_special_function_cut_signature_recursion_codex_20260628.py",
    "04-computation/lrc14_owner_cut_resurrection_calculus_codex_20260628.py",
    "04-computation/lrc14_owner_cut_dual_certificate_synthesis_codex_20260628.py",
    "04-computation/lrc14_charal_owner_cut_recursion_prototype_codex_20260628.py",
    "04-computation/lrc14_transcendence_cut_chiral_synthesis_codex_20260628.py",
    "04-computation/lrc14_offgrid_resonance_transparency_rprime_closure_codex_20260628.py",
    "04-computation/lrc14_two_adic_offgrid_relocation_codex_20260628.py",
    "04-computation/lrc14_quniform_topology_arithmetic_break_codex_20260628.py",
    "04-computation/lrc14_covering_floor_duality_transfer_codex_20260628.py",
    "04-computation/lrc14_two_branch_obstruction_helly_codex_20260628.py",
    "04-computation/lrc14_additive_energy_sheet_sidecar_codex_20260628.py",
    "04-computation/lrc14_one_branch_mirror_endpoint_support_codex_20260628.py",
    "04-computation/lrc14_two_branch_wall_signature_atlas_codex_20260628.py",
    "04-computation/lrc14_two_adic_descent_loss_ledger_codex_20260628.py",
    "04-computation/lrc14_component_spine_certificate_codex_20260628.py",
    "04-computation/lrc14_euler_mascheroni_harmonic_intercept_codex_20260628.py",
    "04-computation/lrc14_canonical_corridor_fence_codex_20260628.py",
    "04-computation/lrc14_euler_mascheroni_wall_budget_codex_20260628.py",
    "04-computation/lrc14_euler_mascheroni_endpoint_spine_ledger_codex_20260628.py",
    "04-computation/lrc14_overlap_tax_harmonic_sieve_remainder_codex_20260628.py",
    "04-computation/lrc14_two_adic_branch_cover_certificate_codex_20260628.py",
    "04-computation/lrc14_overlap_menger_cut_certificate_codex_20260628.py",
    "07-reflections/the-one-inequality-that-completes-lrc14-and-why-resonance-dissolves.md",
    "07-reflections/the-galois-group-of-the-apex-prime-and-where-its-symmetry-breaks.md",
    "07-reflections/lrc14-owner-cut-resurrection-calculus-codex-20260628.md",
    "07-reflections/lrc14-owner-cut-dual-certificate-synthesis-codex-20260628.md",
    "07-reflections/lrc14-charal-owner-cut-recursion-prototype-codex-20260628.md",
    "07-reflections/lrc14-owner-cut-chiral-transcendence-synthesis-codex-20260628.md",
    "07-reflections/lrc14-offgrid-resonance-transparency-rprime-closure-codex-20260628.md",
    "07-reflections/lrc14-two-adic-offgrid-relocation-codex-20260628.md",
    "07-reflections/lrc14-quniform-topology-arithmetic-break-codex-20260628.md",
    "07-reflections/lrc14-covering-floor-duality-transfer-codex-20260628.md",
    "07-reflections/lrc14-two-branch-obstruction-helly-codex-20260628.md",
    "07-reflections/lrc14-additive-energy-sheet-sidecar-codex-20260628.md",
    "07-reflections/lrc14-one-branch-mirror-endpoint-support-codex-20260628.md",
    "07-reflections/lrc14-two-branch-wall-signature-atlas-codex-20260628.md",
    "07-reflections/lrc14-two-adic-descent-loss-ledger-codex-20260628.md",
    "07-reflections/lrc14-component-spine-certificate-codex-20260628.md",
    "07-reflections/lrc14-euler-mascheroni-harmonic-intercept-codex-20260628.md",
    "07-reflections/lrc14-canonical-corridor-fence-codex-20260628.md",
    "07-reflections/lrc14-euler-mascheroni-wall-budget-codex-20260628.md",
    "07-reflections/lrc14-euler-mascheroni-endpoint-spine-finite-part-codex-20260628.md",
    "07-reflections/lrc14-overlap-tax-harmonic-sieve-remainder-codex-20260628.md",
    "07-reflections/lrc14-two-adic-branch-cover-certificate-codex-20260628.md",
    "07-reflections/lrc14-overlap-menger-cut-certificate-codex-20260628.md",
    "agents/broadcast/MSG-1407-from-kind-pasteur-2026-06-28-kps-s258-the-critical-path-ma.md",
    "agents/broadcast/MSG-1408-from-kind-pasteur-2026-06-28-kps-s259-the-covering-floor-i.md",
    "comms/POKE-COORDINATION.md",
    "poke-forum/post_1782662454798.md",
]


OBLIGATION_WEIGHT = {
    "owner_support_exactness": 16,
    "endpoint_cut_resurrection": 15,
    "height_flex_leak": 12,
    "residue_breakpoint": 11,
    "boundary_order": 10,
    "branch_sheet_address": 9,
    "p_adic_stability": 8,
    "colored_variance_floor": 8,
    "finite_checkability": 8,
    "quotient_legality": 8,
    "bulk_prime_floor": 5,
    "factorization_gate": 5,
    "transcendence_guard": 5,
}


@dataclass(frozen=True)
class Motif:
    name: str
    family: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    corpus_terms: tuple[str, ...]
    theorem_use: str
    guardrail: str
    priority: int


MOTIFS = [
    Motif(
        name="menger_endpoint_cut_resurrection",
        family="min-cut / max-flow proof carrier",
        preserves=frozenset(
            {
                "owner_support_exactness",
                "endpoint_cut_resurrection",
                "finite_checkability",
                "quotient_legality",
                "residue_breakpoint",
            }
        ),
        destroys=frozenset({"colored_variance_floor", "bulk_prime_floor"}),
        corpus_terms=("Menger", "deletion cut", "owner support", "endpoint owner", "resurrection"),
        theorem_use=(
            "Turn the HYP-3406 owner-support repair into a min-cut theorem: "
            "a quotient may forget endpoint owners only when a bounded owner "
            "transversal makes every cut-code bucket theorem-exit pure."
        ),
        guardrail="a cut table is only useful if it proves a max-flow/min-cut dual, not just a label",
        priority=24,
    ),
    Motif(
        name="schwarz_christoffel_owner_polygon",
        family="boundary prevertex / polygonal side order carrier",
        preserves=frozenset(
            {
                "owner_support_exactness",
                "endpoint_cut_resurrection",
                "boundary_order",
                "finite_checkability",
                "quotient_legality",
                "height_flex_leak",
            }
        ),
        destroys=frozenset({"bulk_prime_floor"}),
        corpus_terms=("Schwarz", "Christoffel", "prevertex", "polygon", "boundary", "owner strip"),
        theorem_use=(
            "Read each actual-packet fiber as a polygonal boundary word.  "
            "The theorem exit must be invariant under legal prevertex moves, "
            "or the first changed side is the endpoint-owner debt."
        ),
        guardrail="conformal language cannot replace exact endpoint intervals or side order",
        priority=18,
    ),
    Motif(
        name="bring_radical_branch_kernel",
        family="quintic branch / monodromy address carrier",
        preserves=frozenset(
            {
                "branch_sheet_address",
                "height_flex_leak",
                "quotient_legality",
                "finite_checkability",
                "residue_breakpoint",
            }
        ),
        destroys=frozenset({"owner_support_exactness", "endpoint_cut_resurrection"}),
        corpus_terms=("Bring", "radical", "quintic", "resolvent", "branch", "monodromy"),
        theorem_use=(
            "Use Bring-radical thinking only as a warning about branch sheets: "
            "do not scalarize a mixed fiber unless the branch address and lost "
            "endpoint owner are carried."
        ),
        guardrail="generic quintic language is a danger signal, not a solver for the LRC14 packet",
        priority=10,
    ),
    Motif(
        name="barban_davenport_halberstam_variance_gate",
        family="variance-in-arithmetic-progressions sidecar",
        preserves=frozenset(
            {
                "colored_variance_floor",
                "residue_breakpoint",
                "bulk_prime_floor",
                "finite_checkability",
                "quotient_legality",
            }
        ),
        destroys=frozenset({"owner_support_exactness", "endpoint_cut_resurrection"}),
        corpus_terms=("Barban", "Davenport", "Halberstam", "variance", "colored", "discrepancy"),
        theorem_use=(
            "Use BDH-style variance as the colored CRT half-boundary floor "
            "after the endpoint cut is named; it should bound fiber spread, "
            "not choose the owner side."
        ),
        guardrail="variance cannot identify which endpoint-owner channel was forgotten",
        priority=9,
    ),
    Motif(
        name="krasner_hlw_no_free_slider_guard",
        family="p-adic stability / exponential independence guard",
        preserves=frozenset(
            {
                "p_adic_stability",
                "branch_sheet_address",
                "transcendence_guard",
                "height_flex_leak",
                "quotient_legality",
            }
        ),
        destroys=frozenset({"owner_support_exactness", "colored_variance_floor"}),
        corpus_terms=("Krasner", "Lindemann", "Weierstrass", "Hensel", "HLW", "no-free-slider"),
        theorem_use=(
            "Use Krasner stability to forbid hidden branch sliding under "
            "same-residue height moves, and use HLW only as a no-free-slider "
            "guard for exponential sidecars."
        ),
        guardrail=(
            "Lindemann-Weierstrass gives transcendence/linear independence "
            "statements, not the algebraic independence claimed by the forum note"
        ),
        priority=8,
    ),
    Motif(
        name="sophie_germain_quartic_factor_gate",
        family="quartic factorization / two-petal gate",
        preserves=frozenset(
            {
                "factorization_gate",
                "height_flex_leak",
                "branch_sheet_address",
                "finite_checkability",
            }
        ),
        destroys=frozenset({"owner_support_exactness", "bulk_prime_floor"}),
        corpus_terms=("Sophie", "Germain", "quartic", "factor", "phi4", "Lee-Yang"),
        theorem_use=(
            "Use a^4+4b^4 factorization as a two-petal sidecar for quartic "
            "height walls before invoking a higher branch equation."
        ),
        guardrail="factorization only helps after the exact owner or height coordinate is named",
        priority=5,
    ),
    Motif(
        name="ramanujan_soldner_mertens_bulk_zero",
        family="logarithmic-integral / prime-density threshold sidecar",
        preserves=frozenset(
            {
                "bulk_prime_floor",
                "colored_variance_floor",
                "finite_checkability",
            }
        ),
        destroys=frozenset({"owner_support_exactness", "endpoint_cut_resurrection", "height_flex_leak"}),
        corpus_terms=("Ramanujan", "Soldner", "Meissel", "Mertens", "logarithmic integral"),
        theorem_use=(
            "Treat the Ramanujan-Soldner and Meissel-Mertens constants as "
            "bulk normalization warnings for prime/log-density tails, not as "
            "local packet separators."
        ),
        guardrail="log-prime constants do not see the HYP-3406 owner-support leak",
        priority=-4,
    ),
]


@dataclass(frozen=True)
class PacketRow:
    name: str
    kernel: str
    residue_word: str
    v2_word: tuple[int, ...]
    height_word: str
    unit_slot: tuple[int, int, int]
    owner_support: tuple[str, ...]
    boundary_sector: str
    branch_sheet: str

    @property
    def coarse_base(self) -> tuple[str, str, str]:
        return ("eq14", "six_unit_boundary", "expanded_residue_failure")

    @property
    def endpoint_cut_word(self) -> tuple[str, ...]:
        # The cut word keeps just the owner channels, not theorem-exit labels.
        return tuple(sorted(channel.split(":")[0] for channel in self.owner_support))

    @property
    def charal_signature(self) -> tuple[object, ...]:
        # Charal = chiral/character/arc-lamination signature.  It is deliberately
        # recursive: residue branch, boundary side order, branch sheet, then cut.
        return (
            self.residue_word,
            self.boundary_sector,
            self.branch_sheet,
            self.endpoint_cut_word,
        )


ROWS = [
    PacketRow(
        "single swap 1->26",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 1, 2, 2, 3),
        "H13_to_26",
        (1, 2, 2),
        ("12:g2", "13:g1"),
        "right_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 1->40",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 2, 2, 3, 3),
        "H13_to_40",
        (1, 2, 2),
        ("12:g2", "13:g1", "2:g2"),
        "right_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 1->54",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 1, 2, 2, 3),
        "H13_to_26",
        (1, 2, 2),
        ("12:g2", "13:g1"),
        "right_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 3->26",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 1, 2, 2, 3),
        "H13_to_26",
        (2, 1, 2),
        ("11:g1", "12:g2", "13:g1", "6:g2"),
        "middle_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 3->40",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 2, 2, 3, 3),
        "H13_to_40",
        (2, 1, 2),
        ("11:g1", "12:g2", "13:g1", "6:g2"),
        "middle_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 3->54",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 1, 2, 2, 3),
        "H13_to_26",
        (2, 1, 2),
        ("11:g1", "12:g2", "13:g1", "6:g2"),
        "middle_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 5->26",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 1, 2, 2, 3),
        "H13_to_26",
        (2, 2, 1),
        ("10:g2", "11:g1", "12:g2", "13:g1", "9:g1"),
        "left_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 5->40",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 2, 2, 3, 3),
        "H13_to_40",
        (2, 2, 1),
        ("10:g2", "11:g1", "12:g2", "13:g1", "9:g1"),
        "left_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 5->54",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 1, 2, 2, 3),
        "H13_to_26",
        (2, 2, 1),
        ("10:g2", "11:g1", "12:g2", "9:g1"),
        "left_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 9->54",
        "positive-Haar-open",
        "R13_to_26_family",
        (0, 1, 1, 1, 1, 2, 2, 3),
        "H13_to_26",
        (2, 2, 1),
        ("10:g2", "11:g1", "12:g2", "13:g1", "5:g1", "7:g7", "8:g2"),
        "left_owner_cover",
        "covering_sheet",
    ),
    PacketRow(
        "petal 13->26",
        "unit-petal-named",
        "R13_to_26_family",
        (0, 1, 1, 1, 1, 2, 2, 3),
        "H13_to_26",
        (1, 2, 2),
        ("12:g2", "1:g1"),
        "right_owner_petal",
        "unit_petal_sheet",
    ),
    PacketRow(
        "GW-shell alias 12->132",
        "positive-Haar-open",
        "R12_to_132_family",
        (0, 1, 1, 1, 2, 2, 3),
        "H12_to_132",
        (2, 2, 2),
        ("11:g1", "13:g1", "5:g1", "6:g2", "7:g7"),
        "gw_cover_shell",
        "covering_sheet",
    ),
    PacketRow(
        "single swap 12->48",
        "positive-Haar-open",
        "R12_to_132_family",
        (0, 1, 1, 1, 2, 3, 4),
        "H12_to_48",
        (2, 2, 2),
        ("5:g1", "6:g2", "7:g7"),
        "gw_cover_shell",
        "covering_sheet",
    ),
    PacketRow(
        "P10+GW",
        "unit-petal-named",
        "R12_to_132_family",
        (0, 1, 1, 2, 2, 3, 3),
        "H10_plus_GW",
        (2, 2, 2),
        ("6:g2", "7:g7"),
        "gw_unit_petal",
        "unit_petal_sheet",
    ),
]


SIGNATURE_LEVELS = [
    ("C0_coarse", ("coarse_base",)),
    ("C1_residue_branch", ("coarse_base", "residue_word")),
    ("C2_height_or_v2", ("coarse_base", "residue_word", "height_word")),
    ("C3_unit_slot", ("coarse_base", "residue_word", "height_word", "unit_slot")),
    ("C4_endpoint_cut", ("coarse_base", "residue_word", "endpoint_cut_word")),
    ("C5_charal_signature", ("coarse_base", "charal_signature")),
]


def corpus_text() -> str:
    chunks: list[str] = []
    for rel in KEY_FILES:
        path = ROOT / rel
        if path.exists():
            chunks.append(path.read_text(encoding="utf-8", errors="ignore").lower())
    return "\n".join(chunks)


def corpus_hits(motif: Motif, text: str) -> int:
    return sum(text.count(term.lower()) for term in motif.corpus_terms)


def motif_score(motif: Motif, text: str) -> int:
    base = sum(OBLIGATION_WEIGHT[key] for key in motif.preserves)
    penalty = sum(OBLIGATION_WEIGHT[key] for key in motif.destroys)
    hit_bonus = min(corpus_hits(motif, text), 100) // 4
    return base - penalty + motif.priority + hit_bonus


def value(row: PacketRow, field: str) -> object:
    attr = getattr(row, field)
    return attr() if callable(attr) else attr


def key(row: PacketRow, fields: tuple[str, ...]) -> tuple[object, ...]:
    out: list[object] = []
    for field in fields:
        field_value = value(row, field)
        if isinstance(field_value, tuple):
            out.append(field_value)
        else:
            out.append(field_value)
    return tuple(out)


def fibers(rows: list[PacketRow], fields: tuple[str, ...]) -> dict[tuple[object, ...], list[PacketRow]]:
    grouped: dict[tuple[object, ...], list[PacketRow]] = defaultdict(list)
    for row in rows:
        grouped[key(row, fields)].append(row)
    return grouped


def mixed_count(rows: list[PacketRow], fields: tuple[str, ...]) -> int:
    return sum(
        1
        for fiber in fibers(rows, fields).values()
        if len({row.kernel for row in fiber}) > 1
    )


def largest_mixed(rows: list[PacketRow], fields: tuple[str, ...]) -> int:
    mixed_sizes = [
        len(fiber)
        for fiber in fibers(rows, fields).values()
        if len({row.kernel for row in fiber}) > 1
    ]
    return max(mixed_sizes) if mixed_sizes else 0


def minimal_repairs(rows: list[PacketRow], baseline: tuple[str, ...], candidates: list[str]) -> list[tuple[str, ...]]:
    good: list[tuple[str, ...]] = []
    for size in range(1, len(candidates) + 1):
        for combo in combinations(candidates, size):
            fields = baseline + combo
            if mixed_count(rows, fields) == 0:
                good.append(combo)
        if good:
            return good
    return []


def tournament_edges(ranked: list[Motif]) -> set[tuple[str, str]]:
    return {(a.name, b.name) for a, b in combinations(ranked, 2)}


def directed_triangles(names: list[str], edges: set[tuple[str, str]]) -> int:
    total = 0
    for a, b, c in combinations(names, 3):
        total += int((a, b) in edges and (b, c) in edges and (c, a) in edges)
        total += int((a, c) in edges and (c, b) in edges and (b, a) in edges)
    return total


def hamiltonian_path_count(names: list[str], edges: set[tuple[str, str]]) -> int:
    return int(all((a, b) in edges for a, b in zip(names, names[1:])))


def main() -> None:
    text = corpus_text()
    ranked = sorted(MOTIFS, key=lambda motif: (-motif_score(motif, text), motif.name))
    names = [motif.name for motif in ranked]
    edges = tournament_edges(ranked)
    scores = {motif.name: motif_score(motif, text) for motif in MOTIFS}

    print("HYP-3440 BRING / SCHWARZ-CHRISTOFFEL / MENGER CHARAL SIGNATURE SCOUT")
    print("=" * 78)
    print("status=synthesis / recursive-sidecar router; not an LRC14 proof")
    print("new_handles=HYP-3440,T1401,LTI-401,LTT-301")
    print("base=HYP-3406 enlarged-bank residue failure and owner-support repair")
    print(
        "upstream=HYP-3407/HYP-3408/HYP-3409/HYP-3410/HYP-3411/"
        "HYP-3412/HYP-3413/HYP-3414/HYP-3416/HYP-3417/HYP-3419/"
        "HYP-3420/HYP-3421/HYP-3422/HYP-3423/HYP-3424/HYP-3425/HYP-3426/HYP-3427/HYP-3428/HYP-3429/"
        "HYP-3430-HYP-3437 boundary, guardrail, recursive-sidecar, charal-recursion, "
        "Galois, special-function, GW-gate, owner-cut, quotient-ladder, "
        "owner-current, cut-decision-tree, chiral-owner, off-grid-transparency, "
        "two-adic-relocation, q-uniform-topology guardrail, floor-transfer, "
        "Helly-obstruction, energy-sheet packet, one-branch mirror, "
        "wall-signature, two-adic-loss, endpoint-spine, scalar-firewall, "
        "corridor-fence, harmonic-budget, overlap-tax, branch-cover, "
        "bad-core extractor, and overlap Menger-cut atlases"
    )
    print(
        "critical_path=HYP-3415/HYP-3418: non-covering is closed by q-witness; "
        "covering needs the single decorrelation floor inequality, now a "
        "2-adic even-speed descent problem rather than an apex-7 census problem"
    )
    print(
        "topology_guardrail=HYP-3423: SC/Menger boundary language may certify "
        "q-uniform order or owner-cut legality, but q-specific magnitude and "
        "floor claims must restore HYP-3413 arithmetic, HYP-3417 owner currents, "
        "or HYP-3418/HYP-3421/HYP-3422 two-adic floor data"
    )
    print(
        "incoming_extension=HYP-3424 makes floor transfer the parent router; "
        "HYP-3425 Helly supplies the branch-relocation interval target; "
        "parallel HYP-3425 additive energy requires energy-plus-sheet packets. "
        "HYP-3426 reduces branch choice to one-branch endpoint-support triples; "
        "HYP-3427 upgrades positive windows to wall signatures; HYP-3428 "
        "records halving/descent loss classes; HYP-3429 compresses survivor "
        "windows to endpoint spines; HYP-3430-HYP-3434 block scalar harmonic "
        "compression; HYP-3435 turns branch survival into endpoint-gate "
        "certificates; HYP-3437 reserves the overlap-tax Menger-cut extractor. "
        "HYP-3440 endpoint cuts may feed those packets only as owner/branch/wall debt."
    )
    print()

    print("SELECTED ANGLES")
    for motif in ranked[:2]:
        print(f"  {motif.name} score={scores[motif.name]}")
        print(f"    family={motif.family}")
        print(f"    preserves={','.join(sorted(motif.preserves))}")
        print(f"    destroys={','.join(sorted(motif.destroys))}")
        print(f"    theorem_use={motif.theorem_use}")
        print(f"    guardrail={motif.guardrail}")
    print()

    print("FULL MOTIF RANKING")
    for motif in ranked:
        print(
            f"  {motif.name:42s} score={scores[motif.name]:4d} "
            f"corpus_hits={corpus_hits(motif, text):4d}"
        )
    print()

    print("RECURSIVE CHARAL SIGNATURE LADDER")
    for label, fields in SIGNATURE_LEVELS:
        print(
            f"  {label:22s} fibers={len(fibers(ROWS, fields)):2d} "
            f"mixed={mixed_count(ROWS, fields):1d} "
            f"largest_mixed={largest_mixed(ROWS, fields):2d} "
            f"fields={fields}"
        )
    print()

    print("MINIMAL REPAIRS OVER C1_residue_branch")
    candidates = [
        "v2_word",
        "height_word",
        "unit_slot",
        "owner_support",
        "endpoint_cut_word",
        "boundary_sector",
        "branch_sheet",
        "charal_signature",
    ]
    for repair in minimal_repairs(ROWS, ("coarse_base", "residue_word"), candidates):
        print(f"  {repair}")
    print()

    print("MIXED FIBER DIAGNOSTICS AT C2_height_or_v2")
    c2_fields = ("coarse_base", "residue_word", "height_word")
    for fiber_key, fiber in fibers(ROWS, c2_fields).items():
        kernels = Counter(row.kernel for row in fiber)
        if len(kernels) > 1:
            print(f"  fiber_key={fiber_key} size={len(fiber)} kernels={dict(kernels)}")
            for row in fiber:
                print(
                    "    "
                    f"{row.kernel:19s} {row.name:24s} "
                    f"unit_slot={row.unit_slot} cut={row.endpoint_cut_word}"
                )
    print()

    print("REQUESTED MOTIF INTERPRETATION")
    print("  Bring radical: branch-sheet warning for illegal quintic scalarization.")
    print("  Schwarz-Christoffel: owner supports are polygon sides / prevertex order.")
    print("  Barban-Davenport-Halberstam: variance gate for colored CRT bulk fibers.")
    print("  Menger cuts: exact carrier for endpoint-owner resurrection.")
    print("  Ramanujan-Soldner + Meissel-Mertens: prime/log bulk normalization only.")
    print("  Sophie Germain: quartic two-petal factor gate for height-wall sidecars.")
    print("  Krasner: p-adic branch stability for same-residue height moves.")
    print("  Hermite-Lindemann-Weierstrass: no-free-slider guard, not an owner separator.")
    print()

    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof motifs / recursive sidecars, not runners or arcs")
    print(
        "  pairwise_observable=retained current proof obligations minus destroyed "
        "owner/height/variance/stability coordinates"
    )
    print("  switch/gauge=higher motif score; ties use motif name")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={directed_triangles(names, edges)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(names, edges)}")
    print("  priority_path=" + " -> ".join(names))
    print()

    print("CONCLUSION")
    print(
        "  The best HYP-3440 route is to keep the lower-payload motif router "
        "subordinate to HYP-3419/HYP-3420's concrete owner-cut evidence, "
        "HYP-3421/HYP-3422's off-grid/two-adic floor route, HYP-3423's "
        "topology/arithmetic guardrail, HYP-3424's floor-transfer router, "
        "the paired HYP-3425 Helly / energy-sheet sidecars, HYP-3426's "
        "one-branch mirror endpoint-support audit, HYP-3427's wall-signature atlas, "
        "HYP-3428's two-adic descent loss ledger, HYP-3429's endpoint-spine "
        "certificate, HYP-3435's branch-cover certificate, and HYP-3437's "
        "overlap-tax Menger-cut target: "
        "formalize owner-support exactness as a "
        "bounded owner-transversal / Menger min-cut / Schwarz-Christoffel "
        "boundary-order theorem, then feed that certificate into the recursive "
        "quotient ladder.  HYP-3414 blocks a universal singleton-owner shortcut, "
        "and HYP-3417 gives margin-1 owner-current witnesses on the current "
        "fibers.  This remains off the LRC14 critical path unless it helps the "
        "HYP-3415/HYP-3418 covering floor: the completion theorem is the uniform "
        "R' > 0 decorrelation inequality, sharpened to a 2-adic even-speed "
        "descent route.  Bring radicals, Krasner stability, and HLW enter as "
        "branch/no-free-slider guards after the endpoint cut is retained; BDH "
        "and Meissel-Mertens belong to named bulk/floor, wall-signature, or energy-plus-sheet "
        "sidecars, not to the local HYP-3406 leak."
    )


if __name__ == "__main__":
    main()
