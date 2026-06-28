#!/usr/bin/env python3
"""HYP-3407 scout: boundary-uniformization cut-stability atlas for LRC14.

This is a research router, not a proof.  It integrates the newest
HYP-3405/HYP-3406 frontier with the user's requested lenses:

  * Bring radical / quintic branch data;
  * Schwarz-Christoffel boundary polygons;
  * Barban-Davenport-Halberstam mean-square residue control;
  * Menger cuts and endpoint-owner resurrection;
  * recursive chiral signatures;
  * Ramanujan-Soldner, Sophie Germain, Hermite-Lindemann-Weierstrass,
    Krasner, and Meissel-Mertens guardrails.

Tournament Analysis uses proof carriers as vertices.  The chosen observable is
the amount of HYP-3405/HYP-3406 payload a carrier preserves, minus the amount
of load-bearing sidecar data it would forget.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True)
class Carrier:
    code: str
    name: str
    move: str
    preserves: tuple[str, ...]
    forgets_if_naive: tuple[str, ...]
    first_test: str
    theorem_pull: str
    risk: str
    scores: dict[str, int]

    @property
    def total(self) -> int:
        weights = {
            "h3405_fit": 4,
            "h3406_fit": 5,
            "finite_test": 3,
            "global_path": 3,
            "sidecar_legality": 3,
            "analytic_or_local_control": 2,
            "recursion_payload": 2,
            "risk_penalty": -3,
        }
        return sum(weights[key] * value for key, value in self.scores.items())


def carriers() -> list[Carrier]:
    return [
        Carrier(
            "C01",
            "Boundary-uniformization Menger zipper",
            (
                "Model the forbidden arcs as a Schwarz-Christoffel contact polygon, "
                "then apply endpoint-owner min-cuts to decide which boundary "
                "channels must be resurrected when residue/height data collide."
            ),
            (
                "AP/GW boundary atoms",
                "strict-open intervals",
                "endpoint-owner support",
                "Menger cut value",
            ),
            ("branch sheet", "prime-channel variance"),
            (
                "On the HYP-3406 petal 13->26 vs single-swap 26/40/54 fiber, "
                "build the owner-support bipartite graph and compute the minimum "
                "endpoint cut separating unit-petal-named from positive-Haar-open."
            ),
            (
                "Every expanded-bank residue+height collision has an owner-support "
                "Menger cut; the cut either resurrects the missing owner sidecar, "
                "lands on AP/GW, or names a finite residual family."
            ),
            "Needs exact owner graph extraction for more HYP-2963 rows.",
            {
                "h3405_fit": 3,
                "h3406_fit": 4,
                "finite_test": 4,
                "global_path": 3,
                "sidecar_legality": 4,
                "analytic_or_local_control": 1,
                "recursion_payload": 2,
                "risk_penalty": 1,
            },
        ),
        Carrier(
            "C02",
            "Krasner collar stability",
            (
                "Treat AP-collar rows as local p-adic disks around contact roots. "
                "Krasner-style stability says a row cannot change theorem exit "
                "while all load-bearing roots stay inside the same labelled disks."
            ),
            (
                "unit-height lift",
                "nonunit-height packet",
                "local root disk",
                "sidecar repair matrix",
            ),
            ("global owner support", "SC accessory parameters"),
            (
                "Attach p-adic disk labels to the HYP-3405 AP versus 13->27 "
                "fiber and verify that the unit-height lift (13,0)->(13,1) is "
                "exactly the first disk-exit coordinate."
            ),
            (
                "AP-collar strict-open rows remain strict-open under local disk "
                "stability; boundary-vs-strict changes require a named disk exit, "
                "which HYP-3405 identifies as unit height on the first fiber."
            ),
            "Local p-adic stability must be tied to real safe-interval exits.",
            {
                "h3405_fit": 4,
                "h3406_fit": 2,
                "finite_test": 4,
                "global_path": 2,
                "sidecar_legality": 4,
                "analytic_or_local_control": 4,
                "recursion_payload": 1,
                "risk_penalty": 1,
            },
        ),
        Carrier(
            "C03",
            "BDH-Mertens owner discrepancy",
            (
                "Use Barban-Davenport-Halberstam as the model for a mean-square "
                "owner-residue discrepancy bound, with Meissel-Mertens prime-channel "
                "normalization tracking the cost of enlarging residue banks."
            ),
            (
                "residue word",
                "owner-support word",
                "mean-square discrepancy",
                "exceptional fibers",
            ),
            ("single worst fiber", "unit-height lift"),
            (
                "For the HYP-3406 enlarged banks, compute an L2 discrepancy "
                "between residue-only fiber exits and owner-support repaired exits; "
                "test whether the exceptional mass is confined to the two known "
                "families."
            ),
            (
                "Residue-only forgetting is allowed only with a BDH-style "
                "mean-square bound plus an explicit exceptional-fiber ledger; "
                "the first exceptional families must be height and owner leaks."
            ),
            "Average control cannot replace finite exceptional classification.",
            {
                "h3405_fit": 1,
                "h3406_fit": 4,
                "finite_test": 3,
                "global_path": 4,
                "sidecar_legality": 3,
                "analytic_or_local_control": 4,
                "recursion_payload": 1,
                "risk_penalty": 2,
            },
        ),
        Carrier(
            "C04",
            "Recursive chiral signature deck",
            (
                "Extend the HYP-3123/HYP-3124 cross-sector orientation lesson: "
                "recursively delete endpoint owners and record mirror/converse "
                "signature words so owner support cannot collapse with its mirror."
            ),
            (
                "endpoint-owner support",
                "cross-sector orientation",
                "tail/tip child packets",
                "chiral guard word",
            ),
            ("analytic density", "p-adic disk size"),
            (
                "On the petal 13->26 fiber, compute owner-support child decks after "
                "deleting each active endpoint owner and check whether the recursive "
                "chiral word separates petal-named from positive-open rows."
            ),
            (
                "Owner support becomes recursive: if a quotient forgets endpoint "
                "orientation, the chiral deck resurrects it or produces a named "
                "mirror-collapse debt."
            ),
            "Could overfit unless tested on a larger HYP-2963 bank.",
            {
                "h3405_fit": 1,
                "h3406_fit": 4,
                "finite_test": 3,
                "global_path": 3,
                "sidecar_legality": 4,
                "analytic_or_local_control": 1,
                "recursion_payload": 4,
                "risk_penalty": 1,
            },
        ),
        Carrier(
            "C05",
            "Bring branch-sheet packet",
            (
                "Use the Bring radical only as a branch-locus template: proof rows "
                "live on labelled sheets, and forbidden quotient moves are the ones "
                "that identify sheets with different theorem exits."
            ),
            (
                "branch sheet",
                "monodromy word",
                "resolvent sidecar",
                "theorem exit",
            ),
            ("owner support", "strict-open mass"),
            (
                "Represent HYP-3405 boundary rows and HYP-3406 owner leaks as sheet "
                "changes in a small branch graph; test whether the known leaks are "
                "exactly the first nontrivial monodromy edges."
            ),
            (
                "A legal quotient cannot cross a branch cut unless it keeps the "
                "sheet label, proves sheet-invariance, or discharges by AP/GW, "
                "strict-open mass, owner cut, or state-lift debt."
            ),
            "A quintic metaphor is dangerous unless reduced to finite sheet labels.",
            {
                "h3405_fit": 2,
                "h3406_fit": 2,
                "finite_test": 2,
                "global_path": 2,
                "sidecar_legality": 4,
                "analytic_or_local_control": 1,
                "recursion_payload": 3,
                "risk_penalty": 2,
            },
        ),
        Carrier(
            "C06",
            "Sophie-Germain quartic split",
            (
                "Use a^4+4b^4=(a^2-2ab+2b^2)(a^2+2ab+2b^2) as a warning that "
                "quartic/Fejer moments may be secretly two owner channels, not one "
                "scalar mass."
            ),
            (
                "quartic moment",
                "two quadratic owner factors",
                "sign/chiral split",
            ),
            ("exact period", "mean-square sieve"),
            (
                "For HYP-3406 owner-support words, split fourth-moment owner "
                "payloads into plus/minus quadratic factors and test whether one "
                "factor is the petal-vs-positive separator."
            ),
            (
                "Every quartic certificate must either factor into legal owner "
                "channels or retain a named sign/chiral sidecar."
            ),
            "Likely a local algebraic guardrail, not a terminal proof route.",
            {
                "h3405_fit": 1,
                "h3406_fit": 3,
                "finite_test": 3,
                "global_path": 1,
                "sidecar_legality": 4,
                "analytic_or_local_control": 1,
                "recursion_payload": 2,
                "risk_penalty": 1,
            },
        ),
        Carrier(
            "C07",
            "HLW exponential-independence guardrail",
            (
                "Hermite-Lindemann-Weierstrass warns that independent exponential "
                "or Fourier channels cannot be collapsed to a small algebraic "
                "scalar unless the period lattice and relations are retained."
            ),
            (
                "Fourier mode",
                "period lattice",
                "linear relation status",
                "exact-period boundary",
            ),
            ("owner support", "finite cut"),
            (
                "Audit any analytic-lifting ledger that separates HYP-3406 rows: "
                "does it prove relation invariance, or is it only a row hash over "
                "independent modes?"
            ),
            (
                "Analytic exponential packets are legal only when their period "
                "lattice, relation rank, and theorem exit are carried or discharged."
            ),
            "Mostly a guardrail against false simplification.",
            {
                "h3405_fit": 1,
                "h3406_fit": 1,
                "finite_test": 2,
                "global_path": 2,
                "sidecar_legality": 4,
                "analytic_or_local_control": 3,
                "recursion_payload": 1,
                "risk_penalty": 2,
            },
        ),
        Carrier(
            "C08",
            "Ramanujan-Soldner critical normalizer",
            (
                "Use the Ramanujan-Soldner zero of li only as a sign-change "
                "normalizer: it can mark when an averaged prime-channel correction "
                "changes sign, not as a standalone LRC invariant."
            ),
            (
                "critical zero",
                "smoothing sign",
                "prime-channel normalization",
            ),
            ("finite owner cut", "height lift"),
            (
                "Add a critical-normalization column to the BDH-Mertens discrepancy "
                "toy table and verify that it does not separate known fibers unless "
                "owner or height sidecars are also present."
            ),
            (
                "Named constants are allowed as calibration columns only after the "
                "finite LRC predicate they preserve is stated."
            ),
            "High numerology risk; useful mainly as a calibration discipline.",
            {
                "h3405_fit": 0,
                "h3406_fit": 1,
                "finite_test": 2,
                "global_path": 1,
                "sidecar_legality": 3,
                "analytic_or_local_control": 2,
                "recursion_payload": 0,
                "risk_penalty": 2,
            },
        ),
        Carrier(
            "C09",
            "Schwarz-Christoffel accessory audit",
            (
                "Separate true boundary angles from accessory parameters.  If an "
                "SC polygon model changes endpoint-owner support through an "
                "accessory parameter, that parameter must be a sidecar."
            ),
            (
                "prevertex order",
                "interior angle",
                "endpoint-owner channel",
                "accessory parameter",
            ),
            ("prime variance", "Bring sheet"),
            (
                "Build the contact polygon for AP, GW 12->24, 12->36, and 13->27; "
                "identify which angle data are identical and which accessory field "
                "separates boundary from strict-open."
            ),
            (
                "Boundary angle quotients are legal only with accessory-parameter "
                "debt discharged by exact rational strict-open intervals or owner "
                "cuts."
            ),
            "SC maps can hide hard data in accessory parameters.",
            {
                "h3405_fit": 4,
                "h3406_fit": 2,
                "finite_test": 3,
                "global_path": 2,
                "sidecar_legality": 4,
                "analytic_or_local_control": 2,
                "recursion_payload": 1,
                "risk_penalty": 1,
            },
        ),
        Carrier(
            "C10",
            "Meissel-Mertens prime-channel budget",
            (
                "Use sum_p<=x 1/p = log log x + B1 + o(1) as a channel-rank "
                "budget: each new prime or modulus layer costs a controlled "
                "loglog amount and must buy a real sidecar."
            ),
            (
                "prime-channel count",
                "modulus budget",
                "sieve normalization",
            ),
            ("finite exceptional fibers", "SC geometry"),
            (
                "Compare sidecar payload cost for residue, v2, height, owner "
                "support, and exact period against the mixed-fiber reduction they "
                "actually purchase on HYP-3406 banks."
            ),
            (
                "Sieve layers are legal only when their prime-channel budget buys "
                "a named reduction in mixed theorem-exit fibers."
            ),
            "Budgeting does not by itself prove positivity.",
            {
                "h3405_fit": 0,
                "h3406_fit": 3,
                "finite_test": 4,
                "global_path": 2,
                "sidecar_legality": 3,
                "analytic_or_local_control": 3,
                "recursion_payload": 0,
                "risk_penalty": 1,
            },
        ),
    ]


def edge(a: Carrier, b: Carrier) -> tuple[str, str]:
    if a.total != b.total:
        return (a.code, b.code) if a.total > b.total else (b.code, a.code)
    # Tie Hamiltonian path is lexical in carrier code.
    return (a.code, b.code) if a.code < b.code else (b.code, a.code)


def tournament_fingerprint(items: list[Carrier]) -> dict[str, object]:
    codes = [item.code for item in items]
    edges = {frozenset((a.code, b.code)): edge(a, b) for a, b in combinations(items, 2)}
    cycles = 0
    for x, y, z in combinations(codes, 3):
        directed = {(x, y): False, (y, x): False, (x, z): False, (z, x): False, (y, z): False, (z, y): False}
        for a, b in ((x, y), (x, z), (y, z)):
            u, v = edges[frozenset((a, b))]
            directed[(u, v)] = True
        if (directed[(x, y)] and directed[(y, z)] and directed[(z, x)]) or (
            directed[(x, z)] and directed[(z, y)] and directed[(y, x)]
        ):
            cycles += 1
    ordered = sorted(items, key=lambda item: (-item.total, item.code))
    return {
        "vertices": len(items),
        "score_hist": dict(sorted(Counter(item.total for item in items).items())),
        "directed_3cycles": cycles,
        "hamiltonian_path_count": 1 if cycles == 0 else "not counted",
        "priority_path": " -> ".join(item.code for item in ordered),
    }


def generated_tests(items: list[Carrier]) -> list[str]:
    axes = [
        ("branch", "sheet label", "Bring branch graph"),
        ("boundary", "endpoint owner", "SC polygon/accessory audit"),
        ("cut", "min cut", "Menger owner graph"),
        ("local", "p-adic disk", "Krasner collar"),
        ("analytic", "mean-square deficit", "BDH-Mertens ledger"),
        ("recursive", "chiral child deck", "tail/tip owner recursion"),
    ]
    debts = [
        "unit-height lift",
        "endpoint-owner support",
        "exact-period boundary",
        "off-grid floor",
        "state-lift/H7 label",
    ]
    out = []
    for i, (axis, sidecar, test) in enumerate(axes, start=1):
        debt = debts[(i - 1) % len(debts)]
        out.append(
            f"G{i:02d} {axis}: retain {sidecar}; first failure should be {debt}; "
            f"finite test = {test} on HYP-3405/HYP-3406 named fibers."
        )
    out.append(
        "G07 zipper theorem: compose boundary polygon -> owner min-cut -> "
        "Krasner local disk -> BDH exceptional ledger; any unclosed row must "
        "emit a named chiral or state-lift debt."
    )
    return out


def main() -> None:
    items = carriers()
    ordered = sorted(items, key=lambda item: (-item.total, item.code))
    print("HYP-3407 BOUNDARY-UNIFORMIZATION CUT-STABILITY ATLAS")
    print("=" * 78)
    print("status=synthesis / executable proof-carrier router; not an LRC14 proof")
    print("frontier_inputs=HYP-3405 AP-collar certificate + HYP-3406 owner-support repair")
    print("new_lenses=Bring radical, Schwarz-Christoffel, BDH, Menger cuts, chiral recursion")
    print()

    print("CURRENT HARD FACTS")
    print("  HYP-3405: AP and GW 12->24 are the only boundary-tight AP-collar rows.")
    print("  HYP-3405: every other AP-collar row through speed 84 is strict-open.")
    print("  HYP-3405: the first AP-collar missing coordinate is unit-height (13,0)->(13,1).")
    print("  HYP-3406: residue-only fails on enlarged banks.")
    print("  HYP-3406: height/v2 repairs the first leak but not the stronger owner leak.")
    print("  HYP-3406: residue + owner_support is exact through the scanned 872-row bank.")
    print()

    print("RANKED CARRIERS")
    for item in ordered:
        print(f"  {item.code} score={item.total:3d}  {item.name}")
        print(f"      move: {item.move}")
        print(f"      first_test: {item.first_test}")
    print()

    print("PROCEDURALLY GENERATED NEXT TESTS")
    for line in generated_tests(items):
        print(f"  {line}")
    print()

    print("BOUNDARY-UNIFORMIZATION THEOREM TARGET")
    print("  For every primitive expanded-bank packet after q-witness and AP/GW exits:")
    print("    either residue+owner_support is theorem-exit exact,")
    print("    or the first failure is witnessed by one of:")
    print("      unit-height disk exit, endpoint-owner Menger cut, SC accessory debt,")
    print("      exact-period/BDH exceptional fiber, recursive chiral mirror debt,")
    print("      state-lift/H7 label, or a newly named finite residual.")
    print("  A quotient is proof-legal only if it preserves that sidecar, proves it")
    print("  irrelevant on the fiber, or resurrects it by an explicit cut/stability map.")
    print()

    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof carriers, not runners, arcs, raw constants, or raw residues")
    print("  pairwise_observable=retained HYP-3405/HYP-3406 payload minus forgotten sidecar debt")
    print("  switch_gauge=A -> B iff A has larger weighted carrier score; ties use code order")
    fingerprint = tournament_fingerprint(items)
    for key, value in fingerprint.items():
        print(f"  {key}={value}")
    print()

    print("ASSUMPTION CHALLENGE")
    print("  Alternate vertices considered:")
    print("    runners, gaps, fixed circle sections, section boundaries, wall crossings,")
    print("    residues, cover arcs, Fourier modes, matroid circuits, endpoint cuts,")
    print("    SC prevertices, p-adic root disks, branch sheets, prime moduli,")
    print("    chiral child decks, and proof obligations.")
    print("  Chosen vertices:")
    print("    proof carriers, because the preserved predicate is theorem-exit exactness.")
    print("  Preserved predicate:")
    print("    boundary-tight vs strict-open vs positive-Haar-open vs unit-petal/K33/H7 debt.")
    print("  Destroyed by naive quotients:")
    print("    unit height, endpoint owner, accessory parameter, exact period, off-grid floor,")
    print("    branch sheet, chiral orientation, and exceptional-fiber identity.")
    print()

    print("CONCLUSION")
    print("  The requested lenses converge on one disciplined move:")
    print("    do not ask a quotient to remember less unless a labelled packet theorem")
    print("    says exactly how the forgotten coordinate returns.")
    print("  The highest-leverage next computation is the owner-support Menger graph on")
    print("  the HYP-3406 petal/single-swap fibers, paired with the HYP-3405 Krasner-style")
    print("  unit-height collar audit.")


if __name__ == "__main__":
    main()
