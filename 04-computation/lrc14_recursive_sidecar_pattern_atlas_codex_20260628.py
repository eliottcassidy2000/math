#!/usr/bin/env python3
"""HYP-3409: recursive sidecar pattern atlas for LRC14.

This is a proof-pattern router, not a proof.  It extracts the recursive shape
shared by HYP-3405, HYP-3406, HYP-3407, and HYP-3408:

    legal quotient -> mixed theorem-exit fiber -> first missing sidecar
    -> repaired quotient -> next quotient.

Tournament Analysis uses recursion operators / proof obligations as vertices,
not runners, arcs, constants, or raw residues.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


WEIGHTS = {
    "exit_purity": 5,
    "sidecar_precision": 4,
    "finite_testability": 4,
    "theorem_pull": 3,
    "recursive_depth": 2,
    "risk_penalty": -3,
}


@dataclass(frozen=True)
class Pattern:
    code: str
    name: str
    recursive_move: str
    preserves: tuple[str, ...]
    destroys_if_naive: tuple[str, ...]
    first_test: str
    theorem_pull: str
    challenged_assumption: str
    scores: dict[str, int]

    @property
    def total(self) -> int:
        return sum(WEIGHTS[key] * value for key, value in self.scores.items())


def patterns() -> list[Pattern]:
    return [
        Pattern(
            code="R00",
            name="mixed-fiber resurrection loop",
            recursive_move=(
                "Treat every compression failure as a typed sidecar request: "
                "quotient rows, inspect mixed theorem-exit fibers, restore the "
                "first destroyed coordinate, then quotient again."
            ),
            preserves=("theorem exit", "boundary/strict status", "named debt"),
            destroys_if_naive=("unit-height flex", "endpoint owner", "exact period"),
            first_test=(
                "Replay the HYP-3405 AP-vs-13->27 leak and HYP-3406 enlarged "
                "owner leaks under the same quotient/fiber/repaired-quotient API."
            ),
            theorem_pull=(
                "A finite lemma can be stated as: every mixed fiber either gains "
                "one declared sidecar or exits through AP/GW, strict-open mass, "
                "state-lift debt, or a named residual."
            ),
            challenged_assumption=(
                "The proof object is not a single best invariant; it is a stack "
                "of legal forgetful maps with local repair obligations."
            ),
            scores={
                "exit_purity": 5,
                "sidecar_precision": 5,
                "finite_testability": 5,
                "theorem_pull": 5,
                "recursive_depth": 4,
                "risk_penalty": 0,
            },
        ),
        Pattern(
            code="R01",
            name="owner-cut recursion",
            recursive_move=(
                "When residue/height packets still mix exits, split the fiber by "
                "endpoint-owner support and then recurse on owner deletion, "
                "Menger cut, or tail-tip child decks."
            ),
            preserves=("endpoint owner support", "exit purity", "cut boundary"),
            destroys_if_naive=("which endpoint carries the obstruction", "tail/tip side"),
            first_test=(
                "Build the exact Menger-style owner graph for HYP-3406's "
                "petal 13->26 and petal 10->20 leaks."
            ),
            theorem_pull=(
                "Candidate theorem: after residue and height sidecars, every "
                "persistent enlarged-bank leak is an owner cut, an owner-child "
                "deck, or a named finite residual."
            ),
            challenged_assumption=(
                "Owner labels are not cosmetic row metadata; they are the next "
                "recursive boundary coordinate after height stops helping."
            ),
            scores={
                "exit_purity": 5,
                "sidecar_precision": 5,
                "finite_testability": 4,
                "theorem_pull": 5,
                "recursive_depth": 4,
                "risk_penalty": 0,
            },
        ),
        Pattern(
            code="R02",
            name="collar-to-bank lift",
            recursive_move=(
                "Use the AP-collar unit-height leak as the local base case and "
                "the expanded bank as the induction stress test."
            ),
            preserves=("strict-open witness mass", "boundary atoms", "height sidecar"),
            destroys_if_naive=("global owner support", "larger-bank collision class"),
            first_test=(
                "Attach each HYP-3405 AP-collar sidecar to the first HYP-3406 "
                "bank where the analogous sidecar fails or becomes insufficient."
            ),
            theorem_pull=(
                "A concrete finite lemma becomes a chamber theorem only when "
                "the AP/GW atoms and strict-open witnesses survive enlargement."
            ),
            challenged_assumption=(
                "Local AP-collar exactness is not automatically global; it must "
                "state which sidecar is inherited and which sidecar is newly lost."
            ),
            scores={
                "exit_purity": 5,
                "sidecar_precision": 4,
                "finite_testability": 5,
                "theorem_pull": 4,
                "recursive_depth": 4,
                "risk_penalty": 0,
            },
        ),
        Pattern(
            code="R03",
            name="height-then-owner escalation",
            recursive_move=(
                "Order sidecars by first failure: residue leaks to height/v2; "
                "height-persistent leaks escalate to owner support."
            ),
            preserves=("residue word", "v2/height word", "owner support"),
            destroys_if_naive=("sidecar order", "which repair is minimal"),
            first_test=(
                "Track the first failure of residue, residue+v2, "
                "residue+height, and residue+owner_support while extending the "
                "HYP-3406 scan."
            ),
            theorem_pull=(
                "The finite lemma should prove a sidecar priority chain, not a "
                "flat catalogue of packets."
            ),
            challenged_assumption=(
                "Height and owner are not interchangeable labels; the current "
                "data suggests a recursive escalation order."
            ),
            scores={
                "exit_purity": 4,
                "sidecar_precision": 5,
                "finite_testability": 4,
                "theorem_pull": 4,
                "recursive_depth": 4,
                "risk_penalty": 0,
            },
        ),
        Pattern(
            code="R04",
            name="finite chamber terminal router",
            recursive_move=(
                "Make every sidecar recursion terminate at an allowed theorem "
                "exit: AP/GW boundary, strict-open mass, q-witness, state-lift/H7, "
                "off-grid floor, exact-period/BDH exception, or named residual."
            ),
            preserves=("legal terminal exits", "finite residual labels"),
            destroys_if_naive=("termination proof", "exception ledger"),
            first_test=(
                "Add a terminal-exit column to each HYP-3406 mixed-fiber family "
                "before testing stronger sidecars."
            ),
            theorem_pull=(
                "This is the bridge from finite computation to a rigorous finite "
                "lemma: no recursive branch is allowed to end as an unnamed mix."
            ),
            challenged_assumption=(
                "A high-scoring sidecar is not enough; every recursion branch "
                "needs a terminal theorem label."
            ),
            scores={
                "exit_purity": 4,
                "sidecar_precision": 4,
                "finite_testability": 5,
                "theorem_pull": 4,
                "recursive_depth": 3,
                "risk_penalty": 0,
            },
        ),
        Pattern(
            code="R05",
            name="chiral child-deck recursion",
            recursive_move=(
                "If owner support is not yet enough, attach orientation/mirror "
                "child decks and recurse on the tail-tip edge witness."
            ),
            preserves=("owner side", "mirror/chiral side", "tail-tip orientation"),
            destroys_if_naive=("left/right boundary memory", "child gluing class"),
            first_test=(
                "For each owner leak, compare owner-deleted, tip-extended, and "
                "mirror-swapped child packets for theorem-exit purity."
            ),
            theorem_pull=(
                "This imports HYP-3124-style edge recursion only after the "
                "owner cut has named the boundary side."
            ),
            challenged_assumption=(
                "The next vertex set need not be runners or arcs; it can be "
                "child proof obligations generated by one mixed fiber."
            ),
            scores={
                "exit_purity": 4,
                "sidecar_precision": 4,
                "finite_testability": 3,
                "theorem_pull": 4,
                "recursive_depth": 5,
                "risk_penalty": 1,
            },
        ),
        Pattern(
            code="R06",
            name="local stability gate",
            recursive_move=(
                "Use Krasner/contact-root language only as a gate: local lifts "
                "are stable if the contact/root and owner packet is unchanged."
            ),
            preserves=("contact roots", "p-adic height", "owner support"),
            destroys_if_naive=("exit status under same-residue lift",),
            first_test=(
                "Record contact-root packets for same-residue strict moves "
                "12->26, 2->16, 13->27 and for tight 12->24."
            ),
            theorem_pull=(
                "A local-stability lemma is legal only when it preserves the "
                "same theorem predicate as the finite packet."
            ),
            challenged_assumption=(
                "P-adic closeness is not the recursive state; stability of the "
                "theorem-facing packet is."
            ),
            scores={
                "exit_purity": 4,
                "sidecar_precision": 4,
                "finite_testability": 4,
                "theorem_pull": 3,
                "recursive_depth": 4,
                "risk_penalty": 1,
            },
        ),
        Pattern(
            code="R07",
            name="quartic factor split",
            recursive_move=(
                "Translate quartic height/flex debt into two quadratic sign or "
                "owner channels before feeding it into the chamber theorem."
            ),
            preserves=("height/flex debt", "quadratic sign channels"),
            destroys_if_naive=("which factor carries the obstruction",),
            first_test=(
                "Apply the Sophie-Germain split to the HYP-3405 unit-height "
                "13->27 leak and both HYP-3406 owner-leak families."
            ),
            theorem_pull=(
                "This is an algebraic subroutine, not a scalar shortcut: it must "
                "return a labelled sidecar usable by the finite lemma."
            ),
            challenged_assumption=(
                "Named algebra is helpful only when it splits a live missing "
                "coordinate into exact packet fields."
            ),
            scores={
                "exit_purity": 3,
                "sidecar_precision": 4,
                "finite_testability": 4,
                "theorem_pull": 3,
                "recursive_depth": 3,
                "risk_penalty": 1,
            },
        ),
        Pattern(
            code="R08",
            name="mean-square exception ledger",
            recursive_move=(
                "Use BDH/Mertens-style averages only after exceptional fibers "
                "are named with owner, period, and smoothing labels."
            ),
            preserves=("exception labels", "period/denominator tail"),
            destroys_if_naive=("prime powers", "exact period", "endpoint owner"),
            first_test=(
                "For any future large-bank average, list the finite exceptional "
                "fibers before reporting a mean-square or entropy number."
            ),
            theorem_pull=(
                "Averaging can support a global tail only after the finite "
                "sidecar recursion has isolated the exceptions."
            ),
            challenged_assumption=(
                "Mean-square control is not a replacement for the finite lemma; "
                "it is a ledger for already-labelled residuals."
            ),
            scores={
                "exit_purity": 3,
                "sidecar_precision": 3,
                "finite_testability": 2,
                "theorem_pull": 3,
                "recursive_depth": 2,
                "risk_penalty": 2,
            },
        ),
        Pattern(
            code="R09",
            name="no-scalar-shadow firewall",
            recursive_move=(
                "Reject scalar shadows until they declare which finite packet "
                "coordinate they preserve or repair."
            ),
            preserves=("exact sidecar obligation", "proof predicate"),
            destroys_if_naive=("row-level witness", "owner/height route"),
            first_test=(
                "When a constant or analytic analogy appears, rewrite it as a "
                "packet column, terminal-exit label, or rejected guardrail."
            ),
            theorem_pull=(
                "This keeps recursive proof pressure on finite packet purity "
                "instead of on evocative scalar similarities."
            ),
            challenged_assumption=(
                "Exotic constants are not vertices of the proof tournament; "
                "translated proof obligations are."
            ),
            scores={
                "exit_purity": 3,
                "sidecar_precision": 2,
                "finite_testability": 3,
                "theorem_pull": 2,
                "recursive_depth": 1,
                "risk_penalty": 2,
            },
        ),
    ]


def ordered_patterns() -> list[Pattern]:
    return sorted(patterns(), key=lambda pattern: (-pattern.total, pattern.code))


def edge(a: Pattern, b: Pattern) -> bool:
    if a.total != b.total:
        return a.total > b.total
    return a.code < b.code


def directed_3cycles(vertices: list[Pattern]) -> int:
    count = 0
    for a, b, c in combinations(vertices, 3):
        ab, bc, ca = edge(a, b), edge(b, c), edge(c, a)
        ba, cb, ac = edge(b, a), edge(c, b), edge(a, c)
        if (ab and bc and ca) or (ba and cb and ac):
            count += 1
    return count


def hamiltonian_path_count(vertices: list[Pattern]) -> int:
    n = len(vertices)
    dp: dict[tuple[int, int], int] = {}
    for idx in range(n):
        dp[(1 << idx, idx)] = 1
    for mask in range(1 << n):
        for last in range(n):
            total = dp.get((mask, last), 0)
            if total == 0:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if edge(vertices[last], vertices[nxt]):
                    dp[(mask | (1 << nxt), nxt)] = (
                        dp.get((mask | (1 << nxt), nxt), 0) + total
                    )
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def print_pattern(pattern: Pattern) -> None:
    print(f"{pattern.code} total={pattern.total:2d}  {pattern.name}")
    print(f"  recursive_move: {pattern.recursive_move}")
    print(f"  preserves: {', '.join(pattern.preserves)}")
    print(f"  destroys_if_naive: {', '.join(pattern.destroys_if_naive)}")
    print(f"  first_test: {pattern.first_test}")
    print(f"  theorem_pull: {pattern.theorem_pull}")
    print(f"  challenged_assumption: {pattern.challenged_assumption}")


def main() -> None:
    ranked = ordered_patterns()

    print("HYP-3409 RECURSIVE SIDECAR PATTERN ATLAS")
    print("status=SYNTHESIS / executable recursion-pattern router; not an LRC14 proof")
    print("source=HYP-3405 + HYP-3406 + HYP-3407 + HYP-3408")
    print()

    print("ABSTRACT RECURSION")
    print("  R_0 = a coarse theorem-facing packet")
    print("  Q_i = a legal forgetful map applied to R_i")
    print("  if theorem_exit is pure on every Q_i fiber: Q_i is legal")
    print("  else: sigma_i = first destroyed sidecar visible in a mixed fiber")
    print("        R_{i+1} = R_i + sigma_i")
    print("        recurse until a terminal theorem exit or named residual appears")
    print()

    print("RANKED RECURSIVE PATTERNS")
    for pattern in ranked:
        print_pattern(pattern)
        print()

    print("CONCRETE TEST IDEAS")
    tests = [
        "Turn HYP-3405 AP-vs-13->27 and HYP-3406 owner leaks into one shared mixed-fiber API.",
        "Extend the HYP-3406 bank beyond (72,20) and record the first failure of residue+owner_support.",
        "Build the endpoint-owner Menger graph for petal 13->26 and petal 10->20 leak families.",
        "Add terminal-exit labels to every unresolved branch: AP/GW, strict-open mass, q-witness, H7/state lift, exact period/BDH exception, off-grid floor, or named residual.",
        "Run child-deck recursion on owner leaks: owner-deleted, tip-extended, and mirror-swapped packets.",
        "Apply the Krasner/contact-root gate and Sophie-Germain quartic split only after exact packet fields are declared.",
    ]
    for idx, test in enumerate(tests, start=1):
        print(f"  {idx}. {test}")
    print()

    print("TOURNAMENT ANALYSIS")
    print("  vertices=recursion operators and proof obligations")
    print("  rejected_vertices=runners, raw arcs, residues, constants, scalar analogies")
    print("  pairwise_observable=weighted theorem-exit purity plus sidecar precision")
    print("  switch=A -> B iff A has larger total score; ties use pattern code")
    print(f"  vertex_count={len(ranked)}")
    print(f"  score_hist={dict(sorted(Counter(p.total for p in ranked).items()))}")
    print(f"  directed_3cycles={directed_3cycles(ranked)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(ranked)}")
    print("  priority_path=" + " -> ".join(p.code for p in ranked))
    print()

    print("ASSUMPTION CHALLENGE")
    print("  considered alternate vertices:")
    print("    runners; gaps; fixed circle sections; section boundaries; wall-crossing events;")
    print("    residues; cover arcs; Fourier modes; matroid circuits; proof obligations.")
    print("  chosen vertices:")
    print("    recursion operators / proof obligations generated by mixed theorem-exit fibers.")
    print("  preserved predicate:")
    print("    theorem-exit purity: boundary-tight, strict-open, positive-Haar-open,")
    print("    unit-petal-named, q-witness, state-lift/H7, AP/GW, or named debt.")
    print("  information destroyed by this quotient:")
    print("    row order, raw runner identity, scalar motif values, and unlabelled analytic shadows.")
    print("  challenged assumption:")
    print("    LRC14 recursion should not be over runners.  It should be over legal")
    print("    forgetful maps whose lost coordinates are repaired exactly when they")
    print("    are needed by the theorem predicate.")


if __name__ == "__main__":
    main()
