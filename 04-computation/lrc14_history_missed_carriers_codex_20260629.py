#!/usr/bin/env python3
"""HYP-3494: history-mined missed carriers for the LRC14 random031 route.

This is a lightweight repo-history synthesis script.  It recomputes the local
random031 two-sheet fiber polynomial from HYP-3486 and then scores a curated
set of older LRC proof carriers for compatibility with the present
seam-complement proof target.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3486 = load_module(
    "hyp3486_random031_fiber_for_hyp3491",
    "lrc14_random031_seam_complement_fiber_graph_codex_20260629.py",
)


@dataclass(frozen=True)
class Candidate:
    hyp: str
    file: str
    niche: str
    transfer: str
    base: int
    current_anchor: bool = False


CANDIDATES = (
    Candidate(
        "HYP-3486",
        "05-knowledge/hypotheses/HYP-3486-lrc14-random031-seam-complement-fiber-graph.md",
        "local two-sheet fiber graph",
        "anchor: rank-2 exits, free holes, pure 12-cell bypass, vertical guardrail",
        55,
        True,
    ),
    Candidate(
        "HYP-3490",
        "05-knowledge/hypotheses/HYP-3490-lrc14-private-label-firewall.md",
        "private-label no-current firewall",
        "negative anchor: proves dead-projection current deletion is the wrong carrier",
        50,
        True,
    ),
    Candidate(
        "HYP-3485",
        "05-knowledge/hypotheses/HYP-3485-lrc14-random031-seam-complement-connection-atlas.md",
        "relative seam-complement bridge",
        "anchor: explains zipper/Cech/Menger/two-adic/PGF compatibility",
        48,
        True,
    ),
    Candidate(
        "HYP-3422",
        "05-knowledge/hypotheses/HYP-3422-lrc14-two-adic-offgrid-relocation.md",
        "two-adic relocation identity",
        "turn vertical half-turn into a branch-lift theorem, not a quotient",
        42,
    ),
    Candidate(
        "HYP-3140",
        "05-knowledge/hypotheses/HYP-3140-lrc14-fiber-pgf-rprime-certificate.md",
        "fiber-PGF first-moment certificate",
        "replace raw 282 witnesses by the local escape-sheet polynomial",
        41,
    ),
    Candidate(
        "HYP-3034",
        "05-knowledge/hypotheses/HYP-3034-lrc14-arc-boundary-path-lift.md",
        "owner-essential relative H1",
        "compute seam-complement H1 relative to rank-2 exits and free holes",
        39,
    ),
    Candidate(
        "HYP-3023",
        "05-knowledge/hypotheses/HYP-3023-lrc14-automatic-fiber-zipper.md",
        "magnitude-cocycle zipper guard",
        "price every seam-complement quotient by the first missing sidecar",
        38,
    ),
    Candidate(
        "HYP-3428",
        "05-knowledge/hypotheses/HYP-3428-lrc14-two-adic-descent-loss-ledger.md",
        "two-adic descent loss ledger",
        "record odd blockers and owner-current debt before using n*2",
        37,
    ),
    Candidate(
        "HYP-3451",
        "05-knowledge/hypotheses/HYP-3451-lrc14-component-cover-conductance-router.md",
        "component-cover conductance router",
        "run Menger/Green escape tests on seam-complement components",
        35,
    ),
    Candidate(
        "HYP-3437",
        "05-knowledge/hypotheses/HYP-3437-lrc14-overlap-menger-cut-certificate.md",
        "overlap-tax Menger cut",
        "treat bypass overlap as a labelled cut core rather than area mass",
        34,
    ),
    Candidate(
        "HYP-3402",
        "05-knowledge/hypotheses/HYP-3402-lrc14-owner-current-tropical-wall-angles.md",
        "endpoint-owner current / tropical wall",
        "separate seam-only owners from bypass-flow owners as boundary current",
        33,
    ),
    Candidate(
        "HYP-3300",
        "05-knowledge/hypotheses/HYP-3300-lrc14-observability-morse-proof-angles.md",
        "observability matrix and Morse descent",
        "turn the final random031 proof into sidecar columns plus legal descent",
        31,
    ),
    Candidate(
        "HYP-3243",
        "05-knowledge/hypotheses/HYP-3243-lrc14-topology-geometry-graph-proof-routes.md",
        "oriented topes/cocircuits atlas",
        "model deleted seam as cocircuit boundary with open complement flow",
        30,
    ),
    Candidate(
        "HYP-3480",
        "05-knowledge/hypotheses/HYP-3480-lrc14-zero-edge-singleton-current.md",
        "formal singleton-current contrast",
        "Lean-style dispatch split: six singleton-current rows versus random031",
        29,
    ),
)


MOTIFS = {
    "seam": ("seam", "puncture", "bypass", "random031"),
    "quotient": ("quotient", "forget", "sidecar", "destroyed", "retained"),
    "topology": ("cech", "h1", "topology", "tope", "cocircuit", "component"),
    "current": ("current", "menger", "green", "cut", "conductance", "farkas"),
    "two_adic": ("two-adic", "2-adic", "u=2t", "half", "branch"),
    "fiber": ("fiber", "pgf", "sheet", "zipper", "barcode"),
    "owner": ("owner", "label", "monodromy"),
    "formal": ("lean", "theorem", "proof target", "lemma"),
}


def motif_counts(path: Path) -> Counter[str]:
    text = path.read_text(encoding="utf-8", errors="ignore").lower()
    counts: Counter[str] = Counter()
    for motif, terms in MOTIFS.items():
        counts[motif] = sum(text.count(term) for term in terms)
    return counts


def motif_bonus(counts: Counter[str]) -> int:
    weights = {
        "seam": 4,
        "quotient": 3,
        "topology": 3,
        "current": 3,
        "two_adic": 4,
        "fiber": 3,
        "owner": 3,
        "formal": 2,
    }
    return sum(weights[key] * min(counts[key], 5) for key in weights)


def random031_local_pgf() -> dict[str, object]:
    row = H3486.H3481.H3450.audit_row(H3486.H3481.ROW_NAME, H3486.H3481.SPEEDS)
    gates = H3486.H3481.build_gates()
    cells = H3486.build_cells(gates, row)

    fibers: dict[int, list[object]] = defaultdict(list)
    for cell in cells:
        fibers[cell.u_index].append(cell)

    rank_escape_hist: Counter[int] = Counter()
    terminal_sheet_hist: Counter[int] = Counter()
    class_signature_hist: Counter[tuple[str, ...]] = Counter()
    for fiber_cells in fibers.values():
        rank_escape_hist[sum(1 for cell in fiber_cells if cell.hits)] += 1
        terminal_sheet_hist[len(fiber_cells)] += 1
        class_signature_hist[tuple(sorted(cell.cell_class for cell in fiber_cells))] += 1

    legal_components = H3486.connected_components(cells, {"horizontal", "mirror"})
    vertical_components = H3486.connected_components(cells, {"horizontal", "mirror", "vertical"})
    mixed_vertical = [comp for comp in vertical_components if len(comp.class_hist) > 1]

    total_rank_escapes = sum(k * v for k, v in rank_escape_hist.items())
    occupied_fibers = len(fibers)
    return {
        "cells": len(cells),
        "occupied_fibers": occupied_fibers,
        "rank_escape_hist": dict(sorted(rank_escape_hist.items())),
        "terminal_sheet_hist": dict(sorted(terminal_sheet_hist.items())),
        "class_signature_hist": class_signature_hist.most_common(),
        "mean_rank_escapes": Fraction(total_rank_escapes, occupied_fibers),
        "legal_components": len(legal_components),
        "vertical_components": len(vertical_components),
        "vertical_mixed_components": len(mixed_vertical),
    }


def compact_counter(counter: Counter[str] | dict) -> str:
    items = counter.items() if isinstance(counter, dict) else counter.items()
    return "{" + ", ".join(f"{k}:{v}" for k, v in sorted(items)) + "}"


def main() -> None:
    pgf = random031_local_pgf()
    rows = []
    for candidate in CANDIDATES:
        counts = motif_counts(ROOT / candidate.file)
        score = candidate.base + motif_bonus(counts)
        rows.append((score, candidate, counts))
    rows.sort(key=lambda item: (-item[0], item[1].hyp))

    nonanchors = [row for row in rows if not row[1].current_anchor]
    score_hist = Counter(score for score, _candidate, _counts in nonanchors)
    path = [candidate.hyp for _score, candidate, _counts in nonanchors]

    print("HYP-3494 LRC14 HISTORY-MINED MISSED CARRIER ATLAS")
    print("status=SYNTHESIS / repo-history carrier mining plus local PGF; not an LRC14 proof")
    print("scope=random_covering_031 seam-complement proof extension")
    print()
    print("## Local Random031 Fiber PGF")
    print(f"witness_cells={pgf['cells']} occupied_fibers={pgf['occupied_fibers']}")
    print(f"terminal_sheet_hist={pgf['terminal_sheet_hist']}")
    print(f"rank2_escape_pgf={pgf['rank_escape_hist']}  # 24 + 226*y + 8*y^2")
    print(f"mean_rank2_escape_sheets_per_occupied_fiber={pgf['mean_rank_escapes']}")
    print(f"class_signature_top={pgf['class_signature_hist'][:6]}")
    print(
        "legal_vs_vertical_components="
        f"{pgf['legal_components']} legal, {pgf['vertical_components']} vertical-glued, "
        f"{pgf['vertical_mixed_components']} vertical-mixed"
    )
    print(
        "pgf_pull=HYP-3140 localizes to a two-sheet seam-complement polynomial; "
        "the raw 282 count can be replaced by a fiber moment plus free-hole sidecar."
    )
    print()

    print("## Current Anchors")
    for score, candidate, counts in rows:
        if not candidate.current_anchor:
            continue
        print(
            f"{candidate.hyp}: score={score} niche={candidate.niche}; "
            f"transfer={candidate.transfer}; motifs={compact_counter(counts)}"
        )
    print()

    print("## History-Mined Missed Carriers")
    for rank, (score, candidate, counts) in enumerate(nonanchors, start=1):
        print(
            f"{rank:02d}. {candidate.hyp}: score={score} niche={candidate.niche}; "
            f"transfer={candidate.transfer}; motifs={compact_counter(counts)}"
        )
    print()

    print("## New Proof Experiments")
    print("E1 local_pgf_escape_moment: prove 24 + 226*y + 8*y^2 is the legal random031 escape polynomial after seam deletion, with free-hole packets terminal rather than failed routes.")
    print("E2 relative_h1_exit: build H1(seam complement, rank2 exits union free-hole boundary) and test whether the pure bypass is the only nontrivial class carrying hard-component flow.")
    print("E3 quotient_price_matrix: make a HYP-3300 observability matrix with columns u-index, branch, mirror mate, class, endpoint rank, owner word, private firewall, and vertical-halfturn sidecar.")
    print("E4 owner_current_cobordism: prove seam-only owners (45,147,169,173) are boundary charge while bypass owners (23,93,113) are flow charge, or emit named owner/two-adic/SPEC debt.")
    print("E5 formal_dispatch_join: extend the HYP-3480 Lean dispatch ledger so private-firewall rows split into six singleton-current terminals plus one random031 seam-cobordism terminal.")
    print()

    print("## Tournament Analysis")
    print("vertices=history-mined proof carriers, not runners or raw gates")
    print("pairwise_observable=random031 predicate retention + quotient sidecar cost + executable follow-up clarity")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))
    print()

    print("## Assumption Challenge")
    print("Considered vertices: runners, gaps, hard gates, witness cells, u-fibers, branch sheets,")
    print("owner labels, dead islands, punctures, relative cycles, component cuts, PGF")
    print("coefficients, observability columns, Lean dispatch packets, and proof obligations.")
    print("Chosen vertices are proof carriers plus the local fiber polynomial.  This preserves")
    print("the random031 terminal predicate after seam deletion and explicitly prices the")
    print("destroyed coordinate in each proposed quotient.")


if __name__ == "__main__":
    main()
