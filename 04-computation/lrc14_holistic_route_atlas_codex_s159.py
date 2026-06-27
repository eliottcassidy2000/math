#!/usr/bin/env python3
"""S159: holistic route atlas for the LRC14 proof history.

This script mines the local research workspace for LRC/LRC14 theorem,
hypothesis, reflection, and forum artifacts.  It is a synthesis aid: it counts
which proof carriers recur, records how the route taxonomy changed over time,
and builds a Tournament Analysis on proof carriers rather than runners.

The pairwise observable for the tournament is a retention vector:

    strict-predicate retention,
    exact arithmetic retention,
    owner/packet label retention,
    adaptability to unbounded families,
    dual-certificate strength,
    executable auditability,
    anti-scalarization guard.

The switch orients carrier A -> B when A wins more coordinates; exact ties use
the listed Hamiltonian path.  This deliberately makes the theorem-use relation
explicit before any quotient is trusted.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import argparse
import re


REPO = Path(__file__).resolve().parents[1]
RESULT = REPO / "05-knowledge" / "results" / "lrc14_holistic_route_atlas_codex_s159.out"

SCAN_DIRS = [
    REPO / "01-canon" / "theorems",
    REPO / "05-knowledge" / "hypotheses",
    REPO / "07-reflections",
    REPO / "poke-forum" / "posts",
]

ROUTES: dict[str, tuple[str, ...]] = {
    "qdiv_small_denominator": ("qdiv", "q-witness", "small-denominator", "multiple of every", "covering set"),
    "ap_gw_boundary": ("ap/gw", "goddyn-wong", "boundary atom", "boundary-only", "12->24"),
    "endpoint_core": ("endpoint", "protection core", "owner-protection", "boundary bridge", "pinch", "winding"),
    "observer_source": ("observer-source", "source-spectrum", "source cone", "source walk", "marked tournament"),
    "raw_tournament_guardrail": ("raw tournament", "isomorphism class", "apex-pressure", "achievable", "a000568"),
    "farey_exact_m": ("farey", "exact m", "q_threshold", "3/41", "17/41", "continued fraction"),
    "singular_series": ("singular series", "sinc", "relation lattice", "theta form", "riesz"),
    "haar_baire": ("haar", "baire", "regular-open", "strict safe", "safe_mu"),
    "c27_unital_k33": ("c27", "unital", "k33", "kuratowski", "state-lift", "d9"),
    "fixed_margin_packet": ("fixed-margin", "labelled packet", "packet-preserving", "johnson", "families"),
    "moon_core_apex": ("moon core", "apex-majority", "14z", "few-apex", "aperture", "comb"),
    "boundary_moment_gk8": ("boundary-moment", "gk8", "l_y", "missed-depth", "sector moment"),
    "lift_packet": ("lift packet", "u=14t", "fourteen lifts", "q/r split"),
    "nork_pinch": ("nork", "pinch-template", "no open residual", "f6"),
    "moment_dual": ("moment dual", "multiplicity", "danger-count", "factorial", "polynomial dual"),
    "twist_ladder": ("twist-ladder", "rational twist", "blocker hypergraph", "denominator wall"),
    "fourier_toeplitz": ("fourier-toeplitz", "toeplitz", "psd", "trigonometric square"),
    "wide_decorrelation": ("wide", "decorrelation", "single-far", "doublet"),
}

PHASES = [
    ("proved endpoint/qdiv base", "THM", 358, 398),
    ("signed/tournament witness experiments", "THM", 407, 491),
    ("singular-series and density side", "THM", 501, 522),
    ("gap-side q-witness covering reduction", "THM", 523, 548),
    ("sector/wide/finite-core machinery", "THM", 555, 566),
    ("apex and tournament-state closure", "THM", 568, 572),
    ("induction/tournament/AP-GW setup", "HYP", 2900, 2924),
    ("Farey/C27/unital/K33 carrier phase", "HYP", 2930, 2947),
    ("Haar/source/gauntlet bridge", "HYP", 2948, 2955),
    ("family-sporadic labelled packets", "HYP", 2956, 2963),
    ("Moon-core, lifts, boundary-moment", "HYP", 2964, 2969),
    ("dual-certificate convergence", "HYP", 2970, 2974),
]


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def iter_markdown_files() -> list[Path]:
    files: list[Path] = []
    for root in SCAN_DIRS:
        if root.exists():
            files.extend(sorted(root.rglob("*.md")))
    return files


def parse_id(path: Path) -> tuple[str, int] | None:
    m = re.search(r"\b(THM|HYP)-(\d+)\b", path.name)
    if m:
        return (m.group(1), int(m.group(2)))
    return None


@dataclass(frozen=True)
class FileHit:
    path: Path
    kind: str
    number: int | None
    routes: frozenset[str]


def scan_corpus() -> tuple[list[FileHit], Counter[str], Counter[tuple[str, str]], Counter[str]]:
    route_counts: Counter[str] = Counter()
    route_pairs: Counter[tuple[str, str]] = Counter()
    phase_counts: Counter[str] = Counter()
    hits: list[FileHit] = []

    for path in iter_markdown_files():
        text = read_text(path)
        lower = text.lower()
        routes = {
            name
            for name, patterns in ROUTES.items()
            if any(pattern in lower for pattern in patterns)
        }
        if not routes and "lrc" not in lower and "lonely runner" not in lower:
            continue

        ident = parse_id(path)
        kind = ident[0] if ident else "DOC"
        number = ident[1] if ident else None
        hits.append(FileHit(path, kind, number, frozenset(routes)))

        for route in routes:
            route_counts[route] += 1
        for a, b in combinations(sorted(routes), 2):
            route_pairs[(a, b)] += 1
        if ident:
            for label, p_kind, lo, hi in PHASES:
                if kind == p_kind and lo <= number <= hi:
                    phase_counts[label] += 1
                    break

    return hits, route_counts, route_pairs, phase_counts


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: tuple[int, int, int, int, int, int, int]
    preserves: str
    destroys: str


CARRIERS = [
    Carrier("source_spectrum_pullback", (5, 4, 5, 4, 2, 3, 5), "observer-source predicate through Farey/Haar/packet labels", "raw row identity and fine wall order"),
    Carrier("endpoint_winding_dual", (5, 5, 5, 3, 5, 4, 5), "open-cover failure as positive-winding endpoint circulation", "bulk Haar geometry away from active SCC"),
    Carrier("fixed_margin_labelled_packet", (5, 4, 5, 4, 3, 4, 5), "family/sporadic status with source and owner labels", "exact representative inside a packet fiber"),
    Carrier("boundary_gap_lift_packet", (5, 5, 5, 3, 3, 5, 5), "strict safe intervals and Q/R lift status", "global source spectrum outside the lift chart"),
    Carrier("twist_ladder_blocker", (5, 4, 3, 5, 4, 5, 4), "finite rational witness or adaptive blocker cover", "continuous intervals between selected twists"),
    Carrier("danger_count_moment_dual", (4, 4, 2, 4, 5, 5, 4), "count-cover obstruction through exact moment inequalities", "endpoint location and packet ownership"),
    Carrier("fourier_toeplitz_psd", (4, 4, 3, 4, 5, 2, 4), "phase-sensitive necessary PSD conditions for nonnegative cover excess", "direct interval ownership until eigenvectors are decoded"),
    Carrier("c27_unital_k33_state", (4, 4, 5, 2, 2, 4, 5), "C27 owner transfer, unit/nonunit fork, K33 state-lift debt", "off-shell magnitude unless Farey data is attached"),
    Carrier("exact_M_Farey", (5, 5, 3, 3, 1, 4, 4), "closed gap value and binding denominator branch", "why the branch is structurally forced"),
    Carrier("haar_baire_boundary", (4, 5, 4, 3, 1, 4, 4), "open-vs-boundary safe set and endpoint debt", "C27/K33 ownership unless reattached"),
    Carrier("qdiv_gate", (5, 5, 2, 4, 1, 5, 4), "small-denominator witness or covering-core necessity", "all structure after the first covering gate"),
    Carrier("singular_series_relation_lattice", (3, 3, 2, 5, 3, 3, 3), "density-side additive relation pressure", "closed-threshold gap and endpoint debt"),
    Carrier("raw_tournament_class", (2, 2, 1, 2, 1, 5, 1), "chosen pairwise ordering only", "qdiv, magnitude, endpoint owners, open mass"),
    Carrier("raw_scalar_mass", (2, 3, 0, 2, 1, 4, 0), "one numeric projection", "almost every theorem-relevant label"),
]


def carrier_edge(a: Carrier, b: Carrier, tie_order: dict[str, int]) -> bool:
    wins = 0
    losses = 0
    for x, y in zip(a.vector, b.vector):
        if x > y:
            wins += 1
        elif x < y:
            losses += 1
    if wins != losses:
        return wins > losses
    return tie_order[a.name] < tie_order[b.name]


def tournament_edges() -> dict[str, set[str]]:
    order = {carrier.name: idx for idx, carrier in enumerate(CARRIERS)}
    edges: dict[str, set[str]] = {carrier.name: set() for carrier in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        if carrier_edge(a, b, order):
            edges[a.name].add(b.name)
        else:
            edges[b.name].add(a.name)
    return edges


def count_directed_triangles(edges: dict[str, set[str]]) -> int:
    count = 0
    names = list(edges)
    for a, b, c in combinations(names, 3):
        sub = [(a, b), (a, c), (b, c)]
        score = Counter()
        for x, y in sub:
            if y in edges[x]:
                score[x] += 1
            else:
                score[y] += 1
        if sorted(score.values()) == [1, 1, 1]:
            count += 1
    return count


def strongly_connected_sizes(edges: dict[str, set[str]]) -> list[int]:
    names = list(edges)
    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    low: dict[str, int] = {}
    sizes: list[int] = []

    def visit(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)

        for w in edges[v]:
            if w not in indices:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])

        if low[v] == indices[v]:
            size = 0
            while True:
                w = stack.pop()
                on_stack.remove(w)
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for name in names:
        if name not in indices:
            visit(name)
    return sorted(sizes, reverse=True)


def count_hamiltonian_paths(edges: dict[str, set[str]]) -> int:
    names = list(edges)
    n = len(names)
    idx = {name: i for i, name in enumerate(names)}
    outmask = [0] * n
    for name, outs in edges.items():
        i = idx[name]
        for out in outs:
            outmask[i] |= 1 << idx[out]

    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp.get((mask, last), 0)
            if not value:
                continue
            available = outmask[last] & ~mask
            nxt = available
            while nxt:
                bit = nxt & -nxt
                j = bit.bit_length() - 1
                dp[(mask | bit, j)] = dp.get((mask | bit, j), 0) + value
                nxt -= bit
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def tournament_summary() -> list[str]:
    edges = tournament_edges()
    scores = {name: len(outs) for name, outs in edges.items()}
    score_hist = Counter(scores.values())
    path = sorted(scores, key=lambda name: (-scores[name], [c.name for c in CARRIERS].index(name)))

    lines = [
        "Tournament Analysis",
        "-------------------",
        "vertices: proof carriers, not runners",
        "pairwise observable: predicate/exactness/labels/adaptability/dual/computation/anti-scalar retention",
        f"score_hist: {dict(sorted(score_hist.items()))}",
        f"directed_3cycles: {count_directed_triangles(edges)}",
        f"SCC_sizes: {strongly_connected_sizes(edges)}",
        f"Hamiltonian_path_count: {count_hamiltonian_paths(edges)}",
        "tie/guiding Hamiltonian path:",
        "  " + " > ".join(path),
        "",
        "carrier retention notes:",
    ]
    for carrier in CARRIERS:
        lines.append(f"- {carrier.name}: preserves {carrier.preserves}; destroys {carrier.destroys}.")
    return lines


def render() -> str:
    hits, route_counts, route_pairs, phase_counts = scan_corpus()
    theorem_hits = [h for h in hits if h.kind == "THM"]
    hyp_hits = [h for h in hits if h.kind == "HYP"]
    forum_hits = [h for h in hits if "poke-forum" in str(h.path)]

    lines: list[str] = []
    lines.append("LRC14 holistic route atlas (codex S159)")
    lines.append("=======================================")
    lines.append("")
    lines.append("Corpus")
    lines.append("------")
    lines.append(f"markdown artifacts scanned with LRC signal: {len(hits)}")
    lines.append(f"canon theorem artifacts: {len(theorem_hits)}")
    lines.append(f"hypothesis artifacts: {len(hyp_hits)}")
    lines.append(f"forum post/comment artifacts: {len(forum_hits)}")
    lines.append("")
    lines.append("Route frequency")
    lines.append("---------------")
    for route, count in route_counts.most_common():
        lines.append(f"{route:28s} {count:4d}")
    lines.append("")
    lines.append("Historical phase density")
    lines.append("------------------------")
    for label, _kind, _lo, _hi in PHASES:
        lines.append(f"{label:42s} {phase_counts[label]:4d}")
    lines.append("")
    lines.append("Top route co-occurrences")
    lines.append("------------------------")
    for (a, b), count in route_pairs.most_common(25):
        lines.append(f"{a} + {b}: {count}")
    lines.append("")
    lines.extend(tournament_summary())
    lines.append("")
    lines.append("Synthesis readout")
    lines.append("-----------------")
    lines.append("1. The repo's earliest durable LRC move was endpoint/qdiv arithmetic, not enumeration.")
    lines.append("2. Raw tournament classes repeatedly served as front filters and falsifiers of lossy quotients.")
    lines.append("3. AP/GW evolved from examples into the unique tested closed-boundary equality stratum.")
    lines.append("4. The qdiv>14 branch became the true strict-counterexample branch after THM-523.")
    lines.append("5. THM-566 refuted fixed bounded-denominator closure and forced adaptive denominators.")
    lines.append("6. The modern object is a labelled source packet: qdiv, exact M/Farey, Haar boundary, C27/K33 owners, lift/moment data.")
    lines.append("7. Current dual routes are converging on the same negation: no open cover can survive all endpoint, moment, twist, and PSD certificates except AP/GW boundary atoms or a state-lift debt.")
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the stored result file")
    args = parser.parse_args()

    text = render()
    print(text, end="")
    if args.write:
        RESULT.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
