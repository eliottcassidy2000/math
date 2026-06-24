#!/usr/bin/env python3
"""S166: Haar zipper cocycle synthesis for the LRC14 proof stack.

This is a proof-interface computation, not a proof of LRC14.

The recent Haar-product scouts identified the local checkerboard

    [[+1, -1],
     [-1, +1]]

as both the two-dimensional Haar product on dyadic children and the elementary
fixed-margin switch in the tournament tiling model.  This script makes the
next abstraction explicit: the mixed Haar coefficient is a local cocycle on a
fixed row/column-margin square.  A handoff is safe only if it preserves that
cocycle, cancels it by a signed discrepancy argument, or routes it to a named
certificate/state-lift label.

Tournament Analysis uses proof carriers as vertices, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations, permutations
from math import gcd
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


HAAR_SWITCH = ((1, -1), (-1, 1))


def row_sums(table: tuple[tuple[int, int], tuple[int, int]]) -> tuple[int, int]:
    return tuple(sum(row) for row in table)


def col_sums(table: tuple[tuple[int, int], tuple[int, int]]) -> tuple[int, int]:
    return tuple(table[0][j] + table[1][j] for j in range(2))


def total(table: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return sum(sum(row) for row in table)


def zipper_cocycle(table: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return sum(table[i][j] * HAAR_SWITCH[i][j] for i in range(2) for j in range(2))


def all_tables_with_total(n: int) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    out = []
    for a in range(n + 1):
        for b in range(n - a + 1):
            for c in range(n - a - b + 1):
                d = n - a - b - c
                out.append(((a, b), (c, d)))
    return out


def margin_key(table: tuple[tuple[int, int], tuple[int, int]]) -> tuple[tuple[int, int], tuple[int, int], int]:
    return (row_sums(table), col_sums(table), total(table))


def fixed_margin_audit(max_total: int = 10) -> dict[str, object]:
    fibers: dict[tuple[tuple[int, int], tuple[int, int], int], list[tuple[tuple[int, int], tuple[int, int]]]] = defaultdict(list)
    augmented: dict[tuple[tuple[int, int], tuple[int, int], int, int], tuple[tuple[int, int], tuple[int, int]]] = {}
    duplicate_augmented = 0
    for n in range(max_total + 1):
        for table in all_tables_with_total(n):
            fibers[margin_key(table)].append(table)
            key = (*margin_key(table), zipper_cocycle(table))
            if key in augmented and augmented[key] != table:
                duplicate_augmented += 1
            augmented[key] = table

    nontrivial = {key: vals for key, vals in fibers.items() if len(vals) > 1}
    step_gcd = 0
    max_span = 0
    example = None
    for key, vals in nontrivial.items():
        zetas = sorted(zipper_cocycle(v) for v in vals)
        diffs = [b - a for a, b in zip(zetas, zetas[1:]) if b != a]
        for diff in diffs:
            step_gcd = gcd(step_gcd, abs(diff))
        span = zetas[-1] - zetas[0]
        if span > max_span:
            max_span = span
            example = (key, vals, zetas)

    return {
        "max_total": max_total,
        "tables": sum(len(all_tables_with_total(n)) for n in range(max_total + 1)),
        "fibers": len(fibers),
        "nontrivial_fibers": len(nontrivial),
        "augmented_keys": len(augmented),
        "duplicate_augmented_keys": duplicate_augmented,
        "zipper_step_gcd": step_gcd,
        "max_zipper_span": max_span,
        "max_span_example": example,
    }


@dataclass(frozen=True, order=True)
class Interval:
    level: int
    index: int

    @property
    def denom(self) -> int:
        return 1 << self.level

    def contains(self, other: "Interval") -> bool:
        return (
            self.index * other.denom <= other.index * self.denom
            and (other.index + 1) * self.denom <= (self.index + 1) * other.denom
        )

    def disjoint(self, other: "Interval") -> bool:
        return (
            (self.index + 1) * other.denom <= other.index * self.denom
            or (other.index + 1) * self.denom <= self.index * other.denom
        )

    def half_sign_on(self, inner: "Interval") -> int:
        lhs = (2 * inner.index + 1) * self.denom
        rhs = (2 * self.index + 1) * inner.denom
        return 1 if lhs < rhs else -1


@dataclass(frozen=True)
class Rect:
    x: Interval
    y: Interval


@dataclass(frozen=True)
class Product1D:
    kind: str
    relation: str
    sign: int


def one_d_product(a: Interval, b: Interval) -> Product1D:
    if a == b:
        return Product1D("indicator", "equal", 1)
    if a.disjoint(b):
        return Product1D("zero", "disjoint", 0)
    if a.contains(b):
        return Product1D("haar", "A_contains_B", a.half_sign_on(b))
    if b.contains(a):
        return Product1D("haar", "B_contains_A", b.half_sign_on(a))
    raise AssertionError("dyadic intervals are nested, equal, or disjoint")


def product_class(px: Product1D, py: Product1D) -> str:
    if px.kind == "zero" or py.kind == "zero":
        return "orthogonal_zero"
    if px.kind == "indicator" and py.kind == "indicator":
        return "same_tile_boundary_atom"
    if px.kind == "indicator" and py.kind == "haar":
        return "vertical_owner_strip"
    if px.kind == "haar" and py.kind == "indicator":
        return "horizontal_owner_strip"
    rels = {px.relation, py.relation}
    if rels == {"A_contains_B"} or rels == {"B_contains_A"}:
        return "nested_refinement"
    if rels == {"A_contains_B", "B_contains_A"}:
        return "cross_zipper_handoff"
    raise AssertionError((px, py))


def intervals(max_depth: int) -> list[Interval]:
    out = []
    for level in range(max_depth + 1):
        for index in range(1 << level):
            out.append(Interval(level, index))
    return out


def rects(max_depth: int) -> list[Rect]:
    ints = intervals(max_depth)
    return [Rect(x, y) for x in ints for y in ints]


def dyadic_product_census(max_depth: int = 4) -> tuple[Counter[str], dict[str, Counter[int]]]:
    class_count: Counter[str] = Counter()
    sign_count: dict[str, Counter[int]] = defaultdict(Counter)
    rs = rects(max_depth)
    for a in rs:
        for b in rs:
            px = one_d_product(a.x, b.x)
            py = one_d_product(a.y, b.y)
            cls = product_class(px, py)
            sign = px.sign * py.sign
            class_count[cls] += 1
            sign_count[cls][sign] += 1
    return class_count, sign_count


CONCEPTS = {
    "haar_product": ["haar product", "haar tile", "mixed haar", "checkerboard"],
    "tournament_tiling": ["tournament tiling", "fixed-margin", "hamiltonian path", "fixed-path"],
    "discrepancy": ["discrepancy", "delta", "component count", "colored resonance"],
    "zipper_handoff": ["zipper", "handoff", "arrow", "certificate handoff"],
    "tope_cocircuit": ["tope", "cocircuit", "endpoint"],
    "exposure_kernel": ["exposure", "unexposed", "source kernel", "no-hidden"],
    "smoothing": ["smoothing", "kaczynski", "kernel", "approach class"],
    "state_lift": ["state-lift", "state lift", "thm-572", "h=7"],
}


DOCS = [
    "05-knowledge/hypotheses/HYP-2989-lrc14-haar-product-discrepancy-tournament-tiling.md",
    "05-knowledge/hypotheses/HYP-2987-lrc14-certificate-handoff-atlas.md",
    "05-knowledge/hypotheses/HYP-2988-lrc14-exposure-poset-proof-pass.md",
    "05-knowledge/hypotheses/HYP-2986-lrc14-tope-wall-cocircuit.md",
    "05-knowledge/hypotheses/HYP-2985-lrc14-admissible-smoothing-dispatcher.md",
    "05-knowledge/hypotheses/HYP-2595-lrc14-colored-resonance-discrepancy.md",
    "05-knowledge/hypotheses/HYP-2594-lrc14-colored-discrepancy-bound.md",
    "05-knowledge/hypotheses/HYP-2450-coefficient-tiling-prime-bridge.md",
    "05-knowledge/hypotheses/HYP-2736-lrc14-torus-line-discrepancy-integer-grid.md",
]


def concept_census() -> tuple[Counter[str], Counter[tuple[str, str]]]:
    hits: Counter[str] = Counter()
    co: Counter[tuple[str, str]] = Counter()
    for rel in DOCS:
        text = (REPO / rel).read_text(encoding="utf-8", errors="ignore").lower()
        present = []
        for name, needles in CONCEPTS.items():
            count = sum(text.count(needle) for needle in needles)
            if count:
                hits[name] += count
                present.append(name)
        for a, b in combinations(sorted(present), 2):
            co[(a, b)] += 1
    return hits, co


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: tuple[int, ...]
    role: str


CARRIERS = [
    Carrier("haar_zipper_cocycle", (5, 5, 5, 5, 4, 5, 5, 4), "local mixed-sign coordinate; complete for 2x2 fixed-margin fibers"),
    Carrier("certificate_handoff_atlas", (5, 4, 4, 5, 5, 4, 5, 5), "global zipper arrows; specifies allowed quotient handoffs"),
    Carrier("exposure_kernel_audit", (5, 4, 4, 4, 5, 5, 4, 5), "tests whether any labelled packet survives all local exits"),
    Carrier("tope_cocircuit_wall", (4, 5, 4, 4, 4, 4, 5, 4), "endpoint-cell geometry; separates open topes from boundary atoms"),
    Carrier("color_resonance_discrepancy", (4, 3, 5, 3, 4, 4, 3, 4), "global CRT/Fourier compatibility for mixed modes"),
    Carrier("admissible_smoothing_clock", (4, 3, 4, 4, 4, 3, 3, 5), "analytic kernel changes with retained approach labels"),
    Carrier("fixed_margin_tiling_shadow", (3, 5, 2, 2, 3, 2, 2, 3), "row/column fiber and tournament tile quotient"),
    Carrier("raw_component_count_K", (2, 1, 1, 1, 1, 1, 1, 1), "safe but lossy continuous boundary count"),
]

TIE_PATH = [carrier.name for carrier in CARRIERS]


def beats(a: Carrier, b: Carrier) -> bool:
    aw = sum(1 for x, y in zip(a.vector, b.vector) if x > y)
    bw = sum(1 for x, y in zip(a.vector, b.vector) if x < y)
    if aw != bw:
        return aw > bw
    return TIE_PATH.index(a.name) < TIE_PATH.index(b.name)


def tournament_adjacency() -> dict[str, set[str]]:
    adj = {carrier.name: set() for carrier in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        if beats(a, b):
            adj[a.name].add(b.name)
        else:
            adj[b.name].add(a.name)
    return adj


def directed_3cycles(adj: dict[str, set[str]]) -> int:
    count = 0
    for a, b, c in combinations(adj, 3):
        if (b in adj[a] and c in adj[b] and a in adj[c]) or (c in adj[a] and b in adj[c] and a in adj[b]):
            count += 1
    return count


def scc_sizes(adj: dict[str, set[str]]) -> list[int]:
    radj = {name: set() for name in adj}
    for a, outs in adj.items():
        for b in outs:
            radj[b].add(a)

    def reach(graph: dict[str, set[str]], start: str) -> set[str]:
        seen = {start}
        q = deque([start])
        while q:
            cur = q.popleft()
            for nxt in graph[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    q.append(nxt)
        return seen

    remaining = set(adj)
    sizes = []
    while remaining:
        start = next(iter(remaining))
        comp = reach(adj, start) & reach(radj, start)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: dict[str, set[str]]) -> list[tuple[str, ...]]:
    paths = []
    names = list(adj)
    for order in permutations(names):
        if all(order[i + 1] in adj[order[i]] for i in range(len(order) - 1)):
            paths.append(order)
    return paths


def main() -> None:
    audit = fixed_margin_audit()
    class_count, sign_count = dyadic_product_census()
    concept_hits, concept_co = concept_census()
    adj = tournament_adjacency()
    paths = hamiltonian_paths(adj)

    print("S166 Haar zipper cocycle / discrepancy handoff synthesis")
    print("=" * 72)
    print()
    print("[1] Local 2x2 fixed-margin zipper cocycle")
    print(f"  switch / Haar product matrix = {HAAR_SWITCH}")
    print("  zeta(T)=sum_ij (-1)^(i+j) T_ij")
    print(f"  audited nonnegative 2x2 tables with total <= {audit['max_total']}")
    print(f"  tables={audit['tables']} margin_fibers={audit['fibers']} nontrivial_fibers={audit['nontrivial_fibers']}")
    print(f"  augmented keys (margins,zeta)={audit['augmented_keys']}")
    print(f"  duplicate augmented keys={audit['duplicate_augmented_keys']}")
    print(f"  zipper step gcd inside nontrivial fibers={audit['zipper_step_gcd']}")
    print(f"  max zipper span inside one margin fiber={audit['max_zipper_span']}")
    key, vals, zetas = audit["max_span_example"]  # type: ignore[misc]
    print(f"  max-span margin key={key}, zeta_values={zetas}")
    print("  conclusion: margins alone collide, but margins plus zeta are a complete local coordinate.")
    print()

    print("[2] Dyadic Haar rectangle product census")
    total_pairs = sum(class_count.values())
    print(f"  depth<=4 rectangles={len(rects(4))}, ordered_products={total_pairs}")
    order = [
        "same_tile_boundary_atom",
        "vertical_owner_strip",
        "horizontal_owner_strip",
        "cross_zipper_handoff",
        "nested_refinement",
        "orthogonal_zero",
    ]
    for cls in order:
        signs = ", ".join(f"{s}:{n}" for s, n in sorted(sign_count[cls].items()))
        print(f"  {cls:26s} {class_count[cls]:7d} share={class_count[cls]/total_pairs:8.5f} signs={signs}")
    print("  nonzero non-atom classes remain sign-balanced before packet labels break symmetry.")
    print()

    print("[3] Recent-work concept co-occurrence")
    print("  documents scanned:")
    for rel in DOCS:
        print(f"    {rel}")
    print("  concept hits:")
    for name, count in concept_hits.most_common():
        print(f"    {name:30s} {count}")
    print("  strongest co-occurrences:")
    for (a, b), count in concept_co.most_common(12):
        print(f"    {a:26s} <-> {b:26s} docs={count}")
    print()

    print("[4] Zipper theorem skeleton")
    teeth = [
        ("Z0 orthogonal cancellation", "disjoint Haar coordinate", "discard safely by discrepancy orthogonality"),
        ("Z1 boundary stop", "same-tile indicator", "AP/GW cocircuit or named boundary atom"),
        ("Z2 owner strip", "one coordinate fixed, one refined", "endpoint owner / Fejer center / Ramanujan label"),
        ("Z3 cross handoff", "opposite coordinate nesting", "zipper arrow between proof clocks"),
        ("Z4 nested descent", "same-direction nesting", "family compression or state-lift descent"),
        ("Z5 residual cocycle", "unpaired zeta survives all exits", "F7 / THM-572 state-lift obligation"),
    ]
    for name, algebra, route in teeth:
        print(f"  {name:24s} | {algebra:34s} | {route}")
    print()
    print("  Candidate zipper-cocycle theorem:")
    print("    On every labelled LRC14 packet fiber, each local mixed Haar cocycle")
    print("    is either paired by color-compatible discrepancy cancellation, stopped")
    print("    at a boundary cocircuit, handed to an owner/period/certificate clock,")
    print("    descended to a smaller packet family, or converted into a named")
    print("    state-lift residual.  If no tooth applies, the packet is the F7 sector.")
    print()

    print("[5] Tournament Analysis on proof carriers")
    print("  vertices: proof carriers, not runners")
    print("  pairwise observable: retained LRC predicate, fixed fiber, mixed sign, endpoint topology, certificate handoff, exposure test, state route, formal check")
    print("  switch/gauge: majority comparison of retention vector; ties use the declared Hamiltonian path")
    for carrier in CARRIERS:
        print(f"  {carrier.name:30s} score={len(adj[carrier.name])} vector={carrier.vector} role={carrier.role}")
    print(f"  score_hist={dict(sorted(Counter(len(v) for v in adj.values()).items()))}")
    print(f"  directed_3cycles={directed_3cycles(adj)}")
    print(f"  scc_sizes={scc_sizes(adj)}")
    print(f"  hamiltonian_path_count={len(paths)}")
    if paths:
        print("  canonical_path=" + " > ".join(paths[0]))
    print()

    print("[6] Assumption challenge")
    print("  alternate vertices considered: runners, dyadic rectangles, endpoint walls,")
    print("  row/column margin fibers, switch moves, color resonances, proof carriers,")
    print("  Fejer atoms, state-lift obligations, and theorem arrows.")
    print("  chosen vertices: proof carriers plus the local zipper cocycle.")
    print("  preserved: strict-witness predicate, fixed packet fiber, mixed Haar sign,")
    print("  endpoint/topology labels, and named certificate exits.")
    print("  destroyed: raw runner identity, continuous component count K, and scalar")
    print("  row/column shadows, unless zeta is retained or discharged.")


if __name__ == "__main__":
    main()
