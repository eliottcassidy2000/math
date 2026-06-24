#!/usr/bin/env python3
"""S165: 2D Haar product rule as tournament-tiling proof algebra.

This is a synthesis scout, not an LRC14 proof.

Recent sessions separated LRC14 rows into open Haar/Baire safe intervals,
boundary-only AP/GW atoms, tope/cocircuit wall packets, Fejer/K33/state-lift
handoffs, and labelled certificate routes.  The prompt asks how discrepancy
theory and the two-dimensional Haar product rule create the same structure as
the tournament tiling model.

Here the model is deliberately finite and exact:

* 1D unnormalised Haar functions satisfy a local product rule:
      h_I h_J is zero, an indicator, or a signed Haar child.
* 2D Haar rectangles multiply coordinatewise.
* Coordinatewise nesting/equality/disjointness induces the same finite set of
  interaction types as fixed-Hamiltonian-path tournament tiles:
      independent zero, same tile, one-coordinate owner strip, nested
      refinement, and cross-coordinate handoff.

The Tournament Analysis below uses product interaction types as vertices.
Pairwise observable: labelled information retained by a Haar-product quotient.
Switch/gauge: orient toward the carrier that remembers more of the LRC proof
packet (bulk locality, endpoint owner, boundary atom, zipper/crossing,
state-lift refinement, quotient guardrail).  Ties use the listed Hamiltonian
path order.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True, order=True)
class Interval:
    level: int
    index: int

    @property
    def left_num(self) -> int:
        return self.index

    @property
    def right_num(self) -> int:
        return self.index + 1

    @property
    def denom(self) -> int:
        return 1 << self.level

    def contains(self, other: "Interval") -> bool:
        return (
            self.left_num * other.denom <= other.left_num * self.denom
            and other.right_num * self.denom <= self.right_num * other.denom
        )

    def disjoint(self, other: "Interval") -> bool:
        return (
            self.right_num * other.denom <= other.left_num * self.denom
            or other.right_num * self.denom <= self.left_num * other.denom
        )

    def half_sign_on(self, inner: "Interval") -> int:
        """Return +1 if inner lies in this interval's left half, else -1."""

        # Compare inner midpoint to this midpoint, avoiding floats.
        # inner_mid = (2*idx+1)/(2*2^level), outer_mid=(2*idx+1)/(2*2^level)
        lhs = (2 * inner.index + 1) * self.denom
        rhs = (2 * self.index + 1) * inner.denom
        return 1 if lhs < rhs else -1

    def text(self) -> str:
        return f"[{self.index}/{self.denom},{self.index + 1}/{self.denom})"


@dataclass(frozen=True)
class Rect:
    x: Interval
    y: Interval

    def text(self) -> str:
        return f"{self.x.text()} x {self.y.text()}"


@dataclass(frozen=True)
class Product1D:
    kind: str
    relation: str
    sign: int
    carrier: Interval | None


def one_d_product(a: Interval, b: Interval) -> Product1D:
    if a == b:
        return Product1D("indicator", "equal", 1, a)
    if a.disjoint(b):
        return Product1D("zero", "disjoint", 0, None)
    if a.contains(b):
        return Product1D("haar", "A_contains_B", a.half_sign_on(b), b)
    if b.contains(a):
        return Product1D("haar", "B_contains_A", b.half_sign_on(a), a)
    raise AssertionError("dyadic intervals should be nested, equal, or disjoint")


@dataclass(frozen=True)
class Product2D:
    class_name: str
    sign: int
    x: Product1D
    y: Product1D


def product_class(px: Product1D, py: Product1D) -> str:
    if px.kind == "zero" or py.kind == "zero":
        return "orthogonal_zero"
    if px.kind == "indicator" and py.kind == "indicator":
        return "same_tile_indicator"
    if px.kind == "indicator" and py.kind == "haar":
        return "vertical_owner_strip"
    if px.kind == "haar" and py.kind == "indicator":
        return "horizontal_owner_strip"

    rels = {px.relation, py.relation}
    if rels == {"A_contains_B"} or rels == {"B_contains_A"}:
        return "nested_refinement"
    if rels == {"A_contains_B", "B_contains_A"}:
        return "cross_handoff"
    raise AssertionError((px, py))


def rect_product(a: Rect, b: Rect) -> Product2D:
    px = one_d_product(a.x, b.x)
    py = one_d_product(a.y, b.y)
    sign = px.sign * py.sign
    return Product2D(product_class(px, py), sign, px, py)


def intervals(max_depth: int) -> list[Interval]:
    out: list[Interval] = []
    for level in range(max_depth + 1):
        for index in range(1 << level):
            out.append(Interval(level, index))
    return out


def rects(max_depth: int) -> list[Rect]:
    ints = intervals(max_depth)
    return [Rect(x, y) for x in ints for y in ints]


def ordered_product_census(max_depth: int) -> tuple[Counter[str], dict[str, Counter[int]]]:
    count: Counter[str] = Counter()
    signs: dict[str, Counter[int]] = defaultdict(Counter)
    rs = rects(max_depth)
    for a in rs:
        for b in rs:
            prod = rect_product(a, b)
            count[prod.class_name] += 1
            signs[prod.class_name][prod.sign] += 1
    return count, signs


CLASS_ORDER = [
    "same_tile_indicator",
    "vertical_owner_strip",
    "horizontal_owner_strip",
    "cross_handoff",
    "nested_refinement",
    "orthogonal_zero",
]


CLASS_READOUT = {
    "same_tile_indicator": (
        "boundary atom",
        "same tile in both coordinates; endpoint/cocircuit atom survives",
        (3, 3, 3, 0, 0, 3),
    ),
    "vertical_owner_strip": (
        "one-owner strip",
        "same x-owner, refined y-scale; one endpoint label survives",
        (3, 2, 2, 0, 1, 3),
    ),
    "horizontal_owner_strip": (
        "one-owner strip",
        "same y-owner, refined x-scale; one endpoint label survives",
        (3, 2, 2, 0, 1, 3),
    ),
    "cross_handoff": (
        "zipper handoff",
        "opposite coordinate nesting; handoff between two wall clocks",
        (3, 1, 3, 3, 2, 2),
    ),
    "nested_refinement": (
        "state lift",
        "same-direction nesting; descend to smaller packet/fiber",
        (3, 1, 1, 1, 3, 2),
    ),
    "orthogonal_zero": (
        "independent cancellation",
        "disjoint in at least one coordinate; discrepancy orthogonality",
        (1, 0, 0, 0, 0, 1),
    ),
}


def dot_score(name: str) -> int:
    # The weights emphasize LRC packet safety over raw scalar mass.
    weights = (4, 5, 5, 4, 4, 6)
    vec = CLASS_READOUT[name][2]
    return sum(w * v for w, v in zip(weights, vec))


def orient(a: str, b: str) -> tuple[str, str]:
    sa = dot_score(a)
    sb = dot_score(b)
    if sa > sb:
        return a, b
    if sb > sa:
        return b, a
    # Tie Hamiltonian path: earlier class in CLASS_ORDER wins.
    return (a, b) if CLASS_ORDER.index(a) < CLASS_ORDER.index(b) else (b, a)


def tournament_edges(vertices: list[str]) -> dict[tuple[str, str], str]:
    edges: dict[tuple[str, str], str] = {}
    for a, b in combinations(vertices, 2):
        winner, loser = orient(a, b)
        edges[(winner, loser)] = "retains-more-labels"
    return edges


def score_histogram(vertices: list[str], edges: dict[tuple[str, str], str]) -> Counter[int]:
    scores: Counter[str] = Counter()
    for a, _ in edges:
        scores[a] += 1
    return Counter(scores[v] for v in vertices)


def has_edge(edges: dict[tuple[str, str], str], a: str, b: str) -> bool:
    return (a, b) in edges


def directed_3cycles(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    total = 0
    for a, b, c in combinations(vertices, 3):
        if (
            has_edge(edges, a, b)
            and has_edge(edges, b, c)
            and has_edge(edges, c, a)
        ) or (
            has_edge(edges, a, c)
            and has_edge(edges, c, b)
            and has_edge(edges, b, a)
        ):
            total += 1
    return total


def hamiltonian_path_count(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    # Small n; dynamic programming over subsets.
    index = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp: dict[tuple[int, int], int] = {}
    for v in vertices:
        dp[(1 << index[v], index[v])] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if not ways:
                continue
            a = vertices[last]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                b = vertices[nxt]
                if has_edge(edges, a, b):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + ways
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def strongly_connected_components(vertices: list[str], edges: dict[tuple[str, str], str]) -> list[list[str]]:
    graph = {v: [] for v in vertices}
    rev = {v: [] for v in vertices}
    for a, b in edges:
        graph[a].append(b)
        rev[b].append(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(comp)
    return comps


def print_census(max_depth: int = 3) -> None:
    counts, signs = ordered_product_census(max_depth)
    total = sum(counts.values())
    print("[1] 2D Haar rectangle product census")
    print(f"  dyadic depth <= {max_depth}; rectangles={(2 ** (max_depth + 1) - 1) ** 2}; ordered pairs={total}")
    print("  product rule: h_R*h_S = (h_I*h_I')*(h_J*h_J')")
    print("  1D factors are zero, indicator, or signed nested Haar.")
    print()
    print(f"  {'class':24s} {'count':>9s} {'share':>9s} {'signs':>22s} readout")
    for name in CLASS_ORDER:
        c = counts[name]
        sign_text = ", ".join(f"{k}:{v}" for k, v in sorted(signs[name].items()))
        print(
            f"  {name:24s} {c:9d} {c / total:9.5f} {sign_text:>22s} "
            f"{CLASS_READOUT[name][1]}"
        )
    print()


def print_dictionary() -> None:
    print("[2] Tournament-tiling dictionary")
    rows = [
        (
            "Haar rectangle R=I x J",
            "fixed-path tile (r,c)",
            "two coordinates carry row/column endpoint clocks",
        ),
        (
            "disjoint coordinate",
            "noninteracting tiles",
            "safe to separate by discrepancy orthogonality",
        ),
        (
            "same coordinate",
            "shared endpoint/owner strip",
            "must retain labelled owner, not just scalar mass",
        ),
        (
            "same tile",
            "boundary cocircuit atom",
            "AP/GW-style endpoint atom; Haar mass can be zero",
        ),
        (
            "same-direction nesting",
            "recursive tile refinement",
            "state-lift/fiber descent, not a global quotient",
        ),
        (
            "opposite-coordinate nesting",
            "zipper/cross handoff",
            "the packet gluing case: one clock refines while another coarsens",
        ),
    ]
    for left, mid, right in rows:
        print(f"  {left:31s} | {mid:29s} | {right}")
    print()


def print_tournament() -> None:
    vertices = CLASS_ORDER[:]
    edges = tournament_edges(vertices)
    print("[3] Proof-carrier Tournament Analysis")
    print("  vertices: Haar-product interaction classes")
    print("  pairwise observable: labelled information retained by quotient")
    print("  switch/gauge: orient toward larger weighted retention vector")
    print("  tie Hamiltonian path:", " -> ".join(vertices))
    print()
    print(f"  {'vertex':24s} {'score':>5s} {'dot':>4s} {'role':20s} vector")
    out_scores = Counter()
    for a, _ in edges:
        out_scores[a] += 1
    for v in vertices:
        role, _, vec = CLASS_READOUT[v]
        print(f"  {v:24s} {out_scores[v]:5d} {dot_score(v):4d} {role:20s} {vec}")
    print()
    print("  score_histogram=", dict(sorted(score_histogram(vertices, edges).items())))
    print("  directed_3cycles=", directed_3cycles(vertices, edges))
    print("  scc_sizes=", [len(c) for c in strongly_connected_components(vertices, edges)])
    print("  hamiltonian_path_count=", hamiltonian_path_count(vertices, edges))
    print()


def print_lrc_target() -> None:
    print("[4] LRC14 proof target suggested by the product algebra")
    print(
        "  Haar-tile discrepancy lemma (target): on each labelled packet fiber, "
        "a non-AP/GW zero-open residual must expose a nonzero 2D Haar tile "
        "coefficient in one of the owner-strip, cross-handoff, or nested-refinement "
        "classes."
    )
    print(
        "  Interpretation: if all such coefficients vanish, the packet has no "
        "local discrepancy address left except same-tile boundary atoms; that "
        "should force the AP/GW boundary skeleton or emit a new THM-572/F7 "
        "state-lift atom."
    )
    print(
        "  Guardrail: the quotient may forget scalar noise, but it may not forget "
        "endpoint owners, C27 transfers, exact-period/Ramanujan labels, or K33 "
        "state-lift debt.  This matches the HYP-2986 handoff atlas."
    )
    print()


def main() -> None:
    print_census(max_depth=3)
    print_dictionary()
    print_tournament()
    print_lrc_target()


if __name__ == "__main__":
    main()
