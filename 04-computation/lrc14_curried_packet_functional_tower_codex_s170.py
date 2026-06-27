#!/usr/bin/env python3
"""Curried packet-functional bridge for LRC14.

This is a synthesis scout, not a proof.  It records the function signatures
that appear when the current LRC14 packet stack, the additive-basis/Fibonacci
work, and the Farey operator lanes are treated as curried maps instead of
scalar invariants.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Iterable


@dataclass(frozen=True)
class Carrier:
    name: str
    signature: str
    partial_evaluation: str
    lrc_use: str
    forgets: str
    retention: tuple[int, ...]


FEATURES = (
    "predicate",
    "fiber_labels",
    "argument_order",
    "dual_certificate",
    "residual_name",
    "farey_scale",
    "additive_carry",
    "anti_scalar",
)


CARRIERS = [
    Carrier(
        "total_curried_packet_evaluator",
        "S -> P(S) -> root -> lane -> fiber -> certificate -> verdict",
        "E(S)(packet)(p/q,e)(farey_lane)(fiber)(cert)",
        "Keeps the full LRC proof object live until the last certificate.",
        "Forgets nothing until a named discharge is returned.",
        (3, 3, 3, 3, 3, 3, 3, 3),
    ),
    Carrier(
        "quotient_cocycle_derivative",
        "Q -> y -> (x,x' in fiber_y) -> lost_coordinate_delta",
        "d_Q(coord)(y)(x,x')",
        "Turns every quotient into an explicit lost-coordinate function.",
        "Forgets the coordinate only after zero/coboundary/dual/residual tests.",
        (3, 3, 3, 2, 3, 2, 2, 3),
    ),
    Carrier(
        "residual_section_router",
        "packet -> section -> grid_exit -> next_certificate_family",
        "R(packet)(section)(owner_strip/cross/nested)",
        "Routes HYP-2963 packets through HYP-2996 residual sections.",
        "Forgets phase location after owner/cross/nested exit is labelled.",
        (3, 3, 2, 2, 3, 2, 1, 3),
    ),
    Carrier(
        "fejer_toeplitz_dual_functional",
        "S -> k -> divisor_fiber -> coefficient -> quadratic_form",
        "Phi(S)(k)(v|k)(atom_bank)",
        "Makes the PSD-cover implication testable on labelled packet fibers.",
        "Forgets endpoint owners unless the packet/family key is retained.",
        (3, 2, 2, 3, 2, 1, 1, 3),
    ),
    Carrier(
        "farey_lane_scheduler",
        "p -> q -> e=14p-q -> lane -> payload",
        "F(p)(q)(e)(root/sum/product/power)",
        "Separates affine-safe, incidence, and stress-test Farey payloads.",
        "Forgets route labels if the payload is used without exact root data.",
        (2, 2, 3, 1, 2, 3, 1, 3),
    ),
    Carrier(
        "pascal_path_rank_section",
        "n -> k -> C(n-k-1,k), then sum_k or choose carry support",
        "A(n)(k)",
        "Keeps the Fibonacci row fiber before Zeckendorf scalarization.",
        "Forgets route/certificate data; only carries additive row support.",
        (2, 2, 2, 1, 1, 1, 3, 3),
    ),
    Carrier(
        "additive_basis_currency",
        "atom_system -> N -> budget -> representation_fiber",
        "B(atoms)(N)(budget)",
        "Classifies a packet as smoothing, bounded arity, or normal form.",
        "Forgets LRC scale unless tethered to the packet root.",
        (2, 2, 2, 1, 1, 1, 2, 2),
    ),
    Carrier(
        "runner_movie_tournament_shadow",
        "S -> t -> pair_observable -> tournament",
        "T(S)(t)(gauge)",
        "Useful exploratory shadow for pairwise structure.",
        "Forgets endpoint debt, exact M, and often the quotient fiber.",
        (1, 1, 1, 0, 0, 0, 0, 2),
    ),
    Carrier(
        "raw_scalar_evaluation",
        "S -> number",
        "m(S), H(T), count(S), product(p,q)",
        "Negative control for quotient guardrails.",
        "Forgets nearly every proof-bearing argument.",
        (1, 0, 0, 0, 0, 0, 0, 0),
    ),
]


def path_rank_row(n: int) -> list[int]:
    """Return [C(n-k-1,k)]_k for the user's Fibonacci indexing."""
    if n <= 0:
        return []
    row = []
    k = 0
    while n - k - 1 >= k:
        row.append(comb(n - k - 1, k))
        k += 1
    return row or [1]


def farey_payload(p: int, q: int, lane: str) -> int | tuple[int, int] | tuple[int, int, int]:
    if lane == "root":
        return (p, q, 14 * p - q)
    if lane == "sum":
        return p + q
    if lane == "product":
        return p * q
    if lane == "power":
        return (q**p, p**q)
    raise ValueError(lane)


def orient(a: Carrier, b: Carrier, tie_order: list[str]) -> tuple[str, str]:
    wins_a = sum(x > y for x, y in zip(a.retention, b.retention))
    wins_b = sum(y > x for x, y in zip(a.retention, b.retention))
    if wins_a > wins_b:
        return a.name, b.name
    if wins_b > wins_a:
        return b.name, a.name
    return (
        (a.name, b.name)
        if tie_order.index(a.name) < tie_order.index(b.name)
        else (b.name, a.name)
    )


def tournament_edges(carriers: list[Carrier]) -> set[tuple[str, str]]:
    order = [c.name for c in carriers]
    edges: set[tuple[str, str]] = set()
    for a, b in combinations(carriers, 2):
        edges.add(orient(a, b, order))
    return edges


def score_hist(carriers: list[Carrier], edges: set[tuple[str, str]]) -> dict[int, int]:
    scores = {c.name: 0 for c in carriers}
    for u, _ in edges:
        scores[u] += 1
    hist: dict[int, int] = {}
    for score in scores.values():
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def directed_3cycles(carriers: list[Carrier], edges: set[tuple[str, str]]) -> int:
    names = [c.name for c in carriers]
    count = 0
    for a, b, c in combinations(names, 3):
        ab = (a, b) in edges
        bc = (b, c) in edges
        ca = (c, a) in edges
        ba = (b, a) in edges
        cb = (c, b) in edges
        ac = (a, c) in edges
        if (ab and bc and ca) or (ba and cb and ac):
            count += 1
    return count


def scc_sizes(carriers: list[Carrier], edges: set[tuple[str, str]]) -> list[int]:
    names = [c.name for c in carriers]
    adj = {n: [] for n in names}
    radj = {n: [] for n in names}
    for u, v in edges:
        adj[u].append(v)
        radj[v].append(u)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes = []

    def rdfs(v: str) -> int:
        seen.add(v)
        total = 1
        for w in radj[v]:
            if w not in seen:
                total += rdfs(w)
        return total

    for name in reversed(order):
        if name not in seen:
            sizes.append(rdfs(name))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(carriers: list[Carrier], edges: set[tuple[str, str]]) -> int:
    names = [c.name for c in carriers]
    idx = {name: i for i, name in enumerate(names)}
    n = len(names)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if (names[last], names[nxt]) in edges:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    return sum(dp[(1 << n) - 1])


def render_row(row: Iterable[int]) -> str:
    return "+".join(str(x) for x in row)


def main() -> None:
    print("LRC14 CURRIED PACKET FUNCTIONAL TOWER")
    print("status: synthesis scout / proof-interface guardrail, not a proof")
    print()

    print("Core total evaluator:")
    print("  E : S -> P(S) -> root -> lane -> fiber -> certificate -> verdict")
    print("  A quotient is a partial evaluation of E plus a declared lost-coordinate function.")
    print()

    print("Quotient audit law:")
    print("  lost_Q(coord)(y)(x,x') = coord(x) - coord(x') on the fiber Q^-1(y)")
    print("  allowed iff lost_Q is zero, reconstructed, coboundary/exact,")
    print("  dual-annihilated, descended to a family, AP/GW boundary equality,")
    print("  or routed to F7/THM-572 residual debt.")
    print()

    print("Fibonacci path-rank partial evaluations:")
    for n in range(2, 10):
        row = path_rank_row(n)
        print(f"  A({n})(k) row={render_row(row):<12} sum={sum(row)}")
    print()

    print("LRC14 unit-excess Farey partial evaluations p -> q -> lane:")
    for p in range(1, 7):
        q = 14 * p - 1
        root = farey_payload(p, q, "root")
        sum_lane = farey_payload(p, q, "sum")
        product_lane = farey_payload(p, q, "product")
        pow_digits = tuple(len(str(x)) for x in farey_payload(p, q, "power"))
        print(
            f"  F({p})({q}): root={root}, sum={sum_lane}, "
            f"product={product_lane}, power_digit_lengths={pow_digits}"
        )
    print()

    print("Named-row curried readout templates:")
    templates = [
        (
            "AP",
            "E(AP)(boundary_packet)(e=0)(root/sum)(same_tile)(no-open)",
            "AP/GW boundary equality atom",
        ),
        (
            "GW",
            "E(GW)(C27_boundary_packet)(e=0)(root/sum)(same_tile)(no-open)",
            "AP/GW boundary equality atom with C27 transfer label",
        ),
        (
            "K33 12->36",
            "E(row)(K33_packet)(p=3)(product)(cross_handoff)(Fejer d159)",
            "strict-safe certificate or named K33/state-lift debt",
        ),
        (
            "P10+GW",
            "E(row)(unit_petal_packet)(root)(sum/product)(owner_strip)(Fejer d280)",
            "strict-safe certificate after petal/GW labels are retained",
        ),
        (
            "covering 12->84",
            "E(row)(covering_packet)(root)(sum)(nested_refinement)(boundary_moment)",
            "positive covering/boundary-moment packet, not scalar all-covered",
        ),
    ]
    for name, call, verdict in templates:
        print(f"  {name}:")
        print(f"    call:    {call}")
        print(f"    verdict: {verdict}")
    print()

    print("Function-carrier retention vectors")
    print("features:", ", ".join(FEATURES))
    for carrier in CARRIERS:
        print(f"  {carrier.name}: {carrier.retention}")
        print(f"    type:    {carrier.signature}")
        print(f"    partial: {carrier.partial_evaluation}")
        print(f"    LRC:     {carrier.lrc_use}")
        print(f"    warning: {carrier.forgets}")
    print()

    edges = tournament_edges(CARRIERS)
    print("Tournament Analysis")
    print("  vertices: curried function carriers / proof obligations, not runners")
    print("  pairwise observable: coordinatewise retention vector over features")
    print("  switch/gauge: majority of strictly larger retention coordinates")
    print("  tie Hamiltonian path:", " > ".join(c.name for c in CARRIERS))
    print(f"  score_hist={score_hist(CARRIERS, edges)}")
    print(f"  directed_3cycles={directed_3cycles(CARRIERS, edges)}")
    print(f"  scc_sizes={scc_sizes(CARRIERS, edges)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(CARRIERS, edges)}")
    print("  top path:")
    for i, carrier in enumerate(CARRIERS, 1):
        print(f"    {i}. {carrier.name}")
    print()

    print("Proof-program takeaway:")
    print("  Treat every LRC map as a curried function.  The moment an argument is")
    print("  precomposed, evaluated, summed, or quotient-collapsed, record the lost")
    print("  coordinate as a function on the remaining fiber.  The proof target is")
    print("  not a new scalar; it is totality of this evaluator: every primitive")
    print("  LRC14 packet returns strict-safe, AP/GW boundary equality, or a named")
    print("  F7/THM-572 residual before the final scalar shadow is reached.")


if __name__ == "__main__":
    main()
