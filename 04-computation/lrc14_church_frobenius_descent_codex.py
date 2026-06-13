#!/usr/bin/env python3
"""Church partial-Frobenius descent as an LRC14 proof skeleton.

codex-2026-06-13

This script is a synthesis atlas, not a new exhaustive search.  It imports the
proof mechanism of Church's arXiv:2508.14876 product-quotient obstruction into
the now-sharper LRC14 Q27 program:

  Church:
    scalar Shioda supersingularity is too weak;
    diagonal forms on every asymmetric partial Frobenius twist are retained;
    non-exceptional rational/elliptic curves descend by a projection-degree drop.

  LRC14:
    scalar "plain q<=27 blocked" is too weak;
    Q27 obligation/support ledgers are retained;
    a non-exceptional blocker should either enter a certified finite atlas or
    descend by losing a resource coordinate.

The output records the exact finite facts already proved in the repo
(HYP-2443/2444/2463/2464/2465) and turns them into a decision-tree target for
the remaining LRC14 proof.

Tournament Analysis uses proof routes as vertices.  Candidate vertices
considered but rejected for this synthesis: runners, individual denominators,
unit twists, deleted core speeds, added speeds, subgroup indices, Frobenius
twists, curve components, and wall events.  Proof-route vertices preserve the
descent/exception structure and deliberately discard raw geometry.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import comb, gcd


N = 14
CORE = tuple(7 * k for k in range(1, 13))
MAX_R = 13 * 84
DIVISORS = (1, 2, 7, 14)
Q27 = tuple(sorted({d * m for d in DIVISORS for m in range(1, 28) if d * m > 1}))
HARD = (260, 351, 442, 611, 702, 793, 962, 1053)
SUPERSINGULAR_PRIME_NORMS = (1091, 2339, 6551, 7643)
SUBGROUP_ORDERS = {"D6": 12, "D7": 14, "A4": 12}


def factor(n: int) -> str:
    out: list[str] = []
    p = 2
    x = n
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            out.append(f"{p}^{e}" if e > 1 else str(p))
        p += 1 if p == 2 else 2
    if x > 1:
        out.append(str(x))
    return " * ".join(out)


def psl2_order(q: int) -> int:
    return q * (q * q - 1) // gcd(2, q - 1)


def antipodal_class(v: int, q: int) -> int | None:
    r = v % q
    if r == 0 or gcd(r, q) != 1:
        return None
    return min(r, (-r) % q)


def q_addresses(q: int) -> tuple[tuple[int, int], ...]:
    return tuple((d, q // d) for d in DIVISORS if q % d == 0 and q // d <= 27)


@dataclass(frozen=True)
class CertifiedBlock:
    name: str
    rows_or_cases: int
    scalar_shadow: str
    retained_channel: str
    certified_exit: str
    source: str


CERTIFIED_BLOCKS = (
    CertifiedBlock(
        "one_stranger_family",
        936,
        "plain shell q<=27 has 8 residual rows",
        "Q27 plus Bprime(any) on S(r)",
        "936/936 Q27 covered; 936/936 Bprime(any) covered",
        "HYP-2444/HYP-2443",
    ),
    CertifiedBlock(
        "hard_packet_replacement_hull",
        sum(comb(8, k) * comb(12, k - 1) for k in range(1, 9)),
        "eight shell-27/13-clock hard residues look stackable",
        "typed Q27 resource ledger",
        "77520/77520 have Q27 witnesses; only 10 miss plain q<=27",
        "HYP-2463",
    ),
    CertifiedBlock(
        "two_stranger_plain_residuals",
        877,
        "one deletion plus two additions can still block every plain shell",
        "13-clock debt, deleted-core address, shell-27 pair class, divisor fiber",
        "877/877 have Q27 witnesses; every residual has a 13-clock speed",
        "HYP-2464",
    ),
    CertifiedBlock(
        "near_core_setcover_cases",
        sum(comb(12, e) for e in range(0, 4)),
        "arbitrary additions might cover all Q27 obligations",
        "primitive Q27 set-cover over safe twists for CORE\\D",
        "all delete counts e<=3 infeasible with budget e+1",
        "HYP-2465",
    ),
)


@dataclass(frozen=True)
class Route:
    name: str
    scalar_shadow: str
    retained_channel: str
    descent_or_exception: str
    exactness: int
    lrc_leverage: int
    descent_drop: int
    computability: int
    risk: int

    def score_tuple(self) -> tuple[int, int, int, int, int]:
        return (
            self.retained_score,
            self.exactness,
            self.lrc_leverage + self.descent_drop,
            self.computability,
            -self.risk,
        )

    @property
    def retained_score(self) -> int:
        if "diagonal forms" in self.retained_channel or "Q27 obligation" in self.retained_channel:
            return 5
        if "resource" in self.retained_channel or "owner" in self.retained_channel:
            return 4
        if "indices" in self.retained_channel:
            return 2
        return 1


ROUTES = (
    Route(
        "church_diagonal_forms_all_twists",
        "Shioda supersingularity",
        "diagonal forms on every asymmetric partial Frobenius twist",
        "finite exceptional curve types or projection-degree drop",
        5,
        4,
        5,
        3,
        1,
    ),
    Route(
        "lrc14_Q27_obligation_setcover",
        "plain blocked/not-blocked denominator shadow",
        "Q27 obligation hypergraph and primitive set-cover number",
        "near-core rows retaining >=9 core speeds are impossible Q27 blockers",
        5,
        5,
        4,
        5,
        1,
    ),
    Route(
        "lrc14_resource_coordinate_compression",
        "hard residues as names",
        "resource vector: 13-clock, deleted-core address, shell-27 class, divisor fiber",
        "compress residuals or open low-clock/fiber witness",
        4,
        5,
        4,
        4,
        2,
    ),
    Route(
        "lrc14_owner_Bprime_exception",
        "row has no small scalar witness",
        "owner-private deletion and Bprime(any runner) support",
        "named exceptional opening rather than generic blocker",
        4,
        5,
        3,
        3,
        2,
    ),
    Route(
        "lrc14_below_nine_core_descent_gap",
        "not near-core after four or more core deletions",
        "support-load ledger over lower retained-core rank",
        "must descend to low clocks, AP/Vstar/2AP, owner, or contradiction",
        1,
        5,
        5,
        2,
        4,
    ),
    Route(
        "hurwitz_1092_91_78_beacons",
        "genus-14 PSL2(F13) arithmetic",
        "indices 91 and 78",
        "search beacons only; no proof without support channel",
        2,
        3,
        1,
        5,
        4,
    ),
    Route(
        "raw_plain_shell_blocking",
        "all q<=27 plain shells blocked",
        "none unless Q27/fiber/owner labels are attached",
        "no descent by itself; many noisy residuals",
        3,
        2,
        1,
        5,
        5,
    ),
    Route(
        "shioda_supersingular_scalar",
        "supersingular cohomological shadow",
        "none by itself",
        "Church counterexample shows scalar gate is insufficient",
        5,
        1,
        0,
        2,
        4,
    ),
)


def route_tournament(routes: tuple[Route, ...]) -> dict[str, object]:
    n = len(routes)
    wins: dict[int, set[int]] = {i: set() for i in range(n)}
    edge_flips_vs_scalar = 0
    for i, j in combinations(range(n), 2):
        a, b = routes[i], routes[j]
        if a.score_tuple() >= b.score_tuple():
            winner, loser = i, j
        else:
            winner, loser = j, i
        wins[winner].add(loser)
        scalar_winner = i if (a.exactness, a.computability) >= (b.exactness, b.computability) else j
        if scalar_winner != winner:
            edge_flips_vs_scalar += 1

    scores = {routes[i].name: len(wins[i]) for i in range(n)}
    cycles: list[tuple[str, str, str]] = []
    for i, j, k in combinations(range(n), 3):
        if j in wins[i] and k in wins[j] and i in wins[k]:
            cycles.append((routes[i].name, routes[j].name, routes[k].name))
        if k in wins[i] and j in wins[k] and i in wins[j]:
            cycles.append((routes[i].name, routes[k].name, routes[j].name))

    graph = {i: wins[i] for i in range(n)}
    reverse: dict[int, set[int]] = defaultdict(set)
    for i, outs in graph.items():
        for j in outs:
            reverse[j].add(i)
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)
    seen.clear()
    sccs: list[list[str]] = []

    def rdfs(v: int, comp: list[int]) -> None:
        seen.add(v)
        comp.append(v)
        for w in reverse[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[int] = []
            rdfs(v, comp)
            sccs.append([routes[i].name for i in sorted(comp)])

    return {
        "scores": scores,
        "score_hist": Counter(scores.values()),
        "cycles": cycles,
        "sccs": sccs,
        "edge_flips_vs_scalar": edge_flips_vs_scalar,
    }


def print_paper_bridge() -> None:
    order = psl2_order(13)
    print("[1] Church paper constants in LRC14 coordinates")
    print(f"  |PSL2(F_13)| = {order} = {factor(order)}")
    print(f"  Hurwitz genus-14 bound 84*(14-1) = {84 * (14 - 1)}")
    print(f"  LRC14 carry window 13*84 = {MAX_R}; same integer = {order == MAX_R}")
    for name, size in SUBGROUP_ORDERS.items():
        idx = order // size
        tags = []
        if idx == comb(14, 2):
            tags.append("C(14,2), LRC q=91 fiber")
        if idx == comb(13, 2):
            tags.append("C(13,2), code72 lambda_5")
        print(f"  subgroup {name:2s}: |H|={size:2d}, index={idx:3d} ({'; '.join(tags)})")
    print("  first supersingular prime norms and residue beacons:")
    for p in SUPERSINGULAR_PRIME_NORMS:
        cls27 = antipodal_class(p, 27)
        print(
            f"    {p:4d}: mod13={p % 13:2d}, mod27={p % 27:2d}, "
            f"pm27={cls27}, mod91={p % 91:2d}"
        )
    print("  caution: the constants are beacons; the proof object is the retained support channel.")
    print()


def print_certified_lrc_blocks() -> None:
    print("[2] Certified LRC14 blocks after the paper import")
    print(f"  Q27 size = {len(Q27)} denominators; Q27={Q27[:12]} ... {Q27[-6:]}")
    for block in CERTIFIED_BLOCKS:
        print(f"  {block.name}:")
        print(f"    rows/cases:       {block.rows_or_cases}")
        print(f"    scalar shadow:    {block.scalar_shadow}")
        print(f"    retained channel: {block.retained_channel}")
        print(f"    certified exit:   {block.certified_exit}")
        print(f"    source:           {block.source}")
    print()


def print_descent_target() -> None:
    print("[3] Candidate Church-style LRC14 descent theorem")
    print("  Suppose S is a primitive 13-speed row with no Q27 witness.")
    print("  The desired proof should force one of the following exits:")
    exits = [
        ("carry_window_near_core", "if S normalizes to 1..1092 and retains >=9 core speeds, HYP-2465 contradicts no-Q27"),
        ("hard_packet_stack", "if S compresses to shell-27/13-clock hard packets, HYP-2463 contradicts no-Q27"),
        ("two_stranger_resource", "if S is a one-delete/two-add plain residual, HYP-2464 gives Q27 plus coordinates"),
        ("owner_exception", "if support concentrates on a runner, Bprime/owner-private deletion should open"),
        ("floor_atom_exception", "AP, Vstar, and nonprimitive 2AP are named wall atoms, not new blockers"),
        ("below_nine_core", "the live gap: rows deleting >=4 core speeds need a support-load descent measure"),
        ("outside_window", "large speeds must trigger divisor/carry/Bprime normalization before the finite atlas applies"),
    ]
    for name, meaning in exits:
        print(f"    {name:23s} -> {meaning}")
    print()
    print("  Proposed descent measure for the live gap:")
    print("    R(S) = (outside_window?, deleted_core_count, Q27_setcover_defect,")
    print("            max_support_load, owner_private_width, low_clock_depth)")
    print("  Generic steps should lower R lexicographically; non-lowering steps must land")
    print("  in one of the named finite exceptions above.")
    print()


def print_route_tournament() -> None:
    result = route_tournament(ROUTES)
    print("[4] Tournament Analysis over proof routes")
    print("  observable=(retained channel, exactness, LRC leverage+descent drop, computability, -risk)")
    print("  tie Hamiltonian path=declaration order")
    print(f"  score histogram={dict(sorted(result['score_hist'].items()))}")
    print(f"  directed 3-cycles={len(result['cycles'])}")
    print(f"  SCCs={result['sccs']}")
    print(f"  edge flips versus scalar-only ranking={result['edge_flips_vs_scalar']}/{comb(len(ROUTES), 2)}")
    print("  route order:")
    for score, name in sorted(((s, n) for n, s in result["scores"].items()), reverse=True):
        print(f"    {score}: {name}")
    print()


def print_next_computation() -> None:
    print("[5] High-leverage next computations")
    moves = [
        (
            "below-nine-core MILP",
            "repeat the HYP-2465 primitive set-cover lower bound for delete_count=4, "
            "but add structured budgets by 13-clock/divisor/owner class instead of raw e+1",
        ),
        (
            "outside-window Bprime normalizer",
            "prove or test that speeds >1092 either dominate an existing core runner, "
            "open Bprime(any), or reduce modulo a carry/fiber without losing blockedness",
        ),
        (
            "support-transport curvature",
            "measure how Q27 obligation covers change under q->2q, q->7q, and 27->9->3; "
            "nonzero defect should be a witness and zero defect should descend",
        ),
        (
            "exception catalogue",
            "make AP/Vstar/2AP, q=91, q=161, owner-private, and low-clock exits a finite typed list, "
            "mirroring Church's finite curve types before descent",
        ),
    ]
    for name, text in moves:
        print(f"  {name}: {text}")


def main() -> None:
    print("Church partial-Frobenius descent x LRC14 Q27")
    print("=" * 72)
    print_paper_bridge()
    print_certified_lrc_blocks()
    print_descent_target()
    print_route_tournament()
    print_next_computation()


if __name__ == "__main__":
    main()
