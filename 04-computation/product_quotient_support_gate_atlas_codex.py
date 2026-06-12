#!/usr/bin/env python3
"""Product-quotient obstruction atlas for arXiv:2508.14876.

codex-2026-06-12

This is not a proof of a new theorem.  It is a structured merge note between
Benjamin Church's product-quotient obstruction paper and the repo's recent
LRC14 / [72,36,16] support-gate investigations.

The goal is to keep three carriers separate:

1. scalar conditions such as Shioda supersingularity or Gleason positivity;
2. retained support channels such as diagonal symmetric forms under all
   partial Frobenius twists, Q27 resource ledgers, or minimum-word designs;
3. exact numerology that may suggest search seeds but is not itself evidence.

Tournament Analysis is included over proof-route vertices, not over curves,
runners, or codewords.  This quotient is deliberately lossy: it records which
ideas currently carry proof leverage across domains.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from math import comb, gcd


ARXIV_ID = "2508.14876"
GENUS = 14
HURWITZ_CONSTANT = 84
PSL_Q = 13
LRC14_SINGLE_STRANGER_CUTOFF = 13 * 84
Q27_RESCUE_WITNESS = 91
TYPE_II_72_LAMBDA5 = 78
SUPERSINGULAR_PRIMES = [1091, 2339, 6551, 7643]
SUBGROUPS = {
    "D6": 12,
    "D7": 14,
    "A4": 12,
}


def factor(n: int) -> str:
    parts: list[str] = []
    p = 2
    x = n
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            parts.append(f"{p}^{e}" if e > 1 else str(p))
        p += 1 if p == 2 else 2
    if x > 1:
        parts.append(str(x))
    return " * ".join(parts)


def psl2_order(q: int) -> int:
    divisor = gcd(2, q - 1)
    return q * (q * q - 1) // divisor


def antipodal_class(r: int, m: int) -> int | None:
    r %= m
    if r == 0 or gcd(r, m) != 1:
        return None
    return min(r, (-r) % m)


def residues(n: int, moduli: list[int]) -> dict[int, str]:
    out: dict[int, str] = {}
    for m in moduli:
        r = n % m
        label = str(r)
        if r == m - 1:
            label += " (-1)"
        cls = antipodal_class(r, m)
        if cls is not None and m in {27, 91}:
            label += f", +/-{cls}"
        out[m] = label
    return out


@dataclass(frozen=True)
class Route:
    name: str
    domain: str
    scalar: str
    retained_channel: str
    status: str
    channel: int
    theorem: int
    repo_bridge: int
    computable: int
    risk: int

    def score_tuple(self) -> tuple[int, int, int, int, int]:
        return (self.channel, self.theorem, self.repo_bridge, self.computable, -self.risk)


ROUTES = [
    Route(
        "diagonal_forms_all_partial_frobenii",
        "product-quotient surfaces",
        "supersingular reductions are scalar/cohomological",
        "diagonal symmetric forms survive every asymmetric Frobenius twist",
        "paper theorem",
        5,
        5,
        5,
        3,
        1,
    ),
    Route(
        "lrc14_Q27_resource_ledger",
        "LRC14",
        "q is blocked / not blocked",
        "shell-27 classes, 13-clock, divisor fibers, Bprime target",
        "exact for one-stranger; open resource theorem",
        5,
        4,
        5,
        5,
        2,
    ),
    Route(
        "order5_marked_fixed_projection",
        "[72,36,16] code",
        "fixed projected Type II [16,8] scalar type",
        "marked tetrad support and heptad/Fano boundary",
        "exact conditional branch reduction",
        5,
        5,
        4,
        5,
        1,
    ),
    Route(
        "code72_minimum_design_support",
        "[72,36,16] code",
        "formal Gleason enumerator is nonnegative",
        "minimum words must realize a 5-(72,16,78) design",
        "open existence / exact scalar and design parameters",
        5,
        3,
        5,
        4,
        3,
    ),
    Route(
        "hurwitz_PSL2_13_indices",
        "numerology bridge",
        "|PSL2(F13)| = 1092",
        "subgroup indices 91 and 78, genus-14 Hurwitz lattice",
        "exact arithmetic; speculative proof leverage",
        3,
        2,
        5,
        5,
        4,
    ),
    Route(
        "random_pentagonal_support_law",
        "Euler partitions",
        "finite-window reciprocal Lyapunov slope",
        "sign support / root placement of sparse pentagonal denominator",
        "computational; analytic theorem open",
        3,
        2,
        4,
        4,
        3,
    ),
    Route(
        "typeII_gleason_scalar_gate",
        "[72,36,16] code",
        "Gleason polynomial positivity and total mass",
        "none by itself; support must be reattached",
        "exact scalar calculation",
        1,
        5,
        4,
        5,
        4,
    ),
    Route(
        "shioda_supersingular_scalar",
        "surfaces over finite fields",
        "Shioda supersingular cohomology",
        "none by itself; Church counterexample shows scalar is insufficient",
        "paper theorem/counterexample context",
        1,
        5,
        4,
        2,
        4,
    ),
]


def route_tournament(routes: list[Route]) -> dict[str, object]:
    n = len(routes)
    wins: dict[int, set[int]] = {i: set() for i in range(n)}
    edge_flips_vs_channel_only = 0
    for i, j in combinations(range(n), 2):
        a, b = routes[i], routes[j]
        if a.score_tuple() > b.score_tuple():
            winner, loser = i, j
        elif b.score_tuple() > a.score_tuple():
            winner, loser = j, i
        else:
            winner, loser = i, j
        wins[winner].add(loser)

        channel_pref = i if (a.channel, -a.risk) >= (b.channel, -b.risk) else j
        if channel_pref != winner:
            edge_flips_vs_channel_only += 1

    scores = {routes[i].name: len(wins[i]) for i in range(n)}
    hist = Counter(scores.values())
    cycles: list[tuple[str, str, str]] = []
    for i, j, k in combinations(range(n), 3):
        if j in wins[i] and k in wins[j] and i in wins[k]:
            cycles.append((routes[i].name, routes[j].name, routes[k].name))
        if k in wins[i] and j in wins[k] and i in wins[j]:
            cycles.append((routes[i].name, routes[k].name, routes[j].name))

    # Kosaraju SCCs.
    graph = {i: wins[i] for i in range(n)}
    rgraph: dict[int, set[int]] = defaultdict(set)
    for i, outs in graph.items():
        for j in outs:
            rgraph[j].add(i)
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
        for w in rgraph[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[int] = []
            rdfs(v, comp)
            sccs.append([routes[i].name for i in comp])

    hp_count = 0
    hp_first: tuple[str, ...] | None = None
    for perm in permutations(range(n)):
        if all(perm[i + 1] in wins[perm[i]] for i in range(n - 1)):
            hp_count += 1
            if hp_first is None:
                hp_first = tuple(routes[i].name for i in perm)

    return {
        "scores": scores,
        "hist": dict(sorted(hist.items())),
        "directed_3cycles": len(cycles),
        "cycle_examples": cycles[:5],
        "sccs": sccs,
        "hamiltonian_paths": hp_count,
        "first_hamiltonian_path": hp_first,
        "edge_flips_vs_channel_only": edge_flips_vs_channel_only,
    }


def print_paper_bridge() -> None:
    group_order = psl2_order(PSL_Q)
    hurwitz_order = HURWITZ_CONSTANT * (GENUS - 1)
    print("[1] Paper facts pulled into repo coordinates")
    print(f"  arXiv:{ARXIV_ID}: product-quotient surfaces over finite fields")
    print(f"  genus g={GENUS}; Hurwitz bound 84*(g-1)={hurwitz_order}")
    print(f"  |PSL2(F_{PSL_Q})|={group_order} = {factor(group_order)}")
    print(f"  LRC14 single-stranger cutoff 13*84={LRC14_SINGLE_STRANGER_CUTOFF}")
    print(f"  same integer? {group_order == hurwitz_order == LRC14_SINGLE_STRANGER_CUTOFF}")
    print("  subgroups used in the examples:")
    for name, order in SUBGROUPS.items():
        index = group_order // order
        extras: list[str] = []
        if index == Q27_RESCUE_WITNESS:
            extras.append("Q27 q=91 rescue witness")
        if index == TYPE_II_72_LAMBDA5:
            extras.append("Type II [72] lambda_5=78")
        if index == comb(14, 2):
            extras.append("C(14,2)")
        if index == comb(13, 2):
            extras.append("C(13,2)")
        print(
            f"    {name}: |H|={order:2d}, index={index:3d}"
            + (f" ({'; '.join(extras)})" if extras else "")
        )
    print(
        f"  further identities: 1092=12*C(14,2)={12*comb(14,2)}, "
        f"14*C(13,2)={14*comb(13,2)}, 3*C(14,3)={3*comb(14,3)}"
    )
    print()


def print_supersingular_residues() -> None:
    print("[2] Explicit supersingular prime residue signatures")
    moduli = [7, 13, 27, 91]
    print("  p      " + "  ".join(f"mod {m:<2d}" for m in moduli))
    for p in SUPERSINGULAR_PRIMES:
        sig = residues(p, moduli)
        print("  " + f"{p:<6d} " + "  ".join(f"{sig[m]:<13s}" for m in moduli))
    all_minus_one_13 = all(p % 13 == 12 for p in SUPERSINGULAR_PRIMES)
    missing_class_hits = [
        p for p in SUPERSINGULAR_PRIMES if antipodal_class(p, 27) == 10
    ]
    print(f"  all listed primes are -1 mod 13: {all_minus_one_13}")
    print(
        "  listed primes in the LRC14 missing +/-10 class mod 27: "
        f"{missing_class_hits}"
    )
    print("  caution: residue echoes are search beacons, not obstruction proofs.")
    print()


def print_gate_rows() -> None:
    print("[3] Scalar-versus-support gate rows")
    print(
        "  route                              channel theorem bridge compute risk  status"
    )
    for r in ROUTES:
        print(
            f"  {r.name:<35s} {r.channel:^7d} {r.theorem:^7d} "
            f"{r.repo_bridge:^6d} {r.computable:^7d} {r.risk:^4d}  {r.status}"
        )
        print(f"    scalar:  {r.scalar}")
        print(f"    support: {r.retained_channel}")
    print()


def print_tournament() -> None:
    fp = route_tournament(ROUTES)
    print("[4] Tournament Analysis over proof-route vertices")
    print(
        "  observable = (retained_channel, theorem_status, repo_bridge, "
        "computability, -risk)"
    )
    print("  tie Hamiltonian path = declaration order in ROUTES")
    print(f"  score histogram: {fp['hist']}")
    print(f"  directed 3-cycles: {fp['directed_3cycles']}")
    print(f"  SCCs: {fp['sccs']}")
    print(f"  Hamiltonian path count: {fp['hamiltonian_paths']}")
    print(f"  first Hamiltonian path: {fp['first_hamiltonian_path']}")
    print(
        "  edge flips versus channel-only ranking: "
        f"{fp['edge_flips_vs_channel_only']}"
    )
    print("  scores:")
    for name, score in sorted(fp["scores"].items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"    {score}: {name}")
    print()


def print_assumption_challenge() -> None:
    print("[5] Assumption challenge and next proof moves")
    print("  alternate vertex sets considered:")
    print(
        "    curves, quotient surfaces, diagonal forms, Frobenius twists, "
        "subgroup indices, supersingular prime residues, LRC runners, gaps, "
        "fixed circle sections, section boundaries, wall-crossing events, "
        "Pisano classes, code supports, matroid circuits, and proof obligations"
    )
    print("  chosen quotient:")
    print(
        "    proof-route vertices.  It preserves whether a scalar obstruction "
        "has a retained side channel and whether the channel is computable."
    )
    print("  information destroyed:")
    print(
        "    actual surface geometry, explicit Magma group data, code support "
        "realizability, and multi-stranger LRC interactions."
    )
    print("  challenged assumption:")
    print(
        "    do not treat supersingularity, Gleason positivity, or q-blocking as "
        "the obstruction.  In all three domains the proof-bearing object is the "
        "extra support channel retained after the scalar quotient."
    )
    print("  concrete next moves:")
    print(
        "    (a) use the D7 index 78 as a design-ledger search seed for the "
        "5-(72,16,78) minimum layer;"
    )
    print(
        "    (b) use the D6/A4 index 91 as the geometric analogue of the LRC14 "
        "q=91 fibered rescue;"
    )
    print(
        "    (c) model LRC14 Q27 descent after Church's partial-Frobenius descent: "
        "each blocked shell should either lower a degree/resource or land in a "
        "finite exceptional type."
    )
    print()


def main() -> None:
    print("Product-quotient support-gate atlas")
    print("arXiv:2508.14876 x LRC14 x Type II [72,36,16]")
    print()
    print_paper_bridge()
    print_supersingular_residues()
    print_gate_rows()
    print_tournament()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
