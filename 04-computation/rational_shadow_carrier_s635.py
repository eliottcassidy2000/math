#!/usr/bin/env python3
"""S635: rationality as a two-shadow carrier obstruction.

The prompt's pi/e thread is corrected and then used as a model:

    If S=x+y and P=xy both descend to a tame field K, then x and y are
    roots of T^2 - S T + P over K.  Thus the unordered hidden pair has
    also descended to K up to a quadratic carrier.

For x=e and y=pi this proves that e+pi and e*pi cannot both be algebraic
(hence cannot both be rational).  It does not prove exactly one is rational;
it proves at least one is transcendental.

The finite-field section is a toy microscope for the same operation:
sum alone and product alone have large fibers, while the joint
(sum, product) carrier reconstructs the unordered pair.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations


def unordered_pairs_mod(p: int) -> list[tuple[int, int]]:
    return [(a, b) for a in range(p) for b in range(a, p)]


def finite_field_shadow_stats(p: int) -> dict[str, object]:
    pairs = unordered_pairs_mod(p)
    fibers: dict[str, dict[object, list[tuple[int, int]]]] = {
        "sum": defaultdict(list),
        "product": defaultdict(list),
        "joint": defaultdict(list),
    }

    for a, b in pairs:
        s = (a + b) % p
        prod = (a * b) % p
        fibers["sum"][s].append((a, b))
        fibers["product"][prod].append((a, b))
        fibers["joint"][(s, prod)].append((a, b))

    def summarize(name: str) -> dict[str, object]:
        sizes = [len(v) for v in fibers[name].values()]
        return {
            "distinct": len(sizes),
            "min_fiber": min(sizes),
            "max_fiber": max(sizes),
            "hist": dict(sorted(Counter(sizes).items())),
            "colliding_keys": sum(1 for size in sizes if size > 1),
        }

    return {
        "p": p,
        "pairs": len(pairs),
        "sum": summarize("sum"),
        "product": summarize("product"),
        "joint": summarize("joint"),
    }


@dataclass(frozen=True)
class Lens:
    name: str
    predicate: int
    reconstructs: int
    obstruction: int
    compression: int
    computable: int
    loss: int


LENSES = [
    Lens("joint_sum_product_carrier", 5, 5, 4, 4, 5, 1),
    Lens("field_descent_obstruction", 5, 4, 5, 3, 5, 1),
    Lens("lrc_relation_lattice", 5, 4, 5, 3, 4, 2),
    Lens("anti_coset_transporter", 4, 4, 4, 4, 4, 1),
    Lens("unit_spine_bulk_split", 4, 3, 4, 4, 4, 2),
    Lens("sequence_shadow_recursion", 3, 3, 4, 5, 5, 2),
    Lens("parity_redei_shadow", 3, 2, 5, 5, 5, 3),
    Lens("product_shadow_only", 2, 1, 3, 5, 5, 4),
    Lens("sum_shadow_only", 2, 1, 2, 5, 5, 5),
]

CRITERIA = [
    "predicate",
    "reconstructs",
    "obstruction",
    "compression",
    "computable",
]

TIE_PATH = [lens.name for lens in LENSES]


def build_lens_tournament() -> dict[str, object]:
    index = {lens.name: i for i, lens in enumerate(LENSES)}
    tie_rank = {name: i for i, name in enumerate(TIE_PATH)}
    edges: dict[tuple[int, int], int] = {}
    wins = Counter()

    for a, b in combinations(LENSES, 2):
        a_votes = 0
        b_votes = 0
        for criterion in CRITERIA:
            av = getattr(a, criterion)
            bv = getattr(b, criterion)
            if av > bv:
                a_votes += 1
            elif bv > av:
                b_votes += 1
        if a.loss < b.loss:
            a_votes += 1
        elif b.loss < a.loss:
            b_votes += 1

        if a_votes > b_votes:
            winner = a
        elif b_votes > a_votes:
            winner = b
        else:
            winner = a if tie_rank[a.name] < tie_rank[b.name] else b

        i, j = index[a.name], index[b.name]
        edges[(i, j)] = 1 if winner is a else 0
        wins[winner.name] += 1

    n = len(LENSES)
    directed_3 = 0
    for a, b, c in combinations(range(n), 3):
        out = Counter()
        for i, j in [(a, b), (a, c), (b, c)]:
            bit = edges[(min(i, j), max(i, j))]
            winner = min(i, j) if bit else max(i, j)
            out[winner] += 1
        if sorted(out.values()) == [1, 1, 1]:
            directed_3 += 1

    full_wins = {lens.name: wins[lens.name] for lens in LENSES}
    return {
        "edges": edges,
        "wins": full_wins,
        "score_hist": dict(sorted(Counter(full_wins.values()).items())),
        "directed_3cycles": directed_3,
        "hamiltonian_paths": hamiltonian_paths(edges, n),
        "sccs": strongly_connected_components(edges, n),
        "ranked": sorted((lens.name for lens in LENSES), key=lambda x: (-wins[x], tie_rank[x])),
    }


def beats(edges: dict[tuple[int, int], int], a: int, b: int) -> bool:
    if a < b:
        return edges[(a, b)] == 1
    return edges[(b, a)] == 0


def hamiltonian_paths(edges: dict[tuple[int, int], int], n: int) -> int:
    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if (prev_mask >> prev) & 1 and beats(edges, prev, last):
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def strongly_connected_components(edges: dict[tuple[int, int], int], n: int) -> list[list[int]]:
    adj = [[] for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if beats(edges, i, j):
            adj[i].append(j)
        else:
            adj[j].append(i)

    index = 0
    stack: list[int] = []
    on_stack = [False] * n
    indices = [-1] * n
    low = [0] * n
    comps: list[list[int]] = []

    def visit(v: int) -> None:
        nonlocal index
        indices[v] = low[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        for w in adj[v]:
            if indices[w] == -1:
                visit(w)
                low[v] = min(low[v], low[w])
            elif on_stack[w]:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            comps.append(sorted(comp))

    for v in range(n):
        if indices[v] == -1:
            visit(v)
    return comps


def fmt_hist(hist: dict[int, int]) -> str:
    return ", ".join(f"{k}:{v}" for k, v in sorted(hist.items()))


def main() -> str:
    lines: list[str] = []
    lines.append("S635 Rational Shadow Carrier")
    lines.append("============================")
    lines.append("")
    lines.append("Correction")
    lines.append("----------")
    lines.append(
        "The usable theorem is not 'exactly one of pi+e and pi*e is rational'."
    )
    lines.append(
        "It is stronger in one direction and weaker in another: pi+e and pi*e "
        "cannot both be algebraic; therefore at least one is transcendental "
        "and hence irrational.  Which one, or whether both are transcendental, "
        "is not known."
    )
    lines.append("")
    lines.append("Algebraic carrier lemma")
    lines.append("-----------------------")
    lines.append(
        "For any field K, if S=x+y and P=xy both lie in K, then x and y are "
        "roots of T^2 - S*T + P in K[T]."
    )
    lines.append(
        "So the pair (x,y) has descended to a quadratic carrier over K.  "
        "Contraposition: if a hidden coordinate is known transcendental over K, "
        "then the two shadows S and P cannot both be K-algebraic/tame."
    )
    lines.append(
        "With K=Qbar and (x,y)=(e,pi), both shadows algebraic would force e "
        "and pi algebraic over Qbar and hence algebraic over Q, contradicting "
        "transcendence."
    )
    lines.append("")
    lines.append("Finite-field microscope")
    lines.append("-----------------------")
    lines.append(
        "Unordered pairs in F_p show the same carrier geometry without analysis:"
    )
    lines.append(
        "sum alone and product alone have fibers, but the joint (sum, product) "
        "shadow is injective because it is exactly the monic quadratic."
    )
    lines.append("")
    lines.append(
        "p  pairs  sum(dist,max,hist)        product(dist,max,hist)      joint max"
    )
    for p in [5, 7, 11, 13, 17]:
        stats = finite_field_shadow_stats(p)
        s = stats["sum"]
        q = stats["product"]
        j = stats["joint"]
        lines.append(
            f"{p:2d} {stats['pairs']:6d} "
            f"({s['distinct']:2d},{s['max_fiber']:2d},{fmt_hist(s['hist']):<6}) "
            f"({q['distinct']:2d},{q['max_fiber']:2d},{fmt_hist(q['hist']):<16}) "
            f"{j['max_fiber']}"
        )
    lines.append("")
    lines.append("Reading")
    lines.append("-------")
    lines.append(
        "Rationality is a field-of-definition quotient, not merely a decimal "
        "flavor.  Odd/even is the same kind of move at the smallest possible "
        "field: a quotient to F_2.  The difference is size, not species."
    )
    lines.append(
        "Addition versus multiplication changes role here: x+y and xy are not "
        "rival operations but complementary elementary symmetric coordinates. "
        "Each one-shadow is lossy; together they reconstruct the unordered pair."
    )
    lines.append(
        "This is the same carrier lesson as S633/S634: a scalar shadow can look "
        "mysterious because the missing side channel is doing the work."
    )
    lines.append("")
    lines.append("Connections")
    lines.append("-----------")
    rows = [
        (
            "pi/e",
            "sum S and product P",
            "both algebraic would descend the hidden pair",
        ),
        (
            "LRC",
            "reset period plus relation lattice",
            "both too tame would miss the resonance carrier",
        ),
        (
            "unit distance",
            "spine edge count plus bulk shell",
            "57 is meaningful as 20+C_hex(3), not as a bare scalar",
        ),
        (
            "SC tournaments",
            "rooted perspective plus converse transporter",
            "fixed/merged/nonfixed layers replace raw complement collapse",
        ),
        (
            "sequence work",
            "value plus shadow family",
            "next-term opacity becomes recursion around the missing channel",
        ),
    ]
    for domain, shadows, obstruction in rows:
        lines.append(f"- {domain}: {shadows}; {obstruction}.")
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append(
        "Vertices are proof-carrier lenses, not constants or runners.  The "
        "pairwise observable is which lens better preserves the target predicate, "
        "reconstructs hidden state, exposes obstruction, compresses search, and "
        "stays computable.  The switch/gauge is majority vote over those criteria, "
        "with lower quotient-loss winning one extra vote; ties follow the listed "
        "Hamiltonian path."
    )
    tour = build_lens_tournament()
    lines.append("Score ranking:")
    for name in tour["ranked"]:
        lines.append(f"  {tour['wins'][name]:2d}  {name}")
    lines.append(f"score histogram: {tour['score_hist']}")
    lines.append(f"directed 3-cycles: {tour['directed_3cycles']}")
    lines.append(f"Hamiltonian paths: {tour['hamiltonian_paths']}")
    readable_sccs = [
        [LENSES[i].name for i in comp]
        for comp in tour["sccs"]
    ]
    lines.append(f"SCCs: {readable_sccs}")
    lines.append("")
    lines.append("Assumption challenge")
    lines.append("--------------------")
    lines.append(
        "For LRC/tournament transfer, do not assume vertices must be runners, "
        "constants, or tournament vertices.  In this lab the useful vertices are "
        "proof obligations and carrier choices.  Preserved predicate: whether a "
        "quotient can still certify the target obstruction.  Destroyed data: "
        "representatives, embeddings, exact times, and individual labels."
    )
    lines.append(
        "The challenged assumption is that rational/irrational is a terminal "
        "classification.  Here it is a warning that a field-of-definition shadow "
        "has lost a side channel."
    )
    lines.append("")
    lines.append("Next algorithmic tests")
    lines.append("----------------------")
    lines.append(
        "1. Add a generic two_shadow_obstruction helper: if two quotient shadows "
        "would reconstruct a forbidden/tame hidden object, at least one shadow "
        "must carry the missing complexity."
    )
    lines.append(
        "2. In LRC, classify speed sets by (reset period, relation lattice rank, "
        "short circuit support) before any dense time search."
    )
    lines.append(
        "3. In unit-distance beams, store edge count as (spine, bulk, symmetry "
        "transporter) rather than a single scalar."
    )
    lines.append(
        "4. In tournament counts, extend fixed/merged/nonfixed layers to H-spectra "
        "and strong-component norm spectra."
    )
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())
