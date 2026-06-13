#!/usr/bin/env python3
"""S640: resonance product-ledger atlas across H=21, LRC, UD, primes, A000568.

The aim is not to prove Goldbach or twin primes.  It is to turn the repo's
recent pattern into a portable proof-design object:

    scalar quotient + retained side channel + residual obstruction.

This script records the shared shape, factors the recurring small integers,
and runs a method-level Tournament Analysis over possible next applications.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations


def factor(n: int) -> str:
    if n == 0:
        return "0"
    if n == 1:
        return "1"
    x = abs(n)
    parts = []
    p = 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            parts.append((p, e))
        p += 1 if p == 2 else 2
    if x > 1:
        parts.append((x, 1))
    return " * ".join(str(p) if e == 1 else f"{p}^{e}" for p, e in parts)


@dataclass(frozen=True)
class ObjectRow:
    name: str
    scalar: str
    integers: tuple[int, ...]
    side_channel: str
    status: str
    proof_instrument: str
    barrier: str


OBJECTS = [
    ObjectRow(
        "H=21 impossibility",
        "Hamiltonian-path count 21",
        (7, 21),
        "strong-component norm, c3/odd-cycle budget, OCF packet",
        "closed: 7 and 21 are permanent H-gaps in repo ledger",
        "finite strong-window closure plus Busch/Moon lower bounds",
        "raw H forgets the strong/OCF side conditions",
    ),
    ObjectRow(
        "LRC n=14",
        "C=2n-1=27 shell clock",
        (14, 27, 3, 9),
        "gcd strata 1/3/9, unit/nonunit shells, carry/owner data",
        "open in full; Res_27 least-positive quotient classified",
        "THM-407 shell fold, Res_27 pinch certificate, carry conservativity",
        "least-positive quotient forgets lift/CRT branch ownership",
    ),
    ObjectRow(
        "unit distance n=21",
        "21 vertices, 57 unit edges",
        (21, 57, 20, 37),
        "20-edge unit spine plus 37=C_hex(3) bulk, traceability, ears",
        "exact stored 57-edge cores traceable; n=22 frontier open",
        "Hamiltonian-spine DP and Moser carrier beam",
        "raw edge count forgets embedding, direction support, spine/bulk split",
    ),
    ObjectRow(
        "A000568 / self-converse",
        "unrooted tournament class count",
        (12, 56, 456, 88),
        "Burnside partitions, fixed/merged/nonfixed, q-deformation, transporters",
        "active method: extend hard sequences by companion shadows",
        "Burnside q-shadow and anti-coset transporter quotient",
        "raw next term forgets fixed-point and orbit side channels",
    ),
    ObjectRow(
        "Schanuel / pi-e",
        "algebraic/rational shadow",
        (2,),
        "algebraic independence, trace/norm/discriminant, transverse shadows",
        "conditional completion; S636 proves at least two of S/P/D transcend",
        "Vieta carrier, Lindemann-Weierstrass fallout, Schanuel conditional",
        "individual transcendence does not control algebraic curves",
    ),
    ObjectRow(
        "twin primes",
        "prime gap 2",
        (2, 246),
        "admissible tuple, singular series, sieve weights, parity barrier",
        "open; bounded gaps known, H_1 <= 246 in Polymath8b ledger",
        "Maynard-Tao/Polymath sieve weights",
        "parity problem: sieve sees almost-primes better than primes",
    ),
    ObjectRow(
        "Goldbach",
        "even N as p+q; odd N as p+p+p",
        (2, 3),
        "major/minor arcs, local congruences, singular series, density",
        "weak/ternary proved by Helfgott; strong/binary open",
        "circle method plus explicit minor-arc bounds",
        "minor arcs and local-to-global transfer are the hard residue side",
    ),
]


@dataclass(frozen=True)
class Route:
    name: str
    preserves_predicate: int
    cross_problem: int
    computable_now: int
    proof_strength: int
    novelty: int
    risk: int


ROUTES = [
    Route("local_obstruction_product_ledger", 5, 5, 4, 4, 4, 1),
    Route("scalar_twins_side_channel_atlas", 5, 5, 5, 3, 4, 1),
    Route("transverse_shadow_fallout_helper", 4, 4, 5, 4, 4, 1),
    Route("discriminant_branch_labels_for_LRC", 5, 4, 4, 4, 5, 2),
    Route("A000568_prime_sieve_shadow", 4, 5, 3, 3, 5, 2),
    Route("unit_distance_ear_singular_series", 4, 4, 3, 3, 5, 3),
    Route("H_spectrum_same_scalar_fibers", 4, 4, 4, 3, 3, 2),
    Route("Schanuel_conditional_transfer", 3, 4, 2, 5, 4, 4),
    Route("raw_scalar_numerology", 1, 2, 5, 1, 1, 5),
]

CRITERIA = [
    "preserves_predicate",
    "cross_problem",
    "computable_now",
    "proof_strength",
    "novelty",
]


def route_tournament():
    idx = {r.name: i for i, r in enumerate(ROUTES)}
    edges = {}
    wins = {r.name: 0 for r in ROUTES}
    for a, b in combinations(ROUTES, 2):
        av = 0
        bv = 0
        for c in CRITERIA:
            if getattr(a, c) > getattr(b, c):
                av += 1
            elif getattr(b, c) > getattr(a, c):
                bv += 1
        if a.risk < b.risk:
            av += 1
        elif b.risk < a.risk:
            bv += 1
        winner = a if av >= bv else b
        i, j = idx[a.name], idx[b.name]
        edges[(i, j)] = 1 if winner is a else 0
        wins[winner.name] += 1
    return edges, wins


def beats(edges, a, b):
    if a < b:
        return edges[(a, b)] == 1
    return edges[(b, a)] == 0


def hamiltonian_count(edges, n):
    @lru_cache(maxsize=None)
    def dp(mask, last):
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


def directed_3cycles(edges, n):
    total = 0
    for a, b, c in combinations(range(n), 3):
        out = Counter()
        for i, j in [(a, b), (a, c), (b, c)]:
            winner = i if beats(edges, i, j) else j
            out[winner] += 1
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def factor_table() -> list[str]:
    all_ints = sorted({x for row in OBJECTS for x in row.integers if x > 1})
    return [f"{x:>4} = {factor(x)}" for x in all_ints]


def main() -> str:
    lines = []
    lines.append("S640 Resonance Product-Ledger Atlas")
    lines.append("====================================")
    lines.append("")
    lines.append("Thesis")
    lines.append("------")
    lines.append(
        "The shared object is not a single theorem.  It is a proof-design "
        "schema: scalar quotient + retained side channel + residual obstruction."
    )
    lines.append(
        "The failures all happen when the scalar is treated as the object.  "
        "Progress happens when the side channel is made first-class."
    )
    lines.append("")
    lines.append("Atlas")
    lines.append("-----")
    for row in OBJECTS:
        lines.append(f"- {row.name}")
        lines.append(f"  scalar: {row.scalar}")
        lines.append(f"  side channel: {row.side_channel}")
        lines.append(f"  status: {row.status}")
        lines.append(f"  instrument: {row.proof_instrument}")
        lines.append(f"  barrier: {row.barrier}")
    lines.append("")
    lines.append("Recurring small integers")
    lines.append("------------------------")
    lines.extend(factor_table())
    lines.append("")
    lines.append("Prime-reading")
    lines.append("-------------")
    lines.append(
        "Prime 2 is the parity/traceability seam: Redei oddness, twin gap 2, "
        "Goldbach parity, LRC even-n first-doubling, and dyadic carry data."
    )
    lines.append(
        "Prime 3 is the first fragmentation prime: LRC n=14 has C=27=3^3, "
        "Goldbach ternary is the solved additive prime theorem, and the "
        "Phi_3/root-of-unity carrier keeps recurring around 7 and 21."
    )
    lines.append(
        "Prime 7 is the forbidden-scalar warning: H=7 and H=21 are impossible "
        "as tournament path counts, while unit-distance n=21 is legal only "
        "because it retains spine/bulk side channels."
    )
    lines.append("")
    lines.append("New applications")
    lines.append("----------------")
    applications = [
        (
            "Singular-series style LRC ledger",
            "For each prime divisor of C=2n-1, record local shell survival "
            "just as twin-prime/Goldbach singular series record local residue "
            "survival.  The n=14 row is the prime-3 test case.",
        ),
        (
            "Same-scalar twin atlas",
            "Systematically search for pairs with the same scalar but different "
            "side channels: same H but different OCF packet, same unit-edge "
            "count but different spines, same C shell but different carry owner.",
        ),
        (
            "Transverse-shadow helper",
            "Generalize S638: once one scalar descends, prove that transverse "
            "shadows become impossible/lonely unless the hidden source descends.",
        ),
        (
            "A000568 sieve shadow",
            "Treat Burnside cycle types like local prime obstructions: fixed, "
            "merged, nonfixed, q-deformed, and transporter shadows are the "
            "tournament analogue of admissible residue data.",
        ),
        (
            "Unit-distance ear singular series",
            "For candidate ears from n=21 cores, record local direction/support "
            "admissibility before embedding search, mirroring prime tuple filters.",
        ),
    ]
    for name, text in applications:
        lines.append(f"- {name}: {text}")
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append(
        "Vertices are next-proof routes, not problem domains.  Pairwise observable: "
        "which route preserves the target predicate, transfers across problems, "
        "is computable now, has proof strength, has novelty, and has lower risk. "
        "Tie Hamiltonian path is the listed route order."
    )
    edges, wins = route_tournament()
    ranked = sorted((r.name for r in ROUTES), key=lambda n: (-wins[n], n))
    for name in ranked:
        lines.append(f"  {wins[name]:2d}  {name}")
    lines.append(f"score histogram: {dict(sorted(Counter(wins.values()).items()))}")
    lines.append(f"directed 3-cycles: {directed_3cycles(edges, len(ROUTES))}")
    lines.append(f"Hamiltonian paths: {hamiltonian_count(edges, len(ROUTES))}")
    lines.append("")
    lines.append("Assumption challenge")
    lines.append("--------------------")
    lines.append(
        "Alternate vertices considered: problem domains, integers, primes, "
        "tournaments, unit-distance points, LRC shells, Burnside partitions, "
        "prime tuples, circle-method arcs, algebraic shadows, and proof routes. "
        "This lab chooses proof routes."
    )
    lines.append(
        "Preserved predicate: whether a proposed quotient still certifies the "
        "domain's hard yes/no property.  Destroyed data: representatives, "
        "embeddings, lift choices, branch labels, local residue weights, and "
        "which transverse shadow carries the obstruction."
    )
    lines.append("")
    lines.append("Concrete next build")
    lines.append("-------------------")
    lines.append(
        "Build a local_obstruction_product_ledger helper with a common schema: "
        "domain, modulus/local prime, forbidden residues, surviving residues, "
        "weight, and lost side channel.  First instances: LRC C=27 shells, "
        "Goldbach mod-p local obstructions, twin-prime admissible tuples, and "
        "unit-distance direction masks."
    )
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())
