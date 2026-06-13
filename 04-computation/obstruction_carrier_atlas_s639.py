#!/usr/bin/env python3
"""S639: shared obstruction-carrier atlas.

This packet ties together the current threads:

* H=21 impossibility for tournament Hamiltonian paths;
* LRC n=14 with C=2n-1=27;
* unit-distance n=21 and the 57=20+37 spine/bulk split;
* A000568 and self-converse/Burnside shadow counts;
* Schanuel/pi-e trace-norm-discriminant carriers;
* twin primes and Goldbach as local-prime side-channel problems.

The goal is not to solve the famous conjectures.  It is to make the common
proof shape explicit enough to generate new finite tests and transfer moves.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import gcd, prod


@dataclass(frozen=True)
class Carrier:
    name: str
    scalar: str
    hidden_witness: str
    side_channels: tuple[str, ...]
    tags: tuple[str, ...]
    finite_probe: str
    open_move: str


CARRIERS = [
    Carrier(
        name="H21_tournament_gap",
        scalar="H=21",
        hidden_witness="strong components plus OCF/odd-cycle packet",
        side_channels=("strongness", "c3", "SCC", "OCF", "Phi3(2)", "self-converse"),
        tags=("finite-window", "tournament", "forbidden-scalar", "root-of-unity", "side-condition", "parity"),
        finite_probe="c3<=10 population extincts at m=13; H=7,21 are permanent gaps",
        open_move="Use H-gap as a quotient-collapse test in non-tournament ledgers.",
    ),
    Carrier(
        name="LRC_n14",
        scalar="n=14, C=27=3^3",
        hidden_witness="lonely time and runner-owner certificate",
        side_channels=("unit/nonunit shells", "gcd strata", "owner", "carry", "lift", "self-converse round classes"),
        tags=("finite-window", "LRC", "local-prime", "root-of-unity", "side-condition", "additive-cover"),
        finite_probe="AP/V* floor rows live on mod-27 shell ledgers; worry set collapses to 64 SC round classes",
        open_move="Run side-channel jackknives over owner/carry/lift to price each residual.",
    ),
    Carrier(
        name="unit_distance_n21",
        scalar="21 vertices, 57 unit edges",
        hidden_witness="Moser slab with unit Hamiltonian spine",
        side_channels=("spine=20", "bulk=37", "direction support", "traceability", "embedding", "Eisenstein shell"),
        tags=("geometry", "root-of-unity", "side-condition", "unit-distance", "spine-bulk", "finite-window"),
        finite_probe="stored exact 57-edge n=21 cores are traceable; 57=20+C_hex(3)",
        open_move="Search scalar twins with same edge count but different spine/direction packets.",
    ),
    Carrier(
        name="A000568_shadow_counts",
        scalar="isomorphism count sequence",
        hidden_witness="Burnside fixed layers and transporter quotients",
        side_channels=("cycle type", "automorphism", "self-converse fixed layer", "q-deformation", "merged/nonfixed pairs"),
        tags=("tournament", "Burnside", "shadow-sequence", "side-condition", "symmetry", "transporters"),
        finite_probe="SC(1..10)=1,1,2,2,8,12,88,176,2752,8784; SC(2m)=A(m,4)",
        open_move="When a term is hard, compute fixed/merged/q/transporter companions first.",
    ),
    Carrier(
        name="Schanuel_pi_e",
        scalar="S=e+pi, P=e*pi",
        hidden_witness="algebraic independence of e and pi",
        side_channels=("trace", "norm", "discriminant", "log(pi)", "power sums", "transverse shadows"),
        tags=("transcendence", "trace-norm", "side-condition", "shadow-sequence", "famous-conjecture"),
        finite_probe="any two of S,P,D reconstruct; algebraic S/P would be lonely by power-sum fallout",
        open_move="Build generic two-root transverse-shadow helper and apply to other carriers.",
    ),
    Carrier(
        name="twin_primes",
        scalar="gap=2",
        hidden_witness="two primes sharing all local residue filters",
        side_channels=("local residues", "singular series", "parity barrier", "level of distribution", "almost-prime shadows"),
        tags=("local-prime", "sieve", "famous-conjecture", "additive-gap", "parity", "side-condition"),
        finite_probe="local survivor density is product over p of (p-2)/(p-1) after removing p=2",
        open_move="Use LRC-style channel jackknife: price each modulus/residue channel instead of one global density.",
    ),
    Carrier(
        name="Goldbach",
        scalar="N=p+q",
        hidden_witness="prime pair representation of an even integer",
        side_channels=("local residues", "major/minor arcs", "parity barrier", "singular series", "representation count"),
        tags=("local-prime", "sieve", "famous-conjecture", "additive-cover", "side-condition", "density"),
        finite_probe="local residues survive iff N mod p is not forcing both primes into zero classes",
        open_move="Try spine/bulk decomposition: forced local pair skeleton plus flexible representation bulk.",
    ),
]


TRANSFER_METHODS = [
    "side_channel_jackknife",
    "transverse_shadow_fallout",
    "finite_window_extinction",
    "spine_bulk_decomposition",
    "Burnside_transporters",
    "local_prime_ledger",
    "quotient_collapse_test",
    "Schanuel_style_completion",
]

METHOD_TAGS = {
    "side_channel_jackknife": {"side-condition", "finite-window", "local-prime"},
    "transverse_shadow_fallout": {"trace-norm", "shadow-sequence", "transcendence"},
    "finite_window_extinction": {"finite-window", "forbidden-scalar", "tournament"},
    "spine_bulk_decomposition": {"spine-bulk", "geometry", "additive-cover"},
    "Burnside_transporters": {"Burnside", "symmetry", "transporters", "tournament"},
    "local_prime_ledger": {"local-prime", "sieve", "parity"},
    "quotient_collapse_test": {"forbidden-scalar", "side-condition", "root-of-unity"},
    "Schanuel_style_completion": {"transcendence", "famous-conjecture", "trace-norm"},
}


def shell_strata(C: int) -> dict[int, list[int]]:
    strata: dict[int, list[int]] = defaultdict(list)
    for r in range(1, C):
        a = min(r, C - r)
        if r != a:
            continue
        strata[gcd(a, C)].append(a)
    return dict(sorted(strata.items()))


def twin_local_survivors(primes: list[int]) -> list[tuple[int, int, int, float]]:
    rows = []
    for k in range(1, len(primes) + 1):
        ps = primes[:k]
        M = prod(ps)
        ok = sum(1 for a in range(M) if gcd(a, M) == 1 and gcd(a + 2, M) == 1)
        rows.append((M, ok, M, ok / M))
    return rows


def goldbach_local_survivors(primes: list[int], N: int) -> list[tuple[int, int, int, float]]:
    rows = []
    for k in range(1, len(primes) + 1):
        ps = primes[:k]
        M = prod(ps)
        ok = sum(1 for a in range(M) if gcd(a, M) == 1 and gcd(N - a, M) == 1)
        rows.append((M, ok, M, ok / M))
    return rows


def carrier_similarity():
    rows = []
    for a, b in combinations(CARRIERS, 2):
        A, B = set(a.tags), set(b.tags)
        shared = sorted(A & B)
        union = A | B
        score = len(shared) / len(union)
        if shared:
            rows.append((score, a.name, b.name, shared))
    return sorted(rows, reverse=True)


def method_scores():
    rows = []
    for method, tags in METHOD_TAGS.items():
        hits = []
        for carrier in CARRIERS:
            overlap = sorted(tags & set(carrier.tags))
            if overlap:
                hits.append((carrier.name, overlap))
        rows.append((len(hits), method, hits))
    return sorted(rows, reverse=True)


@dataclass(frozen=True)
class Route:
    name: str
    finite: int
    transfer: int
    famous_relevance: int
    repo_fit: int
    risk: int


ROUTES = [
    Route("side_channel_jackknife", 5, 5, 4, 5, 1),
    Route("local_prime_ledger", 4, 5, 5, 4, 2),
    Route("transverse_shadow_fallout", 4, 4, 5, 4, 2),
    Route("spine_bulk_decomposition", 5, 4, 3, 5, 1),
    Route("finite_window_extinction", 5, 3, 3, 5, 1),
    Route("Burnside_transporters", 4, 4, 3, 5, 2),
    Route("quotient_collapse_test", 3, 5, 3, 4, 2),
    Route("Schanuel_style_completion", 1, 3, 5, 3, 5),
    Route("raw_scalar_numerology", 1, 1, 1, 2, 5),
]

CRITERIA = ["finite", "transfer", "famous_relevance", "repo_fit"]


def route_tournament():
    idx = {route.name: i for i, route in enumerate(ROUTES)}
    edges: dict[tuple[int, int], int] = {}
    wins = Counter({route.name: 0 for route in ROUTES})
    for a, b in combinations(ROUTES, 2):
        av = 0
        bv = 0
        for criterion in CRITERIA:
            if getattr(a, criterion) > getattr(b, criterion):
                av += 1
            elif getattr(b, criterion) > getattr(a, criterion):
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


def beats(edges: dict[tuple[int, int], int], a: int, b: int) -> bool:
    if a < b:
        return edges[(a, b)] == 1
    return edges[(b, a)] == 0


def hamiltonian_count(edges: dict[tuple[int, int], int], n: int) -> int:
    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if ((prev_mask >> prev) & 1) and beats(edges, prev, last):
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def directed_3cycles(edges: dict[tuple[int, int], int], n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        out = Counter()
        for i, j in [(a, b), (a, c), (b, c)]:
            winner = i if beats(edges, i, j) else j
            out[winner] += 1
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def main() -> str:
    lines: list[str] = []
    lines.append("S639 Obstruction Carrier Atlas")
    lines.append("==============================")
    lines.append("")
    lines.append("Carrier rows")
    lines.append("------------")
    for carrier in CARRIERS:
        lines.append(f"* {carrier.name}: scalar={carrier.scalar}")
        lines.append(f"  hidden={carrier.hidden_witness}")
        lines.append(f"  side_channels={', '.join(carrier.side_channels)}")
        lines.append(f"  finite_probe={carrier.finite_probe}")
        lines.append(f"  open_move={carrier.open_move}")
    lines.append("")
    lines.append("Concrete side-channel ledgers")
    lines.append("-----------------------------")
    strata = shell_strata(27)
    lines.append(f"LRC n=14 shell strata for C=27: {strata}")
    lines.append("Twin-prime local survivor densities over primorial filters:")
    for M, ok, total, density in twin_local_survivors([2, 3, 5, 7, 11]):
        lines.append(f"  mod {M:<5} survivors={ok:<5}/{total:<5} density={density:.6f}")
    lines.append("Goldbach local survivor densities for N=42 over primorial filters:")
    for M, ok, total, density in goldbach_local_survivors([2, 3, 5, 7, 11], 42):
        lines.append(f"  mod {M:<5} survivors={ok:<5}/{total:<5} density={density:.6f}")
    lines.append("")
    lines.append("Highest carrier similarities by shared tags")
    lines.append("-------------------------------------------")
    for score, a, b, shared in carrier_similarity()[:12]:
        lines.append(f"{score:.3f}  {a}  <->  {b}  shared={shared}")
    lines.append("")
    lines.append("Transfer-method coverage")
    lines.append("------------------------")
    for count, method, hits in method_scores():
        hit_summary = "; ".join(f"{name}:{'/'.join(tags)}" for name, tags in hits)
        lines.append(f"{method} hits {count}/{len(CARRIERS)} carriers: {hit_summary}")
    lines.append("")
    lines.append("New hypotheses / applications")
    lines.append("-----------------------------")
    lines.append("1. Side-channel jackknife should be run before density arguments in twin/Goldbach sieves, just as LRC n=14 prices unit/nonunit/carry/owner channels.")
    lines.append("2. Unit-distance scalar twins should be searched by fixing edge count and varying spine/bulk/direction support; the H=21 gap says raw scalar equality is not a proof object.")
    lines.append("3. A000568 hard terms should be extended by transporter and q-shadow companions before chasing the next main term; this is the finite analogue of Schanuel transverse shadows.")
    lines.append("4. H=21 impossibility can act as a quotient-collapse detector: if a model quotient demands H=21, it has forgotten a side condition.")
    lines.append("5. Goldbach representations may benefit from a spine/bulk ledger: local residue skeleton first, flexible representation bulk second, minor arcs as the damage channel.")
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append("Assumption challenge: I considered vertices as problems, primes, scalars, side channels, proof obligations, and transfer methods.  Transfer-method vertices preserve the actionable predicate: which move should be tried next across domains.  They destroy object-specific geometry, so the atlas also prints carrier rows above.")
    edges, wins = route_tournament()
    n = len(ROUTES)
    order = sorted(range(n), key=lambda i: wins[ROUTES[i].name], reverse=True)
    score_hist = Counter(wins.values())
    lines.append(f"vertices={n}")
    lines.append(f"score_hist={dict(sorted(score_hist.items()))}")
    lines.append(f"directed_3cycles={directed_3cycles(edges, n)}")
    lines.append(f"Hamiltonian_paths={hamiltonian_count(edges, n)}")
    lines.append("ranking:")
    for rank, idx in enumerate(order, start=1):
        route = ROUTES[idx]
        lines.append(f"  {rank}. {route.name} (score={wins[route.name]})")
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())
