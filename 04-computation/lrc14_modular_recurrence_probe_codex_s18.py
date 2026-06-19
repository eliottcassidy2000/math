#!/usr/bin/env python3
"""
LRC(14) modular recurrence probe.

Question: mod 30 and mod 6 are useful; what other modular/integer-sequence
addresses are we missing?

This script keeps the object deliberately small and exact.  It treats rational
centers a/b as modular carriers for the small-speed part P subset {1,...,13}.
For a fixed center a/b, a speed v is safe at the center iff

    ||v*a/b|| > 1/14.

For b=2,3 this recovers HYP-2598's universal-center survivor sequence.  Adding
b=5 and b=7 shows how mod 30 and the seven-sector mod-7 layer sit inside a
single squarefree divisor-profile recurrence.  The b>=4 centers are only
address carriers, not universal cluster witnesses: their center gaps are too
small to solve the large-spread cluster by themselves.

Tournament Analysis declaration.
  Vertex set: modular proof obligations, not runners:
    mod6_universal, mod30_address, mod210_address, all_small_denominators,
    support6_coimage_tail, primorial_height_escape, raw_runner_residues.
  Pairwise observable: lexicographic tuple
    (preserves LRC predicate, proof legitimacy, residual compression,
     recurrence/address value).
  Tie path: the listed order is used for exact ties.

Assumption challenge.
  I considered runners, residues, rational centers, prime masks, denominator
  sets, support-six coimage classes, and proof obligations.  Raw runner
  residues preserve too much irrelevant geometry and miss the recurrence.  The
  prime-mask quotient preserves the denominator-survival predicate and exposes
  the hidden inclusion-exclusion sequences, but it destroys cluster-width data;
  this is why b=5/7 are proof addresses rather than standalone witnesses.
"""
from __future__ import annotations

import itertools
import math
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

UNIVERSE = tuple(range(1, 14))
PRIMES = (2, 3, 5, 7)


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def choose(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def dist_num_mod(r: int, b: int) -> int:
    r %= b
    return min(r, b - r)


def center_safe_for_speed(a: int, b: int, v: int) -> bool:
    return 14 * dist_num_mod(a * v, b) > b


def has_safe_center(P: tuple[int, ...], b: int) -> bool:
    """Whether some reduced center a/b is strict-safe for every v in P."""
    for a in range(1, b):
        if math.gcd(a, b) != 1:
            continue
        if all(center_safe_for_speed(a, b, v) for v in P):
            return True
    return False


def survivor_counts_for_denoms(denoms: tuple[int, ...]) -> list[int]:
    counts = [0] * (len(UNIVERSE) + 1)
    for s in range(len(UNIVERSE) + 1):
        for P in itertools.combinations(UNIVERSE, s):
            if any(has_safe_center(P, b) for b in denoms):
                counts[s] += 1
    return counts


def prime_mask(v: int, primes: tuple[int, ...] = PRIMES) -> int:
    mask = 0
    for i, p in enumerate(primes):
        if v % p == 0:
            mask |= 1 << i
    return mask


def residual_counts_hitting_all(primes: tuple[int, ...]) -> list[int]:
    full = (1 << len(primes)) - 1
    counts = [0] * (len(UNIVERSE) + 1)
    for s in range(len(UNIVERSE) + 1):
        for P in itertools.combinations(UNIVERSE, s):
            mask = 0
            for v in P:
                mask |= prime_mask(v, primes)
            if mask == full:
                counts[s] += 1
    return counts


def inclusion_exclusion_survivors(primes: tuple[int, ...]) -> list[int]:
    """Count P that avoid multiples of at least one prime in primes."""
    out: list[int] = []
    idxs = range(len(primes))
    for s in range(len(UNIVERSE) + 1):
        total = 0
        for r in range(1, len(primes) + 1):
            sign = 1 if r % 2 == 1 else -1
            for J in itertools.combinations(idxs, r):
                available = sum(
                    1 for v in UNIVERSE if all(v % primes[j] != 0 for j in J)
                )
                total += sign * choose(available, s)
        out.append(total)
    return out


def polynomial_terms(primes: tuple[int, ...]) -> list[tuple[tuple[int, ...], int]]:
    idxs = range(len(primes))
    terms: list[tuple[tuple[int, ...], int]] = []
    for r in range(1, len(primes) + 1):
        sign = 1 if r % 2 == 1 else -1
        for J in itertools.combinations(idxs, r):
            available = sum(
                1 for v in UNIVERSE if all(v % primes[j] != 0 for j in J)
            )
            terms.append((tuple(primes[j] for j in J), sign * available))
    return terms


def fmt_seq(seq: list[int]) -> str:
    return ", ".join(str(x) for x in seq)


def print_universal_center_recovery() -> None:
    section("RECOVER HYP-2598 AND EXTEND THE PRIME-MASK SIEVE")
    for primes in ((2, 3), (2, 3, 5), (2, 3, 5, 7)):
        survivors = inclusion_exclusion_survivors(primes)
        residual = residual_counts_hitting_all(primes)
        print(f"prime centers {primes}")
        print(f"  survivor sequence: {fmt_seq(survivors)}")
        print(f"  residual hits-all sequence: {fmt_seq(residual)}")
        print("  inclusion-exclusion available-set terms:")
        terms = polynomial_terms(primes)
        pretty = []
        for J, signed_available in terms:
            sign = "+" if signed_available >= 0 else "-"
            pretty.append(f"{sign} C({abs(signed_available)},s) for avoid {J}")
        print("    " + "; ".join(pretty))
    print(
        "\nReadout: HYP-2598 is the (2,3) row.  The next exact sieve rows are "
        "(2,3,5) -> mod 30 and (2,3,5,7) -> mod 210.  The latter is not a "
        "standalone proof, but it is the natural address lattice for combining "
        "parity/triadic centers with the seven-sector support-six tail."
    )


def print_center_denominator_table() -> None:
    section("RATIONAL CENTER DENOMINATORS b <= 30")
    rows = []
    for b in range(2, 31):
        counts = survivor_counts_for_denoms((b,))
        # choose a middle size that appears in the LRC14 residual tables.
        mid5 = counts[5]
        mid8 = counts[8]
        strict_gap = Fraction(1, b)
        universal_cluster = strict_gap > Fraction(2, 7)
        rows.append((b, mid5, mid8, counts, universal_cluster))
    print(
        f"{'b':>3} {'C_b(s=5)':>10} {'C_b(s=8)':>10} "
        f"{'1/b>2/7':>9} {'sequence s=0..13':>44}"
    )
    for b, mid5, mid8, counts, universal in rows:
        if b <= 14 or b in (15, 21, 30):
            print(
                f"{b:>3} {mid5:>10} {mid8:>10} {str(universal):>9} "
                f"{fmt_seq(counts):>44}"
            )
    print(
        "\nOnly b=2,3 are universal cluster centers.  Denominators 5,7,10,14,15,21,30 "
        "still give exact modular addresses for recurring intervals; using them "
        "as if they were universal witnesses is the overclaim to avoid."
    )


def print_prime_mask_inventory() -> None:
    section("DIVISOR-PROFILE MASKS ON {1,...,13}")
    hist = Counter(prime_mask(v) for v in UNIVERSE)
    for mask, count in sorted(hist.items()):
        labels = tuple(p for i, p in enumerate(PRIMES) if mask & (1 << i))
        print(f"  mask {mask:04b} divisors={labels!s:<14} count={count}")
    print(
        "\nThis mask inventory is the finite transfer state.  Adding a center prime "
        "is not a new ad hoc trick; it refines the state by one bit and updates "
        "the survivor polynomial by inclusion-exclusion."
    )


@dataclass(frozen=True)
class Route:
    name: str
    preserves_predicate: int
    proof_legitimacy: int
    residual_compression: int
    recurrence_value: int

    @property
    def score_tuple(self) -> tuple[int, int, int, int]:
        return (
            self.preserves_predicate,
            self.proof_legitimacy,
            self.residual_compression,
            self.recurrence_value,
        )


def tournament(routes: list[Route]) -> None:
    section("TOURNAMENT ANALYSIS: MODULAR PROOF OBLIGATIONS")
    n = len(routes)
    wins = {r.name: 0 for r in routes}
    edges: list[tuple[str, str]] = []
    for i, a in enumerate(routes):
        for j, b in enumerate(routes):
            if i >= j:
                continue
            if a.score_tuple >= b.score_tuple:
                winner, loser = a, b
            else:
                winner, loser = b, a
            wins[winner.name] += 1
            edges.append((winner.name, loser.name))
    print("routes and scores:")
    for r in routes:
        print(f"  {r.name:<28} score={r.score_tuple} wins={wins[r.name]}")
    cycles = 0
    route_names = [r.name for r in routes]
    edge_set = set(edges)
    for a, b, c in itertools.combinations(route_names, 3):
        if (
            ((a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set)
            or ((a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set)
        ):
            cycles += 1
    ordered = sorted(routes, key=lambda r: (wins[r.name], r.score_tuple), reverse=True)
    print(f"score histogram: {sorted(Counter(wins.values()).items())}")
    print(f"directed 3-cycles: {cycles}")
    print("Hamiltonian path:")
    print("  " + " > ".join(r.name for r in ordered))


def main() -> None:
    print("LRC(14) modular recurrence probe - codex S18")
    print_prime_mask_inventory()
    print_universal_center_recovery()
    print_center_denominator_table()
    routes = [
        Route("support6_coimage_tail", 4, 4, 3, 4),
        Route("mod210_address", 3, 2, 4, 4),
        Route("mod30_address", 3, 2, 3, 3),
        Route("mod6_universal", 3, 4, 2, 2),
        Route("primorial_height_escape", 2, 2, 4, 4),
        Route("all_small_denominators", 2, 1, 4, 3),
        Route("raw_runner_residues", 1, 1, 1, 1),
    ]
    tournament(routes)
    section("SESSION READOUT")
    print(
        "Mod 6 is the proven universal-center skeleton.  Mod 30 is the first "
        "non-universal address extension, because denominator 5 joins parity "
        "and triadic data but cannot by itself certify the cluster gap.  The "
        "LRC14 support-six tail forces mod 7, so the natural combined address is "
        "mod 210 or, better, the squarefree divisor-profile mask on {2,3,5,7}."
    )
    print(
        "The hidden recurrence is therefore not a one-dimensional recurrence in "
        "runners.  It is the finite transfer/inclusion-exclusion recurrence of "
        "prime masks, coupled to a signed mod-7 coimage tail.  This explains why "
        "mod 30 looked real and why it was not enough."
    )


if __name__ == "__main__":
    main()
