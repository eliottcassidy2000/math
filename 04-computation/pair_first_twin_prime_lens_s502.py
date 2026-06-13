#!/usr/bin/env python3
"""
pair_first_twin_prime_lens_s502.py

codex-2026-06-01 S502

Think of everything in terms of pairs.

The primitive object is not a vertex or a number by itself.  It is a pair-cell:

    endpoints + relation/label + fiber coordinates.

For tournaments the pair-cell is an unordered pair {i,j} with one orientation
bit.  For prime-pair problems the pair-cell is {a,b} with primality at both
endpoints and labels

    sum = a+b, gap = b-a, midpoint = (a+b)/2.

Twin primes are therefore the fixed-gap row gap=2 in the prime-pair surface.
Goldbach is the fixed-sum column sum=N in the same surface.  This script makes
that duality visible and compares it with the triangular edge board of a
tournament.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt, log


LIMIT = 100_000
MAX_GAP = 60


def prime_sieve(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (((limit - start) // p) + 1)
    return [n for n in range(2, limit + 1) if sieve[n]]


PRIMES = prime_sieve(LIMIT + MAX_GAP)
PRIME_SET = set(PRIMES)


def factor_distinct(n: int) -> tuple[int, ...]:
    factors: list[int] = []
    remaining = n
    for p in PRIMES:
        if p * p > remaining:
            break
        if remaining % p == 0:
            factors.append(p)
            while remaining % p == 0:
                remaining //= p
    if remaining > 1:
        factors.append(remaining)
    return tuple(factors)


def twin_prime_constant_approx() -> float:
    value = 1.0
    for p in PRIMES:
        if p > LIMIT:
            break
        if p > 2:
            value *= 1.0 - 1.0 / ((p - 1) ** 2)
    return value


C2_APPROX = twin_prime_constant_approx()


def singular_gap_multiplier(gap: int) -> float:
    """Relative Hardy-Littlewood multiplier for even prime-pair gap rows."""

    if gap % 2:
        return 0.0
    multiplier = 1.0
    for p in factor_distinct(gap):
        if p > 2:
            multiplier *= (p - 1) / (p - 2)
    return multiplier


def expected_gap_count(gap: int) -> float:
    if gap % 2:
        return 0.0
    # We sum only over odd lower endpoints.  The standard 2*C2 constant is
    # written over all integer starts, so the conditional odd-row density gets
    # the extra parity factor here.
    local = 4.0 * C2_APPROX * singular_gap_multiplier(gap)
    density_sum = 0.0
    for n in range(3, LIMIT - gap + 1, 2):
        density_sum += 1.0 / (log(n) * log(n + gap))
    return local * density_sum


def prime_pair_count_for_gap(gap: int) -> int:
    return sum(1 for p in PRIMES if p <= LIMIT - gap and p + gap in PRIME_SET)


def goldbach_pairs(n: int) -> list[tuple[int, int]]:
    return [(p, n - p) for p in PRIMES if p <= n - p and n - p in PRIME_SET]


def prime_pair_edges(max_value: int, max_gap: int) -> list[tuple[int, int]]:
    edges = []
    for p in PRIMES:
        if p == 2:
            continue
        if p > max_value:
            break
        for q in PRIMES:
            if q <= p or q == 2:
                continue
            if q > max_value or q - p > max_gap:
                break
            edges.append((p, q))
    return edges


def fibs_upto(n: int) -> list[int]:
    fibs = [1, 2]
    while fibs[-1] < n:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs


FIBS = fibs_upto(2 * LIMIT + MAX_GAP)


def zeckendorf_indices(n: int) -> tuple[int, ...]:
    out = []
    remaining = n
    for idx in range(len(FIBS) - 1, -1, -1):
        fib = FIBS[idx]
        if fib <= remaining:
            out.append(idx + 1)
            remaining -= fib
        if remaining == 0:
            break
    return tuple(sorted(out))


@dataclass(frozen=True)
class PairCarry:
    pair: tuple[int, int]
    raw_terms: int
    target_terms: int
    repeats: int
    adjacencies: int
    l1: int
    score: int
    raw_support: tuple[int, ...]
    target_support: tuple[int, ...]


def carry_for_pair(pair: tuple[int, int]) -> PairCarry:
    raw = Counter()
    for endpoint in pair:
        raw.update(zeckendorf_indices(endpoint))
    target = Counter(zeckendorf_indices(sum(pair)))
    all_indices = set(raw) | set(target)
    repeats = sum(max(0, raw[i] - 1) for i in raw)
    raw_support = tuple(sorted(i for i, count in raw.items() if count > 0))
    adjacencies = sum(1 for i in raw_support if raw.get(i + 1, 0) > 0)
    l1 = sum(abs(raw.get(i, 0) - target.get(i, 0)) for i in all_indices)
    return PairCarry(
        pair=pair,
        raw_terms=sum(raw.values()),
        target_terms=sum(target.values()),
        repeats=repeats,
        adjacencies=adjacencies,
        l1=l1,
        score=l1 + repeats + adjacencies,
        raw_support=raw_support,
        target_support=tuple(sorted(target)),
    )


def tournament_range_counts(n: int) -> list[tuple[int, int, float]]:
    total = n * (n - 1) // 2
    return [(gap, n - gap, (n - gap) / total) for gap in range(1, n)]


def print_header(title: str) -> None:
    print()
    print("=" * 112)
    print(title)
    print("=" * 112)


def section_pair_ontology() -> None:
    print_header("Pair-first ontology")
    rows = [
        (
            "tournament",
            "{i,j}",
            "orientation sign s_ij",
            "range j-i, midpoint i+j, flip impact",
        ),
        (
            "twin primes",
            "{p,p+2}",
            "both endpoints prime",
            "fixed gap row g=2",
        ),
        (
            "Goldbach",
            "{p,N-p}",
            "both endpoints prime",
            "fixed sum column s=N",
        ),
        (
            "SC blowup",
            "{u_r,v_s}",
            "lane/cross orientation",
            "four edge-cells between each old pair",
        ),
        (
            "LRC distances",
            "{runner i, runner j}",
            "chosen metric switch",
            "time-varying chord/gap/pressure label",
        ),
    ]
    print(f"{'object':<16} {'pair-cell':<18} {'relation':<26} labels")
    for row in rows:
        print(f"{row[0]:<16} {row[1]:<18} {row[2]:<26} {row[3]}")
    print()
    print("Synthesis:")
    print("  A vertex is a carrier.  The connection between two carriers is the coordinate.")
    print("  Twin primes are the recurrence of the smallest nontrivial prime-pair edge label.")
    print("  Goldbach and twin primes are perpendicular slices of the same prime-pair surface.")


def section_gap_rows() -> None:
    print_header("Twin primes as the gap=2 row of the prime-pair surface")
    twin_count = prime_pair_count_for_gap(2)
    twin_expected = expected_gap_count(2)
    print(f"Prime-pair rows with p,p+gap <= {LIMIT}")
    print(f"C2 approximation from primes <= {LIMIT}: {C2_APPROX:.9f}")
    print(
        f"{'gap':>4} {'count':>7} {'HL-shape':>9} "
        f"{'count/twin':>11} {'shape/twin':>11} {'odd factors':>14}"
    )
    for gap in range(2, MAX_GAP + 1, 2):
        count = prime_pair_count_for_gap(gap)
        expected = expected_gap_count(gap)
        factors = ",".join(str(p) for p in factor_distinct(gap) if p > 2) or "-"
        print(
            f"{gap:4d} {count:7d} {expected:9.1f} "
            f"{count / twin_count:11.3f} {expected / twin_expected:11.3f} "
            f"{factors:>14}"
        )
    print()
    print("Read this as edge-row geometry, not isolated prime numerology:")
    print("  multiplying the gap by a new odd prime removes fewer residue classes,")
    print("  so the whole row gets a Hardy-Littlewood singular-series boost.")


def section_sum_gap_duality() -> None:
    print_header("Goldbach/twin duality in (sum, gap) coordinates")
    print("For same-parity endpoints, the coordinate change is invertible:")
    print("  endpoints = ((sum-gap)/2, (sum+gap)/2)")
    print("  twin primes: fixed gap=2")
    print("  Goldbach:    fixed even sum=N")
    print()
    print(f"{'N':>7} {'Goldbach pairs':>15} {'gap range':>15} {'best carry pair':>18} {'carry':>7}")
    for n in (100, 500, 1_000, 5_000, 10_000, 50_000, 100_000):
        pairs = goldbach_pairs(n)
        gaps = [q - p for p, q in pairs]
        best = min((carry_for_pair(pair) for pair in pairs), key=lambda audit: (audit.score, abs(audit.pair[1] - audit.pair[0]), audit.pair))
        gap_range = f"{min(gaps)}..{max(gaps)}" if gaps else "NONE"
        print(f"{n:7d} {len(pairs):15d} {gap_range:>15} {str(best.pair):>18} {best.score:7d}")
    print()
    print("  Fixed-sum columns usually contain many possible gaps; fixed-gap rows")
    print("  usually contain many possible sums.  The same edge can be read both ways.")


def section_pair_surface_sample() -> None:
    print_header("A small prime-pair surface sample")
    edges = prime_pair_edges(400, 40)
    by_gap = Counter(q - p for p, q in edges)
    by_sum_mod = Counter((p + q) % 30 for p, q in edges)
    by_mid_residue = Counter(((p + q) // 2) % 30 for p, q in edges)
    print(f"Prime pair edges with endpoints <= 400 and gap <= 40: {len(edges)}")
    print("Gap row counts:")
    print("  " + " ".join(f"{gap}:{by_gap[gap]}" for gap in sorted(by_gap)))
    print("Most common sum residues mod 30:")
    print("  " + " ".join(f"{res}:{count}" for res, count in by_sum_mod.most_common(10)))
    print("Most common midpoint residues mod 30:")
    print("  " + " ".join(f"{res}:{count}" for res, count in by_mid_residue.most_common(10)))
    print()
    print("The residue data is a reminder that the pair-cell carries local obstruction")
    print("data.  A prime-pair edge is an arithmetic connection, not just two prime nodes.")


def section_tournament_board() -> None:
    print_header("Tournament edge board as the finite unsieved analogue")
    for n in (14, 18):
        counts = tournament_range_counts(n)
        print(f"n={n}: total pair-cells C(n,2)={n * (n - 1) // 2}")
        print("  range rows: " + " ".join(f"{gap}:{count}" for gap, count, _share in counts))
        print(
            "  entropy by range is triangular before any orientation; "
            "Tournament Analysis then adds one sign bit per cell."
        )
    print()
    print("Analogy:")
    print("  tournament board: all pair-cells exist, then orientation signs decorate them.")
    print("  prime-pair board: all integer pair-cells exist, then the sieve deletes most cells.")
    print("  LRC board: all runner pair-cells exist, then time changes their metric labels.")


def section_sc_twin_gain() -> None:
    print_header("SC blowup: pair cells become twin lanes")
    print("For each old tournament pair {u,v}, SC blowup creates four edge-cells:")
    print("  lane:  u0-v0, u1-v1 follow T")
    print("  cross: u0-v1, u1-v0 follow T^op")
    print()
    print("This is why 'twin' is the right word inside the repo too:")
    print("  the old pair is not discarded; it is doubled into lane/cross structure.")
    print("  Score variation is erased at the vertex level, but remembered in pair labels.")
    print("  That mirrors twin primes: the node property prime is less informative than")
    print("  the repeated connection type gap=2 surviving the sieve.")


def section_carry_debt() -> None:
    print_header("Zeckendorf carry debt on prime-pair edges")
    twin_edges = [(p, p + 2) for p in PRIMES if p <= 1_000 and p + 2 in PRIME_SET]
    selected_twins = [twin_edges[i] for i in (0, 1, 2, 5, 10, 20, len(twin_edges) - 1)]
    print("Twin prime edge audits:")
    print(f"{'pair':>12} {'sum':>6} {'raw support':>22} {'target':>18} {'score':>6}")
    for pair in selected_twins:
        audit = carry_for_pair(pair)
        print(
            f"{str(pair):>12} {sum(pair):6d} "
            f"{str(audit.raw_support):>22} {str(audit.target_support):>18} {audit.score:6d}"
        )
    print()
    print("Best carry Goldbach edge in selected columns:")
    print(f"{'sum':>7} {'pair':>18} {'raw support':>22} {'target':>18} {'score':>6}")
    for n in (100, 500, 1_000, 5_000, 10_000, 50_000):
        audits = [carry_for_pair(pair) for pair in goldbach_pairs(n)]
        best = min(audits, key=lambda audit: (audit.score, abs(audit.pair[1] - audit.pair[0]), audit.pair))
        print(
            f"{n:7d} {str(best.pair):>18} "
            f"{str(best.raw_support):>22} {str(best.target_support):>18} {best.score:6d}"
        )
    print()
    print("This gives a repo-native metric for prime edges:")
    print("  how much local carry work is needed to collapse endpoint supports to the")
    print("  support of the sum.  It is pair data, not endpoint data.")


def section_questions() -> None:
    print_header("New questions suggested by the pair-first lens")
    questions = [
        "Can Tournament Analysis store every pair-cell as (range, midpoint, sign, flip-impact) and only later project to vertex scores?",
        "Can prime-pair work store every edge as (sum, gap, residue singular factor, carry debt) and compare fixed-sum vs fixed-gap fibers?",
        "Is the right Goldbach/twin bridge a TDA of the prime-pair surface, with holes created by local residue deletions and filled by singular-series boosts?",
        "For LRC, do pressure DAGs become clearer if their nodes are runner-pairs or endpoint-runner pairs rather than runners alone?",
        "Can SC blowup be used as the finite tournament model of an arithmetic twin operation: duplicate carriers, keep connection memory, flatten node scores?",
    ]
    for i, question in enumerate(questions, 1):
        print(f"{i}. {question}")


def main() -> None:
    section_pair_ontology()
    section_gap_rows()
    section_sum_gap_duality()
    section_pair_surface_sample()
    section_tournament_board()
    section_sc_twin_gain()
    section_carry_debt()
    section_questions()


if __name__ == "__main__":
    main()
