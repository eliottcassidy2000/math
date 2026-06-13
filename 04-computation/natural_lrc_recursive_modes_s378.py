#!/usr/bin/env python3
"""
natural_lrc_recursive_modes_s378.py

codex-2026-05-31 S378

Bridge the natural-number operation-mode graph thread with the Lonely Runner
frontier.

Natural operation modes:
  addition:       x -> z and y -> z when x + y = z
  multiplication: x -> z and y -> z when x * y = z

After forgetting the second input, addition completes to the total order while
multiplication keeps the sparse divisibility DAG.  Product-sum equations

  x_1 + ... + x_k = x_1 * ... * x_k

are finite critical pairs between those two operation modes.

Lonely Runner analogue:
  for k speeds, the threshold is 1/(k+1).  The initial segment is the
  additive/Dirichlet equality spine; endpoint protection and composite gates
  are controlled by divisibility in k+1.  This script turns that slogan into a
  small feature extractor and recursive n -> n+1 ledger.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()


@dataclass(frozen=True)
class ProductSumMinimum:
    k: int
    product: int
    seed: tuple[int, ...]
    ones: int

    @property
    def seed_len(self) -> int:
        return len(self.seed)

    @property
    def additive_discount(self) -> int:
        return sum(self.seed) - self.seed_len

    def compact(self) -> str:
        if self.ones == 0:
            return "*".join(map(str, self.seed))
        return f"1^{self.ones}+" + "*".join(map(str, self.seed))


@dataclass(frozen=True)
class OperationShadowStats:
    n: int
    additive_edges: int
    multiplicative_edges: int
    multiplicative_density: Fraction
    hasse_edges: int
    richest_add_endpoint: tuple[int, int]
    richest_mul_endpoint: tuple[int, int]


@dataclass(frozen=True)
class SpeedSetFeature:
    label: str
    n: int
    k: int
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap_ratio: Fraction
    components: int
    boundary_witnesses: int
    boundary_modulus: int
    endpoint_count: int
    interval_count: int
    unprotected_count: int
    peel_depth: int
    core_endpoint_count: int
    unit_skeleton: bool
    quotient_layer: str
    speeds_divisible_by_n: int
    nonunit_residues_mod_n: int
    additive_gates: int
    multiplicative_gates: int
    divisor_edges: int
    gcd_value: int


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_ratio(x: Fraction) -> str:
    return f"{float(x):.6f}"


def divisors(n: int) -> list[int]:
    out: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


def tau(n: int) -> int:
    return len(divisors(n))


def phi(n: int) -> int:
    count = 0
    for a in range(1, n + 1):
        if gcd(a, n) == 1:
            count += 1
    return count


def prime_factors(n: int) -> tuple[int, ...]:
    factors: list[int] = []
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            factors.append(d)
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        factors.append(m)
    return tuple(factors)


def is_prime(n: int) -> bool:
    return n >= 2 and phi(n) == n - 1


def prime_sieve(n: int) -> list[bool]:
    primes = [False, False] + [True] * max(0, n - 1)
    for p in range(2, isqrt(n) + 1):
        if primes[p]:
            for m in range(p * p, n + 1, p):
                primes[m] = False
    return primes


def operation_shadow_stats(n: int) -> OperationShadowStats:
    add_edges = n * (n - 1) // 2
    mul_edges = sum(max(0, tau(z) - 2) for z in range(1, n + 1))
    primes = prime_sieve(n)
    hasse_edges = sum(
        1
        for x in range(1, n + 1)
        for p in range(2, n // x + 1)
        if primes[p]
    )
    add_fibers = {
        z: z // 2
        for z in range(2, n + 1)
    }
    mul_fibers = {
        z: sum(1 for d in divisors(z) if d <= z // d and d * d <= z)
        for z in range(2, n + 1)
    }
    richest_add = max(add_fibers.items(), key=lambda item: (item[1], item[0]))
    richest_mul = max(mul_fibers.items(), key=lambda item: (item[1], item[0]))
    return OperationShadowStats(
        n=n,
        additive_edges=add_edges,
        multiplicative_edges=mul_edges,
        multiplicative_density=Fraction(mul_edges, add_edges),
        hasse_edges=hasse_edges,
        richest_add_endpoint=richest_add,
        richest_mul_endpoint=richest_mul,
    )


def enumerate_product_sum_minima(max_k: int) -> dict[int, ProductSumMinimum]:
    max_product = 2 * max_k
    by_k: dict[int, list[ProductSumMinimum]] = defaultdict(list)

    def rec(start: int, product: int, total: int, seed: tuple[int, ...]) -> None:
        if len(seed) >= 2 and product >= total:
            ones = product - total
            k = len(seed) + ones
            if 2 <= k <= max_k:
                by_k[k].append(
                    ProductSumMinimum(
                        k=k,
                        product=product,
                        seed=seed,
                        ones=ones,
                    )
                )

        limit = max_product // product
        for f in range(start, limit + 1):
            rec(f, product * f, total + f, seed + (f,))

    rec(2, 1, 0, tuple())
    minima = {}
    for k in range(2, max_k + 1):
        minima[k] = min(
            by_k[k],
            key=lambda row: (row.product, row.seed_len, row.seed, row.ones),
        )
    return minima


def operation_profile(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    speed_set = set(speeds)
    additive = 0
    multiplicative = 0
    for x, y in combinations(speeds, 2):
        if x + y in speed_set:
            additive += 1
        if x * y in speed_set:
            multiplicative += 1
    for x in speeds:
        if 2 * x in speed_set:
            additive += 1
        if x * x in speed_set:
            multiplicative += 1

    divisor_edges = sum(
        1
        for x in speeds
        for z in speeds
        if x < z and z % x == 0
    )
    return additive, multiplicative, divisor_edges


def classify_report(row) -> str:
    if row.max_gap > 0:
        return "positive_gap"
    if row.boundary_witness_count:
        return "boundary_only"
    return "open_cover_candidate"


def speed_set_feature(label: str, raw_speeds: tuple[int, ...]) -> SpeedSetFeature:
    report = S356.report(label, list(raw_speeds))
    descent = S362.summarize(list(report.speeds))
    n = len(report.speeds) + 1
    g = 0
    for v in report.speeds:
        g = gcd(g, v)
    additive, multiplicative, divisor_edges = operation_profile(report.speeds)
    residues = {v % n for v in report.speeds}
    nonunit_residues = sum(1 for r in residues if r != 0 and gcd(r, n) > 1)
    return SpeedSetFeature(
        label=label,
        n=n,
        k=len(report.speeds),
        speeds=report.speeds,
        classification=classify_report(report),
        forbidden_length=report.forbidden_length,
        max_gap_ratio=(
            Fraction(0) if report.threshold == 0 else report.max_gap / report.threshold
        ),
        components=report.components,
        boundary_witnesses=report.boundary_witness_count,
        boundary_modulus=report.boundary_modulus,
        endpoint_count=descent.endpoint_count,
        interval_count=descent.interval_count,
        unprotected_count=descent.unprotected_count,
        peel_depth=len(descent.peel_layers),
        core_endpoint_count=descent.core_endpoint_count,
        unit_skeleton=descent.unit_skeleton,
        quotient_layer=descent.quotient_layer,
        speeds_divisible_by_n=sum(1 for v in report.speeds if v % n == 0),
        nonunit_residues_mod_n=nonunit_residues,
        additive_gates=additive,
        multiplicative_gates=multiplicative,
        divisor_edges=divisor_edges,
        gcd_value=g,
    )


def initial_segment_features(max_n: int) -> list[SpeedSetFeature]:
    return [
        speed_set_feature(f"initial n={n}", tuple(range(1, n)))
        for n in range(4, max_n + 1)
    ]


def print_operation_shadow_table() -> None:
    print("MODE 1: NATURAL OPERATION SHADOWS")
    print("=" * 86)
    print(
        "Addition becomes the complete transitive order after forgetting the "
        "second input; multiplication remains the divisor skeleton."
    )
    print()
    header = (
        "N   add_edges   mul_edges   mul/add    hasse_edges   "
        "rich_add(z,fib)   rich_mul(z,fib)"
    )
    print(header)
    print("-" * len(header))
    for n in [20, 50, 100, 250, 500, 1000]:
        row = operation_shadow_stats(n)
        print(
            f"{n:>4} {row.additive_edges:>11} {row.multiplicative_edges:>11} "
            f"{fmt_ratio(row.multiplicative_density):>9} {row.hasse_edges:>14} "
            f"{row.richest_add_endpoint!s:>17} {row.richest_mul_endpoint!s:>17}"
        )
    print()


def print_product_sum_lrc_alignment(minima: dict[int, ProductSumMinimum]) -> None:
    print("MODE 2: PRODUCT-SUM MINIMA ALIGNED WITH LRC RUNNER COUNT")
    print("=" * 86)
    print(
        "Here k is both product-sum arity and the number of moving runners; "
        "n=k+1 is the LRC threshold denominator."
    )
    print()
    header = (
        "n  k  phi(n) tau(n)  m(k) seed_len ones discount  minimal_seed"
    )
    print(header)
    print("-" * len(header))
    for n in range(4, 25):
        k = n - 1
        row = minima[k]
        print(
            f"{n:>2} {k:>2} {phi(n):>6} {tau(n):>6} "
            f"{row.product:>5} {row.seed_len:>8} {row.ones:>4} "
            f"{row.additive_discount:>8}  {row.compact()}"
        )
    print()


def print_initial_segment_recursion(rows: list[SpeedSetFeature]) -> None:
    print("MODE 3: INITIAL-SEGMENT LRC RECURSION")
    print("=" * 86)
    print(
        "Initial segments are the Dirichlet equality spine.  The useful n "
        "state is not just tight/untight: it has unit skeleton, divisor gates, "
        "endpoint peeling, and operation-closure coordinates."
    )
    print()
    header = (
        "n  phi  phi/(n-1) tau  endpoints intervals comps bdyW  Q  "
        "peel core addG mulG divE"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        print(
            f"{row.n:>2} {phi(row.n):>4} "
            f"{phi(row.n)/(row.n-1):>9.4f} {tau(row.n):>3} "
            f"{row.endpoint_count:>9} {row.interval_count:>9} "
            f"{row.components:>5} {row.boundary_witnesses:>4} "
            f"{row.boundary_modulus:>3} {row.peel_depth:>5} "
            f"{row.core_endpoint_count:>4} {row.additive_gates:>4} "
            f"{row.multiplicative_gates:>4} {row.divisor_edges:>4}"
        )
    print()


def print_transition_ledger(
    rows: list[SpeedSetFeature],
    minima: dict[int, ProductSumMinimum],
) -> None:
    print("MODE 4: n -> n+1 TRANSITION LEDGER")
    print("=" * 86)
    print(
        "The recursive changes are arithmetical rather than smooth.  Prime n "
        "maximizes the unit skeleton; composite n opens quotient/descent "
        "channels through divisibility."
    )
    print()
    header = "n  type       factors      d_phi d_tau d_peel  dm(k) seed_change"
    print(header)
    print("-" * len(header))
    by_n = {row.n: row for row in rows}
    for n in range(5, max(by_n) + 1):
        prev = by_n[n - 1]
        cur = by_n[n]
        prev_min = minima[n - 2]
        cur_min = minima[n - 1]
        seed_change = cur_min.seed != prev_min.seed
        factors = ".".join(map(str, prime_factors(n))) or "1"
        print(
            f"{n:>2} {'prime' if is_prime(n) else 'comp':<10} {factors:<12} "
            f"{phi(n)-phi(n-1):>5} {tau(n)-tau(n-1):>5} "
            f"{cur.peel_depth-prev.peel_depth:>6} "
            f"{cur_min.product-prev_min.product:>6} "
            f"{'yes' if seed_change else 'no'}"
        )
    print()


def print_hard_family_features() -> None:
    print("MODE 5: FEATURE VECTORS FOR TIGHT AND NEAR-DISPROOF FAMILIES")
    print("=" * 86)
    samples = [
        ("sporadic tight n=5", (1, 3, 4, 7)),
        ("sporadic tight n=6", (1, 3, 4, 5, 9)),
        ("initial n=14", tuple(range(1, 14))),
        (
            "n14 seven-ladder",
            (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91),
        ),
        (
            "n14 2x7 ladder",
            (2, 4, 6, 7, 8, 10, 12, 14, 21, 28, 35, 42, 49),
        ),
        ("initial n=15", tuple(range(1, 15))),
        (
            "n15 3x5 ladder",
            (3, 5, 6, 9, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50),
        ),
        (
            "n15 mixed gates",
            (1, 3, 5, 6, 9, 10, 14, 15, 20, 25, 30, 45, 60, 75),
        ),
    ]
    rows = [speed_set_feature(label, speeds) for label, speeds in samples]

    header = (
        "label                  n class          gap/thresh forb_len comps "
        "bdyW unprot peel core divN nonunit addG mulG divE"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        print(
            f"{row.label:<22} {row.n:>2} {row.classification:<14} "
            f"{fmt_ratio(row.max_gap_ratio):>10} {fmt_frac(row.forbidden_length):>8} "
            f"{row.components:>5} {row.boundary_witnesses:>4} "
            f"{row.unprotected_count:>7} {row.peel_depth:>4} "
            f"{row.core_endpoint_count:>4} {row.speeds_divisible_by_n:>4} "
            f"{row.nonunit_residues_mod_n:>7} {row.additive_gates:>4} "
            f"{row.multiplicative_gates:>4} {row.divisor_edges:>4}"
        )

    print()
    print("Feature notes:")
    for row in rows:
        protection_pressure = Fraction(
            row.endpoint_count - row.unprotected_count,
            max(1, row.endpoint_count),
        )
        divisor_over_add = Fraction(row.divisor_edges, max(1, row.additive_gates))
        print(
            f"  {row.label}: protect_pressure={fmt_ratio(protection_pressure)} "
            f"div/add={fmt_ratio(divisor_over_add)} "
            f"unit_skeleton={row.unit_skeleton} quotient={row.quotient_layer}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 86)
    print(
        "1. The old incomplete-tournament idea is correct only before "
        "forgetting labels.  The additive one-shadow is already the full "
        "order, so its real data are operation fibers and critical pairs."
    )
    print(
        "2. Multiplication is the sparse residue: divisibility survives the "
        "shadow projection.  Product-sum minima are the finite toy model where "
        "additive completion and multiplicative skeleton meet."
    )
    print(
        "3. LRC has the same grammar.  Initial segments are the additive "
        "Dirichlet equality spine; endpoint protection is forced through the "
        "multiplicative/divisibility channel of n=k+1."
    )
    print(
        "4. A useful LRC-TDA feature vector should include: max_gap_ratio, "
        "boundary witness count, endpoint peel depth, unit skeleton size "
        "phi(n), divisor-gate count, nonunit residue count, and operation "
        "closure inside the speed set."
    )
    print(
        "5. The n-recursion is not monotone in n; it jumps at prime/composite "
        "changes in phi(n), tau(n), and factor-packing type.  That is the "
        "natural-number analogue of the repo's tournament recursion story: "
        "the state needs boundary coordinates, not only a scalar invariant."
    )


def main() -> None:
    print("Natural operation modes and LRC recursive metrics (codex-2026-05-31 S378)")
    print()
    minima = enumerate_product_sum_minima(64)
    initial_rows = initial_segment_features(22)
    print_operation_shadow_table()
    print_product_sum_lrc_alignment(minima)
    print_initial_segment_recursion(initial_rows)
    print_transition_ledger(initial_rows, minima)
    print_hard_family_features()
    print_synthesis()


if __name__ == "__main__":
    main()
