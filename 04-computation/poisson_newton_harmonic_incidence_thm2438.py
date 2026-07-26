#!/usr/bin/env python3
"""Exact companion for THM-2438.

All truth-bearing checks use ``require`` so that ``python3`` and
``python3 -O`` execute the same audit.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb, factorial, isqrt


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def forward_difference_at_zero(values: list[int], order: int) -> int:
    row = list(values[: order + 1])
    require(len(row) == order + 1, "insufficient difference row")
    for _ in range(order):
        row = [row[i + 1] - row[i] for i in range(len(row) - 1)]
    return row[0]


def central_trinomial(k: int) -> int:
    return sum(
        factorial(k) // (factorial(j) ** 2 * factorial(k - 2 * j))
        for j in range(k // 2 + 1)
    )


def truncated_egf(values: list[int], degree: int) -> list[Fraction]:
    return [Fraction(values[n], factorial(n)) for n in range(degree + 1)]


def multiply_series(
    left: list[Fraction], right: list[Fraction], degree: int
) -> list[Fraction]:
    out = [Fraction(0) for _ in range(degree + 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j <= degree:
                out[i + j] += a * b
    return out


def factor_pairs(n: int, strict: bool) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for a in range(1, isqrt(n) + 1):
        if n % a:
            continue
        b = n // a
        if a < b or (not strict and a == b):
            pairs.append((a, b))
    return pairs


def fibre_loss(n: int, holes: frozenset[int], strict: bool) -> int:
    ambient = factor_pairs(n, strict)
    retained = [(a, b) for a, b in ambient if a not in holes and b not in holes]
    return len(ambient) - len(retained)


def marked_divisor_incidence(n: int, holes: frozenset[int]) -> int:
    return sum(1 for h in holes if n % h == 0)


def distinct_hole_collisions(n: int, holes: frozenset[int]) -> int:
    return sum(1 for a, b in combinations(sorted(holes), 2) if a * b == n)


def diagonal_hole(n: int, holes: frozenset[int]) -> int:
    r = isqrt(n)
    return int(r * r == n and r in holes)


def additive_pairs(n: int, strict: bool) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for a in range(1, n):
        b = n - a
        if a < b or (not strict and a == b):
            pairs.append((a, b))
    return pairs


def additive_loss(n: int, holes: frozenset[int], strict: bool) -> int:
    ambient = additive_pairs(n, strict)
    retained = [(a, b) for a, b in ambient if a not in holes and b not in holes]
    return len(ambient) - len(retained)


def sieve(limit: int) -> bytearray:
    prime = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        prime[0] = 0
    if limit >= 1:
        prime[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            start = p * p
            prime[start : limit + 1 : p] = b"\x00" * (
                (limit - start) // p + 1
            )
    return prime


def main() -> None:
    max_k = 60
    central = [comb(2 * n, n) for n in range(max_k + 1)]
    bulk = [4**n for n in range(max_k + 1)]
    weak = [(bulk[n] + central[n]) // 2 for n in range(max_k + 1)]
    strict_half = [(bulk[n] - central[n]) // 2 for n in range(max_k + 1)]
    shifted_strict = [strict_half[n + 1] for n in range(max_k)]
    trinomial = [central_trinomial(k) for k in range(max_k + 1)]

    difference_checks = 0
    for k in range(max_k):
        dc = forward_difference_at_zero(central, k)
        db = forward_difference_at_zero(bulk, k)
        dw = forward_difference_at_zero(weak, k)
        ds = forward_difference_at_zero(strict_half, k)
        dshift = forward_difference_at_zero(shifted_strict, k)
        require(dc == trinomial[k], f"central Newton coefficient failed at {k}")
        require(db == 3**k, f"bulk Newton coefficient failed at {k}")
        require(
            dw == (3**k + trinomial[k]) // 2,
            f"weak ternary half failed at {k}",
        )
        require(
            ds == (3**k - trinomial[k]) // 2,
            f"strict ternary half failed at {k}",
        )
        require(
            dshift == (4 * 3**k - trinomial[k] - trinomial[k + 1]) // 2,
            f"shifted strict formula failed at {k}",
        )
        difference_checks += 5

    ternary_checks = 0
    ternary_histogram: list[tuple[int, int, int]] = []
    for k in range(10):
        sums = [sum(word) for word in product((-1, 0, 1), repeat=k)]
        zero = sum(total == 0 for total in sums)
        nonnegative = sum(total >= 0 for total in sums)
        positive = sum(total > 0 for total in sums)
        require(zero == trinomial[k], f"zero-drift word count failed at {k}")
        require(
            nonnegative == (3**k + trinomial[k]) // 2,
            f"nonnegative ternary half failed at {k}",
        )
        require(
            positive == (3**k - trinomial[k]) // 2,
            f"positive ternary half failed at {k}",
        )
        ternary_histogram.append((nonnegative, positive, zero))
        ternary_checks += 3

    egf_degree = 24
    test_sequences = [
        [n**5 - 3 * n**3 + 7 * n + 11 for n in range(egf_degree + 1)],
        [2**n for n in range(egf_degree + 1)],
        [comb(2 * n, n) for n in range(egf_degree + 1)],
        [weak[n] for n in range(egf_degree + 1)],
    ]
    exp_minus = [Fraction((-1) ** n, factorial(n)) for n in range(egf_degree + 1)]
    egf_checks = 0
    for seq_index, seq in enumerate(test_sequences):
        ordinary_egf = truncated_egf(seq, egf_degree)
        transformed = multiply_series(exp_minus, ordinary_egf, egf_degree)
        differences = [
            Fraction(forward_difference_at_zero(seq, k), factorial(k))
            for k in range(egf_degree + 1)
        ]
        require(
            transformed == differences,
            f"Poisson--Newton transform failed for sequence {seq_index}",
        )
        egf_checks += egf_degree + 1

    # Coefficient proof of exp(2x) I_0(2x) and exp(x) I_0(2x).
    bessel_checks = 0
    for n in range(egf_degree + 1):
        central_egf_coefficient = sum(
            Fraction(2 ** (n - 2 * j), factorial(n - 2 * j) * factorial(j) ** 2)
            for j in range(n // 2 + 1)
        )
        ternary_egf_coefficient = sum(
            Fraction(1, factorial(n - 2 * j) * factorial(j) ** 2)
            for j in range(n // 2 + 1)
        )
        require(
            central_egf_coefficient == Fraction(comb(2 * n, n), factorial(n)),
            f"Bessel central coefficient failed at {n}",
        )
        require(
            ternary_egf_coefficient == Fraction(trinomial[n], factorial(n)),
            f"Bessel trinomial coefficient failed at {n}",
        )
        bessel_checks += 2

    hole_universe = tuple(range(1, 8))
    hole_sets = [
        frozenset(combination)
        for size in range(len(hole_universe) + 1)
        for combination in combinations(hole_universe, size)
    ]
    pointwise_checks = 0
    cumulative_checks = 0
    additive_checks = 0
    for holes in hole_sets:
        for n in range(1, 121):
            incidence = marked_divisor_incidence(n, holes)
            collision = distinct_hole_collisions(n, holes)
            diagonal = diagonal_hole(n, holes)
            require(
                fibre_loss(n, holes, False) == incidence - collision,
                f"weak collision tax failed for {holes}, {n}",
            )
            require(
                fibre_loss(n, holes, True) == incidence - collision - diagonal,
                f"strict collision tax failed for {holes}, {n}",
            )
            pointwise_checks += 2

        if holes:
            cutoff = max(2 * max(holes) + 1, 2)
            for n in range(cutoff, cutoff + 25):
                require(
                    additive_loss(n, holes, False) == len(holes),
                    f"weak additive tail failed for {holes}, {n}",
                )
                require(
                    additive_loss(n, holes, True) == len(holes),
                    f"strict additive tail failed for {holes}, {n}",
                )
                additive_checks += 2

            n_limit = max(200, max(holes) ** 2)
            weak_sum = sum(fibre_loss(n, holes, False) for n in range(1, n_limit + 1))
            strict_sum = sum(fibre_loss(n, holes, True) for n in range(1, n_limit + 1))
            floor_sum = sum(n_limit // h for h in holes)
            m = len(holes)
            require(
                weak_sum == floor_sum - m * (m - 1) // 2,
                f"finite weak summatory formula failed for {holes}",
            )
            require(
                strict_sum == floor_sum - m * (m + 1) // 2,
                f"finite strict summatory formula failed for {holes}",
            )
            cumulative_checks += 2

    partial_summation_checks = 0
    support_controls = [
        frozenset(range(1, 80, 3)),
        frozenset(n * (n + 1) // 2 for n in range(1, 18)),
        frozenset({1, 4, 6, 12, 18, 30, 42}),
    ]
    for support in support_controls:
        for n_limit in range(1, 150):
            harmonic = sum(
                (Fraction(1, h) for h in support if h <= n_limit), Fraction(0)
            )
            counts = [
                sum(1 for h in support if h <= n) for n in range(n_limit + 1)
            ]
            abel = Fraction(counts[n_limit], n_limit)
            abel += sum(
                (Fraction(counts[n], n * (n + 1)) for n in range(1, n_limit)),
                Fraction(0),
            )
            require(harmonic == abel, f"partial summation failed at {n_limit}")

            incidence_sum = sum(
                marked_divisor_incidence(n, support) for n in range(1, n_limit + 1)
            )
            floor_identity = sum(n_limit // h for h in support if h <= n_limit)
            require(
                incidence_sum == floor_identity,
                f"divisor-incidence mean identity failed at {n_limit}",
            )
            partial_summation_checks += 2

    triangular_checks = 0
    triangular_partial = Fraction(0)
    for k in range(1, 501):
        triangular = k * (k + 1) // 2
        triangular_partial += Fraction(1, triangular)
        require(
            triangular_partial == Fraction(2 * k, k + 1),
            f"triangular reciprocal telescoping failed at {k}",
        )
        triangular_checks += 1

    prime_limit = 1_000_000
    prime = sieve(prime_limit + 1)
    twin_centres = [
        c for c in range(2, prime_limit + 1) if prime[c - 1] and prime[c + 1]
    ]
    require(twin_centres[:12] == [4, 6, 12, 18, 30, 42, 60, 72, 102, 108, 138, 150],
            "A014574 startup mismatch")
    twin_identity_checks = 0
    for c in twin_centres:
        require(
            Fraction(1, c)
            == Fraction(1, 2) * (Fraction(1, c - 1) + Fraction(1, c + 1))
            - Fraction(1, c * (c * c - 1)),
            f"twin-centre reciprocal identity failed at {c}",
        )
        twin_identity_checks += 1

    # Direct small-N check that reciprocal support is the normalized
    # marked-divisor incidence, before taking any limit.
    twin_small = frozenset(c for c in twin_centres if c <= 5000)
    twin_mean_limit = 20_000
    twin_incidence_sum = sum(
        marked_divisor_incidence(n, twin_small)
        for n in range(1, twin_mean_limit + 1)
    )
    twin_floor_sum = sum(twin_mean_limit // c for c in twin_small)
    require(twin_incidence_sum == twin_floor_sum, "twin divisor-incidence identity")

    print("THM-2438 exact verification")
    print(f"Newton/difference identities: {difference_checks}")
    print(f"ternary word identities: {ternary_checks}")
    print(f"formal EGF coefficients: {egf_checks}")
    print(f"Bessel/central coefficients: {bessel_checks}")
    print(f"operation pointwise identities: {pointwise_checks}")
    print(f"operation cumulative identities: {cumulative_checks}")
    print(f"additive tail identities: {additive_checks}")
    print(f"Abel/incidence identities: {partial_summation_checks}")
    print(f"triangular telescoping identities: {triangular_checks}")
    print(f"A014574 centres through {prime_limit}: {len(twin_centres)}")
    print(f"A014574 reciprocal identities: {twin_identity_checks}")
    print(f"ternary weak/strict/zero k=0..9: {ternary_histogram}")
    print(f"triangular reciprocal partial k=500: {triangular_partial}")
    print(
        "twin-centre marked incidence through "
        f"{twin_mean_limit}: {twin_incidence_sum}/{twin_mean_limit}"
    )
    print("all exact checks passed")


if __name__ == "__main__":
    main()
