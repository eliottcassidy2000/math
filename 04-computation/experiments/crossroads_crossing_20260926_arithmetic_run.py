"""Exact controls for maximal-rise resets and synchronized binary siblings."""
from collections import Counter
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def v2(n):
    require(n > 0, "positive valuation argument")
    return (n & -n).bit_length() - 1


def oddpart(n):
    return n >> v2(n)


def step(n):
    require(n > 0 and n % 2 == 1, "positive odd source")
    return oddpart(3 * n + 1)


def iterate(n, k):
    for _ in range(k):
        n = step(n)
    return n


def compressed(n):
    r = v2(n + 1)
    u = (n + 1) >> r
    return oddpart(3 ** r * u - 1)


def main():
    totals = Counter()
    for n in range(1, 200_000, 2):
        r = v2(n + 1)
        u = (n + 1) >> r
        numerator = 3 ** r * u - 1
        s = v2(numerator)
        z = numerator >> s
        require(iterate(n, r) == z == compressed(n), "first reset")
        other = iterate(2 * n + 1, r)
        require(other == (1 << (s + 1)) * z + 1, "paired exact height")
        require((step(z) == step(other)) == (s == 1), "iff synchronized merge")
        if s >= 2:
            require(step(other) == 3 * (1 << (s - 1)) * z + 1 > step(z),
                    "hostile strict separation")
        if s == 1:
            totals["cheap_resets"] += 1
            # The larger partner reduces to n. The criterion at n itself
            # fails whenever n has a smaller odd binary parent.
            if r >= 2:
                require((3 ** (r - 1) * u - 1) % 4 != 2,
                        "no consecutive tower reductions")
        if (3 ** r * u) % 8 == 3:
            require(compressed(compressed(n)) == compressed(2 * n + 1),
                    "two compressed blocks versus one")
            totals["compressed_two_to_one"] += 1
        totals["sources"] += 1

    mersenne = []
    for r in range(1, 257):
        n = (1 << r) - 1
        s = v2(3 ** r - 1)
        require(s == (1 if r % 2 else 2 + v2(r)), "elementary LTE")
        target = compressed(n)
        require((target < n) == (r in (2, 4, 8)), "complete Mersenne first-reset classification")
        require((target == n) == (r == 1), "only equality in this family")
        if r % 2:
            require(compressed(target) == compressed(2 * n + 1), "Mersenne paired reset")
        if r <= 16:
            mersenne.append([r, n, s, target, (target > n) - (target < n)])

    # Directly recompute the root lane's unwrapped carry table, independently
    # from its iterative carry recurrence, then check the projector sign.
    carries = []
    for word in ("110111", "111011", "111101", "111110"):
        carry = sum((1 << j) * 3 ** word[j + 1:].count("1")
                    for j, bit in enumerate(word) if bit == "1")
        residue = (-pow(3 ** 5, -1, 64) * carry) % 64
        require((3 ** 5 * residue + carry) % 64 == 0, "projector sign")
        carries.append([word, carry, residue])
    require([row[1] for row in carries] == [287, 251, 227, 211], "unwrapped carry values")

    print(json.dumps({
        "scope": "FINITE-EXACT controls; general claims proved in companion note",
        "source_universe": "all 100000 positive odd n below 200000; larger partners below400000",
        "counts": dict(totals),
        "Mersenne_exponents_checked": [1, 256],
        "Mersenne_columns": ["r", "2^r-1", "s", "C(n)", "sign(C(n)-n)"],
        "Mersenne_rows": mersenne,
        "hostile": {"n": 31, "larger_partner": 63, "common_value": iterate(31, 6),
                    "common_value_exceeds_both": iterate(31, 6) > 63},
        "root_unwrapped_carry_audit": carries,
        "theorem_relative_density_among_odds": {"cheap_reset_sources": "1/2", "reducible_larger_partners": "1/4"}
    }, indent=2))


if __name__ == "__main__":
    main()
