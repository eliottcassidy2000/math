"""Exact bounded controls for prime power shells and their carry decoder.

No conditional prime test, floating point, or external dependency.
Prime tests use a sieve of every prime <=100000 and are complete <=10^10.
"""

from collections import Counter
from fractions import Fraction
from math import comb, isqrt


LIMIT = 10 ** 10
flags = bytearray(b"\x01") * 100001
flags[0:2] = b"\x00\x00"
for i in range(2, isqrt(100000) + 1):
    if flags[i]:
        flags[i * i:100001:i] = b"\x00" * (((100000 - i * i) // i) + 1)
PRIMES = [i for i in range(2, 100001) if flags[i]]


def require(test, message):
    if not test:
        raise RuntimeError(message)


def factor_or_prime(n):
    require(2 <= n <= LIMIT, "outside declared primality universe")
    for p in PRIMES:
        if p * p > n:
            return None
        if n % p == 0:
            return p
    require(isqrt(n) <= 100000, "incomplete trial divisor list")
    return None


def next_prime(n):
    q = n + 1
    while factor_or_prime(q) is not None:
        q += 1
    return q


def root_floor(n, exponent):
    lo, hi = 0, 1
    while hi ** exponent <= n:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** exponent <= n:
            lo = mid
        else:
            hi = mid
    return lo


def floor(x):
    return x.numerator // x.denominator


def valuation2(n):
    require(n > 0, "valuation of nonpositive number")
    return (n & -n).bit_length() - 1


def main():
    # Independent divisor-by-divisor primality audit of the sieve.
    for n in range(2, 100001):
        brute = all(n % d for d in range(2, isqrt(n) + 1))
        require(bool(flags[n]) == brute, "sieve independent audit")

    controls = 0
    valuations = Counter()
    parents = [p for p in PRIMES if p <= 97]
    for p in parents:
        for c in range(2, 7):
            if p ** c >= LIMIT - 1000:
                continue
            q = next_prime(p ** c)
            require(p ** c < q < (p + 1) ** c - 1, "nested prime shell")
            require(root_floor(q, c) == p, "exact shell parent decoder")
            d = q - p ** c
            require(d % 2 == int(p == 2), "gap parity seam")
            valuations[c, valuation2(d)] += 1
            controls += 1

    print("FINITE-EXACT; primality universe integers2..10^10")
    print("sieve_verified_by_full_trial_division_integers2..100000")
    print("prime_parents<=97_exponents2..6_under_bound", controls)
    print("gap_valuation_counts", sorted(valuations.items()))
    for c in range(2, 7):
        chain = [2]
        gaps = []
        while chain[-1] ** c < LIMIT - 1000:
            q = next_prime(chain[-1] ** c)
            gaps.append(q - chain[-1] ** c)
            chain.append(q)
        print("finite_greedy_chain", c, chain, "gaps", gaps)

    cubic = [2, 11, 1361, 2521008887]
    for p, q in zip(cubic, cubic[1:]):
        require(next_prime(p ** 3) == q, "displayed cubic prefix")
        require(all(q % d for d in range(2, isqrt(q) + 1)),
                "displayed prime independently trial divided")
    scale = 10 ** 15
    lo = root_floor(cubic[-1] * scale ** 81, 81)
    hi = root_floor((cubic[-1] + 1) * scale ** 81, 81) + 1
    print("any_constant_realizing_cubic_prefix_outer_bracket", lo, hi,
          "denominator", scale)
    children = [q for q in range(9, 27) if factor_or_prime(q) is None]
    require(children == [11, 13, 17, 19, 23], "cube parent2 children")
    print("all_cubic_children_of_2", children)

    # Exact full-binomial carry and halving identities, independently by powers.
    floor_controls = 0
    for p in range(2, 20):
        for c in range(2, 6):
            for numerator in range(64):
                theta = Fraction(numerator, 64)
                x = p + theta
                expansion = sum(comb(c, j) * p ** (c - j) * theta ** j
                                for j in range(c + 1))
                require(expansion == x ** c, "binomial power")
                digit = floor(x ** c) - p ** c
                tail = x ** c - (p ** c + digit)
                require(0 <= tail < 1, "exact residual tail")
                require(digit == floor(expansion - p ** c), "carry digit")
                require(floor(x / 2) == p // 2, "integer halving label")
                require(x / 2 - p // 2 == (p % 2 + theta) / 2,
                        "halving parity remainder")
                floor_controls += 1
    print("full_binomial_and_halving_fraction_controls", floor_controls)
    p, theta = 2, Fraction(9, 10)
    full = (p + theta) ** 3
    truncated = p ** 3 + 3 * p ** 2 * theta + 3 * p * theta ** 2
    require(floor(full) == 24 and floor(truncated) == 23,
            "discarded cubic tail changes floor/primality")
    print("discard_theta_cubed_hostile", truncated, full,
          "floors", floor(truncated), floor(full))

    horizontal = 0
    for p in range(3, 1000, 2):
        for c in range(2, 13):
            actual = valuation2((p + 2) ** c - p ** c)
            expected = 1 if c % 2 else valuation2(c) + valuation2(p + 1) + 1
            require(actual == expected, "horizontal odd-chain valuation")
            horizontal += 1
    print("odd_chain_power_valuation_controls", horizontal)

    algebraic = 0
    for p in parents:
        for t in range(1, 7):
            a = 2 ** t
            require(p ** 3 + 2 ** (3 * t) == (p + a) * (p * p - p * a + a * a),
                    "cubic pure-power gap factorization")
            require(p + a > 1 and p * p - p * a + a * a > 1,
                    "proper cubic factors")
            algebraic += 1
        require(p ** 4 + 4 == (p * p - 2 * p + 2) * (p * p + 2 * p + 2),
                "quartic Sophie Germain gap")
    require(next_prime(7 ** 3) == 347 and 347 - 7 ** 3 == 4,
            "multiples-of-four gap must survive")
    require(next_prime(43 ** 3) == 79531 and 79531 - 43 ** 3 == 24,
            "valuation-three is not the pure-cube-gap obstruction")
    require(next_prime(67 ** 3) == 300779 and 300779 - 67 ** 3 == 16,
            "higher power-of-two gap must survive")
    require(next_prime(2 ** 5) == 37 and 37 - 2 ** 5 == 5,
            "cubic triple is exponent-specific")
    print("cubic_factor_controls", algebraic,
          "quartic_factor_controls", len(parents))
    print("multiple_of_four_positive_control 347=7^3+4;first_prime_above343")
    print("higher_gap_controls 79531=43^3+24;300779=67^3+16;both_greedy")
    print("exponent5_control 37=2^5+5;first_prime_above32")
    print("halving_controls 11->5prime;1361->680composite")
    print("scope:finite_prefixes;no_unconditional_least_Mills_decimal_identification")


if __name__ == "__main__":
    main()
