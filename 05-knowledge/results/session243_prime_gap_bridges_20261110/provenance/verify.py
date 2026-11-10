#!/usr/bin/env python3
"""Exact finite controls for the Session243 prior-arithmetic synthesis."""
import argparse
import json
from math import factorial, gcd, isqrt, prod


def require(ok, message):
    if not ok:
        raise ValueError(message)


def prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def factor(n):
    result = []
    d = 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            result.append((d, e))
        d += 1
    if n > 1:
        result.append((n, 1))
    return result


def first_offset(P, sign, start):
    d = start
    while not prime(P + sign * d):
        d += 1
        require(sign == 1 or P - d >= 2, "negative offset exhausted")
    return d


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fail", choices=["skip-one", "wheel-prime", "brown-packet"])
    args = parser.parse_args()
    primes = [n for n in range(2, 200) if prime(n)]
    rows = []
    for n in range(2, 9):
        P = prod(primes[:n])
        F, L = first_offset(P, 1, 2), first_offset(P, -1, 2)
        f, l = first_offset(P, 1, 1), first_offset(P, -1, 1)
        plus1, minus1 = prime(P + 1), prime(P - 1)
        require(f == (1 if plus1 else F), "plus offset restoration")
        require(l == (1 if minus1 else L), "minus offset restoration")
        interior = [x for x in range(P - l + 1, P + f) if prime(x)]
        require(not interior, "crossing endpoints not consecutive")
        equal = None
        if F == L:
            between = [x for x in range(P - L + 1, P + F) if prime(x)]
            expected = [x for x in [P - 1, P + 1] if prime(x)]
            require(between == expected, "equal offset missing sidecar")
            require((not between) == (not plus1 and not minus1), "iff failed")
            equal = {"raw_gap": F + L, "intervening_primes": between}
        rows.append({"n": n, "P": P, "F": F, "L": L,
                     "plus1_prime": plus1, "minus1_prime": minus1,
                     "nearest_pair": [P - l, P + f], "crossing_gap": f + l,
                     "equal_offset_control": equal})

    r6 = next(r for r in rows if r["n"] == 6)
    require(r6["F"] == r6["L"] == 17, "30030 old offsets")
    require(r6["equal_offset_control"]["intervening_primes"] == [30029], "30030 hostile")
    require(r6["crossing_gap"] == 18, "30030 real gap")
    require(30031 == 59 * 509 and not prime(30031), "Euclid composite control")
    r8 = next(r for r in rows if r["n"] == 8)
    require(r8["F"] == r8["L"] == 23, "positive equal-offset control")
    require(r8["crossing_gap"] == 46, "positive consecutive control")

    obstructions = []
    for d in range(2, 101):
        r = factor(d)[0][0]
        P = 1
        for n, p in enumerate(primes, 1):
            P *= p
            if P % r == 0 and P - d > r:
                break
        require(P % r == d % r == 0, "divisor mismatch")
        require(P - d > r and (P - d) % r == (P + d) % r == 0, "proper divisor")
        require(not prime(P - d) and not prime(P + d), "fixed offset obstruction")
        obstructions.append({"d": d, "r": r, "first_tested_n": n, "P": P})
    require(prime(6 - 3) and gcd(6, 3) == 3, "endpoint-equals-divisor hostile")
    require(prime(30 - 1) and prime(30 + 1), "twin boundary control")

    brown = []
    subset_count = 0
    for n in range(4, 41):
        N = factorial(n)
        M = N // 4
        # Legendre valuation formula provides complete prime-power packets;
        # no trial factorization of a giant factorial is needed.
        packets = []
        for p in primes:
            if p > n:
                break
            power, e = p, 0
            while power <= n:
                e += n // power
                power *= p
            if p == 2:
                e -= 2
            if e:
                packets.append(p ** e)
        require(prod(packets) == M, "factorial packet product")
        divisors = [1]
        for packet in packets:
            divisors += [d * packet for d in divisors]
        subset_count += len(divisors)
        solutions = []
        for x in divisors:
            require(M % x == 0, "not a divisor")
            y = M // x
            require(gcd(x, y) == 1, "split prime packet")
            if x < y and y - x == 1:
                m = x + y
                require(m * m == N + 1, "Brown reconstruction")
                require(gcd(m, N) == 1, "roughness property")
                solutions.append([n, m, x, y])
        s = isqrt(N + 1)
        expected = [s] if s * s == N + 1 else []
        require([row[1] for row in solutions] == expected, "independent square check")
        brown.extend(solutions)
    require(brown == [[4, 5, 2, 3], [5, 11, 5, 6], [7, 71, 35, 36]], "finite Brown list")

    decompositions = []
    for p, q in zip(primes[1:30], primes[2:31]):
        d = (q - p) // 2
        union = set()
        interval_counts = []
        for a in range(p, q, 2):
            part = {x for x in range(a * a + 1, (a + 2) ** 2) if prime(x)}
            require(not union.intersection(part), "overlapping prime sets")
            union.update(part)
            interval_counts.append(len(part))
        direct = {x for x in range(p * p + 1, q * q) if prime(x)}
        require(union == direct and len(interval_counts) == d, "Brocard bracket union")
        require(q * q - p * p == 4 * d * (p + d), "interval length")
        decompositions.append({"p": p, "next_prime": q, "gap": 2 * d,
                               "brackets": d, "prime_count": len(direct)})

    if args.fail == "skip-one":
        require(r6["crossing_gap"] == r6["F"] + r6["L"], "intentional lost offset1")
    if args.fail == "wheel-prime":
        require(prime(30031), "intentional wheel-survivor primality")
    if args.fail == "brown-packet":
        require(gcd(6, 6) == 1, "intentional split prime-power packet")

    print(json.dumps({"status": "PASS", "primorial_rows": rows,
                      "fixed_offsets_checked": len(obstructions),
                      "fixed_offset_examples": obstructions[:10],
                      "Brown_n_universe": [4, 40], "Brown_subsets_checked": subset_count,
                      "Brown_solutions": brown,
                      "prime_square_unions_checked": len(decompositions),
                      "prime_square_examples": decompositions[:8]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
