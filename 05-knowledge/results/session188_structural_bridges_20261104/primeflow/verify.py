#!/usr/bin/env python3
"""Deterministic standard-library checks of exact weighted sigma closure.

No assertion is used: checks are unchanged under python -O.  The direct oracle
uses geometric sums, not the LTE formulas used by the closure implementation.
"""
import itertools
import json
from math import gcd, isqrt, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def valuation(n, r):
    require(n > 0, "valuation is defined here only for positive integers")
    result = 0
    while n % r == 0:
        n //= r
        result += 1
    return result


def order(p, r):
    x = 1
    for d in range(1, r):
        x = x * p % r
        if x == 1:
            return d
    raise RuntimeError("order does not exist")


def geometric_mod(p, n, modulus):
    """Independent binary geometric sum, requiring no division or inverse."""
    power, total = 1, 0
    block_power, block_sum = p % modulus, 1 % modulus
    while n:
        if n & 1:
            total = (total + power * block_sum) % modulus
            power = power * block_power % modulus
        block_sum = block_sum * (1 + block_power) % modulus
        block_power = block_power * block_power % modulus
        n >>= 1
    return total


def direct_profile(p, n, targets):
    value = (pow(p, n) - 1) // (p - 1)
    return tuple(valuation(value, r) for r in targets)


class Source:
    def __init__(self, p, targets, euler=False):
        require(prime(p) and p % 2 == 1, "source must be an odd prime")
        require(not euler or p % 4 == 1, "Euler prime must be 1 modulo 4")
        require(len(set(targets)) == len(targets), "targets must be distinct")
        require(all(prime(r) and r % 2 == 1 and r != p for r in targets),
                "targets must be distinct odd primes other than source")
        self.p, self.targets, self.euler = p, tuple(targets), euler
        self.b = 2 if euler else 1
        self.orders = tuple(order(p, r) for r in targets)
        self.constants = tuple(0 if d == 1 else valuation(pow(p, d)-1, r)
                               for r, d in zip(targets, self.orders))
        self.allowed = tuple(d % (4 if euler else 2) != 0
                             if euler else d % 2 == 1
                             for d in self.orders)

    def role(self, n, physical=True):
        if self.euler:
            return n >= 2 and n % 4 == 2
        return n >= (3 if physical else 1) and n % 2 == 1

    def rank(self, index, k):
        require(k >= 0, "negative demand")
        if k == 0:
            return 1
        r, d, c = (self.targets[index], self.orders[index], self.constants[index])
        return r ** k if d == 1 else d * r ** max(0, k - c)

    def profile(self, n):
        result = []
        for r, d, c in zip(self.targets, self.orders, self.constants):
            result.append(valuation(n, r) if d == 1 else
                          c + valuation(n, r) if n % d == 0 else 0)
        return tuple(result)

    def lower_length(self, lower):
        require(len(lower) == len(self.targets), "wrong demand dimension")
        require(all(isinstance(k, int) and k >= 0 for k in lower), "bad demand")
        if any(k and not allowed for k, allowed in zip(lower, self.allowed)):
            return None
        return lcm(self.b, *(self.rank(i, k) for i, k in enumerate(lower)))

    def closure(self, lower):
        length = self.lower_length(lower)
        return None if length is None else self.profile(length)

    def smallest_box_witness(self, lower, upper):
        require(len(lower) == len(upper) == len(self.targets), "box dimension")
        require(all(0 <= a <= b for a, b in zip(lower, upper)), "invalid box")
        length = self.lower_length(lower)
        if length is None:
            return None
        if any(a > b for a, b in zip(self.profile(length), upper)):
            return None
        if length > 1:
            return length
        forbidden = set()
        for r, d, c, cap in zip(self.targets, self.orders, self.constants, upper):
            if d == 1 and cap == 0:
                forbidden.add(r)
            if d > 1 and d % 2 == 1 and prime(d) and c > cap:
                forbidden.add(d)
        candidate = 3
        while not prime(candidate) or candidate in forbidden:
            candidate += 2
        return candidate


def at_most(a, b):
    return all(x <= y for x, y in zip(a, b))


def check():
    counts = {"rank_recurrences": 0, "direct_profiles": 0,
              "closures": 0, "adjunction_comparisons": 0,
              "monotonicity_covers": 0, "boxes": 0,
              "meet_checks": 0, "unbounded_extension_checks": 0}
    small = [n for n in range(3, 32, 2) if prime(n)]

    # Ranks checked against the first zero of an iterative repunit recurrence.
    for p in small[:7]:
        for r in small[:6]:
            if p == r:
                continue
            source = Source(p, (r,))
            for k in range(1, 4):
                predicted = source.rank(0, k)
                modulus = r ** k
                residue = 0
                first = None
                for n in range(1, predicted + 1):
                    residue = (p * residue + 1) % modulus
                    if residue == 0:
                        first = n
                        break
                require(first == predicted, ("rank", p, r, k, first, predicted))
                counts["rank_recurrences"] += 1

    configs = []
    for p in small:
        targets_pool = [r for r in small[:6] if r != p]
        for euler in (False, True):
            if euler and p % 4 != 1:
                continue
            for targets in itertools.combinations(targets_pool, 3):
                configs.append(Source(p, targets, euler))

    bound = 180
    interval_options = [(lo, hi) for lo in range(3) for hi in range(lo, 3)]
    for source in configs:
        # Cache all physical profiles using direct arbitrary-precision sums.
        physical = [(n, direct_profile(source.p, n, source.targets))
                    for n in range(2, bound + 1) if source.role(n)]
        for n, oracle in physical:
            require(source.profile(n) == oracle, ("profile", source.p, n))
            counts["direct_profiles"] += 1
        with_identity = physical if source.euler else [(1, (0, 0, 0))] + physical

        for lower in itertools.product(range(3), repeat=3):
            length = source.lower_length(lower)
            closed = source.closure(lower)
            if closed is None:
                require(all(not at_most(lower, row) for _, row in physical),
                        "disallowed demand unexpectedly realizable")
                continue
            require(source.role(length, physical=False), "bad role for L")
            require(at_most(lower, closed), "not extensive")
            require(source.closure(closed) == closed, "not idempotent")
            require(source.lower_length(closed) == length, "L not idempotent")
            for r, k in zip(source.targets, closed):
                require(geometric_mod(source.p, length, r ** (k + 1)) != 0,
                        "closure valuation is too small")
                if k:
                    require(geometric_mod(source.p, length, r ** k) == 0,
                            "closure valuation is too large")
            counts["closures"] += 1
            for n, row in with_identity:
                require(at_most(lower, row) == (n % length == 0), "adjunction")
                counts["adjunction_comparisons"] += 1
            for i, allowed in enumerate(source.allowed):
                if not allowed:
                    continue
                bigger = list(lower)
                bigger[i] += 1
                require(at_most(closed, source.closure(tuple(bigger))), "monotone")
                counts["monotonicity_covers"] += 1
            # Fresh-prime extension changes length but not the tested profile.
            fresh = max((*source.targets, *source.orders, 2)) + 1
            while not prime(fresh):
                fresh += 1
            extended = length * fresh
            require(source.role(extended), "extension role")
            for r, k in zip(source.targets, closed):
                require(geometric_mod(source.p, extended, r ** (k + 1)) != 0,
                        "extension changes valuation upward")
                if k:
                    require(geometric_mod(source.p, extended, r ** k) == 0,
                            "extension changes valuation downward")
            counts["unbounded_extension_checks"] += 1

        # Every box with 0<=lower<=upper<=2, compared with all n<=180.
        for specs in itertools.product(interval_options, repeat=3):
            lower = tuple(x[0] for x in specs)
            upper = tuple(x[1] for x in specs)
            predicted = source.smallest_box_witness(lower, upper)
            actual = next((n for n, row in physical
                           if at_most(lower, row) and at_most(row, upper)), None)
            truncated = predicted if predicted is not None and predicted <= bound else None
            require(actual == truncated,
                    ("box", source.p, source.targets, source.euler,
                     lower, upper, predicted, actual))
            if predicted is not None:
                require(source.role(predicted), "nonphysical box witness")
            counts["boxes"] += 1

        # Independent strong-divisibility control implies coordinatewise meet.
        for (n, x), (m, y) in zip(physical[:20], reversed(physical[:20])):
            meet = tuple(min(a, b) for a, b in zip(x, y))
            require(direct_profile(source.p, gcd(n, m), source.targets) == meet,
                    "gcd meet")
            require(source.closure(meet) == meet, "fixed meet")
            counts["meet_checks"] += 1

    weighted = Source(3, (11, 23))
    require(weighted.orders == (5, 11) and weighted.constants == (2, 1), "small data")
    require(weighted.closure((2, 1)) == (3, 1), "weighted join")
    require(weighted.smallest_box_witness((2, 1), (2, 1)) is None, "false box")
    require(weighted.smallest_box_witness((2, 0), (2, 1)) == 5, "first atom")
    require(weighted.smallest_box_witness((0, 1), (2, 1)) == 11, "second atom")
    require(Source(3, (11,)).smallest_box_witness((1,), (1,)) is None, "c=2 control")

    diamond = Source(3, (4561, 368089, 71))
    require(diamond.orders == (15, 21, 35), "diamond orders")
    require(diamond.constants == (1, 1, 1), "diamond base valuations")
    atoms = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    for x, n in zip(atoms, (15, 21, 35)):
        require(diamond.lower_length(x) == n, "diamond atom length")
        require(direct_profile(3, n, diamond.targets) == x, "diamond atom profile")
    for x, y in itertools.combinations(atoms, 2):
        require(diamond.closure(tuple(max(a, b) for a, b in zip(x, y))) == (1, 1, 1),
                "diamond pair join")
        require(tuple(min(a, b) for a, b in zip(x, y)) == (0, 0, 0), "diamond meet")
    require(direct_profile(3, 105, diamond.targets) == (1, 1, 1), "diamond top")
    require(direct_profile(3, 3, diamond.targets) == (0, 0, 0), "physical diamond bottom")

    euler = Source(5, (3, 19), True)
    require(euler.orders == (2, 9), "Euler example orders")
    require(euler.closure((0, 1)) == (3, 1), "Euler overflow")
    require(euler.smallest_box_witness((0, 1), (2, 1)) is None, "Euler rejection")
    require(direct_profile(5, 18, (3, 19)) == (3, 1), "Euler direct control")

    ordinary_local, euler_local = Source(3, (13,)), Source(13, (3,), True)
    require(ordinary_local.smallest_box_witness((1,), (1,)) == 3, "local first row")
    require(euler_local.smallest_box_witness((2,), (2,)) == 18, "local second row")
    require(euler_local.profile(2) == (0,), "global incompatibility control")

    # The ordinary zero case must not silently return exponent zero.
    zero = Source(7, (3, 19))
    smallest = zero.smallest_box_witness((0, 0), (0, 0))
    require(smallest == 5, "zero corner smallest prime")
    require(direct_profile(7, 5, (3, 19)) == (0, 0), "zero corner realization")

    return {
        "status": "PASS", "configurations": len(configs),
        "bounded_length_limit": bound, "counts": counts,
        "weighted_join": {"source": 3, "targets": [11, 23],
                          "atoms": [[2, 0], [0, 1]], "join": [3, 1], "length": 55},
        "nondistributive_diamond": {"source": 3, "targets": list(diamond.targets),
                                   "orders": list(diamond.orders),
                                   "constants": list(diamond.constants),
                                   "physical_lengths": [3, 15, 21, 35, 105]},
        "Euler_overflow": {"source": 5, "targets": [3, 19],
                           "lower": [0, 1], "closure": [3, 1], "length": 18},
        "local_global_control": {"support": [3, 13], "Euler_prime": 13,
                                 "local_lengths": [3, 18],
                                 "column_prescribed_lengths": [3, 2],
                                 "incompatible": True},
        "scope": "Exact single-source partial-target realization, not OPN nonexistence.",
    }


if __name__ == "__main__":
    print(json.dumps(check(), indent=2, sort_keys=True))
