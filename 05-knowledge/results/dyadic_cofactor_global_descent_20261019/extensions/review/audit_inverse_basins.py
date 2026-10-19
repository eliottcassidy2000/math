#!/usr/bin/env python3
"""Independent exact inverse-set controls; supplementary to the proof."""

from collections import defaultdict
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import random


def check(ok, details):
    if not ok:
        raise RuntimeError(details)


def F(n):
    c = n & -n
    return n - n // c + 2 * c - 1


def inverse_compiler(m):
    out = {}
    for a in range(1, m.bit_length() + 1):
        c = 1 << a
        quotient, remainder = divmod(m - 1, c - 1)
        if remainder == 0 and quotient >= 3:
            n = c * (quotient - 2)
            check((n & -n) == c, ("exact_core", m, a, n))
            check(F(n) == m, ("forward_check", m, a, n))
            out[a] = n
    check(len(set(out.values())) == len(out), ("distinctness", m))
    return out


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def proper_divisor_sum(n):
    return sum(d for d in range(1, n) if n % d == 0)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def main():
    source = Path(__file__).resolve().parents[1] / "INVERSE_BASIN_DIVISOR_LATTICES.md"
    result = {"source_sha256": sha256(source.read_bytes()).hexdigest()}
    limit = 100_000
    brute = defaultdict(list)
    # F(n)>=n/2+3, so n<=2m-6 for every positive even preimage of m.
    for n in range(2, 2 * limit - 5, 2):
        m = F(n)
        check(2 * m >= n + 6, ("inverse_search_bound", n))
        if m <= limit:
            brute[m].append(n)
    for m in range(2, limit + 1, 2):
        check(sorted(inverse_compiler(m).values()) == sorted(brute[m]),
              ("complete_inverse_set", m))
    result["complete_general_inverse_even_targets"] = limit // 2
    result["largest_general_target"] = limit
    fixed_counts = {}
    for A in range(1, 513):
        m = (1 << A) * ((1 << (A + 1)) - 1)
        inv = inverse_compiler(m)
        check(list(inv) == divisors(A), ("fixed_divisors", A, inv))
        check(inv[A] == m, ("fixed_self", A))
        if A >= 2:
            p = m - 3
            left, right = (1 << A) + 1, (1 << (A + 1)) - 3
            check(left > 1 and right > 1 and left * right == p,
                  ("core2_composite", A))
        fixed_counts[A] = len(inv)
    result["fixed_exponents_checked"] = len(fixed_counts)
    result["largest_fixed_target_bits"] = m.bit_length()
    for b in range(1, 1025):
        inv = inverse_compiler(1 << b)
        check(list(inv) == divisors(b)[:-1], ("power_divisors", b, inv))
    result["power_target_exponents_checked"] = 1024
    residue_count = 0
    rng = random.Random(20261019)
    pairs = [(A, a) for A in range(1, 257) for a in range(1, 513)]
    pairs += [(rng.randrange(1, 4097), rng.randrange(1, 4097)) for _ in range(10000)]
    for A, a in pairs:
        M = (1 << a) - 1
        target_minus_one = (1 << (2 * A + 1)) - (1 << A) - 1
        r = A % a
        check((target_minus_one % M == 0) == (r == 0), ("necessity", A, a))
        if r:
            if 2 * r + 1 < a:
                residue = (1 << (2 * r + 1)) - (1 << r) - 1
                check(0 < residue < M, ("positive_residue", A, a))
            else:
                residue = (1 << (2 * r + 1 - a)) - (1 << r) - 1
                check(-M < residue < 0, ("negative_residue", A, a))
            check(target_minus_one % M == residue % M, ("residue_equality", A, a))
        residue_count += 1
    result["modular_residue_cases"] = residue_count
    inv496 = inverse_compiler(496)
    check(inv496 == {1: 986, 2: 652, 4: 496}, ("example496", inv496))
    check([a for a, n in inv496.items() if prime(n // (n & -n))] == [2, 4],
          "native_example496")
    check(proper_divisor_sum(652) == proper_divisor_sum(496) == 496,
          "native_divisor_sum496")
    check(proper_divisor_sum(986) == 634 and F(986) == 496, "nonnative_guard986")
    for t in range(1, 513):
        d = 1 << t
        first = d * (d - 1)
        check(F(first) == d * d, ("first_resonance", t))
        initial = d * d * (d - 1)
        second = d * (d * d + d - 1)
        check(F(initial) == second and F(second) == d ** 3,
              ("second_resonance", t))
    result["resonance_link_exponents_checked"] = 512
    result["status"] = "PASS"
    Path(__file__).with_name("inverse_basins_audit_results.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
