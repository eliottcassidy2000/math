#!/usr/bin/env python3
"""Independent exact synthesis check: marked vertices, colored parents, root-switch workload.

No class generation, repository writes, third-party modules, or implicit assertions.
"""
import itertools
import json
import math
from functools import lru_cache


def require(test, message):
    if not test:
        raise RuntimeError(message)


def partitions(n, upper=None):
    if n == 0:
        yield ()
        return
    top = min(n, upper if upper is not None else n)
    top -= (top + 1) % 2
    for size in range(top, 0, -2):
        for multiplicity in range(1, n // size + 1):
            for tail in partitions(n - size * multiplicity, size - 2):
                yield ((size, multiplicity),) + tail


@lru_cache(None)
def terms(n):
    result = []
    for parts in partitions(n):
        e = sum(m * (s - 1) // 2 + s * m * (m - 1) // 2 for s, m in parts)
        e += sum(m * k * math.gcd(s, t) for i, (s, m) in enumerate(parts)
                 for t, k in parts[i + 1:])
        z = math.prod(s ** m * math.factorial(m) for s, m in parts)
        quotient, rem = divmod(math.factorial(n), z)
        require(rem == 0, "Centralizer does not divide factorial")
        result.append((quotient << e, sum(m for _, m in parts),
                       sum(m for s, m in parts if s == 1)))
    return tuple(result)


def colored(n, colors):
    numerator = sum(weight * colors ** cycles for weight, cycles, _ in terms(n))
    quotient, rem = divmod(numerator, math.factorial(n))
    require(rem == 0, "Colored Burnside total is not integral")
    return quotient


def marked(n, k):
    numerator = sum(weight * math.prod(range(fixed - k + 1, fixed + 1))
                    for weight, _, fixed in terms(n) if fixed >= k)
    quotient, rem = divmod(numerator, math.factorial(n))
    require(rem == 0, "Marked Burnside total is not integral")
    return quotient


def literal_marked(n, k):
    """Burnside over S_(n-k), literally following ordered endpoints of each pair."""
    total = 0
    for tail in itertools.permutations(range(k, n)):
        perm = tuple(range(k)) + tail
        unseen = set(itertools.combinations(range(n), 2))
        orbits = 0
        consistent = True
        while unseen and consistent:
            i, j = next(iter(unseen))
            a, b = i, j
            while True:
                unseen.discard(tuple(sorted((a, b))))
                a, b = perm[a], perm[b]
                if (a, b) == (j, i):
                    consistent = False
                    break
                if (a, b) == (i, j):
                    orbits += 1
                    break
        if consistent:
            total += 1 << orbits
    quotient, rem = divmod(total, math.factorial(n - k))
    require(rem == 0, "Literal group sum is not integral")
    return quotient


def main():
    prefix = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056,
              903753248, 154108311168, 48542114686912, 28401423719122304,
              31021002160355166848]
    require([colored(n, 1) for n in range(len(prefix))] == prefix, "A000568 prefix")
    comparisons = 0
    for n in range(41):
        for k in range(min(n, 6) + 1):
            direct = marked(n, k)
            parent = (1 << (k * (k - 1) // 2)) * colored(n - k, 1 << k)
            require(direct == parent, f"Pointing/color mismatch at {(n, k)}")
            comparisons += 1
    literal_checks = 0
    for n in range(8):
        for k in range(n + 1):
            require(marked(n, k) == literal_marked(n, k), f"Literal mismatch {(n, k)}")
            literal_checks += 1
    rows = []
    for n in range(1, 21):
        a, p = colored(n, 1), marked(n, 1)
        homology = p - 2 * a if n % 2 == 0 else None
        if homology is not None:
            require(homology >= 0, "Even root-switch homology negative")
        rows.append({"n": n, "ordinary": a, "pointed": p,
                     "even_root_switch_homology_dimension": homology})
    require(rows[8]["pointed"] == 1716608, "Reported order-nine pointed bank size")
    require(marked(3, 1) == 4 and colored(3, 1) == 2, "Naive unpointing hostile")
    # A tiny relative error is not an exact-integer error bound.
    n = 30
    one_three_cycle = (1 << ((n - 3) * (n - 4) // 2 + n - 2))
    denominator = 3 * math.factorial(n - 3)
    require(one_three_cycle > denominator, "Asymptotic rounding hostile")
    print(json.dumps({"status": "FINITE-EXACT", "direct_vs_colored_checks": comparisons,
                      "literal_group_checks": literal_checks,
                      "bounds": {"n": 40, "ordered_marks": 6, "literal_n": 7},
                      "root_switch_formula": "Inherited Q theorem; dimensions evaluated, Q not built",
                      "n30_single_3cycle_term_greater_than_one": True,
                      "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
