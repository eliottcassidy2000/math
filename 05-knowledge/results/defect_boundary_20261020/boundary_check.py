"""Cycle boundary transfer, signed ordered-list transforms, and live profiles.

Run from repository root:
  python 05-knowledge/results/defect_boundary_20261020/boundary_check.py

No tournament-core enumeration occurs in the main count through n=100.
The small independent control imports the previously audited exact core bank.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import comb, factorial, gcd
from pathlib import Path


def need(test, message):
    if not test:
        raise RuntimeError(message)


def phi(n):
    return sum(gcd(k, n) == 1 for k in range(1, n + 1))


def odd_partitions(maximum, minimum=1):
    """Every odd partition of mass <= maximum exactly once, including empty."""
    def visit(parts, lower, room):
        yield parts
        for d in range(lower, room + 1, 2):
            yield from visit(parts + (d,), d, room - d)
    yield from visit((), minimum, maximum)


def edge_exponent(parts):
    return sum((d - 1) // 2 for d in parts) + sum(
        gcd(parts[i], parts[j]) for i in range(len(parts)) for j in range(i))


def z_value(parts):
    return product_int(d ** m * factorial(m) for d, m in Counter(parts).items())


def product_int(values):
    out = 1
    for value in values:
        out *= value
    return out


def divisibility_profile(parts):
    return {d: sum(p % d == 0 for p in parts)
            for d in range(1, max(parts, default=1) + 1, 2)}


def future_profile_transfer(maximum):
    """Eliminate cycle types descending and forget expired divisibility data."""
    largest = maximum if maximum % 2 else maximum - 1
    active = tuple(range(1, largest + 1, 2))
    states = {(0, (0,) * len(active)): Fraction(1)}
    ways = [1] + [0] * maximum
    history, transitions = [], 0
    phis = {d: phi(d) for d in active}
    for d in range(largest, 2, -2):
        need(active[-1] == d, "active divisor ordering")
        index = {v: j for j, v in enumerate(active)}
        divisors = [v for v in active if d % v == 0]
        smaller_divisors = [index[v] for v in divisors if v != d]
        bases = [Fraction(1 << (k * ((d - 1) // 2) + d * k * (k - 1) // 2),
                          d ** k * factorial(k))
                 for k in range(maximum // d + 1)]
        nxt = defaultdict(Fraction)
        for (mass, profile), coefficient in states.items():
            interaction = sum(phis[v] * profile[index[v]] for v in divisors)
            for k in range((maximum - mass) // d + 1):
                updated = list(profile[:-1])
                for j in smaller_divisors:
                    updated[j] += k
                key = (mass + d * k, tuple(updated))
                nxt[key] += coefficient * bases[k] * (1 << (k * interaction))
                transitions += 1
        for mass in range(d, maximum + 1):
            ways[mass] += ways[mass - d]
        active = active[:-1]
        states = dict(nxt)
        raw = sum(ways)
        need(len(states) <= raw, "aggregation exceeds unmerged multiplicity states")
        history.append((d, len(states), raw))
    need(active == (1,), "final profile is cycle count")
    return {(mass, profile[0]): value for (mass, profile), value in states.items()}, history, transitions


def reconstruct_counts(kernel, maximum):
    result = []
    for n in range(maximum + 1):
        total = Fraction(0)
        for (mass, cycles), value in kernel.items():
            if mass <= n:
                r = n - mass
                total += value * Fraction(1 << (r * (r - 1) // 2 + r * cycles),
                                          factorial(r))
        need(total.denominator == 1, "A(n) integral")
        result.append(total.numerator)
    return result


def lah(m, k):
    if m == k == 0:
        return 1
    if not 1 <= k <= m:
        return 0
    return factorial(m) // factorial(k) * comb(m - 1, k - 1)


def multiplicities(parts, maximum):
    for m in product(*(range(maximum // d + 1) for d in parts)):
        if sum(d * k for d, k in zip(parts, m)) <= maximum:
            yield m


def contraction_indices(m):
    return product(*((range(1, k + 1) if k else (0,)) for k in m))


def fixed_tournaments(ds, m):
    parts = tuple(d for d, k in zip(ds, m) for _ in range(k))
    return 1 << edge_exponent(parts)


def fixed_pairfree(ds, m):
    """Tensor signed Lah transform, normalized as fixed-permutation counts."""
    total = 0
    for k in contraction_indices(m):
        factor = product_int((-d) ** (a - b) * lah(a, b)
                             for d, a, b in zip(ds, m, k))
        total += factor * fixed_tournaments(ds, k)
    return total


def differential_step(ds, target, i, known):
    """Compute c_target from the i-th Möbius-scaled differential equation."""
    m = list(target)
    need(m[i] > 0, "differentiation coordinate is nonzero")
    m[i] -= 1
    m = tuple(m)
    a_i = 1 << ((ds[i] - 1) // 2)
    b = tuple(1 << gcd(ds[i], d) for d in ds)
    total = 0
    for k in contraction_indices(m):
        factor = product_int(lah(a, c) * scale ** c
                             * (d * (scale - 1)) ** (a - c)
                             for d, a, c, scale in zip(ds, m, k, b))
        total += factor * known[k]
    total *= a_i
    total -= 2 * ds[i] * m[i] * known[m]
    if m[i] >= 2:
        prev = list(m)
        prev[i] -= 1
        total -= ds[i] ** 2 * m[i] * (m[i] - 1) * known[tuple(prev)]
    return total


def mobius_map(x, length, scale):
    return scale * x / (1 - length * (scale - 1) * x)


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    print("EXACT: dynamic future-divisor profiles; no individual tournament cores.")
    kernel, history, transitions = future_profile_transfer(100)
    counts = reconstruct_counts(kernel, 100)
    print("selected transfer rows (processed length, retained states, unmerged states):")
    for row in history:
        if row[0] in (99, 51, 31, 21, 15, 11, 9, 7, 5, 3):
            print(row)
    print("peak retained states=", max(row[1] for row in history),
          "peak unmerged states=", max(row[2] for row in history),
          "actual transfer transitions=", transitions)
    print("final kernel states=", len(kernel))
    print("A(0..20)=", counts[:21])
    print("A100 digits=", len(str(counts[100])),
          "decimal_sha256=", sha256(str(counts[100]).encode()).hexdigest())

    # Independent Davis sum using explicit partitions and pairwise gcds.
    brute = [Fraction(0) for _ in range(41)]
    brute_kernel = defaultdict(Fraction)
    partition_controls = 0
    for parts in odd_partitions(40):
        exponent = edge_exponent(parts)
        profile = divisibility_profile(parts)
        dividend = sum(phi(d) * M * M for d, M in profile.items()) - len(parts)
        need(dividend == 2 * exponent, "full gcd/zeta quadratic identity")
        value = Fraction(1 << exponent, z_value(parts))
        brute[sum(parts)] += value
        partition_controls += 1
        if not parts or min(parts) >= 3:
            brute_kernel[(sum(parts), len(parts))] += value
    need(all(brute[n] == counts[n] for n in range(41)), "independent Davis comparison")
    need(all(kernel[key] == value for key, value in brute_kernel.items()),
         "dynamic kernel vs explicit partitions")
    reference = load_module("boundary_prior",
                            Path(__file__).resolve().parents[3]
                            / "04-computation/a000568_core_tail_ladder_codex.py")
    need(all(counts[n] == value for n, value in reference.KNOWN.items()), "known A1..20")
    need(counts[100] == reference.A100, "frozen independently recorded A100")
    print("Independent explicit partition controls through mass40=", partition_controls)
    print("Known A1..20 and frozen A100 matched.")

    # Same mass, length, M3 and current exponent do not determine the M5 interaction.
    left, right = (7, 15), (9, 13)
    summary = lambda p: (sum(p), len(p), sum(d % 3 == 0 for d in p), edge_exponent(p))
    need(summary(left) == summary(right) == (22, 2, 1, 11), "profile hostile common data")
    increments = tuple(edge_exponent(p + (5,)) - edge_exponent(p) for p in (left, right))
    need(increments == (8, 4), "profile hostile missing M5")
    print("Missing-coordinate hostile:", left, right, "common summary", summary(left),
          "adding5 exponent increments", increments)

    # Ordered-list transforms are inverses, including the phase weight d^(m-k).
    matrix_controls = 0
    for d in (1, 3, 5, 9):
        for m in range(13):
            for j in range(m + 1):
                value = sum((-d) ** (m - k) * lah(m, k)
                            * d ** (k - j) * lah(k, j) for k in range(j, m + 1))
                need(value == int(m == j), "signed Lah inverse")
                matrix_controls += 1
    print("Signed ordered-list transform inverse controls=", matrix_controls)

    # Closed operator recurrence on finite odd cycle alphabets. Every available
    # differentiation direction must agree, not merely the chosen recursion.
    operator_controls, direction_controls = 0, 0
    for ds, bound in (((1,), 12), ((3,), 36), ((5,), 60), ((9,), 108),
                      ((3, 5), 45), ((3, 5, 9), 45)):
        states = sorted(multiplicities(ds, bound), key=lambda m: (sum(m), m))
        known = {(0,) * len(ds): 1}
        for m in states[1:]:
            values = [differential_step(ds, m, i, known)
                      for i in range(len(ds)) if m[i]]
            need(len(set(values)) == 1, "mixed differential directions disagree")
            need(values[0] == fixed_pairfree(ds, m), "operator vs signed transform")
            need(values[0] >= 0, "negative fixed-permutation core count")
            known[m] = values[0]
            operator_controls += 1
            direction_controls += len(values)
        if len(ds) == 1:
            print("pure cycle length", ds[0], "pair-free fixed counts m=0..6:",
                  [known[(m,)] for m in range(7)])
    print("Finite-alphabet operator states=", operator_controls,
          "all-direction checks=", direction_controls)
    for d in (1, 3, 5, 9):
        for a, b in ((2, 8), (4, 32), (8, 512)):
            x = Fraction(1, 1000003)
            need(mobius_map(mobius_map(x, d, a), d, b) == mobius_map(x, d, a * b),
                 "Mobius scaling semigroup")
    print("Mobius scaling semigroup controls passed.")

    # Independent exact small core-cycle-index bank, used only as an audit.
    prior = load_module("small_core_audit",
                        Path(__file__).resolve().parents[1]
                        / "weighted_core_transport_20261020/core_count.py")
    small_controls = 0
    for n in range(1, 8):
        atlas = prior.e.census(n)
        core_coeff = defaultdict(Fraction)
        for code, rows in atlas["reps"].items():
            if len(prior.g.clone_chains(rows)) != n:
                continue
            group = prior.automorphisms(rows)
            for p in group:
                core_coeff[prior.cycle_type(p)] += Fraction(1, len(group))
        for parts in odd_partitions(n):
            if sum(parts) != n:
                continue
            ds = tuple(sorted(set(parts)))
            m = tuple(parts.count(d) for d in ds)
            predicted = Fraction(fixed_pairfree(ds, m), z_value(parts))
            need(predicted == core_coeff.get(parts, Fraction(0)),
                 "signed transform vs independent actual core automorphisms")
            small_controls += 1
    print("Actual small core cycle-type controls through7=", small_controls)
    print("ALL CHECKS PASSED (explicit exceptions active under normal and -O Python).")


if __name__ == "__main__":
    main()
