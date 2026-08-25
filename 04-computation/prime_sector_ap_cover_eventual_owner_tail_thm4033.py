#!/usr/bin/env python3
"""Optimization-safe exact audit for the THM-4033 prime-sector finite-owner theorem.

This is regression evidence, not a substitute for the symbolic compactness and
track arguments in PRIME-SECTOR-EVENTUAL-OWNER-THEOREM-PACKET.md.
"""

from fractions import Fraction as Q
from hashlib import sha256
from math import floor, gcd
import sys


sys.stdout.reconfigure(newline="\n")


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, floor(n**0.5) + 1))


def phi(n: int) -> int:
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))


def owners(P: int):
    return tuple(
        (p, q)
        for q in range(1, P)
        for p in range(q)
        if gcd(p, q) == 1
    )


def owner_data(P: int, p: int, q: int):
    data = []
    for s in range(q):
        A = Q(P * p * s, q)
        integer = floor(A)
        data.append((s, integer % P, A - integer))
    return tuple(data)


def cyclic(value: int, modulus: int) -> int:
    result = value % modulus
    require(result != 0, "cyclic distance unexpectedly zero")
    return result


def owner_radii(P: int, p: int, q: int, n: int):
    require(n >= q - 1, "all congruence tracks must be present")
    data = owner_data(P, p, q)
    occupied = {b for _, b, _ in data}
    missing = set(range(P)) - occupied
    require(len(occupied) == q and len(missing) == P - q, "prime grid collision")

    def E(s: int) -> int:
        return n - ((n - s) % q)

    require(all(E(s) > 0 for s, _, _ in data), "all tracks need a positive endpoint")
    plus_radius = max(
        min((Q(cyclic(r - b, P)) - f) / E(s) for s, b, f in data)
        for r in missing
    )
    minus_radius = max(
        min((Q(cyclic(b - r, P) - 1) + f) / E(s)
            for s, b, f in data)
        for r in missing
    )
    require(plus_radius > 0, "positive radius must be nonzero")
    return plus_radius, minus_radius


def owner_leading_radii(P: int, p: int, q: int):
    data = owner_data(P, p, q)
    occupied = {b for _, b, _ in data}
    missing = set(range(P)) - occupied
    plus = max(
        min(Q(cyclic(r - b, P)) - f for _, b, f in data)
        for r in missing
    )
    minus = max(
        min(Q(cyclic(b - r, P) - 1) + f for _, b, f in data)
        for r in missing
    )
    return plus, minus


def track_guard_radii(P: int, p: int, q: int):
    data = owner_data(P, p, q)
    plus = Q(1, q)
    minus = Q(1, q)
    for s, _, f in data:
        if s == 0:
            continue
        require(f > 0, "prime owner has an unexpected nonzero-track wall")
        plus = min(plus, (1 - f) / s)
        minus = min(minus, f / s)
    return plus, minus


def split_circular_arc(P: int, left: Q, right: Q):
    if left < 0:
        return ((Q(0), right), (left + P, Q(P)))
    if right > P:
        return ((left, Q(P)), (Q(0), right - P))
    return ((left, right),)


def owner_arcs_at(P: int, n: int):
    arcs = []
    for p, q in owners(P):
        theta = Q(P * p, q)
        plus, minus = owner_radii(P, p, q, n)
        arcs.append(((p, q), split_circular_arc(P, theta - minus, theta + plus)))
    return tuple(arcs)


def owner_formula_deficit(P: int, m: int) -> Q:
    n = m - 1
    return sum(
        (sum(owner_radii(P, p, q, n), Q(0)) for p, q in owners(P)),
        Q(0),
    ) / P


def closed_constant(P: int) -> Q:
    return sum(
        Q(phi(q) * (2 * (P - q) - 1), P * q)
        for q in range(1, P)
    )


def direct_cover_sequence(P: int, max_m: int):
    max_n = max_m - 1
    walls = {Q(0), Q(1)}
    for e in range(1, max_n + 1):
        walls.update(Q(j, P * e) for j in range(P * e + 1))
    ordered = sorted(walls)
    first_cover_mass = [Q(0) for _ in range(max_m + 1)]
    for left, right in zip(ordered, ordered[1:]):
        x = (left + right) / 2
        seen = {0}
        first = None
        for e in range(1, max_n + 1):
            seen.add(floor(P * e * x) % P)
            if len(seen) == P:
                first = e + 1
                break
        if first is not None:
            first_cover_mass[first] += right - left
    result = []
    cumulative = Q(0)
    for m in range(1, max_m + 1):
        cumulative += first_cover_mass[m]
        result.append(cumulative)
    return tuple(result)


def q_equal_P_boundary_trace(P: int, p: int):
    positive = {e * p % P for e in range(P)}
    negative = {0} | {(e * p - 1) % P for e in range(1, P + 1)}
    return positive, negative


def digest(obj) -> str:
    return sha256(repr(obj).encode("utf-8")).hexdigest()


def main() -> None:
    print("PRIME-SECTOR EVENTUAL OWNER AUDIT")
    primes = tuple(P for P in range(2, 32) if is_prime(P))

    constant_rows = []
    for P in primes:
        total = Q(0)
        by_q = []
        for q in range(1, P):
            expected = (Q(P - q, q), Q(P - q - 1, q))
            values = {owner_leading_radii(P, p, q) for p in range(q) if gcd(p, q) == 1}
            require(values == {expected}, f"leading owner radii changed at P={P}, q={q}")
            by_q.append((q, phi(q), *expected))
            total += phi(q) * sum(expected, Q(0)) / P
        require(total == closed_constant(P), f"closed constant mismatch at P={P}")
        require(all(q_equal_P_boundary_trace(P, p) == (set(range(P)), set(range(P)))
                    for p in range(1, P) if gcd(p, P) == 1),
                f"q=P two-sided boundary robustness trace failed at P={P}")
        constant_rows.append((P, total, tuple(by_q)))
        print(f"P={P}: owners={len(owners(P))}; C_P={total}; q=P_boundary_robust=True")

        n0 = (P * P - 1) // 4 if P > 2 else 1
        for q in range(1, P):
            require((q - 1) * (P - q + 1) <= n0,
                    f"quadratic guard inequality failed at P={P},q={q}")
        for p, q in owners(P):
            plus, minus = owner_radii(P, p, q, n0)
            guard_plus, guard_minus = track_guard_radii(P, p, q)
            require(plus <= guard_plus and minus <= guard_minus,
                    f"owner radius escaped its track guard at P={P},p/q={p}/{q}")
        arcs = owner_arcs_at(P, n0)
        for i, (owner_i, pieces_i) in enumerate(arcs):
            for owner_j, pieces_j in arcs[i + 1:]:
                require(not any(max(a, c) < min(b, d)
                                for a, b in pieces_i for c, d in pieces_j),
                        f"positive owner overlap at P={P}: {owner_i},{owner_j}")

    require(closed_constant(2) == Q(1, 2), "P=2 constant changed")
    require(closed_constant(3) == Q(7, 6), "P=3 constant changed")
    require(closed_constant(7) == Q(127, 35), "P=7 constant changed")
    require(owner_leading_radii(2, 0, 1) == (Q(1), Q(0)), "P=2 seam convention changed")
    require(owner_leading_radii(3, 1, 2) == (Q(1, 2), Q(0)), "P=3 q=P-1 side changed")
    print(f"constant_rows_sha256={digest(tuple(constant_rows))}")

    direct_ranges = {
        2: 20, 3: 24, 5: 32, 7: 30, 11: 42, 13: 54, 17: 84,
        19: 103, 23: 140,
    }
    observed = []
    for P, max_m in direct_ranges.items():
        direct = direct_cover_sequence(P, max_m)
        comparisons = []
        for m in range(P, max_m + 1):
            formula = owner_formula_deficit(P, m)
            comparisons.append((m, 1 - direct[m - 1], formula))
        good = [m for m, actual, formula in comparisons if actual == formula]
        onset = next(
            m for m in good
            if all(actual == formula for mm, actual, formula in comparisons if mm >= m)
        )
        require(all(actual == formula for m, actual, formula in comparisons if m >= onset),
                f"sampled tail mismatch at P={P}")
        expected_onset = 2 if P == 2 else (P * P + 3) // 4
        require(onset == expected_onset, f"sampled sharp-onset pattern changed at P={P}")
        require(all(actual != formula for m, actual, formula in comparisons
                    if P <= m < onset), f"pre-onset hostile pattern changed at P={P}")
        observed.append((P, onset, tuple(comparisons)))
        print(f"P={P}: sampled_owner_formula_onset={onset}; checked_through_m={max_m}; "
              f"direct_sha256={digest(direct)}")

    # Known sharp seven-sector controls keep this general audit tied to the
    # separately proved all-depth packet.
    seven = next(row for row in observed if row[0] == 7)
    require(seven[1] == 13, "P=7 sampled sharp onset changed")
    m12 = next(row for row in seven[2] if row[0] == 12)
    require(m12[1] != m12[2], "P=7,m=12 must remain hostile")
    print(f"small_prime_rows_sha256={digest(tuple(observed))}")

    # First composite hostile to the prime local-track theorem.  At H=4 the
    # q=2 owner has a second boundary-start track; a negative perturbation
    # loses its occupied sector 2.  The prime formula undercounts exactly the
    # missing side 1/(4 E_even).
    composite_direct = direct_cover_sequence(4, 20)
    composite_rows = []
    for m in range(9, 21):
        actual = 1 - composite_direct[m - 1]
        naive = owner_formula_deficit(4, m)
        n = m - 1
        E_even = n if n % 2 == 0 else n - 1
        require(actual - naive == Q(1, 4 * E_even),
                f"H=4 composite hostile changed at m={m}")
        composite_rows.append((m, actual, naive, actual - naive))
    print(f"H=4_composite_hostile_sha256={digest(tuple(composite_rows))}; "
          "naive_prime_formula_fails_through_tail=True")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
