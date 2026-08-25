#!/usr/bin/env python3
"""Independent exact checks for THM-4042's prime-sector clock law.

This file intentionally does not import the theorem's companion.  It
rebuilds owner data, nested asymptotic selectors, cyclic blocks, direct wall
measures, and finite-difference controls from the definitions.
"""

from fractions import Fraction as F
from math import floor, gcd, lcm


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def units(q):
    return tuple(a for a in range(q) if gcd(a, q) == 1)


def owner_data(P, a, q):
    rows = []
    for s in range(q):
        value = F(P * a * s, q)
        integer = floor(value)
        rows.append((s, integer % P, value - integer))
    return tuple(rows)


def nested_leading_winner(P, a, q, side):
    """Return outer leading coefficient and all winning (r,s) pairs.

    This computes max_r min_s on the coefficient before considering any
    horizon phase.  For a positive coefficient the claimed geometry says the
    outer r and inner s are both unique.
    """
    data = owner_data(P, a, q)
    occupied = {b for _, b, _ in data}
    require(len(occupied) == q, (P, a, q, "distinct occupied sectors"))
    rows = []
    for r in set(range(P)) - occupied:
        candidates = []
        for s, b, f in data:
            if side == "plus":
                coefficient = F((r - b) % P) - f
            else:
                coefficient = F((b - r) % P - 1) + f
            candidates.append((coefficient, s))
        inner = min(value for value, _ in candidates)
        rows.append((inner, r, tuple(s for value, s in candidates if value == inner)))
    outer = max(value for value, _, _ in rows)
    winners = tuple((r, tracks) for value, r, tracks in rows if value == outer)
    return outer, winners


def phase_selected_term(P, a, q, phase, side):
    """Direct lexicographic order at infinity of the nested C/(n-c) law."""
    data = owner_data(P, a, q)
    occupied = {b for _, b, _ in data}
    inner = []
    for r in set(range(P)) - occupied:
        terms = []
        for s, b, f in data:
            coefficient = (
                F((r - b) % P) - f
                if side == "plus"
                else F((b - r) % P - 1) + f
            )
            shift = (phase - s) % q
            terms.append((coefficient, shift))
        if any(coefficient == 0 for coefficient, _ in terms):
            inner.append((F(0), 0))
        else:
            inner.append(min(terms))
    return max(inner)


def direct_denominator_rows(P, q):
    result = []
    for phase in range(q):
        row = [F(0) for _ in range(P - 1)]
        for a in units(q):
            for side in ("minus", "plus"):
                coefficient, shift = phase_selected_term(P, a, q, phase, side)
                row[shift] += coefficient / P
        result.append(tuple(row))
    return tuple(result)


def closed_word(P, q):
    if q == 1:
        return (F(2 * P - 3, P),)
    inverse = pow(P, -1, q)
    t = (inverse - 1) % q
    ans = []
    for c in range(q):
        plus = int(gcd(c, q) == 1)
        minus = sum(1 for u in units(q) if t * u % q == c)
        ans.append(F(P - q, P * q) * plus + F(P - q - 1, P * q) * minus)
    return tuple(ans)


def closed_rows(P, q):
    word = closed_word(P, q)
    return tuple(
        tuple(word[(c - phase) % q] if c < q else F(0) for c in range(P - 1))
        for phase in range(q)
    )


def cyclic_period(values):
    length = len(values)
    for d in range(1, length + 1):
        if all(values[i] == values[(i + d) % length] for i in range(length)):
            return d
    raise AssertionError("finite cyclic word has no period")


def predicted_word_period(P, q):
    if q == 1:
        return 1
    answer = 1
    remaining = q
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime:
            prime += 1
            continue
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        valuation = 0
        value = P - 1
        while value % prime == 0:
            value //= prime
            valuation += 1
        answer *= prime ** (1 if q == P - 1 else min(exponent, valuation + 1))
        prime += 1
    if remaining > 1:
        answer *= remaining
    return answer


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def formula_value(P, n):
    total = F(0)
    for q in range(1, P):
        word = closed_word(P, q)
        for c in range(q):
            total += word[(c - n) % q] / (n - c)
    return total


def falling(n, degree):
    answer = 1
    for c in range(degree):
        answer *= n - c
    return answer


def covers(theta, P, m):
    return len({floor(e * theta) % P for e in range(m)}) == P


def direct_deficit(P, m):
    walls = {F(0), F(P)}
    for e in range(1, m):
        walls.update(F(k, e) for k in range(P * e + 1))
    walls = sorted(walls)
    theta_length = sum(
        (
            right - left
            for left, right in zip(walls, walls[1:])
            if not covers((left + right) / 2, P, m)
        ),
        F(0),
    )
    return theta_length / P


def finite_owner_radius(P, a, q, n, side):
    data = []
    occupied = set()
    for s, b, f in owner_data(P, a, q):
        require(s <= n, (P, a, q, n, "track not yet available"))
        E = n - ((n - s) % q)
        require(E > 0, (P, a, q, n, s, "positive terminal time"))
        data.append((b, f, E))
        occupied.add(b)
    radii = []
    for r in set(range(P)) - occupied:
        costs = []
        for b, f, E in data:
            numerator = (
                F((r - b) % P) - f
                if side == "plus"
                else F((b - r) % P - 1) + f
            )
            costs.append(numerator / E)
        radii.append(min(costs))
    return max(radii)


def finite_owner_deficit(P, m):
    n = m - 1
    total = F(0)
    for q in range(1, P):
        for a in units(q):
            total += finite_owner_radius(P, a, q, n, "minus")
            total += finite_owner_radius(P, a, q, n, "plus")
    return total / P


def main():
    winner_checks = 0
    block_checks = 0
    for P in (2, 3, 5, 7, 11, 13, 17, 19):
        for q in range(2, P):
            for a in units(q):
                for side in ("plus", "minus"):
                    coefficient, winners = nested_leading_winner(P, a, q, side)
                    expected_coefficient = F(P - q - (side == "minus"), q)
                    expected_track = (
                        -pow(a, -1, q)
                        if side == "plus"
                        else (1 - pow(P, -1, q)) * pow(a, -1, q)
                    ) % q
                    require(coefficient == expected_coefficient, (P, q, a, side, "coefficient"))
                    require(any(expected_track in tracks for _, tracks in winners), (P, q, a, side, "track"))
                    if coefficient > 0:
                        require(len(winners) == 1 and len(winners[0][1]) == 1, (P, q, a, side, "uniqueness"))
                    winner_checks += 1

        if P <= 13:
            for q in range(1, P):
                require(direct_denominator_rows(P, q) == closed_rows(P, q), (P, q, "word convolution"))
                block_checks += 1

    expected_periods = {2: 1, 3: 2, 5: 6, 7: 60, 11: 420, 13: 27720, 17: 120120}
    period_formula_checks = 0
    for P in (p for p in range(2, 252) if is_prime(p)):
        for q in range(1, P):
            observed = cyclic_period(closed_word(P, q))
            require(observed == predicted_word_period(P, q),
                    (P, q, observed, "p-adic word period"))
            period_formula_checks += 1
    for P, expected in expected_periods.items():
        deltas = tuple(cyclic_period(closed_word(P, q)) for q in range(1, P))
        require(lcm(*deltas) == expected, (P, deltas, "global period"))
        top_expected = 1 if P == 2 else 1
        for prime in range(2, P):
            if (P - 1) % prime == 0 and all(prime % d for d in range(2, int(prime ** 0.5) + 1)):
                top_expected *= prime
        require(deltas[-1] == top_expected, (P, deltas[-1], top_expected, "top radical"))

        # Every affine pole is active already in a translate of the top word.
        top = closed_word(P, P - 1)
        active = {
            c
            for phase in range(P - 1)
            for c in range(P - 1)
            if top[(c - phase) % (P - 1)] > 0
        }
        require(active == set(range(P - 1)), (P, "all poles active"))

        # The (P-1)-st Pi-step difference of QD vanishes identically on the
        # eventual rational phase controller.
        Pi = expected
        for n in range(P, P + 2 * Pi + 1, max(1, Pi // 7)):
            lhs = F(0)
            for j in range(P):
                nn = n + j * Pi
                sign = -1 if (P - 1 - j) % 2 else 1
                # local exact binomial, avoiding an import from math.comb
                choose = 1
                for k in range(1, j + 1):
                    choose = choose * (P - k) // k
                lhs += sign * choose * falling(nn, P - 1) * formula_value(P, nn)
            require(lhs == 0, (P, n, "finite difference"))

    finite_windows = {
        3: (range(3, 9), 3),
        5: (range(5, 13), 7),
        7: (range(7, 20), 13),
    }
    for P, (window, onset) in finite_windows.items():
        bits = tuple(direct_deficit(P, m) == finite_owner_deficit(P, m) for m in window)
        require(bits == tuple(m >= onset for m in window), (P, bits, "finite onset window"))

    direct4 = direct_deficit(4, 20)
    naive4 = finite_owner_deficit(4, 20)
    require(direct4 == F(299, 2907), (direct4, "composite direct"))
    require(naive4 == F(2069, 23256), (naive4, "composite naive"))
    require(direct4 - naive4 == F(1, 72), (direct4 - naive4, "composite gap"))

    print(f"winner_formula_checks={winner_checks}")
    print(f"direct_word_block_checks={block_checks}")
    print(f"p_adic_period_checks={period_formula_checks}")
    print(f"period_profile={expected_periods}")
    print("finite_prime_windows=P3,P5,P7")
    print("composite_P4_m20_gap=1/72")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
