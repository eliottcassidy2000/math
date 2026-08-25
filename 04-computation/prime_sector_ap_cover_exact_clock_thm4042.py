#!/usr/bin/env python3
"""Exact eventual phase audit for THM-4042's prime-sector AP-cover clock.

This script does not assert an onset for the finite-owner formula.  It starts
from the local max-min owner law and selects its eventual rational winner on
each residue phase by the exact order at infinity of C/(n-c).
"""

from fractions import Fraction as Q
from math import floor, gcd, lcm


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def minimal_period(values):
    n = len(values)
    return next(
        d for d in divisors(n)
        if all(values[i] == values[i % d] for i in range(n))
    )


def predicted_word_period(P, q):
    """Closed p-adic formula for the least translation period."""
    if q == 1:
        return 1
    remaining = q
    prime = 2
    answer = 1
    while prime * prime <= remaining:
        if remaining % prime:
            prime += 1
            continue
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if q == P - 1:
            period_exponent = 1
        else:
            valuation = 0
            value = P - 1
            while value % prime == 0:
                value //= prime
                valuation += 1
            period_exponent = min(exponent, valuation + 1)
        answer *= prime ** period_exponent
        prime += 1
    if remaining > 1:
        answer *= remaining
    return answer


def owners_of_denominator(q):
    return tuple(a for a in range(q) if gcd(a, q) == 1)


def eventual_radius_term(P, a, q, phase, side):
    """Return the eventual winning (C,c) for C/(n-c), n=phase mod q."""
    data = []
    occupied = set()
    for s in range(q):
        A = Q(P * a * s, q)
        integer = floor(A)
        frac = A - integer
        residue = integer % P
        occupied.add(residue)
        shift = (phase - s) % q
        data.append((residue, frac, shift))
    missing = sorted(set(range(P)) - occupied)
    require(len(occupied) == q and len(missing) == P - q, (P, a, q, "orbit"))

    inner = []
    for r in missing:
        candidates = []
        for residue, frac, shift in data:
            if side == "plus":
                coefficient = Q((r - residue) % P) - frac
            else:
                coefficient = Q(((residue - r) % P) - 1) + frac
            require(coefficient >= 0, (P, a, q, phase, side, coefficient))
            candidates.append((coefficient, shift))
        # For C>0, eventual order is lexicographic in (C,c), since
        # C/(n-c)=C/n+Cc/n^2+... .  Zero candidates are identical.
        positive = [term for term in candidates if term[0] > 0]
        if len(positive) != len(candidates):
            winner = (Q(0), 0)
        else:
            winner = min(candidates)
        inner.append(winner)
    return max(inner)


def denominator_block(P, q):
    """Aggregate all reduced owners of denominator q, retaining shifts."""
    rows = []
    for phase in range(q):
        row = [Q(0) for _ in range(P - 1)]
        for a in owners_of_denominator(q):
            for side in ("minus", "plus"):
                coefficient, shift = eventual_radius_term(P, a, q, phase, side)
                row[shift] += coefficient / P
        rows.append(tuple(row))
    return tuple(rows)


def explicit_base_word(P, q):
    """Closed word whose cyclic shifts are the denominator-q phase block."""
    if q == 1:
        return (Q(2 * P - 3, P),)
    units = tuple(u for u in range(q) if gcd(u, q) == 1)
    inverse_P = pow(P, -1, q)
    multiplier = (inverse_P - 1) % q
    negative_multiplicity = [0 for _ in range(q)]
    for u in units:
        negative_multiplicity[(multiplier * u) % q] += 1
    plus_weight = Q(P - q, P * q)
    minus_weight = Q(P - q - 1, P * q)
    return tuple(
        plus_weight * (c in units) + minus_weight * negative_multiplicity[c]
        for c in range(q)
    )


def explicit_block(P, q):
    word = explicit_base_word(P, q)
    return tuple(
        tuple(word[(c - phase) % q] if c < q else Q(0) for c in range(P - 1))
        for phase in range(q)
    )


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def asymptotic_radius(P, a, q, side):
    terms = [eventual_radius_term(P, a, q, phase, side)[0] for phase in range(q)]
    require(len(set(terms)) == 1, (P, a, q, side, terms))
    return terms[0]


def global_rows(P, blocks):
    L = lcm(*range(1, P))
    return tuple(
        tuple(
            sum((blocks[q][r % q][c] for q in range(1, P)), Q(0))
            for c in range(P - 1)
        )
        for r in range(L)
    )


def main():
    print("PRIME OWNER ASYMPTOTIC RADII")
    for P in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        total = Q(0)
        for q in range(1, P):
            expected_minus = Q(P - q - 1, q)
            expected_plus = Q(P - q, q)
            for a in owners_of_denominator(q):
                minus = asymptotic_radius(P, a, q, "minus")
                plus = asymptotic_radius(P, a, q, "plus")
                require(minus == expected_minus, (P, a, q, "minus", minus))
                require(plus == expected_plus, (P, a, q, "plus", plus))
                total += minus + plus
        kappa = total / P
        formula = sum(
            (
                Q(len(owners_of_denominator(q)) * (2 * P - 2 * q - 1), q)
                for q in range(1, P)
            ),
            Q(0),
        ) / P
        require(kappa == formula, (P, "leading constant"))
        print(f"P={P};kappa={kappa};radius_formula=PASS")

    print("\nEVENTUAL DENOMINATOR CLOCKS")
    for P in (2, 3, 5, 7, 11, 13):
        blocks = {q: denominator_block(P, q) for q in range(1, P)}
        active_shifts = set()
        maximal_shift_periods = []
        for q, rows in blocks.items():
            require(rows == explicit_block(P, q), (P, q, "closed clock word"))
            vector_period = minimal_period(rows)
            require(vector_period == predicted_word_period(P, q),
                    (P, q, vector_period, "p-adic word period"))
            maximal_shift = tuple(row[q - 1] for row in rows)
            maximal_period = minimal_period(maximal_shift)
            maximal_shift_periods.append(maximal_period)
            active_shifts.update(
                c for row in rows for c, coefficient in enumerate(row) if coefficient > 0
            )
            print(
                f"P={P};q={q};vector_period={vector_period};"
                f"top_shift_period={maximal_period};"
                f"top_shift_support={sum(value > 0 for value in maximal_shift)}"
            )
        L = lcm(*range(1, P))
        rows = global_rows(P, blocks)
        global_period = minimal_period(rows)
        require(active_shifts == set(range(P - 1)), (P, "active shifts", active_shifts))
        print(
            f"P={P};L={L};global_vector_period={global_period};"
            f"phase_bound_sharp={global_period == L};"
            f"all_shifts_active={sorted(active_shifts)};"
            f"top_shift_periods={tuple(maximal_shift_periods)}"
        )

    print("\nCLOSED CLOCK-PERIOD PROFILE")
    for P in tuple(P for P in range(2, 62) if is_prime(P)):
        periods = tuple(minimal_period(explicit_base_word(P, q)) for q in range(1, P))
        predicted = tuple(predicted_word_period(P, q) for q in range(1, P))
        require(periods == predicted, (P, periods, predicted, "closed period profile"))
        exact_phase = lcm(*periods)
        upper_phase = lcm(*range(1, P))
        print(
            f"P={P};Pi={exact_phase};L={upper_phase};L_over_Pi={upper_phase // exact_phase};"
            f"delta={periods}"
        )

    print("RESULT=PASS")


if __name__ == "__main__":
    main()
