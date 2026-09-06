#!/usr/bin/env python3
"""Exact controls for THM-4447's 1/14 composite-clock capacity."""

from fractions import Fraction as Q
from itertools import product
from math import gcd

CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def ceilq(x):
    return -((-x.numerator) // x.denominator)


def floorq(x):
    return x.numerator // x.denominator


def dist(x):
    z = x % 1
    return min(z, 1 - z)


def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def cap_from_g(q, g):
    need(q % g == 0 and g > 0, ("cap typing", q, g))
    m = q // g
    return g * ((m + 6) // 7)


def cap(q, t):
    return cap_from_g(q, gcd(q, abs(t)))


def orbit_count(m, delta):
    return sum(dist(delta + Q(k, m)) < Q(1, 14) for k in range(m))


def orbit_formula(m, delta):
    a = -m * delta - Q(m, 14)
    return ceilq(a + Q(m, 7)) - floorq(a) - 1


def bad_labels(q, t, y):
    return frozenset(
        j for j in range(q) if dist(Q(t, q) * (y + j)) < Q(1, 14)
    )


def event_cells(q, tails):
    walls = {Q(0), Q(1)}
    for t in tails:
        for j in range(q):
            for n in range(-1, abs(t) + 2):
                for sign in (-1, 1):
                    y = Q(q, t) * (Q(n) + Q(sign, 14)) - j
                    if 0 <= y <= 1:
                        walls.add(y)
    walls = sorted(walls)
    samples = set(walls)
    samples.update((a + b) / 2 for a, b in zip(walls, walls[1:]))
    return tuple(sorted(samples))


def best_divisor_certificate(q, gs):
    """Absorb d-divisible tails, then count the residual d-sheet masks."""
    best = None
    for d in divisors(q):
        if d < 2:
            continue
        absorbed = tuple(i for i, g in enumerate(gs) if g % d == 0)
        if len(absorbed) > 2:
            continue
        residual = tuple(i for i in range(3) if i not in absorbed)
        capacity = sum(cap_from_g(d, gcd(gs[i], d)) for i in residual)
        row = (capacity - d, d, absorbed, residual, capacity)
        if best is None or row < best:
            best = row
    return best


def main():
    # Every event wall and intervening chamber for order m <= 140.
    orbit_cases = 0
    low_count_cases = 0
    for m in range(1, 141):
        boundaries = {Q(0), Q(1)}
        for k in range(m):
            boundaries.add((-Q(k, m) + Q(1, 14)) % 1)
            boundaries.add((-Q(k, m) - Q(1, 14)) % 1)
        walls = sorted(boundaries)
        samples = set(walls)
        samples.update((a + b) / 2 for a, b in zip(walls, walls[1:]))
        for delta in samples:
            n = orbit_count(m, delta)
            f = orbit_formula(m, delta)
            top = (m + 6) // 7
            need(n == f, ("orbit formula", m, delta, n, f))
            need(n in (top, top - 1), ("two count levels", m, delta, n))
            need(n <= top, ("sharp cap", m, delta, n))
            orbit_cases += 1
            low_count_cases += n < top

    # Literal q,t,y counts, now including every gcd/3-adic tail type.
    literal_cases = 0
    for q in range(2, 61):
        for t in range(1, 2 * q + 1):
            g = gcd(q, t)
            m = q // g
            for y in (Q(0), Q(1, 14), Q(1, 10), Q(1, 7), Q(5, 22), Q(13, 37)):
                labels = bad_labels(q, t, y)
                need(
                    len(labels) == g * orbit_count(m, Q(t, q) * y),
                    ("gcd orbit multiplicity", q, t, y, labels),
                )
                need(len(labels) <= cap(q, t), ("literal cap", q, t, y, labels))
                literal_cases += 1

    # Exhaust all primitive triples of gcd signatures through q=210.  This
    # intentionally includes tails divisible by three; the theorem is not
    # restricted to ternary units.  Full-row primitivity is gcd(gs)=1.
    unresolved = {}
    signature_cases = 0
    for q in range(2, 211):
        ds = divisors(q)
        bad = []
        for gs in product(ds, repeat=3):
            if gcd(gcd(gs[0], gs[1]), gs[2]) != 1:
                continue
            signature_cases += 1
            cert = best_divisor_certificate(q, gs)
            if cert[0] >= 0:
                bad.append(gs)
        if bad:
            unresolved[q] = tuple(sorted(bad))
    expected = {
        2: ((1, 1, 1), (1, 1, 2), (1, 2, 1), (2, 1, 1)),
        3: ((1, 1, 1),),
        4: ((1, 1, 2), (1, 2, 1), (2, 1, 1)),
    }
    need(unresolved == expected, ("complete small-clock residual", unresolved))

    # Pointwise covers show that the remaining scalar equalities are real
    # method boundaries.  They are pack phases, not non-lonely full rows.
    cover_controls = {
        2: ((1, 5, 7), Q(13, 98)),
        3: ((1, 4, 5), Q(23, 112)),
        4: ((1, 11, 14), Q(86, 539)),
    }
    covers = {}
    for q, (tails, y) in cover_controls.items():
        masks = tuple(bad_labels(q, t, y) for t in tails)
        need(set().union(*masks) == set(range(q)), ("hostile cover", q, tails, y, masks))
        covers[q] = (tails, y, masks)

    # After divisor-two absorption, two odd tails can cover both labels.
    absorbed_controls = {
        2: ((1, 7, 8), Q(13, 98)),
        4: ((1, 11, 14), Q(1, 11)),
    }
    absorbed_records = {}
    for original_q, (tails, y) in absorbed_controls.items():
        masks = tuple(bad_labels(2, t, y) for t in tails)
        even = [i for i, t in enumerate(tails) if t % 2 == 0]
        odd = [i for i, t in enumerate(tails) if t % 2]
        need(
            len(even) == 1 and not masks[even[0]],
            ("absorbed tail is pack-safe", original_q, tails, y, masks),
        )
        need(
            set().union(*(masks[i] for i in odd)) == {0, 1},
            ("two residual tails cover", original_q, tails, y, masks),
        )
        absorbed_records[original_q] = (tails, y, masks)

    # A transition wall can lower the count by one orbit point, hence g
    # labels.  For nonintegral m/7, the same low count persists throughout
    # a whole adjacent chamber; it is not exclusively an endpoint event.
    endpoint_controls = []
    for q, t in ((7, 1), (14, 2), (4, 1), (8, 2), (9, 1)):
        g = gcd(q, t)
        m = q // g
        found = None
        for y in event_cells(q, (t,)):
            labels = bad_labels(q, t, y)
            if len(labels) == cap(q, t) - g:
                found = (y, labels)
                break
        need(found is not None, ("low-count witness", q, t))
        endpoint_controls.append((q, t, g, m, cap(q, t), found))

    print("COMPOSITE_CLOCK_CAPACITY_AUDIT")
    print("CAP B_q(t)=g*ceil(q/(7g)), g=gcd(q,t)")
    print("orbit_cases", orbit_cases, "low_count_cases", low_count_cases)
    print("literal_q_t_y_cases", literal_cases)
    print("signature_cases_q_le_210", signature_cases)
    print("primitive_ten_pack_unresolved", unresolved)
    print("pointwise_cover_controls", covers)
    print("absorbed_divisor_two_controls", absorbed_records)
    print("low_count_transition_controls", endpoint_controls)
    print("checks", CHECKS)
    print("PASS")


if __name__ == "__main__":
    main()
