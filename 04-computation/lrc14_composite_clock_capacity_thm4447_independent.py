#!/usr/bin/env python3
"""Independent exact referee for THM-4447.

This does not import the primary audit or repository geometry code. It audits
the stronger arbitrary-tail form at exact pack clock c=gcd(P).
"""

from fractions import Fraction as F
from itertools import product
from math import gcd


R = F(1, 14)
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def floor_f(x):
    return x.numerator // x.denominator


def ceil_f(x):
    return -floor_f(-x)


def norm(x):
    z = x % 1
    return min(z, 1 - z)


def divs(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def bad(q, t, y):
    return frozenset(j for j in range(q) if norm(F(t, q) * (y + j)) < R)


def exact_count_formula(q, t, y):
    g = gcd(q, t)
    m = q // g
    delta = F(t, q) * y
    a = -m * delta - F(m, 14)
    return g * (ceil_f(a + F(m, 7)) - floor_f(a) - 1)


def capacity(q, t):
    g = gcd(q, t)
    return g * ceil_f(F(q, 7 * g))


def grid_count(m, delta):
    return sum(norm(delta + F(k, m)) < R for k in range(m))


def orbit_formula(m, delta):
    a = -m * delta - F(m, 14)
    return ceil_f(a + F(m, 7)) - floor_f(a) - 1


def translate_samples(m):
    walls = {F(0), F(1)}
    for k in range(m):
        walls.add((-F(k, m) - R) % 1)
        walls.add((-F(k, m) + R) % 1)
    ordered = sorted(walls)
    mids = {(a + b) / 2 for a, b in zip(ordered, ordered[1:])}
    return tuple(sorted(set(ordered) | mids))


def audit_count_formula():
    cells = 0
    resonances = 0
    for m in range(1, 141):
        cap = ceil_f(F(m, 7))
        attained = False
        for delta in translate_samples(m):
            direct = grid_count(m, delta)
            formula = orbit_formula(m, delta)
            check(direct == formula, ("orbit formula", m, delta, direct, formula))
            check(cap - 1 <= direct <= cap, ("two levels", m, delta, direct, cap))

            a = -m * delta - F(m, 14)
            frac_a = a - floor_f(a)
            length = F(m, 7)
            frac_l = length - floor_f(length)
            if frac_l == 0:
                predicted = cap - (frac_a == 0)
                endpoint = frac_a == 0
            else:
                predicted = cap if frac_a > 1 - frac_l else cap - 1
                endpoint = frac_a == 0 or frac_a == 1 - frac_l
            literal_endpoint = any(norm(delta + F(k, m)) == R for k in range(m))
            check(direct == predicted, ("chamber law", m, delta, direct, predicted))
            check(endpoint == literal_endpoint, ("resonance law", m, delta))
            attained |= direct == cap
            resonances += literal_endpoint
            cells += 1
        check(attained, ("capacity unattained", m, cap))

    literal = 0
    ys = (F(0), F(1, 14), F(1, 10), F(1, 7), F(5, 22), F(13, 37), F(86, 539))
    for q in range(2, 85):
        for t in range(1, 2 * q + 2):
            for y in ys:
                direct = len(bad(q, t, y))
                check(direct == exact_count_formula(q, t, y), ("label formula", q, t, y))
                check(direct <= capacity(q, t), ("capacity", q, t, y))
                literal += 1
    return cells, resonances, literal


def divisor_certificate(q, signature):
    candidates = []
    for d in divs(q):
        if d == 1:
            continue
        absorbed = tuple(i for i, g in enumerate(signature) if g % d == 0)
        if len(absorbed) == 3:
            continue
        residual = tuple(i for i in range(3) if i not in absorbed)
        total = 0
        for i in residual:
            gg = gcd(signature[i], d)
            total += gg * ceil_f(F(d, 7 * gg))
        candidates.append((total - d, d, absorbed, total))
    return min(candidates) if candidates else None


def audit_classification():
    unresolved = {}
    signatures_checked = 0
    for q in range(2, 513):
        possible = divs(q)
        local = []
        for sig in product(possible, repeat=3):
            if gcd(gcd(sig[0], sig[1]), sig[2]) != 1:
                continue
            signatures_checked += 1
            cert = divisor_certificate(q, sig)
            if cert is None or cert[0] >= 0:
                local.append(sig)
        if local:
            unresolved[q] = tuple(sorted(local))

    expected = {
        2: ((1, 1, 1), (1, 1, 2), (1, 2, 1), (2, 1, 1)),
        3: ((1, 1, 1),),
        4: ((1, 1, 2), (1, 2, 1), (2, 1, 1)),
    }
    check(unresolved == expected, ("residual classification", unresolved))
    return signatures_checked, unresolved


def audit_divisor_rewrite():
    typed_cases = 0
    for q in range(2, 61):
        body = tuple(range(1, 11)) if q % 2 else tuple(2 * x + 1 for x in range(1, 11))
        pack = {q * c for c in body}
        candidates = [t for t in range(1, 6 * q + 40) if t not in pack]
        triples = []
        for i in range(min(8, len(candidates) - 2)):
            triple = (candidates[i], candidates[i + 1], candidates[-1 - i])
            primitive = gcd(gcd(gcd(*pack), triple[0]), gcd(triple[1], triple[2])) == 1
            if len(set(triple)) == 3 and primitive:
                triples.append(triple)
        for tails in triples:
            physical = pack | set(tails)
            check(len(physical) == 13, ("physical distinctness", q, tails))
            check(gcd(*physical) == 1, ("physical primitivity", q, tails))
            for d in divs(q):
                if d == 1:
                    continue
                absorbed = tuple(t for t in tails if t % d == 0)
                residual = tuple(t for t in tails if t % d)
                reduced = {(q // d) * c for c in body} | {t // d for t in absorbed}
                check(len(absorbed) <= 2, ("three absorbed", q, d, tails))
                check(len(reduced) == 10 + len(absorbed) <= 12, ("pack size", q, d))
                check(
                    {d * r for r in reduced} | set(residual) == physical,
                    ("physical rewrite", q, d, tails),
                )
                typed_cases += 1
    return typed_cases


def audit_equality_and_endpoints():
    examples = {
        (3, (1, 4, 5), F(23, 112)): (
            frozenset({0}),
            frozenset({2}),
            frozenset({1}),
        ),
        (4, (1, 11, 14), F(86, 539)): (
            frozenset({0}),
            frozenset({2}),
            frozenset({1, 3}),
        ),
    }
    records = []
    for (q, tails, y), expected in examples.items():
        masks = tuple(bad(q, t, y) for t in tails)
        check(masks == expected, ("equality masks", q, tails, y, masks))
        check(set().union(*masks) == set(range(q)), ("full cover", q, masks))
        check(sum(capacity(q, t) for t in tails) == q, ("capacity equality", q))
        check(all(len(mask) == capacity(q, t) for t, mask in zip(tails, masks)), ("tight", q))
        check(sum(map(len, masks)) == len(set().union(*masks)), ("partition", q))
        records.append((q, tails, y, masks))

    absorbed_records = []
    absorbed_examples = {
        (2, (1, 7, 8), F(13, 98)),
        (4, (1, 11, 14), F(1, 11)),
    }
    for original_q, tails, y in sorted(absorbed_examples):
        d = 2
        absorbed = tuple(t for t in tails if t % d == 0)
        residual = tuple(t for t in tails if t % d)
        check(len(absorbed) == 1 and len(residual) == 2, ("absorption shape", original_q))
        check(
            all(norm(F(t, d) * (y + j)) >= R for t in absorbed for j in range(d)),
            ("absorbed tail unsafe", original_q),
        )
        masks = tuple(bad(d, t, y) for t in residual)
        check(set().union(*masks) == {0, 1}, ("residual cover", original_q))
        check(sum(capacity(d, t) for t in residual) == d, ("residual equality", original_q))
        absorbed_records.append((original_q, tails, y, absorbed, residual, masks))

    endpoint_records = []
    controls = (
        (7, 1, F(1, 2), True),
        (14, 2, F(1, 2), True),
        (4, 1, F(2, 7), True),
        (8, 2, F(2, 7), True),
        (9, 1, F(0), False),
    )
    for q, t, y, should_resonate in controls:
        g = gcd(q, t)
        mask = bad(q, t, y)
        check(len(mask) == capacity(q, t) - g, ("gcd deficit", q, t, y))
        resonance = any(norm(F(t, q) * (y + j)) == R for j in range(q))
        check(resonance == should_resonate, ("resonance control", q, t, y))
        endpoint_records.append((q, t, y, g, capacity(q, t), mask, resonance))
    return records, absorbed_records, endpoint_records


def main():
    cells, resonances, literal = audit_count_formula()
    typed = audit_divisor_rewrite()
    signatures, unresolved = audit_classification()
    equality, absorbed, endpoints = audit_equality_and_endpoints()
    print("INDEPENDENT_COMPOSITE_CLOCK_REFEREE")
    print("orbit_translate_samples", cells, "resonances", resonances)
    print("literal_label_cases", literal)
    print("typed_divisor_rewrites", typed)
    print("gcd_signatures_q_le_512", signatures)
    print("unresolved", unresolved)
    print("equality", equality)
    print("absorbed_equality", absorbed)
    print("endpoint_deficits", endpoints)
    print("checks", CHECKS)
    print("PASS")


if __name__ == "__main__":
    main()
