#!/usr/bin/env python3
"""Exact hostile checks for THM-2171 (radix quotient-order pumping)."""

from fractions import Fraction
from itertools import combinations, product
from math import gcd


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def state(v, q, j):
    power = q**j
    z = tuple(x // power for x in v)
    r = tuple(x % power for x in v)
    owners = tuple(i for i, x in enumerate(z) if x > 0)
    cuts = tuple(i for i in range(len(v) - 1) if z[i] < z[i + 1])
    return z, r, owners, cuts


def delete_block(v, q, j, k):
    _, rj, _, _ = state(v, q, j)
    zk, _, _, _ = state(v, q, k)
    return tuple(r + q**j * z for r, z in zip(rj, zk))


def circle_distance(x):
    residue = x.numerator % x.denominator
    return Fraction(min(residue, x.denominator - residue), x.denominator)


def exact_maximin(v):
    """Max of the lower envelope, using all affine-cell intersections."""
    times = {Fraction(0)}
    for a in v:
        times.update(Fraction(m, 2 * a) for m in range(2 * a + 1))
    for a, b in combinations(v, 2):
        for den in (a + b, abs(a - b)):
            if den:
                times.update(Fraction(m, den) for m in range(den + 1))
    return max(
        (min(circle_distance(a * t) for a in v), t)
        for t in times
    )


def exhaustive_order_checks():
    transition_checks = 0
    pump_checks = 0
    nontrivial_pumps = 0
    for q in range(2, 6):
        for d in range(2, 6):
            for v in combinations(range(1, 17), d):
                jmax = 0
                while q**jmax <= v[-1]:
                    jmax += 1
                states = [state(v, q, j) for j in range(jmax + 1)]
                for left, right in zip(states, states[1:]):
                    require(set(right[2]) <= set(left[2]), "owner monotonicity")
                    require(set(right[3]) <= set(left[3]), "tie monotonicity")
                    transition_checks += 1
                for j in range(jmax + 1):
                    for k in range(j + 1, jmax + 1):
                        if states[j][2:] != states[k][2:]:
                            continue
                        vp = delete_block(v, q, j, k)
                        require(all(x > 0 for x in vp), "pump positivity")
                        require(
                            all(a < b for a, b in zip(vp, vp[1:])),
                            "pump strict order",
                        )
                        pump_checks += 1
                        nontrivial_pumps += vp != v
    return transition_checks, pump_checks, nontrivial_pumps


def exhaustive_relation_checks():
    checked_relations = 0
    converse_checks = 0
    for q in range(2, 5):
        for v in combinations(range(1, 14), 4):
            jmax = 0
            while q**jmax <= v[-1]:
                jmax += 1
            states = [state(v, q, j) for j in range(jmax + 1)]
            for j in range(jmax + 1):
                for k in range(j + 1, jmax + 1):
                    if states[j][2:] != states[k][2:]:
                        continue
                    vp = delete_block(v, q, j, k)
                    zj, _, _, _ = states[j]
                    zk, _, _, _ = states[k]
                    for a in product(range(-2, 3), repeat=4):
                        if a == (0, 0, 0, 0) or dot(a, v) != 0:
                            continue
                        equal_carry = dot(a, zj) == dot(a, zk)
                        if equal_carry:
                            require(dot(a, vp) == 0, "carry splice")
                            checked_relations += 1
                        if dot(a, vp) == 0:
                            require(equal_carry, "carry splice converse")
                            converse_checks += 1
    return checked_relations, converse_checks


def target_mixed_witness():
    v = (1, 4, 5, 8, 9, 12, 13, 16, 17, 20, 21, 24, 25)
    vp = tuple(range(1, 14))
    r = (0, 2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    s = (0, 3, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0)
    z1, _, o1, t1 = state(v, 2, 1)
    z2, _, o2, t2 = state(v, 2, 2)
    require(delete_block(v, 2, 1, 2) == vp, "target pump")
    require(o1 == o2 == tuple(range(1, 13)), "target owners")
    require(t1 == t2 == (0, 2, 4, 6, 8, 10), "target ties")
    require(
        dot(r, v) == dot(s, v) == dot(r, vp) == dot(s, vp) == 0,
        "target relations",
    )
    require(dot(r, z1) == dot(r, z2) == 0, "target r carries")
    require(dot(s, z1) == dot(s, z2) == 0, "target s carries")
    require(
        (r[3] % 2, s[1] % 2, s[5] % 2) == (1, 1, 1),
        "target mod-2 independence",
    )
    require(
        min(circle_distance(a * Fraction(4, 29)) for a in v)
        == Fraction(3, 29),
        "strict target witness",
    )
    mv, tv = exact_maximin(v)
    mvp, tvp = exact_maximin(vp)
    require(
        (mv, mvp) == (Fraction(3, 29), Fraction(1, 14)),
        "exact target maximin",
    )
    for scale in range(1, 8):
        scaled = tuple(scale * x for x in vp)
        require(exact_maximin(scaled)[0] == mvp, "scale invariance")
        divisor = 0
        for x in scaled:
            divisor = gcd(divisor, x)
        require(tuple(x // divisor for x in scaled) == vp, "gcd normalization")
    return v, vp, r, s, mv, tv, mvp, tvp


def main():
    transitions, pumps, nontrivial = exhaustive_order_checks()
    relation_checks, converse_checks = exhaustive_relation_checks()
    v, vp, r, s, mv, tv, mvp, tvp = target_mixed_witness()
    carry_states = 2729**2
    print("THM-2171 exact hostile audit")
    print(f"monotone owner/tie transitions checked = {transitions}")
    print(f"equal-sidecar block deletions checked = {pumps}")
    print(f"nontrivial positive ordered pumps = {nontrivial}")
    print(f"relation-splice implications checked = {relation_checks}")
    print(f"relation-splice converses checked = {converse_checks}")
    print(f"target witness V = {v}")
    print(f"pumped row V' = {vp}")
    print(f"relations r,s = {r} ; {s}")
    print(f"M(V) = {mv} at t={tv}")
    print(f"M(V') = {mvp} at t={tvp}")
    print(f"rank-two carry states = 2729^2 = {carry_states}")
    print(f"ordered-path state cap = 26*2729^2 = {26 * carry_states}")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
