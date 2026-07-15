#!/usr/bin/env python3
"""Exact verifier for THM-795's lacunary far-flag corollaries.

The script checks the least integer separation factor supplied by the stated
q=333/106 capped-envelope calculation in every far-count stratum f=4,...,13.
It also replays the exact insertion canary showing that safe mass, both
component counts, and the complete component-length multiset do not determine
the next safe mass.

Tournament Analysis is not imposed on the least-constant calculation: its
vertices are ordered rational inequalities, so a binary orientation would add
no invariant.  The companion affine-suspension audit supplies the declared
carrier tournament.  Here the challenged assumption is tested more strongly
by an exact same-fibre/different-output pair: the scalar and length-multiset
carriers tie, while the owner-position/tooth-incidence carrier separates.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations

from lrc14_affine_slope_suspension_codex_S2 import safe_component_counts
from lrc14_certificates import good_intervals


Q = F(333, 106)
A = F(6, 7)


def least_integer(predicate):
    value = 1
    while not predicate(value):
        value += 1
    return value


def geometric_sum(length):
    return sum((A**j for j in range(length)), F(0))


def marked_factor(f):
    inserts = f - 1
    interval = F(f, 98 * (14 - f))
    error = F(2, 7) * geometric_sum(inserts)
    factor = least_integer(lambda r: Q * (r * A**inserts * interval - error) > 1)
    return factor, Q * ((factor - 1) * A**inserts * interval - error), Q * (
        factor * A**inserts * interval - error
    )


def union_factor(f):
    inserts = f - 1
    mass = F(f - 6, 7)
    error = F(2, 7) * geometric_sum(inserts)
    factor = least_integer(lambda r: Q * (r * A**inserts * mass - error) > 1)
    return factor, Q * ((factor - 1) * A**inserts * mass - error), Q * (
        factor * A**inserts * mass - error
    )


def safe_state(speeds):
    intervals = good_intervals(speeds)
    mass = sum((right - left for left, right in intervals), F(0))
    lengths = tuple(sorted(Counter(right - left for left, right in intervals).items()))
    positive, isolated = safe_component_counts(speeds)
    return mass, positive, isolated, lengths


def nine_core_floor():
    best = None
    minimizers = []
    for core in combinations(range(1, 15), 9):
        mass = sum((right - left for left, right in good_intervals(core)), F(0))
        if best is None or mass < best:
            best = mass
            minimizers = [core]
        elif mass == best:
            minimizers.append(core)
    return best, tuple(minimizers)


def main():
    marked = {f: marked_factor(f) for f in range(4, 7)}
    union = {f: union_factor(f) for f in range(7, 13)}

    assert [marked[f][0] for f in range(4, 7)] == [412, 405, 394]
    assert [union[f][0] for f in range(7, 13)] == [27, 17, 14, 13, 13, 13]
    assert all(before <= 1 < after for _, before, after in marked.values())
    assert all(before <= 1 < after for _, before, after in union.values())

    all_far = Q * (
        13 * A**12 - F(2, 7) * sum((A**j for j in range(11)), F(0))
    )
    assert all_far == F(948176665875, 733588221653) > 1

    core_floor, core_minimizers = nine_core_floor()
    assert core_floor == F(10601, 114660)
    assert core_minimizers == ((1, 2, 3, 5, 7, 8, 9, 11, 13),)
    enhanced_before = Q * (18 * A**3 * core_floor - F(254, 343))
    enhanced_after = Q * (19 * A**3 * core_floor - F(254, 343))
    assert enhanced_before == F(55930347, 57900115) < 1
    assert enhanced_after == F(66520746, 57900115) > 1

    base_b = (1, 2, 3, 6, 12)
    base_c = (1, 2, 4, 6, 12)
    state_b = safe_state(base_b)
    state_c = safe_state(base_c)
    assert state_b == state_c
    assert state_b[:3] == (F(4, 7), 12, 0)

    next_b = safe_state(base_b + (59,))
    next_c = safe_state(base_c + (59,))
    assert next_b[0] == F(2425, 4956)
    assert next_c[0] == F(41, 84)
    assert next_b[0] - next_c[0] == F(1, 826)
    assert next_b[1:3] == next_c[1:3] == (44, 0)

    print("THM-795 lacunary far-flag exact certificate")
    print("marked-core factors (f=4..6):")
    for f, (factor, before, after) in marked.items():
        print(f"  f={f}: R={factor}, q-test(R-1)={before}, q-test(R)={after}")
    print("union-mass factors (f=7..12):")
    for f, (factor, before, after) in union.items():
        print(f"  f={f}: R={factor}, q-test(R-1)={before}, q-test(R)={after}")
    print(f"all-far factor: f=13, R=13, q-test={all_far}")
    print("factor vector:", [marked[f][0] for f in range(4, 7)] + [union[f][0] for f in range(7, 13)] + [13])
    print("finite-exact nine-core enhancement:")
    print(f"  min_mu={core_floor}, unique_core={core_minimizers[0]}")
    print(f"  R=18 q-test={enhanced_before}, R=19 q-test={enhanced_after}")
    print("phase-placement no-go:")
    print(f"  common base state: mu={state_b[0]}, r_plus={state_b[1]}, isolated={state_b[2]}")
    print(f"  common length multiset: {state_b[3]}")
    print(f"  after N=59: mu_B={next_b[0]}, mu_C={next_c[0]}, difference={next_b[0]-next_c[0]}")
    print(f"  after N=59: r_plus={next_b[1]}, isolated={next_b[2]}")
    print("PASS: fully lacunary flags are terminal; exact insertion still needs phase incidence.")


if __name__ == "__main__":
    main()
