#!/usr/bin/env python3
"""Concise frozen audit for the AP-cover finite-owner theorem packet.

Uses explicit require() rather than assert so normal and python -O executions
perform the same semantic checks.
"""

from fractions import Fraction as Q
from hashlib import sha256
import sys

from lrc14_ap_cover_components_thm4029 import (
    circular_contains,
    local_missing_sector_radii,
    noncover_components,
)
from lrc14_ap_cover_sequence_engine_thm4029 import (
    a_direct,
    exact_caps,
    farey_cover_distribution,
)
from lrc14_ap_cover_finite_owner_formula_thm4029 import (
    audit_m13_zero_margin_endpoints,
    direct_component_deficit,
    owner_formula_deficit,
    owners,
    prove_local_guards,
    prove_phase_rational_law,
)


sys.stdout.reconfigure(newline="\n")


def require(condition: bool, label: str):
    if not condition:
        raise RuntimeError(label)


def digest(obj) -> str:
    return sha256(repr(obj).encode("utf-8")).hexdigest()


def audit_owner_component_bijection(m: int):
    """Assert the consequence object component by component, not just in total."""
    components = noncover_components(m)
    owner_thetas = tuple(7 * owner for owner in owners())
    require(len(components) == 12, f"m={m}: component count")
    require(len(owner_thetas) == len(set(owner_thetas)) == 12, "owner count")
    seen = set()
    records = []
    for left, right, _ in components:
        hits = tuple(
            theta0 for theta0 in owner_thetas if circular_contains(left, right, theta0)
        )
        require(len(hits) == 1, (m, left, right, "owner/component bijection"))
        theta0 = hits[0]
        lifted = theta0 + (7 if right > 7 and theta0 < left else 0)
        owner = (theta0 % 7) / 7
        minus, plus, _, _ = local_missing_sector_radii(
            owner.numerator, owner.denominator, m
        )
        actual = (lifted - left, right - lifted)
        require(actual == (minus, plus), (m, owner, "individual radii"))
        seen.add(theta0)
        records.append((owner, left, right, minus, plus))
    require(seen == set(owner_thetas), f"m={m}: every owner used once")
    return tuple(sorted(records))


def main():
    print("AP COVER THEOREM-READY AUDIT")

    direct = tuple(a_direct(m) for m in range(1, 31))
    mass, unresolved, witnesses = farey_cover_distribution(200)
    cumulative = Q(0)
    farey = []
    for m in range(1, 201):
        cumulative += mass[m]
        farey.append(cumulative)
    farey = tuple(farey)
    require(direct == farey[:30], "direct x-wall and Farey/Sturmian engines disagree")
    require(all(m in witnesses for m in range(7, 201)), "missing first-cover witness")
    print(f"direct_x_wall_m1_30_sha256={digest(direct)}")
    print(f"farey_cover_time_m1_200_sha256={digest(farey)}")
    print(f"independent_engines_equal_m1_30={direct == farey[:30]}; unresolved_at_200={unresolved}")

    caps, minimizers = exact_caps()
    expected_caps = {
        8: Q(2243, 5880), 9: Q(1979, 4004), 10: Q(55, 91),
        11: Q(66, 91), 12: Q(6, 7), 13: Q(1),
    }
    require(caps == expected_caps, "exhaustive cap census changed")
    thresholds = {}
    crossings = []
    for k in range(8, 14):
        if caps[k] == 1:
            thresholds[k] = None
            continue
        first_above = next(m for m, value in enumerate(farey, 1) if value > caps[k])
        thresholds[k] = first_above - 2
        crossings.append((k, farey[first_above - 2], caps[k], farey[first_above - 1]))
    require(tuple(thresholds[k] for k in range(8, 13)) == (7, 8, 10, 13, 26), "wrong finite N* values")
    require(thresholds[13] is None, "cap-one row should be unbounded")
    print(f"caps={caps}; minimizers_sha256={digest(minimizers)}")
    print(f"crossings={crossings}")
    print("Nstar=(7,8,10,13,26,infinity)")

    d12 = direct_component_deficit(12)
    f12 = owner_formula_deficit(12)
    d13 = direct_component_deficit(13)
    f13 = owner_formula_deficit(13)
    require(d12 != f12, "m=12 must remain the hostile tail-formula failure")
    require(len(noncover_components(12)) == 14, "m=12 component count")
    require(d13 == f13, "sharp m=13 base failed")
    guards13 = prove_local_guards(13)
    guards14 = prove_local_guards(14)
    min13 = min(record[3] for record in guards13)
    min14 = min(record[3] for record in guards14)
    endpoints = audit_m13_zero_margin_endpoints()
    require(min13 == 0 and min14 == Q(1, 13), "unexpected track guard margins")
    require(all(record[-3:] == (False, True, True) for record in endpoints), "m=13 endpoint convention changed")
    owner_components13 = audit_owner_component_bijection(13)
    print(f"m12_immediate_predecessor_hostile=direct:{d12},owner_formula:{f12},components:14")
    print(f"m13_sharp_base=direct:{d13},owner_formula:{f13},components:12")
    print(f"m13_owner_component_records={owner_components13}")
    print(f"track_guards=min_m13:{min13},min_m14:{min14}; endpoint_audit={endpoints}")

    formula_tail = tuple(owner_formula_deficit(m) for m in range(13, 201))
    farey_tail = tuple(1 - farey[m - 1] for m in range(13, 201))
    require(formula_tail == farey_tail, "owner formula and Farey engine disagree on m=13..200")
    print(f"owner_formula_m13_200_sha256={digest(formula_tail)}; farey_equal={formula_tail == farey_tail}")

    laws, constants = prove_phase_rational_law(period=60, minimum_n=12)
    require(len(laws) == 60, "missing residue phase")
    require(constants == {Q(127, 35)}, "phase constants disagree")
    all_terms = [term for _, terms, _ in laws.values() for term in terms]
    require(all(Q(0) <= C and 0 <= c <= 5 for C, c in all_terms), "invalid C/(n-c) term")
    require(all(sum((C for C, _ in terms), Q(0)) == Q(127, 5) for _, terms, _ in laws.values()), "phase coefficient sum changed")
    phase_records = []
    for phase in range(60):
        first_n, terms, _ = laws[phase]
        aggregate = {}
        for coefficient, shift in terms:
            aggregate[shift] = aggregate.get(shift, Q(0)) + coefficient
        aggregate = tuple(
            (coefficient, shift)
            for shift, coefficient in sorted(aggregate.items())
            if coefficient
        )
        require(first_n % 60 == phase, (phase, "first admissible n"))
        for n in (first_n, first_n + 60):
            predicted = sum(
                (coefficient / (n - shift) for coefficient, shift in aggregate),
                Q(0),
            ) / 7
            require(
                predicted == owner_formula_deficit(n + 1),
                (phase, n, "printed phase formula mismatch"),
            )
        phase_records.append((phase, first_n, aggregate))
    print(f"phase_laws=60; phase_law_sha256={digest(laws)}")
    for phase, first_n, aggregate in phase_records:
        print(f"phase={phase:02d};first_n={first_n};terms={aggregate}")
    print("tail_exact=phase-wise sum C/[7(n-c)], n=m-1, 0<=c<=5")
    print("tail_bounds=127/[35(m-1)] <= 1-a(m) <= 127/[35(m-6)] for m>=13")
    print("asymptotic=lim m(1-a(m))=127/35")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
