#!/usr/bin/env python3
"""Exact referee for THM-2645.

The cyclotomic checks use the integral group-ring model
Z[X]/(1+X+...+X^12), so no floating-point tolerance is involved.  Every
logical check uses ``require`` and therefore remains active under ``-O``.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_prime(p):
    if p < 2:
        return False
    return all(p % d for d in range(2, int(p ** 0.5) + 1))


def convolution_profile(left, right, p):
    profile = [0] * p
    for x in left:
        for y in right:
            profile[(x + y) % p] += 1
    return tuple(profile)


def step_class(pair, p):
    a, b = pair
    d = (b - a) % p
    return min(d, (-d) % p)


def character_numerator(a_pair, b_pair, k, p):
    """Coefficient vector for sum_{a,b} zeta^{-k(a+b)}."""
    coeffs = [0] * p
    for a in a_pair:
        for b in b_pair:
            coeffs[(-k * (a + b)) % p] += 1
    return tuple(coeffs)


def cyclotomic_zero(coeffs):
    """Degree-<p vector vanishes at primitive pth root iff constant."""
    return len(set(coeffs)) == 1


def translate_pair(pair, shift, p):
    return tuple(sorted((x + shift) % p for x in pair))


def main():
    p = 13
    require(is_prime(p), "p=13 must be prime")
    pairs = tuple(combinations(range(p), 2))
    require(len(pairs) == 78, "two-subset count")

    class_census = Counter()
    energy_census = Counter()
    raw_energy_census = Counter()
    character_checks = 0
    gauge_checks = 0
    return_gate_checks = 0
    max_oriented_return = 0
    min_return_deficit = None

    for a_pair in pairs:
        a_set = set(a_pair)
        left = tuple(x for x in range(p) if x not in a_set)
        for b_pair in pairs:
            b_set = set(b_pair)
            right = tuple(x for x in range(p) if x not in b_set)
            r_profile = convolution_profile(left, right, p)
            m_profile = convolution_profile(a_pair, b_pair, p)

            require(all(r_profile[c] == 9 + m_profile[c] for c in range(p)),
                    "r=9+m factorization")
            require(sum(r_profile) == 121 and sum(m_profile) == 4,
                    "profile masses")

            matched = step_class(a_pair, p) == step_class(b_pair, p)
            class_census["step_matched" if matched else "step_distinct"] += 1

            if matched:
                require(Counter(m_profile) == Counter({0: 10, 1: 2, 2: 1}),
                        "matched multiplicity profile")
                expected_energy = Fraction(62, 169)
            else:
                require(Counter(m_profile) == Counter({0: 9, 1: 4}),
                        "distinct multiplicity profile")
                expected_energy = Fraction(36, 169)

            mean = Fraction(121, 13)
            centered_energy = sum((Fraction(v) - mean) ** 2 for v in r_profile) / p
            require(centered_energy == expected_energy, "centered energy")

            # THM-2644 cannot fire on this dense multiplicity profile.
            mass = sum(r_profile)
            raw_energy = sum(v * v for v in r_profile)
            defect = mass * mass - raw_energy
            oriented_return = sum(r_profile[c] * r_profile[(-c) % p] for c in range(p))
            require(mass == 121, "multiplicity mass")
            require(raw_energy == (1131 if matched else 1129), "raw energy class")
            require(oriented_return <= raw_energy < defect,
                    "purity/return gate unexpectedly fired")
            raw_energy_census[raw_energy] += 1
            return_gate_checks += 1
            max_oriented_return = max(max_oriented_return, oriented_return)
            gap = defect - oriented_return
            min_return_deficit = gap if min_return_deficit is None else min(min_return_deficit, gap)

            # Exact unnormalized Parseval for m; for k != 0,
            # rhat(k)=P_k/13 where P_k is this two-by-two root sum.
            numerator_energy = p * sum(v * v for v in m_profile) - 16
            require(Fraction(numerator_energy, p * p) == expected_energy,
                    "group-ring Parseval energy")
            energy_census[str(expected_energy)] += 1

            for k in range(1, p):
                numerator = character_numerator(a_pair, b_pair, k, p)
                require(sum(numerator) == 4, "character numerator mass")
                require(not cyclotomic_zero(numerator),
                        "charged character unexpectedly zero")
                character_checks += 1

            # The common-origin gauge (A+t,B-t) preserves the full profile.
            orbit = set()
            for shift in range(p):
                moved_a = translate_pair(a_pair, shift, p)
                moved_b = translate_pair(b_pair, -shift, p)
                orbit.add((moved_a, moved_b))
                require(convolution_profile(moved_a, moved_b, p) == m_profile,
                        "common-origin gauge profile")
                gauge_checks += 1
            require(len(orbit) == p, "common-origin gauge must be free")

    require(class_census == Counter({"step_distinct": 5070, "step_matched": 1014}),
            "difference-class census")
    require(energy_census == Counter({"36/169": 5070, "62/169": 1014}),
            "energy census")
    require(raw_energy_census == Counter({1129: 5070, 1131: 1014}),
            "raw energy census")
    require(return_gate_checks == 6084 and max_oriented_return == 1131,
            "purity/return no-go census")
    require(min_return_deficit == 12379, "purity/return minimum deficit")
    require(Fraction(36, 169) / 12 == Fraction(3, 169),
            "max-mode Cauchy floor")
    require(character_checks == 6084 * 12, "character check count")
    require(gauge_checks == 6084 * 13, "gauge check count")

    print("THM2645 eleven-sheet multiplicity spectrum exact referee")
    print(f"prime={p} two_subsets={len(pairs)} ordered_relation_pairs={len(pairs) ** 2}")
    print(f"charged_character_checks={character_checks} all_nonzero=PASS")
    print("difference_class_census={step_distinct:5070,step_matched:1014}")
    print("centered_energy_census={36/169:5070,62/169:1014}")
    print("max_mode_square_floor=3/169")
    print("purity_return_no_go=M121_Eraw1129or1131_Rmax1131_delta_min13510_min_deficit12379")
    print(f"common_origin_gauge_checks={gauge_checks} orbit_size=13")
    print("normalization=hat_r(k)=P_k/13_for_k_nonzero")
    print("PASS")


if __name__ == "__main__":
    main()
