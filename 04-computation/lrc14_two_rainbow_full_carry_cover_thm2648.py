#!/usr/bin/env python3
"""Exact referee for THM-2648.

Constructs two rainbow matchings for every ordered pair of eleven-sheet
relations on C_13 and verifies that their carry-colour sets cover C_13.
All guards survive optimized Python.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inv(value, p):
    require(value % p != 0, "inverse of zero")
    return pow(value, p - 2, p)


def complement(pair, p):
    missing = set(pair)
    return tuple(x for x in range(p) if x not in missing)


def affine_matching(a_pair, b_pair, orientation, p):
    a0, a1 = a_pair
    b0, b1 = b_pair
    slope = (b1 - b0) * inv(a1 - a0, p) % p
    if orientation == 1:
        return tuple((x, (b0 + slope * (x - a0)) % p) for x in complement(a_pair, p))
    require(orientation == -1, "affine orientation")
    return tuple((x, (b1 - slope * (x - a0)) % p) for x in complement(a_pair, p))


def normalized_template(kind):
    if kind == "parallel":
        targets = (2, 3, 4, 5, 7, 8, 6, 9, 10, 11, 12)
    elif kind == "parallel_inverse":
        targets = (2, 3, 4, 5, 8, 6, 7, 9, 10, 11, 12)
    elif kind == "antiparallel":
        targets = (1, 2, 3, 4, 6, 7, 5, 8, 9, 10, 11)
    elif kind == "antiparallel_inverse":
        targets = (1, 2, 3, 4, 7, 5, 6, 8, 9, 10, 11)
    else:
        raise RuntimeError("unknown template kind")
    return tuple(zip(range(2, 13), targets))


def transported_template(a_pair, b_pair, kind, p):
    a0, a1 = a_pair
    b0, _ = b_pair
    scale = (a1 - a0) % p
    return tuple(
        ((a0 + scale * x) % p, (b0 + scale * y) % p)
        for x, y in normalized_template(kind)
    )


def verify_matching(a_pair, b_pair, matching, p):
    source = complement(a_pair, p)
    target = complement(b_pair, p)
    xs = tuple(x for x, _ in matching)
    ys = tuple(y for _, y in matching)
    colours = tuple((x + y) % p for x, y in matching)
    require(len(matching) == p - 2, "matching edge count")
    require(set(xs) == set(source) and len(set(xs)) == p - 2,
            "matching source bijection")
    require(set(ys) == set(target) and len(set(ys)) == p - 2,
            "matching target bijection")
    require(len(set(colours)) == p - 2, "rainbow colour injectivity")
    holes = frozenset(set(range(p)) - set(colours))
    require(len(holes) == 2, "rainbow hole pair")
    return frozenset(colours), holes


def is_affine(matching, p):
    (x0, y0), (x1, y1) = matching[:2]
    slope = (y1 - y0) * inv(x1 - x0, p) % p
    intercept = (y0 - slope * x0) % p
    return all(y == (slope * x + intercept) % p for x, y in matching)


def minimal_repair_census(a_pair, b_pair, baseline, p):
    """Count compatible repairs moving exactly two or three target values."""
    _, baseline_holes = verify_matching(a_pair, b_pair, baseline, p)
    targets = tuple(y for _, y in baseline)
    sources = tuple(x for x, _ in baseline)
    counts = Counter()

    for i, j in combinations(range(p - 2), 2):
        moved = list(targets)
        moved[i], moved[j] = moved[j], moved[i]
        colours = tuple((sources[r] + moved[r]) % p for r in range(p - 2))
        if len(set(colours)) == p - 2:
            holes = frozenset(set(range(p)) - set(colours))
            if holes.isdisjoint(baseline_holes):
                counts[2] += 1

    for i, j, k in combinations(range(p - 2), 3):
        for cycle in ((j, k, i), (k, i, j)):
            moved = list(targets)
            moved[i], moved[j], moved[k] = (
                targets[cycle[0]], targets[cycle[1]], targets[cycle[2]]
            )
            colours = tuple((sources[r] + moved[r]) % p for r in range(p - 2))
            if len(set(colours)) == p - 2:
                holes = frozenset(set(range(p)) - set(colours))
                if holes.isdisjoint(baseline_holes):
                    counts[3] += 1
    return counts


def main():
    p = 13
    pairs = tuple(combinations(range(p), 2))
    require(len(pairs) == 78, "two-subset count")

    parallel_template = normalized_template("parallel")
    parallel_inverse = normalized_template("parallel_inverse")
    anti_template = normalized_template("antiparallel")
    anti_inverse = normalized_template("antiparallel_inverse")
    c_par, h_par = verify_matching((0, 1), (0, 1), parallel_template, p)
    c_par_inv, h_par_inv = verify_matching((0, 1), (0, 1), parallel_inverse, p)
    c_anti, h_anti = verify_matching((0, 1), (0, 12), anti_template, p)
    c_anti_inv, h_anti_inv = verify_matching((0, 1), (0, 12), anti_inverse, p)
    require(h_par == h_par_inv == frozenset((3, 12)), "parallel template holes")
    require(h_anti == h_anti_inv == frozenset((2, 11)), "antiparallel template holes")
    require(not is_affine(parallel_template, p) and not is_affine(anti_template, p),
            "matched templates must be nonlinear")
    require(not is_affine(parallel_inverse, p) and not is_affine(anti_inverse, p),
            "inverse matched templates must be nonlinear")

    parallel_affine = tuple((x, x) for x in range(2, 13))
    anti_affine = tuple((x, x - 1) for x in range(2, 13))
    repair_parallel = minimal_repair_census((0, 1), (0, 1), parallel_affine, p)
    repair_anti = minimal_repair_census((0, 1), (0, 12), anti_affine, p)
    require(repair_parallel == Counter({3: 2}) and repair_anti == Counter({3: 2}),
            "sharp distance-three repair census")

    collision_1 = tuple(zip(range(2, 13), (2, 4, 11, 5, 6, 7, 8, 9, 3, 10, 12)))
    collision_2 = tuple(zip(range(2, 13), (2, 10, 3, 5, 6, 7, 8, 9, 11, 4, 12)))
    _, collision_h1 = verify_matching((0, 1), (0, 1), collision_1, p)
    _, collision_h2 = verify_matching((0, 1), (0, 1), collision_2, p)
    require(collision_1 != collision_2 and collision_h1 == collision_h2 == frozenset((6, 9)),
            "same-hole nonlinear occurrence collision")

    unrestricted_1 = tuple(zip(range(2, 13), (2, 3, 4, 6, 8, 5, 10, 7, 12, 9, 11)))
    unrestricted_2 = tuple(zip(range(2, 13), (2, 12, 4, 6, 8, 5, 10, 7, 3, 9, 11)))
    _, unrestricted_h1 = verify_matching((0, 1), (0, 1), unrestricted_1, p)
    _, unrestricted_h2 = verify_matching((0, 1), (0, 1), unrestricted_2, p)
    require(unrestricted_h1 == frozenset((0, 2)), "unrestricted first holes")
    require(unrestricted_h2 == frozenset((6, 9)), "unrestricted second holes")
    require(unrestricted_h1.isdisjoint(unrestricted_h2), "unrestricted holes disjoint")
    require(not is_affine(unrestricted_1, p) and not is_affine(unrestricted_2, p),
            "unrestricted optimum must be nonlinear")
    require(len(set(unrestricted_1) & set(unrestricted_2)) == 9,
            "unrestricted optimum overlap")
    require(len(set(unrestricted_1) | set(unrestricted_2)) == 13,
            "unrestricted optimum union")
    unrestricted_anti_1 = tuple((x, y - 1) for x, y in unrestricted_1)
    unrestricted_anti_2 = tuple((x, y - 1) for x, y in unrestricted_2)
    _, unrestricted_anti_h1 = verify_matching((0, 1), (0, 12), unrestricted_anti_1, p)
    _, unrestricted_anti_h2 = verify_matching((0, 1), (0, 12), unrestricted_anti_2, p)
    require(unrestricted_anti_h1 == frozenset((1, 12)),
            "unrestricted antiparallel first holes")
    require(unrestricted_anti_h2 == frozenset((5, 8)),
            "unrestricted antiparallel second holes")
    require(unrestricted_anti_h1.isdisjoint(unrestricted_anti_h2),
            "unrestricted antiparallel holes disjoint")
    require(len(set(unrestricted_anti_1) & set(unrestricted_anti_2)) == 9,
            "unrestricted antiparallel optimum overlap")
    # Distinct permutations cannot differ at exactly one source, so overlap
    # ten (and hence union twelve) is impossible.  The displayed overlap-nine
    # pair attains the unrestricted lower bound thirteen.

    class_census = Counter()
    chart_census = Counter()
    cover_multiplicity = Counter()
    edge_overlap_census = Counter()
    retained_edge_census = Counter()
    cover_character_checks = 0
    chart_swap_character_checks = 0
    three_chart_character_checks = 0
    degenerate_affine_controls = 0
    matching_edges = 0

    for a_pair in pairs:
        a0, a1 = a_pair
        d_a = (a1 - a0) % p
        for b_pair in pairs:
            b0, b1 = b_pair
            d_b = (b1 - b0) % p
            ratio = d_b * inv(d_a, p) % p

            plus = affine_matching(a_pair, b_pair, 1, p)
            minus = affine_matching(a_pair, b_pair, -1, p)
            plus_colours = frozenset((x + y) % p for x, y in plus)
            minus_colours = frozenset((x + y) % p for x, y in minus)

            if ratio not in (1, p - 1):
                class_census["step_distinct"] += 1
                c1, h1 = verify_matching(a_pair, b_pair, plus, p)
                c2, h2 = verify_matching(a_pair, b_pair, minus, p)
                matching1, matching2 = plus, minus
                require(h1 == frozenset(((a0 + b0) % p, (a1 + b1) % p)),
                        "plus affine hole formula")
                require(h2 == frozenset(((a0 + b1) % p, (a1 + b0) % p)),
                        "minus affine hole formula")
                chart_census["affine"] += 2
            elif ratio == 1:
                class_census["step_matched"] += 1
                c1, h1 = verify_matching(a_pair, b_pair, plus, p)
                require(len(plus_colours) == p - 2 and len(minus_colours) == 1,
                        "parallel affine degeneration")
                nonlinear = transported_template(a_pair, b_pair, "parallel", p)
                c2, h2 = verify_matching(a_pair, b_pair, nonlinear, p)
                nonlinear_inverse = transported_template(a_pair, b_pair, "parallel_inverse", p)
                c3, h3 = verify_matching(a_pair, b_pair, nonlinear_inverse, p)
                matching1, matching2 = plus, nonlinear
                base = (a0 + b0) % p
                expected = frozenset((base + d_a * h) % p for h in (3, 12))
                require(h2 == expected, "parallel nonlinear hole transport")
                require(h3 == expected and c3 == c2, "parallel inverse chart transport")
                require(not is_affine(nonlinear, p), "parallel nonlinear transport")
                chart_census["affine"] += 1
                chart_census["nonlinear"] += 1
                degenerate_affine_controls += 1
            else:
                class_census["step_matched"] += 1
                c1, h1 = verify_matching(a_pair, b_pair, minus, p)
                require(len(minus_colours) == p - 2 and len(plus_colours) == 1,
                        "antiparallel affine degeneration")
                nonlinear = transported_template(a_pair, b_pair, "antiparallel", p)
                c2, h2 = verify_matching(a_pair, b_pair, nonlinear, p)
                nonlinear_inverse = transported_template(a_pair, b_pair, "antiparallel_inverse", p)
                c3, h3 = verify_matching(a_pair, b_pair, nonlinear_inverse, p)
                matching1, matching2 = minus, nonlinear
                base = (a0 + b0) % p
                expected = frozenset((base + d_a * h) % p for h in (2, 11))
                require(h2 == expected, "antiparallel nonlinear hole transport")
                require(h3 == expected and c3 == c2, "antiparallel inverse chart transport")
                require(not is_affine(nonlinear, p), "antiparallel nonlinear transport")
                chart_census["affine"] += 1
                chart_census["nonlinear"] += 1
                degenerate_affine_controls += 1

            require(h1.isdisjoint(h2), "two rainbow hole pairs must be disjoint")
            require(c1 | c2 == frozenset(range(p)), "two rainbow charts must cover all carries")
            edge_overlap = len(set(matching1) & set(matching2))
            expected_overlap = 8 if ratio in (1, p - 1) else 1
            require(edge_overlap == expected_overlap, "chart edge-overlap class")
            edge_overlap_census[edge_overlap] += 1
            retained_edge_census[len(set(matching1) | set(matching2))] += 1
            multiplicities = Counter()
            for colour in range(p):
                multiplicity = int(colour in c1) + int(colour in c2)
                require(multiplicity in (1, 2), "cover multiplicity")
                multiplicities[multiplicity] += 1
                cover_multiplicity[multiplicity] += 1
            require(multiplicities == Counter({2: 9, 1: 4}),
                    "per-pair cover multiplicity profile")
            cover_vector = tuple(int(colour in c1) + int(colour in c2)
                                 for colour in range(p))
            signed_vector = tuple(int(colour in c1) - int(colour in c2)
                                  for colour in range(p))
            mean = Fraction(22, 13)
            energy = sum((Fraction(value) - mean) ** 2 for value in cover_vector) / p
            require(energy == Fraction(36, 169), "cover charged energy")
            even_c2_energy = energy / 4
            odd_c2_energy = sum(Fraction(value * value, 4 * p)
                                for value in signed_vector)
            require(even_c2_energy == Fraction(9, 169), "even C2 charged energy")
            require(odd_c2_energy == Fraction(1, 13), "chart-swap charged energy")
            for k in range(1, p):
                # Multiplication by nonzero k permutes the coefficient vector.
                # In Q(zeta_13), a degree-<13 rational vector vanishes only
                # when all thirteen coefficients are equal.
                transformed = tuple(cover_vector[(-k * exponent) % p]
                                    for exponent in range(p))
                require(len(set(transformed)) > 1, "cover character vanished")
                signed_transformed = tuple(signed_vector[(-k * exponent) % p]
                                           for exponent in range(p))
                require(len(set(signed_transformed)) > 1,
                        "chart-swap character vanished")
                cover_character_checks += 1
                chart_swap_character_checks += 1

            if ratio in (1, p - 1):
                three_profile = tuple(
                    int(colour in c1) + int(colour in c2) + int(colour in c3)
                    for colour in range(p)
                )
                require(Counter(three_profile) == Counter({3: 9, 2: 2, 1: 2}),
                        "reflection-stable three-chart profile")
                three_mean = Fraction(33, 13)
                three_energy = sum((Fraction(value) - three_mean) ** 2
                                   for value in three_profile) / p
                require(three_energy == Fraction(94, 169),
                        "three-chart centered energy")
                c3_difference = tuple(int(colour in c1) - int(colour in c2)
                                      for colour in range(p))
                c3_energy = sum(Fraction(value * value, 9 * p)
                                for value in c3_difference)
                require(c3_energy == Fraction(4, 117),
                        "nontrivial C3 chart-character energy")
                for _character in (1, 2):
                    for k in range(1, p):
                        transformed = tuple(c3_difference[(-k * exponent) % p]
                                            for exponent in range(p))
                        require(len(set(transformed)) > 1,
                                "nontrivial C3 carry character vanished")
                        three_chart_character_checks += 1
            matching_edges += 2 * (p - 2)

    require(class_census == Counter({"step_distinct": 5070, "step_matched": 1014}),
            "pair class census")
    require(chart_census == Counter({"affine": 11154, "nonlinear": 1014}),
            "chart type census")
    require(degenerate_affine_controls == 1014, "degenerate affine census")
    require(cover_multiplicity == Counter({2: 54756, 1: 24336}),
            "global cover multiplicity census")
    require(edge_overlap_census == Counter({1: 5070, 8: 1014}),
            "chart edge-overlap census")
    require(retained_edge_census == Counter({21: 5070, 14: 1014}),
            "retained edge census")
    require(matching_edges == 133848, "matching edge census")
    require(cover_character_checks == 73008, "cover character census")
    require(chart_swap_character_checks == 73008, "chart-swap character census")
    require(three_chart_character_checks == 24336, "three-chart character census")

    print("THM2648 TWO-RAINBOW FULL-CARRY COVER EXACT REFEREE")
    print(f"relation_pairs={len(pairs) ** 2} matching_charts={2 * len(pairs) ** 2} matching_edges={matching_edges}")
    print("pair_classes={step_distinct:5070,step_matched:1014}")
    print("atlas_chart_types={affine_all:11154,nonlinear_chosen:1014} degenerate_affine_controls=1014")
    print("normal_form_holes={parallel_nonlinear:(3,12),antiparallel_nonlinear:(2,11)}")
    print("affine_anchored_repair={distance2:0,distance3:2_each_normal_form,union:14}")
    print("unrestricted_optimum={overlap:9,union:13,parallel:(0,2)|(6,9),antiparallel:(1,12)|(5,8)}")
    print("same_hole_nonlinear_collision=parallel_holes_(6,9)")
    print("every_pair_holes=disjoint every_pair_colour_union=13")
    print("cover_multiplicity_census={one:24336,two:54756}")
    print("cover_character_checks=73008_all_nonzero centered_energy=36/169")
    print("chart_swap_checks=73008_all_nonzero C2_energies={even:9/169,odd:1/13}")
    print("matched_three_chart={profile:1^2_2^2_3^9,energy:94/169,C3_each:4/117,checks:24336}")
    print("edge_overlap_census={one:5070,eight:1014} retained_edges={21:5070,14:1014}")
    print("matched_class=one_affine_plus_one_nonlinear")
    print("PASS")


if __name__ == "__main__":
    main()
