#!/usr/bin/env python3
"""Deterministic exact companion for THM-2603.

The finite part is carried out directly in PSL_2(F_13), represented by
determinant-one matrices modulo the central pair {I,-I}.  The final section
checks the integral Barning/Farey hostile without importing any LRC data.
"""

from collections import Counter, deque
from itertools import combinations, product
from math import gcd


P = 13
INF = P
IDENTITY = (1, 0, 0, 1)
MINUS_IDENTITY = (P - 1, 0, 0, P - 1)
CHECKS = 0


def check(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def matrix_det(a: tuple[int, int, int, int]) -> int:
    x, y, z, w = a
    return (x * w - y * z) % P


def matrix_trace(a: tuple[int, int, int, int]) -> int:
    return (a[0] + a[3]) % P


def matrix_mul(
    a: tuple[int, int, int, int], b: tuple[int, int, int, int]
) -> tuple[int, int, int, int]:
    x, y, z, w = a
    r, s, t, u = b
    return (
        (x * r + y * t) % P,
        (x * s + y * u) % P,
        (z * r + w * t) % P,
        (z * s + w * u) % P,
    )


def matrix_inverse(a: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    x, y, z, w = a
    return (w % P, -y % P, -z % P, x % P)


def matrix_power(
    a: tuple[int, int, int, int], exponent: int
) -> tuple[int, int, int, int]:
    out = IDENTITY
    base = a
    while exponent:
        if exponent & 1:
            out = matrix_mul(out, base)
        base = matrix_mul(base, base)
        exponent //= 2
    return out


def psl_canonical(a: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    reduced = tuple(value % P for value in a)
    negative = tuple(-value % P for value in reduced)
    return min(reduced, negative)


def psl_equal(
    a: tuple[int, int, int, int], b: tuple[int, int, int, int]
) -> bool:
    return psl_canonical(a) == psl_canonical(b)


def psl_order(a: tuple[int, int, int, int]) -> int:
    current = psl_canonical(IDENTITY)
    for order in range(1, 200):
        current = psl_canonical(matrix_mul(current, a))
        if current == psl_canonical(IDENTITY):
            return order
    raise AssertionError("projective order exceeded the finite search bound")


def psl_trace_square(a: tuple[int, int, int, int]) -> int:
    """The sign-independent trace invariant of a projective SL_2 class."""

    return ((a[0] + a[3]) % P) ** 2 % P


def generated_group(
    generators: tuple[tuple[int, int, int, int], ...]
) -> set[tuple[int, int, int, int]]:
    identity = psl_canonical(IDENTITY)
    seen = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in generators:
            child = psl_canonical(matrix_mul(current, generator))
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return seen


def projective_action(a: tuple[int, int, int, int], point: int) -> int:
    x, y, z, w = a
    if point == INF:
        numerator, denominator = x, z
    else:
        numerator = (x * point + y) % P
        denominator = (z * point + w) % P
    if denominator == 0:
        return INF
    return numerator * pow(denominator, -1, P) % P


def projective_cycles(a: tuple[int, int, int, int]) -> tuple[tuple[int, ...], ...]:
    seen: set[int] = set()
    cycles = []
    for start in range(P + 1):
        if start in seen:
            continue
        cycle = []
        current = start
        while current not in seen:
            seen.add(current)
            cycle.append(current)
            current = projective_action(a, current)
        cycles.append(tuple(cycle))
    return tuple(cycles)


def cyclic_subgroup(a: tuple[int, int, int, int]) -> set[tuple[int, int, int, int]]:
    return {psl_canonical(matrix_power(a, exponent)) for exponent in range(psl_order(a))}


def conjugacy_class(
    representative: tuple[int, int, int, int],
    group: set[tuple[int, int, int, int]],
) -> set[tuple[int, int, int, int]]:
    return {
        psl_canonical(
            matrix_mul(matrix_mul(element, representative), matrix_inverse(element))
        )
        for element in group
    }


def ordered_conjugate_norm(
    generator: tuple[int, int, int, int],
    root: tuple[int, int, int, int],
    exponent: int,
    forward: bool,
) -> tuple[int, int, int, int]:
    factors = []
    for index in range(7):
        generator_power = matrix_power(generator, index)
        factors.append(
            matrix_mul(
                matrix_mul(generator_power, matrix_power(root, exponent % P)),
                matrix_inverse(generator_power),
            )
        )
    out = IDENTITY
    for factor in factors if forward else reversed(factors):
        out = matrix_mul(out, factor)
    return psl_canonical(out)


def fixed_points(a: tuple[int, int, int, int]) -> tuple[int, ...]:
    return tuple(point for point in range(P + 1) if projective_action(a, point) == point)


def left_coset(
    representative: tuple[int, int, int, int],
    subgroup: set[tuple[int, int, int, int]],
) -> frozenset[tuple[int, int, int, int]]:
    return frozenset(
        psl_canonical(matrix_mul(representative, element)) for element in subgroup
    )


def left_multiply_coset(
    multiplier: tuple[int, int, int, int],
    coset: frozenset[tuple[int, int, int, int]],
) -> frozenset[tuple[int, int, int, int]]:
    return frozenset(
        psl_canonical(matrix_mul(multiplier, representative))
        for representative in coset
    )


def right_multiply_coset(
    coset: frozenset[tuple[int, int, int, int]],
    multiplier: tuple[int, int, int, int],
) -> frozenset[tuple[int, int, int, int]]:
    return frozenset(
        psl_canonical(matrix_mul(representative, multiplier))
        for representative in coset
    )


def permutation_orbits(universe: set, step) -> tuple[tuple, ...]:
    """Exact orbits of a permutation, with a deterministic traversal order."""

    seen = set()
    orbits = []
    for start in sorted(universe, key=repr):
        if start in seen:
            continue
        orbit = []
        current = start
        while current not in seen:
            if current not in universe:
                raise AssertionError("permutation left its declared universe")
            seen.add(current)
            orbit.append(current)
            current = step(current)
        if current != start:
            raise AssertionError("permutation traversal merged distinct cycles")
        orbits.append(tuple(orbit))
    return tuple(orbits)


def integer_matrix_vector(
    a: tuple[tuple[int, int], tuple[int, int]], d: tuple[int, int]
) -> tuple[int, int]:
    return (a[0][0] * d[0] + a[0][1] * d[1], a[1][0] * d[0] + a[1][1] * d[1])


def farey_defect(d: tuple[int, int]) -> int:
    # THM-2056's displayed calibration w=(1,0): F(d)=||d||^2-91 d_1.
    return d[0] * d[0] + d[1] * d[1] - 91 * d[0]


def main() -> None:
    # Integer lifts have determinant one and exact projective torsion 2 and 3.
    x = (0, -1 % P, 1, 0)
    y = (3, -7 % P, 1, -2 % P)
    check(matrix_det(x) == 1, "x is not in SL_2(F_13)")
    check(matrix_det(y) == 1, "y is not in SL_2(F_13)")
    check(matrix_power(x, 2) == MINUS_IDENTITY, "x^2 is not -I")
    check(matrix_power(y, 3) == MINUS_IDENTITY, "y^3 is not -I")
    check(psl_order(x) == 2, "x does not have projective order two")
    check(psl_order(y) == 3, "y does not have projective order three")

    a = matrix_mul(x, y)
    c = matrix_mul(matrix_mul(matrix_mul(x, y), matrix_inverse(x)), matrix_inverse(y))
    d = matrix_mul(matrix_inverse(a), c)
    check(a == (12, 2, 3, 6), "wrong owner matrix A=xy")
    check(c == (5, 9, 9, 6), "wrong root matrix C=[x,y]")
    check(d == (12, 3, 2, 6), "wrong companion matrix A^-1 C")
    check(psl_order(a) == 7, "A does not have order seven")
    check(psl_order(c) == 13, "C does not have order thirteen")
    check(psl_order(d) == 7, "A^-1 C does not have order seven")

    # Enumerate PSL_2(F_13) independently of the chosen generators.
    full_group = {
        psl_canonical((r, s, t, u))
        for r, s, t, u in product(range(P), repeat=4)
        if (r * u - s * t) % P == 1
    }
    generated_xy = generated_group((x, y))
    generated_ac = generated_group((a, c))
    check(len(full_group) == 1092, "wrong order for PSL_2(F_13)")
    check(generated_xy == full_group, "x,y do not generate PSL_2(F_13)")
    check(generated_ac == full_group, "A,C do not generate PSL_2(F_13)")

    a_cycles = projective_cycles(a)
    c_cycles = projective_cycles(c)
    d_cycles = projective_cycles(d)
    check(tuple(map(len, a_cycles)) == (7, 7), "A is not two seven-cycles")
    check(tuple(map(len, c_cycles)) == (13, 1), "C is not 13+1")
    check(tuple(map(len, d_cycles)) == (7, 7), "A^-1 C is not two seven-cycles")
    check(fixed_points(c) == (5,), "C has the wrong projective fixed point")

    # The multiplication order is load-bearing.  This is
    # (A^6 C A^-6)...(A C A^-1)C, not the reverse product.
    descending_norm = IDENTITY
    ascending_norm = IDENTITY
    conjugate_fixed_points = []
    for exponent in range(7):
        a_power = matrix_power(a, exponent)
        conjugate = matrix_mul(matrix_mul(a_power, c), matrix_inverse(a_power))
        descending_norm = matrix_mul(conjugate, descending_norm)
        ascending_norm = matrix_mul(ascending_norm, conjugate)
        points = fixed_points(conjugate)
        check(len(points) == 1, "conjugate root subgroup lost its cusp")
        conjugate_fixed_points.append(points[0])
    check(psl_equal(descending_norm, IDENTITY), "descending norm did not close")
    check(psl_order(ascending_norm) == 13, "reverse norm is not the order-13 hostile")
    check(not psl_equal(matrix_power(c, 7), IDENTITY), "fixed-chart C^7 collapsed")
    check(
        tuple(conjugate_fixed_points) == (5, 11, 13, 4, 10, 7, 8),
        "the moving cusps do not form the expected A-orbit",
    )

    # Conjugate the root cycle to the standard affine translation q -> q+1.
    affine_chart = (1, -1 % P, 3, -2 % P)
    translation = (1, 1, 0, 1)
    affine_c = matrix_mul(matrix_mul(affine_chart, c), matrix_inverse(affine_chart))
    affine_a = matrix_mul(matrix_mul(affine_chart, a), matrix_inverse(affine_chart))
    check(psl_equal(affine_c, translation), "C did not become q -> q+1")
    check(affine_a == (7, 5, 10, 11), "wrong owner action in the affine root chart")
    affine_a_cycles = projective_cycles(affine_a)
    check(
        affine_a_cycles == ((0, 4, 6, 10, 7, 5, 3), (1, 8, 13, 2, 9, 12, 11)),
        "wrong projective owner orbits in the affine chart",
    )
    check(projective_action(affine_a, 8) == INF, "owner action did not leave the affine patch")
    check(projective_action(affine_a, INF) == 2, "owner action did not re-enter the affine patch")

    order_census = Counter(psl_order(element) for element in full_group)
    expected_census = Counter({1: 1, 2: 91, 3: 182, 6: 182, 7: 468, 13: 168})
    check(order_census == expected_census, "wrong PSL_2(F_13) order census")
    check(91 not in order_census, "an element of order 91 appeared")

    # Full ordered norm atlas.  The three g_t represent all three projective
    # conjugacy classes of elements of order seven (trace is projective only
    # through its square).  For every nonzero root exponent, exhaust both
    # orientations of the seven transported parabolic factors.
    trace_generators = {trace: (0, 1, -1 % P, trace) for trace in (3, 5, 6)}
    expected_trace_squares = {3: 9, 5: 12, 6: 10}
    order_seven_elements = {
        element for element in full_group if psl_order(element) == 7
    }
    trace_classes = {}
    for trace, generator in trace_generators.items():
        check(matrix_det(generator) == 1, f"g_{trace} is not determinant one")
        check(psl_order(generator) == 7, f"g_{trace} does not have order seven")
        check(
            psl_trace_square(generator) == expected_trace_squares[trace],
            f"g_{trace} has the wrong projective trace square",
        )
        trace_classes[trace] = conjugacy_class(generator, full_group)
        check(len(trace_classes[trace]) == 156, f"g_{trace} has wrong class size")
    check(
        all(
            trace_classes[left].isdisjoint(trace_classes[right])
            for left, right in combinations((3, 5, 6), 2)
        ),
        "the three order-seven trace classes overlap",
    )
    check(
        set().union(*trace_classes.values()) == order_seven_elements,
        "the trace representatives do not exhaust the order-seven classes",
    )

    expected_forward_closures = {
        3: {6, 8, 9, 10, 11},
        5: {2, 8, 10, 11, 12},
        6: {1, 3, 9, 11, 12},
    }
    expected_reverse_closures = {
        3: {2, 3, 4, 5, 7},
        5: {1, 2, 3, 5, 11},
        6: {1, 2, 4, 10, 12},
    }
    expected_norm_signature = Counter(
        {(1, 4): 5, (13, 4): 2, (3, 1): 2, (6, 3): 2, (2, 0): 1}
    )
    ordered_products = {}
    closure_state_sets = {}
    order_seven_traces = {3, -3 % P, 5, -5 % P, 6, -6 % P}
    for trace, generator in trace_generators.items():
        for orientation, forward in (("F", True), ("R", False)):
            products_by_exponent = {
                exponent: ordered_conjugate_norm(
                    generator, translation, exponent, forward
                )
                for exponent in range(1, 13)
            }
            ordered_products[trace, orientation] = products_by_exponent
            closure_set = {
                exponent
                for exponent, value in products_by_exponent.items()
                if psl_equal(value, IDENTITY)
            }
            closure_state_sets[f"{trace}{orientation}"] = closure_set
            expected = (
                expected_forward_closures[trace]
                if forward
                else expected_reverse_closures[trace]
            )
            check(closure_set == expected, f"wrong {trace}{orientation} closure set")
            signature_census = Counter(
                (psl_order(value), psl_trace_square(value))
                for value in products_by_exponent.values()
            )
            check(
                signature_census == expected_norm_signature,
                f"wrong {trace}{orientation} conjugacy-signature census",
            )
        check(
            all(
                psl_equal(
                    ordered_products[trace, "F"][exponent],
                    matrix_mul(
                        matrix_power(
                            matrix_mul(matrix_power(translation, exponent), generator),
                            7,
                        ),
                        matrix_power(matrix_inverse(generator), 7),
                    ),
                )
                for exponent in range(1, 13)
            ),
            f"the trace-{trace} forward telescope failed",
        )
        check(
            all(
                psl_equal(
                    ordered_products[trace, "R"][exponent],
                    matrix_mul(
                        matrix_power(generator, 7),
                        matrix_power(
                            matrix_mul(
                                matrix_inverse(generator),
                                matrix_power(translation, exponent),
                            ),
                            7,
                        ),
                    ),
                )
                for exponent in range(1, 13)
            ),
            f"the trace-{trace} reverse telescope failed",
        )
        check(
            all(
                matrix_trace(
                    matrix_mul(matrix_power(translation, exponent), generator)
                )
                == (trace - exponent) % P
                and matrix_trace(
                    matrix_mul(matrix_inverse(generator), matrix_power(translation, exponent))
                )
                == (trace + exponent) % P
                for exponent in range(1, 13)
            ),
            f"the trace-{trace} moving-trace formulas failed",
        )
        check(
            all(
                (
                    psl_equal(ordered_products[trace, "F"][exponent], IDENTITY)
                    == ((trace - exponent) % P in order_seven_traces)
                )
                and (
                    psl_equal(ordered_products[trace, "R"][exponent], IDENTITY)
                    == ((trace + exponent) % P in order_seven_traces)
                )
                for exponent in range(1, 13)
            ),
            f"the trace-{trace} closure criterion failed",
        )
        for exponent in range(1, 13):
            reverse_value = ordered_products[trace, "R"][exponent]
            negated_forward = ordered_products[trace, "F"][(-exponent) % P]
            check(
                psl_equal(reverse_value, matrix_inverse(negated_forward)),
                f"orientation inversion failed for trace {trace}, exponent {exponent}",
            )
        check(
            expected_reverse_closures[trace]
            == {(-exponent) % P for exponent in expected_forward_closures[trace]},
            f"reverse closure is not the negative of forward closure for trace {trace}",
        )

    all_nonzero_exponents = set(range(1, 13))
    check(
        set().union(*closure_state_sets.values()) == all_nonzero_exponents,
        "the full trace/orientation atlas misses a nonzero exponent",
    )
    universal_state_covers = []
    for cover_size in range(1, len(closure_state_sets) + 1):
        universal_state_covers = [
            states
            for states in combinations(sorted(closure_state_sets), cover_size)
            if set().union(*(closure_state_sets[state] for state in states))
            == all_nonzero_exponents
        ]
        if universal_state_covers:
            break
    check(cover_size == 3, "the ordered-norm selector cover number is not three")
    check(
        universal_state_covers == [("3F", "3R", "6F"), ("3F", "3R", "6R")],
        "the minimal ordered-norm selector covers are wrong",
    )
    check(
        closure_state_sets["3F"].isdisjoint(closure_state_sets["3R"])
        and closure_state_sets["3F"] | closure_state_sets["3R"]
        == set(range(2, 12)),
        "the two trace-three orientations do not partition exponents 2..11",
    )

    # Exact conjugacy classification of every nonidentity norm output.
    order_sets = {
        order: {element for element in full_group if psl_order(element) == order}
        for order in (2, 3, 6, 13)
    }
    product_values = [
        value
        for products_by_exponent in ordered_products.values()
        for value in products_by_exponent.values()
    ]
    for order in (2, 3, 6):
        representative = next(value for value in product_values if psl_order(value) == order)
        check(
            conjugacy_class(representative, full_group) == order_sets[order],
            f"the order-{order} norm outputs do not lie in the unique class",
        )
    square_unipotent_class = conjugacy_class(translation, full_group)
    nonsquare_unipotent_class = conjugacy_class(matrix_power(translation, 2), full_group)
    check(
        len(square_unipotent_class) == len(nonsquare_unipotent_class) == 84
        and square_unipotent_class.isdisjoint(nonsquare_unipotent_class)
        and square_unipotent_class | nonsquare_unipotent_class == order_sets[13],
        "the two order-thirteen conjugacy classes are wrong",
    )
    check(
        all(
            value in nonsquare_unipotent_class
            for value in product_values
            if psl_order(value) == 13
        ),
        "an ordered norm hit the square unipotent conjugacy class",
    )

    # The original affine pair belongs to the trace-five row, but the
    # simultaneous conjugation rescales U.  This records the exact bridge,
    # not a bare trace analogy.
    trace_five_bridge = (4, 5, 0, 10)
    check(matrix_det(trace_five_bridge) == 1, "the trace-five bridge is singular")
    check(
        psl_equal(
            matrix_mul(
                matrix_mul(trace_five_bridge, affine_a),
                matrix_inverse(trace_five_bridge),
            ),
            trace_generators[5],
        ),
        "the owner generator did not enter the trace-five normal form",
    )
    check(
        psl_equal(
            matrix_mul(
                matrix_mul(trace_five_bridge, translation),
                matrix_inverse(trace_five_bridge),
            ),
            matrix_power(translation, 3),
        ),
        "the trace-five bridge has the wrong root rescaling",
    )

    centralizer_a = {
        element
        for element in full_group
        if psl_equal(matrix_mul(element, a), matrix_mul(a, element))
    }
    centralizer_c = {
        element
        for element in full_group
        if psl_equal(matrix_mul(element, c), matrix_mul(c, element))
    }
    subgroup_a = cyclic_subgroup(a)
    subgroup_c = cyclic_subgroup(c)
    affine_root_subgroup = {
        psl_canonical(
            matrix_mul(matrix_mul(affine_chart, element), matrix_inverse(affine_chart))
        )
        for element in subgroup_c
    }
    check(
        affine_root_subgroup == cyclic_subgroup(translation),
        "the homogeneous root subgroup did not become the affine U",
    )
    check(centralizer_a == subgroup_a, "the C7 centralizer is larger than C7")
    check(centralizer_c == subgroup_c, "the C13 centralizer is larger than C13")
    check(
        centralizer_a & centralizer_c == {psl_canonical(IDENTITY)},
        "the C7 and C13 centralizers meet nontrivially",
    )
    commuting_nontrivial_pairs = sum(
        psl_equal(matrix_mul(left, right), matrix_mul(right, left))
        for left in subgroup_a - {psl_canonical(IDENTITY)}
        for right in subgroup_c - {psl_canonical(IDENTITY)}
    )
    check(commuting_nontrivial_pairs == 0, "a nontrivial C7/C13 pair commuted")

    # The 84-point homogeneous carrier Omega=G/<C>.  The normalizer of the
    # root subgroup is exactly the stabilizer of its unique cusp.  Its
    # quotient by <C> is a six-frame group acting on the right of Omega.
    power_to_exponent = {
        psl_canonical(matrix_power(c, exponent)): exponent for exponent in range(13)
    }
    check(len(power_to_exponent) == 13, "the C13 power table collapsed")
    normalizer_c = {
        element
        for element in full_group
        if psl_canonical(
            matrix_mul(matrix_mul(element, c), matrix_inverse(element))
        )
        in subgroup_c
    }
    cusp_stabilizer = {
        element for element in full_group if projective_action(element, 5) == 5
    }
    check(normalizer_c == cusp_stabilizer, "N_G(C13) is not the cusp stabilizer")
    check(len(normalizer_c) == 78, "the C13 normalizer does not have order 78")
    conjugation_exponents = Counter(
        power_to_exponent[
            psl_canonical(
                matrix_mul(matrix_mul(element, c), matrix_inverse(element))
            )
        ]
        for element in normalizer_c
    )
    square_exponents = (1, 3, 4, 9, 10, 12)
    check(
        conjugation_exponents == Counter({exponent: 13 for exponent in square_exponents}),
        "the normalizer does not act through the six nonzero squares",
    )

    frame = (0, 3, 4, 4)
    check(matrix_det(frame) == 1, "the frame generator is not determinant one")
    check(frame in normalizer_c, "the frame generator does not normalize C13")
    check(psl_order(frame) == 6, "the frame generator does not have order six")
    check(projective_action(frame, 5) == 5, "the frame generator moved the cusp")
    check(
        psl_equal(
            matrix_mul(matrix_mul(frame, c), matrix_inverse(frame)),
            matrix_power(c, 4),
        ),
        "the frame generator does not conjugate C to C^4",
    )
    frame_subgroup = cyclic_subgroup(frame)
    check(
        frame_subgroup & subgroup_c == {psl_canonical(IDENTITY)},
        "the C6 frame subgroup intersects C13 nontrivially",
    )
    check(
        generated_group((c, frame)) == normalizer_c,
        "C13 and the frame generator do not generate the normalizer",
    )

    omega = {left_coset(element, subgroup_c) for element in full_group}
    check(len(omega) == 84, "G/C13 does not have 84 cosets")
    check({len(coset) for coset in omega} == {13}, "a C13 coset has wrong size")
    normalizer_quotient = {left_coset(element, subgroup_c) for element in normalizer_c}
    frame_quotient = {
        left_coset(matrix_power(frame, exponent), subgroup_c) for exponent in range(6)
    }
    check(len(normalizer_quotient) == 6, "N_G(C13)/C13 does not have order six")
    check(frame_quotient == normalizer_quotient, "the chosen frame misses a quotient class")

    coset_base = {}
    for coset in omega:
        bases = {projective_action(representative, 5) for representative in coset}
        check(len(bases) == 1, "G/C13 -> P1(F13) is not well-defined")
        coset_base[coset] = next(iter(bases))
    fibre_sizes = Counter(coset_base.values())
    check(
        fibre_sizes == Counter({point: 6 for point in range(P + 1)}),
        "the projective fibres are not uniformly six-sheeted",
    )

    left_a = lambda coset: left_multiply_coset(a, coset)
    left_c = lambda coset: left_multiply_coset(c, coset)
    right_frame = lambda coset: right_multiply_coset(coset, frame)
    check(all(left_a(coset) in omega for coset in omega), "A did not act on G/C13")
    check(all(left_c(coset) in omega for coset in omega), "C did not act on G/C13")
    check(
        all(right_frame(coset) in omega for coset in omega),
        "the normalizer frame did not act on the right",
    )
    check(
        all(coset_base[right_frame(coset)] == coset_base[coset] for coset in omega),
        "the right frame action moved a projective base point",
    )
    check(
        all(
            left_a(right_frame(coset)) == right_frame(left_a(coset))
            for coset in omega
        ),
        "the left C7 and right C6 actions do not commute",
    )

    a_omega_orbits = permutation_orbits(omega, left_a)
    c_omega_orbits = permutation_orbits(omega, left_c)
    frame_omega_orbits = permutation_orbits(omega, right_frame)
    check(
        Counter(map(len, a_omega_orbits)) == Counter({7: 12}),
        "A does not split G/C13 into twelve seven-cycles",
    )
    check(
        Counter(map(len, c_omega_orbits)) == Counter({1: 6, 13: 6}),
        "C does not split G/C13 into six fixed frames and six 13-cycles",
    )
    check(
        Counter(map(len, frame_omega_orbits)) == Counter({6: 14}),
        "the right frame action is not fourteen six-cycles",
    )
    boundary_fibre = {coset for coset in omega if coset_base[coset] == 5}
    c_fixed_frames = {orbit[0] for orbit in c_omega_orbits if len(orbit) == 1}
    c_long_cells = {
        coset for orbit in c_omega_orbits if len(orbit) == 13 for coset in orbit
    }
    check(c_fixed_frames == boundary_fibre, "the fixed frames are not the cusp fibre")
    check(len(c_long_cells) == 78, "the six affine C13 cycles do not contain 78 frames")
    check(
        c_long_cells == {coset for coset in omega if coset_base[coset] != 5},
        "the C13 cycles are not the preimage of the affine chart",
    )
    check(projective_action(affine_chart, 5) == INF, "the root chart did not send cusp to infinity")

    # The commuting C6 action cyclically permutes both the six boundary
    # frames and the six affine C13 orbits.  On A-orbits it yields two
    # six-cycles, refining 84=12*7 without inventing LRC labels.
    c_long_orbits = {frozenset(orbit) for orbit in c_omega_orbits if len(orbit) == 13}
    a_orbit_set = {frozenset(orbit) for orbit in a_omega_orbits}
    frame_on_orbit = lambda orbit: frozenset(right_frame(coset) for coset in orbit)
    check(
        Counter(map(len, permutation_orbits(boundary_fibre, right_frame)))
        == Counter({6: 1}),
        "the six boundary frames are not one C6 orbit",
    )
    check(
        Counter(map(len, permutation_orbits(c_long_orbits, frame_on_orbit)))
        == Counter({6: 1}),
        "the six affine C13 cycles are not one C6 orbit",
    )
    check(
        Counter(map(len, permutation_orbits(a_orbit_set, frame_on_orbit)))
        == Counter({6: 2}),
        "the twelve C7 cycles are not two C6 frame orbits",
    )

    # The anharmonic S3 is the exact rerooting group of THM-2329.
    inversion = (0, 5, 5, 0)  # slope g -> 1/g; 5^2=-1 mod 13.
    tricycle = (0, -1 % P, 1, 1)  # slope g -> -1/(g+1).
    anharmonic = generated_group((inversion, tricycle))
    check(len(anharmonic) == 6, "the anharmonic group is not S3")
    check(psl_order(inversion) == 2, "wrong anharmonic involution")
    check(psl_order(tricycle) == 3, "wrong anharmonic tricycle")

    orbit_partition = []
    orbit_stabilizers = []
    orbit_seen: set[int] = set()
    for point in range(P + 1):
        if point in orbit_seen:
            continue
        orbit = tuple(sorted({projective_action(element, point) for element in anharmonic}))
        orbit_seen.update(orbit)
        stabilizer_size = sum(projective_action(element, point) == point for element in anharmonic)
        orbit_partition.append(orbit)
        orbit_stabilizers.append(stabilizer_size)
    check(
        tuple(orbit_partition)
        == ((0, 12, 13), (1, 6, 11), (2, 4, 5, 7, 8, 10), (3, 9)),
        "wrong anharmonic orbit partition",
    )
    check(tuple(orbit_stabilizers) == (2, 2, 1, 3), "wrong S3 stabilizers")
    check(all((root * root + root + 1) % P == 0 for root in (3, 9)), "wrong C3 roots")

    # The Barning branches are a free ternary semigroup mechanism, not C3
    # torsion.  Their mod-two image has only identity and coordinate swap.
    barning = (
        ((0, 1), (-1, 2)),
        ((0, 1), (1, 2)),
        ((1, 0), (2, 1)),
    )
    determinants = tuple(
        branch[0][0] * branch[1][1] - branch[0][1] * branch[1][0]
        for branch in barning
    )
    mod_two = tuple(
        tuple(tuple(value % 2 for value in row) for row in branch) for branch in barning
    )
    check(determinants == (1, -1, 1), "wrong Barning determinants")
    check(
        mod_two
        == (
            ((0, 1), (1, 0)),
            ((0, 1), (1, 0)),
            ((1, 0), (0, 1)),
        ),
        "the Barning mod-two image was not C2",
    )
    # Growth witnesses rule out order three over the positive rational cone.
    seed = (1, 2)
    first_children = tuple(integer_matrix_vector(branch, seed) for branch in barning)
    check(first_children == ((2, 3), (2, 5), (1, 4)), "wrong Barning children")
    check(all(gcd(*child) == 1 for child in first_children), "primitivity was lost")
    third_iterate = seed
    for _ in range(3):
        third_iterate = integer_matrix_vector(barning[2], third_iterate)
    check(third_iterate == (1, 8), "the third branch behaved like order-three torsion")

    # THM-2056 hostile: Farey/unimodular structure survives, but the fixed
    # Euclidean owner defect changes sign in both directions under branches.
    parent = (1, 10)
    child = integer_matrix_vector(barning[0], parent)
    grandchild = integer_matrix_vector(barning[2], child)
    check((child, grandchild) == ((10, 19), (10, 39)), "wrong hostile orbit")
    hostile_defects = (farey_defect(parent), farey_defect(child), farey_defect(grandchild))
    check(hostile_defects == (10, -449, 711), "wrong Farey-fan hostile defects")
    check(gcd(*parent) == gcd(*child) == gcd(*grandchild) == 1, "hostile lost visibility")

    print("THM-2603 Hurwitz/projective root-owner atlas exact companion")
    print("free_factors=x:order2,y:order3 generated_group=PSL2(F13):1092")
    print(f"owner_A=xy:{a} cycles={a_cycles}")
    print(f"root_C=[x,y]:{c} cycles={c_cycles}")
    print(f"companion_A^-1C:{d} cycles={d_cycles}")
    print(
        "nonabelian_norm=PASS order=6,5,4,3,2,1,0 "
        f"moving_cusps={tuple(conjugate_fixed_points)} reverse_order=13"
    )
    print(f"affine_chart=C:q->q+1 owner_cycles={affine_a_cycles} owner_8->inf->2")
    print("element_orders=1:1,2:91,3:182,6:182,7:468,13:168 order91=ABSENT")
    print("order7_trace_square_classes=t3:9x156,t5:12x156,t6:10x156 exhaustive=468")
    print("norm_closure_t3=F:6,8,9,10,11 R:2,3,4,5,7 union:2,3,4,5,6,7,8,9,10,11")
    print("norm_closure_t5=F:2,8,10,11,12 R:1,2,3,5,11 union:1,2,3,5,8,10,11,12")
    print("norm_closure_t6=F:1,3,9,11,12 R:1,2,4,10,12 union:1,2,3,4,9,10,11,12")
    print("norm_closure_total_union=1,2,3,4,5,6,7,8,9,10,11,12 orientation_R(a)=F(-a)^-1")
    print(
        "norm_signature_per_state="
        "identity:5,unipotent_nonsquare:2,order3:2,order6:2,order2:1"
    )
    print("norm_selector_cover_number=3 minimal=3F+3R+6F;3F+3R+6R trace5=REDUNDANT")
    print("norm_trace_orientation_physical_selector=NOT_CONSTRUCTED")
    print("centralizers=C7:7,C13:13 commuting_nontrivial_C7xC13_pairs=0")
    print(
        "root_normalizer=78=C13:13*C6:6 "
        "frame_conjugation_exponents=1,3,4,9,10,12"
    )
    print(
        "homogeneous_carrier=G/H~G/U:84 projective_fibres=14x6 "
        "left_A_orbits=12x7 right_frame_orbits=14x6"
    )
    print(
        "left_C_orbits=6x1+6x13 boundary_frames=6 affine_frames=78 "
        "frame_on_A_orbits=2x6"
    )
    print("left_C7_right_C6=COMMUTE physical_equivariant_bijection=NOT_CONSTRUCTED")
    print(
        "anharmonic_S3_orbits="
        "boundary3:(0,12,inf);harmonic3:(1,6,11);"
        "generic6:(2,4,5,7,8,10);equianharmonic2:(3,9)"
    )
    print("anharmonic_stabilizers=2,2,1,3 transverse_split=3+6+2")
    print("barning_mod2=swap,swap,identity ternary_arity_is_not_C3")
    print(f"farey_hostile_defects={hostile_defects} safe->unsafe->safe")
    print("physical_common_LRC_carrier=NOT_CONSTRUCTED")
    print(f"ALL_CHECKS=PASS checks={CHECKS}")


if __name__ == "__main__":
    main()
