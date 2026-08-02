#!/usr/bin/env python3
"""Exact companion for THM-3147's endpoint profile-jet observer.

The script works on the complete THM-3136 active-prefix bank through degree
nine.  It checks the one- and two-variable profile kernels under every
successive pole subtraction, identifies the rank-tail facets seen by the
length jet, and verifies that one singleton-count variable detects precisely
the remaining finite nonprincipal-required failures.
"""

import argparse
from collections import Counter, defaultdict, deque
from contextlib import redirect_stdout
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
MAXIMUM_DEGREE = 9


def load(name, filename):
    spec = spec_from_file_location(name, HERE / filename)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


a = load(
    "thm3136_companion",
    "gmc_one_sided_fixed_reference_hasse_no_go_thm3136.py",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def shifted_phi_all(coefficients, base_powers, removed):
    """Evaluate every monomial type through degree nine after removing R."""

    atom_count = len(coefficients)
    removed_powers = tuple(
        sum(root ** degree for root in removed)
        for degree in range(1, MAXIMUM_DEGREE + 1)
    )
    power_arrays = tuple(
        tuple(base_powers[atom][degree] - removed_powers[degree]
              for atom in range(atom_count))
        for degree in range(MAXIMUM_DEGREE)
    )
    values = {(): (1,) * atom_count}
    phi = {(): sum(coefficients)}
    for shape in a.ALL_SHAPES:
        exponent, remainder, merged_terms, multiplicity = a.RECURRENCE[shape]
        array = [
            power_arrays[exponent - 1][atom] * values[remainder][atom]
            for atom in range(atom_count)
        ]
        for merged, coefficient in merged_terms:
            merged_array = values[merged]
            for atom in range(atom_count):
                array[atom] -= coefficient * merged_array[atom]
        require(all(value % multiplicity == 0 for value in array),
                "profile monomial recurrence lost integrality")
        array = tuple(value // multiplicity for value in array)
        values[shape] = array
        phi[shape] = sum(
            coefficient * value
            for coefficient, value in zip(coefficients, array)
        )
    return phi


def two_variable_profile(phi):
    """Return coefficients of sum Phi(m_lambda) t^N u^ell w^m1."""

    out = defaultdict(int)
    out[(0, 0, 0)] = phi[()]
    for shape in a.ALL_SHAPES:
        out[(sum(shape), len(shape), shape.count(1))] += phi[shape]
    return dict(out)


def profile_value(profile, degree, length, singleton):
    if degree < 0 or length < 0 or singleton < 0:
        return 0
    return profile.get((degree, length, singleton), 0)


def endpoint_length_jet(profile):
    """Taylor coefficients at u=1 after setting w=1."""

    answer = {}
    for degree in range(MAXIMUM_DEGREE + 1):
        length_coefficients = [
            sum(profile_value(profile, degree, length, singleton)
                for singleton in range(length + 1))
            for length in range(degree + 1)
        ]
        for order in range(degree + 1):
            answer[(degree, order)] = sum(
                comb(length, order) * length_coefficients[length]
                for length in range(order, degree + 1)
            )
    return answer


def verify_profile_kernel(before, after, pole):
    """Check the exact length and length-singleton subtraction kernels."""

    before_jet = endpoint_length_jet(before)
    after_jet = endpoint_length_jet(after)
    length_checks = bivariate_checks = 0
    for degree in range(MAXIMUM_DEGREE + 1):
        for order in range(degree + 1):
            lhs = after_jet[(degree, order)]
            rhs = (
                before_jet.get((degree, order), 0)
                - pole * before_jet.get((degree - 1, order), 0)
                - pole * after_jet.get((degree - 1, order - 1), 0)
            )
            require(lhs == rhs, "length endpoint-jet kernel drift")
            length_checks += 1

        for length in range(degree + 1):
            for singleton in range(length + 1):
                lhs = (
                    profile_value(after, degree, length, singleton)
                    - pole * profile_value(
                        after, degree - 1, length, singleton)
                    + pole * profile_value(
                        after, degree - 1, length - 1, singleton - 1)
                    + pole ** 2 * profile_value(
                        after, degree - 2, length - 1, singleton)
                    - pole ** 2 * profile_value(
                        after, degree - 2, length - 1, singleton - 1)
                )
                rhs = (
                    profile_value(before, degree, length, singleton)
                    - pole * profile_value(
                        before, degree - 1, length, singleton)
                )
                require(lhs == rhs,
                        "length-singleton profile kernel drift")
                bivariate_checks += 1
    return length_checks, bivariate_checks


def connected_subset(vertices, edges):
    vertices = set(vertices)
    if not vertices:
        return False
    adjacency = {vertex: set() for vertex in vertices}
    for fine, coarse in edges:
        if fine in vertices and coarse in vertices:
            adjacency[fine].add(coarse)
            adjacency[coarse].add(fine)
    seen = {next(iter(vertices))}
    queue = deque(seen)
    while queue:
        vertex = queue.popleft()
        for neighbor in adjacency[vertex] - seen:
            seen.add(neighbor)
            queue.append(neighbor)
    return seen == vertices


def verify_observer_facets():
    """Check finite facet geometry for all rank tails and the N=7 sidecar."""

    facet_checks = 0
    for degree in range(2, MAXIMUM_DEGREE + 1):
        types, edges, _, closures = a.POSETS.get(
            degree, (None, None, None, None))
        if types is None:
            types = tuple(a.s.partitions(degree))
            edges = a.s.hasse_edges(types)
            closures = a.s.upper_closures(types, edges)
        universe = set(range(len(types)))
        for maximum_length in range(1, degree):
            upset = {
                index for index, shape in enumerate(types)
                if len(shape) <= maximum_length
            }
            require(all(closures[index] <= upset for index in upset),
                    "rank tail is not an upset")
            require(connected_subset(upset, edges)
                    and connected_subset(universe - upset, edges),
                    "rank-tail facet connectivity drift")
            facet_checks += 1

    degree = 7
    types, edges, _, closures = a.POSETS[degree]
    universe = set(range(len(types)))
    singleton_upset = {
        index for index, shape in enumerate(types)
        if len(shape) <= 3
        or (len(shape) == 4 and shape.count(1) >= 2)
    }
    require(all(closures[index] <= singleton_upset
                for index in singleton_upset),
            "singleton-refined set is not an upset")
    require(connected_subset(singleton_upset, edges)
            and connected_subset(universe - singleton_upset, edges),
            "singleton-refined facet connectivity drift")
    minimal = a.upset_minimal_shapes(
        frozenset(singleton_upset), types, edges)
    require(minimal == ((4, 1, 1, 1), (3, 2, 1, 1)),
            "singleton-refined minimal antichain drift")
    return facet_checks, minimal


def current_vector(phi, q_values, degree):
    shapes = tuple(a.g.partitions(degree))
    phi_h = sum(phi[shape] for shape in shapes)
    q_h = sum(q_values[shape] for shape in shapes)
    vector = {
        shape: Fraction(phi_h * q_values[shape] - phi[shape] * q_h)
        for shape in shapes
    }
    require(sum(vector.values()) == 0,
            "profile observer current lost zero mass")
    return vector, phi_h, q_h


def rank_tail_from_endpoint(length_coefficients, maximum_length):
    degree = len(length_coefficients) - 1
    endpoint = [
        sum(comb(length, order) * length_coefficients[length]
            for length in range(order, degree + 1))
        for order in range(degree + 1)
    ]
    require(endpoint[0] == 0, "endpoint current is not zero mass")
    reconstructed = sum(
        (-1) ** (order - maximum_length)
        * comb(order - 1, maximum_length) * endpoint[order]
        for order in range(maximum_length + 1, degree + 1)
    )
    direct = sum(length_coefficients[:maximum_length + 1])
    require(reconstructed == direct, "rank-tail endpoint inversion drift")
    return direct, tuple(endpoint)


def run_audit():
    facet_checks, singleton_antichain = verify_observer_facets()
    prefix_cases = transform_steps = 0
    length_kernel_checks = bivariate_kernel_checks = 0
    stats = {degree: Counter() for degree in range(5, 10)}
    first_length = first_singleton = None
    scalar_jet_hostile = None
    digest = sha256()

    for support_a, support_b in a.f.support_universe():
        for invariant in (0, 1):
            numerator, _, remaining, _, _, _ = a.f.reduced_row_fraction(
                invariant, support_a, support_b)
            polynomial = tuple(numerator[5:])
            poles = tuple(sorted(remaining.elements(), reverse=True))
            coefficients = tuple(
                coefficient for coefficient, _ in a.g.BANKS[invariant])
            atom_roots = tuple(
                tuple(a.g.residual_roots(
                    invariant, rows, support_a, support_b))
                for _, rows in a.g.BANKS[invariant]
            )
            base_powers = tuple(
                tuple(sum(root ** degree for root in roots)
                      for degree in range(1, MAXIMUM_DEGREE + 1))
                for roots in atom_roots
            )
            q_roots = tuple(a.g.residual_roots(
                invariant, a.dominant_rows(invariant), support_a, support_b))
            q_values = a.monomial_values(q_roots, ())
            previous_profile = None

            for prefix_length in range(len(polynomial)):
                prefix_cases += 1
                removed = poles[:prefix_length]
                phi = shifted_phi_all(coefficients, base_powers, removed)
                profile = two_variable_profile(phi)
                if previous_profile is not None:
                    length_count, bivariate_count = verify_profile_kernel(
                        previous_profile, profile, poles[prefix_length - 1])
                    length_kernel_checks += length_count
                    bivariate_kernel_checks += bivariate_count
                    transform_steps += 1
                previous_profile = profile

                for degree in range(5, 10):
                    vector, phi_h, _ = current_vector(
                        phi, q_values, degree)
                    types, edges, upsets, closures = a.POSETS[degree]
                    stats[degree]["total"] += 1
                    principal_bad = any(
                        sum(vector[types[index]] for index in principal) < 0
                        for principal in closures
                    )
                    stats[degree]["principal_bad"] += principal_bad

                    negative_tails = tuple(
                        (maximum_length, sum(
                            value for shape, value in vector.items()
                            if len(shape) <= maximum_length))
                        for maximum_length in range(1, degree)
                        if sum(value for shape, value in vector.items()
                               if len(shape) <= maximum_length) < 0
                    )
                    stats[degree]["ranktail_bad"] += bool(negative_tails)
                    if not principal_bad and negative_tails:
                        stats[degree]["length_observer"] += 1
                        require(all(maximum_length == degree - 2
                                    for maximum_length, _ in negative_tails),
                                "unexpected principal-pass rank tail")
                        record = (
                            degree, support_a, support_b, invariant,
                            prefix_length, removed, negative_tails, vector,
                        )
                        if first_length is None:
                            first_length = record

                    singleton_mass = None
                    if degree == 7:
                        singleton_mass = sum(
                            value for shape, value in vector.items()
                            if len(shape) <= 3
                            or (len(shape) == 4 and shape.count(1) >= 2)
                        )
                        if (not principal_bad and not negative_tails
                                and singleton_mass < 0):
                            stats[degree]["singleton_observer"] += 1
                            cut, upset = a.upset_minimum(
                                vector, types, upsets)
                            minimal = a.upset_minimal_shapes(
                                upset, types, edges)
                            require(
                                cut == singleton_mass
                                and minimal
                                == ((4, 1, 1, 1), (3, 2, 1, 1)),
                                "singleton hostile minimum upset drift",
                            )
                            stats[degree]["singleton_minimum_checked"] += 1
                            record = (
                                degree, support_a, support_b, invariant,
                                prefix_length, removed, singleton_mass, vector,
                            )
                            if first_singleton is None:
                                first_singleton = record

                    digest.update(repr((
                        degree, support_a, support_b, invariant,
                        prefix_length, phi_h, principal_bad,
                        negative_tails, singleton_mass,
                    )).encode())

                if (support_a, support_b, invariant, prefix_length) == (
                        1, 2, 0, 1):
                    degree = 6
                    length_coefficients = tuple(
                        sum(profile_value(profile, degree, length, singleton)
                            for singleton in range(length + 1))
                        for length in range(degree + 1)
                    )
                    endpoint = tuple(
                        sum(comb(length, order)
                            * length_coefficients[length]
                            for length in range(order, degree + 1))
                        for order in range(degree + 1)
                    )
                    require(length_coefficients
                            == (0, 0, 0, 20, 368, 264, -40)
                            and endpoint
                            == (612, 2612, 4308, 3332, 1088, 24, -40),
                            "scalar-positive endpoint-jet hostile drift")
                    scalar_jet_hostile = (length_coefficients, endpoint)

    require(prefix_cases == 8241, "active-prefix count drift")
    require(transform_steps == 8011, "successive-prefix transform count drift")
    require(length_kernel_checks == 55 * transform_steps,
            "length-kernel check count drift")
    require(bivariate_kernel_checks == 220 * transform_steps,
            "bivariate-kernel check count drift")
    expected = {
        5: (0, 0, 0, 0),
        6: (5641, 5641, 0, 0),
        7: (5386, 5278, 19, 16),
        8: (6510, 6451, 0, 0),
        9: (6767, 6647, 8, 0),
    }
    for degree, row in expected.items():
        actual = (
            stats[degree]["principal_bad"],
            stats[degree]["ranktail_bad"],
            stats[degree]["length_observer"],
            stats[degree]["singleton_observer"],
        )
        require(actual == row, f"degree-{degree} observer census drift")
    require(sum(stats[n]["length_observer"] for n in stats) == 27
            and sum(stats[n]["singleton_observer"] for n in stats) == 16,
            "27+16 observer split drift")
    require(first_length is not None and first_singleton is not None
            and scalar_jet_hostile is not None,
            "distinguished observer controls missing")

    length_degree, length_a, length_b, length_i, length_j, length_r, \
        length_tails, length_vector = first_length
    length_coefficients = tuple(
        sum(value for shape, value in length_vector.items()
            if len(shape) == length)
        for length in range(length_degree + 1)
    )
    length_mass, length_endpoint = rank_tail_from_endpoint(
        length_coefficients, length_degree - 2)
    require(
        (length_degree, length_a, length_b, length_i, length_j, length_r)
        == (7, 1, 2, 0, 1, (5,))
        and length_mass == -1151498720
        and (-length_endpoint[6] + 6 * length_endpoint[7])
        == -1151498720,
        "first rank-tail observer witness drift",
    )

    singleton_degree, singleton_a, singleton_b, singleton_i, singleton_j, \
        singleton_r, singleton_mass, singleton_vector = first_singleton
    require(
        (singleton_degree, singleton_a, singleton_b,
         singleton_i, singleton_j)
        == (7, 2, 10, 1, 24)
        and singleton_mass == -2816353101792484700160000,
        "first singleton-refined observer witness drift",
    )
    singleton_length_profile = tuple(
        sum(value for shape, value in singleton_vector.items()
            if len(shape) == length)
        for length in range(8)
    )
    require(all(sum(singleton_length_profile[:maximum_length + 1]) >= 0
                for maximum_length in range(1, 7)),
            "length-invisible hostile acquired a bad rank tail")

    stopping = a.all_degree_stopping_control()
    require(stopping[6][-1] == 1760080,
            "THM-3136 partial-reference boundary drift")

    print("THM-3147 length-singleton endpoint-jet facet observer")
    print("length_profile=Z_X(t,u)=prod_x(1+(u-1)xt)/(1-xt)")
    print("length_pole_kernel=(1-Mt)/(1-(1-u)Mt)")
    print("endpoint_recurrence=A'_(n,j)=A_(n,j)-M*A_(n-1,j)-M*A'_(n-1,j-1)")
    print("singleton_profile=prod_x(1+(uw-1)xt+u(1-w)x^2t^2)/(1-xt)")
    print("singleton_pole_kernel=(1-Mt)/(1+(uw-1)Mt+u(1-w)M^2t^2)")
    print(
        f"prefixes={prefix_cases}:transforms={transform_steps}:"
        f"length_kernel_checks={length_kernel_checks}:"
        f"singleton_kernel_checks={bivariate_kernel_checks}"
    )
    print(f"facet_connectivity_checks={facet_checks}:singleton_antichain={singleton_antichain}")
    for degree in range(5, 10):
        row = stats[degree]
        print(
            f"degree=N{degree}:principal_bad={row['principal_bad']}:"
            f"ranktail_bad={row['ranktail_bad']}:"
            f"principal_pass_length_caught={row['length_observer']}:"
            f"principal_pass_singleton_caught={row['singleton_observer']}:"
            f"singleton_minima_checked={row['singleton_minimum_checked']}"
        )
    print("nonprincipal_completion=length:27+singleton:16=43")
    print(
        "first_length_witness="
        f"N{length_degree}:support=({length_a},{length_b}):I{length_i+1}:"
        f"j{length_j}:R={length_r}:mass={length_mass}"
    )
    print(
        "first_length_endpoint="
        f"A6={length_endpoint[6]}:A7={length_endpoint[7]}:"
        f"-A6+6A7={length_mass}"
    )
    print(
        "first_singleton_witness="
        f"N{singleton_degree}:support=({singleton_a},{singleton_b}):"
        f"I{singleton_i+1}:j{singleton_j}:R_first_last="
        f"({singleton_r[0]},{singleton_r[-1]}):mass={singleton_mass}"
    )
    print(
        "length_only_kernel_hostile="
        "delta_(4,1,1,1)-delta_(2,2,2,1):"
        "same_length_profile;singleton_facet_separates"
    )
    print(
        "scalar_endpoint_jet_hostile="
        f"L6={scalar_jet_hostile[0]}:A6={scalar_jet_hostile[1]}"
    )
    print(
        "partial_reference_boundary="
        f"N6_bottom_Bernstein_min={stopping[6][-1]}"
    )
    print("scope=finite_profile_observer;not_full_Hasse_or_original_response")
    print(f"semantic_sha256={digest.hexdigest()}")
    print("all_exact_checks=PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is None:
        run_audit()
    else:
        with args.output.open("w", encoding="utf-8", newline="\n") as handle:
            with redirect_stdout(handle):
                run_audit()


if __name__ == "__main__":
    main()
