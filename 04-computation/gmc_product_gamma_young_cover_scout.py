#!/usr/bin/env python3
"""Exact Young-lattice scout for the anchored product-Gamma response banks.

This is evidence for a prospective all-degree extension of THM-3110, not a
proof.  It imports the exact atom construction from the THM-3110 companion,
then tests the dominant-normalized response

    R_i(lambda) = Phi_i(s_lambda)/(B_i c_i s_lambda(Q_i))

on a stated finite universe.  Here Q_i is the residual alphabet of the
distinguished positive atom, c_i is its coefficient, and B_i is the forced
positive polynomial divisor from the THM-3110 companion.  It also probes two
possible proof mechanisms:

* monotonicity of each individual negative-atom Schur ratio; and
* a fixed capacitated transport from negative atoms to non-dominant positive
  atoms using the stronger pairwise decrement order.

It also separates three tempting Young-lattice reductions: minimization of the
response itself on the one-row ray, minimization of every cover increment by
the one-row cover, and monotonicity by added-box content.  These are logically
different, and the latter two acquire exact hostiles in the census.

Only Python integers and fractions are used.
"""

from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


COMPANION = Path(__file__).with_name(
    "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
)
SPEC = spec_from_file_location("thm3110_companion", COMPANION)
g = module_from_spec(SPEC)
SPEC.loader.exec_module(g)


# Rows are written as exponent pairs alpha*a+beta*b.  Sorting matches the
# canonical bank representation in the imported companion.
DOMINANT_ROWS = (
    tuple(sorted(((0, 3), (2, 0), (2, 0), (1, 0)), reverse=True)),
    tuple(sorted(((0, 3), (1, 1), (2, 0), (1, 0), (1, 0)),
                 reverse=True)),
)
DOMINANT_COEFFICIENTS = (1, 2)
MAXIMUM_DEGREE = 16
GENERIC_COVER_SOURCE_DEGREES = (5, 6, 7)
ROW_TAIL_MAXIMUM_DEGREE = 128


def young_covers(partition):
    """Return each partition obtained by adding one Young-diagram box."""

    out = []
    for row in range(len(partition)):
        if row == 0 or partition[row] < partition[row - 1]:
            out.append(partition[:row] + (partition[row] + 1,)
                       + partition[row + 1:])
    out.append(partition + (1,))
    # The construction is duplicate-free, but retain a guard against drift.
    require(len(out) == len(set(out)), "duplicate Young cover")
    return tuple(out)


def added_box_profile(mu, nu):
    """Return (row, column, coarm, coleg, content) for a Young cover."""

    padded = mu + (0,) * (len(nu) - len(mu))
    changed = tuple(index for index, (left, right)
                    in enumerate(zip(padded, nu)) if right != left)
    require(len(changed) == 1, "cover does not add exactly one row-part")
    row_index = changed[0]
    require(nu[row_index] == padded[row_index] + 1,
            "cover does not add exactly one box")
    row = row_index + 1
    column = nu[row_index]
    coarm = column - 1
    coleg = row - 1
    return row, column, coarm, coleg, coarm - coleg


def chamber(a, b):
    return "tight" if b < 2 * a else "wide"


def interval_pole_factor(length, maximum):
    """Return prod_(r=0)^(length-1) (1-r/maximum)^(-1)."""

    require(0 <= length <= maximum, "invalid pole-factor interval")
    out = Fraction(1)
    for root in range(length):
        out /= 1 - Fraction(root, maximum)
    return out


def simple_top_pole_constant(roots, maximum):
    """Leading constant of h_n(roots)/maximum^n at a simple top root."""

    if maximum not in roots:
        return Fraction(0)
    require(roots.count(maximum) == 1, "top pole is not simple")
    out = Fraction(1)
    skipped = False
    for root in roots:
        if root == maximum and not skipped:
            skipped = True
            continue
        out /= 1 - Fraction(root, maximum)
    return out


def support_universe():
    return tuple(
        (a, b)
        for a in range(1, 11)
        for b in range(a + 1, min(3 * a + 4, 21) + 1)
    )


def maximum_transport(invariant, edges):
    """Maximum negative-mass to non-dominant-positive-mass transport."""

    bank = g.BANKS[invariant]
    negative = tuple(index for index, (coefficient, _) in enumerate(bank)
                     if coefficient < 0)
    positive = tuple(index for index, (coefficient, rows) in enumerate(bank)
                     if coefficient > 0
                     and rows != DOMINANT_ROWS[invariant])
    source = len(negative) + len(positive)
    target = source + 1
    flow = g.Dinic(target + 1)
    for local, index in enumerate(negative):
        flow.add_edge(source, local, -bank[index][0])
    for local, index in enumerate(positive):
        flow.add_edge(len(negative) + local, target, bank[index][0])
    negative_local = {index: local for local, index in enumerate(negative)}
    positive_local = {index: local for local, index in enumerate(positive)}
    total_mass = sum(bank[index][0] for index in positive)
    for negative_index, positive_index in edges:
        flow.add_edge(
            negative_local[negative_index],
            len(negative) + positive_local[positive_index],
            total_mass,
        )
    return total_mass, flow.maximum_flow(source, target)


def generic_cover_newton_audit():
    """Prove the first three cover ranks on both full integer chambers.

    For |mu|=n and nu covering mu, the cross-multiplied numerator divided by
    B_i has chamber degree at most

        (2(n+1)-9) + 2n = 4n-7.

    The first summand uses the THM-3110 response degree bound and the second
    uses the elementary degree-2n bound for s_mu(Q_i).  Exact Newton recovery,
    one excess degree shell, and three off-grid points certify the polynomial.
    """

    expected = {
        5: (13, 7980, 7980, 0, Fraction(2620)),
        6: (17, 20520, 20518, 2, Fraction(6000)),
        7: (21, 45540, 45534, 6, Fraction(532800)),
    }
    results = {}
    zero_records = []
    maximum_target_degree = max(GENERIC_COVER_SOURCE_DEGREES) + 1

    for invariant in (0, 1):
        for name in ("tight", "wide"):
            audit_points = set()
            for source_degree in GENERIC_COVER_SOURCE_DEGREES:
                polynomial_degree = 4 * source_degree - 7
                audit_points.update({
                    (first, second)
                    for first in range(polynomial_degree + 2)
                    for second in range(polynomial_degree + 2 - first)
                })
                audit_points.update({
                    (polynomial_degree + 2, polynomial_degree + 1),
                    (polynomial_degree + 1, polynomial_degree + 2),
                    (polynomial_degree + 2, 2),
                })

            cached = {}
            for first, second in sorted(audit_points):
                a, b = g.chamber_parameters(name, first, second)
                cached[first, second] = (
                    a,
                    b,
                    g.complete_and_elementary(
                        g.residual_roots(
                            invariant, DOMINANT_ROWS[invariant], a, b),
                        maximum_target_degree,
                    ),
                    g.atom_symmetric_data(
                        invariant, a, b, maximum_target_degree),
                )

            for source_degree in GENERIC_COVER_SOURCE_DEGREES:
                polynomial_degree = 4 * source_degree - 7
                key = source_degree
                result = results.setdefault(
                    key, {"slots": 0, "positive": 0, "zero": 0,
                          "minimum": None})
                for mu in g.partitions(source_degree):
                    for nu in young_covers(mu):
                        values = {}
                        for first in range(polynomial_degree + 2):
                            for second in range(
                                    polynomial_degree + 2 - first):
                                a, b, q_data, data = cached[first, second]
                                cross = (
                                    g.phi_from_data(data, nu)
                                    * g.schur_value(q_data, mu)
                                    - g.phi_from_data(data, mu)
                                    * g.schur_value(q_data, nu)
                                )
                                values[first, second] = Fraction(
                                    cross, g.forced_divisor(invariant, a, b))

                        coefficients = g.newton_triangle(
                            values, polynomial_degree)
                        for index, coefficient in coefficients.items():
                            total_order = sum(index)
                            if total_order <= polynomial_degree:
                                result["slots"] += 1
                                require(coefficient >= 0,
                                        "negative generic cover coefficient")
                                if coefficient == 0:
                                    result["zero"] += 1
                                    zero_records.append((
                                        source_degree, invariant + 1, name,
                                        mu, nu, index,
                                    ))
                                else:
                                    result["positive"] += 1
                                    current = result["minimum"]
                                    if current is None or coefficient < current:
                                        result["minimum"] = coefficient
                            elif total_order == polynomial_degree + 1:
                                require(coefficient == 0,
                                        "generic cover degree guard failed")

                        for first, second in (
                            (polynomial_degree + 2, polynomial_degree + 1),
                            (polynomial_degree + 1, polynomial_degree + 2),
                            (polynomial_degree + 2, 2),
                        ):
                            a, b, q_data, data = cached[first, second]
                            exact = Fraction(
                                (g.phi_from_data(data, nu)
                                 * g.schur_value(q_data, mu)
                                 - g.phi_from_data(data, mu)
                                 * g.schur_value(q_data, nu)),
                                g.forced_divisor(invariant, a, b),
                            )
                            reconstructed = g.evaluate_newton(
                                coefficients, first, second,
                                polynomial_degree)
                            require(exact == reconstructed,
                                    "generic cover off-grid guard failed")

    lines = []
    for source_degree in GENERIC_COVER_SOURCE_DEGREES:
        polynomial_degree, slots, positive, zero, minimum = expected[
            source_degree]
        result = results[source_degree]
        require((result["slots"], result["positive"], result["zero"],
                 result["minimum"]) == (slots, positive, zero, minimum),
                "generic cover census drift")
        lines.append(
            f"generic_cover={source_degree}->{source_degree + 1}:"
            f"degree={polynomial_degree}:slots={slots}:positive={positive}:"
            f"zero={zero}:negative=0:min={minimum}"
        )
    require(len(zero_records) == 8, "generic cover zero-record drift")
    return tuple(lines), tuple(zero_records)


def row_tail_audit(supports):
    """Stress the complete-homogeneous ray far beyond the Schur census."""

    checks = positive = below_limit = increasing = 0
    minimum_cross = None
    for a, b in supports:
        for invariant, bank in enumerate(g.BANKS):
            q_complete = g.complete_and_elementary(
                g.residual_roots(
                    invariant, DOMINANT_ROWS[invariant], a, b),
                ROW_TAIL_MAXIMUM_DEGREE,
            )[0]
            atom_complete = tuple(
                g.complete_and_elementary(
                    g.residual_roots(invariant, rows, a, b),
                    ROW_TAIL_MAXIMUM_DEGREE,
                )[0]
                for _, rows in bank
            )
            phi = tuple(
                sum(coefficient * values[degree]
                    for (coefficient, _), values in zip(bank, atom_complete))
                for degree in range(ROW_TAIL_MAXIMUM_DEGREE + 1)
            )
            normalization = (
                g.forced_divisor(invariant, a, b)
                * DOMINANT_COEFFICIENTS[invariant]
            )
            maximum_root = 3 * b - 1
            p_a = interval_pole_factor(a, maximum_root)
            p_2a = interval_pole_factor(2 * a, maximum_root)
            x = Fraction(p_a * p_a, p_2a)
            if invariant == 0:
                limit = Fraction(
                    (1 - x) ** 2,
                    g.forced_divisor(invariant, a, b),
                )
            else:
                y = Fraction(
                    p_2a * interval_pole_factor(b, maximum_root),
                    interval_pole_factor(a + b, maximum_root) * p_a,
                )
                limit = Fraction(
                    (1 - x) * (1 - y),
                    g.forced_divisor(invariant, a, b),
                )

            for degree in range(5, ROW_TAIL_MAXIMUM_DEGREE + 1):
                checks += 1
                require(phi[degree] > 0, "row-tail positivity failed")
                positive += 1
                response = Fraction(
                    phi[degree], normalization * q_complete[degree])
                require(response < limit,
                        "row response crossed its top-pole limit")
                below_limit += 1
                if degree > 5:
                    cross = (
                        phi[degree] * q_complete[degree - 1]
                        - phi[degree - 1] * q_complete[degree]
                    )
                    require(cross > 0, "row-tail monotonicity failed")
                    increasing += 1
                    record = (
                        cross, invariant + 1, a, b, degree - 1, degree)
                    if minimum_cross is None or cross < minimum_cross[0]:
                        minimum_cross = record

    require((checks, positive, below_limit, increasing)
            == (28520, 28520, 28520, 28290),
            "row-tail census drift")
    require(minimum_cross == (29433312, 1, 1, 2, 5, 6),
            "row-tail minimum-cross drift")
    return (
        f"one_row_tail=degrees5..{ROW_TAIL_MAXIMUM_DEGREE}:"
        f"checks={checks}:positive={positive}:increasing={increasing}:"
        f"below_limit={below_limit}:minimum_cross={minimum_cross}"
    )


def main():
    supports = support_universe()
    require(len(supports) == 115, "support-universe count drift")
    generic_lines, generic_zero_records = generic_cover_newton_audit()
    row_tail_line = row_tail_audit(supports)

    # Candidate decrement-order edges are intersected over every cover in each
    # chamber.  An edge N->P survives if d_N >= d_P everywhere tested, where
    # d_S=s_mu(S)/s_mu(Q)-s_nu(S)/s_nu(Q).
    transport_edges = {}
    ordered_transport_edges = {}
    for invariant, bank in enumerate(g.BANKS):
        negative = tuple(index for index, (coefficient, _) in enumerate(bank)
                         if coefficient < 0)
        positive = tuple(index for index, (coefficient, rows) in enumerate(bank)
                         if coefficient > 0
                         and rows != DOMINANT_ROWS[invariant])
        for name in ("tight", "wide"):
            transport_edges[invariant, name] = {
                (negative_index, positive_index)
                for negative_index in negative
                for positive_index in positive
            }
            ordered_transport_edges[invariant, name] = set(
                transport_edges[invariant, name])

    cover_count = 0
    individual_negative_count = 0
    individual_negative_first_hostile = None
    minimum_delta = None
    minimum_ratio = None
    grouped_minimum_delta = {}
    grouped_minimum_ratio = {}
    box_profiles = set()
    row_shape_checks = 0
    row_shape_equality_count = 0
    row_shape_first_hostile = None
    row_shape_minimum_margin = None
    row_shape_minimum_positive_margin = None
    row_response_minimum = None
    row_response_degree_minimum = {}
    row_cover_checks = 0
    row_cover_violation_count = 0
    row_cover_first_hostile = None
    row_cover_minimum_gap = None
    content_universal_first_hostile = None
    content_minima_first_hostile = None
    top_pole_formula_checks = 0

    for a, b in supports:
        name = chamber(a, b)
        for invariant, bank in enumerate(g.BANKS):
            dominant_index = next(
                index for index, (coefficient, rows) in enumerate(bank)
                if (coefficient, rows) == (
                    DOMINANT_COEFFICIENTS[invariant],
                    DOMINANT_ROWS[invariant],
                )
            )
            q_data = g.complete_and_elementary(
                g.residual_roots(
                    invariant, DOMINANT_ROWS[invariant], a, b),
                MAXIMUM_DEGREE,
            )
            atom_data = tuple(
                g.complete_and_elementary(
                    g.residual_roots(invariant, rows, a, b),
                    MAXIMUM_DEGREE,
                )
                for _, rows in bank
            )
            atom_roots = tuple(
                g.residual_roots(invariant, rows, a, b)
                for _, rows in bank
            )

            # The root-ordered candidate is stronger than empirical decrement
            # order: it additionally requires S_negative<=S_positive
            # coordinatewise throughout the chamber sample.
            ordered_active = ordered_transport_edges[invariant, name]
            root_order_failures = {
                edge for edge in ordered_active
                if any(left > right for left, right in zip(
                    atom_roots[edge[0]], atom_roots[edge[1]]))
            }
            ordered_active.difference_update(root_order_failures)

            shapes = tuple(
                partition
                for degree in range(5, MAXIMUM_DEGREE + 1)
                for partition in g.partitions(degree)
            )
            q_values = {
                partition: g.schur_value(q_data, partition)
                for partition in shapes
            }
            atom_values = tuple(
                {
                    partition: g.schur_value(data, partition)
                    for partition in shapes
                }
                for data in atom_data
            )
            phi_values = {
                partition: sum(
                    coefficient * atom_values[index][partition]
                    for index, (coefficient, _) in enumerate(bank)
                )
                for partition in shapes
            }

            normalization = (
                g.forced_divisor(invariant, a, b)
                * DOMINANT_COEFFICIENTS[invariant]
            )

            maximum_root = 3 * b - 1
            q_pole = simple_top_pole_constant(
                g.residual_roots(
                    invariant, DOMINANT_ROWS[invariant], a, b),
                maximum_root,
            )
            response_pole = sum(
                coefficient * simple_top_pole_constant(roots, maximum_root)
                for (coefficient, _), roots in zip(bank, atom_roots)
            )
            pole_limit = Fraction(response_pole, normalization * q_pole)
            p_a = interval_pole_factor(a, maximum_root)
            p_2a = interval_pole_factor(2 * a, maximum_root)
            x = Fraction(p_a * p_a, p_2a)
            if invariant == 0:
                expected_pole_limit = Fraction(
                    (1 - x) ** 2,
                    g.forced_divisor(invariant, a, b),
                )
            else:
                p_b = interval_pole_factor(b, maximum_root)
                p_ab = interval_pole_factor(a + b, maximum_root)
                y = Fraction(p_2a * p_b, p_ab * p_a)
                expected_pole_limit = Fraction(
                    (1 - x) * (1 - y),
                    g.forced_divisor(invariant, a, b),
                )
            require(0 < x < 1, "first pole cross-ratio left (0,1)")
            if invariant == 1:
                require(0 < y < 1, "second pole cross-ratio left (0,1)")
            require(pole_limit == expected_pole_limit,
                    "row top-pole factorization failed")
            top_pole_formula_checks += 1

            # This is weaker than cover monotonicity and is tested separately:
            # at fixed degree, does the normalized response attain its minimum
            # on the one-row shape?
            for degree in range(5, MAXIMUM_DEGREE + 1):
                row_shape = (degree,)
                row_q = q_values[row_shape]
                row_phi = phi_values[row_shape]
                require(row_q > 0, "dominant row alphabet vanished")
                row_response = Fraction(
                    row_phi, normalization * row_q)
                row_response_record = (
                    row_response, invariant + 1, a, b, degree)
                if (row_response_minimum is None
                        or row_response < row_response_minimum[0]):
                    row_response_minimum = row_response_record
                if (degree not in row_response_degree_minimum
                        or row_response
                        < row_response_degree_minimum[degree][0]):
                    row_response_degree_minimum[degree] = row_response_record
                for partition in g.partitions(degree):
                    partition_q = q_values[partition]
                    if partition_q == 0:
                        continue
                    margin = Fraction(
                        phi_values[partition] * row_q
                        - row_phi * partition_q,
                        normalization * partition_q * row_q,
                    )
                    row_shape_checks += 1
                    record = (
                        margin, invariant + 1, a, b, degree, partition,
                    )
                    if margin == 0:
                        row_shape_equality_count += 1
                    elif (row_shape_minimum_positive_margin is None
                          or margin < row_shape_minimum_positive_margin[0]):
                        row_shape_minimum_positive_margin = record
                    if (row_shape_minimum_margin is None
                            or margin < row_shape_minimum_margin[0]):
                        row_shape_minimum_margin = record
                    if margin < 0:
                        key = degree, invariant, a, b, partition
                        if (row_shape_first_hostile is None
                                or key < row_shape_first_hostile[0]):
                            row_shape_first_hostile = key, record

            for degree in range(5, MAXIMUM_DEGREE):
                row_mu = (degree,)
                row_nu = (degree + 1,)
                row_cross = (
                    phi_values[row_nu] * q_values[row_mu]
                    - phi_values[row_mu] * q_values[row_nu]
                )
                row_delta = Fraction(
                    row_cross,
                    normalization * q_values[row_mu] * q_values[row_nu],
                )
                content_buckets = {}
                for mu in g.partitions(degree):
                    q_mu = q_values[mu]
                    if q_mu == 0:
                        continue
                    phi_mu = phi_values[mu]
                    require(phi_mu > 0, "tested source response is not positive")
                    for nu in young_covers(mu):
                        q_nu = q_values[nu]
                        if q_nu == 0:
                            continue
                        phi_nu = phi_values[nu]
                        require(phi_nu > 0,
                                "tested target response is not positive")
                        cross = phi_nu * q_mu - phi_mu * q_nu
                        require(cross > 0,
                                "dominant-normalized cover monotonicity failed")
                        cover_count += 1
                        profile = added_box_profile(mu, nu)
                        box_profiles.add(profile)

                        delta = Fraction(
                            cross,
                            normalization * q_mu * q_nu,
                        )
                        content_buckets.setdefault(profile[-1], []).append(
                            (delta, mu, nu, profile)
                        )
                        row_cover_checks += 1
                        row_gap = delta - row_delta
                        row_record = (
                            row_gap, invariant + 1, a, b, degree, mu, nu,
                            row_delta, delta, profile,
                        )
                        if (row_cover_minimum_gap is None
                                or row_gap < row_cover_minimum_gap[0]):
                            row_cover_minimum_gap = row_record
                        if row_gap < 0:
                            row_cover_violation_count += 1
                            key = degree, invariant, a, b, mu, nu
                            if (row_cover_first_hostile is None
                                    or key < row_cover_first_hostile[0]):
                                row_cover_first_hostile = key, row_record
                        ratio = Fraction(phi_nu * q_mu, phi_mu * q_nu)
                        delta_record = (
                            delta, invariant + 1, a, b, degree, mu, nu)
                        ratio_record = (
                            ratio, invariant + 1, a, b, degree, mu, nu)
                        if minimum_delta is None or delta < minimum_delta[0]:
                            minimum_delta = delta_record
                        if minimum_ratio is None or ratio < minimum_ratio[0]:
                            minimum_ratio = ratio_record
                        group = invariant, name
                        if (group not in grouped_minimum_delta
                                or delta < grouped_minimum_delta[group][0]):
                            grouped_minimum_delta[group] = delta_record
                        if (group not in grouped_minimum_ratio
                                or ratio < grouped_minimum_ratio[group][0]):
                            grouped_minimum_ratio[group] = ratio_record

                        # Every decrement has the same positive denominator
                        # q_mu*q_nu, so pairwise comparisons use integer
                        # numerators only.
                        decrements = tuple(
                            atom_values[index][mu] * q_nu
                            - atom_values[index][nu] * q_mu
                            for index in range(len(bank))
                        )
                        require(decrements[dominant_index] == 0,
                                "dominant decrement drift")
                        for index, (coefficient, rows) in enumerate(bank):
                            if coefficient >= 0:
                                continue
                            individual_negative_count += 1
                            if (decrements[index] < 0
                                    and individual_negative_first_hostile is None):
                                individual_negative_first_hostile = (
                                    invariant + 1, a, b, degree, mu, nu,
                                    coefficient, rows, decrements[index],
                                )

                        active = transport_edges[invariant, name]
                        failed = {
                            edge for edge in active
                            if decrements[edge[0]] < decrements[edge[1]]
                        }
                        active.difference_update(failed)
                        ordered_transport_edges[invariant, name].difference_update(
                            failed)

                # A content-prefix proof would require smaller-content boxes
                # to have no smaller increment than larger-content boxes.  We
                # test both the universal interval order and the weaker order
                # of the minimum increment in each content class.
                contents = sorted(content_buckets)
                if content_universal_first_hostile is None:
                    for lower_index, lower_content in enumerate(contents):
                        lower_minimum = min(content_buckets[lower_content])
                        for upper_content in contents[lower_index + 1:]:
                            upper_maximum = max(content_buckets[upper_content])
                            if lower_minimum[0] < upper_maximum[0]:
                                content_universal_first_hostile = (
                                    invariant + 1, a, b, degree,
                                    lower_content, upper_content,
                                    lower_minimum, upper_maximum,
                                )
                                break
                        if content_universal_first_hostile is not None:
                            break
                if content_minima_first_hostile is None:
                    minima = tuple(
                        (content, min(content_buckets[content]))
                        for content in contents
                    )
                    for (lower_content, lower_minimum), (
                            upper_content, upper_minimum) in zip(
                                minima, minima[1:]):
                        if lower_minimum[0] < upper_minimum[0]:
                            content_minima_first_hostile = (
                                invariant + 1, a, b, degree,
                                lower_content, upper_content,
                                lower_minimum, upper_minimum,
                            )
                            break

    require(cover_count == 555662, "Young-cover census drift")
    require(minimum_delta[0] == Fraction(
        125872476615673643,
        522634777460648197207357773696,
    ), "minimum cover increment drift")
    require(minimum_ratio[0] == Fraction(53329143091, 53305526863),
            "minimum cover ratio drift")

    print("status=HYPOTHESIS_EXACT_FINITE_EVIDENCE")
    for line in generic_lines:
        print(line)
    print("generic_cover_structural_zeros=" + repr(generic_zero_records))
    print("supports=a:1..10;b:a+1..min(3a+4,21);count=115")
    print("invariants=I1,I2;source_degrees=5..15;target_degrees=6..16")
    print(f"young_covers={cover_count}:strict_pass={cover_count}:equal=0:fail=0")
    print("R_i(lambda)=Phi_i(s_lambda)/(B_i*c_i*s_lambda(Q_i));c_1=1;c_2=2")
    print("B_1=a^4*b^2*(b-a)^3;B_2=a^3*b^2*(b-a)^4")
    print("Q_1=rows(3b,2a,2a,a);Q_2=rows(3b,a+b,2a,a,a)")
    print("minimum_delta=" + repr(minimum_delta))
    print("minimum_delta_box="
          + repr(added_box_profile(minimum_delta[-2], minimum_delta[-1])))
    print("minimum_ratio=" + repr(minimum_ratio))
    print("minimum_ratio_box="
          + repr(added_box_profile(minimum_ratio[-2], minimum_ratio[-1])))
    print(f"added_box_profiles={len(box_profiles)}:"
          f"coarm=0..{max(profile[2] for profile in box_profiles)}:"
          f"coleg=0..{max(profile[3] for profile in box_profiles)}:"
          f"content={min(profile[4] for profile in box_profiles)}.."
          f"{max(profile[4] for profile in box_profiles)}")
    if row_shape_first_hostile is None:
        print(
            f"one_row_response_minimum={row_shape_checks}:PASS:"
            f"equal={row_shape_equality_count}:"
            f"minimum_positive_margin={row_shape_minimum_positive_margin}"
        )
    else:
        print(
            f"one_row_response_minimum={row_shape_checks}:FAIL:"
            f"first={row_shape_first_hostile[1]}:"
            f"minimum_margin={row_shape_minimum_margin}:"
            f"equal={row_shape_equality_count}:"
            f"minimum_positive_margin={row_shape_minimum_positive_margin}"
        )
    print("one_row_response_global_minimum=" + repr(row_response_minimum))
    print("one_row_response_degree_minimizers=" + repr(tuple(
        (degree,) + row_response_degree_minimum[degree][1:]
        for degree in range(5, MAXIMUM_DEGREE + 1)
    )))
    print(f"one_row_top_pole_formula={top_pole_formula_checks}:PASS")
    print(row_tail_line)
    print(
        f"one_row_cover_minimum={row_cover_checks}:"
        f"violations={row_cover_violation_count}:"
        f"first={None if row_cover_first_hostile is None else row_cover_first_hostile[1]}:"
        f"minimum_gap={row_cover_minimum_gap}"
    )
    print("content_universal_order=FAIL:first="
          + repr(content_universal_first_hostile))
    print("content_minima_order=FAIL:first="
          + repr(content_minima_first_hostile))
    for invariant in (0, 1):
        for name in ("tight", "wide"):
            delta_record = grouped_minimum_delta[invariant, name]
            ratio_record = grouped_minimum_ratio[invariant, name]
            print(
                f"group_extrema=I{invariant + 1}:{name}:"
                f"minimum_delta={delta_record}:"
                f"box={added_box_profile(delta_record[-2], delta_record[-1])}:"
                f"minimum_ratio={ratio_record}:"
                f"box={added_box_profile(ratio_record[-2], ratio_record[-1])}"
            )
    print(f"individual_negative_tests={individual_negative_count}")
    if individual_negative_first_hostile is None:
        print("individual_negative_ratio_monotonicity=PASS")
    else:
        print("individual_negative_ratio_monotonicity=FAIL:first="
              + repr(individual_negative_first_hostile))

    for invariant in (0, 1):
        for name in ("tight", "wide"):
            edges = transport_edges[invariant, name]
            demand, flow = maximum_transport(invariant, edges)
            print(
                f"fixed_decrement_transport=I{invariant + 1}:{name}:"
                f"edges={len(edges)}:demand={demand}:flow={flow}:"
                f"residual={demand - flow}"
            )
            ordered_edges = ordered_transport_edges[invariant, name]
            ordered_demand, ordered_flow = maximum_transport(
                invariant, ordered_edges)
            print(
                f"root_ordered_decrement_transport=I{invariant + 1}:{name}:"
                f"edges={len(ordered_edges)}:demand={ordered_demand}:"
                f"flow={ordered_flow}:residual={ordered_demand - ordered_flow}"
            )

    print("implication=cover_theorem+degree5_base=>all_degree_schur_positivity")


if __name__ == "__main__":
    main()
