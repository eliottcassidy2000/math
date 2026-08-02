#!/usr/bin/env python3
"""Exact companion for THM-3136's one-sided fixed-reference Hasse no-go.

Only the signed bank functional is plethystically shifted by the active
reduced-pole prefix R_j.  The distinguished residual alphabet Q is kept
physical and unshifted:

    H_mu = Phi^j(h_N) m_mu(Q) - Phi^j(m_mu) h_N(Q).

The script scans the complete THM-3120 support/prefix bank in degrees 5..9,
checks every coarsening upset and the equivalent max-flow certificate, and
compares its raw labelled-deletion image with the co-shifted current ruled
out by THM-3128.
"""

import argparse
from collections import Counter
from contextlib import redirect_stdout
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import comb, factorial
from pathlib import Path


HERE = Path(__file__).resolve().parents[1] / "04-computation"
MAXIMUM_DEGREE = 9


def load(name, filename):
    spec = spec_from_file_location(name, HERE / filename)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


h = load(
    "thm3128_companion",
    "gmc_active_prefix_labelled_deletion_no_preimage_thm3128.py",
)
s = load(
    "thm3127_companion",
    "gmc_partition_refinement_strassen_thm3127.py",
)
g = h.g
f = h.f

ALL_SHAPES = tuple(sorted(
    (shape
     for degree in range(1, MAXIMUM_DEGREE + 1)
     for shape in g.partitions(degree)),
    key=lambda shape: (sum(shape), len(shape), shape),
))
RECURRENCE = {}
for shape in ALL_SHAPES:
    exponent = shape[-1]
    remainder = shape[:-1]
    merged_terms = []
    for old_exponent in set(remainder):
        merged = list(remainder)
        merged.remove(old_exponent)
        merged.append(old_exponent + exponent)
        merged = tuple(sorted(merged, reverse=True))
        merged_terms.append((merged, Counter(merged)[old_exponent + exponent]))
    RECURRENCE[shape] = (
        exponent, remainder, tuple(merged_terms), Counter(shape)[exponent])


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dominant_rows(invariant):
    if invariant == 0:
        rows = ((0, 3), (2, 0), (2, 0), (1, 0))
    else:
        rows = ((0, 3), (1, 1), (2, 0), (1, 0), (1, 0))
    return tuple(sorted(rows, reverse=True))


def monomial_values(roots, removed):
    return h.all_monomial_values(
        roots, removed, maximum_degree=MAXIMUM_DEGREE)


def elementary_polynomial(roots):
    answer = [1]
    for root in roots:
        out = [0] * (len(answer) + 1)
        for degree, value in enumerate(answer):
            out[degree] += value
            out[degree + 1] += root * value
        answer = out
    return answer


def all_degree_stopping_control():
    """Reconstruct the scalar-positive infinite bottom-facet hostile."""

    invariant, a, b, removed_pole = 0, 1, 2, 5
    numerator, denominator, remaining, _, _, _ = f.reduced_row_fraction(
        invariant, a, b)
    polynomial = tuple(numerator[5:])
    poles = tuple(sorted(remaining.elements(), reverse=True))
    require(
        polynomial == (36, -108, -72)
        and poles == (5, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "selected reduced row fraction drift",
    )
    tail_flag, _ = f.pole_flag(polynomial, poles[1:])
    require(tail_flag == (36, 144, 72),
            "one-sided scalar re-flag drift")
    require(all(root > 0 for root in poles[1:])
            and len(poles) - 1 >= len(polynomial),
            "positive suffix-denominator hypotheses drift")

    signed_elementary = [0]
    for coefficient, rows in g.BANKS[invariant]:
        roots = g.residual_roots(invariant, rows, a, b)
        signed_elementary = f.add_scaled(
            signed_elementary, elementary_polynomial(roots), coefficient)
    require(tuple(signed_elementary) == (0, 0, 0, 0, 0, 16, 40),
            "unshifted elementary response drift")

    shifted_elementary = []
    for degree in range(21):
        numerator_value = (
            signed_elementary[degree]
            if degree < len(signed_elementary) else 0
        )
        shifted_elementary.append(
            numerator_value
            - removed_pole * (shifted_elementary[-1] if degree else 0)
        )
    require(shifted_elementary[5] == 16, "shifted e5 drift")
    require(all(
        shifted_elementary[degree]
        == -40 * (-5) ** (degree - 6)
        for degree in range(6, 21)
    ), "shifted elementary geometric tail drift")

    q_roots = tuple(g.residual_roots(
        invariant, dominant_rows(invariant), a, b))
    require(q_roots == (1, 1, 2, 3, 4, 5),
            "fixed distinguished alphabet drift")
    q_elementary = tuple(elementary_polynomial(q_roots))
    require(q_elementary == (1, 16, 100, 310, 499, 394, 120),
            "fixed-Q elementary polynomial drift")

    tail_denominator = f.divide_linear_if_exact(denominator, removed_pole)
    require(tail_denominator is not None,
            "selected pole did not divide row denominator")
    scalar = f.series_coefficients(numerator, tail_denominator, 20)
    require(all(scalar[degree] > 0 for degree in range(5, 21)),
            "finite scalar-positive control drift")
    q_denominator = [1]
    for root in q_roots:
        q_denominator = f.multiply_linear(q_denominator, root)
    q_complete = f.series_coefficients([1], q_denominator, 20)
    h6_bottom = (
        scalar[6] * q_elementary[6]
        - shifted_elementary[6] * q_complete[6]
    )
    require((scalar[6], q_complete[6], h6_bottom)
            == (612, 297412, 11969920),
            "degree-six bottom coordinate drift")
    for degree in range(8, 21, 2):
        bottom = -shifted_elementary[degree] * q_complete[degree]
        require(bottom == 40 * 5 ** (degree - 6) * q_complete[degree] > 0,
                "even elementary-tail hostile drift")

    # Interpolate plethystically between the fixed reference Q and the
    # co-shifted reference Q-5 by Q_z=Q-[5z], 0<=z<=1.  Then
    # H_(Q_z)=H_Q(1-5zt) and E_(Q_z)=E_Q/(1+5zt).  In degree six the
    # bottom coordinate is an exact polynomial in z.  Its Bernstein
    # coefficients certify positivity on the entire interval at once.
    qz_e6 = tuple(
        q_elementary[6 - degree] * (-5) ** degree
        for degree in range(7)
    )
    qz_h6 = (q_complete[6], -5 * q_complete[5])
    partial_bottom = tuple(
        scalar[6] * qz_e6[degree]
        - shifted_elementary[6]
        * (qz_h6[degree] if degree < len(qz_h6) else 0)
        for degree in range(7)
    )
    require(
        partial_bottom
        == (11969920, -11342040, 7634700, -23715000,
            38250000, -30600000, 9562500),
        "partial-reference power polynomial drift",
    )
    bernstein = tuple(
        sum(
            Fraction(partial_bottom[degree] * comb(index, degree),
                     comb(6, degree))
            for degree in range(index + 1)
        )
        for index in range(7)
    )
    require(
        bernstein
        == (11969920, 10079580, 8698220, 6640090,
            5269440, 3400520, 1760080)
        and all(bernstein[index] > bernstein[index + 1]
                for index in range(6)),
        "partial-reference Bernstein certificate drift",
    )
    return (
        polynomial, poles, tail_flag, q_roots, h6_bottom,
        partial_bottom, bernstein,
    )


def efficient_upsets(types, closures):
    """Enumerate upsets without iterating through all 2^p(N) subsets."""

    q = len(types)
    order = tuple(sorted(range(q), key=lambda i: (len(types[i]), types[i])))
    position = {vertex: i for i, vertex in enumerate(order)}
    require(
        all(position[coarse] <= position[fine]
            for fine in range(q) for coarse in closures[fine]),
        "coarsest-first order drift",
    )
    answer = []

    def recurse(k, chosen):
        if k == q:
            answer.append(frozenset(chosen))
            return
        vertex = order[k]
        recurse(k + 1, chosen)
        if all(coarse == vertex or coarse in chosen
               for coarse in closures[vertex]):
            chosen.add(vertex)
            recurse(k + 1, chosen)
            chosen.remove(vertex)

    recurse(0, set())
    return tuple(answer)


POSETS = {}
EXPECTED_UPSETS = {5: 10, 6: 27, 7: 47, 8: 168, 9: 573}
for degree in range(5, MAXIMUM_DEGREE + 1):
    types = tuple(s.partitions(degree))
    edges = s.hasse_edges(types)
    closures = s.upper_closures(types, edges)
    upsets = efficient_upsets(types, closures)
    require(len(upsets) == EXPECTED_UPSETS[degree],
            f"upset census drift at N={degree}")
    POSETS[degree] = (types, edges, upsets, closures)


def upset_minimum(vector, types, upsets):
    values = tuple(vector[shape] for shape in types)
    records = tuple((sum(values[i] for i in upset), upset)
                    for upset in upsets)
    return min(records, key=lambda item: (item[0], len(item[1]), tuple(item[1])))


def upset_minimal_shapes(upset, types, edges):
    chosen = set(upset)
    return tuple(
        types[vertex]
        for vertex in sorted(chosen)
        if not any(fine in chosen and coarse == vertex
                   for fine, coarse in edges)
    )


def weight(shape):
    answer = 1
    for part in shape:
        answer *= factorial(part)
    return answer


def lowered(shape, part):
    target = list(shape)
    target.remove(part)
    if part > 1:
        target.append(part - 1)
    return tuple(sorted(target, reverse=True))


def raw_labelled_delete(vector):
    """Apply B=W_(N-1) A W_N^-1 to raw Young-gap coefficients."""

    factorial_target = Counter()
    for shape, raw in vector.items():
        coordinate = Fraction(raw, weight(shape))
        for part, multiplicity in Counter(shape).items():
            factorial_target[lowered(shape, part)] += (
                coordinate * part * multiplicity
            )
    return {
        shape: value * weight(shape)
        for shape, value in factorial_target.items()
        if value
    }


def prefix_phi(coefficients, base_powers, removed):
    """Evaluate the whole signed bank in one vectorized Newton recurrence."""

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
    phi = {}
    for shape in ALL_SHAPES:
        exponent, remainder, merged_terms, multiplicity = RECURRENCE[shape]
        array = [
            power_arrays[exponent - 1][atom] * values[remainder][atom]
            for atom in range(atom_count)
        ]
        for merged, coefficient in merged_terms:
            merged_array = values[merged]
            for atom in range(atom_count):
                array[atom] -= coefficient * merged_array[atom]
        require(all(value % multiplicity == 0 for value in array),
                "vectorized monomial recurrence lost integrality")
        array = tuple(value // multiplicity for value in array)
        values[shape] = array
        if sum(shape) >= 5:
            phi[shape] = sum(
                coefficient * value
                for coefficient, value in zip(coefficients, array)
            )
    return phi


def coefficient_vectors(phi, q_values):
    answer = {}
    for degree in range(5, MAXIMUM_DEGREE + 1):
        shapes = tuple(g.partitions(degree))
        phi_h = sum(phi[shape] for shape in shapes)
        q_h = sum(q_values[shape] for shape in shapes)
        require(phi[(degree,)] == 0,
                "common-prefix row marginal cancellation failed")
        vector = {
            shape: Fraction(phi_h * q_values[shape] - phi[shape] * q_h)
            for shape in shapes
        }
        require(sum(vector.values()) == 0,
                "fixed-reference current lost zero mass")
        answer[degree] = (vector, phi_h, q_h)
    return answer


def format_vector(vector):
    return ";".join(
        f"{shape}:{value}"
        for shape, value in sorted(
            vector.items(), key=lambda item: (len(item[0]), item[0]))
        if value
    )


def run_audit():
    stopping = all_degree_stopping_control()
    prefix_cases = current_cases = upset_checks = flow_checks = 0
    positive_scalar = zero_scalar = negative_scalar = 0
    fixed_passes = 0
    first_fixed_failure = None
    first_interior_failure = None
    first_nonprincipal_required = None
    minimum_cut = None
    hostile_comparison = None
    degree_stats = {degree: Counter() for degree in range(5, MAXIMUM_DEGREE + 1)}
    digest = sha256()

    for support_index, (a, b) in enumerate(f.support_universe(), 1):
        for invariant in (0, 1):
            numerator, _, remaining, _, _, _ = f.reduced_row_fraction(
                invariant, a, b)
            polynomial = tuple(numerator[5:])
            poles = tuple(sorted(remaining.elements(), reverse=True))
            degree_flag = len(polynomial) - 1
            require(degree_flag <= len(poles), "active pole prefix overflow")
            flag, _ = f.pole_flag(polynomial, poles)
            require(len(flag) == degree_flag + 1 and all(x > 0 for x in flag),
                    "active scalar Newton flag compatibility failed")
            q_roots = tuple(g.residual_roots(
                invariant, dominant_rows(invariant), a, b))
            q_fixed = monomial_values(q_roots, ())
            coefficients = tuple(
                coefficient for coefficient, _ in g.BANKS[invariant])
            atom_roots = tuple(
                tuple(g.residual_roots(invariant, rows, a, b))
                for _, rows in g.BANKS[invariant]
            )
            base_powers = tuple(
                tuple(sum(root ** degree for root in roots)
                      for degree in range(1, MAXIMUM_DEGREE + 1))
                for roots in atom_roots
            )
            for prefix_length in range(degree_flag + 1):
                prefix_cases += 1
                removed = poles[:prefix_length]
                phi = prefix_phi(coefficients, base_powers, removed)
                fixed = coefficient_vectors(phi, q_fixed)
                for degree in range(5, MAXIMUM_DEGREE + 1):
                    current_cases += 1
                    vector, phi_h, q_h = fixed[degree]
                    if phi_h > 0:
                        positive_scalar += 1
                    elif phi_h == 0:
                        zero_scalar += 1
                    else:
                        negative_scalar += 1

                    types, edges, upsets, closures = POSETS[degree]
                    cut, upset = upset_minimum(vector, types, upsets)
                    upset_checks += len(upsets)
                    flow = s.maxflow_feasible(
                        tuple(vector[shape] for shape in types), edges)
                    flow_checks += 1
                    passed = cut >= 0
                    digest.update(
                        repr((degree, a, b, invariant, prefix_length,
                              phi_h, q_h, tuple(vector[shape] for shape in types),
                              cut)).encode()
                    )
                    stats = degree_stats[degree]
                    stats["total"] += 1
                    stats["scalar_positive"] += phi_h > 0
                    stats["scalar_negative"] += phi_h < 0
                    stats["pass"] += passed
                    stats["positive_scalar_pass"] += passed and phi_h > 0
                    stats["negative_scalar_pass"] += passed and phi_h < 0
                    top_bad = vector[(degree,)] < 0
                    bottom_bad = vector[(1,) * degree] > 0
                    principal_bad = []
                    for root, principal in enumerate(closures):
                        mass = sum(vector[types[i]] for i in principal)
                        if mass < 0:
                            principal_bad.append((mass, root, principal))
                    extreme_bad = top_bad or bottom_bad
                    stats["top_bad"] += top_bad
                    stats["bottom_complement_bad"] += bottom_bad
                    stats["extreme_union_bad"] += extreme_bad
                    stats["principal_bad"] += bool(principal_bad)
                    stats["interior_fail"] += (not passed and not extreme_bad)
                    stats["interior_principal_caught"] += (
                        not passed and not extreme_bad and bool(principal_bad))
                    stats["interior_nonprincipal_required"] += (
                        not passed and not extreme_bad and not principal_bad)
                    stats["nonprincipal_required_including_bottom"] += (
                        not passed and not principal_bad)
                    require(flow == passed,
                            "upset/max-flow duality disagreed on fixed reference")
                    if passed:
                        fixed_passes += 1
                    elif first_fixed_failure is None:
                        first_fixed_failure = (
                            degree, a, b, invariant, prefix_length, removed,
                            phi_h, q_h, cut,
                            tuple(types[i] for i in sorted(upset)), vector,
                        )
                    if not passed and not extreme_bad:
                        record = (
                            degree, a, b, invariant, prefix_length, removed,
                            phi_h, cut,
                            tuple(types[i] for i in sorted(upset)),
                            tuple((mass, types[root])
                                  for mass, root, _ in principal_bad),
                            vector,
                        )
                        key = (
                            degree, a, b, invariant, prefix_length,
                            len(upset), tuple(types[i] for i in sorted(upset)),
                        )
                        if (first_interior_failure is None
                                or key < first_interior_failure[0]):
                            first_interior_failure = (key, record)
                        if not principal_bad and (
                                first_nonprincipal_required is None
                                or key < first_nonprincipal_required[0]):
                            first_nonprincipal_required = (key, record)
                    if minimum_cut is None or cut < minimum_cut[0]:
                        minimum_cut = (
                            cut, degree, a, b, invariant, prefix_length,
                            tuple(types[i] for i in sorted(upset)),
                        )

                    # Load-bearing comparison with THM-3128's first active
                    # co-shifted hostile.  The fixed reference must escape by
                    # changing the labelled-deletion image, not by moving in
                    # its forbidden affine fibre.
                    if (degree, a, b, invariant, prefix_length) == (5, 1, 3, 1, 5):
                        q_coshift = monomial_values(q_roots, removed)
                        old_vector = coefficient_vectors(phi, q_coshift)[degree][0]
                        expected_old = dict(zip(
                            types,
                            map(Fraction, (
                                -2324160, 544320, 2237760, -656640,
                                -915840, 972000, 142560,
                            )),
                        ))
                        require(old_vector == expected_old,
                                "THM-3128 co-shift hostile drift")
                        expected_fixed = dict(zip(
                            types,
                            map(Fraction, (
                                89261280, 381574080, 306501120, 716713920,
                                621743040, -721199520, -1394593920,
                            )),
                        ))
                        require(vector == expected_fixed,
                                "fixed-reference repair vector drift")
                        old_image = raw_labelled_delete(old_vector)
                        fixed_image = raw_labelled_delete(vector)
                        top = (4,)
                        old_top = old_image.get(top, 0)
                        fixed_top = fixed_image.get(top, 0)
                        require(old_top == -1779840 and fixed_top == 470835360,
                                "top deletion coordinate repair drift")
                        require(old_image != fixed_image,
                                "fixed reference remained in forbidden deletion fibre")
                        old_cut, old_upset = upset_minimum(old_vector, types, upsets)
                        require(old_cut < 0 and not s.maxflow_feasible(
                            tuple(old_vector[shape] for shape in types), edges),
                            "THM-3128 hostile unexpectedly passed")
                        hostile_comparison = (
                            old_cut, cut, old_top, fixed_top,
                            tuple(types[i] for i in sorted(old_upset)),
                        )

    require(prefix_cases == 8241, "active-prefix census drift")
    require(current_cases == 41205, "N5..9 current census drift")
    require(hostile_comparison is not None,
            "THM-3128 fixed-reference comparison was not reached")
    require((positive_scalar, zero_scalar, negative_scalar, fixed_passes)
            == (33162, 0, 8043, 16858),
            "global fixed-reference census drift")
    expected_degree_stats = {
        5: (8241, 8241, 8241, 8241, 0, 0, 0, 0, 0),
        6: (8241, 2600, 5298, 2600, 2943, 5641, 5641, 0, 0),
        7: (8241, 2820, 7041, 2820, 1200, 1492, 2692, 2729, 35),
        8: (8241, 1731, 5638, 1731, 2603, 6268, 6276, 234, 0),
        9: (8241, 1466, 6944, 1466, 1297, 1983, 3240, 3535, 8),
    }
    for degree, expected in expected_degree_stats.items():
        stats = degree_stats[degree]
        actual = (
            stats["total"], stats["pass"], stats["scalar_positive"],
            stats["positive_scalar_pass"], stats["top_bad"],
            stats["bottom_complement_bad"], stats["extreme_union_bad"],
            stats["interior_fail"], stats["interior_nonprincipal_required"],
        )
        require(actual == expected, f"degree-{degree} facet census drift")
    require(sum(degree_stats[n]["interior_principal_caught"]
                for n in expected_degree_stats) == 6455,
            "interior-principal count drift")
    require(sum(degree_stats[n]["interior_nonprincipal_required"]
                for n in expected_degree_stats) == 43,
            "nonprincipal-required count drift")
    require(first_nonprincipal_required is not None,
            "first nonprincipal witness disappeared")
    nonprincipal_record = first_nonprincipal_required[1]
    np_degree, _, _, _, _, _, _, np_mass, np_upset = nonprincipal_record[:9]
    np_types, np_edges, _, _ = POSETS[np_degree]
    np_indices = frozenset(np_types.index(shape) for shape in np_upset)
    require(
        np_mass == -1151498720
        and upset_minimal_shapes(np_indices, np_types, np_edges)
        == ((3, 1, 1, 1, 1), (2, 2, 1, 1, 1)),
        "first nonprincipal upset/minimal antichain drift",
    )
    print("THM-3136 one-sided fixed-reference elementary-tail Hasse no-go")
    print(
        "all_degree_control="
        f"P:{stopping[0]};poles:{stopping[1]};tail_flag:{stopping[2]};"
        f"Q:{stopping[3]};H_1^6:{stopping[4]}"
    )
    print(
        "partial_reference_N6="
        f"power:{stopping[5]};bernstein:{stopping[6]};"
        f"uniform_bottom>={stopping[6][-1]}"
    )
    print("elementary_response=(16*t^5+40*t^6)/(1+5*t)")
    print("infinite_hostile=scalar_positive_all_N>=5;bottom_facet_fails_even_N>=6")
    print("current=Phi^R(h_N)m_mu(Q)-Phi^R(m_mu)h_N(Q)")
    print("holotopy=bank_prefix_moves;distinguished_Q_fixed")
    print(f"active_prefixes={prefix_cases}:N5..9_currents={current_cases}")
    print(f"scalar_PhiR_hN=positive:{positive_scalar}:zero:{zero_scalar}:negative:{negative_scalar}")
    print(f"upset_sums_checked={upset_checks}:maxflows_checked={flow_checks}")
    print(f"fixed_reference_passes={fixed_passes}/{current_cases}")
    print(f"thm3128_hostile_comparison={hostile_comparison}")
    print("thm3128_deletion_image=fixed_reference_changes_image:PASS")
    print(
        "global_minimum_cut="
        f"mass:{minimum_cut[0]}:N{minimum_cut[1]}:"
        f"support:({minimum_cut[2]},{minimum_cut[3]}):"
        f"I{minimum_cut[4]+1}:j{minimum_cut[5]}:upset_size:{len(minimum_cut[6])}"
    )
    for degree in range(5, MAXIMUM_DEGREE + 1):
        stats = degree_stats[degree]
        print(
            f"degree=N{degree}:total={stats['total']}:pass={stats['pass']}:"
            f"scalar_positive={stats['scalar_positive']}:"
            f"positive_scalar_pass={stats['positive_scalar_pass']}:"
            f"scalar_negative={stats['scalar_negative']}:"
            f"negative_scalar_pass={stats['negative_scalar_pass']}:"
            f"top_bad={stats['top_bad']}:"
            f"bottom_complement_bad={stats['bottom_complement_bad']}:"
            f"extreme_union_bad={stats['extreme_union_bad']}:"
            f"principal_bad={stats['principal_bad']}:"
            f"interior_fail={stats['interior_fail']}:"
            f"interior_principal={stats['interior_principal_caught']}:"
            f"interior_nonprincipal={stats['interior_nonprincipal_required']}:"
            f"allprincipal_pass_but_fail={stats['nonprincipal_required_including_bottom']}"
        )
    if first_fixed_failure is None:
        print("fixed_reference_hasse=ALL_PASS")
    else:
        print(f"fixed_reference_first_failure={first_fixed_failure[:-1]}")
        print(f"fixed_reference_failure_vector={format_vector(first_fixed_failure[-1])}")
    if first_interior_failure is None:
        print("first_interior_failure=NONE")
    else:
        print(f"first_interior_failure={first_interior_failure[1][:-1]}")
        print(f"first_interior_vector={format_vector(first_interior_failure[1][-1])}")
    if first_nonprincipal_required is None:
        print("first_nonprincipal_required=NONE")
    else:
        print(
            "first_nonprincipal_required="
            f"{first_nonprincipal_required[1][:-1]}"
        )
        print(
            "first_nonprincipal_vector="
            f"{format_vector(first_nonprincipal_required[1][-1])}"
        )
    print("scope=finite_exact_fixed_reference_current;no_nonrow_reconstruction_claim")
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
