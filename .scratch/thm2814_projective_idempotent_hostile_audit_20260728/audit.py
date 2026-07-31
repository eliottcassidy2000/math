#!/usr/bin/env python3
"""Independent hostile audit of the proposed THM-2814 classification.

This companion deliberately does not import the proposed theorem companion.
It checks:

* the complete row/column-torus orbit classification on the unit locus;
* the distinction between a fixed common coefficient line and a gauge orbit;
* the exact scope boundary over rings and for unrelated vertex gauges;
* the field-idempotent no-go and linear provenance identity;
* the orientation and projective loss in the THM-2779 Weyl application;
* the honest-character boundary and quotient loss in THM-2807.
"""

from collections import Counter
from itertools import product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mobius(vector, modulus):
    return (vector[0] - vector[1] - vector[2] + vector[3]) % modulus


def cross_ratio(vector, modulus):
    v00, v10, v01, v11 = vector
    denominator = v10 * v01 % modulus
    require(denominator != 0, "cross-ratio denominator vanished")
    return v00 * v11 * pow(denominator, -1, modulus) % modulus


def row_column_transform(vector, gauge, modulus):
    """Apply (r0,r1,c0,c1) to (v00,v10,v01,v11)."""
    v00, v10, v01, v11 = vector
    r0, r1, c0, c1 = gauge
    return (
        r0 * c0 * v00 % modulus,
        r1 * c0 * v10 % modulus,
        r0 * c1 * v01 % modulus,
        r1 * c1 * v11 % modulus,
    )


def normalize_unit_square(vector, modulus):
    """An independently chosen gauge taking a unit square to (1,1,1,kappa)."""
    v00, v10, v01, _v11 = vector
    gauge = (
        1,
        v00 * pow(v10, -1, modulus) % modulus,
        pow(v00, -1, modulus),
        pow(v01, -1, modulus),
    )
    normalized = row_column_transform(vector, gauge, modulus)
    return normalized, gauge


def field_orbit_audit(modulus):
    units = tuple(range(1, modulus))
    fibres = {kappa: set() for kappa in units}
    normalization_checks = 0

    for vector in product(units, repeat=4):
        kappa = cross_ratio(vector, modulus)
        fibres[kappa].add(vector)
        normalized, gauge = normalize_unit_square(vector, modulus)
        require(
            normalized == (1, 1, 1, kappa),
            f"normal form failed over F_{modulus}",
        )
        require(
            cross_ratio(
                row_column_transform(vector, gauge, modulus),
                modulus,
            )
            == kappa,
            f"normalizing gauge changed kappa over F_{modulus}",
        )
        normalization_checks += 1

    expected_orbit_size = (modulus - 1) ** 3
    stabilizer_sizes = []
    for kappa in units:
        representative = (1, 1, 1, kappa)

        # Fix r0=1 to select one representative of the scalar kernel in the
        # four-parameter gauge group.
        effective_orbit = {
            row_column_transform(
                representative,
                (1, r1, c0, c1),
                modulus,
            )
            for r1, c0, c1 in product(units, repeat=3)
        }
        require(
            effective_orbit == fibres[kappa],
            f"kappa was not a complete orbit invariant over F_{modulus}",
        )
        require(
            len(effective_orbit) == expected_orbit_size,
            f"effective orbit size changed over F_{modulus}",
        )

        stabilizer = tuple(
            gauge
            for gauge in product(units, repeat=4)
            if row_column_transform(representative, gauge, modulus)
            == representative
        )
        expected_kernel = {
            (
                scalar,
                scalar,
                pow(scalar, -1, modulus),
                pow(scalar, -1, modulus),
            )
            for scalar in units
        }
        require(
            set(stabilizer) == expected_kernel,
            f"full-gauge stabilizer was not the scalar kernel over F_{modulus}",
        )
        stabilizer_sizes.append(len(stabilizer))

    kappa_one_mu = Counter(
        mobius(vector, modulus) for vector in fibres[1]
    )
    require(
        kappa_one_mu[0] > 0
        and sum(
            count
            for value, count in kappa_one_mu.items()
            if value != 0
        )
        > 0,
        f"kappa=1 orbit did not contain both additive behaviours over F_{modulus}",
    )

    branch_a_successes = tuple(
        (alpha, beta)
        for alpha, beta in product(units, repeat=2)
        if alpha != 1
        and beta != 1
        and mobius(
            (1, alpha, beta, alpha * beta % modulus),
            modulus,
        )
        != 0
    )
    require(
        len(branch_a_successes) == (modulus - 2) ** 2,
        f"fixed-line Segre census changed over F_{modulus}",
    )

    return {
        "field": modulus,
        "vectors": normalization_checks,
        "orbits": len(fibres),
        "orbit_size": expected_orbit_size,
        "stabilizer_sizes": tuple(sorted(set(stabilizer_sizes))),
        "kappa1_mu_histogram": tuple(sorted(kappa_one_mu.items())),
        "branch_a_nonzero_mu": len(branch_a_successes),
    }


def branch_and_scope_hostiles():
    modulus = 5

    # Same kappa=1 orbit: a physically meaningful fixed-line difference can
    # be nonzero, while row/column normalization makes the raw Mobius scalar
    # zero.
    segre = (1, 2, 2, 4)
    segre_normalized, segre_gauge = normalize_unit_square(segre, modulus)
    require(cross_ratio(segre, modulus) == 1, "Segre hostile gained holonomy")
    require(mobius(segre, modulus) == 1, "Segre hostile lost additive face")
    require(
        segre_normalized == (1, 1, 1, 1)
        and mobius(segre_normalized, modulus) == 0,
        "row/column gauge failed to erase chart-dependent Segre mu",
    )

    # Even for kappa != 1, raw mu=0/nonzero is not an orbit invariant.
    curved_zero = (1, 2, 2, 3)
    curved_normalized, curved_gauge = normalize_unit_square(
        curved_zero,
        modulus,
    )
    require(
        cross_ratio(curved_zero, modulus) == 2
        and mobius(curved_zero, modulus) == 0,
        "curved zero-face hostile changed",
    )
    require(
        curved_normalized == (1, 1, 1, 2)
        and mobius(curved_normalized, modulus) == 1,
        "normalized curved defect changed",
    )

    # Arbitrary independent vertex rephasing would destroy kappa.  The
    # row/column factorization of the gauge is therefore load-bearing.
    flat = (1, 1, 1, 1)
    arbitrary_vertex_rescale = (1, 1, 1, 2)
    require(
        cross_ratio(flat, modulus) == 1
        and cross_ratio(arbitrary_vertex_rescale, modulus) == 2,
        "unrelated-vertex gauge hostile changed",
    )

    # Over Z, equal formal cross-ratio does not classify unit-gauge orbits:
    # the only units are +/-1, so entry magnitudes cannot be normalized.
    integer_a = (1, 2, 3, 6)
    integer_b = (1, 1, 1, 1)
    integer_unit_orbit = {
        (
            r0 * c0 * integer_a[0],
            r1 * c0 * integer_a[1],
            r0 * c1 * integer_a[2],
            r1 * c1 * integer_a[3],
        )
        for r0, r1, c0, c1 in product((-1, 1), repeat=4)
    }
    require(integer_b not in integer_unit_orbit, "integer nonunit hostile failed")
    require(
        integer_a[0] * integer_a[3]
        == integer_a[1] * integer_a[2],
        "integer formal kappa changed",
    )

    # In Z/6, nonzero does not imply invertible, and this denominator is zero.
    z6_vector = (1, 2, 3, 1)
    z6_denominator = z6_vector[1] * z6_vector[2] % 6
    require(
        z6_denominator == 0 and gcd(z6_denominator, 6) != 1,
        "nonzero/nonunit ring hostile changed",
    )

    return {
        "segre": (
            segre,
            mobius(segre, modulus),
            segre_normalized,
            mobius(segre_normalized, modulus),
            segre_gauge,
        ),
        "curved_zero": (
            curved_zero,
            cross_ratio(curved_zero, modulus),
            mobius(curved_zero, modulus),
            curved_normalized,
            mobius(curved_normalized, modulus),
            curved_gauge,
        ),
        "unrelated_vertex": (
            cross_ratio(flat, modulus),
            cross_ratio(arbitrary_vertex_rescale, modulus),
        ),
        "integer_equal_formal_kappa_distinct_unit_orbits": (
            integer_a,
            integer_b,
            len(integer_unit_orbit),
        ),
        "z6_nonzero_denominator": z6_denominator,
    }


def idempotent_and_provenance_audit():
    modulus = 13
    idempotents = tuple(
        value
        for value in range(modulus)
        if value * value % modulus == value
    )
    require(idempotents == (0, 1), "field idempotents changed")

    common_faces = []
    for alpha, beta in product(idempotents, repeat=2):
        square = (1, alpha, beta, alpha * beta % modulus)
        if all(square):
            common_faces.append(mobius(square, modulus))
    require(common_faces == [0], "field-idempotent common atom gained mu")

    checks = 0
    atoms = 5
    weight_values = (-1, 0, 2)
    for left_word in product((0, 1), repeat=atoms):
        for right_word in product((0, 1), repeat=atoms):
            for weights in product(weight_values, repeat=atoms):
                bare = sum(weights)
                left = sum(
                    weights[index] * left_word[index]
                    for index in range(atoms)
                )
                right = sum(
                    weights[index] * right_word[index]
                    for index in range(atoms)
                )
                both = sum(
                    weights[index]
                    * left_word[index]
                    * right_word[index]
                    for index in range(atoms)
                )
                face = bare - left - right + both
                absent = sum(
                    weights[index]
                    * (1 - left_word[index])
                    * (1 - right_word[index])
                    for index in range(atoms)
                )
                require(face == absent, "linear provenance identity failed")
                checks += 1
    require(
        checks == 2**atoms * 2**atoms * len(weight_values) ** atoms,
        "provenance universe changed",
    )

    # The field hypothesis is essential.  In F2^4, nontrivial idempotents
    # can have nonzero product and a disjoint nonzero complement atom.
    one = (1, 1, 1, 1)
    alpha = (1, 1, 0, 0)
    beta = (1, 0, 1, 0)

    def component_product(left, right):
        return tuple(x * y % 2 for x, y in zip(left, right))

    alpha_beta = component_product(alpha, beta)
    ring_mu = tuple(
        (one[index] - alpha[index] - beta[index] + alpha_beta[index]) % 2
        for index in range(4)
    )
    require(
        all(any(value for value in vector) for vector in (one, alpha, beta, alpha_beta))
        and any(ring_mu),
        "product-ring field-scope hostile failed",
    )

    return {
        "field_idempotents": idempotents,
        "common_faces": tuple(common_faces),
        "linear_provenance_checks": checks,
        "product_ring_square": (one, alpha, beta, alpha_beta),
        "product_ring_mu": ring_mu,
    }


FIELD = 53
ORDER = 13


def primitive_thirteenth_root():
    for candidate in range(2, FIELD):
        if pow(candidate, ORDER, FIELD) != 1:
            continue
        if all(
            pow(candidate, exponent, FIELD) != 1
            for exponent in range(1, ORDER)
        ):
            return candidate
    raise RuntimeError("no primitive thirteenth root in F_53")


ZETA = primitive_thirteenth_root()


def thm2779_application_audit():
    ratios = []
    for root in range(ORDER):
        # Standard composition: (TM)e_r=T(M e_r).
        tm_target = (root + 1) % ORDER
        tm_coefficient = pow(ZETA, -root, FIELD)
        mt_target = (root + 1) % ORDER
        mt_coefficient = pow(ZETA, -(root + 1), FIELD)
        require(tm_target == mt_target, "Weyl paths changed endpoints")
        ratios.append(
            tm_coefficient * pow(mt_coefficient, -1, FIELD) % FIELD
        )
    require(
        set(ratios) == {ZETA},
        "oriented Weyl path ratio was not zeta",
    )

    reverse_ratio = pow(ZETA, -1, FIELD)
    cocycle_normal_form = (1, 1, 1, ZETA)
    require(
        cross_ratio(cocycle_normal_form, FIELD) == ZETA,
        "Weyl cocycle normal form changed",
    )

    one_dimensional_solutions = tuple(
        (left, right)
        for left, right in product(range(1, FIELD), repeat=2)
        if left * right % FIELD
        == ZETA * right * left % FIELD
    )
    require(
        not one_dimensional_solutions,
        "a nonzero one-dimensional Weyl pair appeared",
    )

    # On basis lines M is the identity permutation and T is a cycle.  The
    # scalar commutator is therefore killed both on those lines and in PGL.
    t_lines = tuple((root + 1) % ORDER for root in range(ORDER))
    m_lines = tuple(range(ORDER))
    tm_lines = tuple(t_lines[m_lines[root]] for root in range(ORDER))
    mt_lines = tuple(m_lines[t_lines[root]] for root in range(ORDER))
    require(tm_lines == mt_lines, "projective basis-line shadows did not commute")

    return {
        "field": FIELD,
        "zeta": ZETA,
        "oriented_TM_over_MT": tuple(sorted(set(ratios))),
        "reverse_orientation": reverse_ratio,
        "normal_form": cocycle_normal_form,
        "normal_form_mu": mobius(cocycle_normal_form, FIELD),
        "one_dimensional_solutions": len(one_dimensional_solutions),
        "basis_line_TM_equals_MT": tm_lines == mt_lines,
    }


def thm2807_application_audit():
    n0 = 3454614
    nplus = 3454627
    na = 4143978
    quotient = 169

    pure_gap = nplus - n0
    vertical_gap = na - nplus
    diagonal_gap = na - n0
    require(
        (pure_gap, vertical_gap, diagonal_gap)
        == (13, 169 * 4079, 689364),
        "THM-2807 address gaps changed",
    )
    require(
        pure_gap + vertical_gap == diagonal_gap,
        "translation triangle stopped closing",
    )

    formal_pure = (1, 0)
    formal_vertical = (0, 4079)
    formal_diagonal = (1, 4079)
    require(
        (
            formal_pure[0] + formal_vertical[0],
            formal_pure[1] + formal_vertical[1],
        )
        == formal_diagonal,
        "formal translation boundary changed",
    )

    # Honest character with chi(Z1)=chi(Z2)=zeta.
    pure_phase = ZETA
    vertical_phase = pow(ZETA, 4079, FIELD)
    diagonal_phase = pow(ZETA, 4080, FIELD)
    boundary_holonomy = (
        pure_phase
        * vertical_phase
        * pow(diagonal_phase, -1, FIELD)
        % FIELD
    )
    require(boundary_holonomy == 1, "honest character gained triangle holonomy")
    require(
        vertical_phase == pow(ZETA, 10, FIELD),
        "vertical residue phase changed",
    )

    # This four-corner completion is an algebraic character control, not a
    # claim that THM-2807 supplies its missing vertical-only physical vertex.
    character_completion = (
        1,
        pure_phase,
        vertical_phase,
        diagonal_phase,
    )
    require(
        cross_ratio(character_completion, FIELD) == 1,
        "honest translation character gained square curvature",
    )
    require(
        mobius(character_completion, FIELD)
        == (1 - pure_phase) * (1 - vertical_phase) % FIELD
        != 0,
        "flat-holonomy character completion lost its additive hostile",
    )

    quotient_addresses = (n0 % quotient, nplus % quotient, na % quotient)
    require(
        quotient_addresses == (85, 98, 98),
        "mod-169 collapse changed",
    )

    # Exhaust all F_53-valued characters of Z/169.  Every one kills the
    # vertical edge because it is the quotient identity.
    quotient_character_generators = tuple(
        value
        for value in range(1, FIELD)
        if pow(value, quotient, FIELD) == 1
    )
    require(
        len(quotient_character_generators) == ORDER,
        "quotient character census changed",
    )
    quotient_vertical_phases = {
        pow(generator, vertical_gap, FIELD)
        for generator in quotient_character_generators
    }
    quotient_pure_phases = tuple(
        pow(generator, pure_gap, FIELD)
        for generator in quotient_character_generators
    )
    quotient_diagonal_phases = tuple(
        pow(generator, diagonal_gap, FIELD)
        for generator in quotient_character_generators
    )
    require(
        quotient_vertical_phases == {1}
        and quotient_pure_phases == quotient_diagonal_phases,
        "a quotient character retained the forgotten vertical coordinate",
    )

    return {
        "addresses": (n0, nplus, na),
        "gaps": (pure_gap, vertical_gap, diagonal_gap),
        "formal_edges": (formal_pure, formal_vertical, formal_diagonal),
        "vertical_residue": 4079 % ORDER,
        "character_phases": (
            pure_phase,
            vertical_phase,
            diagonal_phase,
        ),
        "boundary_holonomy": boundary_holonomy,
        "character_completion": character_completion,
        "completion_kappa": cross_ratio(character_completion, FIELD),
        "completion_mu": mobius(character_completion, FIELD),
        "quotient_addresses": quotient_addresses,
        "quotient_characters": len(quotient_character_generators),
        "quotient_vertical_phases": tuple(sorted(quotient_vertical_phases)),
    }


def main():
    orbit_audits = tuple(
        field_orbit_audit(modulus)
        for modulus in (3, 5, 7, 13)
    )
    branch_scope = branch_and_scope_hostiles()
    provenance = idempotent_and_provenance_audit()
    heisenberg = thm2779_application_audit()
    triangle = thm2807_application_audit()

    print("THM-2814 PROJECTIVE/IDEMPOTENT HOSTILE AUDIT")
    print("status=FINITE-EXACT + independent algebraic proof audit")
    for audit in orbit_audits:
        print(
            "row_column_orbits="
            f"F{audit['field']};"
            f"unit_squares={audit['vectors']};"
            f"orbits={audit['orbits']};"
            f"orbit_size={audit['orbit_size']};"
            f"full_stabilizer_sizes={audit['stabilizer_sizes']};"
            f"kappa1_mu_histogram={audit['kappa1_mu_histogram']};"
            f"fixed_line_nonzero_mu={audit['branch_a_nonzero_mu']}"
        )
    print(
        "branch_A_vs_B="
        f"segre={branch_scope['segre']};"
        f"curved_zero={branch_scope['curved_zero']};"
        "verdict=mu, including its vanishing, is chart-dependent under "
        "row/column gauge; kappa is complete on the unit locus"
    )
    print(
        "scope_hostiles="
        f"arbitrary_vertex_kappas={branch_scope['unrelated_vertex']};"
        "integer_equal_formal_kappa_distinct_unit_orbits="
        f"{branch_scope['integer_equal_formal_kappa_distinct_unit_orbits']};"
        f"z6_denominator={branch_scope['z6_nonzero_denominator']};"
        "repair=state field or four unit entries, and require a "
        "row/column-factorized line gauge"
    )
    print(
        "idempotent_provenance="
        f"F13_idempotents={provenance['field_idempotents']};"
        f"common_faces={provenance['common_faces']};"
        f"linear_checks={provenance['linear_provenance_checks']};"
        f"product_ring_square={provenance['product_ring_square']};"
        f"product_ring_mu={provenance['product_ring_mu']};"
        "verdict=field no-go and linear joint-absent provenance pass; "
        "field scope is essential"
    )
    print(
        "thm2779_application="
        f"field={heisenberg['field']};zeta={heisenberg['zeta']};"
        f"TM_over_MT={heisenberg['oriented_TM_over_MT']};"
        f"reverse={heisenberg['reverse_orientation']};"
        f"normal_form={heisenberg['normal_form']};"
        f"normal_form_mu={heisenberg['normal_form_mu']};"
        f"one_dimensional_solutions={heisenberg['one_dimensional_solutions']};"
        f"projective_basis_lines_commute={heisenberg['basis_line_TM_equals_MT']};"
        "verdict=correct abstract lifted-projective multiplier with stated "
        "orientation; not four physical scalar toggles"
    )
    print(
        "thm2807_application="
        f"addresses={triangle['addresses']};gaps={triangle['gaps']};"
        f"formal_edges={triangle['formal_edges']};"
        f"vertical_residue={triangle['vertical_residue']};"
        f"character_phases={triangle['character_phases']};"
        f"boundary_holonomy={triangle['boundary_holonomy']};"
        f"completion={triangle['character_completion']};"
        f"completion_kappa={triangle['completion_kappa']};"
        f"completion_mu={triangle['completion_mu']};"
        f"quotient_addresses={triangle['quotient_addresses']};"
        f"quotient_characters={triangle['quotient_characters']};"
        f"quotient_vertical_phases={triangle['quotient_vertical_phases']};"
        "verdict=honest full triangle is flat; residue 10 cannot survive "
        "through the mod-169 quotient character"
    )
    print(
        "audit_verdict=ACCEPT_AFTER_SCOPE_REPAIR: the two-branch correction "
        "is mathematically necessary and the THM-2779/2807 applications are "
        "properly bounded.  Promote only on the field/unit locus with a "
        "row/column-factorized gauge, call kappa intrinsic, and call "
        "(kappa-1)w a normalized/covariant defect rather than a "
        "gauge-invariant scalar."
    )


if __name__ == "__main__":
    main()
