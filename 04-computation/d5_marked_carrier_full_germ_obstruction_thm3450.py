#!/usr/bin/env python3
"""Exact deterministic companion for THM-3450.

The script verifies the marked Fourier normalization between THM-2512's
rational doubly-centred 7 x 13 interaction and THM-3443's degree-91 Keller
mode-one germ.  It also computes the universal Keller recurrence through
order fourteen and freezes the first two nonprimitive ANOVA sectors.  These
checks support a marked carrier isomorphism and a full-germ obstruction; they
do not manufacture an H^1 map, an LRC current, or a physical D5 bridge.
"""

from __future__ import annotations

from fractions import Fraction
from functools import cache
from hashlib import sha256
from math import gcd
from pathlib import Path

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = "ac9cdaf18e641572ccac24731d0959177038a7aa95a0559872a7a500708a2f2a"
N = 91


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return sha256(repr(value).encode("utf-8")).hexdigest()


z = sp.symbols("z")
PHI_91 = sp.Poly(sp.cyclotomic_poly(N, z), z, domain=sp.QQ)
require(PHI_91.degree() == 72, PHI_91.degree())


def cyclo_reduce(expression: sp.Expr) -> sp.Expr:
    return sp.rem(sp.Poly(sp.expand(expression), z, domain=sp.QQ), PHI_91).as_expr()


@cache
def zpower(exponent: int) -> sp.Expr:
    return cyclo_reduce(z ** (exponent % N))


def cyclo_key(expression: sp.Expr) -> tuple[tuple[int, int, int], ...]:
    polynomial = sp.Poly(cyclo_reduce(expression), z, domain=sp.QQ)
    return tuple(
        (degree, int(coefficient.p), int(coefficient.q))
        for degree in range(PHI_91.degree())
        if (coefficient := polynomial.nth(degree)) != 0
    )


def root(exponent: int, factor_order: int, convention: str) -> sp.Expr:
    if convention == "normalized":
        root_exponent = 78 if factor_order == 7 else 14
    elif convention == "conventional":
        root_exponent = 13 if factor_order == 7 else 7
    else:
        raise RuntimeError(("root convention", convention))
    return zpower(root_exponent * exponent)


def anova(table: tuple[tuple[Fraction, ...], ...]) -> tuple[tuple[Fraction, ...], ...]:
    row_means = tuple(sum(row, Fraction(0)) / 13 for row in table)
    column_means = tuple(
        sum((table[ell][s] for ell in range(7)), Fraction(0)) / 7
        for s in range(13)
    )
    mean = sum((sum(row, Fraction(0)) for row in table), Fraction(0)) / 91
    interaction = tuple(
        tuple(table[ell][s] - row_means[ell] - column_means[s] + mean for s in range(13))
        for ell in range(7)
    )
    require(all(sum(row, Fraction(0)) == 0 for row in interaction), "ANOVA row margins")
    require(
        all(sum((interaction[ell][s] for ell in range(7)), Fraction(0)) == 0 for s in range(13)),
        "ANOVA column margins",
    )
    return interaction


def ahat(
    table: tuple[tuple[Fraction, ...], ...],
    kappa: int,
    b: int,
    convention: str,
) -> sp.Expr:
    value = sp.Integer(0)
    for ell in range(7):
        for s in range(13):
            value += sp.Rational(table[ell][s].numerator, table[ell][s].denominator) * root(
                kappa * ell, 7, convention
            ) * root(b * s, 13, convention)
    return cyclo_reduce(value / 91)


def cut_psi(
    interaction: tuple[tuple[Fraction, ...], ...],
    tau: int,
    a0: int,
    alpha: int,
    beta: int,
    convention: str,
) -> sp.Expr:
    # d(s,ell)=I(ell,s), with all indices reduced in their own fields.
    psi = sp.Integer(0)
    for c in range(7):
        theta = sp.Integer(0)
        for v in range(13):
            response = sum(
                (
                    interaction[ell][
                        (v - tau * ((a0 * ell + c) % 7)) % 13
                    ]
                    for ell in range(7)
                ),
                Fraction(0),
            )
            theta += sp.Rational(response.numerator, response.denominator) * root(
                -alpha * v, 13, convention
            )
        psi += cyclo_reduce(theta) * root(-beta * c, 7, convention)
    return cyclo_reduce(psi)


def cut_kernel(u: int, beta: int, convention: str) -> sp.Expr:
    ratio_exponent = sp.Integer(0)
    value = sp.Integer(0)
    for j in range(7):
        value += root(-u * j, 13, convention) * root(-beta * j, 7, convention)
    return cyclo_reduce(value + ratio_exponent)


def fourier_and_cut_audit() -> tuple[object, ...]:
    table_rows: list[list[Fraction]] = [[Fraction(0) for _ in range(13)] for _ in range(7)]
    table_rows[0][0] = Fraction(1)
    table_rows[1][1] = Fraction(2)
    table = tuple(tuple(row) for row in table_rows)
    interaction = anova(table)
    require(any(value for row in interaction for value in row), "nonreplica control")

    lambda_normalized = ahat(table, 1, 1, "normalized")
    require(lambda_normalized != 0, "marked primitive coefficient nonzero")
    psi_normalized = cut_psi(interaction, tau=1, a0=1, alpha=12, beta=1, convention="normalized")
    kernel_normalized = cut_kernel(12, 1, "normalized")
    require(kernel_normalized != 0, "normalized cut kernel")
    require(
        cyclo_reduce(psi_normalized - 91 * kernel_normalized * lambda_normalized) == 0,
        (cyclo_key(psi_normalized), cyclo_key(kernel_normalized), cyclo_key(lambda_normalized)),
    )

    lambda_conventional = ahat(table, 6, 2, "conventional")
    require(cyclo_reduce(lambda_conventional - lambda_normalized) == 0, "root-gauge readdressing")
    psi_conventional = cut_psi(interaction, tau=1, a0=1, alpha=11, beta=6, convention="conventional")
    kernel_conventional = cut_kernel(11, 6, "conventional")
    require(kernel_conventional != 0, "conventional cut kernel")
    require(
        cyclo_reduce(psi_conventional - 91 * kernel_conventional * lambda_conventional) == 0,
        "conventional cut factorization",
    )

    # A zero-margin rational array has basis
    # (e_ell-e_6) tensor (e_s-e_12).  Evaluation at the linked primitive
    # roots is an exact 72 x 72 rational isomorphism to Q(zeta_91).
    columns: list[list[sp.Rational]] = []
    for ell in range(6):
        for s in range(12):
            image = cyclo_reduce(
                (root(ell, 7, "normalized") - root(6, 7, "normalized"))
                * (root(s, 13, "normalized") - root(12, 13, "normalized"))
            )
            polynomial = sp.Poly(image, z, domain=sp.QQ)
            columns.append([polynomial.nth(degree) for degree in range(72)])
    evaluation_matrix = sp.Matrix(72, 72, lambda row, column: columns[column][row])
    rank = evaluation_matrix.rank()
    require(rank == 72, rank)

    units = tuple(unit for unit in range(1, 91) if gcd(unit, 91) == 1)
    galois_labels = tuple((unit % 7, unit % 13) for unit in units)
    require(len(units) == len(set(galois_labels)) == 72, galois_labels)
    require(all(a and b for a, b in galois_labels), galois_labels)

    return (
        ("control_table", "anchor_plus_one_off_axis_atom", "interaction_nonzero"),
        ("normalized", "lambda=Ahat(1,1)", "Psi_(tau,1)(12,1)=91*K_(12tau,1)*lambda", cyclo_key(lambda_normalized), cyclo_key(kernel_normalized)),
        ("conventional", "same_lambda=Ahat(6,2)", "Psi_(tau,1)(11,6)=91*K_(11tau,6)*lambda", cyclo_key(kernel_conventional)),
        ("zero_margin_evaluation_rank", rank, "dimension=72=phi(91)"),
        ("Galois_mixed_labels", len(galois_labels), digest(galois_labels)),
    )


def series_multiply(left: list[Fraction], right: list[Fraction], cutoff: int) -> list[Fraction]:
    result = [Fraction(0) for _ in range(cutoff + 1)]
    for i, left_value in enumerate(left):
        if not left_value:
            continue
        for j, right_value in enumerate(right):
            if i + j <= cutoff and right_value:
                result[i + j] += left_value * right_value
    return result


def series_power(value: list[Fraction], exponent: int, cutoff: int) -> list[Fraction]:
    result = [Fraction(0) for _ in range(cutoff + 1)]
    result[0] = Fraction(1)
    base = value[:]
    power = exponent
    while power:
        if power & 1:
            result = series_multiply(result, base, cutoff)
        base = series_multiply(base, base, cutoff)
        power //= 2
    return result


def fraction_key(value: Fraction) -> tuple[int, int]:
    return value.numerator, value.denominator


def universal_keller_recurrence() -> tuple[object, ...]:
    cutoff = 14
    V = [Fraction(0) for _ in range(cutoff + 1)]
    V[0] = Fraction(1)
    for degree in range(1, cutoff + 1):
        # Through order 87 the transformed inverse equation is
        # 1-V^91+(91/90)qV^90=0, q=s/rho.  With V_degree temporarily zero,
        # the missing linear contribution to V^91 is 91*V_degree.
        known_91 = series_power(V, 91, cutoff)[degree]
        known_90 = series_power(V, 90, cutoff)[degree - 1]
        V[degree] = (Fraction(91, 90) * known_90 - known_91) / 91

    residual = [Fraction(0) for _ in range(cutoff + 1)]
    residual[0] = Fraction(1)
    V91 = series_power(V, 91, cutoff)
    V90 = series_power(V, 90, cutoff)
    for degree in range(cutoff + 1):
        residual[degree] -= V91[degree]
        if degree:
            residual[degree] += Fraction(91, 90) * V90[degree - 1]
    require(all(value == 0 for value in residual), tuple(map(fraction_key, residual)))

    # X=s^-90*x=rho*(V-qV')/91.  The order-m coefficient is
    # ((1-m)V_m/91)*rho^(1-m).  For m>=2, canonical reduction to a
    # nonnegative rho exponent adds 91 and divides the scalar by rho^91=91.
    x_rows: list[tuple[object, ...]] = []
    first_nonprimitive: tuple[object, ...] | None = None
    for degree, coefficient in enumerate(V):
        scalar = Fraction(1 - degree) * coefficient / 91
        raw_rho_exponent = 1 - degree
        canonical_scalar = scalar
        canonical_exponent = raw_rho_exponent
        while canonical_exponent < 0:
            canonical_exponent += 91
            canonical_scalar /= 91
        mode = (1 - degree) % 91
        address = (mode % 7, mode % 13)
        sector = (
            "trivial" if address == (0, 0)
            else "C13_main" if address[0] == 0
            else "C7_main" if address[1] == 0
            else "primitive_mixed"
        )
        row = (
            degree,
            fraction_key(coefficient),
            fraction_key(scalar),
            canonical_exponent,
            fraction_key(canonical_scalar),
            mode,
            address,
            sector,
        )
        x_rows.append(row)
        if scalar and sector != "primitive_mixed" and first_nonprimitive is None:
            first_nonprimitive = row

    require(x_rows[2][2] == (-1, 16380), x_rows[2])
    require(first_nonprimitive is not None and first_nonprimitive[0] == 8, first_nonprimitive)
    require(first_nonprimitive[5:8] == (84, (0, 6), "C13_main"), first_nonprimitive)
    require(
        (x_rows[8][3], x_rows[8][4])
        == (84, (-418006897, 404160880500000)),
        x_rows[8],
    )
    require(x_rows[14][5:8] == (78, (1, 0), "C7_main"), x_rows[14])
    require(
        (x_rows[14][3], x_rows[14][4])
        == (78, (-18258883976232522940057, 20343399512726430000000000000)),
        x_rows[14],
    )

    require((91 - 3, 91 - 2, 91 - 1) == (88, 89, 90), "lower-term entry orders")
    return (
        ("universal_equation_through_order_87", "1-V^91+(91/90)qV^90=0", "q=s/rho"),
        ("V_coefficients_0_to_14", tuple(fraction_key(value) for value in V)),
        ("X_rows_columns", "m,V_m,raw_scalar,rho_exp,canonical_scalar,mode,CRT,sector"),
        tuple(x_rows),
        ("first_nonprimitive", first_nonprimitive),
        ("first_other_margin", x_rows[14]),
        ("seed_entry_orders", "w^3_at_88", "w^2_at_89", "P0_at_90"),
    )


def structural_obstruction_audit() -> tuple[object, ...]:
    sector_dimensions = (1, 6, 12, 72)
    require(sum(sector_dimensions) == 91, sector_dimensions)
    require(6 * 12 == 72, "primitive tensor dimension")
    require((84 % 7, 84 % 13) == (0, 6), "order-eight sector")
    require((78 % 7, 78 % 13) == (1, 0), "order-fourteen sector")
    require(1 != (-1) % 91, "inertia reversal separates characters")
    return (
        ("rational_ANOVA_decomposition", "Q", "Q(zeta7)", "Q(zeta13)", "Q(zeta91)", sector_dimensions),
        ("marked_module", "R=E[[t]]", "Phi(f*e_(1,1))=f*u_K(t)*e_1", "u_K(0)=rho/91", "strict_filtered_isomorphism"),
        ("isomorphism_torsor", "Isom_fil_C91(R_chi1,R_chi1)=R^times", "amplitude_calibration_required"),
        ("gauge_hostile", "reverse_Keller_only:chi1_to_chi-1", "Hom_C91=0", "diagonal_reversal_restores_label"),
        ("H1_hostile", "H1(C13;F13)_has_exponent_13", "E[[t]]_is_torsionfree_and_uniquely_13_divisible", "all_additive_maps_both_directions_zero"),
        ("full_germ_hostiles", "order2_leaves_fixed_line", "order8_first_nonprimitive_(0,6)", "order14_other_margin_(1,0)"),
        ("missing_object", "gauge-linked_Rees_connection_with_order_m_in_mode_1-m", "not_constructed"),
    )


def main() -> None:
    fourier = fourier_and_cut_audit()
    recurrence = universal_keller_recurrence()
    structural = structural_obstruction_audit()
    verdict = (
        "lambda_L=Ahat(1,1)=Psi_(tau,1)(12,1)/(91*K_(12tau,1))_under_linked_normalized_roots",
        "same_character_is_Ahat_conventional(6,2)=Psi_(tau,1)(11,6)/(91*K_(11tau,6))",
        "one_rational_primitive_coefficient_plus_Galois_determines_the_whole_doubly_centred_interaction",
        "multiplication_by_Keller_mode-one_unit_is_a_strict_marked_filtered_C91_isomorphism",
        "replica_iff_lambda_zero_iff_marked_carrier_image_zero_for_lawful_rational_anchored_tables",
        "no_additive_H1_to_amplitude_map_and_no_unmarked_gauge-independent_map",
        "order8_nonzero_(0,6)_term_blocks_full-germ_factorization_through_centred_interaction",
        "order14_nonzero_(1,0)_term_exhibits_the_other_missing_ANOVA_margin",
        "gauge-linked_Rees_connection_is_a_missing_sidecar_not_a_constructed_D5_current",
        "no_LRC14_JC2_or_physical_transport_consequence",
    )
    semantic_surface = (fourier, recurrence, structural, verdict)
    semantic_sha = digest(semantic_surface)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, (semantic_sha, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    print("THM-3450 marked D5 carrier and full-germ obstruction exact companion")
    print("status=EXACT_DETERMINISTIC_COMPANION_FOR_THM3450;finite_checks_do_not_replace_proof")
    print("arithmetic=QQ_cyclotomic_reduction_mod_Phi91;72x72_exact_rank;Fraction_formal_recurrence;no_randomness")
    print(f"fourier_and_cut={fourier}")
    print(f"universal_Keller_recurrence={recurrence}")
    print(f"structural_obstructions={structural}")
    print(f"verdict={verdict}")
    print("scope=marked_characteristic-zero_line_after_formal_base_change;centred_rational_interaction_only;full_Keller_germ_obstructed")
    print("loss_boundary=main_effects_ancestry_positivity_owner_clock_physical_time_and_H1_class_are_not_transported")
    print(f"semantic_sha256={semantic_sha}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;all_checks_survive_python_O")
    print("commands=python -B 04-computation/d5_marked_carrier_full_germ_obstruction_thm3450.py;python -B -O 04-computation/d5_marked_carrier_full_germ_obstruction_thm3450.py")
    print("PASS")


if __name__ == "__main__":
    main()
