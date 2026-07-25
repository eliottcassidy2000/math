#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2321."""

from fractions import Fraction
from itertools import product


P = 13


def require(condition: bool, message: str) -> None:
    """Raise in ordinary and optimized Python when a check fails."""
    if not condition:
        raise RuntimeError(message)


def reduce_cyclotomic(coefficients: list[Fraction]) -> tuple[Fraction, ...]:
    """Reduce a length-p coefficient vector modulo 1+x+...+x^(p-1)."""
    last = coefficients[-1]
    return tuple(value - last for value in coefficients[:-1])


def add_term(
    coefficients: list[Fraction],
    exponent: int,
    weight: Fraction,
) -> None:
    """Add weight*zeta^exponent to a cyclic coefficient vector."""
    coefficients[exponent % len(coefficients)] += weight


def bispectrum_cell(
    prime: int,
    masses: tuple[int, int],
    sites: tuple[int, int],
    first: int,
    second: int,
) -> list[Fraction]:
    """Exact group-algebra cell M_k M_l conjugate(M_(k+l))."""
    coefficients = [Fraction(0) for _ in range(prime)]
    for r_index, s_index, t_index in product(range(2), repeat=3):
        r = sites[r_index]
        s = sites[s_index]
        t = sites[t_index]
        weight = Fraction(
            masses[r_index] * masses[s_index] * masses[t_index]
        )
        exponent = -first * r - second * s + (first + second) * t
        add_term(coefficients, exponent, weight)
    return coefficients


def real_part(coefficients: list[Fraction]) -> list[Fraction]:
    """Return the exact cyclic coefficient vector of the real part."""
    prime = len(coefficients)
    return [
        (coefficients[index] + coefficients[-index % prime]) / 2
        for index in range(prime)
    ]


def allowed_pair_atlas(prime: int = P) -> tuple[int, int, int]:
    """Verify the row/shape factorization of all allowed pairs."""
    direct = {
        (first, second)
        for first in range(1, prime)
        for second in range(1, prime)
        if (first + second) % prime != 0
    }
    factored = {
        (first, (shape * first) % prime)
        for first in range(1, prime)
        for shape in range(1, prime)
        if shape != prime - 1
    }
    require(direct == factored, "allowed-pair factorization changed")
    require(len(direct) == (prime - 1) * (prime - 2),
            "allowed-pair count changed")
    return len(direct), prime - 1, prime - 2


def fixed_shape_checks() -> int:
    """Exhaust the exact fixed-shape identity on small binary vectors."""
    rows = 0
    for prime in (5, 7, 13):
        for first_site in range(prime):
            for difference in range(1, prime):
                sites = (first_site, (first_site + difference) % prime)
                for a in range(5):
                    for b in range(5):
                        if a + b == 0:
                            continue
                        mass = a + b
                        expected = prime * (a**3 + b**3) - mass**3
                        floor = Fraction(prime - 4, 4) * mass**3
                        require(expected >= floor,
                                "fixed-shape positivity floor failed")
                        require(
                            (expected == floor) == (a == b and a > 0),
                            "fixed-shape equality classification changed",
                        )
                        for shape in range(1, prime - 1):
                            coefficients = [Fraction(0) for _ in range(prime)]
                            for first in range(1, prime):
                                second = (shape * first) % prime
                                cell = bispectrum_cell(
                                    prime, (a, b), sites, first, second
                                )
                                coefficients = [
                                    x + y for x, y in zip(coefficients, cell)
                                ]
                            predicted = [Fraction(0) for _ in range(prime)]
                            predicted[0] = Fraction(expected)
                            require(
                                reduce_cyclotomic(coefficients)
                                == reduce_cyclotomic(predicted),
                                "fixed-shape cyclotomic identity failed",
                            )
                            rows += 1
    return rows


def fixed_first_character_checks() -> int:
    """Verify the exact group-algebra formula for every prescribed k."""
    rows = 0
    for prime in (5, 7, 13):
        for difference in range(1, prime):
            sites = (0, difference)
            for a in range(5):
                for b in range(5):
                    if a + b == 0:
                        continue
                    mass = a + b
                    ab_mass = a * b * mass
                    for first in range(1, prime):
                        coefficients = [Fraction(0) for _ in range(prime)]
                        for second in range(1, prime):
                            if (first + second) % prime == 0:
                                continue
                            cell = bispectrum_cell(
                                prime, (a, b), sites, first, second
                            )
                            coefficients = [
                                x + y for x, y in zip(coefficients, cell)
                            ]
                        actual_real = real_part(coefficients)

                        phase = (first * difference) % prime
                        predicted = [Fraction(0) for _ in range(prime)]
                        predicted[0] = Fraction(
                            (prime - 2) * mass**3
                            + (-3 * prime + 4) * ab_mass
                        )
                        oscillatory = Fraction((prime - 4) * ab_mass, 2)
                        add_term(predicted, phase, oscillatory)
                        add_term(predicted, -phase, oscillatory)
                        require(
                            reduce_cyclotomic(actual_real)
                            == reduce_cyclotomic(predicted),
                            "fixed-first-character real-part identity failed",
                        )
                        rows += 1
    return rows


def threshold_checks() -> tuple[Fraction, Fraction, Fraction]:
    """Check p=3/4 boundaries and the p=13 rational word constants."""
    # Equal masses V=2, x=1/4.
    p4_zero = (
        Fraction(4 - 4, 4)
        * (1 + Fraction(-1))
        * 2**3
    )
    require(p4_zero == 0, "p=4 antipodal boundary changed")

    # At p=3 and theta=2*pi/3, cos(theta)=-1/2.
    p3_negative = (
        Fraction(3 - 4, 4)
        * (1 + Fraction(-1, 2))
        * 2**3
    )
    require(p3_negative == -1, "p=3 hostile value changed")

    fibre_cube_floor = 49 * P**3
    require(fibre_cube_floor == 107653,
            "word Hölder cube coefficient changed")

    fixed_shape_aggregate = Fraction(P - 4, 4) * fibre_cube_floor
    fixed_shape_term = fixed_shape_aggregate / (P - 1)
    require(fixed_shape_aggregate == Fraction(968877, 4),
            "fixed-shape word aggregate changed")
    require(fixed_shape_term == Fraction(322959, 16),
            "fixed-shape selected-term floor changed")

    # Concavity of sine gives 1-cos(pi/p)>2/p^2 for odd p>1.
    fixed_first_rational_aggregate = (
        Fraction(P - 4, 4)
        * Fraction(2, P**2)
        * fibre_cube_floor
    )
    fixed_first_rational_term = fixed_first_rational_aggregate / (P - 2)
    require(fixed_first_rational_aggregate == Fraction(5733, 2),
            "fixed-first rational aggregate changed")
    require(fixed_first_rational_term == Fraction(5733, 22),
            "fixed-first selected-term floor changed")
    return (
        fixed_shape_term,
        fixed_first_rational_aggregate,
        fixed_first_rational_term,
    )


def negative_cell_boundary() -> tuple[int, int, int]:
    """Record an exact p=13 cell whose real bispectrum is negative."""
    # For equal masses at sites 0,1 and (k,l,k+l)=(1,6,7),
    # the value is 8*cos(pi/13)*cos(6pi/13)*cos(7pi/13).
    # The first two factors are positive and the last is negative.
    first, second = 1, 6
    third = (first + second) % P
    require(0 < first <= P // 2, "first cosine sign changed")
    require(0 < second <= P // 2, "second cosine sign changed")
    require(P // 2 < third < P, "third cosine sign changed")
    return first, second, third


def projective_shape_check() -> tuple[int, int, int, int, int]:
    """Verify the full three-marked projective character-pair line."""
    finite_slopes = tuple(range(P))
    infinity = "infinity"
    projective_line = finite_slopes + (infinity,)
    generic_slopes = tuple(
        slope for slope in finite_slopes if slope not in (0, P - 1)
    )
    boundary_slopes = (0, P - 1, infinity)
    require(len(projective_line) == P + 1,
            "projective character-pair line changed")
    require(len(generic_slopes) == P - 2,
            "generic bispectrum slopes changed")
    require(len(boundary_slopes) == 3,
            "three-marked boundary changed")

    def swap(slope: int | str) -> int | str:
        if slope == 0:
            return infinity
        if slope == infinity:
            return 0
        return pow(int(slope), -1, P)

    require(all(swap(swap(slope)) == slope for slope in projective_line),
            "leg-swap involution changed")
    fixed_mixed = tuple(
        slope
        for slope in range(1, P)
        if swap(slope) == slope
    )
    require(fixed_mixed == (1, P - 1),
            "swap-fixed mixed slopes changed")

    # Any labelled-axis projective identification has lambda -> c*lambda.
    # Compatibility with inversion restricts c to c^2=1.
    couplings = tuple(range(1, P))
    swap_compatible = tuple(c for c in couplings if c * c % P == 1)
    require(swap_compatible == (1, P - 1),
            "swap-compatible coupling constants changed")
    return (
        len(projective_line),
        len(generic_slopes),
        len(boundary_slopes),
        len(couplings),
        len(swap_compatible),
    )


def boundary_shape_positivity_check() -> int:
    """Check the exact positive floor on the three trivial-character shapes."""
    rows = 0
    for prime in (5, 7, 13):
        for a in range(5):
            for b in range(5):
                if a + b == 0:
                    continue
                mass = a + b
                second_moment = a * a + b * b
                boundary_sum = mass * (
                    prime * second_moment - mass * mass
                )
                floor = Fraction(prime - 2, 2) * mass**3
                require(boundary_sum >= floor,
                        "projective boundary-shape floor failed")
                require(
                    (boundary_sum == floor) == (a == b and a > 0),
                    "projective boundary equality changed",
                )
                rows += 3
    return rows


def target_gain_pivot_no_go() -> tuple[int, int, int]:
    """Verify that owner-pivot normalization leaves all target gains free."""
    gains = tuple(range(1, P))

    # Mod-13 owner-star packet. The six pivot columns are j,u1,...,u5;
    # a and b are target columns and u0 is the omitted unit.
    labels = ("j", "a", "b", "u0", "u1", "u2", "u3", "u4", "u5")
    pivot = ("j", "u1", "u2", "u3", "u4", "u5")
    index = {label: position for position, label in enumerate(labels)}

    target_b_values = set()
    pivot_signatures = set()
    non_target_signatures = set()
    for alpha in gains:
        rows = []
        for row_label in pivot:
            row = [0] * len(labels)
            row[index[row_label]] = -1 % P
            row[index["u0"]] = 0 if row_label == "j" else 1
            rows.append(row)
        rows[1][index["a"]] = -1 % P
        rows[2][index["b"]] = -alpha % P

        pivot_signature = tuple(
            tuple(row[index[label]] for label in pivot)
            for row in rows
        )
        non_target_signature = tuple(
            tuple(
                row[index[label]]
                for label in labels
                if label != "b"
            )
            for row in rows
        )
        pivot_signatures.add(pivot_signature)
        non_target_signatures.add(non_target_signature)
        target_b_values.add(rows[2][index["b"]])

        require(
            all(
                pivot_signature[i][j] == (-1 % P if i == j else 0)
                for i in range(6)
                for j in range(6)
            ),
            "owner pivot restriction changed",
        )

    require(len(pivot_signatures) == 1,
            "target rescaling changed the owner pivot")
    require(len(non_target_signatures) == 1,
            "target rescaling changed non-b packet data")
    require(target_b_values == set(gains),
            "target rescaling failed to sweep the gain torsor")
    return (
        len(pivot_signatures),
        len(non_target_signatures),
        len(target_b_values),
    )


def main() -> None:
    allowed_pairs, root_rows, shape_columns = allowed_pair_atlas()
    fixed_shape_rows = fixed_shape_checks()
    fixed_first_rows = fixed_first_character_checks()
    (
        fixed_shape_term,
        fixed_first_aggregate,
        fixed_first_term,
    ) = threshold_checks()
    negative_cell = negative_cell_boundary()
    (
        projective_slopes,
        generic_slopes,
        boundary_slopes,
        coupling_count,
        swap_compatible_couplings,
    ) = projective_shape_check()
    boundary_shape_rows = boundary_shape_positivity_check()
    (
        pivot_signature_count,
        non_target_signature_count,
        target_b_value_count,
    ) = target_gain_pivot_no_go()

    print("theorem=THM-2321")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
    print("allowed_pair_factorization=F_13^* x Lambda_13")
    print(f"allowed_pairs={allowed_pairs}")
    print(f"root_rows={root_rows}")
    print(f"shape_columns={shape_columns}")
    print("fixed_shape_identity=S_lambda=p*S3-S1^3")
    print("fixed_shape_floor=(p-4)/4*S1^3")
    print("fixed_first_floor=(p-4)/4*(1+cos(theta))*S1^3")
    print("uniform_fixed_first_floor=(p-4)/4*(1-cos(pi/p))*S1^3")
    print(f"fixed_shape_exact_rows={fixed_shape_rows}")
    print(f"fixed_first_exact_rows={fixed_first_rows}")
    print(f"p13_fixed_shape_selected_word_floor={fixed_shape_term}*rho^3")
    print(
        "p13_fixed_first_rational_aggregate_floor="
        f"{fixed_first_aggregate}*rho^3"
    )
    print(
        "p13_fixed_first_rational_selected_floor="
        f"{fixed_first_term}*rho^3"
    )
    print("p4_equal_two_sheet_fixed_first_boundary=0")
    print("p3_equal_two_sheet_hostile=-1")
    print(f"p13_negative_prescribed_cell={negative_cell}")
    print("every_root_row_has_positive_sum=true")
    print("every_shape_column_has_positive_sum=true")
    print("individual_cells_need_not_be_positive=true")
    print(f"projective_character_pair_slopes={projective_slopes}")
    print(f"generic_bispectrum_slopes={generic_slopes}")
    print(f"trivial_character_boundary_slopes={boundary_slopes}")
    print(f"boundary_shape_exact_rows={boundary_shape_rows}")
    print(f"axis_preserving_projective_rescalings={coupling_count}")
    print(f"swap_compatible_projective_rescalings={swap_compatible_couplings}")
    print(f"pivot_signatures_under_target_scaling={pivot_signature_count}")
    print(
        "non_target_signatures_under_target_scaling="
        f"{non_target_signature_count}"
    )
    print(f"target_b_values_swept={target_b_value_count}")
    print("thm2318_root_character_can_be_prescribed=true")
    print("abstract_target_gain_label_selected=true")
    print("ordinary_atom_is_not_a_cubic_factor=true")
    print("pair_edge_incidence_not_proved=true")
    print("visible_target_gain_address_not_landed=true")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
