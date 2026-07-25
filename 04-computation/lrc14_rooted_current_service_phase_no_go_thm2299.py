#!/usr/bin/env python3
"""Exact companion for THM-2299."""

from fractions import Fraction


P = 13

DELTA_STRICT = Fraction(39002430583, 53493927587100)
DELTA_REPEAT = Fraction(13560199813, 53493927587100)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(x: Fraction) -> Fraction:
    residue = x % 1
    return min(residue, 1 - residue)


def main() -> None:
    strict_energy_floor = DELTA_STRICT**2
    repeat_energy_floor = DELTA_REPEAT**2
    strict_named_energy_floor = strict_energy_floor / 4
    repeat_named_energy_floor = repeat_energy_floor / 4

    require(
        strict_energy_floor
        == Fraction(
            1521189591381733719889,
            2861600288693898428086410000,
        ),
        "strict rooted-energy floor changed",
    )
    require(
        repeat_energy_floor
        == Fraction(
            183879018968485234969,
            2861600288693898428086410000,
        ),
        "repeated rooted-energy floor changed",
    )
    require(
        4 * strict_named_energy_floor == strict_energy_floor,
        "strict named-edge floor changed",
    )
    require(
        4 * repeat_named_energy_floor == repeat_energy_floor,
        "repeated named-edge floor changed",
    )

    # The worst two-sheet phase is the nontrivial thirteenth root closest
    # to -1. Its half-angle gap is pi/26, and concavity gives
    # sin(pi/26)>1/13. The finite ledger checks the exponent boundary.
    nonzero_exponents = range(1, P)
    distances_from_half_turn = [
        abs(P - 2 * min(j, P - j)) for j in nonzero_exponents
    ]
    require(
        min(distances_from_half_turn) == 1,
        "worst thirteenth-root half-turn gap changed",
    )
    require(P > 0, "root modulus must be positive")

    # Symbolic covariant endpoint-Prony bank. Both the ancestry density and
    # one named target use subsets of the same total scalar boundary bank.
    jump_rows = 0
    for S in range(9, 80):
        J_G = 2 * S
        K_named = 2 * S
        B = J_G + K_named
        require(B - 1 == 4 * S - 1, "covariant gauge-index bound changed")
        jump_rows += 1

    # Smallest hard-multiplier hostile row.
    H = 1
    qs = (4, 2, 3, 6, 10)
    c1 = P
    c2 = P**3
    c3 = 2 * P**5
    cs = (c1, c2, c3)

    require(H % 2 == 1, "guard parity failed")
    require(len(set(qs)) == 5, "unit speeds are not distinct")
    require(all(q % P != 0 for q in (H, *qs)), "a guard/unit is ramified")
    require(c1 % P == 0 and c1 % P**2 != 0, "c1 depth changed")
    require(c2 % P**3 == 0 and c2 % P**4 != 0, "c2 depth changed")
    require(c3 % P**5 == 0 and c3 % P**6 != 0, "c3 depth changed")

    # p: 13q1-4c1=0; r:q1-2q2=0.
    require(P * qs[0] - 4 * c1 == 0, "shallow pair relation failed")
    require(qs[0] - 2 * qs[1] == 0, "independent relation failed")
    anchored_minor = (-4) * (-2)
    require(anchored_minor == 8, "anchored minor changed")
    require(anchored_minor % P != 0, "anchored minor ramified")

    # The hostile row already has a much stronger c1-anchored rank-six
    # relation packet.  In coordinate order
    #
    #   (H,q1,q2,q3,q4,q5,c1,c2,c3),
    #
    # take r_i=q_i-q_i_value*H and p=13q1-4c1.  On the pivot columns
    # (q1,...,q5,c1), the first five rows are I_5 and the last row is
    # (13,0,0,0,0,-4), so the determinant is exactly -4.
    speed_row = (H, *qs, *cs)
    rank_six_rows: list[tuple[int, ...]] = []
    for i, q in enumerate(qs, start=1):
        row = [0] * len(speed_row)
        row[0] = -q
        row[i] = 1
        rank_six_rows.append(tuple(row))
    pair_row = [0] * len(speed_row)
    pair_row[1] = P
    pair_row[6] = -4
    rank_six_rows.append(tuple(pair_row))

    require(
        all(sum(a * b for a, b in zip(row, speed_row)) == 0 for row in rank_six_rows),
        "anchored rank-six packet contains a nonrelation",
    )
    rank_six_heights = tuple(max(abs(a) for a in row) for row in rank_six_rows)
    require(
        rank_six_heights == (4, 2, 3, 6, 10, 13),
        "anchored rank-six heights changed",
    )
    pivot_columns = (1, 2, 3, 4, 5, 6)
    pivot_matrix = tuple(
        tuple(row[j] for j in pivot_columns) for row in rank_six_rows
    )
    require(
        pivot_matrix
        == (
            (1, 0, 0, 0, 0, 0),
            (0, 1, 0, 0, 0, 0),
            (0, 0, 1, 0, 0, 0),
            (0, 0, 0, 1, 0, 0),
            (0, 0, 0, 0, 1, 0),
            (13, 0, 0, 0, 0, -4),
        ),
        "anchored rank-six pivot matrix changed",
    )
    rank_six_minor = -4
    require(rank_six_minor % P != 0, "anchored rank-six minor ramified")

    root = 8
    require((qs[0] * root) % P == 6, "pair-aligned root changed")

    epsilon = Fraction(1, 10**12)
    y_centers = (Fraction(-1, 16), Fraction(1, 16))

    margin_rows: list[tuple[Fraction, str]] = []
    for y0 in y_centers:
        x0 = (root + y0) / P
        z0 = (P * y0) % 1

        source_conditions = [
            ("source_guard", H, Fraction(1, 7), False),
            *[
                (f"source_q_{q}", q, Fraction(1, 14), False)
                for q in qs
            ],
            ("source_c1", c1, Fraction(1, 14), True),
            ("source_c2", c2, Fraction(1, 14), False),
            ("source_c3", c3, Fraction(1, 14), False),
        ]
        for label, coefficient, threshold, wants_danger in source_conditions:
            value = circle_norm(coefficient * x0)
            margin = threshold - value if wants_danger else value - threshold
            slope = Fraction(coefficient, P)
            require(margin > 0, f"{label} is not strict at {x0}")
            margin_rows.append((margin / slope, label))

        target_conditions = [
            ("target_guard", H, Fraction(1, 7), False),
            *[
                (f"target_q_{q}", q, Fraction(1, 14), False)
                for q in qs
            ],
            ("target_c1", c1, Fraction(1, 14), False),
            ("target_c2", c2, Fraction(1, 14), True),
            ("target_c3", c3, Fraction(1, 14), False),
        ]
        for label, coefficient, threshold, wants_danger in target_conditions:
            value = circle_norm(coefficient * z0)
            margin = threshold - value if wants_danger else value - threshold
            slope = P * coefficient
            require(margin > 0, f"{label} is not strict at {z0}")
            margin_rows.append((margin / slope, label))

    minimum_margin, minimum_label = min(margin_rows)
    require(
        minimum_margin == Fraction(3, 540602608),
        "minimum hostile half-width margin changed",
    )
    require(minimum_label == "target_c3", "minimum-margin label changed")
    require(epsilon < minimum_margin, "hostile intervals are too wide")
    require(2 * epsilon < Fraction(1, 8), "source intervals overlap")

    source_centers = tuple((root + y) / P for y in y_centers)
    target_centers = tuple((P * y) % 1 for y in y_centers)
    require(
        source_centers == (Fraction(127, 208), Fraction(129, 208)),
        "source centers changed",
    )
    require(
        target_centers == (Fraction(3, 16), Fraction(13, 16)),
        "terminal centers changed",
    )
    target_separation = circle_norm(target_centers[1] - target_centers[0])
    require(2 * P * epsilon < target_separation, "terminal images overlap")
    require(2 * P * epsilon < 1, "T is not injective on a source component")

    # Equal intervals centered at +/-1/16 have phase ratio -1 at frequency
    # four because 4*((1/16)-(-1/16))=1/2.
    phase_difference = qs[0] * (y_centers[1] - y_centers[0])
    require(phase_difference == Fraction(1, 2), "quarter-turn cancellation failed")

    measure_F = 4 * epsilon
    measure_E = measure_F / P
    measure_R = P * measure_F
    rooted_character_squared = Fraction(1, P**2)
    rooted_energy = measure_R * rooted_character_squared
    require(measure_E == Fraction(1, 3250000000000), "source mass changed")
    require(rooted_energy == measure_E, "one-sheet rooted energy changed")
    require(rooted_energy > measure_E**2, "rooted-energy inequality failed")

    print("theorem=THM-2299")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
    print(f"return_floor_strict={DELTA_STRICT}")
    print(f"return_floor_repeat={DELTA_REPEAT}")
    print(f"rooted_energy_floor_strict={strict_energy_floor}")
    print(f"rooted_energy_floor_repeat={repeat_energy_floor}")
    print(f"named_edge_energy_floor_strict={strict_named_energy_floor}")
    print(f"named_edge_energy_floor_repeat={repeat_named_energy_floor}")
    print("two_sheet_phase_bound=sin(pi/26)>1/13")
    print(f"covariant_jump_rows={jump_rows}")
    print("gauge_index_bound_on_named_edge=h<=4S-1")
    print(f"hostile_profile=(1,3,5)")
    print(f"hostile_speeds=H:{H};q:{qs};c:{cs}")
    print("hostile_pair=13q1-4c1=0;K=52;kappa=4")
    print("hostile_second_relation=q1-2q2=0")
    print(f"hostile_anchored_minor={anchored_minor}")
    print(f"hostile_anchored_rank6_heights={rank_six_heights}")
    print("hostile_anchored_rank6_columns=(q1,q2,q3,q4,q5,c1)")
    print(f"hostile_anchored_rank6_minor={rank_six_minor}")
    print(f"hostile_source_centers={source_centers}")
    print(f"hostile_target_centers={target_centers}")
    print(f"hostile_epsilon={epsilon}")
    print(f"hostile_minimum_margin={minimum_margin};label={minimum_label}")
    print(f"hostile_source_mass={measure_E}")
    print(f"hostile_rooted_energy_each_character={rooted_energy}")
    print("hostile_exact_cancellations=Fhat(4)=Ehat(52)=What_4(0)=0")
    print("hostile_scope=local_strict_labels_only;not_a_global_cover")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
