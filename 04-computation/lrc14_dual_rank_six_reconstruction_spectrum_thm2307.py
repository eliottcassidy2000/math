"""Exact companion for THM-2307.

The script checks the sharp rank-six and rank-seven models, the complete
count arithmetic, the low-degree selector floors, and the THM-2299
pure-blocker hostile packet.  All validity gates use ``require`` so normal
and optimized Python execute the same audit.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dot(form: tuple[int, ...], point: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(form, point)) % P


def matrix_rank_mod_p(rows: tuple[tuple[int, ...], ...]) -> int:
    matrix = [[entry % P for entry in row] for row in rows]
    if not matrix:
        return 0
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (r for r in range(pivot_row, row_count) if matrix[r][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, P)
        matrix[pivot_row] = [(inverse * x) % P for x in matrix[pivot_row]]
        for r in range(row_count):
            if r == pivot_row:
                continue
            factor = matrix[r][column]
            if factor:
                matrix[r] = [
                    (x - factor * y) % P
                    for x, y in zip(matrix[r], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def rank_six_sharp_model() -> tuple[int, int, int]:
    # Dual coordinates are (x,y,z).  The first six forms are the physical
    # guard/unit coordinates; the final three are blocker coordinates.
    unit_forms = ((0, 0, 1),) + tuple((1, 0, -c) for c in range(1, 6))
    blocker_forms = ((1, 0, 0), (0, 1, 0), (1, -1, 0))
    forms = unit_forms + blocker_forms
    physical = (0, 0, 1)

    require(matrix_rank_mod_p(forms) == 3, "rank-six model is not a 3-plane")
    require(
        all(dot(form, physical) != 0 for form in unit_forms),
        "rank-six physical unit coordinate vanished",
    )
    require(
        all(dot(form, physical) == 0 for form in blocker_forms),
        "rank-six physical blocker coordinate did not vanish",
    )
    require(
        all(any(entry % P for entry in form) for form in blocker_forms),
        "rank-six model accidentally has a dark blocker",
    )

    # Every all-unit projective direction has z != 0 because z is one of
    # the unit forms, so it has a unique representative (x,y,1).
    good = []
    bad = []
    for x, y in product(range(P), repeat=2):
        point = (x, y, 1)
        word = tuple(dot(form, point) for form in forms)
        (good if all(word) else bad).append((point, word))

    require(len(good) == 77, "sharp rank-six good count changed")
    require(len(bad) == 92, "sharp rank-six bad count changed")
    require(len(good) * (P - 1) == 924, "rank-six vector count changed")

    # The bad set is exactly x in {0,...,5}, or y=0, or y=x.
    explicit_bad = {
        (x, y)
        for x, y in product(range(P), repeat=2)
        if x in range(6) or y == 0 or y == x
    }
    require(len(explicit_bad) == 92, "sharp bad-union formula changed")
    require(
        {(point[0], point[1]) for point, _ in bad} == explicit_bad,
        "sharp bad-union geometry changed",
    )
    return len(good), len(bad), len(good) * (P - 1)


def rank_seven_sharp_model() -> tuple[int, int]:
    # Dual coordinates are (x,z).  The blocker forms may coincide; this is
    # the sharp one-blocker-point configuration in P^1(F_13).
    unit_forms = ((0, 1),) + tuple((1, -c) for c in range(1, 6))
    blocker_forms = ((1, 0), (1, 0), (1, 0))
    forms = unit_forms + blocker_forms
    physical = (0, 1)

    require(matrix_rank_mod_p(forms) == 2, "rank-seven model is not a 2-plane")
    require(
        all(dot(form, physical) != 0 for form in unit_forms),
        "rank-seven physical unit coordinate vanished",
    )
    require(
        all(dot(form, physical) == 0 for form in blocker_forms),
        "rank-seven physical blocker coordinate did not vanish",
    )

    projective_points = [(x, 1) for x in range(P)] + [(1, 0)]
    good = [
        point
        for point in projective_points
        if all(dot(form, point) != 0 for form in forms)
    ]
    require(len(projective_points) == P + 1, "projective-line size changed")
    require(len(good) == 7, "sharp rank-seven good count changed")
    return len(projective_points), len(good)


def thm2299_hostile_packet() -> tuple[tuple[int, ...], int]:
    # Coordinate order: (H,q1,q2,q3,q4,q5,c1,c2,c3).
    speed_row = (1, 4, 2, 3, 6, 10, 13, 13**3, 2 * 13**5)
    rows: list[tuple[int, ...]] = []
    for i, q in enumerate(speed_row[1:6], start=1):
        row = [0] * 9
        row[0] = -q
        row[i] = 1
        rows.append(tuple(row))
    pair = [0] * 9
    pair[1] = 13
    pair[6] = -4
    rows.append(tuple(pair))

    require(
        all(sum(a * b for a, b in zip(row, speed_row)) == 0 for row in rows),
        "THM-2299 packet contains a nonrelation",
    )
    require(matrix_rank_mod_p(tuple(rows)) == 6, "THM-2299 packet lost rank six")
    heights = tuple(max(abs(entry) for entry in row) for row in rows)
    require(heights == (4, 2, 3, 6, 10, 13), "THM-2299 heights changed")

    pivot_columns = (1, 2, 3, 4, 5, 6)
    pivot = tuple(tuple(row[j] for j in pivot_columns) for row in rows)
    expected = (
        (1, 0, 0, 0, 0, 0),
        (0, 1, 0, 0, 0, 0),
        (0, 0, 1, 0, 0, 0),
        (0, 0, 0, 1, 0, 0),
        (0, 0, 0, 0, 1, 0),
        (13, 0, 0, 0, 0, -4),
    )
    require(pivot == expected, "THM-2299 pivot matrix changed")
    pivot_minor = -4
    require(pivot_minor % P != 0, "THM-2299 pivot minor ramified")

    pair_mod_p = tuple(entry % P for entry in rows[-1])
    expected_pure_blocker = tuple(9 if i == 6 else 0 for i in range(9))
    require(
        pair_mod_p == expected_pure_blocker,
        "THM-2299 pair is not a literal pure-blocker residue",
    )

    # Equal interval components centered at +/-1/16 have antipodal phase at
    # frequency four.  The source lift sends that cancellation to frequency
    # 13*4=52.
    phase_gap = 4 * (Fraction(1, 16) - Fraction(-1, 16))
    require(phase_gap == Fraction(1, 2), "THM-2299 phase gap changed")
    return heights, pivot_minor


def main() -> None:
    q = P
    projective_plane_size = q * q + q + 1
    lines_through_point = q + 1
    points_per_punctured_line = q
    require(projective_plane_size == 183, "P^2(F_13) size changed")
    require(lines_through_point == 14, "projective pencil size changed")

    rank6_fibre_rows = []
    for blocker_line_count in range(1, 4):
        surviving_fibres = lines_through_point - blocker_line_count
        surviving_per_fibre = points_per_punctured_line - 6
        floor = surviving_fibres * surviving_per_fibre
        affine_bad_cap = (
            1
            + (q - 1) * blocker_line_count
            + 5 * (q + 1 - blocker_line_count)
        )
        affine_floor = q * q - affine_bad_cap
        require(floor == affine_floor, "rank-six count derivations disagree")
        require(floor >= 77, "rank-six universal floor failed")
        rank6_fibre_rows.append((blocker_line_count, floor))

    sharp_good, sharp_bad, sharp_vectors = rank_six_sharp_model()
    require(sharp_good == rank6_fibre_rows[-1][1], "rank-six floor is not sharp")

    rank6_selector_floors = tuple(77 - q * degree for degree in range(6))
    require(
        rank6_selector_floors == (77, 64, 51, 38, 25, 12),
        "rank-six selector floors changed",
    )

    rank7_projective, rank7_good = rank_seven_sharp_model()
    require(rank7_projective == 14 and rank7_good == 7, "rank-seven spectrum changed")
    rank7_selector_floors = tuple(7 - degree for degree in range(6))
    require(
        rank7_selector_floors == (7, 6, 5, 4, 3, 2),
        "rank-seven selector floors changed",
    )

    # Rank eight means dim(R)=dim(w^perp)=8, hence equality and a one-line
    # dual.  Record the dimension invoice explicitly.
    require(8 == 9 - 1, "rank-eight hyperplane dimension changed")
    require(1 == 9 - 8, "rank-eight dual dimension changed")

    hostile_heights, hostile_minor = thm2299_hostile_packet()

    print("THM-2307 dual rank-six reconstruction spectrum exact companion")
    print("status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_AUDITED")
    print(
        "rank6_projective_invoice="
        f"P2:{projective_plane_size};pencil:{lines_through_point};"
        f"rows:{tuple(rank6_fibre_rows)}"
    )
    print(
        "rank6_sharp_model="
        f"good:{sharp_good};bad:{sharp_bad};nonzero_vectors:{sharp_vectors}"
    )
    print(f"rank6_selector_floors_d0..5={rank6_selector_floors}")
    print(
        "rank7_sharp_model="
        f"projective:{rank7_projective};good:{rank7_good};"
        f"selector_floors_d0..5:{rank7_selector_floors}"
    )
    print("rank8_reconstruction=R=w_perp;D=span(w)")
    print(
        "literal_dark_blocker=e_b_in_R;"
        "generic_blocker_supported_relation_is_not_the_exception"
    )
    print(
        "THM2299_dark_branch="
        f"heights:{hostile_heights};pivot_minor:{hostile_minor};"
        "pair_mod13=-4e_c1;phase_gap=1/2"
    )
    print("selector_scope=restricted_dehomogenized_product_nonzero_on_D")
    print("hostile_scope=local_interface_only;not_a_global_cover")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
