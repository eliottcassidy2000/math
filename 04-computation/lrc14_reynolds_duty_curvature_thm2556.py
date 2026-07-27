#!/usr/bin/env python3
"""Exact companion for THM-2556: Reynolds duty and mixed curvature."""

from __future__ import annotations

from fractions import Fraction


P = 13
VARIABLES = 6  # a0,...,a5
CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def affine_solution_count(equations: list[tuple[list[int], int]]) -> int:
    """Count solutions of coeffs.variable=rhs over F_13."""
    if not equations:
        return P**VARIABLES
    matrix = [[c % P for c in coeffs] + [rhs % P]
              for coeffs, rhs in equations]
    row = 0
    for col in range(VARIABLES):
        pivot = next((i for i in range(row, len(matrix))
                      if matrix[i][col] % P), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        inv = pow(matrix[row][col], -1, P)
        matrix[row] = [(v * inv) % P for v in matrix[row]]
        for i in range(len(matrix)):
            if i != row and matrix[i][col]:
                scale = matrix[i][col]
                matrix[i] = [(matrix[i][j] - scale * matrix[row][j]) % P
                             for j in range(VARIABLES + 1)]
        row += 1
        if row == len(matrix):
            break
    for values in matrix:
        if all(values[j] == 0 for j in range(VARIABLES)) and values[-1] != 0:
            return 0
    rank = sum(any(values[j] != 0 for j in range(VARIABLES))
               for values in matrix)
    return P ** (VARIABLES - rank)


def hyperplanes(x: int, y: int) -> list[tuple[list[int], int]]:
    """The nine zero-coordinate equations from packet chart (13)."""
    rows: list[tuple[list[int], int]] = []
    # -a1,...,-a5=0.
    for variable in range(1, 6):
        coeffs = [0] * VARIABLES
        coeffs[variable] = -1
        rows.append((coeffs, 0))
    # a1+...+a5=0.
    coeffs = [0] * VARIABLES
    for variable in range(1, 6):
        coeffs[variable] = 1
    rows.append((coeffs, 0))
    # -a0=0.
    coeffs = [0] * VARIABLES
    coeffs[0] = -1
    rows.append((coeffs, 0))
    # x-a1=0 and y-a2=0.
    coeffs = [0] * VARIABLES
    coeffs[1] = -1
    rows.append((coeffs, -x))
    coeffs = [0] * VARIABLES
    coeffs[2] = -1
    rows.append((coeffs, -y))
    require(len(rows) == 9, "packet must have nine coordinate walls")
    return rows


def inclusion_exclusion_count(x: int, y: int) -> int:
    walls = hyperplanes(x, y)
    total = 0
    for mask in range(1 << len(walls)):
        equations = [walls[i] for i in range(len(walls)) if mask >> i & 1]
        term = affine_solution_count(equations)
        total += -term if mask.bit_count() & 1 else term
    return total


def closed_count(x: int, y: int) -> int:
    return (2_316_060
            + 210_552 * (int(x == 0) + int(y == 0))
            + 12 * int((x + y) % P == 0)
            + 19_128 * int(x == 0 and y == 0))


def fourier_bins(values: dict[tuple[int, int], int], u: int, v: int) -> list[int]:
    bins = [0] * P
    for (x, y), value in values.items():
        bins[(-u * x - v * y) % P] += value
    return bins


def cyclotomic_equals_integer(bins: list[int], integer: int) -> bool:
    difference = bins.copy()
    difference[0] -= integer
    # Phi_13 is the only rational relation among 1,zeta,...,zeta^12.
    return len(set(difference)) == 1


def cyclotomic_zero(bins: list[int]) -> bool:
    return len(set(bins)) == 1


def main() -> None:
    values: dict[tuple[int, int], int] = {}
    ie_checks = 0
    for x in range(P):
        for y in range(P):
            closed = closed_count(x, y)
            direct = inclusion_exclusion_count(x, y)
            require(direct == closed, "closed count disagrees with inclusion-exclusion")
            values[(x, y)] = closed
            ie_checks += 1

    histogram: dict[int, int] = {}
    for value in values.values():
        histogram[value] = histogram.get(value, 0) + 1
    require(histogram == {2_756_304: 1, 2_526_612: 24,
                          2_316_072: 12, 2_316_060: 132},
            "canonical duty histogram failed")

    n7 = 6 * (6**8 + 6) // 7
    fibre_size = 7**8 * 13**6
    require(n7 == 1_439_676, "septimal unit count failed")
    require(fibre_size == 27_825_593_350_009, "target fibre size failed")

    # Exact l2 operator norm of the within-fibre covariance residual.
    unit_counts = {q: n7 * count for q, count in values.items()}
    require(all(0 <= count <= fibre_size for count in unit_counts.values()),
            "unit count left its target fibre")
    variance_numerators = {
        q: count * (fibre_size - count)
        for q, count in unit_counts.items()
    }
    maximum_variance = max(variance_numerators.values())
    maximizing_fibres = [q for q, value in variance_numerators.items()
                         if value == maximum_variance]
    require(maximizing_fibres == [(0, 0)],
            "canonical covariance norm has the wrong maximizing fibre")
    residual_norm_sq = Fraction(maximum_variance, fibre_size)
    require(residual_norm_sq
            == Fraction(1_932_053_149_688_864_170_956_480,
                        567_869_252_041),
            "exact covariance operator norm failed")

    fourier_checks = 0
    all_nontrivial = 0
    for u in range(P):
        for v in range(P):
            bins = fourier_bins(values, u, v)
            if (u, v) == (0, 0):
                predicted = 396_907_776
            else:
                predicted = (19_128
                             + 2_737_176 * (int(u == 0) + int(v == 0))
                             + 156 * int(u == v))
                require(predicted > 0, "nontrivial duty Fourier coefficient vanished")
                all_nontrivial += 1
            require(cyclotomic_equals_integer(bins, predicted),
                    "duty Fourier formula failed")
            fourier_checks += 1
    require(all_nontrivial == 168, "duty spectrum is not full")

    stabilizer_checks = 0
    curvature_character_checks = 0
    curvature_nonzero = 0
    energy_classes: dict[str, set[int]] = {"axis": set(), "anti": set(), "generic": set()}
    for gx in range(P):
        for gy in range(P):
            if (gx, gy) == (0, 0):
                continue
            shifted_equal = all(values[(x, y)]
                                == values[((x - gx) % P, (y - gy) % P)]
                                for x in range(P) for y in range(P))
            require(not shifted_equal, "duty profile has a nontrivial stabilizer")
            stabilizer_checks += 1

            defect: dict[tuple[int, int], int] = {}
            for x in range(P):
                for y in range(P):
                    defect[(x, y)] = (7 * values[(x, y)]
                                      - sum(values[((x - j * gx) % P,
                                                    (y - j * gy) % P)]
                                            for j in range(7)))
            energy = sum(value * value for value in defect.values())
            if gx == 0 or gy == 0:
                energy_classes["axis"].add(energy)
            elif (gx + gy) % P == 0:
                energy_classes["anti"].add(energy)
            else:
                energy_classes["generic"].add(energy)

            live_for_gain = 0
            for u in range(P):
                for v in range(P):
                    bins = fourier_bins(defect, u, v)
                    is_zero = cyclotomic_zero(bins)
                    orthogonal = (u * gx + v * gy) % P == 0
                    require(is_zero == orthogonal,
                            "curvature zero set is not the orthogonal line")
                    if not is_zero:
                        live_for_gain += 1
                        curvature_nonzero += 1
                    curvature_character_checks += 1
            require(live_for_gain == 156, "gain does not have 156 live characters")

    require(energy_classes == {
        "axis": {24_559_042_191_264},
        "anti": {49_102_678_687_104},
        "generic": {49_102_698_046_752},
    }, "curvature energy classes failed")
    require(stabilizer_checks == 168, "wrong stabilizer scan size")
    require(curvature_character_checks == 28_392,
            "wrong gain/character check count")
    require(curvature_nonzero == 168 * 156,
            "wrong live curvature character count")

    # Toy carrier: Q(x,y)=x and S=(F13*)^2.
    raw_neutral_rank = sum(1 for x in range(P) for y in range(P)
                           if int(x != 0 and y != 0)
                           != int(x != 0 and (y - 1) % P != 0))
    require(raw_neutral_rank == 24, "raw neutral unit-wall rank failed")
    nu = [0] + [12] * (P - 1)  # common denominator 13
    reynolds_target_rank = sum(1 for x in range(P)
                               if nu[x] != nu[(x - 1) % P])
    require(reynolds_target_rank == 2, "Reynolds target commutator rank failed")
    toy_fourier_colours = 0
    for u in range(1, P):
        bins = [0] * P
        for x in range(P):
            bins[(-u * x) % P] += nu[x]
        require(not cyclotomic_zero(bins), "toy Reynolds colour vanished")
        toy_fourier_colours += 1
    require(toy_fourier_colours == 12, "wrong toy target-colour count")

    # c_host=1_(y=0): U is one, J is zero, and target x-shifts preserve it.
    hostile_U = [sum(1 for y in range(P) if y == 0) for _x in range(P)]
    hostile_J = [sum(1 for y in range(P) if y == 0 and _x != 0 and y != 0)
                 for _x in range(P)]
    require(hostile_U == [1] * P and hostile_J == [0] * P,
            "toy nonunit-section hostile failed")

    print("THM-2556 exact Reynolds-duty curvature referee")
    print(f"target_fibres={ie_checks} inclusion_exclusion_subsets_each=512")
    print("duty_histogram=2756304:1,2526612:24,2316072:12,2316060:132")
    print(f"N7={n7} fibre_size={fibre_size}")
    print(f"residual_norm_sq={residual_norm_sq} residual_max_fibre=0,0")
    print(f"duty_fourier_checks={fourier_checks} nontrivial_support={all_nontrivial}")
    print(f"stabilizer_checks={stabilizer_checks}")
    print(f"curvature_character_checks={curvature_character_checks} live={curvature_nonzero} live_per_gain=156")
    print("energy_axis=24559042191264 energy_antidiagonal=49102678687104 energy_generic=49102698046752")
    print(f"toy_raw_neutral_rank={raw_neutral_rank} toy_reynolds_target_rank={reynolds_target_rank} toy_target_colours={toy_fourier_colours}")
    print(f"explicit_require_checks={CHECKS}")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
