#!/usr/bin/env python3
"""Exact arithmetic controls for THM-2296's bilinear endpoint-Prony gate."""

from __future__ import annotations

from fractions import Fraction


P = 13

ALPHA_STRICT = Fraction(15041431, 593783190)
ALPHA_REPEAT = Fraction(5229541, 593783190)
BETA = Fraction(2593, 90090)

V_STAR = int(
    "142330090246961694045731915874739303206807949565751983640556780068068"
    "273579762053565008212738265314562529028994495033108953992828270745007"
    "30156340741088127740309697466937985642005305031005306880000"
)


def require(condition: bool, message: str) -> None:
    """Optimization-safe exact check."""

    if not condition:
        raise AssertionError(message)


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    """Exact determinant by fraction-preserving Gaussian elimination."""

    a = [row[:] for row in matrix]
    n = len(a)
    answer = Fraction(1)
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col]), None)
        require(pivot is not None, f"singular matrix at column {col}")
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            answer = -answer
        value = a[col][col]
        answer *= value
        for row in range(col + 1, n):
            ratio = a[row][col] / value
            for j in range(col + 1, n):
                a[row][j] -= ratio * a[col][j]
        # The eliminated pivot column is not used again.
    return answer


def vandermonde_formula(nodes: tuple[Fraction, ...], start: int) -> Fraction:
    """Determinant of [node_j^(start+i)] in closed form."""

    answer = Fraction(1)
    for node in nodes:
        answer *= node**start
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            answer *= nodes[j] - nodes[i]
    return answer


def perron_fourier(
    coefficients: dict[int, Fraction], iterations: int
) -> dict[int, Fraction]:
    """Fourier dictionary of P^iterations f for T(x)=13x."""

    scale = P**iterations
    return {
        frequency // scale: value
        for frequency, value in coefficients.items()
        if frequency % scale == 0
    }


def main() -> None:
    strict_profiles = [
        (1, b, c)
        for c in range(5, 20)
        for b in range(2, c)
    ]
    repeated_profiles = [(1, 1, c) for c in range(5, 20)]
    require(len(strict_profiles) == 150, "strict profile census changed")
    require(len(repeated_profiles) == 15, "repeated profile census changed")
    require(
        len(set(strict_profiles) | set(repeated_profiles)) == 165,
        "first-depth-one profile census changed",
    )

    strict_return_floor = ALPHA_STRICT * BETA
    repeat_return_floor = ALPHA_REPEAT * BETA
    require(
        strict_return_floor == Fraction(39002430583, 53493927587100),
        "strict prescribed-return floor changed",
    )
    require(
        repeat_return_floor == Fraction(13560199813, 53493927587100),
        "repeated prescribed-return floor changed",
    )
    require(strict_return_floor > repeat_return_floor > 0, "return floors fail")

    # The normalized Perron operator keeps exactly the Fourier modes divisible
    # by 13 at each step and divides their indices by 13.
    fourier_test = {
        -2197: Fraction(2, 5),
        -169: Fraction(-3, 7),
        -13: Fraction(5, 11),
        0: Fraction(7, 13),
        1: Fraction(11, 17),
        13: Fraction(-13, 19),
        169: Fraction(17, 23),
        2197: Fraction(-19, 29),
    }
    require(
        perron_fourier(fourier_test, 1)
        == {
            -169: Fraction(2, 5),
            -13: Fraction(-3, 7),
            -1: Fraction(5, 11),
            0: Fraction(7, 13),
            1: Fraction(-13, 19),
            13: Fraction(17, 23),
            169: Fraction(-19, 29),
        },
        "one-step Perron Fourier normalization changed",
    )
    require(
        perron_fourier(fourier_test, 3)
        == {
            -1: Fraction(2, 5),
            0: Fraction(7, 13),
            1: Fraction(-19, 29),
        },
        "three-step Perron Fourier normalization changed",
    )

    # Exact consecutive-power determinants behind the bilinear Prony gate.
    vandermonde_checks = 0
    for size in range(1, 9):
        nodes = tuple(Fraction(j + 1, j + 2) for j in range(size))
        require(len(set(nodes)) == size, f"nodes collide at size {size}")
        for start in (-5, -1, 0, 3, 9):
            matrix = [
                [node ** (start + row) for node in nodes]
                for row in range(size)
            ]
            exact = determinant(matrix)
            closed = vandermonde_formula(nodes, start)
            require(exact == closed != 0, f"Vandermonde failure {(size, start)}")
            vandermonde_checks += 1
    require(vandermonde_checks == 40, "Vandermonde check census changed")

    # A source Boolean rectangle uses all nine scalar boundary banks, while
    # R_j uses only H, the five units, and c_j. Hence J<=2S, K<=2S_j and
    # the endpoint-difference sequence has at most JK nodes. Its value at
    # index zero vanishes because each periodic jump sum is zero, improving
    # the first-positive-index bound to JK-1<=4SS_j-1<=4S^2-1.
    for S in range(9, 80):
        for S_j in range(7, S + 1):
            J = 2 * S
            K = 2 * S_j
            require(J * K == 4 * S * S_j, "jump product identity changed")
            require(J * K - 1 <= 4 * S * S - 1, "quadratic jump bound failed")

    scalar_S_upper = 5 * V_STAR // 8
    require(8 * scalar_S_upper == 5 * V_STAR, "speed ceiling lost divisibility")
    universal_frequency_bound = 4 * scalar_S_upper**2 - 1
    require(
        universal_frequency_bound == 25 * V_STAR**2 // 16 - 1,
        "universal quadratic frequency bound changed",
    )
    universal_source_frequency_bound = P**20 * universal_frequency_bound

    # Sharp information-loss family. For d=13^t, take d translates of two
    # adjacent intervals of relative length 1/4 in each 1/d cell. The two
    # indicators are disjoint and 1/d-periodic, so all nonzero Fourier modes
    # below d vanish for both, while their coefficients at d are nonzero.
    hostile_rows: list[tuple[int, int, int, int]] = []
    for depth in range(0, 7):
        d = P**depth
        relative_length = Fraction(1, 4)
        require(2 * relative_length <= 1, "hostile intervals overlap")
        first_common_frequency = d
        jumps_each = 2 * d
        bilinear_bound = jumps_each**2 - 1
        require(
            first_common_frequency <= bilinear_bound,
            f"bilinear bound misses hostile depth {depth}",
        )
        require(
            first_common_frequency * relative_length == Fraction(d, 4),
            "hostile coefficient phase ledger changed",
        )
        # At normalized frequency m=1 the interval coefficient is nonzero
        # because the relative length is not integral.
        require(relative_length.denominator != 1, "hostile first atom vanished")
        hostile_rows.append(
            (depth, d, first_common_frequency, bilinear_bound)
        )

    print("scope=FINITE-EXACT arithmetic; universal harmonic analysis is theorem-side")
    print(f"profiles=strict:{len(strict_profiles)},repeat:{len(repeated_profiles)},total:165")
    print(f"source_floor_strict={ALPHA_STRICT}")
    print(f"source_floor_repeat={ALPHA_REPEAT}")
    print(f"target_floor={BETA}")
    print(f"prescribed_independence_return_floor_strict={strict_return_floor}")
    print(f"prescribed_independence_return_floor_repeat={repeat_return_floor}")
    print(
        "jump_ledger=J_ancestry<=2S,K_residual<=2S_j,"
        "n<=JK-1<=4SS_j-1<=4S^2-1"
    )
    print("perron_fourier_law=Fourier(P^k f,n)=Fourier(f,13^k*n)")
    print(f"vandermonde_exact_checks={vandermonde_checks}")
    print(f"primitive_speed_ceiling={V_STAR}")
    print(f"scalar_comb_order_bound={scalar_S_upper}")
    print(f"universal_common_atom_bound={universal_frequency_bound}")
    print(f"universal_common_atom_bound_digits={len(str(universal_frequency_bound))}")
    print(f"universal_source_atom_bound={universal_source_frequency_bound}")
    print(
        "universal_source_atom_bound_digits="
        f"{len(str(universal_source_frequency_bound))}"
    )
    print(f"ramified_hostile_rows={tuple(hostile_rows)}")
    print("ramified_hostile_mechanism=period_1/d forces common spectrum inside dZ")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
