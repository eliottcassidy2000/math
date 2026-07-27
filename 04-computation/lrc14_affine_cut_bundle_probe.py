#!/usr/bin/env python3
"""Exact affine-cut-bundle referee for the truncated 13 x 7 Radon fold.

The program checks the full source-index equivariance, the induced DFT phase
law, the 42-cut incidence invoice, and the kappa=2 carry hostile.  It is
independent of the THM-2436 atlas engine and edits no theorem prose.
"""

from __future__ import annotations

from itertools import product


P = 13
Q = 7
N = P * Q
ZERO = (0,) * 12


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def add(x: tuple[int, ...], y: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(x, y))


def mul(x: tuple[int, ...], y: tuple[int, ...]) -> tuple[int, ...]:
    """Multiply in Z[z]/(1+z+...+z^12)."""
    raw = [0] * 23
    for i, a in enumerate(x):
        for j, b in enumerate(y):
            raw[i + j] += a * b
    for exponent in range(22, 11, -1):
        coefficient = raw[exponent]
        if coefficient:
            for j in range(12):
                raw[exponent - 12 + j] -= coefficient
            raw[exponent] = 0
    return tuple(raw[:12])


def zeta_pow(exponent: int) -> tuple[int, ...]:
    exponent %= P
    if exponent == 12:
        return (-1,) * 12
    answer = [0] * 12
    answer[exponent] = 1
    return tuple(answer)


def dft(vector: list[int], alpha: int) -> tuple[int, ...]:
    answer = ZERO
    for v, coefficient in enumerate(vector):
        if coefficient:
            answer = add(
                answer,
                tuple(coefficient * a for a in zeta_pow(-alpha * v)),
            )
    return answer


def divisors(integer: int) -> list[int]:
    return [candidate for candidate in range(1, integer + 1) if integer % candidate == 0]


def poly_div_exact(numerator: list[int], denominator: list[int]) -> list[int]:
    remainder = numerator[:]
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    while len(remainder) >= len(denominator):
        coefficient = remainder[-1]
        shift = len(remainder) - len(denominator)
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            remainder[index + shift] -= coefficient * value
        while remainder and remainder[-1] == 0:
            remainder.pop()
    require(not remainder, "cyclotomic division was not exact")
    return quotient


def cyclotomic(integer: int, cache: dict[int, list[int]]) -> list[int]:
    if integer in cache:
        return cache[integer]
    polynomial = [-1] + [0] * (integer - 1) + [1]
    for divisor in divisors(integer)[:-1]:
        polynomial = poly_div_exact(polynomial, cyclotomic(divisor, cache))
    cache[integer] = polynomial
    return polynomial


PHI91 = cyclotomic(N, {1: [-1, 1]})
DEGREE91 = len(PHI91) - 1


def clean(polynomial: dict[int, int]) -> dict[int, int]:
    return {exponent % N: value for exponent, value in polynomial.items() if value}


def cyclic_from_exponents(exponents: list[int]) -> dict[int, int]:
    answer: dict[int, int] = {}
    for exponent in exponents:
        exponent %= N
        answer[exponent] = answer.get(exponent, 0) + 1
    return clean(answer)


def cyclic_multiply(
    left: dict[int, int], right: dict[int, int]
) -> dict[int, int]:
    answer: dict[int, int] = {}
    for left_exponent, left_value in left.items():
        for right_exponent, right_value in right.items():
            exponent = (left_exponent + right_exponent) % N
            answer[exponent] = answer.get(exponent, 0) + left_value * right_value
    return clean(answer)


def cyclic_shift(polynomial: dict[int, int], exponent: int) -> dict[int, int]:
    return clean({(key + exponent) % N: value for key, value in polynomial.items()})


def cyclotomic_remainder(polynomial: dict[int, int]) -> tuple[int, ...]:
    coefficients = [0] * N
    for exponent, value in polynomial.items():
        coefficients[exponent % N] += value
    for exponent in range(N - 1, DEGREE91 - 1, -1):
        leading = coefficients[exponent]
        if leading:
            shift = exponent - DEGREE91
            for index, value in enumerate(PHI91):
                coefficients[index + shift] -= leading * value
    require(not any(coefficients[DEGREE91:]), "Phi91 reduction failed")
    return tuple(coefficients[:DEGREE91])


def kernel91(alpha_tau: int, beta: int) -> dict[int, int]:
    return cyclic_from_exponents(
        [(-7 * alpha_tau * s - 13 * beta * s) % N for s in range(Q)]
    )


def mixed91(
    defect: list[list[int]], alpha: int, beta: int
) -> dict[int, int]:
    answer: dict[int, int] = {}
    for h, r in product(range(P), range(Q)):
        coefficient = defect[h][r]
        if coefficient:
            exponent = (-7 * alpha * h - 13 * beta * r) % N
            answer[exponent] = answer.get(exponent, 0) + coefficient
    return clean(answer)


def cut_phase91(
    defect: list[list[int]], tau: int, a: int, alpha: int, beta: int
) -> dict[int, int]:
    """The c-DFT of the v-DFT, embedded through zeta13=w^7, xi7=w^13."""
    answer: dict[int, int] = {}
    for c, h, r in product(range(Q), range(P), range(Q)):
        coefficient = defect[h][r]
        if coefficient:
            s = (a * r + c) % Q
            exponent = (-7 * alpha * (h + tau * s) - 13 * beta * c) % N
            answer[exponent] = answer.get(exponent, 0) + coefficient
    return clean(answer)


def check_cut_phase_intertwiner(defect: list[list[int]]) -> tuple[int, int, int]:
    coefficient_rows = 0
    nonzero_values = 0
    zero_character_checks = 0
    for alpha, beta, tau, a in product(
        range(1, P), range(1, Q), range(1, P), range(1, Q)
    ):
        K = kernel91(alpha * tau, beta)
        lambda_exponent = (-7 * alpha * tau - 13 * beta) % N
        require(lambda_exponent != 0, "geometric ratio became one")
        require((7 * lambda_exponent) % N != 0, "geometric ratio has seventh power one")
        require(any(cyclotomic_remainder(K)), "geometric kernel vanished")

        for h, r in product(range(P), range(Q)):
            direct = sorted(
                (
                    -7 * alpha * (h + tau * ((a * r + c) % Q))
                    - 13 * beta * c
                )
                % N
                for c in range(Q)
            )
            expected = sorted(
                (
                    -7 * alpha * tau * s
                    - 13 * beta * s
                    - 7 * alpha * h
                    + 13 * beta * a * r
                )
                % N
                for s in range(Q)
            )
            require(direct == expected, "cut-phase coefficient sign failed")
            coefficient_rows += 1

        direct_value = cut_phase91(defect, tau, a, alpha, beta)
        expected_value = cyclic_multiply(K, mixed91(defect, alpha, -beta * a))
        require(direct_value == expected_value, "cut-phase factorization failed")
        require(any(cyclotomic_remainder(direct_value)), "cut-phase value vanished")
        nonzero_values += 1

    for alpha, tau, a in product(range(1, P), range(1, P), range(1, Q)):
        require(not cut_phase91(defect, tau, a, alpha, 0), "zero cut character survived")
        zero_character_checks += 1
    return coefficient_rows, nonzero_values, zero_character_checks


def affine_cut(
    defect: list[list[int]], tau: int, a: int, c: int
) -> list[int]:
    require(1 <= a < Q, "a must be a unit modulo seven")
    return [
        sum(
            defect[(v - tau * ((a * r + c) % Q)) % P][r]
            for r in range(Q)
        )
        for v in range(P)
    ]


def pullback(
    defect: list[list[int]], A: int, H: int, B: int, C: int
) -> list[list[int]]:
    Ainv = pow(A, -1, P)
    Binv = pow(B, -1, Q)
    return [
        [
            defect[(Ainv * (h - H)) % P][(Binv * (r - C)) % Q]
            for r in range(Q)
        ]
        for h in range(P)
    ]


def two_row_defect() -> list[list[int]]:
    defect = [[0] * Q for _ in range(P)]
    defect[0][5] = 1
    defect[0][3] = -1
    defect[1][5] = 1
    defect[1][4] = -1
    require(all(sum(row) == 0 for row in defect), "control is not row-zero")
    return defect


def check_source_index_equivariance() -> int:
    """Check the identity after r=B*s+C, with u=v-H.

    This exhausts every A,B,C,tau,a,c,u,s.  Since H occurs only through u,
    it covers all H and v without repeating each common translation 13 times.
    """
    checks = 0
    for A in range(1, P):
        Ainv = pow(A, -1, P)
        for B in range(1, Q):
            Binv = pow(B, -1, Q)
            for C in range(Q):
                mapped_indices = {
                    (
                        (Ainv * tau) % P,
                        (a * B) % Q,
                        (a * C + c) % Q,
                    )
                    for tau in range(P)
                    for a in range(1, Q)
                    for c in range(Q)
                }
                require(len(mapped_indices) == P * (Q - 1) * Q, "index map folded")
                require(
                    sum(1 for tau, _, _ in mapped_indices if tau == 0)
                    == (Q - 1) * Q,
                    "tau=0 bank was not preserved",
                )
                for tau in range(P):
                    mapped_tau = (Ainv * tau) % P
                    require((tau == 0) == (mapped_tau == 0), "slope type changed")
                    for a in range(1, Q):
                        mapped_a = (a * B) % Q
                        for c in range(Q):
                            mapped_c = (a * C + c) % Q
                            for u in range(P):
                                mapped_u = (Ainv * u) % P
                                for s in range(Q):
                                    r = (B * s + C) % Q
                                    require(
                                        (Binv * (r - C)) % Q == s,
                                        "row substitution has inverse error",
                                    )
                                    left_h = (
                                        Ainv * (u - tau * ((a * r + c) % Q))
                                    ) % P
                                    right_h = (
                                        mapped_u
                                        - mapped_tau
                                        * ((mapped_a * s + mapped_c) % Q)
                                    ) % P
                                    require(left_h == right_h, "source index mismatch")
                                    checks += 1
    return checks


def check_dft_phase_law() -> int:
    """Check the output-affine Fourier law on every delta basis vector."""
    checks = 0
    for A in range(1, P):
        Ainv = pow(A, -1, P)
        for H in range(P):
            for alpha in range(P):
                phase = zeta_pow(-alpha * H)
                for source in range(P):
                    vector = [0] * P
                    vector[source] = 1
                    moved = [vector[(Ainv * (v - H)) % P] for v in range(P)]
                    left = dft(moved, alpha)
                    right = mul(phase, dft(vector, (alpha * A) % P))
                    require(left == right, "DFT phase law failed")
                    checks += 1
    return checks


def main() -> None:
    cuts = {
        tuple((a * r + c) % Q for r in range(Q))
        for a in range(1, Q)
        for c in range(Q)
    }
    require(len(cuts) == 42, "affine cut bank is not free")
    require(all(sorted(cut) == list(range(Q)) for cut in cuts), "cut not bijective")
    good_invoice = len(cuts) * (P - Q + 1)
    nonzero_components = len(cuts) * (P - 1)
    require((good_invoice, nonzero_components) == (294, 504), "invoice drifted")
    print("affine cuts / nonzero-slope components:", len(cuts), nonzero_components)
    print("universal good-component invoice:", good_invoice, "/", nonzero_components)

    orbit = {
        (
            pow(A, -1, P),
            B,
            C,
            (-pow(A, -1, P) * H) % P,
        )
        for A in range(1, P)
        for H in range(P)
        for B in range(1, Q)
        for C in range(Q)
    }
    # This is the orbit of the seed (tau,a,c,v)=(1,1,0,0).
    require(len(orbit) == 6_552, "affine gauge action is not transitive")
    require(len(orbit) == (P - 1) * (Q - 1) * Q * P, "index size drifted")
    print("affine gauge orbit / nonzero evaluation-index set:", len(orbit), 6_552)
    print("stabilizer size: 1 (the action is regular)")

    multiplicities = []
    for h in range(P):
        for r in range(Q):
            multiplicity = 0
            for tau in range(1, P):
                for a in range(1, Q):
                    for c in range(Q):
                        multiplicity += sum(
                            h == (v - tau * ((a * r + c) % Q)) % P
                            for v in range(P)
                        )
            multiplicities.append(multiplicity)
    require(set(multiplicities) == {504}, "uniform collapse is not uniform")
    print("uniform invariant linear-collapse coefficient per source:", 504)

    source_checks = check_source_index_equivariance()
    require(source_checks == 25_041_744, "source-check universe drifted")
    print("exact source-index identities:", source_checks)
    print("all 504 nonzero components permute; all 42 tau=0 components stay zero-type")

    phase_checks = check_dft_phase_law()
    require(phase_checks == 26_364, "phase-check universe drifted")
    print("exact DFT basis phase identities:", phase_checks)

    defect = two_row_defect()
    factor_rows, cut_values, zero_cut_values = check_cut_phase_intertwiner(defect)
    require(factor_rows == 471_744, "factorization universe drifted")
    require(cut_values == 5_184, "cut-phase value universe drifted")
    require(zero_cut_values == 864, "zero-character universe drifted")
    print("exact cut-phase factorization coefficient rows:", factor_rows)
    print(
        "nonzero mixed cut-phase values / zero-character checks:",
        cut_values,
        zero_cut_values,
    )
    component_counts = []
    for a in range(1, Q):
        for c in range(Q):
            component_counts.append(
                sum(
                    affine_cut(defect, tau, a, c) != [0] * P
                    for tau in range(1, P)
                )
            )
    require(component_counts == [12] * 42, "two-row control lost a component")
    print("THM-2506 two-row control good components:", sum(component_counts), "/ 504")

    bank = [
        value
        for tau in range(1, P)
        for a in range(1, Q)
        for c in range(Q)
        for value in affine_cut(defect, tau, a, c)
    ]
    require(sum(bank) == 0, "uniform invariant linear collapse did not vanish")
    energy = sum(value * value for value in bank)
    require(energy > 0, "quadratic energy failed positive control")
    print("uniform invariant linear collapse on row-zero control:", sum(bank))

    moved = pullback(defect, A=1, H=2, B=1, C=2)
    left = affine_cut(moved, tau=1, a=1, c=0)
    matched = affine_cut(defect, tau=1, a=1, c=2)
    right = [matched[(v - 2) % P] for v in range(P)]
    require(left == right, "kappa=2 hostile did not become a component permutation")
    original = affine_cut(defect, tau=1, a=1, c=0)
    require(
        not any(left == original[shift:] + original[:shift] for shift in range(P)),
        "single-cut hostile disappeared",
    )
    print("kappa=2,tau=1: cut (1,0) maps to cut (1,2) plus v -> v-2: PASS")
    print("single-cut output was not any translate of the old cut (1,0): PASS")
    cut_phase_covariance = 0
    for alpha, beta, tau, a in product(
        range(1, P), range(1, Q), range(1, P), range(1, Q)
    ):
        left_phase = cut_phase91(moved, tau, a, alpha, beta)
        right_phase = cyclic_shift(
            cut_phase91(defect, tau, a, alpha, beta),
            -7 * alpha * 2 + 13 * beta * a * 2,
        )
        require(left_phase == right_phase, "mixed cut-phase covariance failed")
        cut_phase_covariance += 1
    require(cut_phase_covariance == 5_184, "mixed covariance universe drifted")
    print("kappa=2 mixed cut-phase covariance values:", cut_phase_covariance)
    moved_bank = [
        value
        for tau in range(1, P)
        for a in range(1, Q)
        for c in range(Q)
        for value in affine_cut(moved, tau, a, c)
    ]
    moved_energy = sum(value * value for value in moved_bank)
    require(moved_energy == energy, "quadratic bank energy is not gauge invariant")
    print("quadratic bank energy before / after kappa=2:", energy, moved_energy)
    print("quadratic energy is invariant but phase-forgetting: VERIFIED BOUNDARY")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
