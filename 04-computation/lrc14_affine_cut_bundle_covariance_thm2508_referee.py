#!/usr/bin/env python3
"""Independent exact referee for the THM-2508 affine-cut bundle.

The calculation is over F_547, which contains roots of unity of exact orders
13 and 7.  It verifies the operator identity coefficientwise, not only on a
sample defect.  The concrete two-row THM-2506 hostile is then used as a
positive control for all primitive colours and all affine cuts.
"""

from math import isqrt


P13 = 13
P7 = 7
MOD = 547  # 547 - 1 = 6 * 91


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in range(2, isqrt(n) + 1):
        if n % p == 0:
            return False
    return True


def primitive_root(prime: int) -> int:
    factors = (2, 3, 7, 13)
    for g in range(2, prime):
        if all(pow(g, (prime - 1) // q, prime) != 1 for q in factors):
            return g
    raise RuntimeError("no primitive root")


require(is_prime(MOD), "modulus is not prime")
GEN = primitive_root(MOD)
ZETA = pow(GEN, (MOD - 1) // P13, MOD)
XI = pow(GEN, (MOD - 1) // P7, MOD)
require(pow(ZETA, P13, MOD) == 1 and ZETA != 1, "bad 13th root")
require(pow(XI, P7, MOD) == 1 and XI != 1, "bad 7th root")


def root_power(root: int, exponent: int, order: int) -> int:
    return pow(root, exponent % order, MOD)


def mixed_dft(defect: list[list[int]], alpha: int, beta: int) -> int:
    return sum(
        defect[h][r]
        * root_power(ZETA, -alpha * h, P13)
        * root_power(XI, -beta * r, P7)
        for h in range(P13)
        for r in range(P7)
    ) % MOD


def radon(defect: list[list[int]], tau: int, a: int, c: int) -> list[int]:
    return [
        sum(defect[(v - tau * ((a * r + c) % P7)) % P13][r]
            for r in range(P7))
        for v in range(P13)
    ]


def dft13(vector: list[int], alpha: int) -> int:
    return sum(
        vector[v] * root_power(ZETA, -alpha * v, P13)
        for v in range(P13)
    ) % MOD


def cut_coefficient(
    defect: list[list[int]], tau: int, a: int, alpha: int, beta: int
) -> int:
    return sum(
        root_power(XI, -beta * c, P7)
        * dft13(radon(defect, tau, a, c), alpha)
        for c in range(P7)
    ) % MOD


def kernel(slope: int, beta: int) -> int:
    return sum(
        root_power(XI, -beta * j, P7)
        * root_power(ZETA, -slope * j, P13)
        for j in range(P7)
    ) % MOD


def translate(defect: list[list[int]], A: int, C: int) -> list[list[int]]:
    return [
        [defect[(h - A) % P13][(r - C) % P7] for r in range(P7)]
        for h in range(P13)
    ]


def two_row_hostile() -> list[list[int]]:
    defect = [[0] * P7 for _ in range(P13)]
    defect[0][5] = 1
    defect[0][3] = -1
    defect[1][5] = 1
    defect[1][4] = -1
    require(all(sum(row) == 0 for row in defect), "hostile is not row-zero")
    return defect


def main() -> None:
    # Coefficientwise proof of
    # C_(tau,a)(alpha,beta)=K(alpha*tau,beta)D~(alpha,-beta*a).
    factorization_rows = 0
    nonzero_kernels = 0
    for alpha in range(1, P13):
        for beta in range(1, P7):
            for tau in range(1, P13):
                for a in range(1, P7):
                    K = kernel(alpha * tau, beta)
                    require(K != 0, "geometric kernel vanished")
                    nonzero_kernels += 1
                    for h in range(P13):
                        for r in range(P7):
                            direct = sum(
                                root_power(XI, -beta * c, P7)
                                * root_power(
                                    ZETA,
                                    -alpha * (h + tau * ((a * r + c) % P7)),
                                    P13,
                                )
                                for c in range(P7)
                            ) % MOD
                            expected = (
                                K
                                * root_power(ZETA, -alpha * h, P13)
                                * root_power(XI, beta * a * r, P7)
                            ) % MOD
                            require(direct == expected, "operator factorization failed")
                            factorization_rows += 1

    defect = two_row_hostile()
    primitive = [
        mixed_dft(defect, alpha, beta)
        for alpha in range(1, P13)
        for beta in range(1, P7)
    ]
    require(all(value != 0 for value in primitive), "primitive hostile colour vanished")

    cut_values = 0
    for alpha in range(1, P13):
        for beta in range(1, P7):
            for tau in range(1, P13):
                for a in range(1, P7):
                    direct = cut_coefficient(defect, tau, a, alpha, beta)
                    expected = (
                        kernel(alpha * tau, beta)
                        * mixed_dft(defect, alpha, -beta * a)
                    ) % MOD
                    require(direct == expected and direct != 0,
                            "primitive cut coefficient failed")
                    cut_values += 1

    # The zero cut-character cannot descend to a pure F_13 current.
    zero_cut_values = 0
    for alpha in range(1, P13):
        for tau in range(1, P13):
            for a in range(1, P7):
                require(cut_coefficient(defect, tau, a, alpha, 0) == 0,
                        "zero cut-character survived row-sum law")
                zero_cut_values += 1

    # Raw covariance under every CRT translation.  The cut coordinate absorbs
    # the representative carry by c -> c+a*C.
    covariance_rows = 0
    for A in range(P13):
        for C in range(P7):
            moved = translate(defect, A, C)
            for tau in range(1, P13):
                for a in range(1, P7):
                    for c in range(P7):
                        left = radon(moved, tau, a, c)
                        right0 = radon(defect, tau, a, (c + a * C) % P7)
                        right = [right0[(v - A) % P13] for v in range(P13)]
                        require(left == right, "raw cut-bundle covariance failed")
                        covariance_rows += P13

    print("THM-2508 affine-cut bundle independent exact referee: PASS")
    print(f"field=F_{MOD}; generator={GEN}; ord(zeta)={P13}; ord(xi)={P7}")
    print(f"operator_factorization_coefficients={factorization_rows}")
    print(f"nonzero_geometric_kernels={nonzero_kernels}")
    print(f"hostile_primitive_colours={len(primitive)}; nonzero_cut_coefficients={cut_values}")
    print(f"zero_cut_character_checks={zero_cut_values}; all_zero=PASS")
    print(f"raw_translation_covariance_components={covariance_rows}")
    print("mixed cut character is necessary; pure F_13 descent remains zero")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
