#!/usr/bin/env python3
"""Independent exact referee for the reserved THM-2512 candidate.

This companion checks only the algebraic live-table transplant and a finite
Boolean-before-DFT model of the forward temporal composition.  It does not
edit or promote the reserved theorem.
"""

from fractions import Fraction
from itertools import product
from math import isqrt


P = 13
Q = 7
MOD = 547  # MOD - 1 = 6 * 7 * 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_prime(integer: int) -> bool:
    if integer < 2:
        return False
    return all(integer % divisor for divisor in range(2, isqrt(integer) + 1))


def primitive_root(prime: int) -> int:
    factors = (2, 3, 7, 13)
    for candidate in range(2, prime):
        if all(pow(candidate, (prime - 1) // factor, prime) != 1 for factor in factors):
            return candidate
    raise RuntimeError("primitive root not found")


require(is_prime(MOD), "finite-field modulus is not prime")
GENERATOR = primitive_root(MOD)
ZETA = pow(GENERATOR, (MOD - 1) // P, MOD)
XI = pow(GENERATOR, (MOD - 1) // Q, MOD)
require(pow(ZETA, P, MOD) == 1 and ZETA != 1, "bad thirteenth root")
require(pow(XI, Q, MOD) == 1 and XI != 1, "bad seventh root")


def root_power(root: int, exponent: int, order: int) -> int:
    return pow(root, exponent % order, MOD)


def fraction_mod(value: Fraction) -> int:
    return value.numerator * pow(value.denominator, -1, MOD) % MOD


def anova(table: list[list[Fraction]]) -> list[list[Fraction]]:
    """Return I_(ell,s), with both marginal means zero."""

    row_means = [sum(row, Fraction(0)) / P for row in table]
    column_means = [
        sum((table[ell][s] for ell in range(Q)), Fraction(0)) / Q
        for s in range(P)
    ]
    grand = sum((sum(row, Fraction(0)) for row in table), Fraction(0)) / (P * Q)
    answer = [
        [table[ell][s] - row_means[ell] - column_means[s] + grand for s in range(P)]
        for ell in range(Q)
    ]
    require(
        all(sum(answer[ell], Fraction(0)) == 0 for ell in range(Q)),
        "ANOVA row mean survived",
    )
    require(
        all(sum((answer[ell][s] for ell in range(Q)), Fraction(0)) == 0 for s in range(P)),
        "ANOVA column mean survived",
    )
    return answer


def transpose_interaction(interaction: list[list[Fraction]]) -> list[list[Fraction]]:
    return [[interaction[ell][s] for ell in range(Q)] for s in range(P)]


def replica_table(add_defect: bool) -> list[list[Fraction]]:
    owner_mass = Fraction(5)
    word = [Fraction(s * s) for s in range(P)]
    table = [
        [word[s] + (owner_mass if ell == 0 else 0) for s in range(P)]
        for ell in range(Q)
    ]
    if add_defect:
        table[1][1] += 1
    require(
        all(table[ell][0] == owner_mass * (ell == 0) for ell in range(Q)),
        "owner anchor drifted",
    )
    require(all(value >= 0 for row in table for value in row), "control is not nonnegative")
    return table


def rectangles(table: list[list[Fraction]]) -> tuple[Fraction, ...]:
    return tuple(
        table[0][s] - table[ell][s] - table[0][0] + table[ell][0]
        for ell in range(1, Q)
        for s in range(1, P)
    )


def mixed_dft(defect: list[list[Fraction]], alpha: int, gamma: int) -> int:
    return sum(
        fraction_mod(defect[s][ell])
        * root_power(ZETA, -alpha * s, P)
        * root_power(XI, -gamma * ell, Q)
        for s in range(P)
        for ell in range(Q)
    ) % MOD


def owner_dft(table: list[list[Fraction]], kappa: int, b: int) -> int:
    unnormalized = sum(
        fraction_mod(table[ell][s])
        * root_power(XI, kappa * ell, Q)
        * root_power(ZETA, b * s, P)
        for ell in range(Q)
        for s in range(P)
    ) % MOD
    return unnormalized * pow(P * Q, -1, MOD) % MOD


def radon(
    defect: list[list[Fraction]], tau: int, a: int, c: int
) -> list[Fraction]:
    return [
        sum(
            (defect[(v - tau * ((a * ell + c) % Q)) % P][ell] for ell in range(Q)),
            Fraction(0),
        )
        for v in range(P)
    ]


def dft13(vector: list[Fraction], alpha: int) -> int:
    return sum(
        fraction_mod(vector[v]) * root_power(ZETA, -alpha * v, P)
        for v in range(P)
    ) % MOD


def cut_phase(
    defect: list[list[Fraction]], tau: int, a: int, alpha: int, beta: int
) -> int:
    return sum(
        root_power(XI, -beta * c, Q) * dft13(radon(defect, tau, a, c), alpha)
        for c in range(Q)
    ) % MOD


def geometric_kernel(alpha_tau: int, beta: int) -> int:
    return sum(
        root_power(ZETA, -alpha_tau * j, P) * root_power(XI, -beta * j, Q)
        for j in range(Q)
    ) % MOD


def anova_mod(table: list[list[int]]) -> list[list[int]]:
    inv_p = pow(P, -1, MOD)
    inv_q = pow(Q, -1, MOD)
    inv_pq = pow(P * Q, -1, MOD)
    row_means = [sum(row) % MOD * inv_p % MOD for row in table]
    column_means = [sum(table[ell][s] for ell in range(Q)) % MOD * inv_q % MOD for s in range(P)]
    grand = sum(map(sum, table)) % MOD * inv_pq % MOD
    return [
        [(table[ell][s] - row_means[ell] - column_means[s] + grand) % MOD for s in range(P)]
        for ell in range(Q)
    ]


def cut_phase_mod(
    transposed: list[list[int]], tau: int, a: int, alpha: int, beta: int
) -> int:
    answer = 0
    for c, s, ell in product(range(Q), range(P), range(Q)):
        v = (s + tau * ((a * ell + c) % Q)) % P
        answer += (
            transposed[s][ell]
            * root_power(XI, -beta * c, Q)
            * root_power(ZETA, -alpha * v, P)
        )
    return answer % MOD


def selected_cut_functional_mod(table: list[list[int]]) -> int:
    interaction = anova_mod(table)
    transposed = [[interaction[ell][s] for ell in range(Q)] for s in range(P)]
    return cut_phase_mod(transposed, tau=1, a=1, alpha=1, beta=1)


def main() -> None:
    replica = replica_table(add_defect=False)
    defect_table = replica_table(add_defect=True)
    replica_interaction = anova(replica)
    defect_interaction = anova(defect_table)
    replica_d = transpose_interaction(replica_interaction)
    defect_d = transpose_interaction(defect_interaction)

    require(not any(any(value for value in row) for row in replica_interaction), "replica interaction survived")
    require(any(any(value for value in row) for row in defect_interaction), "defect interaction vanished")
    require(not any(rectangles(replica)), "replica rectangle survived")
    require(any(rectangles(defect_table)), "defect rectangle vanished")

    # Rectangle(A)=rectangle(I) is checked coefficientwise on the 91 basis tables.
    rectangle_basis_checks = 0
    for owner_ell, target_s in product(range(Q), range(P)):
        basis = [[Fraction(0) for _ in range(P)] for _ in range(Q)]
        basis[owner_ell][target_s] = 1
        interaction = anova(basis)
        require(rectangles(basis) == rectangles(interaction), "rectangle/ANOVA operator mismatch")
        rectangle_basis_checks += len(rectangles(basis))

    # Canon signs: d~(alpha,gamma)=91 Ahat(-gamma,-alpha).
    dft_sign_checks = 0
    for alpha in range(1, P):
        for gamma in range(1, Q):
            left = mixed_dft(defect_d, alpha, gamma)
            right = (P * Q) * owner_dft(defect_table, -gamma, -alpha) % MOD
            require(left == right, "mixed DFT transpose sign failed")
            require(left != 0, "one-entry rectangle lost a mixed colour")
            dft_sign_checks += 1

    # Direct cut DFT and the THM-2508 factorization on every primitive index.
    factorization_checks = 0
    replica_zero_checks = 0
    for alpha, beta, tau, a in product(range(1, P), range(1, Q), range(1, P), range(1, Q)):
        kernel = geometric_kernel(alpha * tau, beta)
        require(kernel != 0, "geometric cut kernel vanished")
        direct = cut_phase(defect_d, tau, a, alpha, beta)
        via_mixed = kernel * mixed_dft(defect_d, alpha, -beta * a) % MOD
        via_owner = kernel * (P * Q) * owner_dft(defect_table, beta * a, -alpha) % MOD
        require(direct == via_mixed == via_owner, "cut/owner factorization failed")
        require(direct != 0, "defect primitive cut coefficient vanished")
        require(cut_phase(replica_d, tau, a, alpha, beta) == 0, "replica cut coefficient survived")
        factorization_checks += 1
        replica_zero_checks += 1

    # Boolean-before-DFT finite temporal toy.  Raw entries stay nonnegative;
    # ANOVA and the selected cut character are applied only after integration.
    # The physical grid is (old root branch w, slow cell u, future cell v).
    slow_size = 3
    future_size = 5
    physical_integrals = [[Fraction(0) for _ in range(P)] for _ in range(Q)]
    disintegrated_integrals = [
        [[Fraction(0) for _ in range(P)] for _ in range(Q)] for _ in range(P)
    ]
    for ell, s in product(range(Q), range(P)):
        raw_value = defect_table[ell][s]
        physical_total = Fraction(0)
        for w, u, v in product(range(P), range(slow_size), range(future_size)):
            carrier = Fraction((w + 1) * (u + 1), 14)  # positive, mean one in (w,u)
            physical_total += raw_value * carrier
        physical_integrals[ell][s] = physical_total / (P * slow_size * future_size)
        for w in range(P):
            base_total = Fraction(0)
            for u, v in product(range(slow_size), range(future_size)):
                carrier = Fraction((w + 1) * (u + 1), 14)
                base_total += raw_value * carrier / P
            disintegrated_integrals[w][ell][s] = base_total / (slow_size * future_size)
        require(
            sum((disintegrated_integrals[w][ell][s] for w in range(P)), Fraction(0))
            == physical_integrals[ell][s]
            == raw_value,
            "old-root disintegration did not recover the raw integral",
        )

    raw_table_mod = [[fraction_mod(value) for value in row] for row in physical_integrals]
    selected = selected_cut_functional_mod(raw_table_mod)
    require(selected != 0, "selected live mixed cut functional vanished")

    # One future Boolean root-shift cell (shift=1) occupies two of five fast
    # cells.  DFT in the shift is taken only after multiplying raw components.
    rho = Fraction(2, future_size)
    rho_mod = fraction_mod(rho)
    raw_future_product_checks = 0
    for w, ell, s in product(range(P), range(Q), range(P)):
        raw_value = defect_table[ell][s]
        product_total = Fraction(0)
        for u, v in product(range(slow_size), range(future_size)):
            carrier = Fraction((w + 1) * (u + 1), 14)
            collision_cell = int(v < 2)
            product_total += raw_value * carrier * collision_cell / P
        product_mean = product_total / (slow_size * future_size)
        require(
            product_mean == rho * disintegrated_integrals[w][ell][s],
            "raw Boolean future cell did not mix componentwise",
        )
        raw_future_product_checks += 1
    future_colour_checks = 0
    for colour in range(1, P):
        collision_mean = rho_mod * root_power(ZETA, -colour, P) % MOD
        scaled_raw = [
            [(collision_mean * value) % MOD for value in row]
            for row in raw_table_mod
        ]
        observed = selected_cut_functional_mod(scaled_raw)
        require(observed == collision_mean * selected % MOD, "componentwise future multiplication lost functional")
        require(observed != 0, "future colour killed selected functional")
        future_colour_checks += 1

    # Rebase hostile from THM-2478: same future base, two ancestry sheets,
    # old deep phases 0 and 1/13.  At radius 1/14 they have opposite truth.
    depth = 1
    delay = 2
    deep_speed = P**depth
    future_base = Fraction(0)
    sheet_zero = 0
    sheet_one = 1
    phase_zero = Fraction(deep_speed * (future_base + sheet_zero), P**delay)
    phase_one = Fraction(deep_speed * (future_base + sheet_one), P**delay)

    def dangerous(phase: Fraction) -> bool:
        reduced = phase - (phase.numerator // phase.denominator)
        distance = min(reduced, 1 - reduced)
        return distance < Fraction(1, 14)

    require(phase_zero == 0 and phase_one == Fraction(1, P), "rebase phases drifted")
    require(dangerous(phase_zero) and not dangerous(phase_one), "rebase hostile lost opposite truth")

    print("THM-2512 lawful interaction/cut transplant independent exact referee: PASS")
    print(f"field=F_{MOD}; generator={GENERATOR}; ord(zeta)={P}; ord(xi)={Q}")
    print(f"anova_rectangle_basis_checks={rectangle_basis_checks}")
    print(f"mixed_dft_sign_checks={dft_sign_checks}")
    print(f"primitive_cut_factorizations={factorization_checks}; replica_zero_checks={replica_zero_checks}")
    print("replica_interaction=ZERO; one_rectangle_defect=ALL_5184_PRIMITIVE_CUTS_NONZERO")
    print(f"old_root_disintegration_components={P * Q * P}; recovery_checks={Q * P}")
    print(
        f"raw_future_product_checks={raw_future_product_checks}; "
        f"boolean_before_dft_future_colours={future_colour_checks}; "
        "selected_functional_retained=PASS"
    )
    print("rebase_hostile_phases=0,1/13; danger_truth=1,0")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
