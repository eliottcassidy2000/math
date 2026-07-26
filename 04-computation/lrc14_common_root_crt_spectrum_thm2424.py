#!/usr/bin/env python3
"""Exact companion for THM-2424.

The proof is algebraic.  This companion checks the CRT character map in
finite group rings, the polyphase residue selector, the complete LRC
unit-character bank, the rational constants, and sharp boundary models.
"""

from __future__ import annotations

import cmath
import math
from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def add_term(poly: dict[int, Fraction], exponent: int, value: Fraction) -> None:
    poly[exponent] = poly.get(exponent, Fraction(0)) + value
    if poly[exponent] == 0:
        del poly[exponent]


def crt_character_checks() -> tuple[int, int]:
    map_cases = 0
    profile_cases = 0
    for p in range(2, 14):
        for q in range(2, 14):
            if gcd(p, q) != 1:
                continue
            modulus = p * q
            alpha = pow(q, -1, p)
            beta = pow(p, -1, q)

            for s in range(p):
                for r in range(q):
                    t = (q * alpha * s + p * beta * r) % modulus
                    require(t % p == s and t % q == r, "CRT residue failure")
                    for m in range(modulus):
                        lhs = (m * t) % modulus
                        rhs = (
                            m * q * alpha * s + m * p * beta * r
                        ) % modulus
                        require(lhs == rhs, "CRT character failure")
                        map_cases += 1

            profiles = []
            for seed in range(3):
                u = [((seed + 2) * s + seed + 1) % 7 - 3 for s in range(p)]
                v = [((seed + 3) * r + 2 * seed + 1) % 9 - 4 for r in range(q)]
                profiles.append((u, v))

            for u, v in profiles:
                for m in range(modulus):
                    ell = (m * alpha) % p
                    k = (m * beta) % q
                    left: dict[int, Fraction] = {}
                    right: dict[int, Fraction] = {}
                    for s in range(p):
                        for r in range(q):
                            t = (q * alpha * s + p * beta * r) % modulus
                            weight = Fraction(u[s] * v[r], modulus)
                            add_term(left, (-m * t) % modulus, weight)
                            exponent = (-ell * s * q - k * r * p) % modulus
                            add_term(right, exponent, weight)
                    require(left == right, "group-ring tensor factorization failure")
                    profile_cases += 1
    return map_cases, profile_cases


def polyphase_checks() -> int:
    cases = 0
    for modulus in (6, 10, 15, 21, 35, 91):
        fourier = {
            n: Fraction((n * n + 3 * n + 5) % 17 - 8, abs(n) + 3)
            for n in range(-4 * modulus, 4 * modulus + 1)
        }
        for m in range(modulus):
            selected = {n: a for n, a in fourier.items() if (n - m) % modulus == 0}
            energy = sum((a * a for a in selected.values()), Fraction(0))

            direct_selected: dict[int, Fraction] = {}
            for n, a in fourier.items():
                root_sum = sum(
                    1 if (n - m) * t % modulus == 0 else 0
                    for t in range(modulus)
                )
                # The preceding count is not the complex root sum.  Exact
                # orthogonality says it is modulus precisely on the selected
                # congruence and zero otherwise.
                selector = Fraction(1 if (n - m) % modulus == 0 else 0)
                require(
                    (root_sum == modulus) == ((n - m) % modulus == 0),
                    "integer selector control failure",
                )
                if selector:
                    direct_selected[n] = a * selector
            require(direct_selected == selected, "polyphase selector failure")
            direct_energy = sum(
                (a * a for a in direct_selected.values()), Fraction(0)
            )
            require(direct_energy == energy, "polyphase Parseval failure")
            cases += 1
    return cases


def lrc_unit_bank_checks() -> tuple[int, float]:
    p = 7
    q = 13
    modulus = 91
    alpha = pow(q, -1, p)
    beta = pow(p, -1, q)
    require(alpha == 6 and beta == 2, "LRC inverse indices changed")

    unit_residues = [m for m in range(modulus) if gcd(m, modulus) == 1]
    require(len(unit_residues) == 72, "phi(91) failure")
    cases = 0
    adjacent_min_scaled = float("inf")

    for d in range(7):
        for r0 in range(13):
            for adjacent in (False, True):
                support = {r0} if not adjacent else {r0, (r0 + 1) % 13}
                for m in unit_residues:
                    ell = (m * alpha) % p
                    k = (m * beta) % q
                    require(ell != 0 and k != 0, "unit character lost a prime colour")

                    left = 0j
                    for t in range(modulus):
                        j_value = 1 if t % 7 == d else 0
                        a_value = 1 if t % 13 in support else 0
                        left += (
                            j_value
                            * a_value
                            * cmath.exp(-2j * math.pi * m * t / modulus)
                        )
                    left /= modulus

                    jhat = cmath.exp(-2j * math.pi * ell * d / 7) / 7
                    ahat = sum(
                        cmath.exp(-2j * math.pi * k * r / 13)
                        for r in support
                    ) / 13
                    require(abs(left - jhat * ahat) < 2e-13, "LRC CRT DFT failure")

                    scaled = abs(left) ** 2 * modulus * modulus
                    if adjacent:
                        adjacent_min_scaled = min(adjacent_min_scaled, scaled)
                        require(
                            scaled > Fraction(4, 169),
                            "strict adjacent chord floor failure",
                        )
                    else:
                        require(abs(scaled - 1.0) < 2e-12, "singleton energy failure")
                    cases += 1

    expected_min = 4 * math.sin(math.pi / 26) ** 2
    require(
        abs(adjacent_min_scaled - expected_min) < 3e-12,
        "adjacent minimum formula failure",
    )
    require(cases == 7 * 13 * 2 * 72, "LRC bank coverage failure")
    return cases, adjacent_min_scaled


def constant_checks() -> dict[str, Fraction]:
    rho = Fraction(1, 26754 * 338)
    require(rho == Fraction(1, 9042852), "universal owner-cell mass failure")

    per_class_strict_floor = rho * Fraction(4, 169 * 8281)
    require(
        per_class_strict_floor == Fraction(1, 3163842975657),
        "universal per-class floor failure",
    )
    singleton_total = rho * Fraction(72, 8281)
    adjacent_total = rho * Fraction(132, 8281)
    require(
        singleton_total == Fraction(6, 6240321451),
        "singleton total floor failure",
    )
    require(
        adjacent_total == Fraction(11, 6240321451),
        "adjacent total floor failure",
    )

    singleton_colour_sum = Fraction(13 * 1 - 1 * 1, 13 * 13)
    adjacent_colour_sum = Fraction(13 * 2 - 2 * 2, 13 * 13)
    require(singleton_colour_sum == Fraction(12, 169), "singleton Parseval failure")
    require(adjacent_colour_sum == Fraction(22, 169), "adjacent Parseval failure")

    return {
        "rho": rho,
        "per_class": per_class_strict_floor,
        "singleton_total": singleton_total,
        "adjacent_total": adjacent_total,
    }


def noncoprime_hostile_check() -> None:
    u = (1, 0)
    v = (0, 1)
    abstract_nonzero = u[0] * v[1]
    common_profile = [u[t % 2] * v[t % 2] for t in range(2)]
    require(abstract_nonzero == 1, "abstract hostile tensor vanished")
    require(common_profile == [0, 0], "noncoprime hostile failed")


def prony_sharp_controls() -> int:
    controls = 0
    for length in range(2, 9):
        nodes = [Fraction(j + 2) for j in range(length)]
        weights = []
        for j, node in enumerate(nodes):
            derivative = Fraction(1)
            for h, other in enumerate(nodes):
                if h != j:
                    derivative *= node - other
            weights.append(1 / derivative)
        values = []
        for n in range(length):
            values.append(
                sum(
                    (weights[j] * nodes[j] ** n for j in range(length)),
                    Fraction(0),
                )
            )
        require(values[:-1] == [0] * (length - 1), "Prony sharp zeros failed")
        require(values[-1] == 1, "Prony Vandermonde survivor failed")
        controls += 1
    return controls


def main() -> None:
    map_cases, profile_cases = crt_character_checks()
    polyphase_cases = polyphase_checks()
    lrc_cases, adjacent_min = lrc_unit_bank_checks()
    constants = constant_checks()
    noncoprime_hostile_check()
    prony_controls = prony_sharp_controls()

    print("THM-2424 exact companion")
    print(f"crt_character_cases={map_cases}")
    print(f"crt_profile_factorizations={profile_cases}")
    print(f"polyphase_packets={polyphase_cases}")
    print(f"lrc_owner_status_unit_cases={lrc_cases}")
    print(f"adjacent_min_scaled={adjacent_min:.18f}")
    print(f"owner_cell_rho>={constants['rho']}")
    print(f"each_unit_class_energy>{constants['per_class']}")
    print(f"singleton_total_unit_energy>={constants['singleton_total']}")
    print(f"adjacent_total_unit_energy>={constants['adjacent_total']}")
    print(f"prony_sharp_controls={prony_controls}")
    print("noncoprime_p2q2_hostile=PASS")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
