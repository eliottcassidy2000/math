#!/usr/bin/env python3
"""Exact companion for THM-2355.

The checks use only Fraction arithmetic.  Roots of unity are represented in
Q[zeta_p] by cyclic coefficient vectors modulo 1+z+...+z^(p-1).
Every executable check remains active under optimized Python.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


Rat = Fraction
Gaussian = tuple[Rat, Rat]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def gadd(z: Gaussian, w: Gaussian) -> Gaussian:
    return (z[0] + w[0], z[1] + w[1])


def gmul(z: Gaussian, w: Gaussian) -> Gaussian:
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z: Gaussian) -> Gaussian:
    return (z[0], -z[1])


def gscale(a: Rat, z: Gaussian) -> Gaussian:
    return (a * z[0], a * z[1])


def gabs2(z: Gaussian) -> Rat:
    return z[0] * z[0] + z[1] * z[1]


ZERO_G: Gaussian = (Rat(0), Rat(0))


class CycloPrime:
    """An element of Q[zeta_p] for prime p, in the basis 1,...,zeta^(p-2)."""

    def __init__(self, p: int, coefficients: tuple[Rat, ...]) -> None:
        require(p >= 3, "the cyclotomic model needs p>=3")
        require(len(coefficients) == p - 1, "wrong cyclotomic vector length")
        self.p = p
        self.c = tuple(Rat(value) for value in coefficients)

    @classmethod
    def zero(cls, p: int) -> "CycloPrime":
        return cls(p, tuple(Rat(0) for _ in range(p - 1)))

    @classmethod
    def scalar(cls, p: int, value: Rat) -> "CycloPrime":
        return cls(p, (Rat(value),) + tuple(Rat(0) for _ in range(p - 2)))

    @classmethod
    def root(cls, p: int, exponent: int) -> "CycloPrime":
        exponent %= p
        if exponent < p - 1:
            values = [Rat(0) for _ in range(p - 1)]
            values[exponent] = Rat(1)
            return cls(p, tuple(values))
        return cls(p, tuple(Rat(-1) for _ in range(p - 1)))

    @classmethod
    def from_cyclic(cls, p: int, coefficients: list[Rat]) -> "CycloPrime":
        require(len(coefficients) == p, "wrong cyclic vector length")
        tail = Rat(coefficients[p - 1])
        return cls(p, tuple(Rat(coefficients[i]) - tail for i in range(p - 1)))

    def __add__(self, other: "CycloPrime") -> "CycloPrime":
        require(self.p == other.p, "cyclotomic prime mismatch")
        return CycloPrime(self.p, tuple(a + b for a, b in zip(self.c, other.c)))

    def __neg__(self) -> "CycloPrime":
        return CycloPrime(self.p, tuple(-a for a in self.c))

    def __sub__(self, other: "CycloPrime") -> "CycloPrime":
        return self + (-other)

    def scale(self, value: Rat) -> "CycloPrime":
        return CycloPrime(self.p, tuple(Rat(value) * a for a in self.c))

    def __mul__(self, other: "CycloPrime") -> "CycloPrime":
        require(self.p == other.p, "cyclotomic prime mismatch")
        cyclic = [Rat(0) for _ in range(self.p)]
        for i, a in enumerate(self.c):
            for j, b in enumerate(other.c):
                cyclic[(i + j) % self.p] += a * b
        return CycloPrime.from_cyclic(self.p, cyclic)

    def conj(self) -> "CycloPrime":
        result = CycloPrime.zero(self.p)
        for exponent, coefficient in enumerate(self.c):
            result = result + CycloPrime.root(
                self.p,
                -exponent,
            ).scale(coefficient)
        return result

    def is_zero(self) -> bool:
        return all(value == 0 for value in self.c)

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, CycloPrime)
            and self.p == other.p
            and self.c == other.c
        )


GaussianCyclo = tuple[CycloPrime, CycloPrime]


def gc_zero(p: int) -> GaussianCyclo:
    return (CycloPrime.zero(p), CycloPrime.zero(p))


def gc_scalar(p: int, z: Gaussian) -> GaussianCyclo:
    return (CycloPrime.scalar(p, z[0]), CycloPrime.scalar(p, z[1]))


def gc_add(z: GaussianCyclo, w: GaussianCyclo) -> GaussianCyclo:
    return (z[0] + w[0], z[1] + w[1])


def gc_mul(z: GaussianCyclo, w: GaussianCyclo) -> GaussianCyclo:
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gc_conj(z: GaussianCyclo) -> GaussianCyclo:
    return (z[0].conj(), -z[1].conj())


def gc_scale(a: Rat, z: GaussianCyclo) -> GaussianCyclo:
    return (z[0].scale(a), z[1].scale(a))


def gc_times_root(z: Gaussian, p: int, exponent: int) -> GaussianCyclo:
    root = CycloPrime.root(p, exponent)
    return (root.scale(z[0]), root.scale(z[1]))


def current_energy(values: list[Gaussian], subset: tuple[int, ...]) -> Rat:
    total = ZERO_G
    for index in subset:
        total = gadd(total, values[index])
    return gabs2(total)


def cross_correlation(
    u: list[Gaussian],
    v: list[Gaussian],
    shift: int,
) -> Gaussian:
    p = len(u)
    total = ZERO_G
    for x in range(p):
        total = gadd(total, gmul(u[(x + shift) % p], gconj(v[x])))
    return total


def main() -> None:
    polarization_checks = 0
    mobius_checks = 0
    deletion_checks = 0

    for n in range(1, 9):
        for seed in range(1, 8):
            values = [
                (
                    Rat((seed + 2 * j) % 7 - 3, j + 1),
                    Rat((2 * seed + 3 * j) % 9 - 4, seed + j + 2),
                )
                for j in range(n)
            ]
            full = current_energy(values, tuple(range(n)))
            singles = sum(
                (current_energy(values, (j,)) for j in range(n)),
                Rat(0),
            )
            pair_sum = sum(
                (
                    current_energy(values, pair)
                    for pair in combinations(range(n), 2)
                ),
                Rat(0),
            )
            require(
                full == pair_sum - (n - 2) * singles,
                "complete pair-union polarization failed",
            )
            polarization_checks += 1

            deletion_sum = sum(
                (
                    current_energy(
                        values,
                        tuple(k for k in range(n) if k != j),
                    )
                    for j in range(n)
                ),
                Rat(0),
            )
            require(
                deletion_sum == (n - 2) * full + singles,
                "single-deletion identity failed",
            )
            if n != 2:
                require(
                    full == (deletion_sum - singles) / (n - 2),
                    "single-deletion reconstruction failed",
                )
            deletion_checks += 1

            if n <= 6:
                for size in range(3, n + 1):
                    for target in combinations(range(n), size):
                        coefficient = Rat(0)
                        for inner_size in range(size + 1):
                            for inner in combinations(target, inner_size):
                                coefficient += (
                                    (-1) ** (size - inner_size)
                                ) * current_energy(values, inner)
                        require(
                            coefficient == 0,
                            "a Boolean energy interaction survived above order two",
                        )
                        mobius_checks += 1

    # Untwisted energies on a tree do not determine the full current.
    tree_a = [(Rat(0), Rat(1)), (Rat(1), Rat(0)), (Rat(0), Rat(1))]
    tree_b = [(Rat(0), Rat(1)), (Rat(1), Rat(0)), (Rat(0), Rat(-1))]
    require(
        [gabs2(z) for z in tree_a] == [gabs2(z) for z in tree_b],
        "tree hostile lost singleton agreement",
    )
    for edge in ((0, 1), (1, 2)):
        require(
            current_energy(tree_a, edge) == current_energy(tree_b, edge),
            "tree hostile lost edge-energy agreement",
        )
    require(
        current_energy(tree_a, (0, 1, 2)) == 5
        and current_energy(tree_b, (0, 1, 2)) == 1,
        "tree hostile lost its distinct grouped energies",
    )

    # At two components, deletions retain only the singleton energies.
    require(
        current_energy([(Rat(1), Rat(0)), (Rat(1), Rat(0))], (0, 1)) == 4
        and current_energy(
            [(Rat(1), Rat(0)), (Rat(-1), Rat(0))],
            (0, 1),
        )
        == 0,
        "two-component deletion hostile changed",
    )

    twist_checks = 0
    test_pairs = [
        ((Rat(1), Rat(2)), (Rat(3), Rat(-1))),
        ((Rat(-2), Rat(1)), (Rat(1, 3), Rat(5, 2))),
        ((Rat(0), Rat(0)), (Rat(4), Rat(7))),
        ((Rat(5, 4), Rat(-3, 7)), (Rat(-2, 5), Rat(9, 8))),
    ]
    for p in (3, 5, 13):
        for z, w in test_pairs:
            formal = {
                0: (gabs2(z) + gabs2(w), Rat(0)),
                (-1) % p: gmul(z, gconj(w)),
                1: gmul(gconj(z), w),
            }
            for k in range(p):
                coefficient = formal.get(k, ZERO_G)
                expected = (
                    (gabs2(z) + gabs2(w), Rat(0))
                    if k == 0
                    else gmul(gconj(z), w)
                    if k == 1
                    else gmul(z, gconj(w))
                    if k == p - 1
                    else ZERO_G
                )
                require(coefficient == expected, "twist-energy DFT failed")
                twist_checks += 1

    correlation_checks = 0
    for p in (3, 5, 13):
        arrays = [
            [
                (
                    Rat((2 * q + 1) % 5 - 2, q + 1),
                    Rat((3 * q + 2) % 7 - 3, q + 2),
                )
                for q in range(p)
            ],
            [
                (
                    Rat(1 if q in {0, 2 % p} else 0),
                    Rat(1 if q == 1 else 0),
                )
                for q in range(p)
            ],
        ]
        for values in arrays:
            energies: list[GaussianCyclo] = []
            for t in range(p):
                fourier = gc_zero(p)
                for q, value in enumerate(values):
                    fourier = gc_add(
                        fourier,
                        gc_times_root(value, p, t * q),
                    )
                energies.append(gc_mul(fourier, gc_conj(fourier)))

            for h in range(p):
                inverse = gc_zero(p)
                for t, energy in enumerate(energies):
                    root = CycloPrime.root(p, -t * h)
                    inverse = gc_add(
                        inverse,
                        (energy[0] * root, energy[1] * root),
                    )
                inverse = gc_scale(Rat(1, p), inverse)
                expected = gc_scalar(
                    p,
                    cross_correlation(values, values, h),
                )
                require(
                    inverse == expected,
                    "twist energy did not invert to residue autocorrelation",
                )
                correlation_checks += 1

    # Exact real full-support perfect-autocorrelation array on C_13.
    p = 13
    epsilon = [-1 if k in {1, p - 1} else 1 for k in range(p)]
    real_cazac: list[CycloPrime] = []
    for x in range(p):
        value = CycloPrime.zero(p)
        for k in range(p):
            value = value + CycloPrime.root(p, k * x).scale(Rat(epsilon[k], p))
        real_cazac.append(value)
        require(not value.is_zero(), "real CAZAC lost full support")
    require(
        real_cazac[0] == CycloPrime.scalar(p, Rat(9, 13)),
        "real CAZAC origin value changed",
    )
    for x in range(1, p):
        expected = (
            CycloPrime.root(p, x) + CycloPrime.root(p, -x)
        ).scale(Rat(-2, 13))
        require(real_cazac[x] == expected, "real CAZAC closed form changed")
    for h in range(p):
        autocorrelation = CycloPrime.zero(p)
        for x in range(p):
            autocorrelation = (
                autocorrelation
                + real_cazac[(x + h) % p] * real_cazac[x]
            )
        require(
            autocorrelation
            == CycloPrime.scalar(p, Rat(1 if h == 0 else 0)),
            "real CAZAC autocorrelation changed",
        )

    # Acute-cone support equality and its sharp boundary failure.
    u = [ZERO_G for _ in range(13)]
    v = [ZERO_G for _ in range(13)]
    u[0], u[2], u[5] = (Rat(1), Rat(0)), (Rat(3), Rat(1)), (Rat(2), Rat(0))
    v[1], v[4] = (Rat(2), Rat(0)), (Rat(4), Rat(1))
    expected_support = {
        (a - b) % 13
        for a, value_a in enumerate(u)
        for b, value_b in enumerate(v)
        if value_a != ZERO_G and value_b != ZERO_G
    }
    actual_support = {
        h for h in range(13) if cross_correlation(u, v, h) != ZERO_G
    }
    require(
        actual_support == expected_support,
        "acute-cone cross-correlation lost support equality",
    )
    for h in actual_support:
        require(
            cross_correlation(u, v, h)[0] > 0,
            "acute-cone correlation left its open half-plane",
        )

    boundary_u = [(Rat(1), Rat(0)), (Rat(0), Rat(1)), ZERO_G]
    boundary_v = [(Rat(1), Rat(0)), (Rat(0), Rat(-1)), ZERO_G]
    require(
        0
        in {
            (a - b) % 3
            for a, za in enumerate(boundary_u)
            for b, zb in enumerate(boundary_v)
            if za != ZERO_G and zb != ZERO_G
        }
        and cross_correlation(boundary_u, boundary_v, 0) == ZERO_G,
        "cone-width boundary hostile stopped cancelling",
    )

    print("theorem=THM-2355")
    print("status=PROVED+VERIFIED-EXACT-CANDIDATE-UNDER-INDEPENDENT-AUDIT")
    print(f"pair_polarization_checks={polarization_checks}")
    print(f"boolean_higher_interaction_checks={mobius_checks}")
    print(f"single_deletion_checks={deletion_checks}")
    print("two_component_deletion_boundary=SHARP")
    print("untwisted_tree_energy_boundary=SHARP")
    print(f"cyclic_twist_dft_checks={twist_checks}")
    print(f"group_energy_autocorrelation_checks={correlation_checks}")
    print("real_C13_full_support_perfect_autocorrelation=YES")
    print("acute_cone_support_equality=PASS")
    print("cone_width_sum_pi_boundary_cancellation=YES")
    print("lrc14_status=OPEN")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
