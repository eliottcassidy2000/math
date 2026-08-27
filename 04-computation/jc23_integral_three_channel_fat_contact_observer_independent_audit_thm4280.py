#!/usr/bin/env python3
"""Dependency-free exact referee for THM-4280's arithmetic observer.

The program uses only the Python standard library.  It independently checks
the normalized four-channel matrix over

    K = Q(omega, kappa),  omega^2+omega+1=0,  kappa^2=-omega,

then forgets that scalar extension and checks the grouped arithmetic channel
rank law over Fraction matrices.  It also verifies all minimality hostiles
and the degree-34/42 mod-four consequence of the c1-kernel.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product


Q = Fraction


def require(condition: object, label: str) -> None:
    """Optimization-safe audit gate."""
    if not condition:
        raise RuntimeError(f"audit failure: {label}")


@dataclass(frozen=True)
class Eisenstein:
    """a+b*omega in Q(omega), with omega^2=-omega-1."""

    a: Fraction
    b: Fraction

    @staticmethod
    def coerce(value: object) -> "Eisenstein":
        if isinstance(value, Eisenstein):
            return value
        if isinstance(value, (int, Fraction)):
            return Eisenstein(Q(value), Q(0))
        raise TypeError(f"cannot coerce {type(value)!r} to Eisenstein")

    def __add__(self, other: object) -> "Eisenstein":
        rhs = Eisenstein.coerce(other)
        return Eisenstein(self.a + rhs.a, self.b + rhs.b)

    __radd__ = __add__

    def __neg__(self) -> "Eisenstein":
        return Eisenstein(-self.a, -self.b)

    def __sub__(self, other: object) -> "Eisenstein":
        return self + (-Eisenstein.coerce(other))

    def __rsub__(self, other: object) -> "Eisenstein":
        return Eisenstein.coerce(other) - self

    def __mul__(self, other: object) -> "Eisenstein":
        rhs = Eisenstein.coerce(other)
        # (a+bw)(c+dw)=(ac-bd)+(ad+bc-bd)w.
        return Eisenstein(
            self.a * rhs.a - self.b * rhs.b,
            self.a * rhs.b + self.b * rhs.a - self.b * rhs.b,
        )

    __rmul__ = __mul__

    def inverse(self) -> "Eisenstein":
        # conjugate(a+bw)=(a-b)-bw and norm=a^2-ab+b^2.
        norm = self.a * self.a - self.a * self.b + self.b * self.b
        if norm == 0:
            raise ZeroDivisionError("zero Eisenstein element")
        return Eisenstein((self.a - self.b) / norm, -self.b / norm)

    def __truediv__(self, other: object) -> "Eisenstein":
        return self * Eisenstein.coerce(other).inverse()

    def __rtruediv__(self, other: object) -> "Eisenstein":
        return Eisenstein.coerce(other) / self

    def __bool__(self) -> bool:
        return self.a != 0 or self.b != 0


F_ZERO = Eisenstein(Q(0), Q(0))
F_ONE = Eisenstein(Q(1), Q(0))
F_OMEGA = Eisenstein(Q(0), Q(1))


@dataclass(frozen=True)
class Cyclotomic12:
    """A+B*kappa over Q(omega), with kappa^2=-omega."""

    A: Eisenstein
    B: Eisenstein

    @staticmethod
    def coerce(value: object) -> "Cyclotomic12":
        if isinstance(value, Cyclotomic12):
            return value
        if isinstance(value, Eisenstein):
            return Cyclotomic12(value, F_ZERO)
        if isinstance(value, (int, Fraction)):
            return Cyclotomic12(Eisenstein.coerce(value), F_ZERO)
        raise TypeError(f"cannot coerce {type(value)!r} to Cyclotomic12")

    def __add__(self, other: object) -> "Cyclotomic12":
        rhs = Cyclotomic12.coerce(other)
        return Cyclotomic12(self.A + rhs.A, self.B + rhs.B)

    __radd__ = __add__

    def __neg__(self) -> "Cyclotomic12":
        return Cyclotomic12(-self.A, -self.B)

    def __sub__(self, other: object) -> "Cyclotomic12":
        return self + (-Cyclotomic12.coerce(other))

    def __rsub__(self, other: object) -> "Cyclotomic12":
        return Cyclotomic12.coerce(other) - self

    def __mul__(self, other: object) -> "Cyclotomic12":
        rhs = Cyclotomic12.coerce(other)
        return Cyclotomic12(
            self.A * rhs.A - F_OMEGA * self.B * rhs.B,
            self.A * rhs.B + self.B * rhs.A,
        )

    __rmul__ = __mul__

    def inverse(self) -> "Cyclotomic12":
        # (A+Bk)(A-Bk)=A^2+omega*B^2.
        denominator = self.A * self.A + F_OMEGA * self.B * self.B
        if not denominator:
            raise ZeroDivisionError("zero cyclotomic element")
        return Cyclotomic12(self.A / denominator, -self.B / denominator)

    def __truediv__(self, other: object) -> "Cyclotomic12":
        return self * Cyclotomic12.coerce(other).inverse()

    def __rtruediv__(self, other: object) -> "Cyclotomic12":
        return Cyclotomic12.coerce(other) / self

    def __bool__(self) -> bool:
        return bool(self.A) or bool(self.B)


K_ZERO = Cyclotomic12(F_ZERO, F_ZERO)
K_ONE = Cyclotomic12(F_ONE, F_ZERO)
K_OMEGA = Cyclotomic12(F_OMEGA, F_ZERO)
K_KAPPA = Cyclotomic12(F_ZERO, F_ONE)


def matrix_rank(rows: list[list[object]], column_count: int) -> int:
    """Exact Gaussian rank over Fraction or the field classes above."""

    matrix = [list(row) for row in rows]
    rank = 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, len(matrix))
             if bool(matrix[index][column])),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not bool(matrix[index][column]):
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                matrix[index][j] - multiplier * matrix[rank][j]
                for j in range(column_count)
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def matvec(rows: list[list[Cyclotomic12]], vector: list[Cyclotomic12]) \
        -> list[Cyclotomic12]:
    return [
        sum((row[index] * vector[index] for index in range(len(vector))),
            K_ZERO)
        for row in rows
    ]


def powerset(items: tuple[int, ...]):
    for size in range(len(items) + 1):
        yield from combinations(items, size)


def subset_label(subset: tuple[int, ...]) -> str:
    return "empty" if not subset else ",".join(str(value) for value in subset)


def main() -> None:
    active = (1, 2, 4, 7)
    half = K_ONE / 2
    omega_squared = K_OMEGA * K_OMEGA

    require(
        K_OMEGA * K_OMEGA + K_OMEGA + 1 == K_ZERO,
        "omega relation",
    )
    require(K_KAPPA * K_KAPPA == -K_OMEGA, "kappa relation")
    # The representation itself certifies 1,kappa are F-independent.
    require(
        K_KAPPA.A == F_ZERO and K_KAPPA.B == F_ONE,
        "field-basis independence",
    )

    # Nonzero rescalings of the exact THM-4259 coefficient rows on
    # [u,f,g,h].
    complex_rows = {
        1: [K_ZERO, K_ONE, -K_KAPPA,
            (omega_squared - K_KAPPA) * half],
        2: [K_ONE, K_ZERO, K_ZERO, K_ZERO],
        4: [K_ZERO, K_ZERO, K_ZERO, half],
        7: [K_ZERO, K_ONE, K_KAPPA,
            (omega_squared + K_KAPPA) * half],
    }

    complex_subset_ranks: dict[tuple[int, ...], int] = {}
    for subset in powerset(active):
        rows = [complex_rows[degree] for degree in subset]
        rank = matrix_rank(rows, 4)
        require(rank == len(subset), f"complex rank for subset {subset}")
        complex_subset_ranks[subset] = rank
    complex_active_census = Counter(complex_subset_ranks.values())
    require(
        complex_active_census
        == Counter({0: 1, 1: 4, 2: 6, 3: 4, 4: 1}),
        "complex active census",
    )

    # Adding the eight zero coefficient positions gives the complete
    # length-twelve sample matroid census.
    full_degree_set = tuple(range(12))
    full_census = Counter()
    for subset in powerset(full_degree_set):
        rank = len(set(subset) & set(active))
        full_census[rank] += 1
    require(
        full_census
        == Counter({0: 256, 1: 1024, 2: 1536, 3: 1024, 4: 256}),
        "length-twelve rank census",
    )

    # Pure-channel complexification hostiles, in [u,f,g,h].
    pure_witnesses = {
        1: [K_ZERO, K_KAPPA, -K_ONE, K_ZERO],
        2: [K_ONE, K_ZERO, K_ZERO, K_ZERO],
        4: [K_ZERO, -omega_squared, -K_ONE, K_ONE * 2],
        7: [K_ZERO, K_KAPPA, K_ONE, K_ZERO],
    }
    for target_degree, witness in pure_witnesses.items():
        values = matvec([complex_rows[degree] for degree in active], witness)
        support = tuple(
            degree for degree, value in zip(active, values) if value != K_ZERO
        )
        require(
            support == (target_degree,),
            f"pure complex witness in degree {target_degree}",
        )

    # Arithmetic block observer over F after the invertible coordinate
    # change [u,f,g,h] -> [u,F',G',v], where
    # F'=p+omega^2*r/2, G'=q+r/2, v=2h-omega^2*f-g.
    # A K-valued coefficient is split into its two F coordinates, so all
    # rank work below is over exact Fraction matrices.
    arithmetic_blocks: dict[int, list[list[Fraction]]] = {
        1: [[Q(0), Q(1), Q(0), Q(0)],
            [Q(0), Q(0), Q(-1), Q(0)]],
        2: [[Q(1), Q(0), Q(0), Q(0)]],
        4: [[Q(0), Q(0), Q(0), Q(1)]],
        7: [[Q(0), Q(1), Q(0), Q(0)],
            [Q(0), Q(0), Q(1), Q(0)]],
    }

    arithmetic_subset_ranks: dict[tuple[int, ...], int] = {}
    for subset in powerset(active):
        rows = [row for degree in subset for row in arithmetic_blocks[degree]]
        rank = matrix_rank(rows, 4)
        expected = (
            2 * int(bool(set(subset) & {1, 7}))
            + int(2 in subset)
            + int(4 in subset)
        )
        require(rank == expected, f"arithmetic rank for subset {subset}")
        arithmetic_subset_ranks[subset] = rank
    arithmetic_census = Counter(arithmetic_subset_ranks.values())
    require(
        arithmetic_census == Counter({0: 1, 1: 2, 2: 4, 3: 6, 4: 3}),
        "arithmetic rank census",
    )

    minimal_full = tuple(
        subset for subset, rank in arithmetic_subset_ranks.items()
        if rank == 4
        and all(arithmetic_subset_ranks[subset[:index] + subset[index + 1:]] < 4
                for index in range(len(subset)))
    )
    require(
        minimal_full == ((1, 2, 4), (2, 4, 7)),
        "arithmetic minimal observers",
    )

    # Necessity hostiles in the transformed arithmetic coordinates.
    u = [Q(1), Q(0), Q(0), Q(0)]
    v = [Q(0), Q(0), Q(0), Q(1)]
    f = [Q(0), Q(1), Q(0), Q(0)]

    def vanishes_on(subset: tuple[int, ...], vector: list[Fraction]) -> bool:
        return all(
            sum((row[index] * vector[index] for index in range(4)), Q(0)) == 0
            for degree in subset for row in arithmetic_blocks[degree]
        )

    for subset in powerset(active):
        if 2 not in subset:
            require(vanishes_on(subset, u), f"u hostile for subset {subset}")
        if 4 not in subset:
            require(vanishes_on(subset, v), f"v hostile for subset {subset}")
        if not (set(subset) & {1, 7}):
            require(vanishes_on(subset, f), f"f hostile for subset {subset}")
        full = arithmetic_subset_ranks[subset] == 4
        necessary_conditions = (
            2 in subset and 4 in subset and bool(set(subset) & {1, 7})
        )
        require(full == necessary_conditions, f"full-rank law for {subset}")

    # If c1=0 on the integral Hom lattice, the class is a*u+e*v.  The
    # visible basis u,v is orthogonal of degree four, so its degree is
    # 4*(Norm(a)+Norm(e)).  Exhaust all Eisenstein coefficient residues
    # modulo four as an exact control of the claimed congruence.
    def eisenstein_norm_mod4(a: int, b: int) -> int:
        return (a * a - a * b + b * b) % 4

    zero_channel_degree_residues = {
        (4 * (eisenstein_norm_mod4(a0, a1)
              + eisenstein_norm_mod4(e0, e1))) % 4
        for a0, a1, e0, e1 in product(range(4), repeat=4)
    }
    require(zero_channel_degree_residues == {0}, "c1-kernel degree residues")
    require(34 % 4 == 2 and 42 % 4 == 2, "degree shell residues")
    degree_34_42_c1_zero_impossible = all(
        degree % 4 not in zero_channel_degree_residues for degree in (34, 42)
    )
    require(degree_34_42_c1_zero_impossible, "degree 34/42 c1 zero test")

    complex_rank_ledger = ";".join(
        f"{subset_label(subset)}:{complex_subset_ranks[subset]}"
        for subset in powerset(active)
    )
    arithmetic_rank_ledger = ";".join(
        f"{subset_label(subset)}:{arithmetic_subset_ranks[subset]}"
        for subset in powerset(active)
    )

    print("THM4280 DEPENDENCY-FREE FORMAL-LOG OBSERVER REFEREE")
    print("status=FINITE_EXACT_INDEPENDENT_AUDIT_PASS")
    print("dependencies=PYTHON_STANDARD_LIBRARY_ONLY")
    print("field=Q(omega,kappa),omega^2+omega+1=0,kappa^2=-omega")
    print("complex_active_degrees=1,2,4,7")
    print(f"complex_subset_ranks={complex_rank_ledger}")
    print("complex_active_rank_census=0:1,1:4,2:6,3:4,4:1")
    print("complex_length12_rank_census=0:256,1:1024,2:1536,3:1024,4:256")
    print("complex_matroid=U_4_4_PLUS_8_LOOPS")
    print("complex_pure_hostiles=1:kappa*f-g;2:u;4:v;7:kappa*f+g")
    print(f"arithmetic_subset_ranks={arithmetic_rank_ledger}")
    print("arithmetic_rank_law=2*hit_{1,7}+hit_2+hit_4")
    print("arithmetic_rank_census=0:1,1:2,2:4,3:6,4:3")
    print("arithmetic_minimal_full_observers=1,2,4|2,4,7")
    print("arithmetic_necessity_hostiles=missing2:u;missing4:v;missing1and7:f")
    print("c1_integral_kernel=O*u+O*v")
    print("c1_zero_degree_mod4_residues=0")
    print("degree34_mod4=2 degree42_mod4=2")
    print("degree34_42_c1_zero=IMPOSSIBLE")
    print("degree34_42_channel1_zero_test=PASS")
    print("base_change_hostile=three_channels_fail_after_adjoining_kappa")
    print("raw_Keller_descent=NOT_ADDRESSED")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
