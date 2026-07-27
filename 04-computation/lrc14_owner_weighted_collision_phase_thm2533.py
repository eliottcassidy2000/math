#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2533.

The checks are deliberately split by type.

* Q(zeta_13) arithmetic verifies the C_13 projector module and the
  owner-multiplier identity without floating point.
* Prime-field controls verify the grouped-jump Prony block bound and every
  consecutive Fourier minor used by the guard argument.
* Rational interval arithmetic verifies the three/four guard-root law and
  the minimal four-mode toothpick stalk.
* A sharp cyclotomic BV fibre is completed to a nonnegative rational
  7-by-13 table.  It satisfies every common stationary unit-safe support
  condition but is not asserted to be a lawful anchored THM-2449 table.

Only the Python standard library is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


P13 = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def valuation(number: int, prime: int) -> int:
    """Return the p-adic valuation of a nonzero integer."""

    require(number != 0, "valuation of zero is not used in this referee")
    answer = 0
    remaining = abs(number)
    while remaining % prime == 0:
        answer += 1
        remaining //= prime
    return answer


def interval_step_fourier_nonzero(frequency: int, denominator: int) -> bool:
    """Nonvanishing test for 1_[0,1/denominator] at an integer frequency."""

    return frequency == 0 or frequency % denominator != 0


@dataclass(frozen=True)
class Cyclo13:
    """An exact element of Q[zeta_13] in the basis 1,zeta,...,zeta^11."""

    coefficients: tuple[Fraction, ...]

    def __post_init__(self) -> None:
        require(len(self.coefficients) == 12, "bad cyclotomic dimension")

    @classmethod
    def zero(cls) -> "Cyclo13":
        return cls((Fraction(0),) * 12)

    @classmethod
    def one(cls) -> "Cyclo13":
        return cls((Fraction(1),) + (Fraction(0),) * 11)

    @classmethod
    def rational(cls, value: int | Fraction) -> "Cyclo13":
        return cls((Fraction(value),) + (Fraction(0),) * 11)

    @classmethod
    def from_cyclic(cls, values: tuple[Fraction, ...]) -> "Cyclo13":
        """Reduce a length-13 cyclic polynomial modulo Phi_13."""

        require(len(values) == 13, "bad cyclic representative")
        tail = values[12]
        return cls(tuple(values[index] - tail for index in range(12)))

    @classmethod
    def zeta_power(cls, exponent: int) -> "Cyclo13":
        values = [Fraction(0)] * 13
        values[exponent % 13] = Fraction(1)
        return cls.from_cyclic(tuple(values))

    def __add__(self, other: "Cyclo13") -> "Cyclo13":
        return Cyclo13(tuple(a + b for a, b in zip(self.coefficients, other.coefficients)))

    def __sub__(self, other: "Cyclo13") -> "Cyclo13":
        return Cyclo13(tuple(a - b for a, b in zip(self.coefficients, other.coefficients)))

    def __neg__(self) -> "Cyclo13":
        return Cyclo13(tuple(-value for value in self.coefficients))

    def __mul__(self, other: "Cyclo13") -> "Cyclo13":
        values = [Fraction(0)] * 13
        for left_index, left in enumerate(self.coefficients):
            for right_index, right in enumerate(other.coefficients):
                values[(left_index + right_index) % 13] += left * right
        return Cyclo13.from_cyclic(tuple(values))

    def scale(self, scalar: int | Fraction) -> "Cyclo13":
        scalar = Fraction(scalar)
        return Cyclo13(tuple(scalar * value for value in self.coefficients))

    def __truediv__(self, scalar: int | Fraction) -> "Cyclo13":
        return self.scale(Fraction(1, 1) / Fraction(scalar))

    def __pow__(self, exponent: int) -> "Cyclo13":
        require(exponent >= 0, "negative cyclotomic power")
        answer = Cyclo13.one()
        base = self
        power = exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power //= 2
        return answer

    def sigma(self, unit: int) -> "Cyclo13":
        require(unit % 13 != 0, "Galois multiplier is not a unit")
        values = [Fraction(0)] * 13
        for exponent, coefficient in enumerate(self.coefficients):
            values[(unit * exponent) % 13] += coefficient
        return Cyclo13.from_cyclic(tuple(values))

    def conjugate(self) -> "Cyclo13":
        return self.sigma(-1)

    def is_zero(self) -> bool:
        return all(value == 0 for value in self.coefficients)


def cyclo_sum(values: tuple[Cyclo13, ...] | list[Cyclo13]) -> Cyclo13:
    answer = Cyclo13.zero()
    for value in values:
        answer = answer + value
    return answer


def translate(vector: tuple[Cyclo13, ...], displacement: int) -> tuple[Cyclo13, ...]:
    require(len(vector) == 13, "translation needs one C_13 fibre")
    return tuple(vector[(index + displacement) % 13] for index in range(13))


def projector(vector: tuple[Cyclo13, ...], frequency: int) -> tuple[Cyclo13, ...]:
    require(len(vector) == 13, "projector needs one C_13 fibre")
    return tuple(
        cyclo_sum(
            [
                Cyclo13.zeta_power(-frequency * shift)
                * vector[(index + shift) % 13]
                for shift in range(13)
            ]
        )
        / 13
        for index in range(13)
    )


def fibre_norm_squared(vector: tuple[Cyclo13, ...]) -> Cyclo13:
    return cyclo_sum([value * value.conjugate() for value in vector])


def polynomial_multiply(
    left: list[Cyclo13], right: list[Cyclo13]
) -> list[Cyclo13]:
    answer = [Cyclo13.zero()] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            answer[left_degree + right_degree] = (
                answer[left_degree + right_degree]
                + left_coefficient * right_coefficient
            )
    return answer


def polynomial_evaluate(polynomial: list[Cyclo13], value: Cyclo13) -> Cyclo13:
    answer = Cyclo13.zero()
    for coefficient in reversed(polynomial):
        answer = answer * value + coefficient
    return answer


def prime_factors(number: int) -> tuple[int, ...]:
    factors = []
    candidate = 2
    remaining = number
    while candidate * candidate <= remaining:
        if remaining % candidate == 0:
            factors.append(candidate)
            while remaining % candidate == 0:
                remaining //= candidate
        candidate += 1
    if remaining > 1:
        factors.append(remaining)
    return tuple(factors)


def primitive_root(prime: int) -> int:
    factors = prime_factors(prime - 1)
    for candidate in range(2, prime):
        if all(pow(candidate, (prime - 1) // factor, prime) != 1 for factor in factors):
            return candidate
    raise AssertionError("primitive root not found")


def root_of_order(prime: int, order: int) -> int:
    require((prime - 1) % order == 0, "requested order is absent")
    root = pow(primitive_root(prime), (prime - 1) // order, prime)
    require(pow(root, order, prime) == 1, "bad root order")
    for divisor in prime_factors(order):
        require(pow(root, order // divisor, prime) != 1, "root order collapsed")
    return root


def determinant_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    require(all(len(row) == len(matrix) for row in matrix), "determinant is not square")
    work = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, len(work)):
            multiplier = work[row][column] * inverse % prime
            for entry in range(column, len(work)):
                work[row][entry] = (
                    work[row][entry] - multiplier * work[column][entry]
                ) % prime
    return determinant % prime


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    if not matrix:
        return 0
    work = [[entry % prime for entry in row] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (candidate for candidate in range(row, len(work)) if work[candidate][column]),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, prime)
        work[row] = [entry * inverse % prime for entry in work[row]]
        for other in range(len(work)):
            if other == row:
                continue
            multiplier = work[other][column]
            work[other] = [
                (left - multiplier * right) % prime
                for left, right in zip(work[other], work[row])
            ]
        row += 1
        if row == len(work):
            break
    return row


def fraction_mod(value: Fraction, prime: int) -> int:
    return value.numerator * pow(value.denominator, -1, prime) % prime


def evaluate_mod(value: Cyclo13, root: int, prime: int) -> int:
    return sum(
        fraction_mod(coefficient, prime) * pow(root, exponent, prime)
        for exponent, coefficient in enumerate(value.coefficients)
    ) % prime


def torus_distance(value: Fraction) -> Fraction:
    residue = value - (value.numerator // value.denominator)
    return min(residue, 1 - residue)


def guard_positions(phase: Fraction) -> tuple[int, ...]:
    return tuple(
        shift
        for shift in range(13)
        if torus_distance(phase + Fraction(shift, 13)) <= Fraction(1, 7)
    )


def danger_roots(speed: int, base: Fraction, radius: Fraction) -> tuple[int, ...]:
    """Dangerous roots x=(base+u)/13 for one integer speed."""

    return tuple(
        shift
        for shift in range(13)
        if torus_distance(Fraction(speed, 13) * (base + shift)) < radius
    )


def cyclically_consecutive(positions: tuple[int, ...]) -> bool:
    values = set(positions)
    return any(
        values == {(start + offset) % 13 for offset in range(len(positions))}
        for start in range(13)
    )


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def main() -> None:
    # ------------------------------------------------------------------
    # Exact C_13 projector module and a Boolean periodic owner multiplier.
    # ------------------------------------------------------------------
    fibres = tuple(
        tuple(
            Cyclo13.rational(((fibre + 2) * (root + 3) + root * root) % 19 - 9)
            for root in range(13)
        )
        for fibre in range(4)
    )
    projector_checks = 0
    for vector in fibres:
        projections = tuple(projector(vector, frequency) for frequency in range(13))
        require(
            tuple(cyclo_sum([projection[root] for projection in projections]) for root in range(13))
            == vector,
            "projectors do not reconstruct the fibre",
        )
        projector_checks += 13
        for left in range(13):
            for right in range(13):
                composed = projector(projections[right], left)
                expected = projections[left] if left == right else (Cyclo13.zero(),) * 13
                require(composed == expected, "projectors are not orthogonal idempotents")
                projector_checks += 1
            for displacement in range(13):
                expected = tuple(
                    Cyclo13.zeta_power(left * displacement) * value
                    for value in projections[left]
                )
                require(
                    translate(projections[left], displacement) == expected,
                    "projected fibre has the wrong translation character",
                )
                projector_checks += 1

    owner = (1, 0, 1, 1)
    multiplier_checks = 0
    for fibre, vector in enumerate(fibres):
        weighted = tuple(value.scale(owner[fibre]) for value in vector)
        for frequency in range(13):
            left = projector(weighted, frequency)
            right = tuple(
                value.scale(owner[fibre]) for value in projector(vector, frequency)
            )
            require(left == right, "periodic owner does not commute with projector")
            require(
                fibre_norm_squared(left)
                == fibre_norm_squared(projector(vector, frequency)).scale(owner[fibre]),
                "Boolean owner square identity failed",
            )
            multiplier_checks += 2

    one_branch = tuple(
        Cyclo13.rational(int(root == 4)) for root in range(13)
    )
    eventwise_colours = tuple(
        frequency
        for frequency in range(1, 13)
        if not fibre_norm_squared(projector(one_branch, frequency)).is_zero()
    )
    require(eventwise_colours == tuple(range(1, 13)), "one branch lost a root colour")

    # ---------------------------------------------------------------
    # Grouped-jump Prony: exact block nonsingularity and sharp controls.
    # ---------------------------------------------------------------
    prime = 911  # 910=2*5*7*13.
    root_7 = root_of_order(prime, 7)
    prony_vandermonde_checks = 0
    node_bank = tuple(pow(root_7, index, prime) for index in range(7))
    for width in range(1, 7):
        for nodes in combinations(node_bank, width):
            for start in range(7):
                matrix = tuple(
                    tuple(pow(node, start + row, prime) for node in nodes)
                    for row in range(width)
                )
                require(determinant_mod(matrix, prime) != 0, "Prony block vanished")
                prony_vandermonde_checks += 1

    root_5 = root_of_order(prime, 5)
    sharp_nodes = tuple(pow(root_5, index, prime) for index in range(5))
    sharp_amplitudes = sharp_nodes
    sharp_sequence = tuple(
        sum(amplitude * pow(node, height, prime) for amplitude, node in zip(sharp_amplitudes, sharp_nodes))
        % prime
        for height in range(5)
    )
    require(sharp_sequence[:4] == (0, 0, 0, 0), "sharp Prony zeros disappeared")
    require(sharp_sequence[4] != 0, "sharp Prony endpoint vanished")

    raw_hostile = ((root_5, 1), (root_5, -1))
    grouped = {}
    for node, amplitude in raw_hostile:
        grouped[node] = (grouped.get(node, 0) + amplitude) % prime
    grouped = {node: amplitude for node, amplitude in grouped.items() if amplitude}
    require(not grouped, "coincident grouped jumps failed to cancel")
    require(
        all(sum(amplitude * pow(node, height, prime) for node, amplitude in raw_hostile) % prime == 0 for height in range(12)),
        "raw-jump cancellation hostile is not identically zero",
    )

    prony_frequency_bounds = 0
    for colour in range(1, 13):
        for grouped_nodes in range(1, 21):
            require(
                colour + 13 * (grouped_nodes - 1) <= 13 * grouped_nodes - 1,
                "Prony frequency bound failed",
            )
            prony_frequency_bounds += 1

    # -------------------------------------------------------------------
    # THM-2349 bookkeeping with the future-owner mark retained in f=W.
    # -------------------------------------------------------------------
    # The theorem is invoked with e=F and f=W=gF.  Pointwise Boolean
    # multiplication gives 0 <= W <= F, so every f-cell is an e-cell, while
    # the extra factor fhat(X) in the conclusion remembers which endpoint is
    # future-owned.  The interval steps below are only a finite control for
    # that marked/unmarked distinction; the D_c/D_d support geometry remains
    # the inherited hypothesis of THM-2349.
    carrier_cells = tuple(range(12))
    F_cells = {cell for cell in carrier_cells if cell <= 8}
    future_owner_cells = {cell for cell in carrier_cells if cell % 3 != 2}
    W_cells = F_cells & future_owner_cells
    require(W_cells, "future-owner restriction became empty")
    require(W_cells <= F_cells, "W=gF is not a substep of F")
    require(
        all(int(cell in W_cells) == int(cell in F_cells) * int(cell in future_owner_cells)
            for cell in carrier_cells),
        "Boolean product W=gF failed pointwise",
    )

    c3 = 13 * 13
    F_denominator = 97
    W_denominator = 2 * F_denominator
    require(Fraction(1, W_denominator) <= Fraction(1, F_denominator), "W support escaped F")
    marked_restart_checks = 0
    for colour in range(1, 13):
        X = 13 * colour
        Y = X + c3
        multiplier = (Y - X) // c3
        require(Y == X + multiplier * c3, "restart edge is not a c3 edge")
        require(gcd(multiplier, 91) == 1, "restart edge is not a 91-unit")
        require(valuation(X, 13) == valuation(Y, 13) == 1, "restart left grade one")
        require(X // 13 % 13 == colour, "marked endpoint has the wrong root residue")
        require(Y // 13 % 13 == colour, "other endpoint has the wrong root residue")
        require(
            interval_step_fourier_nonzero(X, W_denominator),
            "future-owner-marked Fourier atom vanished",
        )
        require(
            interval_step_fourier_nonzero(X, F_denominator)
            and interval_step_fourier_nonzero(Y, F_denominator),
            "carrier Fourier atom vanished",
        )
        require(multiplier % 7 and multiplier % 13, "deep-comb factor vanished")
        marked_restart_checks += 1

    # If one forgets fhat(X), a 91-unit edge can be present without touching
    # the future-owned atom.  At positions 0,1,14 on one c3-line, 0--1 is a
    # unit edge, whereas the marked position 14 meets them by multipliers 14
    # and 13.  Thus the unmarked conclusion cannot recover future ownership.
    hostile_positions = (0, 1, 14)
    hostile_marked = {14}
    hostile_unit_edges = tuple(
        (left, right)
        for left, right in combinations(hostile_positions, 2)
        if gcd(abs(right - left), 91) == 1
    )
    hostile_marked_unit_edges = tuple(
        edge for edge in hostile_unit_edges if hostile_marked.intersection(edge)
    )
    require(hostile_unit_edges == ((0, 1),), "hostile unit-edge pattern changed")
    require(not hostile_marked_unit_edges, "hostile edge accidentally retained the mark")

    # ----------------------------------------------------------------
    # Guard blocks, consecutive Fourier minors, and the k=4 stalk split.
    # ----------------------------------------------------------------
    guard_counts = {3: 0, 4: 0}
    guard_phase_checks = 0
    phase_bank = []
    for denominator in (1, 2, 7, 13, 14, 26, 91, 182):
        for numerator in range(denominator):
            phase = Fraction(numerator, denominator)
            positions = guard_positions(phase)
            require(len(positions) in (3, 4), "guard has the wrong root count")
            require(cyclically_consecutive(positions), "guard roots are not consecutive")
            base_is_danger = torus_distance(13 * phase) < Fraction(1, 7)
            require(
                (len(positions) == 3) == base_is_danger,
                "guard count disagrees with the inverse-root law",
            )
            guard_counts[len(positions)] += 1
            guard_phase_checks += 1
            phase_bank.append((phase, positions))

    root_13 = root_of_order(prime, 13)
    consecutive_minor_checks = 0
    for width in (1, 2, 3, 4):
        for frequencies in combinations(range(13), width):
            for start in range(13):
                matrix = tuple(
                    tuple(
                        pow(root_13, frequency * ((start + row) % 13), prime)
                        for frequency in frequencies
                    )
                    for row in range(width)
                )
                require(
                    determinant_mod(matrix, prime) != 0,
                    "consecutive Fourier minor vanished",
                )
                consecutive_minor_checks += 1

    stalk_rank_checks = 0
    for frequencies in combinations(range(13), 4):
        for start in range(13):
            three_rows = tuple(
                tuple(
                    pow(root_13, frequency * ((start + row) % 13), prime)
                    for frequency in frequencies
                )
                for row in range(3)
            )
            require(rank_mod(three_rows, prime) == 3, "three-row stalk rank failed")
            stalk_rank_checks += 1

    guard_safe_and_next_danger = Fraction(2, 7) * Fraction(10, 13)
    guard_safe_and_next_safe = Fraction(5, 7) * Fraction(9, 13)
    require(guard_safe_and_next_danger == Fraction(20, 91), "bad stalk measure")
    require(guard_safe_and_next_safe == Fraction(45, 91), "bad transverse measure")
    require(
        guard_safe_and_next_danger + guard_safe_and_next_safe == Fraction(5, 7),
        "guard stalks do not partition the safe carrier",
    )

    # ------------------------------------------------------------
    # Three distinct gains give three disjoint 72-element orbits.
    # ------------------------------------------------------------
    gain_triple_checks = 0
    for gains in combinations(range(1, 13), 3):
        slices = []
        for gain in gains:
            orbit = {
                (source, target, gain * target % 13)
                for source in range(1, 7)
                for target in range(1, 13)
            }
            require(len(orbit) == 72, "fixed-gain orbit has the wrong size")
            slices.append(orbit)
        require(len(set().union(*slices)) == 216, "three gain slices overlap")
        require(
            all(
                collision * pow(target, -1, 13) % 13 in gains
                for source, target, collision in set().union(*slices)
            ),
            "gain ratio changed inside an orbit",
        )
        gain_triple_checks += 1

    # --------------------------------------------------------------------
    # Sharp common-safe four-mode cyclotomic BV fibre and rational table.
    # --------------------------------------------------------------------
    guard_speed = 5
    epsilon = Fraction(1, 10000)
    expected_guard = tuple(
        shift
        for shift in range(13)
        if (guard_speed * shift) % 13 in (0, 1, 12)
    )
    common_safe_checks = 0
    for base in (-epsilon, Fraction(0), epsilon):
        require(
            danger_roots(guard_speed, base, Fraction(1, 7)) == expected_guard,
            "guard word is not stable on the hostile interval",
        )
        common_safe_checks += 1
        for speed in range(1, 13):
            require(
                danger_roots(speed, base, Fraction(1, 14)) == (0,),
                "ordinary unit danger is not the common singleton",
            )
            common_safe_checks += 1

    polynomial = [Cyclo13.one()]
    for root_exponent in (-1, 0, 1):
        polynomial = polynomial_multiply(
            polynomial,
            [-Cyclo13.zeta_power(root_exponent), Cyclo13.one()],
        )
    require(len(polynomial) == 4, "sharp polynomial has the wrong degree")
    require(all(not coefficient.is_zero() for coefficient in polynomial), "sharp coefficient vanished")

    hostile_fibre = tuple(
        polynomial_evaluate(
            polynomial,
            Cyclo13.zeta_power(guard_speed * root),
        )
        for root in range(13)
    )
    hostile_zero_set = tuple(root for root, value in enumerate(hostile_fibre) if value.is_zero())
    require(hostile_zero_set == expected_guard, "sharp fibre misses the guard block")
    require(len(hostile_zero_set) == 3, "sharp fibre has extra zeros")

    hostile_projections = tuple(projector(hostile_fibre, frequency) for frequency in range(13))
    hostile_active = tuple(
        frequency
        for frequency, projection in enumerate(hostile_projections)
        if any(not value.is_zero() for value in projection)
    )
    expected_active = tuple(sorted({0, guard_speed, 2 * guard_speed % 13, 3 * guard_speed % 13}))
    require(hostile_active == expected_active, "sharp fibre does not have four modes")
    require(0 in hostile_active, "sharp fibre lost its mean")
    hostile_gains = tuple(frequency for frequency in hostile_active if frequency)
    require(len(hostile_gains) == 3, "sharp fibre does not have three gains")

    safe_roots = tuple(root for root in range(13) if root not in expected_guard)
    signed_coefficients = []
    for root in range(13):
        coefficients = list(hostile_fibre[root].coefficients) + [Fraction(0)]
        signed_coefficients.append(tuple(Fraction(91) * value for value in coefficients))
    worst_negative = max(
        [-value for root in safe_roots for value in signed_coefficients[root]]
        + [Fraction(0)]
    )
    background = ceil_fraction(worst_negative) + 1

    table = {}
    for root in range(13):
        safe = root in safe_roots
        for source in range(7):
            for target in range(13):
                value = Fraction(background if safe else 0)
                if source == 0:
                    value += signed_coefficients[root][target]
                require(value >= 0, "rational hostile table is not nonnegative")
                table[(root, source, target)] = value

    root_7_table = root_of_order(prime, 7)
    inverse_91 = pow(91, -1, prime)
    rational_table_checks = 0
    table_channels = {}
    for source_frequency in range(1, 7):
        for target_frequency in range(1, 13):
            for root in range(13):
                transformed = 0
                for source in range(7):
                    for target in range(13):
                        transformed += (
                            fraction_mod(table[(root, source, target)], prime)
                            * pow(root_7_table, source_frequency * source, prime)
                            * pow(root_13, target_frequency * target, prime)
                        )
                transformed = transformed * inverse_91 % prime
                expected = evaluate_mod(
                    hostile_fibre[root].sigma(target_frequency), root_13, prime
                )
                require(transformed == expected, "rational table Galois completion failed")
                table_channels[(source_frequency, target_frequency, root)] = transformed
                rational_table_checks += 1

    hostile_modes = set()
    inverse_13 = pow(13, -1, prime)
    for source_frequency in range(1, 7):
        for target_frequency in range(1, 13):
            channel = tuple(
                table_channels[(source_frequency, target_frequency, root)]
                for root in range(13)
            )
            active = []
            for collision_frequency in range(13):
                coefficient = sum(
                    pow(root_13, -collision_frequency * root, prime) * channel[root]
                    for root in range(13)
                ) * inverse_13 % prime
                if coefficient:
                    active.append(collision_frequency)
            expected = sorted(
                {
                    0,
                    target_frequency * guard_speed % 13,
                    target_frequency * 2 * guard_speed % 13,
                    target_frequency * 3 * guard_speed % 13,
                }
            )
            require(active == expected, "rational table changed the four-mode support")
            for collision_frequency in active:
                if collision_frequency:
                    hostile_modes.add(
                        (source_frequency, target_frequency, collision_frequency)
                    )
            rational_table_checks += 13
    require(len(hostile_modes) == 216, "sharp rational table lost a gain slice")
    require(
        {
            collision * pow(target, -1, 13) % 13
            for source, target, collision in hostile_modes
        }
        == set(hostile_gains),
        "sharp rational table changed its gain set",
    )

    maximum_table_entry = max(table.values())
    require(all(value.denominator == 1 for value in table.values()), "table is not integral")

    print("THM-2533 owner-weighted collision-phase referee")
    print(
        f"projector-module identities={projector_checks}; "
        f"Boolean owner identities={multiplier_checks}; eventwise colours={len(eventwise_colours)}"
    )
    print(
        f"grouped-jump Prony minors={prony_vandermonde_checks}; "
        f"frequency bounds={prony_frequency_bounds}; sharp sequence={sharp_sequence}; "
        "coincident-jump hostile groups to L=0"
    )
    print(
        f"THM-2349 future-owner restart residues={marked_restart_checks}; "
        f"marked 91-unit c3 edges={marked_restart_checks}; "
        f"unmarked-only hostile edges={len(hostile_unit_edges)}"
    )
    print(
        f"guard phases={guard_phase_checks} counts={guard_counts}; "
        f"consecutive Fourier minors={consecutive_minor_checks}; "
        f"four-mode stalk ranks={stalk_rank_checks}"
    )
    print(
        f"minimal k=4 stalk measures: C_H&T^-1E_H={guard_safe_and_next_danger}, "
        f"C_H&T^-1C_H={guard_safe_and_next_safe}"
    )
    print(
        f"gain triples={gain_triple_checks}; each gives 3*72=216 modes; "
        f"common-safe interval checks={common_safe_checks}"
    )
    print(
        f"sharp cyclotomic fibre: zeros={hostile_zero_set}; active={hostile_active}; "
        f"nonzero gains={hostile_gains}"
    )
    print(
        f"rational-table completion checks={rational_table_checks}; background={background}; "
        f"max_entry={maximum_table_entry}; exact mixed modes={len(hostile_modes)}"
    )
    print(
        "VERIFIED: periodic Boolean owners commute with root projectors; grouped jumps give "
        "the sharp Prony block bound; the inherited THM-2349 restart retains a future-owner "
        "mark on every nonzero root residue; guard support forces at least four root modes "
        "and three gain slices, while equality is confined to the one-step guard stalk and "
        "is sharp even for a nonnegative rational common-safe table (not a lawful anchored "
        "table)."
    )


if __name__ == "__main__":
    main()
