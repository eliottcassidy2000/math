#!/usr/bin/env python3
"""Exact finite controls for THM-2365.

The analytic endpoint and BV arguments are proved in the theorem text.  This
companion checks the complete F_13^3 shift geometry, the target-action
projection and finite-difference identity, the sharp line-energy inequality,
the circulant/non-circulant controls, and rational cyclotomic rigidity.
"""

from fractions import Fraction
from itertools import product


P = 13
CELLS = P**3
TARGETS = P**2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def norm_sq(value: tuple[int, int]) -> int:
    """Squared modulus of a Gaussian integer."""
    return value[0] * value[0] + value[1] * value[1]


def add(
    left: tuple[int, int], right: tuple[int, int]
) -> tuple[int, int]:
    return left[0] + right[0], left[1] + right[1]


def neg(value: tuple[int, int]) -> tuple[int, int]:
    return -value[0], -value[1]


def orbit_projection(
    table: dict[tuple[int, int, int], Fraction]
) -> tuple[dict[tuple[int, int, int], Fraction], dict[int, Fraction]]:
    """Project onto H(r,s,t)=h(r-t)."""
    orbit_means: dict[int, Fraction] = {}
    for difference in range(P):
        orbit_means[difference] = (
            sum(
                table[((t + difference) % P, s, t)]
                for s, t in product(range(P), repeat=2)
            )
            / (P * P)
        )
    projected = {
        (r, s, t): orbit_means[(r - t) % P]
        for r, s, t in product(range(P), repeat=3)
    }
    return projected, orbit_means


def explicit_group_average(
    table: dict[tuple[int, int, int], Fraction]
) -> dict[tuple[int, int, int], Fraction]:
    """Average H(r+v,s+u,t+v) over the target action."""
    return {
        (r, s, t): (
            sum(
                table[
                    (
                        (r + v) % P,
                        (s + u) % P,
                        (t + v) % P,
                    )
                ]
                for u, v in product(range(P), repeat=2)
            )
            / (P * P)
        )
        for r, s, t in product(range(P), repeat=3)
    }


def normalized_drift(
    table: dict[tuple[int, int, int], Fraction],
    projected: dict[tuple[int, int, int], Fraction],
) -> Fraction:
    return (
        sum(
            (table[index] - projected[index]) ** 2
            for index in table
        )
        / CELLS
    )


def average_difference_energy(
    table: dict[tuple[int, int, int], Fraction]
) -> Fraction:
    """Return (2|G|)^-1 sum_g ||H-T_g H||_2^2."""
    total = Fraction(0)
    for u, v in product(range(P), repeat=2):
        total += (
            sum(
                (
                    table[(r, s, t)]
                    - table[
                        (
                            (r + v) % P,
                            (s + u) % P,
                            (t + v) % P,
                        )
                    ]
                )
                ** 2
                for r, s, t in product(range(P), repeat=3)
            )
            / CELLS
        )
    return total / (2 * P * P)


def cyclotomic_transform(
    values: tuple[Fraction, ...], colour: int
) -> tuple[Fraction, ...]:
    """Represent sum_r values[r] zeta^(colour*r) in Q[zeta_13]."""
    require(len(values) == P, "wrong cyclotomic vector length")
    coefficients = [Fraction(0) for _ in range(P)]
    for residue, value in enumerate(values):
        coefficients[(colour * residue) % P] += value
    # Reduce zeta^12=-(1+zeta+...+zeta^11).
    return tuple(
        coefficients[index] - coefficients[P - 1]
        for index in range(P - 1)
    )


# A nonnegative circulant hostile with h(0)=0.
profile = tuple(
    Fraction(0) if difference == 0
    else Fraction((difference * difference + 3 * difference + 5) % 19 + 1, 23)
    for difference in range(P)
)
circulant = {
    (r, s, t): profile[(r - t) % P]
    for r, s, t in product(range(P), repeat=3)
}
require(
    all(circulant[(t, s, t)] == 0 for s, t in product(range(P), repeat=2)),
    "circulant hostile lost its diagonal zero",
)
circulant_projection, orbit_means = orbit_projection(circulant)
require(circulant_projection == circulant, "circulant projection changed")
require(tuple(orbit_means[d] for d in range(P)) == profile, "wrong orbit means")
circulant_drift = normalized_drift(circulant, circulant_projection)
require(circulant_drift == 0, "circulant hostile acquired drift")

# Break one off-diagonal cell while preserving nonnegativity and the diagonal.
noncirculant = dict(circulant)
noncirculant[(1, 0, 0)] += 1
require(
    all(
        noncirculant[(t, s, t)] == 0
        for s, t in product(range(P), repeat=2)
    ),
    "noncirculant control changed the diagonal",
)
projected, noncirculant_means = orbit_projection(noncirculant)
explicit = explicit_group_average(noncirculant)
require(projected == explicit, "orbit and group-average projections disagree")
noncirculant_drift = normalized_drift(noncirculant, projected)
require(noncirculant_drift > 0, "noncirculant control has zero drift")
require(
    average_difference_energy(noncirculant) == noncirculant_drift,
    "finite-difference drift identity failed",
)
require(
    any(
        sum(noncirculant[(r, s, t)] for r in range(P))
        != sum(noncirculant[(r, 0, 0)] for r in range(P))
        for s, t in product(range(P), repeat=2)
    ),
    "row-sum sufficient drift control disappeared",
)

# The diagonal-plane line sum is the 2D transform of H(t,s,t), hence zero.
diagonal_cells = 0
for s, t in product(range(P), repeat=2):
    require(noncirculant[(t, s, t)] == 0, "diagonal target line survived")
    diagonal_cells += 1
require(diagonal_cells == P * P, "wrong diagonal-plane size")

# Check the sharp 1/13 energy inequality on every nonzero target line.
line_controls = 0
for first, second in product(range(P), repeat=2):
    if first == second == 0:
        continue
    entries: dict[int, tuple[int, int]] = {}
    running = (0, 0)
    for colour in range(1, P):
        value = (
            ((first + 2) * colour + 3 * second) % 17 - 8,
            ((second + 5) * colour + first) % 19 - 9,
        )
        entries[colour] = value
        running = add(running, value)
    entries[0] = neg(running)
    require(
        sum(value[0] for value in entries.values()) == 0
        and sum(value[1] for value in entries.values()) == 0,
        "target-line cancellation failed",
    )
    total_energy = sum(norm_sq(value) for value in entries.values())
    unit_energy = sum(
        norm_sq(entries[colour]) for colour in range(1, P)
    )
    require(
        P * unit_energy >= total_energy,
        "unit-deep 1/13 energy inequality failed",
    )
    line_controls += 1
require(line_controls == TARGETS - 1, "wrong nonzero target count")

sharp_entries = {0: (-12, 0)}
sharp_entries.update({colour: (1, 0) for colour in range(1, P)})
sharp_total = sum(norm_sq(value) for value in sharp_entries.values())
sharp_unit = sum(
    norm_sq(sharp_entries[colour]) for colour in range(1, P)
)
sharp_ratio = Fraction(sharp_unit, sharp_total)
require(sharp_ratio == Fraction(1, P), "sharp line ratio changed")

# Rational prime-cyclotomic rigidity: every nonconstant rational profile has
# all twelve nontrivial transforms nonzero.
cyclotomic_checks = 0
for seed in range(64):
    values = tuple(
        Fraction(0) if residue == 0
        else Fraction(
            (seed + 3) * (residue + 1) % 29 + 1,
            (seed + residue) % 11 + 1,
        )
        for residue in range(P)
    )
    require(any(values), "cyclotomic control is zero")
    for colour in range(1, P):
        require(
            any(cyclotomic_transform(values, colour)),
            "nonconstant rational profile lost a nonzero colour",
        )
        cyclotomic_checks += 1
require(cyclotomic_checks == 64 * (P - 1), "wrong rigidity census")

profile_colours = sum(
    any(cyclotomic_transform(profile, colour))
    for colour in range(1, P)
)
require(profile_colours == P - 1, "circulant inverse line lost a colour")

# The two BV Fourier factors contribute 1/(4*pi^2), and the two-sided
# Basel sum contributes pi^2/3.
bv_constant = Fraction(1, 4) * Fraction(1, 3)
require(bv_constant == Fraction(1, 12), "BV covariance constant changed")

print("theorem=THM-2365")
print("status=PROVED-CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print(f"target_shift_cells={CELLS}")
print(f"target_action_orbits={P}")
print(f"cells_per_orbit={P * P}")
print(f"diagonal_plane_cells={diagonal_cells}")
print(f"circulant_hostile_drift={circulant_drift}")
print(f"noncirculant_control_drift={noncirculant_drift}")
print("finite_difference_identity=PASS")
print(f"nonzero_target_line_controls={line_controls}")
print(f"sharp_unit_deep_energy_ratio={sharp_ratio}")
print(f"rational_cyclotomic_checks={cyclotomic_checks}")
print(f"circulant_inverse_line_nonzero_colours={profile_colours}")
print(f"bv_covariance_constant={bv_constant}")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
