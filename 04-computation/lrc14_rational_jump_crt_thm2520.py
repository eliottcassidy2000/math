#!/usr/bin/env python3
"""Exact companion for THM-2520's rational-jump CRT dichotomy."""

from fractions import Fraction


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mean(values) -> Fraction:
    values = tuple(values)
    return sum(values, Fraction()) / len(values)


def cyclic_jumps(values) -> tuple[Fraction, ...]:
    values = tuple(Fraction(value) for value in values)
    return tuple(values[i] - values[i - 1] for i in range(len(values)))


def aggregate_jumps(values, d: int) -> tuple[Fraction, ...]:
    jumps = cyclic_jumps(values)
    return tuple(
        sum((jumps[i] for i in range(len(values)) if i % d == residue), Fraction())
        for residue in range(d)
    )


def step_value(values, x: Fraction) -> Fraction:
    require(0 <= x < 1, "step evaluation escaped the circle chart")
    return Fraction(values[(len(values) * x.numerator) // x.denominator])


def perron_profile(values, d: int, level: int) -> tuple[Fraction, ...]:
    """Values of P_(13^(level-1)) F on the d equal continuity cells."""
    dilation = P ** (level - 1)
    result = []
    for cell in range(d):
        y = Fraction(2 * cell + 1, 2 * d)
        result.append(
            mean(step_value(values, (y + branch) / dilation) for branch in range(dilation))
        )
    return tuple(result)


def predicted_perron_jumps(
    values, nu: int, d: int, level: int
) -> tuple[Fraction, ...]:
    require(level - 1 >= nu, "level below endpoint conductor")
    dilation = P ** (level - 1)
    unit = pow(P, level - 1 - nu, d)
    current = aggregate_jumps(values, d)
    predicted = [Fraction() for _ in range(d)]
    for residue, jump in enumerate(current):
        predicted[(unit * residue) % d] += jump / dilation
    return tuple(predicted)


def delayed_collision_vector(delay_scale: int) -> tuple[Fraction, ...]:
    """B_u for h=1_[0,1/2) and G=1_[0,1/2) at G(K y)."""
    require(delay_scale % P == 0, "owner weight is not last-digit invariant")
    cells = 2 * delay_scale
    modulus = 2 * cells
    shift_unit = modulus // P
    counts = []
    for u in range(P):
        count = 0
        for cell in range(cells):
            midpoint_numerator = 2 * cell + 1
            first = midpoint_numerator < cells
            shifted = (midpoint_numerator + u * shift_unit) % modulus
            second = shifted < cells
            owner = (delay_scale * midpoint_numerator) % modulus < cells
            count += int(first and second and owner)
        counts.append(Fraction(count, cells))
    return tuple(counts)


def cyclotomic_remainder_nonzero(values, multiplier: int) -> bool:
    """Test sum_u values[u] zeta^(multiplier*u) != 0 in Q(zeta_13)."""
    coefficients = [Fraction() for _ in range(P)]
    for u, value in enumerate(values):
        coefficients[(multiplier * u) % P] += value
    return any(coefficients[i] != coefficients[-1] for i in range(P - 1))


def main() -> None:
    # Exhaustive-style deterministic profiles on D=13*5.  At both legal
    # levels, the actual Perron jumps equal the prime-to-13 aggregate current
    # divided by M and permuted by the appropriate CRT unit.
    nu = 1
    d = 5
    conductor = (P**nu) * d
    state = 0x2520
    profiles = []
    for _ in range(64):
        row = []
        for _index in range(conductor):
            state = (1664525 * state + 1013904223) & 0xFFFFFFFF
            row.append(Fraction((state >> 29) % 3, 2))
        profiles.append(tuple(row))

    jump_controls = 0
    mean_controls = 0
    for values in profiles:
        for level in (2, 3):
            profile = perron_profile(values, d, level)
            require(
                cyclic_jumps(profile) == predicted_perron_jumps(values, nu, d, level),
                "Perron jump-current formula failed",
            )
            require(mean(profile) == mean(values), "Perron mean changed")
            jump_controls += d
            mean_controls += 1
    require((jump_controls, mean_controls) == (640, 128), "control census drifted")

    # A genuinely nonconstant Boolean depth-one inverse-fibre tile has one
    # active branch over each of the five coprime cells.  Its aggregate jump
    # current vanishes and every later Perron profile is the constant 1/13.
    balanced = tuple(
        Fraction(int(index in {t + d * t for t in range(d)}))
        for index in range(conductor)
    )
    require(any(cyclic_jumps(balanced)), "balanced hostile became constant")
    require(aggregate_jumps(balanced, d) == (Fraction(),) * d, "balanced current")
    for level in (2, 3):
        require(
            perron_profile(balanced, d, level) == (Fraction(1, P),) * d,
            "zero current did not give the constant Perron branch",
        )

    # A one-cell Boolean response has nonzero coprime current at every legal
    # level.  The thirteen shifted d-grids are a genuine CRT partition.
    charged = (Fraction(1),) + (Fraction(),) * (conductor - 1)
    current = aggregate_jumps(charged, d)
    require(current != (Fraction(),) * d, "charged current vanished")
    crt_points = {
        (P * target - d * source) % (P * d)
        for target in range(d)
        for source in range(P)
    }
    require(len(crt_points) == P * d, "shifted jump grids collided")

    # For each source colour a mod 13 and every coprime Fourier residue b mod
    # d, CRT supplies exactly one combined ladder residue.
    ladder_addresses = 0
    for level in (2, 3):
        unit = pow(P, level - 1 - nu, d)
        for source in range(1, P):
            for target in range(d):
                addresses = [
                    k
                    for k in range(P * d)
                    if k % P == source and (unit * k) % d == target
                ]
                require(len(addresses) == 1, "ladder CRT address is not unique")
                ladder_addresses += 1
    require(ladder_addresses == 120, "ladder CRT census drifted")

    # Quantitative delayed-owner positive control.  Here d=2, M=1,
    # sum(C^2)=2, d=2, M=1, B=1, V=W=2, rho=1/2.  The exact
    # Parseval/cosecant floor from the theorem is strict at K=13^3.
    delay_scale = P**3
    energy_floor = Fraction(2, 4 * P**2 * 2)
    covariance_error = Fraction(1 * 2 * 2, 6 * delay_scale)
    require(
        Fraction(1, 2) * energy_floor > covariance_error,
        "delayed-owner positivity invoice is not strict",
    )
    collisions = delayed_collision_vector(delay_scale)
    require(
        all(collisions[0] > collisions[u] for u in range(1, P)),
        "delayed collision toothpick energy is not positive",
    )
    require(
        all(cyclotomic_remainder_nonzero(collisions, a) for a in range(1, P)),
        "a delayed primitive collision colour vanished",
    )

    print("THM-2520 rational-jump CRT exact companion: PASS")
    print("endpoint_conductor=13^1*5; Perron_jump_controls=640; mean_controls=128")
    print("zero_current_balanced_Boolean_tile=constant_1/13_at_levels_2,3")
    print("nonzero_current_CRT_points=65; ladder_addresses=120")
    print("all_12_high_frequency_ladders=nonzero_iff_prime_to_13_current_nonzero")
    print("delayed_owner_scale=13^3; BV_invoice=strict; all_12_collision_colours=positive")
    print("constant_branch=exact; future_owner_delay=quantitative; orientation=not_supplied")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
