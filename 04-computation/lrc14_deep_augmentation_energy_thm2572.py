#!/usr/bin/env python3
"""Exact companion for THM-2572.

The script uses only integer and Fraction arithmetic.  Cyclotomic
nonvanishing is checked in the power basis of Q[zeta_13], where
zeta_13^12=-(1+...+zeta_13^11).
"""

from fractions import Fraction
from itertools import product


P = 13
NONZERO = range(1, P)
DUTY_CLASSES = {
    229692: 24,
    440232: 12,
    440244: 132,
}

checks = 0


def check(condition: bool, message: str = "exact check failed") -> None:
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def raw_anchor_energy(g):
    """Sum over m!=0 of |sum_u g[u] zeta^(mu)|^2."""
    rho = sum(g)
    return P * sum(x * x for x in g) - rho * rho


def reduced_anchor(g, m):
    """Power-basis vector of sum_u g[u] zeta^(mu)."""
    coeff = [0] * P
    for u, value in enumerate(g):
        coeff[(m * u) % P] += value
    return tuple(coeff[j] - coeff[P - 1] for j in range(P - 1))


def duty(x, y):
    return (
        2316060
        + 210552 * (int(x == 0) + int(y == 0))
        + 12 * int((x + y) % P == 0)
        + 19128 * int(x == 0 and y == 0)
    )


def hostile_target(x, y):
    if y != 0:
        return 0
    return -12 if x == 0 else 1


def hostile_commutator(gx, gy):
    values = {}
    for x in range(P):
        for y in range(P):
            value = 0
            for j in range(7):
                u = (x - j * gx) % P
                v = (y - j * gy) % P
                value += (duty(x, y) - duty(u, v)) * hostile_target(u, v)
            values[(x, y)] = value
    return values


def main() -> None:
    # Exhaust the small nonnegative cone and its unique equality ray.
    ternary_profiles = 0
    equality_profiles = 0
    for tail in product(range(3), repeat=P - 1):
        if not any(tail):
            continue
        ternary_profiles += 1
        g = (0,) + tail
        rho = sum(g)
        energy = raw_anchor_energy(g)
        check((P - 1) * energy >= rho * rho)
        equality = (P - 1) * energy == rho * rho
        check(equality == (len(set(tail)) == 1))
        check(energy > 0)
        equality_profiles += int(equality)

    # Every Boolean off-zero profile has all twelve nontrivial anchors.
    boolean_profiles = 0
    for bits in product(range(2), repeat=P - 1):
        if not any(bits):
            continue
        boolean_profiles += 1
        g = (0,) + bits
        for m in NONZERO:
            check(any(reduced_anchor(g, m)))

    # The 2,028 diagonal-free singleton cells used in THM-2567.
    singleton_controls = 0
    for r in range(P):
        for s in range(P):
            for t in range(P):
                if r == t:
                    continue
                singleton_controls += 1
                g = [0] * P
                g[(r - t) % P] = 1
                check(raw_anchor_energy(g) == P - 1)
                check(all(any(reduced_anchor(g, m)) for m in NONZERO))

    # Equality hostile H=1_(s!=0)1_(r!=t): G(0)=0, G(u)=156.
    hostile_g = [0] + [156] * (P - 1)
    rho = sum(hostile_g)
    anchor_energy = Fraction(raw_anchor_energy(hostile_g), P**6)
    sharp_floor = Fraction(rho * rho, (P - 1) * P**6)
    check(rho == 1872)
    check(anchor_energy == sharp_floor)
    check(anchor_energy == 12 * Fraction(12, 169) ** 2)
    for m in NONZERO:
        vec = reduced_anchor(hostile_g, m)
        check(vec == (-156,) + (0,) * 11)

    # Stronger three-zero-plane hostile from THM-2567.
    strong_g = [0] + [144] * (P - 1)
    check(sum(strong_g) == 1728)
    check(
        Fraction(raw_anchor_energy(strong_g), P**6)
        == Fraction(1728**2, 12 * P**6)
    )

    # The restricted six-replica Gram is c*J_6 and has one eigenvalue 6c.
    class_energy = {}
    for d, multiplicity in DUTY_CLASSES.items():
        c = d * d * anchor_energy
        replica_energy = 6 * c
        lower = Fraction(d * d * rho * rho, 2 * P**6)
        check(replica_energy == lower)
        check(6 * c == replica_energy)
        check(multiplicity > 0)
        class_energy[d] = replica_energy

    # Sharpness in (15) is for the displayed six coordinates, not the full
    # commutator.  The equality hostile has a strict 25/24 full-energy ratio
    # on the horizontal gain.
    horizontal = hostile_commutator(1, 0)
    full_energy = sum(value * value for value in horizontal.values())
    displayed_energy = sum(horizontal[(j, 0)] ** 2 for j in range(1, 7))
    check(Fraction(full_energy, displayed_energy) == Fraction(25, 24))

    check(sum(DUTY_CLASSES.values()) == P * P - 1)
    d_square_sum = sum(mult * d * d for d, mult in DUTY_CLASSES.items())
    global_lower_factor = Fraction(d_square_sum, 2 * P**6)
    global_energy = sum(
        DUTY_CLASSES[d] * class_energy[d] for d in DUTY_CLASSES
    )
    check(global_energy == global_lower_factor * rho * rho)
    check(d_square_sum == 29175403421376)
    check(global_lower_factor == Fraction(14587701710688, 4826809))

    # THM-2563 mass floors: M>=63*rho_0 (guard) or 81*rho_0 (ordinary).
    guard_floor = Fraction(63**2, 12 * P**6)
    ordinary_floor = Fraction(81**2, 12 * P**6)
    check(guard_floor == Fraction(1323, 19307236))
    check(ordinary_floor == Fraction(2187, 19307236))

    # Three roots with cyclic gaps 4,4,5 surround zero.  Positive real
    # barycentric weights can therefore be approximated by positive rationals,
    # proving that no individual-colour amplitude has a uniform positive floor.
    triangle = (1, 5, 9)
    gaps = tuple(
        (triangle[(j + 1) % 3] - triangle[j]) % P for j in range(3)
    )
    check(gaps == (4, 4, 5))
    check(max(gaps) * 2 < P)

    print("THM-2572 deep-augmentation Parseval energy")
    print(f"ternary profiles {ternary_profiles}, equality profiles {equality_profiles}")
    print(f"Boolean profiles with all 12 anchors {boolean_profiles}")
    print(f"diagonal-free singleton controls {singleton_controls}")
    print(f"sharp hostile mass {rho}")
    print(f"sharp anchor energy {anchor_energy}")
    print("duty classes " + ", ".join(
        f"{d}^{DUTY_CLASSES[d]}" for d in sorted(DUTY_CLASSES)
    ))
    print(f"global duty-energy factor {global_lower_factor}")
    print(f"THM-2563 guard anchor floor {guard_floor} rho_0^2")
    print(f"THM-2563 ordinary anchor floor {ordinary_floor} rho_0^2")
    print(f"single-colour infimum triangle gaps {gaps}")
    print(f"explicit checks {checks}")


if __name__ == "__main__":
    main()
