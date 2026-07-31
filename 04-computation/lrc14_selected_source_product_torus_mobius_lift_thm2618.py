#!/usr/bin/env python3
"""Exact referee for THM-2618.

The script is dependency-free and deliberately uses ``require`` rather than
``assert`` so that its normal and ``python -O`` executions have identical
semantics.  All cyclic labels are elements of Z/13Z.
"""

from fractions import Fraction


P = 13
FULL = (1 << P) - 1
TAU = 1
PHI13_BITS = (1 << 12) - 1

CHECKS = 0


def require(condition, message):
    """Raise on a failed exact check, including under ``python -O``."""

    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def bit(mask, r):
    return (mask >> (r % P)) & 1


def cyclic_word(mask, alpha, tau=TAU):
    """The thirteen-bit word read from ``alpha`` in direction ``tau``."""

    value = 0
    for j in range(P):
        value = (value << 1) | bit(mask, alpha + j * tau)
    return value


def selected_packet(mask, tau=TAU):
    """Return (alpha,q,occupied tail,empty head), or None for constants."""

    if mask == 0 or mask == FULL:
        return None
    words = [cyclic_word(mask, alpha, tau) for alpha in range(P)]
    alpha = max(range(P), key=words.__getitem__)
    require(words.count(words[alpha]) == 1, "nonconstant rotation tie")
    q = next(j for j in range(1, P) if bit(mask, alpha + j * tau) == 0)
    tail = (alpha + (q - 1) * tau) % P
    head = (tail + tau) % P
    require(bit(mask, tail) == 1, "selected tail is not occupied")
    require(bit(mask, head) == 0, "selected head is not empty")
    return alpha, q, tail, head


def translate_mask(mask, t):
    """Translate the support L to L+t."""

    out = 0
    for r in range(P):
        if bit(mask, r):
            out |= 1 << ((r + t) % P)
    return out


def mobius_transform(values):
    """Boolean-lattice Mobius transform indexed by subset masks."""

    coeff = values[:]
    for i in range(P):
        step = 1 << i
        for mask in range(1 << P):
            if mask & step:
                coeff[mask] -= coeff[mask ^ step]
    return coeff


def zeta_transform(values):
    """Inverse of ``mobius_transform`` on the Boolean lattice."""

    out = values[:]
    for i in range(P):
        step = 1 << i
        for mask in range(1 << P):
            if mask & step:
                out[mask] += out[mask ^ step]
    return out


def zeta_power_mod_phi13(exponent):
    """X^exponent in F_2[X]/(1+X+...+X^12), in the degree <12 basis."""

    exponent %= P
    if exponent < 12:
        return 1 << exponent
    # X^12 = 1+X+...+X^11 in characteristic two.
    return PHI13_BITS


def cyclotomic_profile_sum(mask, alpha):
    """Sum_{s in mask} zeta^(-alpha*s) in the exact Phi_13 quotient."""

    total = 0
    for s in range(P):
        if bit(mask, s):
            total ^= zeta_power_mod_phi13(-alpha * s)
    return total


def mod_one(x):
    return x - (x.numerator // x.denominator)


def circular_distance(x):
    y = mod_one(x)
    return min(y, 1 - y)


def ordinary_danger(theta):
    return circular_distance(theta) < Fraction(1, 14)


def selector_and_mobius_checks():
    packets = [selected_packet(mask) for mask in range(1 << P)]
    require(packets[0] is None, "zero word must be unselected")
    require(packets[FULL] is None, "one word must be unselected")

    truth = []
    nonconstant_selected = 0
    for rho in range(P):
        row = [int(packets[mask] is not None and packets[mask][2] == rho)
               for mask in range(1 << P)]
        coeff = mobius_transform(row)
        require(coeff[0] == 0, "empty Mobius coefficient must vanish")
        require(zeta_transform(coeff) == row, "Mobius inversion failure")
        truth.append(row)

    for mask in range(1 << P):
        packet = packets[mask]
        selected_count = sum(truth[rho][mask] for rho in range(P))
        if packet is None:
            require(selected_count == 0, "constant word acquired a selector")
            continue
        nonconstant_selected += 1
        require(selected_count == 1, "nonconstant word lacks a unique selector")
        require(truth[packet[2]][mask] == 1, "selector truth-table mismatch")
        require(bit(mask, packet[2]) == 1, "selector violates source occupancy")

        for t in range(P):
            moved = translate_mask(mask, t)
            moved_packet = packets[moved]
            require(moved_packet is not None, "translation made a word constant")
            require(moved_packet[2] == (packet[2] + t) % P,
                    "selected tail is not root-covariant")
            require(moved_packet[3] == (packet[3] + t) % P,
                    "selected head is not root-covariant")

    require(nonconstant_selected == 8190, "wrong nonconstant selector count")
    return len(packets), nonconstant_selected


def profile_checks():
    profile_count = 0
    colour_count = 0
    for mask in range(1, FULL):
        profile_count += 1
        orbit = [translate_mask(mask, t) for t in range(P)]
        require(len(set(orbit)) == P, "proper profile has a translation stabilizer")
        for alpha in range(1, P):
            require(cyclotomic_profile_sum(mask, alpha) != 0,
                    "proper profile lost a primitive cyclotomic colour")
            colour_count += 1

    require(profile_count == 8190, "wrong proper-profile count")
    require(colour_count == 98280, "wrong profile-colour count")

    # A representative atom orbit checks (28)--(31) combinatorially.
    representative = (1 << 0) | (1 << 2) | (1 << 5)
    orbit = [translate_mask(representative, t) for t in range(P)]
    require(orbit[0] == representative, "orbit lost its base atom")
    require(len(set(orbit)) == P, "atom orbit is not regular")
    for t in range(P):
        for u in range(P):
            require((orbit[t] == orbit[u]) == (t == u),
                    "distinct profile atoms are not orthogonal labels")
    for alpha in range(1, P):
        require(cyclotomic_profile_sum(representative, alpha) != 0,
                "representative profile character vanished")

    return profile_count, colour_count


def ordinary_graft_wall_atlas_checks():
    walls = set()
    for s in range(P):
        centre = Fraction(s, P)
        walls.add(mod_one(centre - Fraction(1, 14)))
        walls.add(mod_one(centre + Fraction(1, 14)))
    walls = sorted(walls)
    require(len(walls) == 26, "ordinary-graft wall atlas must have 26 walls")

    chambers = []
    for i, left in enumerate(walls):
        right = walls[(i + 1) % len(walls)]
        if i + 1 == len(walls):
            right += 1
        chambers.append(mod_one((left + right) / 2))

    danger_histogram = {1: 0, 2: 0}
    for theta in walls + chambers:
        dangerous = [s for s in range(P)
                     if ordinary_danger(theta - Fraction(s, P))]
        require(len(dangerous) in (1, 2), "translated grid has wrong danger count")
        danger_histogram[len(dangerous)] += 1
        if not ordinary_danger(theta):
            require(any(s != 0 for s in dangerous),
                    "safe base state has no translated ordinary-graft hole")

    require(set(danger_histogram) == {1, 2}, "danger histogram lost a count class")
    require(all(danger_histogram[n] > 0 for n in (1, 2)),
            "both one-hole and two-hole chambers must occur")
    return len(walls), len(chambers), danger_histogram


def finite_hostile_checks():
    # Equal masses on a regular C_13 orbit: every primitive mean vanishes.
    for alpha in range(1, P):
        full_sum = 0
        for t in range(P):
            full_sum ^= zeta_power_mod_phi13(-alpha * t)
        require(full_sum == 0, "regular-orbit primitive mean did not vanish")

    identity = [[int(x == y) for y in range(P)] for x in range(P)]
    shift = [[int(y == (x + 1) % P) for y in range(P)] for x in range(P)]
    for matrix in (identity, shift):
        require([sum(row) for row in matrix] == [1] * P,
                "permutation square has wrong row marginals")
        require([sum(matrix[x][y] for x in range(P)) for y in range(P)]
                == [1] * P, "permutation square has wrong column marginals")

    identity_diagonal = [identity[s][s] for s in range(P)]
    shift_diagonal = [shift[s][s] for s in range(P)]
    require(identity_diagonal == [1] * P,
            "identity square lost its constant positive diagonal")
    require(shift_diagonal == [0] * P,
            "unit-shift square acquired a diagonal entry")
    for alpha in range(1, P):
        diagonal_colour = 0
        for s, value in enumerate(identity_diagonal):
            if value & 1:
                diagonal_colour ^= zeta_power_mod_phi13(-alpha * s)
        require(diagonal_colour == 0,
                "constant identity diagonal acquired a primitive colour")

    return P - 1, 2


def main():
    word_count, selected_count = selector_and_mobius_checks()
    profile_count, colour_count = profile_checks()
    wall_count, chamber_count, danger_histogram = ordinary_graft_wall_atlas_checks()
    hostile_colours, square_count = finite_hostile_checks()

    print("THM-2618 selected-source product-torus Mobius lift exact referee")
    print(f"boolean_words={word_count}")
    print(f"nonconstant_selected_words={selected_count}")
    print("mobius_truth_tables=13 inverse_reconstructions=13")
    print(f"nonempty_proper_profiles={profile_count}")
    print(f"primitive_profile_coefficients={colour_count}")
    print("cyclotomic_certificate=F2[X]/(1+X+...+X^12)")
    print(f"ordinary_graft_walls={wall_count} chambers={chamber_count}")
    print("ordinary_graft_danger_histogram="
          + ",".join(f"{key}:{danger_histogram[key]}" for key in sorted(danger_histogram)))
    print("profile_orbit_size=13 stabilizer_size=1")
    print(f"equal_mass_primitive_zero_colours={hostile_colours}")
    print(f"permutation_square_controls={square_count}")
    print(f"exact_require_calls={CHECKS}")
    print("status=PASS")


if __name__ == "__main__":
    main()
