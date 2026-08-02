#!/usr/bin/env python3
"""Exact companion for THM-3084's cyclic-root/Kummer eigenline theorem."""

from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13
MODULUS = 53
ZETA = 16


def inverse(value):
    return pow(value, MODULUS - 2, MODULUS)


def shift(vector, amount=1):
    amount %= len(vector)
    return tuple(vector[(index - amount) % len(vector)] for index in range(len(vector)))


def scalar_multiple(left, right):
    """Return the scalar when left=lambda*right over F_53, or None."""
    require(all(value % MODULUS for value in right), "torus vector required")
    scalar = left[0] * inverse(right[0]) % MODULUS
    return scalar if all(a % MODULUS == scalar * b % MODULUS for a, b in zip(left, right)) else None


def norm(vector):
    answer = 1
    for value in vector:
        answer = answer * value % MODULUS
    return answer


require(pow(ZETA, P, MODULUS) == 1, "zeta order divides 13")
require(all(pow(ZETA, exponent, MODULUS) != 1 for exponent in range(1, P)), "zeta primitive")


# Every nontrivial character line identifies the cyclic root torsor with the
# diagonal mu_13 fibre, and the equal-weight norm is c^13.
character_cells = 0
for character in range(1, P):
    c = character + 1
    vector = tuple(c * pow(ZETA, character * index, MODULUS) % MODULUS for index in range(P))
    require(norm(vector) == pow(c, P, MODULUS), "character-line norm")
    for amount in range(P):
        expected = pow(ZETA, -character * amount, MODULUS)
        moved = shift(vector, amount)
        require(moved == tuple(expected * value % MODULUS for value in vector), "shift eigenvalue")
        require(scalar_multiple(moved, vector) == expected, "fibre scalar")
        character_cells += 1
require(character_cells == 156, "character-cell count")


# Positive {1,2}^13 profiles align only on the two constant vectors.
positive_aligned = 0
positive_failed = 0
for vector in product((1, 2), repeat=P):
    scalar = scalar_multiple(shift(vector), vector)
    if scalar is None:
        positive_failed += 1
    else:
        require(scalar == 1 and len(set(vector)) == 1, "positive eigenline must be constant")
        positive_aligned += 1
require((positive_aligned, positive_failed) == (2, 8190), "positive-cone census")


# Deterministic generic torus profiles leave their augmentation fibre.
generic_hostiles = 0
for seed in range(100):
    exceptional = 2 + seed % 51
    vector = (exceptional,) + (1,) * (P - 1)
    require(scalar_multiple(shift(vector), vector) is None, "generic non-eigen hostile")
    require(norm(shift(vector)) == norm(vector), "norm remains shift invariant")
    generic_hostiles += 1


# A fixed singleton root, projected first to a nonzero character, gives the
# same character-line phase pointer.  Do not sum roots before this check.
singleton_cells = 0
for root in range(P):
    for character in range(1, P):
        coefficient = pow(ZETA, -character * root, MODULUS)
        require(coefficient != 0, "singleton charged coefficient")
        next_root_coefficient = pow(ZETA, -character * ((root + 1) % P), MODULUS)
        require(
            next_root_coefficient == pow(ZETA, -character, MODULUS) * coefficient % MODULUS,
            "singleton pointer phase",
        )
        singleton_cells += 1
require(singleton_cells == 156, "singleton-cell count")


# V4 needs two independent characters; either one alone has two-point fibres.
v4 = tuple(product((0, 1), repeat=2))
two_character_codes = {(x, y): ((-1) ** x, (-1) ** y) for x, y in v4}
require(len(set(two_character_codes.values())) == 4, "two V4 characters reconstruct")
for coordinate in (0, 1):
    fibres = {}
    for point in v4:
        value = two_character_codes[point][coordinate]
        fibres.setdefault(value, []).append(point)
    require(sorted(map(len, fibres.values())) == [2, 2], "single V4 character fibre")


print("THM-3084 CYCLIC ROOT POINTER AND KUMMER EIGENLINE")
print("field=F53 zeta13=16 order=13")
print(f"character_shift_norm_cells={character_cells}")
print("criterion=F(Px)=F(x) iff Px=lambda*x iff x_j=c*lambda^(-j)")
print(f"positive_binary_profiles=8192 aligned={positive_aligned} failed={positive_failed}")
print(f"generic_non_eigen_hostiles={generic_hostiles} norm_still_invariant=PASS")
print(f"singleton_pointer_character_cells={singleton_cells}")
print("V4=two_characters_bijective;single_character_fibres=2+2")
print("LRC=point_fixed_THM2537_selector_before_character_projection")
print("boundary=cross_packet_transport;zeros;composite_or_unequal_weight")
print("all_exact_checks=PASS")
