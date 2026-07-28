#!/usr/bin/env python3
"""Dependency-free exact referee for the THM-2716 C4 transporter."""

from itertools import product
from math import gcd


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def target(arrow):
    source, degree = arrow
    return (source + degree) % 2


def inverse(arrow):
    return target(arrow), (-arrow[1]) % 4


def compose(first, second):
    check(target(first) == second[0], "noncomposable arrows")
    return first[0], (first[1] + second[1]) % 4


# The action groupoid C4 \ltimes (C4/<B^2>).
arrows = {(source, degree) for source in range(2) for degree in range(4)}
identities = {(0, 0), (1, 0)}
generators = identities | {(0, 1), (0, 3)}
closure = set(generators)
while True:
    old = set(closure)
    closure |= {inverse(arrow) for arrow in old}
    closure |= {
        compose(first, second)
        for first, second in product(old, repeat=2)
        if target(first) == second[0]
    }
    if closure == old:
        break
check(closure == arrows, "two odd arms do not generate all eight arrows")
check(all(inverse(inverse(arrow)) == arrow for arrow in arrows),
      "inverse involution failed")
for arrow in arrows:
    source = arrow[0]
    check(compose((source, 0), arrow) == arrow, "left identity failed")
    check(compose(arrow, (target(arrow), 0)) == arrow, "right identity failed")
for first, second, third in product(arrows, repeat=3):
    if target(first) == second[0] and target(second) == third[0]:
        check(compose(compose(first, second), third)
              == compose(first, compose(second, third)),
              "associativity failed")

hom_degrees = {
    (source, target_object): tuple(
        degree for degree in range(4)
        if target((source, degree)) == target_object
    )
    for source, target_object in product(range(2), repeat=2)
}
check(set(hom_degrees.values()) == {(0, 2), (1, 3)},
      "homogeneous fibres are not the even/odd C2 torsors")
cross = hom_degrees[(0, 1)]
post_j = tuple(compose((0, degree), (1, 2))[1] for degree in cross)
check(cross == (1, 3) and post_j == (3, 1), "target deck swap failed")

# A functorial section of the pair groupoid chooses one odd arm; the inverse
# cross arrow and identities are forced.  Deck translation swaps the choices.
sections = tuple(
    (degree, inverse((0, degree))[1])
    for degree in cross
)
check(sections == ((1, 3), (3, 1)), "pair-groupoid sections failed")
deck_sections = tuple(((forward + 2) % 4, (backward + 2) % 4)
                      for forward, backward in sections)
check(deck_sections == sections[::-1], "deck action does not swap sections")
check(not any(section == deck_sections[index]
              for index, section in enumerate(sections)),
      "an invariant section unexpectedly exists")

# The physical macro and debt cycles have opposite parity types.  Exactly the
# two odd offsets give an equivariant type alignment E->A and G->B.
macro_phases = [4]
debt_phases = [7]
for _ in range(3):
    macro_phases.append((13 * macro_phases[-1]) % 17)
    debt_phases.append((13 * debt_phases[-1]) % 170)
check(tuple(macro_phases) == (4, 1, 13, 16), "macro C4 phases failed")
check(tuple(debt_phases) == (7, 91, 163, 79), "debt C4 phases failed")
macro_types = ("E", "G", "E", "G")
debt_types = ("B", "A", "B", "A")
alignments = tuple(
    offset for offset in range(4)
    if all((macro_types[index] == "E") ==
           (debt_types[(index + offset) % 4] == "A")
           for index in range(4))
)
check(alignments == (1, 3), "macro/debt alignment torsor failed")

# Linearization of B+B^3 is the all-ones incidence matrix.  It acts by 2 on
# the invariant line and by 0 on the sign line.
norm_matrix = ((1, 1), (1, 1))


def matvec(matrix, vector):
    return tuple(sum(row[index] * vector[index] for index in range(len(vector)))
                 for row in matrix)


check(matvec(norm_matrix, (1, 1)) == (2, 2), "norm-two line failed")
check(matvec(norm_matrix, (1, -1)) == (0, 0), "sign-line kernel failed")

# On the independent half sheet, p-epsilon is invariant.  The two odd arms
# toggle epsilon, so their sheet-character sums are +2 and -2.
half_objects = tuple(product(range(2), repeat=2))
half_components = tuple(
    tuple(obj for obj in half_objects if (obj[0] - obj[1]) % 2 == component)
    for component in range(2)
)
check(tuple(map(len, half_components)) == (2, 2), "half-sheet split failed")
half_boundaries = tuple(
    sum((-1) ** (character * (degree % 2)) for degree in cross)
    for character in range(2)
)
check(half_boundaries == (2, -2), "half-character boundaries failed")

# THM-2707 lift translations lie in the odd cyclic group Z/13^6.  Counts of
# n-torsion in a cyclic group are gcd(n,R), hence only zero has order 2 or 4.
R = 13 ** 6
check(gcd(2, R) == gcd(4, R) == 1, "odd-deck obstruction failed")

# In F_13[C7], good-prime Frobenius is basis inversion.  Its nontrivial
# exponent orbits over F_13 are the three inverse pairs; no nonzero or unit
# aggregate can be killed by iterating this automorphism.
frobenius_exponents = tuple((13 * exponent) % 7 for exponent in range(7))
inverse_exponents = tuple((-exponent) % 7 for exponent in range(7))
frobenius_pairs = tuple(tuple(sorted((exponent, (-exponent) % 7)))
                          for exponent in range(1, 4))
check(frobenius_exponents == inverse_exponents,
      "13-Frobenius is not C7 inversion")
check(frobenius_pairs == ((1, 6), (2, 5), (3, 4)),
      "inverse-character pair descent failed")

# Exact THM-2706 displayed-bank census.  Unequal cardinality rules out an
# involution exchanging either pair of full selected banks.
raw_banks = (5850, 4958)
safe_banks = (5850, 4386)
check(raw_banks[0] != raw_banks[1] and safe_banks[0] != safe_banks[1],
      "displayed forward/reverse banks unexpectedly pair")

print("C4 ARM TRANSPORTER EXACT REFEREE")
print(f"objects=2; arrows={len(arrows)}; hom_degrees={hom_degrees}")
print(f"cross_hom={cross}; target_deck_swap={post_j}; sections={sections}")
print(f"macro_phases={tuple(macro_phases)}; debt_phases={tuple(debt_phases)}; alignments={alignments}")
print(f"norm_matrix={norm_matrix}; invariant_image={matvec(norm_matrix, (1, 1))}; sign_image={matvec(norm_matrix, (1, -1))}")
print(f"half_components={half_components}; half_character_boundaries={half_boundaries}")
print(f"thm2707_R={R}; order2_translation_count={gcd(2, R)}; order4_translation_count={gcd(4, R)}")
print(f"F13_C7_frobenius={frobenius_exponents}; inverse_pairs={frobenius_pairs}")
print(f"displayed_transit_banks_raw={raw_banks}; safe={safe_banks}; exchange_involution=false")
print("scope=relative_Z2_degree and finite C4 support/incidence groupoid; no physical current")
print("verdict=PASS")
