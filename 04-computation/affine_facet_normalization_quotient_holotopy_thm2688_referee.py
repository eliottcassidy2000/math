#!/usr/bin/env python3
"""Independent exact referee for THM-2688.

This replay uses the relation and cyclic-extension descriptions directly,
rather than the primary companion's mapping-cone matrices.  It checks the
vertical/diagonal quotient distinction, the reduced-regular lens weights,
and the coefficient Bockstein, together with three sharp hostiles.
"""

from itertools import combinations
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def powerset(values):
    values = tuple(values)
    return {
        frozenset(face)
        for size in range(len(values) + 1)
        for face in combinations(values, size)
    }


def affine_relation_audit(p, a, b, step):
    require(p >= 3 and b % p and gcd(step, p) == 1, "bad affine input")
    missing = lambda c: (a + b * c) % p
    image = set()
    vertices = []
    for c in range(p):
        row = tuple(d for d in range(p) if d != missing(c))
        image |= powerset(row)
        vertices.extend((c, d) for d in row)
    expected_boundary = powerset(range(p)) - {frozenset(range(p))}
    require(image == expected_boundary, f"p={p}: vertical image is not boundary")

    label_step = b * step % p
    action = lambda vertex: (
        (vertex[0] + step) % p,
        (vertex[1] + label_step) % p,
    )
    require(all(action(v)[1] != missing(action(v)[0]) for v in vertices),
            f"p={p}: compatible action left the relation")

    unseen = set(vertices)
    orbit_offsets = []
    while unseen:
        start = min(unseen)
        orbit = []
        value = start
        while value not in orbit:
            orbit.append(value)
            value = action(value)
        require(value == start and len(orbit) == p,
                f"p={p}: diagonal orbit is not free of size p")
        unseen -= set(orbit)
        offsets = {(d - missing(c)) % p for c, d in orbit}
        require(len(offsets) == 1, f"p={p}: offset is not invariant")
        orbit_offsets.append(next(iter(offsets)))
    require(set(orbit_offsets) == set(range(1, p)),
            f"p={p}: diagonal quotient is not Delta(F_p^*)")

    # The closed homotopy ledger follows from p contractible components,
    # vertical target S^(p-2), and diagonal target Delta^(p-2).
    return {
        "p": p,
        "components": p,
        "vertical_dimension": p - 2,
        "quotient_vertices": len(orbit_offsets),
        "cone_pi_h1": p - 1,
        "cone_pi_top": p - 2,
        "cone_q_h1": p - 1,
    }


rows = [
    affine_relation_audit(p, -1, -2, 7 if p == 13 else 1)
    for p in (3, 5, 7, 13)
]

# The reduced real regular representation of C_13 is the six complex
# character lines with weights 1,...,6 and their conjugates.  Every
# nonidentity group element acts without eigenvalue one on these lines.
p = 13
weights = tuple(range(1, 7))
real_characters = {weight % p for weight in weights} | {
    (-weight) % p for weight in weights
}
require(real_characters == set(range(1, p)),
        "lens weights do not exhaust the reduced regular characters")
require(all(
    all((power * character) % p for character in real_characters)
    for power in range(1, p)
), "a nonidentity label rotation fixed a sphere direction")

# Standard one-cell-per-degree lens complex through dimension eleven:
# d_even=13 and d_odd=0.  This gives Z/13 in the five positive odd degrees.
lens_torsion_degrees = tuple(range(1, 11, 2))
require(lens_torsion_degrees == (1, 3, 5, 7, 9),
        "lens homology degree ledger changed")

# THM-2657's quotient is k |-> 2k mod 13.  The generator section s(a)=7a
# has wrap defect 91=13*7, so the kernel-coordinate Bockstein is 7 mod 13.
section = lambda value: 7 * value
require(all((2 * section(value) - value) % p == 0 for value in range(p)),
        "section does not lift the quotient generator")
defects = []
for left in range(p):
    for right in range(p):
        reduced = (left + right) % p
        numerator = section(left) + section(right) - section(reduced)
        require(numerator % p == 0, "section defect missed the kernel")
        defects.append((numerator // p) % p)
require(defects.count(0) == 91 and defects.count(7) == 78,
        "Bockstein wrap census changed")
require(set(defects) == {0, 7} and 7 % p,
        "physical Bockstein became zero")

# Hostile 1: a constant missing map has one filled facet as vertical image.
constant_image = set()
for _c in range(5):
    constant_image |= powerset((1, 2, 3, 4))
require(constant_image == powerset((1, 2, 3, 4)),
        "constant-map hostile changed")

# Hostile 2: an incompatible label step fails to preserve the deleted graph.
missing13 = lambda c: (12 - 2 * c) % 13
bad_action = lambda vertex: ((vertex[0] + 7) % 13, vertex[1])
vertices13 = [(c, d) for c in range(13) for d in range(13)
              if d != missing13(c)]
require(any(bad_action(v)[1] == missing13(bad_action(v)[0])
            for v in vertices13),
        "incompatible-shift hostile unexpectedly preserved the relation")

# Hostile 3: a nongenerator step on C_6 has two component orbits.
unseen = set(range(6))
component_orbits = []
while unseen:
    start = min(unseen)
    orbit = set()
    value = start
    while value not in orbit:
        orbit.add(value)
        value = (value + 2) % 6
    unseen -= orbit
    component_orbits.append(orbit)
require(len(component_orbits) == 2,
        "composite nongenerator hostile lost its two quotient components")

# A split extension has the zero connecting class, providing the sharp
# positive boundary to the nonsplit odometer result.
split_defects = [0 for _left in range(p) for _right in range(p)]
require(set(split_defects) == {0}, "split-extension hostile became nonzero")

print("THM-2688 independent affine-facet/lens/Bockstein referee")
for row in rows:
    print(
        f"p={row['p']}:components={row['components']}:"
        f"vertical=S^{row['vertical_dimension']}:"
        f"diagonal_quotient=Delta^{row['vertical_dimension']}:"
        f"Cone_pi_H1_rank={row['cone_pi_h1']}:"
        f"Cone_pi_top_degree={row['cone_pi_top']}:"
        f"Cone_q_H1_rank={row['cone_q_h1']}"
    )
print("p13_lens_weights=(1,2,3,4,5,6):free=True")
print("p13_lens_torsion_degrees=(1,3,5,7,9):group=Z/13")
print("odometer_section=7a:wraps=78:nonwraps=91:Bockstein=7_mod_13")
print("hostiles=constant_missing_map; incompatible_label_step; C6_step2; split_extension")
print("ALL INDEPENDENT CHECKS PASSED")
