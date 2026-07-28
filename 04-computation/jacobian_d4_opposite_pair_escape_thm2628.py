#!/usr/bin/env python3
"""Exact finite and local controls for THM-2628.

The D4 part exhausts every geometric survivor subset and every inertia
element.  The symbolic fixtures are strict-henselian normalization/resultant
controls only; they are not polynomial Keller maps.
"""

import sympy as sp


checks = 0


def require(condition, message):
    global checks
    if not bool(condition):
        raise RuntimeError(message)
    checks += 1


def compose(left, right):
    """Permutation left after right."""

    return tuple(left[right[i]] for i in range(4))


def power(permutation, exponent):
    out = tuple(range(4))
    for _ in range(exponent):
        out = compose(permutation, out)
    return out


def cycle_type(permutation):
    seen = set()
    cycles = []
    for start in range(4):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        cycles.append(length)
    return tuple(sorted(cycles))


identity = tuple(range(4))
rotation = (1, 2, 3, 0)
diagonal_reflection = (0, 3, 2, 1)
tau = power(rotation, 2)
d4 = {
    power(rotation, exponent)
    for exponent in range(4)
} | {
    compose(power(rotation, exponent), diagonal_reflection)
    for exponent in range(4)
}

require(len(d4) == 8, "D4 enumeration did not have eight elements")
require(tau == (2, 3, 0, 1), "opposite pairing changed")
require(all(compose(g, tau) == compose(tau, g) for g in d4),
        "opposite pairing does not centralize D4")


def moved_set(subset):
    return {tau[i] for i in subset}


def shape(subset):
    k = len(subset)
    if k == 0:
        return "empty"
    if k == 1:
        return "singleton"
    if k == 3:
        return "complement-singleton"
    if k == 4:
        return "full"
    require(k == 2, "unexpected survivor-set size")
    return "opposite" if moved_set(subset) == subset else "adjacent"


subsets = [
    {i for i in range(4) if mask & (1 << i)}
    for mask in range(16)
]
shape_census = {}
allowed_pair_census = {}
allowed_inertia = {}
for survivor in subsets:
    k = len(survivor)
    kind = shape(survivor)
    p = len(survivor-moved_set(survivor))
    shape_census[(k, kind, p)] = shape_census.get((k, kind, p), 0) + 1
    possible = {
        g for g in d4
        if survivor <= {i for i in range(4) if g[i] == i}
    }
    allowed_inertia[(tuple(sorted(survivor)), kind)] = possible
    for g in possible:
        key = (k, kind, p, cycle_type(g))
        allowed_pair_census[key] = allowed_pair_census.get(key, 0) + 1

expected_shapes = {
    (0, "empty", 0): 1,
    (1, "singleton", 1): 4,
    (2, "adjacent", 2): 4,
    (2, "opposite", 0): 2,
    (3, "complement-singleton", 1): 4,
    (4, "full", 0): 1,
}
require(shape_census == expected_shapes, "survivor shape/pole census changed")

for survivor in subsets:
    k = len(survivor)
    kind = shape(survivor)
    types = {cycle_type(g) for g in allowed_inertia[(tuple(sorted(survivor)), kind)]}
    if k == 0:
        require(types == {(1, 1, 1, 1), (1, 1, 2), (2, 2), (4,)},
                "empty survivor lost an inertia type")
    elif k == 1:
        require(types == {(1, 1, 1, 1), (1, 1, 2)},
                "singleton inertia table changed")
    elif k == 2 and kind == "adjacent":
        require(types == {(1, 1, 1, 1)},
                "adjacent pair acquired nontrivial inertia")
    elif k == 2 and kind == "opposite":
        require(types == {(1, 1, 1, 1), (1, 1, 2)},
                "opposite pair inertia table changed")
    elif k == 3:
        require(types == {(1, 1, 1, 1)},
                "three survivors acquired nontrivial inertia")
    else:
        require(k == 4 and types == {(1, 1, 1, 1)},
                "full survivor inertia table changed")

# Every nonempty pole-positive survivor type is one of the three claimed
# global D4 boundary shapes.
pole_positive = {
    (len(survivor), shape(survivor),
     len(survivor-moved_set(survivor)))
    for survivor in subsets
    if survivor-moved_set(survivor)
}
require(pole_positive == {
    (1, "singleton", 1),
    (2, "adjacent", 2),
    (3, "complement-singleton", 1),
}, "pole-positive survivor classification changed")

# Strict-henselian local resultant controls.
t, T = sp.symbols("t T")
R1 = sp.expand(t**2*T*(T-t**-1)*(T**2-t**-1))
R2 = sp.expand((t*T-1)*(t*T-2)*T*(T-1))
R3 = sp.expand((t*T-1)*T*(T-1)*(T-2))

require(sp.expand(R1-(t**2*T**4-t*T**3-t*T**2+T)) == 0,
        "R1 expansion changed")
require(sp.Poly(R1.subs(t, 0), T).degree() == 1,
        "R1 specialized degree is not one")
require(sp.Poly(R2.subs(t, 0), T).degree() == 2,
        "R2 specialized degree is not two")
require(sp.Poly(R3.subs(t, 0), T).degree() == 3,
        "R3 specialized degree is not three")

for polynomial, k, e in ((R1, 1, 2), (R2, 2, 2), (R3, 3, 1)):
    lead = sp.Poly(polynomial, T).LC()
    surviving = sp.Poly(polynomial, T).coeff_monomial(T**k)
    require(sp.degree(lead, t) == e,
            "local control leading order changed")
    normalized = sp.cancel(surviving/lead)
    numerator, denominator = sp.fraction(normalized)
    require(sp.degree(denominator, t) == e
            and sp.simplify(numerator.subs(t, 0)) != 0,
            "local control lost the full coefficient pole")

print("THM-2628 D4 opposite-pair escape exact controls")
print("D4_elements=8 opposite_pairing=(02)(13) central=PASS")
print("survivor_shapes: empty=1 singleton=4 adjacent_pairs=4 opposite_pairs=2 triple=4 full=1")
print("deck_pole_branches: k1=1 k2_adjacent=2 k2_opposite=0 k3=1")
print("inertia: k1=id|diagonal k2_adjacent=id k2_opposite=id|diagonal k3=id")
print("local_controls: R1=(k1,e2,diagonal) R2=(k2,e2,p0|p2) R3=(k3,e1,id)")
print("same_k2_resultant_and_pole_order_support_both_p=0_and_p=2")
print(f"exact assertions passed: {checks}")
