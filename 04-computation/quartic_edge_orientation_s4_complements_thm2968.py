"""Exact companion for THM-2968.

Enumerate the edge and oriented-four-cycle actions of S4 in one common
signed-pair frame, verify their pointwise central-sign twist, decompose the
two six-dimensional permutation characters and their fifth exterior powers,
and check the THM-2769 specialization of the THM-2864 discriminant formulas.
"""

from itertools import permutations

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


S4 = list(permutations(range(4)))


def compose_perm(g, h):
    return tuple(g[h[i]] for i in range(len(g)))


def parity(p):
    return sum(
        p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))
    ) % 2


def cycle_type(p):
    used = set()
    parts = []
    for i in range(len(p)):
        if i in used:
            continue
        j = i
        length = 0
        while j not in used:
            used.add(j)
            length += 1
            j = p[j]
        parts.append(length)
    return tuple(sorted(parts, reverse=True))


def rotate_to_zero(cycle):
    j = cycle.index(0)
    return cycle[j:] + cycle[:j]


def inverse_cycle(cycle):
    return rotate_to_zero((cycle[0], cycle[3], cycle[2], cycle[1]))


def relabel_cycle(g, cycle):
    return rotate_to_zero(tuple(g[i] for i in cycle))


def relabel_edge(g, edge):
    return tuple(sorted((g[edge[0]], g[edge[1]])))


def signed_pair_action(pairs, relabel, g):
    lookup = {}
    for j, pair in enumerate(pairs):
        lookup[pair[0]] = (j, 0)
        lookup[pair[1]] = (j, 1)
    require(len(lookup) == 6, "signed-pair labels are not six distinct objects")
    sigma = [None] * 3
    delta = [None] * 3
    for i, pair in enumerate(pairs):
        target0 = relabel(g, pair[0])
        require(target0 in lookup, "signed-pair relabelling left the carrier")
        j, bit = lookup[target0]
        require(
            relabel(g, pair[1]) == pairs[j][bit ^ 1],
            "signed-pair action is inconsistent",
        )
        sigma[i] = j
        delta[i] = bit
    require(sorted(sigma) == [0, 1, 2], "block action is not a permutation")
    return tuple(delta), tuple(sigma)


def compose_signed(a, b):
    # First b, then a, for e_(i,x) -> e_(sigma(i),x+delta_i).
    delta_a, sigma_a = a
    delta_b, sigma_b = b
    delta = tuple(delta_b[i] ^ delta_a[sigma_b[i]] for i in range(3))
    sigma = compose_perm(sigma_a, sigma_b)
    return delta, sigma


def central_twist(a, exponent):
    require(exponent in (0, 1), "central-twist exponent is not a bit")
    delta, sigma = a
    return tuple(bit ^ exponent for bit in delta), sigma


def generated_subgroup(generators):
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    frontier = [identity]
    while frontier:
        g = frontier.pop()
        for h in generators:
            for product in (compose_perm(g, h), compose_perm(h, g)):
                if product not in seen:
                    seen.add(product)
                    frontier.append(product)
    return seen


# The pair order is the common matching quotient frame.  The orientation
# order and within-pair gauges are deliberately frozen: with these exact
# labels rho_or(g)=z^parity(g) rho_edge(g) pointwise on all 24 elements.
edge_pairs = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)
orientation_pairs = (
    ((0, 2, 1, 3), (0, 3, 1, 2)),
    ((0, 3, 2, 1), (0, 1, 2, 3)),
    ((0, 1, 3, 2), (0, 2, 3, 1)),
)

cycles = {rotate_to_zero(c) for c in permutations(range(4))}
orientation_objects = {c for pair in orientation_pairs for c in pair}
require(len(cycles) == 6, "oriented four-cycle census changed")
require(orientation_objects == cycles, "orientation-pair labels changed")
require(
    all(inverse_cycle(pair[0]) == pair[1] for pair in orientation_pairs),
    "an orientation pair is not an inverse pair",
)

S3 = list(permutations(range(3)))
W = {
    (delta, sigma)
    for delta in (
        (a, b, c) for a in range(2) for b in range(2) for c in range(2)
    )
    for sigma in S3
}
H0 = {
    (delta, sigma)
    for delta, sigma in W
    if sum(delta) % 2 == 0
}
H1 = {
    (delta, sigma)
    for delta, sigma in W
    if sum(delta) % 2 == parity(sigma)
}
require(len(W) == 48, "signed-pair group census changed")
require(len(H0) == len(H1) == 24, "complement census changed")


def action_data(name, pairs, relabel):
    image = set()
    table = {}
    ordered = tuple(x for pair in pairs for x in pair)
    lookup = {x: j for j, x in enumerate(ordered)}
    require(len(lookup) == 6, f"{name}: ordered carrier has repeated labels")
    for g in S4:
        delta, sigma = signed_pair_action(pairs, relabel, g)
        p = sum(delta) % 2
        s = parity(sigma)
        six_perm = tuple(lookup[relabel(g, x)] for x in ordered)
        six_type = cycle_type(six_perm)
        sheet_type = cycle_type(g)
        table.setdefault(sheet_type, set()).add(six_type)
        require(parity(six_perm) == p, f"{name}: ambient sign != flip parity")
        require(s == parity(g), f"{name}: matching sign != quartic sign")
        image.add((delta, sigma))
    require(len(image) == 24, f"{name}: action is not faithful")
    require(
        all(len(types) == 1 for types in table.values()),
        f"{name}: a quartic conjugacy class split",
    )
    return image, {key: next(iter(values)) for key, values in table.items()}


edge_image, edge_table = action_data("edge", edge_pairs, relabel_edge)
orientation_image, orientation_table = action_data(
    "orientation", orientation_pairs, relabel_cycle
)
require(edge_image == H0, "edge lift is not exactly H0")
require(orientation_image == H1, "orientation lift is not exactly H1")

for g in S4:
    edge_rho = signed_pair_action(edge_pairs, relabel_edge, g)
    orientation_rho = signed_pair_action(orientation_pairs, relabel_cycle, g)
    require(
        orientation_rho == central_twist(edge_rho, parity(g)),
        f"pointwise central-sign twist failed at {g}",
    )
    for h in S4:
        gh = compose_perm(g, h)
        require(
            signed_pair_action(edge_pairs, relabel_edge, gh)
            == compose_signed(
                signed_pair_action(edge_pairs, relabel_edge, g),
                signed_pair_action(edge_pairs, relabel_edge, h),
            ),
            "edge map is not a homomorphism",
        )
        require(
            signed_pair_action(orientation_pairs, relabel_cycle, gh)
            == compose_signed(
                signed_pair_action(orientation_pairs, relabel_cycle, g),
                signed_pair_action(orientation_pairs, relabel_cycle, h),
            ),
            "orientation map is not a homomorphism",
        )

global_flip = ((1, 1, 1), (0, 1, 2))
require(global_flip not in edge_image, "edge image contains the central flip")
require(
    global_flip not in orientation_image,
    "orientation image contains the central flip",
)

class_order = ((1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,))
expected_edge = {
    (1, 1, 1, 1): (1, 1, 1, 1, 1, 1),
    (2, 1, 1): (2, 2, 1, 1),
    (2, 2): (2, 2, 1, 1),
    (3, 1): (3, 3),
    (4,): (4, 2),
}
expected_orientation = {
    (1, 1, 1, 1): (1, 1, 1, 1, 1, 1),
    (2, 1, 1): (2, 2, 2),
    (2, 2): (2, 2, 1, 1),
    (3, 1): (3, 3),
    (4,): (4, 1, 1),
}
require(edge_table == expected_edge, "edge cycle table changed")
require(orientation_table == expected_orientation, "orientation cycle table changed")

class_sizes = (1, 6, 3, 8, 6)
irreducibles = {
    "1": (1, 1, 1, 1, 1),
    "sgn": (1, -1, 1, 1, -1),
    "2": (2, 0, 2, -1, 0),
    "3": (3, 1, -1, 0, -1),
    "3sgn": (3, -1, -1, 0, 1),
}


def decompose_permutation(table, det_character):
    chi_u = tuple(table[key].count(1) for key in class_order)
    chi_l5 = tuple(det_character[i] * chi_u[i] for i in range(5))

    def decompose(chi):
        result = {}
        for name, irr in irreducibles.items():
            numerator = sum(
                class_sizes[i] * chi[i] * irr[i] for i in range(5)
            )
            require(
                numerator % 24 == 0,
                f"nonintegral character multiplicity for {name}",
            )
            result[name] = numerator // 24
        return result

    return decompose(chi_u), decompose(chi_l5)


edge_u, edge_l5 = decompose_permutation(edge_table, (1, 1, 1, 1, 1))
orientation_u, orientation_l5 = decompose_permutation(
    orientation_table, irreducibles["sgn"]
)
require(
    edge_u == {"1": 1, "sgn": 0, "2": 1, "3": 1, "3sgn": 0},
    "edge U decomposition changed",
)
require(edge_l5 == edge_u, "edge Lambda5 decomposition changed")
require(
    orientation_u == {"1": 1, "sgn": 0, "2": 1, "3": 0, "3sgn": 1},
    "orientation U decomposition changed",
)
require(
    orientation_l5 == {"1": 0, "sgn": 1, "2": 1, "3": 1, "3sgn": 0},
    "orientation Lambda5 decomposition changed",
)

# The odd generator is a sheet transposition already specified at the S4
# level.  In the frozen common frame it is exactly THM-2965's s0/s1 control.
odd_transposition = (0, 1, 3, 2)
edge_odd = signed_pair_action(edge_pairs, relabel_edge, odd_transposition)
orientation_odd = signed_pair_action(
    orientation_pairs, relabel_cycle, odd_transposition
)
require(
    edge_odd == ((0, 0, 0), (0, 2, 1)),
    "edge transposition control changed",
)
require(
    orientation_odd == ((1, 1, 1), (0, 2, 1)),
    "orientation transposition control changed",
)

ternary_sheet_generator = (2, 0, 1, 3)
expected_ternary = ((0, 1, 1), (1, 2, 0))
edge_ternary = signed_pair_action(
    edge_pairs, relabel_edge, ternary_sheet_generator
)
orientation_ternary = signed_pair_action(
    orientation_pairs, relabel_cycle, ternary_sheet_generator
)
require(cycle_type(ternary_sheet_generator) == (3, 1), "ternary control is not order three")
require(edge_ternary == expected_ternary, "edge ternary control changed")
require(
    orientation_ternary == expected_ternary,
    "orientation ternary control changed",
)
require(
    generated_subgroup((odd_transposition, ternary_sheet_generator)) == set(S4),
    "binary--ternary controls do not generate S4",
)

# Exact THM-2769 specialization of the THM-2864 discriminant gates.
t = sp.symbols("t")
p = -2
q = -8 * t
r = 1 - 4 * t
Delta = sp.expand(
    256 * r**3
    - 128 * p**2 * r**2
    + 144 * p * q**2 * r
    - 27 * q**4
    + 16 * p**4 * r
    - 4 * p**3 * q**2
)
J_or = sp.expand(
    1024 * r**3 + 768 * p**2 * r**2 - 288 * p * q**2 * r + 27 * q**4
)
expected_delta = -2**12 * t**2 * (27 * t**2 - 14 * t + 3)
expected_jor = 2**12 * (t - 1) * (27 * t**3 - 25 * t**2 + 8 * t - 1)
require(
    sp.expand(Delta - expected_delta) == 0,
    "THM-2769 quartic discriminant changed",
)
require(
    sp.expand(J_or - expected_jor) == 0,
    "THM-2769 orientation wall changed",
)
require(
    not sp.Poly(q * Delta * J_or, t).is_zero,
    "physical common separability open is empty",
)

edge_disc = sp.expand(2**6 * q**2 * Delta**2)
edge_square = sp.expand((2**3 * q * Delta) ** 2)
orientation_disc = sp.expand(2**66 * q**12 * J_or**4 * Delta**3)
orientation_square_factor = sp.expand(2**33 * q**6 * J_or**2 * Delta)
require(edge_disc == edge_square, "edge discriminant square identity changed")
require(
    orientation_disc == sp.expand(Delta * orientation_square_factor**2),
    "orientation discriminant square class changed",
)

edge_table_ordered = {key: edge_table[key] for key in class_order}
orientation_table_ordered = {key: orientation_table[key] for key in class_order}

print("THM-2968 quartic signed-pair complement bridge")
print(f"S4={len(S4)} W={len(W)} H0={len(H0)} H1={len(H1)}")
print(f"orientation_pairs={orientation_pairs}")
print("pointwise_rho_or=z^sgn*rho_edge checks=24")
print(f"edge_image={len(edge_image)} class=H0 condition=p=0")
print(f"orientation_image={len(orientation_image)} class=H1 condition=p=s")
print(f"edge_cycle_table={edge_table_ordered}")
print(f"orientation_cycle_table={orientation_table_ordered}")
print("edge_U=1+2+3 edge_Lambda5=1+2+3")
print("orientation_U=1+2+3sgn orientation_Lambda5=sgn+2+3")
print(f"odd_sheet_transposition={odd_transposition}")
print("edge_transposition_control=(000,(0,2,1))")
print("orientation_transposition_control=(111,(0,2,1))")
print(f"ternary_sheet_generator={ternary_sheet_generator}")
print("common_ternary_control=(011,(1,2,0))")
print("binary_ternary_generate_S4=24")
print(f"THM2769_Delta={sp.factor(Delta)}")
print(f"THM2769_J_or={sp.factor(J_or)}")
print("edge_discriminant_squareclass=1")
print("orientation_discriminant_squareclass=Delta")
print("THM2769_generic_S4_witness=imported")
print("physical_generic_S4_realizes=H0+H1")
print("quartic_edge_orientation_s4_complements_thm2968=PASS")
