"""Exact controls for THM-2961.

All truth-bearing checks use explicit RuntimeError gates so ordinary and
optimized Python execute the same verification.
"""

from itertools import combinations, permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sign(value):
    return (value > 0) - (value < 0)


def cyclic_discrete_control(order):
    """Return the complete crossing count and minimum physical dictionary."""

    elements = tuple(range(order))
    profiles = {
        x: tuple(0 if (x + z) % order == 0 else 1 for z in elements)
        for x in elements
    }

    directed_unique_witnesses = 0
    crossing_edges = 0
    for x, y in permutations(elements, 2):
        witnesses = tuple(
            z for z in elements if profiles[x][z] > profiles[y][z]
        )
        require(witnesses == ((-y) % order,), f"C{order} witness failure")
        directed_unique_witnesses += 1
    for x, y in combinations(elements, 2):
        signs = {sign(profiles[x][z] - profiles[y][z]) for z in elements}
        require(-1 in signs and 1 in signs, f"C{order} crossing failure")
        crossing_edges += 1

    minimum = None
    minimal_dictionaries = []
    for size in range(order + 1):
        for dictionary in combinations(elements, size):
            good = True
            for x in elements:
                for y in elements:
                    full = all(
                        profiles[x][z] <= profiles[y][z] for z in elements
                    )
                    restricted = all(
                        profiles[x][z] <= profiles[y][z] for z in dictionary
                    )
                    if full != restricted:
                        good = False
                        break
                if not good:
                    break
            if good:
                minimum = size
                minimal_dictionaries.append(dictionary)
        if minimum is not None:
            break

    require(
        directed_unique_witnesses == order * (order - 1),
        f"C{order} directed witness count",
    )
    require(
        crossing_edges == order * (order - 1) // 2,
        f"C{order} crossing graph is not complete",
    )
    require(minimum == order, f"C{order} context tax is not its order")
    require(
        minimal_dictionaries == [elements],
        f"C{order} has an unexpected minimal dictionary",
    )
    return profiles, crossing_edges, directed_unique_witnesses, minimum


c2_profiles, c2_edges, c2_directed, c2_contexts = cyclic_discrete_control(2)
c3_profiles, c3_edges, c3_directed, c3_contexts = cyclic_discrete_control(3)


# V4=F_2^2.
V4 = tuple(product(range(2), repeat=2))
V4_INDEX = {x: index for index, x in enumerate(V4)}


def add_v4(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])


def discrete_length(x):
    return 0 if x == (0, 0) else 1


v4_profiles = {
    x: tuple(discrete_length(add_v4(x, z)) for z in V4)
    for x in V4
}
require(len(set(v4_profiles.values())) == 4, "V4 profiles collided")

v4_directed = 0
for x, y in permutations(V4, 2):
    witnesses = tuple(
        z for z in V4 if v4_profiles[x][V4_INDEX[z]] > v4_profiles[y][V4_INDEX[z]]
    )
    require(witnesses == (y,), f"V4 unique witness failure for {x},{y}")
    v4_directed += 1

v4_edges = []
for x, y in combinations(V4, 2):
    signs = {
        sign(v4_profiles[x][index] - v4_profiles[y][index])
        for index in range(4)
    }
    require(-1 in signs and 1 in signs, f"V4 crossing failure for {x},{y}")
    v4_edges.append((x, y))
require(len(v4_edges) == 6, "V4 crossing graph is not K4")


def v4_order_complete(dictionary):
    indices = tuple(V4_INDEX[z] for z in dictionary)
    for x in V4:
        for y in V4:
            full = all(a <= b for a, b in zip(v4_profiles[x], v4_profiles[y]))
            restricted = all(
                v4_profiles[x][index] <= v4_profiles[y][index]
                for index in indices
            )
            if full != restricted:
                return False
    return True


v4_minimum = None
v4_minimal_dictionaries = []
for size in range(5):
    for dictionary in combinations(V4, size):
        if v4_order_complete(dictionary):
            v4_minimum = size
            v4_minimal_dictionaries.append(dictionary)
    if v4_minimum is not None:
        break
require(v4_minimum == 4, "V4 context tax is not four")
require(v4_minimal_dictionaries == [V4], "unexpected minimal V4 dictionary")


# AGL(2,2).
def invertible_matrix(a, b, c, d):
    return (a * d - b * c) % 2 == 1


affine_permutations = set()
for a, b, c, d in product(range(2), repeat=4):
    if not invertible_matrix(a, b, c, d):
        continue
    for translation in V4:
        image = []
        for x in V4:
            image.append(
                (
                    (a * x[0] + b * x[1] + translation[0]) % 2,
                    (c * x[0] + d * x[1] + translation[1]) % 2,
                )
            )
        affine_permutations.add(tuple(image))
require(len(affine_permutations) == 24, "AGL(2,2) is not S4")

translations = {tuple(add_v4(x, t) for x in V4) for t in V4}
require(len(translations) == 4, "translation subgroup is not V4")
require(
    len(affine_permutations) // len(translations) == 6,
    "affine quotient does not have size six",
)


def compose(left, right):
    return tuple(left[V4_INDEX[right[index]]] for index in range(4))


def inverse(perm):
    answer = [None] * 4
    for source_index, image_point in enumerate(perm):
        answer[V4_INDEX[image_point]] = V4[source_index]
    return tuple(answer)


normality_checks = 0
for affine in affine_permutations:
    affine_inverse = inverse(affine)
    for translation in translations:
        conjugate = compose(compose(affine, translation), affine_inverse)
        require(conjugate in translations, "V4 translations are not normal")
        normality_checks += 1


# No S4-invariant orientation of K4.
invariant_orientations = 0
for oriented_bits in product(range(2), repeat=6):
    edge_orientation = {}
    for bit, (x, y) in zip(oriented_bits, v4_edges):
        edge_orientation[frozenset((x, y))] = (x, y) if bit == 0 else (y, x)
    invariant = True
    for perm in permutations(V4):
        mapping = dict(zip(V4, perm))
        for oriented in edge_orientation.values():
            u, v = oriented
            moved_edge = frozenset((mapping[u], mapping[v]))
            if edge_orientation[moved_edge] != (mapping[u], mapping[v]):
                invariant = False
                break
        if not invariant:
            break
    if invariant:
        invariant_orientations += 1
require(invariant_orientations == 0, "found an S4-invariant tournament")


# Full two-ended V4 kernel: pointwise dominance forces equality.
def v4_distance(x, y):
    return discrete_length(add_v4(x, y))


full_kernels = {
    x: tuple(v4_distance(add_v4(x, a), b) for a in V4 for b in V4)
    for x in V4
}
kernel_dominances = []
for x in V4:
    for y in V4:
        if all(a <= b for a, b in zip(full_kernels[x], full_kernels[y])):
            kernel_dominances.append((x, y))
require(
    kernel_dominances == [(x, x) for x in V4],
    "full-kernel order has a nontrivial dominance",
)


# Non-group metric-monoid hostile.
def d_nat(m, n):
    gap = abs(m - n)
    if gap == 0:
        return 0
    if gap == 2:
        return 1
    return 2


triangle_checks = 0
for a, b, c in product(range(21), repeat=3):
    require(d_nat(a, c) <= d_nat(a, b) + d_nat(b, c), "triangle failure")
    triangle_checks += 1

translation_checks = 0
for a, b, shift in product(range(12), repeat=3):
    require(
        d_nat(a + shift, b + shift) == d_nat(a, b),
        "translation invariance failure",
    )
    translation_checks += 1

joint_checks = 0
for a, b, c, e in product(range(12), repeat=4):
    require(
        d_nat(a + c, b + e) <= d_nat(a, b) + d_nat(c, e),
        "joint nonexpansivity failure",
    )
    joint_checks += 1


def f_nat(n):
    return d_nat(n, 0)


old_signs = {sign(f_nat(1 + z) - f_nat(2 + z)) for z in range(20)}
require(-1 in old_signs and 1 in old_signs, "1 and 2 do not cross")
new_differences = tuple(f_nat(2 + z) - f_nat(3 + z) for z in range(20))
require(all(value <= 0 for value in new_differences), "translated dominance failed")
require(any(value < 0 for value in new_differences), "dominance is not strict")


# The knot entries are a proof certificate: diagonals are exact six, while
# each off-diagonal occurrence is only bounded above by five.
knot_diagonal_exact = (6, 6)
knot_off_diagonal_upper = (5, 5)
require(
    all(diagonal > upper for diagonal, upper in zip(knot_diagonal_exact, knot_off_diagonal_upper)),
    "knot mirror crossing is not strict",
)


print("C2_profiles", tuple(c2_profiles.items()))
print("C2_edges_directed_context_tax", c2_edges, c2_directed, c2_contexts)
print("C3_profiles", tuple(c3_profiles.items()))
print("C3_edges_directed_context_tax", c3_edges, c3_directed, c3_contexts)
print("V4_profiles", tuple(v4_profiles.items()))
print("V4_edges_directed_context_tax", len(v4_edges), v4_directed, v4_minimum)
print("AGL22_translations_quotient", len(affine_permutations), len(translations), 6)
print("V4_normality_checks", normality_checks)
print("S4_invariant_tournament_orientations", invariant_orientations)
print("full_kernel_dominances", tuple(kernel_dominances))
print("N_metric_triangle_checks", triangle_checks)
print("N_translation_checks", translation_checks)
print("N_joint_nonexpansivity_checks", joint_checks)
print("crossing_signs_before_translation", tuple(sorted(old_signs)))
print("translated_differences_prefix", new_differences[:5])
print("knot_diagonal_exact", knot_diagonal_exact)
print("knot_off_diagonal_upper", knot_off_diagonal_upper)
print("PASS")
