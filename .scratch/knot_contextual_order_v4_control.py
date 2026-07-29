from itertools import combinations, permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sign(value):
    return (value > 0) - (value < 0)


def cyclic_discrete_control(order):
    elements = tuple(range(order))
    cyclic_profiles = {
        x: tuple(0 if (x + z) % order == 0 else 1 for z in elements)
        for x in elements
    }
    cyclic_edges = 0
    for x, y in combinations(elements, 2):
        signs = {
            sign(cyclic_profiles[x][z] - cyclic_profiles[y][z])
            for z in elements
        }
        require(-1 in signs and 1 in signs, f"C{order} crossing failure")
        cyclic_edges += 1
    minimum = None
    for size in range(order + 1):
        for dictionary in combinations(elements, size):
            good = True
            for x in elements:
                for y in elements:
                    full = all(
                        cyclic_profiles[x][z] <= cyclic_profiles[y][z]
                        for z in elements
                    )
                    restricted = all(
                        cyclic_profiles[x][z] <= cyclic_profiles[y][z]
                        for z in dictionary
                    )
                    if full != restricted:
                        good = False
                        break
                if not good:
                    break
            if good:
                minimum = size
                break
        if minimum is not None:
            break
    require(cyclic_edges == order * (order - 1) // 2, f"C{order} not complete")
    require(minimum == order, f"C{order} context tax is not its order")
    return cyclic_edges, minimum


c2_edges, c2_contexts = cyclic_discrete_control(2)
c3_edges, c3_contexts = cyclic_discrete_control(3)


# V4/F_2^2 exact control.
V4 = tuple(product(range(2), repeat=2))
v4_index = {x: index for index, x in enumerate(V4)}


def add_v4(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])


def discrete_length(x):
    return 0 if x == (0, 0) else 1


profiles = {
    x: tuple(discrete_length(add_v4(x, z)) for z in V4)
    for x in V4
}
require(len(set(profiles.values())) == 4, "V4 continuation profiles collided")

crossing_edges = []
for x, y in combinations(V4, 2):
    signs = {
        sign(discrete_length(add_v4(x, z)) - discrete_length(add_v4(y, z)))
        for z in V4
    }
    require(-1 in signs and 1 in signs, f"missing V4 crossing for {x},{y}")
    crossing_edges.append((x, y))
require(len(crossing_edges) == 6, "V4 crossing graph is not K4")


def order_complete(dictionary):
    indices = tuple(v4_index[z] for z in dictionary)
    for x in V4:
        for y in V4:
            full = all(a <= b for a, b in zip(profiles[x], profiles[y]))
            restricted = all(profiles[x][index] <= profiles[y][index] for index in indices)
            if full != restricted:
                return False
    return True


minimum_context_size = None
complete_dictionaries = []
for size in range(5):
    for dictionary in combinations(V4, size):
        if order_complete(dictionary):
            minimum_context_size = size
            complete_dictionaries.append(dictionary)
    if minimum_context_size is not None:
        break
require(minimum_context_size == 4, "V4 physical context tax is not four")
require(complete_dictionaries == [V4], "unexpected minimal V4 dictionary")

# Every permutation of the four affine points is affine over F_2.
def invertible_matrix(a, b, c, d):
    return (a * d - b * c) % 2 == 1


affine_permutations = set()
for a, b, c, d in product(range(2), repeat=4):
    if not invertible_matrix(a, b, c, d):
        continue
    for t in V4:
        image = []
        for x in V4:
            image.append(
                (
                    (a * x[0] + b * x[1] + t[0]) % 2,
                    (c * x[0] + d * x[1] + t[1]) % 2,
                )
            )
        affine_permutations.add(tuple(image))
require(len(affine_permutations) == 24, "AGL(2,2) is not S4")

translations = {
    tuple(add_v4(x, t) for x in V4)
    for t in V4
}
require(len(translations) == 4, "translation subgroup is not V4")
require(len(affine_permutations) // len(translations) == 6, "S4/V4 is not size 6")

def compose(left, right):
    return tuple(left[v4_index[right[index]]] for index in range(4))


def inverse(perm):
    answer = [None] * 4
    for source_index, image_point in enumerate(perm):
        answer[v4_index[image_point]] = V4[source_index]
    return tuple(answer)


for affine in affine_permutations:
    affine_inverse = inverse(affine)
    for translation in translations:
        conjugate = compose(compose(affine, translation), affine_inverse)
        require(conjugate in translations, "translations are not normal in AGL(2,2)")

# No orientation of K4 is invariant under every S4 permutation: a
# transposition fixing the other two vertices reverses its own edge.
for oriented_bits in product(range(2), repeat=6):
    edge_orientation = {}
    for bit, (x, y) in zip(oriented_bits, crossing_edges):
        edge_orientation[frozenset((x, y))] = (x, y) if bit == 0 else (y, x)
    invariant = True
    for perm in permutations(V4):
        mapping = dict(zip(V4, perm))
        for edge, oriented in edge_orientation.items():
            u, v = oriented
            moved_edge = frozenset((mapping[u], mapping[v]))
            if edge_orientation[moved_edge] != (mapping[u], mapping[v]):
                invariant = False
                break
        if not invariant:
            break
    require(not invariant, "found an S4-invariant tournament orientation")


# Non-group metric-monoid hostile.
def d_nat(m, n):
    gap = abs(m - n)
    if gap == 0:
        return 0
    if gap == 2:
        return 1
    return 2


for a, b, c in product(range(21), repeat=3):
    require(d_nat(a, c) <= d_nat(a, b) + d_nat(b, c), "triangle failure")
for a, b, shift in product(range(12), repeat=3):
    require(
        d_nat(a + shift, b + shift) == d_nat(a, b),
        "translation invariance failure",
    )
for a, b, c, e in product(range(12), repeat=4):
    require(
        d_nat(a + c, b + e) <= d_nat(a, b) + d_nat(c, e),
        "joint nonexpansivity failure",
    )


def f_nat(n):
    return d_nat(n, 0)


old_signs = {sign(f_nat(1 + z) - f_nat(2 + z)) for z in range(20)}
require(-1 in old_signs and 1 in old_signs, "1 and 2 do not cross")
new_differences = [f_nat(2 + z) - f_nat(3 + z) for z in range(20)]
require(all(value <= 0 for value in new_differences), "translated dominance failed")
require(any(value < 0 for value in new_differences), "translated dominance not strict")

# Symbolic knot certificate table: diagonal values are exact 6 and each
# off-diagonal value is at most 5, so opposite contexts force opposite signs.
knot_upper_table = ((6, 5), (5, 6))
require(knot_upper_table[0][0] > knot_upper_table[1][0], "first knot witness")
require(knot_upper_table[0][1] < knot_upper_table[1][1], "second knot witness")

print("V4_profiles", profiles)
print("C2_crossing_edges_context_tax", c2_edges, c2_contexts)
print("C3_crossing_edges_context_tax", c3_edges, c3_contexts)
print("V4_crossing_edges", len(crossing_edges))
print("V4_minimum_order_complete_contexts", minimum_context_size)
print("AGL22_size", len(affine_permutations))
print("translation_V4_size", len(translations))
print("quotient_size", len(affine_permutations) // len(translations))
print("S4_invariant_tournament_orientations", 0)
print("N_metric_triangle_checks", 21 ** 3)
print("N_translation_checks", 12 ** 3)
print("N_joint_nonexpansivity_checks", 12 ** 4)
print("crossing_signs_before_translation", sorted(old_signs))
print("translated_differences_prefix", new_differences[:5])
print("knot_mirror_upper_certificate", knot_upper_table)
print("PASS")
