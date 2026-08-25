#!/usr/bin/env python3
"""Clean-room boundary-DP and order-six audit for THM-4123."""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb, factorial, gcd


EXPECTED_SEMANTIC = "0f227c83e27b0f2beda7d432c46c2821318cf5cd1277c4c6bafd6e1de67963ec"
EXPECTED_EXTENDED = "b0ee0ce581dccf3af6ba937473c45717b6132f7a396a2dee36be8b62f0efd178"


def require(condition, label):
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def digest(value):
    return sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def fraction_pair(value):
    value = Fraction(value)
    return (value.numerator, value.denominator)


def gcd_all(values):
    return reduce(gcd, (abs(value) for value in values), 0)


def integer_ceiling(value):
    value = Fraction(value)
    return -((-value.numerator) // value.denominator)


def odd_ceiling(value):
    ceiling = integer_ceiling(value)
    return ceiling if ceiling % 2 else ceiling + 1


def coset_ceiling(value, anchor, modulus):
    value = Fraction(value)
    if modulus == 0:
        require(value == anchor, "zero layer lattice means a constant layer")
        return anchor
    return anchor + modulus * integer_ceiling((value - anchor) / modulus)


def decode(code, order):
    adjacency = [0] * order
    bit = 0
    for left in range(order):
        for right in range(left + 1, order):
            if (code >> bit) & 1:
                adjacency[left] |= 1 << right
            else:
                adjacency[right] |= 1 << left
            bit += 1
    return tuple(adjacency)


def is_strong(adjacency):
    order = len(adjacency)
    full = (1 << order) - 1
    for root in range(order):
        reached = 1 << root
        frontier = reached
        while frontier:
            following = 0
            todo = frontier
            while todo:
                bit = todo & -todo
                todo ^= bit
                following |= adjacency[bit.bit_length() - 1]
            following &= ~reached
            reached |= following
            frontier = following
        if reached != full:
            return False
    return True


def one_defect_count(adjacency):
    total = 0
    for word in permutations(range(len(adjacency))):
        defects = sum(
            not ((adjacency[word[index]] >> word[index + 1]) & 1)
            for index in range(len(word) - 1)
        )
        total += defects == 1
    return total


def boundary(adjacency):
    order = len(adjacency)
    size = 1 << order
    full = size - 1
    starts = [[0] * order for _ in range(size)]
    ends = [[0] * order for _ in range(size)]
    for vertex in range(order):
        starts[1 << vertex][vertex] = 1
        ends[1 << vertex][vertex] = 1
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        for vertex in range(order):
            bit = 1 << vertex
            if not mask & bit:
                continue
            rest = mask ^ bit
            todo = rest
            while todo:
                other_bit = todo & -todo
                todo ^= other_bit
                other = other_bit.bit_length() - 1
                if adjacency[vertex] & other_bit:
                    starts[mask][vertex] += starts[rest][other]
                if adjacency[other] & bit:
                    ends[mask][vertex] += ends[rest][other]

    transitions = [[0] * order for _ in range(order)]
    for left_mask in range(1, full):
        right_mask = full ^ left_mask
        left_todo = left_mask
        while left_todo:
            left_bit = left_todo & -left_todo
            left_todo ^= left_bit
            left = left_bit.bit_length() - 1
            if not ends[left_mask][left]:
                continue
            right_todo = right_mask
            while right_todo:
                right_bit = right_todo & -right_todo
                right_todo ^= right_bit
                right = right_bit.bit_length() - 1
                transitions[left][right] += (
                    ends[left_mask][left] * starts[right_mask][right]
                )
    return (
        tuple(starts[full]), tuple(ends[full]),
        tuple(tuple(row) for row in transitions),
    )


def boundary_responses(data):
    starts, ends, transitions = data
    order = len(starts)
    values = []
    for mask in range(1 << order):
        value = sum(
            starts[vertex] if mask & (1 << vertex) else ends[vertex]
            for vertex in range(order)
        )
        for left in range(order):
            if mask & (1 << left):
                continue
            for right in range(order):
                if mask & (1 << right):
                    value += transitions[left][right]
        values.append(value)
    return tuple(values)


def literal_add_ear(adjacency, cut):
    order = len(adjacency)
    child = list(adjacency) + [0]
    for vertex in range(order):
        if cut & (1 << vertex):
            child[order] |= 1 << vertex
        else:
            child[vertex] |= 1 << order
    return tuple(child)


def literal_hamilton_count(adjacency):
    return sum(
        all(adjacency[word[index]] & (1 << word[index + 1])
            for index in range(len(word) - 1))
        for word in permutations(range(len(adjacency)))
    )


def literal_responses(adjacency):
    return tuple(
        literal_hamilton_count(literal_add_ear(adjacency, cut))
        for cut in range(1 << len(adjacency))
    )


def coefficients(values, order):
    base = values[0]
    linear = tuple(values[1 << vertex] - base for vertex in range(order))
    curvature = [[0] * order for _ in range(order)]
    for left, right in combinations(range(order), 2):
        value = (
            values[1 << left] + values[1 << right] - base
            - values[(1 << left) | (1 << right)]
        )
        curvature[left][right] = curvature[right][left] = value
    flat = linear + tuple(
        curvature[left][right]
        for left, right in combinations(range(order), 2)
    )
    return linear, tuple(tuple(row) for row in curvature), gcd_all(flat)


def layer_data(values, linear, curvature, size):
    order = len(linear)
    masks = tuple(mask for mask in range(1 << order) if mask.bit_count() == size)
    anchor = values[masks[0]]
    difference_gcd = gcd_all(values[mask] - anchor for mask in masks)
    exchange_increments = []
    if 1 <= size <= order - 1:
        for left, right in combinations(range(order), 2):
            others = tuple(
                vertex for vertex in range(order) if vertex not in (left, right)
            )
            for context in combinations(others, size - 1):
                context_mask = sum(1 << vertex for vertex in context)
                actual = (
                    values[context_mask | (1 << right)]
                    - values[context_mask | (1 << left)]
                )
                predicted = linear[right] - linear[left] - sum(
                    curvature[right][vertex] - curvature[left][vertex]
                    for vertex in context
                )
                require(actual == predicted, "clean-room exchange identity")
                exchange_increments.append(actual)
    require(gcd_all(exchange_increments) == difference_gcd,
            "clean-room Johnson lattice equality")
    image = tuple(sorted({values[mask] for mask in masks}))
    mean = Fraction(sum(values[mask] for mask in masks), len(masks))
    floor = coset_ceiling(mean, anchor, difference_gcd)
    require(max(image) >= floor, "clean-room layer coset rounding")
    return {
        "size": size,
        "mean": fraction_pair(mean),
        "maximum": max(image),
        "anchor": anchor,
        "lattice": difference_gcd,
        "lattice_floor": floor,
        "image": image,
    }


def analyze(adjacency):
    order = len(adjacency)
    values = boundary_responses(boundary(adjacency))
    base = values[0]
    defect_one = one_defect_count(adjacency)
    linear, curvature, global_lattice = coefficients(values, order)
    require(values[-1] == base, "clean-room constant ears")
    require(all(value % 2 for value in values), "clean-room oddness")
    weight = Fraction(sum(
        curvature[left][right]
        for left, right in combinations(range(order), 2)
    ), 2)
    require(2 * weight == (order - 1) * base + defect_one,
            "clean-room cut-weight identity")

    layers = []
    for size in range(order + 1):
        record = layer_data(values, linear, curvature, size)
        multiplicity = (
            comb(order - 2, size - 1) if 1 <= size <= order - 1 else 0
        )
        predicted_sum = (
            comb(order, size) * base
            + multiplicity * ((order - 1) * base + defect_one)
        )
        actual_sum = Fraction(*record["mean"]) * comb(order, size)
        require(actual_sum == predicted_sum, "clean-room slice sum")
        beta = (
            Fraction(2 * size * (order - size), order * (order - 1))
            if 1 <= size <= order - 1 else Fraction(0)
        )
        require(record["mean"] == fraction_pair(base + beta * weight),
                "clean-room slice mean")
        if record["lattice"]:
            require(record["lattice"] % global_lattice == 0,
                    "clean-room global-to-layer lattice divisibility")
        layers.append(record)

    balanced_size = order // 2
    balanced = layers[balanced_size]
    q = order * order // 4
    balanced_bound = (
        (1 + Fraction(q, order)) * base
        + Fraction(q, order * (order - 1)) * defect_one
    )
    require(balanced["mean"] == fraction_pair(balanced_bound),
            "clean-room balanced mean")
    for size in {order // 2, (order + 1) // 2}:
        require(layers[size]["mean"] == fraction_pair(balanced_bound),
                "clean-room central-layer equality")
        require(layers[size]["maximum"] >= odd_ceiling(balanced_bound),
                "clean-room balanced odd-ceiling")

    reduced = (
        Fraction((order + 1) * (order + 2), 4 * order) * base
        + Fraction(order * (order - 1) + 2, 4 * order * (order - 1))
        * defect_one
    )
    index = order - 2 if order % 2 == 0 else order - 3
    difference = Fraction(index, 4 * order) * (
        base + Fraction(defect_one, order - 1)
    )
    require(balanced_bound - reduced == difference,
            "clean-room THM-4115 comparison")

    edge_weights = tuple(
        Fraction(curvature[left][right], 2)
        for left, right in combinations(range(order), 2)
    )
    fields = tuple(
        Fraction(linear[vertex])
        - sum(Fraction(curvature[vertex][other], 2)
              for other in range(order) if other != vertex)
        for vertex in range(order)
    )
    energy = sum(field * field for field in fields) + sum(
        edge * edge for edge in edge_weights
    )
    variance_floor = base + weight / 2 + energy / (2 * weight)

    return {
        "order": order,
        "H": base,
        "F1": defect_one,
        "W": fraction_pair(weight),
        "global_lattice": global_lattice,
        "balanced_bound": fraction_pair(balanced_bound),
        "balanced": balanced,
        "layer_lattices": tuple(layer["lattice"] for layer in layers),
        "layer_maxima": tuple(layer["maximum"] for layer in layers),
        "variance_floor": fraction_pair(variance_floor),
        "full_maximum": max(values),
    }


def compact_named(record):
    balanced = record["balanced"]
    return {
        "order": record["order"],
        "H_F1_W": (record["H"], record["F1"], record["W"]),
        "global_lattice": record["global_lattice"],
        "balanced": (
            balanced["size"], balanced["mean"], balanced["maximum"],
            balanced["anchor"], balanced["lattice"],
            balanced["lattice_floor"], balanced["image"],
        ),
        "layer_lattices": record["layer_lattices"],
        "layer_maxima": record["layer_maxima"],
        "variance_floor": record["variance_floor"],
        "full_maximum": record["full_maximum"],
    }


def main():
    order_rows = []
    totals = [0, 0, 0, 0]
    for order in range(2, 6):
        parent_count = 1 << (order * (order - 1) // 2)
        strong_count = 0
        histogram = {}
        strong_histogram = {}
        equality_count = strong_equality_count = 0
        for code in range(parent_count):
            adjacency = decode(code, order)
            record = analyze(adjacency)
            lattice = record["balanced"]["lattice"]
            histogram[lattice] = histogram.get(lattice, 0) + 1
            equality = (
                record["balanced"]["maximum"]
                == odd_ceiling(Fraction(*record["balanced_bound"]))
            )
            equality_count += equality
            if is_strong(adjacency):
                strong_count += 1
                strong_histogram[lattice] = strong_histogram.get(lattice, 0) + 1
                strong_equality_count += equality
        row = {
            "order": order,
            "parents": parent_count,
            "strong": strong_count,
            "ears": parent_count * (1 << order),
            "layers": parent_count * (order + 1),
            "ordering_checks": parent_count * factorial(order),
            "balanced_lattice_histogram": tuple(sorted(histogram.items())),
            "strong_balanced_lattice_histogram": tuple(sorted(strong_histogram.items())),
            "oddceil_equality_all_strong": (equality_count, strong_equality_count),
        }
        order_rows.append(row)
        totals[0] += row["parents"]
        totals[1] += row["ears"]
        totals[2] += row["layers"]
        totals[3] += row["ordering_checks"]

    require(tuple(totals) == (1098, 33864, 6502, 124468),
            "clean-room census totals")
    require(tuple((row["parents"], row["strong"]) for row in order_rows)
            == ((2, 0), (8, 2), (64, 24), (1024, 544)),
            "clean-room parent census")
    require(tuple(row["balanced_lattice_histogram"] for row in order_rows) == (
        ((2, 2),), ((0, 2), (2, 6)), ((2, 48), (8, 16)), ((2, 1024),)
    ), "clean-room balanced lattice census")
    require(tuple(row["strong_balanced_lattice_histogram"] for row in order_rows) == (
        (), ((0, 2),), ((2, 24),), ((2, 544),)
    ), "clean-room strong balanced lattice census")
    require(tuple(row["oddceil_equality_all_strong"] for row in order_rows) == (
        (2, 0), (2, 2), (0, 0), (24, 24)
    ), "clean-room equality census")

    named_codes = (
        ("c3", 2, 3), ("code8", 8, 5), ("code759", 759, 5),
        ("regular_c5", 76, 5), ("code40", 40, 5),
        ("code20", 20, 6), ("code30571", 30571, 6),
    )
    named = {
        name: compact_named(analyze(decode(code, order)))
        for name, code, order in named_codes
    }
    for name, code, order in (
        ("c3", 2, 3), ("code20", 20, 6), ("code30571", 30571, 6)
    ):
        boundary_values = boundary_responses(boundary(decode(code, order)))
        require(boundary_values == literal_responses(decode(code, order)),
                f"literal child-ordering replay for {name}")

    require(named["c3"]["balanced"][1:3] == ((5, 1), 5),
            "clean-room C3 equality")
    require(named["code8"]["H_F1_W"] == named["code759"]["H_F1_W"]
            == (9, 30, (33, 1)), "clean-room equal-mean hostile")
    require(named["code8"]["balanced"][-1] != named["code759"]["balanced"][-1],
            "clean-room different balanced images")
    require(named["regular_c5"]["balanced"][1:3] == ((42, 1), 43),
            "clean-room sharp odd-ceiling")
    require(Fraction(*named["code8"]["variance_floor"])
            > Fraction(*named["code8"]["balanced"][1])
            and Fraction(*named["code40"]["variance_floor"])
            < Fraction(*named["code40"]["balanced"][1]),
            "clean-room variance incomparability")
    require(named["code20"]["layer_maxima"] == (19, 91, 133, 131, 123, 91, 19),
            "clean-room code20 layer maxima")
    require(named["code30571"]["layer_lattices"] == (0, 18, 6, 6, 6, 18, 0),
            "clean-room code30571 layer lattices")

    core_ledger = {
        "formula": "sum_|S|=m F=C(n,m)H+C(n-2,m-1)((n-1)H+F1)",
        "balanced_bound": "(1+floor(n^2/4)/n)H+floor(n^2/4)F1/(n(n-1))",
        "layer_lattice": "gcd_exchange_increments_equals_gcd_layer_differences",
        "order_rows": order_rows,
        "totals": tuple(totals),
        "named": named,
        "scope": "one large balanced child and exact layer coset; no interval or maximizing-layer theorem",
    }
    semantic = digest(core_ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "shared semantic digest")

    strong6 = no_balanced_global = 0
    first_failure = None
    balanced_lattice_histogram = {}
    for code in range(1 << 15):
        adjacency = decode(code, 6)
        if not is_strong(adjacency):
            continue
        strong6 += 1
        values = boundary_responses(boundary(adjacency))
        balanced_masks = tuple(mask for mask in range(64) if mask.bit_count() == 3)
        balanced_maximum = max(values[mask] for mask in balanced_masks)
        global_maximum = max(values)
        linear, curvature, _ = coefficients(values, 6)
        balanced_lattice = layer_data(values, linear, curvature, 3)["lattice"]
        balanced_lattice_histogram[balanced_lattice] = (
            balanced_lattice_histogram.get(balanced_lattice, 0) + 1
        )
        if balanced_maximum < global_maximum:
            no_balanced_global += 1
            if first_failure is None:
                maximum_masks = tuple(
                    mask for mask, value in enumerate(values) if value == global_maximum
                )
                layer_maxima = tuple(
                    max(values[mask] for mask in range(64) if mask.bit_count() == size)
                    for size in range(7)
                )
                first_failure = (
                    code, global_maximum, balanced_maximum, maximum_masks, layer_maxima
                )

    order6 = {
        "strong": strong6,
        "without_balanced_global_maximizer": no_balanced_global,
        "first_failure": first_failure,
        "balanced_lattice_histogram": tuple(sorted(balanced_lattice_histogram.items())),
    }
    require(strong6 == 22320, "order-six strong labelled census")
    require(no_balanced_global == 1920, "order-six balanced-global hostile count")
    require(first_failure == (20, 133, 131, (48,), (19, 91, 133, 131, 123, 91, 19)),
            "first balanced-global failure")
    extended = digest({"core_semantic": semantic, "order6": order6})
    if EXPECTED_EXTENDED is not None:
        require(extended == EXPECTED_EXTENDED, "extended semantic digest")

    print("status=PASS")
    print("implementation=clean-room Start/End/Q boundary plus literal ordering controls")
    print("slice_sum=C(n,m)H+C(n-2,m-1)((n-1)H+F1)")
    print("balanced_bound=(1+floor(n^2/4)/n)H+floor(n^2/4)F1/(n(n-1))")
    print(f"order_rows={order_rows}")
    print(f"totals_parents_ears_layers_orderings={tuple(totals)}")
    print(f"literal_named_replays={(('c3', 2), ('code20', 20), ('code30571', 30571))}")
    print(f"order6_census={order6}")
    print("scope=large balanced child and layer coset only; interval and global-max location remain OPEN")
    print(f"semantic_sha256={semantic}")
    print(f"extended_semantic_sha256={extended}")


if __name__ == "__main__":
    main()
