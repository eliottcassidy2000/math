#!/usr/bin/env python3
"""Subset-DP referee for the balanced slice and layer lattice in THM-4123."""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb, factorial, gcd


EXPECTED_SEMANTIC = "0f227c83e27b0f2beda7d432c46c2821318cf5cd1277c4c6bafd6e1de67963ec"


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
        old = -1
        while reached != old:
            old = reached
            for vertex in range(order):
                if (reached >> vertex) & 1:
                    reached |= adjacency[vertex]
        if reached != full:
            return False
    return True


def hamilton_count(adjacency):
    order = len(adjacency)
    full = (1 << order) - 1
    paths = [[0] * order for _ in range(1 << order)]
    for vertex in range(order):
        paths[1 << vertex][vertex] = 1
    for mask in range(1, full + 1):
        for last in range(order):
            count = paths[mask][last]
            if not count:
                continue
            available = adjacency[last] & (full ^ mask)
            while available:
                bit = available & -available
                available ^= bit
                paths[mask | bit][bit.bit_length() - 1] += count
    return sum(paths[full])


def one_defect_count(adjacency):
    total = 0
    for word in permutations(range(len(adjacency))):
        defects = sum(
            not ((adjacency[word[index]] >> word[index + 1]) & 1)
            for index in range(len(word) - 1)
        )
        total += defects == 1
    return total


def add_ear(adjacency, cut):
    order = len(adjacency)
    child = list(adjacency) + [0]
    for vertex in range(order):
        if (cut >> vertex) & 1:
            child[order] |= 1 << vertex
        else:
            child[vertex] |= 1 << order
    return tuple(child)


def responses(adjacency):
    return tuple(
        hamilton_count(add_ear(adjacency, cut))
        for cut in range(1 << len(adjacency))
    )


def coefficients(values, order):
    base = values[0]
    linear = tuple(values[1 << vertex] - base for vertex in range(order))
    curvature = [[0] * order for _ in range(order)]
    for left in range(order):
        for right in range(left + 1, order):
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
                require(actual == predicted, "exchange coefficient identity")
                exchange_increments.append(actual)
    exchange_gcd = gcd_all(exchange_increments)
    require(exchange_gcd == difference_gcd, "Johnson exchange lattice equality")
    image = tuple(sorted({values[mask] for mask in masks}))
    mean = Fraction(sum(values[mask] for mask in masks), len(masks))
    lattice_floor = coset_ceiling(mean, anchor, difference_gcd)
    require(max(image) >= lattice_floor, "layer-coset mean rounding")
    if difference_gcd:
        require(difference_gcd % 2 == 0, "odd response layer lattice is even")
        require(all((value - anchor) % difference_gcd == 0 for value in image),
                "single layer response coset")
    else:
        require(len(image) == 1, "zero lattice constant image")
    return {
        "size": size,
        "mean": fraction_pair(mean),
        "maximum": max(image),
        "anchor": anchor,
        "lattice": difference_gcd,
        "lattice_floor": lattice_floor,
        "image": image,
    }


def analyze(adjacency):
    order = len(adjacency)
    values = responses(adjacency)
    base = values[0]
    defect_one = one_defect_count(adjacency)
    linear, curvature, global_lattice = coefficients(values, order)
    require(values[-1] == base, "constant source and sink ears")
    require(all(value % 2 for value in values), "Redei oddness control")

    weight = Fraction(sum(
        curvature[left][right]
        for left, right in combinations(range(order), 2)
    ), 2)
    require(weight.denominator == 1, "integral total cut weight")
    require(2 * weight == (order - 1) * base + defect_one,
            "one-defect cut-weight identity")

    layers = []
    for size in range(order + 1):
        record = layer_data(values, linear, curvature, size)
        repair_multiplicity = (
            comb(order - 2, size - 1) if 1 <= size <= order - 1 else 0
        )
        predicted_sum = (
            comb(order, size) * base
            + repair_multiplicity * ((order - 1) * base + defect_one)
        )
        actual_sum = Fraction(*record["mean"]) * comb(order, size)
        require(actual_sum == predicted_sum, "fixed-cardinality sum identity")
        beta = (
            Fraction(2 * size * (order - size), order * (order - 1))
            if 1 <= size <= order - 1 else Fraction(0)
        )
        predicted_mean = base + beta * weight
        require(record["mean"] == fraction_pair(predicted_mean),
                "fixed-cardinality cut-weight mean")
        if record["lattice"]:
            require(record["lattice"] % global_lattice == 0,
                    "global lattice divides each layer lattice")
        layers.append(record)

    balanced_size = order // 2
    balanced = layers[balanced_size]
    q = order * order // 4
    balanced_bound = (
        (1 + Fraction(q, order)) * base
        + Fraction(q, order * (order - 1)) * defect_one
    )
    require(balanced["mean"] == fraction_pair(balanced_bound),
            "balanced closed-form mean")
    for size in {order // 2, (order + 1) // 2}:
        require(layers[size]["mean"] == fraction_pair(balanced_bound),
                "each central layer has the balanced mean")
        require(layers[size]["maximum"] >= odd_ceiling(balanced_bound),
                "balanced odd-ceiling lower bound")

    reduced_variance_bound = (
        Fraction((order + 1) * (order + 2), 4 * order) * base
        + Fraction(order * (order - 1) + 2, 4 * order * (order - 1))
        * defect_one
    )
    comparison_index = order - 2 if order % 2 == 0 else order - 3
    expected_difference = Fraction(comparison_index, 4 * order) * (
        base + Fraction(defect_one, order - 1)
    )
    require(balanced_bound - reduced_variance_bound == expected_difference,
            "coefficientwise comparison with THM-4115")
    require(expected_difference >= 0, "balanced comparison nonnegative")

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
    total_parents = total_ears = total_layers = total_orderings = 0
    for order in range(2, 6):
        parent_count = 1 << (order * (order - 1) // 2)
        strong_count = 0
        balanced_histogram = {}
        strong_histogram = {}
        equality_count = strong_equality_count = 0
        for code in range(parent_count):
            adjacency = decode(code, order)
            record = analyze(adjacency)
            lattice = record["balanced"]["lattice"]
            balanced_histogram[lattice] = balanced_histogram.get(lattice, 0) + 1
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
            "balanced_lattice_histogram": tuple(sorted(balanced_histogram.items())),
            "strong_balanced_lattice_histogram": tuple(sorted(strong_histogram.items())),
            "oddceil_equality_all_strong": (equality_count, strong_equality_count),
        }
        order_rows.append(row)
        total_parents += row["parents"]
        total_ears += row["ears"]
        total_layers += row["layers"]
        total_orderings += row["ordering_checks"]

    require(tuple((row["parents"], row["strong"]) for row in order_rows)
            == ((2, 0), (8, 2), (64, 24), (1024, 544)),
            "labelled parent and strong census")
    require(tuple(row["balanced_lattice_histogram"] for row in order_rows) == (
        ((2, 2),), ((0, 2), (2, 6)), ((2, 48), (8, 16)), ((2, 1024),)
    ), "balanced layer lattice histograms")
    require(tuple(row["strong_balanced_lattice_histogram"] for row in order_rows) == (
        (), ((0, 2),), ((2, 24),), ((2, 544),)
    ), "strong balanced layer lattice histograms")
    require(tuple(row["oddceil_equality_all_strong"] for row in order_rows) == (
        (2, 0), (2, 2), (0, 0), (24, 24)
    ), "balanced odd-ceiling equality census")
    require((total_parents, total_ears, total_layers, total_orderings)
            == (1098, 33864, 6502, 124468), "primary census totals")

    named = {
        name: compact_named(analyze(decode(code, order)))
        for name, code, order in (
            ("c3", 2, 3),
            ("code8", 8, 5),
            ("code759", 759, 5),
            ("regular_c5", 76, 5),
            ("code40", 40, 5),
            ("code20", 20, 6),
            ("code30571", 30571, 6),
        )
    }
    require(named["c3"]["balanced"][1:3] == ((5, 1), 5),
            "C3 balanced equality")
    require(named["code8"]["H_F1_W"] == named["code759"]["H_F1_W"]
            == (9, 30, (33, 1)), "equal-mean hostile invariants")
    require(named["code8"]["balanced"][1] == named["code759"]["balanced"][1]
            == (144, 5), "equal balanced hostile means")
    require(named["code8"]["balanced"][-1] != named["code759"]["balanced"][-1],
            "equal mean does not determine balanced image")
    require(named["regular_c5"]["balanced"][1:3] == ((42, 1), 43),
            "regular C5 makes odd-ceiling sharp")
    require(Fraction(*named["code8"]["variance_floor"])
            > Fraction(*named["code8"]["balanced"][1]),
            "code8 favors exact full-cube variance floor")
    require(Fraction(*named["code40"]["variance_floor"])
            < Fraction(*named["code40"]["balanced"][1]),
            "code40 favors balanced mean floor")
    require(named["code20"]["layer_maxima"] == (19, 91, 133, 131, 123, 91, 19),
            "code20 nonbalanced global maximum")
    require(named["code20"]["full_maximum"] == 133
            > named["code20"]["balanced"][2] == 131,
            "balanced maximum need not be global")
    require(named["code30571"]["layer_lattices"] == (0, 18, 6, 6, 6, 18, 0),
            "code30571 strict layer lattice refinement")
    require(named["code30571"]["balanced"][1:6]
            == ((738, 5), 189, 159, 6, 153),
            "code30571 balanced lattice rounding")

    ledger = {
        "formula": "sum_|S|=m F=C(n,m)H+C(n-2,m-1)((n-1)H+F1)",
        "balanced_bound": "(1+floor(n^2/4)/n)H+floor(n^2/4)F1/(n(n-1))",
        "layer_lattice": "gcd_exchange_increments_equals_gcd_layer_differences",
        "order_rows": order_rows,
        "totals": (total_parents, total_ears, total_layers, total_orderings),
        "named": named,
        "scope": "one large balanced child and exact layer coset; no interval or maximizing-layer theorem",
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")

    print("status=PASS")
    print("implementation=subset Hamiltonian DP plus literal one-defect orderings")
    print("slice_sum=C(n,m)H+C(n-2,m-1)((n-1)H+F1)")
    print("balanced_bound=(1+floor(n^2/4)/n)H+floor(n^2/4)F1/(n(n-1))")
    print(f"order_rows={order_rows}")
    print(f"totals_parents_ears_layers_orderings={ledger['totals']}")
    print(f"C3_balanced={named['c3']['balanced']}")
    print(f"equal_mean_hostiles_code8_code759={(named['code8']['balanced'], named['code759']['balanced'])}")
    print(f"variance_incomparability_code8_code40={(named['code8']['variance_floor'], named['code8']['balanced'][1], named['code40']['variance_floor'], named['code40']['balanced'][1])}")
    print(f"regular_C5_balanced={named['regular_c5']['balanced']}")
    print(f"code20_layer_maxima={named['code20']['layer_maxima']};balanced_not_global=PASS")
    print(f"code30571_balanced={named['code30571']['balanced']};layer_lattices={named['code30571']['layer_lattices']}")
    print("scope=large balanced child and layer coset only; interval and global-max location remain OPEN")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
