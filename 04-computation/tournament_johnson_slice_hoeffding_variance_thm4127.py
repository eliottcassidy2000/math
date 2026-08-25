#!/usr/bin/env python3
"""Exact Johnson-slice variance and support referee for THM-4127."""

from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
from math import comb, gcd


EXPECTED_SEMANTIC = "5a6ab550db55c636209cba0281d91dd9c2c6161cb5c532f14d78733754fdfd67"


def require(condition, label):
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def digest(value):
    return sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def pair(value):
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
        require(value == anchor, "constant layer support floor")
        return anchor
    return anchor + modulus * integer_ceiling((value - anchor) / modulus)


def decode(code, order):
    adjacency = [0] * order
    bit_index = 0
    for left in range(order):
        for right in range(left + 1, order):
            if (code >> bit_index) & 1:
                adjacency[left] |= 1 << right
            else:
                adjacency[right] |= 1 << left
            bit_index += 1
    return tuple(adjacency)


def is_strong(adjacency):
    order = len(adjacency)
    full = (1 << order) - 1
    for root in range(order):
        reached = 1 << root
        frontier = reached
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            fresh = adjacency[vertex] & ~reached
            reached |= fresh
            frontier |= fresh
        if reached != full:
            return False
    return True


def boundary_data(adjacency):
    """Start/end path counts and ordered exposed-gap counts."""
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

    exposed = [[0] * order for _ in range(order)]
    for left_mask in range(1, full):
        right_mask = full ^ left_mask
        left_todo = left_mask
        while left_todo:
            left_bit = left_todo & -left_todo
            left_todo ^= left_bit
            left = left_bit.bit_length() - 1
            left_count = ends[left_mask][left]
            if not left_count:
                continue
            right_todo = right_mask
            while right_todo:
                right_bit = right_todo & -right_todo
                right_todo ^= right_bit
                right = right_bit.bit_length() - 1
                exposed[left][right] += left_count * starts[right_mask][right]
    return tuple(starts[full]), tuple(ends[full]), tuple(tuple(row) for row in exposed)


def responses_from_boundary(data):
    starts, ends, exposed = data
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
                    value += exposed[left][right]
        values.append(value)
    return tuple(values)


def add_ear(adjacency, cut):
    order = len(adjacency)
    child = list(adjacency) + [0]
    for vertex in range(order):
        if cut & (1 << vertex):
            child[order] |= 1 << vertex
        else:
            child[vertex] |= 1 << order
    return tuple(child)


def hamilton_count(adjacency):
    order = len(adjacency)
    full = (1 << order) - 1
    paths = [[0] * order for _ in range(1 << order)]
    for vertex in range(order):
        paths[1 << vertex][vertex] = 1
    for mask in range(1, full + 1):
        for last, count in enumerate(paths[mask]):
            if not count:
                continue
            available = adjacency[last] & (full ^ mask)
            while available:
                bit = available & -available
                available ^= bit
                paths[mask | bit][bit.bit_length() - 1] += count
    return sum(paths[full])


def literal_responses(adjacency):
    return tuple(
        hamilton_count(add_ear(adjacency, cut))
        for cut in range(1 << len(adjacency))
    )


def cut_coefficients(values, order):
    H = values[0]
    linear = tuple(values[1 << vertex] - H for vertex in range(order))
    weights = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    for left, right in combinations(range(order), 2):
        curvature = (
            values[1 << left] + values[1 << right] - H
            - values[(1 << left) | (1 << right)]
        )
        require(curvature % 2 == 0, "integral symmetric cut weight")
        weight = Fraction(curvature, 2)
        require(weight >= 0, "nonnegative cut weight")
        weights[left][right] = weights[right][left] = weight
    degrees = tuple(sum(weights[vertex]) for vertex in range(order))
    field = tuple(linear[vertex] - degrees[vertex] for vertex in range(order))
    require(sum(field) == 0, "zero-sum orientation field")
    W = sum(weights[left][right] for left, right in combinations(range(order), 2))
    require(W > 0, "positive cut weight")
    for mask, actual in enumerate(values):
        predicted = H + sum(
            field[vertex] for vertex in range(order) if mask & (1 << vertex)
        ) + sum(
            weights[left][right]
            for left, right in combinations(range(order), 2)
            if bool(mask & (1 << left)) != bool(mask & (1 << right))
        )
        require(predicted == actual, "complete cut-field reconstruction")
    return H, field, tuple(tuple(row) for row in weights), degrees, W


def hoeffding_coordinates(field, weights, degrees, W):
    order = len(field)
    if order == 2:
        return None, None
    average_edge = Fraction(2 * W, order * (order - 1))
    u = tuple(
        Fraction(degrees[vertex] - Fraction(2 * W, order), order - 2)
        for vertex in range(order)
    )
    require(sum(u) == 0, "zero-sum edge-star coordinate")
    residual = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    for left, right in combinations(range(order), 2):
        value = weights[left][right] - average_edge - u[left] - u[right]
        residual[left][right] = residual[right][left] = value
    require(all(sum(residual[vertex]) == 0 for vertex in range(order)),
            "zero-row-sum harmonic edge residual")
    q2 = sum(
        weights[left][right] ** 2
        for left, right in combinations(range(order), 2)
    )
    u2 = sum(value * value for value in u)
    z2 = sum(
        residual[left][right] ** 2
        for left, right in combinations(range(order), 2)
    )
    require(
        q2 == Fraction(2 * W * W, order * (order - 1)) + (order - 2) * u2 + z2,
        "orthogonal edge-energy decomposition",
    )
    if order == 3:
        require(z2 == 0, "order-three harmonic edge space vanishes")
    return u, tuple(tuple(row) for row in residual)


def analyze(adjacency):
    order = len(adjacency)
    values = responses_from_boundary(boundary_data(adjacency))
    require(values[0] == values[-1], "constant ears")
    require(all(value >= values[0] for value in values), "ear lower support")
    require(all(value % 2 for value in values), "Redei odd response")
    H, field, weights, degrees, W = cut_coefficients(values, order)
    coordinates = hoeffding_coordinates(field, weights, degrees, W)
    u, residual = coordinates if coordinates != (None, None) else (None, None)

    h2 = sum(value * value for value in field)
    d2 = sum(value * value for value in degrees)
    q2 = sum(
        weights[left][right] ** 2
        for left, right in combinations(range(order), 2)
    )
    chd = sum(field[vertex] * degrees[vertex] for vertex in range(order))
    z2 = (
        sum(residual[left][right] ** 2
            for left, right in combinations(range(order), 2))
        if residual is not None else Fraction(0)
    )
    cube_mean = Fraction(sum(values), len(values))
    cube_variance = sum((Fraction(value) - cube_mean) ** 2 for value in values) / len(values)
    require(cube_mean == H + W / 2, "full-cube mean")
    require(cube_variance == (h2 + q2) / 4, "full-cube Walsh variance")
    cube_floor = cube_mean + cube_variance / (cube_mean - H)

    disjoint_energy = sum(
        weights[a][b] * weights[c][d]
        for a, b in combinations(range(order), 2)
        for c, d in combinations(range(order), 2)
        if (a, b) < (c, d) and len({a, b, c, d}) == 4
    )
    require(2 * disjoint_energy == W * W + q2 - d2,
            "disjoint-edge energy identity")

    layers = []
    for size in range(order + 1):
        masks = tuple(mask for mask in range(1 << order) if mask.bit_count() == size)
        layer_values = tuple(values[mask] for mask in masks)
        mean = Fraction(sum(layer_values), len(layer_values))
        variance = sum((Fraction(value) - mean) ** 2 for value in layer_values) / len(layer_values)
        maximum = max(layer_values)
        anchor = layer_values[0]
        lattice = gcd_all(value - anchor for value in layer_values)
        require(lattice == 0 or lattice % 2 == 0, "odd layer response coset")

        if size in (0, order):
            require((mean, variance, maximum, lattice) == (H, 0, H, 0),
                    "endpoint layer boundary")
            layers.append({
                "size": size, "mean": pair(mean), "variance": pair(variance),
                "maximum": maximum, "lattice": lattice,
            })
            continue

        complement = order - size
        alpha = Fraction(size * complement, order * (order - 1))
        require(mean == H + 2 * alpha * W, "Johnson-slice mean")
        if order == 2:
            predicted_variance = h2 / 2
            g = field
        else:
            g = tuple(
                field[vertex] + (order - 2 * size) * u[vertex]
                for vertex in range(order)
            )
            for mask, actual in zip(masks, layer_values):
                predicted = mean + sum(
                    g[vertex] for vertex in range(order) if mask & (1 << vertex)
                ) - 2 * sum(
                    residual[left][right]
                    for left, right in combinations(range(order), 2)
                    if mask & (1 << left) and mask & (1 << right)
                )
                require(predicted == actual, "pointwise Johnson-Hoeffding decomposition")
            if order == 3:
                predicted_variance = alpha * sum(value * value for value in g)
            else:
                c2 = Fraction(
                    size * (size - 1) * complement * (complement - 1),
                    order * (order - 1) * (order - 2) * (order - 3),
                )
                predicted_variance = (
                    alpha * sum(value * value for value in g) + 4 * c2 * z2
                )
                delta = 4 * c2
                raw_variance = (
                    alpha * h2
                    + Fraction(2 * alpha * (order - 2 * size), order - 2) * chd
                    + delta * q2 + (alpha - delta) * d2
                    + (delta - 4 * alpha * alpha) * W * W
                )
                require(raw_variance == predicted_variance,
                        "raw aggregate slice variance")
        require(variance == predicted_variance, "exact Johnson-slice variance")

        support_floor = mean + variance / (mean - H)
        require(maximum >= support_floor, "slice support inequality")
        support_equality = maximum == support_floor
        endpoint_support = all(value in (H, maximum) for value in layer_values)
        require(support_equality == endpoint_support,
                "support inequality equality classification")
        if order >= 4:
            closed_floor = (
                mean + sum(value * value for value in g) / (2 * W)
                + Fraction(
                    2 * (size - 1) * (complement - 1),
                    (order - 2) * (order - 3) * W,
                ) * z2
            )
            require(closed_floor == support_floor, "closed Johnson support floor")
        odd_floor = odd_ceiling(support_floor)
        lattice_floor = coset_ceiling(support_floor, anchor, lattice)
        require(maximum >= lattice_floor >= odd_floor,
                "Johnson-coset support rounding")
        layers.append({
            "size": size,
            "mean": pair(mean),
            "variance": pair(variance),
            "maximum": maximum,
            "lattice": lattice,
            "support_floor": pair(support_floor),
            "odd_floor": odd_floor,
            "coset_floor": lattice_floor,
            "support_equality": support_equality,
        })

    weighted_slice_variance = sum(
        Fraction(comb(order, size), 1 << order) * Fraction(*layers[size]["variance"])
        for size in range(order + 1)
    )
    require(
        cube_variance == weighted_slice_variance + Fraction(W * W, 2 * order * (order - 1)),
        "exact total-variance relation",
    )

    central_sizes = tuple(sorted({order // 2, (order + 1) // 2}))
    central = tuple(layers[size] for size in central_sizes)
    best_central_floor = max(Fraction(*layer["support_floor"]) for layer in central)
    best_central_odd = max(layer["odd_floor"] for layer in central)
    best_central_coset = max(layer["coset_floor"] for layer in central)
    central_maximum = max(layer["maximum"] for layer in central)
    balanced_mean = Fraction(*central[0]["mean"])
    require(all(Fraction(*layer["mean"]) == balanced_mean for layer in central),
            "equal central means")

    if order == 2:
        require(best_central_floor == cube_floor, "order-two slice/cube boundary")
    elif order == 3:
        require(z2 == 0, "order-three degree-two boundary")
    elif order % 2 == 0:
        layer = central[0]
        difference = Fraction(*layer["support_floor"]) - cube_floor
        require(difference == Fraction(disjoint_energy, W * (order - 3)),
                "even central-to-cube surplus")
    else:
        inner = sum(field[vertex] * u[vertex] for vertex in range(order))
        lower_difference = Fraction(*central[0]["support_floor"]) - cube_floor
        upper_difference = Fraction(*central[1]["support_floor"]) - cube_floor
        require(
            lower_difference == Fraction(disjoint_energy, W * (order - 2)) + inner / W,
            "odd lower-central surplus",
        )
        require(
            upper_difference == Fraction(disjoint_energy, W * (order - 2)) - inner / W,
            "odd upper-central surplus",
        )
        require(
            best_central_floor - cube_floor
            == Fraction(disjoint_energy, W * (order - 2)) + abs(inner) / W,
            "best odd-central dominance",
        )

    return {
        "order": order,
        "H": H,
        "W": pair(W),
        "cube_floor": pair(cube_floor),
        "cube_odd_floor": odd_ceiling(cube_floor),
        "balanced_mean": pair(balanced_mean),
        "central": central,
        "best_central_floor": pair(best_central_floor),
        "best_central_odd": best_central_odd,
        "best_central_coset": best_central_coset,
        "central_maximum": central_maximum,
        "global_maximum": max(values),
        "disjoint_energy": pair(disjoint_energy),
        "field_degree_inner": pair(chd),
    }


def compact_named(record):
    lower = record["central"][0]
    return {
        "order": record["order"],
        "H_W": (record["H"], record["W"]),
        "lower_central": (
            lower["size"], lower["mean"], lower["variance"],
            lower["support_floor"], lower["odd_floor"], lower["coset_floor"],
            lower["maximum"], lower["lattice"],
        ),
        "upper_variance": record["central"][-1]["variance"],
        "cube": (record["cube_floor"], record["cube_odd_floor"]),
        "central_global_maximum": (
            record["central_maximum"], record["global_maximum"]
        ),
    }


def main():
    require(literal_responses((0,)) == (1, 1), "literal order-one constant boundary")
    order_rows = []
    scanned = profiled = 0
    for order in range(2, 7):
        tournament_count = 1 << (order * (order - 1) // 2)
        strong_count = 0
        profiled_rational_strict = profiled_rational_equal = 0
        rational_strict = rational_equal = 0
        rounded_strict = rounded_equal = 0
        mean_odd_improve = mean_odd_tie = 0
        mean_coset_improve = mean_coset_tie = 0
        floor_attained = 0
        for code in range(tournament_count):
            adjacency = decode(code, order)
            scanned += 1
            strong = is_strong(adjacency)
            if order == 6 and not strong:
                continue
            record = analyze(adjacency)
            profiled += 1
            profiled_floor = Fraction(*record["best_central_floor"])
            profiled_cube = Fraction(*record["cube_floor"])
            profiled_rational_strict += profiled_floor > profiled_cube
            profiled_rational_equal += profiled_floor == profiled_cube
            if not strong:
                continue
            strong_count += 1
            central_floor = Fraction(*record["best_central_floor"])
            cube_floor = Fraction(*record["cube_floor"])
            rational_strict += central_floor > cube_floor
            rational_equal += central_floor == cube_floor
            rounded_strict += record["best_central_coset"] > record["cube_odd_floor"]
            rounded_equal += record["best_central_coset"] == record["cube_odd_floor"]
            # Recompute the mean rounding from any actual layer value; maximum
            # is a valid anchor for the same coset and avoids storing masks.
            mean_coset = max(
                coset_ceiling(Fraction(*layer["mean"]), layer["maximum"], layer["lattice"])
                for layer in record["central"]
            )
            mean_odd = odd_ceiling(Fraction(*record["balanced_mean"]))
            mean_odd_improve += record["best_central_odd"] > mean_odd
            mean_odd_tie += record["best_central_odd"] == mean_odd
            mean_coset_improve += record["best_central_coset"] > mean_coset
            mean_coset_tie += record["best_central_coset"] == mean_coset
            floor_attained += record["central_maximum"] == record["best_central_coset"]

        row = {
            "order": order,
            "labelled": tournament_count,
            "strong": strong_count,
            "profiled": tournament_count if order < 6 else strong_count,
            "profiled_central_vs_cube_rational": (
                profiled_rational_strict, profiled_rational_equal
            ),
            "strong_central_vs_cube_rational": (rational_strict, rational_equal),
            "strong_central_vs_cube_rounded": (rounded_strict, rounded_equal),
            "strong_slice_vs_mean_odd": (mean_odd_improve, mean_odd_tie),
            "strong_slice_vs_mean_coset": (mean_coset_improve, mean_coset_tie),
            "strong_slice_floor_attained": floor_attained,
        }
        order_rows.append(row)

    require(scanned == 33866 and profiled == 1098 + 22320,
            "labelled scan and exact-profile totals")
    require(tuple((row["labelled"], row["strong"]) for row in order_rows)
            == ((2, 0), (8, 2), (64, 24), (1024, 544), (32768, 22320)),
            "labelled and strong tournament census")
    require(order_rows[0]["profiled_central_vs_cube_rational"] == (0, 2)
            and order_rows[1]["profiled_central_vs_cube_rational"] == (0, 8),
            "all-labelled order-two/three slice-cube equality")
    require(order_rows[1]["strong_central_vs_cube_rational"] == (0, 2),
            "order-three exact equality")
    require(order_rows[2]["strong_central_vs_cube_rational"] == (24, 0)
            and order_rows[2]["strong_central_vs_cube_rounded"] == (0, 24),
            "order-four rational gain and rounded tie")
    require(order_rows[3]["strong_central_vs_cube_rational"] == (544, 0)
            and order_rows[3]["strong_central_vs_cube_rounded"] == (544, 0)
            and order_rows[3]["strong_slice_vs_mean_odd"] == (480, 64)
            and order_rows[3]["strong_slice_vs_mean_coset"] == (480, 64)
            and order_rows[3]["strong_slice_floor_attained"] == 24,
            "order-five strict and sharp census")
    require(order_rows[4]["strong_central_vs_cube_rational"] == (22320, 0)
            and order_rows[4]["strong_central_vs_cube_rounded"] == (22320, 0)
            and order_rows[4]["strong_slice_vs_mean_odd"] == (22080, 240)
            and order_rows[4]["strong_slice_vs_mean_coset"] == (21840, 480),
            "order-six strict census")

    named = {
        name: compact_named(analyze(decode(code, order)))
        for name, code, order in (
            ("c3", 2, 3), ("code8", 8, 5), ("code40", 40, 5),
            ("regular_c5", 76, 5), ("code759", 759, 5),
            ("code20", 20, 6), ("code30571", 30571, 6),
        )
    }
    expected_named_core = {
        "c3": ((5, 1), (0, 1), (5, 1), 5, 5, 5, (5, 1), 5, (5, 5)),
        "code8": ((144, 5), (1469, 25), (3145, 99), 33, 33, 41, (994, 33), 31, (41, 41)),
        "code40": ((192, 5), (261, 25), (505, 13), 39, 39, 43, (479, 13), 37, (43, 43)),
        "regular_c5": ((42, 1), (1, 1), (1135, 27), 43, 43, 43, (358, 9), 41, (43, 43)),
        "code759": ((144, 5), (1529, 25), (287, 9), 33, 33, 43, (333, 11), 31, (43, 43)),
        "code20": ((422, 5), (17991, 25), (10399, 109), 97, 97, 131, (9636, 109), 89, (131, 133)),
        "code30571": ((738, 5), (4401, 25), (2837, 19), 151, 153, 189, (2619, 19), 139, (189, 189)),
    }
    for name, expected in expected_named_core.items():
        lower = named[name]["lower_central"]
        actual = (
            lower[1], lower[2], lower[3], lower[4], lower[5], lower[6],
            named[name]["cube"][0], named[name]["cube"][1],
            named[name]["central_global_maximum"],
        )
        require(actual == expected, f"named central control {name}")

    for name, code, order in (
        ("c3", 2, 3), ("code20", 20, 6), ("code30571", 30571, 6)
    ):
        adjacency = decode(code, order)
        require(responses_from_boundary(boundary_data(adjacency)) == literal_responses(adjacency),
                f"literal child replay {name}")

    hostile_record = analyze(decode(11, 5))
    odd_layer_hostile = (
        11, hostile_record["H"], hostile_record["W"],
        hostile_record["disjoint_energy"], hostile_record["field_degree_inner"],
        hostile_record["central"][0]["support_floor"],
        hostile_record["cube_floor"],
        hostile_record["central"][1]["support_floor"],
    )
    require(odd_layer_hostile == (
        11, 3, (21, 1), (60, 1), (-66, 1),
        (421, 21), (141, 7), (155, 7),
    ), "one odd central layer can lose to the cube floor")

    first_odd_split = None
    preceding_strong_ties = 0
    for code in range(35):
        adjacency = decode(code, 7)
        if not is_strong(adjacency):
            continue
        record = analyze(adjacency)
        lower_variance = record["central"][0]["variance"]
        upper_variance = record["central"][1]["variance"]
        if lower_variance != upper_variance:
            first_odd_split = (
                code, preceding_strong_ties, record["balanced_mean"],
                lower_variance, upper_variance,
            )
            break
        preceding_strong_ties += 1
    require(first_odd_split == (34, 1, (1209, 7), (234944, 49), (1205856, 245)),
            "first strong odd central-variance split")

    ledger = {
        "hoeffding": "F=H+2mqW/(n(n-1))+sum_S(h+(n-2m)u)-2sum_pairs_S z",
        "variance": "alpha*||h+(n-2m)u||^2+4*c2*||z||^2",
        "support": "J_m=mu_m+Var_m/(mu_m-H) with Johnson-coset rounding",
        "central_dominance": "even:D4/(W(n-3));odd_best:D4/(W(n-2))+abs(<h,u>)/W",
        "total_variance": "V_cube=sum_m binomial_weight*V_m+W^2/(2n(n-1))",
        "order_rows": order_rows,
        "named": named,
        "odd_layer_hostile": odd_layer_hostile,
        "first_odd_split": first_odd_split,
        "scope": "central support floors; no slice interval or balanced global maximizer",
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")

    print("status=PASS")
    print("implementation=Start/End/exposed-gap DP plus exact Johnson-Hoeffding algebra")
    print("variance=alpha*||h+(n-2m)u||^2+4*c2*||z||^2")
    print("central_dominance=even:D4/(W(n-3));odd_best:D4/(W(n-2))+abs(<h,u>)/W")
    print(f"order_rows={order_rows}")
    print(f"totals_scanned_profiled={(scanned, profiled)}")
    print(f"named_controls={named}")
    print(f"odd_layer_hostile={odd_layer_hostile}")
    print(f"first_odd_variance_split={first_odd_split}")
    print("scope=central support floors only; slice intervals and balanced global maximizers remain OPEN")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
