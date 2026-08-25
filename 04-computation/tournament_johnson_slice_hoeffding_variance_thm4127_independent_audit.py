#!/usr/bin/env python3
"""Clean-room boundary-state and incidence audit for THM-4127."""

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
    return reduce(gcd, (abs(int(value)) for value in values), 0)


def ceiling(value):
    value = Fraction(value)
    return -((-value.numerator) // value.denominator)


def odd_ceiling(value):
    result = ceiling(value)
    return result if result % 2 else result + 1


def coset_ceiling(value, anchor, modulus):
    if modulus == 0:
        require(Fraction(value) == anchor, "clean-room constant slice")
        return anchor
    return anchor + modulus * ceiling((Fraction(value) - anchor) / modulus)


def tournament(order, code):
    rows = [0] * order
    place = 0
    for small in range(order):
        for large in range(small + 1, order):
            if code & (1 << place):
                rows[small] |= 1 << large
            else:
                rows[large] |= 1 << small
            place += 1
    return tuple(rows)


def reverse_tournament(rows):
    order = len(rows)
    reverse = [0] * order
    for source in range(order):
        targets = rows[source]
        while targets:
            bit = targets & -targets
            targets ^= bit
            reverse[bit.bit_length() - 1] |= 1 << source
    return tuple(reverse)


def reach(rows, root):
    seen = 1 << root
    active = seen
    while active:
        bit = active & -active
        active ^= bit
        fresh = rows[bit.bit_length() - 1] & ~seen
        seen |= fresh
        active |= fresh
    return seen


def is_strong(rows):
    full = (1 << len(rows)) - 1
    return reach(rows, 0) == full and reach(reverse_tournament(rows), 0) == full


def end_table(rows):
    order = len(rows)
    table = [[0] * order for _ in range(1 << order)]
    for vertex in range(order):
        table[1 << vertex][vertex] = 1
    for mask in range(1, 1 << order):
        for last, count in enumerate(table[mask]):
            if not count:
                continue
            available = rows[last] & ~mask
            while available:
                bit = available & -available
                available ^= bit
                table[mask | bit][bit.bit_length() - 1] += count
    return table


def clean_room_cut_field(rows):
    order = len(rows)
    full = (1 << order) - 1
    forward = end_table(rows)
    backward = end_table(reverse_tournament(rows))
    ends = tuple(forward[full])
    starts = tuple(backward[full])
    H = sum(ends)
    require(H == sum(starts), "clean-room start/end Hamilton count")

    exposed = [[0] * order for _ in range(order)]
    for left in range(order):
        for right in range(order):
            if left == right:
                continue
            free = full ^ (1 << left) ^ (1 << right)
            subset = free
            while True:
                left_mask = subset | (1 << left)
                right_mask = full ^ left_mask
                exposed[left][right] += (
                    forward[left_mask][left] * backward[right_mask][right]
                )
                if subset == 0:
                    break
                subset = (subset - 1) & free

    weights = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    for left, right in combinations(range(order), 2):
        symmetric = exposed[left][right] + exposed[right][left]
        require(symmetric % 2 == 0, "clean-room integral edge weight")
        value = Fraction(symmetric, 2)
        weights[left][right] = weights[right][left] = value
    field = []
    for vertex in range(order):
        outgoing = sum(exposed[vertex])
        incoming = sum(exposed[other][vertex] for other in range(order))
        require((incoming - outgoing) % 2 == 0,
                "clean-room integral orientation field")
        field.append(starts[vertex] - ends[vertex] + Fraction(incoming - outgoing, 2))
    require(sum(field) == 0, "clean-room zero field sum")
    return H, tuple(field), tuple(tuple(row) for row in weights)


def polynomial_values(H, field, weights):
    order = len(field)
    values = []
    for mask in range(1 << order):
        value = H
        value += sum(field[i] for i in range(order) if mask & (1 << i))
        value += sum(
            weights[i][j]
            for i, j in combinations(range(order), 2)
            if bool(mask & (1 << i)) != bool(mask & (1 << j))
        )
        require(value.denominator == 1, "clean-room integral response")
        values.append(value.numerator)
    return tuple(values)


def add_ear(rows, mask):
    order = len(rows)
    child = list(rows) + [0]
    for vertex in range(order):
        if mask & (1 << vertex):
            child[order] |= 1 << vertex
        else:
            child[vertex] |= 1 << order
    return tuple(child)


def literal_values(rows):
    return tuple(
        sum(end_table(add_ear(rows, mask))[-1])
        for mask in range(1 << len(rows))
    )


def audit(rows):
    order = len(rows)
    H, field, weights = clean_room_cut_field(rows)
    values = polynomial_values(H, field, weights)
    require(values[0] == values[-1] == H, "clean-room constant ears")
    require(min(values) == H, "clean-room lower support")
    require(all(value % 2 for value in values), "clean-room Redei parity")

    edges = tuple(combinations(range(order), 2))
    W = sum(weights[i][j] for i, j in edges)
    require(W > 0, "clean-room positive weight")
    degrees = tuple(sum(weights[i]) for i in range(order))
    h2 = sum(value * value for value in field)
    d2 = sum(value * value for value in degrees)
    q2 = sum(weights[i][j] ** 2 for i, j in edges)
    chd = sum(field[i] * degrees[i] for i in range(order))

    if order >= 3:
        edge_average = Fraction(2 * W, order * (order - 1))
        u = tuple(
            (degrees[i] - Fraction(2 * W, order)) / (order - 2)
            for i in range(order)
        )
        z = [[Fraction(0) for _ in range(order)] for _ in range(order)]
        for i, j in edges:
            z[i][j] = z[j][i] = weights[i][j] - edge_average - u[i] - u[j]
        require(sum(u) == 0 and all(sum(row) == 0 for row in z),
                "clean-room harmonic constraints")
        z2 = sum(z[i][j] ** 2 for i, j in edges)
        require(q2 == Fraction(2 * W * W, order * (order - 1))
                + (order - 2) * sum(value * value for value in u) + z2,
                "clean-room edge ANOVA energy")
        if order == 3:
            require(z2 == 0, "clean-room order-three edge boundary")
    else:
        u = None
        z = None
        z2 = Fraction(0)

    cube_mean = Fraction(sum(values), len(values))
    cube_variance = sum((Fraction(value) - cube_mean) ** 2 for value in values) / len(values)
    require(cube_mean == H + W / 2 and cube_variance == (h2 + q2) / 4,
            "clean-room full-cube moments")
    cube_floor = cube_mean + cube_variance / (cube_mean - H)

    edge_list = list(edges)
    D4 = sum(
        weights[edge_list[first][0]][edge_list[first][1]]
        * weights[edge_list[second][0]][edge_list[second][1]]
        for first in range(len(edge_list))
        for second in range(first + 1, len(edge_list))
        if set(edge_list[first]).isdisjoint(edge_list[second])
    )
    require(2 * D4 == W * W + q2 - d2, "clean-room D4 identity")

    layers = []
    for size in range(order + 1):
        masks = tuple(mask for mask in range(1 << order) if mask.bit_count() == size)
        image = tuple(values[mask] for mask in masks)
        mean = Fraction(sum(image), len(image))
        variance = sum((Fraction(value) - mean) ** 2 for value in image) / len(image)
        maximum = max(image)
        anchor = image[0]
        lattice = gcd_all(value - anchor for value in image)
        if size in (0, order):
            require((mean, variance, maximum, lattice) == (H, 0, H, 0),
                    "clean-room endpoint slice")
            layers.append({
                "size": size, "mean": pair(mean), "variance": pair(variance),
                "maximum": maximum, "lattice": lattice,
            })
            continue

        q = order - size
        alpha = Fraction(size * q, order * (order - 1))
        require(mean == H + 2 * alpha * W, "clean-room slice mean")
        if order == 2:
            predicted = h2 / 2
            g = field
        else:
            g = tuple(field[i] + (order - 2 * size) * u[i] for i in range(order))
            for mask, actual in zip(masks, image):
                reconstructed = mean
                reconstructed += sum(g[i] for i in range(order) if mask & (1 << i))
                reconstructed -= 2 * sum(
                    z[i][j] for i, j in edges
                    if mask & (1 << i) and mask & (1 << j)
                )
                require(reconstructed == actual, "clean-room pointwise Hoeffding law")
            if order == 3:
                predicted = alpha * sum(value * value for value in g)
            else:
                c2 = Fraction(
                    size * (size - 1) * q * (q - 1),
                    order * (order - 1) * (order - 2) * (order - 3),
                )
                predicted = alpha * sum(value * value for value in g) + 4 * c2 * z2
                delta = 4 * c2
                raw = (
                    alpha * h2
                    + Fraction(2 * alpha * (order - 2 * size), order - 2) * chd
                    + delta * q2 + (alpha - delta) * d2
                    + (delta - 4 * alpha * alpha) * W * W
                )
                require(raw == predicted, "clean-room raw incidence variance")
        require(variance == predicted, "clean-room slice variance")

        floor = mean + variance / (mean - H)
        require(maximum >= floor, "clean-room slice support floor")
        equality = maximum == floor
        require(equality == all(value in (H, maximum) for value in image),
                "clean-room support equality case")
        if order >= 4:
            closed = (
                mean + sum(value * value for value in g) / (2 * W)
                + Fraction(2 * (size - 1) * (q - 1), (order - 2) * (order - 3))
                * z2 / W
            )
            require(closed == floor, "clean-room closed support formula")
        odd_floor = odd_ceiling(floor)
        coset_floor = coset_ceiling(floor, anchor, lattice)
        require(maximum >= coset_floor >= odd_floor, "clean-room coset rounding")
        layers.append({
            "size": size, "mean": pair(mean), "variance": pair(variance),
            "maximum": maximum, "lattice": lattice,
            "support_floor": pair(floor), "odd_floor": odd_floor,
            "coset_floor": coset_floor, "support_equality": equality,
        })

    weighted_within = sum(
        Fraction(comb(order, size), 1 << order) * Fraction(*layers[size]["variance"])
        for size in range(order + 1)
    )
    require(cube_variance == weighted_within + Fraction(W * W, 2 * order * (order - 1)),
            "clean-room law of total variance")

    central_indices = tuple(sorted({order // 2, (order + 1) // 2}))
    central = tuple(layers[index] for index in central_indices)
    mean = Fraction(*central[0]["mean"])
    require(all(Fraction(*layer["mean"]) == mean for layer in central),
            "clean-room common central mean")
    best_floor = max(Fraction(*layer["support_floor"]) for layer in central)
    if order == 2:
        require(best_floor == cube_floor, "clean-room order-two boundary")
    elif order >= 4 and order % 2 == 0:
        require(Fraction(*central[0]["support_floor"]) - cube_floor
                == D4 / (W * (order - 3)), "clean-room even dominance")
    elif order >= 5:
        hu = sum(field[i] * u[i] for i in range(order))
        require(Fraction(*central[0]["support_floor"]) - cube_floor
                == D4 / (W * (order - 2)) + hu / W,
                "clean-room odd lower surplus")
        require(Fraction(*central[1]["support_floor"]) - cube_floor
                == D4 / (W * (order - 2)) - hu / W,
                "clean-room odd upper surplus")
        require(best_floor - cube_floor == D4 / (W * (order - 2)) + abs(hu) / W,
                "clean-room best odd dominance")

    return {
        "order": order,
        "H": H,
        "W": pair(W),
        "cube_floor": pair(cube_floor),
        "cube_odd_floor": odd_ceiling(cube_floor),
        "balanced_mean": pair(mean),
        "central": central,
        "best_central_floor": pair(best_floor),
        "best_central_odd": max(layer["odd_floor"] for layer in central),
        "best_central_coset": max(layer["coset_floor"] for layer in central),
        "central_maximum": max(layer["maximum"] for layer in central),
        "global_maximum": max(values),
        "disjoint_energy": pair(D4),
        "field_degree_inner": pair(chd),
    }


def compact(record):
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
        "central_global_maximum": (record["central_maximum"], record["global_maximum"]),
    }


def main():
    require(literal_values((0,)) == (1, 1),
            "clean-room literal order-one boundary")
    order_rows = []
    scanned = profiled = 0
    for order in range(2, 7):
        labelled = 1 << (order * (order - 1) // 2)
        strong_count = 0
        profiled_rational = [0, 0]
        rational = [0, 0]
        rounded = [0, 0]
        mean_odd = [0, 0]
        mean_coset = [0, 0]
        attained = 0
        for code in range(labelled):
            rows = tournament(order, code)
            scanned += 1
            strong = is_strong(rows)
            if order == 6 and not strong:
                continue
            record = audit(rows)
            profiled += 1
            profiled_floor = Fraction(*record["best_central_floor"])
            profiled_cube = Fraction(*record["cube_floor"])
            profiled_rational[0] += profiled_floor > profiled_cube
            profiled_rational[1] += profiled_floor == profiled_cube
            if not strong:
                continue
            strong_count += 1
            central_floor = Fraction(*record["best_central_floor"])
            cube_floor = Fraction(*record["cube_floor"])
            rational[0] += central_floor > cube_floor
            rational[1] += central_floor == cube_floor
            rounded[0] += record["best_central_coset"] > record["cube_odd_floor"]
            rounded[1] += record["best_central_coset"] == record["cube_odd_floor"]
            old_odd = odd_ceiling(Fraction(*record["balanced_mean"]))
            old_coset = max(
                coset_ceiling(Fraction(*layer["mean"]), layer["maximum"], layer["lattice"])
                for layer in record["central"]
            )
            mean_odd[0] += record["best_central_odd"] > old_odd
            mean_odd[1] += record["best_central_odd"] == old_odd
            mean_coset[0] += record["best_central_coset"] > old_coset
            mean_coset[1] += record["best_central_coset"] == old_coset
            attained += record["central_maximum"] == record["best_central_coset"]
        order_rows.append({
            "order": order, "labelled": labelled, "strong": strong_count,
            "profiled": labelled if order < 6 else strong_count,
            "profiled_central_vs_cube_rational": tuple(profiled_rational),
            "strong_central_vs_cube_rational": tuple(rational),
            "strong_central_vs_cube_rounded": tuple(rounded),
            "strong_slice_vs_mean_odd": tuple(mean_odd),
            "strong_slice_vs_mean_coset": tuple(mean_coset),
            "strong_slice_floor_attained": attained,
        })

    require((scanned, profiled) == (33866, 23418), "clean-room census totals")
    require(tuple((row["labelled"], row["strong"]) for row in order_rows)
            == ((2, 0), (8, 2), (64, 24), (1024, 544), (32768, 22320)),
            "clean-room labelled/strong census")
    require(order_rows[0]["profiled_central_vs_cube_rational"] == (0, 2)
            and order_rows[1]["profiled_central_vs_cube_rational"] == (0, 8),
            "clean-room all-labelled small-order equality")
    require(order_rows[1]["strong_central_vs_cube_rational"] == (0, 2),
            "clean-room order-three equality")
    require(order_rows[2]["strong_central_vs_cube_rational"] == (24, 0)
            and order_rows[2]["strong_central_vs_cube_rounded"] == (0, 24),
            "clean-room order-four census")
    require(order_rows[3]["strong_central_vs_cube_rational"] == (544, 0)
            and order_rows[3]["strong_central_vs_cube_rounded"] == (544, 0)
            and order_rows[3]["strong_slice_vs_mean_odd"] == (480, 64)
            and order_rows[3]["strong_slice_vs_mean_coset"] == (480, 64)
            and order_rows[3]["strong_slice_floor_attained"] == 24,
            "clean-room order-five census")
    require(order_rows[4]["strong_central_vs_cube_rational"] == (22320, 0)
            and order_rows[4]["strong_central_vs_cube_rounded"] == (22320, 0)
            and order_rows[4]["strong_slice_vs_mean_odd"] == (22080, 240)
            and order_rows[4]["strong_slice_vs_mean_coset"] == (21840, 480),
            "clean-room order-six census")

    named_specs = (
        ("c3", 2, 3), ("code8", 8, 5), ("code40", 40, 5),
        ("regular_c5", 76, 5), ("code759", 759, 5),
        ("code20", 20, 6), ("code30571", 30571, 6),
    )
    named = {name: compact(audit(tournament(order, code)))
             for name, code, order in named_specs}
    for name, code, order in (
        ("c3", 2, 3), ("code20", 20, 6), ("code30571", 30571, 6)
    ):
        rows = tournament(order, code)
        H, field, weights = clean_room_cut_field(rows)
        require(polynomial_values(H, field, weights) == literal_values(rows),
                f"clean-room literal child replay {name}")

    hostile_record = audit(tournament(5, 11))
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
    ), "clean-room one-sided odd central hostile")

    first_odd_split = None
    preceding_strong_ties = 0
    for code in range(35):
        rows = tournament(7, code)
        if not is_strong(rows):
            continue
        record = audit(rows)
        low = record["central"][0]["variance"]
        high = record["central"][1]["variance"]
        if low != high:
            first_odd_split = (
                code, preceding_strong_ties, record["balanced_mean"], low, high
            )
            break
        preceding_strong_ties += 1
    require(first_odd_split == (34, 1, (1209, 7), (234944, 49), (1205856, 245)),
            "clean-room first odd variance split")

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
        require(semantic == EXPECTED_SEMANTIC, "shared frozen semantic digest")

    print("status=PASS")
    print("implementation=clean-room forward/reverse path tables plus raw slice incidences")
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
