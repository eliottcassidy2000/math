#!/usr/bin/env python3
"""Clean-room permutation audit for THM-4128."""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from json import dumps
from math import gcd


EXPECTED_SEMANTIC = "d4dcc84bddb0e17848c97f578eaf5554191dea2d6c22841c1c9c14a80420b88f"


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def digest(value):
    return sha256(
        dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def pair(value):
    value = Fraction(value)
    return value.numerator, value.denominator


def gcd_all(values):
    result = 0
    for value in values:
        result = gcd(result, abs(value))
    return result


def odd_ceiling(value):
    value = Fraction(value)
    integer = -(-value.numerator // value.denominator)
    return integer if integer % 2 else integer + 1


def coset_ceiling(value, anchor, modulus):
    value = Fraction(value)
    if modulus == 0:
        require(value == anchor, "constant layer floor")
        return anchor
    quotient = (value - anchor) / modulus
    return anchor + modulus * (-(-quotient.numerator // quotient.denominator))


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
        seen = 1 << root
        frontier = seen
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            fresh = adjacency[vertex] & ~seen
            seen |= fresh
            frontier |= fresh
        if seen != full:
            return False
    return True


def defect_exposure(adjacency):
    """Enumerate zero/one-defect orderings, independently of boundary DPs."""
    order = len(adjacency)
    capacity = [[0] * order for _ in range(order)]
    H = 0
    for ordering in permutations(range(order)):
        bad = []
        for left, right in zip(ordering, ordering[1:]):
            if not adjacency[left] & (1 << right):
                bad.append((left, right))
                if len(bad) == 2:
                    break
        if not bad:
            H += 1
            for left, right in zip(ordering, ordering[1:]):
                capacity[left][right] += 1
        elif len(bad) == 1:
            left, right = bad[0]
            require(adjacency[right] & (1 << left), "tournament bad-gap reversal")
            capacity[right][left] += 1
    require(H > 0, "Redei path existence")
    for left, right in combinations(range(order), 2):
        require(not (capacity[left][right] and capacity[right][left]),
                "exposure follows the tournament orientation")
    return H, tuple(tuple(row) for row in capacity)


def directed_cut_values(H, capacity):
    order = len(capacity)
    return tuple(
        H + sum(
            capacity[tail][head]
            for tail in range(order) if mask & (1 << tail)
            for head in range(order) if not mask & (1 << head)
        )
        for mask in range(1 << order)
    )


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


def add_ear(adjacency, cut):
    order = len(adjacency)
    child = list(adjacency) + [0]
    for vertex in range(order):
        if cut & (1 << vertex):
            child[order] |= 1 << vertex
        else:
            child[vertex] |= 1 << order
    return tuple(child)


def literal_responses(adjacency):
    return tuple(
        hamilton_count(add_ear(adjacency, cut))
        for cut in range(1 << len(adjacency))
    )


def analyze(adjacency):
    order = len(adjacency)
    require(order >= 4, "audit order range")
    H, capacity = defect_exposure(adjacency)
    values = directed_cut_values(H, capacity)
    require(values[0] == values[-1] == H, "constant ears")
    require(all(value >= H and value % 2 for value in values), "odd lower support")

    weights = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    for left, right in combinations(range(order), 2):
        weight = Fraction(capacity[left][right] + capacity[right][left], 2)
        weights[left][right] = weights[right][left] = weight
    degrees = tuple(sum(weights[vertex]) for vertex in range(order))
    field = tuple(
        sum(
            weights[vertex][other]
            if adjacency[vertex] & (1 << other)
            else -weights[vertex][other]
            for other in range(order) if other != vertex
        )
        for vertex in range(order)
    )
    require(sum(field) == 0, "zero-sum directed exposure field")
    W = sum(weights[left][right] for left, right in combinations(range(order), 2))

    for mask, actual in enumerate(values):
        quadratic = H + sum(
            field[vertex] for vertex in range(order) if mask & (1 << vertex)
        ) + sum(
            weights[left][right]
            for left, right in combinations(range(order), 2)
            if bool(mask & (1 << left)) != bool(mask & (1 << right))
        )
        require(quadratic == actual, "directed-to-symmetric cut conversion")

    average_edge = Fraction(2 * W, order * (order - 1))
    u = tuple(
        Fraction(degrees[vertex] - Fraction(2 * W, order), order - 2)
        for vertex in range(order)
    )
    residual = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    for left, right in combinations(range(order), 2):
        value = weights[left][right] - average_edge - u[left] - u[right]
        residual[left][right] = residual[right][left] = value
    require(sum(u) == 0 and all(sum(row) == 0 for row in residual),
            "independent Hoeffding zero sums")

    h2 = sum(value * value for value in field)
    u2 = sum(value * value for value in u)
    z2 = sum(
        residual[left][right] ** 2
        for left, right in combinations(range(order), 2)
    )
    C_hd = sum(field[vertex] * degrees[vertex] for vertex in range(order))
    require(sum(field[vertex] * u[vertex] for vertex in range(order))
            == Fraction(C_hd, order - 2), "raw tilt identity")

    edges = tuple(combinations(range(order), 2))
    D4 = sum(
        weights[a][b] * weights[c][d]
        for index, (a, b) in enumerate(edges)
        for c, d in edges[index + 1:]
        if len({a, b, c, d}) == 4
    )
    require(D4 > 0, "positive disjoint exposure")
    gamma = (
        Fraction(W * W, order * (order - 1))
        + Fraction(z2, (order - 2) * (order - 3)) - u2
    )
    require(gamma == Fraction(2 * D4, (order - 2) * (order - 3)),
            "independent curvature collapse")
    theta = Fraction((order - 3) * C_hd, 2 * D4)
    J0 = (
        H + Fraction(order * W, 2 * (order - 1)) + Fraction(h2, 2 * W)
        + Fraction((order - 2) * z2, 2 * (order - 3) * W)
    )
    grid = tuple(range(-(order - 2), order - 1, 2))
    distance = min(abs(Fraction(t) - theta) for t in grid)
    predicted_t = tuple(t for t in grid if abs(Fraction(t) - theta) == distance)

    layers = []
    for size in range(1, order):
        t = order - 2 * size
        layer_values = tuple(
            values[mask] for mask in range(1 << order) if mask.bit_count() == size
        )
        mean = Fraction(sum(layer_values), len(layer_values))
        variance = sum(
            (Fraction(value) - mean) ** 2 for value in layer_values
        ) / len(layer_values)
        floor = mean + variance / (mean - H)
        envelope = (
            J0 + Fraction(C_hd * t, W * (order - 2))
            - Fraction(D4 * t * t, W * (order - 2) * (order - 3))
        )
        require(floor == envelope, "raw-incidence envelope")
        anchor = layer_values[0]
        lattice = gcd_all(value - anchor for value in layer_values)
        odd_floor = odd_ceiling(floor)
        coset_floor = coset_ceiling(floor, anchor, lattice)
        require(max(layer_values) >= coset_floor >= odd_floor,
                "independent rounded support")
        layers.append({
            "size": size,
            "t": t,
            "floor": floor,
            "odd_floor": odd_floor,
            "coset_floor": coset_floor,
            "maximum": max(layer_values),
            "lattice": lattice,
        })

    rational_maximum = max(layer["floor"] for layer in layers)
    actual_t = tuple(sorted(
        layer["t"] for layer in layers if layer["floor"] == rational_maximum
    ))
    require(actual_t == predicted_t, "independent nearest-grid optimizer")
    require(max(layer["odd_floor"] for layer in layers)
            == odd_ceiling(rational_maximum), "independent odd monotonicity")
    central_t = (0,) if order % 2 == 0 else (-1, 1)
    central_optimal = not set(actual_t).isdisjoint(central_t)
    criterion = (
        (order - 3) * abs(C_hd) <= (2 if order % 2 == 0 else 4) * D4
    )
    require(central_optimal == criterion, "independent centrality criterion")
    best_candidate_coset = max(
        layer["coset_floor"] for layer in layers if layer["t"] in actual_t
    )
    return {
        "order": order,
        "strong": is_strong(adjacency),
        "H": H,
        "W": W,
        "D4": D4,
        "C_hd": C_hd,
        "theta": theta,
        "J0": J0,
        "predicted_t": predicted_t,
        "central_optimal": central_optimal,
        "best_odd": odd_ceiling(rational_maximum),
        "best_candidate_coset": best_candidate_coset,
        "best_coset": max(layer["coset_floor"] for layer in layers),
        "global_maximum": max(values),
        "layers": tuple(layers),
    }


def compact(record):
    return {
        "order": record["order"],
        "strong": record["strong"],
        "H_W_D4_Chd": (
            record["H"], pair(record["W"]), pair(record["D4"]),
            pair(record["C_hd"]),
        ),
        "theta": pair(record["theta"]),
        "predicted_t": record["predicted_t"],
        "central_optimal": record["central_optimal"],
        "layers": tuple(
            (
                layer["size"], layer["t"], pair(layer["floor"]),
                layer["odd_floor"], layer["coset_floor"],
                layer["maximum"], layer["lattice"],
            )
            for layer in record["layers"]
        ),
        "global_maximum": record["global_maximum"],
    }


def main():
    order_rows = []
    total_scanned = 0
    first_central_failure = None
    first_coset_reordering = None
    for order in range(4, 7):
        labelled = 1 << (order * (order - 1) // 2)
        strong_count = 0
        central_pass = central_fail = 0
        strong_central_pass = strong_central_fail = 0
        coset_reorder = strong_coset_reorder = 0
        optimizer_histogram = {}
        for code in range(labelled):
            record = analyze(decode(code, order))
            total_scanned += 1
            strong_count += record["strong"]
            central_pass += record["central_optimal"]
            central_fail += not record["central_optimal"]
            if record["strong"]:
                strong_central_pass += record["central_optimal"]
                strong_central_fail += not record["central_optimal"]
            optimizer_histogram[record["predicted_t"]] = (
                optimizer_histogram.get(record["predicted_t"], 0) + 1
            )
            reordered = record["best_coset"] > record["best_candidate_coset"]
            coset_reorder += reordered
            strong_coset_reorder += reordered and record["strong"]
            if first_central_failure is None and not record["central_optimal"]:
                first_central_failure = (order, code, compact(record))
            if first_coset_reordering is None and reordered:
                first_coset_reordering = (order, code, compact(record))
        order_rows.append({
            "order": order,
            "labelled": labelled,
            "strong": strong_count,
            "central_pass_fail": (central_pass, central_fail),
            "strong_central_pass_fail": (
                strong_central_pass, strong_central_fail,
            ),
            "coset_reorder": coset_reorder,
            "strong_coset_reorder": strong_coset_reorder,
            "optimizer_histogram": tuple(sorted(optimizer_histogram.items())),
        })

    require(total_scanned == 64 + 1024 + 32768, "independent census size")
    named = {
        name: compact(analyze(decode(code, order)))
        for name, code, order in (
            ("code2", 2, 6),
            ("code11", 11, 5),
            ("code20", 20, 6),
            ("code30571", 30571, 6),
        )
    }
    require(first_central_failure[:2] == (6, 2), "first noncentral code")
    require(first_coset_reordering[:2] == (6, 140), "first coset reorder code")
    require(named["code2"]["theta"] == (141, 125), "code2 theta")
    require(named["code20"]["layers"][2][5] == 131
            and named["code20"]["global_maximum"] == 133,
            "code20 location hostile")
    for code, order in ((2, 6), (20, 6), (30571, 6)):
        adjacency = decode(code, order)
        H, capacity = defect_exposure(adjacency)
        require(directed_cut_values(H, capacity) == literal_responses(adjacency),
                f"clean-room literal replay {order}:{code}")

    ledger = {
        "theorem": "THM-4128",
        "envelope": "J(t)=J0+C_hd*t/(W(n-2))-D4*t^2/(W(n-2)(n-3))",
        "vertex": "theta=(n-3)C_hd/(2D4); nearest parity-grid t",
        "curvature": "Gamma=2D4/((n-2)(n-3))>0",
        "centrality": "even:(n-3)|C_hd|<=2D4;odd:(n-3)|C_hd|<=4D4",
        "order_rows": order_rows,
        "first_central_failure": first_central_failure,
        "first_coset_reordering": first_coset_reordering,
        "named": named,
        "scope": (
            "rational and ordinary-odd support-floor optimizer; exact layer "
            "cosets may reorder; no actual maximizing layer or interval"
        ),
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")

    print("status=ACCEPT")
    print("implementation=permutation defect<=1 exposure enumeration plus direct directed cuts")
    print("envelope=J(t)=J0+C_hd*t/(W(n-2))-D4*t^2/(W(n-2)(n-3))")
    print("vertex=theta=(n-3)C_hd/(2D4); nearest admissible parity-grid point")
    print("centrality=even:(n-3)|C_hd|<=2D4;odd:(n-3)|C_hd|<=4D4")
    print(f"order_rows={order_rows}")
    print(f"first_central_failure={first_central_failure}")
    print(f"first_coset_reordering={first_coset_reordering}")
    print(f"named={named}")
    print("scope=rational/ordinary-odd floors; exact cosets may reorder; actual maximizing layers remain undetermined")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
