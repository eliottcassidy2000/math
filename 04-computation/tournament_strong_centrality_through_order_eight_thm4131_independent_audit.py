#!/usr/bin/env python3
"""Pure-standard-library independent audit for THM-4131.

Tournament isomorphism classes are generated recursively by vertex extension
and invariant canonical certificates.  Ear exposures are computed without
the primary exposed-gap table: for every tournament arc, two contracted
adjacent blocks are counted, the good block and its reversed one-defect
block.  The resulting directed capacities reconstruct every ear response.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from json import dumps
from math import gcd


EXPECTED_SEMANTIC = "8e5ef116e577c3d3ab5dd4bea50953581a457ffbaae4f80504060c8793fed578"
EXPECTED_ALL_COUNTS = {4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
EXPECTED_STRONG_COUNTS = {4: 1, 5: 6, 6: 35, 7: 353, 8: 6008}


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def digest(value):
    return sha256(
        dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def raw_digest(lines):
    return sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def pair(value):
    value = Fraction(value)
    return value.numerator, value.denominator


def gcd_all(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def coset_ceiling(value, anchor, modulus):
    value = Fraction(value)
    if modulus == 0:
        require(value == anchor, "constant-layer floor")
        return anchor
    quotient = (value - anchor) / modulus
    return anchor + modulus * (-(-quotient.numerator // quotient.denominator))


def decode_bits(bits, order):
    require(len(bits) == order * (order - 1) // 2, "certificate bit length")
    adjacency = [0] * order
    cursor = 0
    for left in range(order):
        for right in range(left + 1, order):
            if bits[cursor] == "1":
                adjacency[left] |= 1 << right
            elif bits[cursor] == "0":
                adjacency[right] |= 1 << left
            else:
                raise RuntimeError("nonbinary canonical certificate")
            cursor += 1
    return tuple(adjacency)


def decode_code(code, order):
    bits = []
    for cursor in range(order * (order - 1) // 2):
        bits.append("1" if code & (1 << cursor) else "0")
    return decode_bits("".join(bits), order)


def bit_certificate(adjacency, ordering):
    code = 0
    order = len(adjacency)
    for left_index in range(order):
        left = ordering[left_index]
        for right_index in range(left_index + 1, order):
            code <<= 1
            right = ordering[right_index]
            if adjacency[left] & (1 << right):
                code |= 1
    return code


def stable_vertex_colors(adjacency):
    """Invariant directed color refinement, beginning with outdegree."""
    order = len(adjacency)
    initial = tuple(adjacency[vertex].bit_count() for vertex in range(order))
    palette = {value: index for index, value in enumerate(sorted(set(initial)))}
    colors = tuple(palette[value] for value in initial)
    while True:
        color_count = max(colors) + 1
        signatures = []
        for vertex in range(order):
            out_counts = [0] * color_count
            for other in range(order):
                if adjacency[vertex] & (1 << other):
                    out_counts[colors[other]] += 1
            signatures.append((colors[vertex], tuple(out_counts)))
        signature_palette = {
            value: index for index, value in enumerate(sorted(set(signatures)))
        }
        refined = tuple(signature_palette[value] for value in signatures)
        if len(set(refined)) == len(set(colors)):
            return refined
        colors = refined


@lru_cache(maxsize=None)
def canonical_code(adjacency):
    """Minimum certificate among permutations preserving invariant colors."""
    order = len(adjacency)
    colors = stable_vertex_colors(adjacency)
    blocks = tuple(
        tuple(vertex for vertex in range(order) if colors[vertex] == color)
        for color in range(max(colors) + 1)
    )
    pools = tuple(tuple(permutations(block)) for block in blocks)
    best = None
    for choices in product(*pools):
        ordering = tuple(vertex for block in choices for vertex in block)
        candidate = bit_certificate(adjacency, ordering)
        if best is None or candidate < best:
            best = candidate
    require(best is not None, "nonempty canonical-label search")
    return best


def format_code(code, order):
    width = order * (order - 1) // 2
    return format(code, f"0{width}b")


def extend_by_vertex(adjacency, old_to_new_mask):
    old_order = len(adjacency)
    new_bit = 1 << old_order
    child = list(adjacency) + [0]
    for vertex in range(old_order):
        if old_to_new_mask & (1 << vertex):
            child[vertex] |= new_bit
        else:
            child[old_order] |= 1 << vertex
    return tuple(child)


def canonical_augmentation_through_eight():
    """Deletion-closed vertex augmentation with canonical deduplication."""
    representatives = {"": (0,)}
    by_order = {1: representatives}
    for order in range(2, 9):
        children = {}
        for parent in representatives.values():
            for mask in range(1 << (order - 1)):
                child = extend_by_vertex(parent, mask)
                code = canonical_code(child)
                label = format_code(code, order)
                if label not in children:
                    children[label] = decode_bits(label, order)
        representatives = dict(sorted(children.items()))
        by_order[order] = representatives
    for order, expected in EXPECTED_ALL_COUNTS.items():
        require(len(by_order[order]) == expected,
                f"pure canonical augmentation count at order {order}")
    return by_order


def reachable(adjacency, root, reverse=False):
    order = len(adjacency)
    reached = 1 << root
    frontier = reached
    while frontier:
        bit = frontier & -frontier
        frontier ^= bit
        vertex = bit.bit_length() - 1
        if reverse:
            neighborhood = 0
            for other in range(order):
                if adjacency[other] & bit:
                    neighborhood |= 1 << other
        else:
            neighborhood = adjacency[vertex]
        fresh = neighborhood & ~reached
        reached |= fresh
        frontier |= fresh
    return reached


def is_strong(adjacency):
    full = (1 << len(adjacency)) - 1
    return reachable(adjacency, 0) == full and reachable(adjacency, 0, True) == full


def forward_reverse_paths(adjacency):
    """Hamilton-path counts by last and by first vertex on every subset."""
    order = len(adjacency)
    size = 1 << order
    ending = [[0] * order for _ in range(size)]
    starting = [[0] * order for _ in range(size)]
    for vertex in range(order):
        ending[1 << vertex][vertex] = 1
        starting[1 << vertex][vertex] = 1
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        vertices = mask
        while vertices:
            bit = vertices & -vertices
            vertices ^= bit
            vertex = bit.bit_length() - 1
            rest = mask ^ bit
            others = rest
            while others:
                other_bit = others & -others
                others ^= other_bit
                other = other_bit.bit_length() - 1
                if adjacency[other] & bit:
                    ending[mask][vertex] += ending[rest][other]
                if adjacency[vertex] & other_bit:
                    starting[mask][vertex] += starting[rest][other]
    return ending, starting


def contracted_block_exposure(adjacency):
    """Count good and reversed one-defect adjacent blocks for every arc."""
    order = len(adjacency)
    size = 1 << order
    full = size - 1
    ending, starting = forward_reverse_paths(adjacency)
    H = sum(ending[full])
    require(H > 0, "Hamilton path existence")

    before_entry = [[0] * order for _ in range(size)]
    after_exit = [[0] * order for _ in range(size)]
    for vertex in range(order):
        before_entry[0][vertex] = 1
        after_exit[0][vertex] = 1
    for mask in range(1, size):
        for boundary in range(order):
            before_entry[mask][boundary] = sum(
                ending[mask][last]
                for last in range(order)
                if adjacency[last] & (1 << boundary)
            )
            after_exit[mask][boundary] = sum(
                starting[mask][first]
                for first in range(order)
                if adjacency[boundary] & (1 << first)
            )

    def block_count(entry, exit_vertex):
        remaining = full ^ (1 << entry) ^ (1 << exit_vertex)
        answer = 0
        left = remaining
        while True:
            right = remaining ^ left
            answer += before_entry[left][entry] * after_exit[right][exit_vertex]
            if left == 0:
                break
            left = (left - 1) & remaining
        return answer

    capacity = [[0] * order for _ in range(order)]
    for tail in range(order):
        heads = adjacency[tail]
        while heads:
            head_bit = heads & -heads
            heads ^= head_bit
            head = head_bit.bit_length() - 1
            good_block = block_count(tail, head)
            reversed_one_defect_block = block_count(head, tail)
            capacity[tail][head] = good_block + reversed_one_defect_block
    return H, tuple(tuple(row) for row in capacity)


def directed_cut_responses(H, capacity):
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
    ending, unused_starting = forward_reverse_paths(adjacency)
    del unused_starting
    return sum(ending[-1])


def add_ear(adjacency, cut):
    order = len(adjacency)
    child = list(adjacency) + [0]
    ear_bit = 1 << order
    for vertex in range(order):
        if cut & (1 << vertex):
            child[order] |= 1 << vertex
        else:
            child[vertex] |= ear_bit
    return tuple(child)


def literal_responses(adjacency):
    return tuple(
        hamilton_count(add_ear(adjacency, cut))
        for cut in range(1 << len(adjacency))
    )


def analyze(adjacency):
    order = len(adjacency)
    H, capacity = contracted_block_exposure(adjacency)
    values = directed_cut_responses(H, capacity)
    require(values[0] == values[-1] == H, "constant-ear equality")
    require(all(value >= H and value % 2 == 1 for value in values),
            "Redei lower support and parity")

    weights = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    for tail in range(order):
        for head in range(order):
            if capacity[tail][head]:
                require(adjacency[tail] & (1 << head), "capacity follows arc")
                weights[tail][head] = weights[head][tail] = Fraction(
                    capacity[tail][head], 2
                )
    degrees = tuple(sum(row) for row in weights)
    field = tuple(
        sum(
            weights[vertex][other]
            if adjacency[vertex] & (1 << other)
            else -weights[vertex][other]
            for other in range(order) if other != vertex
        )
        for vertex in range(order)
    )
    require(sum(field) == 0, "zero-sum contracted-block field")
    W = sum(weights[left][right] for left, right in combinations(range(order), 2))
    edges = tuple(combinations(range(order), 2))
    D4 = sum(
        weights[a][b] * weights[c][d]
        for index, (a, b) in enumerate(edges)
        for c, d in edges[index + 1:]
        if len({a, b, c, d}) == 4
    )
    C_hd = sum(field[vertex] * degrees[vertex] for vertex in range(order))
    require(W > 0 and D4 > 0, "positive contracted-block exposure masses")
    theta = Fraction((order - 3) * C_hd, 2 * D4)

    layers = []
    for size in range(1, order):
        t = order - 2 * size
        layer_values = tuple(
            values[mask] for mask in range(1 << order) if mask.bit_count() == size
        )
        mean = Fraction(sum(layer_values), len(layer_values))
        variance = sum((Fraction(value) - mean) ** 2 for value in layer_values)
        variance /= len(layer_values)
        require(mean > H, "positive independent layer mean gap")
        rational_floor = mean + variance / (mean - H)
        anchor = layer_values[-1]
        lattice = gcd_all(value - anchor for value in layer_values)
        exact_floor = coset_ceiling(rational_floor, anchor, lattice)
        require(max(layer_values) >= exact_floor, "independent coset support")
        layers.append({
            "size": size,
            "t": t,
            "floor": rational_floor,
            "coset_floor": exact_floor,
            "maximum": max(layer_values),
            "lattice": lattice,
            "values": tuple(sorted(layer_values)),
        })

    grid = tuple(range(-(order - 2), order - 1, 2))
    distance = min(abs(Fraction(t) - theta) for t in grid)
    predicted_rational_t = tuple(t for t in grid if abs(Fraction(t) - theta) == distance)
    best_rational = max(layer["floor"] for layer in layers)
    rational_t = tuple(sorted(layer["t"] for layer in layers
                              if layer["floor"] == best_rational))
    require(rational_t == predicted_rational_t, "independent nearest-grid rule")
    best_coset = max(layer["coset_floor"] for layer in layers)
    coset_t = tuple(sorted(layer["t"] for layer in layers
                          if layer["coset_floor"] == best_coset))
    best_actual = max(layer["maximum"] for layer in layers)
    actual_t = tuple(sorted(layer["t"] for layer in layers
                           if layer["maximum"] == best_actual))
    central = {0} if order % 2 == 0 else {-1, 1}
    rational_central = not central.isdisjoint(rational_t)
    coset_central = not central.isdisjoint(coset_t)
    actual_central = not central.isdisjoint(actual_t)
    normalized_tilt = Fraction(
        (order - 3) * abs(C_hd), (2 if order % 2 == 0 else 4) * D4
    )
    require(rational_central == (normalized_tilt <= 1),
            "independent normalized centrality criterion")
    return {
        "order": order,
        "H": H,
        "W": W,
        "D4": D4,
        "C_hd": C_hd,
        "theta": theta,
        "normalized_tilt": normalized_tilt,
        "rational_t": rational_t,
        "coset_t": coset_t,
        "actual_t": actual_t,
        "rational_central": rational_central,
        "coset_central": coset_central,
        "actual_central": actual_central,
        "layers": tuple(layers),
        "values": values,
    }


def profile_signature(record):
    return (
        record["order"], record["H"], pair(record["W"]), pair(record["D4"]),
        pair(record["C_hd"]), pair(record["theta"]),
        record["rational_t"], record["coset_t"], record["actual_t"],
        tuple(
            (
                layer["size"], layer["t"], pair(layer["floor"]),
                layer["coset_floor"], layer["maximum"], layer["lattice"],
                layer["values"],
            )
            for layer in record["layers"]
        ),
    )


def compact(record):
    return {
        "order": record["order"],
        "H_W_D4_Chd": (
            record["H"], pair(record["W"]), pair(record["D4"]),
            pair(record["C_hd"]),
        ),
        "theta": pair(record["theta"]),
        "normalized_tilt": pair(record["normalized_tilt"]),
        "rational_t": record["rational_t"],
        "coset_t": record["coset_t"],
        "actual_t": record["actual_t"],
        "layers": tuple(
            (
                layer["size"], layer["t"], pair(layer["floor"]),
                layer["coset_floor"], layer["maximum"], layer["lattice"],
            )
            for layer in record["layers"]
        ),
    }


def increment(histogram, key):
    histogram[key] = histogram.get(key, 0) + 1


def summarize_order(order, records):
    require(len(records) == EXPECTED_STRONG_COUNTS[order], "strong class row size")
    rational_histogram = {}
    coset_histogram = {}
    actual_histogram = {}
    worst = max(record["normalized_tilt"] for record in records)
    worst_records = tuple(record for record in records
                          if record["normalized_tilt"] == worst)
    central = {0} if order % 2 == 0 else {-1, 1}
    strict_coset_margins = []
    for record in records:
        increment(rational_histogram, record["rational_t"])
        increment(coset_histogram, record["coset_t"])
        increment(actual_histogram, record["actual_t"])
        best_central = max(
            layer["coset_floor"] for layer in record["layers"]
            if layer["t"] in central
        )
        best_noncentral = max(
            layer["coset_floor"] for layer in record["layers"]
            if layer["t"] not in central
        )
        strict_coset_margins.append(best_central - best_noncentral)
    representative = min(worst_records, key=profile_signature)
    return {
        "order": order,
        "all_classes": EXPECTED_ALL_COUNTS[order],
        "strong_classes": len(records),
        "rational_central_pass_fail": (
            sum(record["rational_central"] for record in records),
            sum(not record["rational_central"] for record in records),
        ),
        "coset_central_pass_fail": (
            sum(record["coset_central"] for record in records),
            sum(not record["coset_central"] for record in records),
        ),
        "actual_central_pass_fail": (
            sum(record["actual_central"] for record in records),
            sum(not record["actual_central"] for record in records),
        ),
        "coset_reorders_rational": sum(
            record["coset_t"] != record["rational_t"] for record in records
        ),
        "rational_optimizer_histogram": tuple(sorted(rational_histogram.items())),
        "coset_optimizer_histogram": tuple(sorted(coset_histogram.items())),
        "actual_optimizer_histogram": tuple(sorted(actual_histogram.items())),
        "worst_normalized_tilt": pair(worst),
        "worst_tilt_multiplicity": len(worst_records),
        "worst_tilt_packet": compact(representative),
        "theta_zero_count": sum(record["theta"] == 0 for record in records),
        "minimum_strict_coset_margin": min(strict_coset_margins),
        "profile_fingerprint": digest(sorted(profile_signature(record)
                                             for record in records)),
    }


def named_controls():
    records = {
        name: analyze(decode_code(code, order))
        for name, code, order in (
            ("code2_n6", 2, 6),
            ("code140_n6", 140, 6),
            ("code20_n6", 20, 6),
        )
    }
    require(records["code2_n6"]["rational_t"] == (2,),
            "code 2 rational noncentral control")
    require(records["code140_n6"]["coset_t"] == (4,),
            "code 140 exact-coset reorder control")
    require(records["code20_n6"]["rational_t"] == (0,)
            and 0 in records["code20_n6"]["coset_t"]
            and 0 not in records["code20_n6"]["actual_t"],
            "code 20 support/actual boundary")
    for name, code in (("code2_n6", 2), ("code140_n6", 140), ("code20_n6", 20)):
        adjacency = decode_code(code, 6)
        require(records[name]["values"] == literal_responses(adjacency),
                f"literal child replay {name}")
    return {name: compact(record) for name, record in records.items()}


def main():
    universes = canonical_augmentation_through_eight()
    rows = []
    canonical_fingerprints = {}
    strong_canonical_fingerprints = {}
    combined_strong_labels = []
    for order in range(4, 9):
        all_records = universes[order]
        canonical_fingerprints[order] = raw_digest(tuple(all_records))
        strong_items = tuple(
            (label, adjacency) for label, adjacency in all_records.items()
            if is_strong(adjacency)
        )
        strong_adjacencies = tuple(adjacency for unused_label, adjacency in strong_items)
        strong_labels = tuple(label for label, unused_adjacency in strong_items)
        strong_canonical_fingerprints[order] = raw_digest(strong_labels)
        combined_strong_labels.extend(f"{order}:{label}" for label in strong_labels)
        require(len(strong_adjacencies) == EXPECTED_STRONG_COUNTS[order],
                f"independent strong filter count at order {order}")
        records = tuple(analyze(adjacency) for adjacency in strong_adjacencies)
        rows.append(summarize_order(order, records))

    require(tuple(row["actual_central_pass_fail"][1] for row in rows)
            == (0, 0, 4, 0, 1702), "actual-central hostile boundary counts")
    require(tuple(row["theta_zero_count"] for row in rows)
            == (1, 6, 11, 79, 162), "zero-tilt class counts")
    require(tuple(row["minimum_strict_coset_margin"] for row in rows)
            == (2, 6, 6, 22, 16), "strict central coset margins")
    require(all(row["rational_central_pass_fail"][1] == 0 for row in rows),
            "independent strong rational centrality")
    require(all(row["coset_central_pass_fail"][1] == 0 for row in rows),
            "independent strong exact-coset centrality")
    controls = named_controls()
    ledger = {
        "theorem": "THM-4131",
        "universe": "one isomorphism class of every strong tournament, orders 4..8",
        "all_class_counts": tuple(EXPECTED_ALL_COUNTS.items()),
        "strong_class_counts": tuple(EXPECTED_STRONG_COUNTS.items()),
        "order_rows": rows,
        "controls": controls,
        "scope": (
            "finite exact rational and exact-coset support centrality; actual response "
            "maxima may be noncentral; no all-order strong theorem or interval"
        ),
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")

    print("status=ACCEPT")
    print("implementation=pure canonical augmentation plus contracted good/reversed-block DP")
    print(f"canonical_universe_sha256={tuple(canonical_fingerprints.items())}")
    print(f"strong_canonical_universe_sha256={tuple(strong_canonical_fingerprints.items())}")
    print(f"combined_strong_canonical_sha256={raw_digest(tuple(combined_strong_labels))}")
    print(f"order_rows={rows}")
    print(f"controls={controls}")
    print("scope=FINITE-EXACT orders 4..8; support floors central, actual maxima need not be")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
