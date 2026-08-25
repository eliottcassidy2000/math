#!/usr/bin/env python3
"""Nauty/Start-End exact audit for THM-4131.

The theorem universe is the isomorphism-class stream produced by
``gentourng -q -c n`` for 4 <= n <= 8.  Ear responses are reconstructed by
a Start/End/exposed-gap subset DP.  Nothing in this file imports either the
THM-4128 audit or the independent THM-4131 implementation.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from json import dumps
from math import gcd
import os
import shutil
import subprocess


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


def file_digest(path):
    hasher = sha256()
    with open(path, "rb") as source:
        while True:
            block = source.read(1 << 20)
            if not block:
                break
            hasher.update(block)
    return hasher.hexdigest()


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
    require(len(bits) == order * (order - 1) // 2, "triangle bit length")
    adjacency = [0] * order
    cursor = 0
    for left in range(order):
        for right in range(left + 1, order):
            if bits[cursor] == "1":
                adjacency[left] |= 1 << right
            elif bits[cursor] == "0":
                adjacency[right] |= 1 << left
            else:
                raise RuntimeError("nonbinary gentourng record")
            cursor += 1
    return tuple(adjacency)


def decode_code(code, order):
    bits = []
    for cursor in range(order * (order - 1) // 2):
        bits.append("1" if code & (1 << cursor) else "0")
    return decode_bits("".join(bits), order)


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


def find_gentourng():
    configured = os.environ.get("THM4131_GENTOURNG")
    candidates = (
        configured,
        shutil.which("gentourng"),
        "/opt/homebrew/bin/gentourng",
        "/usr/local/bin/gentourng",
    )
    for candidate in candidates:
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    raise RuntimeError("gentourng executable not found")


def nauty_strong_classes(executable, order):
    process = subprocess.run(
        [executable, "-q", "-c", str(order)],
        check=False,
        capture_output=True,
        text=True,
    )
    require(process.returncode == 0, f"gentourng exit at order {order}")
    require(not process.stderr.strip(), f"gentourng stderr at order {order}")
    records = tuple(line.strip() for line in process.stdout.splitlines() if line.strip())
    require(len(records) == EXPECTED_STRONG_COUNTS[order],
            f"gentourng strong count at order {order}")
    require(len(set(records)) == len(records), f"gentourng duplicate at order {order}")
    return records


def start_end_exposed(adjacency):
    """Return full-path starts, ends, and ordered exposed-gap counts."""
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
        todo = mask
        while todo:
            bit = todo & -todo
            todo ^= bit
            vertex = bit.bit_length() - 1
            rest = mask ^ bit
            candidates = rest
            while candidates:
                other_bit = candidates & -candidates
                candidates ^= other_bit
                other = other_bit.bit_length() - 1
                if adjacency[vertex] & other_bit:
                    starts[mask][vertex] += starts[rest][other]
                if adjacency[other] & bit:
                    ends[mask][vertex] += ends[rest][other]

    exposed = [[0] * order for _ in range(order)]
    for left_mask in range(1, full):
        right_mask = full ^ left_mask
        left_vertices = left_mask
        while left_vertices:
            left_bit = left_vertices & -left_vertices
            left_vertices ^= left_bit
            left = left_bit.bit_length() - 1
            left_count = ends[left_mask][left]
            if not left_count:
                continue
            right_vertices = right_mask
            while right_vertices:
                right_bit = right_vertices & -right_vertices
                right_vertices ^= right_bit
                right = right_bit.bit_length() - 1
                exposed[left][right] += left_count * starts[right_mask][right]
    return tuple(starts[full]), tuple(ends[full]), tuple(tuple(row) for row in exposed)


def response_vector(adjacency):
    starts, ends, exposed = start_end_exposed(adjacency)
    order = len(adjacency)
    values = []
    for mask in range(1 << order):
        value = 0
        for vertex in range(order):
            value += starts[vertex] if mask & (1 << vertex) else ends[vertex]
        for left in range(order):
            if mask & (1 << left):
                continue
            for right in range(order):
                if mask & (1 << right):
                    value += exposed[left][right]
        values.append(value)
    return tuple(values)


def exposure_packet(adjacency, values):
    order = len(adjacency)
    H = values[0]
    weights = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    for left, right in combinations(range(order), 2):
        curvature = (
            values[1 << left] + values[1 << right] - H
            - values[(1 << left) | (1 << right)]
        )
        require(curvature >= 0 and curvature % 2 == 0, "exposure curvature")
        weights[left][right] = weights[right][left] = Fraction(curvature, 2)
    degrees = tuple(sum(row) for row in weights)
    field = tuple(values[1 << vertex] - H - degrees[vertex]
                  for vertex in range(order))
    require(sum(field) == 0, "zero-sum exposure field")
    for vertex in range(order):
        divergence = sum(
            weights[vertex][other]
            if adjacency[vertex] & (1 << other)
            else -weights[vertex][other]
            for other in range(order) if other != vertex
        )
        require(divergence == field[vertex], "directed exposure divergence")
    W = sum(weights[left][right] for left, right in combinations(range(order), 2))
    edges = tuple(combinations(range(order), 2))
    D4 = sum(
        weights[a][b] * weights[c][d]
        for index, (a, b) in enumerate(edges)
        for c, d in edges[index + 1:]
        if len({a, b, c, d}) == 4
    )
    C_hd = sum(field[vertex] * degrees[vertex] for vertex in range(order))
    require(W > 0 and D4 > 0, "positive exposure masses")
    return H, W, D4, C_hd


def analyze(adjacency):
    order = len(adjacency)
    values = response_vector(adjacency)
    H, W, D4, C_hd = exposure_packet(adjacency, values)
    require(values[-1] == H, "constant-ear equality")
    require(all(value >= H and value % 2 == 1 for value in values),
            "Redei lower support and parity")
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
        require(mean > H, "positive nonconstant-layer mean gap")
        rational_floor = mean + variance / (mean - H)
        anchor = layer_values[0]
        lattice = gcd_all(value - anchor for value in layer_values)
        exact_floor = coset_ceiling(rational_floor, anchor, lattice)
        require(max(layer_values) >= exact_floor, "exact-coset support floor")
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
    require(rational_t == predicted_rational_t, "THM-4128 nearest-grid rule")
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
            "THM-4128 normalized centrality criterion")
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


def summarize_order(order, records, labels):
    require(len(records) == EXPECTED_STRONG_COUNTS[order], "strong class row size")
    rational_histogram = {}
    coset_histogram = {}
    actual_histogram = {}
    worst = max(record["normalized_tilt"] for record in records)
    worst_indices = tuple(index for index, record in enumerate(records)
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
    representative_index = min(worst_indices, key=lambda index: profile_signature(records[index]))
    representative = compact(records[representative_index])
    row = {
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
        "worst_tilt_multiplicity": len(worst_indices),
        "worst_tilt_packet": representative,
        "theta_zero_count": sum(record["theta"] == 0 for record in records),
        "minimum_strict_coset_margin": min(strict_coset_margins),
        "profile_fingerprint": digest(sorted(profile_signature(record)
                                             for record in records)),
    }
    return row


def named_controls():
    controls = {
        name: compact(analyze(decode_code(code, order)))
        for name, code, order in (
            ("code2_n6", 2, 6),
            ("code140_n6", 140, 6),
            ("code20_n6", 20, 6),
        )
    }
    require(controls["code2_n6"]["rational_t"] == (2,),
            "code 2 rational noncentral control")
    require(controls["code140_n6"]["coset_t"] == (4,),
            "code 140 exact-coset reorder control")
    require(controls["code20_n6"]["rational_t"] == (0,)
            and 0 in controls["code20_n6"]["coset_t"]
            and 0 not in controls["code20_n6"]["actual_t"],
            "code 20 support/actual boundary")
    return controls


def main():
    executable = find_gentourng()
    rows = []
    nauty_fingerprints = {}
    for order in range(4, 9):
        labels = nauty_strong_classes(executable, order)
        nauty_fingerprints[order] = raw_digest(labels)
        records = []
        for label in labels:
            adjacency = decode_bits(label, order)
            require(is_strong(adjacency), f"gentourng -c strong record {order}:{label}")
            records.append(analyze(adjacency))
        rows.append(summarize_order(order, tuple(records), labels))

    require(tuple(row["actual_central_pass_fail"][1] for row in rows)
            == (0, 0, 4, 0, 1702), "actual-central hostile boundary counts")
    require(tuple(row["theta_zero_count"] for row in rows)
            == (1, 6, 11, 79, 162), "zero-tilt class counts")
    require(tuple(row["minimum_strict_coset_margin"] for row in rows)
            == (2, 6, 6, 22, 16), "strict central coset margins")
    require(all(row["rational_central_pass_fail"][1] == 0 for row in rows),
            "strong rational centrality through order eight")
    require(all(row["coset_central_pass_fail"][1] == 0 for row in rows),
            "strong exact-coset centrality through order eight")
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

    print("status=PASS")
    print("implementation=nauty gentourng -q -c plus Start/End/exposed-gap DP")
    print("nauty_version_contract=audited with Homebrew nauty 2.9.3; THM4131_GENTOURNG overrides executable")
    print(f"gentourng_binary_sha256={file_digest(executable)}")
    print(f"nauty_strong_stream_sha256={tuple(nauty_fingerprints.items())}")
    print(f"order_rows={rows}")
    print(f"controls={controls}")
    print("scope=FINITE-EXACT orders 4..8; support floors central, actual maxima need not be")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
