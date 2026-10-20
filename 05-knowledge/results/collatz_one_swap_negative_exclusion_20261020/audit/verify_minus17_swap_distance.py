#!/usr/bin/env python3
"""Independent exact cyclic-swap distance audit, using bit-mask graph states."""

from collections import Counter, deque
from hashlib import sha256
from math import comb
from pathlib import Path
import json

CHECKS = 0


def require(condition, detail):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(detail)


def bits_to_mask(bits):
    return sum(bit << i for i, bit in enumerate(bits))


def mask_to_bits(mask, q):
    return [(mask >> i) & 1 for i in range(q)]


def word(mask, q):
    return "".join(map(str, mask_to_bits(mask, q)))


def carry(bits):
    output, power3 = 0, 1
    for i in range(len(bits) - 1, -1, -1):
        if bits[i]:
            output += (1 << i) * power3
            power3 *= 3
    return output


# Derive the target cycle rather than accepting its supplied parity word.
state, orbit, parity = -17, [], []
while state not in orbit:
    orbit.append(state)
    parity.append(state % 2)
    state = (3 * state + 1) // 2 if state % 2 else state // 2
require(state == -17 and len(orbit) == 11, "-17 orbit closure")
p, q = sum(parity), len(parity)
source = [((i + 1) * p) // q - (i * p) // q for i in range(q)]
start = bits_to_mask(source)

# The graph includes every edge, including q-1 -> 0, in both directions.
distance, parent, queue = {start: 0}, {}, deque([start])
edge_occurrences = 0
while queue:
    current = queue.popleft()
    for edge in range(q):
        other = (edge + 1) % q
        if ((current >> edge) & 1) == ((current >> other) & 1):
            continue
        changed = current ^ (1 << edge) ^ (1 << other)
        edge_occurrences += 1
        require(changed.bit_count() == p, "weight-preserving edge")
        if changed not in distance:
            distance[changed] = distance[current] + 1
            parent[changed] = (current, edge)
            queue.append(changed)
require(len(distance) == comb(q, p) == 330, "full fixed-weight graph reached")

targets = []
for rotation in range(q):
    target = parity[rotation:] + parity[:rotation]
    target_mask = bits_to_mask(target)
    prefix, imbalances = 0, []
    for a, b in zip(source, target):
        prefix += a - b
        imbalances.append(prefix)
    require(prefix == 0, "same weight for flow bound")
    median = sorted(imbalances)[q // 2]
    flow_lower = sum(abs(value - median) for value in imbalances)
    require(distance[target_mask] == flow_lower, ("flow bound matches BFS", rotation))
    targets.append({
        "rotation": rotation,
        "word": "".join(map(str, target)),
        "bfs_distance": distance[target_mask],
        "flow_lower_bound": flow_lower,
        "prefix_imbalances": imbalances,
        "minimizing_circulation": -median,
    })

minimum = min(row["bfs_distance"] for row in targets)
require(minimum == 3, "exact necklace distance is three")
target = min(targets, key=lambda row: row["bfs_distance"])
target_mask = bits_to_mask(list(map(int, target["word"])))
cursor, reverse_edges, reverse_path = target_mask, [], [target_mask]
while cursor != start:
    cursor, edge = parent[cursor]
    reverse_edges.append(edge)
    reverse_path.append(cursor)
path = [word(mask, q) for mask in reversed(reverse_path)]
path_edges = list(reversed(reverse_edges))
require(len(path_edges) == 3, "independent shortest path length")

# Separately validate the supplied constructive witness.
claimed_edges = [0, 4, 3]
claimed_words = ["".join(map(str, source))]
walk = list(source)
for edge in claimed_edges:
    other = (edge + 1) % q
    require(walk[edge] != walk[other], ("witness uses legal swap", edge))
    walk[edge], walk[other] = walk[other], walk[edge]
    claimed_words.append("".join(map(str, walk)))
require(claimed_words == [
    "01011011011", "10011011011", "10010111011", "10001111011"
], "supplied path reproduced")
require(walk == parity[7:] + parity[:7], "supplied path reaches rotation seven")
clock, b = 2**q - 3**p, carry(walk)
require(clock == -139 and b == 12649 and b // clock == -91,
        "target carry and negative fixed point")
x = b // clock
for bit in walk:
    require(x % 2 == bit, "target parity consistency")
    x = (3 * x + 1) // 2 if bit else x // 2
require(x == -91, "target cycle closure")

layers = dict(sorted(Counter(distance.values()).items()))
require(layers == {0: 1, 1: 8, 2: 29, 3: 60, 4: 82, 5: 90, 6: 43, 7: 14, 8: 3},
        "distance histogram")
require([row["bfs_distance"] for row in targets] == [6, 5, 5, 6, 4, 6, 5, 3, 5, 4, 4],
        "complete rotation distance list")
print(json.dumps({
    "status": "PASS: exact cyclic adjacent-swap necklace distance is three",
    "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "explicit_checks": CHECKS,
    "p": p, "q": q,
    "full_graph_vertices": len(distance),
    "directed_edge_occurrences": edge_occurrences,
    "distance_layers": layers,
    "source": "".join(map(str, source)),
    "target_cycle_word": "".join(map(str, parity)),
    "rotations": targets,
    "minimum_necklace_distance": minimum,
    "independent_shortest_path": path,
    "independent_shortest_path_swap_starts": path_edges,
    "supplied_path": claimed_words,
    "supplied_path_swap_starts": claimed_edges,
    "target_clock": clock, "target_carry": b, "target_fixed_point": b // clock,
}, indent=2, sort_keys=True))
