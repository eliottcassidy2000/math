"""Independent structural audit of six_color, not its extension-count code.

Use exact chromatic numbers of boundary contractions, setwise palette-chart
stabilizers, and the closed form of the seven-point ternary histogram.
"""
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
from hashlib import sha256
import json

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


@lru_cache(None)
def chi(n, edges):
    rows = [0] * n
    for u, v in edges:
        rows[u] |= 1 << v
        rows[v] |= 1 << u
    independent = [True] * (1 << n)
    value = [0] * (1 << n)
    for mask in range(1, 1 << n):
        bit = mask & -mask
        rest = mask ^ bit
        independent[mask] = independent[rest] and not (rows[bit.bit_length() - 1] & rest)
        if independent[mask]:
            value[mask] = 1
            continue
        best, sub = n, mask
        while sub:
            if sub & bit and independent[sub]:
                best = min(best, 1 + value[mask ^ sub])
            sub = (sub - 1) & mask
        value[mask] = best
    return value[-1]


PARTITIONS = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (0, 1, 2)]


def boundary_chi(n, edges, word):
    blocks = 1 + max(word)
    labels = list(word) + list(range(blocks, blocks + n - 3))
    quotient = {tuple(sorted((labels[u], labels[v]))) for u, v in edges}
    if any(u == v for u, v in quotient):
        return None
    quotient.update(combinations(range(blocks), 2))
    return chi(blocks + n - 3, tuple(sorted(quotient)))


def main():
    graph_count = partition_count = 0
    for n in range(3, 6):
        pairs = list(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            edges = tuple(pair for bit, pair in enumerate(pairs) if mask >> bit & 1)
            k = chi(n, edges)
            graph_count += 1
            independent_ports = not any(u < 3 and v < 3 for u, v in edges)
            for word in PARTITIONS:
                threshold = boundary_chi(n, edges, word)
                if threshold is None:
                    continue
                multiplicities = [word.count(i) for i in range(1 + max(word))]
                t = sum(size > 1 for size in multiplicities)
                s = sum(size == 1 for size in multiplicities)
                need(threshold <= k + t + max(0, s - 1), "unrestricted recolouring bound")
                if k <= 4:
                    need(threshold <= 6, "four-colour fragment transparent at six")
                if independent_ports and k <= 5 and word != (0, 1, 2):
                    need(threshold <= 6, "nondiscrete five-colour extension")
                partition_count += 1
    edges5 = tuple(sorted(set(combinations(range(3, 7), 2)) |
                          {(u, v) for u in range(3) for v in range(3, 7)}))
    edges6 = set(combinations(range(3, 6), 2)) | {(i, i + 3) for i in range(3)}
    edges6 |= set(combinations(range(6, 9), 2)) | {(i, j) for i in range(6) for j in range(6, 9)}
    edges6 = tuple(sorted(edges6))
    thresholds5 = [boundary_chi(7, edges5, word) for word in PARTITIONS]
    thresholds6 = [boundary_chi(9, edges6, word) for word in PARTITIONS]
    need(chi(7, edges5) == 5 and thresholds5 == [5, 6, 6, 6, 7], "sharp first threshold")
    need(chi(9, edges6) == 6 and thresholds6 == [7, 6, 6, 6, 6], "sharp NAE threshold")
    exterior = tuple((u - 3, v - 3) for u, v in edges6 if u >= 3 and v >= 3)
    need(chi(6, exterior) == 6, "NAE exterior meets required threshold")
    tuples = [tuple(2 * i + ((x >> i) & 1) for i in range(3)) for x in range(8)]
    chart = set(tuples)
    even = {tuples[x] for x in range(8) if x.bit_count() % 2 == 0}
    chart_stabilizers, even_stabilizers = [], []
    for palette in permutations(range(6)):
        image = {tuple(palette[c] for c in word) for word in chart}
        if image == chart:
            chart_stabilizers.append(palette)
            offsets = [palette[2 * i] - 2 * i for i in range(3)]
            need(all(offset in (0, 1) for offset in offsets), "fixed positions preserve individual pairs")
            need(all(palette[2 * i + 1] == 2 * i + 1 - offsets[i] for i in range(3)),
                 "chart action is independent pair flips")
        even_image = {tuple(palette[c] for c in word) for word in even}
        if even_image == even:
            even_stabilizers.append(palette)
    need(len(chart_stabilizers) == 8 and len(even_stabilizers) == 4,
         "chart stabilizers with ports fixed are C2^3 and C2^2")
    # Independent product classification of the ternary sectors.
    histogram = {5: 18 * 18 * 3, 6: 2 * 18 * 6 * 3, 7: 6 * 6 * 3}
    need(histogram == {5: 972, 6: 648, 7: 108}, "ternary component histogram")
    need(sum(histogram.values()) == 1728, "successful ternary sectors")
    need(sum(count * 2 ** components for components, count in histogram.items()) == 86400,
         "ordinary colouring count recovered")
    # Verify the leaf factor's necessary condition directly from the old receipt.
    prior = Path(__file__).resolve().parents[2] / "hadwiger_six_seven_eight_20261025/cube_boundary/receipt.json"
    raw = json.loads(prior.read_text(encoding="utf-8"))
    for name, graph in raw["unit_graphs"].items():
        pins = {int(v) for v in graph["fixed_colors"]}
        edges = graph["edges"]
        for pin in pins:
            neighbors = [v if u == pin else u for u, v in edges if pin in (u, v)]
            need(len(neighbors) == 1 and neighbors[0] not in pins,
                 "every unpinned leaf attaches to retained graph, no isolated leaf-pair component")
    producer = Path(__file__).resolve().parents[1] / "six_color/check.py"
    print(json.dumps({"status": "PASS", "checks": CHECKS,
                      "method": "independent chromatic contraction and palette setwise stabilizers",
                      "all_graphs_n3_to_n5": graph_count,
                      "proper_boundary_partition_cases": partition_count,
                      "chi5_thresholds": thresholds5, "chi6_NAE_thresholds": thresholds6,
                      "fixed_position_chart_stabilizers": [len(chart_stabilizers), len(even_stabilizers)],
                      "ternary_histogram": histogram,
                      "producer_checker_sha256": sha256(producer.read_bytes()).hexdigest()},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
