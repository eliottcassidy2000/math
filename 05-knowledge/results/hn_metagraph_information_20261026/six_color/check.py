"""Exact unprescribed six-colour boundary and paired-frame controls.

Run: python check.py > receipt.json
     python -O check.py > optimized.json
No third-party packages; all checks remain active under -O.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

CHECKS = 0


def require(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)


def canon(word):
    names = {}
    return tuple(names.setdefault(a, len(names)) for a in word)


PARTITIONS = sorted({canon(w) for w in product(range(3), repeat=3)})


def count_extensions(n, edges, ports, word, palette=6):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    colors = [-1] * n
    for v, a in zip(ports, word):
        colors[v] = a
    if any(colors[u] >= 0 and colors[u] == colors[v] for u, v in edges):
        return 0

    def rec():
        remaining = [v for v in range(n) if colors[v] < 0]
        if not remaining:
            return 1
        v = max(remaining, key=lambda u: (sum(colors[z] >= 0 for z in adj[u]), len(adj[u])))
        forbidden = {colors[u] for u in adj[v]}
        total = 0
        for a in range(palette):
            if a not in forbidden:
                colors[v] = a
                total += rec()
        colors[v] = -1
        return total

    return rec()


def norm(p, q):
    a, b = p[0] - q[0], p[1] - q[1]
    return a*a + a*b + b*b


def unit_edges(points):
    require(len(set(points)) == len(points), "distinct exact points")
    return [(i, j) for i, j in combinations(range(len(points)), 2) if norm(points[i], points[j]) == 1]


def bipartite_components(n, edges):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    side = {}
    components = 0
    for start in range(n):
        if start in side:
            continue
        components += 1
        side[start] = 0
        stack = [start]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v in side:
                    if side[v] == side[u]:
                        return None
                else:
                    side[v] = 1 - side[u]
                    stack.append(v)
    return components


def recolor_control():
    # Every 3-port equality partition and every initial k-colouring of the
    # three ports. Build the densest graph consistent with that colouring
    # and with the required independent boundary blocks. Add one interior
    # vertex in each old colour class, so all permitted old-colour conflicts
    # are exposed. This is a constructive-test universe, not all graphs.
    count = 0
    for k in range(1, 7):
        for old_boundary in product(range(k), repeat=3):
            old = list(old_boundary) + list(range(k))
            for partition in PARTITIONS:
                blocks = [[i for i, label in enumerate(partition) if label == b]
                          for b in sorted(set(partition))]
                edges = [(u, v) for u, v in combinations(range(k+3), 2)
                         if old[u] != old[v]
                         and not (u < 3 and v < 3 and partition[u] == partition[v])]
                new = list(old)
                fresh = k
                # Recolour every nonsingleton block, then keep the first
                # singleton in each old colour and separate its repetitions.
                singles = []
                for block in blocks:
                    if len(block) > 1:
                        for v in block:
                            new[v] = fresh
                        fresh += 1
                    else:
                        singles.append(block[0])
                seen = set()
                for v in singles:
                    if old[v] in seen:
                        new[v] = fresh
                        fresh += 1
                    seen.add(old[v])
                t = sum(len(b) > 1 for b in blocks)
                s = len(singles)
                bound = k + t + max(0, s-1)
                require(canon(new[:3]) == partition, "requested boundary equality partition")
                require(all(new[u] != new[v] for u, v in edges), "recolouring remains proper")
                require(max(new) < bound, "fresh-colour bound")
                count += 1
    return count


def frame_control():
    cube = list(range(8))
    even = [0, 3, 5, 6]
    tuples = [tuple(2*j + ((x >> j) & 1) for j in range(3)) for x in cube]
    require(len({canon(w) for w in tuples}) == 1, "all endpoint triples have one S6 orbit")
    # All 720 palette permutations, not an assumed transitivity statement.
    images = {tuple(p[a] for a in tuples[0]) for p in permutations(range(6))}
    require(len(images) == 120 and all(w in images for w in tuples), "injective-triple orbit has size120")
    diff_partitions = {}
    for x, y in product(cube, repeat=2):
        key = canon(tuples[x] + tuples[y])
        d = x ^ y
        require(key not in diff_partitions or diff_partitions[key] == d, "orbit determines difference")
        diff_partitions[key] = d
    require(len(diff_partitions) == 8, "eight two-frame equality partitions")
    commuting = []
    for perm in permutations(even):
        f = dict(zip(even, perm))
        if all(f[x ^ b] == (f[x] ^ b) for x, b in product(even, repeat=2)):
            commuting.append(perm)
            require(all(f[x] == (x ^ f[0]) for x in even), "centralizer consists of translations")
    require(len(commuting) == 4, "four of24 deterministic S4 tables survive")
    relations = 0
    for mask in range(1 << 16):
        relation = {(even[i], even[j]) for i in range(4) for j in range(4)
                    if mask >> (4*i+j) & 1}
        if all(((x ^ b, y ^ b) in relation) for x, y in relation for b in even):
            differences = {x ^ y for x, y in relation}
            require(relation == {(x, y) for x, y in product(even, repeat=2) if x ^ y in differences},
                    "all invariant relations are difference kernels")
            relations += 1
    require(relations == 16, "16 invariant completion relations")
    return {"palette_permutations": 720, "injective_triple_orbit": 120,
            "two_frame_orbits": 8, "four_state_relations_tested": 65536,
            "translation_invariant_relations": relations,
            "S4_permutations_tested": 24, "equivariant_permutations": len(commuting)}


def sharp_abstract_controls():
    # K4 joined to three independent ports: chi5, omits only all-distinct.
    edges5 = list(combinations(range(3, 7), 2)) + [(p, v) for p in range(3) for v in range(3, 7)]
    weights5 = [count_extensions(7, edges5, (0, 1, 2), w) for w in PARTITIONS]
    require(weights5 == [120, 24, 24, 24, 0], "sharp chi5 three-port control")
    # Join K3 to the inherited triangle-with-three-pendants NAE gadget.
    edges6 = list(combinations(range(3, 6), 2)) + [(i, i+3) for i in range(3)]
    edges6 += list(combinations(range(6, 9), 2)) + [(a, b) for a in range(6) for b in range(6, 9)]
    weights6 = [count_extensions(9, edges6, (0, 1, 2), w) for w in PARTITIONS]
    require(weights6 == [0, 48, 48, 48, 12], "sharp chi6 NAE control")
    # Four equidistant points give a rank3 Gram matrix based at one point.
    gram = [[Q(1) if i == j else Q(1, 2) for j in range(3)] for i in range(3)]
    determinant = sum((-1 if sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)) % 2 else 1)
                      * gram[0][p[0]] * gram[1][p[1]] * gram[2][p[2]]
                      for p in permutations(range(3)))
    require(determinant == Q(1, 2), "K4 metric obstruction")
    return {"chi5_weights": weights5, "chi6_NAE_weights": weights6,
            "K4_Gram_determinant": str(determinant), "unit_realizable": False}


def metric_holonomy_control():
    points = [(Q(0), Q(0)), (Q(1), Q(0)), (Q(0), Q(1)),
              (Q(10), Q(0)), (Q(11), Q(0)), (Q(10), Q(1)), (Q(100), Q(0))]
    edges = unit_edges(points)
    require(edges == list(combinations(range(3), 2)) + list(combinations(range(3, 6), 2)),
            "complete unit graph is two disjoint triangles plus one isolate")
    coloring = [0, 1, 2, 0, 1, 2, 0]
    require(all(coloring[u] != coloring[v] for u, v in edges), "explicit genuine three-colouring")
    # x_j=(bit at A_j, bit at B_j, bit at common Z); gain flips first two.
    gain = 3
    require(gain.bit_count() % 2 == 0 and gain ^ gain ^ gain == gain, "even nonzero cycle gain")
    solutions = [xs for xs in product(range(8), repeat=3)
                 if all(xs[(j+1) % 3] == (xs[j] ^ gain) for j in range(3))]
    require(not solutions, "conditional frame CSP is inconsistent")
    q_success = 0
    total_colorings = 0
    component_histogram = Counter()
    for q in product(range(3), repeat=7):
        same = [(u, v) for u, v in edges if q[u] == q[v]]
        components = bipartite_components(7, same)
        if components is not None:
            q_success += 1
            component_histogram[components] += 1
            total_colorings += 2**components
    require(q_success == 1728, "all ternary pair maps: 1728 successful")
    require(total_colorings == 86400, "sum over q equals genuine six-colouring count")
    require(total_colorings == count_extensions(7, edges, (), ()), "independent ordinary colour backtracking")
    weights = [count_extensions(7, edges, (0, 3, 6), w) for w in PARTITIONS]
    require(weights == [400]*5, "full unprescribed one-frame boundary weights")
    return {"points_in_1_alpha_basis": [[str(a), str(b)] for a, b in points], "edges": edges,
            "three_coloring": coloring, "gain_bits": [1, 1, 0], "frame_CSP_solutions": len(solutions),
            "ternary_maps_tested": 3**7, "ternary_maps_successful": q_success,
            "component_histogram": dict(sorted(component_histogram.items())),
            "six_colorings": total_colorings, "three_port_weights": weights}


def unpin_previous_unit_gadgets():
    here = Path(__file__).resolve()
    prior = here.parents[2] / "hadwiger_six_seven_eight_20261025" / "cube_boundary" / "receipt.json"
    data = json.loads(prior.read_text(encoding="utf-8"))
    output = {"input_receipt_sha256": sha256(prior.read_bytes()).hexdigest()}
    base_points = [(Q(0), Q(0)), (Q(1), Q(0)), (Q(0), Q(1)),
                   (Q(0), Q(-1)), (Q(2), Q(0)), (Q(0), Q(2))]
    base_edges = unit_edges(base_points)
    base_weights = [count_extensions(6, base_edges, (3, 4, 5), w) for w in PARTITIONS]
    require(base_weights == [60, 68, 68, 68, 71], "unprescribed six-colour base weights")
    for name, packet in data["unit_graphs"].items():
        pts = [tuple(Q(a) for a in point) for point in packet["points"]]
        edges = unit_edges(pts)
        require(edges == [tuple(edge) for edge in packet["edges"]], "all old metric edges independently recomputed")
        degrees = Counter(v for edge in edges for v in edge)
        pins = {int(v) for v in packet["fixed_colors"]}
        require(all(degrees[v] == 1 for v in pins), "old fixed vertices really are leaves")
        remainder = [v for v in range(len(pts)) if v not in pins]
        rem_edges = [(u, v) for u, v in edges if u not in pins and v not in pins]
        if name == "F":
            require(len(pins) == 25 and len(remainder) == 7 and rem_edges == base_edges,
                    "F is six-vertex base plus isolated centre after deleting leaves")
            weights = [6 * 5**25 * w for w in base_weights]
        else:
            require(name == "G" and len(pins) == 15 and len(remainder) == 4 and not rem_edges,
                    "G consists of three port stars and a separate star")
            weights = [6 * 5**15]*5
        require(all(w > 0 for w in weights), "every unprescribed boundary orbit survives")
        output[name] = {"vertices": len(pts), "edges": len(edges), "unprescribed_leaves": len(pins),
                        "weights_by_partition": weights}
    output["base_weights"] = base_weights
    return output


def main():
    output = {"status": "FINITE-EXACT", "partitions": PARTITIONS,
              "constructive_recolor_cases": recolor_control(), "frame_symmetry": frame_control(),
              "abstract_sharp_controls": sharp_abstract_controls(),
              "metric_holonomy_hostile": metric_holonomy_control(),
              "unpin_previous_gadgets": unpin_previous_unit_gadgets()}
    output["checks"] = CHECKS
    output["source_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(output, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
