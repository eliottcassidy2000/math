"""Exact controls for the trivalent blue-periphery black-hub theorem.

Run from repository root with python or python -O.  No removable assertions.
The all-order proof is in blue.md; finite controls test the construction.
"""
import csv
from collections import Counter, defaultdict, deque
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OLD_PATH = ROOT / "05-knowledge/results/metagraph_growth_20261020/growth_check.py"
SPEC = spec_from_file_location("inherited_growth", OLD_PATH)
G = module_from_spec(SPEC)
SPEC.loader.exec_module(G)
arc, need = G.arc, G.need
C5_PATHS = [(0, 1, 2, 3, 4), (0, 2, 3, 4, 1), (0, 2, 4, 1, 3)]


def cyclic(n):
    return tuple(sum(1 << v for v in range(n)
                     if 1 <= (v - u) % n <= n // 2) for u in range(n))


def ordered_join(parts):
    n = sum(map(len, parts))
    rows, offset = [], 0
    for a in parts:
        later = sum(1 << v for v in range(offset + len(a), n))
        rows.extend((row << offset) | later for row in a)
        offset += len(a)
    return tuple(rows)


def is_path(a, p):
    return sorted(p) == list(range(len(a))) and all(
        arc(a, u, v) for u, v in zip(p, p[1:]))


def blue(a, p):
    b = G.induced(a, p)
    return b == G.transpose(b)


def flip(a, u, v):
    b = list(a)
    b[u] ^= 1 << v
    b[v] ^= 1 << u
    return tuple(b)


def complement_along(a, p):
    need(is_path(a, p), "complement needs an actual Hamiltonian path")
    path_edges = {tuple(sorted(e)) for e in zip(p, p[1:])}
    b = a
    for u, v in G.pairs(len(a)):
        if (u, v) not in path_edges:
            b = flip(b, u, v)
    return b


def reversal(p):
    result = [None] * len(p)
    for u, v in zip(p, reversed(p)):
        result[u] = v
    return result


def hamilton_paths(a):
    n = len(a)
    answer = []

    def visit(p, unused):
        if not unused:
            answer.append(p)
            return
        for v in range(n):
            if unused >> v & 1 and (not p or arc(a, p[-1], v)):
                visit(p + (v,), unused ^ (1 << v))

    visit((), (1 << n) - 1)
    return answer


def strong(a):
    n = len(a)

    def reach(reverse):
        seen, todo = {0}, [0]
        while todo:
            u = todo.pop()
            for v in range(n):
                if v not in seen and (arc(a, v, u) if reverse else arc(a, u, v)):
                    seen.add(v)
                    todo.append(v)
        return len(seen) == n

    return reach(False) and reach(True)


def blue_tilings(n):
    pairs = G.free_pairs(n)
    index = {e: i for i, e in enumerate(pairs)}
    seen, orbits = set(), []
    for e in pairs:
        if e in seen:
            continue
        reflected = (n - 1 - e[1], n - 1 - e[0])
        orbit = {e, reflected}
        seen.update(orbit)
        orbits.append(sum(1 << index[x] for x in orbit))
    for bits in range(1 << len(orbits)):
        code = sum(mask for i, mask in enumerate(orbits) if bits >> i & 1)
        yield G.tournament(n, code, fixed_path=True)


def nonstrong_blue_sources(n):
    """Unique first strong component, arbitrary symmetric middle, dual sink."""
    for size in range(1, n // 2 + 1):
        for code in range(1 << len(G.free_pairs(size))):
            first = G.tournament(size, code, fixed_path=True)
            if not strong(first):
                continue
            for middle in blue_tilings(n - 2 * size):
                yield ordered_join([first, middle, G.transpose(first)])


def words(weight):
    if weight == 0:
        yield ()
    for size in (1, 3):
        if size <= weight:
            for rest in words(weight - size):
                yield (size,) + rest


def check_word(word):
    a = ordered_join([(0,) if size == 1 else cyclic(3) for size in word])
    r = len(a)
    u = ordered_join([a, cyclic(5), G.transpose(a)])
    n = len(u)
    hub = flip(flip(u, 0, 1), 0, n - 1)
    center_profiles = []
    for middle in C5_PATHS:
        p = tuple(range(r)) + tuple(r + v for v in middle) + tuple(range(r + 5, n))
        need(is_path(u, p) and blue(u, p), "pure-blue source control")
        s = complement_along(u, p)
        need(blue(s, p), "tau preserves blue")
        need(arc(s, p[-1], p[0]), "partner Hamiltonian cycle")
        q = p[1:] + p[:1]
        need(is_path(s, q), "cycle rotation is a path")
        need(not blue(s, q), "rotated partner must be black")
        need(complement_along(s, q) == hub, "common hub in original labels")
        first, second = reversal(p), reversal(q)
        need(all(second[first[p[i]]] == p[(i + 2) % n] for i in range(n)),
             "two positional reversals compose to shift two")
        degrees = [row.bit_count() for row in s]
        need(set(degrees[r:r + 5]) == {r + 2}, "central score shell")
        need(max(degrees[:r]) <= r and min(degrees[r + 5:]) >= r + 4,
             "outer score separation")
        core = G.induced(s, tuple(range(r, r + 5)))
        center_profiles.append(G.hp_subset_dp(core))
    need(center_profiles == [11, 15, 13], "three distinct intrinsic center types")


def graph_components(vertices, edges):
    adjacency = defaultdict(set)
    for u, v, count in edges:
        if count:
            adjacency[u].add(v)
            adjacency[v].add(u)
    unseen, result = set(vertices), []
    while unseen:
        start = next(iter(unseen))
        component, todo = {start}, [start]
        unseen.remove(start)
        while todo:
            for v in adjacency[todo.pop()]:
                if v in unseen:
                    unseen.remove(v)
                    component.add(v)
                    todo.append(v)
        result.append(component)
    return result


def load_census(n):
    stem = ROOT / "04-computation/experiments"
    with (stem / f"tiling_line_metagraph_census_20261020_classes_n{n}.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    node = {int(r["class"]): min(int(r["class"]), int(r["sigmaClass"])) for r in rows}
    data = {}
    for r in rows:
        v = node[int(r["class"])]
        if v not in data:
            scores = list(map(int, r["scores"].split(",")))
            triangles = n * (n - 1) * (n - 2) // 6 - sum(d * (d - 1) // 2 for d in scores)
            data[v] = {"mass": 0, "blue": 0, "type": r["type"],
                       "H": int(r["H"]), "triangles": triangles}
        data[v]["mass"] += int(r["tc"])
        data[v]["blue"] += int(r["bm"])
    counts = {"blue": Counter(), "black": Counter()}
    with (stem / f"tiling_line_metagraph_census_20261020_lines_n{n}.tsv").open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pair = tuple(sorted((node[int(row["classA"])], node[int(row["classB"])])))
            for colour in counts:
                counts[colour][pair] += int(row[colour])
    return data, counts


def degree(counter):
    out = Counter()
    for (u, v), count in counter.items():
        out[u] += count
        out[v] += count
    return out


def bridges(vertices, counters, support=False):
    totals = sum(counters, Counter())
    edges = [(u, v, count) for (u, v), count in totals.items() if count]
    old = len(graph_components(vertices, edges))
    answer = set()
    for u, v, count in edges:
        if u == v or (not support and count != 1):
            continue
        remaining = [(x, y, k) for x, y, k in edges if (x, y) != (u, v)]
        if len(graph_components(vertices, remaining)) > old:
            answer.add((u, v))
    return answer


def main():
    print("PROVED construction / FINITE-EXACT controls: trivalent blue edges have black detours.")
    c5 = cyclic(5)
    # Actual normalized paths in the three partner presentations.
    qs = [(3, 1, 4, 2, 0), (1, 2, 3, 4, 0), (2, 4, 1, 3, 0)]
    target_orders = [(0, 1, 2, 3, 4), (0, 2, 3, 1, 4), (1, 3, 0, 2, 4)]
    common = (26, 20, 25, 18, 0)
    for p, q, order in zip(C5_PATHS, qs, target_orders):
        s = G.tau(G.induced(c5, p))
        need(is_path(s, q) and not blue(s, q), "order-five black witness")
        v = G.tau(G.induced(s, q))
        need(G.induced(v, order) == common, "order-five identical hub witness")
        need(G.hp_subset_dp(v) == 5, "order-five hub path count")
        print("n=5: source path", p, "; H(partner)=", G.hp_subset_dp(s),
              "; black rerouting", q, "; hub relabelling", order)
    need(common[-1] == 0, "small hub has a sink")
    pinned = []
    for p in C5_PATHS:
        s = G.tau(G.induced(c5, p))
        pinned.append(sum(is_path(s, q) for q in permutations(range(5))
                          if q[0] == 0 and q[-1] == 4))
    need(pinned == [1, 3, 1], "fixed-port hostile")
    print("Pinned-endpoint HP counts:", pinned, "(two branches cannot reroute locally).")

    cases = 0
    by_weight = []
    for weight in range(1, 13):
        count = 0
        for word in words(weight):
            check_word(word)
            count += 1
            cases += 1
        by_weight.append(count)
    print("All A-words of weight 1..12:", by_weight, "; total", cases)
    print("All three labelled common-hub routes checked per word; total", 3 * cases,
          "; tournament orders 7..29.")
    for n in (3, 5, 7, 9):
        a = cyclic(n)
        p = tuple(range(n))
        need(blue(a, p), "regular blue base")
        need(all(blue(a, p[i:] + p[:i]) for i in range(n)), "regular rotation hostile")
    print("Regular hostile: all cycle cuts remain blue in cyclic controls n=3,5,7,9.")

    small_centers = {}
    for n in range(3, 8):
        mixed_sources, boundary_failures, edge_failures, rotation_cases = 0, 0, 0, 0
        for u in blue_tilings(n):
            paths = hamilton_paths(u)
            black_paths = [p for p in paths if not blue(u, p)]
            if n <= 5 and black_paths:
                small_centers.setdefault(G.canonical(u), (u, paths))
            if strong(u):
                continue
            if n >= 4:
                for p in paths:
                    s = complement_along(u, p)
                    q = p[1:] + p[:1]
                    need(is_path(s, q) and not blue(s, q), "nonstrong rotation lemma")
                    need(u[p[0]].bit_count() + u[p[1]].bit_count() > n,
                         "strict source-degree inequality")
                    rotation_cases += 1
            if black_paths:
                mixed_sources += 1
                if not any((p[0], p[1], p[-1]) == (0, 1, n - 1) for p in black_paths):
                    boundary_failures += 1
                witness = None
                for p in black_paths:
                    if (p[0], p[-1]) != (0, n - 1):
                        continue
                    for k in range(1, n):
                        if (k - 1, k) not in tuple(zip(p, p[1:])):
                            continue
                        end_sum = (u[k - 1].bit_count() + u[k].bit_count()
                                   - int(k == 1) + int(k == n - 1))
                        if end_sum != n - 1:
                            witness = p, k
                            break
                    if witness:
                        break
                if witness is None:
                    edge_failures += 1
                else:
                    p, k = witness
                    q0 = tuple(range(k, n)) + tuple(range(k))
                    cut = p.index(k)
                    q1 = p[cut:] + p[:cut]
                    s0 = G.tau(u)
                    s1 = complement_along(u, p)
                    hub0 = complement_along(s0, q0)
                    hub1 = complement_along(s1, q1)
                    need(hub0 == hub1 and not blue(s0, q0) and not blue(s1, q1),
                         "generalized shared-edge black hub")
        print("All blue presentations n=%d: mixed nonstrong=%d; anchored failures=%d;"
              " shared-edge failures=%d; all-path rotation controls=%d" %
              (n, mixed_sources, boundary_failures, edge_failures, rotation_cases))
        if n < 6:
            need(mixed_sources == 0, "minimal mixed nonstrong order")
        if n == 6:
            need(boundary_failures > 0, "first boundary-method hostile")

    hostile = ordered_join([(0,), G.tournament(4, 0, True), (0,)])
    hp = hamilton_paths(hostile)
    need(hostile == (62, 36, 40, 50, 38, 0), "anchored hostile rows")
    black_hp = [p for p in hp if not blue(hostile, p)]
    need(len(hp) == 5 and len(black_hp) == 2, "anchored hostile mixed profile")
    need([p for p in hp if (p[0], p[1], p[-1]) == (0, 1, 5)] == [tuple(range(6))],
         "no anchored alternate path")
    print("Minimal anchored hostile 1->R4->1, black paths:", black_hp)
    internal_black = (0, 4, 1, 2, 3, 5)
    internal_q0 = (2, 3, 4, 5, 0, 1)
    internal_q1 = (2, 3, 5, 0, 4, 1)
    internal_s1 = complement_along(hostile, internal_black)
    internal_hub = flip(flip(hostile, 0, 5), 1, 2)
    need(complement_along(G.tau(hostile), internal_q0) == internal_hub,
         "internal cut first hub")
    need(complement_along(internal_s1, internal_q1) == internal_hub,
         "internal cut second hub")
    need(not blue(G.tau(hostile), internal_q0) and not blue(internal_s1, internal_q1),
         "internal cut is black on both sides")
    need(hostile[1].bit_count() + hostile[2].bit_count() == 4,
         "internal degree-imbalance certificate")
    print("Same hostile repaired by internal shared edge (1,2), degree sum 4 != 5;",
          "hub rows", internal_hub)
    for n in range(3, 8):
        need(set(nonstrong_blue_sources(n)) ==
             {u for u in blue_tilings(n) if not strong(u)},
             "canonical SCC source generator agrees with complete small blue cube")
    for n in (8, 9):
        sources = list(nonstrong_blue_sources(n))
        mixed, anchored_fail, edge_fail = 0, 0, 0
        for u in sources:
            black_paths = [p for p in hamilton_paths(u) if not blue(u, p)]
            if not black_paths:
                continue
            mixed += 1
            if not any((p[0], p[1], p[-1]) == (0, 1, n - 1) for p in black_paths):
                anchored_fail += 1
            good = False
            for p in black_paths:
                if (p[0], p[-1]) != (0, n - 1):
                    continue
                for u0, v0 in zip(p, p[1:]):
                    if v0 != u0 + 1:
                        continue
                    end_sum = (u[u0].bit_count() + u[v0].bit_count()
                               - int(u0 == 0) + int(v0 == n - 1))
                    if end_sum != n - 1:
                        good = True
                        break
                if good:
                    break
            if not good:
                edge_fail += 1
        need((len(sources), mixed, anchored_fail, edge_fail) ==
             ((70, 67, 5, 3) if n == 8 else (519, 511, 9, 3)),
             "small endpoint-release boundary census")
        print("Canonical SCC sources n=%d: all=%d mixed=%d anchored failures=%d shared-edge failures=%d" %
              (n, len(sources), mixed, anchored_fail, edge_fail))
    r4 = G.tournament(4, 0, True)
    for middle in [(), (0,), cyclic(3), cyclic(5)]:
        u = ordered_join([r4, middle, G.transpose(r4)])
        paths = hamilton_paths(u)
        anchored = [p for p in paths if (p[0], p[-1]) == (0, len(u) - 1)]
        need(bool(anchored) and all(blue(u, p) for p in anchored),
             "two-endpoint obstruction family")
        need(any(not blue(u, p) for p in paths), "obstruction family is mixed")
    print("R4 -> Z -> R4^op endpoint obstruction: every fixed-endpoint path is blue;",
          "Z=empty,1,C3,C5 checked.")
    target = G.tau(hostile)

    def merged_key(t):
        return min(G.canonical(t), G.canonical(G.transpose(t)))

    start, finish = merged_key(hostile), merged_key(target)
    representatives, parent = {start: hostile}, {start: None}
    todo = deque([start])
    while todo and finish not in representatives:
        current = todo.popleft()
        a = representatives[current]
        for p in hamilton_paths(a):
            if blue(a, p):
                continue
            b = G.tau(G.induced(a, p))
            nxt = merged_key(b)
            if nxt not in representatives:
                representatives[nxt] = b
                parent[nxt] = (current, p)
                todo.append(nxt)
    need(finish in representatives, "hostile black detour exists")
    route, current = [], finish
    while parent[current] is not None:
        previous, p = parent[current]
        route.append((previous, p, current))
        current = previous
    route.reverse()
    need(len(route) == 3, "shortest hostile detour length")
    prescribed = [(0, 2, 3, 4, 1, 5), (0, 1, 3, 4, 5, 2), (0, 1, 2, 4, 5, 3)]
    a = hostile
    for p in prescribed:
        need(is_path(a, p) and not blue(a, p), "explicit hostile detour black step")
        a = G.tau(G.induced(a, p))
        print("n=6 detour path", p, "-> rows", a)
    final_order = (0, 1, 4, 5, 2, 3)
    need(G.induced(a, final_order) == target, "exact terminal isomorphism")
    print("n=6 shortest black detour length=3; terminal relabelling", final_order)

    padded_cases = 0
    for weight in range(2, 9):
        for word in words(weight):
            a = ordered_join([(0,) if x == 1 else cyclic(3) for x in word])
            r = len(a)
            for center, center_paths in small_centers.values():
                m = len(center)
                u = ordered_join([a, center, G.transpose(a)])
                n = len(u)
                center_black = next(p for p in center_paths if not blue(center, p))
                p_black = tuple(range(r)) + tuple(r + v for v in center_black) + tuple(range(r + m, n))
                need(not blue(u, p_black), "padded black path")
                hub = flip(flip(u, 0, 1), 0, n - 1)
                for center_blue in center_paths:
                    if not blue(center, center_blue):
                        continue
                    p_blue = tuple(range(r)) + tuple(r + v for v in center_blue) + tuple(range(r + m, n))
                    need(blue(u, p_blue), "padded blue source path")
                    for p in (p_blue, p_black):
                        s = complement_along(u, p)
                        q = p[1:] + p[:1]
                        need(is_path(s, q) and not blue(s, q), "padded black rotation")
                        need(complement_along(s, q) == hub, "padded common hub")
                    padded_cases += 1
    print("Padded mixed-center three-black-line bypass controls:", padded_cases,
          "; center classes", len(small_centers), "; A weight 2..8.")

    for n in range(3, 8):
        data, counts = load_census(n)
        black_neighbors = defaultdict(set)
        blue_neighbors = defaultdict(set)
        for colour, target in [("black", black_neighbors), ("blue", blue_neighbors)]:
            for (u, v), count in counts[colour].items():
                if count:
                    target[u].add(v)
                    target[v].add(u)
        branches = [v for v in data if data[v]["type"] == "PB" and data[v]["mass"] == 3]
        for v in branches:
            neighbors = blue_neighbors[v]
            need(len(neighbors) == 3, "three distinct branch neighbors")
            hubs = set.intersection(*(black_neighbors[w] for w in neighbors))
            need(bool(hubs), "census black hub")
            print("Census n=%d branch=%d neighbors=%s common black hubs=%s" %
                  (n, v, sorted(neighbors), sorted(hubs)))
        true_bridges = bridges(data, [counts["blue"], counts["black"]])
        support_bridges = bridges(data, [counts["blue"], counts["black"]], support=True)
        need(not any(v in edge for v in branches for edge in support_bridges),
             "branch has no support bridge")
        print("Census n=%d: multigraph bridges=%d; simple-support bridges=%d" %
              (n, len(true_bridges), len(support_bridges)))

    # Same n=5 nodes, masses, blue degrees and all pure-blue incidences,
    # but a deliberately unrealized mixed-locus rewiring.
    data, actual = load_census(5)
    fake_blue = Counter({(0, 11): 1, (2, 9): 1, (2, 6): 1, (3, 6): 1,
                         (5, 6): 1, (0, 2): 1, (0, 1): 1, (3, 3): 1})
    fake_black = Counter({(0, 1): 6, (1, 7): 2, (2, 4): 4,
                          (3, 4): 2, (4, 5): 4, (2, 3): 6})
    # Normalize key order before comparing the actual pure-blue incidences.
    fake_blue = Counter({tuple(sorted(k)): v for k, v in fake_blue.items()})
    need(degree(fake_blue) == degree(actual["blue"]), "same blue degree at every actual node")
    need(degree(fake_black) == degree(actual["black"]), "same black degree at every actual node")
    need(sum(fake_blue.values()) == 8 and sum(fake_black.values()) == 24, "n5 global counts")
    need(all(v % 2 == 0 for v in fake_black.values()), "doubled black lines")
    pure = {v for v in data if data[v]["type"] == "PB"}
    actual_pure = Counter({k: v for k, v in actual["blue"].items()
                           if v and (k[0] in pure or k[1] in pure)})
    fake_pure = Counter({k: v for k, v in fake_blue.items()
                         if v and (k[0] in pure or k[1] in pure)})
    need(actual_pure == fake_pure, "all actual pure-blue incidence retained")
    need(all(fake_black[tuple(sorted((4, v)))] for v in (2, 3, 5)),
         "trivalent neighbors retain a common black hub")
    fake_bridges = bridges(data, [fake_blue, fake_black])
    need(fake_bridges == {(0, 2), (0, 11), (2, 9)}, "unrealized mixed bridge")
    need(data[0]["type"] == data[2]["type"] == "MX", "bridge endpoints mixed")
    print("UNREALIZED n5 hostile: exact fibres, coloured degrees, pure-blue incidences and")
    print("common black hub retained, but mixed edge (0,2) is a multigraph bridge.")
    def current_row(counter, vertex):
        row = Counter()
        for (u, v), count in counter.items():
            if u == vertex:
                row[data[v]["triangles"] - data[vertex]["triangles"]] += count
            if v == vertex:
                row[data[u]["triangles"] - data[vertex]["triangles"]] += count
        return row
    need(current_row(actual["blue"], 0) == Counter({-3: 1, 1: 2}),
         "actual node-zero current inventory")
    need(current_row(fake_blue, 0) == Counter({-3: 1, 1: 1, 0: 1}),
         "fake node-zero current violation")
    need(current_row(actual["blue"], 3) == Counter({-1: 2, 1: 1}),
         "actual node-three current inventory")
    need(current_row(fake_blue, 3) == Counter({0: 2, 1: 1}),
         "fake node-three current violation")
    print("Incoming endpoint-current gate REJECTS the fake: node0 blue row")
    print("z^-3+2z becomes z^-3+z+1; node3 row 2z^-1+z becomes 2+z.")
    print("No actual tournament bridge is claimed. All controls passed.")


if __name__ == "__main__":
    main()
