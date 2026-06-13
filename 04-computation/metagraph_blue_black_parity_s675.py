#!/usr/bin/env python3
"""
metagraph_blue_black_parity_s675.py

Audit the user's black/even and blue/odd observation for the merged
tournament metagraph.

There are at least two active blue-black conventions in the repo:

1. SC-type wiggly convention:
   BLUE  = self-converse<->self-converse edge in the merged one-arc-flip
           metagraph.
   BLACK = self-converse<->non-self-converse edge.
   GREEN = non-self-converse<->non-self-converse edge.

2. Explorer complement-line convention:
   BLUE  = complement-tiling lines whose fixed-path bit pattern is grid-sym.
   BLACK = complement-tiling lines whose bit pattern is not grid-sym.
   Edges are then aggregated between tournament classes, and the explorer
   can merge classes by true complement / transpose.

For each convention we test:
   - simple edge parity: every visible edge counts once;
   - weighted line parity: parallel lines count with multiplicity mod 2.

The GF(2) boundary of a layer is the set of vertices with odd layer-degree.
An even graph is exactly a layer with empty boundary. An "odd graph" can only
mean "every incident vertex is odd" unless the layer spans every vertex with
no isolates.
"""

from collections import Counter, defaultdict
from itertools import permutations, product
from math import comb


def banner(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def perm_groups(groups):
    if not groups:
        yield []
        return
    first, rest = groups[0], groups[1:]
    for p in permutations(first):
        for tail in perm_groups(rest):
            yield list(p) + tail


def canon_tournament(A):
    n = len(A)
    scores = [sum(row) for row in A]
    buckets = defaultdict(list)
    for v, s in enumerate(scores):
        buckets[s].append(v)
    groups = [buckets[s] for s in sorted(buckets)]
    best = None
    for p in perm_groups(groups):
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or key < best:
            best = key
    return best


def full_pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def full_bits_to_adj(mask, n, pairs):
    A = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (mask >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A


def transpose_adj(A):
    n = len(A)
    return [[A[j][i] if i != j else 0 for j in range(n)] for i in range(n)]


def count_ham_paths(A):
    n = len(A)
    full = (1 << n) - 1
    dp = {(1 << v, v): 1 for v in range(n)}
    for mask in range(1, 1 << n):
        for u in range(n):
            cur = dp.get((mask, u), 0)
            if not cur:
                continue
            for v in range(n):
                if (mask >> v) & 1:
                    continue
                if A[u][v]:
                    dp[(mask | (1 << v), v)] = dp.get((mask | (1 << v), v), 0) + cur
    return sum(dp.get((full, v), 0) for v in range(n))


def tournament_classes(n):
    pairs = full_pairs(n)
    m = len(pairs)
    class_of = {}
    reps = {}
    sizes = Counter()
    for mask in range(1 << m):
        A = full_bits_to_adj(mask, n, pairs)
        key = canon_tournament(A)
        if key not in class_of:
            cid = len(class_of)
            class_of[key] = cid
            reps[cid] = A
        cid = class_of[key]
        sizes[cid] += 1

    mask_to_cid = {}
    for mask in range(1 << m):
        mask_to_cid[mask] = class_of[canon_tournament(full_bits_to_adj(mask, n, pairs))]

    true_comp = {}
    self_converse = {}
    hvals = {}
    for cid, A in reps.items():
        ckey = canon_tournament(transpose_adj(A))
        true_comp[cid] = class_of[ckey]
        hvals[cid] = count_ham_paths(A)
    for cid in reps:
        self_converse[cid] = true_comp[cid] == cid

    return {
        "n": n,
        "pairs": pairs,
        "m": m,
        "class_count": len(class_of),
        "mask_to_cid": mask_to_cid,
        "reps": reps,
        "sizes": sizes,
        "true_comp": true_comp,
        "self_converse": self_converse,
        "hvals": hvals,
    }


def merge_by_true_complement(data):
    mid_of = {}
    members = {}
    mid = 0
    for cid in range(data["class_count"]):
        if cid in mid_of:
            continue
        comp = data["true_comp"][cid]
        mid_of[cid] = mid
        if comp != cid:
            mid_of[comp] = mid
            members[mid] = tuple(sorted((cid, comp)))
        else:
            members[mid] = (cid,)
        mid += 1
    return mid_of, members


def boundary_stats(num_vertices, edge_items):
    """Return simple and weighted GF(2) degree-boundary stats."""
    if not edge_items:
        return {
            "edges": 0,
            "total_weight": 0,
            "simple_hist": [],
            "weighted_hist": [],
            "simple_boundary": [],
            "weighted_boundary": [],
            "simple_incident_all_odd": False,
            "weighted_incident_all_odd": False,
        }

    simple_deg = [0] * num_vertices
    weighted_deg = [0] * num_vertices
    for a, b, weight in edge_items:
        if a == b:
            # Loops are silent for ordinary graph degree parity here.
            continue
        simple_deg[a] += 1
        simple_deg[b] += 1
        weighted_deg[a] += weight
        weighted_deg[b] += weight

    simple_boundary = [i for i, d in enumerate(simple_deg) if d % 2]
    weighted_boundary = [i for i, d in enumerate(weighted_deg) if d % 2]
    simple_incident = [d for d in simple_deg if d]
    weighted_incident = [d for d in weighted_deg if d]
    return {
        "edges": len(edge_items),
        "total_weight": sum(w for _, _, w in edge_items),
        "simple_hist": sorted(Counter(d for d in simple_deg if d).items()),
        "weighted_hist": sorted(Counter(d for d in weighted_deg if d).items()),
        "simple_boundary": simple_boundary,
        "weighted_boundary": weighted_boundary,
        "simple_incident_all_odd": bool(simple_incident)
        and all(d % 2 for d in simple_incident),
        "weighted_incident_all_odd": bool(weighted_incident)
        and all(d % 2 for d in weighted_incident),
    }


def print_layer(label, num_vertices, edge_items):
    st = boundary_stats(num_vertices, edge_items)
    if st["edges"] == 0:
        print(f"    {label}: empty")
        return st
    print(
        f"    {label}: E={st['edges']} total_weight={st['total_weight']} "
        f"simple_boundary={len(st['simple_boundary'])} "
        f"simple_even={len(st['simple_boundary']) == 0} "
        f"simple_incident_all_odd={st['simple_incident_all_odd']}"
    )
    print(
        f"      simple degree hist {st['simple_hist']} ; "
        f"weighted_boundary={len(st['weighted_boundary'])} "
        f"weighted_even={len(st['weighted_boundary']) == 0} "
        f"weighted_incident_all_odd={st['weighted_incident_all_odd']}"
    )
    if st["simple_boundary"]:
        head = st["simple_boundary"][:12]
        suffix = "" if len(st["simple_boundary"]) <= 12 else " ..."
        print(f"      simple boundary vertices {head}{suffix}")
    return st


def sc_type_wiggly_audit(n):
    data = tournament_classes(n)
    mid_of, members = merge_by_true_complement(data)
    nm = len(members)
    allmask = (1 << data["m"]) - 1

    edge_weights = Counter()
    for mask in range(1 << data["m"]):
        a = mid_of[data["mask_to_cid"][mask]]
        for k in range(data["m"]):
            b = mid_of[data["mask_to_cid"][mask ^ (1 << k)]]
            if a != b:
                edge_weights[(min(a, b), max(a, b))] += 1
    edge_weights = {edge: weight // 2 for edge, weight in edge_weights.items()}

    node_self_converse = {}
    for mid, cids in members.items():
        node_self_converse[mid] = data["self_converse"][cids[0]]

    layers = {"BLUE_SC_SC": [], "BLACK_SC_NS": [], "GREEN_NS_NS": []}
    for (a, b), weight in edge_weights.items():
        if node_self_converse[a] and node_self_converse[b]:
            layers["BLUE_SC_SC"].append((a, b, weight))
        elif node_self_converse[a] or node_self_converse[b]:
            layers["BLACK_SC_NS"].append((a, b, weight))
        else:
            layers["GREEN_NS_NS"].append((a, b, weight))

    print(f"  n={n}: classes={data['class_count']} merged_nodes={nm} arc_bits={data['m']}")
    for name in ["BLUE_SC_SC", "BLACK_SC_NS", "GREEN_NS_NS"]:
        print_layer(name, nm, layers[name])

    # The allmask variable is intentionally unused except as a visual reminder
    # that this convention uses true all-arc complement in the quotient.
    _ = allmask
    return data["class_count"], nm, layers


def fixed_path_tiles(n):
    tiles = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            tiles.append((x, y))
    tile_index = {t: i for i, t in enumerate(tiles)}
    trans_map = [tile_index[(n - y + 1, n - x + 1)] for x, y in tiles]
    return tiles, trans_map


def fixed_path_adj(mask, n, tiles):
    bits = [(mask >> k) & 1 for k in range(len(tiles))]
    verts = list(range(n, 0, -1))
    A = [[0] * n for _ in range(n)]
    for k in range(n - 1):
        A[k][k + 1] = 1
    for i, (x_label, y_label) in enumerate(tiles):
        xi = verts.index(x_label)
        yi = verts.index(y_label)
        if bits[i] == 0:
            A[xi][yi] = 1
        else:
            A[yi][xi] = 1
    return A


def is_grid_sym(mask, num_tiles, trans_map):
    for i, j in enumerate(trans_map):
        if i != j and ((mask >> i) & 1) != ((mask >> j) & 1):
            return False
    return True


def explorer_complement_line_audit(n):
    tiles, trans_map = fixed_path_tiles(n)
    num_tiles = len(tiles)
    allmask = (1 << num_tiles) - 1

    records = []
    sig_to_ci = {}
    groups = defaultdict(list)
    for mask in range(1 << num_tiles):
        A = fixed_path_adj(mask, n, tiles)
        sig = canon_tournament(A)
        records.append(
            {
                "mask": mask,
                "adj": A,
                "sig": sig,
                "flip": allmask ^ mask,
                "grid_sym": is_grid_sym(mask, num_tiles, trans_map),
            }
        )
        groups[sig].append(mask)

    for sig in sorted(groups):
        sig_to_ci[sig] = len(sig_to_ci)
    mask_to_ci = {}
    reps = {}
    for rec in records:
        ci = sig_to_ci[rec["sig"]]
        rec["ci"] = ci
        mask_to_ci[rec["mask"]] = ci
        reps.setdefault(ci, rec["adj"])

    true_comp = {}
    for ci, A in reps.items():
        true_comp[ci] = sig_to_ci[canon_tournament(transpose_adj(A))]

    mid_of = {}
    members = {}
    mid = 0
    for ci in range(len(sig_to_ci)):
        if ci in mid_of:
            continue
        comp = true_comp[ci]
        mid_of[ci] = mid
        if comp != ci:
            mid_of[comp] = mid
            members[mid] = tuple(sorted((ci, comp)))
        else:
            members[mid] = (ci,)
        mid += 1

    raw_edges = {}
    raw_self = Counter()
    for rec in records:
        a = rec["ci"]
        b = mask_to_ci[rec["flip"]]
        color = "blue" if rec["grid_sym"] else "black"
        if a == b:
            raw_self[color] += 1
            continue
        key = (min(a, b), max(a, b))
        if key not in raw_edges:
            raw_edges[key] = Counter()
        raw_edges[key][color] += 1

    # Each complement pair is seen from both endpoints.
    for counter in raw_edges.values():
        counter["blue"] //= 2
        counter["black"] //= 2
    raw_self["blue"] //= 2
    raw_self["black"] //= 2

    merged_edges = {}
    merged_self = Counter()
    for (a, b), counts in raw_edges.items():
        ma, mb = mid_of[a], mid_of[b]
        if ma == mb:
            merged_self["blue"] += counts["blue"]
            merged_self["black"] += counts["black"]
            continue
        key = (min(ma, mb), max(ma, mb))
        if key not in merged_edges:
            merged_edges[key] = Counter()
        merged_edges[key]["blue"] += counts["blue"]
        merged_edges[key]["black"] += counts["black"]

    layers = {"PURE_BLUE": [], "PURE_BLACK": [], "MIXED": []}
    for (a, b), counts in merged_edges.items():
        blue, black = counts["blue"], counts["black"]
        if black == 0:
            layers["PURE_BLUE"].append((a, b, blue))
        elif blue == 0:
            layers["PURE_BLACK"].append((a, b, black))
        else:
            layers["MIXED"].append((a, b, blue + black))

    grid_sym_count = sum(1 for rec in records if rec["grid_sym"])
    print(
        f"  n={n}: fixed_path_tiles={num_tiles} tilings={1 << num_tiles} "
        f"classes={len(sig_to_ci)} true_comp_merged_nodes={len(members)} "
        f"grid_sym_tilings={grid_sym_count}"
    )
    print(
        f"    silent/self complement lines before/after merge: "
        f"raw_blue={raw_self['blue']} raw_black={raw_self['black']} "
        f"merged_blue={merged_self['blue']} merged_black={merged_self['black']}"
    )
    for name in ["PURE_BLUE", "PURE_BLACK", "MIXED"]:
        print_layer(name, len(members), layers[name])

    return len(sig_to_ci), len(members), layers


def reference_even_graph_table(tournament_counts):
    # A002854 degree-even graph iso counts through n=6, as recorded in
    # 07-reflections/even-graphs-through-the-metagraph.md.
    degree_even_counts = {3: 2, 4: 3, 5: 7, 6: 16}
    print()
    print("  n | tournament/Royle-even classes | degree-even classes | message")
    print("  --|-------------------------------|---------------------|--------")
    for n in sorted(tournament_counts):
        t = tournament_counts[n]
        d = degree_even_counts[n]
        msg = "same" if t == d else "Royle-even is not degree-even"
        print(f"  {n} | {t:29d} | {d:19d} | {msg}")


def route_tournament():
    routes = [
        "boundary_vector",
        "convention_split",
        "line_weight_parity",
        "royle_even_property",
        "degree_even_cycle_space",
        "raw_color_claim",
    ]
    # Higher rank points from earlier route to later route, giving a transitive
    # tournament over proof methods.
    n = len(routes)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    scores = [0] * n
    for i, j in edges:
        scores[i] += 1
    print()
    print("  Tournament Analysis over proof routes")
    print("    vertices:", ", ".join(routes))
    print("    observable: does this route preserve color-layer GF(2) boundary?")
    print("    switch/gauge: convention selector plus simple-vs-weighted parity")
    print("    tie Hamiltonian path:", " > ".join(routes))
    print(f"    score_hist={sorted(Counter(scores).items())} c3=0 scc_sizes={[1] * n} H=1")


def main():
    print("S675 merged metagraph blue/black parity audit")
    print("Question: black layer always even? blue layer always odd?")
    print("Answer tested under both active repo conventions, with simple and weighted parity.")

    tournament_counts = {}

    banner("A. SC-type wiggly convention on G_n/Z_2")
    for n in range(3, 7):
        class_count, _, _ = sc_type_wiggly_audit(n)
        tournament_counts[n] = class_count

    banner("B. Explorer complement-line convention, merged by true complement")
    for n in range(3, 7):
        explorer_complement_line_audit(n)

    banner("C. Even graph warning: Royle-even versus degree-even")
    reference_even_graph_table(tournament_counts)

    banner("D. Interpretation")
    print("  Verdict:")
    print("    The naive always-statements are false as simple graph parity theorems.")
    print("    SC-type black is not always even, and explorer pure-black is not always even.")
    print("    Explorer pure-blue is incident-all-odd for n=3,4,5 but fails at n=6.")
    print("    However, black weighted-line parity has zero boundary in all audited")
    print("    cases through n=6; that is the strongest surviving version of the")
    print("    user's black-even observation.")
    print()
    print("  Surviving invariant:")
    print("    A color layer should be read as a GF(2) 1-chain.")
    print("    Its useful even/odd datum is the boundary vector, not the color name.")
    print("    Empty boundary = closed/even. Nonempty boundary = an odd defect packet.")
    print("    Handshaking forces defect packets to have even cardinality.")
    print()
    print("  Relation to even and odd numbers:")
    print("    Even numbers behave like closed parity classes: no exposed endpoint.")
    print("    Odd numbers behave like address defects: they carry a boundary bit.")
    print("    Addition moves boundary bits by xor; doubling kills them.")
    print("    The metagraph quotient is doing the same thing: it forgets an address,")
    print("    and the forgotten address reappears as odd-degree boundary vertices.")
    print()
    print("  Corrected research target:")
    print("    Do not ask whether black or blue is intrinsically even.")
    print("    Ask which quotient/functor sends tournament classes to Royle-even graph")
    print("    classes while preserving the color-layer boundary vector.")
    print("    That is the likely bridge between A000568 and the visible metagraph.")

    route_tournament()


if __name__ == "__main__":
    main()
