"""Exact pair-tile, medial-map, flag, and small-tournament controls.

Run: python -X utf8 -B 04-computation/experiments/medial_pair_tiles_20261004.py
Only the standard library is used. Checks remain active with Python -O.
"""

from collections import Counter
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial
import json


def need(condition, message):
    if not condition:
        raise ValueError(message)


def graph_from_relation(objects, relation):
    adjacency = [set() for _ in objects]
    for i, j in combinations(range(len(objects)), 2):
        if relation(objects[i], objects[j]):
            adjacency[i].add(j)
            adjacency[j].add(i)
    return adjacency


def graph_edges(adjacency):
    return {(i, j) for i, neighbours in enumerate(adjacency) for j in neighbours if i < j}


def automorphisms(adjacency, edge_colors=None):
    """Exhaustive adjacency-preserving backtracking, with no known-group filter."""
    n = len(adjacency)
    matrix = [[-1] * n for _ in range(n)]
    for i, neighbours in enumerate(adjacency):
        for j in neighbours:
            matrix[i][j] = 0 if edge_colors is None else edge_colors[tuple(sorted((i, j)))]
    signatures = [tuple(sorted(Counter(matrix[i][j] for j in adjacency[i]).items())) for i in range(n)]
    domains = [[j for j in range(n) if signatures[j] == signatures[i]] for i in range(n)]
    assigned = {}
    used = set()

    def recurse():
        if len(assigned) == n:
            yield tuple(assigned[i] for i in range(n))
            return
        source = max((i for i in range(n) if i not in assigned),
                     key=lambda i: (sum(j in assigned for j in adjacency[i]), len(adjacency[i]), -i))
        for target in domains[source]:
            if target in used:
                continue
            if any(matrix[source][i] != matrix[target][j] for i, j in assigned.items()):
                continue
            assigned[source] = target
            used.add(target)
            yield from recurse()
            used.remove(target)
            del assigned[source]
    yield from recurse()


def pair_tiles(n, diagonal=False):
    iterator = combinations_with_replacement if diagonal else combinations
    return tuple(iterator(range(n), 2))


def pair_overlap_graph(n, diagonal=False):
    tiles = pair_tiles(n, diagonal)
    return tiles, graph_from_relation(tiles, lambda e, f: bool(set(e) & set(f)))


def medial_from_incidence(vertex_edge_sets, face_edge_sets, edge_count):
    return graph_from_relation(range(edge_count), lambda e, f:
        any(e in star and f in star for star in vertex_edge_sets)
        and any(e in boundary and f in boundary for boundary in face_edge_sets))


def map_data(name, faces):
    """Simple closed polyhedral maps with cyclic face boundaries."""
    faces = tuple(tuple(face) for face in faces)
    vertices = tuple(sorted(set().union(*map(set, faces))))
    edges = tuple(sorted({tuple(sorted((face[i], face[(i + 1) % len(face)])))
                          for face in faces for i in range(len(face))}))
    edge_index = {e: i for i, e in enumerate(edges)}
    vertex_edge_sets = tuple(frozenset(i for i, edge in enumerate(edges) if v in edge) for v in vertices)
    face_edge_sets = tuple(frozenset(edge_index[tuple(sorted((face[i], face[(i + 1) % len(face)])))]
                                     for i in range(len(face))) for face in faces)
    need(all(sum(e in f for f in face_edge_sets) == 2 for e in range(len(edges))), "two faces per edge")
    primal = graph_from_relation(vertices, lambda u, v: tuple(sorted((u, v))) in edge_index)
    medial = medial_from_incidence(vertex_edge_sets, face_edge_sets, len(edges))
    dual_medial = medial_from_incidence(face_edge_sets, vertex_edge_sets, len(edges))
    need(medial == dual_medial, "medial agrees exactly under dual incidence swap")
    flags = tuple((v, e, f) for f, boundary in enumerate(face_edge_sets)
                  for e in sorted(boundary) for v in edges[e])
    flag_graph = graph_from_relation(flags, lambda a, b: sum(x != y for x, y in zip(a, b)) == 1)
    edge_colors = {ij: next(k for k in range(3) if flags[ij[0]][k] != flags[ij[1]][k])
                   for ij in graph_edges(flag_graph)}
    need(len(flags) == 4 * len(edges), "four flags per edge")
    need(all(len(row) == 3 for row in flag_graph), "three flag neighbours")
    need(all({edge_colors[tuple(sorted((i, j)))] for j in row} == {0, 1, 2}
             for i, row in enumerate(flag_graph)), "one neighbour per flag color")
    return dict(name=name, vertices=vertices, edges=edges, faces=faces, primal=primal,
                vertex_edge_sets=vertex_edge_sets, face_edge_sets=face_edge_sets,
                medial=medial, flags=flags, flag_graph=flag_graph, edge_colors=edge_colors)


def components_with_colors(adjacency, colors, allowed):
    unseen = set(range(len(adjacency)))
    components = set()
    while unseen:
        component = {min(unseen)}
        pending = list(component)
        while pending:
            i = pending.pop()
            for j in adjacency[i]:
                if colors[tuple(sorted((i, j)))] in allowed and j not in component:
                    component.add(j)
                    pending.append(j)
        unseen -= component
        components.add(frozenset(component))
    return components


def audit_map(data):
    flags = data["flags"]
    recovered = []
    for retained_type in range(3):
        actual = components_with_colors(data["flag_graph"], data["edge_colors"], {0, 1, 2} - {retained_type})
        expected = {frozenset(i for i, flag in enumerate(flags) if flag[retained_type] == value)
                    for value in {f[retained_type] for f in flags}}
        need(actual == expected, "typed flag components recover original cells")
        recovered.append(len(actual))
    map_aut = list(automorphisms(data["flag_graph"], data["edge_colors"]))
    medial_aut = list(automorphisms(data["medial"]))
    black = set(data["vertex_edge_sets"])
    white = set(data["face_edge_sets"])
    color_actions = Counter()
    for permutation in medial_aut:
        black_image = {frozenset(permutation[e] for e in f) for f in black}
        white_image = {frozenset(permutation[e] for e in f) for f in white}
        if (black_image, white_image) == (black, white):
            color_actions["preserve"] += 1
        elif (black_image, white_image) == (white, black):
            color_actions["exchange"] += 1
        else:
            color_actions["fails_embedded_faces"] += 1
    need(color_actions["preserve"] == len(map_aut), "colored medial reconstruction")
    need(color_actions["fails_embedded_faces"] == 0, "abstract medial automorphisms respect these test embeddings")
    return dict(name=data["name"], VEF=recovered, flags=len(flags),
                map_automorphisms=len(map_aut), medial_automorphisms=len(medial_aut),
                medial_face_actions=dict(color_actions))


def cube_faces():
    faces = []
    for fixed in range(3):
        free = [i for i in range(3) if i != fixed]
        for value in (0, 1):
            face = []
            for x, y in ((0, 0), (1, 0), (1, 1), (0, 1)):
                bits = [0, 0, 0]
                bits[fixed], bits[free[0]], bits[free[1]] = value, x, y
                face.append(sum(bit << i for i, bit in enumerate(bits)))
            faces.append(tuple(face))
    return faces


def tournament_relabel(n, mask, permutation):
    pairs = pair_tiles(n)
    index = {pair: i for i, pair in enumerate(pairs)}
    result = 0
    for bit, (i, j) in enumerate(pairs):
        u, v = permutation[i], permutation[j]
        direction = (mask >> bit) & 1
        if u > v:
            u, v, direction = v, u, 1 - direction
        result |= direction << index[(u, v)]
    return result


def canonical_tournament(n, mask):
    return min(tournament_relabel(n, mask, p) for p in permutations(range(n)))


def tournament_scores(n, mask):
    scores = [0] * n
    for bit, (i, j) in enumerate(pair_tiles(n)):
        scores[i if (mask >> bit) & 1 else j] += 1
    return tuple(sorted(scores))


def order_tournament(order):
    position = {v: i for i, v in enumerate(order)}
    return sum(int(position[i] < position[j]) << bit for bit, (i, j) in enumerate(pair_tiles(len(order))))


def main():
    pair_controls = []
    for n in range(2, 7):
        tiles, graph = pair_overlap_graph(n, diagonal=True)
        automorphism_count = sum(1 for _ in automorphisms(graph))
        need(automorphism_count == factorial(n), "diagonal completion has exactly vertex relabelings")
        need({len(graph[i]) for i, e in enumerate(tiles) if e[0] == e[1]} == {n - 1}, "diagonal degree")
        need({len(graph[i]) for i, e in enumerate(tiles) if e[0] != e[1]} == {2 * n - 2}, "ordinary pair degree")
        pair_controls.append(dict(n=n, ordered=n*n, strict=n*(n-1)//2, with_diagonal=len(tiles),
                                  completed_automorphisms=automorphism_count))
    off_diagonal = {}
    for n in (3, 4, 5, 6):
        _, graph = pair_overlap_graph(n)
        count = sum(1 for _ in automorphisms(graph))
        need(count == (48 if n == 4 else factorial(n)), "Johnson graph symmetry boundary")
        off_diagonal[n] = count

    tetrahedron = map_data("tetrahedron", combinations(range(4), 3))
    cube = map_data("cube", cube_faces())
    pyramid = map_data("square pyramid", [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)])
    map_controls = [audit_map(data) for data in (tetrahedron, cube, pyramid)]
    need([(r["map_automorphisms"], r["medial_automorphisms"]) for r in map_controls] == [(24,48),(48,48),(8,16)], "positive and self-duality controls")
    _, tetra_pairs = pair_overlap_graph(4)
    need(tetra_pairs == tetrahedron["medial"], "four-vertex pair overlap is the tetrahedral medial")
    edges = tetrahedron["edges"]
    complement = tuple(edges.index(tuple(sorted(set(range(4)) - set(e)))) for e in edges)
    need(all((j in tetra_pairs[i]) == (complement[j] in tetra_pairs[complement[i]])
             for i in range(6) for j in range(6)), "complement is an extra medial automorphism")
    star0 = frozenset(i for i, e in enumerate(edges) if 0 in e)
    complement_star = frozenset(complement[i] for i in star0)
    need(complement_star in tetrahedron["face_edge_sets"] and complement_star not in tetrahedron["vertex_edge_sets"], "extra symmetry exchanges vertex star and face triangle")
    pyramid_line = graph_from_relation(pyramid["edges"], lambda e,f: bool(set(e)&set(f)))
    opposite_spokes = tuple(pyramid["edges"].index(e) for e in ((0,4),(2,4)))
    need(opposite_spokes[1] in pyramid_line[opposite_spokes[0]]
         and opposite_spokes[1] not in pyramid["medial"][opposite_spokes[0]], "line graph forgets the face-corner restriction")
    need((len(graph_edges(pyramid_line)),len(graph_edges(pyramid["medial"]))) == (18,16), "smallest polyhedral line/medial hostile")

    metric_vertices = ((0,0,0), (1,0,0), (0,2,0), (0,0,3))
    distance = lambda i,j: sum((metric_vertices[i][k] - metric_vertices[j][k])**2 for k in range(3))
    metric_aut = [p for p in permutations(range(4))
                  if all(distance(i,j) == distance(p[i],p[j]) for i,j in combinations(range(4),2))]
    need(len(metric_aut) == 1, "smallest polyhedral geometric symmetry hostile")

    n5_tiles, n5_overlap = pair_overlap_graph(5)
    table = Counter("equal" if e == f else "overlap" if set(e)&set(f) else "disjoint"
                    for e, f in product(n5_tiles, repeat=2))
    unordered = Counter("overlap" if set(e)&set(f) else "disjoint" for e,f in combinations(n5_tiles,2))
    need(table == {"equal":10,"overlap":60,"disjoint":30}, "100-cell orbital partition")
    need(unordered == {"overlap":30,"disjoint":15}, "45 unordered distinct pairs")
    need(all(len(row) == 6 for row in n5_overlap), "Johnson degree6")
    petersen = graph_from_relation(n5_tiles, lambda e,f: not(set(e)&set(f)))
    need(all(len(row) == 3 for row in petersen) and len(graph_edges(petersen)) == 15, "Petersen pair model")
    need(all(sum(v in e for e in n5_tiles) == 4 for v in range(5)), "five stars each have four pairs")

    n4_orbits = Counter(canonical_tournament(4, mask) for mask in range(64))
    need(sorted(n4_orbits.values()) == [8,8,24,24], "four tournament classes")
    n4_classes = [dict(scores=tournament_scores(4, mask), labelled_count=count) for mask,count in n4_orbits.items()]
    flag_index = {flag: i for i,flag in enumerate(tetrahedron["flags"])}
    face_index = {frozenset(face):i for i,face in enumerate(tetrahedron["faces"])}
    edge_index = {e:i for i,e in enumerate(edges)}
    transitive_to_flag = {}
    for order in permutations(range(4)):
        flag = (order[0], edge_index[tuple(sorted(order[:2]))], face_index[frozenset(order[:3])])
        transitive_to_flag[order_tournament(order)] = flag_index[flag]
    need(len(transitive_to_flag) == 24 and len(set(transitive_to_flag.values())) == 24, "transitive tournament to flag bijection")
    for a,b in combinations(transitive_to_flag,2):
        need(((a^b).bit_count() == 1) == (transitive_to_flag[b] in tetrahedron["flag_graph"][transitive_to_flag[a]]), "arc reversal equals one flag adjacency")

    free_pairs = [(i,j) for i,j in n5_tiles if j-i >= 2]
    fixed_path_pairs = [(i,j) for i,j in n5_tiles if j-i == 1]
    path_masks = [sum(((bits >> k)&1) << n5_tiles.index(e) for k,e in enumerate(free_pairs)) for bits in range(64)]
    n5_path_classes = {canonical_tournament(5,mask) for mask in path_masks}
    need((len(free_pairs),len(fixed_path_pairs),len(n5_path_classes)) == (6,4,12), "five-vertex tiling is a presentation of twelve classes")
    need(1*4 == 2*2 and len(set((1,4))) != len(set((2,2))), "product value loses diagonal type")

    report = dict(status="PROVED elementary constructions; FINITE-EXACT controls; no inferred F(G)/E(G) definition",
                  pair_tiles=pair_controls, off_diagonal_automorphisms=off_diagonal,
                  maps=map_controls,
                  smallest_polyhedral_hostile=dict(vertices=metric_vertices,
                      squared_edge_lengths=sorted(distance(i,j) for i,j in combinations(range(4),2)),
                      geometric_symmetry_order=1, abstract_map_order=24, abstract_medial_order=48,
                      complement_sends_star0_to_face_edges=[edges[i] for i in sorted(complement_star)]),
                  line_medial_hostile=dict(map="square pyramid", line_edges=18,medial_edges=16,
                      adjacent_only_in_line_graph=[[0,4],[2,4]],
                      n5_pair_graph_degree=6,medial_map_degree=4),
                  n5=dict(pair_vertices=10,stars=5,star_size=4,stars_per_pair=2,
                          overlap_degree=6,disjoint_degree=3,ordered_pair_table=dict(table),
                          unordered_distinct_pairs=dict(unordered),Petersen_edges=15,
                          complete_tournament_bits=10,fixed_path_bits=4,free_tiles=6,optional_loop_tiles=5,
                          fifteen_positions="6 free + 4 path + 5 loops",path_presentations=64,path_isomorphism_classes=12),
                  n4_tournament_classes=n4_classes, transitive_flag_bijection=dict(vertices=24,edges=36),
                  arithmetic_loss_hostile="1*4=2*2=4; off-diagonal pair and diagonal double remain different cells")
    print(json.dumps(report,indent=2))
    print("PASS: all checks remain active under -O")


if __name__ == "__main__":
    main()
