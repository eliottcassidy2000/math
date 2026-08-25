#!/usr/bin/env python3
"""Independent hostile audit for THM-4069.

This version does not build the cycle space from a basis.  It scans all
subgraphs of K_n, retains the Eulerian ones, scans all (n-1)-edge subsets for
spanning trees, and recognizes simple cycles by their degree/connectivity
profile.  The overlap with the primary audit is therefore only the theorem's
definitions.
"""

from itertools import combinations, permutations


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def edges_of_complete_graph(n):
    return list(combinations(range(n), 2))


def degrees(mask, edges, n):
    result = [0 for _ in range(n)]
    for position, (left, right) in enumerate(edges):
        if (mask >> position) & 1:
            result[left] += 1
            result[right] += 1
    return result


def connected_support(mask, edges, n):
    present = [index for index, degree in enumerate(degrees(mask, edges, n))
               if degree]
    if not present:
        return False
    adjacency = [[] for _ in range(n)]
    for position, (left, right) in enumerate(edges):
        if (mask >> position) & 1:
            adjacency[left].append(right)
            adjacency[right].append(left)
    reached = {present[0]}
    frontier = [present[0]]
    while frontier:
        vertex = frontier.pop()
        for neighbor in adjacency[vertex]:
            if neighbor not in reached:
                reached.add(neighbor)
                frontier.append(neighbor)
    return reached == set(present)


def eulerian_masks(n):
    edges = edges_of_complete_graph(n)
    return [mask for mask in range(1 << len(edges))
            if all(degree % 2 == 0 for degree in degrees(mask, edges, n))]


def cycle_masks_by_recognition(n):
    edges = edges_of_complete_graph(n)
    cycles = []
    for mask in range(1, 1 << len(edges)):
        degree_list = degrees(mask, edges, n)
        if sum(degree != 0 for degree in degree_list) < 3:
            continue
        if all(degree in (0, 2) for degree in degree_list) \
                and connected_support(mask, edges, n):
            cycles.append(mask)
    return set(cycles)


def tree_masks_by_subset_scan(n):
    edges = edges_of_complete_graph(n)
    trees = []
    for chosen in combinations(range(len(edges)), n - 1):
        parent = list(range(n))

        def find(vertex):
            while parent[vertex] != vertex:
                parent[vertex] = parent[parent[vertex]]
                vertex = parent[vertex]
            return vertex

        acyclic = True
        mask = 0
        for position in chosen:
            left, right = edges[position]
            root_left = find(left)
            root_right = find(right)
            if root_left == root_right:
                acyclic = False
                break
            parent[root_left] = root_right
            mask |= 1 << position
        if acyclic and len({find(vertex) for vertex in range(n)}) == 1:
            trees.append(mask)
    return trees


def fundamental_cycles_from_tree_mask(n, tree_mask):
    edges = edges_of_complete_graph(n)
    edge_position = {edge: position for position, edge in enumerate(edges)}
    adjacency = [[] for _ in range(n)]
    for position, (left, right) in enumerate(edges):
        if (tree_mask >> position) & 1:
            adjacency[left].append(right)
            adjacency[right].append(left)
    cycles = []
    for chord_position, (source, target) in enumerate(edges):
        if (tree_mask >> chord_position) & 1:
            continue
        predecessor = {source: None}
        queue = [source]
        for vertex in queue:
            if vertex == target:
                break
            for neighbor in adjacency[vertex]:
                if neighbor not in predecessor:
                    predecessor[neighbor] = vertex
                    queue.append(neighbor)
        check(target in predecessor, "subset called a tree is disconnected")
        cycle = 1 << chord_position
        vertex = target
        while predecessor[vertex] is not None:
            previous = predecessor[vertex]
            cycle |= 1 << edge_position[tuple(sorted((vertex, previous)))]
            vertex = previous
        cycles.append(cycle)
    return cycles


def canonical_word(mask, n):
    edges = edges_of_complete_graph(n)
    edge_position = {edge: position for position, edge in enumerate(edges)}
    best = None
    for ordering in permutations(range(n)):
        candidate = 0
        for position, (left, right) in enumerate(edges):
            if (mask >> position) & 1:
                image = tuple(sorted((ordering[left], ordering[right])))
                candidate |= 1 << edge_position[image]
        if best is None or candidate < best:
            best = candidate
    return best


def quotient_edges(states, generators, n):
    canonical = {state: canonical_word(state, n) for state in states}
    vertices = sorted(set(canonical.values()))
    vertex_number = {word: number for number, word in enumerate(vertices)}
    edges = set()
    for state in states:
        source = vertex_number[canonical[state]]
        for generator in generators:
            target_state = state ^ generator
            check(target_state in canonical, "cycle generator left Eulerian universe")
            target = vertex_number[canonical[target_state]]
            if source != target:
                edges.add(tuple(sorted((source, target))))
    return vertices, edges


def tree_is_star(tree_mask, n):
    edge_list = edges_of_complete_graph(n)
    return max(degrees(tree_mask, edge_list, n)) == n - 1


def graph_name(mask, n):
    degree_list = tuple(sorted(degrees(mask, edges_of_complete_graph(n), n)))
    edge_count = mask.bit_count()
    names = {
        (0, (0, 0, 0, 0)): "empty",
        (3, (0, 2, 2, 2)): "triangle+isolated",
        (4, (2, 2, 2, 2)): "four-cycle",
    }
    return names.get((edge_count, degree_list),
                     "edges=%d,degrees=%s" % (edge_count, degree_list))


def main():
    print("THM-4069 independent subset-scan audit")
    print("universe=all edge subsets; tree_source=(n-1)-subset scan")
    print()

    print("cycle-space and generator-union census")
    print("  n all_subgraphs Eulerian iso_classes trees stars odd_bases simple_cycles union path_E envelope_E")
    for n in range(3, 7):
        all_edge_count = len(edges_of_complete_graph(n))
        states = eulerian_masks(n)
        canonical_classes = {canonical_word(state, n) for state in states}
        trees = tree_masks_by_subset_scan(n)
        expected_trees = n ** (n - 2)
        check(len(trees) == expected_trees, "Cayley tree count failed")
        recognized_cycles = cycle_masks_by_recognition(n)
        generator_union = set()
        star_count = 0
        odd_basis_count = 0
        for tree in trees:
            generators = fundamental_cycles_from_tree_mask(n, tree)
            generator_union.update(generators)
            star = tree_is_star(tree, n)
            odd_basis = all(generator.bit_count() % 2 == 1
                            for generator in generators)
            check(star == odd_basis,
                  "nonstar has all-odd basis or star has even generator")
            star_count += int(star)
            odd_basis_count += int(odd_basis)
        check(generator_union == recognized_cycles,
              "tree-generator union differs from recognized simple cycles")
        check(len(states) == 1 << ((n - 1) * (n - 2) // 2),
              "Eulerian subset count failed")
        edge_position = {edge: position
                         for position, edge in enumerate(edges_of_complete_graph(n))}
        path_tree = sum(1 << edge_position[(vertex, vertex + 1)]
                        for vertex in range(n - 1))
        path_generators = fundamental_cycles_from_tree_mask(n, path_tree)
        path_quotient = quotient_edges(states, path_generators, n)[1]
        envelope_quotient = quotient_edges(states, recognized_cycles, n)[1]
        check(path_quotient == envelope_quotient,
              "path quotient did not equal the canonical envelope")
        print("  %d %d %d %d %d %d %d %d %d %d %d"
              % (n, 1 << all_edge_count, len(states), len(canonical_classes),
                 len(trees), star_count, odd_basis_count,
                 len(recognized_cycles), len(generator_union),
                 len(path_quotient), len(envelope_quotient)))
    print()

    n = 4
    states = eulerian_masks(n)
    edge_list = edges_of_complete_graph(n)
    edge_position = {edge: position for position, edge in enumerate(edge_list)}
    star_mask = sum(1 << edge_position[(0, vertex)] for vertex in range(1, n))
    path_mask = sum(1 << edge_position[(vertex, vertex + 1)]
                    for vertex in range(n - 1))
    star_generators = fundamental_cycles_from_tree_mask(n, star_mask)
    path_generators = fundamental_cycles_from_tree_mask(n, path_mask)
    cycles = cycle_masks_by_recognition(n)
    star_vertices, star_edges = quotient_edges(states, star_generators, n)
    path_vertices, path_edges = quotient_edges(states, path_generators, n)
    envelope_vertices, envelope_edges = quotient_edges(states, cycles, n)
    check(star_vertices == path_vertices == envelope_vertices,
          "n=4 vertex quotients disagree")
    check(len(star_edges) == 2 and len(path_edges) == 3,
          "n=4 P3/K3 witness failed")
    check(path_edges == envelope_edges,
          "n=4 path quotient is not the all-cycle envelope")

    labels = {number: graph_name(mask, n)
              for number, mask in enumerate(star_vertices)}

    def labeled(edge_set):
        return sorted((labels[left], labels[right]) for left, right in edge_set)

    print("n=4 direct quotient witness")
    print("  vertices=%s" % [labels[index] for index in range(len(labels))])
    print("  star=%s" % labeled(star_edges))
    print("  path=%s" % labeled(path_edges))
    print("  envelope=%s" % labeled(envelope_edges))
    print()
    print("PASS: independent edge-subset and spanning-tree scans agree with THM-4069")


if __name__ == "__main__":
    main()
