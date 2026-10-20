#!/usr/bin/env python3
"""Exact eight-port matching atlas and cubic closure audit.

Universe: all 105 perfect matchings of F2^3, every ternary word on eight
ports, and all 105 closures of each of two specified eight-poles.
No random or sampled step; independent colour-count and F2^2-flow engines.
Run: python3 04-computation/boundary_support_20261020_seven.py
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import prod
import networkx as nx
import sympy as sp


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def matchings(vertices):
    if not vertices:
        yield ()
        return
    a = vertices[0]
    for j in range(1, len(vertices)):
        b = vertices[j]
        for tail in matchings(vertices[1:j] + vertices[j + 1:]):
            yield ((a, b),) + tail


ATLAS = tuple(matchings(tuple(range(8))))
KINDS = tuple(len({a ^ b for a, b in m}) for m in ATLAS)


def accepts(matching, word):
    return all(word[a] == word[b] for a, b in matching)


def component_count(first, second):
    graph = nx.Graph()
    graph.add_nodes_from(range(8))
    graph.add_edges_from(first + second)
    return nx.number_connected_components(graph)


def modular_rank(matrix, prime=101):
    """Independent finite-field elimination gives a rational-rank lower bound."""
    rows = [[int(x) % prime for x in row] for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(x * inverse) % prime for x in rows[rank]]
        for i in range(rank + 1, len(rows)):
            scale = rows[i][column]
            if scale:
                rows[i] = [(a - scale * b) % prime for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def graph_flow_count(vertices, edges):
    """All ordered binary even-subgraph pairs covering E, multiedge safe."""
    incidence = [0] * vertices
    for i, (a, b) in enumerate(edges):
        check(a != b, "loop not supported")
        incidence[a] ^= 1 << i
        incidence[b] ^= 1 << i
    rows = list(incidence)
    pivots = []
    r = 0
    for j in range(len(edges)):
        pivot = next((k for k in range(r, len(rows)) if rows[k] >> j & 1), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for k in range(len(rows)):
            if k != r and rows[k] >> j & 1:
                rows[k] ^= rows[r]
        pivots.append(j)
        r += 1
    basis = []
    for j in set(range(len(edges))) - set(pivots):
        vector = 1 << j
        for k, p in enumerate(pivots):
            if rows[k] >> j & 1:
                vector |= 1 << p
        basis.append(vector)
    vectors = [0]
    for vector in basis:
        vectors += [v ^ vector for v in vectors]
    full = (1 << len(edges)) - 1
    return sum((a | b) == full for a in vectors for b in vectors)


def pole_boundary_counts(vertices, internal, port_vertices):
    """Backtrack all proper partial colourings; retain ordered boundary."""
    edges = internal + tuple((v, -1) for v in port_vertices)
    masks = [0] * vertices
    assignment = [-1] * len(edges)
    counts = Counter()

    def recurse(left):
        if not left:
            counts[tuple(assignment[len(internal):])] += 1
            return
        best = None
        options = None
        for i in left:
            a, b = edges[i]
            blocked = masks[a] | (masks[b] if b >= 0 else 0)
            legal = tuple(c for c in range(3) if not (blocked >> c & 1))
            if not legal:
                return
            if options is None or len(legal) < len(options):
                best, options = i, legal
                if len(legal) == 1:
                    break
        a, b = edges[best]
        remaining = tuple(i for i in left if i != best)
        for c in options:
            masks[a] ^= 1 << c
            if b >= 0:
                masks[b] ^= 1 << c
            assignment[best] = c
            recurse(remaining)
            masks[a] ^= 1 << c
            if b >= 0:
                masks[b] ^= 1 << c
        assignment[best] = -1

    recurse(tuple(range(len(edges))))
    return counts


def pole_audit(name, source):
    deleted = {(i, i + 5) for i in range(4)}
    internal = tuple(sorted(tuple(sorted(e)) for e in source.edges()
                            if tuple(sorted(e)) not in deleted))
    ports = (0, 1, 2, 3, 5, 6, 7, 8)
    boundary = pole_boundary_counts(10, internal, ports)
    hist = {k: Counter() for k in (1, 2, 4)}
    simple_hist = {k: Counter() for k in (1, 2, 4)}
    uncolourable = []
    for matching, kind in zip(ATLAS, KINDS):
        count = sum(multiplicity for word, multiplicity in boundary.items()
                    if accepts(matching, word))
        edges = internal + tuple((ports[a], ports[b]) for a, b in matching)
        check(count == graph_flow_count(10, edges), "independent flow mismatch")
        hist[kind][count] += 1
        simple = len({tuple(sorted(e)) for e in edges}) == len(edges)
        if simple:
            simple_hist[kind][count] += 1
        if count == 0:
            g = nx.MultiGraph()
            g.add_nodes_from(range(10))
            g.add_edges_from(edges)
            uncolourable.append((kind, simple, len(list(nx.bridges(g))),
                                simple and nx.is_isomorphic(nx.Graph(g), nx.petersen_graph())))
    print("pole", name, "boundary_support", len(boundary),
          "boundary_multiplicities", sorted(Counter(boundary.values()).items()))
    print("closures_by_directions", {k: sorted(h.items()) for k, h in hist.items()})
    print("simple_closures_by_directions",
          {k: sorted(h.items()) for k, h in simple_hist.items()})
    print("uncolourable_kind_simple_bridges_isPetersen", sorted(Counter(uncolourable).items()))
    return boundary


def tournament_cube_audit():
    """The same eight-state matching object, with actual n=4 class fibres."""
    pairs = tuple(combinations(range(4), 2))
    tiles = ((0, 2), (0, 3), (1, 3))

    def bits(mask):
        return {e: 1 if e[1] == e[0] + 1 else (mask >> tiles.index(e) & 1)
                for e in pairs}

    def canonical(mask):
        word = bits(mask)
        codes = []
        for p in permutations(range(4)):
            code = 0
            for i, (a, b) in enumerate(pairs):
                u, v = p[a], p[b]
                bit = word[tuple(sorted((u, v)))] ^ (u > v)
                code |= bit << i
            codes.append(code)
        return min(codes)

    labels = tuple(canonical(mask) for mask in range(8))
    classes = sorted(set(labels))
    pi = tuple(classes.index(label) for label in labels)
    sigma = tuple(((x & 1) << 2) | (x & 2) | ((x & 4) >> 2) for x in range(8))
    converse = {pi[x]: pi[sigma[x]] for x in range(8)}
    check(all(converse[pi[x]] == pi[sigma[x]] for x in range(8)), "converse descends")
    merged_classes = sorted({min(c, converse[c]) for c in set(pi)})
    merged_pi = tuple(merged_classes.index(min(c, converse[c])) for c in pi)
    wiring_hist = Counter()
    merged_wiring_hist = Counter()
    commuting = Counter()
    admissible = Counter()
    line_wirings = {}
    for matching, kind in zip(ATLAS, KINDS):
        wiring = tuple(sorted(tuple(sorted((pi[a], pi[b]))) for a, b in matching))
        wiring_hist[wiring] += 1
        merged_wiring = tuple(sorted(tuple(sorted((merged_pi[a], merged_pi[b])))
                                     for a, b in matching))
        merged_wiring_hist[merged_wiring] += 1
        permutation = {}
        for a, b in matching:
            permutation[a], permutation[b] = b, a
        commutes = all(permutation[sigma[x]] == sigma[permutation[x]] for x in range(8))
        mixed_free = all(permutation[sigma[x]] != x for x in range(8))
        if commutes:
            commuting[kind] += 1
        if commutes and mixed_free:
            admissible[kind] += 1
            coloured_wiring = tuple(sorted((int(sigma[a] == a),
                                           *sorted((pi[a], pi[b])))
                                          for a, b in matching))
            line_wirings[coloured_wiring] = line_wirings.get(coloured_wiring, 0) + 1
    check(sorted(Counter(pi).values()) == [1, 1, 1, 5], "n4 class fibres")
    check(sorted(wiring_hist.values()) == [15, 15, 15, 60], "n4 pushforward wiring")
    check(sorted(Counter(merged_pi).values()) == [1, 2, 5], "n4 merged class fibres")
    check(sorted(merged_wiring_hist.values()) == [15, 30, 60], "n4 merged pushforward wiring")
    check(sum(commuting.values()) == 9 and sum(admissible.values()) == 6, "involution constraints")
    check(len(line_wirings) == 1, "eligible marked pairings should have same n4 wiring")
    print("tournament_n4_tile_order", tiles, "class_map", pi, "sigma", sigma)
    print("tournament_n4_all105_pushforward_wiring_hist", sorted(wiring_hist.items()))
    print("tournament_n4_converse_class_map", sorted(converse.items()),
          "merged_class_map", merged_pi)
    print("tournament_n4_all105_merged_pushforward_wiring_hist", sorted(merged_wiring_hist.items()))
    print("tournament_n4_sigma_commuting_by_directions", sorted(commuting.items()))
    print("tournament_n4_sigma_commuting_mixedfree_by_directions", sorted(admissible.items()))
    print("tournament_n4_eligible_coloured_wirings", sorted(line_wirings.items()))


def main():
    print("FINITE-EXACT: all 105 matching closures, all 6561 boundary words")
    check(len(ATLAS) == 105, "wrong matching count")
    check(Counter(KINDS) == {1: 7, 2: 42, 4: 56}, "wrong direction census")
    print("direction_census", sorted(Counter(KINDS).items()))
    four_sets = Counter(tuple(sorted({a ^ b for a, b in m}))
                        for m, k in zip(ATLAS, KINDS) if k == 4)
    check(len(four_sets) == 7 and set(four_sets.values()) == {8}, "four-direction fibres")
    print("four_direction_fibres", sorted(four_sets.items()))

    parity_words = []
    acceptance_hist = Counter()
    blind_hist = Counter()
    witness = None
    for word in product(range(3), repeat=8):
        sizes = Counter(word)
        if any(sizes[c] % 2 for c in range(3)):
            continue
        parity_words.append(word)
        matching_count = sum(accepts(m, word) for m in ATLAS)
        formula = prod(prod(range(1, sizes[c], 2)) for c in range(3))
        check(matching_count == formula, "monochromatic matching count")
        parallel = sum(accepts(m, word) for m, k in zip(ATLAS, KINDS) if k == 1)
        periods = sum(all(word[x] == word[x ^ d] for x in range(8)) for d in range(1, 8))
        check(parallel == periods, "parallel support versus period")
        shape = tuple(sorted(sizes.values(), reverse=True))
        acceptance_hist[(shape, matching_count, parallel)] += 1
        if parallel == 0:
            blind_hist[shape] += 1
            if witness is None:
                witness = word
    check(len(parity_words) == 1641, "wrong parity word count")
    print("parity_word_matching_parallel_hist", sorted(acceptance_hist.items()))
    print("seven_parallel_blind_words", sum(blind_hist.values()), sorted(blind_hist.items()))
    print("parallel_blind_word", witness, "accepted_matching",
          next(m for m in ATLAS if accepts(m, witness)))

    gram = sp.Matrix([[3 ** component_count(m, n) for n in ATLAS] for m in ATLAS])
    # The Gram matrix can be checked without constructing a 6561 x 105 matrix.
    columns = [sum(1 << j for j, word in enumerate(parity_words) if accepts(m, word))
               for m in ATLAS]
    for i, m in enumerate(ATLAS):
        check(gram[i, i] == 81, "wrong diagonal")
        for j, n in enumerate(ATLAS):
            check(gram[i, j] == (columns[i] & columns[j]).bit_count(), "Gram word count")
    print("gram_entry_hist", sorted(Counter(gram).items()))
    exact_rank = gram.to_DM().rank()
    mod_rank = modular_rank(gram.tolist())
    check(exact_rank == mod_rank == 91, "rank91 independent elimination")
    polynomial = gram.charpoly()
    t = polynomial.gen
    expected = t ** 14 * (t - 945) * (t - 210) ** 20 * (t - 120) ** 14 * (t - 30) ** 56
    check(sp.Poly(expected, t) == polynomial.as_poly(), "Gram characteristic polynomial")
    print("gram_rank_exact_and_mod101", exact_rank, mod_rank)
    print("gram_characteristic_polynomial", sp.factor(polynomial.as_expr()))
    print("gram_rank_by_atlas_class",
          {k: gram.extract([i for i, v in enumerate(KINDS) if v == k],
                           [i for i, v in enumerate(KINDS) if v == k]).to_DM().rank()
           for k in (1, 2, 4)})
    print("gram_rank_by_pairs_of_classes",
          {str(ks): gram.extract([i for i, v in enumerate(KINDS) if v in ks],
                                [i for i, v in enumerate(KINDS) if v in ks]).to_DM().rank()
           for ks in ((1, 2), (1, 4), (2, 4))})
    coefficient = {(2, 3): 2, (2, 2): -2, (2, 1): 0,
                   (4, 2): -1, (4, 1): 1}
    for i, matching in enumerate(ATLAS):
        if KINDS[i] != 1:
            continue
        row = [coefficient.get((k, component_count(matching, n)), 0)
               for n, k in zip(ATLAS, KINDS)]
        row[i] = -4
        check(sp.Matrix([row]) * gram == sp.zeros(1, 105),
              "parallel kernel recovery identity")
        # An independent direct word check, rather than using only Gram-nullity.
        for word in parity_words:
            check(sum(co * accepts(n, word) for co, n in zip(row, ATLAS)) == 0,
                  "parallel identity direct words")
    print("parallel_recovery_coefficients_for_4B", sorted(coefficient.items()),
          "verified_all_7_parallel_kernels")
    index = {frozenset(frozenset(pair) for pair in m): i for i, m in enumerate(ATLAS)}
    relations = []
    for left in combinations(range(1, 8), 3):
        left = (0,) + left
        right = tuple(x for x in range(8) if x not in left)
        row = [0] * 105
        for perm in permutations(right):
            inversions = sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4))
            m = frozenset(frozenset(pair) for pair in zip(left, perm))
            row[index[m]] = (-1) ** inversions
        relations.append(row)
    relation_matrix = sp.Matrix(relations)
    check(relation_matrix * gram == sp.zeros(35, 105), "determinantal relation")
    print("four_by_four_determinant_relations_count_rank", len(relations), relation_matrix.to_DM().rank())

    petersen = nx.petersen_graph()
    prism = nx.Graph()
    prism.add_edges_from((i, (i + 1) % 5) for i in range(5))
    prism.add_edges_from((i + 5, (i + 1) % 5 + 5) for i in range(5))
    prism.add_edges_from((i, i + 5) for i in range(5))
    pole_audit("Petersen_cut_four_spokes", petersen)
    pole_audit("pentagonal_prism_cut_four_spokes", prism)
    tournament_cube_audit()
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
