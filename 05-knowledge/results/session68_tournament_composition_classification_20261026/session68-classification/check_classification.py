#!/usr/bin/env python3
"""Standalone exact controls for decorated modular classification.

No repo input, third-party package, assert statement, or network required.
The exhaustive universes are recorded in the emitted receipt.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial
import hashlib
import json

CHECKS = Counter()


def require(condition, tag):
    CHECKS[tag] += 1
    if not condition:
        raise RuntimeError(tag)


def adj_from_mask(n, mask):
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = (mask >> k) & 1
        a[j][i] = 1 - a[i][j]
    return tuple(map(tuple, a))


def sub(a, vertices):
    return tuple(tuple(a[i][j] for j in vertices) for i in vertices)


def opposite(a):
    return tuple(tuple(0 if i == j else 1-a[i][j] for j in range(len(a)))
                 for i in range(len(a)))


def compose(q, blocks):
    address = [(i, j) for i, block in enumerate(blocks)
               for j in range(len(block))]
    return tuple(tuple(blocks[i][j][v] if i == u else q[i][u]
                       for u, v in address) for i, j in address)


C3 = ((0, 1, 0), (0, 0, 1), (1, 0, 0))
TT2 = ((0, 1), (0, 0))
ONE = ((0,),)
A = compose(TT2, [ONE, C3])
AO = opposite(A)


def automorphisms(a):
    n = len(a)
    degree = list(map(sum, a))
    candidates = [[v for v in range(n) if degree[v] == degree[u]]
                  for u in range(n)]
    order = sorted(range(n), key=lambda u: (len(candidates[u]), u))
    out, mapping, used = [], {}, set()

    def dfs(k):
        if k == n:
            out.append(tuple(mapping[i] for i in range(n)))
            return
        u = order[k]
        for v in candidates[u]:
            if v not in used and all(a[u][x] == a[v][y]
                                     for x, y in mapping.items()):
                mapping[u] = v
                used.add(v)
                dfs(k + 1)
                used.remove(v)
                del mapping[u]
    dfs(0)
    return out


def brute_code(a):
    n = len(a)
    return min(tuple(a[p[i]][p[j]] for i in range(n) for j in range(i+1, n))
               for p in permutations(range(n)))


@lru_cache(None)
def canonical_tree(a):
    """Exponential independent strong-module detector plus decorated node code."""
    n = len(a)
    if n == 1:
        return ("v",), 1
    whole = (1 << n) - 1
    modules = []
    for m in range(1, whole + 1):
        inside = [i for i in range(n) if m >> i & 1]
        if all(len({a[v][u] for u in inside}) == 1
               for v in range(n) if not (m >> v & 1)):
            modules.append(m)
    strong = [m for m in modules if all(not (m & t) or m & t == m
                                       or m & t == t for t in modules)]
    proper = [m for m in strong if m != whole]
    children = [m for m in proper if not any(m != t and m & t == m
                                            for t in proper)]
    vertices = [[i for i in range(n) if m >> i & 1] for m in children]
    require(sum(children) == whole, "children_partition")
    q = tuple(tuple(0 if i == j else a[v[0]][w[0]]
                    for j, w in enumerate(vertices)) for i, v in enumerate(vertices))
    child_data = [canonical_tree(sub(a, v)) for v in vertices]
    codes = [x[0] for x in child_data]
    internal = 1
    for _, order in child_data:
        internal *= order
    k = len(children)
    qautos = automorphisms(q)
    stabilizer = [g for g in qautos if all(codes[i] == codes[g[i]]
                                          for i in range(k))]
    scores = list(map(sum, q))
    if sorted(scores) == list(range(k)):
        order = sorted(range(k), key=lambda i: -scores[i])
        require(len(stabilizer) == 1, "linear_has_no_child_permutations")
        code = ("L", tuple(codes[i] for i in order))
    else:
        require(k >= 3, "prime_arity")
        for m in range(1, (1 << k) - 1):
            if m.bit_count() < 2:
                continue
            ii = [i for i in range(k) if m >> i & 1]
            require(any(len({q[v][u] for u in ii}) == 2
                        for v in range(k) if not (m >> v & 1)), "quotient_prime")
        options = []
        for p in permutations(range(k)):
            bits = tuple(q[p[i]][p[j]] for i in range(k) for j in range(i+1, k))
            options.append((bits, tuple(codes[p[i]] for i in range(k))))
        code = ("P", min(options))
    return code, internal * len(stabilizer)


def c3set(a):
    return tuple(t for t in combinations(range(len(a)), 3)
                 if all(sum(a[i][j] for j in t) == 1 for i in t))


def hypergraph_automorphism_count(n, edges):
    """Independent complete backtracking; degree/codegree only prune candidates."""
    edges = set(edges)
    degrees = [sum(i in t for t in edges) for i in range(n)]
    pairs = [[sum(i in t and j in t for t in edges) for j in range(n)]
             for i in range(n)]
    candidates = [[j for j in range(n) if degrees[i] == degrees[j]]
                  for i in range(n)]
    order = sorted(range(n), key=lambda i: (len(candidates[i]), i))
    mapping, used = {}, set()

    def dfs(k):
        if k == n:
            return 1
        i = order[k]
        total = 0
        for j in candidates[i]:
            if j in used or any(pairs[i][x] != pairs[j][y] for x, y in mapping.items()):
                continue
            if any((tuple(sorted((i, x, y))) in edges) !=
                   (tuple(sorted((j, mapping[x], mapping[y]))) in edges)
                   for x, y in combinations(mapping, 2)):
                continue
            mapping[i] = j
            used.add(j)
            total += dfs(k+1)
            used.remove(j)
            del mapping[i]
        return total
    return dfs(0)


def all_subset_hamilton(a):
    n = len(a)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    h = [0] * size
    for m in range(1, size):
        for j in range(n):
            if not (m >> j & 1):
                continue
            rem = m ^ (1 << j)
            dp[m][j] = 1 if not rem else sum(dp[rem][i] for i in range(n)
                                             if rem >> i & 1 and a[i][j])
        h[m] = sum(dp[m])
    return h


def path_cover_profile(a):
    n = len(a)
    h = all_subset_hamilton(a)
    pc = [[0] * (n+1) for _ in range(1 << n)]
    pc[0][0] = 1
    for m in range(1, 1 << n):
        first = m & -m
        s = m
        while s:
            if s & first:
                r = m ^ s
                for d in range(1, m.bit_count()+1):
                    pc[m][d] += h[s] * pc[r][d-1]
            s = (s-1) & m
    return pc[-1][1:]


def polynomial_next(p):
    cube = [0] * (3 * (len(p)-1) + 1)
    for i, a in enumerate(p):
        for j, b in enumerate(p):
            for k, c in enumerate(p):
                cube[i+j+k] += a*b*c
    out = [0] * (len(cube)+1)
    for i, value in enumerate(cube):
        diag = p[i//3] if i % 3 == 0 else 0
        require((value-diag) % 3 == 0, "polynomial_integrality")
        out[i] += (value-diag)//3
    for i, value in enumerate(p):
        out[3*i+1] += value
    return out


def word_code(word):
    if len(word) == 1:
        return word[0], 0
    step = len(word)//3
    children = [word_code(word[i*step:(i+1)*step]) for i in range(3)]
    codes = tuple(x[0] for x in children)
    key = min(codes[i:]+codes[:i] for i in range(3))
    return key, sum(x[1] for x in children) + int(codes[0] == codes[1] == codes[2])


def cyclic_types(previous):
    keys = list(previous)
    out = {}
    for a, b, c in product(keys, repeat=3):
        t = (a, b, c)
        key = min(t, (b, c, a), (c, a, b))
        exponent = previous[a]+previous[b]+previous[c]+int(a == b == c)
        require(key not in out or out[key] == exponent, "type_exponent_consistency")
        out[key] = exponent
    return out


def w2_permutations():
    # Independent literal imprimitive action on nine addresses.
    return [tuple(3*((i//3+top) % 3)+(i % 3+bottom[i//3]) % 3
                  for i in range(9))
            for top in range(3) for bottom in product(range(3), repeat=3)]


def main():
    receipts = {"schema": "session68-classification-v1"}
    classes = []
    total = 0
    for n in range(1, 6):
        by_tree, by_brute = {}, {}
        for mask in range(1 << (n*(n-1)//2)):
            a = adj_from_mask(n, mask)
            code, order = canonical_tree(a)
            direct = brute_code(a)
            require(code not in by_tree or by_tree[code] == direct, "tree_implies_iso")
            require(direct not in by_brute or by_brute[direct] == code, "iso_implies_tree")
            by_tree[code] = direct
            by_brute[direct] = code
            require(order == len(automorphisms(a)), "recursive_aut_vs_arc_backtracking")
            total += 1
        classes.append(len(by_tree))
    receipts["modular_census"] = {"labelled_inputs": total, "classes_n1_to_n5": classes}

    family = []
    reference_c3, reference_pc = None, None
    for bits in product(range(2), repeat=3):
        a = compose(C3, [AO if bit else A for bit in bits])
        c3 = c3set(a)
        pc = path_cover_profile(a)
        group = automorphisms(a)
        expected = 81 if len(set(bits)) == 1 else 27
        require(len(group) == expected, "twelve_vertex_aut")
        if reference_c3 is None:
            reference_c3, reference_pc = c3, pc
        require(c3 == reference_c3, "same_labelled_c3")
        require(pc == reference_pc, "same_full_path_cover_profile")
        require(pc[0] == 462915, "twelve_vertex_hamilton")
        require(all(factorial(d)*count % 81 == 0 for d, count in enumerate(pc, 1)),
                "max_aut_divides_ordered_path_covers")
        family.append({"bits": bits, "aut": len(group), "ordinary_path_orbits": pc[0]//len(group)})
    receipts["twelve_vertex_palette"] = {"cyclic_triples": len(reference_c3),
                                          "path_cover_profile": reference_pc,
                                          "members": family}

    # Full labelled C3 fibre, not merely the eight terminal-color presentations.
    realized, type_counts = set(), Counter()
    for top_reverse in range(2):
        for line_bits in product(range(2), repeat=3):
            for triangle_bits in product(range(2), repeat=3):
                blocks = []
                for line, triangle in zip(line_bits, triangle_bits):
                    block = [list(row) for row in A]
                    for j in range(1, 4):
                        block[0][j], block[j][0] = 1-line, line
                    for i in range(1, 4):
                        for j in range(1, 4):
                            if i != j:
                                block[i][j] = C3[i-1][j-1] ^ triangle
                    blocks.append(tuple(map(tuple, block)))
                a = compose(opposite(C3) if top_reverse else C3, blocks)
                require(c3set(a) == reference_c3, "full_fibre_labelled_c3")
                realized.add(a)
                # Explicit upper-slot and terminal-triangle gauges into the palette.
                slot = (0, 2, 1) if top_reverse else (0, 1, 2)
                target_colors = [0]*3
                g = []
                for i in range(3):
                    target_colors[slot[i]] = line_bits[i]
                    local = (0, 2, 1, 3) if triangle_bits[i] != line_bits[i] else (0, 1, 2, 3)
                    g.extend(4*slot[i]+local[j] for j in range(4))
                target = compose(C3, [AO if b else A for b in target_colors])
                require(all(a[i][j] == target[g[i]][g[j]] for i in range(12) for j in range(12)),
                        "full_fibre_gauge_into_palette")
                require({tuple(sorted(g[v] for v in edge)) for edge in reference_c3} == set(reference_c3),
                        "normalizing_gauge_preserves_hypergraph")
                type_counts[sum(line_bits)] += 1
    require(len(realized) == 128, "full_fibre_128_distinct_realizations")
    require(dict(type_counts) == {0: 16, 1: 48, 2: 48, 3: 16}, "full_fibre_type_orbit_sizes")
    hypergraph_aut = hypergraph_automorphism_count(12, reference_c3)
    require(hypergraph_aut == 1296, "full_hypergraph_aut_backtracking")
    receipts["full_twelve_vertex_realization_fibre"] = {
        "labelled_realizations": len(realized), "ordinary_types": 4,
        "hypergraph_automorphism_order": hypergraph_aut,
        "orbit_sizes_by_number_of_B_blocks": dict(type_counts),
        "type_separation": "Distinct B-block counts give distinct degree multisets.",
        "exhaustiveness": "Cited independent node-choice theorem; all 2^7 choices enumerated."}

    # A deliberately noncanonical partition: its block-preserving group is not full Aut.
    nine = compose(C3, [C3]*3)
    q = compose(C3, [C3, ONE, ONE])
    alternative = compose(q, [ONE, ONE, ONE, C3, C3])
    require(alternative == nine, "noncanonical_partition_same_tournament")
    full_group = automorphisms(nine)
    partition = [{0}, {1}, {2}, {3,4,5}, {6,7,8}]
    preserving = [g for g in full_group if all({g[v] for v in block} in partition
                                               for block in partition)]
    require(len(full_group) == 81 and len(preserving) == 27, "noncanonical_partition_hostile")
    receipts["partition_hostile"] = {"order": 9, "full_aut": len(full_group),
                                      "presentation_block_preserving_aut": len(preserving)}

    carousel9 = tuple(tuple(int(0 < (j-i) % 9 <= 4) for j in range(9)) for i in range(9))
    # Primeness checked by every nontrivial proper vertex subset, not a group-order inference.
    for mask in range(1, (1 << 9)-1):
        if mask.bit_count() < 2:
            continue
        inside = [i for i in range(9) if mask >> i & 1]
        require(any(len({carousel9[v][u] for u in inside}) == 2
                    for v in range(9) if not (mask >> v & 1)), "carousel9_prime_subset_control")
    require(len(automorphisms(carousel9)) == 9, "carousel9_aut")
    require(len(c3set(carousel9)) == len(c3set(nine)) == 30, "regular9_same_c3_count_different_hierarchy")
    receipts["prime_not_simple_group_control"] = {
        "carousel_order": 9, "carousel_automorphism_order": 9,
        "iterated_C3_automorphism_order": 81, "common_cyclic_triple_count": 30}

    polys = [[2]]
    for d in range(1, 6):
        polys.append(polynomial_next(polys[-1]))
    for d, p in enumerate(polys):
        m, internal = 3**d, (3**d-1)//2
        require(len(p) == internal+1, "full_exponent_range")
        require(all(c >= 2 and c % 2 == 0 for c in p), "complement_paired_full_support")
        require(p[0] == p[-1] == 2, "exactly_two_extreme_types")
        mass = sum(Fraction(c, 3**e) for e, c in enumerate(p)) * 3**internal
        require(mass == 2**m, "orbit_mass")
        if d:
            require(sum(p) == (sum(polys[d-1])**3+2*sum(polys[d-1]))//3,
                    "necklace_class_count")
    receipts["polynomials"] = [{"depth": d, "terminal_blocks": 3**d,
                                "coefficients": p, "class_count": sum(p)}
                               for d, p in enumerate(polys)]

    types = {0: 0, 1: 0}
    for d in range(4):
        histogram = Counter(types.values())
        require([histogram[e] for e in range(len(polys[d]))] == polys[d],
                "explicit_cyclic_tuple_enumeration")
        if d < 3:
            types = cyclic_types(types)

    group = w2_permutations()
    require(len(set(group)) == 81, "literal_w2_order")
    all_codes, orbit_reps = {}, set()
    for word in product(range(2), repeat=9):
        images = {tuple(word[g[i]] for i in range(9)) for g in group}
        stabilizer = sum(tuple(word[g[i]] for i in range(9)) == word for g in group)
        key, exponent = word_code(word)
        require(stabilizer == 3**exponent, "literal_action_vs_recursive_stabilizer")
        require(key not in all_codes or all_codes[key] == min(images), "word_code_orbit")
        all_codes[key] = min(images)
        orbit_reps.add(min(images))
    require(len(all_codes) == len(orbit_reps) == 24, "all_512_colorings_classified")
    receipts["independent_color_action"] = {"labelled_colorings": 512, "permutations": 81,
                                            "isomorphism_types": 24}
    receipts["checks"] = dict(sorted(CHECKS.items()))
    receipts["total_checks"] = sum(CHECKS.values())
    print(json.dumps(receipts, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
