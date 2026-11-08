#!/usr/bin/env python3
"""Independent literal-graph audit. Standard library; no producer imports/asserts.

Universe: all labelled A of orders 1..4, all root subsets, all labelled B of
orders 1..3; all labelled A of order 5/root subsets for coefficient and
strongness criteria; named direct fixed-menu families through 13 vertices.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations, combinations, product
from math import comb, factorial
from fractions import Fraction
import hashlib
import json

checks = Counter()

def require(condition, label, detail=None):
    checks[label] += 1
    if not condition:
        raise RuntimeError(json.dumps({'failure': label, 'detail': detail}, sort_keys=True))

def tournaments(n):
    pairs = tuple(combinations(range(n), 2))
    for code in range(1 << len(pairs)):
        rows = [0] * n
        for bit, (i, j) in enumerate(pairs):
            if code >> bit & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        yield tuple(rows)

def transitive(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))

def join(parts):
    rows = []
    total = sum(map(len, parts))
    offset = 0
    for part in parts:
        later = ((1 << total) - 1) ^ ((1 << (offset + len(part))) - 1)
        rows.extend((row << offset) | later for row in part)
        offset += len(part)
    return tuple(rows)

def graft(a, root, b):
    na, nb = len(a), len(b)
    bmask = ((1 << nb) - 1) << na
    return tuple(row | (0 if root >> i & 1 else bmask)
                 for i, row in enumerate(a)) + tuple((row << na) | root for row in b)

@lru_cache(None)
def hamilton(rows):
    """Literal vertex-level subset/last-vertex count, not a graft formula."""
    n = len(rows)
    if n == 0:
        return 1
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for used in range(1, size):
        for last, value in enumerate(dp[used]):
            if value:
                targets = rows[last] & (size - 1 - used)
                while targets:
                    low = targets & -targets
                    nxt = low.bit_length() - 1
                    dp[used | low][nxt] += value
                    targets -= low
    return sum(dp[-1])

@lru_cache(None)
def ordered_profile(rows):
    """Permutation/backward-adjacency cuts, distinct from subset HP count."""
    n = len(rows)
    out = [0] * (n + 1)
    if n == 0:
        return (1,)
    for p in permutations(range(n)):
        backward = sum(not (rows[x] >> y & 1) for x, y in zip(p, p[1:]))
        forward = n - 1 - backward
        for optional in range(forward + 1):
            out[1 + backward + optional] += comb(forward, optional)
    return tuple(out)

def coefficients(rows, root):
    """Erase B; mark mandatory broken A arcs and optional legal B slots."""
    n = len(rows)
    out = [0] * (n + 2)
    for p in permutations(range(n)):
        optional = int(bool(root >> p[0] & 1)) + int(not (root >> p[-1] & 1))
        mandatory = 0
        for x, y in zip(p, p[1:]):
            eligible = not (root >> x & 1) and bool(root >> y & 1)
            forward = bool(rows[x] >> y & 1)
            if not forward:
                if not eligible:
                    break
                mandatory += 1
            elif eligible:
                optional += 1
        else:
            for j in range(optional + 1):
                out[mandatory + j] += comb(optional, j)
    return tuple(out)

def skeleton_coefficients(rows, root):
    """Second coefficient engine: count quotient words with a repeatable B."""
    n = len(rows)
    full = (1 << n) - 1
    b_symbol = n
    @lru_cache(None)
    def rec(used, last, count):
        out = [0] * (n + 2)
        if used == full:
            out[count] += 1
        if count < n + 1 and last != b_symbol:
            if last == -1 or not (root >> last & 1):
                tail = rec(used, b_symbol, count + 1)
                out = [x + y for x, y in zip(out, tail)]
        for v in range(n):
            if used >> v & 1:
                continue
            legal = (last == -1 or (last == b_symbol and root >> v & 1)
                     or (last not in (-1, b_symbol) and rows[last] >> v & 1))
            if legal:
                tail = rec(used | (1 << v), v, count)
                out = [x + y for x, y in zip(out, tail)]
        return tuple(out)
    return rec(0, -1, 0)

def reachability(rows):
    reach = [row | 1 << v for v, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return reach

def strong(rows):
    full = (1 << len(rows)) - 1
    return all(row == full for row in reachability(rows))

def strong_condition(rows, root):
    n = len(rows)
    full = (1 << n) - 1
    reach = reachability(rows)
    source = sum(1 << v for v in range(n) if reach[v] == full)
    sink = sum(1 << v for v in range(n) if all(row >> v & 1 for row in reach))
    return bool(source & root) and bool(sink & (full ^ root))

def q_value(rows, j):
    return sum(comb(j, ell) * f for ell, f in enumerate(ordered_profile(rows)) if ell <= j)

def transformed(coeff, q):
    limit = len(coeff) - 1
    lam = [sum((-1) ** (ell - j) * comb(ell, j) * coeff[ell]
               for ell in range(max(1, j), limit + 1)) for j in range(limit + 1)]
    return sum(lam[j] * q(j) for j in range(limit + 1))

def literal_rooted_state(rows, root):
    """All-subset HP and root-start tables, without any graft transfer law."""
    n = len(rows)
    size = 1 << n
    full = size - 1
    dp = [[0] * n for _ in range(size)]
    start = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    start[1 << root][root] = 1
    for used in range(1, size):
        for last in range(n):
            value, fixed = dp[used][last], start[used][last]
            if value or fixed:
                targets = rows[last] & (full ^ used)
                while targets:
                    low = targets & -targets
                    nxt = low.bit_length() - 1
                    dp[used | low][nxt] += value
                    start[used | low][nxt] += fixed
                    targets -= low
    h = [sum(v) for v in dp]
    h[0] = 1
    sr = [sum(v) for v in start]
    u = sum(h[sub] * sr[full ^ sub] for sub in range(size) if not (sub >> root & 1))
    return (h[full], dp[full][root], u, h[full ^ (1 << root)])

def scc_sizes_in_order(rows):
    reach = reachability(rows)
    unseen = set(range(len(rows)))
    components = []
    while unseen:
        v = min(unseen)
        part = {w for w in unseen if reach[v] >> w & 1 and reach[w] >> v & 1}
        unseen -= part
        components.append(part)
    fixed_components = tuple(components)
    components.sort(key=lambda part: -sum(bool(rows[next(iter(part))] >> next(iter(other)) & 1)
                                          for other in fixed_components if other != part))
    return [len(part) for part in components]

def rational_rank(matrix):
    rows = [[Fraction(x) for x in row] for row in matrix]
    pivot = 0
    for col in range(len(rows[0])):
        found = next((i for i in range(pivot, len(rows)) if rows[i][col]), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        scale = rows[pivot][col]
        rows[pivot] = [x / scale for x in rows[pivot]]
        for i in range(len(rows)):
            if i != pivot:
                scale = rows[i][col]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot

def main():
    bank_b = [b for n in range(1, 4) for b in tournaments(n)]
    literal_digest = hashlib.sha256()
    for n in range(1, 6):
        for a in tournaments(n):
            for root in range(1 << n):
                coeff = coefficients(a, root)
                require(coeff[0] == hamilton(a), 'constant_is_HA')
                bound = min(root.bit_count(), n - root.bit_count()) + 1
                require(not any(coeff[bound + 1:]), 'two_sided_degree_bound')
                if n <= 4:
                    require(coeff == skeleton_coefficients(a, root), 'independent_skeleton_coefficients')
                compatible = strong_condition(a, root)
                require(strong(graft(a, root, (0,))) == compatible, 'source_sink_iff_singleton')
                if compatible:
                    require(coeff[2] >= 1, 'strong_compatible_c2_positive')
                if n <= 4:
                    for b in bank_b:
                        literal = hamilton(graft(a, root, b))
                        profile = ordered_profile(b)
                        predicted = sum(c * f for c, f in zip(coeff[1:], profile[1:]))
                        require(literal == predicted, 'literal_graph_expansion', [a, root, b, literal, predicted])
                        require(strong(graft(a, root, b)) == compatible, 'source_sink_iff_arbitrary_B')
                        if compatible:
                            require(literal >= 2 ** len(b) + 1, 'literal_sharp_strong_growth_floor')
                        literal_digest.update(f'{a}:{root}:{b}:{literal}\n'.encode())
                        require(literal == transformed(coeff, lambda j: q_value(b, j)), 'newton_reconstruction')
    # The full ordered profile counts empty-free bins, not arbitrary colourings.
    for n in range(1, 6):
        for b in tournaments(n):
            profile = ordered_profile(b)
            require(profile[1] == hamilton(b), 'profile_H_endpoint')
            if n >= 2:
                require(profile[2] >= 2 ** n - 2, 'universal_two_path_floor')
    # Sharp degree for each possible pair of positive colour class sizes up to 5.
    for r in range(1, 6):
        for c in range(1, 6):
            # A total order with first R and last C and maximal C->R transitions.
            word = [1]
            rr, cc = r - 1, c
            while rr and cc > 1:
                word += [0, 1]
                rr -= 1
                cc -= 1
            word += [1] * rr + [0] * cc
            root = sum(bit << i for i, bit in enumerate(word))
            # Skeleton DP avoids factorial permutation enumeration here.
            coeff = skeleton_coefficients(transitive(r + c), root)
            bound = min(r, c) + 1
            require(coeff[bound] > 0, 'degree_sharpness')
            require(not any(coeff[bound + 1:]), 'degree_sharpness_no_overflow')
    shell = transitive(4)
    root = 0b0101
    coeff = coefficients(shell, root)
    require(coeff == (1, 9, 10, 2, 0, 0), 'three_mode_shell_coefficients')
    triangle = (2, 4, 1)
    menu = (triangle, transitive(1), transitive(2))
    direct_modes = []
    for m in range(1, 4):
        b = join([triangle] * m)
        literal = hamilton(graft(shell, root, b))
        expected = -5 * 3 ** m + 4 * 12 ** m + 2 * 33 ** m
        require(literal == expected, 'literal_three_mode_family')
        require(strong(graft(shell, root, b)), 'literal_three_mode_strong')
        direct_modes.append({'m': m, 'vertices': 4 + 3 * m, 'H': literal})
    order_blind = []
    for word in product(range(len(menu)), repeat=3):
        parts = [menu[x] for x in word]
        b = join(parts)
        literal = hamilton(graft(shell, root, b))
        expected = transformed(coeff, lambda j: __import__('math').prod(q_value(part, j) for part in parts))
        require(literal == expected, 'literal_variable_menu_transform')
        sorted_literal = hamilton(graft(shell, root, join(sorted(parts))))
        require(literal == sorted_literal, 'literal_menu_order_blindness')
        order_blind.append(literal)
    # Sharp universal growth, with literal equality graphs of unbounded-form family.
    for n in range(1, 11):
        graph = graft(transitive(2), 1, transitive(n))
        require(hamilton(graph) == 2 ** n + 1, 'sharp_growth_equality_family')
        require(strong(graph), 'sharp_growth_equality_strong')
    # Full rooted-state collisions: actual graphs, no transfer-matrix imports.
    f4 = (8, 1, 3, 6)
    require(strong(f4) and hamilton(f4) == 5, 'collision_F4_control')
    require(q_value(f4, 2) == 32, 'collision_F4_Q2_control')
    collision_rows = []
    pairs = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (2, 1)]
    for k, ell in pairs:
        h, q = 3 ** k * 5 ** ell, 12 ** k * 32 ** ell
        expected = (3 * h + 12 * q, 3 * q, 6 * h + 14 * q - int(k + ell == 0), h + 4 * q)
        signatures = set()
        for word in sorted(set(permutations([0] * k + [1] * ell))):
            parts = [triangle if x == 0 else f4 for x in word] + [triangle]
            b = join(parts)
            graph = graft(transitive(2), 1, b)
            newroot = 2 + len(b) - 3
            state = literal_rooted_state(graph, newroot)
            require(state == expected, 'literal_four_state_collision', [k, ell, word, state, expected])
            require(strong(graph), 'literal_four_state_collision_strong')
            out = [row.bit_count() for row in graph]
            require([i for i, degree in enumerate(out) if degree == 1] == [0], 'collision_unique_outdegree_one')
            require([i for i, degree in enumerate(out) if degree == len(graph) - 2] == [1], 'collision_unique_indegree_one')
            signature = tuple(scc_sizes_in_order(b))
            require(signature == tuple([3 if x == 0 else 4 for x in word] + [3]), 'collision_recovered_SCC_word')
            signatures.add(signature)
        require(len(signatures) == comb(k + ell, k), 'collision_nonisomorphism_certificate_count')
        collision_rows.append({'k': k, 'ell': ell, 'state': expected, 'distinct_SCC_words': len(signatures)})
    # Minimal linear gadget summary: physical contexts and transitive gadgets.
    contexts = [(transitive(2), 1), (transitive(3), 1)]
    contexts += [(transitive(2 * m), (1 << m) - 1) for m in range(2, 5)]
    response = []
    for i, (a, roots) in enumerate(contexts):
        coeffs = skeleton_coefficients(a, roots)
        if i >= 2:
            m = i
            require(coeffs[m + 1] == factorial(m) ** 2, 'balanced_shell_top_coefficient')
        row = []
        for n in range(1, 6):
            literal = hamilton(graft(a, roots, transitive(n)))
            expected = transformed(coeffs, lambda j: j ** n)
            require(literal == expected, 'literal_linear_observer_response')
            row.append(literal)
        response.append(row)
    observer_ranks = []
    for r in range(1, 5):
        d = r + 1
        matrix = [row[:d] for row in response[:d]]
        rank = rational_rank(matrix)
        require(rank == d, 'minimal_linear_gadget_summary_rank')
        observer_ranks.append({'r': r, 'dimension': d, 'literal_response_rank': rank})
    # Hostile controls reject the tempting wrong shortcuts explicitly.
    hostile = {}
    proper = graft(shell, root, triangle)
    actual = hamilton(proper)
    missing_factorial = sum(coeff[k] * (ordered_profile(triangle)[k] // __import__('math').factorial(k))
                            for k in range(1, 4))
    require(actual != missing_factorial, 'hostile_unordered_profile_without_factorial')
    hostile['omit_ordering_factorial'] = {'correct': actual, 'wrong': missing_factorial}
    wrong_empty = sum(coeff[k] * q_value(triangle, k) for k in range(1, 4))
    require(actual != wrong_empty, 'hostile_empty_bins_in_positive_expansion')
    hostile['allow_empty_bins'] = {'correct': actual, 'wrong': wrong_empty}
    require(coeff[3] > 0, 'hostile_degree_root_count_only')
    hostile['degree_r_not_r_plus_one'] = {'r': 2, 'c3': coeff[3]}
    bad_shell = transitive(3)
    bad_root = 0b010
    require(not strong(graft(bad_shell, bad_root, triangle)), 'hostile_nonempty_proper_roots_not_enough')
    hostile['nonempty_proper_R_without_SCC_condition'] = {'A': list(bad_shell), 'R_mask': bad_root}
    require(hamilton(triangle) == hamilton(join([triangle, (0,)])), 'hostile_same_H_B_endpoint')
    one = hamilton(graft(shell, root, triangle))
    two = hamilton(graft(shell, root, join([triangle, (0,)])))
    require(one != two, 'hostile_scalar_H_B_insufficient')
    hostile['H_B_alone_insufficient'] = {'H_B': 3, 'graft_H': [one, two]}
    output = {
        'status': 'PASS',
        'scope': 'Finite-exact literal graph and independent coefficient audit; asymptotic density claim requires written proof, not this computation.',
        'checks': dict(sorted(checks.items())),
        'literal_count_digest_sha256': literal_digest.hexdigest(),
        'three_mode_family': direct_modes,
        'menu_words': len(order_blind),
        'minimal_linear_gadget_observer': observer_ranks,
        'physical_four_state_collision_families': collision_rows,
        'hostile_controls': hostile,
        'optimization_safe': 'All checks use explicit exceptions; no assert statements or producer imports.'
    }
    print(json.dumps(output, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
