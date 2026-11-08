#!/usr/bin/env python3
"""Independent finite controls for the fixed finite graft-menu obstruction.

Only the standard library is used. No producer code or inherited verifier is
imported. Every guard raises explicitly, including under Python optimization.
The all-order argument is reviewed in AUDIT.md; this script is a finite test.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd, lcm
import json

RAYS = ((1, 1, 1, 1), (1, 0, 2, 1), (1, 0, 1, 0), (0, 0, 1, 0))
O = ((1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 2, 0), (0, 0, 0, 1))
Z = ((1, -1, 1, 0), (0, 0, 0, 1), (0, 0, 2, -1), (0, 0, 0, 1))
J = ((1, -1, 1, 0), (1, -1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 0))


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def apply(matrix, state):
    return tuple(sum(a * b for a, b in zip(row, state)) for row in matrix)


def cone_coordinates(state):
    D, E, U, d = state
    return E, d - E, D - d, U - D - d + E


def check_physical_constraints(state):
    D, E, U, d = state
    coordinates = cone_coordinates(state)
    need(min(coordinates) >= 0 and D >= 1 and d >= 1, ('cone', state))
    need(D % 2 == d % 2 == 1 and (U - E) % 2 == 0, ('parity', state))
    recovered = tuple(sum(coordinates[j] * RAYS[j][i] for j in range(4)) for i in range(4))
    need(recovered == state, ('ray expansion', state))


def all_tournaments(order):
    arcs = tuple(combinations(range(order), 2))
    for bits in range(1 << len(arcs)):
        rows = [0] * order
        for index, (a, b) in enumerate(arcs):
            winner, loser = (a, b) if bits & (1 << index) else (b, a)
            rows[winner] += 1 << loser
        yield tuple(rows)


def valid_path(rows, vertices):
    return all(rows[a] & (1 << b) for a, b in zip(vertices, vertices[1:]))


def permutation_profile(rows, root):
    """Count paths and two ordered bins directly, with no transfer formulas."""
    vertices = tuple(range(len(rows)))
    D = E = U = 0
    for word in permutations(vertices):
        good = valid_path(rows, word)
        D += good
        E += good and word[-1] == root
        cut = word.index(root)
        U += valid_path(rows, word[:cut]) and valid_path(rows, word[cut:])
    d = sum(valid_path(rows, word) for word in permutations(v for v in vertices if v != root))
    return D, E, U, d


def recursive_profile(rows, root):
    """Deletion recursion independently counts induced-subset path endpoints."""
    full = (1 << len(rows)) - 1

    @lru_cache(None)
    def ending(mask, end):
        if mask == 1 << end:
            return 1
        previous = mask ^ (1 << end)
        return sum(ending(previous, v) for v in range(len(rows))
                   if previous & (1 << v) and rows[v] & (1 << end))

    @lru_cache(None)
    def starting(mask, start):
        if mask == 1 << start:
            return 1
        later = mask ^ (1 << start)
        return sum(starting(later, v) for v in range(len(rows))
                   if later & (1 << v) and rows[start] & (1 << v))

    @lru_cache(None)
    def H(mask):
        if not mask:
            return 1
        return sum(ending(mask, v) for v in range(len(rows)) if mask & (1 << v))

    U = sum(starting(mask, root) * H(full ^ mask)
            for mask in range(1, full + 1) if mask & (1 << root))
    return H(full), ending(full, root), U, H(full ^ (1 << root))


def literal_graft(source, source_root, gadget, gadget_root):
    n = len(source)
    rows = list(source) + [row << n for row in gadget]
    for a in range(n):
        for b in range(len(gadget)):
            if a == source_root:
                rows[n + b] |= 1 << a
            else:
                rows[a] |= 1 << (n + b)
    return tuple(rows), n + gadget_root


def literal_O(source, root):
    n = len(source)
    return tuple(row | (1 << n) for row in source) + (0,), root


def literal_Z(source, root):
    rows, _ = literal_graft(source, root, (0,), 0)
    return rows, root


def bin_statistics(rows, root):
    """Enumerate permutations with every bin cut, including empty ordinary bins."""
    n = len(rows)
    H = Q = s = e = alpha = delta = T = 0
    for word in permutations(range(n)):
        if valid_path(rows, word):
            H += 1
            s += word[0] == root
            e += word[-1] == root
        for cut in range(n + 1):
            first, second = word[:cut], word[cut:]
            if valid_path(rows, first) and valid_path(rows, second):
                Q += 1
                alpha += bool(first) and first[0] == root
                delta += bool(first) and first[-1] == root
            if first and first[0] == root and valid_path(rows, first):
                for second_cut in range(cut, n + 1):
                    T += valid_path(rows, word[cut:second_cut]) and valid_path(rows, word[second_cut:])
    q = b = 0
    for word in permutations(v for v in range(n) if v != root):
        q += valid_path(rows, word)
        for cut in range(n):
            b += valid_path(rows, word[:cut]) and valid_path(rows, word[cut:])
    return H, Q, s, e, alpha, delta, T, q, b


def graft_matrix(stats, singleton):
    H, Q, s, e, alpha, delta, T, q, b = stats
    return ((H, -H, Q - H, 2 * H - Q),
            (e, -e, delta - e, e - delta),
            (alpha, s - alpha, 2 * T - alpha, 3 * alpha - s - 2 * T),
            (1, 0, 0, 0) if singleton else (q, -q, b - q, 2 * q - b))


def sharp_bound(matrix):
    images = [apply(matrix, ray) for ray in RAYS]
    need(all(image[0] > 0 for image in images), 'positive denominator on every ray')
    return max(Fraction(image[2], image[0]) for image in images)


def transitive(order):
    return tuple(sum(1 << j for j in range(i + 1, order)) for i in range(order))


def cyclic_block(m):
    rows = list(transitive(m)) + [0, 0]
    for v in range(m):
        rows[v] |= 1 << m
        rows[m + 1] |= 1 << v
    rows[m] |= 1 << (m + 1)
    return tuple(rows), m


def main():
    profile_count = sink_count = 0
    states = set()
    small_sources = []
    gadgets = []
    for n in range(1, 6):
        for rows in all_tournaments(n):
            for root in range(n):
                state = permutation_profile(rows, root)
                need(recursive_profile(rows, root) == state, ('independent path counts', rows, root))
                check_physical_constraints(state)
                profile_count += 1
                states.add(state)
                if not rows[root]:
                    need(state == (state[0],) * 4, ('sink equality', rows, root))
                    sink_count += 1
                if n <= 3:
                    small_sources.append((rows, root, state))
                if n <= 4:
                    stats = bin_statistics(rows, root)
                    H, Q, s, e, alpha, delta, T, q, b = stats
                    need((H, e, alpha, q) == state, ('profile and bin statistics', rows, root))
                    matrix = graft_matrix(stats, n == 1)
                    denominator = (H, H + Q, Q, Q - H)
                    numerator = (2 * alpha, 2 * alpha + 2 * T - s, 2 * T, 2 * T - alpha)
                    need(Q >= 2 * H > 0, ('strict denominator', rows, root))
                    for ray, den, num in zip(RAYS, denominator, numerator):
                        image = apply(matrix, ray)
                        need(image[0] == den and image[2] == num, ('ray values', rows, root, ray))
                    kappa = sharp_bound(matrix)
                    gadgets.append((rows, root, stats, matrix, kappa))
    need(profile_count == 5405, 'complete rooted universe through order five')
    need(len(small_sources) == 29 and len(gadgets) == 285, 'complete small rooted universes')

    ratio_checks = literal_graft_checks = operation_checks = 0
    for rows, root, stats, matrix, kappa in gadgets:
        for state in states:
            image = apply(matrix, state)
            check_physical_constraints(image)
            need(image[0] >= state[0], ('induced addition', state, rows, root))
            need(image[2] <= kappa * image[0], ('ratio', state, rows, root))
            ratio_checks += 1
        for source, source_root, state in small_sources:
            graph, mark = literal_graft(source, source_root, rows, root)
            actual = recursive_profile(graph, mark)
            need(actual == apply(matrix, state), ('literal graft', source, source_root, rows, root))
            need(actual[0] >= state[0] and actual[2] <= kappa * actual[0], 'literal graft bound')
            literal_graft_checks += 1
        # Applying any fixed gadget to a singleton gives the first-ray ratio.
        singleton_output = apply(matrix, (1, 1, 1, 1))
        need(singleton_output[0] == stats[0] and singleton_output[2] == 2 * stats[4], 'singleton boundary')
    for source, root, state in small_sources:
        for matrix, graph_op in ((O, literal_O), (Z, literal_Z)):
            graph, mark = graph_op(source, root)
            need(recursive_profile(graph, mark) == apply(matrix, state), 'literal O/Z')
            operation_checks += 1
        graph, mark = literal_O(source, root)
        graph, mark = literal_graft(graph, mark, (0,), 0)
        graph, mark = literal_graft(graph, mark, (0,), 0)
        need(recursive_profile(graph, mark)[0] == state[0] + 4 * state[2], 'literal OJJ readout')
        operation_checks += 1

    need(graft_matrix(bin_statistics((0,), 0), True) == J, 'singleton fourth row')
    need(sharp_bound(J) == sharp_bound(Z) == 2, 'sharp singleton and Z constants')
    need(apply(O, RAYS[3])[:3] == (0, 0, 2), 'O has no finite ray ratio')
    approximating_family_checks = 0
    for m in range(1, 8):
        rows, root = cyclic_block(m)
        need(recursive_profile(rows, root) == (2 ** m + 1, 1, 2 ** m + 1, 1), 'third ray family')
        need(recursive_profile(transitive(m + 1), 0) == (1, 0, 2 ** m, 1), 'fourth ray family')
        approximating_family_checks += 2
    # Exact convex combinations independently track the two limiting ratios.
    limit_checks = 0
    for _, _, _, matrix, kappa in gadgets:
        for j in (2, 3):
            limiting_image = apply(matrix, RAYS[j])
            target = Fraction(limiting_image[2], limiting_image[0])
            previous_gap = None
            for m in (1, 2, 4, 8, 16, 32):
                state = ((2 ** m + 1, 1, 2 ** m + 1, 1) if j == 2 else (1, 0, 2 ** m, 1))
                image = apply(matrix, state)
                actual = Fraction(image[2], image[0])
                gap = abs(actual - target)
                need(actual <= kappa and (previous_gap is None or gap <= previous_gap), 'ray approximation')
                previous_gap = gap
                limit_checks += 1
    for n in range(1, 10):
        # Larger orders are literal path counts rather than factorial bin enumeration.
        state = recursive_profile(transitive(n), 0)
        need(state == ((1, 1, 1, 1) if n == 1 else (1, 0, 2 ** (n - 1), 1)), 'transitive profile')
        output, mark = literal_graft((0,), 0, transitive(n), 0)
        actual = recursive_profile(output, mark)
        need(actual[0] == 1 and actual[2] == 2 ** n, 'no menu-independent ratio constant')

    B4 = (8, 1, 3, 6)
    selected = [((0,), 0), ((2, 0), 0), ((2, 0), 1), ((2, 4, 1), 0), (B4, 0), (B4, 1)]
    menu = [('O', O), ('Z', Z)]
    menu_bounds = {'Z': '2'}
    for index, (rows, root) in enumerate(selected):
        matrix = graft_matrix(bin_statistics(rows, root), len(rows) == 1)
        menu.append(('G' + str(index), matrix))
        menu_bounds['G' + str(index)] = str(sharp_bound(matrix))
    kappa = max(sharp_bound(matrix) for name, matrix in menu if name != 'O')
    need(kappa == 4, 'declared mixed-menu ratio')
    seeds = sorted(set(state for _, _, state in small_sources))
    sink_graph = (10, 12, 9, 0)  # Directed triangle followed by a global sink.
    sink_state = recursive_profile(sink_graph, 3)
    need(sink_state == (3, 3, 3, 3), 'nontrivial sink seed')
    seeds.extend([sink_state, recursive_profile(transitive(9), 0)])
    seeds = sorted(set(seeds))
    C = 25
    core_limit = (kappa * C).__floor__()
    need(any(state[0] <= C and state[2] > core_limit for state in seeds), 'large-U seed exception present')
    B_C = C * core_limit + sum(state[0] <= C for state in seeds)
    word_count = bounded_count = stationary_count = large_U_seed_tails = 0
    observed = set()
    max_depth = 5
    for seed in seeds:
        stack = [(seed, (seed[0], seed[2]), 0, False, '')]
        while stack:
            state, core, tail, had_non_O, word = stack.pop()
            D, E, U, d = state
            check_physical_constraints(state)
            need((D, U) == (core[0], core[1] * 2 ** tail), ('tail decomposition', seed, word))
            h, g = D + 2 * U, 2 * U
            need(h > g > 0 and 3 * g >= 2 * h and g % 2 == 0, 'endpoint chamber')
            endpoint = apply(J, apply(O, state))
            need(endpoint == (h, D, h, D), 'OJ endpoint')
            N = apply(J, endpoint)[0]
            need(N == D + 4 * U, 'OJJ scalar')
            if D <= C:
                bounded_count += 1
                observed.add(N)
                if had_non_O:
                    need(1 <= core[1] <= core_limit, ('finite core', seed, word, state, core))
                else:
                    need(core == (seed[0], seed[2]), 'all-O seed core')
                    large_U_seed_tails += core[1] > core_limit
            word_count += 1
            if len(word) < max_depth:
                for index, (name, matrix) in enumerate(menu):
                    after = apply(matrix, state)
                    need(after[0] >= D, 'monotone prefix pruning')
                    if name == 'O':
                        stack.append((after, core, tail + 1, had_non_O, word + str(index)))
                    else:
                        need(after[2] <= kappa * after[0], 'menu bound')
                        stack.append((after, (after[0], after[2]), 0, True, word + str(index)))
    expected_words = len(seeds) * sum(len(menu) ** depth for depth in range(max_depth + 1))
    need(word_count == expected_words, 'exhaustive word count')
    need(large_U_seed_tails > 0, 'large-U all-O tails tested')
    for X in (1, 4, 5, 17, 100, 1000, 10 ** 6, 10 ** 12):
        actual_count = sum(N <= X for N in observed)
        need(actual_count <= B_C * X.bit_length(), ('uniform count bound', X))
    for seed in seeds:
        for k in (0, 1, 2, 10, 100, 1000):
            state = (seed[0], 0 if k else seed[1], seed[2] * 2 ** k, seed[3])
            need(apply(J, apply(J, apply(O, state)))[0] == seed[0] + 4 * seed[2] * 2 ** k, 'long O tail')
    stationary = sink_state
    for _ in range(10000):
        stationary = apply(Z, stationary)
        need(stationary == sink_state, 'unbounded stationary Z run')
        stationary_count += 1

    congruence_checks = discarded_line_checks = 0
    for b in range(2, 13):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            for c in range(-20, 21):
                for N in range(-40, 101):
                    modulus = 2 * (a + b)
                    equivalent = (a * N - c) % modulus == 0
                    quotient_integer = (a * N - c) % (a + b) == 0
                    candidate_even_g = quotient_integer and ((a * N - c) // (a + b)) % 2 == 0
                    need(equivalent == candidate_even_g, 'exact even-g congruence')
                    if equivalent:
                        g = (a * N - c) // (a + b)
                        h = N - g
                        need(a * h - b * g == c and g % 2 == 0, 'line reconstruction')
                        if h > g > 0 and 3 * g >= 2 * h and 3 * a < 2 * b:
                            Hc = (-3 * c) // (2 * b - 3 * a)
                            need(h <= Hc and N <= 2 * Hc - 1, 'discarded line finite ceiling')
                            discarded_line_checks += 1
                    if abs(c) < a:
                        need((a - c) % modulus != 0, 'small-intercept common gap')
                    congruence_checks += 1
    # An explicit common-gap bank and its progression are independent controls.
    lines = [(3, 4, 2), (3, 4, -2), (4, 5, 3), (2, 3, 0)]
    modulus = lcm(*(2 * (a + b) for a, b, _ in lines))
    need(modulus % 2 == 0, 'common modulus even')
    for N in range(1, 100001, modulus):
        need(N % 2 == 1 and all((a * N - c) % (2 * (a + b)) != 0 for a, b, c in lines), 'common residue gap')
    hostile_count = 0
    for N in range(1, 2000, 2):
        matching = [c for c in (0, -2, -4, -6, -8) if (2 * N - c) % 10 == 0]
        need(len(matching) == 1, 'five-line residue partition')
        if N >= 17:
            c = matching[0]
            g, h = (2 * N - c) // 5, (3 * N + c) // 5
            need(h > g > 0 and g % 2 == 0 and 3 * g >= 2 * h, 'five-line positive chamber')
        hostile_count += 1
    for X in range(1, 100):
        for m in (2, 10, 14, 90):
            for r in range(1, m + 1, 2):
                count = sum(N % m == r % m for N in range(1, X + 1))
                need(count == max(0, 1 + (X - r) // m), 'exact progression count')

    print(json.dumps({
        'status': 'PASS',
        'implementation': 'standalone independent permutation/bin counts and deletion recursion; no imports from producer',
        'rooted_profile_presentations_through_order_5': profile_count,
        'distinct_source_profiles_through_order_5': len(states),
        'sink_root_presentations': sink_count,
        'rooted_gadgets_through_order_4': len(gadgets),
        'exact_ray_and_source_ratio_checks': ratio_checks,
        'literal_graft_pairs_all_roots_source_through_3_gadget_through_4': literal_graft_checks,
        'literal_O_Z_and_OJJ_checks': operation_checks,
        'literal_extreme_ray_family_checks': approximating_family_checks,
        'exact_limiting_ratio_checks': limit_checks,
        'mixed_menu_bound': str(kappa),
        'mixed_menu_operation_bounds': menu_bounds,
        'mixed_menu_alphabet_size': len(menu),
        'physical_seed_states': seeds,
        'word_depth': max_depth,
        'all_word_presentations': word_count,
        'bounded_intercept_word_presentations': bounded_count,
        'distinct_bounded_counts_in_word_universe': len(observed),
        'boundary_cap': C,
        'finite_core_U_cap': core_limit,
        'envelope_pair_bound': B_C,
        'all_O_large_U_seed_controls': large_U_seed_tails,
        'stationary_Z_steps': stationary_count,
        'largest_exact_O_tail_exponent': 1000,
        'even_g_line_congruence_checks': congruence_checks,
        'discarded_line_positive_candidate_checks': discarded_line_checks,
        'five_line_no_common_gap_odd_count_controls': hostile_count,
        'scope': 'Finite controls only. The proof, not enumeration, establishes arbitrary gadget orders and finite word lengths. No unrestricted finite-line or scalar-coverage conclusion.'
    }, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
