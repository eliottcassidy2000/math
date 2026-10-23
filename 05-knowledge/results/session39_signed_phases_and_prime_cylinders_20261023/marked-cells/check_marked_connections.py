"""Exact checks of cyclic-three-cell affine connections. Standard library only."""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def compose(f, g):
    """f after g, each encoded (epsilon, delta) over F3."""
    return (f[0] * g[0], (f[0] * g[1] + f[1]) % 3)


def inverse(f):
    return (f[0], (-f[0] * f[1]) % 3)


def apply(f, x):
    return (f[0] * x + f[1]) % 3


def gauge(edges, shifts):
    return {e: (f[0], (f[1] + shifts[e[1]] - f[0] * shifts[e[0]]) % 3)
            for e, f in edges.items()}


def connection_key(edges):
    return tuple((u, v, e, d) for (u, v), (e, d) in sorted(edges.items()))


def connection_canon(n, edges):
    return min(connection_key(gauge(edges, a)) for a in product(range(3), repeat=n))


def holonomies(n, edges, tree):
    adjacency = [[] for _ in range(n)]
    for u, v in tree:
        f = edges[u, v]
        adjacency[u].append((v, f))
        adjacency[v].append((u, inverse(f)))
    transport = [None] * n
    transport[0] = (1, 0)
    queue = [0]
    for u in queue:
        for v, f in adjacency[u]:
            if transport[v] is None:
                transport[v] = compose(f, transport[u])
                queue.append(v)
    need(all(f is not None for f in transport), 'tree does not connect every cell')
    return tuple(compose(inverse(transport[v]), compose(f, transport[u]))
                 for (u, v), f in sorted(edges.items()) if (u, v) not in tree)


def holonomy_key(n, edges, tree):
    hs = holonomies(n, edges, tree)
    conjugates = [tuple((e, (d + (1 - e) * a) % 3) for e, d in hs)
                  for a in range(3)]
    return (tuple(f[0] for _, f in sorted(edges.items())), min(conjugates))


def sections(n, edges):
    return [xs for xs in product(range(3), repeat=n)
            if all(xs[v] == apply(f, xs[u]) for (u, v), f in edges.items())]


def tournament(n, edges, polarity=None):
    """Named oriented C3 cells; absent blocks uniformly forward for i<j.

    polarity[e]=1 uses the selected matching as forward arcs, 2 its complement.
    """
    polarity = polarity or {}
    A = [[0] * (3 * n) for _ in range(3 * n)]
    for i in range(n):
        for x in range(3):
            A[3 * i + x][3 * i + (x + 1) % 3] = 1
    for i, j in combinations(range(n), 2):
        for x in range(3):
            for y in range(3):
                z = int(y == apply(edges[i, j], x)) if (i, j) in edges else 1
                if (i, j) in edges and polarity.get((i, j), 1) == 2:
                    z = 1 - z
                A[3 * i + x][3 * j + y] = z
                A[3 * j + y][3 * i + x] = 1 - z
    return A


def opposite(A):
    return [list(row) for row in zip(*A)]


def incidence(A):
    n = len(A)
    c = [0] * n
    edge = [[0] * n for _ in range(n)]
    for u, v, w in combinations(range(n), 3):
        if A[u][v] == A[v][w] == A[w][u]:
            for x in (u, v, w):
                c[x] += 1
            for x, y in ((u, v), (u, w), (v, w)):
                edge[x][y] += 1
                edge[y][x] += 1
    return c, edge


def signature(A, c=None, edge=None):
    if c is None:
        c, edge = incidence(A)
    return tuple(sorted(tuple(sorted((c[v], edge[u][v])
                                   for v in range(len(A)) if A[u][v]))
                        for u in range(len(A))))


def ternary(A, c, edge):
    return tuple(sorted(tuple(sum(1 for v in range(len(A))
                                   if A[u][v] and (c[v] + edge[u][v]) % 3 == z)
                              for z in range(3)) for u in range(len(A))))


def hamiltonian_count(A):
    n = len(A)
    outs = [sum(A[u][v] << v for v in range(n)) for u in range(n)]
    dp = [[0] * n for _ in range(1 << n)]
    for u in range(n):
        dp[1 << u][u] = 1
    for mask in range(1, 1 << n):
        for u, value in enumerate(dp[mask]):
            if not value:
                continue
            choices = outs[u] & ~mask
            while choices:
                bit = choices & -choices
                choices -= bit
                v = bit.bit_length() - 1
                dp[mask | bit][v] += value
    return sum(dp[-1])


def bitrows(A):
    return [sum(z << v for v, z in enumerate(row)) for row in A]


def verify_arc_gauge(A, B, shifts):
    permutation = [3 * i + (x + shifts[i]) % 3
                   for i in range(len(shifts)) for x in range(3)]
    return all(A[u][v] == B[permutation[u]][permutation[v]]
               for u in range(len(A)) for v in range(len(A)))


def three_cell_census():
    support = [(0, 1), (0, 2), (1, 2)]
    tree = {(0, 1), (1, 2)}
    by_canon = defaultdict(list)
    by_holo = defaultdict(set)
    for signs in product((-1, 1), repeat=3):
        for phases in product(range(3), repeat=3):
            edges = dict(zip(support, zip(signs, phases)))
            A = tournament(3, edges)
            need([sum(A[3 * i]) for i in range(3)] == [3, 4, 5], 'intrinsic cell scores')
            key = connection_canon(3, edges)
            hkey = holonomy_key(3, edges, tree)
            by_holo[hkey].add(key)
            hs = holonomies(3, edges, tree)
            roots = [x for x in range(3) if all(apply(h, x) == x for h in hs)]
            need(len(sections(3, edges)) == len(roots), 'point section criterion')
            need(all(verify_arc_gauge(A, tournament(3, gauge(edges, a)), a)
                     for a in product(range(3), repeat=3)), 'connection gauge not an arc isomorphism')
            by_canon[key].append((edges, len(roots)))
    need(len(by_canon) == 16, 'three-cell orbit census')
    need(len(by_holo) == 16 and all(len(v) == 1 for v in by_holo.values()),
         'holonomy classification differs from brute rotation orbits')
    records = []
    for key, members in sorted(by_canon.items()):
        edges, point_sections = members[0]
        A = tournament(3, edges)
        aut = sum(gauge(edges, a) == edges for a in product(range(3), repeat=3))
        need(len(members) * aut == 27, 'orbit-stabilizer failure')
        records.append({'canonical_edges': key, 'family_members': len(members),
                        'automorphisms': aut, 'point_sections': point_sections,
                        'hamiltonian_paths': hamiltonian_count(A),
                        'holonomies': holonomies(3, edges, tree)})
    need(Counter(r['automorphisms'] for r in records) == {1: 4, 3: 12}, 'automorphism histogram')
    return {'configurations': 216, 'ordinary_isotypes': 16, 'records': records}


def joint_cycle_census():
    support = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]
    tree = {(0, 1), (1, 2), (2, 3)}
    by_holo = {}
    orbit_sets = set()
    total = 0
    # Canonical orbit enumeration is independent of the tree construction.
    for signs in product((-1, 1), repeat=5):
        for phases in product(range(3), repeat=5):
            edges = dict(zip(support, zip(signs, phases)))
            key = connection_canon(4, edges)
            hkey = holonomy_key(4, edges, tree)
            if hkey in by_holo:
                need(by_holo[hkey] == key, 'joint holonomies lose an orbit')
            by_holo[hkey] = key
            orbit_sets.add(key)
            hs = holonomies(4, edges, tree)
            expected = sum(all(apply(h, x) == x for h in hs) for x in range(3))
            need(len(sections(4, edges)) == expected, 'joint point sections differ')
            total += 1
    need(len(by_holo) == len(orbit_sets), 'joint holonomy over-distinguishes orbits')
    hostile = {(0, 1): (1, 0), (1, 2): (1, 0), (2, 3): (1, 0),
               (0, 2): (-1, 0), (0, 3): (-1, 1)}
    hs = holonomies(4, hostile, tree)
    need(hs == ((-1, 0), (-1, 1)), 'joint reflection hostile changed')
    need(all(sum(apply(h, x) == x for x in range(3)) == 1 for h in hs), 'individual reflection has no point')
    need(not sections(4, hostile), 'inconsistent roots nevertheless have a section')
    return {'configurations': total, 'rotation_orbits': len(orbit_sets),
            'individual_fixed_points_do_not_glue': {'edges': connection_key(hostile),
                'based_holonomies': hs, 'individual_root_fixed_points': [[0], [2]],
                'joint_sections': 0}}


def cycle_edges(m, delta=0, signs=None):
    support = [(i, i + 1) for i in range(m - 1)] + [(0, m - 1)]
    support.sort()
    signs = signs or [1] * len(support)
    return {e: (s, delta if e == (0, m - 1) else 0)
            for e, s in zip(support, signs)}


def hostile_families():
    output = []
    for m in (3, 4, 5):
        records = []
        signatures = []
        converse_signatures = []
        for delta in range(3):
            A = tournament(m, cycle_edges(m, delta))
            c, edge = incidence(A)
            E = signature(A, c, edge)
            E_star = signature(opposite(A), c, edge)
            signatures.append(E)
            converse_signatures.append(E_star)
            K = ternary(A, c, edge)
            records.append({'closing_phase': delta, 'bitrows': bitrows(A),
                'cell_scores': [sum(A[3 * i]) for i in range(m)],
                'cell_triangle_incidences': c[::3],
                'hamiltonian_paths': hamiltonian_count(A), 'ternary_signature': K})
            need(all(c[3 * i:3 * i + 3] == [c[3 * i]] * 3 for i in range(m)), 'triangle cell not uniform')
            # Session38's oriented color difference is an exact vertex gradient.
            for u, v in combinations(range(3 * m), 2):
                if not A[u][v]:
                    u, v = v, u
                need(((c[v] + edge[u][v]) - (c[u] + edge[u][v])) % 3 == (c[v] - c[u]) % 3,
                     'ternary gradient identity')
        expected_H = {3: [2697, 2625, 2721], 4: [401085, 401397, 401721],
                      5: [45132777, 45150825, 45102765]}[m]
        need([r['hamiltonian_paths'] for r in records] == expected_H, 'hostile Hamiltonian counts')
        both_equal = len(set(signatures)) == len(set(converse_signatures)) == 1
        need(both_equal == (m == 5), 'exact-E blindness boundary changed')
        row = {'cells': m, 'records': records, 'exact_E_equal_in_both_orientations': both_equal}
        if m == 5:
            row['common_exact_E'] = signatures[0]
            row['common_converse_exact_E'] = converse_signatures[0]
            need(len({tuple(map(tuple, r['ternary_signature'])) for r in records}) == 1, 'ternary differs despite E equality')
        output.append(row)
    return output


def girth_five_checks():
    support = sorted(cycle_edges(5))
    count = 0
    for signs in product((-1, 1), repeat=5):
        base = dict(zip(support, zip(signs, (0,) * 5)))
        A0 = tournament(5, base)
        c0, e0 = incidence(A0)
        expected = (signature(A0, c0, e0), signature(opposite(A0), c0, e0))
        for phases in product(range(3), repeat=5):
            A = tournament(5, dict(zip(support, zip(signs, phases))))
            c, e = incidence(A)
            need((signature(A, c, e), signature(opposite(A), c, e)) == expected,
                 'girth-five signed phase blindness fails')
            count += 1
    controls = 0
    for bits in product((1, 2), repeat=5):
        polarity = dict(zip(support, bits))
        for signs in ((1,) * 5, (-1,) * 5, (1, -1, 1, -1, 1)):
            base = dict(zip(support, zip(signs, (0,) * 5)))
            A0 = tournament(5, base, polarity)
            c0, e0 = incidence(A0)
            expected = (signature(A0, c0, e0), signature(opposite(A0), c0, e0))
            for phases in ((0, 1, 2, 1, 2), (2, 0, 1, 2, 1)):
                A = tournament(5, dict(zip(support, zip(signs, phases))), polarity)
                c, e = incidence(A)
                need((signature(A, c, e), signature(opposite(A), c, e)) == expected,
                     'mixed polarity girth-five blindness fails')
                controls += 1
    return {'all_signed_phase_assignments_on_C5': count,
            'additional_mixed_polarity_controls': controls}


def marked_motif_histogram(A, cells, k):
    """Canonicalize only within occupied named cells, not between cells."""
    result = Counter()
    for subset in combinations(range(len(A)), k):
        blocks = [[v for v in subset if v // 3 == i] for i in range(cells)]
        occupancy = tuple(map(len, blocks))
        code = None
        for blocks_order in product(*(permutations(block) for block in blocks)):
            order = sum((tuple(b) for b in blocks_order), ())
            value = 0
            for i, j in combinations(range(k), 2):
                value = 2 * value + A[order[i]][order[j]]
            code = value if code is None else min(code, value)
        result[occupancy, code] += 1
    return result


def marked_refinement(A, root):
    n = len(A)
    colors = [u // 3 + 1 for u in range(n)]
    colors[root] = 0
    rounds = 0
    while len(set(colors)) < n:
        old_classes = sorted(set(colors))
        signatures = [(colors[u], tuple(sum(A[u][v] for v in range(n)
                                             if colors[v] == color)
                                        for color in old_classes)) for u in range(n)]
        palette = {s: i for i, s in enumerate(sorted(set(signatures)))}
        next_colors = [palette[s] for s in signatures]
        rounds += 1
        need(len(set(next_colors)) > len(set(colors)), 'marked refinement stalled')
        colors = next_colors
    return rounds


def strengthened_cycle_controls():
    output = []
    for m in range(3, 9):
        graphs = [tournament(m, cycle_edges(m, delta)) for delta in range(3)]
        Es = []
        all_rounds = []
        detectors = []
        for delta, A in enumerate(graphs):
            c, edge = incidence(A)
            Es.append(signature(A, c, edge))
            predicted_scores = [3 * m - 6] + [3 * m - 2 - 3 * i for i in range(1, m - 1)] + [5]
            need([sum(A[3 * i]) for i in range(m)] == predicted_scores, 'general intrinsic score formula')
            need(len(set(predicted_scores)) == m, 'cells not intrinsically distinguished')
            R = [3 * (m - 1 - i) + (-x) % 3 for i in range(m) for x in range(3)]
            need(all(A[u][v] == A[R[v]][R[u]] for u in range(3 * m) for v in range(3 * m)),
                 'literal converse map failed')
            for u in range(3 * m):
                seen = {u}
                todo = [u]
                for v in todo:
                    for w in range(3 * m):
                        if A[v][w] and w not in seen:
                            seen.add(w)
                            todo.append(w)
                need(len(seen) == 3 * m, 'cycle tournament not strong')
            rounds = [marked_refinement(A, root) for root in range(3 * m)]
            need(max(rounds) <= 1 + m // 2, 'one-mark propagation exceeded proved bound')
            all_rounds.append(Counter(rounds))
            Z = 0
            Q = 0
            for xs in product(range(3), repeat=m):
                path = all(xs[i + 1] == xs[i] for i in range(m - 1))
                Z += path and xs[-1] == (xs[0] + delta) % 3
                # a=xs[0], b=a+1 is the directed internal C3 arc.
                Q += path and xs[-1] == (xs[0] + 1 + delta) % 3
            need(Z == (3 if delta == 0 else 0), 'matching transversal detector')
            need(Q == (3 if delta == 2 else 0), 'one-internal-arc chirality detector')
            detectors.append({'phase': delta, 'matching_transversals': Z, 'split_root_chiral_motifs': Q})
        need((len(set(Es)) == 1) == (m >= 5), 'girth threshold controls differ')
        row = {'cells': m, 'self_converse_and_strong_all_phases': True,
               'all_actual_vertex_marks_tested_per_phase': 3 * m,
               'refinement_round_histograms': [dict(sorted(C.items())) for C in all_rounds],
               'marked_detectors': detectors,
               'exact_E_equal_all_phases': len(set(Es)) == 1}
        if m <= 5:
            motif_records = []
            for k in range(1, m + 2):
                hist = [marked_motif_histogram(A, m, k) for A in graphs]
                all_equal = hist[0] == hist[1] == hist[2]
                signs_equal = hist[1] == hist[2]
                need(all_equal == (k < m), 'marked zero/nonzero motif threshold')
                need(signs_equal == (k <= m), 'marked plus/minus motif threshold')
                motif_records.append({'vertices': k, 'distinct_marked_motifs': [len(h) for h in hist],
                                      'all_phases_equal': all_equal, 'nonzero_phases_equal': signs_equal})
            row['complete_cell_marked_induced_motif_histogram_comparisons'] = motif_records
        output.append(row)
    return output


def equitable_block_controls():
    regular = []
    for mask in range(512):
        B = [[(mask >> (3 * i + j)) & 1 for j in range(3)] for i in range(3)]
        rows = [sum(row) for row in B]
        cols = [sum(B[i][j] for i in range(3)) for j in range(3)]
        if len(set(rows)) != 1 or len(set(cols)) != 1:
            continue
        q = rows[0]
        need(cols[0] == q, 'equitable block double count')
        if q in (1, 2):
            P = B if q == 1 else [[1 - x for x in row] for row in B]
            pi = [row.index(1) for row in P]
            matches = [(e, d) for e in (-1, 1) for d in range(3)
                       if all((e * x + d) % 3 == pi[x] for x in range(3))]
            need(len(matches) == 1, 'selected matching lacks a unique affine description')
        regular.append(q)
    need(Counter(regular) == {0: 1, 1: 6, 2: 6, 3: 1}, 'regular binary block classification')
    return {'all_binary_3_by_3_blocks': 512, 'regular_blocks': len(regular),
            'counts_by_degree': dict(sorted(Counter(regular).items()))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = {'scope': 'Exact finite controls of named cyclic-three-cell connection theorems; not an unrestricted tournament wiring theorem.',
              'equitable_block_controls': equitable_block_controls(),
              'three_cell_census': three_cell_census(),
              'two_cycle_joint_census': joint_cycle_census(),
              'explicit_phase_hostiles': hostile_families(),
              'girth_five_controls': girth_five_checks(),
              'strengthened_cycle_controls': strengthened_cycle_controls(),
              'checks_before_receipt_comparison': CHECKS}
    result = json.loads(json.dumps(result))
    path = HERE / 'marked_connections_receipt.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2) + '\n')
    else:
        need(json.loads(path.read_text()) == result, 'frozen receipt differs')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
