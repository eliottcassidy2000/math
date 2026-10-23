"""Exact controls for arbitrary regular-five blocks and shared chart potentials.

No assertions: checks remain active under python -O. All cells and chart
coordinates are named; there is no independent re-rooting at adjacent blocks.
"""
from collections import Counter
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
        raise ValueError(message)


def shift(block, a, b):
    """Target block after source/target-owner translations by a and b."""
    p = len(block)
    mask = (1 << p) - 1
    return tuple(((block[(i-a) % p] << b) |
                  (block[(i-a) % p] >> (p-b))) & mask for i in range(p))


def stabilizer(block):
    p = len(block)
    return tuple((a, b) for a in range(p) for b in range(p)
                 if shift(block, a, b) == block)


def direct_relation(source, target):
    p = len(source)
    return tuple((a, b) for a in range(p) for b in range(p)
                 if all(((source[i] >> j) & 1) ==
                        ((target[(i+a) % p] >> ((j+b) % p)) & 1)
                        for i in range(p) for j in range(p)))


def regular_blocks(p, d):
    """Backtrack rows subject to every column having the same declared degree."""
    rows = [sum(1 << j for j in js) for js in combinations(range(p), d)]
    columns = [0] * p

    def rec(prefix):
        if len(prefix) == p:
            if all(c == d for c in columns):
                yield tuple(prefix)
            return
        remaining = p-len(prefix)-1
        for row in rows:
            nxt = [columns[j] + ((row >> j) & 1) for j in range(p)]
            if any(c > d or c+remaining < d for c in nxt):
                continue
            previous = columns[:]
            columns[:] = nxt
            yield from rec(prefix+[row])
            columns[:] = previous
    yield from rec([])


def affine_pattern(p, s, subset):
    return tuple(sum(1 << j for j in range(p) if (j-s*i) % p in subset)
                 for i in range(p))


def chosen_matching(block):
    """Actual owner map: singleton or winning vertex of selected target pair.

    Used only on blocks with nontrivial rotation stabilizer. A stored edge
    is traversed backwards using this map's inverse, not a fresh selection.
    """
    selected = block if block[0].bit_count() <= 2 else tuple(31-row for row in block)
    images = []
    for row in selected:
        entries = [j for j in range(5) if row >> j & 1]
        if len(entries) == 1:
            images.append(entries[0])
        else:
            a, b = entries
            images.append(a if (b-a) % 5 in (1, 2) else b)
    return tuple(images)


def selector_bank(block):
    selected = block if block[0].bit_count() <= 2 else tuple(31-row for row in block)
    winner = chosen_matching(block)
    if selected[0].bit_count() == 1:
        return (winner,)
    loser = tuple(next(j for j in range(5) if row >> j & 1 and j != winner[i])
                  for i, row in enumerate(selected))
    return (winner, loser)


def converse_block(block):
    p = len(block)
    return tuple(sum((1-((block[(-i) % p] >> ((-j) % p)) & 1)) << j
                     for j in range(p)) for i in range(p))


def solve(n, relations, p=5):
    """Root-domain propagation on a connected graph of partial bijections."""
    neighbors = [[] for _ in range(n)]
    for (i, j), pairs in relations.items():
        f = dict(pairs)
        g = {b: a for a, b in pairs}
        if len(f) != len(pairs) or len(g) != len(pairs):
            raise ValueError('relation is not a partial bijection')
        neighbors[i].append((j, f))
        neighbors[j].append((i, g))
    solutions = []
    for root in range(p):
        values = [None] * n
        values[0] = root
        queue = [0]
        good = True
        for i in queue:
            for j, f in neighbors[i]:
                if values[i] not in f:
                    good = False
                    break
                value = f[values[i]]
                if values[j] is None:
                    values[j] = value
                    queue.append(j)
                elif values[j] != value:
                    good = False
                    break
            if not good:
                break
        if good and all(v is not None for v in values):
            solutions.append(tuple(values))
    return solutions


def tournament(blocks, p=5):
    count = max(max(k) for k in blocks)+1
    out = [[0] * (count*p) for _ in range(count*p)]
    for k in range(count):
        for i in range(p):
            for j in range(p):
                out[k*p+i][k*p+j] = int(0 < (j-i) % p <= (p-1)//2)
    for (a, b), rows in blocks.items():
        for i in range(p):
            for j in range(p):
                value = (rows[i] >> j) & 1
                out[a*p+i][b*p+j] = value
                out[b*p+j][a*p+i] = 1-value
    return out


def hamiltonian_path(out):
    path = []
    for v in range(len(out)):
        position = next((i for i, u in enumerate(path) if out[v][u]), len(path))
        path.insert(position, v)
    need(all(out[u][v] for u, v in zip(path, path[1:])), 'path insertion failed')
    return path


def complement_along_path(out, path):
    kept = set(zip(path, path[1:]))
    return [[0 if i == j else out[i][j] if (i, j) in kept or (j, i) in kept
             else 1-out[i][j] for j in range(len(out))] for i in range(len(out))]


def refine(out, partition):
    current = tuple(tuple(sorted(C)) for C in partition)
    while True:
        following = []
        for C in current:
            buckets = {}
            for v in C:
                signature = tuple(sum(out[v][u] for u in D) for D in current)
                buckets.setdefault(signature, []).append(v)
            following.extend(tuple(buckets[key]) for key in sorted(buckets))
        following = tuple(following)
        if following == current:
            return current
        current = following


def main():
    p = 5
    regular = [[int(0 < (j-i) % p <= 2) for j in range(p)] for i in range(p)]
    aut = [q for q in permutations(range(p)) if all(
        regular[i][j] == regular[q[i]][q[j]] for i in range(p) for j in range(p))]
    need(len(aut) == 5, 'internal automorphism group')
    need(set(aut) == {tuple((i+a) % p for i in range(p)) for a in range(p)}, 'nonrotational automorphism')

    degrees = []
    total_relations = 0
    refinement_tests = 0
    two_subset_flags = Counter()
    displacement_flag_tests = 0
    symmetric_owner_map_tests = 0
    duplicate_example = (3, 3, 12, 24, 20)
    all_blocks = {}
    for d in range(1, p):
        blocks = list(regular_blocks(p, d))
        all_blocks[d] = blocks
        reps = {}
        hist = Counter()
        duplicates = 0
        expected_symmetric = {affine_pattern(p, s, set(S))
                              for s in range(1, p)
                              for S in combinations(range(p), d)}
        observed_symmetric = set()
        for B in blocks:
            H = stabilizer(B)
            need(len(H) in (1, 5), 'wrong stabilizer size')
            need(len({a for a, b in H}) == len(H) == len({b for a, b in H}), 'chart relation not functional')
            orbit = {shift(B, a, b) for a in range(p) for b in range(p)}
            key = min(orbit)
            reps[key] = len(H)
            need(len(H)*len(orbit) == 25, 'orbit-stabilizer')
            hist[len(H)] += 1
            if len(set(B)) < p:
                duplicates += 1
                need(len(H) == 1, 'duplicate-row block unexpectedly movable')
            if len(H) == p:
                observed_symmetric.add(B)
                s = next(b for a, b in H if a == 1)
                need(H == tuple((a, s*a % p) for a in range(p)), 'nonlinear chart gain')
                need(B == affine_pattern(p, s, {j for j in range(p) if B[0] >> j & 1}), 'wrong orbit pattern')
                q = chosen_matching(B)
                need(len(set(q)) == 5, 'symmetric block selector is not bijective')
                need(all(q[i] == (s*i+q[0]) % 5 for i in range(5)), 'selector did not retain actual stabilizer slope')
                opposite_bank = set(selector_bank(converse_block(B)))
                expected_bank = {tuple((-r[(-i) % 5]) % 5 for i in range(5)) for r in selector_bank(B)}
                need(opposite_bank == expected_bank, 'selected-owner bank not converse-closed')
                symmetric_owner_map_tests += 1
            if d in (2, 3):
                selected = B if d == 2 else tuple(31-row for row in B)
                transpose = tuple(sum(((selected[i] >> j) & 1) << i for i in range(5)) for j in range(5))
                def pair_flag(row):
                    u, v = [i for i in range(5) if row >> i & 1]
                    return int((u-v) % 5 in (1, 4))
                row_flags = tuple(map(pair_flag, selected))
                column_flags = tuple(map(pair_flag, transpose))
                row_constant = len(set(row_flags)) == 1
                column_constant = len(set(column_flags)) == 1
                need((row_constant and column_constant) == (len(H) == 5), 'pair flags miss rigid degree-two block')
                two_subset_flags[d, len(H), row_constant, column_constant] += 1
                if len(H) == 1:
                    split = [[i for i in range(5) if row_flags[i] == bit] for bit in range(2)]
                    split += [[5+j for j in range(5) if column_flags[j] == bit] for bit in range(2)]
                    split = [C for C in split if C]
                    need(all(len(C) == 1 for C in refine(tournament({(0, 1): B}), split)), 'intrinsic pair flag failed zero-owner refinement')
            if d in (1, 4):
                q = chosen_matching(B)
                word = tuple((q[(i+1) % 5]-q[i]) % 5 for i in range(5))
                need((len(set(word)) == 1) == (len(H) == 5), 'displacement flag missed rigid matching')
                if len(H) == 1:
                    split = [[i for i in range(5) if word[i] == value] for value in sorted(set(word))]
                    split.append(list(range(5, 10)))
                    need(all(len(C) == 1 for C in refine(tournament({(0, 1): B}), split)), 'displacement flag failed zero-owner refinement')
                displacement_flag_tests += 1
            A = tournament({(0, 1): B})
            for mark in range(5):
                partition = [(mark,), tuple(i for i in range(5) if i != mark), tuple(range(5, 10))]
                result = refine(A, partition)
                need(all(len(C) == 1 for C in result), 'one owner failed ordinary five-cell refinement')
                refinement_tests += 1
        need(observed_symmetric == expected_symmetric, 'slope-pattern classification')
        # Independently compare every local source/target orbit representative pair,
        # including empty relations between genuinely different attachment types.
        for B in reps:
            for C in reps:
                R = direct_relation(B, C)
                need(bool(R) == (B == C), 'ordinary orbit classifier')
                total_relations += 1
            for a, b in product(range(p), repeat=2):
                C = shift(B, a, b)
                R = direct_relation(B, C)
                H = stabilizer(B)
                expected = tuple(sorted(((a+x) % p, (b+y) % p) for x, y in H))
                need(R == expected, 'actual chart coset changed owners')
                total_relations += 1
        degrees.append({'degree': d, 'labelled_blocks': len(blocks),
                        'stabilizer_sizes': dict(sorted(hist.items())),
                        'ordinary_orbits': len(reps),
                        'orbits_by_stabilizer': dict(sorted(Counter(reps.values()).items())),
                        'blocks_with_repeated_rows': duplicates})
    need([r['labelled_blocks'] for r in degrees] == [120, 2040, 2040, 120], 'regular block counts')
    need([r['ordinary_orbits'] for r in degrees] == [8, 88, 88, 8], 'orbit counts')
    need(duplicate_example in all_blocks[2], 'repeated row witness missing')
    need(stabilizer(duplicate_example) == ((0, 0),), 'repeated row witness not rigid')

    # Full affine/singleton/empty triangle bank: each entry is a possible chart
    # relation, rather than an independent vertex-owner matching section.
    relation_bank = [()] + [((a, b),) for a, b in product(range(p), repeat=2)]
    relation_bank += [tuple((x, (s*x+t) % p) for x in range(p))
                      for s in range(1, p) for t in range(p)]
    need(len(relation_bank) == 46, 'relation alphabet')
    sections = Counter()
    for R, S, T in product(relation_bank, repeat=3):
        got = solve(3, {(0, 1): R, (1, 2): S, (0, 2): T})
        # Independent join, with no spanning-tree or propagation logic.
        truth = sorted((a, b, c) for a, b in R for bb, c in S
                       if b == bb and (a, c) in T)
        need(got == truth, 'shared triangle chart mismatch')
        need(len(got) in (0, 1, p), 'unexpected solution count')
        sections[len(got)] += 1

    I = affine_pattern(p, 1, {0})
    D = affine_pattern(p, 2, {0})
    U = affine_pattern(p, 1, {1})
    systems = {
        'balanced': {(0, 1): I, (1, 2): I, (0, 2): I},
        'nontrivial_gain': {(0, 1): I, (1, 2): I, (0, 2): D},
        'shifted_closing_block': {(0, 1): I, (1, 2): I, (0, 2): U},
        'repeated_row_pin': {(0, 1): duplicate_example, (1, 2): I, (0, 2): I},
    }
    fixtures = {}
    for name, B in systems.items():
        R = {e: direct_relation(v, v) for e, v in B.items()}
        solutions = solve(3, R)
        A = tournament(B)
        path = hamiltonian_path(A)
        K = complement_along_path(A, path)
        for solution in solutions:
            f = [5*i+(v+solution[i]) % 5 for i in range(3) for v in range(5)]
            need(all(A[i][j] == A[f[i]][f[j]] for i in range(15) for j in range(15)), 'chart solution not actual tournament automorphism')
            transported_path = [f[v] for v in path]
            transported_K = complement_along_path(A, transported_path)
            need(all(K[i][j] == transported_K[f[i]][f[j]] for i in range(15) for j in range(15)), 'literal path-complement transport lost ownership')
        fixtures[name] = {'blocks': {str(e): list(v) for e, v in B.items()},
                          'cell_degrees': [sum(A[5*i]) for i in range(3)],
                          'automorphism_shifts': solutions,
                          'extra_root_owners': int(len(solutions) > 1),
                          'actual_path': path}
        for mark in range(15):
            partition = [(mark,)] + [tuple(v for v in range(5*k, 5*k+5) if v != mark)
                                     for k in range(3)]
            need(all(len(C) == 1 for C in refine(A, partition)), 'one owner failed network refinement')
            refinement_tests += 1
    need(len(fixtures['balanced']['automorphism_shifts']) == 5, 'balanced symmetry lost')
    need(len(fixtures['nontrivial_gain']['automorphism_shifts']) == 1, 'cycle failed to rigidify')
    need(len(fixtures['repeated_row_pin']['automorphism_shifts']) == 1, 'repeated-row pin failed')
    cross_relations = {e: direct_relation(systems['balanced'][e], systems['shifted_closing_block'][e])
                       for e in systems['balanced']}
    need(all(len(R) == 5 for R in cross_relations.values()), 'local orbit hostile lost local feasibility')
    need(solve(3, cross_relations) == [], 'lost shared closing phase')
    # Every actual affine-matching triangle. Owner-section holonomy is a
    # different predicate from chart automorphisms; translations may be nonflat
    # while the full common rotation group remains genuine.
    affine_maps = [tuple((s*i+t) % 5 for i in range(5)) for s in range(1, 5) for t in range(5)]
    owner_triangle_types = Counter()
    for q01, q12, q02 in product(affine_maps, repeat=3):
        inverse02 = {y: x for x, y in enumerate(q02)}
        loop = tuple(inverse02[q12[q01[x]]] for x in range(5))
        gain = (loop[1]-loop[0]) % 5
        fixed = [x for x in range(5) if loop[x] == x]
        B = {(0, 1): tuple(1 << x for x in q01),
             (1, 2): tuple(1 << x for x in q12),
             (0, 2): tuple(1 << x for x in q02)}
        groups = solve(3, {e: stabilizer(v) for e, v in B.items()})
        if gain != 1:
            need(len(fixed) == 1 and len(groups) == 1, 'nontrivial gain must select one fixed owner and destroy symmetry')
            A = tournament(B)
            mark = fixed[0]
            partition = [(mark,), [x for x in range(5) if x != mark], list(range(5, 10)), list(range(10, 15))]
            need(all(len(C) == 1 for C in refine(A, partition)), 'intrinsic cycle fixed owner failed to split component')
            owner_triangle_types['rigid_unique_intrinsic_owner'] += 1
        else:
            need(len(groups) == 5, 'translation holonomy incorrectly destroyed chart symmetry')
            if fixed:
                need(len(fixed) == 5, 'identity owner loop fixed wrong number')
                owner_triangle_types['balanced_flat_five_owner_sections'] += 1
            else:
                owner_triangle_types['balanced_nonflat_no_owner_section_but_five_charts'] += 1
    # Converse is obtained by negating the internal coordinates and complementing
    # each block; it preserves every self-chart equation and solution count.
    for name, B in systems.items():
        C = {e: converse_block(v) for e, v in B.items()}
        got = solve(3, {e: direct_relation(v, v) for e, v in C.items()})
        wanted = sorted(tuple((-x) % p for x in v) for v in fixtures[name]['automorphism_shifts'])
        need(got == wanted, 'converse changed shared charts')
    edge = affine_pattern(5, 1, {0, 1})
    transposed = tuple(sum(((edge[i] >> j) & 1) << i for i in range(5)) for j in range(5))
    forward = chosen_matching(edge)
    fresh_reverse = chosen_matching(transposed)
    inverse = tuple(forward.index(i) for i in range(5))
    need(inverse != fresh_reverse, 'reverse-selection hostile vanished')
    need(tuple(fresh_reverse[forward[i]] for i in range(5)) == (4, 0, 1, 2, 3), 'wrong synthetic backtrack translation')

    return {'status': 'FINITE-EXACT', 'universe': 'all regular 5-by-5 blocks of degrees1..4, named R5 cells',
            'degree_census': degrees, 'total_labelled_blocks': sum(r['labelled_blocks'] for r in degrees),
            'total_ordinary_orbits': sum(r['ordinary_orbits'] for r in degrees),
            'direct_source_target_relation_tests': total_relations,
            'one_owner_ordinary_refinement_tests': refinement_tests,
            'displacement_flag_block_tests': displacement_flag_tests,
            'symmetric_intrinsic_matching_tests': symmetric_owner_map_tests,
            'two_subset_intrinsic_flag_census': [dict(degree=d, stabilizer_size=h,
                  row_flag_constant=r, column_flag_constant=c, blocks=count)
                  for (d, h, r, c), count in sorted(two_subset_flags.items())],
            'repeated_row_rigid_witness': list(duplicate_example),
            'triangle_relation_alphabet': len(relation_bank),
            'triangle_systems': len(relation_bank)**3,
            'triangle_counts_by_sections': dict(sorted(sections.items())),
            'actual_affine_matching_triangles': len(affine_maps)**3,
            'actual_owner_holonomy_types': dict(sorted(owner_triangle_types.items())),
            'independent_reverse_selection_hostile': {'block': list(edge), 'forward': forward,
                  'required_inverse': inverse, 'fresh_transpose_winner': fresh_reverse},
            'actual_tournament_fixtures': fixtures,
            'pairwise_feasible_but_globally_incompatible': {str(e): list(R) for e, R in cross_relations.items()},
            'checks': CHECKS}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    output = json.loads(json.dumps(main()))
    destination = HERE/'receipt.json'
    if args.write:
        destination.write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
    else:
        need(json.loads(destination.read_text(encoding='utf-8')) == output, 'receipt mismatch')
    print(json.dumps({'checks': output['checks'], 'blocks': output['total_labelled_blocks'],
                      'orbits': output['total_ordinary_orbits'],
                      'triangle_systems': output['triangle_systems'],
                      'triangle_counts': output['triangle_counts_by_sections']}, sort_keys=True))
