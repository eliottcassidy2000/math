"""Exact controls for the marked-matrix / spanning-tree / subset-fill bridge."""
from collections import Counter
from itertools import permutations, product, combinations
from pathlib import Path
import hashlib
import json
import math
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def mul(a, b):
    c = Counter()
    for s, x in a.items():
        for t, y in b.items():
            if not s & t:
                c[s | t] += x * y
    return {s: x for s, x in c.items() if x}

def det(matrix):
    n = len(matrix)
    answer = Counter()
    for p in permutations(range(n)):
        term = {0: (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))}
        for i in range(n):
            term = mul(term, matrix[i][p[i]])
        answer.update(term)
    return {s: x for s, x in answer.items() if x}

def path_table(rows):
    n = len(rows)
    table = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        table[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            rest = mask ^ (1 << v)
            if rest:
                table[mask][v] = sum(table[rest][u] for u in range(n) if rows[u] >> v & 1)
    return [1] + [sum(table[s]) for s in range(1, 1 << n)]

def literal_endpoints(rows):
    n = len(rows)
    out = [[0] * n for _ in range(n)]
    for p in permutations(range(n)):
        if all(rows[p[i]] >> p[i + 1] & 1 for i in range(n - 1)):
            out[p[0]][p[-1]] += 1
    return out

def marked_endpoints(rows):
    # Square-zero matrix multiplication, retaining the start mark.
    n = len(rows)
    result = []
    for r in range(n):
        state = [{1 << r: 1} if v == r else {} for v in range(n)]
        for _ in range(n - 1):
            nxt = [Counter() for _ in range(n)]
            for v in range(n):
                for w in range(n):
                    if rows[v] >> w & 1:
                        nxt[w].update(mul(state[v], {1 << w: 1}))
            state = nxt
        result.append([state[s].get((1 << n) - 1, 0) for s in range(n)])
    return result

def tree_endpoints(rows):
    # Edge i->j has weight x_i.  Principal cofactor is the weighted out-tree sum.
    n = len(rows)
    matrix = [[{} for _ in range(n)] for _ in range(n)]
    for i in range(n):
        matrix[i][i] = {1 << u: 1 for u in range(n) if rows[u] >> i & 1}
        for j in range(n):
            if rows[i] >> j & 1:
                matrix[i][j] = {1 << i: -1}
    result = []
    for r in range(n):
        keep = [v for v in range(n) if v != r]
        polynomial = det([[matrix[i][j] for j in keep] for i in keep])
        result.append([polynomial.get(((1 << n) - 1) ^ (1 << s), 0) for s in range(n)])
    return result

def forest_coefficients(rows):
    # Principal minors, graded by deleted root-set size, expand det(zI+L).
    n = len(rows)
    matrix = [[{} for _ in range(n)] for _ in range(n)]
    for i in range(n):
        matrix[i][i] = {1 << u: 1 for u in range(n) if rows[u] >> i & 1}
        for j in range(n):
            if rows[i] >> j & 1:
                matrix[i][j] = {1 << i: -1}
    result = [0] * (n + 1)
    for roots in range(1 << n):
        keep = [v for v in range(n) if not roots >> v & 1]
        polynomial = det([[matrix[i][j] for j in keep] for i in keep])
        result[roots.bit_count()] += sum(polynomial.values())
    return result

def literal_cover_coefficients(rows):
    n = len(rows)
    ordered = [0] * (n + 1)
    for p in permutations(range(n)):
        for cuts in range(1 << (n - 1)):
            if all(cuts >> i & 1 or rows[p[i]] >> p[i + 1] & 1 for i in range(n - 1)):
                ordered[1 + cuts.bit_count()] += 1
    need(all(ordered[k] % math.factorial(k) == 0 for k in range(n + 1)), 'free component ordering')
    return [ordered[k] // math.factorial(k) for k in range(n + 1)]

def direct_bin_count(rows, bins):
    h = path_table(rows)
    total = 0
    for assignment in product(range(bins), repeat=len(rows)):
        masks = [0] * bins
        for v, color in enumerate(assignment):
            masks[color] |= 1 << v
        total += math.prod(h[mask] for mask in masks)
    return total

def odd_cycles(rows):
    n = len(rows)
    cycles = []
    for size in range(3, n + 1, 2):
        for vertices in combinations(range(n), size):
            r = vertices[0]
            for tail in permutations(vertices[1:]):
                p = (r,) + tail
                if all(rows[p[i]] >> p[(i + 1) % size] & 1 for i in range(size)):
                    cycles.append(sum(1 << v for v in vertices))
    return cycles

def fill_from_cycles(n, cycles):
    p = {0: 1}
    for v in range(n):
        p = mul(p, {0: 1, 1 << v: 1})
    for cycle in cycles:
        p = mul(p, {0: 1, cycle: 2})
    return p

def tournament(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a, b = (i, j) if code >> k & 1 else (j, i)
        rows[a] |= 1 << b
    return rows

def digraph(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate((i, j) for i in range(n) for j in range(n) if i != j):
        if code >> k & 1:
            rows[i] |= 1 << j
    return rows

def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]

def power(rows, exponent):
    n = len(rows)
    a = [[int(rows[i] >> j & 1) for j in range(n)] for i in range(n)]
    p = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(exponent):
        p = matmul(p, a)
    return p

def brute_trees_by_leaves(rows):
    n = len(rows)
    counts = Counter()
    for root in range(n):
        vertices = [v for v in range(n) if v != root]
        menus = [[u for u in range(n) if rows[u] >> v & 1] for v in vertices]
        for choices in product(*menus):
            parent = dict(zip(vertices, choices))
            valid = True
            for v in vertices:
                seen = set()
                while v != root and v not in seen:
                    seen.add(v)
                    v = parent[v]
                if v != root:
                    valid = False
                    break
            if valid:
                degrees = Counter(choices)
                counts[sum(degrees[v] == 0 for v in range(n))] += 1
    return dict(sorted(counts.items()))

digraph_count = 0
for n in range(1, 5):
    for code in range(1 << (n * (n - 1))):
        rows = digraph(n, code)
        literal = literal_endpoints(rows)
        need(literal == marked_endpoints(rows), ('marked walk', n, code))
        need(literal == tree_endpoints(rows), ('marked tree', n, code))
        covers = forest_coefficients(rows)
        need(covers == literal_cover_coefficients(rows), ('marked forest', n, code))
        for bins in range(4):
            from_forest = sum(math.prod(range(bins - k + 1, bins + 1)) * covers[k]
                              for k in range(min(bins, n) + 1))
            need(from_forest == direct_bin_count(rows, bins), ('falling-factorial bin law', n, code, bins))
        digraph_count += 1

tournament_count = 0
subset_cells = 0
birth_checks = 0
spectra = {}
for n in range(1, 6):
    values = Counter()
    for code in range(1 << (n * (n - 1) // 2)):
        rows = tournament(n, code)
        h = path_table(rows)
        cycles = odd_cycles(rows)
        p = fill_from_cycles(n, cycles)
        need([p[s] for s in range(1 << n)] == h, ('all subset fill', n, code))
        need(all(x > 0 and x % 2 for x in h), ('positive odd cells', n, code))
        if n >= 2:
            old = rows[:-1]
            old = [row & ((1 << (n - 1)) - 1) for row in old]
            oldp = fill_from_cycles(n - 1, odd_cycles(old))
            newfactor = {0: 1, 1 << (n - 1): 1}
            for c in cycles:
                if c >> (n - 1) & 1:
                    newfactor[c] = newfactor.get(c, 0) + 2
            need(mul(oldp, newfactor) == p, ('exact birth factor', n, code))
            birth_checks += 1
        tournament_count += 1
        subset_cells += 1 << n
        values[h[-1]] += 1
    spectra[n] = dict(sorted(values.items()))

# Directed books: 0->1->B->0, with B transitive.
books = []
for m in range(1, 9):
    n = m + 2
    rows = [2, ((1 << n) - 1) ^ 3] + [1 | sum(1 << j for j in range(i + 1, n)) for i in range(2, n)]
    need(path_table(rows)[-1] == 2 ** m + 1, ('book formula', m))
    books.append({'pages': m, 'rows': rows, 'H': 2 ** m + 1})
fake = fill_from_cycles(5, [7, 11, 19])
need(fake[31] == 7 and all(v > 0 and v % 2 for v in fake.values()), 'hostile fake fill')
book = books[2]['rows']
need(Counter(odd_cycles(book)) == Counter({7: 1, 11: 1, 19: 1, 31: 1}), 'forced five cycle')

# A missing arc destroys7/21 exclusion already at very small orders.
missing = {}
for n in (5, 6):
    for code in range(1 << (n * (n - 1) // 2)):
        if len(missing) == 2:
            break
        rows = tournament(n, code)
        for i in range(n):
            for j in range(n):
                if rows[i] >> j & 1:
                    changed = rows.copy()
                    changed[i] ^= 1 << j
                    value = path_table(changed)[-1]
                    if value in (7, 21) and value not in missing:
                        need(sum(map(sum, literal_endpoints(changed))) == value, 'missing arc independent')
                        opposite = changed.copy()
                        opposite[j] |= 1 << i
                        missing[value] = {'n': n, 'tournament_rows': rows, 'deleted_arc': [i, j], 'rows': changed,
                                          'H_before': path_table(rows)[-1], 'H_after': value,
                                          'H_opposite_completion': path_table(opposite)[-1]}
    if len(missing) == 2:
        break
need(set(missing) == {7, 21}, 'missing-pair witnesses')

tt3 = [6, 4, 0]
star4 = [14, 0, 0, 0]
cycle3_padded = [2, 4, 1, 0, 0, 0]
need(sum(map(sum, power(tt3, 3))) == 0 and path_table(tt3)[-1] == 1, 'nilpotent hostile')
need(sum(map(sum, power(tt3, 2))) == 1, 'correct path exponent')
need(path_table(star4)[-1] == 0 and brute_trees_by_leaves(star4) == {3: 1}, 'tree not path')
need(sum(power(cycle3_padded, 6)[i][i] for i in range(6)) == 3 and path_table(cycle3_padded)[-1] == 0,
     'closed walk not spanning cycle')

sources = [
    '01-canon/theorems/THM-002-ocf.md',
    '01-canon/theorems/THM-029-alpha1-gap.md',
    '01-canon/theorems/THM-1460-arborescences-are-the-determinantal-relaxation-of-H.md',
    '01-canon/theorems/THM-1745-leaf-graded-arborescence-filtration-721-shadow.md',
    '01-canon/theorems/THM-1870-the-cycle-count-spectral-boundary-is-the-hamiltonian-length.md',
    '01-canon/theorems/THM-4986-the-7-21-omission-proof-has-an-exact-spanning-tree-twin-block-floor-crevasse-lemma-and-the-first-creator-contrast.md',
    '05-knowledge/hypotheses/HYP-9640-chan-pak-as-a-matrix-semigroup-orbit-and-one-defect-digraph-spectra.md',
    '05-knowledge/results/squarefree_cycle_inventory_20261019/SQUAREFREE_CYCLE_BRIDGES.md',
]
report = {
    'status': 'FINITE-EXACT controls of stated all-order identities; no coverage proof',
    'all_loopless_labelled_digraphs_through_order4': digraph_count,
    'weighted_forest_and_bin_controls': {'graphs': digraph_count, 'bin_values': [0, 1, 2, 3]},
    'all_labelled_tournaments_through_order5': tournament_count,
    'subset_cells': subset_cells, 'birth_factor_checks': birth_checks,
    'labelled_H_histograms': spectra,
    'directed_books': books,
    'fake_three_page_fill': {'cycle_supports': [7, 11, 19], 'all32_cells_positive_odd': True, 'H': fake[31],
                             'real_completion_H': path_table(book)[-1]},
    'one_missing_arc_witnesses': missing,
    'matrix_hostiles': {'TT3_A3_zero_H1': True, 'three_star_one_tree_H0': True,
                        'padded_C3_trace_A6_3_H0': True},
    'source_pins': [{'path': p, 'sha256': hashlib.sha256((ROOT / p).read_bytes()).hexdigest()} for p in sources],
}
suffix = '.optimized' if sys.flags.optimize else ''
(HERE / ('receipt' + suffix + '.json')).write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps({k: report[k] for k in ('all_loopless_labelled_digraphs_through_order4',
      'all_labelled_tournaments_through_order5', 'subset_cells', 'birth_factor_checks', 'one_missing_arc_witnesses')}))
