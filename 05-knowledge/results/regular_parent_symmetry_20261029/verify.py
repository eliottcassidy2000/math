"""Exact native tests for regular-parent growth and residual symmetry laws."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial, prod
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
GATES = 0


def need(value, why):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(why)


def iso_maps(A, B, stop=False):
    n = len(A)
    if n != len(B):
        return []
    da, db = [v.bit_count() for v in A], [v.bit_count() for v in B]
    if sorted(da) != sorted(db):
        return []
    chosen, used, out = [], set(), []
    def visit():
        u = len(chosen)
        if u == n:
            out.append(tuple(chosen))
            return stop
        for v in range(n):
            if v in used or da[u] != db[v]:
                continue
            if any(((A[u] >> j) & 1) != ((B[v] >> chosen[j]) & 1) for j in range(u)):
                continue
            chosen.append(v); used.add(v)
            finished = visit()
            chosen.pop(); used.remove(v)
            if finished:
                return True
        return False
    visit()
    return out


def from_mask(n, mask):
    rows = [0]*n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if (mask >> bit) & 1:
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return rows


def regular_cards(rows):
    n = len(rows)
    if n % 2:
        return []
    m = (n-2)//2
    degrees = [r.bit_count() for r in rows]
    return [v for v in range(n) if all(degrees[u]-((rows[u] >> v) & 1) == m
                                      for u in range(n) if u != v)]


def add_vertex(rows, mask):
    n = len(rows)
    # mask is the set of OLD vertices beaten by the new vertex.
    return [row | (0 if (mask >> u) & 1 else 1 << n) for u, row in enumerate(rows)] + [mask]


def permute_mask(mask, p):
    return sum(((mask >> u) & 1) << p[u] for u in range(len(p)))


def carousel(n):
    return [sum(1 << v for v in range(n) if 0 < (v-u) % n <= (n-1)//2) for u in range(n)]


def paley7():
    return [sum(1 << v for v in range(7) if (v-u) % 7 in [1, 2, 4]) for u in range(7)]


def lex_product(A, B):
    a, b = len(A), len(B)
    return [sum(1 << (i*b+t) for t in range(b) if (B[s] >> t) & 1)
            | sum(((1 << b)-1) << (j*b) for j in range(a) if (A[i] >> j) & 1)
            for i in range(a) for s in range(b)]


def regular_representatives(n):
    """Generate all labelled regular tournaments by degree-pruned edge recursion."""
    edges = list(combinations(range(n), 2))
    target = (n-1)//2
    degrees, left, rows = [0]*n, [n-1]*n, [0]*n
    representatives, count = [], 0
    def visit(i):
        nonlocal count
        if i == len(edges):
            count += 1
            if not any(iso_maps(rows, R, True) for R in representatives):
                representatives.append(rows[:])
            return
        u, v = edges[i]
        left[u] -= 1; left[v] -= 1
        for a, b in [(u, v), (v, u)]:
            degrees[a] += 1
            rows[a] |= 1 << b
            if all(degrees[x] <= target <= degrees[x]+left[x] for x in [u, v]):
                visit(i+1)
            rows[a] ^= 1 << b
            degrees[a] -= 1
        left[u] += 1; left[v] += 1
    visit(0)
    return representatives, count


def extension_checks():
    total_card_histogram = {}
    for n in [2, 4, 6]:
        histogram = Counter()
        for mask in range(1 << (n*(n-1)//2)):
            A = from_mask(n, mask)
            cards = regular_cards(A)
            need(len(cards) <= 2, 'at most two regular cards')
            if len(cards) == 2:
                u, v = cards
                need(all(((A[w] >> u) & 1) == ((A[w] >> v) & 1)
                         for w in range(n) if w not in cards), 'two regular cards are a clone pair')
            histogram[len(cards)] += 1
        total_card_histogram[n] = dict(sorted(histogram.items()))
    records = []
    for n in [1, 3, 5, 7, 9]:
        if n == 9:
            parents, labelled = [lex_product(carousel(3), carousel(3))], None
        else:
            parents, labelled = regular_representatives(n)
        for parent in parents:
            need(all(row.bit_count() == (n-1)//2 for row in parent), 'parent is regular')
            group = iso_maps(parent, parent)
            cut_orbits, cut_stabs = {}, Counter()
            vertex_orbits, vertex_stabs = set(), Counter()
            for v in range(n):
                key = min(p[v] for p in group)
                if key not in vertex_orbits:
                    vertex_orbits.add(key)
                    vertex_stabs[sum(p[v] == v for p in group)] += 1
            children, child_stabs = [], Counter()
            rigid_labels = 0
            for mask in range(1 << n):
                key = min(permute_mask(mask, p) for p in group)
                stab = sum(permute_mask(mask, p) == mask for p in group)
                if key not in cut_orbits:
                    cut_orbits[key] = mask
                    cut_stabs[stab] += 1
                child = add_vertex(parent, mask)
                actual = iso_maps(child, child)
                need(all(p[n] == n for p in actual), 'every native child automorphism fixes new vertex')
                need(len(actual) == stab, 'full child group equals subset stabilizer')
                rigid_labels += int(stab == 1)
                if not any(iso_maps(child, old, True) for old in children):
                    children.append(child)
                    child_stabs[stab] += 1
            predicted = cut_stabs.copy()
            predicted.subtract(vertex_stabs)
            predicted = +predicted
            need(predicted == child_stabs, 'stabilizer-refined unmarked child count')
            need(len(children) == len(cut_orbits)-len(vertex_orbits), 'subset orbit minus vertex orbit')
            rejected_clones = {row | (1 << u) for u, row in enumerate(parent)}
            accepted_masks = [mask for mask in cut_orbits.values() if mask not in rejected_clones]
            need(len(accepted_masks) == len(children), 'clone-filter canonical augmentation count')
            need(Counter(sum(permute_mask(mask, p) == mask for p in group) for mask in accepted_masks)
                 == child_stabs, 'clone-filter canonical augmentation groups')
            records.append({'n': n, 'parent_rows': parent, 'parent_aut_order': len(group),
                            'regular_labelled_universe': labelled, 'subset_orbits': len(cut_orbits),
                            'vertex_orbits': len(vertex_orbits), 'child_classes': len(children),
                            'child_group_order_inventory': dict(sorted(child_stabs.items())),
                            'rigid_attachment_masks': rigid_labels, 'canonical_child_masks': accepted_masks})
    subclass = []
    for n in [1, 3, 5, 7]:
        inventory = Counter()
        for row in records:
            if row['n'] == n:
                inventory.update(row['child_group_order_inventory'])
        mass = sum(factorial(n+1)//s*c for s, c in inventory.items())
        labelled_regular = next(x['regular_labelled_universe'] for x in records if x['n'] == n)
        need(mass == (n+1)*(2**n-n)*labelled_regular, 'all-size labelled mass identity')
        exactly_two = (n+1)*n*labelled_regular
        exactly_one = (n+1)*(2**n-2*n)*labelled_regular
        if n+1 <= 6:
            hist = total_card_histogram[n+1]
            need(mass == hist.get(1, 0)+hist.get(2, 0), 'independent labelled regular-card mass')
            need(exactly_one == hist.get(1, 0) and exactly_two == hist.get(2, 0),
                 'all-size labelled deletion-multiplicity identities')
        subclass.append({'child_n': n+1, 'classes': sum(inventory.values()),
                         'group_order_inventory': dict(sorted(inventory.items())),
                         'labelled_mass': mass, 'exactly_one_regular_card': exactly_one,
                         'exactly_two_regular_cards': exactly_two})
    return {'all_even_labelled_through_6': total_card_histogram, 'parents': records,
            'regular_card_subclass': subclass}


def cycle_partition(p):
    seen, blocks = set(), []
    for u in range(len(p)):
        if u in seen:
            continue
        block, v = [], u
        while v not in seen:
            seen.add(v); block.append(v); v = p[v]
        blocks.append(tuple(sorted(block)))
    return tuple(sorted(blocks))


def max_transitive(rows):
    n = len(rows)
    for size in range(n, 0, -1):
        for subset in combinations(range(n), size):
            # A tournament is transitive iff its induced scores are 0,...,size-1.
            if sorted(sum((rows[u] >> v) & 1 for v in subset) for u in subset) == list(range(size)):
                return size


def seed_checks():
    result = []
    for name, rows in [('C3', carousel(3)), ('R5', carousel(5)), ('P7', paley7())]:
        need(not nontrivial_modules(rows), 'base substitution tournament is prime')
        n = len(rows)
        group = iso_maps(rows, rows)
        nontrivial = [p for p in group if p != tuple(range(n))]
        index = min(n-len(cycle_partition(p)) for p in nontrivial)
        minimal = {cycle_partition(p) for p in nontrivial if n-len(cycle_partition(p)) == index}
        generic = Counter()
        for partition in minimal:
            col = [next(j for j, block in enumerate(partition) if u in block) for u in range(n)]
            order = sum(all(col[p[u]] == col[u] for u in range(n)) for p in group)
            generic[order] += 1
        census = {}
        for q in [2, 3, 4]:
            counts = Counter()
            for col in product(range(q), repeat=n):
                counts[sum(all(col[p[u]] == col[u] for u in range(n)) for p in group)] += 1
            expected = ({1: q**3-q, 3: q} if name == 'C3' else
                        {1: q**5-q, 5: q} if name == 'R5' else
                        {1: q**7-7*q**3+6*q, 3: 7*(q**3-q), 21: q})
            need(dict(counts) == expected, 'complete seed stabilizer distribution')
            census[q] = dict(sorted(counts.items()))
        result.append({'seed': name, 'b': n, 'h': len(group), 'd': index,
                       'minimal_partition_groups': dict(sorted(generic.items())),
                       'max_transitive': max_transitive(rows), 'coloring_census': census})
    return result


def regular_orbits(name, q):
    return ((q**3-q)//3 if name == 'C3' else (q**5-q)//5 if name == 'R5'
            else (q**7-7*q**3+6*q)//21)


def interval(x, scale=10**12):
    lower = x.numerator*scale//x.denominator
    return {'lower': lower, 'upper': lower+1, 'denominator': scale}


def threshold_checks():
    out = []
    for name, b, h, d in [('C3', 3, 3, 2), ('R5', 5, 5, 4), ('P7', 7, 21, 4)]:
        for q in [2, 4, 8, 16, 32]:
            free, group_order = q, 1
            for k in range(1, 5):
                free = regular_orbits(name, free)
                group_order = h*group_order**b
                n = b**k
                probability = Fraction(group_order*free, q**n)
                seed_good = Fraction(h*regular_orbits(name, q), q**b)
                bottom = seed_good**(n//b)
                error = Fraction(0)
                for level in range(2, k+1):
                    m = b**(level-1)
                    g = h**((m-1)//(b-1))
                    error += Fraction(n, b*m)*(h-1)*Fraction(g, q**m)**d
                need(0 <= bottom-probability <= error, 'higher-level union bound vs exact all-depth recurrence')
                if name == 'C3':
                    explicit = Fraction(n, q**6) / (1-Fraction(3, q*q))
                    need(bottom-probability <= explicit, 'explicit triangle bottom coupling bound')
                out.append({'seed': name, 'q': q, 'depth': k, 'n': n,
                            'critical_parameter_n_over_qd': str(Fraction(n, q**d)),
                            'rigidity_probability_interval': interval(probability),
                            'bottom_probability_interval': interval(bottom),
                            'higher_symmetry_bound_interval': interval(min(error, Fraction(1)))})
    return out


def hostile_checks():
    result, triples = [], []
    for flip in [False, True]:
        rows = [0]*12
        for u, v in combinations(range(12), 2):
            b, s = divmod(u, 4); c, t = divmod(v, 4)
            arc = (c-b) % 3 == 1 if b != c else s == 0 or (s > 0 and (t-s) % 3 == 1)
            if flip and b == c == 2:
                arc = not arc
            rows[u if arc else v] |= 1 << (v if arc else u)
        group = iso_maps(rows, rows)
        need(len(group) == (27 if flip else 81), 'hostile actual group order')
        triples.append([t for t in combinations(range(12), 3)
                        if all(sum((rows[u] >> v) & 1 for v in t) == 1 for u in t)])
        regular = sum(all(permute_mask(mask, p) != mask for p in group[1:])
                      for mask in range(1 << 12))
        # iso_maps enumerates identity first, but check that dependency explicitly.
        need(group[0] == tuple(range(12)), 'identity is first group element')
        need(regular == (1728 if flip else 1620), 'binary regular masks distinguish hostile')
        result.append({'flip': flip, 'group_order': len(group), 'regular_binary_masks': regular})
    need(triples[0] == triples[1] and len(triples[0]) == 67, 'labelled cyclic triples coincide')
    return {'tournaments': result, 'common_cyclic_triples': 67}


def nontrivial_modules(rows):
    n = len(rows)
    modules = []
    for size in range(2, n):
        for subset in combinations(range(n), size):
            mask = sum(1 << u for u in subset)
            if all((rows[u] & mask) in [0, mask] for u in range(n) if u not in subset):
                modules.append(list(subset))
    return modules


def induced(rows, vertices):
    return [sum(((rows[u] >> v) & 1) << j for j, v in enumerate(vertices)) for u in vertices]


def forbidden_prime_checks():
    p = 11
    qr = {i*i % p for i in range(1, p)}
    rows = [sum(1 << v for v in range(p) if (v-u) % p in qr) for u in range(p)]
    need(not nontrivial_modules(rows), 'Paley eleven is prime')
    cards = []
    for v in range(p):
        card = induced(rows, [u for u in range(p) if u != v])
        mods = nontrivial_modules(card)
        need(not mods, 'every Paley-eleven card is prime of order ten')
        cards.append({'deleted': v, 'nontrivial_modules': len(mods)})
    seed = paley7()
    group = iso_maps(seed, seed)
    p7_cards = [len(nontrivial_modules(induced(seed, [u for u in range(7) if u != v])))
                for v in range(7)]
    need(p7_cards == [0]*7, 'all Paley-seven cards are prime of order six')
    need(not nontrivial_modules(seed), 'Paley-seven substitution is prime')
    need(sum(permute_mask(3, g) == 3 for g in group) == 1, 'cut {0,1} has trivial stabilizer')
    need(max_transitive(seed) == 3, 'Paley seven transitive maximum')
    need(max_transitive(induced(seed, [2, 3, 4])) == 3, 'zero-color quotient transitive triple')
    coloring, tower = [0], [0]
    construction = []
    for k in range(1, 5):
        coloring = [1-c if u in [0, 1] else c for u in range(7) for c in coloring]
        n = 7**k
        need(sum(coloring) == (n-3**k)//2, 'explicit cut weight recurrence')
        zero_transitive_words = list(product([2, 3, 4], repeat=k))
        indices = [sum(x*7**(k-1-j) for j, x in enumerate(word)) for word in zero_transitive_words]
        need(all(coloring[i] == 0 for i in indices), 'explicit all-zero transitive product')
        if k <= 2:
            tower = lex_product(seed, tower)
            mask = sum(c << i for i, c in enumerate(coloring))
            child = add_vertex(tower, mask)
            full = iso_maps(child, child)
            need(len(full) == 1, 'native full rigid Paley child at size 8 or 50')
            need(max_transitive(induced(child, indices+[n])) == 3**k+1,
                 'native transitive witness meets all-size upper bound')
        construction.append({'depth': k, 'parent_n': n, 'cut_weight': sum(coloring),
                             'complement_weight': n-sum(coloring),
                             'child_max_transitive': 3**k+1})
    c5 = carousel(5)
    r5_construction = []
    coloring, tower = [0], [0]
    for k in range(1, 5):
        coloring = [1-c if u == 0 else c for u in range(5) for c in coloring]
        n = 5**k
        need(sum(coloring) == (n-3**k)//2, 'carousel-five parity weight recurrence')
        if k <= 2:
            tower = lex_product(c5, tower)
            child = add_vertex(tower, sum(c << i for i, c in enumerate(coloring)))
            need(len(iso_maps(child, child)) == 1, 'native rigid carousel child at size 6 or 26')
            indices = [sum(x*5**(k-1-j) for j, x in enumerate(word))
                       for word in product([1, 2, 3], repeat=k)]
            need(all(coloring[i] == 0 for i in indices), 'carousel all-zero transitive product')
            need(max_transitive(induced(child, indices+[n])) == 3**k+1,
                 'carousel transitive witness attains all-size bound')
        r5_construction.append({'depth': k, 'child_n': n+1, 'cut_weight': sum(coloring),
                                'child_max_transitive': 3**k+1})
    c3 = carousel(3)
    triangle_tower = lex_product(c3, c3)
    triangle_colors = [int(u == 0) ^ int(v == 0) for u in range(3) for v in range(3)]
    triangle_child = add_vertex(triangle_tower, sum(c << i for i, c in enumerate(triangle_colors)))
    bad_witness = [0, 2, 4, 6, 9]
    need(bool(iso_maps(induced(triangle_child, bad_witness), c5)), 'rigid attachment creates forbidden R5')
    need(len(iso_maps(triangle_child, triangle_child)) == 1, 'guard hostile is genuinely rigid')
    need(bool(nontrivial_modules(induced(c5, [1, 2, 3, 4]))), 'R5 deletion fails prime guard')
    growth_vertices = [[0, 1, 3], [0, 1, 2, 3], list(range(5))]
    growth = [len(iso_maps(induced(c5, V), induced(c5, V))) for V in growth_vertices]
    need(growth == [3, 1, 5], 'two successive additions can introduce a new symmetry prime')
    return {'P11_rows': rows, 'P11_nontrivial_modules': 0, 'P11_cards': cards,
            'P7_card_nontrivial_module_counts': p7_cards,
            'rigid_P7_free_construction': r5_construction,
            'failed_guard_hostile': {'child_rows': triangle_child, 'induced_R5_vertices': bad_witness},
            'unprotected_two_step_growth_orders': growth,
            'rigid_P11_free_construction': construction}


def polynomial_product(a, b):
    out = Counter()
    for i, x in a.items():
        for j, y in b.items():
            out[i+j] += x*y
    return out


def ternary_inventory():
    polynomial = Counter({0: 2})
    out = []
    for k in range(5):
        n = 3**k
        point_exponent = (n-1)//2-k
        child = polynomial.copy()
        child[point_exponent] -= 1
        need(min(child.values()) >= 0, 'regular-parent polynomial subtraction is nonnegative')
        need(child[0] == (1 if k <= 1 else 2), 'rigid unmarked ternary children')
        out.append({'depth': k, 'parent_n': n, 'subset_orbits': sum(polynomial.values()),
                    'child_classes': sum(child.values()), 'point_stabilizer_exponent': point_exponent,
                    'child_stabilizer_exponent_inventory': dict(sorted((j, c) for j, c in child.items() if c))})
        cube = polynomial_product(polynomial_product(polynomial, polynomial), polynomial)
        for j, c in polynomial.items():
            cube[3*j] -= c
        need(all(c % 3 == 0 for c in cube.values()), 'nonconstant triple orbit division is integral')
        next_polynomial = Counter({j: c//3 for j, c in cube.items() if c})
        for j, c in polynomial.items():
            next_polynomial[3*j+1] += c
        polynomial = next_polynomial
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'certificate.json')
    args = parser.parse_args()
    result = {'extensions': extension_checks(), 'seeds': seed_checks(),
              'threshold_controls': threshold_checks(), 'hostile': hostile_checks(),
              'ternary_unmarked_children': ternary_inventory(),
              'forbidden_prime_guard': forbidden_prime_checks()}
    for row in result['ternary_unmarked_children'][:3]:
        parent = next(x for x in result['extensions']['parents'] if x['n'] == row['parent_n'])
        predicted = {3**j: c for j, c in row['child_stabilizer_exponent_inventory'].items()}
        need(predicted == parent['child_group_order_inventory'], 'polynomial vs native unmarked child groups')
    result['gates'] = GATES
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'gates': GATES, 'parent_cases': len(result['extensions']['parents']),
                      'output': str(args.output), 'sha256': hashlib.sha256(args.output.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
