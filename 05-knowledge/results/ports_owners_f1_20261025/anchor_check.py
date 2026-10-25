"""Exact canonical rank-raising tournament forest, through order seven.

The imported helper generates isomorphism classes only. Anchor ownership,
rank, deletion, pair predicates and record reconstruction are implemented
locally. No general deterministic transducer is materialized here.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import comb
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
DEP = HERE.parent / 'metagraph_growth_20261020' / 'enumeration.py'
spec = spec_from_file_location('anchor_class_helper', DEP)
E = module_from_spec(spec)
spec.loader.exec_module(E)


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def arc(rows, u, v):
    return (rows[u] >> v) & 1


def induced(rows, vertices):
    return tuple(sum(arc(rows, u, v) << j for j, v in enumerate(vertices))
                 for u in vertices)


def transitive(rows, vertices):
    order = tuple(sorted(vertices, key=lambda u: -sum(arc(rows, u, v) for v in vertices)))
    return order if all(arc(rows, u, v) for u, v in combinations(order, 2)) else None


@lru_cache(None)
def rank(rows):
    n = len(rows)
    for size in range(n, -1, -1):
        if any(transitive(rows, vs) is not None for vs in combinations(range(n), size)):
            return n - size
    raise RuntimeError('empty set is transitive')


def triangle_rank(rows):
    n = len(rows)
    triangles = [set(vs) for vs in combinations(range(n), 3)
                 if all(sum(arc(rows, u, v) for v in vs) == 1 for u in vs)]
    for r in range(n + 1):
        if any(all(set(xs) & tri for tri in triangles)
               for xs in combinations(range(n), r)):
            return r
    raise RuntimeError('all vertices hit all triangles')


def encode(rows, anchors):
    base = transitive(rows, tuple(v for v in range(len(rows)) if v not in anchors))
    if base is None:
        return None
    header = sum(arc(rows, anchors[i], anchors[j]) << p
                 for p, (i, j) in enumerate(combinations(range(len(anchors)), 2)))
    word = tuple(sum(arc(rows, x, v) << i for i, x in enumerate(anchors)) for v in base)
    return (header, word), base


def decode(a, record):
    header, word = record
    n = a + len(word)
    rows = [0] * n
    for p, (i, j) in enumerate(combinations(range(a), 2)):
        u, v = (i, j) if (header >> p) & 1 else (j, i)
        rows[u] |= 1 << v
    for p, letter in enumerate(word, a):
        for i in range(a):
            u, v = (i, p) if (letter >> i) & 1 else (p, i)
            rows[u] |= 1 << v
        rows[p] |= sum(1 << q for q in range(p + 1, n))
    return tuple(rows)


@lru_cache(None)
def owner(rows, a):
    choices = [value[0] for anchors in permutations(range(len(rows)), a)
               if (value := encode(rows, anchors)) is not None]
    return min(choices) if choices else None


def project(a, record):
    """Remove the final anchor from an (a+1)-anchor record."""
    rows = decode(a + 1, record)
    old_pairs = list(combinations(range(a + 1), 2))
    header = sum(((record[0] >> old_pairs.index((i, j))) & 1) << p
                 for p, (i, j) in enumerate(combinations(range(a), 2)))
    word = tuple(letter & ((1 << a) - 1) for letter in record[1])
    need(decode(a, (header, word)) == induced(rows, tuple(range(a)) + tuple(range(a + 1, len(rows)))),
         'projection is literal vertex deletion')
    return header, word


def homogeneous_pairs(rows, vertices=None):
    vertices = tuple(range(len(rows))) if vertices is None else tuple(vertices)
    return {tuple(sorted((u, v))) for u, v in combinations(vertices, 2)
            if all(arc(rows, u, w) == arc(rows, v, w) for w in vertices if w not in (u, v))}


def disagreement(rows, u, v):
    return {w for w in range(len(rows)) if w not in (u, v)
            and arc(rows, u, w) != arc(rows, v, w)}


def chain_sizes(rows):
    neighbours = defaultdict(set)
    for u, v in homogeneous_pairs(rows):
        neighbours[u].add(v)
        neighbours[v].add(u)
    seen, result = set(), []
    for u in range(len(rows)):
        if u in seen:
            continue
        stack, comp = [u], set()
        while stack:
            v = stack.pop()
            if v in comp:
                continue
            comp.add(v)
            stack.extend(neighbours[v] - comp)
        need(transitive(rows, tuple(comp)) is not None, 'homogeneous pair component is transitive')
        seen.update(comp)
        result.append(len(comp))
    return sorted(result)


def pair_free_record(a, record):
    header, word = record
    if any(x == y for x, y in zip(word, word[1:])):
        return False
    anchors = decode(a, (header, ()))
    for i, j in combinations(range(a), 2):
        if all(arc(anchors, i, k) == arc(anchors, j, k) for k in range(a) if k not in (i, j)):
            if all(((b >> i) & 1) == ((b >> j) & 1) for b in word):
                return False
    for i in range(a):
        for p, b in enumerate(word):
            if (all(not ((c >> i) & 1) for c in word[:p]) and
                    all((c >> i) & 1 for c in word[p + 1:]) and
                    all(((b >> j) & 1) == 1 - arc(anchors, i, j)
                        for j in range(a) if j != i)):
                return False
    return True


def canonical_augmentation(rows):
    r = rank(rows)
    need(r > 0, 'augmentation requires positive feedback rank')
    a = r - 1
    candidates = []
    eligible = []
    for x in range(len(rows)):
        survivors = tuple(v for v in range(len(rows)) if v != x)
        parent = induced(rows, survivors)
        if rank(parent) != a:
            continue
        eligible.append(x)
        want = owner(parent, a)
        for aa in permutations(survivors, a):
            code = encode(rows, aa + (x,))
            if code is not None and project(a, code[0]) == want:
                candidates.append((code[0], x, aa, code[1]))
    need(candidates, 'every positive rank class has an admissible augmentation')
    best = min(c[0] for c in candidates)
    selected = [c for c in candidates if c[0] == best]
    parents = {project(a, c[0]) for c in selected}
    need(len(parents) == 1, 'automorphic witnesses have one parent record')
    parent_record = parents.pop()
    need(rank(decode(a, parent_record)) == a, 'rank drops exactly one')
    need(len(best[1]) == len(parent_record[1]), 'base word length preserved')
    need(E.canonical(decode(r, best)) == E.canonical(rows), 'all child arcs reconstructed')
    for code, x, aa, base in selected:
        need(decode(r, code) == induced(rows, aa + (x,) + base), 'literal child port permutation')
        need(decode(a, parent_record) == induced(rows, aa + base), 'literal parent port permutation')
        need(pair_free_record(r, code) == (not homogeneous_pairs(rows)), 'child clone filter exact')
        parent_rows = decode(a, parent_record)
        need(pair_free_record(a, parent_record) == (not homogeneous_pairs(parent_rows)), 'parent clone filter exact')
        survivors = tuple(v for v in range(len(rows)) if v != x)
        expected_pairs = {pair for pair in combinations(survivors, 2)
                          if disagreement(rows, *pair) <= {x}}
        need(homogeneous_pairs(rows, survivors) == expected_pairs, 'singleton-disagreement clone update')
    return best, parent_record, selected, eligible


def main():
    layers = {0: {()}, 1: {(0,)}}
    for n in range(2, 8):
        codes = {E.canonical(E.extension(rows, mask))
                 for rows in layers[n - 1] for mask in range(1 << (n - 1))}
        layers[n] = {E.decode(n, code) for code in codes}
    sizes = [len(layers[n]) for n in range(8)]
    need(sizes == [1, 1, 1, 2, 4, 12, 56, 456], 'complete class universe')
    rank_counts, forest, outdegrees = {}, {}, defaultdict(Counter)
    ordinary_owner_hostile = None
    clone_hostile = None
    child_count = 0
    normalizations, reverse_normalizations = {}, {}
    for n in range(8):
        rank_counts[n] = Counter()
        for rows in sorted(layers[n]):
            r = rank(rows)
            need(r == triangle_rank(rows), 'maximum-transitive versus triangle-hitting rank')
            alpha = n - r
            need(n <= (1 << alpha) - 1 if n else alpha == 0, 'elementary Ramsey bound')
            rank_counts[n][r] += 1
            if not r:
                continue
            child_count += 1
            child, parent, selected, eligible = canonical_augmentation(rows)
            p_rows = decode(r - 1, parent)
            p_code = E.canonical(p_rows)
            forest[(n, E.canonical(rows))] = (n - 1, p_code)
            outdegrees[(n - 1, r - 1)][p_code] += 1
            ordinary = owner(rows, r)
            g_key, l_key = (n, r, child), (n, r, ordinary)
            need(g_key not in normalizations and l_key not in reverse_normalizations,
                 'augmentation/ordinary normalization is a bijection on classes')
            normalizations[g_key] = l_key
            reverse_normalizations[l_key] = g_key
            if ordinary_owner_hostile is None and project(r - 1, ordinary) != owner(decode(r - 1, project(r - 1, ordinary)), r - 1):
                ordinary_owner_hostile = {'n': n, 'rank': r, 'rows': rows,
                    'ordinary_owner': ordinary, 'ordinary_parent_projection': project(r - 1, ordinary),
                    'augmentation_owner': child, 'canonical_parent': parent}
            if clone_hostile is None:
                for _, x, _, _ in selected:
                    survivors = [v for v in range(n) if v != x]
                    newpairs = homogeneous_pairs(rows, survivors) - homogeneous_pairs(rows)
                    if newpairs:
                        clone_hostile = {'rows': rows, 'deleted': x, 'new_parent_pairs': sorted(newpairs)}
                        break
    for (n, code), (pn, pcode) in forest.items():
        r = rank(E.decode(n, code))
        node = (n, code)
        depth = 0
        while node in forest:
            node = forest[node]
            depth += 1
        need(depth == r and node[0] == n - r, 'unique transitive root at fixed alpha')
    for n in range(1, 7):
        parent_code = E.canonical(tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n)))
        expected = (1 << n) - n - comb(n - 1, 2) - 1
        need(outdegrees[(n, 0)][parent_code] == expected, 'closed one-anchor first transition')
    histograms = {}
    for (n, r), counts in sorted(outdegrees.items()):
        parents = [rows for rows in layers[n] if rank(rows) == r]
        histogram = Counter(counts[E.canonical(rows)] for rows in parents)
        need(sum(k * v for k, v in histogram.items()) == rank_counts[n + 1][r + 1], 'canonical transition total')
        histograms[f'{n},{r}'] = dict(sorted(histogram.items()))
    need(ordinary_owner_hostile is not None, 'ordinary-owner projection hostile found')
    column = []
    ids = {}
    for n in range(3, 8):
        for rows in sorted(layers[n]):
            if n - rank(rows) != 3:
                continue
            key = (n, E.canonical(rows))
            ids[key] = f'{n}.{sum(v.startswith(str(n)+".") for v in ids.values()) + 1}'
            h = sum(all(arc(rows, u, v) for u, v in zip(p, p[1:]))
                    for p in permutations(range(n)))
            need(h == E.ham_count(rows), 'permutation HP count versus subset DP')
            column.append({'id': ids[key], 'parent': ids.get(forest.get(key)),
                           'rows': rows, 'H': h, 'clone_sizes': chain_sizes(rows)})
    need(len(column) == 9, 'complete alpha=3 column has nine classes')
    print('Class universe:', sizes)
    print('Feedback-rank rows:', {n: dict(sorted(c.items())) for n, c in rank_counts.items()})
    print('Canonical augmentation children:', child_count, 'PASS')
    print('Parent outdegree histograms (n,rank):', histograms)
    print('Ordinary-owner projection hostile:', json.dumps(ordinary_owner_hostile, sort_keys=True))
    print('Deletion-induced clone merger:', json.dumps(clone_hostile, sort_keys=True))
    print('Complete alpha=3 forest:', json.dumps(column, sort_keys=True))
    print('Imported class helper SHA256:', sha256(DEP.read_bytes().replace(b'\r\n', b'\n')).hexdigest())
    print('PASS: exact ranks, literal ports, canonical ownership, clone tests and forest roots')
    print('Scope: complete through n=7; general finite-state transducer proved separately, not materialized')


if __name__ == '__main__':
    main()
