"""Independent target-core audit: no producer imports or path-word canonicalizer."""
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
TARGETS = HERE.parent/'targets'


def need(x, why):
    if not x:
        raise RuntimeError(why)


def decode(n, code):
    a = [0]*n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if code >> k & 1 else (j, i)
        a[u] |= 1 << v
    return tuple(a)


def restrict(a, vertices):
    return tuple(sum((1 << j) for j, w in enumerate(vertices) if a[v] >> w & 1)
                 for v in vertices)


def counts(a):
    n = len(a)
    incoming = tuple(sum((1 << j) for j in range(n) if a[j] >> i & 1) for i in range(n))
    @lru_cache(None)
    def ending(mask, last):
        old = mask ^ (1 << last)
        if not old:
            return 1
        return sum(ending(old, j) for j in range(n) if old & incoming[last] & (1 << j))
    def h(mask):
        return 1 if not mask else sum(ending(mask, j) for j in range(n) if mask >> j & 1)
    full = (1 << n)-1
    return h(full), tuple(h(full ^ (1 << j)) for j in range(n))


def sccs(a):
    # Independent forward searches rather than producer's transitive closure.
    reachable = []
    for start in range(len(a)):
        seen, todo = {start}, [start]
        while todo:
            v = todo.pop()
            for w in range(len(a)):
                if a[v] >> w & 1 and w not in seen:
                    seen.add(w)
                    todo.append(w)
        reachable.append(seen)
    parts = {tuple(j for j in range(len(a)) if j in reachable[i] and i in reachable[j])
             for i in range(len(a))}
    return sorted(parts, key=lambda p: -sum(bool(a[p[0]] >> q[0] & 1) for q in parts if p != q))


def core(a):
    parts = sccs(a)
    active = [i for i, p in enumerate(parts) if len(p) > 1]
    need(active, 'source nontransitive')
    lo, hi = active[0], active[-1]
    vertices = tuple(v for part in parts[lo:hi+1] for v in part)
    return restrict(a, vertices), lo, len(parts)-1-hi


@lru_cache(None)
def key(a):
    # Degree-ordered canonical form: all permutations inside each degree block.
    # Isomorphisms preserve degrees, so this is complete, not just a signature.
    groups = defaultdict(list)
    for i, r in enumerate(a):
        groups[r.bit_count()].append(i)
    orderings = product(*(permutations(groups[d]) for d in sorted(groups)))
    return len(a), min(sum(((a[p[i]] >> p[j]) & 1) << k
                          for k, (i, j) in enumerate(combinations(range(len(a)), 2)))
                       for blocks in orderings for p in [sum(blocks, ())])


def join(a, b):
    n = len(a)
    return tuple(r | (((1 << len(b))-1) << n) for r in a)+tuple(r << n for r in b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='TARGET-AUDIT.json')
    args = parser.parse_args()
    report = json.loads((TARGETS/'receipt.json').read_text())
    bank_path = ROOT/report['bank']
    need(sha256(bank_path.read_bytes()).hexdigest() == report['bank_sha256'], 'bank pin')
    bank = json.loads(bank_path.read_text())
    pairs = {(55, 111), (55, 167), (55, 169), (105, 317), (105, 319), (35, 107), (35, 109)}
    targets = {b for a, b in pairs}
    found = defaultdict(lambda: defaultdict(set))
    classes = deletions = hits = 0
    for ns, rows in bank['orders'].items():
        n = int(ns)
        for code, h, *_ in rows:
            if h not in targets:
                continue
            a = decode(n, int(code))
            actual_h, ds = counts(a)
            need(actual_h == h and len(sccs(a)) == 1, 'selected strong child')
            classes += 1
            for root, old_h in enumerate(ds):
                deletions += 1
                if (old_h, h) not in pairs:
                    continue
                source = restrict(a, tuple(v for v in range(n) if v != root))
                middle, left, right = core(source)
                found[old_h, h][key(middle)].add((left, right))
                hits += 1
    for old_h, family in report['families'].items():
        h = int(old_h)
        for target in family['targets']:
            supplied = {}
            for k in family['kernels']:
                pads = k['branches'][str(target)]['all_core_padding_pairs']
                if pads:
                    supplied[key(tuple(k['rows']))] = set(map(tuple, pads))
            need(supplied == found[h, target], ('exact cores and every padding pair', h, target))
    intersection55 = set(found[55, 167]) & set(found[55, 169])
    intersection105 = set(found[105, 317]) & set(found[105, 319])
    need(len(intersection55) == 2 and not intersection105, 'common-source verdicts')
    need(set(found[55, 111]) == intersection55, '111 branch same source cores')
    need(not found[35, 107] and not found[35, 109], 'negative controls')
    # Internal singleton must survive; exterior singleton must disappear.
    c3, q4 = (2, 4, 1), (10, 4, 1, 6)
    internal = join(join(c3, (0,)), q4)
    middle, left, right = core(internal)
    need(len(middle) == 8 and (left, right) == (0, 0), 'retain internal singleton')
    padded = join(join(join((0,), (0,)), internal), (0,))
    padded_middle, left, right = core(padded)
    need(key(padded_middle) == key(internal) and (left, right) == (2, 1), 'strip exterior only')
    # Explicit common-source cuts from proof, checked by recursive predecessor DP.
    controls = []
    examples = [((70, 12, 8, 17, 39, 79, 30), (119, 5), (167, 169)),
                ((78, 68, 8, 18, 39, 79, 28), (16, 31), (167, 169))]
    for a, cuts, hs in examples:
        need(counts(a)[0] == 55, 'literal source55')
        for cut, h in zip(cuts, hs):
            n = len(a)
            child = tuple(r | (0 if cut >> i & 1 else 1 << n) for i, r in enumerate(a))+(cut,)
            need(counts(child)[0] == h, 'literal alternative branch')
            controls.append(h)
    receipt = {'status': 'PASS; independent source-kernel and padding audit',
               'selected_children': classes, 'vertex_deletions': deletions, 'relevant_deletions': hits,
               'branch_core_counts': {f'{a}->{b}': len(found[a, b]) for a, b in sorted(pairs)},
               'common55': len(intersection55), 'common105': len(intersection105),
               'literal_branch_controls': controls, 'internal_singleton_control': 'PASS',
               'canonicalization': 'minimum over all within-degree permutations; no producer imports',
               'bank_sha256': report['bank_sha256']}
    (HERE/args.output).write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
