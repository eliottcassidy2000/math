#!/usr/bin/env python3
"""Count-dipping source reset; complete fixed-order H105 replacement controls.

Universe: every unlabelled H105 source of orders8 and9, from the inherited
complete strong bank and SCC multiplication. Reducible nonsingleton factors
would have to be3*35, requiring at least10 vertices. No all-order optimality
claim is made for the observed replacement bottleneck.
"""
from collections import defaultdict, deque
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import heapq
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BANK = ROOT / '04-computation/experiments/h_spectrum_multiplicities_20261021/audit_S50/inventory_375.json'
TARGETS = HERE.parent / 'targets/receipt.json'

def require(q, message):
    if not q:
        raise RuntimeError(message)

def induced(rows, keep):
    return tuple(sum(((rows[u] >> v) & 1) << j for j, v in enumerate(keep)) for u in keep)

def decode(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if code >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)

def valid(rows):
    n = len(rows)
    return all(not (rows[i] >> i & 1) and rows[i] < 1 << n for i in range(n)) and all(
        ((rows[i] >> j & 1) + (rows[j] >> i & 1)) == 1 for i, j in combinations(range(n), 2))

@lru_cache(None)
def hp(rows):
    if not rows:
        return 1
    @lru_cache(None)
    def finish(last, unused):
        if not unused:
            return 1
        return sum(finish(v, unused ^ (1 << v)) for v in range(len(rows))
                   if unused >> v & 1 and rows[last] >> v & 1)
    return sum(finish(v, ((1 << len(rows))-1) ^ (1 << v)) for v in range(len(rows)))

def words(rows):
    def walk(p, remaining):
        if not remaining:
            yield p
        else:
            for v in remaining:
                if not p or rows[p[-1]] >> v & 1:
                    yield from walk(p + (v,), tuple(w for w in remaining if w != v))
    yield from walk((), tuple(range(len(rows))))

def code_at(rows, p):
    return sum(((rows[p[i]] >> p[j]) & 1) << k
               for k, (i, j) in enumerate(combinations(range(len(p)), 2)))

@lru_cache(None)
def canon(rows):
    # Restrict to Hamilton orders, a nonempty isomorphism-invariant set.
    return min((code_at(rows, p), p) for p in words(rows))

def join(a, b):
    return tuple(r | (((1 << len(b))-1) << len(a)) for r in a) + tuple(r << len(a) for r in b)

def chain(n):
    return tuple(sum(1 << j for j in range(i+1, n)) for i in range(n))

def birth(rows, cut):
    n = len(rows)
    return tuple(r | ((1-(cut >> i & 1)) << n) for i, r in enumerate(rows)) + (cut,)

def fixed_order(strong, n):
    result = []
    for key, core in sorted(strong.items()):
        if len(core) <= n:
            extra = n-len(core)
            for a in range(extra+1):
                rows = join(join(chain(a), core), chain(extra-a))
                result.append({'key': key, 'prefix': a, 'rows': rows})
    require(len({canon(t['rows'])[0] for t in result}) == len(result), 'no duplicate padded sources')
    return result

def replacement_graph(graphs):
    cards = defaultdict(list)
    for i, source in enumerate(graphs):
        rows = source['rows']
        for v in range(len(rows)):
            card = induced(rows, tuple(w for w in range(len(rows)) if w != v))
            code, _ = canon(card)
            cards[code].append((i, v, hp(card)))
    adj = defaultdict(list)
    for ls in cards.values():
        for (i, v, h), (j, w, k) in combinations(ls, 2):
            require(h == k, 'isomorphic cards have same H')
            if i != j:
                adj[i].append((j, h, v, w))
                adj[j].append((i, h, w, v))
    return adj

def bottleneck(graphs, adj, start_keys, end_keys):
    starts = [i for i, t in enumerate(graphs) if t['key'] in start_keys]
    ends = {i for i, t in enumerate(graphs) if t['key'] in end_keys}
    best = {i: 105 for i in starts}
    queue = [(-105, i) for i in starts]
    heapq.heapify(queue)
    pred = {i: None for i in starts}
    while queue:
        negative, i = heapq.heappop(queue)
        h = -negative
        if h != best[i]:
            continue
        if i in ends:
            path = []
            while pred[i] is not None:
                j, hcard, v, w = pred[i]
                path.append((j, i, hcard, v, w))
                i = j
            return h, list(reversed(path))
        for j, hcard, v, w in adj[i]:
            new = min(h, hcard)
            if new > best.get(j, -1):
                best[j] = new
                pred[j] = (i, hcard, v, w)
                heapq.heappush(queue, (-new, j))
    return None, []

def replacement_witness(left, right, v, w):
    left_keep = tuple(i for i in range(len(left)) if i != v)
    right_keep = tuple(i for i in range(len(right)) if i != w)
    lcard, rcard = induced(left, left_keep), induced(right, right_keep)
    lc, lp = canon(lcard)
    rc, rp = canon(rcard)
    require(lc == rc, 'common card')
    mapping = {left_keep[lp[i]]: right_keep[rp[i]] for i in range(len(lp))}
    require(all((left[a] >> b & 1) == (right[mapping[a]] >> mapping[b] & 1)
                for a in mapping for b in mapping), 'named surviving arcs')
    cut = sum((right[w] >> mapping[a] & 1) << i for i, a in enumerate(left_keep))
    actual = birth(lcard, cut)
    target_order = tuple(mapping[a] for a in left_keep) + (w,)
    require(actual == induced(right, target_order), 'actual new vertex reconstructs target')
    return {'left_rows': left, 'right_rows': right, 'delete_left': v, 'new_right': w,
            'survivor_map': sorted(mapping.items()), 'card_rows': lcard,
            'card_H': hp(lcard), 'new_vertex_cut_in_card_order': cut}

def literal(rows):
    return sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(len(p)-1))
               for p in permutations(range(len(rows))))

def main():
    bank = json.loads(BANK.read_text())
    strong = {}
    for ns, entries in bank['orders'].items():
        n = int(ns)
        for code, h, _, _ in entries:
            if h == 105:
                rows = decode(n, int(code))
                require(valid(rows) and hp(rows) == 105, 'actual strong source')
                c, p = canon(rows)
                strong[n, c] = induced(rows, p)
    require(len(strong) == 18, 'complete strong H105 bank')
    targets = json.loads(TARGETS.read_text())['families']['105']['kernels']
    branch_keys = []
    for N in (317, 319):
        branch_keys.append({tuple(t['key']) for t in targets
                            if t['branches'][str(N)]['all_core_padding_pairs']})
    require(not (branch_keys[0] & branch_keys[1]), 'inherited all-order branch disjointness')
    controls = []
    chosen = None
    for n, expected in ((8, 19), (9, 37)):
        graphs = fixed_order(strong, n)
        require(len(graphs) == expected, 'complete fixed-order source universe')
        adj = replacement_graph(graphs)
        floor, path = bottleneck(graphs, adj, *branch_keys)
        require(floor == 9 and len(path) == 2, 'finite replacement bottleneck and witness')
        controls.append({'order': n, 'classes': len(graphs), 'maximum_minimum_card_H': floor})
        if n == 8:
            chosen = graphs, path
    graphs, path = chosen
    witnesses = [replacement_witness(graphs[i]['rows'], graphs[j]['rows'], v, w)
                 for i, j, h, v, w in path]
    require([w['card_H'] for w in witnesses] == [19, 9], 'displayed reset trace')
    start = witnesses[0]['left_rows']
    end = witnesses[-1]['right_rows']
    branch_cuts = {}
    for N, rows in ((317, start), (319, end)):
        cuts = [cut for cut in range(1 << len(rows)) if hp(birth(rows, cut)) == N]
        require(bool(cuts), 'actual branch cut')
        branch_cuts[N] = cuts
        require(literal(birth(rows, cuts[0])) == N, 'independent literal branch count')
    all_named = {i: i for i in range(len(start))}
    for w in witnesses:
        mapping = dict(w['survivor_map'])
        all_named = {old: mapping[current] for old, current in all_named.items() if current in mapping}
        require(literal(w['card_rows']) == w['card_H'], 'independent literal card count')
        require(literal(w['right_rows']) == 105, 'independent literal source count')
    require(len(all_named) == 6, 'six original vertices survive the two replacements')
    ports = (1, 2)  # Both start inside the original strong H105 component.
    joint = []
    current_ports = ports
    for i, rows in enumerate([start] + [w['right_rows'] for w in witnesses]):
        if i:
            transport = dict(witnesses[i-1]['survivor_map'])
            current_ports = tuple(transport[p] for p in current_ports)
        cube = [hp(induced(rows, tuple(v for v in range(len(rows))
                   if v not in {current_ports[j] for j in range(2) if mask >> j & 1})))
                for mask in range(4)]
        joint.append({'ports': current_ports, 'deletion_cube': cube,
                      'mod49': [v % 49 for v in cube], 'mod343': [v % 343 for v in cube]})
    require(joint[0]['deletion_cube'] != joint[-1]['deletion_cube'],
            'retained names do not imply preserved joint response')
    result = {'status': 'PASS; explicit reset and finite exact fixed-order bottlenecks',
              'strong_H105_sources': len(strong), 'fixed_order_controls': controls,
              'count_sequence': [105] + [h for w in witnesses for h in (w['card_H'], 105)],
              'replacements': witnesses,
              'branch_cuts': branch_cuts, 'original_to_final_survivors': sorted(all_named.items()),
              'retained_named_pairs': list(combinations(sorted(all_named), 2)),
              'named_pair_joint_deletion_controls': joint,
              'input_sha256': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                               for p in (BANK, TARGETS)},
              'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
