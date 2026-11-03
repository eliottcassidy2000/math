#!/usr/bin/env python3
"""Small actual H31 branch and neutral-native certificates; no producer imports.

The only carrier input is the pinned universal two-port receipt. All named
graph operations, path counts, inverse ancestry and terminal programs below
are reconstructed independently. No global support census is performed.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import permutations
from pathlib import Path
from math import gcd
import hashlib
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
INPUT = RESULTS / 'adaptive_state_induction_20261103/two_port_envelope.json'
PIN = 'fc14da2ead29fa9c9fe5218ad96c87bce606c0cf41f5f0d5a4325917ebaeeafe'
C3 = (2, 4, 1)


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def validate(T):
    n = len(T)
    need(all(0 <= row < 1 << n and not row >> i & 1 for i, row in enumerate(T)), 'row bounds')
    need(all(((T[i] >> j) ^ (T[j] >> i)) & 1 for i in range(n) for j in range(i)), 'tournament arcs')


@lru_cache(None)
def H(T):
    n = len(T)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for S in range(1, 1 << n):
        for last, value in enumerate(dp[S]):
            if value:
                choices = T[last] & ~S
                while choices:
                    bit = choices & -choices
                    choices -= bit
                    dp[S | bit][bit.bit_length() - 1] += value
    return sum(dp[-1])


def paths(T):
    n = len(T)
    result = []
    def visit(word, left):
        if not left:
            result.append(word)
            return
        choices = left if not word else left & T[word[-1]]
        while choices:
            bit = choices & -choices
            choices -= bit
            visit(word + (bit.bit_length() - 1,), left ^ bit)
    visit((), (1 << n) - 1)
    need(len(result) == H(T), 'independent path DFS versus subset DP')
    return tuple(result)


def permutation_count(T):
    return sum(all(T[a] >> b & 1 for a, b in zip(p, p[1:]))
               for p in permutations(range(len(T))))


def ear(T, incoming):
    n = len(T)
    return tuple(row | ((1 << n) if incoming >> i & 1 else 0)
                 for i, row in enumerate(T)) + (((1 << n) - 1) ^ incoming,)


def birth(T, root):
    return ear(T, 1 << root)


def apply(T, ports, word):
    ports = list(ports)
    for code in word:
        mask, drop = divmod(code, 3)
        n = len(T)
        T = ear(T, sum(1 << ports[i] for i in range(2) if mask >> i & 1))
        if drop < 2:
            ports[drop] = n
    return T, tuple(ports)


def relabel(T, names):
    U = [0] * len(T)
    for i, row in enumerate(T):
        U[names[i]] = sum(1 << names[j] for j in range(len(T)) if row >> j & 1)
    return tuple(U)


def native(T, p):
    need(sorted(p) == list(range(len(T))) and
         all(T[a] >> b & 1 for a, b in zip(p, p[1:])), 'legal native path')
    kept = {frozenset(pair) for pair in zip(p, p[1:])}
    U = list(T)
    for i in range(len(T)):
        for j in range(i):
            if frozenset((i, j)) not in kept:
                U[i] ^= 1 << j
                U[j] ^= 1 << i
    U = tuple(U)
    validate(U)
    return U


def strong(T, vertices=None):
    V = tuple(range(len(T))) if vertices is None else tuple(vertices)
    mask = sum(1 << v for v in V)
    for start in V:
        reached = 1 << start
        while True:
            old = reached
            for v in V:
                if reached >> v & 1:
                    reached |= T[v] & mask
            if old == reached:
                break
        if reached != mask:
            return False
    return bool(V)


def ancestry(T):
    validate(T)
    if not strong(T):
        return {'strong': False, 'terminal_programs': []}
    V, inverse = list(range(len(T))), []
    while len(V) > 3:
        candidates = []
        for x in V:
            incoming = [v for v in V if T[v] >> x & 1]
            if len(incoming) == 1 and strong(T, [v for v in V if v != x]):
                candidates.append((x, incoming[0]))
        need(len(candidates) <= 1, 'unique strong inverse')
        if not candidates:
            break
        x, root = candidates[0]
        V.remove(x)
        inverse.append((x, root))
    births = tuple(reversed(inverse))
    out = {'strong': True, 'core_order': len(V), 'core_vertices': V,
           'births': births, 'terminal_programs': []}
    if len(V) != 3:
        return out
    zero = births[0][1] if births else V[0]
    one = next(v for v in V if T[zero] >> v & 1)
    two = next(v for v in V if v not in (zero, one))
    names = (zero, one, two) + tuple(x for x, _ in births)
    index = {v: i for i, v in enumerate(names)}
    history = tuple(index[r] for _, r in births)
    states = {(u, v): ((), (u, v)) for u in V for v in V if T[u] >> v & 1}
    for x, root in births:
        nxt = {}
        for ports, (word, initial) in states.items():
            if root not in ports:
                continue
            mask = 1 << ports.index(root)
            for drop in range(3):
                pp = list(ports)
                if drop < 2:
                    pp[drop] = x
                nxt.setdefault(tuple(pp), (word + (3 * mask + drop,), initial))
        states = nxt
    rebuilt = C3
    for root in history:
        rebuilt = birth(rebuilt, root)
    need(relabel(rebuilt, names) == T, 'literal inverse-ancestry reconstruction')
    for ports, (word, initial) in states.items():
        replay, pair = apply(C3, tuple(index[v] for v in initial), word)
        need(relabel(replay, names) == T and tuple(names[v] for v in pair) == ports,
             'literal terminal program and port transport')
    out.update(history=history, canonical_names=names,
               terminal_programs=[{'ports': ports, 'word': word, 'initial_ports': initial}
                                  for ports, (word, initial) in sorted(states.items())])
    return out


def isomorphisms(T, U):
    if len(T) != len(U):
        return []
    return [p for p in permutations(range(len(T)))
            if all(T[i].bit_count() == U[p[i]].bit_count() for i in range(len(T)))
            and all(bool(T[i] >> j & 1) == bool(U[p[i]] >> p[j] & 1)
                    for i in range(len(T)) for j in range(i))]


def advance(v, M):
    w = [0] * len(v)
    for i, j, c in M:
        w[j] += c * v[i]
    return tuple(w)


def output(v, O):
    return sum(x * y for x, y in zip(v, O))


def vecword(v, word, matrices):
    for c in word:
        v = advance(v, matrices[c])
    return v


def closure(seed, alphabet, cap, matrices, observer):
    vectors, words, ids, edges = [seed], [()], {seed: 0}, []
    for v, word in zip(vectors, words):
        row = []
        for c in alphabet:
            u = advance(v, matrices[c])
            h = output(u, observer)
            need(h > output(v, observer), 'proper-cut strict growth on actual source')
            if h > cap:
                continue
            if u not in ids:
                ids[u] = len(vectors)
                vectors.append(u)
                words.append(word + (c,))
            row.append((ids[u], c))
        edges.append(row)
    return vectors, words, edges


def support(vectors, words, edges, targets, observer):
    counts = [output(v, observer) for v in vectors]
    masks, witnesses = {}, {}
    for i in sorted(range(len(vectors)), key=lambda i: counts[i], reverse=True):
        mask = 0
        for j, target in enumerate(targets):
            if counts[i] == target:
                mask |= 1 << j
                witnesses[i, j] = ()
            else:
                for child, c in edges[i]:
                    if masks[child] >> j & 1:
                        mask |= 1 << j
                        witnesses[i, j] = (c,) + witnesses[child, j]
                        break
        masks[i] = mask
    return masks, witnesses


def free_roots(T, cap):
    nodes = [(T, ())]
    targets = {}
    for U, word in nodes:
        h = H(U)
        if h in (93, 95, 97):
            targets.setdefault(h, word)
        for root in range(len(U)):
            V = birth(U, root)
            need(H(V) > h, 'all-root growth')
            if H(V) <= cap:
                nodes.append((V, word + (root,)))
    return {'nodes': len(nodes), 'maximum_births': max(len(w) for _, w in nodes),
            'target_root_words': targets}


def main():
    need(hashlib.sha256(INPUT.read_bytes()).hexdigest() == PIN, 'carrier input pin')
    data = json.loads(INPUT.read_text())
    matrices, observer = data['matrices'], data['observer']
    triangle = advance(tuple(data['seed']), matrices[8])
    need(output(triangle, observer) == 3, 'TT2 to triangle embedding')
    old_letters, new_letters = tuple(range(3, 9)), tuple(range(3, 12))

    # All actual adaptive ordinary classes H<=31, with first C3 root fixed0.
    graphs = [(C3, ())]
    for T, history in graphs:
        for root in ((0,) if not history else range(len(T))):
            U = birth(T, root)
            if H(U) <= 31:
                graphs.append((U, history + (root,)))
    fibre_graphs = [(T, word) for T, word in graphs if H(T) == 31]
    need(len(graphs) == 31 and len(fibre_graphs) == 3, 'complete small ordinary universe')
    fibre = []
    source_classes = []
    for T, history in fibre_graphs:
        a = ancestry(T)
        need(tuple(a['history']) == history, 'actual canonical history')
        need(len(isomorphisms(T, T)) == 1, 'source rigidity')
        source_classes.append({'rows': T, 'history': history, 'ancestry': a,
                               'free_root_descendants': free_roots(T, 97)})
        for terminal in a['terminal_programs']:
            v = vecword(triangle, terminal['word'], matrices)
            need(output(v, observer) == H(T), 'terminal literal source count')
            need(apply(C3, terminal['initial_ports'], terminal['word']) == (T, terminal['ports']),
                 'chronological target source replay')
            fibre.append({'rows': T, 'ports': terminal['ports'], 'source_program': terminal,
                          'history': history, 'vector': v})
    need(len(fibre) == len({row['vector'] for row in fibre}) == 37, '37 actual marked observers')
    old_vectors, old_words, _ = closure(triangle, old_letters, 31, matrices, observer)
    need(len(old_vectors) == 318 and {v for v in old_vectors if output(v, observer) == 31}
         == {row['vector'] for row in fibre}, 'independent source-image agreement')

    histograms = {}
    for name, alphabet in (('six_adaptive', old_letters), ('nine_proper_cuts', new_letters)):
        hist = Counter()
        for row in fibre:
            vv, ww, ee = closure(row['vector'], alphabet, 97, matrices, observer)
            masks, witnesses = support(vv, ww, ee, (93, 95, 97), observer)
            hist[masks[0]] += 1
            row[name] = {'mask': masks[0], 'descendant_states': len(vv),
                         'maximum_births': max(map(len, ww)),
                         'branches': {str(h): witnesses[0, i] for i, h in enumerate((93, 95, 97))
                                      if masks[0] >> i & 1}}
        histograms[name] = dict(sorted(hist.items()))
    need(histograms['six_adaptive'] == {0: 6, 2: 3, 4: 10, 5: 18}, 'inherited fibre histogram')
    need(histograms['nine_proper_cuts'] == {0: 6, 2: 3, 4: 10, 5: 16, 7: 2}, 'mask3 exact enlargement')

    # The two mask7 observers are opposite ordered slots of one actual source.
    good = [row for row in fibre if row['nine_proper_cuts']['mask'] == 7]
    need({tuple(row['history']) for row in good} == {(0, 3, 0, 5)} and
         {tuple(row['ports']) for row in good} == {(3, 5), (5, 3)}, 'one physical branching pair')
    source_word = (4, 8, 3, 5)
    T, ports = apply(C3, (0, 1), source_word)
    v = vecword(triangle, source_word, matrices)
    selected = {}
    for h, word in ((93, (3, 3)), (95, (10,)), (97, (6,))):
        U, pp = apply(T, ports, word)
        a = ancestry(U)
        need(H(U) == h and strong(U), 'literal first branch')
        marked = [p for p in a['terminal_programs'] if tuple(p['ports']) == pp]
        need(marked, 'supplied marked reentry')
        if len(U) == 8:
            need(permutation_count(U) == h, '8! literal branch control')
        selected[h] = {'word': word, 'rows': U, 'ports': pp, 'H': h,
                       'ancestry': a, 'marked_reentry': marked[0]}
    bad95, badports = apply(T, ports, (9,))
    need(bad95 == selected[95]['rows'] and
         not any(tuple(p['ports']) == badports for p in ancestry(bad95)['terminal_programs']),
         'same H95 graph but wrong retained pair fails marked reentry')

    vv, ww, ee = closure(v, new_letters, 295, matrices, observer)
    need(len(vv) == 195 and max(map(len, ww)) == 4, 'targeted next-level universe')
    next_level = {}
    for h in (93, 95, 97):
        masks, witnesses = support(vv, ww, ee, (3*h, 3*h+2, 3*h+4), observer)
        ids = [i for i, value in enumerate(vv) if output(value, observer) == h]
        histogram = dict(sorted(Counter(masks[i] for i in ids).items()))
        need(all(masks[i] not in (6, 7) for i in ids), 'no next simultaneous positive branches')
        next_level[h] = {'states': len(ids), 'mask_histogram': histogram,
                         'sources': [{'source_word': ww[i], 'mask': masks[i],
                                      'branches': {str(3*h+2*j): witnesses[i, j] for j in range(3)
                                                   if masks[i] >> j & 1}} for i in ids]}
    literal_dag = []
    for i, (value, word) in enumerate(zip(vv, ww)):
        U, pp = apply(T, ports, word)
        need(H(U) == output(value, observer) and strong(U), 'all195 actual graph recounts')
        literal_dag.append({'word': word, 'H': H(U), 'rows': U, 'ports': pp,
                            'vector': value, 'edges': ee[i]})

    # An H-neutral native reset at the six-vertex H31 source.
    A = fibre_graphs[0][0]
    need(A == (10, 36, 1, 22, 7, 29), 'named native source')
    fixed_ports = (5, 1)
    path = (0, 3, 4, 1, 5, 2)
    U = native(A, path)
    need(U == (60, 41, 26, 48, 34, 4) and H(U) == 31, 'named neutral endpoint')
    before95, pp95 = apply(A, fixed_ports, (3, 3))
    after97, pp97 = apply(U, fixed_ports, (7, 7))
    need(H(before95) == 95 and H(after97) == 97, 'source-preserving named branching')
    need(permutation_count(before95) == 95 and permutation_count(after97) == 97,
         'independent8! native-branch controls')
    native_graphs, native_ids, native_words, native_edges = [A], {A: 0}, [()], []
    for V, history in zip(native_graphs, native_words):
        row = []
        for p in paths(V):
            W = native(V, p)
            if H(W) != 31:
                continue
            if W not in native_ids:
                native_ids[W] = len(native_graphs)
                native_graphs.append(W)
                native_words.append(history + (p,))
            row.append({'to': native_ids[W], 'path': p})
        native_edges.append(row)
    need(len(native_graphs) == 6 and all(len(e) == 2 for e in native_edges), 'exact named neutral hexagon')
    native_records = []
    for i, V in enumerate(native_graphs):
        ia, iu = isomorphisms(A, V), isomorphisms(U, V)
        need(len(ia) + len(iu) == 1, 'exact two rigid ordinary classes')
        native_records.append({'rows': V, 'ancestry': ancestry(V),
                               'ordinary_class': 'A' if ia else 'U',
                               'isomorphism_from_class_representative': (ia or iu)[0],
                               'path_from_A': native_words[i], 'edges': native_edges[i]})
    source_native = []
    for V, history in fibre_graphs:
        hh, neutral, self_count = Counter(), [], 0
        for p in paths(V):
            W = native(V, p)
            hh[H(W)] += 1
            if H(W) == 31:
                neutral.append({'path': p, 'rows': W, 'ancestry': ancestry(W)})
                self_count += bool(isomorphisms(V, W))
        need(self_count == 0, 'allthree class-staying groups trivial')
        source_native.append({'history': history, 'native_H_histogram': dict(sorted(hh.items())),
                              'neutral_edges': neutral, 'self_native_generators': self_count})
    stabilized = []
    for incoming in (0, (1 << len(A)) - 1):
        V = ear(A, incoming)
        histogram = Counter(H(native(V, p)) for p in paths(V))
        need(H(V) == 31 and 31 not in histogram, 'one source/sink stabilization has no neutral first native step')
        stabilized.append({'incoming': incoming, 'rows': V, 'native_H_histogram': dict(sorted(histogram.items()))})

    def queries(V):
        return tuple(H(apply(V, fixed_ports, tuple(3*mask+2 for mask in word))[0])
                     for word in data['rank']['suffixes'])
    X_A, X_U = queries(A), queries(U)
    need(X_A != X_U and gcd(*X_A) == gcd(*X_U) == 1, 'same H/content loses branch response')
    result = {
        'status': 'PASS', 'scope': 'H31 actual source/branch tests; no global spectrum or closure claim',
        'source_sha256': {str(INPUT.relative_to(RESULTS)): PIN},
        'command': '3*incoming_mask+discard; positive polarity; masks1,2 are adaptive',
        'all_adaptive_H_le_31': {'ordinary_classes': len(graphs), 'maximum_order': max(len(T) for T, _ in graphs),
                               'carrier_states': len(old_vectors), 'H31_sources': source_classes},
        'H31_marked_fibre': fibre, 'H31_branch_histograms': histograms,
        'mask3_first_level': {'source_rows': T, 'ports': ports, 'source_word': source_word,
                              'branches': selected, 'wrong_95_pair': badports},
        'nine_command_next_level': {'states': len(vv), 'maximum_births': max(map(len, ww)),
                                    'fibres': next_level, 'literal_DAG': literal_dag},
        'native_fixed_pair_repair': {'source': A, 'ports': fixed_ports, 'path': path, 'target': U,
                                     'source_queries': X_A, 'target_queries': X_U,
                                     'branch95': {'word': (3, 3), 'rows': before95, 'ports': pp95},
                                     'branch97_after_native': {'word': (7, 7), 'rows': after97, 'ports': pp97},
                                     'free_root_target_tree': free_roots(U, 97)},
        'named_native_component': native_records, 'three_source_native_controls': source_native,
        'one_stabilization_hostiles': stabilized,
    }
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
