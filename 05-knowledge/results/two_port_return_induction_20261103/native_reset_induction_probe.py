#!/usr/bin/env python3
"""Literal named native returns, unique ancestry, and two-port queries.

No repository imports. All construction clocks and ports are physical names.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import permutations, product
from math import gcd, isqrt
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def validate(T):
    n = len(T)
    need(all(not (T[i] >> i & 1) and T[i] < (1 << n) for i in range(n)), 'rows')
    need(all(((T[i] >> j) ^ (T[j] >> i)) & 1 for i in range(n) for j in range(i)), 'tournament')


@lru_cache(None)
def H(T):
    n = len(T)
    D = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        D[1 << i][i] = 1
    for S in range(1, 1 << n):
        for v, count in enumerate(D[S]):
            if count:
                rest = T[v] & ~S
                while rest:
                    bit = rest & -rest
                    rest -= bit
                    D[S | bit][bit.bit_length() - 1] += count
    return sum(D[-1])


def paths(T):
    return tuple(p for p in permutations(range(len(T)))
                 if all(T[a] >> b & 1 for a, b in zip(p, p[1:])))


def native(T, p):
    need(sorted(p) == list(range(len(T))) and all(T[a] >> b & 1 for a, b in zip(p, p[1:])), 'native guard')
    edges = {frozenset((a, b)) for a, b in zip(p, p[1:])}
    R = list(T)
    for i in range(len(T)):
        for j in range(i):
            if frozenset((i, j)) not in edges:
                R[i] ^= 1 << j
                R[j] ^= 1 << i
    return tuple(R)


def ear(T, incoming):
    n = len(T)
    return tuple(row | ((1 << n) if incoming >> i & 1 else 0)
                 for i, row in enumerate(T)) + (((1 << n) - 1) ^ incoming,)


def birth(T, root):
    return ear(T, 1 << root)


def strong(T, vertices=None):
    V = tuple(range(len(T))) if vertices is None else tuple(vertices)
    if not V:
        return False
    target = sum(1 << v for v in V)
    for start in V:
        reached = 1 << start
        while True:
            old = reached
            for v in V:
                if reached >> v & 1:
                    reached |= T[v] & target
            if reached == old:
                break
        if reached != target:
            return False
    return True


def ancestry(T):
    validate(T)
    if not strong(T):
        return {'strong': False, 'in_two_port_language': False}
    V = list(range(len(T)))
    inverse = []
    while len(V) > 3:
        candidates = []
        for x in V:
            ins = [u for u in V if T[u] >> x & 1]
            if len(ins) == 1 and strong(T, [u for u in V if u != x]):
                candidates.append((x, ins[0]))
        need(len(candidates) <= 1, 'unique strong adaptive inverse')
        if not candidates:
            break
        x, root = candidates[0]
        V.remove(x)
        inverse.append((x, root))
    births = tuple(reversed(inverse))
    result = {'strong': True, 'core_vertices': V, 'core_order': len(V),
              'births': births, 'in_two_port_language': False}
    if len(V) != 3:
        return result
    order = tuple(V) + tuple(x for x, _ in births)
    chron = {v: i for i, v in enumerate(order)}
    root_word = tuple(chron[r] for _, r in births)
    width = max((len({r for r in root_word[t - 3:] if r < t})
                 for t in range(3, len(T))), default=0)
    states = {(u, v): ([], (u, v)) for u in V for v in V if T[u] >> v & 1}
    for x, root in births:
        new = {}
        for ports, (program, initial_ports) in states.items():
            if root not in ports:
                continue
            slot = ports.index(root)
            for drop in range(3):
                pp = list(ports)
                if drop < 2:
                    pp[drop] = x
                new.setdefault(tuple(pp), (program + [3 * slot + drop], initial_ports))
        states = new
    need(bool(states) == (width <= 2), 'independent width iff versus literal slot programs')
    for final, (program, initial) in states.items():
        replay = list(initial)
        for (x, root), code in zip(births, program):
            slot, drop = divmod(code, 3)
            need(replay[slot] == root, 'named program root replay')
            if drop < 2:
                replay[drop] = x
        need(tuple(replay) == final, 'named program terminal replay')
    result.update(root_word=root_word, chronological_names=order, width=width,
                  in_two_port_language=bool(states),
                  terminal_programs=[{'ports': p, 'commands': w, 'initial_ports': initial}
                                     for p, (w, initial) in sorted(states.items())])
    return result


QUERIES = ((), (1,), (2,), (3,), (0, 1), (1, 1), (2, 1), (3, 1),
           (0, 2), (1, 2), (2, 2), (3, 2), (2, 0, 1), (3, 0, 1),
           (0, 2, 1), (2, 2, 1), (3, 2, 1))


def query(T, ports, word):
    for mask in word:
        incoming = sum(1 << ports[i] for i in range(2) if mask >> i & 1)
        T = ear(T, incoming)
    return H(T)


def X(T, ports):
    return tuple(query(T, ports, word) for word in QUERIES)


def paired_queries(T, U, ports):
    before, after = X(T, ports), X(U, ports)
    return {'ports': ports, 'before': before, 'after': after,
            'gcd_before': gcd(*before), 'gcd_after': gcd(*after),
            'adaptive_child_counts_before': [H(birth(T, v)) for v in ports],
            'adaptive_child_counts_after': [H(birth(U, v)) for v in ports]}


def graph_record(T):
    return {'rows': T, 'H': H(T), 'ancestry': ancestry(T)}


def history_graph(word):
    T = (2, 4, 1)
    for root in word:
        T = birth(T, root)
    return T


def ordinary_history_key(T, a=None):
    a = ancestry(T) if a is None else a
    if a.get('core_order') != 3:
        return None
    births = a['births']
    if not births:
        return ()
    zero = births[0][1]
    core = a['core_vertices']
    one = next(v for v in core if T[zero] >> v & 1)
    two = next(v for v in core if v not in (zero, one))
    mapping = {v: i for i, v in enumerate((zero, one, two) + tuple(x for x, _ in births))}
    return tuple(mapping[root] for _, root in births)


def canonical_names(T, a=None):
    a = ancestry(T) if a is None else a
    births = a['births']
    need(a.get('core_order') == 3 and births, 'rigid C3 descendant')
    zero = births[0][1]
    one = next(v for v in a['core_vertices'] if T[zero] >> v & 1)
    two = next(v for v in a['core_vertices'] if v not in (zero, one))
    return (zero, one, two) + tuple(x for x, _ in births)


def procedural_sources(n):
    if n == 3:
        return [((), (2, 4, 1))]
    out = []
    for suffix in product(*(range(i) for i in range(4, n))):
        word = (0,) + suffix
        T = history_graph(word)
        if ancestry(T)['in_two_port_language']:
            out.append((word, T))
    return out


def neutral_procedural():
    counts = {}
    witnesses = []
    for n in range(4, 8):
        sources = procedural_sources(n)
        counts[n] = {'width_two_classes': len(sources), 'native_paths': 0,
                     'neutral_paths': 0, 'neutral_to_width_two': 0, 'nonisomorphic_targets': 0}
        seen = set()
        for word, T in sources:
            h = H(T)
            a = ancestry(T)
            source_ports = {tuple(row['ports']) for row in a['terminal_programs']}
            pp = paths(T)
            need(len(pp) == h, 'source H independent literal paths')
            counts[n]['native_paths'] += len(pp)
            for p in pp:
                U = native(T, p)
                if H(U) != h:
                    continue
                counts[n]['neutral_paths'] += 1
                b = ancestry(U)
                if not b['in_two_port_language']:
                    continue
                counts[n]['neutral_to_width_two'] += 1
                key = ordinary_history_key(U, b)
                if key == word or (word, key) in seen:
                    continue
                seen.add((word, key))
                counts[n]['nonisomorphic_targets'] += 1
                target_ports = {tuple(row['ports']) for row in b['terminal_programs']}
                common = sorted(source_ports & target_ports)
                witnesses.append({'source_word': word, 'source': graph_record(T), 'path': p,
                                  'endpoint_word': key, 'endpoint': graph_record(U),
                                  'common_legal_ordered_ports': common,
                                  'chosen_queries': None if not common else paired_queries(T, U, common[0])})
    return {'scope': 'all width-two C3 root histories through seven vertices, first root fixed0',
            'counts': counts, 'nonisomorphic_witnesses': witnesses}


def known_returns():
    controls = [
        ('H93', (110, 28, 104, 112, 37, 2, 50),
         ((0, 2, 5, 1, 3, 6, 4), (6, 4, 3, 0, 2, 5, 1)), 33, 283),
        ('H95', (110, 92, 104, 48, 37, 66, 24),
         ((1, 2, 3, 4, 0, 5, 6), (3, 1, 2, 4, 6, 0, 5)), 44, 287),
        ('H37', (40, 17, 3, 6, 13, 30),
         ((4, 2, 0, 5, 3, 1), (2, 0, 5, 3, 4, 1),
          (2, 0, 5, 4, 1, 3), (0, 1, 5, 4, 2, 3)), None, None),
    ]
    rows = []
    for label, T, pp, outgoing, target in controls:
        stages = [graph_record(T)]
        U = T
        for p in pp:
            U = native(U, p)
            stages.append(graph_record(U))
        need(H(T) == H(U), 'return H')
        row = {'label': label, 'paths': pp, 'stages': stages,
               'query_pairs': [paired_queries(T, U, ports)
                               for ports in ((0, 1), (2, 4))]}
        if target is not None:
            incoming = ((1 << len(U)) - 1) ^ outgoing
            child = ear(U, incoming)
            need(H(child) == target == len(paths(child)), 'literal repaired ear')
            row['repaired_ear'] = {'outgoing_mask': outgoing, 'target': target,
                                   'incoming_old_names': [i for i in range(len(U)) if incoming >> i & 1],
                                   'child': graph_record(child)}
        rows.append(row)
    return rows


def constant_to_chain():
    rows = []
    outputs = [3, 5, 9]
    for m in range(3, 11):
        outputs.append(sum(outputs[-3:]))
    for m in range(1, 11):
        T = history_graph((0,) * m)
        n = len(T)
        p = (3, 2, 0) + tuple(range(n - 1, 3, -1)) + (1,)
        U = native(T, p)
        Q = (0,) + tuple(range(n - 1, 3, -1)) + (1, 3, 2)
        need(all(U[Q[i]] >> Q[(i+1) % 3] & 1 for i in range(3)), 'mapped seed cycle')
        need(all(bool(U[Q[j]] >> Q[i] & 1) == (j == i-1)
                 for i in range(3, n) for j in range(i)), 'literal all-pair newborn-chain map')
        a, b = ancestry(T), ancestry(U)
        ports = (2, 0)
        need(all(any(tuple(x['ports']) == ports for x in aa['terminal_programs']) for aa in (a, b)),
             'both named endpoints reenter')
        need(H(T) == 2 ** (m+1) + 1 and H(U) == outputs[m], 'two exact growth laws')
        row = {'m': m, 'source_rows': T, 'path': p, 'endpoint_rows': U,
               'endpoint_chronological_names': Q, 'source_H': H(T), 'endpoint_H': H(U),
               'source_program': [x for x in a['terminal_programs'] if tuple(x['ports']) == ports][0],
               'endpoint_program': [x for x in b['terminal_programs'] if tuple(x['ports']) == ports][0]}
        if m <= 4:
            row['queries'] = paired_queries(T, U, ports)
        rows.append(row)
    relation = [rows[3]['queries']['before'][i] - 6*rows[2]['queries']['before'][i]
                + 11*rows[1]['queries']['before'][i] - 6*rows[0]['queries']['before'][i] for i in range(17)]
    defect = rows[3]['endpoint_H'] - 6*rows[2]['endpoint_H'] + 11*rows[1]['endpoint_H'] - 6*rows[0]['endpoint_H']
    need(relation == [0]*17 and defect == -2, 'actual four-source linear obstruction')
    T, U = tuple(rows[1]['source_rows']), tuple(rows[1]['endpoint_rows'])
    need([H(birth(T, v)) for v in (2, 0)] == [23, 17], 'H9 old menu')
    need([H(birth(U, v)) for v in (2, 0)] == [17, 31], 'H9 reset menu')
    child = birth(U, 0)
    need(H(child) == len(paths(child)) == 31, 'literal H31 repair')
    need(ancestry(child)['in_two_port_language'], 'repaired child language')
    return {'rows': rows, 'query_vector_relation': relation, 'native_output_defect': defect,
            'H9_repair_child': graph_record(child)}


def module_controls():
    counts = {'adaptive_sources': 0, 'proper_nonsingleton_modules': 0}
    for n in range(3, 8):
        words = [()] if n == 3 else [(0,) + s for s in product(*(range(i) for i in range(4, n)))]
        for word in words:
            T = history_graph(word)
            counts['adaptive_sources'] += 1
            for mask in range(1, (1 << n)-1):
                V = [i for i in range(n) if mask >> i & 1]
                if len(V) < 2 or any((T[v] & mask) not in (0, mask) for v in range(n) if not mask >> v & 1):
                    continue
                counts['proper_nonsingleton_modules'] += 1
                need(all(not (T[x] >> y & 1 and T[y] >> z & 1 and T[z] >> x & 1)
                         for x, y, z in permutations(V, 3)), 'every proper module transitive')
    return counts


def self_native_group(T):
    n = len(T)
    key = ordinary_history_key(T)
    source_names = canonical_names(T)
    generators = []
    for p in paths(T):
        U = native(T, p)
        a = ancestry(U)
        if ordinary_history_key(U, a) != key:
            continue
        target_names = canonical_names(U, a)
        f = [None] * n
        for u, v in zip(source_names, target_names):
            f[u] = v
        f = tuple(f)
        need(all(bool(T[i] >> j & 1) == bool(U[f[i]] >> f[j] & 1)
                 for i in range(n) for j in range(n)), 'literal unique isomorphism')
        generators.append({'path': p, 'permutation': f})
    group = {tuple(range(n)): []}
    todo = list(group)
    while todo:
        f = todo.pop()
        for i, row in enumerate(generators):
            g = row['permutation']
            fg = tuple(f[g[v]] for v in range(n))
            if fg not in group:
                group[fg] = group[f] + [i]
                todo.append(fg)
    pairs = {(i, j) for i in range(n) for j in range(n) if i != j}
    orbits = []
    while pairs:
        i, j = min(pairs)
        orbit = sorted({(f[i], f[j]) for f in group})
        orbits.append(orbit)
        pairs.difference_update(orbit)
    return {'source': graph_record(T), 'generators': generators, 'group_order': len(group),
            'elements': [{'permutation': f, 'generator_word': w} for f, w in sorted(group.items())],
            'ordered_pair_orbits': orbits}


def factor_escape():
    T = history_graph((0, 3, 3))
    ports = (3, 5)
    pp = ((0, 3, 5, 4, 1, 2), (4, 1, 0, 2, 3, 5))
    U = native(T, pp[0])
    V = native(U, pp[1])
    need(H(T) == H(U) == H(V) == 15, 'neutral factor escape')
    stages = []
    for W in (T, U, V):
        a = ancestry(W)
        x = X(W, ports)
        marked = any(tuple(r['ports']) == ports for r in a.get('terminal_programs', []))
        stages.append({'graph': graph_record(W), 'ports': ports, 'X': x, 'D': gcd(*x),
                       'marked_reentry': marked})
    need([s['D'] for s in stages] == [3, 1, 1], 'factor divisor escape')
    need([s['marked_reentry'] for s in stages] == [True, False, True], 'marked reentry guard')
    child = birth(V, 5)
    need(H(child) == len(paths(child)) == 41, 'literal prime41')
    need(ancestry(child)['in_two_port_language'], 'prime child width-two')
    group = self_native_group(T)
    need(group['group_order'] == 4 and len(group['generators']) == 2, 'exact C4')
    orbit = next(o for o in group['ordered_pair_orbits'] if ports in o)
    a = ancestry(T)
    legal = {tuple(r['ports']) for r in a['terminal_programs']}
    need(sorted(set(orbit) & legal) == [(3, 4), (3, 5)], 'exact marked subset of orbit')
    return {'paths': pp, 'stages': stages, 'prime_child': graph_record(child),
            'self_native_group': group, 'marked_port_orbit': orbit,
            'legal_pairs_in_orbit': sorted(set(orbit) & legal)}


def module_reversal_family():
    rows = []
    for k in range(2, 7):
        T = history_graph((0,) + (3,) * k)
        n = len(T)
        ports = (3, n-1)
        V = list(T)
        for i in range(4, n):
            for j in range(4, i):
                V[i] ^= 1 << j
                V[j] ^= 1 << i
        V = tuple(V)
        a = ancestry(V)
        need(any(tuple(p['ports']) == ports for p in a['terminal_programs']), 'reversal marked reentry')
        xx = paired_queries(T, V, ports)
        h = 3 * (2**k + 1)
        expected = 6 * 2**k + 2 * 3**k - 1
        prime = lambda value: value >= 2 and all(value % d for d in range(2, isqrt(value)+1))
        need(H(T) == H(V) == h, 'module reversal H')
        need(xx['gcd_before'] == 3 and xx['gcd_after'] == 1, 'all retained queried divisors')
        need(H(birth(V, n-1)) == expected and H(birth(V, 3)) == 3*(2**(k+1)+1), 'oldest port formula')
        need(gcd(h, H(birth(V, 3)), expected) == 1, 'three-query divisor certificate')
        need(prime(expected) if k <= 5 else expected == 7*263 and prime(263), 'bounded prime-lane controls')
        p = (0, 3) + tuple(range(n-1, 3, -1)) + (1, 2)
        q = tuple(range(n-2, 3, -1)) + (1, 0, 2, 3, n-1)
        W = native(native(T, p), q)
        zz = paired_queries(T, W, ports)
        need(H(W) == h, 'latest pair transposition endpoint H')
        need(zz['gcd_after'] == (1 if k == 2 else 3), 'latest two is not oldest for k>2')
        rows.append({'k': k, 'source_rows': T, 'reversed_rows': V, 'ports': ports, 'queries': xx,
                     'H': h, 'oldest_child_H': expected, 'endpoint_ancestry': a,
                     'prime_lane': {'prime': prime(expected), 'seven_times_prime': expected % 7 == 0 and prime(expected//7)},
                     'literal_leaf_module': list(range(4, n)),
                     'stabilized_clock': {'births': k*(k-1)//2, 'native_moves': k*(k-1),
                                          'deletions': k*(k-1)//2, 'peak_order': n+1},
                     'latest_pair_only_control': {'paths': (p, q), 'middle_H': H(native(T,p)),
                                                  'endpoint_rows': W, 'queries': zz}})
    return rows


def main():
    print(json.dumps({'scope': 'selected named native returns; exact adaptive ancestry and actual static queries',
                      'static_queries': QUERIES, 'known_returns': known_returns(),
                      'neutral_procedural': neutral_procedural(),
                      'constant_to_chain': constant_to_chain(),
                      'proper_module_controls': module_controls(),
                      'factor_escape': factor_escape(),
                      'module_reversal_family': module_reversal_family(),
                      'self_native_hostiles': [self_native_group(history_graph(w))
                                               for w in ((0,), (0,0), (0,3), (0,3,3,3), (0,3,3,3,3))]},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
