#!/usr/bin/env python3
"""Exact two-live-port adaptive language; bounded target scans, not a tail proof."""
from collections import Counter, deque
from functools import lru_cache
from fractions import Fraction
from itertools import permutations
from math import isqrt
import argparse
import hashlib
import json


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def birth(T, v):
    n = len(T)
    return tuple(row | ((1 << n) if i == v else 0)
                 for i, row in enumerate(T)) + (((1 << n)-1) ^ (1 << v),)


def count(T):
    n = len(T)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v, value in enumerate(dp[mask]):
            if not value:
                continue
            choices = T[v] & ~mask
            while choices:
                bit = choices & -choices
                choices -= bit
                dp[mask | bit][bit.bit_length()-1] += value
    return sum(dp[-1])


def literal_vector(T, ports):
    covers = set()
    for word in permutations(range(len(T))):
        for cuts in range(1 << (len(T)-1)):
            parts, start = [], 0
            for i in range(len(T)-1):
                if cuts >> i & 1:
                    parts.append(word[start:i+1])
                    start = i+1
                elif not (T[word[i]] >> word[i+1] & 1):
                    break
            else:
                parts.append(word[start:])
                covers.add(tuple(sorted(parts)))
    slot = {v: i for i, v in enumerate(ports)}
    vec = Counter()
    for cover in covers:
        state = tuple(sorted((slot.get(p[0], -1), slot.get(p[-1], -1)) for p in cover))
        if sum(e == -1 for _, e in state) <= 1:
            vec[state] += 1
    return vec


@lru_cache(None)
def kernel(state, root, drop):
    """Two slots. drop0/1 replaces that old slot; drop2 forgets newborn."""
    rename = {0: 0, 1: 1}
    if drop < 2:
        del rename[drop]
        rename[2] = drop
    incoming = [None] + [i for i, (_, e) in enumerate(state) if e == root]
    outgoing = [None] + [i for i, (s, _) in enumerate(state) if s != root]
    result = Counter()
    for i in incoming:
        for j in outgoing:
            if i is not None and i == j:
                continue
            parts = [p for k, p in enumerate(state) if k != i and k != j]
            parts.append((2 if i is None else state[i][0],
                          2 if j is None else state[j][1]))
            parts = [(rename.get(s, -1), rename.get(e, -1)) for s, e in parts]
            if sum(e == -1 for _, e in parts) <= 1:
                result[tuple(sorted(parts))] += 1
    return tuple(sorted(result.items()))


def compile_carrier():
    literal = literal_vector((2, 4, 1), (0, 1))
    seen = set(literal)
    todo = deque(sorted(seen))
    while todo:
        state = todo.popleft()
        for root in range(2):
            for drop in range(3):
                for target, _ in kernel(state, root, drop):
                    if target not in seen:
                        seen.add(target)
                        todo.append(target)
    states = tuple(sorted(seen))
    index = {s: i for i, s in enumerate(states)}
    matrices = tuple(tuple((i, index[t], c) for i, s in enumerate(states)
                          for t, c in kernel(s, r, d))
                     for r in range(2) for d in range(3))
    seed = tuple(literal[s] for s in states)
    observer = tuple(int(len(s) == 1) for s in states)
    return states, matrices, seed, observer


def advance(v, matrix):
    out = [0]*len(v)
    for i, j, c in matrix:
        out[j] += v[i]*c
    return tuple(out)


def output(v, observer):
    return sum(x*y for x, y in zip(v, observer))


def decode(word):
    T, ports, roots = (2, 4, 1), [0, 1], []
    for letter in word:
        root, drop = divmod(letter, 3)
        roots.append(ports[root])
        n = len(T)
        T = birth(T, ports[root])
        if drop < 2:
            ports[drop] = n
    return T, tuple(ports), tuple(roots)


def future_width(roots):
    return max((len({r for r in roots[i:] if r < 3+i})
                for i in range(len(roots))), default=0)


def add_basis(v, echelon):
    row = [Fraction(x) for x in v]
    for pivot, b in sorted(echelon.items()):
        c = row[pivot]
        if c:
            row = [x-c*y for x, y in zip(row, b)]
    pivot = next((i for i, x in enumerate(row) if x), None)
    if pivot is None:
        return False
    c = row[pivot]
    echelon[pivot] = tuple(x/c for x in row)
    return True


def rank_certificate(matrices, seed, observer, observation_letters=None):
    def closure(start, mats):
        echelon, rows, words = {}, [start], [()]
        need(add_basis(start, echelon), 'nonzero initial row')
        i = 0
        while i < len(rows):
            for a, M in enumerate(mats):
                v = advance(rows[i], M)
                if add_basis(v, echelon):
                    rows.append(v)
                    words.append(words[i]+(a,))
            i += 1
        return rows, words
    reach, prefixes = closure(seed, matrices)
    dual = tuple(tuple((j, i, c) for i, j, c in M) for M in matrices)
    selected_dual = dual if observation_letters is None else tuple(dual[i] for i in observation_letters)
    observe, reversed_suffixes = closure(observer, selected_dual)
    hankel = [tuple(sum(a*b for a, b in zip(r, c)) for c in observe) for r in reach]
    basis = {}
    for row in hankel:
        add_basis(row, basis)
    # Express every one-letter continuation query in the exact query basis.
    # Track the original-basis coefficients while reducing, independently of
    # the reachability rows. These are rational linear charts, not new sources.
    query_echelon = {}
    dimension = len(observe)
    for k, original in enumerate(observe):
        row = list(map(Fraction, original))
        coefficients = [Fraction(int(i == k)) for i in range(dimension)]
        for pivot, (b, c) in sorted(query_echelon.items()):
            factor = row[pivot]
            row = [x-factor*y for x, y in zip(row, b)]
            coefficients = [x-factor*y for x, y in zip(coefficients, c)]
        pivot = next(i for i, x in enumerate(row) if x)
        factor = row[pivot]
        query_echelon[pivot] = (tuple(x/factor for x in row),
                                tuple(x/factor for x in coefficients))
    def coordinates(target):
        row = list(map(Fraction, target))
        coefficients = [Fraction(0)]*dimension
        for pivot, (b, c) in sorted(query_echelon.items()):
            factor = row[pivot]
            row = [x-factor*y for x, y in zip(row, b)]
            coefficients = [x+factor*y for x, y in zip(coefficients, c)]
        need(not any(row), 'closed continuation query space')
        need(all(sum(coefficients[i]*observe[i][j] for i in range(dimension)) == target[j]
                 for j in range(len(seed))), 'full carrier intertwining')
        return coefficients
    query_matrices = []
    for M in dual:
        columns = [coordinates(advance(q, M)) for q in observe]
        query_matrices.append([[str(columns[j][i]) for j in range(dimension)]
                               for i in range(dimension)])
    result = {'reachable_rank_Q': len(reach), 'observable_rank_Q': len(observe),
            'hankel_rank_Q': len(basis), 'prefixes': prefixes,
            'suffixes': [tuple(reversed(w)) for w in reversed_suffixes],
            'hankel': hankel, 'query_matrices': query_matrices,
            'query_seed': [sum(x*y for x, y in zip(seed, q)) for q in observe]}
    if observation_letters is not None:
        result['observation_letter_map'] = observation_letters
    return result


def small_audit(states, matrices, seed, observer):
    todo = [((), seed)]
    checked = literal = 0
    while todo:
        word, vec = todo.pop()
        T, ports, roots = decode(word)
        need(output(vec, observer) == count(T), 'subset Hamilton check')
        need(future_width(roots) <= 2, 'future-live width')
        checked += 1
        if len(T) <= 5:
            actual = literal_vector(T, ports)
            need(vec == tuple(actual[s] for s in states), 'literal full cover vector')
            need(not (set(actual)-set(states)), 'complete support')
            literal += 1
        if len(T) < 6:
            for a, M in enumerate(matrices):
                todo.append((word+(a,), advance(vec, M)))
    return {'named_programs_through_order6': checked,
            'literal_cover_vectors_through_order5': literal}


def saturation(limit, matrices, seed, observer):
    seen, todo = {seed: ()}, deque([seed])
    witnesses = {3: ()}
    lengths = Counter()
    tested = 0
    while todo:
        v = todo.popleft()
        word = seen[v]
        h = output(v, observer)
        lengths[len(word)] += 1
        for a, M in enumerate(matrices):
            w = advance(v, M)
            hw = output(w, observer)
            tested += 1
            need(hw > h and hw % 2 == 1, 'actual strict odd growth')
            if hw > limit or w in seen:
                continue
            seen[w] = word+(a,)
            todo.append(w)
            witnesses.setdefault(hw, seen[w])
    missing = [h for h in range(3, limit+1, 2) if h not in witnesses and h not in (7, 21)]
    digest = hashlib.sha256()
    for v in sorted(seen):
        digest.update((','.join(map(str, v))+'\n').encode())
    # Every scalar witness is replayed on its actual named graph independently.
    for h, word in witnesses.items():
        T, _, roots = decode(word)
        need(count(T) == h and future_width(roots) <= 2, 'literal-source witness replay')
    return {'cutoff': limit, 'distinct_vectors': len(seen), 'tested_transitions': tested,
            'state_sha256': digest.hexdigest(), 'shortest_program_lengths': dict(sorted(lengths.items())),
            'missing_allowed': missing,
            'witness_programs': {str(h): witnesses[h] for h in sorted(witnesses)}}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=999)
    args = p.parse_args()
    need(args.limit >= 127, 'target bound at least127')
    states, matrices, seed, observer = compile_carrier()
    result = {'scope': 'two retained named ports, arbitrary root degree; no all-target theorem',
              'letter': '3*root_slot+discard; discard2 forgets newborn',
              'states': states, 'matrices': matrices, 'seed': seed, 'observer': observer,
              'carrier_dimension': len(states), 'audit': small_audit(states, matrices, seed, observer),
              'rank': rank_certificate(matrices, seed, observer),
              'coverage': saturation(args.limit, matrices, seed, observer)}
    word = result['coverage']['witness_programs']['127']
    T, ports, roots = decode(word)
    direct = sum(all(T[p[i]] >> p[i+1] & 1 for i in range(len(T)-1))
                 for p in permutations(range(len(T))))
    need(direct == 127, 'independent literal127')
    result['H127'] = {'program': word, 'roots': roots, 'rows': T,
                      'retained_ports': ports, 'literal_H': direct}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
