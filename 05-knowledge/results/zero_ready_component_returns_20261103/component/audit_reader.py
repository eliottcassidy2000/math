#!/usr/bin/env python3
"""Audit the reader against literal DP and independently enumerated short words."""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'inverse/reader.py'
spec = importlib.util.spec_from_file_location('reader_under_audit', SOURCE)
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def counts(a):
    n = len(a)
    ends = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for v in range(n):
        ends[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not mask >> v & 1:
                continue
            for w in range(n):
                if not mask >> w & 1 and a[v] >> w & 1:
                    ends[mask | (1 << w)][w] += ends[mask][v]
        h[mask] = sum(ends[mask])
    return h, ends


def fields(a, root):
    h, end = counts(a)
    reverse = tuple(sum(1 << j for j in range(len(a)) if a[j] >> i & 1)
                    for i in range(len(a)))
    hr, start = counts(reverse)
    need(h == hr, 'independent reversal')
    p = [0] * len(h)
    for mask in range(len(h)):
        sub = mask
        while True:
            p[mask] += h[sub] * h[mask ^ sub]
            if sub == 0:
                break
            sub = (sub - 1) & mask
    full = len(h)-1
    return (h[full], end[full][root],
            sum(start[m][root] * h[full ^ m] for m in range(len(h))),
            sum(end[m][root] * h[full ^ m] for m in range(len(h))),
            p[full ^ (1 << root)], h[full ^ (1 << root)], p[full],
            sum(start[m][root] * p[full ^ m] for m in range(len(h))))


def append_word(a, root, word):
    a = list(a)
    for bit in word:
        new = len(a)
        cut = 0 if bit == '1' else 1 << root
        a = [row + ((0 if cut >> i & 1 else 1) << new) for i, row in enumerate(a)]
        a.append(cut)
    return tuple(a)


def connected(a):
    for source in range(len(a)):
        reached = {source}
        todo = [source]
        while todo:
            u = todo.pop()
            for v in range(len(a)):
                if a[u] >> v & 1 and v not in reached:
                    reached.add(v)
                    todo.append(v)
        if len(reached) != len(a):
            return False
    return True


def standard_h(word):
    """Explicit pair sum, no reader recurrence: each earlier1 contributes at a0."""
    return 1 + sum((1 << (i-j-1)) + int(i == j+1)
                   for i, bit in enumerate(word) if bit == '0'
                   for j in range(i) if word[j] == '1')


def closed(seed, word):
    h, e, t, f, p, q, _, _ = seed
    r = sum(1 << i for i, bit in enumerate(word) if bit == '0')
    return (h + q*(standard_h(word)-1) + (t-q)*r + (q-e)*(r & 1),
            f + p*r, r)


def all_words(maximum):
    for length in range(maximum+1):
        yield from (''.join(bits) for bits in itertools.product('01', repeat=length))


def main():
    cap, maximum = 1001, 10
    seed_graphs = [((2, 0), 0), ((2, 4, 1), 0), ((2, 12, 9, 1), 3)]
    seeds = [(a, root, fields(a, root)) for a, root in seed_graphs]
    literal = pair_sum_controls = scalar_controls = inversions = negative = positive = 0
    tried = 0
    for word in all_words(7):
        literal_h = counts(append_word((0,), 0, word))[0][-1]
        need(literal_h == standard_h(word) == reader.word_h(word), 'one-apex pair sum')
        pair_sum_controls += 1
    records = []
    for a, root, seed in seeds:
        need(reader.all_fields(a)[root] == seed, 'independent seed fields')
        strong_seed = connected(a)
        for word in all_words(4):
            graph = append_word(a, root, word)
            actual = fields(graph, root)
            need(reader.graph(a, root, word) == graph, 'literal named reader graph')
            need(reader.run(seed, word) == actual, 'literal eight-field update')
            h, f, r = closed(seed, word)
            need((h, f) == (actual[0], actual[3]), 'literal closed formula')
            need(connected(graph) == (strong_seed if not word else word[-1] == '0'),
                 'order>=2 seed strong-ending criterion')
            literal += 1
        expected = {target: set() for target in range(1, cap+1)}
        max_possible_r = 0
        h0, e, t, f0, p, q, _, _ = seed
        clock, multiplier = 1, 2
        while h0 + multiplier*f0 <= cap:
            k = t-q+multiplier*p
            max_possible_r = max(max_possible_r, (cap-h0-multiplier*f0)//k)
            clock += 1
            multiplier *= 2
        need(max_possible_r < 1 << maximum, 'brute word cap contains every inverse candidate')
        for word in all_words(maximum):
            h, f, r = closed(seed, word)
            out = reader.run(seed, word)
            need((h, f) == (out[0], out[3]), 'closed formula vs complete state replay')
            need(0 <= standard_h(word)-1 <= r, 'interval bound for every brute word')
            scalar_controls += 1
            if (not word and not strong_seed) or (word and word[-1] != '0'):
                continue
            if word:
                need(reader.decode_r(r) == word, 'reversed zero-mask bijection')
            clock, multiplier = 1, 2
            while h + multiplier*f <= cap:
                target = h + multiplier*f
                expected[target].add((clock, word))
                if word:
                    residual = target-h0-multiplier*f0
                    k, d = t-q+multiplier*p, q-e
                    lower = max(1, (residual-d+k+q-1)//(k+q))
                    upper = residual//k
                    need(lower <= r <= upper, 'independent inverse interval')
                clock += 1
                multiplier *= 2
        seed_positive = seed_negative = 0
        for target in range(1, cap+1):
            found, count = reader.inverse(seed, target, include_empty=strong_seed)
            tried += count
            observed = {(hit['prefix'], hit['word']) for hit in found}
            need(len(observed) == len(found), 'inverse emits no duplicate words')
            need(observed == expected[target], ('complete brute inverse match', a, root, target))
            inversions += 1
            if observed:
                positive += 1
                seed_positive += 1
            else:
                negative += 1
                seed_negative += 1
        records.append({'seed_rows': a, 'root': root, 'fields': seed,
                        'include_empty': strong_seed, 'max_required_r': max_possible_r,
                        'positive_targets': seed_positive, 'negative_targets': seed_negative})
    # The arithmetic remains valid for a singleton, but the strong-ending
    # interpretation requires an explicit order>=2 boundary.
    singleton = fields((0,), 0)
    hostile = append_word((0,), 0, '0')
    need(not connected(hostile) and closed(singleton, '0')[:2] == (1, 2),
         'singleton/zero hostile boundary')
    return {'status': 'PASS', 'literal_full_field_graphs': literal,
            'independent_one_apex_pair_sum_controls': pair_sum_controls,
            'word_maximum': maximum, 'all_word_scalar_controls': scalar_controls,
            'target_cap': cap, 'target_seed_inverse_checks': inversions,
            'positive_target_checks': positive, 'negative_target_checks': negative,
            'inverse_candidate_integers': tried, 'seeds': records,
            'singleton_zero_boundary': {'rows': hostile, 'strong': False, 'h': 1, 'F': 2},
            'reader_source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'auditor_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope': 'Independent finite checks of the narrow closed formula and inverse. The interval proof requires 0<=H_word-1<=r and the stated positive seed bounds. Strong-ending interpretation additionally assumes seed order>=2 and the retained mark reaches every seed vertex.'}


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
