#!/usr/bin/env python3
"""Twelve-command envelope: arbitrary two-port cuts and controlled forgetting.
No bounded-H saturation: zero-growth commands are permitted here.
"""
from pathlib import Path
import importlib.util
from collections import Counter, deque
from functools import lru_cache
import json

sp = importlib.util.spec_from_file_location('live', Path(__file__).with_name('live_port.py'))
live = importlib.util.module_from_spec(sp)
sp.loader.exec_module(live)


@lru_cache(None)
def kernel(state, mask, drop):
    incoming = {i for i in range(2) if mask >> i & 1}
    rename = {0: 0, 1: 1}
    if drop < 2:
        del rename[drop]
        rename[2] = drop
    ins = [None] + [i for i, (_, e) in enumerate(state) if e in incoming]
    outs = [None] + [i for i, (s, _) in enumerate(state) if s not in incoming]
    result = Counter()
    for i in ins:
        for j in outs:
            if i is not None and i == j:
                continue
            parts = [p for k, p in enumerate(state) if k != i and k != j]
            parts.append((2 if i is None else state[i][0],
                          2 if j is None else state[j][1]))
            parts = [(rename.get(s, -1), rename.get(e, -1)) for s, e in parts]
            if sum(e == -1 for _, e in parts) <= 1:
                result[tuple(sorted(parts))] += 1
    return tuple(sorted(result.items()))


def compile_carrier(letters):
    initial = live.literal_vector((2, 0), (0, 1))
    seen, todo = set(initial), deque(sorted(initial))
    while todo:
        s = todo.popleft()
        for mask, drop in letters:
            for t, _ in kernel(s, mask, drop):
                if t not in seen:
                    seen.add(t)
                    todo.append(t)
    states = tuple(sorted(seen))
    index = {s: i for i, s in enumerate(states)}
    matrices = tuple(tuple((i, index[t], c) for i, s in enumerate(states)
                          for t, c in kernel(s, mask, drop)) for mask, drop in letters)
    seed = tuple(initial[s] for s in states)
    observer = tuple(int(len(s) == 1) for s in states)
    return states, matrices, seed, observer


def decode(word):
    T, ports = (2, 0), [0, 1]
    for letter in word:
        mask, drop = divmod(letter, 3)
        n = len(T)
        old_to_new = {ports[i] for i in range(2) if mask >> i & 1}
        T = tuple(row | ((1 << n) if i in old_to_new else 0)
                  for i, row in enumerate(T)) + (sum(1 << i for i in range(n) if i not in old_to_new),)
        if drop < 2:
            ports[drop] = n
    return T, tuple(ports)


def audit(states, matrices, seed, observer):
    todo = [((), seed)]
    counts = Counter()
    while todo:
        word, vec = todo.pop()
        T, ports = decode(word)
        live.need(live.output(vec, observer) == live.count(T), 'literal graph Hamilton count')
        counts['named_graph_programs'] += 1
        if len(word) <= 2:
            literal = live.literal_vector(T, ports)
            live.need(vec == tuple(literal[s] for s in states), 'literal full cover distribution')
            live.need(not(set(literal)-set(states)), 'full cover support')
            counts['literal_cover_vectors'] += 1
        if len(word) < 3:
            for a, M in enumerate(matrices):
                todo.append((word+(a,), live.advance(vec, M)))
    # The triangle and the six adaptive commands embed literally.
    T, ports = decode((8,))
    live.need(T == (2, 4, 1) and ports == (0, 1), 'triangle seed embedding')
    return dict(counts)


def main():
    letters = [(mask, drop) for mask in range(4) for drop in range(3)]
    states, matrices, seed, observer = compile_carrier(letters)
    static_states, static_matrices, static_seed, static_observer = compile_carrier([(mask, 2) for mask in range(4)])
    live.need(states == static_states and seed == static_seed and observer == static_observer,
              'common positive carrier and initial source')
    # Query words use the four static letters only; nevertheless every new
    # port-replacement command is checked on the FULL25-dimensional carrier.
    certificate = live.rank_certificate(matrices, seed, observer, [2, 5, 8, 11])
    static_certificate = live.rank_certificate(static_matrices, static_seed, static_observer)
    live.need(certificate['hankel_rank_Q'] == static_certificate['hankel_rank_Q'] == 17,
              'common exact observer')
    print(json.dumps({'scope': 'twelve-command two-live-port envelope; not all tournaments or a target induction',
                      'letter': '3*incoming_mask+discard; suffix queries encode incoming_mask alone',
                      'states': states, 'matrices': matrices, 'seed': seed, 'observer': observer,
                      'carrier_dimension': len(states), 'rank': certificate,
                      'static_ranks': {k: static_certificate[k] for k in ['reachable_rank_Q', 'observable_rank_Q', 'hankel_rank_Q']},
                      'audit': audit(states, matrices, seed, observer)}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
