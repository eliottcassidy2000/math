#!/usr/bin/env python3
"""Literal positive profiles and the variable-order scope control for H=97309."""
import json
from functools import lru_cache


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def decode(n, code):
    choices = [0] * n
    for size in range(n, 1, -1):
        choices[size - 1], code = code % size, code // size
    need(code == 0, 'code range')
    graph = [0] * n
    for v in range(1, n):
        for u in range(v):
            if choices[v] == u + 1:
                graph[v] |= 1 << u
            else:
                graph[u] |= 1 << v
    return tuple(graph)


def history(graph):
    depths, choices, code = [0], [0], 0
    for v in range(1, len(graph)):
        backward = graph[v] & ((1 << v) - 1)
        need(backward.bit_count() <= 1, 'not a chronological forest')
        choice = backward.bit_length()
        choices.append(choice)
        depths.append(depths[choice - 1] + 1 if choice else 0)
        code = code * (v + 1) + choice
    return {'chronological_code': code, 'choices_including_vertex_zero': choices,
            'vertex_depths': depths, 'maximum_history_depth': max(depths)}


def profiles(graph):
    n = len(graph)
    incoming = tuple(sum(1 << u for u in range(n) if graph[u] >> v & 1)
                     for v in range(n))

    @lru_cache(None)
    def endpoint(mask, v):
        if mask == 1 << v:
            return 1
        rest = mask ^ (1 << v)
        return sum(endpoint(rest, u) for u in range(n) if rest & incoming[v] & (1 << u))

    @lru_cache(None)
    def start(mask, v):
        if mask == 1 << v:
            return 1
        rest = mask ^ (1 << v)
        return sum(start(rest, u) for u in range(n) if rest & graph[v] & (1 << u))

    @lru_cache(None)
    def count(mask):
        return sum(endpoint(mask, v) for v in range(n) if mask >> v & 1) if mask else 1

    full = (1 << n) - 1
    h = count(full)
    ends = [endpoint(full, v) for v in range(n)]
    starts = [start(full, v) for v in range(n)]
    need(sum(starts) == sum(ends) == h, 'endpoint/start totals')
    deletions = [count(full ^ (1 << v)) for v in range(n)]
    covers = [sum(start(s, v) * count(full ^ s) for s in range(1, full + 1)
                  if s >> v & 1) for v in range(n)]
    return {'H': h, 'starts': starts, 'ends': ends,
            'vertex_deletion_H': deletions, 'root_two_path_cover_U': covers,
            'singleton_child_H_including_sink': [h] + [h + u - e for u, e in zip(covers, ends)]}


def main():
    graph = decode(12, 25112126)
    graph_history = history(graph)
    profile = profiles(graph)
    need(profile['H'] == 97309 and graph_history['maximum_history_depth'] == 3,
         'order12 positive witness')
    mask, forced = (1 << 12) - 1, []
    while mask:
        candidates = [v for v in range(12) if mask >> v & 1
                      and (graph[v] & mask).bit_count() <= 1]
        if len(candidates) != 1:
            break
        v = candidates[0]
        parent = (graph[v] & mask).bit_length() - 1
        forced.append([v, parent])
        mask ^= 1 << v
    need(forced[:8] == [[11, 1], [10, 3], [9, 2], [8, 6], [7, 0], [6, 4], [5, 0], [4, 0]],
         'forced elimination certificate')
    factors, joined, offset = [], [], 0
    for n, code, value in [(6, 50, 31), (7, 296, 43), (7, 58, 73)]:
        factor = decode(n, code)
        description = history(factor)
        factor_profile = profiles(factor)
        need(description['maximum_history_depth'] == 1 and factor_profile['H'] == value,
             'depth-one factor witness')
        factors.append({'n': n, 'H': value, 'outgoing_masks': factor,
                        **description, 'profile': factor_profile})
        for v in range(offset):
            joined[v] |= ((1 << n) - 1) << offset
        joined.extend(row << offset for row in factor)
        offset += n
    join_history = history(joined)
    need(len(joined) == 20 and join_history['maximum_history_depth'] == 1,
         'joined depth-one history')
    boundaries = [0, 6, 13, 20]
    for i in range(3):
        for j in range(i + 1, 3):
            for v in range(boundaries[i], boundaries[i + 1]):
                for u in range(boundaries[j], boundaries[j + 1]):
                    need(joined[v] >> u & 1 and not joined[u] >> v & 1, 'ordinal cross edge')
    need(31 * 43 * 73 == 97309, 'factor product')
    return {'status': 'PASS', 'order12_witness': {'n': 12, 'outgoing_masks': graph,
            **graph_history, 'profile': profile, 'forced_reverse_deletions': forced},
            'depth_one_factors': factors, 'order20_depth_one_witness': {
                'n': 20, 'H': 97309, 'outgoing_masks': joined, **join_history,
                'count_method': 'Exact ordinal-join multiplication of three independently literal-DP-counted factors; every cross edge verified.'},
            'scope': 'Order12 depth-three positive witness and exact depth-one order20 construction. The full same-order shallow exclusion is established separately by the exhaustive shallow census, not by these positive profiles.'}


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
