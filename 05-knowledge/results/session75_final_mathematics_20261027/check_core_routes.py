#!/usr/bin/env python3
"""Independent stdlib checker for the explicit two-black-move E_t core route.
No imported repository code; no assert statements; normal and -O are equivalent.
Bounded controls only: not an exceptional-classification census.
"""
import itertools
import copy


def need(condition, message):
    if not condition:
        raise ValueError(message)


def compositions(total):
    if total == 0:
        yield ()
    for part in (1, 3):
        if part <= total:
            for rest in compositions(total - part):
                yield (part,) + rest


def set_arc(matrix, u, v):
    matrix[u][v], matrix[v][u] = True, False


def build(parts, signs, t):
    """Root=0, block vertices ascend, middles ascend, sigma(v)=n-1-v."""
    need(len(parts) >= 2 and parts[-1] == 3, 'invalid block word')
    need(all(p in (1, 3) for p in parts), 'invalid block size')
    need(len(parts) == len(signs) and signs[0] == 1 and signs[-1] == 0,
         'invalid root signs')
    need(t >= 1, 'invalid t')
    a = 1 + sum(parts)
    n = 2 * a + t - 1
    matrix = [[False] * n for _ in range(n)]
    # The ascending label order is the layer order; initialize all pairs downward.
    for v in range(n):
        for u in range(v):
            set_arc(matrix, v, u)
    blocks, offset = [], 1
    for size in parts:
        blocks.append(tuple(range(offset, offset + size)))
        offset += size
    for block, sign in zip(blocks, signs):
        for v in block:
            set_arc(matrix, 0 if sign else v, v if sign else 0)
        if len(block) == 3:
            for u, v in zip(block, block[1:] + block[:1]):
                set_arc(matrix, u, v)
    for i, block in enumerate(blocks):
        for later in blocks[i + 1:]:
            for u in block:
                for v in later:
                    set_arc(matrix, u, v)
    for u in range(a):
        for v in range(u + 1, a):
            if matrix[u][v]:
                set_arc(matrix, n - 1 - v, n - 1 - u)
            else:
                set_arc(matrix, n - 1 - u, n - 1 - v)
    middle = list(range(a, a + t - 1))
    jumps = [0] + middle + [n - 1]
    for u, v in zip(jumps, jumps[1:]):
        set_arc(matrix, u, v)
    return matrix, blocks, middle


def is_tournament(matrix):
    return all(not matrix[u][u] for u in range(len(matrix))) and all(
        matrix[u][v] != matrix[v][u]
        for u in range(len(matrix)) for v in range(u + 1, len(matrix)))


def path(matrix, word, expected=None):
    if expected is None:
        expected = set(range(len(matrix)))
    need(len(word) == len(expected) and set(word) == set(expected), 'wrong path vertex set')
    for u, v in zip(word, word[1:]):
        need(matrix[u][v], 'illegal path edge %s->%s' % (u, v))


def cycle(matrix, word, expected=None):
    path(matrix, word, expected)
    need(matrix[word[-1]][word[0]], 'illegal cycle closing edge')


def converse(matrix):
    n = len(matrix)
    return [[u != v and not matrix[u][v] for v in range(n)] for u in range(n)]


def reverse_path_arcs(matrix, word):
    path(matrix, word)
    result = [row[:] for row in matrix]
    for u, v in zip(word, word[1:]):
        result[u][v], result[v][u] = False, True
    return result


def native_k(matrix, word):
    """Independent literal K: preserve path edges, reverse every other pair."""
    path(matrix, word)
    preserve = {frozenset((u, v)) for u, v in zip(word, word[1:])}
    n = len(matrix)
    result = [[False] * n for _ in range(n)]
    for u in range(n):
        for v in range(u + 1, n):
            result[u][v] = matrix[u][v] if frozenset((u, v)) in preserve else not matrix[u][v]
            result[v][u] = not result[u][v]
    return result


def blue(matrix, word):
    reversal = dict(zip(word, reversed(word)))
    return all(matrix[u][v] == matrix[reversal[v]][reversal[u]]
               for u in range(len(matrix)) for v in range(u + 1, len(matrix)))


def strong(matrix, vertices=None):
    verts = set(range(len(matrix))) if vertices is None else set(vertices)
    if len(verts) <= 1:
        return True
    start = min(verts)
    for backward in (False, True):
        reached, queue = {start}, [start]
        while queue:
            u = queue.pop()
            for v in verts - reached:
                if matrix[v][u] if backward else matrix[u][v]:
                    reached.add(v)
                    queue.append(v)
        if reached != verts:
            return False
    return True


def rotations(block):
    return [block[i:] + block[:i] for i in range(len(block))]


def flattened(block_paths):
    return [vertex for block_path in block_paths for vertex in block_path]


def all_choices(blocks):
    for rp in rotations(blocks[-1]):
        for ap in itertools.product(*(rotations(block) for block in blocks[:-1])):
            for tp in itertools.product(*(rotations(block) for block in blocks)):
                yield rp, flattened(ap), flattened(tp)


def selected_choices(blocks):
    for rshift in range(3):
        for ashift in range(3):
            for ushift in range(3):
                rp = rotations(blocks[-1])[rshift]
                ap = [rotations(b)[ashift % len(b)] for b in blocks[:-1]]
                tp = [rotations(b)[(ushift + j) % len(b)] for j, b in enumerate(blocks)]
                yield rp, flattened(ap), flattened(tp)


def route(matrix, middle, rp, earlier, top_bottom_path):
    n = len(matrix)
    r, r1, r2 = rp
    x, a1, ak = 0, earlier[0], earlier[-1]
    # Every top HP from sigma(x) is the sigma-reverse of a bottom HP ending x.
    top = [n - 1] + [n - 1 - v for v in reversed(top_bottom_path)]
    p = [r, r1, r2, x] + middle + top + earlier
    q = [r, x, a1] + top[::-1] + middle[::-1] + earlier[:0:-1] + [r2, r1]
    final_cycle = [x, r2] + earlier[1:] + middle + top + [a1]
    return dict(P=p, Q=q, native_second=q[::-1], final_cycle=final_cycle,
                r=r, r1=r1, r2=r2, x=x, A=earlier, U=top, C=middle,
                attachments=[[a1, r1], [r1, x]])


def verify(matrix, data):
    n = len(matrix)
    p, q = data['P'], data['Q']
    r, r1, r2, a1 = data['r'], data['r1'], data['r2'], data['A'][0]
    score0 = [sum(row) for row in matrix]
    need(is_tournament(matrix), 'input not a tournament')
    need(min(score0 + [n - 1 - d for d in score0]) == 2, 'initial delta != 2')
    cycle(matrix, p)
    need(strong(matrix), 'initial tournament not strong')
    need(score0[r] == 2 and score0[data['A'][-1]] <= 5, 'P endpoint score claims fail')
    need(score0[p[0]] + score0[p[-1]] <= 7 < n - 1, 'P score blackness fails')
    need(not blue(matrix, p), 'P blue')
    t1 = reverse_path_arcs(matrix, p)
    cycle(t1, q)
    need(strong(t1), 'first R intermediate not strong')
    score1 = [sum(row) for row in t1]
    need(score1[r] == 1 and score1[r1] == 2, 'Q endpoint scores fail')
    need(score1[q[0]] + score1[q[-1]] == 3 != n - 1, 'Q score blackness fails')
    need(not blue(t1, q), 'Q blue')
    t2 = reverse_path_arcs(t1, q)
    need(sum(t2[r]) == 0, 'r is not sink')
    need(all(t2[v][r] for v in range(n) if v != r), 'sink direction fails')
    cycle(t2, data['final_cycle'], set(range(n)) - {r, r1})
    for u, v in data['attachments']:
        need(t2[u][v], 'attachment fails')
    need(data['attachments'] == [[a1, r1], [r1, 0]], 'wrong attachment names')
    need(strong(t2, set(range(n)) - {r}), 'deleted-r core not strong')
    need(not strong(t2), 'terminal unexpectedly strong')
    # Independently implement literal native K at both named words.
    k1 = native_k(matrix, p)
    need(k1 == converse(t1), 'first K/R parity mismatch')
    need(strong(k1), 'first native K intermediate not strong')
    native_second = data['native_second']
    need(native_second == q[::-1], 'native second word not literal reverse(Q)')
    path(k1, native_second)
    need(not blue(k1, native_second), 'native second word blue')
    k2 = native_k(k1, native_second)
    need(k2 == t2, 'two-step K endpoint differs from R endpoint')
    # Explicit edge-origin audit: no ambiguous unnamed transition is skipped.
    ep = {frozenset(e) for e in zip(p, p[1:])}
    eq = {frozenset(e) for e in zip(q, q[1:])}
    for u in range(n):
        for v in range(u + 1, n):
            e = frozenset((u, v))
            need(t2[u][v] == (matrix[u][v] ^ (e in ep) ^ (e in eq)), 'edge parity ledger fails')
    return {'n': n, 'P_endpoint_scores': [score0[p[0]], score0[p[-1]]],
            'Q_endpoint_scores': [score1[q[0]], score1[q[-1]]],
            'singleton': r, 'singleton_type': 'sink', 'core_order': n - 1,
            'initial_strong': True, 'first_R_strong': True, 'first_K_strong': True,
            'deleted_singleton_core_strong': True, 'literal_K_endpoint_equal': True}


def negative_controls():
    matrix, blocks, middle = build((1, 3), (1, 0), 1)
    rp, earlier, top_bottom = next(all_choices(blocks))
    good = route(matrix, middle, rp, earlier, top_bottom)
    verify(matrix, good)
    controls = []
    for name in ('illegal_P', 'wrong_native_word', 'wrong_final_cycle', 'wrong_attachment'):
        bad = copy.deepcopy(good)
        if name == 'illegal_P':
            bad['P'][0], bad['P'][1] = bad['P'][1], bad['P'][0]
        elif name == 'wrong_native_word':
            bad['native_second'] = bad['Q'][:]
        elif name == 'wrong_final_cycle':
            bad['final_cycle'][0], bad['final_cycle'][1] = bad['final_cycle'][1], bad['final_cycle'][0]
        elif name == 'wrong_attachment':
            bad['attachments'][1] = [0, bad['r1']]
        try:
            verify(matrix, bad)
        except ValueError as exc:
            controls.append({'control': name, 'rejected': True, 'reason': str(exc)})
        else:
            raise ValueError('negative control accepted: ' + name)
    return controls


def main():
    checked = 0
    for parts, signs in [((1,3),(1,0)), ((3,3),(1,0)),
                          ((1,1,1,3),(1,0,1,0)), ((3,3,3),(1,0,0))]:
        for t in (1,2,6):
            matrix, blocks, middle = build(parts, signs, t)
            for rp, earlier, top_bottom in all_choices(blocks):
                verify(matrix, route(matrix, middle, rp, earlier, top_bottom))
                checked += 1
    need(len(negative_controls()) == 4, 'missing hostile controls')
    print(f'PASS: {checked} named two-black routes; four hostile controls rejected.')

if __name__ == '__main__':
    main()
