"""Independent finite replay of two pinned concurrent packets; no repo imports.

Construct the P7 marker witness literally, exhaust S7, count Hamilton paths
two ways, and check both nonabelian frame conventions on all 21^3 four-site
configurations. Also replay the five supplied order-eight separating contexts,
not the exhaustive horizon-seven partition which selected those five pairs.
All guards remain active under python -O. Stdlib only; emits JSON to stdout.
"""
from collections import deque
from itertools import combinations, permutations, product
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def paley7():
    return tuple(tuple((j-i) % 7 in (1, 2, 4) for j in range(7))
                 for i in range(7))


def preserves(a, b, p):
    return all(a[i][j] == b[p[i]][p[j]]
               for i in range(len(a)) for j in range(i+1, len(a)))


def canonical_tuple(group, ports):
    return min(tuple(g[x] for x in ports) for g in group)


def marker_graph(t, ports, order):
    n = len(t)
    a = [list(row) + [False]*3 for row in t] + [[False]*(n+3) for _ in ports]
    for j, port in enumerate(ports):
        for v in range(n):
            a[n+j][v] = v == port
            a[v][n+j] = v != port
    for i, j in combinations(range(3), 2):
        a[n+i][n+j] = order.index(i) < order.index(j)
        a[n+j][n+i] = not a[n+i][n+j]
    return tuple(tuple(row) for row in a)


def hamilton_dp(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            if count:
                for w in range(n):
                    if not (mask >> w & 1) and a[v][w]:
                        dp[mask | 1 << w][w] += count
    return sum(dp[-1])


def hamilton_dfs(a):
    n = len(a)
    def visit(path, unused):
        if not unused:
            return 1
        return sum(visit(path + (v,), unused - {v})
                   for v in unused if a[path[-1]][v])
    return sum(visit((v,), set(range(n)) - {v}) for v in range(n))


def compose(g, h):
    """g after h, as literal permutations."""
    return tuple(g[h[v]] for v in range(len(g)))


def inverse(g):
    return tuple(g.index(v) for v in range(len(g)))


def frame_check(group):
    ident = tuple(range(7))
    edge_list = ((0, 1), (0, 2), (2, 3))
    left_records, right_records = set(), set()
    pair_checks = 0
    gauge = next(g for g in group if g != ident)
    for tail in product(group, repeat=3):
        charts = (ident,) + tail
        invs = tuple(inverse(g) for g in charts)
        left = {(u, v): compose(invs[u], charts[v]) for u, v in edge_list}
        right = {(u, v): compose(charts[v], invs[u]) for u, v in edge_list}
        left_records.add(tuple(left[e] for e in edge_list))
        right_records.add(tuple(right[e] for e in edge_list))
        for records in (left, right):
            for u, v in edge_list:
                records[v, u] = inverse(records[u, v])
        for start in range(4):
            q = deque([(start, ident, ident)])
            seen = {start}
            while q:
                u, a, b = q.popleft()
                need(a == compose(invs[start], charts[u]), 'left-gauge path product')
                need(b == compose(charts[u], invs[start]), 'right-gauge path product')
                pair_checks += 1
                for v in range(4):
                    if (u, v) in left and v not in seen:
                        seen.add(v)
                        q.append((v, compose(a, left[u, v]), compose(right[u, v], b)))
        lc = tuple(compose(gauge, g) for g in charts)
        rc = tuple(compose(g, gauge) for g in charts)
        for u, v in edge_list:
            need(compose(inverse(lc[u]), lc[v]) == left[u, v], 'left gauge invariance')
            need(compose(rc[v], inverse(rc[u])) == right[u, v], 'right gauge invariance')
    need(len(left_records) == len(right_records) == 21**3, 'all frame alignments distinct')
    noncommuting = next((g, h) for g, h in product(group, repeat=2)
                       if compose(g, h) != compose(h, g))
    return {'sites': 4, 'relative_configurations': 21**3,
            'tree': edge_list, 'path_product_checks_including_identities': pair_checks,
            'left_gauge_gain': 'inverse(g_i) * g_j',
            'right_gauge_gain': 'g_j * inverse(g_i)',
            'noncommuting_witness': noncommuting}


def presentation(n, mask):
    a = [[False]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            forward = True if j == i+1 else bool(mask >> k & 1)
            if j != i+1:
                k += 1
            a[i][j], a[j][i] = forward, not forward
    return tuple(tuple(row) for row in a)


def free_mask(a):
    bits = [a[i][j] for i in range(len(a)) for j in range(i+2, len(a))]
    return sum(bit << k for k, bit in enumerate(bits))


def apply_word(mask, word):
    a = presentation(6, mask)
    for op in word:
        n = len(a)
        if op == 'complement':
            a = presentation(n, free_mask(a) ^ ((1 << ((n-1)*(n-2)//2))-1))
        elif op == 'singleton-right':
            a = tuple(tuple(row) + (True,) for row in a) + ((False,)*(n+1),)
        elif op == 'singleton-left':
            a = ((False,) + (True,)*n,) + tuple((False,) + tuple(row) for row in a)
        else:
            raise ValueError(op)
    return a


def residual_check():
    pairs = ((100, 689), (225, 564), (326, 697), (334, 923), (459, 798))
    expected = ((12544, 301120), (45120, 268544), (301632, 74112),
                (301120, 12544), (268544, 45120))
    results = []
    for i, (pair, masks) in enumerate(zip(pairs, expected)):
        word = (('complement',) if i < 2 else ()) + (
            'singleton-right', 'singleton-left', 'complement')
        a, b = (apply_word(x, word) for x in pair)
        need((free_mask(a), free_mask(b)) == masks, 'literal separating context masks')
        degree_a, degree_b = tuple(map(sum, a)), tuple(map(sum, b))
        candidates, isos = 0, 0
        for p in permutations(range(8)):
            if all(degree_a[v] == degree_b[p[v]] for v in range(8)):
                candidates += 1
                isos += preserves(a, b, p)
        need(isos == 0, 'order-eight pair is nonisomorphic')
        results.append({'pair': pair, 'word': word, 'terminal_masks': masks,
                        'permutations_exhausted': 40320, 'degree_compatible': candidates,
                        'isomorphisms': isos})
    return results


def main():
    p7 = paley7()
    group = tuple(g for g in permutations(range(7)) if preserves(p7, p7, g))
    need(len(group) == 21, 'P7 full automorphism group')
    affine = {tuple((pow(2, a, 7)*v+b) % 7 for v in range(7))
              for a in range(3) for b in range(7)}
    need(set(group) == affine, 'literal affine group')
    singleton_stabilizers = [sum(g[v] == v for g in group) for v in range(7)]
    need(singleton_stabilizers == [3]*7, 'no one-port free frame')
    frame = (0, 1)
    need(sum(all(g[v] == v for v in frame) for g in group) == 1, 'two-port free frame')
    need(len({tuple(g[v] for v in frame) for g in group}) == 21, 'free frame orbit')
    a, b = (0, 1, 3), (0, 1, 5)
    pair_keys = [[canonical_tuple(group, (ports[i], ports[j]))
                  for i, j in permutations(range(3), 2)] for ports in (a, b)]
    need(pair_keys[0] == pair_keys[1], 'all six ordered pair observations coincide')
    joint_keys = [canonical_tuple(group, ports) for ports in (a, b)]
    need(joint_keys[0] != joint_keys[1], 'joint triple observations differ')
    common = [sum(all(p7[u][v] for u in ports) for v in range(7)) for ports in (a, b)]
    need(common == [0, 1], 'common outneighbors')
    markers = []
    for order in permutations(range(3)):
        graphs = [marker_graph(p7, ports, order) for ports in (a, b)]
        dp = list(map(hamilton_dp, graphs))
        dfs = list(map(hamilton_dfs, graphs))
        need(dp == dfs, 'independent subset DP and unmemorized path DFS')
        need(dp in ([4281, 4315], [4287, 4325]), 'marker expected counts')
        need(dp[0] != dp[1], 'uncoloured attachment separates')
        markers.append({'marker_order': order, 'Hamilton_DP': dp, 'Hamilton_DFS': dfs})
    need(markers[0]['Hamilton_DP'] == [4281, 4315], 'specified marker order')
    output = {'status': 'PASS', 'imports_from_source_packets': False,
              'P7': {'group_order': len(group), 'S7_permutations_exhausted': 5040,
                     'minimum_free_frame': frame, 'minimum_free_frame_size': 2,
                     'single_port_stabilizer_sizes': singleton_stabilizers,
                     'triple_ports': [a, b], 'all_ordered_pair_keys': pair_keys,
                     'joint_keys': joint_keys, 'common_outneighbors': common,
                     'markers': markers, 'four_site_frame_replay': frame_check(group)},
              'five_supplied_residual_contexts': residual_check(),
              'not_replayed': ['complete native complement/reset closure at n=7',
                               'complete horizon-7 observational partition',
                               'order-nine class generators or comparative timings',
                               'all other strong-atom frame-size claims']}
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
