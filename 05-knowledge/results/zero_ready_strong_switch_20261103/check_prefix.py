"""Independent small-universe audit of PREFIX.md; standard library only."""
import itertools
import json
from pathlib import Path


def profile(rows):
    n = len(rows)
    incoming = [sum(1 << u for u in range(n) if rows[u] >> v & 1)
                for v in range(n)]
    end = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                prev = mask ^ (1 << v)
                if not prev:
                    end[mask][v] = 1
                else:
                    end[mask][v] = sum(end[prev][u] for u in range(n)
                                       if prev & incoming[v] & (1 << u))
        h[mask] = sum(end[mask])
    full = (1 << n) - 1
    p = sum(h[s] * h[full ^ s] for s in range(1 << n))
    triangles = sum(h[(1 << a) | (1 << b) | (1 << c)] == 3
                    for a, b, c in itertools.combinations(range(n), 3))
    return h[full], p, triangles


def main():
    counts = []
    for n in range(1, 6):
        edges = list(itertools.combinations(range(n), 2))
        equality = [0, 0]
        for mask in range(1 << len(edges)):
            rows = [0] * n
            for i, (a, b) in enumerate(edges):
                if mask >> i & 1:
                    rows[a] |= 1 << b
                else:
                    rows[b] |= 1 << a
            h, p, tri = profile(rows)
            assert p >= 1 << n
            assert p >= 2 * h + (1 << n) - 2
            assert (p == 1 << n) == (tri == 0)
            equality[0] += p == 1 << n
            if tri:
                assert p >= 3 * (1 << (n - 1))
                assert (p == 3 * (1 << (n - 1))) == (tri == 1)
                equality[1] += p == 3 * (1 << (n - 1))
        counts.append({'n': n, 'labelled_tournaments': 1 << len(edges),
                       'transitive_equality': equality[0],
                       'single_triangle_equality': equality[1]})
    packet = Path(__file__).resolve().parent
    inherited = json.loads((packet.parent / 'triangle_return_170807_20261103'
                            / 'zero_ready.json').read_text(encoding='utf-8'))
    source = next(x for x in inherited if x['H'] == 170807)
    assert (source['h_B'], source['F_w']) == (2195, 42153)
    h, f = source['h_B'], source['F_w']
    values = [h + 2*f, h + 4*f, h + 8*f, 3*h + 12*f]
    assert values == [86501, 170807, 339419, 512421]
    assert h + 16*f == 676643 > 3*source['H'] + 4
    assert all(3*source['H'] + c not in values for c in (2, 4))
    receipt = {'status': 'FINITE-EXACT', 'universe': counts,
               'graphs_checked': sum(x['labelled_tournaments'] for x in counts),
               'source_H': source['H'], 'fixed_suffix': [h, f],
               'prefix_order_at_most_3_values': values,
               'prefix_order_at_least_4_lower_bound': h + 16*f}
    (packet / 'prefix_check.json').write_text(json.dumps(receipt, indent=2) + '\n',
                                             encoding='utf-8')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
