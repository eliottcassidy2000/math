"""Independent explicit-hostile audit for the parent odd-cell network packet."""
from itertools import product
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def require(b, message):
    if not b:
        raise ValueError(message)


def automorphisms(A):
    n = len(A)
    result = []
    # Backtracking correspondence search with both incidence orientations.
    def grow(mapping, available):
        u = len(mapping)
        if u == n:
            result.append(tuple(mapping))
            return
        for v in available:
            if all(A[u][i] == A[v][j] and A[i][u] == A[j][v]
                   for i, j in enumerate(mapping)):
                grow(mapping + [v], available - {v})
    grow([], set(range(n)))
    return result


def eleven():
    A = [[0] * 11 for _ in range(11)]
    for i in range(3):
        A[i][(i + 1) % 3] = 1
    # Four-vertex strong tournaments, independently expanded in both parts.
    for offset in (3, 7):
        for i, j in ((0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (3, 0)):
            A[offset + i][offset + j] = 1
    for i, j in product(range(4), range(3)):
        A[3 + i][j] = 1
        A[j][7 + i] = 1
    selected = {(2, 0), (3, 1)}
    for i, j in product(range(4), repeat=2):
        A[3 + i][7 + j] = int((i, j) in selected)
        A[7 + j][3 + i] = 1 - A[3 + i][7 + j]
    B = [[0] * 11 for _ in range(11)]
    for i, j in product(range(3), repeat=2):
        B[i][j] = 1
    for i, j in product(range(8), repeat=2):
        B[3 + i][3 + j] = int((j - i) % 8 in (0, 1, 2))
    require(all(A[i][j] + A[j][i] == (i != j) for i, j in product(range(11), repeat=2)), 'R11 tournament')
    require([sum(row) for row in A] == [5] * 11, 'R11 degree')
    require([sum(row) for row in B] == [3] * 11, 'R11 block row degree')
    require([sum(B[i][j] for i in range(11)) for j in range(11)] == [3] * 11, 'R11 block column degree')
    G = automorphisms(A)
    K = [g for g in G if all(B[i][j] == B[i][g[j]] for i, j in product(range(11), repeat=2))]
    require(len(G) == len(K) == 3, 'R11 actual internal/kernel group')
    require(all(sum(i == g[i] for i in range(11)) in (8, 11) for g in K), 'R11 eight fixed owners')
    return A, B, {'internal_group_order': len(G), 'one_sided_kernel_order': len(K),
                  'nonidentity_fixed_owners': list(range(3, 11))}


def nine():
    A = [[0] * 9 for _ in range(9)]
    for u, v in product(range(9), repeat=2):
        a, x = divmod(u, 3)
        b, y = divmod(v, 3)
        A[u][v] = int((y - x) % 3 == 1) if a == b else int((b - a) % 3 == 1)
    require(all(sum(row) == 4 for row in A), 'R9 degree')
    for u in range(9):
        for v in range(u):
            common = sum(A[u][w] * A[v][w] for w in range(9))
            require(common == (3 if u // 3 == v // 3 else 1), 'R9 intrinsic module detector')
    G = automorphisms(A)
    require(len(G) == 81, 'R9 internal automorphisms')
    B = [[int(u // 3 == v // 3) for v in range(9)] for u in range(9)]
    pairs = [(g, h) for g, h in product(G, repeat=2)
             if all(B[u][v] == B[g[u]][h[v]] for u, v in product(range(9), repeat=2))]
    require(len(pairs) == 2187, 'R9 pair group')
    require(all(sum(left == g for left, right in pairs) == 27 for g in G), 'R9 target fibre')
    # Each of the six triples supports its own nontrivial fixed-elsewhere
    # rotation, proving an arbitrary set of five marks cannot be a base.
    full_maps = [g + tuple(9 + x for x in h) for g, h in pairs]
    for start in range(0, 18, 3):
        rotation = list(range(18))
        for k in range(3):
            rotation[start + k] = start + (k + 1) % 3
        require(tuple(rotation) in full_maps, 'R9 independent triple rotation')
    return {'internal_group_order': 81, 'pair_group_order': 2187,
            'target_fibre_above_every_source_chart': 27, 'minimum_pair_owner_base': 6}


def main():
    A, B, R11 = eleven()
    R9 = nine()
    parent = json.loads((HERE.parent / 'networks.json').read_text(encoding='utf8'))
    expected = parent['prime_eleven_boundary']
    require([''.join(map(str, row)) for row in A] == expected['internal_rows'], 'independently built R11 differs')
    require([''.join(map(str, row)) for row in B] == expected['block_rows'], 'independently built block differs')
    print(json.dumps({'R9': R9, 'R11': R11, 'parent_literal_fixture_agreement': True}, sort_keys=True))


if __name__ == '__main__':
    main()
