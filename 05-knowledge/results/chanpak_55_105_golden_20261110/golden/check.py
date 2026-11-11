"""Exact, dependency-free audit of scoped golden-bridge divisibility.

Run normally and with -O. No checks depend on assert. Atom matrices are
independently recompiled from literal Boolean tournaments and subset DP.
The all-order join theorem is inherited; large witnesses use that theorem.
"""
from collections import deque
from functools import reduce
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
M = {
    's': ((0, 1, 1), (0, 0, 1), (1, 0, 1)),
    'c': ((1, 0, 4), (1, 0, 3), (3, 1, 6)),
    'q0': ((1, 0, 8), (1, 0, 7), (7, 1, 18)),
    'q1': ((3, 0, 6), (3, 0, 5), (9, 3, 14)),
    'q2': ((4, 1, 7), (3, 1, 6), (6, 1, 14)),
    'q4': ((3, 0, 12), (3, 0, 9), (5, 1, 14)),
    'q5': ((4, 1, 9), (4, 0, 7), (7, 2, 11)),
}
I = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
WORDS = {
    3: ('q5', 's', 'q2', 'q5', 's', 's', 'q5'),
    5: ('q2', 'c', 'c', 'q2', 's', 's', 's', 's', 'q2'),
    7: ('q5', 'q0', 'q2', 'q5', 'q2', 'q5'),
}
checks = 0


def need(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)


def mul(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def vec(v, a):
    return tuple(sum(v[k]*a[k][j] for k in range(3)) for j in range(3))


def det(a):
    return sum((-1 if sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)) % 2 else 1)
               *a[0][p[0]]*a[1][p[1]]*a[2][p[2]] for p in permutations(range(3)))


def smith(a):
    d1 = reduce(gcd, (abs(x) for row in a for x in row))
    d12 = reduce(gcd, (abs(a[i][j]*a[k][l]-a[i][l]*a[k][j])
                       for i, k in combinations(range(3), 2)
                       for j, l in combinations(range(3), 2)))
    return (d1, d12//d1, abs(det(a))//d12)


def valuation(x, p):
    if not x:
        return None
    k = 0
    while x % p == 0:
        x //= p
        k += 1
    return k


def rref(a, p):
    a = [list(row) for row in a]
    r = 0
    for c in range(3):
        pivot = next((k for k in range(r, len(a)) if a[k][c] % p), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        z = pow(a[r][c] % p, -1, p)
        a[r] = [x*z % p for x in a[r]]
        for k in range(len(a)):
            if k != r:
                z = a[k][c] % p
                a[k] = [(x-z*y) % p for x, y in zip(a[k], a[r])]
        r += 1
    return tuple(tuple(row) for row in a[:r])


def subspace_graph(p):
    queue = deque([I])
    seen = {I: ()}
    while queue:
        space = queue.popleft()
        for name, matrix in M.items():
            child = rref([vec(row, matrix) for row in space], p)
            if child not in seen:
                seen[child] = seen[space]+(name,)
                queue.append(child)
    return {'states': len(seen),
            'ranks': [sum(len(a) == k for a in seen) for k in range(4)],
            'shortest_letter_word_to_zero': seen.get(())}


def chart(n, bits):
    rows = [0]*n
    k = 0
    for i, j in combinations(range(n), 2):
        forward = j == i+1 or bits >> k & 1
        if j != i+1:
            k += 1
        u, v = (i, j) if forward else (j, i)
        rows[u] |= 1 << v
    return tuple(rows)


def atom(name):
    n, bits = (1, 0) if name == 's' else (3, 0) if name == 'c' else (4, int(name[1:]))
    rows = list(chart(n, bits))
    for i, j in combinations(range(n), 2):
        if j != i+1:
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return tuple(sum((1 << (n-1-j)) for j in range(n) if rows[n-1-i] >> j & 1)
                 for i in range(n))


def join(rows, block):
    m, n = len(rows), len(block)
    if not m:
        return block
    out = [r | (((1 << n)-1) << m) for r in rows]+[r << m for r in block]
    out[m-1] &= ~(1 << m)
    out[m] |= 1 << (m-1)
    return tuple(out)


def literal(word):
    rows = ()
    for name in word:
        rows = join(rows, atom(name))
    return rows


def state(rows):
    n = len(rows)
    lim = 1 << n
    end = [[0]*n for _ in range(lim)]
    start = [[0]*n for _ in range(lim)]
    for i in range(n):
        end[1 << i][i] = 1
    start[1 << (n-1)][n-1] = 1
    for mask in range(1, lim):
        for i in range(n):
            if not (end[mask][i] or start[mask][i]):
                continue
            options = rows[i] & ~mask
            while options:
                bit = options & -options
                options -= bit
                j = bit.bit_length()-1
                end[mask | bit][j] += end[mask][i]
                start[mask | bit][j] += start[mask][i]
    total = [sum(a) for a in end]
    total[0] = 1
    h, b = total[-1], end[-1][-1]
    c = sum(sum(start[mask])*total[(lim-1) ^ mask]
            for mask in range(lim) if mask >> (n-1) & 1)
    return (h-b, b, c)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='CHECK.json')
    args = parser.parse_args()
    # A different path to each atom: actual 1/2/3-vertex old sources, subset DP.
    prefixes = [literal(('s',)*i) for i in (1, 2, 3)]
    base = tuple(state(rows) for rows in prefixes)
    need(det(base) == 2, 'actual state basis index two')
    for name, matrix in M.items():
        y1, y2, y3 = [state(join(rows, atom(name))) for rows in prefixes]
        c2 = tuple(2*b+a-c for a, b, c in zip(y1, y2, y3))
        need(all(x % 2 == 0 for x in c2), 'integer recovery')
        cc = tuple(x//2 for x in c2)
        aa = tuple(x-2*z for x, z in zip(y2, cc))
        bb = tuple(x-z for x, z in zip(y1, cc))
        need((aa, bb, cc) == matrix, ('literal atom matrix', name))
    data = {}
    for p, word in WORDS.items():
        matrix = reduce(mul, (M[x] for x in word), I)
        need(all(x % p == 0 for row in matrix for x in row), ('divisible', p))
        normalized = tuple(tuple(x//p for x in row) for row in matrix)
        need(gcd(det(normalized), p) == 1, ('primitive contraction', p))
        image = {tuple(x % p for x in vec(v, normalized)) for v in product(range(p), repeat=3)}
        need(len(image) == p**3, ('full normalized phase survives', p))
        powers = []
        power = I
        for k in range(1, 5):
            power = mul(power, matrix)
            d = reduce(gcd, (abs(x) for row in power for x in row))
            need(valuation(d, p) == k, ('exact valuation depth', p, k))
            powers.append({'depth': k, 'entry_gcd_valuation': valuation(d, p)})
        rows = literal(word)
        n = len(rows)
        need(all(not (r >> i & 1) for i, r in enumerate(rows)), 'no loops')
        need(all(((rows[i] >> j & 1)+(rows[j] >> i & 1)) == 1
                 for i, j in combinations(range(n), 2)), 'tournament')
        data[str(p)] = {'word': word, 'vertex_cost': n, 'matrix': matrix,
                        'determinant': det(matrix), 'smith': smith(matrix),
                        'normalized_matrix': normalized,
                        'normalized_determinant': det(normalized),
                        'normalized_smith': smith(normalized),
                        'literal_block_rows': rows, 'first_port': 0, 'last_port': n-1,
                        'powers': powers}
    short3 = mul(M['q4'], M['q1'])
    need(state(literal(('q4', 'q1'))) == short3[0], 'literal eight-vertex composite')
    # Complete subspace searches: 28,64,116 are *all* subspaces at these primes.
    graphs = {str(p): subspace_graph(p) for p in (3, 5, 7, 11)}
    for p in (3, 5, 7):
        need(graphs[str(p)]['states'] == 2*p*p+2*p+4, ('full Grassmann reachability', p))
    need(graphs['11']['states'] == 1, 'eleven-unit hostile')
    # Retain the two inherited actual H55 source witnesses, no conflation.
    twins = []
    w7 = data['7']['matrix']
    for word, expected in [(('s', 's', 'q1', 's'), (34, 21, 61)),
                           (('s', 'q4', 's', 's'), (32, 23, 63))]:
        rows = literal(word)
        v = state(rows)
        need(v == expected, 'inherited H55 witness')
        output = vec(v, w7)
        twins.append({'word': word, 'rows': rows, 'state': v,
                      'after_W7': output,
                      'normalized_phase_mod7': [x//7 % 7 for x in output]})
    need(twins[0]['normalized_phase_mod7'] != twins[1]['normalized_phase_mod7'],
         'zero raw phase does not identify H55 twins')
    # Integral phase transport at higher precision, not only the residue field.
    for p in (3, 5, 7):
        u = data[str(p)]['normalized_matrix']
        for k in (1, 2, 3):
            need(gcd(det(u), p**k) == 1, ('phase invertible at every tested precision', p, k))
    receipt = {'status': 'PASS', 'checks': checks,
               'universe': 'seven literal golden atoms; exact 3-coordinate single-bridge observer',
               'determinants': {name: det(a) for name, a in M.items()},
               'primitive_contractions': data,
               'short_three_contraction': {'word': ['q4', 'q1'], 'matrix': short3,
                    'smith': smith(short3), 'vertex_cost': 8},
               'complete_subspace_graphs': graphs, 'H55_twins': twins}
    (HERE/args.output).write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'checks': checks,
                      'subspace_counts': {p: d['states'] for p, d in graphs.items()}}))


if __name__ == '__main__':
    main()
