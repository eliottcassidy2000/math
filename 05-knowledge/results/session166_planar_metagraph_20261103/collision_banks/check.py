#!/usr/bin/env python3
"""Self-contained exact controls for sharp physical interpolation banks."""
from fractions import Fraction as F
from itertools import product
from math import comb
import hashlib
import json

checks = 0

def require(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)


def step(v, a):
    s, u, d = v
    return (u, s + u + d, s) if a == 'C' else (s + d, 2*u, 0)


def state(word):
    v = (1, 2, 0)
    for a in word:
        v = step(v, a)
    return v


def cinv(v):
    s, u, d = v
    return (d, s, u - d - s)


def ell(v, k):
    for _ in range(k):
        v = cinv(v)
    return v[2]


def allwords(depth):
    return [''.join(p) for n in range(depth + 1) for p in product('CD', repeat=n)]


def hbank(d):
    return ['D'*j + 'C'*k for k in range(d + 1) for j in range(d - k + 1)]


def bbank(d):
    return hbank(d) + ['D'*m + 'C' + 'D'*n + 'C'*k
                      for k in range(d + 1) for n in range(1, d - k + 1)
                      for m in range(d - k - n + 1)]


def monomials(d, homogeneous):
    return [(i, j, k) for i in range(d + 1) for j in range(d + 1 - i)
            for k in range(d + 1 - i - j)
            if not homogeneous or i + j + k == d]


def evmatrix(nodes, mons):
    return [[s**i * u**j * d**k for i, j, k in mons] for s, u, d in nodes]


def rank_q(matrix):
    a = [[F(x) for x in row] for row in matrix]
    if not a:
        return 0
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        pivot = a[r][c]
        a[r] = [x/pivot for x in a[r]]
        for i in range(r + 1, len(a)):
            q = a[i][c]
            if q:
                a[i] = [x - q*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def rank_mod(matrix, p=1000003):
    a = [[x % p for x in row] for row in matrix]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [(x*inv) % p for x in a[r]]
        for i in range(r + 1, len(a)):
            q = a[i][c]
            if q:
                a[i] = [(x - q*y) % p for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def determinant(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    prev, sign = 1, 1
    for k in range(n - 1):
        p = next((i for i in range(k, n) if a[i][k]), None)
        if p is None:
            return 0
        if p != k:
            a[k], a[p] = a[p], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = pivot*a[i][j] - a[i][k]*a[k][j]
                require(num % prev == 0, 'Bareiss exact division')
                a[i][j] = num//prev
            a[i][k] = 0
        prev = pivot
    return sign*a[-1][-1]


def add_guard(adj, q, letter):
    n = len(adj)
    out = adj[:]
    out[q] |= 1 << n
    out.append(((1 << n) - 1) ^ (1 << q))
    return out, n if letter == 'C' else q


def strong(adj):
    n = len(adj)
    for start in range(n):
        seen = 1 << start
        old = -1
        while seen != old:
            old = seen
            for i in range(n):
                if seen >> i & 1:
                    seen |= adj[i]
        if seen != (1 << n) - 1:
            return False
    return True


def paths(adj, allowed, start=None):
    n = len(adj)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        if allowed >> v & 1 and (start is None or start == v):
            dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        if mask & ~allowed:
            continue
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            choices = adj[last] & allowed & ~mask
            while choices:
                bit = choices & -choices
                choices -= bit
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += count
    return sum(dp[allowed])


def physical_state(word):
    adj, q = [2, 4, 1], 0
    for a in word:
        adj, q = add_guard(adj, q, a)
    n = len(adj)
    require(strong(adj), 'literal strongness')
    pred = [i for i in range(n) if adj[i] >> q & 1]
    require(len(pred) == 1, 'literal unique predecessor')
    full = (1 << n) - 1
    h = paths(adj, full)
    s = paths(adj, full, q)
    d = paths(adj, full ^ (1 << q), pred[0])
    return (s, h-s, d)


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def period_matrix(ts, relations):
    matrix = [[F(0) for _ in range(3*len(ts))] for _ in range(3*len(relations))]
    for r, relation in enumerate(relations):
        for i, t in enumerate(ts):
            for a in range(3):
                basis = tuple(int(j == a) for j in range(3))
                col = cross(t, basis)
                for b in range(3):
                    matrix[3*r + b][3*i + a] = relation[i]*col[b]
    return matrix


def motion_matrix(ts):
    # Columns: 3 common motions and one branch tangent for each direction.
    m = len(ts)
    return [[int(a == c) for c in range(3)] +
            [t[a] if i == j else 0 for j in range(m)]
            for i, t in enumerate(ts) for a in range(3)]


def apply_mat(a, v):
    return [sum(x*y for x, y in zip(row, v)) for row in a]


def lower_general(v, d):
    w = v
    for _ in range(d - 1):
        w = cinv(w)
    ans = w[0] - w[2] - 1
    for j in range(d - 1):
        ans *= ell(v, j)
    return ans


def lower_homogeneous(v, d):
    ans = 1
    for j in range(d):
        ans *= ell(v, j)
    return ans


def main():
    literal_words = allwords(5)
    for word in literal_words:
        require(physical_state(word) == state(word), 'literal count transfer: '+word)
    require(len(literal_words) == 63, 'literal universe')

    bank_rows = []
    # An exact full-rank reduction modulo p is a valid characteristic-zero lower bound.
    for d in range(9):
        for homogeneous, build in [(True, hbank), (False, bbank)]:
            words = build(d)
            nodes = [state(w) for w in words]
            mons = monomials(d, homogeneous)
            expected = comb(d + 2, 2) if homogeneous else comb(d + 3, 3)
            require(len(words) == expected == len(mons), 'bank cardinality')
            require(len(set(words)) == len(set(nodes)) == expected, 'distinct bank states')
            radius = d if homogeneous else (d + 1 if d else 0)
            require(max(map(len, words)) == radius, 'sharp bank radius')
            matrix = evmatrix(nodes, mons)
            require(rank_mod(matrix) == expected, 'modular full-rank certificate')
            if d <= 3:
                require(rank_q(matrix) == expected, 'independent rational rank')
            bank_rows.append({'degree': d, 'homogeneous': homogeneous,
                              'states': expected, 'radius': radius})
        if d:
            for word in allwords(d - 1):
                require(lower_homogeneous(state(word), d) == 0, 'homogeneous radius obstruction')
            require(any(lower_homogeneous(state(w), d) for w in hbank(d)),
                    'nonzero homogeneous obstruction is detected at sharp radius')
            for word in allwords(d):
                require(lower_general(state(word), d) == 0, 'total-degree radius obstruction')
            require(any(lower_general(state(w), d) for w in bbank(d)),
                    'nonzero total-degree obstruction is detected at sharp radius')

    # Separate geometry controls for the line factorization proof.
    for k in range(9):
        for j in range(k):
            for n in range(9):
                v = state('D'*n + 'C'*k)
                require(ell(v, j) > 0, 'new positive plane avoids every earlier plane')
        require(ell(state('C'*k), k) == 0, 'own plane')

    words = ['', 'D', 'DD', 'C', 'DC', 'CC']
    ts = [state(w) for w in words]
    mons = [(2,0,0), (0,2,0), (0,0,2), (1,1,0), (1,0,1), (0,1,1)]
    eq = evmatrix(ts, mons)
    require(determinant(eq) == 192, 'six-direction quadratic determinant')
    require(rank_mod(eq, 3) < 6, 'characteristic-three hostile')
    relations = [(2,-3,1,0,0,0), (-3,1,0,-1,1,0), (2,-1,0,-2,0,1)]
    for relation in relations:
        require(all(sum(relation[i]*ts[i][j] for i in range(6)) == 0 for j in range(3)),
                'exact tangent relation')
    require(rank_q(relations) == 3, 'relation basis')
    pmat = period_matrix(ts, relations)
    gmat = motion_matrix(ts)
    require(rank_q(pmat) == 9, 'complete wedge period rank')
    require(rank_q(gmat) == 9, 'motion and tangent gauge rank')
    for j in range(9):
        require(not any(apply_mat(pmat, [row[j] for row in gmat])), 'gauge contained in kernel')

    ns = [(0,F(-1,2),0), (0,F(-3,2),0), (0,F(-7,2),0),
          (0,F(1,4),F(1,4)), (0,F(1,8),F(1,8))]
    mat = [[0,0,F(1,2)], [0,0,F(-1,2)], [F(1,2),F(-1,2),1]]
    for t, n in zip(ts, ns):
        require(list(cross(t,n)) == apply_mat(mat, t), 'hostile linear wedge response')
    p5 = period_matrix(ts[:5], [r[:5] for r in relations[:2]])
    flat = [x for n in ns for x in n]
    require(not any(apply_mat(p5, flat)), 'five-direction hostile invisible')
    g5 = motion_matrix(ts[:5])
    require(rank_q([row+[x] for row,x in zip(g5,flat)]) > rank_q(g5),
            'five-direction hostile is not common motion')
    require(ts[5][2]*(ts[5][0]-ts[5][1]+ts[5][2]) == -2, 'sixth direction sees conic')

    print(json.dumps({'status':'PASS', 'always_active_checks':checks,
                      'literal_words':len(literal_words), 'literal_max_vertices':8,
                      'bank_ranks_prime':1000003, 'bank_controls':bank_rows,
                      'quadratic_determinant':192, 'complete_period_rank':9,
                      'hostile':'five directions have invisible non-common motion',
                      'source_sha256':hashlib.sha256(open(__file__,'rb').read()).hexdigest()},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
