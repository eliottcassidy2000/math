"""Exact finite controls; only the Python standard library. No assertions."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, factorial
import hashlib
import json
from pathlib import Path

SOURCE = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
checks = 0


def require(ok, context):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(context)


def tournament(n, code):
    rows = [0] * n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if (code >> bit) & 1:
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return tuple(rows)


def path_counts(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for mask in range(1, 1 << n):
        for u in range(n):
            if not (mask >> u) & 1:
                continue
            rest = mask ^ (1 << u)
            dp[mask][u] = (sum(dp[rest][v] for v in range(n)
                                  if (rows[v] >> u) & 1) if rest else 1)
        h[mask] = sum(dp[mask])
    return h


def profile(rows):
    h = path_counts(rows)
    n = len(rows)
    q = [[0] * (n + 1) for _ in h]
    q[0][0] = 1
    for mask in range(1, len(h)):
        sub = mask
        while sub:
            rest = mask ^ sub
            for k in range(1, mask.bit_count() + 1):
                q[mask][k] += h[sub] * q[rest][k - 1]
            sub = (sub - 1) & mask
    return tuple(q[-1][1:])


def poly(q):
    result = [F(0)] * (len(q) + 1)
    falling = [F(1)]
    for k, value in enumerate(q, 1):
        nxt = [F(0)] * (len(falling) + 1)
        for i, c in enumerate(falling):
            nxt[i] -= (k - 1) * c
            nxt[i + 1] += c
        falling = nxt
        for i, c in enumerate(falling):
            result[i] += value * c / factorial(k)
    return tuple(result)


def rank(vectors):
    basis = {}
    for vector in vectors:
        row = list(map(F, vector))
        for col in sorted(basis):
            if row[col]:
                scale = row[col]
                row = [a - scale * b for a, b in zip(row, basis[col])]
        pivot = next((i for i, x in enumerate(row) if x), None)
        if pivot is not None:
            scale = row[pivot]
            basis[pivot] = [x / scale for x in row]
    return len(basis)


def solve(matrix, rhs):
    n = len(rhs)
    a = [list(map(F, row)) + [F(b)] for row, b in zip(matrix, rhs)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if a[i][col])
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [x / scale for x in a[col]]
        for i in range(n):
            if i != col:
                scale = a[i][col]
                a[i] = [x - scale * y for x, y in zip(a[i], a[col])]
    return [row[-1] for row in a]


@lru_cache(None)
def skeleton(r):
    m = 2 * r

    def count(k):
        @lru_cache(None)
        def rec(mask, used, last):
            if mask == (1 << m) - 1 and used == k:
                return 1
            total = 0
            # last=-2 is the start; last=-1 is X. Xs are indistinguishable.
            if used < k and (last == -2 or (last >= 0 and last % 2 == 1)):
                total += rec(mask, used + 1, -1)
            for v in range(m):
                if (mask >> v) & 1:
                    continue
                legal = last == -2 or (last == -1 and v % 2 == 0) or (last >= 0 and last < v)
                if legal:
                    total += rec(mask | (1 << v), used, v)
            return total
        return rec(0, 0, -2)
    return tuple(count(k) for k in range(1, r + 2))


def wrap(rows, r):
    n = len(rows)
    out = list(rows) + [0] * (2 * r)
    for j in range(2 * r):
        v = n + j
        for u in range(n):
            a, b = (u, v) if j % 2 == 0 else (v, u)
            out[a] |= 1 << b
        for k in range(j + 1, 2 * r):
            out[v] |= 1 << (n + k)
    return tuple(out)


def recover(n, observations):
    d = (n + 1) // 2
    require(len(observations) == d - 1, 'probe number')
    q = []
    for r, obs in enumerate(observations):
        c = skeleton(r)
        val = F(obs - sum(x * y for x, y in zip(c, q)), c[-1])
        q.append(val)
    exponents = list(range(1 if n % 2 else 2, n, 2))
    values = [sum(q[k - 1] * comb(t, k) for k in range(1, t + 1)) - t ** n
              for t in range(1, d)]
    coeff = solve([[t ** e for e in exponents] for t in range(1, d)], values)
    result = [F(0)] * (n + 1)
    result[n] = F(1)
    for e, v in zip(exponents, coeff):
        result[e] = v
    return tuple(result)


def witness(n, j):
    rows = [sum(1 << v for v in range(u + 1, n)) for u in range(n)]
    if j:
        end = 2 * j
        rows[0] ^= 1 << end
        rows[end] |= 1
    return tuple(rows)


def predicted(n, j):
    result = [0] * (n + 1)
    result[n] = 1
    if j:
        for k in range(1, 2 * j, 2):
            result[n - k - 1] += 2 * comb(2 * j - 1, k)
    return tuple(result)


def gate_step(row, upper):
    out = {}
    for j, weight in row.items():
        if upper:
            entries = ((j - 1, -2 * j * (j*j - j + 1)),
                       (j, j**4), (j + 1, 2 * j * (j*j + j + 1)))
        else:
            entries = ((j - 1, -j * (j - 1) * (2*j - 1)),
                       (j, j*j * (j*j - 2)), (j + 1, j * (j + 1) * (2*j + 1)))
        for k, coeff in entries:
            if k > 0 and coeff:
                out[k] = out.get(k, 0) + weight * coeff
    return {k: v for k, v in out.items() if v}


def gate_ranks():
    checked = 0
    for n in range(1, 13):
        exps = list(range(1 if n % 2 else 2, n + 1, 2))
        d = len(exps)
        layer = [{1: 1}]
        all_rows = layer[:]
        repeated = layer[:]
        for depth in range(7):
            matrix = [[sum(weight * j**e for j, weight in row.items()) for e in exps]
                      for row in all_rows]
            single_word = [[sum(weight * j**e for j, weight in row.items()) for e in exps]
                           for row in repeated]
            require(rank(matrix) == rank(single_word) == min(depth + 1, d), ('gate rank', n, depth))
            require(rank([row[:-1] for row in matrix]) == rank([row[:-1] for row in single_word]) == min(depth + 1, d - 1),
                    ('gate affine rank', n, depth))
            checked += 1
            layer = [gate_step(row, upper) for row in layer for upper in (False, True)]
            all_rows.extend(layer)
            repeated.append(gate_step(repeated[-1], False))
    return checked


def main():
    census = []
    probes = 0
    for r, expected in enumerate(((1,), (3, 1), (9, 10, 2), (31, 83, 50, 6), (117, 684, 908, 336, 24))):
        require(skeleton(r) == expected, ('skeleton', r))
    for n in range(1, 6):
        vectors = []
        d = (n + 1) // 2
        for code in range(1 << (n * (n - 1) // 2)):
            rows = tournament(n, code)
            q = profile(rows)
            coefficients = poly(q)
            require(coefficients[0] == 0 and coefficients[-1] == 1, 'monic')
            require(all(not c for i, c in enumerate(coefficients) if i % 2 != n % 2), 'parity')
            observations = []
            for r in range(d - 1):
                actual = path_counts(wrap(rows, r))[-1]
                require(actual == sum(c * v for c, v in zip(skeleton(r), q)), 'actual context')
                observations.append(actual)
                probes += 1
            require(recover(n, observations) == coefficients, 'affine recovery')
            vectors.append(q)
        linear = rank(vectors)
        affine = rank([[a - b for a, b in zip(v, vectors[0])] for v in vectors])
        require((linear, affine) == (d, d - 1), ('ranks', n))
        census.append({'n': n, 'labelled': len(vectors), 'rank': linear, 'affine_rank': affine})
    basis_receipt = []
    for n in range(1, 12):
        basis = []
        d = (n + 1) // 2
        for j in range(d):
            q = profile(witness(n, j))
            coefficients = poly(q)
            require(coefficients == predicted(n, j), ('backedge identity', n, j))
            observations = [sum(c * v for c, v in zip(skeleton(r), q)) for r in range(d - 1)]
            require(recover(n, observations) == coefficients, ('basis recovery', n, j))
            basis.append(coefficients)
        require(rank(basis) == d, ('basis rank', n))
        basis_receipt.append([n, d])
    # Inherited hostile: same H, distinct next coordinate, realized wrapper separation.
    a, b = (16, 9, 3, 5, 14), (24, 17, 3, 6, 12)
    qa, qb = profile(a), profile(b)
    require(qa[0] == qb[0] == 15 and (qa[1], qb[1]) == (78, 90), 'scalar hostile')
    require((path_counts(wrap(a, 1))[-1], path_counts(wrap(b, 1))[-1]) == (123, 135), 'wrapper hostile')
    # Boundary: a reversed two-vertex tournament still has only one Hamilton path.
    require(profile((0, 1)) == (1, 2), 'two-vertex boundary')
    gate_universes = gate_ranks()
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == SOURCE, 'source changed during run')
    print(json.dumps({'status': 'PASS', 'source_sha256': SOURCE, 'checks': checks,
                      'census': census, 'actual_probe_evaluations': probes,
                      'basis_ranks_through_11': basis_receipt,
                      'scalar_hostile': [123, 135], 'gate_order_depth_universes': gate_universes}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
