#!/usr/bin/env python3
"""Same-agent separate audit: literal path permutations and dense F3 matrices.

Does not import check.py or any inherited research implementation. Reads its
output only as candidate data, then verifies orders and minimality directly.
"""
from itertools import combinations, permutations, product
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def matmul(a, b):
    cols = list(zip(*b))
    return [[sum(x*y for x, y in zip(row, col)) % 3 for col in cols] for row in a]


def identity(d):
    return [[int(i == j) for j in range(d)] for i in range(d)]


def mpow(a, n):
    out = identity(len(a))
    while n:
        if n % 2:
            out = matmul(out, a)
        a = matmul(a, a)
        n //= 2
    return out


def rank(a):
    a = [list(row) for row in a]
    r = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(r, len(a)) if a[j][col] % 3), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        if a[r][col] % 3 == 2:
            a[r] = [2*x % 3 for x in a[r]]
        for j in range(r+1, len(a)):
            q = a[j][col] % 3
            a[j] = [(x-q*y) % 3 for x, y in zip(a[j], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def matrix(n):
    out = [[0]*(2*n) for _ in range(2*n)]
    for i in range(n):
        out[i][n+i] = 1
        out[n+i][i] = 1
        out[n+i][n+i] += 1
        out[n+i][n+(i+1) % n] += 1
    return out


def primes(n):
    d, answer = 2, []
    while d*d <= n:
        if n % d == 0:
            answer.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        answer.append(n)
    return answer


def count_literal(mask, cut):
    r = [[False]*7 for _ in range(7)]
    for k, (i, j) in enumerate(combinations(range(6), 2)):
        r[i][j] = bool(mask >> k & 1)
        r[j][i] = not r[i][j]
    for i in range(6):
        r[i][6] = bool(cut >> i & 1)
        r[6][i] = not r[i][6]
    return sum(all(r[p[i]][p[i+1]] for i in range(6)) for p in permutations(range(7)))


def main():
    candidate = json.loads((HERE/'results.json').read_text())
    seeds = (153, 8345, 8337, 8369, 177, 185)
    birth_counts = [count_literal(m, 2) for m in seeds]
    require(birth_counts == [109, 115, 101, 77, 73, 81], 'literal reader path enumeration')
    orders = []
    for row in candidate['arbitrary_widths']:
        n, order = row['n'], row['order']
        a, eye = matrix(n), identity(2*n)
        require(mpow(a, order) == eye, 'matrix order upper bound')
        for p in primes(order):
            require(mpow(a, order//p) != eye, 'matrix order minimality')
        # Independently check cyclicity and one-site observability.
        power, observable, reachable = eye, [], []
        for _ in range(2*n):
            observable.append(power[0])
            reachable.append([row[0] for row in power])
            power = matmul(power, a)
        require(rank(observable) == rank(reachable) == 2*n, 'both cyclic vectors')
        orders.append([n, order])
    fixed_spaces = []
    for n in (1, 3, 9, 27):
        a, eye = matrix(n), identity(2*n)
        for t in (1, 2, 4):
            power = mpow(a, t)
            require(rank([[x-y for x, y in zip(r, s)] for r, s in zip(power, eye)]) == 2*n,
                    'no nonzero shorter residue period')
        q = 1
        while q <= n:
            power = mpow(a, 8*q)
            nullity = 2*n-rank([[x-y for x, y in zip(r, s)] for r, s in zip(power, eye)])
            require(nullity == 2*q, 'all-state period subspace')
            fixed_spaces.append([n, 8*q, nullity])
            q *= 3
    # A second route to the history census: dense linear output rows, then
    # the proved signed carry identity, instead of the main gate simulator.
    history_census = []
    sigma = lambda x: ((x % 3)+1) % 3-1
    for n in (1, 2, 3, 4):
        a, power, obs = matrix(n), identity(2*n), []
        for _ in range(2*n+1):
            obs.append(power[0])
            power = matmul(power, a)
        histories = [set() for _ in range(2*n)]
        for v in product(range(3), repeat=2*n):
            ys = [sum(x*y for x, y in zip(row, v)) % 3 for row in obs]
            for bit in (0, 1):
                h = [ys[0]+3*bit]
                histories[0].add(tuple(h))
                for t in range(2*n-1):
                    h.append((h[-1]+sigma(ys[t+2]-ys[t])+sigma(ys[t+1])-sigma(ys[t+2])) % 6)
                    histories[t+1].add(tuple(h))
        counts = [len(h) for h in histories]
        expected = [6 if length == 1 else (2*9**n if length == 2*n else 10*3**(length-1))
                    for length in range(1, 2*n+1)]
        require(counts == expected, 'carry reconstruction census')
        history_census.append({'n': n, 'histories_by_length': counts})
    require(history_census == candidate['physical_boundary_histories'], 'two history implementations')
    result = {'auditor': 'same agent; distinct permutation/matrix algorithms; no research imports',
              'birth_reader_permutation_counts': birth_counts, 'matrix_orders': orders,
              'fixed_space_dimensions': fixed_spaces,
              'carry_history_census': history_census,
              'boundaries': 'does not independently reprove the inherited seed cover or gate-lift theorem'}
    (HERE/'audit_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
