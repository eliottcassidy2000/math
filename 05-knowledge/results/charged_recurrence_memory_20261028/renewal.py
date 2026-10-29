"""Actual fixed-path complement: class compression, leakage and first-return memory."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json
import charged
from memory import rank_mod

ports = charged.paley.ports
oracle = ports.oracle
HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def presentation(n, bits):
    rows = [0] * n
    free = 0
    for i in range(n):
        for j in range(i + 1, n):
            forward = True if j == i + 1 else bool(bits >> free & 1)
            if j > i + 1:
                free += 1
            a, b = (i, j) if forward else (j, i)
            rows[a] |= 1 << b
    return tuple(rows)


def classify_all(n):
    buckets = defaultdict(list)
    reps, ids = [], []
    for x in range(1 << ((n - 1) * (n - 2) // 2)):
        T = presentation(n, x)
        signature = oracle.colours(T)[1]
        found = next((i for i in buckets[signature] if oracle.isomorphisms(T, reps[i])), None)
        if found is None:
            found = len(reps)
            buckets[signature].append(found)
            reps.append(T)
        ids.append(found)
    return reps, ids


def multiply(a, b, mod=None):
    out = [[0] * len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for k, v in enumerate(row):
            if v:
                for j, w in enumerate(b[k]):
                    if w:
                        out[i][j] += v * w
        if mod:
            out[i] = [x % mod for x in out[i]]
    return out


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def subtract(a, b, mod=None):
    return [[(x - y) % mod if mod else x - y for x, y in zip(u, v)] for u, v in zip(a, b)]


def components(K):
    seen = set()
    c = b = 0
    for start in range(len(K)):
        if start in seen:
            continue
        c += 1
        colour = {start: 0}
        todo = [start]
        bipartite = True
        while todo:
            i = todo.pop()
            seen.add(i)
            for j, value in enumerate(K[i]):
                if not value:
                    continue
                if j in colour:
                    bipartite &= colour[j] != colour[i]
                else:
                    colour[j] = 1 - colour[i]
                    todo.append(j)
        b += int(bipartite)
    return c, b


def census(n):
    reps, ids = classify_all(n)
    r, N = len(reps), len(ids)
    need(r == oracle.burnside(n), 'all ordinary tournament classes recovered from presentations')
    mass = Counter(ids)
    K = [[0] * r for _ in range(r)]
    for x, i in enumerate(ids):
        K[i][ids[(N - 1) ^ x]] += 1
    need(K == list(map(list, zip(*K))), 'actual complement incidence symmetry')
    for i, T in enumerate(reps):
        need(mass[i] * len(oracle.isomorphisms(T, T, True)) == oracle.hamilton(T),
             'presentation mass equals Hamilton paths divided by automorphisms')
    mod = 1000000007
    M = [[v * pow(mass[i], -1, mod) % mod for v in row] for i, row in enumerate(K)]
    I = identity(r)
    leakage = subtract(I, multiply(M, M, mod), mod)
    c, b = components(K)
    rank = rank_mod(leakage, mod)
    need(rank == r - c - b, 'exact leakage rank via graph upper bound and modular lower bound')
    pairs = [(i, j) for i in range(r) for j in range(r) if K[i][j]]
    incidence = [[int(a == i) for a in range(r)] + [int(a == j) for a in range(r)] for i, j in pairs]
    linear = rank_mod(incidence, mod)
    need(linear == 2 * r - c - b, 'minimal invariant linear closure rank')
    kernels = [leakage]
    for _ in range(7):
        kernels.append([[(-x) % mod for x in row] for row in multiply(M, kernels[-1], mod)])
    F = [I]
    for t in range(8):
        value = multiply(M, F[t], mod)
        for j in range(t):
            term = multiply(kernels[j], F[t - 1 - j], mod)
            value = [[(a + b) % mod for a, b in zip(u, v)] for u, v in zip(value, term)]
        expected = M if (t + 1) % 2 else I
        need(value == expected, 'memory recurrence recovers retained-path complement powers')
        F.append(value)
    # Direct rational hidden-state iteration, independent of closed kernel formula.
    if n == 4:
        rational = [[Q(v, mass[i]) for v in row] for i, row in enumerate(K)]
        def project(V):
            sums = [[sum(V[x][j] for x in range(N) if ids[x] == i) / mass[i]
                     for j in range(r)] for i in range(r)]
            return [sums[ids[x]][:] for x in range(N)]
        J = [[Q(int(ids[x] == j)) for j in range(r)] for x in range(N)]
        UJ = [J[(N - 1) ^ x] for x in range(N)]
        hidden = subtract(UJ, project(UJ))
        expected = subtract(I, multiply(rational, rational))
        for _ in range(6):
            shifted = [hidden[(N - 1) ^ x] for x in range(N)]
            observed = project(shifted)
            actual = [observed[ids.index(i)] for i in range(r)]
            need(actual == expected, 'literal hidden transition kernel versus formula over rationals')
            hidden = subtract(shifted, observed)
            expected = [[-x for x in row] for row in multiply(rational, expected)]
    witness = next(({'class': i, 'presentations': [x, y],
                     'complement_classes': [ids[(N - 1) ^ x], ids[(N - 1) ^ y]]}
                    for i in range(r) for x in range(N) if ids[x] == i
                    for y in range(x + 1, N) if ids[y] == i and ids[(N - 1) ^ x] != ids[(N - 1) ^ y]), None)
    return {'order': n, 'presentations': N, 'ordinary_classes': r, 'components': c,
            'bipartite_components': b, 'leakage_rank': rank,
            'minimal_complement_only_linear_dimension': linear,
            'minimal_complement_only_set_carrier': len(pairs),
            'fibre_masses': [mass[i] for i in range(r)], 'incidence': K,
            'nonlumpability_witness': witness}


def main(write=True):
    global CHECKS
    CHECKS = 0
    rows = [census(n) for n in range(3, 7)]
    left = ports.frames.p.controls.build_hostile(False)
    right = ports.frames.p.controls.build_hostile(True)
    def cyclic(T):
        return [v for v in combinations(range(len(T)), 3)
                if all((T[v[i]] >> v[(i + 1) % 3]) & 1 for i in range(3))
                or all((T[v[(i + 1) % 3]] >> v[i]) & 1 for i in range(3))]
    need(cyclic(left) == cyclic(right), '81/27 identical labelled cyclic triples')
    auts = [len(oracle.isomorphisms(T, T, True)) for T in (left, right)]
    need(auts == [81, 27], '81/27 full transporter control')
    out = {'status': 'FINITE-EXACT ordinary-class universe n=3..6; scoped complement-only closure',
           'checks': CHECKS, 'rows': rows, 'hostile_automorphism_orders': auts,
           'hostile_labelled_cyclic_triples': len(cyclic(left))}
    if write:
        (HERE / 'renewal.json').write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'checks': CHECKS, 'rows': [(r['order'], r['ordinary_classes'], r['leakage_rank'],
                      r['minimal_complement_only_linear_dimension'], r['minimal_complement_only_set_carrier']) for r in rows],
                     'hostile': auts}))
    return out


if __name__ == '__main__':
    main()
