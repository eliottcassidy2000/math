#!/usr/bin/env python3
"""Two-polarity endpoint covers: exact polynomial carrier and small controls.

No repository imports. The only inherited receipt read supplies old query words.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import json


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def skeletons(k):
    parts = [(a, b) for a in range(-1, k) for b in range(-1, k)
             if (a, b) != (-1, -1)]
    result = []
    for size in range(k+1):
        for subset in combinations(parts, size):
            used = [v for part in subset for v in set(part) if v >= 0]
            if len(used) == len(set(used)):
                result.append(tuple(sorted(subset)))
    return tuple(sorted(result))


def split(state):
    return tuple(p for p in state if p != (-1, -1)), state.count((-1, -1))


def kernel(skeleton, count, k, polarity, mask, discard):
    state = list(skeleton) + [(-1, -1)]*count
    incoming = [None] + [i for i, (_, e) in enumerate(state)
                        if (polarity == 1 if e == -1 else bool(mask >> e & 1))]
    outgoing = [None] + [i for i, (s, _) in enumerate(state)
                        if (polarity == 0 if s == -1 else not (mask >> s & 1))]
    rename = {i: i for i in range(k)}
    if discard < k:
        del rename[discard]
        rename[k] = discard
    answer = Counter()
    for i in incoming:
        for j in outgoing:
            if i is not None and i == j:
                continue
            parts = [p for t, p in enumerate(state) if t != i and t != j]
            parts.append((k if i is None else state[i][0],
                          k if j is None else state[j][1]))
            parts = tuple(sorted((rename.get(a, -1), rename.get(b, -1)) for a, b in parts))
            answer[split(parts)] += 1
    return answer


def operator(k, command):
    result = []
    for state in skeletons(k):
        A = kernel(state, 0, k, *command)
        B = kernel(state, 1, k, *command)
        for (target, degree), value in A.items():
            B[target, degree+1] -= value
        B = Counter({key: value for key, value in B.items() if value})
        need(all(v > 0 for v in B.values()), 'nonnegative first derivative coefficients')
        result.append((A, B))
    return result


def apply(vector, k, op):
    result = Counter()
    for state, (A, B) in zip(skeletons(k), op):
        for (old, degree), value in vector.items():
            if old != state:
                continue
            for (target, power), weight in A.items():
                result[target, degree+power] += value*weight
            if degree:
                for (target, power), weight in B.items():
                    result[target, degree-1+power] += value*degree*weight
    return +result


def literal(T, ports):
    n = len(T)
    slot = {v: i for i, v in enumerate(ports)}
    covers = set()
    for word in permutations(range(n)):
        for cuts in range(1 << (n-1)):
            parts, start = [], 0
            for i in range(n-1):
                if cuts >> i & 1:
                    parts.append(word[start:i+1]); start = i+1
                elif not (T[word[i]] >> word[i+1] & 1):
                    break
            else:
                parts.append(word[start:]); covers.add(tuple(sorted(parts)))
    result = Counter()
    for cover in covers:
        state = tuple(sorted((slot.get(p[0], -1), slot.get(p[-1], -1)) for p in cover))
        result[split(state)] += 1
    return result


def birth(T, ports, polarity, mask, discard):
    n = len(T); slot = {v: i for i, v in enumerate(ports)}
    old_to_new = [bool(mask >> slot[i] & 1) if i in slot else polarity == 1 for i in range(n)]
    child = tuple(row | ((1 << n) if old_to_new[i] else 0) for i, row in enumerate(T))
    child += (sum((1 << i) for i in range(n) if not old_to_new[i]),)
    ports = list(ports)
    if discard < len(ports):
        ports[discard] = n
    return child, tuple(ports)


@lru_cache(None)
def H(T):
    n = len(T); dp = [[0]*n for _ in range(1 << n)]
    for i in range(n): dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            if count:
                choices = T[v] & ~mask
                while choices:
                    bit = choices & -choices; choices -= bit
                    dp[mask | bit][bit.bit_length()-1] += count
    return sum(dp[-1])


def scalar(vector):
    return sum(value for (state, count), value in vector.items() if len(state)+count == 1)


def transitive(n):
    return tuple(sum(1 << j for j in range(i+1, n)) for i in range(n))


def clone_cycle(m):
    # a=0 -> the transitive module -> b=1 -> a.
    return ((1 << (m+2))-4, 1) + tuple(2 + sum(1 << j for j in range(i+1, m+2))
                                             for i in range(2, m+2))


def stirling(n, k):
    return sum((-1)**(k-j)*__import__('math').comb(k, j)*j**n for j in range(k+1)) // factorial(k)


def solve(M, b):
    A = [list(map(Fraction, row))+[Fraction(y)] for row, y in zip(M, b)]
    for j in range(len(b)):
        pivot = next(i for i in range(j, len(b)) if A[i][j])
        A[j], A[pivot] = A[pivot], A[j]
        d = A[j][j]; A[j] = [x/d for x in A[j]]
        for i in range(len(b)):
            if i != j:
                d = A[i][j]; A[i] = [x-d*y for x, y in zip(A[i], A[j])]
    return [row[-1] for row in A]


def main():
    commands = list(product(range(2), range(4), range(3)))
    ops = {c: operator(2, c) for c in commands}
    identities = 0
    for command, op in ops.items():
        for state in skeletons(2):
            for m in range(6):
                need(apply(Counter({(state, m): 1}), 2, op) == kernel(state, m, 2, *command),
                     'universal multiplicity control')
                identities += 1
    sources = []
    for n in (2, 3):
        edges = list(combinations(range(n), 2))
        for bits in range(1 << len(edges)):
            rows = [0]*n
            for i, (a, b) in enumerate(edges):
                if bits >> i & 1: rows[a] |= 1 << b
                else: rows[b] |= 1 << a
            sources.append(tuple(rows))
    sources += [transitive(4), (10, 4, 1, 6)]
    replays = 0
    for T in sources:
        v = literal(T, (0, 1))
        need(scalar(v) == H(T), 'literal initial Hamilton count')
        for command in commands:
            U, ports = birth(T, (0, 1), *command)
            actual = literal(U, ports); expected = apply(v, 2, ops[command])
            need(actual == expected, 'literal complete path-cover update')
            need(scalar(expected) == H(U), 'independent Hamilton DP')
            replays += 1
    # Old seventeen-query equality survives arbitrary positive-polarity words:
    # the extra terminal singleton is an ordinal sink block.
    old = json.loads((Path(__file__).resolve().parents[1]/'adaptive_state_induction_20261103/live_port.json').read_text())
    X = []
    for n in (2, 3):
        values = []
        for word in old['rank']['suffixes']:
            T, ports = transitive(n), (0, 1)
            for letter in word:
                root, discard = divmod(letter, 3)
                T, ports = birth(T, ports, 0, 1 << root, discard)
            values.append(H(T))
        X.append(values)
    need(X[0] == X[1], 'same old seventeen queries')
    hostile = [H(birth(transitive(n), (0, 1), 1, 1, 2)[0]) for n in (2, 3)]
    need(hostile == [1, 3], 'negative-polarity hostile distinguishes sources')
    # Strong clone sources; alternating outdegree-one/indegree-one births.
    families = []
    for r in range(4):
        values = []
        for m in range(1, 7):
            T, ports = clone_cycle(m), (0,)
            for _ in range(r):
                T, ports = birth(T, ports, 1, 0, 0)
                T, ports = birth(T, ports, 0, 1, 0)
            values.append(H(T))
        d = r+2
        coeff = solve([[factorial(j)*stirling(m, j) for j in range(1, d+1)]
                       for m in range(1, d+1)], values[:d])
        need(all(x.denominator == 1 and x >= 0 for x in coeff) and coeff[-1] > 0,
             'positive top run coefficient')
        need(all(sum(c*factorial(j)*stirling(m, j) for j, c in enumerate(coeff, 1)) == h
                 for m, h in enumerate(values, 1)), 'remaining clone controls')
        families.append(dict(rounds=r, max_runs=d, run_coefficients=[int(x) for x in coeff],
                             counts_m1_through6=values))
    serialized = []
    states = skeletons(2); index = {state: i for i, state in enumerate(states)}
    for command, op in ops.items():
        blocks = []
        for i, (A, B) in enumerate(op):
            for derivative, matrix in enumerate((A, B)):
                for (state, degree), value in sorted(matrix.items()):
                    blocks.append([i, index[state], derivative, degree, value])
        serialized.append(dict(command=command, entries=blocks))
    print(json.dumps(dict(status='PASS', skeleton_counts=[len(skeletons(k)) for k in (1, 2)],
                         skeletons=states, operators=serialized, multiplicity_checks=identities,
                         literal_sources=len(sources), literal_transitions=replays,
                         old_query_twins=X, dual_birth_counts=hostile,
                         strong_clone_families=families), indent=2, sort_keys=True))


if __name__ == '__main__': main()
