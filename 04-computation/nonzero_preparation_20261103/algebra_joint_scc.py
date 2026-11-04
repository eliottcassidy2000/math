#!/usr/bin/env python3
"""Independent terminal-SCC and full joint-zero sidecar.

Uses literal permutations, not the DFS counter in algebra_joint_guard.py.
No imports of producer code or receipts.
"""
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json

GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def build(n, code):
    rows = [0] * n
    for i, (u, v) in enumerate(combinations(range(n), 2)):
        if code >> i & 1:
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return tuple(rows)


def reachability(rows, mask):
    reach = [(row | (1 << i)) & mask if mask >> i & 1 else 0 for i, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    terminal = {v for v in range(len(rows)) if mask >> v & 1 and
                all(reach[u] >> v & 1 for u in range(len(rows)) if mask >> u & 1)}
    initial = {v for v in range(len(rows)) if mask >> v & 1 and reach[v] == mask}
    return initial, terminal


def paths(rows):
    n, size = len(rows), 1 << len(rows)
    P = [0] * size
    P[0] = 1
    A = [[0] * size for _ in rows]
    E = [[0] * size for _ in rows]
    Q = [[[0] * size for _ in rows] for _ in rows]
    for k in range(1, n + 1):
        for path in permutations(range(n), k):
            if all(rows[u] >> v & 1 for u, v in zip(path, path[1:])):
                mask = sum(1 << v for v in path)
                P[mask] += 1
                A[path[0]][mask] += 1
                E[path[-1]][mask] += 1
                Q[path[0]][path[-1]][mask] += 1
    return P, A, E, Q


def multiply(A, B):
    C = [0] * len(A)
    for mask in range(len(A)):
        s = mask
        while True:
            C[mask] += A[s] * B[mask ^ s]
            if not s:
                break
            s = (s - 1) & mask
    return C


def coeff(A, B, mask):
    total, s = 0, mask
    while True:
        total += A[s] * B[mask ^ s]
        if not s:
            return total
        s = (s - 1) & mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = ap.parse_args()
    pair_cases = marked_cases = sole_exit_cases = 0
    census = []
    for n in range(2, 6):
        pcount = mcount = 0
        for code in range(1 << (n * (n - 1) // 2)):
            rows = build(n, code)
            full = (1 << n) - 1
            initial, terminal = reachability(rows, full)
            P, A, E, Q = paths(rows)
            for b, w in combinations(range(n), 2):
                J = coeff(E[b], E[w], full)
                need((J > 0) == bool({b, w} & terminal), 'two ending paths iff a terminal endpoint')
                pcount += 1
            if n < 3 or len(terminal) != n:
                continue
            P2 = multiply(P, P)
            positive_by_root = [0] * n
            for a in range(n):
                exits = [v for v in range(n) if rows[a] >> v & 1]
                sigma = coeff(P, A[a], full) - P[full]
                U = full ^ (1 << a)
                I, C = reachability(rows, U)
                if len(exits) == 1:
                    need(exits[0] in I, 'strong sole-exit port lies in deletion initial SCC')
                    need((sigma > 0) == (len(C) == n - 1), 'positive one-exit defect iff strong deletion')
                    sole_exit_cases += 1
                for w in range(n):
                    F = coeff(E[w], P, full)
                    G = coeff(multiply(E[w], P), A[a], full) + coeff(P2, Q[a][w], full)
                    tau = G - F
                    predicted = (a != w and len(exits) == 1 and
                                 ((len(C) == n - 1 and w == exits[0]) or
                                  (len(C) < n - 1 and w not in C)))
                    need((tau == 0) == predicted, 'complete SCC classification including all other degrees and a=w')
                    positive_by_root[w] += tau == 0 and sigma > 0
                    mcount += 1
            need(all(v <= 1 for v in positive_by_root), 'unique positive auxiliary at each fixed root')
        pair_cases += pcount
        marked_cases += mcount
        census.append({'n': n, 'ending_pairs': pcount, 'strong_marked_pairs': mcount})
    need(pair_cases == 10650 and marked_cases == 14002, 'complete finite universes')
    result = {'gates': GATES, 'ending_pairs': pair_cases, 'strong_marked_pairs': marked_cases,
              'sole_exit_auxiliaries': sole_exit_cases, 'census': census,
              'scope': 'finite controls of analytic all-order terminal-SCC and unique-positive-auxiliary corollaries'}
    raw = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    args.output.write_bytes(raw)
    print(json.dumps({'gates':GATES,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))


if __name__ == '__main__':
    main()
