#!/usr/bin/env python3
"""Full marked joint guard, with a literal omitted-exit path witness.

Standalone finite controls. No producer imports or response-table inputs.
"""
from itertools import combinations
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


def graph(n, code):
    rows = [0] * n
    for i, (u, v) in enumerate(combinations(range(n), 2)):
        if code >> i & 1:
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return tuple(rows)


def tables(rows):
    """Enumerate directed simple paths by DFS, retaining both endpoints."""
    n, size = len(rows), 1 << len(rows)
    P = [0] * size
    P[0] = 1
    A = [[0] * size for _ in rows]
    E = [[0] * size for _ in rows]
    Q = [[[0] * size for _ in rows] for _ in rows]
    def visit(start, end, mask):
        P[mask] += 1
        A[start][mask] += 1
        E[end][mask] += 1
        Q[start][end][mask] += 1
        for v in range(n):
            if not (mask >> v & 1) and rows[end] >> v & 1:
                visit(start, v, mask | (1 << v))
    for v in range(n):
        visit(v, v, 1 << v)
    return P, A, E, Q


def conv(A, B):
    C = [0] * len(A)
    for m in range(len(A)):
        s = m
        while True:
            C[m] += A[s] * B[m ^ s]
            if s == 0:
                break
            s = (s - 1) & m
    return C


def coefficient(A, B, mask):
    s, out = mask, 0
    while True:
        out += A[s] * B[mask ^ s]
        if s == 0:
            return out
        s = (s - 1) & mask


def shortest(rows, a, w):
    todo = [(a,)]
    seen = {a}
    for path in todo:
        if path[-1] == w:
            return path
        for v in range(len(rows)):
            if v not in seen and rows[path[-1]] >> v & 1:
                seen.add(v)
                todo.append(path + (v,))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = ap.parse_args()
    census, positive, wrong_exit_zero, scalar_two_exit = [], [], None, None
    for n in range(1, 6):
        graphs = strongs = cases = witnesses = fullguards = 0
        for code in range(1 << (n * (n - 1) // 2)):
            rows = graph(n, code)
            graphs += 1
            P, A, E, Q = tables(rows)
            full = (1 << n) - 1
            paths = {(a, w): shortest(rows, a, w) for a in range(n) for w in range(n)}
            strong = n >= 3 and all(paths.values())
            strongs += strong
            P2 = conv(P, P)
            WP = [conv(E[w], P) for w in range(n)]
            F = [coefficient(E[w], P, full) for w in range(n)]
            for a in range(n):
                exits = [b for b in range(n) if rows[a] >> b & 1]
                EB = [sum(E[b][s] for b in exits) for s in range(full + 1)]
                EBA = conv(EB, A[a])
                EBP = conv(EB, P)
                sigma = coefficient(P, A[a], full) - P[full]
                for w in range(n):
                    cases += 1
                    G = coefficient(WP[w], A[a], full) + coefficient(P2, Q[a][w], full)
                    tau = G - F[w]
                    debt = coefficient(EBA, E[w], full) + coefficient(EBP, Q[a][w], full)
                    need(tau == debt, 'exact failed-join decomposition, including coincident marks')
                    path = paths[a, w]
                    if path and ((a != w and len(exits) >= 2) or (a == w and exits)):
                        omitted = [b for b in exits if b not in path]
                        need(bool(omitted), 'shortest-path omitted exit')
                        b = omitted[0]
                        used = sum(1 << v for v in path) | (1 << b)
                        need(P[full ^ used] >= 1 and Q[a][w][used ^ (1 << b)] >= 1,
                             'literal positive K-surplus cover')
                        need(tau >= 1, 'strict joint defect from explicit cover')
                        witnesses += 1
                    if strong:
                        need(tau >= len(exits) - (a != w), 'all-strong outdegree lower bound')
                        right = (a != w and exits == [w] and E[w][full ^ (1 << a)] > 0)
                        need((tau == 0 and sigma > 0) == right, 'complete positive joint guard locus')
                        if tau == 0 and sigma in (1, 2):
                            fullguards += 1
                            if not any(r['sigma'] == sigma for r in positive):
                                positive.append({'rows': rows, 'a': a, 'w': w, 'H': P[full], 'F': F[w], 'sigma': sigma, 'tau': tau})
                        if tau == 0 and len(exits) == 1 and exits[0] != w and wrong_exit_zero is None:
                            wrong_exit_zero = {'rows': rows, 'a': a, 'w': w, 'exit': exits[0], 'sigma': sigma, 'tau': tau}
                        if sigma == 2 and len(exits) == 2 and scalar_two_exit is None:
                            scalar_two_exit = {'rows': rows, 'a': a, 'w': w, 'sigma': sigma, 'tau': tau}
        need(cases == n*n*2 ** (n*(n-1)//2), 'full labelled marked universe')
        census.append({'n': n, 'graphs': graphs, 'strong_graphs': strongs, 'marked_cases': cases,
                       'positive_cover_witnesses': witnesses, 'digit_guards': fullguards})
    need({r['sigma'] for r in positive} == {1,2}, 'both actual positive guards')
    need(wrong_exit_zero is not None and wrong_exit_zero['sigma'] == 0, 'positivity hypothesis indispensable')
    need(scalar_two_exit is not None and scalar_two_exit['tau'] > 0, 'scalar two-exit ready hostile')
    # Literal C2 at a=0 of (2,4,1), retaining old root w=1. The new
    # transported auxiliary3 exits to old0, not to retained1.
    transported = (18, 28, 25, 1, 8)
    P, A, E, Q = tables(transported)
    full, a, w = 31, 3, 1
    F = coefficient(E[w], P, full)
    G = coefficient(conv(E[w], P), A[a], full) + coefficient(conv(P, P), Q[a][w], full)
    need(coefficient(P, A[a], full) == P[full] and G == F, 'actual transported zero auxiliary')
    need(transported[a] == 1 and a != w, 'transported exit differs from retained root')
    result = {'scope': 'complete all-strong theorem has an analytic proof; finite controls cover all labelled tournaments through5 and all ordered mark pairs',
              'gates': GATES, 'census': census, 'positive': positive,
              'wrong_exit_zero': wrong_exit_zero, 'scalar_two_exit': scalar_two_exit,
              'transported_zero': {'rows': transported, 'a': a, 'w': w, 'exit':0, 'H':P[full], 'F':F, 'sigma':0, 'tau':0}}
    raw = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    args.output.write_bytes(raw)
    print(json.dumps({'gates': GATES, 'cases': sum(r['marked_cases'] for r in census), 'bytes':len(raw),
                      'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))


if __name__ == '__main__':
    main()
