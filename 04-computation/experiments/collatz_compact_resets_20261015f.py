#!/usr/bin/env python3
"""Exact first-reset packets with compact native graph and integer clocks.

No ROOT oracle, guessed sign, or floating-point acceptance. The graph only
contains explicitly observed shortcut edges. Compression remembers full clocks.
"""
import argparse
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def need(test, message):
    if not test:
        raise AssertionError(message)


def T(n):
    return (3*n+1)//2 if n & 1 else n//2


def ceilroot(n, d):
    v = math.isqrt(n//d)
    return v if d*v*v == n else v+1


class Graph:
    def __init__(self, M):
        self.M = M
        self.seen = np.zeros(M+1, dtype=np.bool_)
        self.seen[2::2] = True
        self.outside = set()
        self.edge_count = M//2
        self.endpoint = np.arange(M+1, dtype=np.uint64)
        self.clock = np.zeros(M+1, dtype=np.uint32)
        power = 2
        while power <= M:
            self.endpoint[power::power] //= 2
            self.clock[power::power] += 1
            power *= 2

    def has(self, n):
        return bool(self.seen[n]) if n <= self.M else n in self.outside

    def observe(self, n):
        need(n != 1, 'no edge out of stopped ROOT')
        if not self.has(n):
            if n <= self.M:
                self.seen[n] = True
            else:
                self.outside.add(n)
            self.edge_count += 1

    def extend(self, cols, resets=1):
        histogram = Counter()
        for source in cols:
            x = int(source)
            if x == 1:
                continue
            need(x & 1, 'every frontier terminal is odd')
            for _ in range(resets):
                ones = 0
                while x != 1:
                    self.observe(x)
                    x = T(x)
                    a = 1
                    while not (x & 1):
                        self.observe(x)
                        x //= 2
                        a += 1
                    if a >= 2:
                        histogram[ones] += 1
                        break
                    ones += 1
        return dict(sorted(histogram.items()))

    def close(self, cols, inverse):
        # Old source paths remain valid because edges are only added. They are
        # shortcuts with their actual lengths, not fresh inferred graph edges.
        memo_end = np.zeros(self.M+1, dtype=np.uint64)
        memo_clock = np.zeros(self.M+1, dtype=np.uint32)
        memo_end[1] = 1
        memo_out = {}

        def get(n):
            if n <= self.M:
                e = int(memo_end[n])
                return (e, int(memo_clock[n])) if e else None
            return memo_out.get(n)

        def put(n, e, d):
            need(e < 2**64 and d < 2**32, 'memo integer storage bounds')
            if n <= self.M:
                memo_end[n], memo_clock[n] = e, d
            else:
                memo_out[n] = (e, d)

        target = np.zeros(len(cols), dtype=np.uint64)
        extra = np.zeros(len(cols), dtype=np.uint32)
        for ci, source in enumerate(cols):
            x = int(source)
            path, where = [], set()
            while get(x) is None:
                need(x not in where, 'no witnessed nonROOT cycle')
                where.add(x)
                if not self.has(x):
                    put(x, x, 0)
                    break
                if x <= self.M and self.clock[x]:
                    y, cost = int(self.endpoint[x]), int(self.clock[x])
                else:
                    y, cost = T(x), 1
                path.append((x, cost))
                x = y
            e, d = get(x)
            for x, cost in reversed(path):
                d += cost
                put(x, e, d)
            target[ci], extra[ci] = get(int(source))
        need(int(extra.max())+int(self.clock.max()) < 2**32,
             'checked integer storage bounds')
        self.endpoint[1:] = target[inverse]
        self.clock[1:] += extra[inverse]

    def replay(self, sources):
        for n in sources:
            x = n
            for _ in range(int(self.clock[n])):
                need(x != 1 and self.has(x), 'replayed edge is observed, ROOT-stopped')
                x = T(x)
            need(x == int(self.endpoint[n]), 'literal endpoint and integer clock')

    def internal_clock_audit(self):
        checked = 0
        for lo in range(2, self.M+1, 1 << 19):
            n = np.arange(lo, min(self.M+1, lo+(1 << 19)), dtype=np.uint64)
            y = np.where(n & 1, (3*n+1)//2, n//2)
            good = y <= self.M
            n, y = n[good], y[good]
            need(bool(np.all(self.seen[n])), 'every original source edge was observed')
            need(bool(np.all(self.endpoint[n] == self.endpoint[y])), 'all internal edge endpoints')
            need(bool(np.all(self.clock[n] == self.clock[y]+1)), 'all internal exact clock increments')
            checked += len(n)
        return checked


def row_norms(k, graph, inverse, column_count, matrix=False):
    X, b = 1 << k, 1 << (k//2)
    need(32*b*X < 2**63, 'signed int64 event/row bound')
    # inverse[n-1] corresponds to original source n; zero is only a sentinel.
    target = np.empty(8*X+1, dtype=np.int64)
    target[0], target[1:] = 0, inverse
    value = np.full(8*X+1, -1, dtype=np.int64)
    value[2::3] = 2
    value[(graph.clock & 1) != 0] *= -1
    value[0] = 0
    delta = np.zeros(column_count, dtype=np.int64)
    np.add.at(delta, target[X+1:4*X+1], value[X+1:4*X+1])
    norms, matrices = [], []
    for a in range(X, 2*X, b):
        row = b*delta
        nn = np.arange(a, a+b, dtype=np.int64)
        remaining = a+b-1-nn
        for j in (1, 2, 3, 4):
            ids = 4*nn+j
            np.add.at(row, target[ids], value[ids]*remaining)
            np.add.at(delta, target[ids], value[ids])
        ids = nn+1
        np.add.at(row, target[ids], -value[ids]*remaining)
        np.add.at(delta, target[ids], -value[ids])
        norms.append(int(np.sum(np.abs(row))))
        if matrix:
            matrices.append(row)
    return norms, matrices


def run(k, verify=False, resets=1, full=False):
    X, b = 1 << k, 1 << (k//2)
    graph = Graph(8*X)
    cols, inverse = np.unique(graph.endpoint[1:], return_inverse=True)
    # Halving-only coefficients have constant sign at each odd terminal as N
    # varies. Therefore the norm of a bin is the sum of the exact L_0(N).
    nn = np.arange(X, 2*X, dtype=np.int64)
    l0 = 4*nn-3*((2*nn+3)//6)+(4*nn+1)//6-(2*nn+1)//6
    norms = [int(v) for v in l0.reshape(-1, b).sum(axis=1)]
    if verify:
        direct_norms, _ = row_norms(k, graph, inverse, len(cols))
        need(norms == direct_norms, 'halving closed form agrees with event sweep')
    del nn, l0
    initial = sum(v*v for v in norms)
    need(initial <= 64*b*X**3, 'inherited initial envelope')
    bound = F(ceilroot(initial, b))
    stages = []
    for j in range(1, k+9):
        histogram = graph.extend(cols, resets)
        graph.close(cols, inverse)
        cols, inverse = np.unique(graph.endpoint[1:], return_inverse=True)
        newer, matrix = row_norms(k, graph, inverse, len(cols), verify)
        need(all(v <= u for u, v in zip(norms, newer)), 'native row monotonicity')
        rs = [max(0, 10*v-7*u) for u, v in zip(norms, newer)]
        r2, q2 = sum(v*v for v in rs), sum(v*v for v in newer)
        bound = F(7, 10)*bound+F(3, 10)*ceilroot(r2, 9*b)
        # Repeatable independent literal-edge replay, including fine annular
        # endpoints and a source grid. Every source retains its entire clock.
        graph.replay(sorted({1, X+1, 2*X-1, 4*X+3, 8*X} |
                            set(range(1, 8*X+1, max(1, X//32)))))
        checked_clocks = graph.internal_clock_audit()
        if verify:
            from collatz_gram_slack_20261015b import close_graph
            from collatz_odd_ray_overlap_20261015d import rows
            edges = {n: T(n) for n in range(2, 8*X+1) if graph.has(n)}
            edges.update({n: T(n) for n in graph.outside})
            cache, cycles = close_graph(edges, range(1, 8*X+1))
            need(not cycles, 'independent graph closure cycle audit')
            need(all(cache[n] == (int(graph.endpoint[n]),
                                  -1 if int(graph.clock[n]) & 1 else 1,
                                  int(graph.clock[n])) for n in cache if n <= 8*X),
                 'independent full-source clock closure')
            direct, direct_cols = rows(k, cache)
            pos = {e: i for i, e in enumerate(direct_cols)}
            need(np.array_equal(np.array(matrix), direct[:, [pos[int(e)] for e in cols]]),
                 'independent original c_I(n) coefficients')
        record = dict(stage=j, residual_square=[r2, 9*b],
                      row_envelope=[q2, b], residual_over_X=math.sqrt(r2/(9*b))/X,
                      Q_over_X=math.sqrt(q2/b)/X, pass_seven=r2 <= 441*b*X*X,
                      row_norms=newer, residual_numerators=rs, terminals=len(cols),
                      edges=graph.edge_count, outside_edges=len(graph.outside),
                      unresolved_sources=int(np.count_nonzero(graph.endpoint[1:] != 1)),
                      internal_clock_edges_checked=checked_clocks,
                      max_native_clock=int(graph.clock.max()), reset_ones_histogram=histogram,
                      discounted_upper_bound=[bound.numerator, bound.denominator])
        stages.append(record)
        print(json.dumps({key: record[key] for key in ('stage', 'residual_over_X',
                         'Q_over_X', 'terminals', 'max_native_clock')} | {'k': k}), flush=True)
        norms = newer
        if q2 < 49*b*X*X and not full:
            break
    result = dict(k=k, resets_per_stage=resets, initial_envelope=[initial, b],
                  source_interval=[1, 8*X], bin_width=b, stages=stages,
                  identity_padding=k+8-len(stages), all_stage_residuals_below_7X=
                  all(s['pass_seven'] for s in stages),
                  discounted_bound_below_7p5X=2*bound < 15*X,
                  independent_full_source_check=verify)
    data = json.dumps(result, indent=2, sort_keys=True)+'\n'
    path = ROOT/f'05-knowledge/results/collatz_compact_resets_20261015f_k{k}_r{resets}.json'
    path.write_text(data)
    print(json.dumps(dict(output=str(path), sha256=hashlib.sha256(data.encode()).hexdigest())), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, nargs='+', default=[8, 9, 19, 20])
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--resets', type=int, default=1)
    parser.add_argument('--full', action='store_true')
    args = parser.parse_args()
    need(args.resets > 0, 'positive reset count')
    for k in args.k:
        run(k, args.verify, args.resets, args.full)
