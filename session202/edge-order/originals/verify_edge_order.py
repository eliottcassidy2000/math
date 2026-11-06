#!/usr/bin/env python3
import hashlib
import itertools
import json
from pathlib import Path


def need(p, msg):
    if not p:
        raise RuntimeError(msg)


def transitive(n):
    return [sum(1 << j for j in range(i + 1, n)) for i in range(n)]


def all_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        rows = [0] * n
        for k, (u, v) in enumerate(pairs):
            if bits & (1 << k):
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
        yield rows


def build(P, Q):
    # Boundary x=0, y=1, b=2, followed by P, then Q.
    p, q = len(P), len(Q)
    n = 3 + p + q
    rows = [0] * n
    def edge(a, b):
        rows[a] |= 1 << b
    edge(0, 1)
    edge(1, 2)
    edge(2, 0)
    for i in range(p):
        v = 3 + i
        edge(0, v); edge(v, 1); edge(v, 2)
        for j in range(p):
            if P[i] & (1 << j):
                edge(v, 3 + j)
    for i in range(q):
        v = 3 + p + i
        edge(0, v); edge(v, 1); edge(2, v)
        for j in range(q):
            if Q[i] & (1 << j):
                edge(v, 3 + p + j)
    for i in range(p):
        for j in range(q):
            edge(3 + i, 3 + p + j)
    return rows


def all_subset_counts(rows, force_forward=False, cycle=False):
    n, N = len(rows), 1 << len(rows)
    dp = [[0] * n for _ in range(N)]
    for v in range(n):
        if cycle and v != 0:
            continue
        if force_forward and v == 1:
            continue
        dp[1 << v][v] = 1
    for S in range(1, N):
        available = (N - 1) ^ S
        if force_forward and not (S & 1):
            available &= ~2
        for v, z in enumerate(dp[S]):
            if not z:
                continue
            choices = rows[v] & available
            while choices:
                bit = choices & -choices
                w = bit.bit_length() - 1
                dp[S | bit][w] += z
                choices -= bit
    if cycle:
        return sum(dp[-1][v] for v in range(n) if rows[v] & 1)
    out = [sum(row) for row in dp]
    out[0] = 1
    return out


def module_stats(rows):
    h = all_subset_counts(rows)
    S = (1 << len(rows)) - 1
    c = sum(h[A] * h[S ^ A] for A in range(S + 1))
    return h[S], c


def run():
    controls = []
    for p in range(1, 5):
        for q in range(1, 5):
            rows = build(transitive(p), transitive(q))
            H = all_subset_counts(rows)[-1]
            f = all_subset_counts(rows, force_forward=True)[-1]
            C = all_subset_counts(rows, cycle=True)
            need(f == 2**p + 2**q + 1, 'forward sector')
            need(H - f == 2**(p + q), 'reverse sector')
            need(H == (2**p + 1) * (2**q + 1), 'total count')
            need(C == 1, 'unique cycle')
            controls.append(dict(p=p, q=q, order=len(rows), H=H, f=f, r=H-f, C=C))
    modules = [T for n in range(1, 4) for T in all_tournaments(n)]
    stats = [module_stats(T) for T in modules]
    tests = 0
    for P, (hp, cp) in zip(modules, stats):
        for Q, (hq, cq) in zip(modules, stats):
            rows = build(P, Q)
            H = all_subset_counts(rows)[-1]
            f = all_subset_counts(rows, force_forward=True)[-1]
            need(f == cp*hq + hp*cq + hp*hq, 'general forward formula')
            need(H-f == cp*cq, 'general reverse formula')
            need(H == (cp+hp)*(cq+hq), 'general total formula')
            tests += 1
    proof = Path(__file__).with_name('EDGE_ORDER.md')
    return dict(status='PASS', scope='Literal constructions, not eventual coverage',
                proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),
                transitive_controls=controls, arbitrary_module_pairs=tests,
                arbitrary_modules=len(modules))


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
