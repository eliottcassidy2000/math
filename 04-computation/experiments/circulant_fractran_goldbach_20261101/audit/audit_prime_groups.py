#!/usr/bin/env python3
"""Independent audit of the prime-groups "reading 3" table and the interleaving obstruction
(prime_groups_third_indexing.out; results note section 2.5).

Q_n = square-sum graph on {1..n}, i ~ j iff i + j is a perfect square (i != j).
Unlike the session script (which HARD-CODES the admissible orders 15, 16, 17, 23, 25, ...), this script
DECIDES Hamiltonian-path and Hamiltonian-cycle existence for every n <= NMAX by exhaustive backtracking with
sound pruning, so the admissible list itself is recomputed.

usage (from the repo root):
  python3 04-computation/experiments/circulant_fractran_goldbach_20261101/audit/audit_prime_groups.py [NMAX]
"""
import math
import sys
import time

import networkx as nx

FAILS, MISMATCH = [], []


def report(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        FAILS.append(name)


def claim(name, ok, detail=''):
    print(f'  [{"CLAIM OK" if ok else "CLAIM WRONG"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        MISMATCH.append(name)


def Q(n):
    adj = {v: set() for v in range(1, n + 1)}
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            r = math.isqrt(i + j)
            if r * r == i + j:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def unvisited_connected(adj, used, extra):
    """are the unvisited vertices (plus `extra`, the current end) connected?"""
    rest = [v for v in adj if v not in used]
    if not rest:
        return True
    start = extra if extra is not None else rest[0]
    seen = {start}
    stack = [start]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen and (w not in used):
                seen.add(w)
                stack.append(w)
    return all(v in seen for v in rest)


def hamiltonian(adj, cycle):
    """exact decision by backtracking (returns True/False). Pruning: (i) the unvisited vertices plus the current
    end stay connected; (ii) an unvisited vertex with no unvisited/endpoint neighbour left is fatal; (iii) at most one
    unvisited vertex may be left with a single remaining option (it must be the final vertex) for paths; for cycles the
    final vertex must also be adjacent to the start."""
    n = len(adj)
    if n == 1:
        return not cycle
    degs = {v: len(adj[v]) for v in adj}
    if cycle and min(degs.values()) < 2:
        return False
    if not cycle and sum(1 for v in adj if degs[v] <= 1) > 2:
        return False
    if any(d == 0 for d in degs.values()):
        return False
    order = sorted(adj, key=lambda v: degs[v])
    starts = order[:1] if cycle else order
    leaves = [v for v in adj if degs[v] == 1]
    if not cycle and leaves:
        starts = leaves[:1]      # a leaf must be an endpoint; paths are reversible, so start at one of them
    sys.setrecursionlimit(10000)
    for s in starts:
        used = {s}
        path = [s]

        def rec():
            v = path[-1]
            if len(path) == n:
                return (not cycle) or (s in adj[v])
            # pruning on unvisited vertices
            single = 0
            for u in adj:
                if u in used:
                    continue
                avail = sum(1 for w in adj[u] if w not in used or w == v or (cycle and w == s))
                need = 2
                if avail == 0:
                    return False
                if avail < need:
                    single += 1
            if single > (0 if cycle else 1):      # paths: only the final vertex may have one option left
                return False
            if not unvisited_connected(adj, used, v):
                return False
            cand = sorted((w for w in adj[v] if w not in used),
                          key=lambda w: sum(1 for x in adj[w] if x not in used))
            for w in cand:
                used.add(w)
                path.append(w)
                if rec():
                    return True
                path.pop()
                used.discard(w)
            return False
        if rec():
            return True
    return False


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 45
    t0 = time.time()
    data = {}
    for n in range(1, nmax + 1):
        adj = Q(n)
        G = nx.Graph()
        G.add_nodes_from(adj)
        G.add_edges_from((u, w) for u in adj for w in adj[u] if u < w)
        conn = nx.is_connected(G)
        bic = n >= 3 and nx.is_biconnected(G)
        hp = hamiltonian(adj, cycle=False)
        hc = hamiltonian(adj, cycle=True) if n >= 3 else False
        deg1 = sum(1 for v in adj if len(adj[v]) == 1)
        data[n] = dict(conn=conn, bic=bic, hp=hp, hc=hc, deg1=deg1, mindeg=min(len(adj[v]) for v in adj))
    print(f'  decided Hamiltonian path/cycle for Q_1..Q_{nmax} in {time.time() - t0:.0f}s')
    hp_orders = [n for n in range(1, nmax + 1) if data[n]['hp'] and n > 1]
    hc_orders = [n for n in range(3, nmax + 1) if data[n]['hc']]
    print(f'  Hamiltonian-path orders (n >= 2): {hp_orders}')
    print(f'  Hamiltonian-cycle orders: {hc_orders}')
    adm = hp_orders[:25]
    claim('first 25 Hamiltonian-path-admissible orders = 15,16,17,23,25,...,45 (session hard-coded this list)',
          adm == [15, 16, 17, 23] + list(range(25, 46)), str(adm))
    report('2-connected iff n >= 31 (session132) on 3..NMAX', all(data[n]['bic'] == (n >= 31) for n in range(3, nmax + 1)))
    report('Hamiltonian cycle iff n >= 32 on 3..NMAX; Q_31 2-connected without one',
           hc_orders == list(range(32, nmax + 1)) and data[31]['bic'] and not data[31]['hc'])
    primes = [p for p in range(2, 100) if all(p % d for d in range(2, p))]
    groups = {'A': [2, 3, 5, 7, 11, 13], 'B': [17, 19, 23, 31, 47], 'C': [29, 41, 59, 71], 'D': [37, 43, 53, 61, 71],
              'E': [73, 79, 83, 89, 97]}
    idx = {g: [primes.index(p) + 1 for p in ps] for g, ps in groups.items()}
    claim('group indices B {7,8,9,11,15}, C {10,13,17,20}, D {12,14,16,18,20}',
          idx['B'] == [7, 8, 9, 11, 15] and idx['C'] == [10, 13, 17, 20] and idx['D'] == [12, 14, 16, 18, 20], str(idx))
    rows = []
    for g, ps in groups.items():
        for p, i in zip(ps, idx[g]):
            n = adm[i - 1]
            d = data[n]
            rows.append((g, p, n, d))
            print(f'   {g} {p:2d} -> Q_{n}: connected={d["conn"]}, 2-connected={d["bic"]}, Ham path={d["hp"]}, '
                  f'Ham cycle={d["hc"]}, deg1={d["deg1"]}, mindeg={d["mindeg"]}')
    ok_A = all(d['hp'] and not d['hc'] and not d['bic'] for g, p, n, d in rows if g == 'A')
    ok_DE = all(n >= 32 and d['hc'] and d['bic'] for g, p, n, d in rows if g in 'DE')
    B_hc = {d['hc'] for g, p, n, d in rows if g == 'B'}
    C_hc = {d['hc'] for g, p, n, d in rows if g == 'C'}
    claim('A -> Q_15..Q_26 all with a path, no cycle, not 2-connected; D, E -> orders >= 32 with cycles, 2-connected; '
          'B and C straddle the cycle threshold', ok_A and ok_DE and B_hc == {True, False} and C_hc == {True, False})
    # interleaving: any property monotone in n is a threshold in the index under an increasing indexing
    def interleave(X, Y):
        return min(X) < max(Y) and min(Y) < max(X)
    report('B, C, D index sets pairwise interleave (so no threshold-in-index property is constant on two of them '
           'with different values)', interleave(idx['B'], idx['C']) and interleave(idx['B'], idx['D'])
           and interleave(idx['C'], idx['D']))
    # A, E against the rest: thresholds DO separate them (scope of the obstruction)
    print(f'     note: A (indices 1..6) and E (21..25) are separated from B, C, D by thresholds; e.g. Ham cycle '
          f'(threshold 32) is False on all of A and True on all of D, E')
    print()
    print('SESSION CLAIMS FOUND WRONG: ' + ('; '.join(MISMATCH) if MISMATCH else 'none'))
    print('AUDIT CHECK FAILURES: ' + ('; '.join(FAILS) if FAILS else 'none'))
    return 1 if FAILS else 0


if __name__ == '__main__':
    sys.exit(main())
