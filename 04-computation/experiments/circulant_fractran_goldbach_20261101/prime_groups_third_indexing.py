#!/usr/bin/env python3
"""mac-mini-2026-11-01-fractranTM: the owner's five prime groups against square-sum graphs, the third indexing left
open by session132 (prime_graphs/PROOFS.md: readings 1 and 2 done), and the interleaving obstruction.

Owner's groups (as written; 71 appears twice and 67 is missing, preserved as in session131/132):
  A {2,3,5,7,11,13}  B {17,19,23,31,47}  C {29,41,59,71}  D {37,43,53,61,71}  E {73,79,83,89,97}
Monster-level reading (CITED, session131 moonshine/README.md): A = primes of |Fi22|, A u B = primes of the Baby
Monster, A u B u C = Ogg's supersingular primes; D (with 67 for 71) and E are the non-supersingular primes < 100.

Reading 3: map the k-th prime p_k to the k-th Hamiltonian-path-admissible order a_k of the square-sum graph Q_n
(a = 15, 16, 17, 23, 25, 26, ...; OEIS A090461, Gerbicz's tail CITED), and test whether any of the computed
properties of Q_{a_k} (connected, 2-connected, Hamiltonian path, Hamiltonian cycle, number of degree-1 vertices,
minimum degree) is constant on each group and separates the groups.

Interleaving obstruction (PROVED, elementary): the index sets of B, C, D in the prime order are
B = {7,8,9,11,15}, C = {10,13,17,20}, D = {12,14,16,18,20 (or 19 with 67)}; they interleave, so no property that
depends monotonically on the index (a threshold "k >= k0") can separate B, C, D under readings 1-3, which are all
increasing maps; any separating invariant must be non-monotone in the size, i.e. arithmetic (as the Monster levels
are: divisibility of group orders).

PRE-REGISTERED (before the first run): (R3a) no tested property is constant on every group and different between
some two groups B, C, D; (R3b) "Hamiltonian cycle" is false on all of A (orders 15..26 < 32).
usage: python3 prime_groups_third_indexing.py
"""
import sys
import time

GROUPS = {'A': [2, 3, 5, 7, 11, 13], 'B': [17, 19, 23, 31, 47], 'C': [29, 41, 59, 71], 'D': [37, 43, 53, 61, 71],
          'E': [73, 79, 83, 89, 97]}
FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def graph(n):
    sq = {k * k for k in range(2, 2 * n)}
    return {v: [w for w in range(1, n + 1) if w != v and v + w in sq] for v in range(1, n + 1)}


def connected(G, removed=()):
    vs = [v for v in G if v not in removed]
    if not vs:
        return True
    seen = {vs[0]}
    stack = [vs[0]]
    while stack:
        v = stack.pop()
        for w in G[v]:
            if w not in removed and w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(vs)


def ham(G, cycle, limit=20.0):
    n = len(G)
    t0 = time.time()
    starts = sorted(G, key=lambda v: len(G[v]))
    for s in (starts if not cycle else starts[:1]):
        path = [s]
        used = {s}

        def rec():
            if time.time() - t0 > limit:
                raise TimeoutError
            if len(path) == n:
                return (not cycle) or (s in G[path[-1]])
            v = path[-1]
            for w in sorted((w for w in G[v] if w not in used), key=lambda w: sum(1 for x in G[w] if x not in used)):
                used.add(w)
                path.append(w)
                if rec():
                    return True
                path.pop()
                used.discard(w)
            return False
        try:
            if rec():
                return True
        except TimeoutError:
            return None
    return False


def main():
    P = primes_upto(100)
    assert len(P) == 25
    adm = [n for n in range(2, 60) if n in (15, 16, 17, 23) or n >= 25][:25]
    print(f'first 25 admissible orders a_k: {adm}')
    idx = {g: [P.index(p) + 1 for p in ps] for g, ps in GROUPS.items()}
    print(f'group index sets in the prime order: {idx}')
    props = {}
    for k, n in enumerate(adm, 1):
        G = graph(n)
        deg = [len(G[v]) for v in G]
        two_conn = connected(G) and all(connected(G, removed=(v,)) for v in G)
        hp = ham(G, cycle=False)
        hc = ham(G, cycle=True)
        props[k] = {'order': n, 'connected': connected(G), '2-connected': two_conn, 'Ham path': hp,
                    'Ham cycle': hc, 'deg1': deg.count(1), 'mindeg': min(deg)}
    names = ['connected', '2-connected', 'Ham path', 'Ham cycle', 'deg1', 'mindeg']
    print('reading 3 table (group: prime -> order: properties)')
    for g, ps in GROUPS.items():
        for p in ps:
            k = P.index(p) + 1
            pr = props[k]
            print(f'  {g} {p:>2} -> Q_{pr["order"]:<2}: ' + ', '.join(f'{nm}={pr[nm]}' for nm in names))
    sep = []
    for nm in names:
        vals = {g: {props[P.index(p) + 1][nm] for p in ps} for g, ps in GROUPS.items()}
        const = all(len(v) == 1 for v in vals.values())
        if const:
            distinct = len({next(iter(v)) for g, v in vals.items() if g in 'BCD'}) > 1
            if distinct:
                sep.append(nm)
        print(f'  {nm}: values per group {vals}')
    check('R3a: no tested property is constant on every group and separates two of B, C, D', not sep, f'{sep}')
    check('R3b: Hamiltonian cycle false on all of A', all(props[P.index(p) + 1]['Ham cycle'] is False
                                                           for p in GROUPS['A']))
    inter = (min(idx['C']) < max(idx['B']) and min(idx['D']) < max(idx['B']) and min(idx['D']) < max(idx['C']))
    check('interleaving: B, C, D index ranges overlap pairwise (so no threshold classifier separates them)', inter)
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
