#!/usr/bin/env python3
"""Exhaust finite simple labelled unicyclic graphs through order 6.
Compares leaf-peeling/deletion theorem against independent permutation search.
No assert; reject malformed and non-unicyclic graphs; check hostile controls.
"""
from collections import Counter
from itertools import combinations, permutations
from math import isqrt


def require(cond, message):
    if not cond:
        raise ValueError(message)


def validate(n, edges):
    require(type(n) is int and n >= 3, 'n must be an integer >= 3')
    E = set()
    for e in edges:
        require(isinstance(e, (tuple, list)) and len(e) == 2, 'edge must have two ends')
        u, v = e
        require(type(u) is int and type(v) is int and 0 <= u < n and 0 <= v < n and u != v, 'invalid edge ends')
        e = tuple(sorted((u,v)))
        require(e not in E, 'repeated edge')
        E.add(e)
    require(len(E) == n, 'unicyclic graph must have n edges')
    adj = [set() for _ in range(n)]
    for u,v in E:
        adj[u].add(v); adj[v].add(u)
    seen = {0}; stack = [0]
    while stack:
        for w in adj[stack.pop()] - seen:
            seen.add(w); stack.append(w)
    require(len(seen) == n, 'graph disconnected')
    return E, adj


def theorem(n, edges):
    E, adj = validate(n, edges)
    # A unicyclic graph's 2-core is exactly its cycle.
    degree = list(map(len, adj)); alive = set(range(n)); stack = [v for v in alive if degree[v] == 1]
    while stack:
        v = stack.pop()
        if v not in alive:
            continue
        alive.remove(v)
        for w in adj[v] & alive:
            degree[w] -= 1
            if degree[w] == 1: stack.append(w)
    C = {e for e in E if e[0] in alive and e[1] in alive}
    require(len(C) == len(alive) >= 3 and all(degree[v] == 2 for v in alive), 'unexpected cycle core')
    S = {v for v in alive if len(adj[v]) == 3}
    bad = any(len(adj[v]) >= 4 if v in alive else len(adj[v]) >= 3 for v in range(n))
    D = set() if bad else {e for e in C if S <= set(e)}
    closed_count = 0 if bad else (len(C) if len(S) == 0 else 2 if len(S) == 1 else 1 if len(S) == 2 and tuple(sorted(S)) in C else 0)
    require(len(D) == closed_count, 'closed count mismatch')
    return E, C, D, S, bad


def exhaustive_paths(n, E):
    # Generate paths by permutations, no appeal to bridge/cycle criterion.
    paths = set()
    for p in permutations(range(n)):
        if p[0] > p[-1]:
            continue
        path_edges = frozenset(tuple(sorted((p[i],p[i+1]))) for i in range(n-1))
        if path_edges <= E:
            paths.add(path_edges)
    return paths


def components(n, E):
    adj = [set() for _ in range(n)]
    for u,v in E: adj[u].add(v); adj[v].add(u)
    unseen = set(range(n)); count = 0
    while unseen:
        count += 1; stack = [unseen.pop()]
        while stack:
            for w in adj[stack.pop()] & unseen:
                unseen.remove(w); stack.append(w)
    return count


def verify(n, edges):
    E,C,D,S,bad = theorem(n,edges)
    # Independently count components after every edge-deletion set.
    edge_list = sorted(E)
    for r in range(len(E)+1):
        for removed_tuple in combinations(edge_list,r):
            removed = set(removed_tuple)
            expected = 1 + r - bool(removed & C)
            require(components(n, E-removed) == expected, 'deletion component count mismatch')
            if r == 2:
                require((expected == 3) == (removed <= E-C), 'three-way cut classification failed')
    paths = exhaustive_paths(n,E)
    predicted = {frozenset(E-{e}) for e in D}
    require(paths == predicted, 'path edge-set mismatch')
    if paths:
        common = set.intersection(*(set(p) for p in paths))
        union = set.union(*(set(p) for p in paths))
        require(common == E-D, 'forced edge characterization failed')
        require(E-union == (D if len(D)==1 else set()), 'forbidden edge characterization failed')
    return len(paths),len(C),len(S),bad


def main():
    for n in range(3,7):
        possible = list(combinations(range(n),2)); counts = Counter(); total = 0
        for edges in combinations(possible,n):
            try:
                validate(n,edges)
            except ValueError as exc:
                if str(exc) != 'graph disconnected': raise
                continue
            h,c,s,bad = verify(n,edges)
            counts[h] += 1; total += 1
        print(f'n={n} connected unicyclic labelled={total}; path counts={dict(sorted(counts.items()))}')
    # The rejects must stay rejects under python -O.
    invalid = [(2, [(0,1)]), (3,[(0,1),(1,2)]), (3,[(0,0),(0,1),(1,2)]),
               (3,[(0,1),(1,0),(1,2)]), (3,[(0,1),(1,2),(0,3)]),
               (6,[(0,1),(1,2),(0,2),(3,4),(4,5),(3,5)])]
    for n,edges in invalid:
        try: theorem(n,edges)
        except ValueError: pass
        else: raise RuntimeError('malformed input accepted')
    # Triangle + branch at an off-cycle vertex, and triangle + two leaves at one cycle vertex.
    controls = [(6,[(0,1),(1,2),(0,2),(0,3),(3,4),(3,5)]),
                (5,[(0,1),(1,2),(0,2),(0,3),(0,4)])]
    for n,edges in controls:
        h,_,_,bad = verify(n,edges)
        require(bad and h==0, 'branch control unexpectedly passed')
    print('all edge-deletion subsets checked; two-edge 3-component cuts are exactly bridge pairs')
    print(f'hostile controls: {len(invalid)} malformed/disconnected rejected; {len(controls)} branch obstructions verified')
    for n in (15,16,17):
        edges = [(u-1,v-1) for u,v in combinations(range(1,n+1),2) if isqrt(u+v)**2 == u+v]
        E,C,D,S,bad=theorem(n,edges)
        require(not bad and len(D)==1, 'square-sum graph not uniquely traceable')
        require((3,4) in E-C, '4-5 not a bridge')
        require((0,2) in D, '1-3 not the sole omitted cycle edge')
        require((1,3) not in E, '2-4 unexpectedly present')
        print(f'G{n}: cycle={sorted((u+1,v+1) for u,v in C)}, deletion={sorted((u+1,v+1) for u,v in D)}; 4-5 forced, 1-3 forbidden, 2-4 absent')

if __name__ == '__main__': main()
