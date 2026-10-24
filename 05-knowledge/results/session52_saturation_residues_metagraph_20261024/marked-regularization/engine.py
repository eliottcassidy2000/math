#!/usr/bin/env python3
"""Dependency-free exact tournament operations for the Type-I marked lift.

Vertices of D(T): U_i=i, V_i=n+i, z=2*n. A[i][j]=1 means i -> j.
No assertions are used: checks remain active under python -O.
"""
from collections import Counter
from itertools import combinations, permutations


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def tournament(n, mask):
    a = [[0] * n for _ in range(n)]
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = (mask >> bit) & 1
        a[j][i] = 1 - a[i][j]
    return tuple(map(tuple, a))


def validate(a):
    n = len(a)
    need(all(len(row) == n for row in a), 'nonsquare tournament')
    need(all(a[i][i] == 0 for i in range(n)), 'nonzero diagonal')
    need(all(a[i][j] in (0, 1) and a[i][j] + a[j][i] == 1
             for i in range(n) for j in range(i + 1, n)), 'not a tournament')
    return True


def converse(a):
    return tuple(zip(*a))


def induced(a, vertices):
    return tuple(tuple(a[i][j] for j in vertices) for i in vertices)


def regularize(a):
    n = len(a)
    r = [[0] * (2*n+1) for _ in range(2*n+1)]
    for i in range(n):
        r[2*n][i] = r[n+i][2*n] = 1
        for j in range(n):
            r[i][j] = r[n+i][n+j] = a[i][j]
            r[i][n+j] = int(i == j) + a[j][i]
            r[n+j][i] = 1 - r[i][n+j]
    return tuple(map(tuple, r))


def regularize_arcs(a):
    """Independent edge-by-edge transcription of the published recipe."""
    n = len(a)
    arcs = {(2*n, i) for i in range(n)} | {(n+i, 2*n) for i in range(n)}
    arcs |= {(i, n+i) for i in range(n)}
    for i in range(n):
        for j in range(n):
            if a[i][j]:
                arcs.update(((i,j), (n+i,n+j), (j,n+i), (n+j,i)))
    return tuple(tuple(int((i,j) in arcs) for j in range(2*n+1))
                 for i in range(2*n+1))


def scores(a):
    return tuple(map(sum, a))


def unordered_edges(n):
    return {frozenset(p) for p in combinations(range(n), 2)}


def reverse_edges(a, edges):
    r = [list(row) for row in a]
    for edge in edges:
        i, j = tuple(edge)
        r[i][j], r[j][i] = r[j][i], r[i][j]
    return tuple(map(tuple, r))


def changed_edges(a, b):
    return {frozenset((i,j)) for i in range(len(a)) for j in range(i+1,len(a))
            if a[i][j] != b[i][j]}


def block(n, i, j):
    return {frozenset(e) for e in ((i,j),(n+i,n+j),(i,n+j),(j,n+i))}


def directed_cycle(a, cycle):
    return len(set(cycle)) == len(cycle) and all(
        a[cycle[k]][cycle[(k+1) % len(cycle)]] for k in range(len(cycle)))


def hamilton_paths(a):
    return [p for p in permutations(range(len(a)))
            if all(a[p[k]][p[k+1]] for k in range(len(a)-1))]


def hamilton_path_count_dp(a):
    n = len(a)
    d = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1, 1 << n):
        for end in range(n):
            count = d.get((mask,end), 0)
            if count:
                for nxt in range(n):
                    if not mask >> nxt & 1 and a[end][nxt]:
                        key = (mask | (1 << nxt), nxt)
                        d[key] = d.get(key, 0) + count
    return sum(d.get(((1 << n)-1, i), 0) for i in range(n))


def lifted_cycle(path):
    n = len(path)
    return (2*n,) + tuple(v for i in reversed(path) for v in (i,n+i))


def source_from_cycle(cycle):
    n = (len(cycle)-1)//2
    need(cycle[0] == 2*n, 'cycle has wrong distinguished root')
    need(len(set(cycle)) == 2*n+1, 'cycle repeats vertices')
    order = cycle[1::2]
    need(set(order) == set(range(n)), 'cycle does not alternate from U')
    need(cycle[2::2] == tuple(i+n for i in order), 'cycle violates matching')
    return tuple(reversed(order))


def path_complement(a, path):
    fixed = {frozenset(path[k:k+2]) for k in range(len(path)-1)}
    return reverse_edges(a, unordered_edges(len(a)) - fixed)


def cycle_complement(a, cycle):
    fixed = {frozenset((cycle[k],cycle[(k+1) % len(cycle)]))
             for k in range(len(cycle))}
    return reverse_edges(a, unordered_edges(len(a)) - fixed)


def common_out(a, i, j):
    return sum(x and y for x,y in zip(a[i], a[j]))


def universal_common_out_vertices(a):
    return [i for i in range(len(a)) if all(
        common_out(a,i,j) > 0 for j in range(len(a)) if i != j)]


def is_map(a, b, p):
    n = len(a)
    return len(b) == n and len(p) == n and set(p) == set(range(n)) and all(
        a[i][j] == b[p[i]][p[j]] for i in range(n) for j in range(n))


def vertex_signatures(a):
    n = len(a)
    co = tuple(tuple(common_out(a,i,j) for j in range(n)) for i in range(n))
    ci = tuple(tuple(sum(a[k][i] and a[k][j] for k in range(n))
                     for j in range(n)) for i in range(n))
    sig = [tuple(sorted((a[i][j],co[i][j],ci[i][j])
                        for j in range(n) if j != i)) for i in range(n)]
    return sig, co


def isomorphisms(a, b, fixed=None, stop_after=None):
    """Independent exact backtracking; pruning uses only isomorphism invariants.

    Unlike the marked theorem, this routine has no knowledge of U,V,z or D.
    Every surviving bijection is checked on the full adjacency matrix.
    """
    if len(a) != len(b):
        return []
    n = len(a)
    sa, ca = vertex_signatures(a)
    sb, cb = vertex_signatures(b)
    if Counter(sa) != Counter(sb):
        return []
    p = dict(fixed or {})
    if len(set(p.values())) != len(p) or any(sa[i] != sb[j] for i,j in p.items()):
        return []
    if any(a[i][k] != b[j][l] for i,j in p.items() for k,l in p.items()):
        return []
    answer = []
    def search():
        if stop_after is not None and len(answer) >= stop_after:
            return
        if len(p) == n:
            q = tuple(p[i] for i in range(n))
            need(is_map(a,b,q), 'backtracking admitted a non-isomorphism')
            answer.append(q)
            return
        used = set(p.values())
        best, choices = None, None
        for i in range(n):
            if i in p:
                continue
            candidates = [j for j in range(n) if j not in used and sa[i] == sb[j]
                          and all(a[i][k] == b[j][l] and ca[i][k] == cb[j][l]
                                  for k,l in p.items())]
            if not candidates:
                return
            if choices is None or len(candidates) < len(choices):
                best, choices = i, candidates
        for j in choices:
            p[best] = j
            search()
            del p[best]
    search()
    return answer


def isomorphic(a, b):
    return bool(isomorphisms(a,b,stop_after=1))


def brute_automorphisms(a):
    return [p for p in permutations(range(len(a))) if is_map(a,a,p)]


def diagonal_lift(p):
    n = len(p)
    return tuple(p) + tuple(n+i for i in p) + (2*n,)


def order_join(*blocks):
    n = sum(map(len, blocks))
    a = [[0]*n for _ in range(n)]
    off = 0
    for b in blocks:
        for i in range(len(b)):
            for j in range(len(b)):
                a[off+i][off+j] = b[i][j]
            for j in range(off+len(b),n):
                a[off+i][j] = 1
        off += len(b)
    return tuple(map(tuple,a))


def recenter_source_map(q):
    """D(1+Q) -> D(Q+1), with Q at old indices 1..n-1/new 0..n-2."""
    n = len(q)+1
    p = [None]*(2*n+1)
    p[0], p[n], p[2*n] = 2*n, n-1, 2*n-1
    for i in range(1,n):
        p[i], p[n+i] = i-1, n+i-1
    return tuple(p)


def scc_sizes_in_order(a):
    """Transitive-closure SCCs, then topological order by dominance."""
    n = len(a)
    reach = [[bool(a[i][j] or i==j) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    remaining = set(range(n)); comps=[]
    while remaining:
        i=min(remaining)
        c={j for j in remaining if reach[i][j] and reach[j][i]}
        comps.append(c); remaining-=c
    frozen = tuple(comps)
    comps.sort(key=lambda c: -sum(a[min(c)][min(d)] for d in frozen if c != d))
    return tuple(map(len,comps))


def matrix_lines(a):
    return [''.join(map(str,row)) for row in a]
