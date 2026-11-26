"""Exact replay for prime-window matchings and their induced Moore control.

    python3 prime_contacts.py
    python3 prime_contacts.py ../moore57_quadratic_tensors_20261125
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import factorial, isqrt, prod
from pathlib import Path
import sys

P = tuple(x for x in range(2, 102) if all(x % d for d in range(2, isqrt(x) + 1)))
CASES = (P[1:-1], P[:-1], P[1:], P)
C6 = (3, 61, 83, 17, 47, 97)
FORK = (5, 31, 59, 11, 89)


def square_graph(labels):
    return {p: {q for q in labels if p != q and isqrt(p + q)**2 == p + q}
            for p in labels}


def components(g):
    todo = set(g)
    while todo:
        seen = {min(todo)}
        pending = list(seen)
        while pending:
            w = pending.pop()
            pending.extend(g[w] - seen)
            seen.update(g[w])
        todo -= seen
        yield seen


def cover_count(g, W):
    W = sorted(W)
    ix = {q: i for i, q in enumerate(W)}
    active = set(W)
    masks = sorted((sum(1 << ix[q] for q in g[p] & active) for p in W),
                   key=int.bit_count)

    @lru_cache(None)
    def count(i, used):
        if i == len(W):
            return 1
        ans, free = 0, masks[i] & ~used
        while free:
            bit = free & -free
            free ^= bit
            ans += count(i+1, used | bit)
        return ans

    return count(0, 0)


def forests(g, W):
    edges = sorted((p, q) for p in W for q in g[p] if p < q)
    for n_drop in range(len(edges)+1):
        ways, witness = 0, None
        for drop in combinations(range(len(edges)), n_drop):
            keep = set(range(len(edges))) - set(drop)
            deg = Counter()
            parts = {p: {p} for p in W}
            for i in sorted(keep):
                p, q = edges[i]
                a, b = parts[p], parts[q]
                deg[p] += 1
                deg[q] += 1
                if a is b or max(deg[p], deg[q]) > 2:
                    break
                a |= b
                for v in b:
                    parts[v] = a
            else:
                ways += 1
                if witness is None:
                    witness = [edges[i] for i in sorted(keep)]
        if ways:
            return len(W) - len(edges) + n_drop, ways, witness
    raise AssertionError


def list_paths(g):
    ans = []
    for W in components(g):
        p = min(p for p in W if len(g[p]) <= 1)
        path = [p]
        while g[p] - set(path):
            p, = g[p] - set(path)
            path.append(p)
        assert set(path) == W
        ans.append(path)
    return ans


def prime_checks():
    pi0 = (79, 97, 31, 29, 89, 23, 19, 17, 13, 7, 5, 67, 59, 101,
           53, 47, 41, 83, 37, 73, 71, 2, 61, 11, 3, 43)
    perm = dict(zip(P, pi0))
    assert len(P) == len(set(pi0)) == 26 and all(perm[perm[p]] == p for p in P)
    assert sorted(p for p in P if perm[p] not in square_graph(P)[p]) == [37, 67]
    for labels, expected in zip(CASES, [(4,7), (3,6), (3,7), (2,6)]):
        g = square_graph(labels)
        isolates = sorted(p for p in g if not g[p])
        active = set(g) - set(isolates)
        count = cover_count(g, active)
        assert count == 1
        defects = len(isolates)
        cover = [forests(g, W) for W in components(g)]
        linear = sum(t[0] for t in cover)
        assert (defects, linear) == expected
        pathg = {p:set() for p in g}
        for _, _, witness in cover:
            for p, q in witness:
                pathg[p].add(q)
                pathg[q].add(p)
        optimal = sum((-1)**j * (factorial(defects)//factorial(j))
                      for j in range(defects+1))
        print("primes", (labels[0],labels[-1]), "length", len(g),
              "edges", sum(map(len,g.values()))//2,
              "isolates", isolates, "square permutation defects", defects,
              "optimal derangements", optimal, "spanning path defects", linear-1,
              "optimal square forests", prod(t[1] for t in cover))
        print("  paths", list_paths(pathg))
    h = square_graph(C6 + FORK)
    pairs = set(map(frozenset, zip(C6, C6[1:]+C6[:1])))
    pairs.update(map(frozenset, [(5,31),(5,59),(5,11),(11,89)]))
    assert {frozenset((p,q)) for p in h for q in h[p]} == pairs
    # Note shortest square-sum 4-cycles, with and without 2, by largest label.
    primes = [x for x in range(2,194) if all(x%d for d in range(2,isqrt(x)+1))]
    g = square_graph(primes)
    first = {}
    for x,y in combinations(primes,2):
        for a,b in combinations(sorted(g[x]&g[y]),2):
            v = sorted((x,y,a,b))
            key = "through 2" if v[0] == 2 else "odd only"
            first[key] = min(first.get(key, v[::-1]), v[::-1])
    assert {k:v[0] for k,v in first.items()} == {"through 2":167,"odd only":193}
    print("first square-sum C4 labels by largest prime", first)


def hoffman_singleton():
    H = [set() for _ in range(50)]

    def add(v,w):
        H[v].add(w)
        H[w].add(v)

    for i in range(5):
        for x in range(5):
            add(5*i+x, 5*i+(x+1)%5)
            add(25+5*i+x, 25+5*i+(x+2)%5)
            for j in range(5):
                add(5*i+x, 25+5*j+(x+i*j)%5)
    assert all(len(V)==7 for V in H)
    assert all(len(H[v]&H[w]) == (w not in H[v]) for v,w in combinations(range(50),2))
    return H


def hexagons(H):
    def grow(path):
        if len(path)==6:
            if path[1]<path[-1] and path[0] in H[path[-1]]:
                yield path
        else:
            for x in sorted(H[path[-1]]):
                if x>path[0] and x not in path:
                    yield from grow((*path,x))

    for v in range(len(H)):
        yield from grow((v,))


def small_moore_check():
    H = hoffman_singleton()
    hist, holonomy = Counter(), Counter()
    loops = [(0,2,4),(1,3,5),(0,2,5,3),(0,4,1,3),(2,4,1,5)]
    for C in hexagons(H):
        U = set(range(50)) - set(C) - set().union(*(H[v] for v in C))
        g = {v:H[v]&U for v in U}
        shapes = tuple(sorted(tuple(sorted(len(g[v]) for v in W))
                              for W in components(g)))
        D = [next(iter(H[a]&H[b])) for a,b in zip(C[:3],C[3:])]
        apex, = set.intersection(*(H[d] for d in D))
        assert apex in U
        hist[shapes] += 1
        V = [set(w for w in H[v] - set(C) if len(H[w]&set(C))==1) for v in C]
        assert list(map(len,V)) == [4]*6
        types = []
        for loop in loops:
            p = {}
            for x in V[loop[0]]:
                y = x
                for target in loop[1:]+loop[:1]:
                    y, = H[y]&V[target]
                p[x] = y
            assert set(p) == set(p.values())
            todo = set(p)
            lengths = []
            while todo:
                y = root = min(todo)
                length = 0
                while y in todo:
                    todo.remove(y)
                    length += 1
                    y = p[y]
                assert y == root
                lengths.append(length)
            types.append(tuple(sorted(lengths)))
        holonomy[tuple(types)] += 1
    expected = {((1,1,1,1,4), (2,)*6, (2,)*6):5250}
    assert hist == expected
    assert holonomy == {((4,),(4,),(2,2),(2,2),(2,2)):5250}
    print("Hoffman–Singleton prism loop cycle types",dict(holonomy))
    print("Hoffman–Singleton hexagon antineighborhoods", dict(hist))
    Q = [0,5,7,4,6,9,2,3,30,29,8,15,13,11,1,26,46,25,16,14,10,35,27,32,31,12]
    assert len(set(Q))==26
    assert all(Q[j] in H[Q[i]] for i,j in combinations(range(26),2)
               if isqrt(P[i]+P[j])**2==P[i]+P[j])
    extra = sum(Q[j] in H[Q[i]] and isqrt(P[i]+P[j])**2!=P[i]+P[j]
                for i,j in combinations(range(26),2))
    assert extra > 0
    print("A non-induced 26-prime embedding in Hoffman–Singleton, added edges", extra)


def chart_check(directory):
    sys.path.insert(0,str(directory.resolve()))
    from planted_stars import orbit_chart

    chart, stars = orbit_chart('period_fourteen_site.txt', 4, ((10,45,2),(27,49,0)))
    M = {(a,x,b):y for (a,b),p in chart.items() for x,y in enumerate(p)}
    for (a,x), ys in stars.items():
        for b,y in ys:
            for k,v in [((a,x,b),y),((b,y,a),x)]:
                assert M.get(k,v) == v
                M[k] = v

    def leaf(a,x):
        return 58+56*a+x

    H = [set() for _ in range(3250)]

    def join(v,w):
        H[v].add(w)
        H[w].add(v)

    for a in range(57):
        join(0,1+a)
        for x in range(56):
            join(1+a,leaf(a,x))
    for (a,x,b),y in M.items():
        join(leaf(a,x),leaf(b,y))
    assert sum(map(len,H))//2==13231
    for v, N in enumerate(H):
        reached = N | {v}
        for w in N:
            S = H[w] - {v}
            assert not reached & S
            reached |= S

    hist, triple_hist = Counter(), Counter()
    pet = set()
    for a in range(57):
        for x in range(56):
            for b,c in combinations((i for i in range(57) if (a,x,i) in M),2):
                u,v = M[a,x,b], M[a,x,c]
                t,s = M.get((c,v,b)), M.get((b,u,c))
                z,w = M.get((b,t,a)), M.get((c,s,a))
                kind = ('unknown' if z is None or w is None
                        else 'concurrent' if z==w else 'separated')
                hist[kind] += 1
                if {a,b,c} == {0,28,56}:
                    triple_hist[kind] += 1
                if kind == 'concurrent':
                    pet.add(frozenset([0,1+b,leaf(b,u),leaf(a,x),leaf(c,v),1+c,
                                       1+a,leaf(b,t),leaf(c,s),leaf(a,z)]))
    assert hist == {'unknown':262598,'concurrent':1770,'separated':27168}
    assert triple_hist=={'concurrent':168} and len(pet)==295
    print("root hexagons", dict(sorted(hist.items())), "on rows 0,28,56",
          dict(triple_hist), "distinct rooted Petersen subgraphs",len(pet))

    # A nonconcurrent hexagon carrying the 11-prime induced pattern in ALL completions.
    Q = (0,2,115,58,1654,29)
    F = (60,1542,2019,3196,1124)
    D = (1,118,1651)
    assert all(H[v]&H[w]=={d} for v,w,d in zip(Q[:3],Q[3:],D))
    assert H[D[0]] & H[D[1]] == {61} and H[D[0]] & H[D[2]] == {111}
    mark = dict(zip(C6+FORK,Q+F))
    assert len(set(mark.values())) == len(mark)
    for a,b in combinations(mark,2):
        v,w = mark[a],mark[b]
        if isqrt(a+b)**2 == a+b:
            assert v in H[w]
        else:
            assert v not in H[w]
            assert (len(H[v])==57 or len(H[w])==57 or H[v]&H[w]
                    or any(H[t]&H[w] for t in H[v]))
    print("forced induced 11-prime pattern", mark)
    print("forced nonconcurrent Petersen test", Q, D, (61,111))


if __name__=='__main__':
    prime_checks()
    small_moore_check()
    if len(sys.argv) > 1:
        chart_check(Path(sys.argv[1]))
