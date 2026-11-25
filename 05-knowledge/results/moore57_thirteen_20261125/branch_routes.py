"""Route Q25 join Q32 through arbitrary 56-slot perfect matchings at the 57 branches."""
from itertools import combinations
from math import isqrt
from random import Random


def square_edges(n, shift=0):
    return [(i+shift, j+shift) for i, j in combinations(range(1,n+1),2)
            if isqrt(i+j)**2 == i+j]


LEFT = list(range(1,26))
RIGHT = list(range(26,58))
IN_LEFT = square_edges(25)
IN_RIGHT = square_edges(32, 25)
CROSS = [(u,v) for u in LEFT for v in RIGHT]
EDGES = IN_LEFT + CROSS + IN_RIGHT
Q = 56


def matching_row(lists):
    """Distinct representatives for a family of allowed slots on one row."""
    owner = {}

    def augment(v, seen):
        for c in lists[v]:
            if c in seen:
                continue
            seen.add(c)
            if c not in owner or augment(owner[c], seen):
                owner[c] = v
                return True
        return False

    for v in lists:
        assert augment(v, set())
    return {v:c for c,v in owner.items()}


def route(match):
    used = {u:set() for u in LEFT+RIGHT}
    routes = {}

    def put(u,v,c):
        s = match[u,v][c]
        assert c not in used[u] and s not in used[v]
        used[u].add(c)
        used[v].add(s)
        routes[min(u,v),max(u,v)] = (c,s) if u<v else (s,c)

    def available(u,v):
        return [c for c,s in enumerate(match[u,v]) if c not in used[u] and s not in used[v]]

    # At each of the last two left rows, prearrange one redundant exclusion.
    for last, mate, early, right in ((8,1,2,26),(18,7,3,27)):
        put(last,mate,0)
        same = match[last,right][0]
        put(early,right,match[right,early][same])

    for u,v in IN_LEFT:
        if (u,v) not in routes:
            put(u,v,available(u,v)[0])
    deg = {u:sum(u in e for e in IN_LEFT) for u in LEFT}
    order = [2,3] + sorted(set(LEFT)-{2,3,8,18},key=lambda u:(-deg[u],u)) + [8,18]
    for index,u in enumerate(order,1):
        lists = {v:available(u,v) for v in RIGHT if (u,v) not in routes}
        if index <= 23:
            assert min(map(len,lists.values())) >= len(lists)
        else:
            assert min(map(len,lists.values())) >= 31 and max(map(len,lists.values())) >= 32
        for v,c in matching_row(lists).items():
            put(u,v,c)

    deg = {v:sum(v in e for e in IN_RIGHT) for v in RIGHT}
    for u,v in sorted(IN_RIGHT,key=lambda e: -sum(deg[i] for i in e)):
        put(u,v,available(u,v)[0])
    assert set(routes)==set(EDGES) and sum(map(len,used.values()))==2*878
    return routes


def synthetic_matchings(seed):
    rng = Random(seed)
    match = {}
    for u,v in EDGES:
        p = list(range(Q))
        if seed is not None:
            rng.shuffle(p)
        inv = [0]*Q
        for i,j in enumerate(p):
            inv[j] = i
        match[u,v], match[v,u] = p, inv
    return match


if __name__ == "__main__":
    assert list(map(len,(IN_LEFT,CROSS,IN_RIGHT)))==[32,800,46]
    for seed in (None,13,25,32,57):
        assert len(route(synthetic_matchings(seed)))==878
    print("All 878 routes verified for identical matchings and four shuffled fixtures.")
