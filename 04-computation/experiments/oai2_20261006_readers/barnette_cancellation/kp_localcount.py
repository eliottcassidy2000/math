import itertools, sys
import knight_paths2 as K
n = 6
eng = K.Engine(n)
def count(pts):
    path = [n*p[0] + p[1] for p in pts]
    EP = set(eng.eid[(min(u, v), max(u, v))] for u, v in zip(path, path[1:]))
    blocking = []
    for v in range(eng.nV):
        inc = [e for e in eng.inc[v] if e not in EP]
        for D in itertools.combinations(inc, 6):
            if K.find_hc(eng.E, eng.nV, set(D), list(EP)) is None:
                blocking.append((v, D))
    kinds = {}
    I = set(path[1:-1]); a, b = path[0], path[-1]
    for v, D in blocking:
        keep = [w for w in eng.adj[v] if eng.eid[(min(v, w), max(v, w))] not in D]
        if v in path: kind = "at-path-square"
        elif a in keep and b in keep: kind = "closing"
        elif any(w in I for w in keep) or len([w for w in keep if w not in I]) <= 1: kind = "starve"
        else: kind = "other"
        kinds[kind] = kinds.get(kind, 0) + 1
    print(pts, "local blocking 6-sets:", len(blocking), kinds, flush=True)
for pts in ([(0,0),(1,2),(2,4)], [(0,0),(1,2),(5,3)], [(0,0),(1,2),(2,4),(3,0)], [(0,0),(1,2),(3,3),(1,4)]):
    count(pts)
