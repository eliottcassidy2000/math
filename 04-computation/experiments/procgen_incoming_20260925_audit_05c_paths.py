import networkx as nx
from math import isqrt
def Q(N):
    g = nx.Graph(); g.add_nodes_from(range(1, N+1))
    for x in range(1, N+1):
        for y in range(x+1, N+1):
            s = x+y
            if isqrt(s)**2 == s: g.add_edge(x, y)
    return g
q25 = Q(25)
ok, cert = nx.check_planarity(q25, counterexample=True)
A, B = [3, 5, 12], [4, 11, 13]
branch = set(A+B)
# walk each branch-to-branch path in the witness
paths = []
for a in A:
    for nb in cert.neighbors(a):
        path = [a, nb]
        prev, cur = a, nb
        while cur not in branch:
            nxt = [w for w in cert.neighbors(cur) if w != prev]
            assert len(nxt) == 1
            prev, cur = cur, nxt[0]
            path.append(cur)
        paths.append(path)
for p in sorted(paths):
    print(p, "sums:", [x+y for x, y in zip(p, p[1:])])
# own verification that these 9 paths form a K33 subdivision in Q25
ends = sorted((min(p[0],p[-1]), max(p[0],p[-1])) for p in paths)
print("endpoint pairs:", ends)
inter = [v for p in paths for v in p[1:-1]]
print("interiors disjoint:", len(inter) == len(set(inter)), "avoid branch:", not (set(inter) & branch))
print("all edges square sums:", all(isqrt(x+y)**2 == x+y for p in paths for x, y in zip(p, p[1:])))
