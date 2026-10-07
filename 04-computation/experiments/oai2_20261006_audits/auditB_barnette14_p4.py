import itertools
import networkx as nx
Q = [(a, b) for a in range(8) for b in range(8) if a < b and bin(a ^ b).count('1') == 1]
G = nx.Graph()
G.add_edges_from((('o', a), ('o', b)) for a, b in Q if a != 0 and b != 0)
G.add_edges_from((('i', a), ('i', b)) for a, b in Q if a != 0 and b != 0)
for x, y in zip([1, 2, 4], [1, 2, 4]):
    G.add_edge(('o', x), ('i', y))
H = nx.convert_node_labels_to_integers(G)
n = H.number_of_nodes(); adj = [list(H[v]) for v in range(n)]
# all Hamiltonian cycles (as edge sets)
cycles = []
def rec(path, seen):
    v = path[-1]
    if len(path) == n:
        if path[0] in adj[v] and path[1] < path[-1]:
            cycles.append(frozenset(frozenset((path[i], path[(i+1) % n])) for i in range(n)))
        return
    for w in adj[v]:
        if not seen[w]:
            seen[w] = True; path.append(w); rec(path, seen); path.pop(); seen[w] = False
seen = [False]*n; seen[0] = True; rec([0], seen)
print("Hamiltonian cycles:", len(cycles))
bad = 0; tot = 0
for v in range(n):
    for a in adj[v]:
        for b in adj[a]:
            if b == v: continue
            for c in adj[b]:
                if c in (v, a): continue
                if v > c: continue   # undirected path counted once
                tot += 1
                P = {frozenset((v, a)), frozenset((a, b)), frozenset((b, c))}
                if not any(P <= C for C in cycles): bad += 1
print("3-edge paths:", tot, " in no Hamiltonian cycle:", bad)
