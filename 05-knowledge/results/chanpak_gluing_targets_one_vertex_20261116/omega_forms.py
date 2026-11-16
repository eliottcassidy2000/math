"""All graphs G (up to isomorphism, no restriction to conflict graphs) with I(G,2) in {7, 21},
and every way to write G as a join G1 v G2 (split of its co-components), recording the pair
(I(G1,2), I(G2,2)).  I(G,2) = sum_k alpha_k 2^k, alpha_k = number of independent k-sets.
Method: I(G,2) = 1 + 2|V| + 4 alpha_2 + ... <= 21 forces |V| <= 10 and alpha_2 <= 4, so the
complement Gc has <= 4 edges; enumerate all labelled complements on |V| <= 10 vertices."""
import itertools, networkx as nx
def I2(G):
    # independence polynomial at 2 by brute force over independent sets (|V| <= 10)
    nodes = list(G.nodes()); tot = 0
    for r in range(len(nodes) + 1):
        for S in itertools.combinations(nodes, r):
            if all(not G.has_edge(a, b) for a, b in itertools.combinations(S, 2)): tot += 2 ** r
    return tot
found = {7: [], 21: []}
for v in range(1, 11):
    pairs = list(itertools.combinations(range(v), 2))
    for m in range(0, 5):
        for E in itertools.combinations(pairs, m):
            Gc = nx.Graph(); Gc.add_nodes_from(range(v)); Gc.add_edges_from(E)
            G = nx.complement(Gc)
            val = I2(G)
            if val in found and not any(nx.is_isomorphic(G, H) for H in found[val]):
                found[val].append(G)
def joins(G):
    Gc = nx.complement(G); comps = [set(c) for c in nx.connected_components(Gc)]
    out = set()
    for r in range(1, len(comps)):
        for S in itertools.combinations(range(len(comps)), r):
            A = set().union(*(comps[i] for i in S)); B = set(G.nodes()) - A
            a, b = I2(G.subgraph(A)), I2(G.subgraph(B)); out.add(tuple(sorted((a, b))))
    return sorted(out)
for val in (7, 21):
    print(f"I(G,2) = {val}: {len(found[val])} graphs")
    for G in found[val]:
        Gc = nx.complement(G)
        alpha = [sum(1 for S in itertools.combinations(G.nodes(), k) if all(not G.has_edge(a, b) for a, b in itertools.combinations(S, 2))) for k in range(1, 4)]
        co = nx.number_connected_components(Gc); cc = nx.number_connected_components(G)
        kind = 'join' if co > 1 else ('disjoint union' if cc > 1 else 'prime (neither)')
        print(f"  |V|={G.number_of_nodes()} alpha={alpha} complement edges={sorted(Gc.edges())} -> {kind}; join splits {joins(G)}")
