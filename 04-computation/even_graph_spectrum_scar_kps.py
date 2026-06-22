"""
The even-graph spectrum analogue of the H=7 scar (kind-pasteur S31h).

Owner ask: does the EVEN-graph spectrum have a unique 'forbidden' class mirroring K_3
(the tournament's unique permanent gap 7 = I(K_3,2), K_3 forbidden as a CONFLICT graph)?

KEY OBSERVATION to test: K_3 (the triangle) is FORBIDDEN as a conflict graph Omega(T)
(=> H=7 impossible) BUT K_3 IS itself a valid EVEN graph (all degrees 2).  More generally
K_r is an even graph iff r is ODD.  So the even-graph independence spectrum should REALIZE
the clique values 2r+1 (incl. 7 = I(K_3,2)) that the tournament spectrum forbids.

=> Hypothesis: the even-graph I(G,2)-spectrum is SMOOTH where the tournament H-spectrum has
its 7-scar.  The odd-clique obstruction is orientation-born; forgetting orientation heals it.
The even graph's own imperfection is the odd HOLE (C_5, E_7), the SPGT dual of the odd clique.
"""
import networkx as nx
from itertools import combinations

def I_at_2(G):
    nodes = list(G.nodes()); adj = {v:set(G.neighbors(v)) for v in nodes}
    tot = 0
    for r in range(len(nodes)+1):
        for S in combinations(nodes, r):
            if all(b not in adj[a] for a,b in combinations(S,2)):
                tot += 2**r
    return tot

atlas = nx.graph_atlas_g()
print("EVEN-graph I(G,2) spectrum (connected, all-even-degree graphs <=7 vertices):")
even_spec = set()
clique_even = {}
for G in atlas:
    n = G.number_of_nodes()
    if n < 1: continue
    if not nx.is_connected(G): continue
    if any(d % 2 for _, d in G.degree()): continue   # even graph: all degrees even
    val = I_at_2(G)
    even_spec.add(val)
    # is it a clique K_r? (r odd => even graph)
    if G.number_of_edges() == n*(n-1)//2:
        clique_even[n] = val
print(f"  even-graph spectrum (<=45): {sorted(v for v in even_spec if v<=45)}")
print(f"  odd-clique even graphs K_r -> I(K_r,2)=2r+1: {clique_even}")
print(f"  is 7 in the even-graph spectrum? {7 in even_spec}  (=> via K_3, the triangle)")

# tournament strong-atom H-spectrum (validated, m<=7) for contrast
atoms = sorted(set([3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45]+
    [25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,101,103,105]))
LIM=45
reach=[False]*(LIM+1); reach[1]=True; ch=True
while ch:
    ch=False
    for h in range(1,LIM+1):
        if reach[h]:
            for a in atoms:
                if h*a<=LIM and not reach[h*a]: reach[h*a]=True; ch=True
H_spec = [h for h in range(3,LIM+1,2) if reach[h]]
print(f"\nTOURNAMENT H-spectrum (achievable, <=45): {H_spec}")
print(f"  is 7 in the tournament spectrum? {reach[7]}  (NO -- K_3 forbidden as conflict graph)")

print("\n=== CONTRAST ===")
print(f"  7 in EVEN-graph spectrum: {7 in even_spec};   7 in TOURNAMENT spectrum: {reach[7]}")
print(f"  21 in EVEN-graph spectrum: {21 in even_spec};  21 in TOURNAMENT spectrum: {reach[21]}")
gaps_even = [v for v in range(3,46,2) if v not in even_spec]
print(f"  ODD gaps in the even-graph spectrum <=45: {gaps_even}")
