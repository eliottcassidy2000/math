import importlib.util
spec = importlib.util.spec_from_file_location('m','/Users/e/Documents/GitHub/math/04-computation/hampath_cycle_correspondence_mac-mini-2026-06-15-S6.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

# KEY: T_G's edges are ALWAYS oriented consistently with the vertex order (acyclic on G's edges).
# So T_G is strong ONLY via the non-edges (which become backward arcs). 
# Claim to test: exists-order-T_G-strong  <=>  G has a Ham cycle in its COMPLEMENT? or
#   <=> G is NOT complete and the non-edges + acyclic give strong?
# Let's directly correlate exists-order-strong with: "complement has a Ham path/cycle",
# "G not complete", etc.

import networkx as nx
for n in range(3,8):
    reps = m.all_noniso_graphs(n)
    # tabulate exists-strong vs (G complete?), (complement connected?), (complement Ham path?)
    rows={}
    for es in reps:
        adj=m.edges_to_adj(n,es)
        exS=m.exists_order_tournament_strong(n,adj)
        full = n*(n-1)//2
        complete = (len(es)==full)
        # complement
        G=nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(es)
        Gc=nx.complement(G)
        comp_connected = nx.is_connected(Gc) if n>0 else True
        comp_adj=[[Gc.has_edge(i,j) for j in range(n)] for i in range(n)]
        comp_hp=m.has_hamiltonian_path(n,comp_adj)
        key=(int(exS),int(complete),int(comp_connected),int(comp_hp))
        rows[key]=rows.get(key,0)+1
    print(f"n={n}: (exStrong, Gcomplete, compConnected, compHamPath) -> count")
    for k in sorted(rows): print("   ",k,rows[k])
    # test implication: exStrong <=> NOT complete?  and exStrong <=> compConnected?
    fp_complete=sum(v for k,v in rows.items() if k[0]==1 and k[1]==1)
    # exStrong <=> compConnected
    bad_cc=sum(v for k,v in rows.items() if k[0]!=k[2])
    print(f"    exStrong-and-complete count: {fp_complete}  | exStrong != compConnected count: {bad_cc}")
