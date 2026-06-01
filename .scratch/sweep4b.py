from lrc_coloring_core import *
from planarity import is_planar
import collections

# For each set: count nonplanar cells; check correspondence nonplanar <=> omega>=5.
# Also: is nonplanarity EXACTLY the t-near-0 K_n cells for n>=5, or do mid-cells go nonplanar?
print("PLANARITY (real DMP test):")
for speeds in [(1,2,3,4),(1,2,3,5),(2,3,5,7),(1,2,4,8),
               (1,2,3,4,5),(1,3,4,5,9),(2,3,4,5,6),
               (1,2,3,4,5,6),(1,3,5,7,9,11),(1,2,4,8,16,32)]:
    n=len(speeds)+1
    crit=critical_times(speeds,n); mids=cell_midpoints(crit)
    np_cells=[]; mism=0
    for t in mids:
        N,edges,pts=danger_graph(speeds,n,t)
        pl=is_planar(N,edges)
        w=clique_number_arc(speeds,n,t)
        if not pl:
            np_cells.append((t,w))
        # K5 subgraph (omega>=5) forces nonplanar; check equivalence
        if (not pl) != (w>=5):
            mism+=1
    print(f"{speeds} n={n}: nonplanar {len(np_cells)}/{len(mids)} ; nonplanar<=>omega>=5 mismatches: {mism}")
    if np_cells:
        print("     nonplanar cells (t, omega):", [(str(t),w) for t,w in np_cells])
