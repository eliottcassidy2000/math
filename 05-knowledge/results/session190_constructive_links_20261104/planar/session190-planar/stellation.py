#!/usr/bin/env python3
"""Eleven-vertex face-stellation hostile, direct vertex-order counts."""
import json
from verify import cycles,hp_witness
base_edges={(0,1),(1,2),(0,2)}|{(u,v) for u in (0,1,2) for v in (3,4)}
faces=[(3,0,1),(3,1,2),(3,2,0),(4,0,1),(4,1,2),(4,2,0)]
E=base_edges|{tuple(sorted((5+i,v))) for i,f in enumerate(faces) for v in f};V=set(range(11))
cycle=[0,5,3,6,1,8,4,9,2,7]
if len(set(cycle))!=10 or set(cycle)!=V-{10} or not all(tuple(sorted((a,b))) in E for a,b in zip(cycle,cycle[1:]+cycle[:1])):raise RuntimeError('cycle witness')
G={v:{w for e in E if v in e for w in e if w!=v} for v in V}
row={'vertices':sorted(V),'edges':sorted(E),'base_faces':faces,'independent_set':list(range(5,11)),'hamilton_cycle_count':cycles(V,E),'hamilton_path_witness':hp_witness(G),'delete10_cycle_witness':cycle,'deletion_cycle_counts':{v:cycles(V-{v},{e for e in E if v not in e}) for v in V}}
if row['hamilton_cycle_count']!=0 or set(v for v,c in row['deletion_cycle_counts'].items() if c)>set(range(5,11)):raise RuntimeError('nonham certificate')
if not all(row['deletion_cycle_counts'][v]>0 for v in range(5,11)):raise RuntimeError('new vertex deletion')
print(json.dumps(row,indent=2,sort_keys=True))
