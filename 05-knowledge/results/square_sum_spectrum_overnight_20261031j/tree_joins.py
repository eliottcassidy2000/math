#!/usr/bin/env python3
"""Port-linear tree joins; owner, scale and branch coordinates are retained."""
from dataclasses import dataclass
from itertools import product
from math import isqrt
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031i'))
from elliptic_ports import join,Closed,Atlas,condition_cube,edge,need,digest
from pythagorean_glue import from_cube

@dataclass
class Packet:
    graph:Closed
    ports:dict
    owner_roots:dict
    source_bits:list
    spent:list

def source(owner,G,ports):
    for name,(e,d) in ports.items():
        Closed(G.vertices,G.base,G.moves,e,d).check_port()
    need(len({e for e,d in ports.values()})==len(ports),'one name per physical port')
    return Packet(G,{(owner,name):(e,d) for name,(e,d) in ports.items()},
                  {owner:1},[(owner,len(G.moves))],[])

def attach(A,pa,B,pb):
    need(not (A.owner_roots.keys()&B.owner_roots.keys()),'tree edge joins distinct components')
    need(pa in A.ports and pb in B.ports,'unspent registered ports')
    ea,da=A.ports[pa];eb,db=B.ports[pb]
    need(da==db,'pairwise common square class')
    GA=Closed(A.graph.vertices,A.graph.base,A.graph.moves,ea,da)
    GB=Closed(B.graph.vertices,B.graph.base,B.graph.moves,eb,db)
    G,receipt=join(GA,GB)
    if receipt['method']=='matched norm conic':
        roots=[isqrt(s)*a for s,a in zip(receipt['frame']['source_scales'],receipt['input_alignment_roots'])]
    else:roots=[receipt['left_scale_root'],receipt['right_scale_root']]
    ports={};owners={}
    for side,used,r in [(A,pa,roots[0]),(B,pb,roots[1])]:
        owners.update({o:r*s for o,s in side.owner_roots.items()})
        for key,(e,d) in side.ports.items():
            if key==used:continue
            mapped=edge(*(r*r*v for v in e));ports[key]=(mapped,d)
            Closed(G.vertices,G.base,G.moves,mapped,d).check_port()
    spent=A.spent+B.spent+[(pa,pb)]
    result=Packet(G,ports,owners,A.source_bits+B.source_bits,spent)
    need(len(result.graph.moves)==sum(b for o,b in result.source_bits),'exact additive source clock')
    return result,{'class':da,'method':receipt['method'],'scale_root_bits':[r.bit_length() for r in roots],
                   'output_vertices':len(G.vertices),'output_bits':len(G.moves),
                   'remaining_ports':len(ports),'spent_pairs':len(spent)}

def recover(P,originals,bits):
    es=P.graph.edges(bits);offset=0;owner_sets=[]
    for owner,count in P.source_bits:
        G=originals[owner];s=P.owner_roots[owner]**2;V={s*v for v in G.vertices};owner_sets.append(V)
        inherited={edge(u//s,v//s) for u,v in es if u in V and v in V}
        restored={e for key,(e,d) in G.all_ports.items() if any((owner,key) in pair for pair in P.spent)}
        need(inherited|restored==G.edges(bits[offset:offset+count]),'exact original cycle after restoring spent ports')
        offset+=count
    need(sum(map(len,owner_sets))==len(set().union(*owner_sets))==len(P.graph.vertices),'exact disjoint owners')
    return digest(P.graph.walk(bits))

def components(E,n):
    adj={v:[] for v in range(n)}
    for a,b in E:adj[a].append(b);adj[b].append(a)
    need(all(len(ns)==2 for ns in adj.values()),'cyclic hostile has degree two')
    seen=set();lengths=[]
    for v in range(n):
        if v in seen:continue
        todo=[v];seen.add(v);count=0
        while todo:
            z=todo.pop();count+=1
            for w in adj[z]:
                if w not in seen:seen.add(w);todo.append(w)
        lengths.append(count)
    return sorted(lengths)

def main():
    atlas=Atlas();G,pins=condition_cube(atlas.cycle(2812,'A'),(160,1440),10)
    ports={'one':((9,16),1),'five':((20,80),5),'ten':((160,1440),10)}
    G.all_ports=ports
    B=from_cube(atlas.cycle(2861,'A'),(196,2304));B.all_ports={'one':((196,2304),1)}
    originals={0:G,1:B,2:G,3:G}
    P=[source(i,originals[i],originals[i].all_ports) for i in range(4)]
    C,r1=attach(P[0],(0,'one'),P[1],(1,'one'))
    C,r2=attach(C,(0,'five'),P[2],(2,'five'))
    C,r3=attach(C,(2,'ten'),P[3],(3,'ten'))
    hashes=[]
    for j in range(8):
        bits=''.join(str(hashlib.sha256(f'{j}:{b}'.encode()).digest()[0]%2) for b in range(len(C.graph.moves)))
        hashes.append(recover(C,originals,bits))
    need(len(C.graph.vertices)==11297 and len(C.graph.moves)==97,'mixed tree exact product')
    reused=False;cycle=False
    try:attach(C,(0,'five'),source(4,G,ports),(4,'five'))
    except RuntimeError:reused=True
    try:attach(C,(0,'five'),C,(3,'five'))
    except RuntimeError:cycle=True
    need(reused and cycle,'linear port and tree ownership guards')
    # Euler-brick triangle: every port is forced and each pair shares a vertex.
    V={44**2,117**2,240**2};E={edge(a,b) for a in V for b in V if a<b}
    triangle=Closed(V,E,[],min(E));triangle.check_port()
    tp={str(j):(e,1) for j,e in enumerate(sorted(E))};triangle.all_ports=tp
    center=source(0,triangle,tp);star=[]
    for j in range(3):
        center,r=attach(center,(0,str(j)),source(j+1,triangle,tp),(j+1,str((j+1)%3)))
        star.append(r)
    star_hash=recover(center,{i:triangle for i in range(4)},'')
    # Abstract cyclic incidence control: arithmetic feasibility is not asserted.
    cyclic=[]
    for bits in product((0,1),repeat=3):
        edges={edge(3*j,3*j+2) for j in range(3)}
        for j,flip in enumerate(bits):
            a=[3*j+1,3*j+2];k=(j+1)%3;b=[3*k,3*k+1]
            if flip:b.reverse()
            edges|={edge(x,y) for x,y in zip(a,b)}
        sizes=components(edges,9)
        need((sizes==[9])==(sum(bits)%2==0),'cyclic wiring retains parity holonomy')
        cyclic.append({'pairing_bits':bits,'component_sizes':sizes})
    print(json.dumps({'status':'passed','native_source_condition':pins,'mixed_class_chain':[r1,r2,r3],
          'native_owner_count':4,'native_cycle_sha256':hashes,'shared_vertex_star':star,'star_cycle_sha256':star_hash,
          'reused_port_rejected':reused,'cycle_edge_rejected':cycle,'cyclic_abstract_hostiles':cyclic,
          'scope':'All-size tree product follows by exact two-edge-cut induction. Native mixed classes1,5,10 and a shared-vertex Euler-brick star are arithmetic controls. The cyclic countermodel is abstract incidence only; it makes no extra arithmetic claim.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
