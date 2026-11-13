#!/usr/bin/env python3
"""Exact signed-support / full-black-core boundary from frozen native matrices."""
from collections import deque
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
GATES=0
def need(ok,msg):
 global GATES
 GATES+=1
 if not ok:raise RuntimeError(str(msg))
def components(vertices,edges):
 unseen=set(vertices);out=[]
 while unseen:
  seed=min(unseen);unseen.remove(seed);seen={seed};queue=deque([seed])
  while queue:
   u=queue.popleft()
   for v in list(unseen):
    if edges.get((u,v),0):unseen.remove(v);seen.add(v);queue.append(v)
  out.append(sorted(seen))
 return out
def connected(vertices,edges):return bool(vertices) and len(components(vertices,edges))==1

def main():
 raw=(HERE/'data.json').read_bytes();data=json.loads(raw);rows=[]
 expected={3:[1,1],4:[1,2],5:[1,1,1,7],6:[1,1,32],7:[1,1,1,1,268]}
 for row in data['rows']:
  n=row['n'];nodes=row['merged_nodes'];codes=[x['code'] for x in nodes];SC={x['code'] for x in nodes if x['self_converse']};ch=[c for c in codes if c not in SC]
  A={(u,v):w for u,v,w in row['merged_A']};B={(u,v):w for u,v,w in row['merged_blue']};K={(u,v):w-B.get((u,v),0) for (u,v),w in A.items()}
  need(all(w>=0 for w in K.values()),'nonnegative black matrix')
  degree={c:sum(K.get((c,v),0) for v in codes) for c in codes};core=[c for c in codes if degree[c]>0]
  full=components(codes,K);sizes=sorted(map(len,full));need(sizes==expected[n],('full black components',n))
  isolated=[cell[0] for cell in full if len(cell)==1]
  need(all(degree[c]==0 for c in isolated),'all observed isolated vertices have zero black degree')
  need(set(isolated)==set(codes)-set(core),'all zero-black-degree vertices isolated')
  C={(u,v):w for u,v,w in row['ordinary_A']};R={x['code']:x['converse'] for x in row['ordinary_nodes']}
  D={(u,v):C.get((u,v),0)-C.get((u,R[v]),0) for u in ch for v in ch}
  signed={(u,v):1 for (u,v),w in D.items() if u!=v and w}
  for u in ch:
   need(u in core,'every chiral node belongs to positive black core')
   for v in ch:
    need(B.get((u,v),0)==0,'chiral block is all black')
    need(K.get((u,v),0)>=2*abs(D[u,v]),'quantitative signed-support embedding')
  sigma=components(ch,signed)
  for block in sigma:need(connected(block,K),'signed component has connected black fibre')
  blocks=sigma+[[c] for c in core if c in SC];owner={u:i for i,g in enumerate(blocks) for u in g}
  need(set(owner)==set(core),'quotient partitions full black core')
  Q={}
  for u in core:
   for v in core:
    if K.get((u,v),0):Q[owner[u],owner[v]]=Q.get((owner[u],owner[v]),0)+K[u,v]
  qverts=list(range(len(blocks)));coreparts=components(core,K);qparts=components(qverts,Q)
  need(len(coreparts)==len(qparts),'contraction preserves component count')
  need(connected(core,K)==connected(qverts,Q),'exact full-core connectivity equivalence')
  for i,g in enumerate(blocks):need(sum(Q.get((i,j),0) for j in qverts)==sum(degree[u] for u in g),'quotient black dart mass')
  scpositive=[c for c in core if c in SC];scparts=components(scpositive,K)
  attachments=[{'self_converse_component':g,'signed_components_touched':sorted({owner[v] for u in g for v in ch if K.get((u,v),0)>0})} for g in scparts]
  if len(sigma)==1:
   need(connected(core,K)==all(x['signed_components_touched'] for x in attachments),'connected Sigma plus SC attachment criterion')
  witness=None
  if n==5:
   need(len(sigma)==2 and all(len(g)==1 for g in sigma) and connected(core,K),'black-connected does not imply signed-connected')
   witness=[2,8,4]
   need(ch==[2,4] and all(K.get((u,v),0)>0 for u,v in zip(witness,witness[1:])),'explicit black path between signed components')
  rows.append({'n':n,'black_component_sizes':sizes,'zero_black_degree_codes':isolated,
               'black_core_size':len(core),'black_core_connected':connected(core,K),
               'signed_components':sigma,'signed_support_connected':connected(ch,signed),
               'quotient_blocks':blocks,'quotient_black_matrix':[[i,j,w] for (i,j),w in sorted(Q.items())],
               'quotient_component_sizes':sorted(map(len,qparts)),
               'self_converse_attachment_components':attachments,'n5_black_path_between_signed_components':witness})
 # Matrix-independent logical hostile: connected signed subgraph alone does
 # not attach a disjoint positive-degree singleton sector. No tournament claim.
 V=[0,1,2,3];K={(0,1):2,(1,0):2,(2,3):2,(3,2):2};signed={(0,1):1,(1,0):1}
 need(connected([0,1],signed) and not connected(V,K),'abstract missing SC attachment hostile')
 return {'status':'FINITE-EXACT; new note RESERVED pending root review','empty_graph_connected':False,
         'rows':rows,'always_active_checks':GATES,'data_sha256':hashlib.sha256(raw).hexdigest(),
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'synthetic_hostile_scope':'abstract black graph with two chiral vertices and two singleton-sector vertices; not claimed tournament-realizable'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
