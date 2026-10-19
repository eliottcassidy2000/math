#!/usr/bin/env python3
from itertools import permutations
from collections import Counter
import json
checks=Counter()
def ck(ok,tag,detail=None):
 checks[tag]+=1
 if not ok:raise RuntimeError((tag,detail))
E=((0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(4,7))
TRANS=[sum(1<<j for j in range(i+1,8))for i in range(8)]
HOST=[210,188,169,81,164,201,22,72]
def histogram(host):
 out=[0]*128
 for p in permutations(range(8)):
  mask=sum((1-((host[p[u]]>>p[v])&1))<<k for k,(u,v)in enumerate(E))
  out[mask]+=1
 return out
for H in(TRANS,HOST):
 for i in range(8):
  ck(not((H[i]>>i)&1),'host_loop')
  for j in range(i+1,8):ck(((H[i]>>j)&1)+((H[j]>>i)&1)==1,'host_orientation')

# Independent injective tree backtracking, prescribed orientation.
def count_embeddings(edges,n,mask,host):
 arcs=[(v,u)if(mask>>k)&1 else(u,v)for k,(u,v)in enumerate(edges)]
 order=sorted(range(n),key=lambda v:-sum(v in e for e in edges))
 assignment=[-1]*n
 def rec(pos,used):
  if pos==n:return 1
  v=order[pos];total=0
  for image in range(len(host)):
   if(used>>image)&1:continue
   good=True
   for a,b in arcs:
    if a==v and assignment[b]>=0 and not((host[image]>>assignment[b])&1):good=False;break
    if b==v and assignment[a]>=0 and not((host[assignment[a]]>>image)&1):good=False;break
   if good:
    assignment[v]=image;total+=rec(pos+1,used|(1<<image));assignment[v]=-1
  return total
 return rec(0,0)

# Independent exact automorphism and oriented-isomorphism census.
edge_set={frozenset(e) for e in E}
automorphisms=[p for p in permutations(range(8)) if {frozenset((p[u],p[v])) for u,v in E}==edge_set]
ck(len(automorphisms)==2,'D_automorphism_group')
def transport(mask,p):
 arcs={(p[v],p[u]) if (mask>>k)&1 else (p[u],p[v]) for k,(u,v) in enumerate(E)}
 return sum(int((v,u) in arcs)<<k for k,(u,v) in enumerate(E))
classes={}
for mask in range(128):
 orbit={transport(mask,p) for p in automorphisms}
 ck(len(orbit)==2,'orientation_orbit_free',mask)
 eps=((mask>>2)^(mask>>4))&1
 ck(all((((m>>2)^(m>>4))&1)==eps for m in orbit),'selector_orbit_invariance',mask)
 classes[min(orbit)]=eps
ck(len(classes)==64 and sum(classes.values())==32,'orientation_class_split')

# Sharpness of the all-order source-free primitive floor at core order seven.
A7=((0,1),(0,2),(0,3),(1,4),(2,5),(4,6))
T7=[sum(1<<j for j in range(i+1,7))for i in range(7)]
H7=[118,44,96,85,102,8,34]
for i in range(7):
 ck(not((H7[i]>>i)&1),'seven_host_loop')
 ck(H7[i].bit_count()!=6,'seven_host_source_free')
 for j in range(i+1,7):ck(((H7[i]>>j)&1)+((H7[j]>>i)&1)==1,'seven_host_orientation')
seven_counts=[]
for host,want in [(T7,66),(H7,81)]:
 value=count_embeddings(A7,7,52,host)
 arcs=[(v,u)if(52>>k)&1 else(u,v)for k,(u,v)in enumerate(A7)]
 literal=sum(all((host[p[u]]>>p[v])&1 for u,v in arcs)for p in permutations(range(7)))
 ck(value==literal==want,'sharp_seven_core_witness')
 seven_counts.append(value)

records=[]
for name,H in [('transitive',TRANS),('inherited_D8_hostile',HOST)]:
 hist=histogram(H)
 defect=(hist[0]+hist[1<<4])%2
 for mask,count in enumerate(hist):
  epsilon=((mask>>2)^(mask>>4))&1
  ck(count%2==(1^(epsilon&defect)),'all_orientation_selector',(name,mask))
  for bit in range(7):
   derivative=(count+hist[mask^(1<<bit)])%2
   ck(derivative==(defect if bit in(2,4)else 0),'all_first_derivatives',(name,mask,bit))
 for mask in(0,5,16,127):
  ck(count_embeddings(E,8,mask,H)==hist[mask],'independent_embedding_count',(name,mask))
 # F is D minus edge(1,5), retaining isolated vertex5.
 F=tuple(e for e in E if e!=(1,5))
 countF=count_embeddings(F,8,0,H)
 ck(countF==hist[0]+hist[16],'retained_isolate_count',name)
 records.append({'host':name,'out_bits':H,'defect_parity':defect,'D_mask5_count':hist[5],
   'D_mask0_count':hist[0],'D_mask16_count':hist[16],'F_reference_count':countF,
   'odd_D_orientations':sum(v%2 for v in hist),'histogram':hist})
ck(records[0]['D_mask5_count']==165,'inherited_transitive_witness')
ck(records[1]['D_mask5_count']==302,'inherited_hostile_witness')
ck(records[0]['defect_parity']==0 and records[1]['defect_parity']==1,'defect_nonconstancy')
print(json.dumps({'status':'PASS','checks':dict(checks),'total_checks':sum(checks.values()),
 'source_edges':E,'mask_convention':'bit1 reverses the listed edge; bit0 keeps it',
 'critical_edge_bits':[2,4],'constant_odd_masks':64,'nonconstant_masks':64,'oriented_isomorphism_classes':len(classes),'constant_odd_isomorphism_classes':32,'nonconstant_isomorphism_classes':32,'undirected_automorphisms':[list(p) for p in automorphisms],
 'sharp_seven_core_witness':{'edges':A7,'orientation_mask':52,'host_out_bits':H7,'counts_transitive_hostile':seven_counts},'records':records,'scope':'Exact controls on two fixed eight-vertex hosts and all128 orientations; universal statement comes from proof.'},indent=2))
