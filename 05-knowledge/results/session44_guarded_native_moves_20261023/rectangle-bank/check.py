#!/usr/bin/env python3
"""Direct certificate checker for universal native rectangle obstructions.

The three image witnesses are checked by literal edge multiplication, not by
trusting the generator, a claimed Smith basis, or vertex-signature completeness.
"""
import argparse,csv,gzip,json,sys
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
HERE=Path(__file__).resolve().parent;COUNT=0

def need(b,why):
 global COUNT
 COUNT+=1
 if not b:raise RuntimeError(why)

def dot(a,b):return sum(x*y for x,y in zip(a,b))
def small_controls():
 cases=0
 for m in [4,8]:
  bitpatterns=range(64) if m==4 else [15,31,63,41,45]
  for bits in bitpatterns:
   A=[[0]*4 for _ in range(4)];adj=[set() for _ in range(4)]
   for k,(i,j) in enumerate(combinations(range(4),2)):
    if not bits>>k&1:continue
    w=1+(k%2);A[i][i]+=w;A[j][j]+=w;A[i][j]-=w;A[j][i]-=w;adj[i].add(j);adj[j].add(i)
   reached={0};todo=[0]
   while todo:
    for v in adj[todo.pop()]:
     if v not in reached:reached.add(v);todo.append(v)
   if len(reached)<4:continue
   vec=list(product(range(m),repeat=4));image={tuple(dot(row,x)%m for row in A) for x in vec}
   H=[x for x in vec if all(dot(row,x)%m==0 for row in A)]
   for c in range(m):
    u=(1,-1,0,0);v=(0,0,1,-1)
    criterion=tuple(c*z%m for z in u) in image and tuple(c*z%m for z in v) in image
    literal=all(all(c*(u[i]*dot(v,x)+v[i]*dot(u,x))%m==0 for i in range(4)) for x in H)
    need(criterion==literal,'rectangle all-bank criterion versus literal finite ring')
    cases+=1
 parity=0
 for k in range(1,7):
  for b in range(k+1):
   for s in range(k):
    for odd in range(1,1<<(k-s),2):
     for a in range(-5,6):
      power=max(0,k-b-s)
      need(((((1<<b)*(1<<s)*odd*a)%(1<<k)==0)==(a%(1<<power)==0)),'valuation divisibility has exact boundary')
      parity+=1
 return dict(rectangle_ring_cases=cases,valuation_cases=parity)

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo',type=Path);p.add_argument('--out',type=Path);a=p.parse_args()
 repo=a.repo or next((q for q in HERE.parents if (q/'05-knowledge/results/metagraph_order8_20261023/merged_edges.tsv.gz').is_file()),None)
 if repo is None:raise RuntimeError('Supply --repo containing the original pinned native atlas.')
 pins=json.loads((HERE/'SOURCE_PINS.json').read_text());blobs={}
 for rec in pins['files']:
  raw=(repo/rec['path']).read_bytes();need(sha256(raw).hexdigest()==rec['sha256'],'source hash '+rec['path']);blobs[rec['key']]=raw
 controls=small_controls();N=3528;root=3527
 edges=[tuple(int(r[k]) for k in ['source','target','blue','black']) for r in csv.DictReader(gzip.decompress(blobs['edges']).decode().splitlines(),delimiter='\t')]
 def apply(x):
  out=[0]*N
  for i,j,b,k in edges:
   z=(b+k)*(x[i]-x[j]);out[i]+=z;out[j]-=z
  return out
 D=json.loads(blobs['electrical'])['order8'];Y=D['radical_lifts_mod4'];X=D['higher_layers'][-1]['lifts'][0];Z=D['eighth_radical_lifts_mod8'][0]
 chars=[[2*v%8 for v in Y[i]] for i in [0,1]]
 for char in chars:need(len(char)==N and char[root]==0 and all(v%8==0 for v in apply(char)),'actual higher-order harmonic separator')
 need(all(v%256==0 for v in apply(X)),'original primitive potential harmonic modulo256')
 need(all(v%8==0 for v in apply(Z)),'actual primitive mod8 character for even-net hostile')
 need((Z[0]-Z[5])%8==7 and 4*(Z[0]-Z[5])%8==4,'even net target amplitudes alone are insufficient')
 charges=[(1238,2246),(881,2762),(2650,2843)]
 witness=json.loads((HERE/'image_witnesses.json').read_text())
 need(witness['input_sha256']==sha256(blobs['edges']).hexdigest(),'generator and checker use identical edge bytes')
 orders=[]
 for i,(pair,entry) in enumerate(zip(charges,witness['solutions'])):
  need(list(pair)==entry['charge'] and entry['coefficient']==4,'image witness names exact signed charge')
  z=entry['potential_mod8'];need(len(z)==N and z[root]==0 and all(type(v) is int and 0<=v<8 for v in z),'full grounded modular image witness')
  image=apply(z);expected=[0]*N;expected[pair[0]]=4;expected[pair[1]]=-4
  need(all((x-y)%8==0 for x,y in zip(image,expected)),'direct entire Laplacian image certificate4d')
  need(entry['integer_correction_div8_nonzero']==[[j,(x-expected[j])//8] for j,x in enumerate(image) if x!=expected[j]],'exact integer correction ledger')
  char=chars[0 if i<2 else 1];difference=(char[pair[0]]-char[pair[1]])%8
  need(difference==2 and 2*difference%8==4,'2d has nonzero pairing with an authenticated harmonic')
  orders.append(dict(charge=pair,exact_order=4,upper_coefficient=4,lower_separator_mod4_lift_index=0 if i<2 else 1,separator_difference=difference,twice_charge_pairing=4,image_potential_sha256=sha256(json.dumps(z,separators=(',',':')).encode()).hexdigest()))
 # Sparse order-two characters cannot detect an even black-pair coefficient.
 sparse=[]
 for node in [1238,2246,881,2762,2650,2843]:
  y=[0]*N;y[node]=4;need(all(v%8==0 for v in apply(y)),'sparse native4e character exists')
  need(all((2*(y[u]-y[v]))%8==0 for u,v in charges),'factor2 kills these sparse order-two separators')
  sparse.append(node)
 # Compare original frozen character and genuinely different full-bank character.
 rect=[]
 for targets in [(881,2762),(2650,2843)]:
  s,r=1238,2246;t,u=targets
  def delta(char,c):return [c*(char[t]-char[u]),-c*(char[t]-char[u]),c*(char[s]-char[r]),-c*(char[s]-char[r])]
  old=delta(X,2);new=delta(chars[0],2)
  need(all(z%8==0 for z in old),'incoming single-character phase guard passes')
  need(any(z%8 for z in new),'universal bank rejects amplitude1')
  need(new[2]%8==new[3]%8==4,'source-side harmonic difference obstructs every target rectangle')
  colored={(i,j):(b,k) for i,j,b,k in edges}
  expected=[4,2,4,4] if targets==(881,2762) else [4,2,4,6]
  observed=[colored[min(a,b),max(a,b)] for a,b in [(s,t),(s,u),(r,t),(r,u)]]
  need(observed==[(0,w) for w in expected],'actual pure-black paired conductances retain factor2')
  rect.append(dict(pure_black_conductances=expected,sources=[s,r],targets=targets,black_pair_coefficient=2,original_frozen_defect=old,other_harmonic_defect=new,other_harmonic_defect_mod8=[z%8 for z in new],exact_allowed_amplitudes='2 divides a',minimum_positive_electrical_amplitude=2))
 # Independently authenticate the scarce actual capacities, not only aggregate weights.
 nodes=list(csv.DictReader(blobs['nodes'].decode().splitlines(),delimiter='\t'));owner=array('H');owner.frombytes(gzip.decompress(blobs['owner']))
 if sys.byteorder!='little':owner.byteswap()
 merged={int(r[k]):int(r['merged_id']) for r in nodes for k in ['class_a','class_b']};chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1
 source=1238;canon=int(nodes[source]['canon_a']);A=[[0]*8 for _ in range(8)]
 for bit,(i,j) in enumerate(combinations(range(8),2)):A[i][j]=(canon>>(27-bit))&1;A[j][i]=1-A[i][j]
 paths=[]
 def dfs(word,remaining):
  if not remaining:paths.append(tuple(word));return
  for v in sorted(remaining):
   if not word or A[word[-1]][v]:dfs(word+[v],remaining-{v})
 dfs([],set(range(8)));need(len(paths)==409,'complete actual source1238 path fibre')
 targetpaths={t:[] for t in [881,2762,2650,2843]}
 for word in paths:
  q=sum(A[word[i]][word[j]]<<bit for bit,(i,j) in enumerate(chords));need(merged[owner[q]]==source,'every source path authenticated in original owner gauge')
  t=merged[owner[full^q]]
  if t in targetpaths:targetpaths[t].append(list(word))
 need(len(targetpaths[2762])==len(targetpaths[2843])==1,'exact distinct-port amplitude cap1 on both incoming rectangles')
 incoming=json.loads(blobs['incoming_receipt'])['paired_selected_port_switches']['rectangles']
 need(all(r['maximum_distinct_actual_observer_amplitude']==1 for r in incoming),'incoming literal capacity scope agrees')
 # Any unit residual-capacity net vector on these four targets is killed,
 # including overlaps/cancellation represented at the NET target level.
 net_cases=0
 for v in product((-1,0,1),repeat=4):
  if sum(v):continue
  target_defect=[4*a%8 for a in v]
  need(all(d==0 for d in target_defect)==(not any(v)),'unit-capacity net parity obstruction')
  net_cases+=1
 out=dict(status='EXACT universal-bank rectangle obstruction',source_commit=pins['commit'],checks=COUNT,small_controls=controls,original_pair_rows=len(edges),charge_orders=orders,rectangles=rect,sparse_order_two_separator_nodes=sparse,source_difference_character=dict(original_mod4_lift_index=0,source_values=[chars[0][1238],chars[0][2246]],support_size=sum(v!=0 for v in chars[0]),not_sparse=True,potential_sha256=sha256(json.dumps(chars[0],separators=(',',':')).encode()).hexdigest()),scarce_native_paths={str(t):p for t,p in targetpaths.items()},net_unit_capacity_cases=net_cases,even_net_insufficiency=dict(net_target_vector=[[0,2],[5,-2]],harmonic_source='original eighth_radical_lifts_mod8[0]',source_defect_mod8=[4,4],scope='signed algebraic hostile, not an authenticated native port replacement'),conclusion='Each incoming nonzero injective rectangle fails the universal mod8 bank. The same source pair has no nonzero incidence perturbation with all net target amplitudes in[-1,1]. Vanishing guards do not certify native wiring or global injectivity.')
 text=json.dumps(out,indent=2,sort_keys=True)+'\n'
 if a.out:a.out.write_text(text)
 else:print(text,end='')
if __name__=='__main__':main()
