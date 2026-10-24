#!/usr/bin/env python3
"""Independent literal-record and topology audit, importing no producer code."""
import json,sys,hashlib,itertools
from pathlib import Path
from collections import Counter
import numpy as np
P=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent/'producer'
OUT=Path(__file__).resolve().parent
checks=Counter()
def need(x,s):
 checks[s]+=1
 if not x:raise RuntimeError(s)
def load(n):return json.loads((P/n).read_text())
n=10; pairs=[(i,j) for i in range(n) for j in range(i+2,n)]; full=(1<<36)-1

for pin in load('SOURCE_PINS.json')['files']:need(hashlib.sha256((P/pin['path']).read_bytes()).hexdigest()==pin['sha256'],'source pin SHA256')
need(hashlib.sha256((P/'sources/regular-eleven-MANIFEST.json').read_bytes()).hexdigest()=='0462baf160610115c8b2436c2c57cce1964d0439ce49e86b9b43ce00f2bda3f2','exact inherited atlas manifest')
need(hashlib.sha256((P/'sources/regular-eleven-eleven-O2-parents.txt').read_bytes()).hexdigest()=='ef22a13423ecf41c90e0d01e3f04bd1a42e173bb519abcad7698fb4aa3979937','exact inherited atlas rows')
D=load('normal.json'); nodes={r['atlas_index']:r for r in D['nodes']}; targets={r['boundary_class_id']:r for r in D['boundary_target_classes']}
def decode(w):
 a=[1<<(i+1) if i<9 else 0 for i in range(10)]
 for k,(u,v) in enumerate(pairs):
  if w>>k&1:a[u]|=1<<v
  else:a[v]|=1<<u
 return a
def norm(a,p):return [sum(((a[u]>>v)&1)<<j for j,v in enumerate(p)) for u in p]
def reverse(a):return [1023^r^(1<<i) for i,r in enumerate(a)]
def colors(a):
 c=[r.bit_count() for r in a]
 while True:
  s=[(c[i],tuple(sum((a[i]>>j)&1 for j in range(10) if c[j]==k) for k in sorted(set(c)))) for i in range(10)]
  names={v:i for i,v in enumerate(sorted(set(s)))};d=[names[v] for v in s]
  if len(set(c))==len(set(d)):return d
  c=d
def mappings(a,b):
 ca,cb=colors(a),colors(b)
 if Counter(ca)!=Counter(cb):return []
 order=sorted(range(10),key=lambda v:(ca.count(ca[v]),v));used=set();p={};out=[]
 def rec(k):
  if k==10:out.append(tuple(p[i] for i in range(10)));return
  v=order[k]
  for w in range(10):
   if w in used or ca[v]!=cb[w]:continue
   if any(((a[v]>>u)&1)!=((b[w]>>p[u])&1) for u in p):continue
   p[v]=w;used.add(w);rec(k+1);used.remove(w);del p[v]
 rec(0);return out
# Enumerate reflection assignments directly in parallel arrays, not Gray traversal.
seen=set();orbits=[]
for k,(i,j) in enumerate(pairs):
 if k in seen:continue
 orb={k,pairs.index((9-j,9-i))};seen|=orb;orbits.append(sorted(orb))
x=np.arange(1<<len(orbits),dtype=np.uint32);degrees=np.tile(np.array([r.bit_count() for r in decode(0)],dtype=np.int16),(len(x),1));words=np.zeros(len(x),dtype=np.uint64)
for bit,orb in enumerate(orbits):
 b=((x>>bit)&1).astype(np.int16)
 for k in orb:
  i,j=pairs[k];degrees[:,i]+=b;degrees[:,j]-=b;words|=b.astype(np.uint64)<<k
cube=set(map(int,words[np.all((degrees==4)|(degrees==5),axis=1)]));ports={r[0]:r for r in D['ports']}
need(len(ports)==len(D['ports']),'unique literal words');need(set(ports)==cube,'independent binary reflection cube exact');need(len(orbits)==20,'reflection dimension')
# Independently find every self-converse atlas class by generic arc-compatible bijections.
tokens=iter((P/'sources/regular-eleven-eleven-O2-parents.txt').read_text().split());count=int(next(tokens));sc=set();atlas=[]
for idx in range(count):
 a=[int(next(tokens)) for _ in range(10)];g=int(next(tokens));group=[tuple(int(next(tokens)) for _ in range(10)) for _ in range(g)];atlas.append(a)
 anti=mappings(a,reverse(a))
 if anti:
  sc.add(idx);need(idx in nodes,'all self converse atlas classes represented');r=nodes[idx];need(r['canonical_rows']==a,'node exact atlas rows');aut=mappings(a,a);need(set(aut)==set(group),'independent full automorphism group');need(len(aut)==r['automorphism_order'],'node automorphism order');c=colors(a);orbits={tuple(sorted({q[v] for q in aut})) for v in range(10)};need(sum(c.count(c[o[0]])-len(o) for o in orbits)==r['refinement_error'],'independent node refinement error');need(set(tuple(q) for q in r['anti_involutions'])=={q for q in anti if all(q[q[v]]==v for v in range(10))},'all anti involutions')
need(next(tokens,None) is None,'atlas consumed');need(count==13333 and sc==set(nodes) and len(sc)==229,'all 229 SC atlas classes exact')
adj=Counter();badj=Counter();internal=Counter();boundary=Counter();wordcount=Counter()
for word,idx,path,target,kind,tid,order in ports.values():
 a=decode(word);u=decode(target);need(sorted(path)==list(range(10)),'source path permutation');need(norm(nodes[idx]['canonical_rows'],path)==a,'literal source path');need(target==full^word,'literal C complement');need(sorted(order)==list(range(10)),'target permutation');wordcount[idx]+=1
 if target in ports:
  need(kind=='internal' and tid==ports[target][1],'internal target owner');need(norm(nodes[tid]['canonical_rows'],order)==u,'internal literal target path');internal[idx]+=1
  if word<target:adj[tuple(sorted((idx,tid)))]+=1
 else:
  need(kind=='boundary' and tid in targets,'boundary target owner');need(norm(u,order)==targets[tid]['canonical_rows'],'boundary literal canonical order');need(sorted(r.bit_count() for r in u)==[3,4,4,4,4,5,5,5,5,6],'boundary scores');boundary[idx]+=1;badj[idx,tid]+=1
need(adj==Counter({(a,b):m for a,b,m in D['adjacency_multiplicities']}),'all adjacency multiplicities');need(badj==Counter({(a,b):m for a,b,m in D['boundary_adjacency']}),'all boundary multiplicities')
need(D['blue_ports']==len(ports) and D['internal_blue_ports']==sum(internal.values()) and D['boundary_blue_edges']==sum(boundary.values()) and D['internal_blue_edges']==sum(adj.values()) and D['boundary_classes']==len(targets),'reported scalar totals')
for i,r in targets.items():need(r['edges_from_balanced_layer']==sum(v for (a,b),v in badj.items() if b==i),'boundary target incidence totals')
need(sum(adj.values())==11452 and sum(m for (a,b),m in adj.items() if a==b)==24,'11452 lines 24 loops');need(sum(boundary.values())==18795 and len(targets)==786,'18795 boundary 786 targets')
for i,r in nodes.items():need((r['blue_presentation_words'],r['internal_degree'],r['boundary_edges'])==(wordcount[i],internal[i],boundary[i]),'node port conservation')
# Independently verify target normal forms unique up to isomorphism within stable-color signature buckets.
buckets={}
for i,r in targets.items():
 a=r['canonical_rows'];c=colors(a);sig=tuple(sorted((c[v],a[v].bit_count(),tuple(sum((a[v]>>j)&1 for j in range(10) if c[j]==k) for k in sorted(set(c)))) for v in range(10)))
 for j in buckets.get(sig,[]):need(not mappings(a,targets[j]['canonical_rows']),'boundary classes pairwise nonisomorphic')
 buckets.setdefault(sig,[]).append(i);aut=mappings(a,a);need(len(aut)==r['automorphism_order'],'target complete group');need(bool(mappings(a,reverse(a))),'boundary self converse');need(all(len({q[v] for q in aut})==c.count(c[v]) for v in range(10)) and r['error']==0,'boundary zero refinement error')
# Small spanning 2-connected certificate: triangle followed by 2-edge ears.
nei={v:set() for v in nodes}
for a,b in adj:
 if a!=b:nei[a].add(b);nei[b].add(a)
positive={i for i,r in nodes.items() if r['refinement_error']>0};partition=Counter()
for (a,b),m in adj.items():partition[int(a in positive)+int(b in positive)]+=m
need(partition=={2:90,1:1857,0:9505},'core edge partition');need(all(nei[v]&positive for v in set(nodes)-positive),'positive core dominates every zero-error node')
start=min(positive);reached={start};todo=[start]
while todo:
 for w in nei[todo.pop()]&positive:
  if w not in reached:reached.add(w);todo.append(w)
need(reached==positive and len(positive)==20,'positive core connected')
triangle=next((a,b,c) for a in sorted(nei) for b in sorted(nei[a]) for c in sorted(nei[a]&nei[b]) if a<b<c)
used=set(triangle);ears=[]
while len(used)<len(nodes):
 candidates=[v for v in sorted(nodes) if v not in used and len(nei[v]&used)>=2];need(bool(candidates),'two-edge ear extension exists');v=candidates[0];a,b=sorted(nei[v]&used)[:2];ears.append([a,v,b]);used.add(v)
cert={'scope':'blue-only induced graph, not full metagraph','triangle':triangle,'two_edge_ears':ears,'vertices':229,'certificate_edges':3+2*len(ears),'conclusion':'connected, no articulation vertices, no bridges'}
(OUT/'topology-certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
# Direct deletion replay of full graph, independent from producer Tarjan.
def reach(omit=None):
 todo=[next(v for v in nei if v!=omit)];seen=set(todo)
 while todo:
  for w in nei[todo.pop()]:
   if w!=omit and w not in seen:seen.add(w);todo.append(w)
 return seen
need(len(reach())==229 and all(len(reach(v))==228 for v in nei),'connected after every vertex deletion');need(D.get('connected_component_count',1)==1 and D.get('component_sizes',[229])==[229] and D['bridge_pairs']==[] and D['articulation_nodes']==[],'reported topology');need(D['components']==[sorted(nodes)],'exact component member list');need(D['cycle_rank']==11452-229+1,'cycle rank')
receipt={'status':'PASS','source_check_sha256':hashlib.sha256((P/'check.py').read_bytes()).hexdigest(),'source_receipt_sha256':hashlib.sha256((P/'normal.json').read_bytes()).hexdigest(),'checks':dict(checks),'SC_classes':229,'ports':len(ports),'internal_lines':sum(adj.values()),'loops':24,'boundary_edges':18795,'boundary_targets':786,'topology_certificate_edges':cert['certificate_edges'],'scope':'complete blue induced multigraph, not full black or full metagraph'}
print(json.dumps(receipt,sort_keys=True,indent=2))
