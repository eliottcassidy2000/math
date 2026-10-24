#!/usr/bin/env python3
"""All nonrigid regular eleven classes from small full-orbit families only."""
import collections as C
import hashlib
import itertools as I
import json
from pathlib import Path
import support as a
HERE=Path(__file__).resolve().parent

def isomorphic(A,B):
 n=len(A);p=[-1]*n
 if sorted(r.bit_count() for r in A)!=sorted(r.bit_count() for r in B):return False
 def go(k,used):
  if k==n:return True
  for w in range(n):
   if used>>w&1 or A[k].bit_count()!=B[w].bit_count():continue
   if any(a.arc(A,k,j)!=a.arc(B,w,p[j]) for j in range(k)):continue
   p[k]=w
   if go(k+1,used|1<<w):return True
  p[k]=-1;return False
 return go(0,0)
def certify(t):
 t=tuple(t);a.valid(t);a.require(len(t)==11 and all(r.bit_count()==5 for r in t),'regular eleven literal degree check')
 a.require(a.strong(t),'regular eleven strong check')
 group=a.automorphisms(t);error,sizes=a.defect(t,group)
 a.require(len(group)%2==1 and all(s%2 for s in sizes),'odd full automorphism orbit sizes')
 a.require(error==11*(len(sizes)-1),'regular eleven actual defect')
 return {'rows':t,'full_automorphism_order':len(group),'vertex_orbit_sizes':sizes,'error':error}
def ensure_distinct(records,label):
 for i,record in enumerate(records):
  for other in records[:i]:a.require(not isomorphic(record['rows'],other['rows']),label)

def orbit_partition_restrictions():
 parts=list(a.oddparts(11));by_q=C.defaultdict(list)
 for p in parts:by_q[len(p)].append(p)
 a.require(sorted(by_q)==[1,3,5,7,9,11],'all odd full orbit count possibilities')
 a.require(len(parts)==12 and [len(by_q[q]) for q in (1,3,5,7,9,11)]==[1,4,3,2,1,1],'all twelve odd orbit partitions by multiplicity')
 # Explicit integer feasibility controls for each discarded orbit partition.
 fixed7=[d for d in range(4) if any(7*b+d==5 for b in (0,1))]
 a.require(not fixed7,'seven orbit plus four fixed cannot be regular')
 fixed53=[d for d in range(3) if any(5*b+3*c+d==5 for b in (0,1) for c in (0,1))]
 a.require(fixed53==[0,2] and not any(sum(d)==3 for d in I.product(fixed53,repeat=3)),'five three and three fixed degree parity obstruction')
 a.require(not any(3*k+d==5 for k in range(4) for d in (0,1)),'three triples and two fixed residue obstruction')
 a.require(not any(9*b+d==5 for b in (0,1) for d in (0,1)),'nine orbit and two fixed degree obstruction')
 a.require(not any(7*b+3*c==5 for b in (0,1) for c in (0,1)),'seven three fixed degree obstruction')
 a.require(not any(1+5*b+d==5 for b in (0,1) for d in range(4)),'five three three coprime crossing obstruction')
 fixed5=[d for d in range(6) if any(5*b+d==5 for b in (0,1))]
 a.require(fixed5==[0,5],'five orbit six fixed extremal degree possibilities')
 a.require(all(not all(r.bit_count() in fixed5 for r in a.graph(6,b)) for b in range(1<<15)),'no six vertex tournament has only source or sink vertices')
 return {'odd_orbit_partitions':parts,'q3_survivor':[1,5,5],'q5_survivors':[],'q7_survivor':[1,1,1,1,1,3,3],'q9_survivor':[1,1,1,1,1,1,1,1,3]}

def three_orbits():
 # Every nontrivial stabilizer is a noncoordinate line in F5^2. Generate
 # its invariant degree-three matrices directly, not the whole regular11 set.
 universe=set()
 for slope in range(1,5):
  for residues in I.combinations(range(5),3):
   universe.add(tuple(sum(1<<y for y in range(5) if (y-slope*x)%5 in residues) for x in range(5)))
 a.require(len(universe)==40,'forty line invariant degree three matrices')
 pending=set(universe);records=[]
 while pending:
  m=min(pending);orbit={a.shift(m,u,v) for u in range(5) for v in range(5)}
  a.require(orbit<=pending and len(orbit)==5,'degree three line matrix orbit size five');pending-=orbit
  a.require(all(r.bit_count()==3 for r in m) and all(sum(r>>j&1 for r in m)==3 for j in range(5)),'degree three matrix regular margins')
  stabilizer=sum(a.shift(m,u,v)==m for u in range(5) for v in range(5))
  a.require(stabilizer==5,'degree three presentation stabilizer five')
  record=certify(a.five_graph(m));record['matrix']=m;records.append(record)
 a.require(C.Counter((r['full_automorphism_order'],r['error']) for r in records)=={(5,22):7,(55,0):1},'seven genuine q3 and one transitive hidden symmetry')
 ensure_distinct(records,'eight q3 candidate ordinary classes distinct')
 return records

def seven_orbits():
 records=[]
 for d,k in ((0,4),(1,3)):
  pending={tuple(s) for s in I.combinations(range(5),k)};subsets=[]
  while pending:
   s=min(pending);orbit={tuple(sorted((x+j)%5 for x in s)) for j in range(5)}
   a.require(orbit<=pending and len(orbit)==5,'fixed C5 attachment rotation orbit');pending-=orbit;subsets.append(s)
  a.require(len(subsets)==(1 if d==0 else 2),'q7 attachment subset orbit count')
  for subset in subsets:
   for sign in ((1,) if d==0 else (1,-1)):
    t=list(a.cyclic(3))+[r<<3 for r in a.cyclic(3)]+[r<<6 for r in a.cyclic(5)]
    for x in range(3):
     for y in range(3):
      u,v=(x,3+y) if d==1 and y==sign*x%3 else (3+y,x);t[u]|=1<<v
     for z in range(5):
      u,v=(x,6+z) if z in subset else (6+z,x);t[u]|=1<<v
      u,v=(6+z,3+x) if z in subset else (3+x,6+z);t[u]|=1<<v
    record=certify(t);record.update({'crossing_degree':d,'matching_sign':sign if d else None,'attachment_subset':subset});records.append(record)
 a.require(C.Counter((r['full_automorphism_order'],r['error'],r['vertex_orbit_sizes']) for r in records)=={(9,66,(1,1,1,1,1,3,3)):1,(3,66,(1,1,1,1,1,3,3)):4},'five genuine q7 classes')
 ensure_distinct(records,'five q7 candidate ordinary classes distinct')
 return records

def nine_orbits():
 fours={a.canonical_small(a.graph(4,b)) for b in range(64)}
 a.require(len(fours)==4,'four ordinary order four tournament types')
 low=[q for q in fours if max(r.bit_count() for r in q)<=2]
 high=[q for q in fours if min(r.bit_count() for r in q)>=1]
 a.require(sorted(sorted(r.bit_count() for r in q) for q in low)==[[0,2,2,2],[1,1,2,2]],'low block only cycle sink or strong four')
 a.require(sorted(sorted(r.bit_count() for r in q) for q in high)==[[1,1,1,3],[1,1,2,2]],'high block only source cycle or strong four')
 strongs=[q for q in fours if a.strong(q)];a.require(len(strongs)==1,'unique strong four internal type');q=strongs[0]
 a.require(len(a.automorphisms(q))==1,'strong four rigid')
 low_targets=[i for i in range(4) if q[i].bit_count()==1];high_targets=[i for i in range(4) if q[i].bit_count()==2];records=[]
 for p in I.permutations(high_targets):
  t=list(a.cyclic(3))+[r<<3 for r in q]+[r<<7 for r in q]
  for x in range(3):
   for v in range(4):t[3+v]|=1<<x;t[x]|=1<<(7+v)
  matching=dict(zip(low_targets,p))
  for u in range(4):
   for v in range(4):
    i,j=(3+u,7+v) if matching.get(u)==v else (7+v,3+u);t[i]|=1<<j
  record=certify(t);record['crossing_matching']=sorted(matching.items());records.append(record)
 a.require(all(r['full_automorphism_order']==3 and r['error']==88 and r['vertex_orbit_sizes']==(1,1,1,1,1,1,1,1,3) for r in records),'two genuine q9 classes')
 a.require(len(records)==2,'two q9 crossing matchings');ensure_distinct(records,'two q9 candidate ordinary classes distinct')
 return records

def transitive_classes():
 universe={tuple(sorted(d if bits>>(d-1)&1 else 11-d for d in range(1,6))) for bits in range(32)}
 a.require(len(universe)==32,'thirty two cyclic degree five connection sets')
 pending=set(universe);records=[]
 while pending:
  s=min(pending);orbit={tuple(sorted(u*d%11 for d in s)) for u in range(1,11)}
  a.require(orbit<=pending,'circulant multiplier orbit partition');pending-=orbit
  t=tuple(sum(1<<((i+d)%11) for d in s) for i in range(11));record=certify(t)
  record.update({'connection_set':s,'multiplier_orbit_size':len(orbit)});records.append(record)
 a.require(C.Counter((r['full_automorphism_order'],r['error'],r['vertex_orbit_sizes']) for r in records)=={(11,0,(11,)):3,(55,0,(11,)):1},'four vertex transitive regular eleven classes')
 ensure_distinct(records,'four multiplier orbits remain distinct under full isomorphism')
 return records

def run():
 pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
 for item in pins['sources']:a.require(hashlib.sha256((HERE/item['archived_path']).read_bytes()).hexdigest()==item['sha256'],'source pin')
 restrictions=orbit_partition_restrictions();q3=three_orbits();q7=seven_orbits();q9=nine_orbits();q1=transitive_classes()
 hidden=next(r for r in q3 if not r['error'])
 matches=[r for r in q1 if isomorphic(r['rows'],hidden['rows'])]
 a.require(len(matches)==1 and matches[0]['full_automorphism_order']==55,'hidden q3 control equals unique transitive Aut55 class')
 selected=q1+[r for r in q3 if r['error']]+q7+q9
 a.require(len(selected)==18,'eighteen nonrigid regular eleven classes')
 ensure_distinct(selected,'all eighteen nonrigid ordinary classes distinct')
 histogram=C.Counter(r['error'] for r in selected)
 a.require(histogram=={0:4,22:7,66:5,88:2},'complete nonrigid regular eleven defect spectrum')
 pointed=sum(len(r['vertex_orbit_sizes']) for r in selected);error=sum(r['error'] for r in selected)
 a.require(pointed==78 and error==660,'nonrigid pointed mass78 and error660')
 a.require(11*len(selected)-pointed==120,'nonrigid unpointing deficit120')
 a.require(19480-120==19360 and 9166-120==9046,'transfer constants from pinned nonregular audit')
 return {'status':'PASS','scope':'Exhaustive nonrigid regular eleven classification from full automorphism orbit partitions. No full regular11 census; rigid class count and almost-regular10 parent count remain uncomputed here.','orbit_partition_restrictions':restrictions,'q3_candidates':q3,'q7_classes':q7,'q9_classes':q9,'q1_classes':q1,'nonrigid_classes':len(selected),'nonrigid_error_histogram':sorted(histogram.items()),'nonrigid_pointed_states':pointed,'nonrigid_error':error,'nonrigid_unpointing_deficit':120,'formulas':{'regular11_class_count':'(P10+120)/11','regular11_error':'10*P10-120','strong11_error':'10*P10+9046','total11_error':'10*P10+19360','total11_class_count':'(U11-10*P10-19360)/11','rigid_regular11_class_count':'(P10-78)/11','P10_definition':'ordinary almost-regular ten-vertex classes with degrees4^5,5^5'},'checks':dict(sorted(a.CHECKS.items()))}
if __name__=='__main__':print(json.dumps(run(),sort_keys=True,indent=2))
