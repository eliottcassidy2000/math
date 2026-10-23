#!/usr/bin/env python3
"""Exact port connection certificates. No Git/network actions; all checks explicit."""
import json, hashlib
from pathlib import Path
from collections import defaultdict,Counter
from itertools import product
import sympy as S
from sympy.polys.matrices import DomainMatrix
from inherited_ports import graph, profiles, components, direct_paths, subsets

def require(x,m):
 if not x: raise RuntimeError(m)

def matchings(v):
 if not v: yield (); return
 a=v[0]
 for b in v[1:]:
  for r in matchings(tuple(x for x in v if x not in (a,b))):
   yield ((a,b),)+r

def rank_exact(a):
 if not a or not a[0]: return 0,[],[]
 dm=DomainMatrix.from_Matrix(S.Matrix(a)).to_field()
 r,p=dm.rref()
 _,q=dm.transpose().rref()
 return len(p),list(p),list(q)

def rank2(a):
 rows=[sum(int(v%2)<<i for i,v in enumerate(row)) for row in a]; piv={}
 for r in rows:
  while r:
   i=r.bit_length()-1
   if i in piv:r^=piv[i]
   else:piv[i]=r;break
 return len(piv)

receipt={'status':'FINITE-EXACT','matching':[]}
for k in range(1,5):
 n=2*k; ms=list(matchings(tuple(range(n))))
 a=[[int(len(components(range(n),list(p)+list(q)))==1) for q in ms] for p in ms]
 rank,pc,pr=rank_exact(a)
 # Every off-diagonal pair of states has a concrete near-perfect continuation.
 separators=[]
 for i,p in enumerate(ms):
  for j in range(i):
   q=ms[j]; ab=next(e for e in p if e not in q)
   # Find a matching R missing a,b making Q a single a-b path.
   r=next(r for r in matchings(tuple(v for v in range(n) if v not in ab))
          if len(components(range(n),list(q)+list(r)))==1)
   require(len(components(range(n),list(p)+list(r)))>1,'separator')
   separators.append([i,j,list(ab),r])
 # Path contexts are each matching with one chosen edge removed.
 # Closure map gives exactly the same acceptance as square connectivity matrix.
 for i,p in enumerate(ms):
  for j,q in enumerate(ms):
   for removed in q:
    r=[e for e in q if e!=removed]
    require(int(len(components(range(n),list(p)+r))==1)==a[i][j], 'closure identity')
 require(len(set(tuple(r) for r in a))==len(ms),'distinct rows')
 null=DomainMatrix.from_Matrix(S.Matrix(a)).to_field().nullspace().to_Matrix()
 require(S.Matrix(a)*null.T==S.zeros(len(ms),null.rows),'kernel certificate')
 sub=S.Matrix(a).extract(pr,pc)
 det=DomainMatrix.from_Matrix(sub).det()
 require(det!=0,'nonzero minor')
 item={'ports':n,'states':len(ms),'rank_Q':rank,'rank_F2':rank2(a),
       'matchings':ms,'matrix':a,'pivot_rows':pr,'pivot_columns':pc,
       'minor_determinant':str(det),'nullspace_rows':[[str(v) for v in row] for row in null.tolist()],
       'pairwise_separator_count':len(separators),'separators':separators}
 receipt['matching'].append(item)
 print('matching',n,len(ms),rank,item['rank_F2'],flush=True)

# Independently enumerate local states and all full-graph HPs, then recover
# each solution's unique local restrictions. All other tensor entries are zero:
# accepted gluing is a HP, and every HP was enumerated.
edges=graph(25); sets={'A':[2,7,9,16,18],'B':[4,5,11,12],'C':[1,3,6,8,10,17]}
core=sorted(set(range(1,26))-set(sum(sets.values(),[])))
ce=[e for e in edges if all(v in core for v in e)]; cs=list(subsets(ce))
blocks={}
for name,vs in sets.items():
 internal=[e for e in edges if all(v in vs for v in e)]
 ports=sorted(e if e[0] in vs else e[::-1] for e in edges if sum(v in vs for v in e)==1)
 states=profiles(vs,internal,ports)
 blocks[name]={'ports':ports,'states':states}
paths=direct_paths(25,edges,18); sols=[]
for p in paths:
 pe={tuple(sorted(e)) for e in zip(p,p[1:])}; ids=[]
 for name,b in blocks.items():
  vs=sets[name]; inside=sorted(e for e in pe if all(v in vs for v in e))
  ports=[i for i,e in enumerate(b['ports']) if tuple(sorted(e)) in pe]
  found=[i for i,s in enumerate(b['states']) if sorted(map(tuple,s['internal_edges']))==inside and s['selected_ports']==ports]
  require(len(found)==1,'unique restriction');ids.append(found[0])
 ci=next(i for i,c in enumerate(cs) if set(c)==set(ce)&pe)
 sols.append(ids+[ci])
require(len(sols)==10,'ten HPs')
receipt['G25']={'solutions':sols,'paths':paths,'core':core,'core_edge_subsets':cs,'blocks':{}}
for i,name in enumerate(blocks):
 rows=[set() for _ in blocks[name]['states']]
 for sol in sols: rows[sol[i]].add(tuple(sol[:i]+sol[i+1:]))
 cols=sorted(set.union(*rows)); a=[[int(c in r) for c in cols] for r in rows]
 groups=defaultdict(list)
 for j,r in enumerate(rows):groups[tuple(sorted(r))].append(j)
 rank,pc,pr=rank_exact(a)
 item={'state_count':len(rows),'total_context_count':16*__import__('math').prod(len(blocks[b]['states']) for b in blocks if b!=name),
 'nonzero_columns':cols,'row_supports':[sorted(r) for r in rows],
 'classes':[{'support':k,'states':v} for k,v in groups.items()],
 'class_count_including_dead':len(groups),'rank_Q':rank,
 'pivot_rows':pr,'pivot_columns':pc,'states':blocks[name]['states'],'ports':blocks[name]['ports']}
 receipt['G25']['blocks'][name]=item
 print('G25',name,len(rows),len(groups),rank,flush=True)
# Tensor quotient can merge all these independently, because equivalence was
# computed against every original other-coordinate value, not just a sum.
classes=[receipt['G25']['blocks'][n]['classes'] for n in blocks]
compressed=0
for gs in product(*classes):
 for ci in range(16):
  rep=[g['states'][0] for g in gs]+[ci]
  if rep in sols:compressed+=__import__('math').prod(len(g['states']) for g in gs)
require(compressed==10,'simultaneous weighted quotient')
receipt['G25']['quotient_weighted_count']=compressed

# Restore endpoint identities as polynomial coefficients, not extra structural states.
endpoint_poly=Counter()
for gs in product(*classes):
 if any(not g['support'] for g in gs):continue
 for name,g in zip(blocks,gs):
  require(len({tuple(blocks[name]['states'][j]['selected_ports']) for j in g['states']})==1,'live class same ports')
 for ci in range(16):
  rep=[g['states'][0] for g in gs]+[ci]
  if rep not in sols:continue
  deg=Counter(v for e in cs[ci] for v in e)
  for name,g in zip(blocks,gs):
   state=blocks[name]['states'][g['states'][0]]
   deg.update(blocks[name]['ports'][j][1] for j in state['selected_ports'])
  core_ends=[v for v in core if deg[v]==1]
  for ids in product(*(g['states'] for g in gs)):
   ends=core_ends+sum([blocks[name]['states'][j]['endpoints'] for name,j in zip(blocks,ids)],[])
   endpoint_poly[tuple(sorted(ends))]+=1
require(endpoint_poly==Counter(tuple(sorted((p[0],p[-1]))) for p in paths),'endpoint polynomial')
receipt['G25']['endpoint_polynomial']=[{'endpoints':e,'coefficient':w} for e,w in sorted(endpoint_poly.items())]

# Exact electrical hostile on the same four labels, unit conductances.
def eliminate(es):
 l=S.zeros(4)
 for u,v in es:l[u,u]+=1;l[v,v]+=1;l[u,v]-=1;l[v,u]-=1
 b=[0,2];i=[1,3];a=l.extract(i,i)
 return l.extract(b,b)-l.extract(b,i)*a.inv()*l.extract(i,b),a.det()
p=[(0,1),(1,2),(2,3)];star=[(0,1),(1,2),(1,3)]
a,da=eliminate(p);b,db=eliminate(star)
require(a==b and da==db==2,'same Schur and determinant')
require(sum(len(direct_paths(4,[(u+1,v+1) for u,v in p],start)) for start in range(1,5))//2==1,'path HP count')
require(sum(len(direct_paths(4,[(u+1,v+1) for u,v in star],start)) for start in range(1,5))==0,'star HP count')
receipt['electrical']={'boundary':[0,2],'path_edges':p,'star_edges':star,'response':[[str(x) for x in row] for row in a.tolist()],'interior_determinant':2,'path_HP_unoriented':1,'star_HP_unoriented':0}
receipt['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('inherited_ports.py')]}
Path(__file__).with_name('certificate.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS all exact checks')
