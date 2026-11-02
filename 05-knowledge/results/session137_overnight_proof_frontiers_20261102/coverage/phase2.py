"""Exact isomorphism-quotient atlas; no inherited code or external class files."""
from explore import tournament,extend,H,strong,literal
from itertools import permutations,product,combinations
from functools import lru_cache
import json,hashlib

def check(x,msg):
 if not x:raise RuntimeError(msg)
def encode(a):
 n=len(a);return sum(((a[i]>>j)&1)<<k for k,(i,j) in enumerate(combinations(range(n),2)))
@lru_cache(maxsize=None)
def canon(n,code):
 a=tournament(n,code);cells={}
 for i in range(n):cells.setdefault(a[i].bit_count(),[]).append(i)
 groups=[list(permutations(cells[d])) for d in sorted(cells)]
 best=None
 for choices in product(*groups):
  p=sum(choices,());c=0;k=0
  for i in range(n):
   for j in range(i+1,n):
    c|=((a[p[i]]>>p[j])&1)<<k;k+=1
  if best is None or c<best:best=c
 return best

def generate(limit):
 rows={0:{0}}
 for n in range(1,limit+1):
  row=set()
  for code in rows[n-1]:
   a=tournament(n-1,code)
   for s in range(1<<(n-1)):row.add(canon(n,encode(extend(a,s))))
  rows[n]=row
  check(len(row)==[1,1,1,2,4,12,56,456][n],'complete class tower control')
  print(json.dumps({'stage':'isomorphism tower','n':n,'classes':len(row)},sort_keys=True),flush=True)
 return rows

def intervals(values):
 out=[]
 for v in sorted(values):
  if out and out[-1][-1]+2==v:out[-1].append(v)
  else:out.append([v])
 return out

def analyse(n,codes):
 vals={c:H(tournament(n,c)) for c in codes};vertices={c for c in codes if strong(tournament(n,c))}
 check(len(vertices)=={5:6,6:35,7:353}[n],'strong class count')
 parent={c:c for c in vertices};forest=[];edgecount=0;sampledp=0
 def root(c):
  while parent[c]!=c:parent[c]=parent[parent[c]];c=parent[c]
  return c
 pairs=list(combinations(range(n),2));ix={ij:k for k,ij in enumerate(pairs)}
 moves=[(1<<k,None,(k,)) for k in range(len(pairs))]
 for r in range(n):
  es=[ix[tuple(sorted((r,j)))] for j in range(n) if j!=r]
  moves += [((1<<k)|(1<<l),r,(k,l)) for k,l in combinations(es,2)]
 witnesses=[]
 for c in sorted(vertices):
  for flip,r,ks in moves:
   d=c^flip;dd=canon(n,d)
   check(dd in vals,'move target exists in complete class tower')
   if dd not in vertices or abs(vals[dd]-vals[c])>2:continue
   edgecount+=1
   if root(c)!=root(dd):
    parent[root(c)]=root(dd)
    check(H(tournament(n,d))==vals[dd],'direct labelled target H')
    check(strong(tournament(n,d)),'direct labelled target strongness')
    forest.append({'source_class_code':c,'target_class_code':dd,'literal_target_code':d,'flipped_pairs':[pairs[k] for k in ks],'source_H':vals[c],'target_H':vals[dd]})
    sampledp+=1
 components={}
 for c in vertices:components.setdefault(root(c),[]).append(c)
 comp=[{'class_count':len(cs),'image':sorted({vals[c] for c in cs}),'class_codes':sorted(cs)} for cs in components.values()]
 comp.sort(key=lambda x:(x['image'],x['class_codes']))
 imax=intervals({vals[c] for c in vertices});images=[x['image'] for x in comp]
 # Distinct components may overlap numerically; compare full images, not only hulls.
 report={'n':n,'strong_classes':len(vertices),'move_presentations':len(vertices)*len(moves),'unit_move_presentations':edgecount,'components':[{'class_count':x['class_count'],'image':x['image']} for x in comp],'global_maximal_odd_intervals':imax,'one_component_per_global_interval':images==imax,'spanning_forest_edges':len(forest)}
 print(json.dumps(report,sort_keys=True,separators=(',',':')),flush=True)
 artifact={'report':report,'components':comp,'spanning_forest':forest,'class_H':{str(c):vals[c] for c in sorted(vertices)}}
 with open('phase2_order%d.json'%n,'w') as f:json.dump(artifact,f,sort_keys=True,separators=(',',':'));f.write('\n')
 return artifact

if __name__=='__main__':
 rows=generate(7);ledger=[analyse(n,rows[n]) for n in [5,6,7]]
 print('PASS semantic_sha256='+hashlib.sha256(json.dumps(ledger,sort_keys=True,separators=(',',':')).encode()).hexdigest())
