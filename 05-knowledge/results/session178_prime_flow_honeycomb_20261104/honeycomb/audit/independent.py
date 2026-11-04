"""Independent direct-coordinate and edge-set audit, active under python -O."""
import json
from pathlib import Path

def check(x,msg):
 if not x:raise RuntimeError(msg)

def nb(v):
 x,y=v
 return {(x-1,y),(x+1,y),(x,y+(-1 if (x+y)&1 else 1))}

def phi(n):
 a=(n&-n).bit_length()-1
 return a,(n>>a)-1

def layout(d):
 # Coordinates from the closed signed-binary expansion, not parent mutation.
 branches={}
 for t in range(d+1):
  for j in range(2**t):
   word=f'{j:0{t}b}' if t else ''
   x=1+sum((1 if bit=='1' else -1)*2**(d-i) for i,bit in enumerate(word))
   branches[word]=(x,2*t)
 paths={}
 for word,(X,Y) in branches.items():
  if not word:continue
  x,y=branches[word[:-1]];end=X-1
  row=sorted(range(min(x,end),max(x,end)+1),reverse=end<x)
  paths[word]=[(u,y) for u in row]+[(end,y+1),(X,y+1),(X,Y)]
 return branches,paths

out={'chart_checks':0,'odd_formula_checks':0,'embeddings':[],'ball_shells':[]}
for n in range(1,200001):
 a,y=phi(n);check((1<<a)*(y+1)==n,'inverse chart')
 if a:
  check(phi(n//2)==(a-1,y),'halving chart');check((a-1,y) in nb((a,y)),'halving adjacency')
 else:
  dest=phi(3*n+1);A,z=dest
  check(A>=1 and z>=0 and z%2==0 and (1<<A)*(z+1)==3*y+4,'odd detector')
  if y%4==0:check(A>=2 and 4*z<=3*y,'contracting formula')
  else:check(A==1 and 2*z==3*y+2,'expanding formula')
  out['odd_formula_checks']+=1
 out['chart_checks']+=1
for d in range(12):
 B,P=layout(d);allb=set(B.values());check(len(allb)==2**(d+1)-1,'branch injectivity')
 ownership={};edges=set();adj={}
 for w,path in P.items():
  check(len(path)==len(set(path)),'simple wire')
  check(path[0]==B[w[:-1]] and path[-1]==B[w],'endpoints')
  D=2**(d-len(w)+1)
  check(len(path)-1==D+(4 if w[-1]=='0' else 2),'wire length')
  for v in path[1:-1]:
   check(v not in allb and v not in ownership,'interior collision');ownership[v]=w
  for u,v in zip(path,path[1:]):
   check(v in nb(u) and u in nb(v),'host adjacency')
   e=tuple(sorted((u,v)));check(e not in edges,'repeated edge');edges.add(e)
   adj.setdefault(u,set()).add(v);adj.setdefault(v,set()).add(u)
 check(len(ownership)==(d+2)*2**(d+1)-4,'interior count')
 for w,v in B.items():
  expected=0 if d==0 else (2 if not w else 1 if len(w)==d else 3)
  check(len(adj.get(v,()))==expected,'branch degree')
 for v in ownership:check(len(adj[v])==2,'wire degree')
 total=len(allb)+len(ownership)
 check(len(edges)==total-1,'tree Euler count')
 # A connected acyclic layout follows from one path per tree edge, plus above
 # disjointness; additionally verify connectedness directly.
 seen={B['']};front=list(seen)
 while front:
  v=front.pop()
  for u in adj.get(v,()):
   if u not in seen:seen.add(u);front.append(u)
 check(len(seen)==total,'connected wire graph')
 out['embeddings'].append({'depth':d,'branches':len(B),'interiors':len(ownership),'edges':len(edges)})
seen={(0,0)};front=set(seen)
for r in range(1,51):
 layer=set().union(*(nb(v) for v in front))-seen
 seen|=layer;front=layer
 check(len(layer)==3*r,'honeycomb shell count')
 check(len(seen)==1+3*r*(r+1)//2,'honeycomb ball formula')
 out['ball_shells'].append(len(layer))
# Inverse tree built as paths, checking distinct labels and ordinary forward edges.
labels={5};front=[5];levels=[1]
for depth in range(1,22):
 new=[]
 for n in front:
  candidates=[2*n]
  if (n-1)%3==0 and ((n-1)//3)%2==1:candidates.append((n-1)//3)
  for m in candidates:
   check(m>0 and (m//2 if m%2==0 else 3*m+1)==n,'inverse arrow')
   check(m not in labels,'inverse tree collision');labels.add(m);new.append(m)
 front=new;levels.append(len(new))
check(len(labels)==704,'inverse ball cardinality')
check(1+3*21*22//2==694<704,'unit dilation obstruction')
check(704<=1+3*42*43//2,'dilation2 not obstructed by volume')
out['inverse_level_counts']=levels;out['inverse_ball']=len(labels);out['status']='PASS'
Path(__file__).with_name('RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
