#!/usr/bin/env python3
"""Exact cyclic-bank incidence controls; no quotient by a guessed graph group."""
from collections import Counter,defaultdict
from itertools import combinations,permutations
from pathlib import Path
import argparse,json

checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)

def compose(p,q):return tuple(p[q[i]] for i in range(len(p)))
def inverse(p):return tuple(p.index(i) for i in range(len(p)))
def powers(p):
 x=tuple(range(len(p)));out=[]
 while True:
  out.append(x);x=compose(p,x)
  if x==out[0]:return out

class Cube:
 def __init__(self,n):
  self.n=n;self.edges=list(combinations(range(n),2));self.ix={e:i for i,e in enumerate(self.edges)}
  self.E=len(self.edges);self.mask=(1<<self.E)-1
 def reduce(self,rows):
  B={}
  for x in rows:
   for p in sorted(B):
    if x>>p&1:x^=B[p]
   v=x&self.mask
   if not v:
    if x:return None
    continue
   p=(v&-v).bit_length()-1
   for j in list(B):
    if B[j]>>p&1:B[j]^=x
   B[p]=x
  return tuple(B[p] for p in sorted(B))
 def fixed(self,g):
  return self.reduce((1<<i)^(1<<self.ix[tuple(sorted((g[a],g[b])))])^((g[a]>g[b])<<self.E)
                     for i,(a,b) in enumerate(self.edges))
 def path(self,word):
  rows=[];keep=0
  for a,b in zip(word,word[1:]):
   i=self.ix[tuple(sorted((a,b)))];keep|=1<<i;rows.append((1<<i)|((a<b)<<self.E))
  return tuple(rows),self.mask^keep
 def translate(self,F,v):return tuple(x^(((x&v).bit_count()%2)<<self.E) for x in F)
 def members(self,F):
  pivots=[((x&self.mask)&-(x&self.mask)).bit_length()-1 for x in F]
  free=[j for j in range(self.E) if j not in pivots];out=set()
  for y in range(1<<len(free)):
   x=sum(((y>>j)&1)<<e for j,e in enumerate(free))
   for e,row in zip(pivots,F):
    x|=(((row>>self.E)^((row&x).bit_count()%2))&1)<<e
   out.add(x)
  return out
 def edge_orbits(self,gens):
  parent=list(range(self.E))
  def root(a):
   while a!=parent[a]:parent[a]=parent[parent[a]];a=parent[a]
   return a
  for g in gens:
   for i,(a,b) in enumerate(self.edges):
    j=self.ix[tuple(sorted((g[a],g[b])))];parent[root(i)]=root(j)
  return len({root(i) for i in range(self.E)})

def all_cyclic(p):
 out={}
 for tail in permutations(range(1,p)):
  w=(0,)+tail;g=[0]*p
  for a,b in zip(w,w[1:]+w[:1]):g[a]=b
  key=min(powers(tuple(g))[1:]);out[key]=key
 return sorted(out)

def cyclic_control(p,full=False):
 C=Cube(p);gens=all_cyclic(p);parents=[];charts=[]
 for g in gens:
  F=C.fixed(g);need(F is not None and C.E-len(F)==(p-1)//2,'cyclic bank dimension')
  group=powers(g);by_parent=set();subcharts=[]
  words=[]
  if full:words=list(permutations(range(p)))
  else:
   for step in group[1:]:
    for start in range(p):
     word=[start]
     while len(word)<p:word.append(step[word[-1]])
     words.append(word)
  for word in words:
    pin,t=C.path(word);source=C.reduce(F+pin)
    if source is None:continue
    chart=C.translate(source,t);members=C.members(chart)
    need(len(members)<=1<<((p-3)//2),'native chart dimension upper bound')
    by_parent|=members;subcharts.append((chart,members))
  parents.append(by_parent);charts.append(subcharts)
 cross_points=set();pairs=0;max_meet=0
 for i,j in combinations(range(len(gens)),2):
  need(C.edge_orbits((gens[i],gens[j]))==1,'distinct prime cyclic groups generate one edge orbit')
  cross_points|=parents[i]&parents[j]
  # Cheap independent incidence oracle: list all chart memberships of each
  # point; no pair of charts from different parents may share two points.
  left=defaultdict(list);right=defaultdict(list)
  for k,(_,S) in enumerate(charts[i]):
   for x in S:left[x].append(k)
  for k,(_,S) in enumerate(charts[j]):
   for x in S:right[x].append(k)
  overlaps=Counter((a,b) for x in left.keys()&right.keys() for a in left[x] for b in right[x])
  need(not overlaps or max(overlaps.values())<=1,'native chart singleton cross-intersections')
  pairs+=len(overlaps);max_meet=max(max_meet,max(overlaps.values(),default=0))
 multiplicity=Counter(x for S in parents for x in S)
 correction=sum(v-1 for v in multiplicity.values() if v>1)
 need(sum(map(len,parents))-correction==len(multiplicity),'point correction union count')
 need({x for x,v in multiplicity.items() if v>1}==cross_points,'complete collision point support')
 return dict(p=p,path_menu='all' if full else 'constant-step',parents=len(gens),charts=sum(map(len,charts)),labelled_union=len(multiplicity),
             parent_mass=sum(map(len,parents)),collision_points=len(cross_points),correction=correction,
             nonempty_cross_chart_pairs=pairs,max_cross_chart_intersection=max_meet,
             multiplicity_histogram=dict(sorted(Counter(multiplicity.values()).items())))

def composite_control(n):
 C=Cube(n);a=tuple((i+1)%n for i in range(n));v=tuple((i+3)%n if i%3==0 else i for i in range(n))
 b=compose(v,compose(a,inverse(v)))
 need(b not in powers(a),'distinct composite cyclic groups')
 F=C.reduce(C.fixed(a)+C.fixed(b));pin,t=C.path(tuple(range(n)));G=C.reduce(F+pin)
 need(G is not None,'composite common legal chart')
 dim=C.E-len(G);members=C.members(C.translate(G,t))
 need(len(members)==1<<dim,'composite native collision dimension')
 need(dim==(n//3-1)//2,'explicit composite dimension formula')
 need(C.edge_orbits((a,b))==(n//3+1)//2,'independent edge-orbit count')
 # Direct adjacency test, independent of the linear fixed-set equations.
 for x in C.members(G):
  for gen in (a,b):
   for i,(u,w) in enumerate(C.edges):
    j=C.ix[tuple(sorted((gen[u],gen[w])))];need(((x>>i)&1)==(((x>>j)&1)^(gen[u]>gen[w])),'direct invariance')
 return dict(n=n,first_generator=a,conjugator=v,second_generator=b,
             generated_edge_orbits=C.edge_orbits((a,b)),native_intersection_dimension=dim,
             native_intersection_size=len(members),two_witness_masks=sorted(members)[:2])

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path(__file__).with_name('intersections.json'));args=parser.parse_args()
 result=dict(prime_controls=[cyclic_control(p) for p in (3,5,7)]+[cyclic_control(5,True)],
             composite_hostiles=[composite_control(n) for n in (9,15,21,27)])
 result['checks']=checks
 result['scope']='all cyclic groups at 3,5,7; all constant-step native charts and every legal path at 5; exact 9,15,21,27 hostiles'
 args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps(result,indent=2))
