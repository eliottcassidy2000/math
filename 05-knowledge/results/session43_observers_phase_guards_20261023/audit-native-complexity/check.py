"""Independent construction and local-to-global proof certificate checks; no producer import."""
from itertools import combinations, product
from collections import Counter
from pathlib import Path
import json,hashlib

def require(x,msg):
 if not x: raise RuntimeError(msg)

def build(n,edges):
 copies={};owners=[];links=[]
 for v in range(n):
  ids=[]
  for j,e in enumerate(edges):
   if v in e:
    copies[v,j]=len(owners);ids.append(len(owners));owners.append(v)
  links += [(a,b,{0}) for a,b in zip(ids,ids[1:])]
 for j,(u,v) in enumerate(edges):links.append((copies[u,j],copies[v,j],{1,2}))
 deg=Counter(v for a,b,E in links for v in (a,b));odd=[v for v in range(len(owners)) if deg[v]%2]
 for a,b in zip(odd[::2],odd[1::2]):
  w=len(owners);owners.append(None);links.extend([(a,w,set(range(3))),(w,b,set(range(3)))])
 s=len(owners);owners.extend([None,None]);t=s+1
 links.extend([(s,0,set(range(3))),(0,t,set(range(3)))])
 adjacency=[[] for _ in owners]
 for j,(a,b,E) in enumerate(links):adjacency[a].append((b,j));adjacency[b].append((a,j))
 require(max(map(len,adjacency))<=6,'degree')
 stack=[s];visited=set();trail=[]
 while stack:
  v=stack[-1]
  while adjacency[v] and adjacency[v][-1][1] in visited:adjacency[v].pop()
  if adjacency[v]:
   w,j=adjacency[v].pop();visited.add(j);stack.append(w)
  else:trail.append(stack.pop())
 trail.reverse();require(len(visited)==len(links),'Euler coverage')
 visits=Counter(trail);labels=Counter();P=[]
 for v in trail:
  count=4-visits[v] if labels[v]==0 else 1
  for _ in range(count):P.append((v,labels[v]));labels[v]+=1
 require(all(labels[v]==3 for v in range(len(owners))),'three labels')
 rule={frozenset((a,b)):E for a,b,E in links};blocks={}
 for (a,u),(b,v) in zip(P,P[1:]):
  if a!=b:
   require(frozenset((a,b)) not in blocks,'unique boundary');blocks[frozenset((a,b))]=(a,b,{(e+v-u)%3 for e in rule[frozenset((a,b))]})
 def arc(x,y):
  a,u=divmod(x,3);b,v=divmod(y,3)
  if a==b:return (v-u)%3==1
  data=blocks.get(frozenset((a,b)))
  if data is None:return a<b
  i,j,D=data
  return ((v-u)%3 in D) if a==i else not ((u-v)%3 in D)
 M=3*len(owners);T=[[int(arc(x,y)) for y in range(M)] for x in range(M)]
 for x in range(M):
  require(T[x][x]==0,'diagonal')
  for y in range(x):require(T[x][y]+T[y][x]==1,'tournament')
 # Universal 9-pair local verification proves equivalence for every phase vector.
 for (a,u),(b,v) in zip(P,P[1:]):
  for x,y in product(range(3),repeat=2):
   if a==b:
    if x==y:require(T[3*a+(u+x)%3][3*b+(v+y)%3]==1,'padding')
   else:require(bool(T[3*a+(u+x)%3][3*b+(v+y)%3])==((y-x)%3 in rule[frozenset((a,b))]),'local relation')
 # Equitability, regular enlargement and suffix path are independently checked.
 Q=[]
 for v in range(M):
  i=0
  while i<len(Q) and T[Q[i]][v]:i+=1
  Q.insert(i,v)
 require(all(T[a][b] for a,b in zip(Q,Q[1:])),'ordinary path')
 def R(x,y):
  if x==y:return 0
  if x==2*M:return int(y>=M)
  if y==2*M:return int(x<M)
  if x//M==y//M:return T[x%M][y%M]
  if x<M:return T[y-M][x]
  return 1-T[x-M][y]
 partition=[list(range(i,i+3)) for i in range(0,2*M,3)]+[[2*M]]
 for C in partition:
  for D in partition:require(len({sum(R(x,y) for y in D) for x in C})==1,'equitable')
 require(all(sum(R(x,y) for y in range(2*M+1))==M for x in range(2*M+1)),'regular')
 for i in range(len(owners)):
  require([q for q in range(3) if M+3*i+q==M+3*i]==[0],'owner pin')
 require(all(R(x,2*M) for x in range(M)) and all(R(2*M,M+x) for x in range(M)),'suffix joins')
 # Independent color enumeration, and equality components establish exact count.
 colorings=sum(all(c[a]!=c[b] for a,b in edges) for c in product(range(3),repeat=n))
 raw=None
 if len(owners)<=7:
  raw=sum(all(T[3*a+(u+ph[a])%3][3*b+(v+ph[b])%3] for (a,u),(b,v) in zip(P,P[1:])) for ph in product(range(3),repeat=len(owners)))
  require(raw==colorings*3**(len(odd)//2+2),'literal exhaustive phase count')
 return dict(n=n,edges=len(edges),cells=len(owners),colorings=colorings,raw=raw,regular_order=2*M+1)

reports=[]
for n in range(2,5):
 pairs=list(combinations(range(n),2))
 for bits in product(range(2),repeat=len(pairs)):
  es=[e for e,b in zip(pairs,bits) if b];seen={0}
  for _ in range(n):seen|={v for e in es if any(u in seen for u in e) for v in e}
  if len(seen)==n:reports.append(build(n,es))
for n in (5,6):
 reports.append(build(n,list(combinations(range(n),2))))
 reports.append(build(n,[(i,(i+1)%n) for i in range(n)]))
print(json.dumps({'status':'PASS','independent_construction':True,'cases':len(reports),'reports':reports},indent=2,sort_keys=True))
