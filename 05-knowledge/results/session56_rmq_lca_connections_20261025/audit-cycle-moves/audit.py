#!/usr/bin/env python3
from itertools import combinations,permutations
from collections import defaultdict,Counter,deque
from pathlib import Path
import json

def require(x,msg):
 if not x: raise RuntimeError(msg)

def induced(n,E):
 S=set(E); ans=[]
 for r in range(3,n+1):
  for V in combinations(range(n),r):
   Q=[e for e in E if e[0] in V and e[1] in V]
   if len(Q)!=r: continue
   nb={v:[] for v in V}
   for u,v in Q:nb[u].append(v);nb[v].append(u)
   if any(len(nb[v])!=2 for v in V):continue
   walk=[V[0],min(nb[V[0]])]
   while len(walk)<r:
    z=next(w for w in nb[walk[-1]] if w!=walk[-2])
    if z in walk:break
    walk.append(z)
   if len(walk)==r and walk[0] in nb[walk[-1]]:ans.append(walk)
 return ans

def data(n,E):
 Cs=induced(n,E);idx={e:i for i,e in enumerate(E)};T=[]
 for C in Cs:
  mask=pat=0
  for u,v in zip(C,C[1:]+C[:1]):
   b=1<<idx[min(u,v),max(u,v)];mask|=b
   if u<v:pat|=b
  T.append((mask,pat))
 fibres=defaultdict(list)
 for m in range(1<<len(E)):
  s=[0]*n
  for i,(u,v) in enumerate(E):s[u if m&(1<<i) else v]+=1
  fibres[tuple(s)].append(m)
 visited=set();components=[];move_edges=0;maxdistance=0
 for m in range(1<<len(E)):
  if m in visited:continue
  q=deque([(m,0)]);visited.add(m);component=[]
  while q:
   x,d=q.popleft();component.append(x);maxdistance=max(maxdistance,d)
   for mask,pat in T:
    if x&mask in (pat,mask^pat):
     y=x^mask
     if x<y:move_edges+=1
     if y not in visited:visited.add(y);q.append((y,d+1))
  components.append(tuple(sorted(component)))
 require(sorted(components)==sorted(tuple(v) for v in fibres.values()),'components differ from fibres')
 return Cs,fibres,move_edges,maxdistance

report={'method':'Independent induced-subgraph 2-regular cycle detection and BFS; explicit RuntimeError checks, unaffected by -O','path_census':[],'generic':[]}
producer=json.loads((Path(__file__).resolve().parent.parent/'cycle-moves/certificate.json').read_text())
for n in range(1,8):
 E=[(u,v) for u,v in combinations(range(n),2) if v-u!=1]
 Cs,F,M,D=data(n,E)
 row={'n':n,'orientations':1<<len(E),'fibres':len(F),'induced_cycle_counts':dict(Counter(map(len,Cs))),'undirected_move_edges':M,'fibre_size_histogram':dict(Counter(map(len,F.values()))),'maximum_BFS_distance_from_representative':D}
 p=producer['census'][n-1]
 for key in ['orientations','fibres','undirected_move_edges']:require(row[key]==p[key],f'producer census {n} {key}')
 for key in ['induced_cycle_counts','fibre_size_histogram']:require({str(k):v for k,v in row[key].items()}==p[key],f'producer census {n} {key}')
 require(all(len(C)<=4 for C in Cs),'path induced length')
 report['path_census'].append(row)
 print('path',n,'PASS',flush=True)
for n in range(1,6):
 pairs=list(combinations(range(n),2));states=witnesses=0
 for mask in range(1<<len(pairs)):
  E=[e for i,e in enumerate(pairs) if mask&(1<<i)]
  Cs,F,M,D=data(n,E);states+=1<<len(E)
  if not Cs:require(all(len(f)==1 for f in F.values()),'forest nonsingleton');continue
  C=max(Cs,key=len);V=set(C);arcs=set(zip(C,C[1:]+C[:1]));m=0
  for i,(u,v) in enumerate(E):
   yes=(u,v) in arcs if u in V and v in V else (u not in V if (u in V)!=(v in V) else u<v)
   if yes:m|=1<<i
  f=next(f for f in F.values() if m in f)
  require(len(f)==2,'sharp fibre not two');require((f[0]^f[1]).bit_count()==len(C),'sharp fibre distance')
  witnesses+=1
 report['generic'].append({'n':n,'graphs':1<<len(pairs),'states':states,'sharpness_witnesses':witnesses})
 print('generic',n,'PASS',flush=True)
# Independently enumerate Hamilton paths in sharp pair.
n=5;E=[(u,v) for u,v in combinations(range(n),2) if v-u!=1];_,F,_,_=data(n,E)
s=tuple(producer['sharpness']['residual_score']);sharp=[]
for m in F[s]:
 arcs={(i,i+1) for i in range(n-1)}|{(u,v) if m&(1<<i) else (v,u) for i,(u,v) in enumerate(E)}
 H=sum(all(e in arcs for e in zip(p,p[1:])) for p in permutations(range(n)))
 require(not any((u,v) in arcs and (v,w) in arcs and (w,u) in arcs and all(tuple(sorted(e)) in E for e in [(u,v),(v,w),(w,u)]) for u,v,w in permutations(range(n),3)),'sharp triangle')
 sharp.append({'mask':m,'Hamilton_paths':H})
require(sorted(x['Hamilton_paths'] for x in sharp)==[11,13],'HP counts')
report['sharpness']=sharp;report['status']='PASS'
Path(__file__).with_name('audit-result'+('-O' if not __debug__ else '')+'.json').write_text(json.dumps(report,indent=2)+'\n')
print('ALL PASS')
