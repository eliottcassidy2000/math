#!/usr/bin/env python3
"""Independent ordinary orbit/physical-path audit and closure controls.
No producer module imports. Source universe n3..6; closure data n3..7.
"""
from collections import Counter,defaultdict,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations
from math import comb,factorial
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent;PACKET=HERE.parent
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)

def refine(cells,matrices):
 n=len(matrices[0]);levels=[cells]
 while True:
  fresh=[]
  for cell in cells:
   groups=defaultdict(list)
   for i in cell:
    sig=tuple(sum(M[i][j] for j in target) for M in matrices for target in cells)
    groups[sig].append(i)
   fresh.extend(sorted(groups.values(),key=lambda z:z[0]))
  fresh.sort(key=lambda z:z[0])
  if fresh==cells:return levels
  need(len(fresh)>len(cells),'strict refinement')
  cells=fresh;levels.append(cells)

def initial(keys):
 d=defaultdict(list)
 for i,key in enumerate(keys):d[key].append(i)
 return sorted(d.values(),key=lambda z:z[0])

def literal(n,source):
 pairs=list(combinations(range(n),2));ix={p:k for k,p in enumerate(pairs)}
 top=(1<<len(pairs))-1;ps=list(permutations(range(n)))
 def arc(t,a,b):return bool((t>>ix[min(a,b),max(a,b)])&1)==(a<b)
 def encode(t,p):return sum(int(arc(t,p[a],p[b]))<<j for j,(a,b) in enumerate(pairs))
 unseen=set(range(top+1));owner={};autos={}
 while unseen:
  t=min(unseen);orbit={encode(t,p) for p in ps}
  need(orbit<=unseen,'overlapping literal orbits')
  for u in orbit:owner[u]=t
  unseen-=orbit;autos[t]=factorial(n)//len(orbit)
 reps=sorted(autos);dual={t:owner[top^t] for t in reps}
 C={};B={};covers={}
 for t in reps:
  cc=Counter();bb=Counter();ff=[0]*(n+1)
  for p in ps:
   bad=sum(not arc(t,p[j],p[j+1]) for j in range(n-1))
   for k in range(bad+1,n+1):ff[k]+=comb(n-1-bad,k-1-bad)
   if bad:continue
   pathmask=sum(1<<ix[min(p[j],p[j+1]),max(p[j],p[j+1])] for j in range(n-1))
   u=owner[top^t^pathmask];cc[u]+=1
   reversal=[0]*n
   for a,b in zip(p,p[::-1]):reversal[a]=b
   if encode(t,reversal)==top^t:bb[u]+=1
  for u,v in cc.items():
   need(v%autos[t]==0,'free path orbit divisor');C[t,u]=v//autos[t]
  for u,v in bb.items():
   need(v%autos[t]==0,'blue free path orbit divisor');B[t,u]=v//autos[t]
  covers[t]=ff[1:]
 for t in reps:
  for u in reps:
   need(C.get((t,u),0)==C.get((u,t),0),'ordinary C symmetry')
   need(B.get((t,u),0)==B.get((u,t),0),'ordinary blue symmetry')
   need(C.get((t,u),0)==C.get((dual[t],dual[u]),0),'converse equivariance')
  need(C.get((t,t),0)%2==0,'ordinary loop doubling')
  if dual[t]!=t:need(not any(k[0]==t for k in B),'blue source must be self-converse')
 def triples(M):return [[a,b,v] for (a,b),v in sorted(M.items()) if v]
 need(triples(C)==source['ordinary_A'],'independent full ordinary adjacency')
 need(triples(B)==source['ordinary_blue'],'independent full ordinary blue')
 src={x['code']:x for x in source['ordinary_nodes']}
 need(reps==sorted(src),'complete ordinary representative list')
 for t in reps:
  need(autos[t]==src[t]['aut'] and dual[t]==src[t]['converse'],'ordinary class metadata')
  need(covers[t]==src[t]['ordered_covers'],'literal ordered covers')
 A=Counter();AB=Counter()
 for (t,u),v in C.items():A[min(t,dual[t]),min(u,dual[u])]+=v
 for (t,u),v in B.items():AB[min(t,dual[t]),min(u,dual[u])]+=v
 need(triples(A)==source['merged_A'] and triples(AB)==source['merged_blue'],'merged projection')
 chirals=[t for t in reps if t<dual[t]];D={(u,v):C.get((u,v),0)-C.get((u,dual[v]),0) for u in chirals for v in chirals}
 seen=set();components=[]
 for u in chirals:
  if u in seen:continue
  todo=[u];seen.add(u);comp=[]
  while todo:
   v=todo.pop();comp.append(v)
   for w in chirals:
    if w not in seen and D[v,w]:seen.add(w);todo.append(w)
  components.append(sorted(comp))
 for u in chirals:
  for v in chirals:
   need(2*(C.get((u,v),0)+C.get((u,dual[v]),0))==A[u,v],'even block sum')
   need(F(A[u,v],4)+F(D[u,v],2)==C.get((u,v),0),'signed block reconstruction')
 return {'n':n,'ordinary':len(reps),'self_converse':sum(t==dual[t] for t in reps),
         'chiral_signed_components':components,'C_sha256':sha256(json.dumps(triples(C)).encode()).hexdigest()}

def main():
 data=json.loads((PACKET/'closure/data.json').read_text());out=[]
 for row in data['rows']:
  nodes=row['merged_nodes'];codes=[x['code'] for x in nodes];idx={c:i for i,c in enumerate(codes)};N=len(nodes)
  A=[[F(0)]*N for _ in nodes];B=[[F(0)]*N for _ in nodes]
  for a,b,v in row['merged_A']:A[idx[a]][idx[b]]=v
  for a,b,v in row['merged_blue']:B[idx[a]][idx[b]]=v
  d=[sum(r) for r in A];P=[[F(x,d[i]) for x in r] for i,r in enumerate(A)];PB=[[F(x,d[i]) for x in r] for i,r in enumerate(B)]
  need(sum(d)==1<<((row['n']-1)*(row['n']-2)//2),'complete fixed-spine mass')
  cells=initial([tuple(x['ordered_covers']) for x in nodes]);results={}
  for label,Ms in [('A',[A]),('P',[P]),('A_colored',[B,[[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]]),('P_colored',[PB,[[a-b for a,b in zip(r,s)] for r,s in zip(P,PB)]])]:
   lev=refine(cells,Ms);counts=[len(c) for c in lev]
   need(counts==row['closure']['Q'][label]['counts'][:-1],('closure agreement',row['n'],label))
   need(counts[-1]==N,'discrete Q closure')
   results[label]=counts
  wl=refine(initial(d),[A]);need(row['n']==3 or len(wl[-1])==N,'degree refinement rigidity')
  result={'n':row['n'],'refinement':results,'weighted_degree_WL':[len(x) for x in wl]}
  if row['n']<=6:result['literal']=literal(row['n'],row)
  if row['n']==5:
   one=initial([(tuple(x['ordered_covers']),x['blue_degree']==1) for x in nodes])
   need(len(one)==N and len(cells)<N,'one-bit minimal closure')
  out.append(result)
 # Arbitrarily long repair depth despite reversibility and two-step detection.
 paths=[]
 for m in range(2,15):
  A=[[int(abs(i-j)==1) for j in range(m+1)] for i in range(m+1)];d=list(map(sum,A))
  P=[[F(x,d[i]) for x in r] for i,r in enumerate(A)]
  levels=refine([[0],list(range(1,m+1))],[P]);need(len(levels)-1==m-1,'unbounded repair depth')
  # From0, actualtwo-stepreturn=1/2; one-block stationaryreturn=1/(2m-1).
  need(sum(P[0][j]*P[j][0] for j in range(m+1))==F(1,2),'path return')
  need(F(1,2)-F(1,2*m-1)>0,'initial second-order defect')
  paths.append({'edges':m,'strict_rounds':len(levels)-1,'counts':[len(v) for v in levels]})
 return {'status':'PASS: independent finite audit; no all-order tournament extrapolation','checks':checks,'rows':out,'reversible_path_controls':paths,'data_sha256':sha256((PACKET/'closure/data.json').read_bytes()).hexdigest(),'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
