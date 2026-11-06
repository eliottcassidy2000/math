import itertools,json,random
import numpy as np

def graph(n):
 pairs=list(itertools.combinations(range(n),2)); ix={e:i for i,e in enumerate(pairs)}; G={t:set() for t in range(1<<len(pairs))}
 def arc(t,a,b):return bool((t>>ix[min(a,b),max(a,b)])&1)==(a<b)
 for t in G:
  for p in itertools.permutations(range(n)):
   if all(arc(t,a,b) for a,b in zip(p,p[1:])):
    keep=sum(1<<ix[min(a,b),max(a,b)] for a,b in zip(p,p[1:]));u=t^(((1<<len(pairs))-1)^keep)
    if t!=u:G[t].add(u)
 return G

def component(G,s,excluded=None):
 seen={s};queue=[s]
 for a in queue:
  for b in sorted(G[a]):
   if tuple(sorted((a,b)))==excluded:continue
   if b not in seen: seen.add(b);queue.append(b)
 return seen
out={};rng=random.Random(208)
for n in (3,4):
 G=graph(n);comps=[];unseen=set(G)
 while unseen:
  C=component(G,min(unseen));unseen-=C;comps.append(C)
 out[n]={'vertices':len(G),'edges':sum(map(len,G.values()))//2,'components':list(map(len,comps)),'leaves':sum(len(e)==1 for e in G.values())}
 for C in comps:
  nodes=sorted(C);E=[(a,b) for a in nodes for b in G[a] if a<b];T={a:set() for a in nodes};seen={nodes[0]};q=[nodes[0]]
  for a in q:
   for b in sorted(G[a]):
    if b not in seen:seen.add(b);q.append(b);T[a].add(b);T[b].add(a)
  caps={}
  for a in nodes:
   for b in T[a]:
    if a<b:
     S=component(T,a,(a,b));caps[a,b]=sum((u in S)!=(v in S) for u,v in E)
  A=np.zeros((len(nodes),len(nodes)),dtype=int); idx={a:i for i,a in enumerate(nodes)}
  for (a,b),w in list(caps.items())+[((a,b),-1) for a,b in E]:
   i,j=idx[a],idx[b];A[i,i]+=w;A[j,j]+=w;A[i,j]-=w;A[j,i]-=w
  vals,vecs=np.linalg.eigh(A)
  for scale in (1,2,4,10,100,1000):
   z=np.rint(scale*vecs[:,0]).astype(int);x={a:int(z[idx[a]]) for a in nodes};eg=sum((x[a]-x[b])**2 for a,b in E);et=sum(c*(x[a]-x[b])**2 for (a,b),c in caps.items())
   if et<eg:
    out[n]['energy_hostile']={'component':nodes,'tree':[[a,b,c] for (a,b),c in sorted(caps.items())],'potential':[x[a] for a in nodes],'graph_energy':eg,'tree_energy':et,'graph_TV':sum(abs(x[a]-x[b]) for a,b in E),'tree_TV':sum(c*abs(x[a]-x[b]) for (a,b),c in caps.items()),'cut_domination':'proved by path-crossing argument; not exhaustively enumerated'};break
  if 'energy_hostile' in out[n]:break
print(json.dumps(out,indent=2))
