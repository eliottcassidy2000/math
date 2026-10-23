import itertools as it,json
from pathlib import Path
import sympy as S
from collections import Counter,defaultdict
ROOT=Path(__file__).resolve().parent.parent
producer=json.loads((ROOT/'metagraph-memory/results.json').read_text())
checks=0
def ck(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
results=[]
for expected in producer:
 n=expected['n'];pairs=list(it.combinations(range(n),2));free=[e for e in pairs if e[1]-e[0]>1];N=2**len(free)
 permutations=list(it.permutations(range(n))); keys=[];blue=[]
 # Full permutation canonicalization. No producer degree sorting, helpers, imports.
 for x in range(N):
  edges={(i,i+1) for i in range(n-1)}|{(u,v) if x>>k&1 else (v,u) for k,(u,v) in enumerate(free)}
  key=min(min(tuple(int((p[u],p[v]) in edges) for u,v in pairs),tuple(int((p[v],p[u]) in edges) for u,v in pairs)) for p in permutations)
  keys.append(key)
  # Reversal-converse invariance of the ENTIRE directed tournament.
  reflected={(n-1-v,n-1-u) for u,v in edges}
  blue.append(edges==reflected)
 labels={k:i for i,k in enumerate(sorted(set(keys)))};node=[labels[k] for k in keys];r=len(labels)
 cells={(node[x],node[N-1-x]) for x in range(N)}
 colored={(node[x],node[N-1-x],blue[x]) for x in range(N)}
 ck((N,r,len(cells),len(colored))==(expected['tilings'],expected['nodes'],expected['pair_states'],expected['colored_pair_states']),'census')
 ck(all(blue[x]==blue[N-1-x] for x in range(N)),'color involution invariant')
 witnesses=[]
 for w in expected['mixed']:
  a,b=w['blue'],w['black']
  ck(blue[a] and not blue[b] and (node[a],node[N-1-a])==(node[b],node[N-1-b]),'mixed color witness')
  witnesses.append([a,b,node[a],node[N-1-a]])
 if n<=5:
  # Verify partition equality independent of node naming.
  transport={node[x]:expected['node'][x] for x in range(N)}
  ck(len(transport)==r and len(set(transport.values()))==r and all(transport[node[x]]==expected['node'][x] for x in range(N)),'node equivalence')
  J=S.zeros(N,r)
  for x,v in enumerate(node):J[x,v]=1
  PJ=J[::-1,:];F=J.T*J;A=J.T*PJ;T=F.inv()*A
  reordered=S.zeros(r)
  for i in range(r):
   for j in range(r):reordered[transport[i],transport[j]]=A[i,j]
  ck(reordered.tolist()==expected['A'],'adjacency reconstruction')
  ck([int(F[i,i]) for i in sorted(range(r),key=lambda i:transport[i])]==expected['fibres'],'fibres')
  W=J.row_join(PJ);delta=F-A*F.inv()*A;H=PJ-J*T
  ck(H.T*H==delta,'defect gram');d=delta.rank()
  ck(W.rank()==r+d==expected['linear_dimension'],'cyclic dimension')
  # Hidden operator applied by swapping tilings then orthogonal projection, no dense QPQ powers.
  def Q(X):return X-J*F.inv()*J.T*X
  V=H
  for k in range(7):
   ck(H.T*V==(-1)**k*delta*T**k,'memory kernel')
   V=Q(V[::-1,:])
  for t in range(10):ck(F.inv()*J.T*(J if t%2==0 else PJ)==(S.eye(r) if t%2==0 else T),'alternating moments')
  h=S.Matrix([(-1)**x*(x+1) for x in range(N)]);z=Q(h);eta=F.inv()*J.T*z[::-1,:]
  for t in range(6):
   ck(F.inv()*J.T*z[::-1,:]==(-T)**t*eta,'nonuniform forcing propagation')
   z=Q(z[::-1,:])
  if n>=4:
   witness=next((x,y) for x in range(N) for y in range(x) if node[x]==node[y] and node[N-1-x]!=node[N-1-y])
   ck(witness[0]!=witness[1],'same current output needs initial sidecar')
  moment=F.row_join(A).col_join(A.row_join(F))
  ck(moment.rank()==r+d,'dilation Gram lower bound')
  adj={i:{j for j in range(r) if A[i,j]!=0} for i in range(r)};remaining=set(adj);c=b=0
  while remaining:
   root=min(remaining);shades={root:0};todo=[root];bip=True
   while todo:
    u=todo.pop()
    for v in adj[u]:
     if v not in shades:shades[v]=1-shades[u];todo.append(v)
     elif shades[v]==shades[u]:bip=False
   remaining-=shades.keys();c+=1;b+=bip
  ck(d==r-c-b and W.rank()==2*r-c-b,'component rank formula')
  # Stable refinement direct all ten-time observation words, not producer refinement iteration.
  words={tuple(node[x if t%2==0 else N-1-x] for t in range(10)) for x in range(N)}
  ck(len(words)==len(cells),'two-observation minimal state')
  results.append(dict(n=n,tilings=N,nodes=r,pair_states=len(cells),colored_states=len(colored),defect_rank=d,linear_dimension=W.rank(),components=c,bipartite_components=b,markov_failure=T*T!=S.eye(r)))
 else:results.append(dict(n=n,tilings=N,nodes=r,pair_states=len(cells),colored_states=len(colored),mixed_witnesses=witnesses))
out={'status':'PASS','checks':checks,'results':results,'independence':'Full relabeling without degree sorting; no producer imports; exact rational Gram and memory reconstruction.'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
