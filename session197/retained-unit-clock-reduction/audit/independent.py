from itertools import permutations,combinations
from functools import lru_cache
import json

def need(x,s):
 if not x:raise RuntimeError(s)
def graphs(n):
 E=list(combinations(range(n),2))
 for z in range(1<<len(E)):
  r=[0]*n
  for k,(a,b) in enumerate(E):r[a if z>>k&1 else b]|=1<<(b if z>>k&1 else a)
  yield tuple(r)
def path(r,p):return all(r[a]>>b&1 for a,b in zip(p,p[1:]))
@lru_cache(None)
def inv(r):
 n=len(r);h=0;e=[0]*n;s=[0]*n;R=[0]*n;K=[[0]*n for _ in r];C=0
 for p in permutations(range(n)):
  if path(r,p):
   h+=1;s[p[0]]+=1;e[p[-1]]+=1
   if n>=3 and p[0]==0 and r[p[-1]]&1:C+=1
  for cut in range(n):
   if path(r,p[:cut]) and path(r,p[cut:]):
    R[p[cut]]+=1
    if cut:K[p[cut-1]][p[cut]]+=1
 return h,e,s,R,K,C
def gadget(L):
 w=[] if L==1 else [1-int(b)for b in bin(L-1)[2:][::-1]];G=[0]*(len(w)+1)
 for i,b in enumerate(w,1):
  G[0 if b else i]|=1<<(i if b else 0)
  for j in range(i+1,len(G)):G[i]|=1<<j
 return tuple(G)
def join(T,entries,G=(0,)):
 n=len(T);r=list(T)+[x<<n for x in G]
 for a in range(n):
  for b in range(len(G)):r[n if b==0 and a in entries else a]|=1<<(a if b==0 and a in entries else n+b)
 return tuple(r)
def profile(r,a):
 h,e,s,R,K,C=inv(r);return h,e[a],R[a]
def permute(r,f):
 out=[0]*len(r)
 for a in range(len(r)):
  for b in range(len(r)):
   if r[a]>>b&1:out[f[a]]|=1<<f[b]
 return tuple(out)
def deleted(T,r):
 V=[x for x in range(len(T))if x!=r]
 return tuple(sum(1<<j for j,y in enumerate(V)if T[x]>>y&1)for x in V)
def dcount(T,r):
 V=[x for x in range(len(T))if x!=r]
 return sum(path(T,p)for p in permutations(V))
def mod(U,r,S):
 n=len(U);T=list(U)+[v<<n for v in S]
 for a in range(n):
  for b in range(len(S)):
   if a==r:T[n+b]|=1<<a
   else:T[a]|=1<<(n+b)
 return tuple(T)
def strong_on(T,V):
 if len(V)==1:return True
 for a in V:
  reach={a};todo=[a]
  for u in todo:
   for v in V:
    if v not in reach and T[u]>>v&1:reach.add(v);todo.append(v)
  if len(reach)!=len(V):return False
 return True
def sink(T,r):
 V=[v for v in range(len(T))if v!=r];out=[]
 for mask in range(1,1<<len(V)):
  S=[v for i,v in enumerate(V)if mask>>i&1]
  if all(T[b]>>a&1 for a in S for b in V if b not in S)and strong_on(T,S):out.append(S)
 need(len(out)==1,'literal terminal SCC uniqueness');return out[0]
mods=states=units=graphmaps=transients=0
for n in range(1,4):
 for U in graphs(n):
  h,e,starts,R,K,C=inv(U)
  for r in range(n):
   d=dcount(U,r)
   for m in range(1,4):
    for S in graphs(m):
     G=inv(S)[0];B=sum(path(S,p[:k])and path(S,p[k:])for p in permutations(range(m))for k in range(m+1))
     T=mod(U,r,S);H,ee,ss,RR,KK,CC=inv(T)
     need((H,ee[r],RR[r],dcount(T,r))==((R[r]-d)*B+(h-e[r]-R[r]+2*d)*G,d*G,(R[r]-d)*B+d*G,d*G),'literal module identities');mods+=1
for n in range(1,6):
 for T in graphs(n):
  h,e,starts,R,K,C=inv(T)
  for r in range(n):
   d=dcount(T,r);Z=mod(T,r,(0,));h1=h-e[r]+R[r]
   beta0=R[r]-h-e[r];beta1=R[r]-h+e[r]-2*d
   need((beta0==beta1)==(e[r]==d),'exact first-step boundary');transients+=1
   if n==1:need(h==1 and e[r]==d==1 and beta0==-1,'singleton convention');continue
   if e[r]!=d:continue
   S=sink(T,r);need(all(T[v]>>r&1 for v in S),'terminal incidence');G=sum(path(T,p)for p in permutations(S));beta=R[r]-h-d
   need(beta%G==0,'offset divisibility');states+=1
   if abs(beta)==1:
    need(len(S)==1 and T[S[0]]==1<<r,'unit terminal singleton');units+=1
    v=S[0];U=T
    for t in range(1,4):
     U=mod(U,r,(0,))
     # Exact clone isomorphism: retained old v followed by each new clone.
     cl=[v]+list(range(n,n+t));outside=[x for x in range(n)if x!=v]
     need(all(U[a]>>b&1 for i,a in enumerate(cl)for b in cl[i+1:]),'transitive clone chain')
     need(all(((U[a]>>x)&1)==((T[v]>>x)&1)and((U[x]>>a)&1)==((T[x]>>v)&1)for a in cl for x in outside),'literal exterior clone map');graphmaps+=1
     if n<=4:need(sum(path(U,p)for p in permutations(range(len(U))))==(h+beta)*2**t-beta,'literal returning H')
print(json.dumps(dict(status='PASS',literal_module_constructions=mods,all_root_initial_boundary_checks=transients,physical_e_equals_d_controls=states,unit_controls=units,exact_clone_graph_maps=graphmaps,method='Literal permutations/cuts, subset-defined sink SCC and exact adjacency; no source DP imports'),sort_keys=True,indent=2))
