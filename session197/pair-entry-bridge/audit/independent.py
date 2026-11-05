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
counts=dict(pair=0,adjacent=0,nonadjacent=0,birth_rows=0)
for n in range(2,4):
 for T in graphs(n):
  h,e,s,R,K,C=inv(T)
  for a,b in combinations(range(n),2):
   for L in range(1,9):
    G=gadget(L);g,eg,sg,Rg,Kg,cg=inv(G);delta=int(len(G)==1)
    V=join(T,{a,b},G);want=(g*h-sg[0]*(e[a]+e[b])+L*(R[a]+R[b])-K[a][b]-K[b][a],h-delta*(e[a]+e[b]),Rg[0]*h+2**(len(G)-1)*(R[a]+R[b]))
    need(profile(V,n)==want,'pair literal profile');counts['pair']+=1
for n in range(3,5):
 for T in graphs(n):
  h,e,s,R,K,C=inv(T)
  if C!=1:continue
  for a in range(n):
   if T[a].bit_count()!=1:continue
   b=T[a].bit_length()-1
   U=join(T,{a,b});hs=h-e[b]+R[b]
   need(profile(U,n)==(hs,h-e[a]-e[b],h+R[a]+R[b]),'first clone literal')
   need(profile(U,a)==(hs,h-e[b],h+R[b]),'last clone literal')
   need(inv(U)[2][n]==2 and inv(U)[2][a]==1,'clone starts');counts['adjacent']+=1
   F=join(U,{n});V=join(join(T,{a}),{n,b});f=list(range(n+2));f[a]=n;f[n]=n+1;f[n+1]=a
   need(permute(V,f)==F,'full marked nonadjacent map')
   need(profile(F,a)==profile(V,n+1) and inv(F)[2][a]==2,'nonadjacent exact root')
   need(profile(F,a)[0]==h+e[a]+R[a]+2*R[b],'nonadjacent H');counts['nonadjacent']+=1
  for a in range(n):
   old=[u for u in range(n)if u!=a];small=tuple(sum(1<<j for j,v in enumerate(old)if T[u]>>v&1)for u in old);d=inv(small)[0]
   Z=join(T,{a});O=join(T,set())
   need(profile(Z,n)==(h-e[a]+R[a],h-e[a],h+R[a]),'J')
   need(profile(Z,a)==(h-e[a]+R[a],d,2*R[a]-d),'Z')
   need(profile(O,a)==(h,0,2*R[a]),'O');counts['birth_rows']+=3
print(json.dumps(dict(status='PASS',checks=counts,method='Literal permutation/cut inventories and independent adjacency construction; no subset DP or imported source counter.'),sort_keys=True,indent=2))
