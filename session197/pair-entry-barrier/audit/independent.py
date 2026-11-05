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
checks=0;excesses=set();physical_roots=0
for n in range(1,6):
 for T in graphs(n):
  h,e,s,R,K,C=inv(T)
  for a in range(n):
   need(R[a]>=h and 0<=e[a]<=h,'physical cone inputs');physical_roots+=1
  if C!=1:continue
  for a in range(n):
   if T[a].bit_count()!=1:continue
   b=T[a].bit_length()-1;U=join(T,{a,b});H,ee,ss,RR,KK,CC=inv(U)
   need(R[b]==2*(h-e[a])-e[b],'literal adjacency identity')
   ex=RR[n]-2*(H-ee[n]);need(ex==R[a]-h+e[b]>=1,'literal first clone excess')
   need(RR[a]<=2*(H-ee[a]),'literal last clone cone');excesses.add(ex);checks+=1
T=(18,20,25,3,8);need(profile(T,0)==(13,2,18),'literal marked inside-cone example')
small=[]
for n in (3,4):
 vals=sorted(set(inv(T)[0]for T in graphs(n)if inv(T)[5]));need(vals==([3]if n==3 else[5]),'literal small strong values');small.append(vals)
print(json.dumps(dict(status='PASS',method='Independent literal permutation and cut inventories; no source DP imports',physical_root_controls=physical_roots,adjacent_clone_pairs=checks,positive_excesses=sorted(excesses),strong_values=small),sort_keys=True,indent=2))
