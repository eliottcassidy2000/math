from itertools import combinations,permutations
from collections import Counter
import json
from pathlib import Path

def need(x,m):
 if not x:raise RuntimeError(m)
def graph(n,mask):
 r=[0]*n
 for k,(u,v) in enumerate(combinations(range(n),2)):
  a,b=(u,v) if mask>>k&1 else(v,u);r[a]|=1<<b
 return r

def cycles(r):
 return sum(all(sum(r[u]>>v&1 for v in t if v!=u)==1 for u in t) for t in combinations(range(len(r)),3))
def partner(r,p):
 n=len(r);s=[((1<<n)-1)^(1<<u)^r[u] for u in range(n)]
 for u,v in zip(p,p[1:]):s[u]^=1<<v;s[v]^=1<<u
 return s

def triangle_palette(bits):
 a=[14,4,8,2];r=[0]*12
 for i in range(12):
  b,u=divmod(i,4);ar=a[u] if not bits[b] else ((15^(1<<u))^a[u]);r[i]=(ar<<(4*b))|(15<<(4*((b+1)%3)))
 return r

def endpoints(r):
 n=len(r);N=1<<n;out=[[0]*n for _ in r]
 for start in range(n):
  dp=[{} for _ in range(N)];dp[1<<start][start]=1
  for mask in range(1,N):
   for u,c in dp[mask].items():
    rest=r[u]&~mask
    while rest:
     b=rest&-rest;rest-=b;v=b.bit_length()-1;target=dp[mask|b];target[v]=target.get(v,0)+c
  for end,c in dp[-1].items():out[start][end]=c
 return out
checks=0;words=0
for n in range(1,6):
 ps=list(permutations(range(n)))
 for mask in range(1<<(n*(n-1)//2)):
  r=graph(n,mask);c=cycles(r);s=[x.bit_count() for x in r]
  for p in ps:
   if not all(r[u]>>v&1 for u,v in zip(p,p[1:])):continue
   words+=1;t=partner(r,p);v=[x.bit_count() for x in t]
   need(v==[n-1-s[i]+(i==p[0])-(i==p[-1]) for i in range(n)],'addressed score transport')
   # n=1 has coincident endpoints, no path edges; formula c3 delta needs n>=2.
   if n>=2:need(cycles(t)-c==s[p[0]]-s[p[-1]]-1,'triangle flux')
   need(partner(t,p)==r,'literal complement involution');checks+=1
rows=[]
for bits,aut in [([0,0,0],81),([1,0,0],27)]:
 r=triangle_palette(bits);mat=endpoints(r);s=[x.bit_count() for x in r];hist=Counter()
 for i in range(12):
  for j in range(12):hist[s[i]-s[j]-1]+=mat[i][j]
 hist={k:v for k,v in sorted(hist.items()) if v}
 need(sum(hist.values())==462915,'H control');need(all(v%aut==0 for v in hist.values()),'free orbit flux divisibility')
 if bits==[0,0,0]:need(all(k%2 for k in hist),'every uniform complement exits C3 count')
 p=[]
 for u in range(12):
  j=next((j for j,v in enumerate(p) if r[u]>>v&1),len(p));p.insert(j,u)
 rows.append({'bits':bits,'aut':aut,'H':sum(hist.values()),'endpoint_matrix':mat,'labelled_flux':hist,'ordinary_tiling_flux':{k:v//aut for k,v in hist.items()},'example_path':p,'example_partner_c3':cycles(partner(r,p))})
print(json.dumps({'status':'PASS','labelled_universe_orders':[1,5],'hamiltonian_words':words,'checks':checks,'controls':rows},sort_keys=True,indent=2))
