#!/usr/bin/env python3
import pathlib,json,itertools as it,math
D=pathlib.Path(__file__).resolve().parent
R=json.loads((D/'result.json').read_text())
def ck(p,m):
 if not p:raise RuntimeError(m)
def hp(n,E,deleted):
 inc=[0]*n
 for j,(u,v) in enumerate(E):
  if not(deleted>>j&1):inc[v]|=1<<u
 end=[0]*(1<<n)
 for v in range(n):end[1<<v]=1<<v
 for S in range(1,1<<n):
  if S&(S-1)==0:continue
  bits=S
  while bits:
   bit=bits&-bits;bits-=bit;v=bit.bit_length()-1
   if end[S^bit]&inc[v]:end[S]|=bit
 return bool(end[-1])
out=[]
for row in R['censuses']:
 n=row['n'];E=[tuple(x) for x in row['edges']];V=set(range(n));r=(n-1)//2;lam=(n-3)//4
 IN={v:{u for u,w in E if w==v} for v in V}
 ck(all(len(IN[v])==r for v in V),'regularity')
 ck(all(len(IN[u]&IN[v])==lam for u,v in it.combinations(V,2)),'double regularity')
 equality=[];expansion=[]
 for t in range(n-1):
  counts=[]
  for T in it.combinations(sorted(V),t):
   U=set().union(*(IN[v] for v in T));counts.append(len(U))
   if t:ck(len(U)*(lam*t+r-lam)>=r*r*t,'Cauchy')
   if len(U)<=t+2:
    for extra in it.combinations(sorted(V-U),t+2-len(U)):
     S=U|set(extra);mask=sum(1<<i for i,(u,v) in enumerate(E) if u in S and v not in T)
     equality.append((mask,tuple(sorted(S)),T))
  expansion.append({'T_size':t,'minimum_union_size':min(counts)})
 ck(len(set(m for m,_,_ in equality))==row['distinct_cost_n_minus_one_Hall_blockers'],'independent Hall count')
 ck({m for m,_,_ in equality}=={x['mask'] for x in row['blockers']},'independent Hall mask equality')
 for b in row['blockers']:
  for cert in b['certificates']:
   S=set(cert['S']);T=set(cert['T']);ck(len(S)-len(T)==2,'deficiency')
   ck(b['mask']==sum(1<<i for i,(u,v) in enumerate(E) if u in S and v not in T),'certificate mask')
 stars=[sum(1<<i for i,(u,v) in enumerate(E) if u==c or v==c) for c in V]
 candidates=sorted(set(stars)|{m for m,_,_ in equality})
 ck(all(not hp(n,E,b) for b in candidates),'HP blocking')
 failures=[]
 for b in candidates:
  for i in range(len(E)):
   if b>>i&1 and not hp(n,E,b^(1<<i)):failures.append((b,i))
 ck(not failures,'inclusion minimality')
 out.append({'n':n,'expansion':expansion,'verified_HP_blockers':len(candidates),'all_verified_inclusion_minimal':True,'scope':'does not certify global minimum or classify exotic blockers'})
print(json.dumps({'status':'PASS','checks':out},indent=2,sort_keys=True))
