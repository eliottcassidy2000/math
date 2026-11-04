import json,math
from pathlib import Path
from fractions import Fraction
P=Path(__file__).parent
A=P.parent/'prime_lanes'
def need(b,m):
 if not b:raise RuntimeError(m)
def allpaths(rows):
 n=len(rows);h=[0]*(1<<n);h[0]=1;e=[[0]*n for _ in h];s=[[0]*n for _ in h]
 def dfs(start,last,mask):
  h[mask]+=1;e[mask][last]+=1;s[mask][start]+=1
  nxt=rows[last]&~mask
  while nxt:
   b=nxt&-nxt;nxt-=b;dfs(start,b.bit_length()-1,mask|b)
 for v in range(n):dfs(v,v,1<<v)
 return h,e,s
# Independently build the eight named suffix vertices, S1,S0,C1,C1.
rows=[0]*8
blocks=[[0],[1],[2,3,4],[5,6,7]]
for i,B in enumerate(blocks):
 for j,C in enumerate(blocks):
  if i<j:
   for b in B:
    for c in C:rows[b]|=1<<c
 if len(B)==3:
  for b,c in zip(B,B[1:]+B[:1]):rows[b]|=1<<c
W={0,2,5};h,e,s=allpaths(rows);N=255
E=[sum(v[w] for w in W) for v in e];S=[sum(v[w] for w in W) for v in s]
z=(h[N],sum(h[x]*h[N^x] for x in range(256)),sum(E[x]*h[N^x] for x in range(256)),sum(h[x]*S[N^x] for x in range(256)),sum(E[x]*S[N^x] for x in range(256)),E[N],S[N])
need(z==(9,576,189,333,100,3,9),'DFS signature')
controls=[]
for t in range(1,5):
 r=[sum(1<<j for j in range(i+1,t+8)) for i in range(t)]+[v<<t for v in rows]+[(1<<t)-1]
 for i in range(8):
  if i in W:r[t+i]|=1<<(t+8)
  else:r[-1]|=1<<(t+i)
 hh,_,_=allpaths(r);need(hh[-1]==189*2**t-91,'full path DFS')
 # Independent Floyd-Warshall reachability.
 reach=[[i==j or bool(r[i]>>j&1) for j in range(len(r))] for i in range(len(r))]
 for k in range(len(r)):
  for i in range(len(r)):
   for j in range(len(r)):reach[i][j]=reach[i][j] or (reach[i][k] and reach[k][j])
 need(all(all(x) for x in reach),'Floyd strong')
 controls.append([t,hh[-1]])
R=json.loads((A/'receipt.json').read_text());roots=[]
for j in range(1,6):
 mod=7**j;period=3*7**(j-1)
 need(pow(2,period,7**(j+1))==(1+7**j)%(7**(j+1)),'exact lift coefficient')
 rr=[t for t in range(period) if (27*pow(2,t,mod)-13)%mod==0]
 need(len(rr)==1,'exhaustive root uniqueness');need(rr[0]==R['seven_adic_roots'][j-1]['root'],'receipt root')
 roots.append([j,period,rr[0]])
densities=[]
for e in range(2,7):
 L0=3*7**(e-2);L1=L0*7
 density=Fraction(1,L0)-Fraction(1,L1)
 need(density==Fraction(2,7**(e-1)),'density')
 densities.append([e,str(density)])
# Deterministic odd trial division, separate from author's all-integer loop.
def factor(n):
 if n%2==0:return 2
 d=3
 while d*d<=n:
  if n%d==0:return d
  d+=2
 return n
prime_clocks=[];composites=[]
for t in range(1,36):
 q=27*2**t-13;f=factor(q)
 if f==q:prime_clocks.append(t)
 else:composites.append([t,f])
need(prime_clocks==[1,4,8,13,19,28,32,35],'prime clocks')
need([x['t'] for x in R['trial_division_prime_witnesses']]==prime_clocks,'receipt clocks')
for w in R['trial_division_prime_witnesses']:need(w['p']==27*2**w['t']-13 and w['H']==7*w['p'] and w['n']==w['t']+9,'witness quantities')
need([x['address'] for x in R['tree'] if x['c']%7==0]==['422'],'unique lane')
out={'status':'PASS','independent_DFS_signature':z,'full_DFS_counts':controls,'exhaustive_residue_roots':roots,'exact_density_fractions':densities,'prime_clocks':prime_clocks,'composite_factor_certificates':composites,'scope':'Proof and independent literal graph/arithmetic audit; no prime infinitude or complete spectrum claim.'}
(P/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
