#!/usr/bin/env python3
import json
from pathlib import Path
P=Path(__file__).resolve().parent

def need(q,why):
 if not q:raise RuntimeError(why)

def join(a,b):
 n=len(a);m=len(b)
 return [x|(((1<<m)-1)<<n)for x in a]+[x<<n for x in b]

def profile(a):
 n=len(a);N=1<<n;F=N-1
 for i,x in enumerate(a):
  need(0<=x<N and not(x>>i&1),'row')
  for j in range(i):need(((x>>j&1)+(a[j]>>i&1))==1,'orientation')
 dp=[[0]*n for _ in range(N)];root=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):dp[1<<v][v]=1
 root[1][0]=1
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   for w in range(n):
    if not S>>w&1 and a[v]>>w&1:
     dp[S|1<<w][w]+=dp[S][v];root[S|1<<w][w]+=root[S][v]
  H[S]=sum(dp[S])
 s=sum(root[F]);L=sum(dp[S][0]*H[F^S]for S in range(N))
 return H[F]+L-s,s

strong={}
for line in (P/'STRONG_PROFILES.tsv').read_text().splitlines():
 B,s,n,*a=map(int,line.split());need(len(a)==n,'source order');strong[B,s]=a
factors={}
for line in (P/'BANK.tsv').read_text().splitlines():
 n,h,code,*a=map(int,line.split())
 if h<=26:factors.setdefault(h,a)
tails={1:[]}
for q in range(2,27):
 for f,a in sorted(factors.items()):
  if q%f==0 and q//f in tails:tails[q]=join(a,tails[q//f]);break
witness={}
for q,t in tails.items():
 candidates=[(q,q,join([0],t))]
 candidates.extend((q*B,q*s,join(a,t))for(B,s),a in strong.items()if q*B<=137 and q*s<=26)
 for B,s,a in candidates:
  if (B,s)not in witness or len(a)<len(witness[B,s]):witness[B,s]=a
expected={tuple(map(int,line.split()))for line in(P/'ALL_PROFILES.tsv').read_text().splitlines()}
need(set(witness)==expected,'complete pair list match')
for p,a in witness.items():need(profile(a)==p,('literal positive witness',p,a))
(P/'ALL_WITNESSES.tsv').write_text(''.join(' '.join(map(str,(B,s,len(a),*a)))+'\n'for(B,s),a in sorted(witness.items())))
report={'status':'PASS','all_literal_profile_witnesses':len(witness),'maximum_witness_order':max(map(len,witness.values())),'tail_counts':sorted(tails),'base137_witnesses':[{'r':r,'j':j,'B':137-2*j,'s':r+3*j,'order':len(witness[137-2*j,r+3*j]),'root':0,'rows':witness[137-2*j,r+3*j]}for r in[1,2]for j in range(9)if(137-2*j,r+3*j)in witness],'scope':'Positive literal witnesses validate data; completeness and negative membership rely on the separate canonical-generation proof.'}
print(json.dumps(report,sort_keys=True,indent=2))
