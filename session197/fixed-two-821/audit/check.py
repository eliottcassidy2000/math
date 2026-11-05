from pathlib import Path
from itertools import combinations
import json,hashlib
P=Path(__file__).resolve().parent
def need(x,s):
 if not x:raise RuntimeError(s)
def walk(a):
 n=len(a);ans=0
 for S in range(1,1<<n):
  v=[int(S>>i&1)for i in range(n)]
  for _ in range(n-1):
   w=[0]*n
   for u,x in enumerate(v):
    if not x:continue
    U=a[u]&S
    while U:
     b=U&-U;U-=b;w[b.bit_length()-1]+=x
   v=w
  ans+=(-1)**(n-S.bit_count())*sum(v)
 return ans
def cycles(a):
 n=len(a);full=(1<<n)-1;ans=0
 def rec(S,v):
  nonlocal ans
  if S==full:ans+=bool(a[v]&1);return
  U=a[v]&~S
  while U:
   b=U&-U;U-=b;rec(S|b,b.bit_length()-1)
 rec(1,0);return ans
def fvs(a):
 n=len(a);largest=0;delete=None
 for S in range(1<<n):
  V=[v for v in range(n)if S>>v&1]
  if sorted((a[v]&S).bit_count()for v in V)==list(range(len(V))) and len(V)>largest:
   largest=len(V);delete=[v for v in range(n)if not S>>v&1]
 return n-largest,delete
def physical(word):
 a=[2,0]
 for mask in word:
  n=len(a);a.append(0)
  for v in range(n):
   if v<2 and mask>>v&1:a[n]|=1<<v
   else:a[v]|=1<<n
 return a
# Direct subset profiles only for the six-vertex seed, not used in the census filter.
def seedprofile(a,r):
 n=len(a);N=1<<n;start=[[0]*n for _ in range(N)];end=[[0]*n for _ in range(N)];h=[1]+[0]*(N-1)
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   U=S^(1<<v)
   if U==0:start[S][v]=end[S][v]=1
   else:
    for w in range(n):
     if U>>w&1:
      if a[v]>>w&1:start[S][v]+=start[U][w]
      else:end[S][v]+=end[U][w]
  h[S]=sum(start[S])
 F=N-1;return [h[F],end[F][r],sum(start[S][r]*h[F^S]for S in range(1,N)if S>>r&1),start[F][r],h[F^(1<<r)]]
records=[];hist={};subsets=0
for line in (P/'lowcycle-821-order13.csv').read_text().splitlines():
 n,c,*a=map(int,line.split(','));need(walk(a)==821 and cycles(a)==c,'independent literal core counts');f,d=fvs(a);need(f>=3,'feedback-two hostile');records.append({'n':n,'C':c,'rows':a,'feedback':f,'deletion':d});hist[str((n,c,f))]=hist.get(str((n,c,f)),0)+1;subsets+=1<<n
need(len(records)==15,'complete hits');need(821**3<5**13,'strong order13 cutoff')
seed=physical([0,0,3,1]);need(seedprofile(seed,0)==[37,5,49,4,5],'exact seed profile')
examples=[]
for k in range(1,7):
 a=physical([0,0,3,1]+[0]*k+[1]);target=37+49*(2**k)
 need(walk(a)==target and cycles(a)==4,'actual binary family');examples.append({'k':k,'H':target,'C':4,'rows':a})
need(walk(physical([0,0,3,1,1]))==81,'k0 boundary is not the same scalar formula')
a=physical([0,0,3,1]+[0]*4+[1]);need(a==[974,2028,2040,2032,2019,1985,1920,1792,1536,1024,1],'literal 821 word')
f,d=fvs(a);need(f==2,'exact positive feedback2')
out={'status':'PASS','strict_cutoff':[821**3,5**13],'complete_lowcycle_bank_classes':80746,'H821_cores':records,'histogram':hist,'core_induced_subset_checks':subsets,'positive_seed':{'rows':seed,'h_e_R_s_d':[37,5,49,4,5]},'positive_H821_feedback':f,'positive_H821_deletion':d,'positive_examples':examples,'k0_H':81}
print(json.dumps(out,sort_keys=True,indent=2))
