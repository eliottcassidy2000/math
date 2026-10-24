#!/usr/bin/env python3
import argparse,itertools,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args();checks=0
def require(x):
 global checks
 checks+=1
 if not x:raise ValueError('control '+str(checks))
def tournament(n,b):
 T=[[0]*n for _ in range(n)]
 for i in range(n-1):T[i][i+1]=1
 free=[(i,j) for i in range(n) for j in range(i+2,n)]
 for k,(i,j) in enumerate(free):T[i][j]=(b>>k)&1;T[j][i]=1-T[i][j]
 return T
def complement(T,path):
 n=len(T);kept=set(zip(path,path[1:]));R=[[0]*n for _ in range(n)]
 require(sorted(path)==list(range(n)))
 require(all(T[i][j] for i,j in kept))
 for i in range(n):
  for j in range(n):
   if i!=j:R[i][j]=T[i][j] if (i,j) in kept or (j,i) in kept else T[j][i]
 return R
def suspend(T):
 n=len(T);S=[r[:]+[0,1] for r in T]+[[1]*n+[0,1],[0]*(n+2)]
 return S
def triangles(T):
 return sum((T[i][j] and T[j][k] and T[k][i]) or (T[j][i] and T[k][j] and T[i][k]) for i,j,k in itertools.combinations(range(len(T)),3))
cases=[]
for bits in (2,3):
 T=tournament(5,bits);K=complement(T,list(range(5)));H=complement(suspend(T),[5,0,1,2,3,4,6]);d=list(map(sum,T));scores=sorted(map(sum,H))
 require(min(d)>0 and max(d)<4)
 require([row[:5] for row in H[:5]]==K)
 require([j for j in range(7) if H[5][j]]==[0])
 require([i for i in range(7) if H[i][6]]==[4])
 require(scores==sorted([1,5]+[5-x for x in d]))
 require([i for i,row in enumerate(H) if sum(row)==1]==[5])
 require([i for i,row in enumerate(H) if sum(row)==5]==[6])
 cases.append({'bits':bits,'T':T,'K':K,'H':H,'T_degrees':d,'H_scores':scores,'H_degree_square_sum':sum(x*x for x in scores)})
f=[1,2,3,0,4];K,J=[c['K'] for c in cases]
require(all(K[i][j]==J[f[i]][f[j]] for i in range(5) for j in range(5)))
require(cases[0]['H_scores']!=cases[1]['H_scores'])
require(cases[0]['H_scores']!=sorted(6-x for x in cases[1]['H_scores']))
require([c['H_degree_square_sum'] for c in cases]==[73,75])
# Independent exhaustive canonicalization also includes converse and owner permutations.
def canonical(T):
 n=len(T)
 return min(tuple(T[p[i]][p[j]] for i in range(n) for j in range(i+1,n)) for p in itertools.permutations(range(n)))
require(canonical(K)==canonical(J))
H0,H1=[c['H'] for c in cases]
require(canonical(H0)!=canonical(H1))
require(canonical(H0)!=canonical(list(map(list,zip(*H1)))))
# Validate suspension identities for all forward-path orders n<=6 (1+2+8+64+1024).
universe=0
for n in range(2,7):
 for b in range(1<<((n-1)*(n-2)//2)):
  T=tournament(n,b);K=complement(T,list(range(n)));H=complement(suspend(T),[n]+list(range(n))+[n+1]);universe+=1
  require([r[:n] for r in H[:n]]==K)
  require(list(map(sum,H))==[n-sum(r) for r in T]+[1,n])
  require([j for j in range(n+2) if H[n][j]]==[0])
  require([i for i in range(n+2) if H[i][n+1]]==[n-1])
  d=list(map(sum,T));ct,ck,ch=map(triangles,(T,K,H))
  require(ch==ct+n)
  require(ck==ct+d[0]-d[-1]-1)
  require(ch-ck==n+1-d[0]+d[-1])
result={'status':'PASS','checks':checks,'cases':cases,'old_target_isomorphism':f,'forward_path_tournaments_checked':universe,'scope':'Explicit order5 to order7 class-descent hostile and suspension identity controls, not harmonic transfer or full metagraph enumeration.'}
s=json.dumps(result,indent=2,sort_keys=True)+'\n'
if a.output:a.output.write_text(s)
else:print(s,end='')
