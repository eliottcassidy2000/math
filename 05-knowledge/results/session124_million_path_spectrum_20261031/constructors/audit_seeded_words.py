from pathlib import Path
DATA_DIR = Path(__file__).resolve().parent
"""Independent literal-permutation seed profiles and direct child subset-DP."""
from itertools import permutations, combinations, product
import json,random

def req(v,m):
 if not v:raise RuntimeError(m)
def adj(n,code):
 a=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=bool(code>>k&1);a[j][i]=not a[i][j]
 return a
def ham(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,1<<n):
  for v in range(n):
   if S>>v&1:
    R=S^(1<<v)
    if R:d[S][v]=sum(d[R][u] for u in range(n) if R>>u&1 and a[u][v])
 return sum(d[-1])
roots={}
for n in range(1,6):
 for code in range(1<<(n*(n-1)//2)):
  a=adj(n,code);paths=[];covers=[]
  for w in permutations(range(n)):
   good=[a[w[i]][w[i+1]] for i in range(n-1)]
   if all(good):paths.append(w)
   for j in range(n+1):
    if all(good[:max(j-1,0)]) and all(good[j:]):covers.append((w[:j],w[j:]))
  for cut in range(1<<n):
   t=sum(not(cut>>w[-1]&1) for w in paths)
   q=sum(bool(Q) and bool(cut>>Q[0]&1) for P,Q in covers)
   F=sum((not P or not(cut>>P[-1]&1)) and (not Q or bool(cut>>Q[0]&1)) for P,Q in covers)
   roots.setdefault((len(paths),F,q,t),(n,code,cut))
raw=json.load(open(DATA_DIR / 'seed_roots.json'))
req(set(roots)=={tuple(x['state'])for x in raw},'seed profiles differ')
print('independent_permutation_seed_states',len(roots))
rng=random.Random(124);cases=0
for state,witness in roots.items():
 n,code,cut=witness
 for rep in range(4):
  word=[rng.randrange(2) for _ in range(rng.randrange(5))]
  h,F,q,t=state
  for b in word:
   if b:F,q,t=F,2*q+h,0
   else:F,q,t=F+q+h-t,2*q,h
  N=n+len(word)+1;a=[[False]*N for _ in range(N)];old=adj(n,code)
  for i in range(n):
   for j in range(n):a[i][j]=old[i][j]
  for i in range(n+len(word)):
   for j in range(max(n,i+1),n+len(word)):a[i][j]=True
  bits=[bool(cut>>i&1)for i in range(n)]+word
  for i,b in enumerate(bits):a[-1][i]=bool(b);a[i][-1]=not b
  req(ham(a)==F,('literal_child',state,word,F,ham(a)))
  cases+=1
print('direct_child_dp_checks',cases)
print('PASS')
