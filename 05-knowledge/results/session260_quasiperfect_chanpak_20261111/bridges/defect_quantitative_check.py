from itertools import combinations,permutations
from pathlib import Path
from random import Random
import argparse,json

def require(q,m):
 if not q:raise RuntimeError(m)
def all_h(out):
 n=len(out);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for S in range(1,1<<n):
  for i,x in enumerate(dp[S]):
   if x:
    todo=out[i]&~S
    while todo:
     bit=todo&-todo;todo-=bit;dp[S|bit][bit.bit_length()-1]+=x
 return [sum(row) if i else 1 for i,row in enumerate(dp)]
def probe(out,kernel=False):
 n=len(out);full=(1<<n)-1;hs=all_h(out);cases=[]
 if kernel:
  paths=[p for p in permutations(range(n)) if all(out[p[i]]>>p[i+1]&1 for i in range(n-1))]
  weights={}
  for u in range(n):
   for v in range(n):
    if out[u]>>v&1:
     L=sum(any(p[i:i+2]==(u,v) for i in range(n-1)) for p in paths)
     flip=out.copy();flip[u]^=1<<v;flip[v]^=1<<u
     R=all_h(flip)[full]-hs[full]+L
     weights[u,v]=L+R
     require(L+R>=2 and (L+R)%2==0,('weight',out,u,v))
 for z in range(n):
  V=full^(1<<z);A=out[z];B=V^A;adj=[0]*n
  for a in range(n):
   if A>>a&1:
    adj[a]=out[a]&B
    for b in range(n):
     if adj[a]>>b&1:adj[b]|=1<<a
  D=0;S=V
  while True:
   D+=all(S>>i&1 or adj[i]&S for i in range(n) if i!=z)
   if not S:break
   S=(S-1)&V
  t=sum((out[a]>>b&1 and out[b]>>c&1 and out[c]>>a&1) or (out[a]>>c&1 and out[c]>>b&1 and out[b]>>a&1) for a,b,c in combinations([i for i in range(n) if i!=z],3))
  if hs[full]-D<hs[V]-1:return {'failure':'old_H_bound','out':out,'z':z,'H':hs[full],'d':D,'oldH':hs[V]}
  require(hs[full]-D>=2*t,('triangle',out,z))
  # All-arc identity for births from this U, checked separately below.
  cases.append((hs[full]-D,hs[V]-1))
 if kernel:
  for A in range(1<<n):
   child=out.copy()+[A]
   for i in range(n):
    if not A>>i&1:child[i]|=1<<n
   inc=all_h(child)[(1<<(n+1))-1]-hs[full]
   expected=sum(w for (u,v),w in weights.items() if A>>u&1 and not A>>v&1)
   require(inc==expected,('kernel',out,A,inc,expected))
 return cases
p=argparse.ArgumentParser();p.add_argument('--max-order',type=int,default=6);args=p.parse_args()
rows=[];rootcases=0;kernelbases=0
for n in range(1,args.max_order+1):
 pairs=list(combinations(range(n),2));graphs=0
 for mask in range(1<<len(pairs)):
  out=[0]*n
  for bit,(i,j) in enumerate(pairs):
   if mask>>bit&1:out[i]|=1<<j
   else:out[j]|=1<<i
  ans=probe(out,kernel=n<=4)
  require(not isinstance(ans,dict),ans);rootcases+=n;graphs+=1;kernelbases+=n<=4
 rows.append({'order':n,'graphs':graphs,'root_cases':graphs*n});print(rows[-1],flush=True)
rng=Random(260);randomcases=0
for n in range(7,11):
 for _ in range(200):
  out=[0]*n
  for i,j in combinations(range(n),2):
   if rng.randrange(2):out[i]|=1<<j
   else:out[j]|=1<<i
  ans=probe(out);require(not isinstance(ans,dict),ans);randomcases+=n
res={'status':'FINITE-EXACT for conjectural inequalities; identity proved in note','exhaustive_rows':rows,'root_cases':rootcases,'kernel_bases':kernelbases,'all_birth_cuts_tested_on_kernel_bases':sum(r['graphs']*2**r['order'] for r in rows if r['order']<=4),'seeded_random_roots_orders_7_to_10':randomcases,'counterexamples':[]}
Path(__file__).with_name('defect_quantitative_receipt.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
