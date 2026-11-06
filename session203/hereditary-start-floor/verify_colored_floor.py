#!/usr/bin/env python3
"""Independent physical construction and rooted path DP. No packet imports."""
import itertools,json,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent
CORE=(102,0,10,83,135,30,182,47)
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def body_from_code(n,code):
 a=[0]*n
 for k,(i,j)in enumerate(itertools.combinations(range(n),2)):
  if code>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a
def build(body,color):
 # Independent tuple labels; do not inherit the producer's inflation indexing.
 labels=[('c',i)for i in range(8)if i!=5]+[('b',i)for i in range(len(body))]
 def arc(u,v):
  if u[0]==v[0]=='b':return bool(body[u[1]]>>v[1]&1)
  if u[0]==v[0]=='c':return bool(CORE[u[1]]>>v[1]&1)
  c,b=(u,v)if u[0]=='c'else(v,u)
  forward=bool(color>>b[1]&1)if c[1]==3 else bool(CORE[c[1]]>>5&1)
  return forward if u[0]=='c'else not forward
 full=[sum(1<<j for j,v in enumerate(labels)if i!=j and arc(u,v))for i,u in enumerate(labels)]
 sink=labels.index(('c',1));require(full[sink]==0,'old sink missing')
 keep=[i for i in range(len(labels))if i!=sink]
 return [sum(1<<j for j,y in enumerate(keep)if full[x]>>y&1)for x in keep],[labels[i]for i in keep]
def start_counts(a):
 n=len(a);N=1<<n;F=N-1;ans=[]
 for root in range(n):
  dp=[[0]*n for _ in range(N)];dp[1<<root][root]=1
  for S in range(N):
   for last in range(n):
    z=dp[S][last]
    if not z:continue
    todo=a[last]&~S
    while todo:
     b=todo&-todo;todo-=b;j=b.bit_length()-1;dp[S|b][j]+=z
  ans.append(sum(dp[F]))
 return ans
def strong(a):
 n=len(a);F=(1<<n)-1
 for root in range(n):
  seen=1<<root;pending=seen
  while pending:
   b=pending&-pending;pending-=b;v=b.bit_length()-1;fresh=a[v]&~seen;seen|=fresh;pending|=fresh
  if seen!=F:return False
 return True
def valid(a):
 for i in range(len(a)):
  require(not(a[i]>>i&1),'loop')
  for j in range(i):require(bool(a[i]>>j&1)!=bool(a[j]>>i&1),'not tournament')
K=[]
for code in range(8):
 body=body_from_code(3,code)
 for color in range(8):
  a,labels=build(body,color);valid(a);require(strong(a),'non-strong kernel');ss=start_counts(a)
  require(min(ss)>=27,'floor below27')
  K.append(dict(body_code=code,body_rows=body,color=color,rows=a,labels=labels,starts=ss,H=sum(ss),minimum=min(ss)))
require(len(K)==64 and min(k['minimum']for k in K)==27,'wrong universe/minimum')
reps=[]
for label,body,colors in [('transitive',[6,4,0],range(8)),('cycle',[2,4,1],[0,1,3,7])]:
 for color in colors:
  a,labels=build(body,color);ss=start_counts(a);reps.append(dict(type=label,color=color,minimum=min(ss),H=sum(ss),rows=a,starts=ss))
# Independent factorial cross-check on one sharp kernel.
a=K[0]['rows'];raw=[0]*len(a)
for w in itertools.permutations(range(len(a))):
 if all(a[u]>>v&1 for u,v in zip(w,w[1:])):raw[w[0]]+=1
require(raw==K[0]['starts'],'permutation control')
payload=dict(kernels=K,representatives=reps)
data=json.dumps(payload,sort_keys=True,indent=2)+'\n'
(P/'COLORED_KERNELS.json').write_text(data)
print(json.dumps(dict(status='PASS',colored_labeled_order_three_bodies=64,kernel_order=9,minimum_positive_start=27,kernels_strong=64,permutation_control_paths=sum(raw),kernel_sha256=hashlib.sha256(data.encode()).hexdigest(),representative_minima=[(r['type'],r['color'],r['minimum'])for r in reps],scope='Actual finite kernels; hereditary implication uses the separate insertion proof.'),sort_keys=True,indent=2))
