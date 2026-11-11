import itertools,json,math
from fractions import Fraction
from collections import Counter

def hp(A,S):
 return sum(all(A[p[i]][p[i+1]] for i in range(len(p)-1)) for p in itertools.permutations(S))
def q(A,t):
 n=len(A);h={s:hp(A,[i for i in range(n) if s>>i&1]) for s in range(1<<n)}
 d={0:1}
 for _ in range(t):
  z={}
  for S in range(1<<n):
   U=S;v=0
   while True:
    v+=h[U]*d.get(S^U,0)
    if U==0:break
    U=(U-1)&S
   z[S]=v
  d=z
 return d.get((1<<n)-1,0)
def base(n,bits):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  u,v=(j,i) if bits>>k&1 else (i,j);A[u][v]=1
 return A
def fan(A,a,k):
 n=len(A);B=[r+[0]*k for r in A]+[[0]*(n+k) for _ in range(k)]
 for i in range(k):
  b=n+i;B[b][a]=1
  for u in range(n):
   if u!=a:B[u][b]=1
  for j in range(i):B[n+j][b]=1
 return B
def induced(A,S):return [[A[i][j] for j in S] for i in S]
def signature(A,a,t):
 n=len(A);E=L=0
 for mask in range(1<<n):
  if not (mask>>a&1):continue
  R=[i for i in range(n) if mask>>i&1];rest=[i for i in range(n) if not(mask>>i&1)]
  w=t*q(induced(A,rest),t-1)
  E+=w*sum(all(A[p[i]][p[i+1]] for i in range(len(p)-1)) for p in itertools.permutations(R) if p[-1]==a)
  for u in range(1<<n):
   if u&~mask or u>>a&1:continue
   U=[i for i in R if u>>i&1];V=[i for i in R if not(u>>i&1)]
   L+=w*hp(A,U)*sum(all(A[p[i]][p[i+1]] for i in range(len(p)-1)) for p in itertools.permutations(V) if p[0]==a)
 D=t*q(induced(A,[i for i in range(n) if i!=a]),t)
 return q(A,t),E,L,D
fans=0
for n in range(1,4):
 for bits in range(1<<(n*(n-1)//2)):
  A=base(n,bits)
  for a in range(n):
   for t in range(1,4):
    vals=[q(fan(A,a,k),t) for k in range(4)]
    Q,E,L,D=signature(A,a,t)
    assert all(vals[k]==(L-D)*(t+1)**k+(Q-E-L+2*D)*t**k+(E-D)*(t-1)**k for k in range(4))
    assert vals[3]==3*t*vals[2]-(3*t*t-1)*vals[1]+(t**3-t)*vals[0]
    u=vals[1]-t*vals[0];v=vals[2]-2*t*vals[1]+t*t*vals[0]
    assert (u+v)%2==0 and (v-u)%2==0
    fans+=1

# Independently evaluate the moment loop-consistency formula, not orbit DP.
n=5;edges=[(i,j) for i in range(n) for j in range(i+2,n)];idx={e:k for k,e in enumerate(edges)}
paths=[]
for p in itertools.permutations(range(n)):
 req={};ok=True
 for u,v in zip(p,p[1:]):
  if abs(u-v)==1:
   if u>v:ok=False;break
  else:req[idx[tuple(sorted((u,v)))]]=int(u>v)
 if ok:paths.append(req)
result={}
for rep in range(4):
 ans=Fraction(0)
 for a in paths:
  for b in paths:
   if any(k in b and b[k]!=v for k,v in a.items()):continue
   req=a|b;adj=[[] for _ in range(n)]
   for k,v in req.items():
    u,w=edges[k];r=v^((rep>>k)&1);adj[u].append((w,r));adj[w].append((u,r))
   signs={};rank=0;ok=True
   for r in range(n):
    if r in signs:continue
    signs[r]=0;todo=[r];size=0
    while todo:
     u=todo.pop();size+=1
     for v,s in adj[u]:
      if v in signs:
       if signs[v]!=(signs[u]^s):ok=False
      else:signs[v]=signs[u]^s;todo.append(v)
    rank+=size-1
   if ok:ans+=Fraction(1,1<<rank)
 result[rep]=str(ans)
assert result=={0:'221/2',1:'227/2',2:'217/2',3:'109'}
births=0
for n in (4,5):
 def ed(m):return [(i,j) for i in range(m) for j in range(i+2,m)]
 def cuts(m):
  return {sum((((x>>i)^(x>>j))&1)<<k for k,(i,j) in enumerate(ed(m))) for x in range(1<<m)}
 C=cuts(n);D=cuts(n+1);old=ed(n);new=ed(n+1)
 groups={}
 for t in range(1<<len(old)):groups.setdefault(min(t^c for c in C),[]).append(t)
 for rep,group in groups.items():
  for a in range(n):
   counts=Counter()
   for t in group:
    A=[[0]*n for _ in range(n)]
    for i in range(n-1):A[i][i+1]=1
    for k,(i,j) in enumerate(old):
     u,v=(j,i) if t>>k&1 else (i,j);A[u][v]=1
    B=fan(A,a,1);order=list(range(a))+[n]+list(range(a,n))
    assert all(B[u][v] for u,v in zip(order,order[1:]))
    b=sum((1-B[order[i]][order[j]])<<k for k,(i,j) in enumerate(new))
    counts[min(b^d for d in D)]+=1
   assert len(counts)==1<<(n-2) and set(counts.values())=={2}
   births+=1
compiler_cases=0
for n in range(1,6):
 old=ed(n);new=ed(n+1);indices={edge:k for k,edge in enumerate(new)}
 for a in range(n):
  order=list(range(a))+[n]+list(range(a,n));pos={v:i for i,v in enumerate(order)}
  def coordinate(u,v):return indices[tuple(sorted((pos[u],pos[v])))]
  constants=0
  for u in range(n):
   if u==a or (a>0 and u==a-1):continue
   if u>a:constants|=1<<coordinate(n,u)
  for z in range(1<<len(new)):
   y=[0]*(n+1)
   for u in range(n):
    if u==a:continue
    key=coordinate(a-1,a) if a>0 and u==a-1 else coordinate(n,u)
    y[u]=((z^constants)>>key)&1
   t=sum((((z>>coordinate(i,j))&1)^y[i]^y[j])<<k for k,(i,j) in enumerate(old))
   born=constants|sum(((t>>k)&1)<<coordinate(i,j) for k,(i,j) in enumerate(old))
   out=born^sum((y[order[i]]^y[order[j]])<<k for k,(i,j) in enumerate(new))
   assert out==z
   compiler_cases+=1
print(json.dumps({'paired_root_inverse_cases':compiler_cases,'birth_sector_cases':births,'fan_recurrence_cases' :fans,'fan_universe':'all labelled n<=3, every root, t=1..3, lengths0..3','moment_order5':result,'status':'PASS'},indent=2))
