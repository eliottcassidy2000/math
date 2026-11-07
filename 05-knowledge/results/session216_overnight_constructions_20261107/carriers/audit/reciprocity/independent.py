from pathlib import Path
import json,itertools,math
P=Path(__file__).parent

def tournaments(n):
 for bits in itertools.product((0,1),repeat=n*(n-1)//2):
  a=[0]*n;j=0
  for u in range(n):
   for v in range(u+1,n):
    if bits[j]:a[u]|=1<<v
    else:a[v]|=1<<u
    j+=1
  yield a

def partitions(vertices):
 if not vertices:yield [];return
 x=vertices[0]
 for rest in partitions(vertices[1:]):
  yield [[x]]+[b[:] for b in rest]
  for j in range(len(rest)):
   yield [([x]+b if j==k else b[:]) for k,b in enumerate(rest)]

def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def ev(p,x):return sum(c*x**i for i,c in enumerate(p))
def reflection(p):
 out=[0]*len(p)
 for i,c in enumerate(p):
  for j in range(i+1):out[j]+=c*math.comb(i,j)*(-1)**j
 return out
seen={};collisions=controls=ringchecks=0;hist={}
for n in range(1,6):
 parts=list(partitions(list(range(n))));N=1<<n;F=N-1
 for a in tournaments(n):
  H=[0]*N;S=[[0]*N for _ in range(n)];E=[[0]*N for _ in range(n)]
  for mask in range(N):
   vs=[j for j in range(n) if mask>>j&1]
   for p in itertools.permutations(vs):
    if all(a[u]>>v&1 for u,v in zip(p,p[1:])):
     H[mask]+=1
     if p:S[p[0]][mask]+=1;E[p[-1]][mask]+=1
  inv=[0]*N;inv[0]=1
  for mask in range(1,N):
   B=mask
   while B:inv[mask]-=H[B]*inv[mask^B];B=(B-1)&mask
   if inv[mask]!=(-1)**mask.bit_count()*H[mask]:raise RuntimeError('P reciprocity')
  for r in range(n):
   FS=[0]*n;GS=[0]*n
   for partition in parts:
    masks=[sum(1<<v for v in b) for b in partition];rootmask=next(m for m in masks if m>>r&1)
    multiplier=math.prod(H[m] for m in masks if m!=rootmask);j=len(partition)-1
    FS[j]+=S[r][rootmask]*multiplier;GS[j]+=E[r][rootmask]*multiplier
   p=[];q=[];fall=[1]
   for j in range(n):
    if j:fall=mul(fall,[-j,1])
    p=add(p,[FS[j]*z for z in fall]);q=add(q,[GS[j]*z for z in fall])
   if reflection(p)!=[(-1)**(n-1)*z for z in q]:raise RuntimeError(('root reciprocity',a,r,p,q))
   for mask in range(N):
    B=mask;actual=0
    while True:
     actual+=S[r][B]*inv[mask^B]
     if B==0:break
     B=(B-1)&mask
    if actual!=(-1)**(mask.bit_count()+1)*E[r][mask]:raise RuntimeError('root ring')
    ringchecks+=1
   s,e=S[r][F],E[r][F]
   R=sum(S[r][B]*H[F^B] for B in range(N));L=sum(E[r][B]*H[F^B] for B in range(N))
   C=[]
   for mask in range(N):
    B=mask;z=0
    while True:
     z+=H[B]*H[mask^B]
     if not B:break
     B=(B-1)&mask
    C.append(z)
   T=sum(S[r][B]*C[F^B] for B in range(N))
   if (ev(p,1),ev(p,2),ev(p,3),ev(q,1),ev(q,2))!=(s,R,T,e,L):raise RuntimeError('profile')
   key=(n%2,e,L,R,T)
   if key in seen:
    if (s-seen[key])%4:raise RuntimeError('same parity collision')
    collisions+=1
   else:seen[key]=s
   controls+=1
print(json.dumps({'rooted_polynomial_controls':controls,'squarefree_root_identity_controls':ringchecks,'same_parity_collisions':collisions,'method':'Unordered vertex-set partitions and direct path enumeration; integer falling-factorial basis; no SymPy'},indent=2))
